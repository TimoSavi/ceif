# ceif Server Automation & Cron Recipes

`ceif` is designed specifically for automated server environments, cron jobs, and Unix shell pipelines. Unlike Python or JVM-based anomaly detection engines, `ceif`:
- Has **zero runtime dependencies** (pure C binary).
- Runs with **constant, bounded memory** via reservoir sampling (will not trigger OOM kills on multi-gigabyte log spikes).
- Natively isolates **multi-tenant categories** (e.g., thousands of separate hosts, ports, or metrics) within a single process.
- Emits standard Unix exit codes (`0` = clean, `2` = anomaly detected).

---

## 1. Quick Reference: The 6 Flags You Need for Cron

| Flag | Purpose | Example |
|:---|:---|:---|
| **`-l FILE`** | Train/learn model from data (`-` for stdin) | `ceif -l baseline.csv -w model.f` |
| **`-r FILE`** | Load pre-trained model | `ceif -r model.f -a current.csv` |
| **`-a FILE`** | Score data for anomalies (`-` for stdin) | `tail -f log.csv | ceif -r model.f -a -` |
| **`-w FILE`** | Save trained model to disk | `ceif -l data.csv -w model.f` |
| **`-z FILE`** | In-place update: loads, trains, saves back | `ceif -z model.f -l latest.csv` |
| **`-O THRESH`** | Outlier trigger threshold (`0.65`, `0.80s`, `99%`) | `ceif -r model.f -a now.csv -O 0.70s` |

---

## 2. Threshold Strategies

`ceif` supports three ways to define what constitutes an "anomaly":

1. **Scaled Score (`-O 0.80s`) — Recommended for most cron jobs**
   Raw scores are normalized $0.0 \dots 1.0$ against the min/max scores seen in that specific category's forest. Setting `-O 0.85s` alerts only when an event enters the top 15% severity band of that forest.
2. **Percentile Cutoff (`-O 99%` or `-O 80%`)**
   - High percentiles (`-O 99%` / `99.5%`): Alerts only on rare, extreme point anomalies.
   - Moderate percentiles (`-O 80%`): Used for **macro population drift detection**, alerting when the fraction of new samples entering the top quintile significantly exceeds baseline expectation (e.g. $\ge 40\%$).
3. **Fixed Score (`-O 0.65`)**
   Direct unscaled isolation score. $0.5$ is typical baseline density; $> 0.65$ indicates significant structural isolation.

---

## 3. Recipe 1: Scheduled Anomaly Scan (Every 5 Minutes)

Scans the last 5 minutes of server metrics or access logs against a pre-trained baseline model.

```bash
#!/bin/bash
# /usr/local/bin/check-server-anomalies.sh
set -euo pipefail

MODEL="/var/lib/ceif/system-metrics.f"
LOG="/var/log/sysstat/metrics-last-5m.csv"
ALERT_EMAIL="noc-alerts@example.com"

# Exit code 2 indicates anomalies were found
if anomalies=$(ceif -r "$MODEL" -a "$LOG" -H -O 0.80s -p "%t ALERT [%C] score=%s metrics=%v"); then
    # Exit 0: Everything normal
    exit 0
elif [ $? -eq 2 ]; then
    # Exit 2: Anomalies detected
    echo "$anomalies" | mail -s "[CEIF ANOMALY] Server Metrics Alert" "$ALERT_EMAIL"
    exit 2
else
    # Exit 1: Fatal execution/parse error
    echo "ceif execution failed" >&2
    exit 1
fi
```

**Crontab entry (`/etc/cron.d/ceif-monitor`):**
```cron
*/5 * * * * root /usr/local/bin/check-server-anomalies.sh >/dev/null 2>&1
```

---

## 4. Recipe 2: Rolling / Self-Updating Model (`-z`)

In dynamic server environments, normal traffic patterns evolve over time. Using `-z` (in-place forest update) automatically ingests the latest batch of samples into the model using reservoir sampling:

```bash
#!/bin/bash
# /usr/local/bin/update-ceif-model.sh
# Ingests the last hour of telemetry into the rolling model without re-reading history.
set -euo pipefail

MODEL="/var/lib/ceif/rolling-traffic.f"
HOURLY_DATA="/var/log/traffic/last-hour.csv"

# Automatically locks, updates reservoir samples, retrains trees, and writes back
flock -x "$MODEL.lock" ceif -z "$MODEL" -l "$HOURLY_DATA" -H -s 256 -t 100
```

**Crontab entry:**
```cron
0 * * * * root /usr/local/bin/update-ceif-model.sh
```

---

## 5. Recipe 3: Multi-Tenant / Per-Device Tracking (`-C`)

If you have metrics from 500 different servers or 2,000 network interfaces in a single CSV file, **do not write bash loops**. `ceif` partitions categories automatically with `-C`:

```text
# Input CSV format: hostname,cpu,memory,iowait,net_err
server-01,12.4,45.1,0.2,0
server-02,89.1,95.0,8.5,12
...
```

**Train individual baseline models per server in one command:**
```bash
ceif -l cluster-baseline.csv -H -C 1 -w /var/lib/ceif/cluster.f
```

**Score the current cluster status:**
```bash
ceif -r /var/lib/ceif/cluster.f -a cluster-now.csv -H -C 1 -O 0.75s      -p "%t HOST=%C SCORE=%s RAW=%v" >> /var/log/ceif/anomalies.log
```
*Each host is scored strictly against its own baseline cluster without cross-contamination.*

---

## 6. Recipe 4: Pruning Retired Hosts / Stale Categories (`-D`)

When servers or containers are decommissioned, their categories remain in the model file unless pruned. The `-D` flag deletes any category not updated within a given time interval:

```bash
# Retrain or refresh model, automatically dropping categories inactive for > 30 days
ceif -z /var/lib/ceif/cluster.f -l recent.csv -D 30D
```
*(Accepts suffix `s` for seconds, `m` for minutes, `H` for hours, `D` for days, `M` for months, `Y` for years).*

---

## 7. Recipe 5: Real-Time Stream Piping (Pipes & `tail -F`)

Monitor incoming log streams in real-time and pipe directly into `logger` (syslog) or an alerting webhook:

```bash
tail -F /var/log/suricata/eve-network.csv |   ceif -r /var/lib/ceif/network.f -a - -O 0.85s -p "%t [SURICATA-ANOMALY] score=%s flow=%v" |   logger -t ceif -p local0.warning
```

---

## 8. Recipe 6: Population Drift & Distribution Shift Detection (`-O 80% -v`)

Standard point-anomaly detection searches for rare outlier events ($< 1\%$ of traffic). In contrast, **population drift detection** evaluates whether the macro-level distribution of telemetry has shifted across an entire population, category, or time window.

### Concept

When evaluating a batch of samples against a baseline model with an 80th-percentile threshold (`-O 80%`), in normal steady-state operations exactly **20%** of samples will exceed the threshold by definition.

If a new batch arrives and suddenly **$\ge 40\%$** of samples score as outliers, the underlying statistical distribution has drifted (e.g., due to software rollouts, changed network routing, customer behavioral shifts, or system degradation):

$$\text{Drift Ratio} = \frac{\%h \text{ (outlier count)}}{\%o \text{ (analyzed rows)}} \ge 0.40 \implies \mathbf{Alert: Population Drift}$$

Unlike univariate statistical tests, `ceif` detects multivariate covariance shifts across all dimensions simultaneously in milliseconds with constant memory.

### Implementation

Combine `-p ""` (to suppress individual sample output) with `-v "%C %h %o %S"` to emit concise aggregate statistics per category:
- `%C`: Category / service name
- `%h`: Count of rows with score $\ge$ threshold
- `%o`: Count of analyzed rows
- `%S`: Mean anomaly score of the batch

```bash
#!/bin/bash
# /usr/local/bin/check-population-drift.sh
set -euo pipefail

MODEL="/var/lib/ceif/baseline.f"
BATCH="/var/log/traffic/last-hour.csv"
DRIFT_THRESHOLD=0.40   # Alert if >= 40% of samples fall into top quintile

drift_alerts=$(ceif -r "$MODEL" -a "$BATCH" -O 80% -p "" -v "%C %h %o %S" | \
awk -v thresh="$DRIFT_THRESHOLD" '
{
    category  = ($1 == "" ? "default" : $1);
    outliers  = $2;
    total     = $3;
    avg_score = $4;
    ratio     = (total > 0 ? (outliers / total) : 0);

    if (ratio >= thresh) {
        printf("DRIFT ALERT [%s]: %.1f%% outliers exceed 80th-percentile (outliers=%d, total=%d, avg_score=%s)\n",
               category, ratio * 100, outliers, total, avg_score);
        drift_found = 1;
    }
}
END {
    if (drift_found) exit 2;
}')

if [ $? -eq 2 ]; then
    echo "$drift_alerts" | mail -s "[CEIF DRIFT] Population Distribution Shift Detected" noc-alerts@example.com
    exit 2
fi
```

---

## 9. Format Directives Quick Reference (`-p STRING` & `-v STRING`)

Customize individual alert outputs (`-p`) or aggregate summary statistics (`-v`):

| Directive | Available In | Output Value | Example |
|:---|:---|:---|:---|
| **`%s`** | `-p`, `-v` | Anomaly score | `0.781204` |
| **`%S`** | `-p`, `-v` | Average anomaly score for analyzed data | `0.654120` |
| **`%C`** | `-p`, `-v` | Category name (from `-C`) | `web-prod-03` |
| **`%h`** | `-p`, `-v` | Count of rows exceeding outlier threshold | `400` |
| **`%o`** | `-p`, `-v` | Count of analyzed rows for forest | `1000` |
| **`%n`** | `-p`, `-v` | Total rows in forest | `1000` |
| **`%l`** | `-p` | Row label (from `-L`) | `request_id_9981` |
| **`%v`** | `-p` | Raw input dimension values | `45.2, 120.4, 0.02` |
| **`%t`** | `-p`, `-v` | Current timestamp (formatted via locale) | `2026-09-08 14:55:00` |
| **`%x`** | `-p`, `-v` | Hex RGB color code (for dashboards) | `#ff2200` |
| **`%m`** | `-p` | Multi-dimension metric summary (with `-j`) | `cpu=98% mem=92%` |

---

## 10. Standard Crontab Best Practices for `ceif`

1. **Always Set `PATH`:** Ensure cron's minimal environment finds your binaries (`PATH=/usr/local/bin:/usr/bin:/bin`).
2. **Use File Locks (`flock`):** When combining `-z` (rolling updates) with `-a` (scoring), use `flock` to prevent concurrent writes to the `.f` model file.
3. **Capture Exit Code `2`:** In monitoring scripts, `exit 2` means anomalies were detected, whereas `exit 1` means syntax or file I/O error. Check `$?` to distinguish between real security/metric alerts and misconfiguration.

---

## 11. Centralized Threshold & Configuration Hierarchy (`-g`)

In complex production environments running multiple cron jobs across different datasets or microservices, avoid hardcoding thresholds across multiple crontabs or scripts. Use a centralized custom config file (`-g`):

```bash
# Both cron jobs evaluate distinct models but apply the centralized threshold defined in cron-rules.rc:
ceif -r /var/lib/ceif/auth.f -g /etc/ceif/cron-rules.rc -a /var/log/auth.csv
ceif -r /var/lib/ceif/db.f   -g /etc/ceif/cron-rules.rc -a /var/log/db.csv
```

### Precedence Hierarchy

Settings are resolved in the following priority order:
1. **`~/.ceifrc`:** User baseline defaults.
2. **Forest Model (`-r` / `-z`):** Settings saved at training time in the model file.
3. **Custom Config (`-g` / `--rc-file`):** Overrides model-stored settings for centralized management.
4. **Direct CLI Options (`-O`, `-t`, `-s`, etc.):** Highest priority; overrides all configuration files.


