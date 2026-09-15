## CEIF Documentation Index

Welcome to the `ceif` documentation suite. Whether you are building from source, calibrating decision boundaries for novel datasets, or deploying automated monitoring pipelines in cron, the following guides cover each aspect:

### Guides

- **[User Manual](manual.md)**  
  The complete command-line reference, covering workflow flags (`-l`, `-a`, `-c`, `-w`, `-z`), format directives (`-p`, `-v`, `-j`), configuration precedence hierarchy (`~/.ceifrc`, `-g`), rc-file parameters, expression syntax (`-Q`), and Unix exit codes.

- **[Building from Source](building.md)**  
  Prerequisites, compiler and build dependencies (including `json-c` installation across RHEL/Rocky, Debian/Ubuntu, Alpine, Arch, and macOS), autotools configuration, and installation.

- **[Tweaking & Testing Guide](tweaking.md)**  
  Practical guide to model calibration and anomaly score maps. Covers score scaling (`-O 0.5s`), percentile thresholds, novelty detection / pure-inlier modeling (`-O 100%`), handling complex topologies (`NEAREST 1`), continuous model updates, stream sampling (`ANALYZE_SAMPLING`), and cluster-based attribute contributions (`%e`).

- **[Server Automation & Cron Recipes](cron-recipes.md)**  
  Production-tested recipes for server automation, scheduled scans, rolling self-updating models (`-z`), automatic reservoir ceiling (`EXTRA_ROWS_FACTOR 3`), multi-tenant category isolation (`-C`), real-time stream piping (`tail -F`), and population distribution shift / drift detection (`-O 80% -v`).
