## Categorized Extended Isolation Forest Tool (`ceif`)

`ceif` is a high-performance command-line utility for anomaly detection and categorization based on the Extended Isolation Forest (EIF) algorithm by [Hariri et al.](https://arxiv.org/abs/1811.02141). 

It is designed for automated production environments, cron jobs, and Unix shell pipelines where diverse datasets must be monitored efficiently without the overhead or dependency footprint of large machine learning frameworks.

### Key Features

* **Multi-Tenant Categorization (`-C`)**: Input data fields can serve as category keys. `ceif` partitions categories automatically, training and evaluating independent forests per category within a single process.
* **Label Tracking (`-L`)**: Designate non-numeric identifiers (e.g., timestamps, UUIDs, hostnames) as labels to identify anomalous records easily in output streams.
* **JSON Model Persistence (`-w`, `-r`, `-z`)**: Models are saved in standard **JSON format by default** (with CSV fallback when compiled without JSON-C), preserving calibration parameters, expressions, and sample reservoirs.
* **Bounded-Memory Reservoir Sampling**: Processes multi-gigabyte files with a constant, configurable sample memory footprint.
* **Continuous In-Place Updates (`-z`) with Automatic Reservoir Ceiling**: Models continuously ingest new batches of samples using reservoir sampling. An automatic ceiling (`EXTRA_ROWS_FACTOR 3`) caps historical extra rows, guaranteeing a minimum 25% acceptance probability for incoming data to prevent long-running models from freezing.
* **Dynamic Input Expressions (`-Q`)**: Transform input dimensions on the fly using arithmetic and mathematical expressions (via `tinyexpr`).
* **Unix Pipeline & Exit Codes**: Emits standard exit codes (`0` = clean, `2` = anomalies detected) for native integration with shell scripts, `cron`, `tail -F`, and `logger`.

---

### Documentation

Comprehensive documentation is available in the [`docs/`](docs/README.md) directory:

* **[User Manual](docs/manual.md)**: Full command-line reference, options, formatting directives, configuration precedence hierarchy (`~/.ceifrc`, `-g`), rc-file parameters, and usage examples.
* **[Building from Source](docs/building.md)**: Prerequisites, build dependencies (`json-c`), package manager installation across distributions, and compilation.
* **[Tweaking & Testing Guide](docs/tweaking.md)**: Calibration guide, heatmaps, score scaling (`-O 0.5s`), percentile thresholds, novelty detection (`-O 100%`), handling complex topologies (`NEAREST 1`), and attribute contribution analysis (`%e`).
* **[Server Automation & Cron Recipes](docs/cron-recipes.md)**: Production cron scripts, scheduled scans, rolling self-updating models (`-z`), multi-tenant tracking (`-C`), stream piping, and population drift detection (`-O 80% -v`).

---

### Algorithm Modifications

`ceif` incorporates key enhancements over the canonical Extended Isolation Forest algorithm:

#### 1. Pairwise Split Point Selection ($p$)
The canonical EIF selects intercept points $p$ uniformly from a rectangular bounding box. When data is uniformly distributed or non-convex, subspaces extending to infinity frequently contain sample points, driving anomaly scores toward ~0.5 across the entire space and causing artifact inliers.

To resolve this, `ceif` uses **Pairwise Split Point Selection**:
1. At each node, two sample points $x_1$ and $x_2$ are randomly selected from the node's sample subset.
2. The split intercept point $p$ is chosen along the line segment between them:
   $$p = x_1 + u \cdot (x_2 - x_1), \quad u \sim \mathcal{U}(-\text{margin}, 1 + \text{margin})$$
3. **Depth-Tapered Margin**: Near the root of the tree, the margin is wider to allow smooth decision boundaries into exterior space. As tree depth increases toward the leaves, the margin contracts toward 0.
4. **Adaptive Density Calibration**: The margin dynamically scales based on the distance between $x_1$ and $x_2$ relative to the average consecutive sample distance in the node. When $x_1$ and $x_2$ span across an empty void or separate clusters, the margin strictly contracts to $[0, 1]$, ensuring the split cleanly bisects the void without bridging artifacts.

#### 2. Nearest Training Point Distance in Leaf Nodes
At leaf nodes, `ceif` evaluates the relative distance between the analyzed point and the nearest training data point in the node (`NEAREST 1`, default). The distance is normalized against the average sample distance within the tree. If the distance is larger than average, the anomaly score is incremented; if smaller, it is reduced. This allows `ceif` to cleanly isolate interior voids, rings, and complex manifold shapes.

---

### Acknowledgements

* Extended Isolation Forest algorithm by Sahand Hariri, Matias Carrasco Kind, and Robert J. Brunner ([arXiv:1811.02141](https://arxiv.org/abs/1811.02141)).
* Expression parsing powered by [tinyexpr](https://github.com/codeplea/tinyexpr) by Lewis Van Winkle.
