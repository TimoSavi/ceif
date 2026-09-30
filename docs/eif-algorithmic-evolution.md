# Beyond EIF: Fundamental Algorithmic Evolution for Isolation Forests

**Author / Maintainer:** Timo Savinen (AI-assisted)

---

## 1. Executive Summary & Context

The Extended Isolation Forest (EIF, Hariri et al., 2019) was proposed to solve the axis-aligned bias of Liu et al.'s original Isolation Forest (2008) by using random linear hyperplanes. However, real-world deployment reveals that both algorithms suffer from severe theoretical limitations:

1. **The Infinite Hyperplane Paradox**: Flat hyperplanes extend to $\pm\infty$, arbitrarily slicing empty space millions of units away and causing radial starburst artifacts.
2. **The Density Bias**: Discrete path length (tree depth) conflates sample count with spatial density. Dense clusters require many cuts to isolate, while diffuse outliers can mimic cluster depths.
3. **Topological Blindness**: Infinite linear cuts struggle to resolve non-convex geometries (cavities, concentric donuts, winding canals) without bridging empty voids.
4. **Dimension Magnitude Distortion**: Pure isotropic Gaussian normals assume unit variance across all features, failing when attributes have unequal physical units.

In `ceif`, practical engineering tackled each of these flaws through dedicated mechanisms:
- `auto_weigth` (dimension scaling)
- Pairwise sample interpolation for $p$ with quadratic depth decay margin
- `NEAREST` distance scaling at leaves
- Exterior Euclidean distance decay in deep space
- Structural Zero Kelvin deepest-leaf minimum score calibration

While highly effective, these mechanisms were initially perceived as "guardrails" around the core EIF tree. The `geif` project was built to test whether these properties could emerge **naturally from first principles** in a next-generation geometric forest. 

Crucially, **empirical testing in `geif` has now validated or invalidated several key theoretical hypotheses**. Notably, while outer stadium decay and continuous metric depth proved highly successful, **data-driven Voronoi bisector splitting did not yield better results** than isotropic Gaussian normals with sample-anchored projections. 

This paper synthesizes these empirical findings, explains why pure Voronoi bisectors degrade ensemble performance, documents the current state-of-the-art architecture, and outlines the next set of fertile geometric theories to test.

---

## 2. Deconstructing the Flaws of Classic EIF

```
Classic EIF Model:
      Random Gaussian Normal n ~ N(0, I)
                 +
      Uniform Intercept p ~ U[min, max]
                 ↓
      Infinite Hyperplane: (x - p) · n = 0
```

| Classic EIF Assumption | Real-World Failure | Why the "Trick" was Needed |
|---|---|---|
| **Hyperplanes extend to $\infty$** | A cut dividing two local points slices through empty space at $X = 80,000$. | Needed exterior distance decay to force outer monotonicity. |
| **Normal vector $n$ is independent of data** | In sparse or manifold data, a random normal cuts across empty voids or misses high-variance axes. | Needed pairwise $(x_2 - x_1)$ interpolation for $p$ and coordinate scaling. |
| **Path length = integer edge hops** | A step of distance 0.001 in a dense core counts the same as a step of distance 10,000 across a void. | Needed `NEAREST` relative distance scaling. |
| **Linear separating primitives** | Surrounding a circular void (donut hole) requires 6–10 straight cuts, bridging clusters in between. | Needed leaf-level density adjustments. |

---

## 3. Four Core Algorithmic Evolutions: Theory vs. Empirical Findings

### Evolution 1: Bounded-Cell Partitioning & Stadium Outer Space
**The Concept:** Instead of an infinite hyperplane in $\mathbb{R}^D$, every tree node represents a **bounded convex polytope** (or oriented bounding box) $C_v \subset \mathbb{R}^D$.

- At the root node: $C_{\text{root}}$ is the calibrated bounding envelope of the training data:
  $$C_{\text{root}} = [\min_j - \text{margin}_j, \; \max_j + \text{margin}_j]$$
- When a node splits: Hyperplane $H$ only partitions the cell $C_v$ into:
  $$C_{\text{left}} = C_v \cap H^-, \quad C_{\text{right}} = C_v \cap H^+$$

#### Empirical Validation in GEIF & CEIF:
In practice, full polytope mesh storage in memory during inference is computationally prohibitive ($O(N \cdot D)$ facets per cell). However, the **mathematical essence of bounded partitioning was successfully achieved** via two unified mechanisms:
1. **Euclidean Stadium Outer Decay:** Any point outside the bounding envelope $C_{\text{root}}$ bypasses internal tree traversal artifacts and decays exponentially toward 1.0 based on its true Euclidean distance $d_{\text{out}}$ to the data hull:
   $$S(x) = 1.0 - (1.0 - S_{\text{edge}}) e^{-\beta \cdot d_{\text{out}}}$$
2. **Zero Kelvin Lower Bound Calibration:** By traversing all trees to derive the theoretical maximum tree depth $\bar{h}_{\max}$, the minimum score is pinned to the structural lower bound $S_{\min} = 2^{-\bar{h}_{\max} / c}$, eliminating artificial score distortion while preserving headroom for inliers.

---

### Evolution 2: Voronoi Bisector Splitting vs. Empirical Findings

#### The Initial Theoretical Hypothesis:
Rather than generating an independent random normal $n \sim \mathcal{N}(0, I)$ and an intercept $p$, generate splits directly as **Voronoi perpendicular bisectors** between random sample pairs $(x_a, x_b) \in S_{\text{node}}$:
$$n = \frac{x_b - x_a}{\|x_b - x_a\|}, \quad p = \frac{x_a + x_b}{2}$$

The theoretical appeal was intuitive: every cut separates at least two real data points, cuts automatically align with the principal variance of local clusters, and dimension scaling seemed naturally absorbed.

#### Empirical Results from GEIF Testing:
Extensive testing across synthetic (`2blob`, `complex2d`) and real-world benchmark datasets demonstrated that **pure Voronoi bisector splitting does NOT yield better results than CEIF's sample-anchored Gaussian hyperplanes**, and in several key metrics performs noticeably worse. In commit `68582b0`, GEIF replaced Voronoi bisectors with isotropic Gaussian normals.

#### Why Voronoi Bisectors Failed (Root Causes):

1. **Ensemble Diversity Collapse (Angular Starvation):**
   The mathematical power of an Isolation Forest relies heavily on **high spherical angular entropy**—having hundreds of trees slicing the feature space from every continuous angle. 
   When split normals are constrained strictly to sample difference vectors $n = x_b - x_a$, the cuts become heavily correlated with the internal chord directions of clusters. Instead of isotropic multi-angle carving, trees generate repetitive, nearly parallel cuts. This loss of ensemble diversity reduces the forest's ability to smoothly approximate curved or non-convex boundaries.

2. **High-Dimensional Chord Degeneracy:**
   In higher dimensions ($D > 3$), due to the concentration of distances, pairwise difference vectors between randomly chosen samples become noisy and poorly conditioned. A single outlier or edge sample in a node frequently gets paired with a cluster core sample, producing an aggressive cut that isolates the outlier prematurely but slices the rest of the cluster at an unfavorable angle.

3. **Polygonal Facet Artifacts vs. Smooth Probability Contours:**
   Because Voronoi bisectors are pinned to discrete sample pairs, decision contours exhibit jagged, polygonal facet edges with sharp vertices. In contrast, Gaussian cuts produce smooth, continuous, probabilistic density isolines when integrated across an ensemble.

#### The Winning Architecture: Data-Anchored Isotropic Gaussian Cuts
The superior solution that emerged from CEIF and GEIF is a hybrid approach:
- **Direction:** Draw $n \sim \mathcal{N}(0, I)$ isotropically (preserving $360^\circ$ continuous angular entropy).
- **Placement:** Anchor $p$ to actual sample projections (using pairwise sample projection intervals $[z_1, z_2]$ and coordinate stretching).
This completely eliminates blind cuts into infinite voids without sacrificing spherical ensemble diversity.

---

### Evolution 3: Continuous / Metric Path Length (Solving the Density Bias)
**The Concept:** In classic iForest, every split edge adds exactly $+1$ to path length, regardless of physical scale. In a continuous geometric forest, edge traversal accumulates a **metric distance weight**:

$$\Delta h = \frac{\text{dist}(x, \text{boundary})}{\text{scale}_{\text{node}}}$$
or
$$\Delta h = \frac{\text{Volume}(C_{\text{child}})}{\text{Volume}(C_{\text{parent}})}$$

#### Empirical Validation in GEIF:
GEIF validated this principle at the leaf level using a **Cauchy-Lorentz relative distance kernel**:
When a query point arrives at a leaf node containing samples $s_1 \dots s_k$, it evaluates its minimum Euclidean distance to the nearest leaf sample:
$$\text{dist}_{\min} = \min_{i} \|x - s_i\|$$
$$\text{rel-dist} = \frac{\text{dist}_{\min}}{\text{dist}_{\text{avg}}} + \text{MIN-REL-DIST}$$
The effective sample count $n$ is modulated by the spatial density:
$$n' = \frac{n}{\text{rel-dist}}$$
and the leaf contributes $c(n')$ to the depth sum. This prevents isolated points falling into large, sparse leaves from mimicking the score of dense cluster cores.

---

### Evolution 4: Quadratic & Spherical Primitives (Native Cavity & Void Resolution)
**The Concept:** Linear hyperplanes are degree-1 polynomials. Non-convex topologies (donuts, crescent canals, interlocking rings) require degree-2 primitives:

A split can choose between:
1. **A Linear Hyperplane**: $(x - p) \cdot n = 0$ (for separating distinct clusters).
2. **A Hyperspherical Bubble**: $\|x - c\|^2 \le R^2$ (for enclosing clusters or isolating central voids).

```
        Linear Split (EIF)                     Spherical Bubble Split
         \                                             . - ~ - .
          \     Data Cluster                         :     c     :  (Cluster or Void
           \                                          .  (R)   .     isolated in 1 cut!)
            \                                            ' - ~ - '
```

- **Algorithmic Implication:**
  - A central cavity (such as the donut hole in `complex2d`) can be isolated in a **single spherical split**, rather than requiring dozens of intersecting planar cuts.
  - Completely eliminates the "bridging effect" where linear planes accidentally join two disconnected clusters across a central void.

---

## 4. Architectural Comparison: EIF vs. CEIF vs. GEIF

| Architectural Dimension | Classic EIF (Hariri 2019) | CEIF (Engineered SOTA) | GEIF (Tested & Validated) |
|---|---|---|---|
| **Split Primitive** | Infinite flat hyperplane | Linear hyperplane + pairwise $p$ | Linear hyperplane + sample-bounded $p$ |
| **Normal Vector Generation** | Isotropic Gaussian $\mathcal{N}(0, I)$ | Stratified / Isotropic + `auto_weigth` | Isotropic Gaussian with coordinate stretching |
| **Voronoi Bisector Cuts** | Not evaluated | Tested & rejected (loss of angular diversity) | Tested & rejected (causes angular collapse) |
| **Void & Cavity Damping** | None (infinite rays) | `NEAREST` distance scaling at leaves | Continuous Cauchy-Lorentz leaf kernel |
| **Outer Space Behavior** | Severe starburst rays & wedges | Euclidean exterior distance decay ($\beta = 0.10$) | Euclidean Stadium decay + continuous $d_{\text{out}}$ |
| **Lower Bound Calibration** | Clamped at arbitrary 0.0 | Structural Zero Kelvin leaf calibration | Theoretical deepest-leaf calibration ($S_{\min}$) |
| **Persistence Footprint** | Complete tree hierarchy | Lightweight sample reservoir / JSON schema | Dynamic in-memory tree reconstruction |

---

## 5. New Hypotheses & Theories to Test with GEIF

With Voronoi bisector splitting disproven and outer decay / leaf damping established as solid foundations, GEIF provides an ideal experimental testbed for several high-potential algorithmic theories:

### Theory 1: Dual-Mode Splitting with Hyperspherical Bubbles (Quadratic Primitives)
- **Problem:** In non-convex geometries (like concentric rings or donut holes in `complex2d`), planar cuts struggle to isolate cavities without creating "bridges" between disconnected clusters.
- **Hypothesis:** Allow trees to randomly select between a linear hyperplane (70%–80%) and a hyperspherical bubble split (20%–30%):
  $$\|x - c\|^2 \le R^2$$
  where center $c$ is a sampled data point and radius $R$ is sampled between intra-node sample distances.
- **Expected Outcome:** Single-cut isolation of circular voids and compact clusters, drastically reducing tree depth and eliminating planar bridging artifacts.

---

### Theory 2: Subspace Feature Bagging (High-Dimensional Sparsity)
- **Problem:** In high dimensions ($D \ge 10$), drawing an isotropic Gaussian normal with non-zero weights in all dimensions dilutes anomalous signals across uninformative noise features (the curse of dimensionality).
- **Hypothesis:** Implement node-level or tree-level subspace sampling: for each split, randomly select a subset of $k \ll D$ active dimensions (e.g., $k = \lceil\sqrt{D}\rceil$ or $k = 2, 3$) and generate the normal vector only in that subspace.
- **Expected Outcome:** Significantly higher detection sensitivity on sparse, high-dimensional datasets and tabular feature spaces with noisy attributes.

---

### Theory 3: Anisotropic / Mahalanobis Ellipsoidal Cuts
- **Problem:** Real-world metrics and telemetry features frequently exhibit strong covariance (e.g. CPU vs. Memory usage, latency vs. throughput). Spherical cuts slice awkwardly across diagonal correlations.
- **Hypothesis:** Instead of isotropic hyperspheres, use local covariance-weighted ellipsoidal splits:
  $$(x - c)^T \Sigma^{-1} (x - c) \le R^2$$
  where $\Sigma$ is approximated using diagonal feature variance or pairwise two-sample difference covariance.
- **Expected Outcome:** Trees adapt tightly to elongated correlation manifolds with far fewer cuts.

---

### Theory 4: Continuous Path Integration (Intermediate Void Traversal Weighting)
- **Problem:** Currently, internal tree edges add a uniform $+1.0$ hop to depth, and continuous metric distance is only evaluated at the final leaf.
- **Hypothesis:** Accumulate distance-weighted increments $\Delta h$ along internal split edges:
  $$\Delta h_i = 1.0 + \alpha \cdot \frac{\text{dist}(x, H)}{\text{margin}_{\text{node}}}$$
  Points traversing through large empty internal gaps between clusters accumulate larger depth reductions immediately, without waiting to reach a leaf.
- **Expected Outcome:** Faster discrimination of internal cavities and multi-modal cluster boundaries.

---

### Theory 5: Dynamic Ensemble Temperature / Multi-Scale Headroom
- **Problem:** The Zero Kelvin calibration ($S_{\min} = 2^{-\bar{h}_{\max}/c}$) establishes an absolute theoretical lower bound, but operational alerting pipelines often require tunable discrimination between core inliers and boundary samples.
- **Hypothesis:** Introduce a continuous temperature / sharpness exponent $\tau$:
  $$S(x) = 2^{-(h(x) / c)^\tau}$$
  Adjusting $\tau$ modulates the steepness of the anomaly knee without altering tree structures or model weights.
- **Expected Outcome:** Flexible operational tuning for production alert thresholds (e.g. centering the alert knee cleanly around $0.65$).

---

### Theory 6: Streaming Reservoir Aging & Time-Decayed Anomaly Scoring
- **Problem:** Production workloads undergo gradual seasonal drift; obsolete patterns in the reservoir can dilute anomaly detection for emerging distributions.
- **Hypothesis:** Augment the reservoir with exponential time decay weights $w_i = e^{-\lambda (t_{\text{now}} - t_i)}$. When evaluating leaf sample density $c(n')$, weight samples by their recency.
- **Expected Outcome:** Continuous self-adapting anomaly detection on streaming, non-stationary time series without retraining from scratch.
