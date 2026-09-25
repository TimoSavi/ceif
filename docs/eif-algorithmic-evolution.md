# Beyond EIF: Fundamental Algorithmic Evolution for Isolation Forests

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
- Exterior distance decay in deep space

While highly effective, these mechanisms operate as "guardrails" around the core EIF tree. This paper examines how the isolation tree algorithm itself can be fundamentally redesigned so that these properties emerge **naturally from first principles**, transforming EIF into a truly native, next-generation geometric anomaly detector.

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
| **Normal vector $n$ is independent of data** | In sparse or manifold data, a random normal cuts across empty voids or misses high-variance axes. | Needed pairwise $(x_2 - x_1)$ interpolation for $p$. |
| **Path length = integer edge hops** | A step of distance 0.001 in a dense core counts the same as a step of distance 10,000 across a void. | Needed `NEAREST` relative distance scaling. |
| **Linear separating primitives** | Surrounding a circular void (donut hole) requires 6–10 straight cuts, bridging clusters in between. | Needed leaf-level density adjustments. |

---

## 3. Four Core Algorithmic Evolutions

### Evolution 1: Bounded-Cell Partitioning (Cell-Bounded Trees)
**The Concept:** Instead of an infinite hyperplane in $\mathbb{R}^D$, every tree node represents a **bounded convex polytope** (or oriented bounding box) $C_v \subset \mathbb{R}^D$.

- At the root node: $C_{\text{root}}$ is the calibrated bounding envelope of the training data:
  $$C_{\text{root}} = [\min_j - \text{margin}_j, \max_j + \text{margin}_j]$$
- When a node splits: Hyperplane $H$ only partitions the cell $C_v$ into:
  $$C_{\text{left}} = C_v \cap H^-, \quad C_{\text{right}} = C_v \cap H^+$$
- **Algorithmic Implication:**
  - If a test point $x \notin C_{\text{root}}$, it is **immediately recognized as an exterior point** at depth 0.
  - No internal hyperplane ever leaks into outer space.
  - The infinite hyperplane artifact is mathematically eliminated at the root level without any secondary patch.

---

### Evolution 2: Voronoi / Random Projection Bisectors (Natural Pairing)
**The Concept:** Rather than generating an independent random normal $n \sim \mathcal{N}(0, I)$ and an independent intercept $p$, generate splits directly from the **data manifold geometry**:

1. At each node, select two distinct sample points $x_a, x_b \in S_{\text{node}}$ with probability proportional to their distance:
   $$P(x_a, x_b) \propto \|x_a - x_b\|^2$$
2. The split normal is naturally the difference vector:
   $$n = \frac{x_b - x_a}{\|x_b - x_a\|}$$
3. The split point is placed along their segment with jitter:
   $$p = \frac{x_a + x_b}{2} + u \cdot n$$

- **Algorithmic Implication:**
  - Splits automatically align with the principal variance of local clusters.
  - Dimension scaling is naturally absorbed because $n$ is formed directly from actual sample differences.
  - Eliminates "blind cuts" that slice empty space without separating actual data.

---

### Evolution 3: Continuous / Metric Path Length (Solving the Density Bias)
**The Concept:** In classic iForest, every split edge adds exactly $+1$ to path length, regardless of physical scale. In a continuous geometric forest, edge traversal accumulates a **metric distance weight**:

$$\Delta h = \frac{\text{dist}(x, \text{boundary})}{\text{node\_scale}}$$
or
$$\Delta h = \frac{\text{Volume}(C_{\text{child}})}{\text{Volume}(C_{\text{parent}})}$$

- **Algorithmic Implication:**
  - Moving across a vast empty void yields a massive $\Delta h$ drop (instantly revealing isolation).
  - Navigating within a compact, dense cluster accumulates small increments.
  - Solves the density-bias problem natively: dense clusters and sparse clusters achieve commensurate anomaly scales without post-hoc leaf calculations.

---

### Evolution 4: Quadratic & Spherical Primitives (Native Cavity & Void Resolution)
**The Concept:** Linear hyperplanes are degree-1 polynomials. Non-convex topologies (donuts, crescent canals, interlocking rings) require degree-2 primitives:

A split can randomly choose between:
1. **A Linear Hyperplane**: $(x - p) \cdot n = 0$ (for separating clusters).
2. **A Hyperspherical Bubble**: $\|x - c\|^2 \le R^2$ (for enclosing clusters or isolating central voids).

```
        Linear Split (EIF)                     Spherical Bubble Split
         \                                             . - ~ - .
          \     Data Cluster                         :     c     :  (Cluster or Void
           \                                          .  (R)   .     isolated in 1 cut!)
            \                                            ' - ~ - '
```

- **Algorithmic Implication:**
  - A central cavity (like the donut hole in `complex2d`) is isolated in a **single spherical split**, rather than requiring dozens of intersecting planar cuts.
  - Completely eliminates the "bridging effect" where linear planes accidentally join two disconnected clusters.

---

## 4. Architectural Comparison: EIF vs CEIF vs Native Geometric Forest

| Feature | Classic EIF (2019) | Current CEIF (Engineered) | Native Geometric Forest (Next-Gen) |
|---|---|---|---|
| **Split Primitive** | Infinite flat hyperplane | Infinite flat hyperplane + pairwise $p$ | Bounded polytope + optional quadratic bubble |
| **Normal Vector $n$** | Isotropic Gaussian $\mathcal{N}(0, I)$ | Stratified / Isotropic + `auto_weigth` | Sample bisector $x_b - x_a$ (data-driven) |
| **Cut Placement $p$** | Uniform $[z_{\min}, z_{\max}]$ | Pairwise sample interpolation + quadratic decay | Midpoint bisector with margin jitter |
| **Inner Density / Voids** | Ignored (pure edge hops) | `NEAREST` distance weighting at leaves | Metric / volume-weighted continuous path length |
| **Outer Space Behavior** | Severe starburst rays & wedges | Exterior exponential distance decay | Cell-bounded root envelope $C_{\text{root}}$ |
| **Scoring Consistency** | Inconsistent across datasets | Min/Max score scaling to $0..1$ | Standardized continuous path-to-volume ratio |

---

## 5. Roadmap: From "Tricks" to Algorithmic Elegance

If we want to evolve `ceif` from an engineered EIF into a pioneering **Geometric Isolation Forest (GIF / CEIF 3.0)**, the natural evolution path is:

1. **Step 1: Unify Inner & Outer Distance into Tree Traversal**
   Instead of checking exterior distance outside `_score()` and `NEAREST` inside leaves, define a continuous traversal metric:
   Every node carries its local bounding radius $R_v$. As a query point travels, distance relative to $R_v$ continuously modulates the traversal step.

2. **Step 2: Voronoi Bisector Cuts as Standard Split Mode**
   Replace independent normal generation with sample-pair difference vectors ($n = x_2 - x_1$). This unifies `generate_p` and `calculate_n` into a single coherent geometric operation.

3. **Step 3: Dual-Mode Splitting (Planes + Spheres)**
   Allow 20–30% of tree splits to be spherical bounds ($\|x - c\| < R$). This provides native, effortless wrapping of cavities, donut holes, and non-linear manifolds with drastically fewer trees.
