## Categorized Extended Isolation Forest Tool

This is a simple command-line program for anomaly detection based on the Extended Isolation Forest algorithm by [Hariri et al.](https://arxiv.org/abs/1811.02141). 
It is intended for environments where there are many diverse simple datasets to be analyzed, but creating a separate program for each dataset is too tedious a task.
It features practical extensions such as:

* Input data fields can be used as category fields. Effectively, this creates an independent forest for each category.
* One input data field can be designated as a label field. The label field provides a unique identifier for each input row, making it easy to identify outlier records (e.g., timestamps or row IDs).
* Forests can be saved to a file for later analysis.
* Sampling is performed using reservoir sampling, allowing very large training datasets to be processed with a bounded memory footprint.
* Existing forests can be continuously updated with new data using reservoir sampling.

See more documentation in the [docs](docs/README.md) directory.

### Algorithm Modifications
#### Selection of Intercept Point ***p***
The original algorithm encounters issues with certain types of datasets due to its selection method for the random intercept point ***p***. 
In the canonical algorithm, intercept point ***p*** is selected uniformly from a rectangular bounding box. If data is uniformly distributed over this rectangular area, all subspaces divided by random hyperplanes will contain sample points.
This causes unbounded subspaces extending to infinity to contain inliers, forcing the anomaly score to approximately 0.5 across the entire space.

To resolve this, the selection of ***p*** in this tool follows these steps:

1. Select a random sample point.
2. Calculate a random adjustment vector ***a*** drawn from a standard normal distribution $\mathcal{N}(0, 1)$. The length &#124;***a***&#124; is proportional to:
  * Tree height (larger at the tree root, decreasing toward leaves)
  * Dimension value range (a wider range results in a larger adjustment)
3. The intercept point ***p*** is calculated by adding ***a*** to the randomly selected sample point.

This has the following effects:

1. There will always be some intercept points ***p*** outside the sample area.
2. Most intercept points ***p*** tend to accumulate where data is already present during early tree construction.

The resulting ***p*** selection area is effectively an enlarged sample boundary rather than an unconstrained bounding box that can create artifact inliers.

#### Nearest Training Point Distance in Leaf Nodes
The relative distance between the analyzed point and the nearest training data point in the node is evaluated at leaf nodes.
The absolute distance is scaled to a relative distance using the average sample distance within the tree. If the relative distance is larger
than average, the score is incremented; if the distance is smaller, the score is reduced.

The scale of dimension attribute values can also be adjusted; see the [tweaking document](docs/tweaking.md) for details.

### Acknowledgements
Thanks to codeplea for tinyexpr (https://github.com/codeplea/tinyexpr).
