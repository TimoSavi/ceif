## Categorized extended isolation forest tool

This is a simple command-line program for anomaly detection based on extended isolation forest by [Hariri et al.](https://arxiv.org/abs/1811.02141). 
It is mentioned for environments were are lot of different simple datasets to be analyzed, but making a separate program for each set is too tedious task.
It has practical extensions like:

* Some input data fields can be used as a category field. Effectively this creates an own forest for each category making each category data independent.
* One input data field can be used as a label field. Label field is an unique label for each input row, making it easy to identify outlier data (e.g. timestamp or row id)
* Forests can be saved to file to be used later in analysis
* Sampling is done using reservoir sampling. This allows very large training data to be used.
* Existing forests can be enhanced by new data. New data is added using reservoir sampling.

See more documents in [docs](docs).

### Algorithm change
The original algorithm has some problems with certain types of datasets. This is due to the selection method of random intercept point ***p***. 
Interception ***p*** is selected from rectangular area and if data is uniformly distributed over rectangular area then all sub-spaces divided by random slopes contain sample points.
This causes all sub-spaces to infinity to have inliers. This gives a anomaly score to be app. 0.5 for the whole space.

To tackle this here the ***p*** selection has following steps:

1. Select a random sample point
2. Calculate a random adjustment vector ***a*** from standard normal distribution [0,1]. The length &#124;***a***&#124; is proportional to:
  * Average sample point distance
  * Tree height (larger at tree root)
  * User given parameter (-R)
3. Interception ***p*** is calculated by adding the ***a*** to randomly selected sample point.

This has following effects:

1. There will always be some ***p***s outside the sample area
2. Most ***p***s tend to accumulated where the data is already at beginning of the building of trees

The ***p*** selection area is effectively an enlarged sample point area and not rectangular area which can cause anomalies.
