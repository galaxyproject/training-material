---
layout: tutorial_hands_on

title: 'Clustering in Machine Learning'
zenodo_link: https://zenodo.org/record/3813447
questions:
- How to use clustering algorithms to categorized data in different clusters
objectives:
- Learn clustering background
- Learn hierarchical clustering algorithm
- Learn k-means clustering algorithm
- Learn DBSCAN clustering algorithm
- Apply clustering algorithms to different datasets
- Learn how to visualize clusters
key_points:
- Using clustering methods, clusters inside a dataset are drawn using hierarchical, k-means and DBSCAN
- For each clustering algorithm, the number of clusters and their respective hyperparameters should be optimised based on the dataset
time_estimation: 2H
contributors:
- khanteymoori
- anuprulez
---

# Introduction
{:.no_toc}

The goal of unsupervised learning is to discover hidden structures or patterns in unlabeled training data. In this tutorial, we will discuss an unsupervised learning task called clustering. Clustering is the grouping of specific objects based on their characteristics and their similarities. These groups are called clusters. A cluster consists of data objects with high inter-similarity and low intra-similarity. It means members of the same cluster, are more similar to each other by a given metric than they are to the members of the other clusters. The goal of clustering is to subdivide a set of items in such a way that similar items fall into the same cluster, whereas dissimilar items fall in different clusters. This brings up two questions: first, how do we decide what is similar; and second, how do we use this information to cluster the items.

Clustering is central to many bioinformatics research. In particular, clustering helps at analyzing unstructured and high-dimensional data in the form of sequences, expressions, texts and images. For example, clustering is often used in gene expression analysis to find groups of genes with similar expression patterns which provides a useful understanding of gene functions, cellular processes, subtypes of cells and gene regulations. Please refer to [ref1](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0171429) and [ref2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5135122/).

We represent an observation/sample as an n-dimensional vector and many such vectors constitute a dataset. To show an example, let us assume that a dataset, shown in Figure 1, contains many samples and each sample has two dimensions each:

>    ![data](images/data_before_clustering.png "Sample data before clustering")

Clustering reveals the following three groups, indicated by different colors:

>    ![data](images/data_after_clustering.png "Sample data after clustering")


Clustering can be divided into two subgroups:

- Hard: Each data point either belongs to a cluster completely or not.
	
- Soft: Instead of putting each data point into a separate cluster, a probability or likelihood of that data point to be in those clusters is assigned. 

The goal of clustering is to determine the internal grouping in a set of unlabeled data. But how to decide what constitutes a good clustering? It can be shown that there is no absolute best criterion which would be independent of the final aim of the clustering approach. Consequently, the users should apply an appropriate criterion based on the problem.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Types of clustering algorithms

Since clustering is a subjective task, there are many algorithms for data clustering. Every method follows a different set of rules for defining similarity. There are many known clustering algorithms. But, only a few of them are popular and can be categorized as follows: 

 - Connectivity models: As the name suggests, these algorithms are based on the notion that the data points closer in data space exhibit more similarity to each other than the data points lying farther away. These models can follow two approaches - in the first approach, they start with classifying all data points into separate clusters and then aggregating them as the distance decreases. In the second approach, all data points are classified as a single cluster and then partitioned as the distance increases. Also, the choice of distance function is subjective. These models are very easy to interpret but lack scalability for handling big datasets. Examples are hierarchical clustering and its variants.

 - Centroid models: These are iterative clustering algorithms in which the notion of similarity is derived by the closeness of a data point to the centroid of the clusters. [k-means](https://en.wikipedia.org/wiki/K-means_clustering) clustering is a popular algorithm that falls into this category. It runs iteratively to find a fixed number of clusters which needs to be specified before running it. Therefore, to fix this number of clusters, it is required to have some prior knowledge about the dataset.

 - Density models: These algorithms search the data space for areas of varied density of data points. They isolate various density regions and assign the data points within these regions in the same cluster. A popular example of a density model is [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN).

 - Distribution models: These clustering algorithms are based on the notion of how probable it is that all data points in a cluster belong to the same distribution such as Gaussian. These algorithms often suffer from overfitting. A popular example is [expectation-maximization](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm) (EM) algorithm which uses multivariate normal distributions.


Now, in this tutorial, we will go through three of the most popular clustering algorithms in detail, hierarchical clustering, k-means, DBSCAN, and a comparison between these methods. Further, we will discuss their parameters and how to apply them to find clusters in the [iris flower dataset](https://en.wikipedia.org/wiki/Iris_flower_data_set) and a few other datasets.


# Clustering distance measures

Since clustering is the grouping of similar objects, a measure that can determine whether two objects are similar or dissimilar is required. There are two main types of measures used to estimate this relation - distance and similarity measures. The notions of distance and similarity are related. The smaller the distance between two objects, the more similar they are to each other. All measures refer to the feature values in some way, but they consider different properties of the feature vector. There is no optimal similarity measure as its usage depends on the problem.

Many clustering algorithms use distance measures to determine the similarity or dissimilarity between any pair of objects. A valid distance measure should be symmetric and obtains its minimum value (usually zero) in case of identical vectors. Clustering requires some methods for computing the distance or (dis)similarity between each pair of observations. The result of this computation is known as a dissimilarity or distance matrix.

The choice of a distance measure is critical in clustering. It defines how the similarity of two elements `(x, y)` is calculated and it will influence the shape of the clusters. The classical distance measures are [euclidean](https://en.wikipedia.org/wiki/Euclidean_distance) and [manhattan](https://en.wikipedia.org/wiki/Taxicab_geometry) distances. For the most common clustering algorithms, the default distance measure is euclidean. If the euclidean distance is chosen, then observations having high magnitudes of their respective features will be clustered together. The same holds for the observations having low magnitudes of their respective features. In Figure 3, we group the cells using euclidean distance and their distance matrix.

![Distances](images/raceid_distance.svg "Euclidean distance between three points (R, P, V) across three features (G1, G2, G3)")


 > ### {% icon question %} Questions
 >
 > 1. Why are there zeroes along the diagonal of the above example distance matrix?
 > 1. Is there any symmetry in this matrix?
 >
 > > ### {% icon solution %} Solution
 > >
 > > 1. The distance between a point to itself is zero.
 > > 1. The distance between point *a* to point *b* is the same as the distance between point *b* to point *a* using the Euclidean distance metric.
 > >
 > {: .solution }
 >
 {: .question }

Other dissimilarity measures exist such as correlation-based distances, which are widely used for gene expression data analyses. Correlation-based distance considers two objects to be similar if their features are highly correlated, even though the observed values may be far apart in terms of euclidean distance. The distance between the two objects is 0 when they are perfectly correlated. [Pearson’s correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) is quite sensitive to outliers. This does not matter when clustering samples, because the correlation is over thousands of genes. During the clustering of genes, it is important to be aware of the possible impact of outliers. This can be mitigated by using [Spearman’s correlation](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient) instead of Pearson’s correlation.

# Different clustering approaches

## Hierarchical clustering

Before using hierarchical clustering, let us first understand the theory behind it. Hierarchical clustering, as the name suggests is an algorithm that builds a hierarchy of clusters. It starts with all the data points assigned to a cluster of their own. Then, the two nearest clusters are merged into the same cluster. In the end, this algorithm terminates when there is only a single cluster left.

Following are the steps that are performed during hierarchical clustering:

1. In the beginning, every data point in the dataset is treated as a cluster which means that we have `N` clusters at the beginning of the algorithm for a dataset of size `N`.

2. The distance between all the points is calculated and two points closest to each other are joined together to a form a cluster. 

3. Next, the point which is closest to the cluster formed in step 2, will be joined to the cluster.

4. Steps 2 and 3 are repeated until one big cluster is formed.

5. Finally, the big cluster is divided into K small clusters with the help of dendrograms. 

Let’s now see how dendrograms help in hierarchical clustering.

>    ![data](images/Hierarchical_clustering_1.png "Hierarchical clustering")

At the bottom, we start with data points, each assigned to separate clusters. Two closest clusters are then merged till we have just one cluster at the top. The height in the dendrogram at which two clusters are merged represents the distance between two clusters in the data space.

The decision of the number of clusters that can best depict different groups can be chosen by observing the dendrogram. The best choice of the number of clusters is the number of vertical lines in the dendrogram cut by a horizontal line that can transverse the maximum distance vertically without intersecting a cluster.

In the above example, the best choice of the number of clusters will be 4 as the red horizontal line in the dendrogram below covers maximum vertical distance AB.
>    ![data](images/Hierarchical_clustering_2.png "Hierarchical clustering")


This algorithm has been implemented above using the bottom-up approach. It is also possible to follow the top-down approach starting with all data points assigned in the same cluster and recursively performing splits till each data point is assigned a separate cluster. The decision of merging two clusters is taken based on the proximity of these clusters.

> ### {% icon comment %} Background of the iris dataset
> The iris flower dataset or Fisher’s iris dataset is a multivariate dataset introduced by the British statistician and biologist Ronald Fisher in his 1936 paper ({% cite Fisher1936 %}).
> Each row of the table represents an iris flower, including its species and dimensions of its botanical parts, sepal and petal, in centimeters.
> For more history of this dataset read here [Wikipedia](https://en.wikipedia.org/wiki/Iris_flower_data_set).
{: .comment}


At the first step, we should upload the iris dataset and two other datasets which will be used at the end of the tutorial.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. **Import** {% icon galaxy-upload %} the file `iris.csv` from [Zenodo](https://zenodo.org/record/3813447/files/iris.csv) or from the data library
>
>    ```
>    https://zenodo.org/record/3813447/files/iris.csv
>    https://zenodo.org/record/3813447/files/circles.csv
>    https://zenodo.org/record/3813447/files/moon.csv
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
>
> 2. **Rename** {% icon galaxy-pencil %} the dataset to `iris`
>
>    {% include snippets/rename_dataset.md %}
>
> 3. Check the **datatype**
>    - Click on the history item to expand it to get more information.
>    - The datatype of the iris dataset should be `csv`.
>    - **Change** {% icon galaxy-pencil %} the datatype *if* it is different than `csv`.
>      - Option 1: Datatypes can be **autodetected**
>      - Option 2: Datatypes can be **manually set**
>
>    {% include snippets/detect_datatype.md datatype="datatypes" %}
>    {% include snippets/change_datatype.md datatype="csv" %}
>
{: .hands_on}

Our objective is to categorise similar flowers in different groups (Figure 6). We know that we have **3** species of iris flowers (versicolor, virginica, setosa) with
**50** samples for each. These species look very much alike as shown in the figure below.

![3 species of iris flowers](images/iris_flowers.png "3 species of iris flowers")

In our dataset, we have the following features measured for each flower: [petal](https://en.wikipedia.org/wiki/Petal) length, petal width, [sepal](https://en.wikipedia.org/wiki/Sepal) length, sepal width

Figure 7 shows the dendrogram of these data.

>    ![data](images/Hierarchical_iris.png "Iris data hierarchical clustering")


We will apply hierarchical clustering to iris dataset to find clusters based on two features (of flowers) - sepal length and width.

> ### {% icon hands_on %} Hands-on: Hierarchical clustering
>
> 1. **Numeric Clustering** {% icon tool %} with the following clustering parameters:
>    - *"Select the format of input data"*: `Tabular Format (tabular,txt)`
>        - {% icon param-file %} *"Data file with numeric values"*: `iris`
>        - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>        - {% icon param-select %} *"Choose how to select data by column"*: `All columns EXCLUDING some by column header name(s)`
>            - {% icon param-text %} *"Type header name(s)"*: `Species`
>        - {% icon param-select %} *"Clustering Algorithm"*: `Hierarchical Agglomerative Clustering`
>        - In *"Advanced option"*    
>            - {% icon param-text %} *"Number of clusters"*: `2`
>            - {% icon param-select %} *"Affinity"*: `Euclidean`
>            - {% icon param-select %} *"Linkage"*: `ward`
>        - In *"Output options"* 
>            - {% icon param-text %} *"width of output"*: `7.0`
>            - {% icon param-text %} *"height of output"*: `5.0`
>            - {% icon param-text %} *"dpi of output"*: `175.0`
> 
> 2. Rename the generated file to `Hierarchical clustering`
{: .hands_on}

If you view the result table, you can see the last column is the label for each cluster and as you see, all the setosa samples are grouped in one cluster and two other species (versicolor and virginica) are grouped in the second cluster. From Figure 3 it is obvious that versicolor and virginica are more similar to each other.

### Visualize hierarchical clustering

The resulting candidate clustering can be visualized using the `Scatterplot with ggplot2` tool. Each sample is color-coded based on its clustering for that sample.
Let's visualize the clustering results to see how groups have been built.

> ### {% icon hands_on %} Hands-on: Visualize hierarchical clustering result
>
> 1. **Scatterplot with ggplot2** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: **Hierarchical clustering**
>    - *"Column to plot on x-axis"*: `1`
>    - *"Column to plot on y-axis"*: `2`
>    - *"Plot title"*: `Hierarchical clustering in iris data`
>    - *"Label for x axis"*: `Sepal length`
>    - *"Label for y axis"*: `Sepal width`
>    - In *"Advanced Options"*:
>        - *"Data point options"*: `User defined point options`
>            - *"relative size of points"*: `2.0`
>        - *"Plotting multiple groups"*: `Plot multiple groups of data on one plot`
>            - *"column differentiating the different groups"*: `6`
>            - *"Color schemes to differentiate your groups"*: `Set 2 - predefined color pallete`
>    - In *"Output options"*:
>        - {% icon param-text %} *"width of output"*: `7.0`
>        - {% icon param-text %} *"height of output"*: `5.0`
>        - {% icon param-text %} *"dpi of output"*: `175.0`
> 
> 2. **View** {% icon galaxy-eye%} the resulting plot
> 3. Rename to `Hierarchical scatter plot`
{: .hands_on}

>    ![data](images/hierarchical_scatter.png "Hierarchical clustering scatter plot")


## K-means clustering

K-means clustering is the most commonly used unsupervised machine learning algorithm for partitioning a given dataset into a set of k clusters, where k represents the number of groups pre-specified by the user.
In k-means clustering, each cluster is represented by its center or centroid which corresponds to the mean of points assigned to the cluster. The basic idea behind k-means clustering consists of defining clusters so that the total intra-cluster variation is minimized.

K-means is popular because of its speed and scalability. There are several k-means algorithms available. The standard algorithm defines the total within-cluster variation as the sum of squared Euclidean distances between items and the corresponding centroid. K is a hyperparameter of the algorithm and the k-means algorithm can be summarized as follows:

1. Specify the number of clusters (k) to be created (to be specified by users).

2. Select k data points or observations randomly from the dataset as the initial cluster centers or means.

3. Assign each data point to their closest centroid, based on the euclidean distance between a data point and its centroid.

4. For each of the k clusters update cluster centroid by calculating the new mean values of all the data points in the cluster.

5. Iteratively minimize the total within the sum of square: iterate steps 3 and 4 until the cluster assignments stop changing or the maximum number of iterations is reached.

The parameters that minimize the cost function are learned through an iterative process of assigning data points to clusters and then moving the clusters. A restriction for the k-means algorithm is that the dataset should be continuous. It won’t work if data is categorical.

> ### {% icon hands_on %} Hands-on: K-means clustering
>
> 1. **Numeric Clustering** {% icon tool %} with the following clustering parameters:
>    - *"Select the format of input data"*: `Tabular Format (tabular,txt)`
>        - {% icon param-file %} *"Data file with numeric values"*: `iris`
>        - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>        - {% icon param-select %} *"Choose how to select data by column"*: `All columns EXCLUDING some by column header name(s)`
>            - {% icon param-text %} *"Type header name(s)"*: `Species`
>        - {% icon param-select %} *"Clustering Algorithm"*: `KMeans`
>        - In *"Advanced option"*    
>            - {% icon param-text %} *"Number of clusters"*: `2`
> 2. Rename the generated file to `k-means clustering`
{: .hands_on}


### Visualize k-means clustering

> ### {% icon hands_on %} Hands-on: Visualize k-means clustering result
>
> 1. **Scatterplot with ggplot2** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: **k-means clustering**
>    - *"Column to plot on x-axis"*: `1`
>    - *"Column to plot on y-axis"*: `2`
>    - *"Plot title"*: `K-means clustering in iris data`
>    - *"Label for x axis"*: `Sepal length`
>    - *"Label for y axis"*: `Sepal width`
>    - In *"Advanced Options"*:
>        - *"Data point options"*: `User defined point options`
>            - *"relative size of points"*: `2.0`
>        - *"Plotting multiple groups"*: `Plot multiple groups of data on one plot`
>            - *"column differentiating the different groups"*: `6`
>            - *"Color schemes to differentiate your groups"*: `Set 2 - predefined color pallete`
>    - In *"Output options"*:
>        - {% icon param-text %} *"width of output"*: `7.0`
>        - {% icon param-text %} *"height of output"*: `5.0`
>        - {% icon param-text %} *"dpi of output"*: `175.0`
> 2. **View** {% icon galaxy-eye%} the resulting plot
> 3. Rename to `k-means scatter plot`
{: .hands_on}

>    ![data](images/k_means_scatter.png "K-means clustering scatter plot")


> ### {% icon question %} Question
>
> How to choose the right number of expected clusters (k)?
>
>
> > ### {% icon solution %} Solution
> >
> > Major difficulty found with k-means is the choice of the number of clusters. Different methods is proposed to solve this problem. 
> > Here, we provide a simple solution. The idea is to compute k-means clustering using different values of clusters k. Next, the within sum of square is drawn according to the number of clusters. The location of a bend (knee) in the plot is generally considered as an indicator of the appropriate number of clusters.
> > ![data](images/number_of_clusters.png "Optimal number of clusters")
> > The plot above represents the variance within the clusters. It decreases as k increases, but it can be seen a bend (or “elbow”) at k = 4. This bend indicates that 
> > additional clusters beyond the fourth have little value. 
> {: .solution}
{: .question}


 > ### {% icon question %} Questions
 >
 > What are the differences between k-means and hierarchical clustering techniques
 >
 > > ### {% icon solution %} Solution
 > >
 > > 1. Hierarchical clustering can’t handle big data well but k-means clustering can. This is because the time complexity of k-means is linear i.e. O(n) while that of hierarchical clustering is quadratic i.e. O(n2).
 > > 
 > > 2. In k-means clustering, since we start with a random choice of clusters, the results produced by running the algorithm multiple times might differ. While results are reproducible in Hierarchical clustering.
 > >
 > > 3. K-means is found to work well when the shape of the clusters is hyperspherical (like circle in 2D, sphere in 3D).
 > >
 > > 4. K-means clustering requires prior knowledge of K i.e. no. of clusters to divide the dataset into. But, you can stop at whatever number of clusters you find appropriate in hierarchical clustering by interpreting the dendrogram
 > >
 > {: .solution }
 >
 {: .question }


## DBSCAN clustering

DBSCAN (Density-based spatial clustering of applications with noise) is a popular clustering algorithm and views clusters as areas of high density separated by areas of low density. Due to this rather generic view, clusters found by DBSCAN can be of any shape, as opposed to k-means which assumes that clusters are convex shaped. The central component to the DBSCAN is the concept of core samples which are present in the areas of high density. A cluster is, therefore, a set of core samples close to one other (measured by some distance measure) and a set of non-core samples that are close to a core sample (but are not core samples themselves). There are two important parameters in DBSCAN algorithm - `min_samples` is the number of samples in a neighborhood for a point to be considered as a core point and `eps` is the maximum distance between two samples for one to be considered as in the neighborhood of the other. Higher the value of `min_samples` or lower the value of eps indicate higher density necessary to form a cluster. DBSCAN does not require one to specify the number of clusters in the data a priori, as opposed to k-means.

> ### {% icon hands_on %} Hands-on: DBSCAN clustering
>
> 1. **Numeric Clustering** {% icon tool %} with the following clustering parameters:
>    - *"Select the format of input data"*: `Tabular Format (tabular,txt)`
>        - {% icon param-file %} *"Data file with numeric values"*: `iris`
>        - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>        - {% icon param-select %} *"Choose how to select data by column"*: `All columns EXCLUDING some by column header name(s)`
>            - {% icon param-text %} *"Type header name(s)"*: `Species`
>        - {% icon param-select %} *"Clustering Algorithm"*: `DBSCAN`
> 2. Rename the generated file to `DBSCAN clustering`
{: .hands_on}


### Visualise DBSCAN clustering

> ### {% icon hands_on %} Hands-on: Visualize DBSCAN clustering result
>
> 1. **Scatterplot with ggplot2** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: **DBSCAN clustering**
>    - *"Column to plot on x-axis"*: `1`
>    - *"Column to plot on y-axis"*: `2`
>    - *"Plot title"*: `DBSCAN Clustering in iris data`
>    - *"Label for x axis"*: `Sepal length`
>    - *"Label for y axis"*: `Sepal width`
>    - In *"Advanced Options"*:
>        - *"Data point options"*: `User defined point options`
>            - *"relative size of points"*: `2.0`
>        - *"Plotting multiple groups"*: `Plot multiple groups of data on one plot`
>            - *"column differentiating the different groups"*: `6`
>            - *"Color schemes to differentiate your groups"*: `Set 2 - predefined color pallete`
>    - In *"Output options"*:
>        - {% icon param-text %} *"width of output"*: `7.0`
>        - {% icon param-text %} *"height of output"*: `5.0`
>        - {% icon param-text %} *"dpi of output"*: `175.0`
> 2. **View** {% icon galaxy-eye%} the resulting plot:
> 3. Rename to `DBSCAN scatter plot`
{: .hands_on}


>    ![data](images/dbscan_scatter.png "DBSCAN clustering scatter plot")


You will also notice that the blue point in the plot is not contained within any cluster. DBSCAN does not necessarily categorize every data point, and is therefore terrific with handling outliers in the dataset. 

> ### {% icon question %} Question
>
> How we can evaluate the clustering results?
>
>
> > ### {% icon solution %} Solution
> >
> > Clustering is an unsupervised learning algorithm; there are no labels or ground truth to compare with the clusters. However, we can still evaluate the performance of the algorithm using intrinsic measures.
> > There is a performance measure for clustering evaluation which is called the [silhouette coefficient](https://en.wikipedia.org/wiki/Silhouette_(clustering)). It is a measure of the compactness and separation of the clusters. 
> > It increases as the quality of the clusters increase; it is large for compact clusters that are far from each other and small for large, overlapping clusters. The silhouette coefficient is calculated per instance; for a set of instances, it is calculated as the mean of the individual samples scores.
> {: .solution}
{: .question}

# Applying clustering algorithms on multiple datasets

You can do the same steps on the other datasets, moon and circles. First, import the data files, [moon.csv](https://zenodo.org/record/3813447/files/moon.csv) and [circles.csv](https://zenodo.org/record/3813447/files/circles.csv) from Zenodo or data library and rename them to `moon` and `circles` respectively. 

## Visualise datasets

> ### {% icon hands_on %} Hands-on: Visualize scatter plot of data
>
> 1. **Scatterplot with ggplot2** {% icon tool %} with the following parameters:
>
>    {% include snippets/select_multiple_datasets.md %}
>    - {% icon param-file %} *"Input tabular dataset"*: `circles` and `moon` as **multiple datasets**
>    - *"Column to plot on x-axis"*: `1`
>    - *"Column to plot on y-axis"*: `2`
>    - *"Plot title"*: `Scatter Plot`
>    - *"Label for x axis"*: `X`
>    - *"Label for y axis"*: `Y`
>    - In *"Output options"*:
>        - {% icon param-text %} *"width of output"*: `7.0`
>        - {% icon param-text %} *"height of output"*: `5.0`
>        - {% icon param-text %} *"dpi of output"*: `175.0`
> 2. **View** {% icon galaxy-eye%} the resulting plots
{: .hands_on}

>    ![data](images/circles_moon_scatter.png "Scatter plot of circles and moon Data")


## Find clusters

Now you can find clusters in these datasets using the aforementioned algorithms.

> ### {% icon hands_on %} Hands-on: Hierarchical clustering of circles and moon datasets
>
> 1. **Numeric Clustering** {% icon tool %} with the following clustering parameters:
>
>    {% include snippets/select_multiple_datasets.md %}
>    - *"Select the format of input data"*: `Tabular Format (tabular,txt)`
>        - {% icon param-file %} *"Data file with numeric values"*: `cirlces` and `moon` as **multiple datasets**
>        - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>        - {% icon param-select %} *"Choose how to select data by column"*: `All`
>        - {% icon param-select %} *"Clustering Algorithm"*: `Hierarchical Agglomerative Clustering`
>        - In *"Advanced option"*    
>            - {% icon param-text %} *"Number of clusters"*: `2`
>            - {% icon param-select %} *"Affinity"*: `Euclidean`
>            - {% icon param-select %} *"Linkage"*: `ward`
> 2. Rename the generated files to `circles hierarchical clustering` and `moon hierarchical clustering` respectively
{: .hands_on}

## Visualise clusters

Then, you can visualize the clustering results using the following steps:

> ### {% icon hands_on %} Hands-on: Visualize hierarchical clustering result on circles and moon datasets.
>
> 1. **Scatterplot with ggplot2** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `circles hierarchical clustering` and `moon hierarchical clustering` as **multiple datasets**
>
>    {% include snippets/select_multiple_datasets.md %}
>    - *"Column to plot on x-axis"*: `1`
>    - *"Column to plot on y-axis"*: `2`
>    - *"Plot title"*: `Hierarchical clustering`
>    - *"Label for x axis"*: `X`
>    - *"Label for y axis"*: `Y`
>    - In *"Advanced Options"*:
>        - *"Data point options"*: `User defined point options`
>            - *"relative size of points"*: `2.0`
>        - *"Plotting multiple groups"*: `Plot multiple groups of data on one plot`
>            - *"column differentiating the different groups"*: `3`
>            - *"Color schemes to differentiate your groups"*: `Set 2 - predefined color pallete`
>    - In *"Output options"*:
>        - {% icon param-text %} *"width of output"*: `7.0`
>        - {% icon param-text %} *"height of output"*: `5.0`
>        - {% icon param-text %} *"dpi of output"*: `175.0`
> 2. **View** {% icon galaxy-eye%} the resulting plots:
{: .hands_on}

In the next steps, you can apply these three algorithms (hierarchical, k-means and DBSCAN) to moon and circles datasets in the same way as explained above. In the k-means algorithm, `k=2` and for the DBSCAN algorithm, the parameters are not the default parameters and should be set as follows: for the circles dataset (`maximum neighborhood distance=0.2` and `minimal core point density=5`) and for the moon datasets (`maximum neighborhood distance=0.3` and `minimal core point density=4`). You can see the scatter plots of the clustering results in Figure 13 and 14.

>    ![data](images/circles_clustering.png "Plot of clustering algorithms on circles dataset")

>    ![data](images/moon_clustering.png "Plot of clustering algorithms on moon dataset")


# Conclusion

In this tutorial, we discussed 3 clustering algorithms which are used to discover structures or patterns in unlabeled data. You learned about the hierarchical, k-means and DBSCAN algorithms. By following steps specified for each clustering tool, we learned how to perform clustering and visualize results using clustering and plotting tools, respectively in Galaxy. There are many other clustering approaches which can be tried out on these datasets to find how they perform and how they compare to the 3 clustering algorithms explained in this tutorial. Different datasets can also be analysed using these algorithms. The clustering algorithms have some parameters which can be altered while performing the analyses to see if they affect the clustering or not. While using clustering algorithms, we need to take care of some important aspects like treating outliers in data and making sure each cluster has sufficient population. Some data pre-processors can also be used to clean the datasets.
