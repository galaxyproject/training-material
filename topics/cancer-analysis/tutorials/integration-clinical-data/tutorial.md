---
layout: tutorial_hands_on

title: "Integration of Clinical Data"
questions: []
objectives: []
key_points: []
time_estimation: 2H
contributions:
  authorship: [larunerdman]
  editing: [hexylena]
  funding: [bioinformatics-ca,erasmusplus]

subtopic: clinical
requirements:
- type: "internal"
  topic_name: data-science
  tutorials:
      - r-basics
      - r-advanced

priority: 11
acronyms:
  IBD: Inflammatory Bowel Disease

notebook:
    language: r
tags:
- R
---

## Loading The Necessary Packages

If you are running this tutorial in Galaxy, consider using the `environment.yml` file to install packages via Conda.

```R
###########
# IBD SNF #
###########

install.packages(c("SNFtool","ExPosition","e1071","RColorBrewer","survival","rms"))
```

Now, make sure you have the following packages loaded:

```R
library(ExPosition) # Chi-squared distance
library(SNFtool) # SNF functions
library(e1071) # hamming.distance()
library(RColorBrewer) # color palettes
library(survival) # survival functions
library(rms) # more advanced survival functions
```

The data we're working is based on {% cite Dhaliwal_2021 %}, however our dataset is
simulated from this data set with cancer outcomes added.
We will perform a similar SNF analysis {% cite Wang_2014 %} to that described in this paper.

```R
## read in data
url <- "https://raw.githubusercontent.com/galaxyproject/training-material/main/topics/cancer-analysis/tutorials/integration-clinical-data/CBW_CAN2021_Module10_LabData.csv"
analysis.df = read.csv(file = url, header = TRUE,as.is = TRUE)
```

Next, let's inspect the data we read in:

```R
head(analysis.df)
str(analysis.df)
```

Here we can see some {IBD} Class labels like `Porto_dx, Porto_ordinal_dx, other_dx, other_ordinal_dx`, as well as survival outcomes like `age` and `can_death`

> <tip-title>R Objects</tip-title>
> There are 3 major datatypes in R that one must be concerned with:
> 1. Vectors, 1 dimensional
>
>    ```
>    my_vec <- c("cat", "dog", "fish", "hamster", "parrot")
>    my_vec[2] == "dog"
>    ```
>
> 2. Data frames, 2 dimensional
>
>    ```
>    my_vec1 <- c("a", "c", "d")
>    my_vec2 <- c(2, 19, 8)
>    my_df <- data.frame(header_name1 = my_vec1, header_name2 = my_vec2)
>    my_df
>    #   header_name1 header_name2
>    # 1            a            2
>    # 2            c           19
>    # 3            d            8
>    colnames(my_df)
>    my_df$header_name1
>    my_df$header_name2
>    ```
>
> 3. Lists, 1 dimensional “bulk” storage
>
>    ```
>    my_list = list("object1" = my_vec, "object2" = my_df)
>    my_list$object1 # Our vector
>    my_list$object2 # Our dataframe
>
>    # Alternative access
>    my_list[['object1']]
>    my_list[['object2']]
>    ```
>
> Summary:
>
> vectors
> :    1 dimensional storage, start with `c(`
>
> data frames
> :    2 dimensional storage, start with `data_frame(`
>
> lists
> :    1 dimensional "bulk" storage, start with `list(`
>
{: .tip}

Now, let's split our dataframe into data-specific data frames. First we'll create column name vectors for each data-type:

```R
# data-type specific column names
hist.features <- c("Granuloma","focal_chronic_duodenitis", "focal_active_colitis",
                   "FEG", "ileitis_mild_cecum", "pattern_involvement_worse.distally")
endosc.features <- c("deep_ulcer_SB", "classic_backwash", "ileal_inflammation", "reverse_gradient",
                     "small._ulcers_SB", "X5_small_ulcers_colon", "deep_ulcer_stomach", "less_5_ulcer_colon",
                     "skip_lesion", "rectal_sparing", "relative_patchiness")
rad.features <- c("SB_thickening")

otherhist.features <-c("basal_plasma_cells", "activity", "gastritis",
                       "duodenitis", "crypt_distortion", "chronic_inflammation")
```

And then create data-type specific data frames.

> <tip-title>What does `names(analysis.df) %in% hist.features` do?</tip-title>
> This creates a vector of TRUE/FALSE for each column name, based on whether or not each name from analysis.df is in hist.features:
>
> ```
> > names(analysis.df) %in% hist.features
>  [1]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE
> [20] FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE
> [39] FALSE
> ```
>
> This is then used to subset the rows.
{: .tip}

We've added a couple `head()` calls in here to preview the data frames as you make them:

```R
## create data-type specific data frames
hist.df <- analysis.df[,names(analysis.df) %in% hist.features]
head(hist.df)
rad.df <-analysis.df[,names(analysis.df) %in% rad.features]
head(rad.df)
endosc.df <-analysis.df[,names(analysis.df) %in% endosc.features]
head(endosc.df)
otherhist.df <-analysis.df[,names(analysis.df) %in% otherhist.features]
head(otherhist.df)
```

Now we'll collect them into a list that we'll use. One list for each distance metric, filled with binary features (yes/no, 1/0). `data.list2` however is multicategorical (values other than 1/0)

```R
## create data list
data.list <- list(hist.df,rad.df,endosc.df)
names(data.list) <- c("hist","rad","endosc")

data.list2 <- list(otherhist.df)
names (data.list2) <-c("otherhist")
```

And with this the data loading and preparation is done!

## Similarity Network Fusion

To perform SNF we need to:

1. Set our model parameters
1. Standardize our datasets
1. Generate matrices representing how different each individual in each data type is from another
1. Use the difference matrices to create 'affinity' or similarity matrices for each data type
1. Fuse data-specific affinity similarity matrices
1. Cluster individuals and evaluate fused clustering

### Setting the Model Parameters

The modeling parameters we need to set are:
- 'Number of neighbors' (size of the group to sparsify matrix using)
- Mu/alpha that handles scaling problems with the similarity metric
- 'Number of iterations'

These values are handily recommended in the SNFtool documentation, either on CRAN, or in RStudio by entering `SNFtool::SNF`, selecting it, and pressing <kbd>F1</kbd>

Please note that these parameters aren't set in stone: evaluate how your clustering changes as you change these!

```R
# First, set all the parameters:
K = 20; # number of neighbors, usually (10~30)
alpha = 0.5; # hyperparameter, usually (0.3~0.8)
T = 20; # Number of Iterations, usually (10~20)
```

### Standardize Our Datasets

Because our data is NOT continuous, we don't need to apply a standard Normal transform to it.

### Generate Distance Matrices

If you have continuous data, it should be separated by from your dichotomous and/or categorical data and normalized using the “standardNormalization” function from the SNFtool package.

We'll be using `lapply` because we want to use the function on a list of each of the genomics datasets. So for every item in `data.list`, we will call the `hamming.distance(as.matrix(x))` function.

```R
## Calculate distance matrices
## (here we calculate Hamming and Chi-Square distances,
## Hamming for binary (0,1 values) and Chi-Square for
## categorical variables)

## If the data are all continuous values, we recommend the users to perform
## standard normalization then using Euclidean distance.

## Calculate the pair-wise distance;
dist.list <- lapply(data.list,
                    function(x){hamming.distance(as.matrix(x))})
str(dist.list)
```

> <tip-title>Anatomy of an R function</tip-title>
> Here is an example function:
>
> ```R
> my_function = function(argument="hello"){return(argument)}
> ```
>
> Part | Meaning
> --- | ---
> `my_function` | The function's name (optional) but can be subsequently used to call `my_function(argument="bye")`
> `function` | Starts a function block
> `argument="hello"` | any arguments your function requires. If values have an `=...`, then that is the optional default value
> `{return(argument)}` | Everything your function does goes inside the `{}`. The return tells your function what to give you back. If not set, will return the last line in function by default
>
> ```R
> my_function() # Returns "hello"
> my_function("yes") # Returns "yes"
> ```
>
> What about this function?
>
> ```R
> my_2nd_function = function(a, b){
> 	a_plus_1 = a+1
> 	b_plus_2 = b+2
> 	my_sum = a_plus_1 + b_plus_2
> 	return(my_sum)
> }
> my_2nd_function(a = 0, b = 0) # What would this return?
> my_2nd_function(a = 3, b = 1) # And this?
> ```
>
{: .tip}

The hamming distance lapply can be understood through this schematic:

![Schematic showing three list objects, hist rad and endosc being transformed. In their initial state they are presented as a series of columns, feat1 and feat1_vec, feat2/feat2_vec and so on, for each of the 3 hist rad and endosc. All of these are transformed via the hammingdistance(as.matrix(x)) function, and the resulting tables are all identical in structure. They are now lists again, except distance matrixes between ind1/ind2/3/4, and ind1/2/3/4. It is not entirely clear where the indicies came from.](./hammingdistance.png)

However for the categorical (non-binary) histology data we'll instead use a `chi2Dist` function, and use `$D` to extract the distance matrix from the output of that function.

```R
dist.list2 <-lapply(data.list2,
                    function(x){chi2Dist((as.matrix(x)))$D})
str(dist.list2)
```

Next we can append each of these into one single large list with everything: `full.dist.list`


```R
# append distance matrix lists
full.dist.list <- append(dist.list,dist.list2)
str(full.dist.list)
# lapply(full.dist.list,names)
```

### Generate "Affinity" / Similarity Matrices

We'll take the list from the last step, and here we'll define a function using the K and alpha we set at the beginning of our script, that will create an 'affinity' or similarity matrix for each distance matrix for each data type.

```R
## next, construct similarity graphs -- we will investigate these more closely later
aff.list <-lapply(full.dist.list,
                  function(x){affinityMatrix(x, K = K,sigma = alpha)})
str(aff.list)
```

This works like the previous `lapply` calls and applies the affinityMatrix function to each element of `full.dist.list`, transforming the distance matrix into an affinity matrix.

### Fuse Affinity Matrices

Now that we've got all of our affinity matricies, we can fuse them together with `SNF`. The SNF function takes our affinity matrices list as the first argument, and the "number of neighbours" and "number of iterations" parameters we set in the beginning.

We'll then also update the row and column names based off of our `studyID` column.

```R
## next, we fuse all the graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(aff.list, K, T)
colnames(W) = rownames(W) = analysis.df$studyID
```

### Cluster fused matrix and evaluate clusters

We will estimate the number of clusters using two different methods

```R
## This function provides two ways to estimate the number of clusters. Note that,
## these two methods cannot guarantee the accuracy of estimated number of
## clusters, but just to offer two insights about the datasets.
estimationResult = estimateNumberOfClustersGivenGraph(W, 2:5)
estimationResult
```

We'll choose the second best option, but you can change this by changing the number below.

```R
# choose
C = estimationResult[[2]] # number of clusters
# we can also set this number manually
# C = 3
```

once you've set the `C` parameter, we can create a vector of cluster assignment:

```R
## you can use spectral clustering to group your data:
group = spectralClustering(W,C) # the final subtypes information
table(group)
displayClusters(W, group)
```

RColorBrewer package has lots of color palettes you can use with heatmaps.
You can see these pallete options using the following command:

```R
# use this to choose your palette name below
display.brewer.all()
```

So now let's generate a nicer heatmap, using one of the brewer palettes:

```R
# create nicer heatmap with brewer palette colors
displayClustersWithHeatmap(W = W,group = group,
                           col = brewer.pal(name = "Spectral",n=10))
```

### Extras: Comparing Clustering to Known Patient Labels

We'll start by creating some labels

```R
## confirmed diagnosis vector
truelabel = analysis.df$other_ordinal_dx

## Porto-defined diagnosis vector
truelabel1 = analysis.df$Porto_ordinal_dx
# grouped by true label=Other_ordinal
```

And then assigning the them to specific colors.

ID | True Label
--- | ---
1 |  Crohn's disease
2 | IBD-U
3 | Ulcerative Colitis

```R
col.vals_oo <- rep(NA,length(truelabel))
col.vals_oo[truelabel == 1] <- "tomato1"
col.vals_oo[truelabel == 2] <- "cornflowerblue"
col.vals_oo[truelabel == 3] <- "gold"

# grouped by true label1=Porto_ordinal
col.vals_po <- rep(NA,length(truelabel1))
col.vals_po[truelabel1 == 1] <- "tomato1"
col.vals_po[truelabel1 == 2] <- "cornflowerblue"
col.vals_po[truelabel1 == 3] <- "gold"

## color top of heatmap by label
displayClustersWithHeatmap(W = W, group = group,
                           col = brewer.pal(name = "Spectral",n=8),
                           ColSideColors = col.vals_oo)

displayClustersWithHeatmap(W = W, group = group,
                           col = brewer.pal(name = "Spectral",n=8),
                           ColSideColors = col.vals_po)
```

Let's test the relationship between our detected clusters and original labels.
```R
#### CREATE 2 BARPLOTS WITH TRUE AND PORTO LABELS
table(truelabel,group)
chisq.test(table(truelabel,group))

table(analysis.df$other_ordinal_dx,group)
chisq.test(table(analysis.df$other_ordinal_dx,group))

barplot(table(analysis.df$other_ordinal_dx,group),
        ylab= "Number of patients", xlab= "SNF groups",
        col=c("tomato1","cornflowerblue", "gold"),
        legend.text=c("Colonic CD","IBD-U",  "UC"))

table(analysis.df$Porto_ordinal_dx,group)
barplot(table(analysis.df$Porto_ordinal_dx,group),
        ylab= "Number of patients",
        xlab= "clustering groups",
        col=c("tomato1","cornflowerblue", "gold"),
        legend.text=c("Colonic CD", "IBDU", "UC"))
```

### Extras: Data-type specific heatmaps

Here we'll create a heatmap specific to one of the datasets; hist, rad, endosc, or otherhist.

```R
# Individual data types
## These similarity graphs have complementary information about clusters. one without the label at top and one with
names(aff.list)

displayClustersWithHeatmap(W = aff.list[["hist"]],
                           group = group,
                           col = brewer.pal(name = "Greens",n=8))

displayClustersWithHeatmap(W = aff.list[["rad"]],group = group, col = brewer.pal(name = "Greens",n=8))
displayClustersWithHeatmap(W = aff.list[["endosc"]],group = group, col = brewer.pal(name = "Greens",n=8))
displayClustersWithHeatmap(W = aff.list[["otherhist"]], group=group, col = brewer.pal(name = "Greens",n=8))
```

## Survival Analysis

To perform our survival analysis, we need to:

1. Create our survival outcome
1. Generate survival descriptive statistics
1. Generate KM curves
1. Test differences in survival by group
1. Plot survival by group
1. Fit Cox proportional hazards model
1. Run model diagnostics

Use `Surv` function to tie the time-to-event (age) with whether an event (cancer death) took placem

```R
analysis.df$survival_outcome <- Surv(time = analysis.df$age,
                                     event = analysis.df$can_death)
```

Every survival analysis uses an object generated by the `Surv` function as an outcome.
The `Surv` function can take the following two forms:

```
Surv(start time, end time, event = 0 or 1)
Surv(time to event, event = 0 or 1)
```

### Generate Survival Descriptive Statistics

Use `survfit` with argument of formula `survival outcome ~ 1`
The `~ 1` indicates we are only interested in survival values themselves.
`log-log` option allows you to fit the CI using the delta method.

```R
# summarizing survival without covariate
survival_fit <- survfit(analysis.df$survival_outcome ~ 1,
                         conf.type = "log-log")
```

### Generating a Basic KM Curve

We'll use the survival fit from the previous step:

```R
plot(survival_fit,col="blue4")
```

### Test Differences in Survival By Group

We'll add the SNF grouped column to our data frame, and create a factor variable from it (needed to be interpreted as a categorical variable).

```R
analysis.df$cluster = group
analysis.df$cluster.fac <- factor(analysis.df$cluster,
                                  levels = 1:max(analysis.df$cluster))
```

We can now test the difference in survival time using peto&peto modification on the
Gehan-Wilcoxon test, using the `survdiff` function


```R
## Cluster Assignment
survdiff(analysis.df$survival_outcome ~ analysis.df$cluster.fac, rho=1)
```

### Plot Survival by Group

We now use the `rms` package here with the `npsurv` and `survplot` functions to generate KM curves with confidence intervals.

FIrst we'll generate stratified KM values with `npsurv`, and then we'll use that to plot our KM curves by strata:

```R

# npsurv (non-parametric survival fit) function is a work around/replacement
#   for survfit since survfit no longer works with survplot which we want to use below

## looking at marginal survival difference by SNF generated cluster
survival.fit.by.snf.group <- npsurv(survival_outcome ~ cluster.fac,
                                    data = analysis.df)
```

We need to pass our `npsurv` result, and then colours for each curve ordered by factor order. Here 'cluster 1' is forestgreen and cluster 2 is darkorchid4. Extra colours are added in case there are more than 2 clusters.

```R
survplot(fit = survival.fit.by.snf.group,
         col=c('forestgreen','darkorchid4', 'dodgerblue','yellow','orange'),
         lwd=2.5,
         col.fill = sapply(c('forestgreen','darkorchid4','dodgerblue','yellow','orange'),
                           function(x){adjustcolor(x, alpha.f = 0.4)}),
         xlab="Age at death")
```

Note the use of `sapply` in the above to adjust colours and make confident interval colours a bit more transparent.


### Cox Proportional Hazards Model

We want to fit a model with our survival outcome and with cluster and age of {IBD} diagnosis as predictors:

Here we use the 'Breslow' method to estimate the baseline hazard function and remains unmodified if there are ties in survival time. The 'Efron' and 'Exact' methods can also be used. In this model, changing the method choice results in a negligible change of results.

```R
## Efron method
coxph.fit <- coxph(survival_outcome ~ cluster.fac + age_at_IBDdiagnosis,
                   data = analysis.df, method = "efron")
summary(coxph.fit)

## Exact method -- commented out because it may run slowly/require a lot of processing
# (coxph.fit <- coxph(survival_outcome ~
#                       cluster.fac + age_at_IBDdiagnosis,
#                     data = analysis.df, method = "exact"))
# summary(coxph.fit)

## Breslow method
coxph.fit <- coxph(survival_outcome ~
                      cluster.fac + age_at_IBDdiagnosis,
                    data = analysis.df,method = "breslow")
summary(coxph.fit)
```

Your 'summary(coxph.fit)' command should output many things including your model coefficients and significance, read it closely and compare the two methods, Efron and Breslow.

Often you will want to save these results outside of R in a table you can open in a spreadsheet program such as excel. To do this, you will need to
1. Extract the coefficients and confidence intervals
2. Combine them in a data frame/matrix
3. Export the results in a .csv file

Results in R are stored in a named list. To see the names for the elements of the list you'd like to extract, use the `str` (structure) function

```R
## Extracting results
str(summary(coxph.fit))
```

Note that we're looking at the structure of the summary object where we saw coefficients and confidence intervals. We don’t have to use the full name of the element in the list, only enough to uniquely identify it.

```R
(coxph.coefs <- summary(coxph.fit)$coef)
(coxph.confint <- summary(coxph.fit)$conf.int)
(coxph.results <- cbind(coxph.coefs,coxph.confint))
```

By surrounding the objects in `()` we can print it automatically, and see that we're going to have more columns than we are necessarily interested in.

To write only the columns we want to the .csv file, we index the column names in our ‘write.csv’ command.
To get the exact names of the columns we’re interested in, use the command ‘colnames’.
We can then select only those columns we're interested in, in our `write.csv` command.

```R
colnames(coxph.results)

write.csv(coxph.results[,c("coef","exp(coef)","se(coef)","Pr(>|z|)","lower .95","upper .95" )],
          file = "Cox-PH-model-results.csv")
```

### Run Model Diagnostics

Diagnostics to check:

1. Proportionality of hazards by covariate
1. Influential observations (outliers)
1. Linearity (functional form) of continuous covariates

#### Proportionality of hazards by covariate

Using cox.zph to test for covariate-specific and global proportional hazards as
well as plotting Schoenfeld residuals to check for non-proportional hazards --
significance implies non-proportionality.

```R
cox.zph(fit = coxph.fit)       # Statistic test for null, linear Schoenfeld residuals
par(mfrow=c(2,1))              # Change plot window to have 2 rows and 1 column
plot(cox.zph(fit = coxph.fit)) # Plot residuals with smoothing splot and ±2 SD error band
par(mfrow=c(1,1))              # Reset plot window.
```

> <tip-title>Recall Schoenfeld Residuals</tip-title>
> Schoenfeld residuals = observed – expected values of the covariates at each failure time. If proportionality of hazards holds, this value should not have a trend over time.
{: .tip}

#### Checking for influential observations (outliers)

Use 'dfbeta' residuals: values indicate change in coefficients when each individual is removed

```R
dfbeta <- residuals(coxph.fit, type = 'dfbeta') ## Dataframe of change in coefficients as each individual removed
par(mfrow=c(2,1))
for(j in 1:2){
  plot(dfbeta[,j],ylab=names(coef(coxph.fit))[j],
       pch=19,col='blue')
  abline(h=0,lty=2)
}
par(mfrow=c(1,1))
```

What we’re looking for here are positive/negative outliers (on the y-axis) since the presence/absence of these individuals in the model will change it drastically.

> <question-title></question-title>
> Do you see any outliers?
> > <solution-title></solution-title>
> > Note that no individuals change any coefficient by more than 0.04.
> > There are no extremely influential points.
> {: .solution}
{: .question}


Checking for linearity in the covariates using plots of martingale residuals
against the individual covariates is the next step.

> <tip-title>Only binary?</tip-title>
> This is not necessary for binary
> variables so we only check it in our age of initial diagnosis covariate
{: .tip}

Note that we only need to check linearity for one of our variables (age of initial diagnosis)

```R
martingale.resids <- residuals(coxph.fit, type = 'martingale')

# Plot
plot(y = martingale.resids,
     x = analysis.df$age_at_IBDdiagnosis,
     ylab = 'Residuals', xlab = 'Age of Initial Diag', pch= 19, col = 'blue')
abline(h=0,lty=2,col='red')
lines(lowess(x = analysis.df$age_at_IBDdiagnosis,
             y = martingale.resids, iter = 0))
```

> <tip-title>Graph interpretation</tip-title>
> Note the way the local linear regression line is on top of the straight fitted line – linearity holds!
{: .tip}

> <tip-title>Recall: Martingale Residuals</tip-title>
> Martingale residuals measure observed (1 or 0) vs expected (0+) event counts at each individual’s failure/censor time. Here we make sure there’s no trend over the covariate’s values.
>
> Note: Martingale residuals occur in the range (-∞, 1]
{: .tip}


## Extra: Graphs investigating different the relationship of indivdiual features relative to SNF features

```R
# check top features by barplot
# non bloody diarrhea
table(analysis.df$non_bloody_diarrhea,group)
barplot(table(analysis.df$non_bloody_diarrhea,group),
        ylab= "Number of patients",
        xlab= "clustering groups",
        main= "Non bloody diarrhea at diagnosis",
        col=c("darkred", "darkgreen"),
        legend.text=c("bloody diarrhea", "non-bloody diarrhea"))



table(analysis.df$reverse_gradient,group)
barplot(table(analysis.df$reverse_gradient,group),
        ylab= "Number of patients",
        xlab= "clustering groups",
        main= "Reverse gradient of inflammation (proximal>distal)",
        col=c("darkred", "darkgreen"),
        legend.text=c("No", "Yes"))


table(analysis.df$gastritis,group)
barplot(table(analysis.df$gastritis,group),
        ylab= "Number of patients",
        xlab= "clustering groups",
        main= "Presence of Gastritis at diagnosis",
        col=c("darkblue", "darkgreen"),
        legend.text=c("no gastritis", "gastritis"))


table(analysis.df$rectal_sparing,group)
barplot(table(analysis.df$rectal_sparing,group),
        ylab= "Number of patients",
        xlab= "clustering groups",
        main= "Rectal sparing",
        col=c("darkred", "darkgreen"),
        legend.text=c("No", "Yes"))


table(analysis.df$focal_active_colitis,group)
barplot(table(analysis.df$focal_active_colitis,group),
        ylab= "Number of patients",
        xlab= "clustering groups",
        main= "Presence of focal Active Colitis",
        col=c("darkred", "darkgreen"),
        legend.text=c("No", "Yes"))

table(analysis.df$activity,group)
barplot(table(analysis.df$activity,group),
        ylab= "Number of patients",
        xlab= "clustering groups",
        main= "Neutrophil Activity",
        col=c("darkred", "darkgreen", "gold", "darkblue", "purple"),
        legend.text=c("none", "focal", "mild", "moderate", "severe"))

table(analysis.df$deep_ulcer_SB,group)
barplot(table(analysis.df$deep_ulcer_SB,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = 'Small bowel Ulceration/Cobblestoning',
        col = c("darkblue", "darkgreen"),
        legend.text=c("no SB", "SB ulceration"))

table(analysis.df$crypt_distortion,group)
barplot(table(analysis.df$crypt_distortion,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = 'Crypt Distortion',
        col = c("darkblue", "darkgreen", "gold"),
        legend.text=c("no distortion", "patchy", "diffuse"))

table(analysis.df$pattern_involvement,group)
barplot(table(analysis.df$pattern_involvement,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = 'Disease worse distally at diagnosis (histological)',
        col= c("darkred", "darkblue"),
        legend.text = c("Yes","No"))


table(analysis.df$chronic_inflammation,group)
barplot(table(analysis.df$chronic_inflammation,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = "Chronic Inflammation",
        col= c("darkblue", "pink"),
        legend.text = c("patchy", "diffuse"))



table(analysis.df$FEG,group)
barplot(table(analysis.df$FEG,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = "Focal Enhanced Gastritis",
        col= c("darkblue", "darkgreen"),
        legend.text = c("no FEG", "FEG"))

table(analysis.df$rectal_sparing,group)
barplot(table(analysis.df$rectal_sparing,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = "Complete Rectal Sparing",
        col= c("darkblue", "darkgreen"),
        legend.text = c("No", "Yes"))

table(analysis.df$focal_chronic_duodenitis,group)
barplot(table(analysis.df$focal_chronic_duodenitis,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = "Focal Chronic Duodenitis",
        col= c("darkblue", "darkgreen"),
        legend.text = c("No", "Yes"))

table(analysis.df$duodenitis,group)
barplot(table(analysis.df$duodenitis,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = "Duodenitis",
        col= c("darkgreen", "darkblue", "purple", "darkred"),
        legend.text = c("None", "Active" ,"Inactive", "Both inactive/active"))

table(analysis.df$basal_plasma_cells,group)
barplot(table(analysis.df$basal_plasma_cells,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = "Basal Plasma Cell",
        col= c("darkgreen", "darkblue", "gold"),
        legend.text = c("none", "few", "plasmacytosis"))

table(analysis.df$Granuloma,group)
barplot(table(analysis.df$Granuloma,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = "Granuloma not associated with crypt rupture",
        col= c("darkgreen", "orange"),
        legend.text = c("No", "Yes"))

table(analysis.df$skip_lesion,group)
barplot(table(analysis.df$skip_lesion,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = "Skip lesions (macroscopic and microscopic)",
        col= c("darkblue", "pink"),
        legend.text = c("No", "Yes"))

table(analysis.df$X5_small_ulcers_colon,group)
barplot(table(analysis.df$X5_small_ulcers_colon,group),
        ylab = "Number of patients",
        xlab = "clustering groups",
        main = "More than 5 small colonic ulcers",
        col= c("darkgreen", "orange"),
        legend.text = c("No", "Yes"))
```
