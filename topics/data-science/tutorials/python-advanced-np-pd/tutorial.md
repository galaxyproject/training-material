---
layout: tutorial_hands_on

title: Advanced Python
level: Intermediate
zenodo_link:  https://zenodo.org/record/3477564
requirements:
 -
   type: internal
   topic_name: data-science
   tutorials:
     - python-basics
follow_up_training:  []

questions:
- How can I analyze data using Python in Galaxy platform?

objectives:
- Use the scientific libraries pandas and numpy to explore tabular datasets

time_estimation:  3H
key_points: []
subtopic: python
contributors:
  - mcmaniou
  - fpsom
  
notebook:
  language: python
---

# Introduction

## Comment 
This tutorial is significantly based on the Carpentries courses [Programming with Python](https://swcarpentry.github.io/python-novice-inflammation/) and [Plotting and Programming in Python](https://swcarpentry.github.io/python-novice-gapminder/). 

## Overview
In this lesson, we will be using Python 3 with some of its most popular scientific libraries. This tutorial assumes that the reader is familiar with the fundamentals of the Python programming language, as well as, how to run Python programs using Galaxy. Otherwise, it is advised to follow the "Introduction to Python" tutorial available in the same platform. We will be using JupyterNotebook, a Python interpreter that comes with everything we need for the lesson. Please note:  JupyterNotebook is only currently available on the [usegalaxy.eu](https://usegalaxy.eu/) and [usegalaxy.org](https://usegalaxy.org/) sites. 


# Analyze data using numpy
NumPy is a python library and it stands for Numerical Python. In general, you should use this library when you want to perform operations and manipulate numerical data, especially if you have matrices or arrays. To tell Python that we’d like to start using NumPy, we need to import it:

```python
import numpy as np
```
A Numpy array contains one or more elements of the same type. To examine the basic functions of the library, we will create an array of random data. These data will correspond to arthritis patients’ inflammation. The rows are the individual patients, and the columns are their daily inflammation measurements. We will use the `random.randint()` function. It has 4 arguments as inputs `randint(low, high=None, size=None, dtype=int)`. `low` nad `high` specify the limits of the random number generator. `size` determines the shape of the array and it can be an integer or a tuple.

```python
np.random.seed(2021)  #create reproducible work
random_data = np.random.randint(1, 25, size=(50,70))
```
If we want to check the data have been loaded, we can print the variable’s value:

```python
print(random_data)
```
```output
[[21 22  1 ...  8 21 12]
 [20 19 15 ...  5 10  7]
 [ 7 24 14 ... 24 19 12]
 ...
 [22 16  3 ...  2 18 21]
 [22 14 23 ...  7 13  4]
 [ 7  7  8 ... 15  8  9]]
```
Now that the data are in memory, we can manipulate them. First, let’s ask what type of thing data refers to:

```python
print(type(random_data))
```
```Output
<class 'numpy.ndarray'>
```

The output tells us that data currently refers to an N-dimensional array, the functionality for which is provided by the NumPy library. These data correspond to arthritis patients’ inflammation. The rows are the individual patients, and the columns are their daily inflammation measurements.

The `type` function will only tell you that a variable is a NumPy array but won’t tell you the type of thing inside the array. We can find out the type of the data contained in the NumPy array.

```python
print(random_data.dtype)
```
```Output
int64
```
This tells us that the NumPy array’s elements are integer numbers.

With the following command, we can see the array’s shape:

```python
print(random_data.shape)
```
```Output
(50, 70)
```
The output tells us that the data array variable contains 50 rows and 70 columns. When we created the variable `random_data` to store our arthritis data, we did not only create the array; we also created information about the array, called members or attributes. This extra information describes `random_data` in the same way an adjective describes a noun. `random_data.shape` is an attribute of `random_data` which describes the dimensions of `random_data`.

If we want to get a single number from the array, we must provide an index in square brackets after the variable name, just as we do in math when referring to an element of a matrix. Our data has two dimensions, so we will need to use two indices to refer to one specific value:

```python
print('first value in data:', random_data[0, 0])
```
```output
first value in data: 21
```

```python
print('middle value in data:', random_data[25, 35])
```
```output
middle value in data: 11
```

The expression random_data[25, 35] accesses the element at row 25, column 35. While this expression may not surprise you, random_data[0, 0] might. Programming languages like Fortran, MATLAB and R start counting at 1 because that’s what human beings have done for thousands of years. Languages in the C family (including C++, Java, Perl, and Python) count from 0 because it represents an offset from the first value in the array (the second value is offset by one index from the first value). As a result, if we have an M×N array in Python, its indices go from 0 to M-1 on the first axis and 0 to N-1 on the second. 

Slicing data
An index like [25, 35] selects a single element of an array, but we can select whole sections as well, using slicing the same way as previously with the strings. For example, we can select the first ten days (columns) of values for the first four patients (rows) like this:

```python
print(random_data[0:4, 0:10])
```
```output
[[21 22  1 14 23 13 22 13 23  7]
 [20 19 15 19 22 19  1 21 24  4]
 [ 7 24 14 11  5  6  3  4  3 17]
 [24  9  4 17 10  4 11  4 13 20]]
```

We don’t have to include the upper and lower bound on the slice. If we don’t include the lower bound, Python uses 0 by default; if we don’t include the upper, the slice runs to the end of the axis, and if we don’t include either (i.e., if we use ‘:’ on its own), the slice includes everything:

```python
small = random_data[:3, 36:]
print('small is:')
print(small)
```
```output
small is:
[[22 23 19 11  5 14 11 20 19 16  9 11 13 22 12 20  2 15 18 19  8 22 20 13
  24 20 10  5 10  8  8  8 21 12]
 [ 7 15 12 10 20  2 12 20 16  7 18 12  1  6 12 17 24  1 18  4  7 13 12 21
   8 10  5  2  8 18 16  5 10  7]
 [23 23  8  5 20 14  3  9 17 13 22  9  2 11 14 23  7 20  8 22  5 21 20  2
  11  4 24 19  2 14 22 24 19 12]]
```

The above example selects rows 0 through 2 and columns 36 through to the end of the array.

## Process the data
NumPy has several useful functions that take an array as input to perform operations on its values. If we want to find the average inflammation for all patients on all days, for example, we can ask NumPy to compute random_data’s mean value:

```python
print(np.mean(random_data))
```
```output
12.567142857142857
```

Let’s use three other NumPy functions to get some descriptive values about the dataset. We’ll also use multiple assignment, a convenient Python feature that will enable us to do this all in one line.

```python
maxval, minval, stdval = np.max(random_data), np.min(random_data), np.std(random_data)

print('maximum inflammation:', maxval)
print('minimum inflammation:', minval)
print('standard deviation:', stdval)
```
```output
maximum inflammation: 24
minimum inflammation: 1
standard deviation: 6.966639309258527
```

How did we know what functions NumPy has and how to use them? If you are working in IPython or in a Jupyter Notebook, there is an easy way to find out. If you type the name of something followed by a dot, then you can use tab completion (e.g. type `np.` and then press Tab) to see a list of all functions and attributes that you can use. After selecting one, you can also add a question mark (e.g. `np.cumprod?`), and IPython will return an explanation of the method! This is the same as doing `help(np.cumprod)`.

When analyzing data, though, we often want to look at variations in statistical values, such as the maximum inflammation per patient or the average inflammation per day. One way to do this is to create a new temporary array of the data we want, then ask it to do the calculation:

```python
patient_0 = random_data[0, :] # 0 on the first axis (rows), everything on the second (columns)
print('maximum inflammation for patient 0:', np.max(patient_0))
```
```output
maximum inflammation for patient 0: 24
```

What if we need the maximum inflammation for each patient over all days (as in the next diagram on the left) or the average for each day (as in the diagram on the right)? As the diagram below shows, we want to perform the operation across an axis:

![Figure 8_Operations_across_axis](../../images/python-advanced-np-pd/Figure8_Operations_across_axis.png)

To support this functionality, most array functions allow us to specify the axis we want to work on. If we ask for the average across axis 0 (rows in our 2D example), we get:

```python
print(np.mean(random_data, axis=0))
```
```output
[12.9  13.4  12.66 12.54 12.76 12.22 10.84 12.32 12.56 13.5  11.36 14.62
 12.18 13.02 13.3  10.6  11.36 10.26 10.22 13.14 14.12 12.5  12.3  12.78
 12.4  12.08 13.42 12.64 11.5  12.08 14.28 13.2  12.48 13.4  12.54 13.84
 12.88 13.74 10.88 11.38 11.7  11.7  12.62 13.8  14.02 12.92 12.98 12.5
 11.52 13.96 11.98 14.5  12.16 11.04 13.96 12.48 11.98 12.32 12.2  13.88
 10.96 13.44 12.26 12.98 11.46 12.52 13.9  13.   12.4  12.36]
```

As a quick check, we can ask this array what its shape is:

```python
print(np.mean(random_data, axis=0).shape)
```
```Output
(70,)
```

The expression (70,) tells us we have an N×1 vector, so this is the average inflammation per day for all patients. If we average across axis 1 (columns in our 2D example), we get the average inflammation per patient across all days.:

```python
print(np.mean(random_data, axis=1))
```
```output
[13.28571429 12.37142857 13.31428571 12.21428571 14.24285714 13.58571429
 11.55714286 13.02857143 12.44285714 12.37142857 13.15714286 12.32857143
 11.37142857 11.94285714 12.37142857 11.81428571 12.48571429 12.47142857
 13.14285714 13.05714286 12.4        13.7        13.28571429 13.01428571
 12.8        12.05714286 11.57142857 12.25714286 12.62857143 12.12857143
 13.94285714 11.88571429 12.02857143 12.98571429 12.15714286 12.24285714
 12.48571429 13.87142857 12.27142857 13.77142857 12.04285714 12.34285714
 13.58571429 13.04285714 12.28571429 11.35714286 10.3        12.72857143
 11.81428571 12.81428571]
```

### Stacking arrays
Arrays can be concatenated and stacked on top of one another, using NumPy’s `vstack` and `hstack` functions for vertical and horizontal stacking, respectively.

```python
import numpy as np

A = np.array([[1,2,3], [4,5,6], [7, 8, 9]])
print('A = ')
print(A)

B = np.hstack([A, A])
print('B = ')
print(B)

C = np.vstack([A, A])
print('C = ')
print(C)
```
```Output
A =
[[1 2 3]
 [4 5 6]
 [7 8 9]]
B =
[[1 2 3 1 2 3]
 [4 5 6 4 5 6]
 [7 8 9 7 8 9]]
C =
[[1 2 3]
 [4 5 6]
 [7 8 9]
 [1 2 3]
 [4 5 6]
 [7 8 9]]
```


### Remove NaN values

Sometimes there are missing values in an array, that could make it difficult to perform operations on it. To remove the `NaN` you must first find their indexes and then replace them. The following example replaces them with `0`.

```python
a = array([[1, 2, 3], [0, 3, NaN]])
print(a)
a[np.isnan(a)] = 0
print(a)
```
```Output
[[ 1.  2.  3.]
 [ 0.  3. nan]]
 
[[1. 2. 3.]
 [0. 3. 0.]]
```


 ### Exercises

1. Write some additional code that slices the first and last columns of A, and stacks them into a 3x2 array. Make sure to print the results to verify your solution.

2. Given the followind array `A`, keep only the elements that are lower that `0.05`.

   ```python
   A = np.array([0.81, 0.025, 0.15, 0.67, 0.01])
   ```


 ### Solutions
1. A ‘gotcha’ with array indexing is that singleton dimensions are dropped by default. That means `A[:, 0]` is a one dimensional array, which won’t stack as desired. To preserve singleton dimensions, the index itself can be a slice or array. For example, `A[:, :1]` returns a two dimensional array with one singleton dimension (i.e. a column vector).

   ```python
   D = np.hstack((A[:, :1], A[:, -1:]))
   print('D = ')
   print(D)
   ```
   ```Output
   D =
   [[1 3]
    [4 6]
    [7 9]]
   ```

2. ```python
   A = A[A<0.05]
   ```

# Use pandas to work with dataframes

Pandas is a widely-used Python library for statistics, particularly on tabular data. If you are familiar with R dataframes, then this is the library that integrates this functionality. A dataframe is a 2-dimensional table with indexes and column names. The indexes indicate the difference in rows, while the column names indicate the difference in columns. You will see later that these two features are useful when you’re manipulating your data. Each column can contain different data types.

Load it with import pandas as `pd`. The alias `pd` is commonly used for pandas.

```python
import pandas as pd
```

There are many ways to create a pandas dataframe. For example you can use a numpy array as input.

```python
data = np.array([['','Col1','Col2'],
                ['Row1',1,2],
                ['Row2',3,4]])
                
print(pd.DataFrame(data=data[1:,1:],
                  index=data[1:,0],
                  columns=data[0,1:]))
```
```output
     Col1 Col2
Row1    1    2
Row2    3    4
```


For the purposes of this tutorial, we will use a file with the annotated differentially expressed genes that was produced in the [Reference-based RNA-Seq data analysis](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html) tutorial.

We can read a tabular file with `pd.read_csv`. The first argument is the filepath of the file to be read. The `sep` argument refers to the symbol used to separate the data into different columns. You can check the rest of the arguments using the `help()` function.

```python
data = pd.read_csv("https://zenodo.org/record/3477564/files/annotatedDEgenes.tabular", sep = "\t")
print(data)
```
```output
          GeneID     Base mean  log2(FC)    StdErr  Wald-Stats        P-value  \
0    FBgn0039155   1086.974295 -4.148450  0.134949  -30.740913  1.617357e-207   
1    FBgn0003360   6409.577128 -2.999777  0.104345  -28.748637  9.419922e-182   
2    FBgn0026562  65114.840564 -2.380164  0.084327  -28.225437  2.850430e-175   
3    FBgn0025111   2192.322369  2.699939  0.097945   27.565978  2.847517e-167   
4    FBgn0029167   5430.067277 -2.105062  0.092547  -22.745964  1.573284e-114   
..           ...           ...       ...       ...         ...            ...   
125  FBgn0035710     26.161771  1.048979  0.232922    4.503559   6.682473e-06   
126  FBgn0035523     70.998197  1.004819  0.223763    4.490561   7.103599e-06   
127  FBgn0038261     44.270577  1.006264  0.224124    4.489756   7.130487e-06   
128  FBgn0039178     23.550056  1.040917  0.232626    4.474631   7.654328e-06   
129  FBgn0034636     24.770519 -1.028531  0.232168   -4.430118   9.418135e-06   

             P-adj Chromosome     Start       End Strand         Feature  \
0    1.387207e-203      chr3R  24141394  24147490      +  protein_coding   
1    4.039734e-178       chrX  10780892  10786958      -  protein_coding   
2    8.149380e-172      chr3R  26869237  26871995      -  protein_coding   
3    6.105789e-164       chrX  10778953  10786907      -  protein_coding   
4    2.698811e-111      chr3L  13846053  13860001      +  protein_coding   
..             ...        ...       ...       ...    ...             ...   
125   1.436480e-04      chr3L   6689326   6703521      -  protein_coding   
126   1.523189e-04      chr3L   4277961   4281585      +  protein_coding   
127   1.525142e-04      chr3R  14798985  14801163      +  protein_coding   
128   1.629061e-04      chr3R  24283954  24288617      +  protein_coding   
129   1.951192e-04      chr2R  21560245  21576035      -  protein_coding   

       Gene name  
0           Kal1  
1           sesB  
2    BM-40-SPARC  
3           Ant2  
4            Hml  
..           ...  
125       SP1173  
126       CG1311  
127      CG14856  
128       CG6356  
129      CG10440  

[130 rows x 13 columns]
```
The columns in a dataframe are the observed variables, and the rows are the observations. Pandas uses backslash `\` to show wrapped lines when output is too wide to fit the screen.

## Explore the data

You can use `index_col` to specify that a column’s values should be used as row headings.

By default row indexes are numbers, but we could use a column of the data. To pass the name of the column to `read_csv`, you can use its `index_col` parameter. Be careful though, because the row indexes must be unique for each row.

```python
data = pd.read_csv("https://zenodo.org/record/3477564/files/annotatedDEgenes.tabular", sep = "\t", index_col = 'GeneID')
print(data)
```
```output
                Base mean  log2(FC)    StdErr  Wald-Stats        P-value  \
GeneID                                                                     
FBgn0039155   1086.974295 -4.148450  0.134949  -30.740913  1.617357e-207   
FBgn0003360   6409.577128 -2.999777  0.104345  -28.748637  9.419922e-182   
FBgn0026562  65114.840564 -2.380164  0.084327  -28.225437  2.850430e-175   
FBgn0025111   2192.322369  2.699939  0.097945   27.565978  2.847517e-167   
FBgn0029167   5430.067277 -2.105062  0.092547  -22.745964  1.573284e-114   
...                   ...       ...       ...         ...            ...   
FBgn0035710     26.161771  1.048979  0.232922    4.503559   6.682473e-06   
FBgn0035523     70.998197  1.004819  0.223763    4.490561   7.103599e-06   
FBgn0038261     44.270577  1.006264  0.224124    4.489756   7.130487e-06   
FBgn0039178     23.550056  1.040917  0.232626    4.474631   7.654328e-06   
FBgn0034636     24.770519 -1.028531  0.232168   -4.430118   9.418135e-06   

                     P-adj Chromosome     Start       End Strand  \
GeneID                                                             
FBgn0039155  1.387207e-203      chr3R  24141394  24147490      +   
FBgn0003360  4.039734e-178       chrX  10780892  10786958      -   
FBgn0026562  8.149380e-172      chr3R  26869237  26871995      -   
FBgn0025111  6.105789e-164       chrX  10778953  10786907      -   
FBgn0029167  2.698811e-111      chr3L  13846053  13860001      +   
...                    ...        ...       ...       ...    ...   
FBgn0035710   1.436480e-04      chr3L   6689326   6703521      -   
FBgn0035523   1.523189e-04      chr3L   4277961   4281585      +   
FBgn0038261   1.525142e-04      chr3R  14798985  14801163      +   
FBgn0039178   1.629061e-04      chr3R  24283954  24288617      +   
FBgn0034636   1.951192e-04      chr2R  21560245  21576035      -   

                    Feature    Gene name  
GeneID                                    
FBgn0039155  protein_coding         Kal1  
FBgn0003360  protein_coding         sesB  
FBgn0026562  protein_coding  BM-40-SPARC  
FBgn0025111  protein_coding         Ant2  
FBgn0029167  protein_coding          Hml  
...                     ...          ...  
FBgn0035710  protein_coding       SP1173  
FBgn0035523  protein_coding       CG1311  
FBgn0038261  protein_coding      CG14856  
FBgn0039178  protein_coding       CG6356  
FBgn0034636  protein_coding      CG10440  

[130 rows x 12 columns]
```

You can use the `DataFrame.info()` method to find out more about a dataframe.

```python
data.info()
```
```output
<class 'pandas.core.frame.DataFrame'>
Index: 130 entries, FBgn0039155 to FBgn0034636
Data columns (total 12 columns):
 #   Column      Non-Null Count  Dtype  
---  ------      --------------  -----  
 0   Base mean   130 non-null    float64
 1   log2(FC)    130 non-null    float64
 2   StdErr      130 non-null    float64
 3   Wald-Stats  130 non-null    float64
 4   P-value     130 non-null    float64
 5   P-adj       130 non-null    float64
 6   Chromosome  130 non-null    object 
 7   Start       130 non-null    int64  
 8   End         130 non-null    int64  
 9   Strand      130 non-null    object 
 10  Feature     130 non-null    object 
 11  Gene name   130 non-null    object 
dtypes: float64(6), int64(2), object(4)
memory usage: 13.2+ KB
```


We learn that this is a DataFrame. It consists of 130 rows and 12 columns. None of the columns contains any missing values. 6 columns contain 64-bit floating point `float64` values, 2 contain 64-bit integer `int64` values and 4 contain character `object` values. It uses 13.2KB of memory.

The `DataFrame.columns` variable stores information about the dataframe’s columns.

Note that this is an attribute, not a method. (It doesn’t have parentheses.) Called a member variable, or just member.

```python
print(data.columns)
```
```output
Index(['Base mean', 'log2(FC)', 'StdErr', 'Wald-Stats', 'P-value', 
	   'P-adj', 'Chromosome', 'Start', 'End', 'Strand', 'Feature', 'Gene name'],
      dtype='object')
```


You could use `DataFrame.T` to transpose a dataframe. The `Transpose` (written `.T`) doesn’t copy the data, just changes the program’s view of it. Like columns, it is a member variable.

```python
print(data.T)
```
```output
GeneID         FBgn0039155     FBgn0003360     FBgn0026562     FBgn0025111  \
Base mean      1086.974295     6409.577128    65114.840564     2192.322369   
log2(FC)          -4.14845       -2.999777       -2.380164        2.699939   
StdErr            0.134949        0.104345        0.084327        0.097945   
Wald-Stats      -30.740913      -28.748637      -28.225437       27.565978   
P-value                0.0             0.0             0.0             0.0   
P-adj                  0.0             0.0             0.0             0.0   
Chromosome           chr3R            chrX           chr3R            chrX   
Start             24141394        10780892        26869237        10778953   
End               24147490        10786958        26871995        10786907   
Strand                   +               -               -               -   
Feature     protein_coding  protein_coding  protein_coding  protein_coding   
Gene name             Kal1            sesB     BM-40-SPARC            Ant2   

GeneID         FBgn0029167     FBgn0039827     FBgn0035085     FBgn0034736  \
Base mean      5430.067277      390.901782      928.263812      330.383023   
log2(FC)         -2.105062       -3.503013       -2.414074       -3.018179   
StdErr            0.092547         0.16003        0.115185        0.158154   
Wald-Stats      -22.745964      -21.889756      -20.958204      -19.083791   
P-value                0.0             0.0             0.0             0.0   
P-adj                  0.0             0.0             0.0             0.0   
Chromosome           chr3L           chr3R           chr2R           chr2R   
Start             13846053        31196915        24945138        22550093   
End               13860001        31203722        24946636        22552113   
Strand                   +               +               +               +   
Feature     protein_coding  protein_coding  protein_coding  protein_coding   
Gene name              Hml          CG1544          CG3770          CG6018   

GeneID     FBgn0264475     FBgn0000071  ...     FBgn0264343     FBgn0038237  \
Base mean   955.454454      468.057926  ...       57.156629       43.409588   
log2(FC)     -2.334486        2.360017  ...        1.036055       -1.105893   
StdErr         0.12423        0.135644  ...        0.216472        0.232871   
Wald-Stats  -18.791643       17.398617  ...        4.786091       -4.748945   
P-value            0.0             0.0  ...        0.000002        0.000002   
P-adj              0.0             0.0  ...        0.000042        0.000049   
Chromosome       chr3L           chr3R  ...           chr2L           chr3R   
Start           820758         6762592  ...         7227733        14513900   
End             821512         6765261  ...         7228641        14558304   
Strand               +               +  ...               +               -   
Feature        lincRNA  protein_coding  ...  protein_coding  protein_coding   
Gene name      CR43883             Ama  ...         CG43799            Pde6   

GeneID         FBgn0020376     FBgn0028939     FBgn0036560     FBgn0035710  \
Base mean         61.54154       31.587685       27.714241       26.161771   
log2(FC)         -1.038158       -1.091024           1.089        1.048979   
StdErr            0.219867        0.232324        0.232753        0.232922   
Wald-Stats       -4.721761       -4.696136        4.678786        4.503559   
P-value           0.000002        0.000003        0.000003        0.000007   
P-adj             0.000055        0.000062        0.000067        0.000144   
Chromosome           chr2L           chr2L           chr3L           chr3L   
Start              4120342        13980125        16035484         6689326   
End                4121627        13983269        16037227         6703521   
Strand                   +               +               +               -   
Feature     protein_coding  protein_coding  protein_coding  protein_coding   
Gene name          Sr-CIII           NimC2          CG5895          SP1173   

GeneID         FBgn0035523     FBgn0038261     FBgn0039178     FBgn0034636  
Base mean        70.998197       44.270577       23.550056       24.770519  
log2(FC)          1.004819        1.006264        1.040917       -1.028531  
StdErr            0.223763        0.224124        0.232626        0.232168  
Wald-Stats        4.490561        4.489756        4.474631       -4.430118  
P-value           0.000007        0.000007        0.000008        0.000009  
P-adj             0.000152        0.000153        0.000163        0.000195  
Chromosome           chr3L           chr3R           chr3R           chr2R  
Start              4277961        14798985        24283954        21560245  
End                4281585        14801163        24288617        21576035  
Strand                   +               +               +               -  
Feature     protein_coding  protein_coding  protein_coding  protein_coding  
Gene name           CG1311         CG14856          CG6356         CG10440  

[12 rows x 130 columns]
```


You can use `DataFrame.describe()` to get summary statistics about the data. `DataFrame.describe()` returns the summary statistics of only the columns that have numerical data.  All other columns are ignored, unless you use the argument `include='all'`. Depending on the data type of each column, the statistics that can't be calculated are replaced with  the value `NaN`.

```python
print(data.describe(include='all'))
```
```output
           Base mean    log2(FC)      StdErr  Wald-Stats        P-value  \
count     130.000000  130.000000  130.000000  130.000000   1.300000e+02   
unique           NaN         NaN         NaN         NaN            NaN   
top              NaN         NaN         NaN         NaN            NaN   
freq             NaN         NaN         NaN         NaN            NaN   
mean     1911.266780   -0.207426    0.164324   -1.901479   4.200552e-07   
std      6888.074171    1.612615    0.042358   11.037352   1.531253e-06   
min        19.150759   -4.148450    0.084327  -30.740913  1.617357e-207   
25%       100.286324   -1.336283    0.128488  -10.004594   4.952408e-31   
50%       237.986359   -1.027150    0.163704   -4.981688   2.109219e-19   
75%       948.656793    1.220304    0.198655    7.692073   6.731953e-12   
max     65114.840564    2.699939    0.232922   27.565978   9.418135e-06   

                P-adj Chromosome         Start           End Strand  \
count    1.300000e+02        130  1.300000e+02  1.300000e+02    130   
unique            NaN          5           NaN           NaN      2   
top               NaN      chr3R           NaN           NaN      +   
freq              NaN         32           NaN           NaN     72   
mean     9.321175e-06        NaN  1.343684e+07  1.344644e+07    NaN   
std      3.278527e-05        NaN  7.970664e+06  7.973420e+06    NaN   
min     1.387207e-203        NaN  1.274480e+05  1.403400e+05    NaN   
25%      1.275950e-28        NaN  7.277516e+06  7.279063e+06    NaN   
50%      2.563262e-17        NaN  1.316155e+07  1.316625e+07    NaN   
75%      3.706020e-10        NaN  1.925043e+07  1.928464e+07    NaN   
max      1.951192e-04        NaN  3.119692e+07  3.120372e+07    NaN   

               Feature Gene name  
count              130       130  
unique               3       130  
top     protein_coding      Sesn  
freq               126         1  
mean               NaN       NaN  
std                NaN       NaN  
min                NaN       NaN  
25%                NaN       NaN  
50%                NaN       NaN  
75%                NaN       NaN  
max                NaN       NaN  
```



### Exercises

1. After reading the data, use `help(data.head)` and `help(data.tail)` to find out what `DataFrame.head` and `DataFrame.tail` do.
   a. What method call will display the first three rows of the data?
   b. What method call will display the last three columns of this data? (Hint: you may need to change your view of the data.)
   
2. As well as the `read_csv` function for reading data from a file, Pandas provides a `to_csv` function to write dataframes to files. Applying what you’ve learned about reading from files, write one of your dataframes to a file called `processed.csv`. You can use `help` to get information on how to use `to_csv`.

### Solutions

1. a. We can check out the first five rows of the data by executing `data.head()` (allowing us to view the head of the DataFrame). We can specify the number of rows we wish to see by specifying the parameter `n` in our call to `data.head()`. To view the first three rows, execute:

    ```python
    data.head(n=3)
    ```

|  | Base mean | log2(FC) | StdErr | Wald-Stats | P-value | P-adj | Chromosome | Start | End | Strand | Feature | Gene name | GeneID |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
| FBgn0039155 | 1086.974295 | -4.148450 | 0.134949 | -30.740913 | 1.617357e-207 | 1.387207e-203 | chr3R | 24141394 | 24147490 | + | protein_coding | Kal1 |
| FBgn0003360 | 6409.577128 | -2.999777 | 0.104345 | -28.748637 | 9.419922e-182 | 4.039734e-178 | chrX | 10780892 | 10786958 | - | protein_coding | sesB |
| FBgn0026562 | 65114.840564 | -2.380164 | 0.084327 | -28.225437 | 2.850430e-175 | 8.149380e-172 | chr3R | 26869237 | 26871995 | - | protein_coding | BM-40-SPARC |

​     b. To check out the last three rows, we would use the command, `data.tail(n=3)`, analogous to `head()` used above. However, here we want to look 	at the last three columns so we need to change our view and then use `tail()`. To do so, we create a new DataFrame in which rows and columns are 	switched:


    data_flipped = data.T

​	We can then view the last three columns of the data by viewing the last three rows of data_flipped:

    data_flipped.tail(n=3)

| GeneID | FBgn0039155 | FBgn0003360 | FBgn0026562 | FBgn0025111 | FBgn0029167 | FBgn0039827 | FBgn0035085 | FBgn0034736 | FBgn0264475 | FBgn0000071 | ... | FBgn0264343 | FBgn0038237 | FBgn0020376 | FBgn0028939 | FBgn0036560 | FBgn0035710 | FBgn0035523 | FBgn0038261 | FBgn0039178 | FBgn0034636 |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
| Strand | + | - | - | - | + | + | + | + | + | + | ... | + | - | + | + | + | - | + | + | + | - | 
| Feature | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding | lincRNA | protein_coding | ... | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding | protein_coding |
| Gene name | Kal1 | sesB | BM-40-SPARC | Ant2 | Hml | CG1544 | CG3770 | CG6018 | CR43883 | Ama | ... | CG43799 | Pde6 | Sr-CIII | NimC2 | CG5895 | SP1173 | CG1311 | CG14856 | CG6356 | CG10440 |

* Note about Pandas DataFrames/Series

A [DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) is a collection of [Series](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.html); The DataFrame is the way Pandas represents a table, and Series is the data-structure Pandas use to represent a column.

Pandas is built on top of the Numpy library, which in practice means that most of the methods defined for Numpy Arrays apply to Pandas Series/DataFrames.

What makes Pandas so attractive is the powerful interface to access individual records of the table, proper handling of missing values, and relational-databases operations between DataFrames.

## Select data
To access a value at the position `[i,j]` of a DataFrame, we have two options, depending on what is the meaning of i in use. Remember that a DataFrame provides an index as a way to identify the rows of the table; a row, then, has a position inside the table as well as a label, which uniquely identifies its entry in the DataFrame.

You can use `DataFrame.iloc[..., ...]` to select values by their (entry) position and basically specify location by numerical index analogously to 2D version of character selection in strings.

```python
print(data.iloc[0, 0])
```
You can also use `DataFrame.loc[..., ...]` to select values by their (entry) label and basically specify location by row name analogously to 2D version of dictionary keys.

```python
print(data.loc["FBgn0039155", "Base mean"])
```
```output
1086.97429520489
```
You can use Python's usual slicing notation, to select all or a subset of rows and/or columns. For example, the following code selects all the columns of the row `"FBgn0039155"`.

```python
print(data.loc["FBgn0039155", :])
```
```output
Base mean        1086.974295
log2(FC)            -4.14845
StdErr              0.134949
Wald-Stats        -30.740913
P-value                  0.0
P-adj                    0.0
Chromosome             chr3R
Start               24141394
End                 24147490
Strand                     +
Feature       protein_coding
Gene name               Kal1
Name: FBgn0039155, dtype: object 
```
Which would get the same result as printing `data.loc["FBgn0039155"]` (without a second index).

You can select multiple columns or rows using `DataFrame.loc` and a named slice or `Dataframe.iloc` and the numbers corresponding to the rows and columns.

```python
print(data.loc['FBgn0003360':'FBgn0029167', 'Base mean':'Wald-Stats'])
print(data.iloc[1:4 , 0:3])
```
```output
                Base mean  log2(FC)    StdErr  Wald-Stats
GeneID                                                   
FBgn0003360   6409.577128 -2.999777  0.104345  -28.748637
FBgn0026562  65114.840564 -2.380164  0.084327  -28.225437
FBgn0025111   2192.322369  2.699939  0.097945   27.565978
FBgn0029167   5430.067277 -2.105062  0.092547  -22.745964

                Base mean  log2(FC)    StdErr
GeneID                                       
FBgn0003360   6409.577128 -2.999777  0.104345
FBgn0026562  65114.840564 -2.380164  0.084327
FBgn0025111   2192.322369  2.699939  0.097945
```

* Note the difference between the 2 outputs.

When choosing or transitioning between `loc` and `iloc`, you should keep in mind that the two methods use slightly different indexing schemes.

`iloc` uses the Python stdlib indexing scheme, where the first element of the range is included and the last one excluded. So `0:10` will select entries `0,...,9`. `loc`, meanwhile, indexes inclusively. So `0:10` will select entries `0,...,10`.

This is particularly confusing when the DataFrame index is a simple numerical list, e.g. `0,...,1000`. In this case `df.iloc[0:1000]` will return 1000 entries, while `df.loc[0:1000]` return 1001 of them! To get 1000 elements using `loc`, you will need to go one lower and ask for `df.loc[0:999]`.

The result of slicing is a new dataframe and can be used in further operations. All the statistical operators that work on entire dataframes work the same way on slices. E.g., calculate max of a slice.

```python
print(data.loc['FBgn0003360':'FBgn0029167', 'Base mean'].max())
```
```output
65114.8405637953
```

## Use conditionals to select data

You can use conditionals to select data. A comparison is applied element by element and returns a similarly-shaped dataframe of `True` and `False`. The last one can be used as a mask to subset the original dataframe. The following example creates a new dataframe consisting only of the columns 'P-adj' and 'Gene name', then keeps the rows that comply with the expression `'P-adj' < 0.000005`

```python
subset = data.loc[:, ['P-adj', 'Gene name']]
print(subset)
```
```output
                     P-adj    Gene name
GeneID                                 
FBgn0039155  1.387207e-203         Kal1
FBgn0003360  4.039734e-178         sesB
FBgn0026562  8.149380e-172  BM-40-SPARC
FBgn0025111  6.105789e-164         Ant2
FBgn0029167  2.698811e-111          Hml
...                    ...          ...
FBgn0035710   1.436480e-04       SP1173
FBgn0035523   1.523189e-04       CG1311
FBgn0038261   1.525142e-04      CG14856
FBgn0039178   1.629061e-04       CG6356
FBgn0034636   1.951192e-04      CG10440

[130 rows x 2 columns]
```
```python
mask = subset.loc[:, 'P-adj'] < 0.000005
new_data = subset[mask]
print(new_data)
```
```output

                     P-adj    Gene name
GeneID                                 
FBgn0039155  1.387207e-203         Kal1
FBgn0003360  4.039734e-178         sesB
FBgn0026562  8.149380e-172  BM-40-SPARC
FBgn0025111  6.105789e-164         Ant2
FBgn0029167  2.698811e-111          Hml
...                    ...          ...
FBgn0039593   1.544012e-07       CG9989
FBgn0265512   3.413938e-07          mlt
FBgn0030326   3.625857e-07       CG2444
FBgn0039485   7.174821e-07      CG17189
FBgn0025836   1.071884e-06     RhoGAP1A

[113 rows x 2 columns]
```

If we have not had specified the column, that the expression should be applied to, then it would have been applied to the entire dataframe. But the dataframe contains different type of data. In that case, an error would occur. 

Consider the following example of a dataframe consisting only of numerical data. The expression and the mask would be normally applied to the data and the mask would return `NaN` for the data that don't comply with the expression. 

```python
subset = data.loc[:, ['StdErr',	'Wald-Stats', 'P-value', 'P-adj']]
mask = subset < 0.05
new_data = subset[mask]
print(new_data)
```
```output
             StdErr  Wald-Stats        P-value          P-adj
GeneID                                                       
FBgn0039155     NaN  -30.740913  1.617357e-207  1.387207e-203
FBgn0003360     NaN  -28.748637  9.419922e-182  4.039734e-178
FBgn0026562     NaN  -28.225437  2.850430e-175  8.149380e-172
FBgn0025111     NaN         NaN  2.847517e-167  6.105789e-164
FBgn0029167     NaN  -22.745964  1.573284e-114  2.698811e-111
...             ...         ...            ...            ...
FBgn0035710     NaN         NaN   6.682473e-06   1.436480e-04
FBgn0035523     NaN         NaN   7.103599e-06   1.523189e-04
FBgn0038261     NaN         NaN   7.130487e-06   1.525142e-04
FBgn0039178     NaN         NaN   7.654328e-06   1.629061e-04
FBgn0034636     NaN   -4.430118   9.418135e-06   1.951192e-04

[130 rows x 4 columns]
```
This is very useful because NaNs are ignored by operations like max, min, average, etc.

```python
print(new_data.describe())
```
```output
       StdErr  Wald-Stats        P-value          P-adj
count     0.0   70.000000   1.300000e+02   1.300000e+02
mean      NaN  -11.035514   4.200552e-07   9.321175e-06
std       NaN    5.695077   1.531253e-06   3.278527e-05
min       NaN  -30.740913  1.617357e-207  1.387207e-203
25%       NaN  -13.059226   4.952408e-31   1.275950e-28
50%       NaN   -9.971703   2.109219e-19   2.563262e-17
75%       NaN   -6.963977   6.731953e-12   3.706020e-10
max       NaN   -4.430118   9.418135e-06   1.951192e-04
```

### Exercises
1. Explain what each line in the following short program does: what is in first, second, etc.?

    ```python
    first = pd.read_csv("https://zenodo.org/record/3477564/files/annotatedDEgenes.tabular", sep = "\t", index_col = 'GeneID')
    second = first[first['log2(FC)'] > 0 ]
    third = second.drop('FBgn0025111')
    fourth = third.drop('StdErr', axis = 1)
    fourth.to_csv('result.csv')
    ```
2.  Explain in simple terms what `idxmin` and `idxmax` do in the short program below. When would you use these methods?

    ```python
    data = pd.read_csv("https://zenodo.org/record/3477564/files/annotatedDEgenes.tabular", sep = "\t", index_col = 'GeneID')
    
    print(data['Base mean'].idxmin())
    print(data['Base mean'].idxmax())
    
    ```
3. Assume Pandas has been imported and the previous dataset has been loaded. Write an expression to select each of the following:

   a. P-value of each gene
   b. all the information of gene `FBgn0039178`
   c. the information of all genes that belong to chromosome `chr3R`

### Solutions
1. Let’s go through this piece of code line by line.

    ```python
    first = pd.read_csv("https://zenodo.org/record/3477564/files/annotatedDEgenes.tabular", sep = "\t", index_col = 'GeneID')
    ```
    This line loads the data into a dataframe called first. The `index_col='GeneID'` parameter selects which column to use as the row labels in the dataframe.

    ```python
    second = first[first['log2(FC)'] > 0 ]
    ```
    This line makes a selection: only those rows of first for which the ‘log2(FC)’ column contains a positive value are extracted. Notice how the Boolean expression inside the brackets is used to select only those rows where the expression is true. 

    ```python
    third = second.drop('FBgn0025111')
    ```
    As the syntax suggests, this line drops the row from second where the label is ‘FBgn0025111’. The resulting dataframe third has one row less than the original dataframe second.

    ```python
    fourth = third.drop('StdErr', axis = 1)
    ```
    Again we apply the drop function, but in this case we are dropping not a row but a whole column. To accomplish this, we need to specify also the axis parameter.

    ```python
    fourth.to_csv('result.csv')
    ```
    The final step is to write the data that we have been working on to a csv file. Pandas makes this easy with the `to_csv()` function. The only required argument to the function is the filename. Note that the file will be written in the directory from which you started the Jupyter or Python session.

2. `idxmin` will return the index value corresponding to the minimum; idxmax will do the same for the maximum value.

    You can use these functions whenever you want to get the row index of the minimum/maximum value and not the actual minimum/maximum value.

    ```output
    FBgn0063667
    FBgn0026562
	```
3. 
    a. 

    ```python
     data['P-value']
	```
    b.
    ```python
       data.loc['FBgn0039178', :]
    ```
    c. 
    ```python
       data[data['Chromosome'] == 'chr3R']
    ```

## Group-by and analyze the data

Many data analysis tasks can be approached using the “split-apply-combine” paradigm: split the data into groups, apply some analysis to each group, and then combine the results.

Pandas makes this very easy through the use of the `groupby()` method, which splits the data into groups. When the data is grouped in this way, the aggregate method `agg()` can be used to apply an aggregating or summary function to each group.

```python
summarised_data = data.groupby('Chromosome').agg({'Base mean':'first', 
                             'log2(FC)': 'max'})
print(summarised_data)
```
```output
              Base mean  log2(FC)
Chromosome                       
chr2L       1269.362573  2.145440
chr2R        928.263812  2.407986
chr3L       5430.067277  2.428607
chr3R       1086.974295  2.360017
chrX        6409.577128  2.699939
```

There are a couple of things that should be noted. The `agg()` method accepts a dictionary as input that specifies the function to be applied to each column. The output is a new dataframe, that each row corresponds to one group. The output dataframe uses the grouping column as index. We could change the last one by simply using the `reset_index()` method.

```python
summarised_data = data.groupby('Chromosome').agg({'Base mean':'first', 
                                                  'log2(FC)': 'max'}).reset_index() 
print(summarised_data)
```
```output
  Chromosome    Base mean  log2(FC)
0      chr2L  1269.362573  2.145440
1      chr2R   928.263812  2.407986
2      chr3L  5430.067277  2.428607
3      chr3R  1086.974295  2.360017
4       chrX  6409.577128  2.699939
```
### Exercises

1. Import the tabular data from the "https://zenodo.org/record/3477564/files/annotatedDEgenes.tabular" link. What are the longest genes in each chromosome?
2. Using the same dataset,  try to find how many genes are found on each strand of each chromosome.


### Solutions
1. 
    ```python
    data = pd.read_csv("https://zenodo.org/record/3477564/files/annotatedDEgenes.tabular", sep = "\t", index_col = 'GeneID')
    
    data['Gene Length'] = data['End'] - data['Start']
    data.groupby('Chromosome').agg(max_length = ('Gene Length', 'max'))
    ```
2. You can group the data according to more than one column.
    ```python
    data.groupby(['Chromosome', 'Strand']).size()
    ```



# Conclusion

This tutorial aims to serve as an introduction to data analysis using the Python programming language.

## Key points

* Python has many libraries offering a variety of capabilities, which makes it popular for beginners, as well as, more experienced users

* You can use scientific libraries like Numpy and Pandas to perform data analysis.



# References

1. Wendi Bacon, Mehmet Tekman, 2021, Trajectory Analysis using Python (Jupyter Notebook) in Galaxy (Galaxy Training Materials). https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-JUPYTER-trajectories/tutorial.html Online; accessed Wed Sep 01 2021
2. Batut et al., 2018, Community-Driven Data Analysis Training for Biology Cell Systems [10.1016/j.cels.2018.05.012
3.  Azalee Bostroem, Trevor Bekolay, and Valentina Staneva (eds): "Software Carpentry: Programming with Python."  Version 2016.06, June
2016, https://github.com/swcarpentry/python-novice-inflammation, 10.5281/zenodo.57492.
4. Allen Lee, Nathan Moore, Sourav Singh, Olav Vahtras (eds). Software Carpentry: Plotting and Programming in Python. http://github.com/swcarpentry/python-novice-plotting, 2018.
5. Python 3.9.7 Documentation, https://docs.python.org/3/tutorial/datastructures.html Online; accessed Wed Sep 01 2021
6. Karlijn Willems, Pandas Tutorial: DataFrames in Python. https://www.datacamp.com/community/tutorials/pandas-tutorial-dataframe-python Online; accessed Wed Sep 01 2021
7. Bérénice Batut, Mallory Freeberg, Mo Heydarian, Anika Erxleben, Pavankumar Videm, Clemens Blank, Maria Doyle, Nicola Soranzo, Peter van Heusden, 2021 Reference-based RNA-Seq data analysis (Galaxy Training Materials). https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html Online; accessed Wed Sep 01 2021
8. Bérénice Batut, Fotis E. Psomopoulos, Toby Hodges, 2021 Advanced R in Galaxy (Galaxy Training Materials). https://training.galaxyproject.org/training-material/topics/introduction/tutorials/r-advanced/tutorial.html Online; accessed Wed Sep 01 2021

