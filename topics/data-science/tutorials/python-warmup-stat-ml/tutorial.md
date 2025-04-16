---
layout: tutorial_hands_on
title: Python - Warm-up for statistics and machine learning
level: Introductory
draft: true
requirements:
 - type: internal
   topic_name: data-science
   tutorials:
   - python-basics
questions:
- to do
objectives:
- to do
key_points:
- to do
time_estimation: 1H
tags:
- elixir
- ai-ml
subtopic: python
contributions:
    authorship:
    - wandrilled
    editing:
    - bebatut
priority: 1
notebook:
  language: python
  pyolite: true
---

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Basic python



```python

X = []

for i in range(10):
    X.append( i**2 )

print(X)
```

    [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]



```python

for x in X:
    print(x)

```

    0
    1
    4
    9
    16
    25
    36
    49
    64
    81



```python
for x in X:
    if x%2 == 1:
        print(x,'is odd')
    else:
        print(x,'is even')
```

    0 is even
    1 is odd
    4 is even
    9 is odd
    16 is even
    25 is odd
    36 is even
    49 is odd
    64 is even
    81 is odd



```python
# list comprehension is a very fine way of compressing all this

X = [ i**2 for i in range(10) ]

Xeven = [ x for x in X if x%2 == 0 ]
Xodd = [ x for x in X if x%2 == 1 ]


print( 'X    ', X )
print( 'Xeven', Xeven )
print( 'Xodd ', Xodd )
```

    X     [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
    Xeven [0, 4, 16, 36, 64]
    Xodd  [1, 9, 25, 49, 81]


[back to the top](#top)

## numpy and vectorized operations


```python
import numpy as np

X_array = np.array(X)

print(X_array)
```

    [ 0  1  4  9 16 25 36 49 64 81]



```python
print(X_array / 2 )
```

    [ 0.   0.5  2.   4.5  8.  12.5 18.  24.5 32.  40.5]



```python
print( np.exp(X_array ) )
print( np.log(X_array ) )
```

    [1.00000000e+00 2.71828183e+00 5.45981500e+01 8.10308393e+03
     8.88611052e+06 7.20048993e+10 4.31123155e+15 1.90734657e+21
     6.23514908e+27 1.50609731e+35]
    [      -inf 0.         1.38629436 2.19722458 2.77258872 3.21887582
     3.58351894 3.8918203  4.15888308 4.39444915]


    /tmp/ipykernel_490123/2855859755.py:2: RuntimeWarning: divide by zero encountered in log
      print( np.log(X_array ) )



```python
print( 'shape' , X_array.shape )
print( 'mean ' , np.mean(X_array) )
print( 'standard deviation' , np.std(X_array) )
```

    shape (10,)
    mean  28.5
    standard deviation 26.852374196707448


### linspace and arange


```python
print( 'linspace 0,2,9 :' , np.linspace(0,2,9) , sep='\t' )
print( 'linspace -0.5,0.5,11 :' , np.linspace(-0.5,0.5,11) , sep='\t' )
print( 'linspace 10,0,11 :' , np.linspace(10,0,11) , sep='\t' )
```

    linspace 0,2,9 :	[0.   0.25 0.5  0.75 1.   1.25 1.5  1.75 2.  ]
    linspace -0.5,0.5,11 :	[-0.5 -0.4 -0.3 -0.2 -0.1  0.   0.1  0.2  0.3  0.4  0.5]
    linspace 10,0,11 :	[10.  9.  8.  7.  6.  5.  4.  3.  2.  1.  0.]



```python
print( "arange 0,2,0.1 :", np.arange(1.5,2,0.1) , sep='\t' )
print( "arange -1,1,0.125 :", np.arange(-1,1,0.125) , sep='\t' )
print( "arange 10,2 :", np.arange(10,2,1) , sep='\t' ) # reverse does not work!
```

    arange 0,2,0.1 :	[1.5 1.6 1.7 1.8 1.9]
    arange -1,1,0.125 :	[-1.    -0.875 -0.75  -0.625 -0.5   -0.375 -0.25  -0.125  0.     0.125
      0.25   0.375  0.5    0.625  0.75   0.875]
    arange 10,2 :	[]


## Basic plotting


```python
import matplotlib.pyplot as plt

plt.plot( [0,1,2,3] , [10,5,7,0.2] )
plt.show()
```



### Adding color, symbols, ...

`matplotlib` offers many options to customize the appearance of your plot.

Here are the (some) common arguments to `plot()` (which can also be applied to many other graphical representations):
 * `color` : could be given as a (red,green,blue) tuple, a [name](https://matplotlib.org/3.1.0/gallery/color/named_colors.html), a hex code, ...  (see [Something better here](https://matplotlib.org/tutorials/colors/colors.html) for all the options)
 * `marker` : symbols for the data point. `'.'` is a point, `'v'` a down triangle, ... see [Something better here](https://matplotlib.org/3.3.3/api/markers_api.html#module-matplotlib.markers) for the list of possibilities.
 * `linestyle` : style of the line. `'-'` is solid, `'--'` is dashed, `''` for no line. See [Something better here](https://matplotlib.org/3.3.3/gallery/lines_bars_and_markers/linestyles.html) for more options
 * `linewidth` : width of the lines
 * `markersize` : size of the markers

You are invited to experiment and explore these options. Here are a few examples:



```python
y1 = [1,2,3,10,5]
y2 = [10,9,7,5.5,6]
y3 = [4,3,1.5,1]

# green, dashed line, with circle markers
plt.plot( y1, color = 'green', marker = 'o', linestyle = '--', linewidth = 2, markersize = 8 )

# blue triangle with no line
plt.plot( y2, color = 'blue', marker = 'v', linestyle = '' , markersize = 16 )

# solid orange line
plt.plot(y3, color = 'orange', marker = '', linestyle = '-', linewidth = 4 )

plt.show()
```
    


Note that:
 * you can call plot several time in a row to make several lines appear (only `plt.show()` causes the figure to appear)
 * the frame of the picture automatically adjust to what it needs to show

### multiple subplots

Now would normally be when we show you how to add labels, titles and legends to figures. 

However, the way `matplotlib` is built, it is actually a bit more efficient to first learn how to create multiple subplots.


Creating multiple plots is possible with the function `plt.subplots()`.
Amon its many arguments, it takes:
 * `nrows` : number of subplot rows
 * `ncols` : number of subplot columns
 * `figsize` : tuple (width,height) of the figure

This function creates a Figure and an Axes object.
The Axes object can be either : 
 * a simple Axe is there is 1 row and 1 columns
 * a list of Axe objects if there is 1 row and multiple columns, or 1 column and multiple rows
 * a list of lists of Axes objects if there is multiple rows and multiple columns



```python
y1 = [1,2,3,10,5]
y2 = [10,9,7,5.5,6]
y3 = [4,3,1.5,1]


# subplots returns a Figure and an Axes object
fig, ax = plt.subplots(nrows=1, ncols=2) # 2 columns and 1 row

# ax is a list with two objects. Each object correspond to 1 subplot

# accessing to the first column ax[0]
ax[0].plot( y1, color = 'green', marker = 'o', linestyle = '--', linewidth = 2, markersize = 8 )

# accessing to the second column ax[1]
ax[1].plot( y2, color = 'blue', marker = 'v', linestyle = '' , markersize = 16 )
ax[1].plot( y3, color = 'orange', marker = '', linestyle = '-' )

plt.show()
```
    


Notice how we call `ax[0].plot(...)` instead of `plt.plot(...)` to specify in which subplots we want to plot.

### multiple subplots - continued

Let's see the same thing with several lines and several columns


```python
y1 = [1,2,3,10,5]
y2 = [10,9,7,5.5,6]
y3 = [4,3,1.5,1]
y4 = [1,2,3,7,5]

# 2 columns and 2 rows, and we also set the figure size
fig, ax = plt.subplots(nrows=2, ncols=2 , figsize = (12,12))

# ax is a list of two lists with two objects each.

# accessing to the first row, first column : ax[0][0]
ax[0][0].plot( y1, color = 'green', marker = 'o', linestyle = '--', linewidth = 2, markersize = 8 )

# accessing to the first row, second column : ax[0][1]
ax[0][1].plot( y2, color = 'blue', marker = 'v', linestyle = '' , markersize = 16 )

# accessing to the second row, first column : ax[1][0]
ax[1][0].plot( y3, color = 'orange', marker = 'x', linestyle = '-' )

# accessing to the first row, second column : ax[1][1]
ax[1][1].plot( y4, color = 'teal', linestyle = '-.' , linewidth=5 )

plt.show()
```
    


### setting up labels

To set the labels at the x-axis, y-axis and title, we use the method of the Axe object:
 * `.set_xlabel(...)`
 * `.set_ylabel(...)`
 * `.set_title(...) `



```python
y1 = [1,2,3,10,5]
y2 = [10,9,7,5.5,6]
y3 = [4,3,1.5,1]

# subplots returns a Figure and an Axes object
fig, ax = plt.subplots(nrows=1, ncols=2 , figsize=(10,5)) # 2 columns and 1 row


# accessing to the first column ax[0]
ax[0].plot( y1, color = 'green', marker = 'o', linestyle = '--', linewidth = 2, markersize = 8 )
ax[0].set_xlabel('x-axis label')
ax[0].set_ylabel('y-axis label')
ax[0].set_title('plot 1')


# accessing to the second column ax[1]
ax[1].plot( y2, color = 'blue', marker = 'v', linestyle = '' , markersize = 16 )
ax[1].plot( y3, color = 'orange', marker = '', linestyle = '-' )
ax[1].set_xlabel('x-axis label')
ax[1].set_ylabel('y-axis label')
ax[1].set_title('plot 2')

plt.show()
```
    


**setting up a legend** 

Each element we add to the figure using `plot()` can be given a label using the `label` argument.
Then, a legend may be added to the figure using the `legend()` method.

This `legend()` method can take a `loc` argument that specifies where it should be plotted. 
Possible values for this argument are: `'best' , 'upper right' , 'upper left' , 'lower left' , 'lower right' , 'right' , 'center left' , 'center right' , 'lower center' , 'upper center' , 'center'` (the default is `best`).



```python

fig, ax = plt.subplots(nrows=1, ncols=1 , figsize=(10,5)) # 2 columns and 1 row

# NB : with 1 col and 1 row, ax is directly the sole subplot we have
#      so to call it we just use ax.plot , ax.set_xlabel , ...

ax.plot( y1, color = 'green', marker = 'o', linestyle = '--', linewidth = 2 , label = 'line A' )
ax.plot( y2, color = 'blue', marker = 'v', linestyle = '' , markersize =  8 , label = 'line B' )
ax.plot( y3, color = 'orange', marker = '', linestyle = '-' , linewidth = 2 , label = 'line C' )

ax.set_xlabel('x-axis label')
ax.set_ylabel('y-axis label')
ax.set_title('plot with a legend')

#adding a legend in the upper right
ax.legend( loc='upper right')

plt.show()

```

    


### additional : writing a figure to a file

Writing a matplotlib figure to a file can be achieved simply by replacing the call to `plt.show()` to `plt.savefig(...)`.

`plt.savefig` takes a number of argument, the most commons are :
 * `fname` : name of the file to write the figure. The extension is used to determine the output format (.pdf,.png, .jpg , .svg ,  ...). Many formats are supported, you can get a list with this command : `plt.gcf().canvas.get_supported_filetypes()`
 * `dpi` : dots per inches , useful to set-up when saving to raster formats (ie., pixel-based such as png or jpeg). The actual size of the image is set using the argument `figsize` of `plt.subplots()`


> <comment-title></comment-title>
>
> in a jupyter notebook the figure will still be shown, whereas in a standard .py script it will not appear on screen.
>
{: .comment}

Here is a demonstration. Apply in on your side and verify that the file `testPlot.png` was created:


```python
import matplotlib.pyplot as plt

y1 = [1,2,3,10,5]
y2 = [10,9,7,5.5,6]
y3 = [4,3,1.5,1]


# subplots returns a Figure and an Axes object
fig, ax = plt.subplots(nrows=1, ncols=2 , figsize = (10,6) ) # 2 columns and 1 row

# ax is a list with two objects. Each object correspond to 1 subplot

# accessing to the first column ax[0]
ax[0].plot( y1, color = 'green', marker = 'o', linestyle = '--', linewidth = 2, markersize = 8 )

# accessing to the second column ax[1]
ax[1].plot( y2, color = 'blue', marker = 'v', linestyle = '' , markersize = 16 )
ax[1].plot( y3, color = 'orange', marker = '', linestyle = '-' )

plt.savefig( 'testPlot.png' , dpi = 90  )
```


## Exercise 00.01 : bringing together numpy and matplotlib

Numpy arrays can be plotted as if they were lists.

1. plot x and y, where:
    * y = 1/(1+exp(-x))
    * x varies between -5 and 5 (plotting around a 100 points should suffice).

2. **Bonus :** plot multiples lines : y = 1/(1+exp(-x*b)) , for the following values of b: 0.5 , 1 , 2 , 4.
    * x still varies between -5 and 5 (plotting around a 100 points should suffice).
    * put a legend in your plot




```python

```

You can load the solution directly in this notebook by uncommenting and running the following line:


```python
# %load  -r -8 solutions/solution_00_01.py
```

bonus question solution:


```python
# %load  -r 9- solutions/solution_00_01.py
```


## Generating random numbers


### the basics


```python
import numpy.random as rd

# random floats between 0 and 1
for i in range(4):
    print( rd.random() )

```

    0.6696103730869407
    0.7426639266737763
    0.6767219223242785
    0.8602105555191791



```python
print( rd.random(size=10) ) # draw directly 10 numbers
```

    [0.37971723 0.80354745 0.4168427  0.70867247 0.17547126 0.43760884
     0.75933345 0.06571168 0.45772397 0.67191214]


### setting the seed: pseudorandomness and reproducibility


```python
rd.seed(42) # setting the seed to 42
print( '1st draw' , rd.random(size=5) )
print( '2nd draw' , rd.random(size=5) )
rd.seed(42)
print( 'after resetting seed' , rd.random(size=5) )
```

    1st draw [0.37454012 0.95071431 0.73199394 0.59865848 0.15601864]
    2nd draw [0.15599452 0.05808361 0.86617615 0.60111501 0.70807258]
    after resetting seed [0.37454012 0.95071431 0.73199394 0.59865848 0.15601864]


### beyond the uniform distribution

numpy offers you quite a large [set of distributions you can draw from](https://docs.scipy.org/doc/numpy-1.15.0/reference/routines.random.html#distributions).

Let's look at the normal distribution:


```python

normalDraw = rd.normal(size = 1000 )

print( 'mean ' , np.mean( normalDraw ) )
print( 'stdev' , np.std( normalDraw ) )
```

    mean  0.025354699638558926
    stdev 1.0003731428167348



```python
normalDraw2 = rd.normal( loc = -2 , scale = 3 , size = 300 ) # loc chnages the location (mean), and scale changes the standard deviation

print( 'mean ' , np.mean( normalDraw2 ) )
print( 'stdev' , np.std( normalDraw2 ) )
```

    mean  -1.9773491637651965
    stdev 2.964622032924749


of course, we could want to plot these drawn numbers:


```python
plt.hist( normalDraw  , alpha = 0.5 , label='loc=0  , scale=1')
plt.hist( normalDraw2 , alpha = 0.5 , label='loc=-2 , scale=3')
plt.legend()
plt.show()
```


## Statistical testing

`numpy.random` let's you draw random numbers ;
`scipy.stats` implements the probability density functions, and Percent point function, as well as the most statistical tests.



```python
import scipy.stats as stats

# plotting the probability density function for 1 of the random draw we just made:

x = np.linspace(-10,10,1001)

normPDF = stats.norm.pdf( x , loc = -2 , scale = 3 )

plt.hist( normalDraw2 , alpha = 0.5 , label='random draw' , density = True) # don't forget density=True
plt.plot(x,normPDF , label='PDF' )
plt.legend()
plt.show()
```


We can also get the expected quantiles of a distribution:


```python
print( '95% quantile of a Chi-square distribution with 3 degrees of freedom:', stats.chi2.ppf(0.95 , df=3))
print( 'fraction of a Chi-square distribution with 3 degrees of freedom above of equal to 5' ,  
      1 - stats.chi2.cdf( 5 , df=3 ) )
```

    95% quantile of a Chi-square distribution with 3 degrees of freedom: 7.814727903251179
    fraction of a Chi-square distribution with 3 degrees of freedom above of equal to 5 0.17179714429673354


And you can apply some classical statistical tests:


```python
# t-test of independance between two random samples:
rd.seed(73)

s1 = rd.normal(size=67)
s2 = rd.normal(size=54 , loc = 0.2)

testStat , pval = stats.ttest_ind(s1,s2 , equal_var=True)  # equal variance : Student's t-test ; unequal : Welch's
#almost all of these stat functions return the same test-statistic , pvalue tuple

print('result of the t-test')
print('\tt:',testStat)
print('\tp-value:',pval)
```

    result of the t-test
        t: 0.26673986193074073
        p-value: 0.7901311339594405


### What is our conclusion for these tests results? What do you think about this?


```python

# Kolmogorov-smirnov test for a chi-square distribution

sample = rd.chisquare(df=13 , size = 43)


# kstest expect as second argument the cdf function of the reference distribution
# this is how to handle the fact that me must set an argument (degree of freedom)
refDistribution = stats.chi2(df=13).cdf

testStat , pval = stats.kstest( sample , refDistribution )
# alternative : 
# testStat , pval = stats.kstest( sample , lambda x : stats.chi2.cdf(x , df=13 ) )

print('result of the Kolmogorov-Smirnov test comparing our sample to a Chi-square distribution with 13 degrees of freedom')
print('\tK:',testStat)
print('\tp-value:',pval)

```

    result of the Kolmogorov-Smirnov test comparing our sample to a Chi-square distribution with 13 degrees of freedom
        K: 0.12249766392962913
        p-value: 0.5003109000967569


If you are interested, this [webpage](https://machinelearningmastery.com/statistical-hypothesis-tests-in-python-cheat-sheet/) references all implemented tests, with examples.

[back to the top](#top)

## Bringing together numpy, numpy.random, and matplotlib

The random generation function return a numpy array, meaning it is fairly trivial to combine it with other arrays:



```python
# combining 

x = np.sort( rd.normal(loc=170 , scale = 23 , size = 100) )

y_theoretical = 0.75 * x + 100 # simple linear relationship : y = a * x + b

measurement_noise = rd.normal(scale = 10 , size = 100) # some noise associated to the measure

y_observed = y_theoretical + measurement_noise # observed = expected + noise

fig,ax = plt.subplots(figsize=(8,8))
plt.plot( x , y_theoretical , label = 'expected' )
plt.plot( x , y_observed , marker = '.' , linestyle='' , alpha = 0.7 , label = 'observed')
plt.legend()
plt.show()
```



## The briefest intro to pandas

`pandas` is a powerful library when doing data analysis, especially in the forms of table.

Basically, it reimplements R data.frame as a DataFrame object and ties together neatly with the libraries we've just seen.



```python
 import pandas as pd
    
df = pd.read_table( 'data/beetle.csv' , sep=',' , index_col=0 ) # pandas automatically detects header.

df.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>dose</th>
      <th>nexp</th>
      <th>ndied</th>
      <th>prop</th>
      <th>nalive</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>49.1</td>
      <td>59</td>
      <td>6</td>
      <td>0.102</td>
      <td>53</td>
    </tr>
    <tr>
      <th>2</th>
      <td>53.0</td>
      <td>60</td>
      <td>13</td>
      <td>0.217</td>
      <td>47</td>
    </tr>
    <tr>
      <th>3</th>
      <td>56.9</td>
      <td>62</td>
      <td>18</td>
      <td>0.290</td>
      <td>44</td>
    </tr>
    <tr>
      <th>4</th>
      <td>60.8</td>
      <td>56</td>
      <td>28</td>
      <td>0.500</td>
      <td>28</td>
    </tr>
    <tr>
      <th>5</th>
      <td>64.8</td>
      <td>63</td>
      <td>52</td>
      <td>0.825</td>
      <td>11</td>
    </tr>
  </tbody>
</table>
</div>




```python
Nrows, Ncols = df.shape
print( 'number of rows:',Nrows, 'number of columns:', Ncols )
print( 'column names' , df.columns )
```

    number of rows: 8 number of columns: 5
    column names Index(['dose', 'nexp', 'ndied', 'prop', 'nalive'], dtype='object')



```python
df.describe()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>dose</th>
      <th>nexp</th>
      <th>ndied</th>
      <th>prop</th>
      <th>nalive</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>count</th>
      <td>8.000000</td>
      <td>8.000000</td>
      <td>8.000000</td>
      <td>8.000000</td>
      <td>8.000000</td>
    </tr>
    <tr>
      <th>mean</th>
      <td>62.800000</td>
      <td>60.125000</td>
      <td>36.375000</td>
      <td>0.602000</td>
      <td>23.750000</td>
    </tr>
    <tr>
      <th>std</th>
      <td>9.599702</td>
      <td>2.232071</td>
      <td>22.557466</td>
      <td>0.367937</td>
      <td>21.985385</td>
    </tr>
    <tr>
      <th>min</th>
      <td>49.100000</td>
      <td>56.000000</td>
      <td>6.000000</td>
      <td>0.102000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>25%</th>
      <td>55.925000</td>
      <td>59.000000</td>
      <td>16.750000</td>
      <td>0.271750</td>
      <td>4.750000</td>
    </tr>
    <tr>
      <th>50%</th>
      <td>62.800000</td>
      <td>60.000000</td>
      <td>40.000000</td>
      <td>0.662500</td>
      <td>19.500000</td>
    </tr>
    <tr>
      <th>75%</th>
      <td>69.675000</td>
      <td>62.000000</td>
      <td>54.750000</td>
      <td>0.919500</td>
      <td>44.750000</td>
    </tr>
    <tr>
      <th>max</th>
      <td>76.500000</td>
      <td>63.000000</td>
      <td>61.000000</td>
      <td>1.000000</td>
      <td>53.000000</td>
    </tr>
  </tbody>
</table>
</div>




```python
# select a single column:
df['dose']
```




    1    49.1
    2    53.0
    3    56.9
    4    60.8
    5    64.8
    6    68.7
    7    72.6
    8    76.5
    Name: dose, dtype: float64




```python
df[ ['ndied','nalive'] ] # select several columns
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ndied</th>
      <th>nalive</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>6</td>
      <td>53</td>
    </tr>
    <tr>
      <th>2</th>
      <td>13</td>
      <td>47</td>
    </tr>
    <tr>
      <th>3</th>
      <td>18</td>
      <td>44</td>
    </tr>
    <tr>
      <th>4</th>
      <td>28</td>
      <td>28</td>
    </tr>
    <tr>
      <th>5</th>
      <td>52</td>
      <td>11</td>
    </tr>
    <tr>
      <th>6</th>
      <td>53</td>
      <td>6</td>
    </tr>
    <tr>
      <th>7</th>
      <td>61</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td>60</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>



### Plotting DataFrame Columns

Because `DataFrame` columns are iterable, they can seamlessly be given as argument to `plot()`.


```python

# plotting the column dose along the x-axis and prop along the y-axis
# I use the + marker, with a teal color.
plt.plot(df['dose'] , df['prop'] , color = 'teal' , linestyle='' , marker = '+' , markersize=10 )
plt.xlabel( 'dose' )
plt.ylabel( 'proportion of dead' )
plt.show()
```
    


DataFrame column can be manipulated like numpy array:


```python

## we can combine columns using normal operators
Odds = df['nalive'] /df['ndied'] # the odds of being alive is nalive / ndead

## adding a new column to the DataFrame is trivial:
df['Odds'] = Odds


## we can also apply numpy function to them
df['logOdds'] = np.log( df['Odds'] )


plt.plot(df['dose'] , df['logOdds'] , color = 'teal' , linestyle='' , marker = '+' , markersize=10 )
plt.xlabel( 'dose' )
plt.ylabel( 'log Odds' )
plt.show()

```



## Exercise 00.02 : tying everything together

1. Read the file `'data/kyphosis.csv'`.
2. how many columns are there ?
3. What is the maximum Age ? 
4. create a new column `Stop` , corresponding to the addition of columns `'Start'` and `'Number'`
5. plot the relationship between `'Age'` and `'Number'` (bonus point : use colors to indicate the presence or absence of kyphosis ).





```python

```

Solutions:


```python
# %load  -r -7 solutions/solution_00_02.py
```


```python
# %load  -r 8-9 solutions/solution_00_02.py
```


```python
# %load  -r 11-12 solutions/solution_00_02.py
```


```python
# %load  -r 14-15 solutions/solution_00_02.py
```


```python
# %load  -r 17-22 solutions/solution_00_02.py
```


```python
# %load  -r 24- solutions/solution_00_02.py
```
