# Galaxy 101-1: The first thing you should try

In this very simple example we will introduce you to bare basics of Galaxy:

* Getting data from UCSC
* Performing simple data manipulation
* Understanding Galaxy's History system

## What are we trying to do?

Suppose you get the following question: ` Mom (or Dad) ... Which coding exon has the highest number of single nucleotide polymorphisms on chromosome 22?`.  You think to yourself "Wow! This is a simple question ... I know exactly where the data is (at UCSC) but how do I actually compute this?" The truth is, there is really no straightforward way of answering this question in a time frame comparable to the attention span of a 7-year-old. Well ... actually there is and it is called Galaxy. So let's try it...

## 0. Organizing your windows and setting up Galaxy account

### 0.0. Getting your display sorted out

To get the most of this tutorial open two browser windows. One you already have (it is this page). To open the other, **right** click [this link](http://usegalaxy.org) and choose "Open in a New Window" (or something similar depending on your operating system and browser):

![open in a new window](http://galaxy.psu.edu/galaxy101/newWindow.png)

Then organize your windows as something like this (depending on the size of your monitor you may or may not be able to organize things this way, but you get the idea):

![Windows side by side](https://galaxyproject.org/galaxy101/twoScreens.png)

### 0.1. Setting up Galaxy account
Go to the **User** link at the top of Galaxy interface and choose **Register** (unless of course you already have an account):

![register](http://galaxy.psu.edu/galaxy101/register.png)

Then enter your information and you're in!

## 1. Getting data from UCSC
### 1.0. Getting coding exons
First thing we will do is to obtain data from UCSC by clicking **Get Data -> UCSC Main**:

![get data from UCSC](http://galaxy.psu.edu/galaxy101/getDataUCSC.png)

You will see UCSC Table Browser interface appearing in your browser window:

![UCSC genes](http://galaxy.psu.edu/galaxy101/ucscGenes.png)

Make sure that your settings are exactly the same as shown on the screen (in particular, **position** should be set to "chr22", **output format** should be set to "BED - browser extensible data", and "Galaxy" should be checked within the **Send output to** option). Click **get output** and you will see the next screen:

![UCSC ganes 2](http://galaxy.psu.edu/galaxy101/ucscGenes2.png)

here make sure **Create one BED record per:** is set to "Coding Exons" and click **Send Query to Galaxy** button. After this you will see your first History Item in Galaxy's right pane. It will go through gray (preparing) and yellow (running) states to become green:

![First history item](http://galaxy.psu.edu/galaxy101/firstHistoryItem.png)

### 1.1. Getting SNPs

Now is the time to obtain SNP data. This is done almost exactly the same way. First thing we will do is to again click on **Get Data -> UCSC Main**:

![get data from UCSC](http://galaxy.psu.edu/galaxy101/getDataUCSC.png)

but now change **group** to "Variation":

![Variation](http://galaxy.psu.edu/galaxy101/variation.png)

so that the whole page looks like this:

![UCSC SNPs](http://galaxy.psu.edu/galaxy101/ucscSNPs.png)

click get output and you should see this:

![UCSC SNPs 2](http://galaxy.psu.edu/galaxy101/ucscSNPs2.png)

where you need to make sure that **Whole Gene** is selected ("Whole Gene" here really means "Whole Feature") and click **Send Query to Galaxy** button. You will get your second item in the history:

![Second history item](http://galaxy.psu.edu/galaxy101/secondHistoryItem.png)

Now we will rename the two history items to "Exons" and "SNPs" by clicking on the Pencil icon adjacent to each item. After changing the name scroll down and click **Save**.  Also we will rename history to "Galaxy 101 (2015)" (or whatever you want) by clicking on **Unnamed history** so everything looks like this:

![Rename](http://galaxy.psu.edu/galaxy101/rename.png)

## 2. Finding Exons with the highest number of SNPs
### 2.0. Joining exons with SNPs
Let's remind ourselves that our objective was to find which exon contains the most SNPs. This first step in answering this question will be joining exons with SNPs (a fancy word for printing exons and SNPs that overlap side by side). This is done using **Operate on Genomics Intervals -> Join** tool:

![Join](http://galaxy.psu.edu/galaxy101/join.png)

make sure your **Exons** are first and **SNPs** are second and click **Execute**. You will get the third history item:

![Third history item](http://galaxy.psu.edu/galaxy101/thirdHistoryItem.png)

which will contain the following data (showing just one row):

`chr22	15528158	15529139	uc011agd.3_cds_0_0_chr22_15528159_f	0	+	chr22	15528374	15528375	rs567527834	0	-`

`chr22	15528158	15529139	uc011agd.3_cds_0_0_chr22_15528159_f	0	+	chr22	15528308	15528309	rs200358901	0	-`

`chr22	15528158	15529139	uc011agd.3_cds_0_0_chr22_15528159_f	0	+	chr22	15528266	15528267	rs545236550	0	-`

`chr22	15528158	15529139	uc011agd.3_cds_0_0_chr22_15528159_f	0	+	chr22	15528178	15528179	rs200562384	0	-`

Let's take a look at this dataset. The first six columns correspond to exons. The last six correspond to SNPs. You can see that exon with ID `uc011agd.3_cds_0_0_chr22_15528159_f` contains many (48) SNPs.

### 2.1. Counting the number of SNPs per exon
Above we've seen that exon `uc011agd.3_cds_0_0_chr22_15528159_f` is repeated many times (48) in the above dataset. Thus we can easily compute the number of SNPs per exon by simply counting the number of repetitions of name for each exon. This can be easily done with the **Join, Subtract, and Group -> Group** tool:

![Group](http://galaxy.psu.edu/galaxy101/group1.png)

choose column 4 by selecting "Column: 4" in **Group by** column. Then click on **Insert Operation** and make sure the interface looks exactly as shown below:

![Group2](http://galaxy.psu.edu/galaxy101/group2.png)

click **Execute**. Your history will look like this:
![Fourth history item](http://galaxy.psu.edu/galaxy101/fourthHistoryItem.png

if you look at the above image you will see that the result of grouping (dataset #4) contains two columns. This first contains the exon name while the second shows the number of times this name has been repeated in dataset #3. 

### 2.2. Sorting exons by SNP count
To see which exon has the highest number of SNPs we can simply sort the dataset #4 on the second column in descending order. This is done with **Filter and Sort -> Sort**:

![Sort](http://galaxy.psu.edu/galaxy101/sort.png)

This will generate the fifth history item:

![Fifth history item](http://galaxy.psu.edu/galaxy101/fifthHistoryItem.png)

and you can now see that the highest number of SNPs per exon is 63. 

### 2.3. Selecting top five
Now let's select top five exons with the highest number of SNPs. For this we will use **Text Manipulation -> Select First** tool:

![Select first](http://galaxy.psu.edu/galaxy101/selectFirst.png)

Clicking **Execute** will produce the sixth history item that will contain just five lines:

![Sixth history item](http://galaxy.psu.edu/galaxy101/sixthHistoryItem.png)

### 2.4. Recovering exon info and displaying data in genome browsers
Now we know that in this dataset the five top exons contain between 26 and 63 SNPs. But what else can we learn about these? To know more we need to get back the positional information (coordinates) of these exons. This information was lost at the grouping step and now all we have is just two columns. To get coordinates back we will match the names of exons in dataset #6 (column 1) against names of the exons in the original dataset #1 (column 4). This can be done with **Join, Subtract and Group -> Compare two Queries** tool (note the settings of the tool in the middle pane):

![Compare](http://galaxy.psu.edu/galaxy101/compare.png)

this adds the seventh dataset to the history:

![Seventh history item](http://galaxy.psu.edu/galaxy101/seventhHistoryItem.png)

The best way to learn about these exons is to look at their genomic surrounding. There is really no better way to do this than using genome browsers. Because this analysis was performed on "standard" human genome (`hg38` in this case), you have two choices - [UCSC Genome Browser](http://genomes.ucsc.edu) and [IGV](https://www.broadinstitute.org/igv/):

http://galaxy.psu.edu/galaxy101/browsers.png

For example, clicking on **display at UCSC main** will show something like this:

http://galaxy.psu.edu/galaxy101/ucsc.png

## 3. Understanding histories
In Galaxy your analyses live in histories such as this one:

http://galaxy.psu.edu/galaxy101/seventhHistoryItem.png

Histories can be very large, you can have as many histories as you want, and all history behavior is controlled by the ![refresh](http://galaxyproject.org/galaxy101/fa-refresh.png), ![cog](http://galaxyproject.org/galaxy101/fa-cog.png), and ![refresh](http://galaxyproject.org/galaxy101/fa-columns.png) buttons on the top of the History pane:

![History detail](http://galaxyproject.org/galaxy101/historyDetail.png)

The ![refresh](http://galaxyproject.org/galaxy101/fa-refresh.png) simply refreshes the history. The ![cog](http://galaxyproject.org/galaxy101/fa-cog.png) button gives you access to myriad of history-specific options:

![History options](http://galaxy.psu.edu/galaxy101/historyOptions.png)

Many of the options here are self-explanatory. If you create a new history, your current history does not disappear. If you would like to list all of your histories just choose **Saved Histories** and you will see a list of all your histories in the center pane:

![History list](http://galaxy.psu.edu/galaxy101/historyList.png)

Yet there is a better way to look for all your histories. This is what the ![refresh](http://galaxyproject.org/galaxy101/fa-columns.png) button is for. It allows you to browse and move datasets across histories:

![History list](http://galaxyproject.org/galaxy101/multiHistoryView.png)

Here, the current history is on the left (**Galaxy 101 (2015)**) and your (or mine in this case) other histories are displayed to the right of the current history. These can be ordered in a variety of ways by clicking the **...** button:

![History list](http://galaxyproject.org/galaxy101/orderingHistories.png)

You can also scroll sideways using trackpad gestures, move datasets across histories by simply clicking and dragging, and search for histories and individual datasets. This interface also allows you to switch to any existing history (i.e., making it current). Click **Done** once you're done.

## :point_right: Continue to [Galaxy 101-2](https://github.com/nekrut/galaxy/wiki/Galaxy101-2)
















