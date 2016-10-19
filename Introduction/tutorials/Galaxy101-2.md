# Galaxy 101-2: Getting to know workflows

## 4. Creating and editing a workflow
### 4.0. Extracting a workflow
Lets take a look at the history again:

![Collapsed history](http://galaxy.psu.edu/galaxy101/historyCollapsed.png)

You can see that this history contains all steps of our analysis. So by building this history we have actually created a complete record of our analysis with Galaxy preserving all parameter settings applied at every step. Wouldn't it be nice to just convert this history into a workflow that we'll be able to execute again and again? This can be done by clicking on the ![cog](http://galaxyproject.org/galaxy101/fa-cog.png) button and selecting **Extract Workflow** option:

![Extract workflow](http://galaxy.psu.edu/galaxy101/extractWorkflow.png)

The center pane will change as shown below and you will be able to choose which steps to include/exclude and how to name the newly created workflow. In this case I named it `galaxy101-2015`:

![Create workflow](http://galaxy.psu.edu/galaxy101/createWorkflow.png)

once you click **Create Workflow** you will get the following message: "Workflow 'galaxy101-2015' created from current history. You can **edit** or **run** the workflow". 

### 4.1. Opening workflow editor

Let's click **edit** (if you click something else and the message in the center page disappears, you can always access all your workflows including the one you just created using the **Workflow** link on top of Galaxy interface). This will open Galaxy's workflow editor (to get this view I clicked the arrow at the lower left corner of the screen, which collapsed the tool pane of the Galaxy interface). It will allow you to examine and change settings of this workflow as shown below. Note that the box corresponding to the **Select First** tool is selected (highlighted with the blueish border) and you can see parameters of this tool on the right pane. This is how you can view and change parameters of all tools involved in the workflow:

![Workflow editor](http://galaxy.psu.edu/galaxy101/wfEditor.png)

### 4.2. Hiding intermediate steps
Among multiple things you can do with workflows I will just mention one. When workflow is executed one is usually interested in the final product and not in the intermediate steps. These steps can be hidden by mousing over a small asterisk in the lower right corner of every tool:

![Hide step](http://galaxy.psu.edu/galaxy101/hideStep.png)

Yet there is a catch. In a newly created workflow all steps are hidden by default and the default behavior of Galaxy is that if all steps of a given workflow are hidden, then nothing gets hidden in the history. This may be counterintuitive, but this is done to decrease the amount of clicking if you do want to hide some steps. So in our case if we want to hide all intermediate steps with the exception of the last one we will click that asterisk in last step of the workflow:

![Last step](http://galaxy.psu.edu/galaxy101/lastStep.png)

Once you do this the representation of the workflow in the bottom right corner of the editor will change with the last step becoming orange. This means that this is the only step, which will generate a dataset visible in the history:

![Workflow Overview](http://galaxy.psu.edu/galaxy101/workflowOverview.png)

### 4.3. Renaming inputs
Right now both inputs to the workflow look exactly the same. This is a problem as will be very confusing which input should be **Exons** and which should be **SNPs**:

![Naming inputs 1](http://galaxy.psu.edu/galaxy101/namingInputs1.png)

One the image above you will see that the top input dataset (the one with the blue border) connects to the **Join tool** first, so it must correspond to the exon data. If you click on this box (in the image above it is already clicked on because it is outlined with the blue border) you will be able to rename the dataset in the right pane:

![Rename Exons](http://galaxy.psu.edu/galaxy101/renameExons.png)

Then click on the second input dataset and rename it "Features" (this would make this workflow a bit more generic, which will be useful later in this tutorial):

![Rename features](http://galaxy.psu.edu/galaxy101/renameFeatures.png)

### 4.4. Renaming outputs
Finally let's rename the workflow's output. For this:

* click on the last dataset (**Compare two Queries**)
* scroll down the rightmost pane and click on ![add action](http://galaxy.psu.edu/galaxy101/addAction.png)
* Type `Top Exons` in the **Rename dataset** text box:

![Top exons](http://galaxy.psu.edu/galaxy101/topExons.png)

### 4.5. Setting parameters "at runtime"
What we are trying to do here is do design a generic workflow. This means that time to time you will need to change parameters within this workflow. For instance, in this tutorial we were selecting 5 exons containing the highest number of SNPs. But what if you need to select 10? Thus it makes sense to leave these types of parameters adjustable. To do this: 

First, select a tool in which you want to set parameters at runtime (`Select first` in this case):

![runtime Tool Selection](http://galaxy.psu.edu/galaxy101/runtimeTool.png)

Next, select parameter you would like to set at runtime. To do this just hover over the ![paramSymbol](http://galaxy.psu.edu/galaxy101/paramSymbol.png) icon so it looks like this:

![runtime Param Selection](http://galaxy.psu.edu/galaxy101/runtimeParam.png)

and click! Your parameter will now be set at runtime.

### 4.6. Save! It is important...
Now let's save the changes we've made by clicking ![cog](http://galaxyproject.org/galaxy101/fa-cog.png) and selecting Save:

![wfSave](http://galaxy.psu.edu/galaxy101/wfSave.png)

## 5. Run workflow on whole genome data
Now that we have a workflow, let's do something grand like, for example, finding exons with the highest number of repetitive elements across the entire human genome. 

### 5.0. Create a new history

First go back into analysis view by clicking **Analyze Data** on top of the Galaxy's interface. Now let's create a new history by clicking ![cog](http://galaxyproject.org/galaxy101/fa-cog.png) and selecting **Create New**:

![Create new history](http://galaxy.psu.edu/galaxy101/createNewHistory.png)

### 5.1. Get Exons
Now let's get coding exons for the entire genome by going to **Get Data -> UCSC Main** and setting up parameters as shown below. Note that this time `region` radio button is set to **genome**:

![Get all genes](http://galaxy.psu.edu/galaxy101/getAllGenes.png)

Click **get output** and you will get the next page (if it looks different from the image below, go back and make sure `output format` is set to **BED - browser extensible format**):

![Get exons](http://galaxy.psu.edu/galaxy101/ucscGenes2.png)

Choose **Coding exons** and click **Send query to Galaxy**.

## 5.2. Get Repeats
Go again to **Get Data -> UCSC Main** and make sure the following settings are selected (in particular `group` = **Repeats** and `track` = **RepeatMasker**):

![All repeats](http://galaxy.psu.edu/galaxy101/allRepeats.png)

Click **get output** and you will get the next page (if it looks different from the image below, go back and make sure `output format` is set to **BED - browser extensible format**):

![All repeats 2](http://galaxy.psu.edu/galaxy101/allRepeats2.png)

Select **Whole gene** and click **Send Query to Galaxy**.

### 5.3. Start the Workflow
At this point you will have two items in your history - one with exons and one with repeats. These datasets are large (especially repeats) and it will take some time for them to become green. Luckily you do not have to wait as Galaxy will automatically start jobs once uploads have ended. So nothing stops us from starting the workflow we have created. First, click on the **Workflow link** at the top of Galaxy interface, mouse over **galaxy101-2015**, click, and select **Run**. Center pane will change to allow you launching the workflow. Select appropriate datasets for `Repeats` and `Exon` inputs as shown below. Now scroll to **Step 6** and will see that we can set up `Select first` parameter at *Runtime* (meaning Now!). So lets put `20` in there (or anything else you want) and scroll further down to click ![Run workflow](http://galaxy.psu.edu/galaxy101/runWorkflowButton.png) to see this:

![Launch workflow](http://galaxy.psu.edu/galaxy101/launchWorkflow.png)

Once workflow has started you will initially be able to see all its steps. Note that you are joining all exons with all repeats, so naturally this will take some time:

![Launched workflow](http://galaxy.psu.edu/galaxy101/launchedWorkflow.png)

### 5.4. Get coffee
As we mentioned above this will take some time, so go get coffee. At last you will see this:

![Final view](http://galaxy.psu.edu/galaxy101/final.png)

## 6. We did not fake this:
The two histories and the workflow described in this page are accessible directly from this page below:

* History [**Galaxy 101 (2015)**](https://usegalaxy.org/u/aun1/h/galaxy-101-2015)
* History [**Exons vs. Repeats**]( https://usegalaxy.org/u/aun1/h/exons-vs-repeats-2015-1)
* Workflow [**Galaxy 101-2015**](https://usegalaxy.org/u/aun1/w/galaxy101-2015-1)

From there you can import histories and workflows to make them your own. For example, to import **Galaxy 101 (2015)** history simply click [this link](https://usegalaxy.org/u/aun1/h/galaxy-101-2015) and select `Import history` link:

![Final view](http://galaxy.psu.edu/galaxy101/importHistory.png)

## 7. If things don't work...
...you need to complain. Use [Galaxy's BioStar Channel](https://usegalaxy.org/biostar/biostar_redirect) to do this. 

## :clap: The End