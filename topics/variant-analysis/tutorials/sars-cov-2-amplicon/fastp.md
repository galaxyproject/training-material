##  Read quality control
Quality control and preprocessing of FASTQ files are essential to providing clean data for downstream analysis. There are many options for quality trimming, for instance [fastp](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234), [trimmomatic](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096) and [trim galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). We will use {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/0.20.1+galaxy0) %} which is an ultra-fast all-in-one FASTQ preprocessor including quality control, adapter trimming and quality filtering.

> #### {% icon hands_on %} Read trimming with `fastp`
> 1. Find the {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/0.20.1+galaxy0) %} tool in the tool bar and click on it. 
> 2. Select *Paired Collection* input
> 3. Your collection (list of pairs) of FASTQ reads should be selected automatically as the input.
> 3. Leave all the other options at their default values and click *Execute* icon to run the tool
>    Once you click *Execute*, new entries are added to your history. Take note of their colors: 
>      * grey indicates the execution is waiting
>    * yellow means the job is running
>    * red means that the job failed (there was an error)
>    * green indicates the run completed successfully.
{: .hands_on}