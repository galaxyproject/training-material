> ### {% icon hands_on %} Hands-on: Visualization of the reads in IGV
>
> 1. Install [IGV](https://software.broadinstitute.org/software/igv/download) (if not already installed)
> 2. Launch IGV on your computer
> 2. Expand the {% icon param-file %} output of **{{ include.tool }}** {% icon tool %}
> 3. Click on the `local` in `display with IGV` to load the reads into the IGV browser
> 4. Zoom on the `{{ include.region_to_zoom }}`
{: .hands_on}

The reads have a direction: they are mapped to the forward or reverse strand, respectively. When hovering over a read, extra information is displayed

> ### {% icon question %} Questions
> 
> 1. What could it mean if a bar in the coverage view is colored?
> 2. What could be the reason why a read is white instead of grey?
> 
> > ### {% icon solution %} Solution
> > 1. If a nucleotide differs from the reference sequence in more than 20% of quality weighted reads, IGV colors the bar in proportion to the read count of each base.
> > 2. They have a mapping quality equal to zero. Interpretation of this mapping quality depends on the mapping aligner as some commonly used aligners use this convention to mark a read with multiple alignments. In such a case, the read also maps to another location with equally good placement. It is also possible the read could not be uniquely placed but the other placements do not necessarily give equally good quality hits.
> {: .solution }
{: .question}

> ### {% icon tip %} Tips for IGV
> 1. Because the number of reads over a region can be quite large, the IGV browser by default only allows to see the reads that fall into a small window. This behaviour can be changed in the IGV from `view > Preferences > Alignments`.
> 2. If the genome of your interest is not there check if it is available via **More...**. If this is not the case, you can add it manually via the menu **Genomes** -> **Load Genome from...**
>
>    ![Select genome in IGV](../../images/igv_select_genome.png "Select genome in IGV")
>
> A general description of the user interface of the IGV browser is available here: [IGV Browser description]({{site.baseurl}}/topics/introduction/tutorials/igv-introduction/tutorial.html)
{: .tip}
