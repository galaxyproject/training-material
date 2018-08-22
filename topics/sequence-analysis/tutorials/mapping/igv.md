> ### {% icon hands_on %} Hands-on: Visualization of the reads in IGV
>
> 1. Install [IGV](https://software.broadinstitute.org/software/igv/download) (if not already installed)
> 1. Expand the {% icon param-file %} output of **{{ include.tool }}** {% icon tool %})
> 1. Click on the `local` in `display with IGV` to load the reads into the IGV browser
> 2. Zoom on the `{{ include.region_to_zoom }}`)
{: .hands_on}

The reads have a direction: they are mapped to the forward or reverse strand, respectively. When hovering over a read, extra information is displayed

> ### {% icon question %} Questions
>
> Some reads have colored lines included. What is this?
>
> > ### {% icon solution %} Solution
> > Try to zoom in in one of those lines and you will see the answer!
> {: .solution }
{: .question}

> ### {% icon comment %} Comments
> Because the number of reads over a region can be quite large, the IGV browser by default only allows to see the reads that fall into a small window. This behaviour can be changed in the IGV from `view > Preferences > Alignments`.
{: .comment}