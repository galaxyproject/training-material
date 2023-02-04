> <hands-on-title> Extract input files </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `AnnData_filtered`
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
> 2. Rename {% icon galaxy-pencil %} the observations annotation `Cell metadata (obs)`
>
> 3. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `AnnData_filtered`
>    - *"What to inspect?"*: `Key-indexed annotation of variables/features (var)`
> 4. Rename {% icon galaxy-pencil %} the annotation of variables `Gene metadata (var)`
>
> 5. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `AnnData_filtered`
>    - *"What to inspect?"*: `The full data matrix`
> 6. Rename {% icon galaxy-pencil %} the output `Expression matrix`
>
> 7. Download the generated files from your history. To do so, just click on the {% icon galaxy-save %} save icon for `Cell metadata (obs)`, `Gene metadata (var)` and `Expression matrix`. We will need those later!
>
{: .hands_on}
