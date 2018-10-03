---
layout: tutorial_hands_on
topic_name: genome-annotation
tutorial_name: annotation-in-apollo
---
> ### Agenda
>
> 1. Prerequisites
> 2. Making an Annotation
> 3. Making the Best Prediction
>
{: .agenda}

# Prerequisites

Before beginning annotation within Galaxy ([CPT Public Galaxy](https://cpt.tamu.edu/galaxy-pub), [CPT TAMU Galaxy](https://cpt.tamu.edu/galaxy)), it is necessary that there is a genome loaded into Apollo. Additionally, the structural annotation *must* be complete, and the functional annotation workflow has **already been run**. The functional annotation workflow opens up the necessary evidence tracks for annotation.

<!-- LINK STRUCTURAL AND FUNCTIONAL ANNOTATION TUTORIALS UPON COMPLETION. -->

# Making an Annotation

Generally, the annotation process is a synthesis between the understanding of phage genomics and the available evidence tracks. The [Center for Phage Technology](https://cpt.tamu.edu) encourages **phage** annotators on the CPT’s Apollo instance to follow some specific conventions (Field -> *Recommended Input*):

> * Name -> *Gene name* (Could be something like **terminase small subunit** or **hypothetical protein**.)
> * Symbol -> *Do not use.*
> * Description -> *Do not use.*
> * DBXRefs -> *Only use if the annotator is experienced; please ensure formatting is correct.*
> * Attributes -> *Do not use.*
> * PubMed IDs -> *Do not use.*
> * Gene Ontology IDs -> *Do not use.*
> * Comments -> *Apply any free-text comments here.* (Could be something like **the e-value(s)** between the annotated gene and homologs or notes to one’s self.)

> ### {% icon comment %} Naming Guidelines
> It is imperative to follow suit with the [UniProt](file:///Users/cptaggies/Downloads/International_Protein_Nomenclature_Guidelines%20(1).pdf) and [NCBI](https://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/) international naming conventions. It allows for standardization and consistency in naming proteins, subsequently aiding data retrieval and improving communication.
{: .comment}

> ### {% icon tip %} Note that…
> Calling genes is [covered in another tutorial](LINK TO TUTORIAL)
{: .tip}

<!-- LINK TUTORIAL UPON COMPLETION -->

# Making the Best Prediction

The [Center for Phage Technology](https://cpt.tamu.edu/) integrates as many data sources as is feasible. Please contact the CPT IT staff (cpt@tamu.edu) if there is another data source needed that is not appearing not currently available in Apollo
