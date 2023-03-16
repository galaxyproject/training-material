> <question-title></question-title>
>
> What is the relation between `{{ include.forward }}` and `{{ include.reverse }}` ?
>
> > <solution-title></solution-title>
> > The data has been sequenced using paired-end sequencing. 
> >
> > The paired-end sequencing is based on the idea that the initial DNA fragments (longer than the actual read length) is sequenced from both sides. This approach results in two reads per fragment, with the first read in forward orientation and the second read in reverse-complement orientation. The distance between both reads is known. Thus, it can be used as an additional piece of information to improve the read mapping.
> > 
> > With paired-end sequencing, each fragment is more covered than with single-end sequencing (only forward orientation sequenced):
> >
> > ```
> >     ------>
> >       ------>
> >         ------>
> >           ------>
> >
> >     ----------------------------- [fragment]
> >
> >     ------>         <------
> >       ------>         <------
> >         ------>         <------
> >           ------>         <------
> > ```
> > 
> > The paired-end sequencing generates then 2 files:
> > - One file with the sequences corresponding to forward orientation of all the fragments
> > - One file with the sequences corresponding to reverse orientation of all the fragments
> >
> > Here `{{ include.forward }}` corresponds to the forward reads and `{{ include.reverse }}` to the reverse reads.
> {: .solution }
{: .question}
