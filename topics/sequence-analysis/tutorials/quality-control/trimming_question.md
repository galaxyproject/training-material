>
>      > <question-title></question-title>
>      >
>      > Why do we run the trimming tool only once on a paired-end dataset and not twice, once for each dataset?
>      >
>      > > <solution-title></solution-title>
>      > >
>      > > The tool can remove sequences if they become too short during the trimming process. For paired-end files it removes entire sequence pairs if one (or both) of the two reads became shorter than the set length cutoff. Reads of a read-pair that are longer than a given threshold but for which the partner read has become too short can optionally be written out to single-end files. This ensures that the information of a read pair is not lost entirely if only one read is of good quality.
>      > >
>      > {: .solution}
>      {: .question}
>