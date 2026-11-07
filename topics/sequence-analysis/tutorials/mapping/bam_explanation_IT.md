Un file BAM ([Binary Alignment Map](https://en.wikipedia.org/wiki/SAM_(file_format)) è un file binario compresso che memorizza le sequenze di letture, se sono state allineate a una sequenza di riferimento (ad esempio un cromosoma) e, in caso affermativo, la posizione sulla sequenza di riferimento in cui sono state allineate.

> <hands-on-title>Ispettare un file BAM/SAM</hands-on-title>
> 
> 1. Ispezionare il {% icon param-file %} output di **{{ include.mapper }}** {% icon tool %}
> 
{: .hands_on}

Un file BAM (o un file SAM, la versione non compressa) consiste in:

- Una sezione di intestazione (le righe che iniziano con `@`) contenente metadati in particolare i nomi dei cromosomi e le lunghezze (righe che iniziano con il simbolo `@SQ`)
- Una sezione di allineamento costituita da una tabella con 11 campi obbligatori e un numero variabile di campi opzionali:

  | Col | Field | Type    | Brief Description                     |
  | --- | ----- | ------- | ------------------------------------- |
  | 1   | QNAME | String  | Query template NAME                   |
  | 2   | FLAG  | Integer | Bitwise FLAG                          |
  | 3   | RNAME | String  | References sequence NAME              |
  | 4   | POS   | Integer | 1- based leftmost mapping POSition    |
  | 5   | MAPQ  | Integer | MAPping Quality                       |
  | 6   | CIGAR | String  | CIGAR String                          |
  | 7   | RNEXT | String  | Ref. name of the mate/next read       |
  | 8   | PNEXT | Integer | Position of the mate/next read        |
  | 9   | TLEN  | Integer | Observed Template LENgth              |
  | 10  | SEQ   | String  | Segment SEQuence                      |
  | 11  | QUAL  | String  | ASCII of Phred-scaled base QUALity+33 |

> <question-title></question-title>
> 
> 1. Quali informazioni si trovano in un file SAM/BAM?
> 2. Quali sono le informazioni aggiuntive rispetto a un file FASTQ?
> 
> > <solution-title></solution-title>
> > 1. Sequenze e informazioni sulla qualità, come un FASTQ
> > 2. Informazioni sulla mappatura, posizione della lettura sul cromosoma, qualità della mappatura, ecc
> {: .solution }
{: .question}
