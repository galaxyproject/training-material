---
layout: tutorial_hands_on
title: Mapping
zenodo_link: https://doi.org/10.5281/zenodo.1324070
questions:
  - Was ist Mapping?
  - Welche zwei Dinge sind entscheidend für ein korrektes Mapping?
  - Was ist eine BAM‑Datei?
objectives:
  - Ein Tool ausführen, um Reads auf ein Referenzgenom zu mappen
  - Erklären, was eine BAM‑Datei ist und was sie enthält
  - Einen Genom‑Browser benutzen, um die Daten zu verstehen
time_estimation: 1h
key_points:
  - Kenne deine Daten!
  - Mapping ist nicht trivial
  - Es gibt viele Mapping‑Algorithmen, abhängig von deinen Daten musst du den passenden wählen
subtopic: basics
requirements:
- type: internal
  topic_name: sequence-analysis
  tutorials:
  - quality-control
follow_up_training:
- type: internal
  topic_name: transcriptomics
  tutorials:
  - ref-based
- type: internal
  topic_name: epigenetics
  tutorials:
  - formation_of_super-structures_on_xi
level: Introductory
edam_ontology:
- topic_0102
contributions:
  authorship:
    - joachimwolff
    - bebatut
    - hexylena
  editing:
    - tflowers15
  translation:
    - Tillsa
    - unode
  funding:
    - unimelb
    - melbournebioinformatics
    - AustralianBioCommons
    - biont
recordings:
- captioners:
  - shiltemann
  date: '2021-02-15'
  galaxy_version: '21.01'
  length: 20M
  youtube_id: 1wm-62E2NkY
  speakers:
  - pvanheus
- youtube_id: B-CFWSkZzRc
  length: 24M
  galaxy_version: 24.1.2.dev0
  date: '2024-09-07'
  speakers:
  - DinithiRajapaksha
  captioners:
  - DinithiRajapaksha
  bot-timestamp: 1725707919
lang: de
tags:
  - english
  - español
  - italiano
translations:
  - en
  - es
  - it
---



Die Sequenzierung erzeugt eine Sammlung von Sequenzen ohne genomischen Kontext. Wir
wissen nicht, zu welchem Teil des Genoms die Sequenzen gehören. Das Mapping der Reads
eines Experiments auf ein Referenzgenom ist ein wichtiger Schritt in der modernen
Genomdatenanalyse. Durch das Mapping werden die Reads einer bestimmten Stelle im Genom
zugeordnet und es können Erkenntnisse wie das Expressionsniveau von Genen gewonnen
werden.

Die Reads sind nicht mit Positionsinformationen versehen, so dass wir nicht wissen, aus
welchem Teil des Genoms sie stammen. Wir müssen die Sequenz des Reads selbst verwenden,
um die entsprechende Region in der Referenzsequenz zu finden. Die Referenzsequenz kann
jedoch recht lang sein (~3 Milliarden Basen beim Menschen), was die Suche nach einer
passenden Region zu einer entmutigenden Aufgabe macht. Da unsere Reads kurz sind, kann
es mehrere, gleich wahrscheinliche Stellen in der Referenzsequenz geben, von denen sie
gelesen worden sein könnten. Dies gilt insbesondere für sich wiederholende Regionen.

Im Prinzip könnten wir eine BLAST-Analyse durchführen, um herauszufinden, wo die
sequenzierten Teile am besten in das bekannte Genom passen. Das müssten wir für jede der
Millionen von Reads in unseren Sequenzierdaten tun. Das Alignment von Millionen kurzer
Sequenzen auf diese Weise kann jedoch einige Wochen dauern. Und wir interessieren uns
nicht für die genaue Übereinstimmung der Basen (Alignment). Was uns interessiert, ist,
"woher diese Reads stammen". Dieser Ansatz wird **Mapping** genannt.

Im Folgenden werden wir einen Datensatz mit dem Mapper **Bowtie2** bearbeiten und die
Daten mit dem Programm **IGV** visualisieren.

> <agenda-title></agenda-title>
> 
> In diesem Tutorial werden wir uns mit folgenden Themen beschäftigen:
> 
> 1. TOC
> {:toc}
> 
{: .agenda}

# Vorbereiten der Daten

> <hands-on-title>Daten-Upload</hands-on-title>
> 
> 1. Erstellen Sie einen neuen Verlauf für dieses Tutorial und geben Sie ihm einen
>    passenden Namen
> 
>    {% snippet faqs/galaxy-de/histories_create_new.md %}
> 
>    {% snippet faqs/galaxy-de/histories_rename.md %}
> 
> 2. Importieren Sie `wt_H3K4me3_read1.fastq.gz` und `wt_H3K4me3_read2.fastq.gz` von
>    [Zenodo](https://zenodo.org/record/1324070) oder aus der Datenbibliothek (fragen
>    Sie Ihren Dozenten)
> 
>    ```
>    https://zenodo.org/record/1324070/files/wt_H3K4me3_read1.fastq.gz
>    https://zenodo.org/record/1324070/files/wt_H3K4me3_read2.fastq.gz
>    ```
> 
>    {% snippet faqs/galaxy-de/datasets_import_via_link.md %}
> 
>    {% snippet faqs/galaxy-de/datasets_import_from_data_library.md %}
> 
>    Standardmäßig nimmt Galaxy den Link als Namen, also benennen Sie sie um.
> 
> 3. Benennen Sie die Dateien in `reads_1` und `reads_2` um
> 
>    {% snippet faqs/galaxy-de/datasets_rename.md %}
>
> 4. Erstelle eine gepaarte Sammlung namens `Paired Reads`
>
>    {% snippet faqs/galaxy-de/collections_build_list_paired.md %}
> 
{: .hands_on}

Wir haben in Galaxy einfach FASTQ-Dateien importiert, die den Paired-End-Daten
entsprechen, die wir direkt von einer Sequenziereinrichtung erhalten haben.

Bei der Sequenzierung können Fehler auftreten, wie z. B. der Aufruf falscher Nukleotide.
Sequenzierungsfehler können die Analyse verfälschen und zu einer Fehlinterpretation der
Daten führen. Der erste Schritt bei jeder Art von Sequenzierungsdaten ist immer die
Überprüfung ihrer Qualität.

Es gibt ein spezielles Tutorial zur [Qualitätskontrolle]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}) von Sequenzierungsdaten. Wir werden die Schritte dort nicht wiederholen. Sie sollten das [tutorial]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}) befolgen und es auf Ihre Daten anwenden, bevor Sie weitermachen.

# Mapping der Reads auf ein Referenzgenom

Beim Read-Mapping werden die Reads an ein Referenzgenom angeglichen. Ein Mapper nimmt
als Eingabe ein Referenzgenom und einen Satz von Reads. Sein Ziel ist es, jeden Read in
der Menge der Reads am Referenzgenom auszurichten, wobei Mismatches, Indels und das
Abschneiden einiger kurzer Fragmente an den beiden Enden der Reads berücksichtigt
werden:

![Erläuterung des Mappings](../../images/mapping/mapping.png "Illustration des Mapping-Prozesses. Die Eingabe besteht aus einer Reihe von Reads und einem Referenzgenom. In der Mitte sind die Ergebnisse des Mappings dargestellt: die Positionen der Reads auf dem Referenzgenom. Der erste Read wird an Position 100 ausgerichtet und das Alignment weist zwei Mismatches auf. Der zweite Read wird an der Position 114 ausgerichtet. Es handelt sich um ein lokales Alignment mit Ausschnitten auf der linken und rechten Seite. Der dritte Read wird an der Position 123 ausgerichtet. Er besteht aus einer 2-Basen-Insertion und einer 1-Base-Deletion.")

Wir benötigen ein Referenzgenom, auf das wir die Reads mappen können.

{% include topics/sequence-analysis/tutorials/mapping/ref_genome_explanation_DE.md answer_3="Diese Daten stammen aus der ChIP-seq von Mäusen, daher werden wir mm10 (*Mus musculus*) verwenden."%}

Derzeit gibt es über 60 verschiedene Mapper, und ihre Zahl wächst. In diesem Tutorial
werden wir [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) verwenden, ein
schnelles und speichereffizientes Open-Source-Tool, das sich besonders gut für das
Alignment von Sequenzierungs-Reads von etwa 50 bis zu 1.000 Basen zu relativ langen
Genomen eignet.

> <hands-on-title>Mapping mit Bowtie2</hands-on-title>
> 1. {% tool [Bowtie2](toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.4.2+galaxy0) %} mit den folgenden Parametern
>    - *"Ist dies eine einzelne oder gepaarte Bibliothek "*: `Paired-end`
>       - {% icon param-file %} *"FASTA/Q-Datei #1"*: `reads_1`
>       - {% icon param-file %} *"FASTA/Q-Datei #2"*: `reads_2`
>       - *"Möchten Sie Paired-End-Optionen festlegen? "*: `No`
> 
>         Sie sollten sich die Parameter dort ansehen, insbesondere die Paarungsorientierung, wenn Sie sie kennen. Sie können die Qualität des Paired-End-Mappings verbessern.
> 
>     - *"Werden Sie ein Referenzgenom aus Ihrer Historie auswählen oder einen eingebauten Index verwenden? "*: `Use a built-in genome index`
>       - *"Referenzgenom auswählen "*: `Mouse (Mus musculus): mm10`
>     - *"Analysemodus auswählen "*: `Default setting only`
> 
>       Sie sollten einen Blick auf die nicht-standard Parameter werfen und versuchen, sie zu verstehen. Sie können einen Einfluss auf das Mapping haben und es verbessern.
> 
>     - *"Speichern Sie die Bowtie2-Mapping-Statistiken in der Historie "*: `Yes`
> 
> 2. Untersuchen Sie die Datei `mapping stats`, indem Sie auf das Symbol {% icon galaxy-eye %} (Auge) Symbol
> 
{: .hands_on}

> <question-title></question-title>
> 
> 1. Welche Informationen werden hier bereitgestellt?
> 2. Wie viele Reads wurden genau 1 Mal gemappt?
> 3. Wie viele Reads wurden mehr als 1 Mal gemappt? Wie ist das möglich? Was sollten wir mit ihnen machen?
> 4. Wie viele Paare von Reads wurden nicht gemappt? Was sind die Ursachen dafür?
> 
> > <solution-title></solution-title>
> > 1. Die hier gegebene Information ist eine quantitative. Wir können sehen, wie viele Sequenzen aufeinander abgestimmt sind. Sie sagt nichts über die Qualität aus.
> > 2. ~90% der Reads wurden genau 1 Mal aligniert
> > 3. ~7 % der Reads wurden >1 Mal übereinstimmend ausgerichtet. Diese werden als
> >    "multi-mapped reads" bezeichnet. Dies kann aufgrund von Wiederholungen im
> >    Referenzgenom geschehen (z. B. mehrere Kopien eines Gens), insbesondere wenn die
> >    Reads klein sind. Es ist schwierig zu entscheiden, woher diese Sequenzen stammen,
> >    und deshalb werden sie von den meisten Pipelines ignoriert. Überprüfen Sie immer
> >    die Statistiken, um sicherzustellen, dass nicht zu viele Informationen in
> >    nachfolgenden Analysen verworfen werden.
> > 4. ~3% der Reads wurden nicht gemappt, weil
> >     - beide Reads des Paares sind aligned, aber ihre Positionen stimmen nicht
> >       mit dem Paar von Reads überein (`aligned discordantly 1 time`)
> >     - Reads dieser Paare sind mehrfach gemappt (`aligned >1 times` in `pairs aligned
> >       0 times concordantly or discordantly`)
> >     - ein Read dieser Paare wird gemappt, aber nicht der gepaarte Read (`aligned
> >       exactly 1 time` in `pairs aligned 0 times concordantly or discordantly`)
> >     - der Rest ist überhaupt nicht gemappt
> {: .solution }
{: .question}

Die Überprüfung der Mapping-Statistiken ist ein wichtiger Schritt, der vor der
Fortsetzung der Analysen durchgeführt werden muss. Es gibt mehrere potenzielle
Fehlerquellen beim Mapping, einschließlich (aber nicht beschränkt auf):

- **Polymerase-Kettenreaktion (PCR)-Artefakte**: Viele
  Hochdurchsatz-Sequenzierungsmethoden (HTS) beinhalten einen oder mehrere PCR-Schritte.
  PCR-Fehler zeigen sich als Mismatches im Alignment, und insbesondere Fehler in frühen
  PCR-Zyklen zeigen sich in mehreren Reads, was fälschlicherweise auf eine genetische
  Variation in der Probe schließen lässt. Ein ähnlicher Fehler sind PCR-Duplikate, bei
  denen dasselbe Lesepaar mehrfach vorkommt und die Berechnung der Abdeckung im
  Alignment verfälscht.
- **Sequenzierungsfehler**: Das Sequenziergerät kann entweder aus physikalischen Gründen
  (z. B. Öl auf einem Illumina-Objektträger) oder aufgrund von Eigenschaften der
  sequenzierten DNA (z. B. Homopolymere) einen fehlerhafte Aussage machen. Da
  Sequenzierfehler oft zufällig sind, können sie beim Variantenaufruf als Singleton
  Reads herausgefiltert werden.
- **Mapping-Fehler**: Der Mapping-Algorithmus kann einen Read an der falschen Stelle in
  der Referenz zuordnen. Dies geschieht häufig im Bereich von Wiederholungen oder
  anderen Regionen mit geringer Komplexität.

Wenn also die Mapping-Statistiken nicht gut sind, sollten Sie die Ursache für diese
Fehler untersuchen, bevor Sie mit Ihren Analysen fortfahren.

Danach sollten Sie einen Blick auf die Reads werfen und die BAM-Datei inspizieren, in
der die Read-Mappings gespeichert sind.

# Inspektion einer BAM-Datei

{% include topics/sequence-analysis/tutorials/mapping/bam_explanation_DE.md mapper="Bowtie2" %}

Die BAM-Datei enthält viele Informationen über jeden Read, insbesondere über die Qualität des Mappings.

> <hands-on-title>Zusammenfassung der Mapping-Qualität</hands-on-title>
> 1. {% tool [Samtools Stats](toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.2+galaxy2) %} mit den folgenden Parametern
>    - {% icon param-file %} *"BAM-Datei "*: `aligned reads` (Ausgabe von **Bowtie2** {% icon tool %})
>    - *"Referenzsequenz verwenden "*: `Locally cached/Use a built-in genome`
>      - *"Genom verwenden "*: `Mouse (Mus musculus): mm10 Full`
> 
> 2. Untersuchen Sie die {% icon param-file %} `Stats`-Datei
> 
{: .hands_on}

> <question-title></question-title>
> 
> 1. Wie hoch ist der Anteil der Mismatches in den gemappten Reads, wenn sie an das Referenzgenom angeglichen werden?
> 2. Was bedeutet die Fehlerrate?
> 3. Was ist die durchschnittliche Qualität? Wie wird sie dargestellt?
> 4. Was ist die durchschnittliche Insertgröße?
> 5. Wie viele Reads haben einen Mapping-Qualitätsscore unter 20?
> 
> > <solution-title></solution-title>
> > 1. Es gibt ~21.900 Mismatches für ~4.753.900 gemappted Basen, was im Durchschnitt
> >    ~0,005 Mismatches pro gemappter Base ergibt.
> > 2. Die Fehlerrate ist der Anteil der Fehlanpassungen pro gemappter Base, also das
> >    unmittelbar zuvor berechnete Verhältnis.
> > 3. Die durchschnittliche Qualität ist der mittlere Qualitätswert der Kartierung. Es
> >    handelt sich um einen Phred-Score, wie er auch in der FASTQ-Datei für jedes
> >    Nukleotid verwendet wird. Hier ist der Score jedoch nicht pro Nukleotid, sondern
> >    pro Read und stellt die Wahrscheinlichkeit der Mapping-Qualität dar.
> > 4. Die Insertgröße ist der Abstand zwischen den beiden Reads in den Paaren.
> > 5. Um die Informationen zu erhalten:
> >      1. {% tool [Filter BAM](toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.5.2+galaxy2) %} mit einem Filter, um nur die Reads mit einer Mappingqualität >= 20 zu behalten
> >      2. {% tool [Samtools Stats](toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.5) %} auf die Ausgabe von **Filter**
> > 
> >    Vor der Filterung: 95.412 Reads und nach der Filterung: 89.664 Reads.
> {: .solution }
{: .question}

# Visualisierung mit einem Genom-Browser

## IGV

Der Integrative Genomics Viewer (IGV) ist ein hochleistungsfähiges Visualisierungstool
für die interaktive Erkundung großer, integrierter Genomdatensätze. Es unterstützt eine
Vielzahl von Datentypen, einschließlich array-basierter und Next-Generation-Sequenzdaten
sowie genomische Annotationen. Im Folgenden werden wir es verwenden, um die gemappten
Reads zu visualisieren.

{% include topics/sequence-analysis/tutorials/mapping/igv_DE.md tool="Bowtie2" region_to_zoom="chr2:98,666,236-98,667,473" %}

## JBrowse

{% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy0) %}
ist ein alternativer, webbasierter Genombrowser. Während IGV eine Software ist, die Sie
herunterladen und ausführen müssen, sind JBrowse-Instanzen online gehostete Websites,
die eine Schnittstelle zum Durchsuchen von Genomdaten bieten. Wir werden es verwenden,
um die gemappten Reads zu visualisieren.

{% include topics/sequence-analysis/tutorials/mapping/jbrowse_DE.md tool="Bowtie2" region_to_zoom="chr2:98,666,236-98,667,473" %}

# Schlussfolgerung

Nach der Qualitätskontrolle ist das Mapping ein wichtiger Schritt bei den meisten
Analysen von Sequenzierungsdaten (RNA-Seq, ChIP-Seq usw.), um festzustellen, wo im Genom
unsere Reads herkommen, und diese Information für nachgeschaltete Analysen zu nutzen.

