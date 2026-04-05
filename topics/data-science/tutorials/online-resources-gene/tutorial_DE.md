---
layout: tutorial_hands_on
title: Lernen über ein Gen über biologische Ressourcen und Formate hinweg
level: Introductory
draft: true
zenodo_link: https://zenodo.org/record/8304465
questions:
  - Wie kann man bioinformatische Ressourcen nutzen, um eine bestimmte Proteinfamilie (Opsine) zu untersuchen?
  - Wie navigiert man im Genome Data Viewer, um Opsine im menschlichen Genom zu finden?
  - Wie identifiziert man Gene, die mit Opsinen assoziiert sind, und analysiert ihre Chromosomenlage?
  - Wie erkundet man Literatur und klinische Kontexte für das Gen OPN1LW?
  - Wie nutzt man Proteinsequenzdateien für Ähnlichkeitssuchen mittels BLAST?
objectives:
  - Ausgehend von einer Textsuche mehrere Webressourcen anzusteuern, um verschiedene Arten von Informationen über ein Gen zu untersuchen, dargestellt in verschiedenen Dateiformaten.
time_estimation: 1H
key_points:
  - Sie können nach Genen und Proteinen mit bestimmten Textbegriffen im NCBI‑Genom suchen.
  - Sobald Sie ein relevantes Gen oder Protein gefunden haben, können Sie dessen Sequenz und Annotation in verschiedenen Formaten von NCBI erhalten.
  - Sie können auch mehr über die Chromosomenlage und die Exon‑Intron‑Zusammensetzung des interessierenden Gens erfahren.
  - NCBI bietet ein BLAST‑Tool, um Ähnlichkeitssuchen mit Sequenzen durchzuführen.
  - Sie können die in diesem Tutorial enthaltenen Ressourcen weiter erkunden, um mehr über genassoziierte Zustände und Varianten zu lernen.
  - Sie können eine FASTA‑Datei mit einer interessierenden Sequenz für BLAST‑Suchen eingeben.
contributions:
  authorship:
  - lisanna
  - bebatut
  - teresa-m
  translation:
  - Tillsa
  - unode
  funding:
  - biont
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


Wenn wir eine bioinformatische Analyse durchführen, z. B. RNA-seq, erhalten wir möglicherweise eine Liste von Gennamen. Dann müssen wir diese Gene untersuchen. Aber wie können wir das tun? Welche Ressourcen stehen dafür zur Verfügung? Und wie kann man durch sie navigieren?

Ziel dieses Tutoriums ist es, sich damit am Beispiel der menschlichen Opsine vertraut zu machen.

Menschliche Opsine befinden sich in den Zellen der Netzhaut. Opsine fangen Licht ein und setzen die Signalfolge in Gang, die zum Sehen führt. Wir werden Fragen über Opsine und Opsin-Gene stellen und dann verschiedene bioinformatische Datenbanken und Ressourcen nutzen, um sie zu beantworten.

> <comment-title></comment-title> Dieses Tutorial ist ein wenig untypisch: Wir werden nicht in Galaxy arbeiten, sondern hauptsächlich außerhalb, indem wir Datenbanken und Tools über ihre eigenen Webschnittstellen navigieren. Das Ziel dieses Tutorials ist es, verschiedene Quellen für biologische Daten in unterschiedlichen Dateiformaten und mit unterschiedlichen Informationen zu veranschaulichen.
{: .comment}

> <agenda-title></agenda-title>
>
> In diesem Tutorium werden wir uns mit:
>
> 1. TOC
> {:toc}
{: .agenda}

# Suche nach humanen Opsinen

Um nach menschlichen Opsinen zu suchen, beginnen wir mit dem [NCBI Genome Data Viewer](https://www.ncbi.nlm.nih.gov/genome/gdv). Der NCBI Genome Data Viewer (GDV) ({% cite rangwala2021accessing %}) ist ein Genom-Browser, der die Erforschung und Analyse von annotierten eukaryotischen Genomassemblies unterstützt. Der GDV-Browser zeigt biologische Informationen an, die einem Genom zugeordnet sind, einschließlich Genannotation, Variationsdaten, BLAST-Alignments und experimentelle Studiendaten aus den Datenbanken NCBI GEO und dbGaP. In den GDV-Versionshinweisen werden neue Funktionen für diesen Browser beschrieben.

> <hands-on-title>Öffnen des NCBI Genome Data Viewer</hands-on-title>
>
> 1. Öffnen Sie den NCBI Genome Data Viewer unter [www.ncbi.nlm.nih.gov/genome/gdv](https://www.ncbi.nlm.nih.gov/genome/gdv/)
>
{: .hands-on}

Die Homepage enthält einen einfachen "Baum des Lebens", in dem der menschliche Knoten hervorgehoben ist, weil er der Standardorganismus für die Suche ist. Wir können das im Feld *Organismen suchen* ändern, aber wir lassen das vorerst, da wir an menschlichen Opsinen interessiert sind.

![Bildschirmfoto der Startseite des Genome Data Viewer, das Wort "opsin" wird in das Suchfeld geschrieben und das Ergebnis wird in der Vorschau angezeigt](./images/GenomeDataViewerpage.png "Genome Data Viewer Startseite")

Das Feld auf der rechten Seite zeigt mehrere Zusammenstellungen des Genoms von Interesse und eine Karte der Chromosomen in diesem Genom. Dort können wir nach Opsinen suchen.

> <hands-on-title>Suche nach OpsinsOpen NCBI Genome Data Viewer</hands-on-title>
>
> 1. Geben Sie `opsin` in das Feld *Suche im Genom* ein
> 2. Klicken Sie auf das Lupensymbol oder drücken Sie <kbd>Eingabe<kbd>
>
{: .hands-on}

Unterhalb des Kastens wird nun eine Tabelle mit Genen, die mit Opsin verwandt sind, zusammen mit ihren Namen und ihrer Lage, d.h. der Chromosomennummer, sowie der Anfangs- und Endposition angezeigt

In der Liste der Gene, die mit dem Suchbegriff Opsin in Verbindung stehen, befinden sich das Rhodopsin-Gen (RHO) und drei Zapfenpigmente, kurz-, mittel- und langwellensensitive Opsine (für die Erkennung von blauem, grünem und rotem Licht). Es gibt noch weitere Einheiten, z. B. eine -LCR (Locus Control Region), mutmaßliche Gene und Rezeptoren.

Mehrere Treffer befinden sich auf dem X-Chromosom, einem der geschlechtsbestimmenden Chromosomen.

> <question-title></question-title>
>
> 1. Wie viele Gene wurden auf dem Chromosom X gefunden?
> 2. Wie viele proteinkodierende Gene gibt es?
>
> > <solution-title></solution-title>
> >
> > 1. Die Treffer in ChrX sind:
> > - OPSIN-LCR
> > - OPN1LW
> > - OP1MW
> > - OPN1MW2
> > - OPN1MW3
> >
> > 2. Wenn wir mit dem Mauszeiger über jedes Gen fahren, öffnet sich ein Kasten, und wir können auf *Details* klicken, um mehr über jedes Gen zu erfahren. Dann erfahren wir, dass das erste (OPSIN-LCR) nicht proteincodierend ist, sondern eine Genregulationsregion darstellt und die anderen proteincodierende Gene sind. Es gibt also 4 proteinkodierende Gene für Opsine auf Chromosom X. Insbesondere enthält Chromosom X ein rotes Pigmentgen (OPN1LW) und drei grüne Pigmentgene (OPN1MW, OPN1MW2 und OPN1MW3 in der Referenzgenomeinheit).
> >
> {: .solution}
{: .question}

Konzentrieren wir uns nun auf ein bestimmtes Opsin, das Gen OPN1LW.

> <hands-on-title>Open Genome Browser für das Gen OPN1LW</hands-on-title>
>
> 1. Klicken Sie auf den blauen Pfeil, der in der Ergebnistabelle erscheint, wenn Sie mit der Maus über die Zeile OPN1LW fahren
>
{: .hands-on}

Sie sollten auf [dieser Seite](https://www.ncbi.nlm.nih.gov/genome/gdv/browser/genome/?id=GCF_000001405.40) gelandet sein, das ist die Genomansicht des Gens OPN1LW.

![Genomdaten-Viewer des Gens OPN1LW, Screenshot der beiden Hauptfenster des Viewers, mit Chromosomen auf der linken Seite und dem Feature-Viewer auf der rechten Seite](./images/GenomeDataViewerofgeneOPN1LW.png "Genomdaten-Viewer des Gens OPN1LW")

Es gibt viele Informationen auf dieser Seite, konzentrieren wir uns momentan auf einen Abschnitt.

1. Der Genomdaten-Viewer, oben, sagt uns, dass wir die Daten des Organismus `Homo sapiens`, der Assemblierung `GRCh38.p14` und insbesondere von `Chr X` (Chromosom X) betrachten. Jede dieser Informationen hat eine eindeutige ID.
2. Das gesamte Chromosom ist direkt darunter dargestellt, und die Positionen entlang der kurzen `p` und langen `q` Arme sind hochgestellt.
3. Ein blauer Kasten unterstreicht, dass wir uns jetzt auf die Region konzentrieren, die dem Gen `OPN1LW` entspricht.

   Es gibt mehrere Möglichkeiten, mit dem Viewer unten zu interagieren. Versuchen Sie zum Beispiel, mit der Maus über die Punkte zu fahren, die Exons in der blauen Box darstellen.

4. In der nachstehenden Grafik ist die Gensequenz als grüne Linie dargestellt, die Exons (proteincodierende Fragmente) sind durch grüne Rechtecke gekennzeichnet.

   Fahren Sie mit der Maus über die grüne Linie, die `NM_020061.6` (unser Gen von Interesse) entspricht, um genauere Informationen zu erhalten.

   > <question-title></question-title>
   >
   > 1. An welcher Stelle befindet sich das OPN1LW-Segment?
   > 2. Wie lang ist das OPN1LW-Segment?
   > 3. Was sind Introns und Exons?
   > 4. Wie viele Exons und Introns befinden sich im OPN1LW-Gen?
   > 5. Wie lang ist der kodierende Bereich insgesamt?
   > 6. Wie ist die Verteilung zwischen kodierenden und nicht kodierenden Regionen? Was bedeutet das in Bezug auf die Biologie?
   > 7. Wie lang ist die Länge des Proteins in Anzahl der Aminosäuren?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. Von 154.144.243 bis 154.159.032
   > > 2. 1.4790 Nukleotide, gefunden bei *Span auf 14790 nt, Nukleotide)*
   > > 3. Eukaryotische Gene sind oft durch nicht-kodierende Regionen unterbrochen, die als Intervallsequenzen oder Introns bezeichnet werden. Die kodierenden Bereiche werden Exons genannt.
   > > 4. Aus diesem Diagramm kann man erkennen, dass das OPN1LW-Gen aus 6 Exons und 5 Introns besteht und dass die Introns viel größer sind als die Exons.
   > > 5. Die CDS-Länge beträgt 1.095 Nukleotide.
   > > 6. Von den 14790 nt des Gens kodieren nur 1095 nt für ein Protein, was bedeutet, dass weniger als 8 % der Basenpaare den Code enthalten. Wenn dieses Gen in Zellen der menschlichen Netzhaut exprimiert wird, wird eine RNA-Kopie des gesamten Gens synthetisiert. Dann werden die Intron-Regionen herausgeschnitten und die Exon-Regionen zusammengefügt, um die reife mRNA zu erzeugen (ein Prozess, der als Spleißen bezeichnet wird), die dann von den Ribosomen bei der Herstellung des roten Opsin-Proteins übersetzt wird. In diesem Fall werden 92 % des ursprünglichen RNA-Transkripts herausgeschnitten, so dass der reine Proteincode übrig bleibt.
   > > 7. Die Länge des resultierenden Proteins beträgt 364 aa, Aminosäuren.
   > {: .solution}
   {: .question}

Aber wie lautet die Sequenz dieses Gens? Es gibt mehrere Möglichkeiten, diese Information zu erhalten, wir werden eine der intuitivsten durchgehen.

> <hands-on-title>Open Genome Browser für das Gen OPN1LW</hands-on-title>
>
> 1. Klicken Sie auf den {% icon tool %} *Werkzeuge* oben rechts in der Box, die das Gen
> 2. Klicken Sie auf *Sequence Text View*
{: .hands-on}

Diese Tafel zeigt die DNA-Sequenz der Introns (in grün), sowie die der Exons (in rosa, einschließlich der übersetzten Proteinsequenz unten).

![Screenshot der Sequenzansicht der NHI-Ressource, Text ist in verschiedenen Farben hervorgehoben](./images/SequenceTextView.png "Sequence Text View")

Diese Sequenzbox zeigt im Moment nicht das gesamte Gen, sondern eine Teilsequenz davon. Mit den Pfeilen *Prev Page* und *Next Page* können Sie sich im genetischen Code stromaufwärts und stromabwärts bewegen oder mit der Schaltfläche *Go To Position* an einer bestimmten Stelle beginnen. Wir schlagen vor, mit dem Anfang des kodierenden Teils des Gens zu beginnen, der sich, wie wir bereits gelernt haben, an Position 154,144,243 befindet.

> <hands-on-title>Zu einer bestimmten Position in der Sequenzansicht gehen</hands-on-title>
>
> 1. Klicken Sie auf *Zur Position gehen*
> 2. Tippen Sie auf `154144243`
>
>    Wir müssen die Kommas entfernen, um den Wert zu validieren
{: .hands-on}

Die hier lila hervorgehobene Sequenz signalisiert eine regulatorische Region.

> <question-title></question-title>
>
> 1. Wie lautet die erste Aminosäure des entstehenden Proteinprodukts?
> 2. Welche ist die letzte?
> 3. Können Sie sich die ersten drei und die letzten drei AAs dieses Proteins merken?
>
> > <solution-title></solution-title>
> >
> > 1. Das entsprechende Protein beginnt mit Methionin, M (das tun sie alle).
> > 2. Das letzte AA des letzten Exons (zu finden auf der 2. Seite) ist Alanin (A). Danach kommt das Stoppcodon TGA, das nicht in eine AA übersetzt wird.
> > 3. Die ersten drei AAs sind: M,A,Q; die letzten drei: S,P,A.
> >
> {: .solution}
{: .question}

Wir können jetzt die *Sequenzansicht* schließen.

Von dieser Ressource können wir auch Dateien in verschiedenen Formaten erhalten, die das Gen beschreiben. Sie sind im Bereich *Download* verfügbar.

1. *Download FASTA* ermöglicht es uns, das einfachste Dateiformat herunterzuladen, um die Nukleotidsequenz des gesamten sichtbaren Bereichs des Genoms darzustellen (länger als nur ein Gen).
2. *Download GenBank flat file* ermöglicht uns den Zugriff auf die auf dieser Seite (und darüber hinaus) verfügbaren Anmerkungen in einem Plaintext-Format.
3. *Download Track Data* ermöglicht es uns, zwei der in den Folien vorgestellten Dateiformate einzusehen: die Formate GFF (GFF3) und BED. Wenn Sie die Tracks ändern, kann es sein, dass die einzelnen Formate nicht mehr verfügbar sind.

# Suche nach weiteren Informationen über unser Gen

Verschaffen wir uns nun einen Überblick über die Informationen, die wir (in der Literatur) über unser Gen haben, indem wir die NCBI-Ressourcen nutzen

> <hands-on-title>Zu einer bestimmten Position in der Sequenzansicht gehen</hands-on-title>
>
> 1. Öffnen Sie die NCBI-Suche unter [www.ncbi.nlm.nih.gov/search](https://www.ncbi.nlm.nih.gov/search/)
> 2. Geben Sie `OPN1LW` in das Suchfeld *Search NCBI* ein
>
{: .hands-on}

![Screenshot der NIH-Ergebnisseite, mit Karten namens Literatur, Gene, Proteine, Genome, Klinisch und PubChem](./images/NIHresults.png "Ergebnisse der Suche nach `OPN1LW` im NCBI").

## Literatur

Beginnen wir mit der Literatur und insbesondere den Ergebnissen von *PubMed* oder *PubMed Central*

> <details-title>Was ist der Unterschied zwischen PubMed und PubMed Central? </details-title>
>
> PubMed ist eine biomedizinische Literaturdatenbank, die die Zusammenfassungen der Veröffentlichungen in der Datenbank enthält.
>
> PubMed Central ist ein Volltext-Repositorium, das den Volltext von Veröffentlichungen in der Datenbank enthält.
>
> Auch wenn die genaue Anzahl der Treffer zeitlich von der obigen Abbildung abweichen kann, sollte jeder Genname in PubMed Central (durchsucht in den Volltexten der Veröffentlichungen) mehr Treffer aufweisen als in PubMed (durchsucht nur in den Abstracts).
>
{: .details}

> <hands-on-title>Öffne PubMed</hands-on-title>
>
> 1. Klicken Sie auf *PubMed* im Feld *Literature*
>
{: .hands-on}

Sie haben in PubMed, einer kostenlosen Datenbank für wissenschaftliche Literatur, die Ergebnisse einer vollständigen Suche nach Artikeln eingegeben, die direkt mit diesem Genlocus in Verbindung stehen.

Wenn Sie auf den Titel eines Artikels klicken, können Sie die Kurzfassung des Artikels sehen. Wenn Sie sich auf einem Universitätscampus befinden, wo es einen Online-Zugang zu bestimmten Zeitschriften gibt, sehen Sie möglicherweise auch Links zu vollständigen Artikeln. PubMed ist Ihr Zugang zu einer großen Vielfalt an wissenschaftlicher Literatur im Bereich der Biowissenschaften. Auf der linken Seite jeder PubMed-Seite finden Sie Links zu einer Beschreibung der Datenbank, zur Hilfe und zu Anleitungen für die Suche.

> <question-title></question-title>
>
> 1. Können Sie erraten, welche Art von Erkrankungen mit diesem Gen verbunden sind?
>
> > <solution-title></solution-title>
> >
> > 1. Wir werden diese Frage später beantworten
> >
> {: .solution}
{: .question}

> <hands-on-title>Zurück zur NCBI-Suchseite</hands-on-title>
>
> 1. Zurück zur [NCBI-Suchseite](https://www.ncbi.nlm.nih.gov/search/all/?term=OPN1LW)
{: .hands-on}

## Klinisch

Konzentrieren wir uns nun auf das Feld *Klinisch*, und zwar speziell auf *OMIM*. OMIM, das Online Mendeliam Inheritance in Man (and woman!), ist ein Katalog menschlicher Gene und genetischer Störungen.

> <hands-on-title>Open OMIM</hands-on-title>
>
> 1. Klicken Sie auf *OMIM* im Feld *Klinisch*
>
{: .hands-on}

Jeder OMIM-Eintrag steht für eine genetische Störung (hier vor allem Arten von Farbenblindheit), die mit Mutationen in diesem Gen verbunden ist.

> <hands-on-title>Lesen Sie so viel, wie es Ihr Interesse erfordert</hands-on-title>
>
> 1. Folgen Sie den Links, um weitere Informationen zu den einzelnen Einträgen zu erhalten
>
{: .hands-on}

> <comment-title>Lesen Sie so viel, wie es Ihr Interesse erfordert</comment-title>
>
> Für weitere Informationen über OMIM selbst klicken Sie bitte auf das OMIM-Logo oben auf der Seite. Über OMIM ist eine Fülle von Informationen zu unzähligen Genen im menschlichen Genom verfügbar, und alle Informationen sind durch Verweise auf die neuesten Forschungsartikel untermauert.
>
{: .comment}

Wie wirken sich Variationen des Gens auf das Proteinprodukt und seine Funktionen aus? Gehen wir zurück zur NIH-Seite und untersuchen wir die Liste der Einzelnukleotid-Polymorphismen (SNPs), die durch genetische Studien in diesem Gen entdeckt wurden.

> <hands-on-title>Öffne dbSNP</hands-on-title>
>
> 1. Zurück zur [NCBI-Suchseite](https://www.ncbi.nlm.nih.gov/search/all/?term=OPN1LW)
> 2. Klicken Sie auf *dbSNP* im Feld *Klinisch*
>
{: .hands-on}

![Screenshot der dbSNPs-Seite über das Gen OPN1LW. Drei Hauptfenster, das linke zum Filtern der Suche anhand von Tags, das mittlere zur Anzeige der Ergebnisse, das rechte für eine detailliertere und programmatische Suche](./images/dbSNPs.png "dbSNP in OPN1LW")

> <question-title></question-title>
>
> 1. Welche klinische Bedeutung haben die rs5986963 und rs5986964 (die ersten beiden Varianten, die zum Zeitpunkt der Erstellung dieses Tutorials aufgelistet waren)?
> 2. Was ist die funktionelle Folge von rs104894912?
> 3. Was ist die funktionelle Folge von rs104894913?
>
> > <solution-title></solution-title>
> >
> > 1. Die klinische Signifikanz ist `benign`, so dass es scheint, dass sie keinen Einfluss auf das endgültige Proteinprodukt haben
> > 2. rs104894912-Mutation führt zu einer `stop_gained`-Variante, die das resultierende Protein zu früh abschneidet und daher `pathogenic` ist
> > 3. rs104894913 Mutation führt zu einem `missense_variant`, auch `pathogenic`.
> >
> {: .solution}
{: .question}

Lassen Sie uns mehr über die Variante rs104894913 herausfinden

> <hands-on-title>Erfahren Sie mehr über eine Variante dbSNP</hands-on-title>
>
> 1. Klicken Sie auf `rs104894913`, um die dazugehörige Seite zu öffnen (https://www.ncbi.nlm.nih.gov/snp/rs104894913)
> 2. Klicken Sie auf *Clinical Significance*
>
>    > <question-title></question-title>
>    >
>    > Welche Art von Erkrankung ist mit der Variante rs104894913 verbunden?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > Der Name der damit verbundenen Krankheit ist "Protan-Defekt". Eine schnelle Internetrecherche mit Ihrer Suchmaschine wird Ihnen zeigen, dass es sich um eine Art von Farbenblindheit handelt.
>    > >
>    > {: .solution}
>    {: .question}
>
> 3. Klicken Sie auf die *Variantendetails*
>
>    > <question-title></question-title>
>    >
>    > 1. Welche Substitution ist mit dieser Variante verbunden?
>    > 2. Welche Auswirkung hat diese Substitution auf das Codon und die Aminosäure?
>    > 3. An welcher Stelle des Proteins befindet sich diese Substitution?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. Die Substitution `NC_000023.10:g.153424319G>A` entspricht dem Wechsel von einem Guanin (G) zu einem Adenin (A)
>    > > 2. Diese Substitution ändert das Codon `GGG`, ein Glycin, in `GAG`, ein Glutathion
>    > > 3. `p.Gly338Glu` bedeutet, dass sich die Substitution an Position 338 des Proteins befindet.
>    > {: .solution}
>    {: .question}
{: .hands-on}

Was bedeutet diese Substitution für das Protein? Schauen wir uns dieses Protein einmal genauer an.

## Protein

> <hands-on-title>Protein öffnen</hands-on-title>
>
> 1. Zurück zur [NCBI-Suchseite](https://www.ncbi.nlm.nih.gov/search/all/?term=OPN1LW)
> 2. Klicken Sie auf *Protein* im Feld *Proteine*
> 3. Klicken Sie auf `OPN1LW – opsin 1, long wave sensitive` in der Box oben
>
{: .hands-on}

![Screenshot der Opsin 1 NIH-Protein-Seite, zwei Hauptfenster. Die linke Seite enthält Informationen über das Gen, die rechte Seite ist ein Inhaltsverzeichnis und eine Reihe von Links zu anderen Ressourcen](./images/Opsin1NIH.png "Opsin 1 NIH protein page")

Auf dieser Seite werden wieder einige Daten präsentiert, die wir kennen (z. B. die Verteilung der Exons entlang der Gensequenz).

> <hands-on-title>Download der Proteinsequenzen</hands-on-title>
>
> 1. Klicken Sie auf *Download Datasets*
> 2. Wählen Sie
>    - `Gene Sequences (FASTA)`
>    - `Transcript sequences (FASTA)`
>    - `Protein sequences (FASTA)`
> 3. Klicken Sie auf die Schaltfläche *Download*
> 4. Öffnen Sie die heruntergeladene ZIP-Datei
{: .hands-on}

> <question-title></question-title>
>
> 1. Was enthält der Ordner?
> 2. Glauben Sie, dass sie gute Datenpraktiken umgesetzt haben?
>
> > <solution-title></solution-title>
> >
> > 1. Der Ordner enthält
> >    - ein Ordner `ncbi_datasets` mit verschiedenen Unterordnern darin, die einige Datendateien (mehrere Formate) führen,
> >    - eine `README.md` (eine Markdown-Datei), die zusammen mit den Daten "reisen" soll und erklärt, wie die Daten gewonnen wurden, wie die Struktur des Unterordners, der die Daten enthält, aussieht und wo man eine ausführliche Dokumentation findet.
> > 2. Es ist auf jeden Fall eine gute Datenverwaltungspraxis, die Benutzer (nicht nur für Ihre Mitarbeiter, sondern auch für Sie selbst in nicht allzu ferner Zukunft, wenn Sie vergessen haben, woher die Datei in Ihrem Download-Ordner stammt) zur Datenquelle und zur Datenstruktur zu führen.
> >
> {: .solution}
{: .question}

# Suche nach Sequenz

Was können wir mit den soeben heruntergeladenen Sequenzen tun? Nehmen wir an, dass wir gerade die Transkripte sequenziert haben, die wir durch ein Experiment isoliert haben - wir kennen also die Sequenz der Entität, die uns interessiert, aber wir wissen nicht, um welche es sich handelt. In diesem Fall müssen wir die gesamte Datenbank mit den der Wissenschaft bekannten Sequenzen durchsuchen und unsere unbekannte Entität mit einem Eintrag abgleichen, der eine Annotation enthält. Tun wir es.

> <hands-on-title>Suche die Proteinsequenz gegen alle Proteinsequenzen</hands-on-title>
>
> 1. Öffnen Sie (mit dem einfachsten Texteditor, den Sie installiert haben) die soeben heruntergeladene Datei `protein.faa`.
> 2. Kopieren Sie den Inhalt
> 3. Öffnen Sie BLAST [blast.ncbi.nlm.nih.gov](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
> 4. Klicken Sie auf `Protein BLAST, protein > protein`
>
>    Wir werden in der Tat eine Proteinsequenz verwenden, um eine Datenbank mit Proteinen zu durchsuchen
>
> 5. Fügen Sie die Proteinsequenz in das große Textfeld ein
> 6. Überprüfen Sie die restlichen Parameter
> 7. Klicken Sie auf die blaue Schaltfläche `BLAST`
{: .hands-on}

Diese Phase wird einige Zeit in Anspruch nehmen, schließlich gibt es irgendwo einen Server, der die Gesamtheit der bekannten Sequenzen mit Ihrem Ziel vergleicht. Wenn die Suche abgeschlossen ist, sollte das Ergebnis etwa so aussehen wie unten:

![Screenshot der BLAST-Ergebnisse, oben eine große Kopfzeile und unten die Ergebnisse in Tabellenform](./images/BLASTresults.png "BLAST results")

> <hands-on-title>Grafische Zusammenfassung der Proteinsequenzen</hands-on-title>
>
> 1. Klicken Sie auf die Registerkarte *Graphic Summary*
{: .hands-on}

Wir haben Zugriff auf einen Kasten mit vielen farbigen Linien. Jede Linie steht für einen Treffer aus Ihrer Blast-Suche. Wenn Sie auf eine rote Linie klicken, gibt der schmale Kasten direkt darüber eine kurze Beschreibung des Treffers.

> <hands-on-title>Beschreibungen der Proteinsequenzen</hands-on-title>
>
> 1. Klicken Sie auf die Registerkarte *Beschreibungen*
{: .hands-on}

> <question-title></question-title>
>
> 1. Was ist der erste Treffer? Haben Sie ihn erwartet?
> 2. Was sind die anderen Treffer und für welche Organismen?
>
> > <solution-title></solution-title>
> >
> > 1. Der erste Treffer ist unser rotes Opsin. Das ist ermutigend, denn die beste Übereinstimmung sollte die Abfragesequenz selbst sein, und Sie haben diese Sequenz aus diesem Geneintrag erhalten.
> > 2. Andere Treffer sind andere Opsine. Dazu gehören Einträge von anderen Primaten (z.B. `Pan troglogytes`).
> >
> {: .solution}
{: .question}

Die Treffer sind für unser rotes Opsin beim Menschen, aber auch für andere Opsine bei anderen Primaten. Das könnte man sich wünschen, wenn man diese Daten zum Beispiel zur Erstellung eines phylogenetischen Baums verwenden will. Wenn wir aber ziemlich sicher sind, dass die Sequenz, die uns interessiert, vom Menschen stammt, können wir die Suche auch nur auf menschliche Sequenzen beschränken.

> <hands-on-title>Filtern einer BLAST-Suche</hands-on-title>
>
> 1. Klicken Sie auf *Suche bearbeiten*
> 2. Geben Sie `Homo sapiens` in das Feld *Organismus* ein
> 3. Klicken Sie auf die blaue Schaltfläche `BLAST`
{: .hands-on}

Mit dieser neuen Suche finden wir die anderen Opsine (grünes, blaues, Stäbchenpigment) in der Liste. Andere Treffer haben eine geringere Anzahl von übereinstimmenden Rückständen. Wenn Sie auf eine der farbigen Linien in der *Grafischen Zusammenfassung* klicken, erhalten Sie weitere Informationen über diesen Treffer und können sehen, wie ähnlich jeder einzelne dem roten Opsin, unserer ursprünglichen Suchsequenz, ist. Je weiter Sie in der Liste nach unten gehen, desto weniger hat jede nachfolgende Sequenz mit dem roten Opsin gemeinsam. Jede Sequenz wird im Vergleich zu rotem Opsin in einem so genannten paarweisen Sequenzabgleich dargestellt. Später werden Sie mehrere Sequenzalignments erstellen, aus denen Sie die Beziehungen zwischen den Genen erkennen können.

> <details-title>Mehr Details zu BLAST-Ergebnissen</details-title>
>
> Die Anzeigen enthalten zwei prominente Maße für die Bedeutung des Treffers:
>
> 1. der BLAST Score - lableled Score (bits)
>
>    Der BLAST-Score gibt die Qualität des besten Alignments zwischen der Abfragesequenz und der gefundenen Sequenz (Treffer) an. Je höher die Punktzahl, desto besser ist das Alignment. Die Punktzahl wird durch Nichtübereinstimmungen und Lücken im besten Alignment verringert. Die Berechnung des Scores ist komplex und beinhaltet eine Substitutionsmatrix, d. h. eine Tabelle, die jedem Paar von ausgerichteten Resten einen Score zuweist. Die am häufigsten verwendete Matrix für das Protein-Alignment ist BLOSUM62.
>
> 2. der Erwartungswert (mit Expect oder E bezeichnet)
>
>     Der Erwartungswert E eines Treffers sagt aus, ob der Treffer wahrscheinlich aus einer zufälligen Ähnlichkeit zwischen Treffer und Abfrage oder aus einer gemeinsamen Abstammung von Treffer und Abfrage resultiert. ()
>
>     > <comment-title>Filtern einer BLAST-Suche</comment-title>
>     >
>     > Wenn E kleiner als $$10\mathrm{e}{-100}$$ ist, wird es manchmal als 0,0 angegeben.
>     {: .comment}
>
>     Der Erwartungswert ist die Anzahl der Treffer, die man rein zufällig erwarten würde, wenn man in einem zufälligen Genom von der Größe des menschlichen Genoms nach seiner Sequenz suchen würde.
>
>      $$E = 25$$ bedeutet, dass man in einem Genom dieser Größe rein zufällig 25 Treffer erwarten könnte. Ein Treffer mit $$E = 25$$ ist also wahrscheinlich ein Zufallstreffer und bedeutet nicht, dass die Treffersequenz eine gemeinsame Abstammung mit Ihrer Suchsequenz hat.
>
>      Erwartungswerte von etwa 0,1 können biologisch signifikant sein oder auch nicht (um dies zu entscheiden, sind weitere Tests erforderlich).
>
>      Aber sehr kleine Werte von E bedeuten, dass der Treffer biologisch signifikant ist. Die Übereinstimmung zwischen Ihrer Suchsequenz und diesem Treffer muss auf eine gemeinsame Abstammung der Sequenzen zurückzuführen sein, denn die Wahrscheinlichkeit ist einfach zu gering, dass die Übereinstimmung zufällig zustande kommt. Zum Beispiel bedeutet $$E = 10\mathrm{e}{-18}$$ für einen Treffer im menschlichen Genom, dass man nur eine zufällige Übereinstimmung in einer Milliarde Milliarden verschiedener Genome von der gleichen Größe wie das menschliche Genom erwarten würde.
>
>      Der Grund, warum wir glauben, dass wir alle von gemeinsamen Vorfahren abstammen, ist, dass eine massive Sequenzähnlichkeit in allen Organismen einfach zu unwahrscheinlich ist, um ein zufälliges Ereignis zu sein. Jede Familie ähnlicher Sequenzen in vielen Organismen muss sich aus einer gemeinsamen Sequenz in einem entfernten Vorfahren entwickelt haben.
>
{: .details}

> <hands-on-title>Entladen </hands-on-title>
>
> 1. Klicken Sie auf die Registerkarte *Beschreibungen*
> 2. Klicken Sie auf einen beliebigen Sequenztreffer
> 3. Klicken Sie auf *Download*
> 4. Wählen Sie `FASTA (aligned sequences)`
{: .hands-on}

Es wird ein neuer, etwas anderer Dateityp heruntergeladen: eine ausgerichtete FASTA. Wenn Sie möchten, können Sie diese vor dem nächsten Abschnitt untersuchen.

Während wir in den vorangegangenen Abschnitten dieses Tutorials ausgiebig die Webschnittstellen der Tools genutzt haben (Genomic Viewer, schnelles Scannen der Literatur, Lesen von Annotationen usw.), ist diese BLAST-Suche ein Beispiel für einen Schritt, den Sie mit Galaxy vollständig automatisieren könnten.

> <hands-on-title> Ähnlichkeitssuche mit BLAST in Galaxy </hands-on-title>
>
> 1. Erstellen Sie einen neuen Verlauf für diese Analyse
>
>    {% snippet faqs/galaxy-de/histories_create_new.md %}
>
> 2. Umbenennen der Historie
>
>    {% snippet faqs/galaxy-de/histories_rename.md %}
>
> 3. Importieren Sie die Proteinsequenz über einen Link aus [Zenodo]({{ page.zenodo_link }}) oder aus gemeinsam genutzten Galaxy-Datenbibliotheken:
>
>    ```text
>    {{ page.zenodo_link }}/files/protein.faa
>    ```
>
>    {% snippet faqs/galaxy-de/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy-de/datasets_import_from_data_library.md %}
>
> 1. {% tool [NCBI BLAST+ blastp](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/2.10.1+galaxy2) %} mit den folgenden Parametern:
>    - *"Protein query sequence(s) "*: `protein.faa`
>    - *"Fachdatenbank/Sequenzen "*: `Locally installed BLAST database`
>    - *"Protein BLAST database "*: `SwissProt`
>
>      Um nur nach annotierten Sequenzen in UniProt zu suchen, müssen wir die neueste Version von `SwissProt` auswählen
>
>    - *"Erwartungswertabschaltung einstellen "*: `0.001`
>    - *"Ausgabeformat "*: `Tabular (extended 25 columns)`
>
{: .hands_on}

> <question-title></question-title>
>
> Glauben Sie, dass wir genau die gleichen Ergebnisse wie bei unserer ursprünglichen Suche nach `opsin` in [www.ncbi.nlm.nih.gov/genome/gdv](https://www.ncbi.nlm.nih.gov/genome/gdv/) sehen? Und warum?
>
> > <solution-title></solution-title>
> >
> > Die Ergebnisse mögen ähnlich sein, aber es gibt durchaus einige Unterschiede. Eine Textsuche unterscheidet sich nicht nur methodisch von einer Sequenzsuche, sondern wir sind in dieser zweiten Runde auch von der Sequenz eines bestimmten Opsins ausgegangen, also von einem Zweig des gesamten Proteinstammbaums. Einige der Familienmitglieder sind sich untereinander sehr ähnlich, so dass diese Art der Suche die gesamte Familie aus einer ziemlich verzerrten Perspektive betrachtet.
> >
> {: .solution}
{: .question}

# Mehr Informationen über unser Protein

Bislang haben wir diese Informationen über Opsine erforscht:
- Wie kann man wissen, welche Proteine eines bestimmten Typs in einem Genom existieren,
- Wie kann man wissen, wo sie sich im Genom befinden?
- Wie erhalte ich mehr Informationen über ein Gen von Interesse,
- wie man ihre Sequenzen in verschiedenen Formaten herunterladen kann,
- wie man diese Dateien für eine Ähnlichkeitssuche verwendet.

Sie sind vielleicht neugierig, wie Sie mehr über die Proteine erfahren können, für die sie kodieren. Wir haben bereits einige Informationen gesammelt (z.B. Krankheiten, die damit verbunden sind), aber in den nächsten Schritten werden wir sie mit Daten über die Proteinstruktur, Lokalisierung, Interaktoren, Funktionen usw. kreuzen.

Das Portal, das man besuchen muss, um alle Informationen über ein Protein zu erhalten, ist [UniProt](https://www.uniprot.org/). Wir können es über eine Textsuche oder über den Gen- oder Proteinnamen durchsuchen. Nehmen wir unser übliches Schlüsselwort `OPN1LW`.

> <hands-on-title>Suche auf UniProt</hands-on-title>
>
> 1. Öffnen [UniProt](https://www.uniprot.org/)
> 2. Geben Sie `OPN1LW` in die Suchleiste ein
> 3. Wählen Sie die Kartenansicht
{: .hands-on}

Der erste Treffer sollte `P04000 · OPSR_HUMAN` sein. Bevor Sie die Seite öffnen, sollten Sie zwei Dinge beachten:

1. Der Name des Proteins `OPSR_HUMAN`, sowie die IDs sind anders als der Genname.
2. Dieser Eintrag hat einen goldenen Stern, was bedeutet, dass er manuell kommentiert und kuratiert wurde.

> <hands-on-title>Öffne ein Ergebnis auf UniProt</hands-on-title>
>
> 1. Klicken Sie auf `P04000 · OPSR_HUMAN`
{: .hands-on}

![Screenshot der Kopfzeile der UniProt-Eingangsseite](./images/UniProt.png "UniProt-Seite")

Dies ist eine lange Seite mit vielen Informationen, wir haben ein [ganzes Tutorial]({% link topics/data-science/tutorials/online-resources-protein/tutorial.md %}) entworfen, um sie durchzugehen.

