---
layout: tutorial_hands_on
topic_name: metagenomics
tutorial_name: general-tutorial
---

# イントロダクション
{:.no_toc}

メタゲノミクスでは、環境内の微生物に関する情報を主に2つの技術で抽出することができます:

- アンプリコンシーケンス（または 16S rRNA/rDNA ）では、生物の rRNA/rDNA のみでシーケンスします
- ショットガンシーケンスでは、環境中の微生物の全ゲノムをシーケンスします

このチュートリアルでは、2種類の解析の一般的な原理とそれぞれの手法の違いを説明します。このような解析をさらに深く知りたい場合は、それぞれの解析の詳細なチュートリアルを確認することをお勧めします。

そのため、私たちは同じ環境から得た2つのデータセットを使用します（1つはアンプリコン用でもう1つはショットガン用）: アルゼンチン Anguil のバルクの土壌で、[project on the Argentinean agricultural pampean soils](https://www.ebi.ac.uk/metagenomics/projects/SRP016633) の研究によって得たもの。このプロジェクトでは、3種類の土地を利用しショットガンとアンプリコンシーケンスを使用して2種類の土壌（バルクと根圏）を解析しました。私たちはアルゼンチン Anguil のバルクの土壌に焦点を当てていきます。

> ### アジェンダ
>
> このチュートリアルでは、以下のことを扱います:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# アンプリコンデータ

アンプリコンシーケンスは特定のゲノム領域の遺伝的な変異を解析する高度に標的化されたアプローチです。
メタゲノミクスの分野において、アンプリコンシーケンスはサンプル中の rRNA データを補足して配列決定することを指しています。
それは細菌や古細菌では 16S で真核生物の場合は 18S です。

> ### {% icon tip %} 背景: 16S リボソーム RNA 遺伝子
> ![The 16S ribosomal RNA gene](../../images/16S_gene.png) <br><br>
>
> 16S rRNA 遺伝子は私たちの目的に非常に適したいくつかの特性を持っています
>
> 1. すべての生物に存在する
> 2. 高度に保存されている領域と高度に可変な領域
> 3. 巨大なリファレンスのデータベース
>
> ![Variable regions](../../images/16S_variableregions.jpg "Variable regions of the 16S rRNA")
>
> 高度に保存された領域を用いて異なった種で標的とする(同じ)遺伝子を特定でき、高度に可変な領域を用いて異なる種が区別できます。
>
{: .tip}

アンプリコンデータでは、サンプル中の配列がどの微生物に由来しているのかを抽出することができます。これは、taxonomic assignation と呼ばれています。
私たちは taxons に配列を割り当ててそれからサンプルでタキソノミーを分類または抽出してみます。

この解析では、[mothur tool suite](https://mothur.org) を使用しますが、使うのはこのツールができることのほんの一部だけです。
より詳しく使用方法を学びたい場合は、[mothur tutorial](../mothur-miseq-sop/tutorial.html) を参照してください。

## データをインポートする

私たちのデータセットはアルゼンチンの2か所の土地から採取した土壌サンプルに由来していて、454 GS FLX Titanium を用いて 16S rDNA V4領域をシーケンスしたものです。このチュートリアルでは、オリジナルの fastq データをダウンサンプリングして fasta に変換されたものを使用しています。オリジナルのデータは EBI Metagenomics で下記の実行番号で取得することができます:

- Pampa の土壌: [SRR531818](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS353016/runs/SRR531818/results/versions/2.0)
- Anguil の土壌: [SRR651839](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS386929/runs/SRR651839/results/versions/2.0)


> ### {% icon hands_on %} ハンズオン: データをアップロードする
>
> 1. [Zenodo](https://zenodo.org/record/815875) または（"Analyses of metagenomics data" 内の）データライブラリからファイルをインポートする
>    - `SRR531818_pampa.fasta`
>    - `SRR651839_anguil.fasta`
>
>    ```
>    https://zenodo.org/record/815875/files/SRR531818_pampa.fasta
>    https://zenodo.org/record/815875/files/SRR651839_anguil.fasta
>    ```
>
>    > ### {% icon tip %} Tip: リンクからデータをインポートする
>    >
>    > * リンクをコピーする
>    > * Galaxy Upload Manager を開く
>    > * **Paste/Fetch Data** を選択する
>    > * テキスト画面にリンクをペーストする
>    > * **Start** を押す
>    {: .tip}
>
>    > ### {% icon tip %} Tip: データライブラリからデータをインポートする
>    >
>    > * 「共有データ」（トップパネル）を選び、それから "Data libraries" を選択する
>    > * Click on "Training data" and then "Analyses of metagenomics data"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
>    デフォルトでは、Galaxy はリンクを名前にしているので、名前を変更します。
>
{: .hands_on}

<!--

Anguil Soil: https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS386929
Pampa Soil https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS353016/runs/SRR531818/results/versions/2.0

Project's data: https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS353016/runs/SRR531818/results/versions/2.0
Project's pipeline: https://www.ebi.ac.uk/metagenomics/pipelines/2.0
Project's QC results: https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS386929/runs/SRR651839/results/versions/2.0
-->

### データセットの準備

mothur でマルチサンプル解析を行いますが、これを行うために、すべてのリードを1つのファイルにまとめて、どのリードがどのサンプルに属しているかを示す、*group file* を作成します。

> ### {% icon hands_on %} ハンズオン: マルチサンプル解析の準備
>
> 1. **Merge.files** {% icon tool %} で次のように設定する
>   - "Merge" → `fasta files`
>   - "Inputs" には2つのサンプル fasta ファイルを選択する
>
> 2. **Make.group** {% icon tool %} で次のように設定する
>   - "Method to create group file" → `Manually specify fasta files and group names`
>   - "Additional": この repeat に2つの要素を追加する
>        - 1つ目の要素
>            - "fasta - Fasta to group" → `SRR531818_pampa` ファイル
>            - "group - Group name" → `pampa`
>        - Second element (click on "Insert Additional")
>            - "fasta - Fasta to group" → `SRR651839_anguil` ファイル
>            - "group - Group name" → `anguil`
>
{: .hands_on}

> ### {% icon tip %} Tip
>
> 私たちは少ない数のサンプルしか持っていないので、手動で指定する方法を使用しました。何百ものサンプルがある場合はこれはすぐに煩わしいものになります。解決方法は？コレクションを使用しましょう！Galaxy のコレクションについて詳しくは、[こちら]() のチュートリアルをご覧ください。
{: .tip}

group ファイルを見てください。これは非常にシンプルなファイルで、2列からなり、1列目はリード名で、2列目はグループ（サンプル）名、今回の場合は `pampa` または `anguil` が入ります。


### 計算のためにファイルの最適化する

私たちは多くの同種の生物をシーケンスしているので、配列の多くは互いに重複していることが予想されます。同じことに莫大な時間を費やすのは計算上無駄なので、 `unique.seqs` コマンドを用いて配列を一意なものにします:

> ### {% icon hands_on %} ハンズオン: 重複している配列を除去する
>
> 1. **Unique.seqs** {% icon tool %} で次のように設定する
>   - "fasta" にはマージした fasta ファイルを選択する
>   - "output format" → `Name File`
>
>    > ### {% icon question %} Question
>    >
>    > 一意な配列はいくつありましたか？除去された重複している配列の数はいくつありましたか？
>    >
>    >    > ### {% icon solution %} Solution
>    >    > 一意な配列は19,502 で重複している配列は498 です。
>    >    >
>    >    > これはこのステップの前の fasta ファイルのライン数と、fasta（もしくは names）アウトプットのライン数を比べることで分かります。
>    >    {: .solution }
>    {: .question}
>
{: .hands_on}

この `Unique.seqs` は2つのファイルをアウトプットとし、1つは一意な配列のみを含んだ fasta ファイルで、*names files* です。
この names file は2列で構成されていて、1列目には一意な配列の配列名があり、２列目には1列目の代表的な配列と同じ配列を持つすべての配列名が書かれています。

```
name          representatives
read_name1    read_name2,read_name,read_name5,read_name11
read_name4    read_name6,read_name,read_name10
read_name7    read_name8
...
```

> ### {% icon hands_on %} ハンズオン: Count sequences
>
> 1. **Count.seqs** {% icon tool %} で次のように設定する
>   - "name" には `Unique.seqs` の name ファイルを選択する
>   - "Use a group file" → `yes`
>   - "group" には `Make.group` の group ファイルを選択する
{: .hands_on}

`Count.seqs` のファイルは、それぞれの一意な配列を代表として表される配列の数を、複数のサンプルにわたって追跡します。このファイルを以降の多くのツールに渡して必要に応じて使用または更新します。

## クオリティコントロール

解析の第一歩はデータの品質をチェックし改善することです。


> ### {% icon comment %} Comment
>
> クオリティコントロールのトピックの詳しい情報については、[こちら]({{site.baseurl}}/topics/sequence-analysis/) のトレーニング資料をご覧ください。
{: .comment}


まずは、データをよく見てみましょう:

> ### {% icon hands_on %} ハンズオン: データを要約する
>
> 1. **Summary.seqs** {% icon tool %} で次のように設定する
>   - "fasta" パラメーターには `Unique.seqs` の fasta ファイルを選択する
>   - "count" には `Count.seqs` の count table を選択する
>   - "output logfile?" → `yes`
>
{: .hands_on}

`summary` のアウトプットファイルはリードごとの情報を表示します。`logfile` のアウトプットにも次のようないくつかの summary の統計が含まれています:

```
              Start    End        NBases     Ambigs   Polymer  NumSeqs
Minimum:      1        80         80         0        3        1
2.5%-tile:    1        104        104        0        3        501
25%-tile:     1        242        242        0        4        5001
Median:       1        245        245        0        4        10001
75%-tile:     1        245        245        0        4        15001
97.5%-tile:   1        247        247        0        6        19501
Maximum:      1        275        275        2        31       20000
Mean:         1        237.519    237.519    0.00495  4.24965
# of unique seqs:   19502
total # of seqs:    20000
```

これは合計で 19,502 の一意な配列があり、80 ～ 275 塩基の長さの異なる配列が主な配列として合計 20,000 配列あることを示しています。また、少なくともいくつかの配列では曖昧な塩基が含まれていることにも注意してください。
さらに、少なくとも1つのリードは 31 塩基の一連のホモポリマーを有していますが、これはおそらく誤りであるのでこのようなリードもフィルタリングしましょう。

もし 20,000 のラウンド数に対して不思議に思っているなら、あなたは正しいです。このチュートリアルではオリジナルのデータセットをサンプルあたり 10,000 のリードにダウンサンプリングして解析のステップにかかる時間を短縮しています。

`Screen.seqs` ツールを使用してデータセット内のリード長、塩基のクオリティ、そしてホモポリマーの最大長をフィルタリングすることができます

次のツールでは曖昧な塩基（`maxambig` パラメーター）、9 塩基以上の連続したホモポリマー（`maxhomop` パラメーター）、そして 275 bp より長い配列や 225 bp より短い配列を除去します。

> ### {% icon hands_on %} ハンズオン: クオリティやリード長で塩基配列をフィルタリングする
>
> 1. **Screen.seqs** {% icon tool %} で次のように設定する
>   - "fasta" には `Unique.seqs` の fasta ファイルを選択する
>   - "minlength" パラメーター → `225`
>   - "maxlength" パラメーター → `275`
>   - "maxambig" パラメーター → `0`
>   - "maxhomop" パラメーター → `8`
>   - "count" には `Count.seqs` の count ファイルを選択する
>
{: .hands_on}

> ### {% icon question %} Question
>
> このスクリーニングステップでいくつのリードが除去されましたか？（ヒント: `Summary.seqs` ツールを再実行してみましょう）
>
> > ### {% icon solution %} Solution
> > 1,804 リードです。
> >
> > これは screen.seqs のアウトプットである bad.accnos のライン数を見るか、スクリーニングステップ前後のサマリーログ間の seqs の合計数を比べることで分かります。
> {: .solution }
{: .question}

## Sequence Alignment

私たちの配列をレファレンスとアライメントすることは OTU の割り当てをよりよく行うのに役立つため [[Schloss et. al.](https://www.ncbi.nlm.nih.gov/pubmed/23018771)] 、16S rRNA の V4 可変領域に配列をアライメントします。このアライメントは Silva のリファレンスデータベースの [[mothur's MiSeq SOP](https://mothur.org/wiki/MiSeq_SOP)] に記載されている通りに作成されています。

> ### {% icon hands_on %} ハンズオン: Align sequences
>
> 1. `silva.v4.fasta` ファイルをヒストリーにインポートする
>
>    ```
>    https://zenodo.org/record/815875/files/silva.v4.fasta
>    ```
>
> 2. **Align.seqs** {% icon tool %} で次のように設定する
>   - "fasta" には `Screen.seqs` のアウトプットである `good.fasta` を選択する
>   - "Select Reference Template from" → `Your history`
>   - "reference" にはリファレンスのファイルである `silva.v4.fasta` を選択する
>   - "flip" → `Yes`
>
>    この手順には数分かかることがあるので、しばらくお待ちください。
>
> 3. **Summary.seqs** {% icon tool %} で次のように設定する
>   - "fasta" パラメーターには `Align.seqs` の aligned アウトプットを選択する
>   - "count" パラメーターには `Screen.seqs` のアウトプットである count_table を選択する
>
{: .hands_on}

アライメントのクオリティを知るために、summary ステップで log のアウトプットを見ることができます:

```
        Start   End NBases  Ambigs  Polymer NumSeqs
Minimum:    2391    10674   9   0   2   1
2.5%-tile:  3080    12071   234 0   4   455
25%-tile:   3080    13424   244 0   4   4545
Median:     3080    13424   245 0   4   9090
75%-tile:   3080    13424   245 0   4   13634
97.5%-tile: 3082    13424   246 0   6   17724
Maximum:    13396   13425   267 0   7   18178
Mean:   3080.6  13380   244.212 0   4.27946
# of unique seqs:   17698
total # of seqs:    18178
```

> ### {% icon question %} Questions
>
> 1. いくつの配列がアライメントされましたか？
> 2. どの位置の間でほとんどのリードがリファレンスとアライメントされていますか？
>
>    > ### {% icon solution %} Solution
>    > 1. 17,698 の配列がアライメントされました
>    > 2. log のアウトプットからほとんどのリードがリファレンスの `3080 から 13424` の位置でうまくアライメントされていることが分かります。これは正確に 16S 遺伝子の V4 標的領域に対応しています。
>    {: .solution }
{: .question}

すべてが同じ領域に重なることを確認するために `Screen.seqs` を再実行して 3,080 より前の位置で開始し、 13,424 以降の位置で終了する配列を取得します。

> ### {% icon hands_on %} ハンズオン: 不適切にアライメントした配列を除去する
>
> 1. **Screen.seqs** {% icon tool %} で次のように設定する
>   - "fasta" → aligned fasta file
>   - "start" → `3080`
>   - "end" → `13424`
>   - "count" には前回の `Screen.seqs` を実行して作成された group ファイルを選択する
>
{: .hands_on}

> ### {% icon question %} Question
> このステップではいくつの配列が除去されましたか？
> > ### {% icon solution %} Solution
> > 4,579 配列が除去されました。これは bad.accnos のアウトプットの配列数です。
> {: .solution }
{: .question}

今私たちは同じアライメント座標で重なっている配列を知っていますが、それらがその領域 **だけ** で重なっていることを確認したいです。そこで配列をフィルタリングして両端のはみ出しているものを除去します。加えて、アライメントの中には、内部ギャップ文字（即ち「-」）を含む列は考慮しませんが、外部ギャップ文字（即ち「.」）だけを含む多くの列があります。これらは情報を失うことなく抜き出すことができます。これはすべて `Filter.seqs` を用いて行います:

> ### {% icon hands_on %} ハンズオン: 配列をフィルタリングする
>
> 1. **Filter.seqs** {% icon tool %} で次のように設定する
>   - "fasta"" には `Sreen.seqs` のアウトプットである `good.fasta` を選択する
>   - "trump" → `.`
{: .hands_on}


## 分類学的な情報を得る

アンプリコンデータ解析における主な質問は以下のとおりです: 環境サンプルにはどのような微生物が存在しているのですか？そしてどれくらいの割合でいるのですか？微生物のコミュニティの構造はどのようになっていますか？

配列を取得してそれらを分類学的に割り当てるというのがこの質問に対する考え方です。そのために、類似性に基づいて配列をグループ化（またはクラスタリング）して Operational Taxonomic Units (OTUs)（単一の「属」または「種」（クラスタリングの閾値に依存する）として扱うことのできる類似の配列のグループ）を定義します。

> ### {% icon tip %} 背景: Operational Taxonomic Units (OTUs)
>
> 16S メタゲノミクスアプローチでは、OTU は 16S rDNA マーカー遺伝子と同種な配列の変異体のクラスターである。これらのクラスターのそれぞれは配列の類似性の閾値に応じて細菌種または属の分類学的な単位を示すことを意図している。典型的には、OTU クラスターでは 16S 遺伝子配列変異体が 97% まで同一だと属レベルで定義されます。98% または 99% 同一だと種の分類まで示唆されます。
>
> ![OTU and cluster with 97% identity threshold](../../images/otu.png "OTU and cluster with 97% identity threshold")
>
> ![OTU graph](../../images/OTU_graph.png "Cladogram of operational taxonomic units (OTUs). Credit: Danzeisen et al. 2013, 10.7717/peerj.237")
>
{: .tip}

まず `Pre.cluster` コマンドを使用して配列をプレクラスタリングし、シーケンス間の違いを 2 塩基まで許容することで、配列の潜在的なエラーを取り除こうと思います。このコマンドはグループごとに配列を分割し、最も多いものから最も少ないものへというように配列の豊富さで並び替えて、他と 2 塩基以上異なるヌクレオチドを特定します。この場合、それらはマージされます。一般的に、100塩基対あたり1塩基の違いを許容することをお勧めします:

> ### {% icon hands_on %} ハンズオン: 配列のプレクラスタリングを行い、望ましくない配列を除去する
>
> 1. **Pre.cluster** {% icon tool %} で次のように設定する
>   - "fasta" には最後に実行した `Filter.seqs` の fasta アウトプットを選択する
>   - "name file or count table" には直近の `Screen.seqs` ステップの count table を選択する
>   - "diffs" → 2
>
>   > ### {% icon question %} Question
>   >
>   >  類似度の高い配列をクラスタリングした後では一意な配列はいくつ残りましたか？
>   > > ### {% icon solution %} Solution
>   > > 10,398 配列
>   > >
>   > > これは fasta のアウトプットのライン数です
>   > {: .solution }
>   {: .question}
>
{: .hands_on}

<!-- optional additional QC: chimera.uchime -->
training set を使用して配列を分類してみようと思います。これは [[mothur's MiSeq SOP](https://mothur.org/wiki/MiSeq_SOP)] でも提供されています。

> ### {% icon hands_on %} ハンズオン: 配列を系統ごとに分類する
>
> 1. Import the `trainset16_022016.pds.fasta` and `trainset16_022016.pds.tax` in your history
>
>    ```
>    https://zenodo.org/record/815875/files/trainset16_022016.pds.fasta
>    https://zenodo.org/record/815875/files/trainset16_022016.pds.tax
>    ```
>
> 2. **Classify.seqs** {% icon tool %} で次のように設定する
>   - "fasta" to the fasta output from `Pre.cluster`
>   - "Select Reference Template from" to `History`
>   - "reference" to `trainset16_022016.pds.fasta` from your history
>   - "Select Taxonomy from" to `History`
>   - "taxonomy" to `trainset16_022016.pds.tax` from your history
>   - "count" to the count table from `Pre.cluster`
>
> This step may take a couple of minutes, now may be a good time to grab a cup of tea :coffee:
>
{: .hands_on}

Have a look at the taxonomy output.

```
name    taxonomy
SRR651839.9109-HXY9DLL01BSHQO-2 Bacteria(99);Verrucomicrobia(99);Spartobacteria(99);unclassified;unclassified;unclassified;
SRR651839.11437-HXY9DLL01AMIPI-2    Bacteria(100);Verrucomicrobia(97);Spartobacteria(97);unclassified;unclassified;unclassified;
SRR651839.15884-HXY9DLL01BDG2F-2    Bacteria(99);unclassified;unclassified;unclassified;unclassified;unclassified;
SRR651839.12048-HXY9DLL01DSOWS-2    Bacteria(100);Verrucomicrobia(100);Spartobacteria(100);unclassified;unclassified;unclassified;
SRR651839.9410-HXY9DLL01BGRG2-2 Bacteria(100);Proteobacteria(100);Alphaproteobacteria(100);Rhizobiales(100);Bradyrhizobiaceae(100);unclassified;
SRR651839.9029-HXY9DLL01E4W9Z-2 Bacteria(100);Verrucomicrobia(100);Spartobacteria(100);unclassified;unclassified;unclassified;
SRR651839.6283-HXY9DLL01CHX25-2 Bacteria(100);Acidobacteria(100);Acidobacteria_Gp4(100);Gp4(100);unclassified;unclassified;
SRR651839.3134-HXY9DLL01DOM67-2 Bacteria(100);Acidobacteria(100);Acidobacteria_Gp6(100);Gp6(100);unclassified;unclassified;
SRR651839.13044-HXY9DLL01ETN7L-2    Bacteria(100);Proteobacteria(100);Alphaproteobacteria(100);Rhizobiales(100);Bradyrhizobiaceae(100);unclassified;
SRR531818.61708-G88ZSJI01AVPPR-2    Bacteria(100);Acidobacteria(99);Acidobacteria_Gp6(99);Gp6(99);unclassified;unclassified;
```

You will see that every read now has a classification.

The next step is then to use this information to know the abundance of the different found taxons. This consists of three steps:
1. first all individual sequences are classified, and get assigned a confidence score (0-100%)
2. next, sequences are grouped at 97% identity threshold (not using taxonomy info)
3. finally, for each cluster, a consensus classification is determined based on the classification of the individual sequences and taking their confidence scores into account

> ### {% icon hands_on %} ハンズオン: Assign sequences to OTUs
>
> 1. **Cluster.split** {% icon tool %} で次のように設定する
>   - "Split by" to `Classification using fasta`
>   - "fasta" to the fasta output from `Pre.cluster`
>   - "taxonomy" to the taxonomy output from `Classify.seqs`
>   - "count" to the count table output from `Pre.cluster`
>   - "Clustering method" to `Average Neighbour`
>   - "cutoff" to `0.15`
>
{: .hands_on}

We obtain a table with the columns being the different identified OTUs, the rows the different distances and the cells the ids of the sequences identified for these OTUs for the different distances.

Next we want to know how many sequences are in each OTU from each group with a distance of 0.03 (97% similarity). We can do this using the `Make.shared` command with the 0.03 cutoff level:

> ### {% icon hands_on %} ハンズオン: Estimate OTU abundance
>
> 2. **Make.shared** {% icon tool %} で次のように設定する
>   - "Select input type" to `OTU list`
>   - "list" to list output from `Cluster.split`
>   - "name file or count table" to the count table from `Pre.cluster`
>   - "label" to `0.03`
>
{: .hands_on}

We probably also want to know the taxonomy for each of our OTUs. We can get the consensus taxonomy for each OTU using the `Classify.otu` command:

> ### {% icon hands_on %} ハンズオン: Classify the OTUs
>
> 3. **Classify.otu** {% icon tool %} で次のように設定する
>   - "list" to output from `Cluster.split`
>   - "count" to the count table from `Pre.cluster`
>   - "Select Taxonomy from" to `History`
>   - "taxonomy" to the taxonomy output from `Classify.seqs`
>   - "label" to `0.03`
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How many OTUs with taxonomic assignation are found for the Anguil sample? And for the Pampa sample?
> 2. What is the annotation of first OTU and its size?
>
>    > ### {% icon solution %} Solution
>    > 1. 2,195 for Anguil and 2,472 for Pampa ("tax.summary")
>    > 2. Otu00001 is associated to 929 sequences and to Bacteria (kingdom), Verrucomicrobia (phylum), Spartobacteria (class) in "taxonomy" file
>    {: .solution }
{: .question}

## Visualization

We have now determined our OTUs and classified them, but looking at a long text file is not very informative.
Let's visualize our data using Krona:

> ### {% icon hands_on %} ハンズオン: Krona
>
>  First we convert our mothur taxonomy file to a format compatible with Krona
>
> - **Taxonomy-to-Krona** {% icon tool %} で次のように設定する
>   - "Taxonomy file" to the taxonomy output from Classify.otu (note: this is a collection input)
>
> - **Krona pie chart** {% icon tool %} で次のように設定する
>   - "Type of input" to `Tabular`
>   - "Input file" to taxonomy output from Taxonomy-to-Krona (collection)
{: .hands_on}

The result is an HTML file with an interactive visualization, for instance try clicking
on one of the rings in the image or playing around with some of the settings.

![Krona output](../../images/krona.png)

This produced a single plot for both your samples, but what if you want to compare
the two samples?

> ### {% icon hands_on %} ハンズオン: Per-sample Krona plots
>
> 1. **Classify.otu** {% icon tool %}
>
>    Hit the rerun button on the `Classify.otu` job in your history and see if you can find settings that will give you per-sample taxonomy data
>
> 2. **Krona** {% icon tool %}
>
>    Now use this new output collection to create per-sample Krona plots
>
{: .hands_on}

In this new Krona output you can switch between the combined plot and the per-sample plots via the selector in the top-left corner.

> ### {% icon question %} Question
> Which soil sample had a higher percentage of Acidobacteria, anguil or pampa? what were the respective percentages?
> > ### {% icon solution %} Solution
> > The anguil sample had a higher proportion of Acidobacteria. The exact percentages can be found by looking at the pie charts at the
> > top right-hand corner after clicking on the label Acidobacteria. For anguil the percentage is 36%, for the pampa sample it is 26%.
> >
> ![krona plot with acidobactaria highlighted](../../images/krona-multisample.png)
>
> {: .solution }
{: .question}


To further explore the community structure, we can visualize it with dedicated tools such as Phinch.

> ### {% icon hands_on %} ハンズオン: Visualization of the community structure with Phinch
>
> 1. **Make.biom** {% icon tool %} で次のように設定する
>   - "shared" to `Make.shared` output
>   - "constaxonomy" to taxonomy output from the first run of `Classify.otu` (collection)
> 1. Expand the dataset and click on the "view biom at phinch" link
>
>     > ### {% icon comment %} Comment
>     >
>     > If this link is not present on your Galaxy, you can download the generated BIOM file and upload directly to Phinch server at [http://phinch.org](http://phinch.org).
>    {: .comment}
>
> 2. Play with the several interactive visualisations:
>
> ![Phinch website interface](../../../../shared/images/phinch_overviewpage.png "Phinch visualizations")
>
{: .hands_on}

Once we have information about the community structure (OTUs with taxonomic structure), we can do more analysis on it: estimation of the diversity of micro-organism, comparison fo diversity between samples, analysis of populations, ... We will not go into detail of such analyses here but you follow our tutorials on amplicon data analyses to learn about them.

# Shotgun metagenomics data

In the previous section, we see how to analyze amplicon data to extract the community structure. Such information can also be extracted from shotgun metagenomic data.

In shotgun data, full genomes of the micro-organisms in the environment are sequenced (not only the 16S or 18S). We can then have access to the rRNA (only a small part of the genomes), but also to the genes of the micro-organisms. Using this information, we can try to answer to questions "What are the micro-organisms doing?" in addition to the question "What micro-organisms are present?".

In this second part, we will use a metagenomic sample of the Pampas Soil ([SRR606451](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS372043/runs/SRR606451/results/versions/2.0)).

## Data upload

> ### {% icon hands_on %} ハンズオン: Data upload
>
> 1. Create a new history
> 2. Import the `SRR606451_pampa` Fasta file from [Zenodo](http://zenodo.org/record/815875) or from the data library (in "Analyses of metagenomics data")
>
>    ```
>    https://zenodo.org/record/815875/files/SRR606451_pampa.fasta
>    ```
>
{: .hands_on}

## Extraction of taxonomic information

As for amplicon data, we can extract taxonomic and community structure information from shotgun data. Different approaches can be used:

- Same approaches as for amplicon data with identification and classification of OTUs

    Such approaches imply a first step of sequence sorting to extract only the 16S and 18S sequences on which the same tools as for amplicon data. However, rRNA sequences represent a low proportion (< 1%) of the shotgun sequences so such an approach is not the most statistically supported

- Assignation of taxonomy on the whole sequences using databases with marker genes

In this tutorial, we use the second approach with MetaPhlAn2. This tools is using a database of ~1M unique clade-specific marker genes (not only the rRNA genes) identified from ~17,000 reference (bacterial, archeal, viral and eukaryotic) genomes.

> ### {% icon hands_on %} ハンズオン: Taxonomic assignation with MetaPhlAn2
>
> 1. **MetaPhlAN2** {% icon tool %} with
>    - "Input file" to the imported file
>    - "Database with clade-specific marker genes" to `locally cached`
>    - "Cached database with clade-specific marker genes" to `MetaPhlAn2 clade-specific marker genes`
>
> This step may take a couple of minutes :coffee:
{: .hands_on}

3 files are generated:

- A tabular file with the community structure

    ```
    #SampleID   Metaphlan2_Analysis
    k__Bacteria 100.0
    k__Bacteria|p__Proteobacteria   86.20712
    k__Bacteria|p__Actinobacteria   13.79288
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria    86.20712
    k__Bacteria|p__Actinobacteria|c__Actinobacteria 13.79288
    ```

    Each line contains a taxa and its relative abundance found for our sample. The file starts with high level taxa (kingdom: `k__`) and go to more precise taxa.


- A BIOM file with the same information as the previous file but in BIOM format

    It can be used then by mothur and other tools requiring community structure information in BIOM format

- A SAM file with the results of the mapping of the sequences on the reference database

> ### {% icon question %} Questions
>
> 1. What is the most precise level we have access to with MetaPhlAn2?
> 2. What are the two orders found in our sample?
> 3. What is the most abundant family in our sample?
>
>    > ### {% icon solution %} Solution
>    > 1. We have access to species level
>    > 2. Pseudomonadales and Solirubrobacterales are found in our sample
>    > 3. The most abundant family is Pseudomonadaceae with 86.21 % of the assigned sequences
>    {: .solution }
{: .question}

Even if the output of MetaPhlAn2 is bit easier to parse than the BIOM file, we want to visualize and explore the community structure with KRONA

> ### {% icon hands_on %} ハンズオン: Interactive visualization with KRONA
>
> 1. **Format MetaPhlAn2 output for Krona** {% icon tool %} with
>    - "Input file" to `Community profile` output of `MetaPhlAn2`
>
> 2. **KRONA pie chart** {% icon tool %} with
>    - "What is the type of your input data" as `MetaPhlan`
>    - "Input file" to the output of `Format MetaPhlAn2`
>
{: .hands_on}

## Extraction of functional information

We would like now to answer the question "What are the micro-organisms doing?" or "Which functions are done by the micro-organisms in the environment?".

In the shotgun data, we have access to the gene sequences from the full genome. We use that to identify the genes, associate them to a function, build pathways, etc to investigate the functional part of the community.

> ### {% icon hands_on %} ハンズオン: Metabolism function identification
>
> 1. **HUMAnN2** {% icon tool %} with
>    - "Input sequence file" to the imported sequence file
>    - "Use of a custom taxonomic profile" to `Yes`
>    - "Taxonomic profile file" to `Community profile` output of `MetaPhlAn2`
>    - "Nucleotide database" to `Locally cached`
>    - "Nucleotide database" to `Full`
>    - "Protein database" to `Locally cached`
>    - "Protein database" to `Full UniRef50`
>    - "Search for uniref50 or uniref90 gene families?" to `uniref50`
>    - "Database to use for pathway computations" to `MetaCyc`
>    - "Advanced Options"
>    - "Remove stratification from output" to `Yes`
>
>    This step is long so we generated the output for you!
>
> 2. Import the 3 files whose the name is starting with "humann2"
>
>    ```
>    https://zenodo.org/record/815875/files/humann2_gene_families_abundance.tsv
>    https://zenodo.org/record/815875/files/humann2_pathways_abundance.tsv
>    https://zenodo.org/record/815875/files/humann2_pathways_coverage.tsv
>    ```
{: .hands_on}

HUMAnN2 generates 3 files

- A file with the abundance of gene families

    Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length. It reflects the relative gene (or transcript) copy number in the community.

    "UNMAPPED" value is the total number of reads which remain unmapped after both alignment steps (nucleotide and translated search). Since other gene features in the table are quantified in RPK units, "UNMAPPED" can be interpreted as a single unknown gene of length 1 kilobase recruiting all reads that failed to map to known sequences.

- A file with the coverage of pathways

    Pathway coverage provides an alternative description of the presence (1) and absence (0) of pathways in a community, independent of their quantitative abundance.

- A file with the abundance of pathways

> ### {% icon question %} Questions
>
> How many gene families and pathways have been identified?
>
>    > ### {% icon solution %} Solution
>    > 44 gene families but no pathways are identified
>    {: .solution }
{: .question}

The RPK for the gene families are quite difficult to interpret in term of relative abundance. We decide then to normalize the values

> ### {% icon hands_on %} ハンズオン: Normalize the gene family abundances
>
> 1. **Renormalize a HUMAnN2 generated table** {% icon tool %} with
>    - "Gene/pathway table" to the gene family table generated with `HUMAnN2`
>    - "Normalization scheme" to `Relative abundance`
>    - "Normalization level" to `Normalization of all levels by community total`
>
>  > ### {% icon question %} Questions
>  >
>  > 1. Which percentage of sequences has not be assigned to a gene family?
>  > 2. What is the most abundant gene family?
>  >
>  >    > ### {% icon solution %} Solution
>  >    > 1. 55% of the sequences has not be assigned to a gene family
>  >    > 2. The most abundant gene family with 25% of sequences is a putative secreted protein
>  >    {: .solution }
>  {: .question}
{: .hands_on}

With the previous analyses, we investigate "Which micro-organims are present in my sample?" and "What function are done by the micro-organisms in my sample?". We can go further in these analyses (for example with combination of functional and taxonomic results). We did not detail that in this tutorial but you can found more analyses in our tutorials on shotgun metagenomic data analyses.

# Conclusion
{:.no_toc}

We can summarize the analyses with amplicon and shotgun metagenomic data:

![Scheme to sum up the analysis](../../images/general-tutorial-scheme.png)

Both analyses are quite complex! However, in this tutorial, we only showed simple cases of metagenomics data analysis with subset of real data.

Check our other tutorials to learn more in details how to analyze metagenomics data.
