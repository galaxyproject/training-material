#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'csv'

renamed_tutorials = {
  '16S Microbial Analysis with Mothur' => '16S Microbial Analysis with mothur (extended)',
  '16S Microbial Analysis with mothur' => '16S Microbial Analysis with mothur (extended)',
  'Quality Control ' => 'Quality Control',
  'RNA-Seq reads to counts' => '1: RNA-Seq reads to counts',
  'RNA-seq counts to genes' => '2: RNA-seq counts to genes',
  'RNA-seq genes to pathways' => '3: RNA-seq genes to pathways',
  'Introduction to Genome Assembly' => 'An Introduction to Genome Assembly',
  'EWAS data analysis of 450k data' => 'Infinium Human Methylation BeadChip',
  'Creating a new tutorial - Defining the technical infrastructure' => 'Tools, Data, and Workflows for tutorials',
  'Creating a new tutorial - Writing content in Markdown' => 'Creating content in Markdown',
  'Running the Galaxy Training material website locally' => 'Running the GTN website locally',
  'Visualizations: JavaScript plugins' => 'JavaScript plugins',
  'Compute and analyze Essential Biodiversity Variables with PAMPA toolsuite' =>
    'Compute and analyze biodiversity metrics with PAMPA toolsuite',
  'Ephemeris for Galaxy Tool Management' => 'Galaxy Tool Management with Ephemeris',
  'Collections: Rule Based Uploader' => 'Rule Based Uploader',
  'Collections: Using dataset collection' => 'Using dataset collections',
  'Data: Downloading and Deleting Data in Galaxy' => 'Downloading and Deleting Data in Galaxy',
  'Histories: Understanding Galaxy history system' => 'Understanding Galaxy history system',
  'Jupyter: Use Jupyter notebooks in Galaxy' => 'Use Jupyter notebooks in Galaxy',
  'Using dataset collection' => 'Using dataset collections',
  'Workflows: Extracting Workflows from Histories' => 'Extracting Workflows from Histories',
  'Workflows: Using Workflow Parameters' => 'Using Workflow Parameters',
  'Exome sequencing data analysis' => 'Exome sequencing data analysis for diagnosing a genetic disease',
  'Galaxy Tool Management' => 'Galaxy Tool Management with Ephemeris',
  'Virtual screening of the SARS-CoV-2 main protease with rDock and pose scoring' =>
    'Virtual screening of the SARS-CoV-2 main protease with rxDock and pose scoring',
  'Refining Genome Annotations with Apollo' => 'Refining Genome Annotations with Apollo (prokaryotes)',
  'Mass spectrometry imaging 1: Loading and exploring MSI data' =>
    'Mass spectrometry imaging: Loading and exploring MSI data',
  'Trajectory Analysis using Python (Jupyter Notebook) in Galaxy' =>
    'Inferring Trajectories using Python (Jupyter Notebook) in Galaxy',
  'Submitting raw sequencing reads to ENA' => 'Submitting sequence data to ENA'
}

moved_topics_by_tutorial = {
  'Formation of the Super-Structures on the Inactive X' => 'Epigenetics',
  'Identification of the binding sites of the Estrogen receptor' => 'Epigenetics',
  'Identification of the binding sites of the T-cell acute lymphocytic leukemia protein 1 (TAL1)' => 'Epigenetics',
  'RAD-Seq Reference-based data analysis' => 'Ecology',
  'RAD-Seq de-novo data analysis' => 'Ecology',
  'RAD-Seq to construct genetic maps' => 'Ecology',
  'Advanced R in Galaxy' => 'Foundations of Data Science',
  'R basics in Galaxy' => 'Foundations of Data Science',
  'Pre-processing of 10X Single-Cell RNA Datasets' => 'Single Cell',
  'Pre-processing of Single-Cell RNA Data' => 'Single Cell',
  'Generating a single cell matrix using Alevin' => 'Single Cell',
  'Downstream Single-cell RNA analysis with RaceID' => 'Single Cell',
  'Clustering 3K PBMCs with Scanpy' => 'Single Cell',
  'Understanding Barcodes' => 'Single Cell',
  'Filter, Plot and Explore Single-cell RNA-seq Data' => 'Single Cell',
  'Trajectory Analysis using Python (Jupyter Notebook) in Galaxy' => 'Single Cell',
  'Inferring Trajectories using Python (Jupyter Notebook) in Galaxy' => 'Single Cell',
  'Bulk RNA Deconvolution with MuSiC' => 'Single Cell'
}

renamed_topics = {
  'User Interface and Features' => 'Using Galaxy and Managing your Data',
  'Data Manipulation' => 'Using Galaxy and Managing your Data',
  'User Interface and Data Manipulation' => 'Using Galaxy and Managing your Data',
  'Assembly) is not working I can do up to multiQC and after unicycler not working' => 'Assembly'
}

ACCEPTABLE_TOPICS = [
  'Assembly',
  'Climate',
  'Computational chemistry',
  'Contributing to the Galaxy Training Material',
  'Development in Galaxy',
  'Ecology',
  'Epigenetics',
  'Foundations of Data Science',
  'Galaxy Server administration',
  'Genome Annotation',
  'Imaging',
  'Introduction to Galaxy Analyses',
  'Metabolomics',
  'Metagenomics',
  'Proteomics',
  'Sequence analysis',
  'Single Cell',
  'Statistics and machine learning',
  'Teaching and Hosting Galaxy training',
  'Transcriptomics',
  'Using Galaxy and Managing your Data',
  'Variant Analysis',
  'Visualisation'
]

DEAD_TUTORIALS = [
  'RNA-seq counts to genes and pathways',
  'Recording Job Metrics'
]

TEST_PHRASES = [
  'TESTING',
  'test',
  'Helena testing'
]

TITLE_TOPIC = /^(?<tutorial>.*) \((?<topic>.*)\)$/

data = CSV.parse($stdin.read, col_sep: "\t", headers: true, quote_char: '|')
data = data.map do |x|
  x[:timestamp] = x['Timestamp']
  x[:rating] = x['How much did you like this tutorial?']
  x[:pro] = x['What did you like?']
  x[:con] = x['What could be improved?']
  x[:tutorial_topic] = x['Tutorial']

  x.delete('Timestamp')
  x.delete('How much did you like this tutorial?')
  x.delete('What did you like?')
  x.delete('What could be improved?')
  x.delete('Tutorial')
  x.delete('Your feedback is always anonymous. Also make it confidential (only visible to admins)?')
  x
end

data = data.reject { |x| x[:timestamp].nil? }
           # remove rows with NaN on note, pro and con
           .reject { |x| (x[:rating].nil? && x[:pro].nil? && x[:con].nil?) }
           # Ensure that they properly have a tutorial name and title
           # If this returns NIL, then we'll ignore those.
           .select { |x| TITLE_TOPIC.match(x[:tutorial_topic]) }

# Mutate
data.each do |x| # Various cleanups/extractions
  # replace NaN in rating by 0, convert to int.
  x[:rating] = x[:rating].to_i

  # Extract dates
  x[:date] = DateTime.parse(x[:timestamp]).strftime('%Y-%m-%d')
  x[:month] = DateTime.parse(x[:timestamp]).strftime('%Y-%m')

  # Extract tutorial/topic
  m = TITLE_TOPIC.match(x[:tutorial_topic])
  x.delete('tutorial_topic')
  x[:topic] = m[:topic]
  x[:tutorial] = m[:tutorial]

  # Replace N/As with nil
  x[:pro] = x[:pro] == 'N/A' || x[:pro] == 'NA' ? nil : x[:pro]
  x[:con] = x[:con] == 'N/A' || x[:con] == 'NA' ? nil : x[:con]

  # Fix renamed topics
  x[:topic] = renamed_topics.fetch(x[:topic], x[:topic])
  x[:tutorial] = renamed_tutorials.fetch(x[:tutorial], x[:tutorial])
  x[:topic] = moved_topics_by_tutorial.fetch(x[:tutorial], x[:topic])
end
    .select { |x| ACCEPTABLE_TOPICS.include? x[:topic] }
    .reject { |x| DEAD_TUTORIALS.include? x[:tutorial] }
    .select do |x|
  !TEST_PHRASES.include? x[:pro] and !TEST_PHRASES.include? x[:con] and !x[:pro].to_s.include? 'Helena Testing'
end

CSV.open('metadata/feedback.csv', 'wb') do |csv|
  csv << [nil, 'note', 'pro', 'con', 'anonymous', 'tutorial', 'topic', 'month', 'date']
  data.each.with_index do |row, index|
    csv << [index, row[:rating], row[:pro], row[:con], nil, row[:tutorial], row[:topic], row[:month], row[:date]]
  end
end

# feedback_over_month = Hash.new
#
# topics = Dir.glob('topics/*/metadata.yaml').map{|f|
#  d = YAML.load(File.open(f).read)
#  {
#    "title" => d['title'],
#    "name" => d['name'],
#  }
# }
#
# pp topics
