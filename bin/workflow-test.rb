#!/usr/bin/env ruby
# frozen_string_literal: true

require 'open3'
require 'json'

GALAXIES = {
  eu: { url: 'https://usegalaxy.eu', key: ENV.fetch('GALAXY_EU_KEY', 'NONE') },
}

def test_workflow(workflow_file, galaxy_id)
  directory = File.dirname(workflow_file)
  workflow_base = File.basename(workflow_file, '.ga')
  workflow_output_json = File.join(directory, "#{workflow_base}.#{galaxy_id}.json")
  galaxy_url = GALAXIES[galaxy_id][:url]
  galaxy_user_key = GALAXIES[galaxy_id][:key]
  cmd = [
    'planemo', '--verbose', 'test',
    '--galaxy_url', galaxy_url,
    '--galaxy_user_key', galaxy_user_key,
    '--no_shed_install',
    '--engine', 'external_galaxy',
    '--polling_backoff', '1',
    '--simultaneous_uploads',
    '--test_output_json', workflow_output_json,
    workflow_file
  ]
  p cmd.join(' ')

  Open3.popen3(*cmd) do |_stdin, stdout, stderr, wait_thr|
    exit_status = wait_thr.value # Process::Status object returned
    File.write("#{directory}/#{workflow_base}.#{galaxy_id}.log", stdout.read)
    File.write("#{directory}/#{workflow_base}.#{galaxy_id}.err", stderr.read)
    puts "#{workflow_file} => #{exit_status} (#{stderr})"
  end
end

# rubocop:disable Layout/LineLength
workflows = %w[
  ./topics/assembly/tutorials/debruijn-graph-assembly/workflows/debruijn-graph.ga
  ./topics/assembly/tutorials/general-introduction/workflows/assembly-general-introduction.ga
  ./topics/assembly/tutorials/unicycler-assembly/workflows/unicycler.ga
  ./topics/computational-chemistry/tutorials/analysis-md-simulations/workflows/main_workflow.ga
  ./topics/computational-chemistry/tutorials/md-simulation-gromacs/workflows/main_workflow.ga
  ./topics/ecology/tutorials/obisindicators/workflows/Obis-indicators.ga
  ./topics/epigenetics/tutorials/tal1-binding-site-identification/workflows/tal1-binding-site-identification-workflow.ga
  ./topics/fair/tutorials/ro-crate-submitting-life-monitor/sort-and-change-case-workflow/sort-and-change-case.ga
  ./topics/genome-annotation/tutorials/annotation-with-prokka/workflows/Galaxy-Workflow-Workflow_constructed_from_history__prokka-workflow_.ga
  ./topics/genome-annotation/tutorials/repeatmasker/workflows/RepeatMasker.ga
  ./topics/introduction/tutorials/galaxy-intro-101/workflows/galaxy-intro-101-workflow.ga
  ./topics/introduction/tutorials/galaxy-intro-101-everyone/workflows/main_workflow.ga
  ./topics/introduction/tutorials/galaxy-intro-peaks2genes/workflows/Galaxy-Introduction-Peaks2Genes-Part-1-Workflow.ga
  ./topics/introduction/tutorials/galaxy-intro-short/workflows/Galaxy-Workflow-galaxy-intro-short.ga
  ./topics/sequence-analysis/tutorials/quality-control/workflows/quality_control.ga
  ./topics/metagenomics/tutorials/beer-data-analysis/workflows/main_workflow.ga
  ./topics/metagenomics/tutorials/general-tutorial/workflows/amplicon.ga
  ./topics/metagenomics/tutorials/pathogen-detection-from-nanopore-foodborne-data/workflows/Nanopore_Datasets_Pathogen_Tracking_among_all_samples.ga
  ./topics/single-cell/tutorials/scatac-preprocessing-tenx/workflows/scATAC-seq-Count-Matrix-Filtering.ga
  ./topics/single-cell/tutorials/scatac-preprocessing-tenx/workflows/scATAC-seq-FASTQ-to-Count-Matrix.ga
  ./topics/single-cell/tutorials/scrna-case_alevin/workflows/Generating-a-single-cell-matrix-using-Alevin-1.9.ga
  ./topics/single-cell/tutorials/scrna-case_monocle3-trajectories/workflows/Trajectory-analysis-using-Monocle3---full-tutorial-workflow.ga
  ./topics/statistics/tutorials/machinelearning/workflows/machine_learning.ga
  ./topics/transcriptomics/tutorials/rna-interactome/workflows/rna-rna-interactome-data-analysis-chira.ga
  ./topics/transcriptomics/tutorials/rna-seq-viz-with-heatmap2/workflows/rna-seq-viz-with-heatmap2.ga
  ./topics/transcriptomics/tutorials/small_ncrna_clustering/workflows/blockclust_clustering.ga
  ./topics/variant-analysis/tutorials/non-dip/workflows/Calling_variants_in_non-diploid_systems.ga
]
# rubocop:enable Layout/LineLength

threads = []

workflows.each do |workflow|
  threads << Thread.new do
    test_workflow(workflow, :eu)
  end
end
threads.map(&:join)
