# frozen_string_literal: true
require 'zip'


module Gtn
  # Parse the git repo to get some facts
  module RoCrate
    GLOBAL_WORKFLOW_OFFSET = 0

    def self.cache
      @@cache ||= Jekyll::Cache.new('RoCrate')
    end

    def self.name2md(site, name)
      return %Q([#{Gtn::Contributors.fetch_name(site, name)}](https://training.galaxyproject.org/training-material/hall-of-fame/#{name}/))
    end

    def self.write(site, dir, associated_material, workflow, url, baseurl, time_based_version: false)
      wfname = workflow['wfname']
      # {"workflow"=>"galaxy-workflow-mouse_novel_peptide_analysis.ga",
      # "tests"=>false,
      # "url"=>
      # "http://0.0.0.0:4002/training-material/topics/.../workflows/galaxy-workflow-mouse_novel_peptide_analysis.ga",
      # "path"=>
      # "topics/proteomics/tutorials/.../galaxy-workflow-mouse_novel_peptide_analysis.ga",
      # "wfid"=>"proteomics-proteogenomics-novel-peptide-analysis",
      # "wfname"=>"galaxy-workflow-mouse_novel_peptide_analysis",
      # "trs_endpoint"=>
      # "http://0.0.0.0:4002/training-material/api/.../versions/galaxy-workflow-mouse_novel_peptide_analysis",
      # "license"=>nil,
      # "creators"=>[],
      # "name"=>"GTN Proteogemics3 Novel Peptide Analysis",
      # "test_results"=>nil,
      # "modified"=>2023-06-07 12:09:36.12 +0200}

      wfdir = File.join(dir, workflow['topic_id'], workflow['tutorial_id'], wfname)
      FileUtils.mkdir_p(wfdir)
      path = File.join(wfdir, 'ro-crate-metadata.json')
      Jekyll.logger.debug "[GTN/API/WFRun] Writing #{path}"
      # We have the `dot` graph code in a variable, we need to pass it to `dot -T png ` on the stdin
      dot_path = File.join(wfdir, "graph.dot")
      File.write(dot_path, workflow['graph_dot'])
      Jekyll.logger.debug "[GTN/API/WFRun] dot -T png #{dot_path} > graph.png"
      `dot -T png '#{dot_path}' > '#{File.join(wfdir, 'graph.png')}'`

      # Our new workflow IDs
      wfurlid = url + baseurl + '/' + workflow['path'].gsub(/.ga$/, '.html')

      uuids = workflow['creators'].map do |c|
        if c.key?('identifier') && !c['identifier'].empty?
          if c['identifier'].start_with?('http')
            c['identifier']
          else
            "https://orcid.org/#{c['identifier']}"
          end
        else
          "##{SecureRandom.uuid}"
        end
      end
      author_uuids = uuids.map { |u| { '@id' => u.to_s } }
      author_linked = workflow['creators'].map.with_index do |c, i|
        {
          '@id' => (uuids[i]).to_s,
          '@type' => c['class'],
          'name' => c['name'],
        }
      end
      wf_ga = JSON.parse(File.read(workflow['path']))

      license = workflow['license'] ? "https://spdx.org/licenses/#{workflow['license']}" : 'https://spdx.org/licenses/CC-BY-4.0'

      version = time_based_version ? Time.now.to_i.to_s : Gtn::ModificationTimes.obtain_modification_count(workflow['path']).to_s + ".#{GLOBAL_WORKFLOW_OFFSET}"

      features = {
        'Includes [Galaxy Workflow Tests](https://training.galaxyproject.org/training-material/faqs/gtn/workflow_run_test.html)' => workflow['tests'],
        'Includes a [Galaxy Workflow Report](https://training.galaxyproject.org/training-material/faqs/galaxy/workflows_report_view.html)' => workflow['features']['report'],
        'Uses [Galaxy Workflow Comments](https://training.galaxyproject.org/training-material/faqs/galaxy/workflows_comments.html)' => workflow['features']['comments'],
        'Uses [subworkflows](https://training.galaxyproject.org/training-material/faqs/galaxy/workflows_subworkflows.html)' => workflow['features']['subworkflows'],
      }

      mat_contribs = [
        ['Workflow Author(s)', workflow['creators'].map { |c| c['name'] }],
        ['Tutorial Author(s)', Gtn::Contributors.get_authors(associated_material).map { |n| name2md(site, n) }],
        ['Tutorial Contributor(s)', Gtn::Contributors.get_non_authors(associated_material).map { |n| name2md(site, n) }],
        ['Funder(s)', Gtn::Contributors.get_funders(associated_material).map { |n| name2md(site, n) }],
      ].reject { |_, v| v.empty? }

      description = %Q(
#{wf_ga['annotation']}

## Associated Tutorial

This workflows is part of the tutorial [#{associated_material['title']}](#{url}#{baseurl}/topics/#{workflow['topic_id']}/tutorials/#{workflow['tutorial_id']}/tutorial.html), available in the [GTN](https://training.galaxyproject.org)

#{"## Features" if features.values.any?}

#{features.select { |_, v| v }.keys.map { |f| "* #{f}" }.join("\n")}

## Thanks to...

#{mat_contribs.map { |k, v| "**#{k}**: #{v.join(', ')}" }.join("\n\n")}

[![gtn star logo followed by the word workflows](https://training.galaxyproject.org/training-material/assets/branding/gtn-workflows.png)](https://training.galaxyproject.org/training-material/)
      ).strip

      crate = {
        '@context' => ['https://w3id.org/ro/crate/1.1/context'],
        '@graph' => [
          {
            '@id': 'ro-crate-metadata.json',
            '@type': 'CreativeWork',
            about: {
              '@id': wfurlid,
            },
            conformsTo: {
              '@id': 'https://w3id.org/ro/crate/1.1'
            },
          },
          {
            '@id': wfurlid,
            '@type': 'Dataset',
            name: workflow['name'],
            description: description,
            version: version,
            license: license,
            datePublished: workflow['modified'].strftime('%Y-%m-%dT%H:%M:%S.%L%:z'),
            # hasPart: [
            #   {
            #     '@id': '#assembly-assembly-quality-control'
            #   }
            # ],
            mainEntity: {
              '@id': "#{wfname}.ga"
            },
            hasPart: [
              {
              '@id': "#{wfname}.ga"
              },
              {
              '@id': "graph.png"
              },
            ]
          },
          {
            '@id': "#{wfname}.ga",
            '@type': %w[
              File
              SoftwareSourceCode
              ComputationalWorkflow
            ],
            author: author_uuids,
            name: workflow['name'],
            programmingLanguage: {
              '@id': 'https://w3id.org/workflowhub/workflow-ro-crate#galaxy'
            },
            image: {
              '@id': 'graph.png'
            }
          },
          {
            '@id': 'graph.png',
            '@type': [
              'File',
              'ImageObject',
              'WorkflowSketch'
            ],
            contentSize: File.size(File.join(wfdir, 'graph.png')),
          },
          {
            '@id': 'https://w3id.org/workflowhub/workflow-ro-crate#galaxy',
            '@type': 'ComputerLanguage',
            identifier: {
              '@id': 'https://galaxyproject.org/'
            },
            name: 'Galaxy',
            url: {
              '@id': 'https://galaxyproject.org/'
            },
            version: '23.1'
          }
        ]
      }
      crate['@graph'] += author_linked
      File.write(path, JSON.pretty_generate(crate))

      wf_ga['tags'].map! { |t| t.gsub('^name:', '').capitalize }
      wf_ga['tags'].push('GTN')
      wf_ga['tags'].push('Galaxy')
      # TODO: Remove this after https://github.com/seek4science/seek/issues/1927
      wf_ga.delete('creator')

      File.write(File.join(wfdir, 'mod.ga'), JSON.pretty_generate(wf_ga))

      zip_path = File.join(wfdir, 'rocrate.zip')
      Jekyll.logger.info "[GTN/API/WFRun] Zipping #{zip_path}"
      if File.exist?(zip_path)
        File.delete(zip_path)
      end
      Zip::File.open(zip_path, create: true) do |zipfile|
        # - The name of the file as it will appear in the archive
        # - The original file, including the path to find it
        zipfile.add('ro-crate-metadata.json', path)
        zipfile.add('graph.png', File.join(wfdir, 'graph.png'))
        zipfile.add("#{wfname}.ga", File.join(wfdir, 'mod.ga'))
      end
    end
  end
end
