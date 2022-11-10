# Plugin to add environment variables to the `site` object in Liquid templates
require 'find'
require 'bibtex'
require 'citeproc/ruby'
require 'csl/styles'
require 'time'

module Jekyll

  class EnvironmentVariablesGenerator < Generator

    def generate(site)

      begin
        git_head = File.open(File.join('.git', 'HEAD')).read.strip.split(' ')[1]
        git_head_ref = File.open(File.join('.git', git_head)).read.strip
      rescue
        git_head_ref = 'none'
      end
      site.config['git_revision'] = git_head_ref
      site.config['gtn_fork'] = ENV['GTN_FORK']
      begin
        tags = `git tag -l`.strip.split.sort
        site.config['git_tags'] = tags
        site.config['git_tags_recent'] = tags.reverse[0..2]
      rescue
        site.config['git_tags'] = []
        site.config['git_tags_recent'] = []
      end

      # Add other environment variables to `site.config` here...

      puts "[GTN/scholar] Creating global bib cache"
      global_bib = BibTeX::Bibliography.new
      bib_paths = [Find.find('./topics'), Find.find('./faqs')].lazy.flat_map(&:lazy)
      bib_paths.each{|path|
        if FileTest.directory?(path)
          if File.basename(path).start_with?('.')
            Find.prune       # Don't look any further into this directory.
          else
            next
          end
        else
          if path =~ /bib$/ then
            for x in BibTeX.open(path)
              x = x.convert_latex
              global_bib << x
            end
          end
        end
      }
      site.config['cached_global_bib'] = global_bib
      puts "[GTN/scholar] Done"
      style = CSL::Style.load("_layouts/g3.csl")
      cp = CiteProc::Processor.new style: style,
                                   format: 'html', locale: 'en'
      cp.import global_bib.to_citeproc
      site.config['cached_citeproc'] = cp

      # Get site age.
      first_commit = Date.parse("2015-06-29")
      today = Date.today()

      site.config['age'] = (today - first_commit).to_f / 365
    end
  end
end
