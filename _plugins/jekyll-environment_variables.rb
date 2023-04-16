# Plugin to add environment variables to the `site` object in Liquid templates
require 'find'
require 'bibtex'
require 'citeproc/ruby'
require 'csl/styles'
require 'time'
require './_plugins/gtn/scholar'

module Jekyll

  class EnvironmentVariablesGenerator < Generator

    ##
    # Environment variables are added to the `site` object in Liquid templates.
    # Here we add the following:
    #  - `site.config['git_revision']` - the current git revision
    #  - `site.config['git_tags']` - an array of all git tags
    #  - `site.config['git_tags_recent']` - an array of the 3 most recent git tags
    #  - `site.config['gtn_fork']` - the fork of the GTN repo
    #  - `site.config['age']` - the age of the site in years
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
      Gtn::Scholar.load_bib(site)

      # Get site age.
      first_commit = Date.parse("2015-06-29")
      today = Date.today()

      site.config['age'] = (today - first_commit).to_f / 365
    end
  end
end
