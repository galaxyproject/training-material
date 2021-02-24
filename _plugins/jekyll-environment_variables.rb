# Plugin to add environment variables to the `site` object in Liquid templates

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

      # Add other environment variables to `site.config` here...
    end

  end

end
