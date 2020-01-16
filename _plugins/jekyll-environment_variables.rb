# Plugin to add environment variables to the `site` object in Liquid templates

module Jekyll

  class EnvironmentVariablesGenerator < Generator

    def generate(site)
      site.config['google_analytics_id'] = ENV['GOOGLE_ANALYTICS_ID']
      # Add other environment variables to `site.config` here...
    end

  end

end