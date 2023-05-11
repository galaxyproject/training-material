require 'jekyll'
# The fastimage gem may not be installed, and that is OK, we will fall back to another method
FASTIMAGE_AVAILABLE = true
begin
  require 'fastimage'
rescue LoadError
  Jekyll.logger.info '[GTN/Images] Could not load fastimage gem, disabling feature (probably due to conda)'
  FASTIMAGE_AVAILABLE = false
end

module Gtn
  module Images
    def self.cache
      @@cache ||= Jekyll::Cache.new('ImageDimensions')
    end

    def self.html_image_dimensions(tuto_dir, url)
      if not FASTIMAGE_AVAILABLE
        return ''
      end

      width, height = self.get_image_dimensions(tuto_dir, url)
      if width && height
        %Q(width="#{width}" height=#{height})
      end
    end

    def self.get_image_dimensions(tuto_dir, url)
      if match = url.match(/^{{\s*site.baseurl\s*}}\/(.*)/)
        self._get_image_dimensions(match[1].strip)
      elsif match = url.match(/{%\s*link\s*(.*)\s*%}/)
        self._get_image_dimensions(match[1].strip)
      elsif !url.match(/https?:\/\//)
        img_path = File.absolute_path(File.join(tuto_dir, url))
        if File.exist?(img_path)
          self._get_image_dimensions(img_path)
        end
      else
        self._get_image_dimensions(img_path)
      end
    end

    def self._get_image_dimensions(path)
      self.cache.getset(path) do
        begin
          FastImage.size(path)
        rescue
          puts "Could not resolve size of #{path}"
        end
      end
    end
  end
end

if $0 == __FILE__
end
