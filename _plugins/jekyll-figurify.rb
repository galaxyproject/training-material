require 'jekyll'

module Jekyll
  class Figurify < Jekyll::Generator

    safe true

    def initialize(config)
      @config = config['figurify'] ||= {}
    end

    def generate(site)
      puts "\n      Applying figurify plugin ✨"
      site.posts.docs.each { |page| figurify page }
      site.pages.each { |page| figurify page }
    end

    private

    def figurify(page)
      num = 0
      page.content = page.content.gsub(/!\[([^\]]*)]\((.+?)\s*(?:"([^"]*)")?\)/) {
        alt = $1
        url = $2
        title = $3
        num += 1

        "<figure>" +
          "<img src=\"#{url}\" alt=\"#{alt}\" title=\"#{title}\">" +
          "<figcaption><span class=\"figcaption-prefix\">#{figcaption_prefix}#{num}:</span> #{title || alt}</figcaption>" +
        "</figure>"
      }
    end

    def figcaption_prefix
      @config['prefix'] || 'Figure '
    end
  end
end
