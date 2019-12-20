require 'jekyll'

module Jekyll
  class Figurify < Jekyll::Generator

    safe true

    def initialize(config)
      @config = config['figurify'] ||= {}
    end

    def generate(site)
      site.pages
        .select { |page| not skip_layout? page.data['layout'] }
        .each { |page| figurify page }
      site.posts.docs
        .select { |post| not skip_layout? post.data['layout'] }
        .each { |post| figurify post }
    end

    private

    def figurify(page)
      num = 0
      page.content = page.content.gsub(/!\[([^\]]*)\]\((.+?)\s*(?:"(.*)")?\)({:(.*)})?/) {
        alt = $1
        url = $2
        title = $3
        style = $5

        if skip_titles?(title) or (title.to_s.empty? and skip_empty?)
          Regexp.last_match
        else
          num += 1

          "<figure id=\"figure-#{num}\">" +
            "<img src=\"#{url}\" alt=\"#{alt}\" #{style}>" +
            "<figcaption><span class=\"figcaption-prefix\">#{figcaption_prefix}#{num}:</span> #{title}</figcaption>" +
          "</figure>"
        end
      }
    end

    def figcaption_prefix
      @config['prefix'] || 'Figure '
    end

    def skip_empty?
      @config['skip_empty'] || false
    end

    def skip_layout?(layout)
      to_skip = @config['skip_layouts'] || []

      if to_skip.empty?
        true
      end

      to_skip.include?(layout)
    end

    def skip_titles?(title)
      to_skip = @config['skip_titles'] || []
      to_skip.include?(title)
    end
  end
end
