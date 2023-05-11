require 'jekyll'

module Jekyll
  class Figurify < Jekyll::Generator
    safe true

    def initialize(config)
      @config = config['figurify'] ||= {}
    end

    def generate(site)
      site.pages
          .select { |page| !skip_layout? page.data['layout'] }
          .each { |page| figurify page, site }
      site.posts.docs
          .select { |post| !skip_layout? post.data['layout'] }
          .each { |post| figurify post, site }
    end

    private

    def figurify(page, site)
      num = 0
      return if page.content.nil?

      tuto_dir = File.dirname(page.path)
      page.content = page.content.gsub(/!\[([^\]]*)\]\((.+?)\s*(?:"(.*)")\)({:(.*)})?/) do
        alt = ::Regexp.last_match(1)
        url = ::Regexp.last_match(2)
        title = ::Regexp.last_match(3)
        style = ::Regexp.last_match(5)

        if skip_titles?(title) or (title.to_s.empty? and skip_empty?)
          Regexp.last_match
        else
          num += 1

          alt.gsub!(/"/, '&quot;')
          alt = "#{alt}. " if alt.strip.length > 0 && !(alt.end_with?('.') || alt.end_with?('!') || alt.end_with?('?'))

          dimensions = Gtn::Images.html_image_dimensions(tuto_dir, url)

          prefix = figcaption_prefix(page, site)
          "<figure id=\"figure-#{num}\">" +
            "<img src=\"#{url}\" alt=\"#{alt}\" #{style} #{dimensions} loading=\"lazy\">" +
            "<figcaption><span class=\"figcaption-prefix\"><strong>#{prefix}#{num}</strong>:</span> #{title}</figcaption>" +
            '</figure>'
        end
      end

      page.content = page.content.gsub(/!\[([^\]]*)\]\((.+?)?\)({:(.*)})?/) do
        alt = ::Regexp.last_match(1)
        url = ::Regexp.last_match(2)
        style = ::Regexp.last_match(4)

        alt.gsub!(/"/, '&quot;')
        alt = "#{alt}. " if alt.strip.length > 0 && !(alt.end_with?('.') || alt.end_with?('!') || alt.end_with?('?'))
        "<img src=\"#{url}\" alt=\"#{alt}\" #{style} loading=\"lazy\">"
      end
    end

    def figcaption_prefix(page, site)
      fig = 'Figure'
      if page['lang']
        lang = page['lang']
        fig = site.data['lang'][lang]['figure']
      end
      @config['prefix'] || (fig + ' ')
    end

    def skip_empty?
      @config['skip_empty'] || false
    end

    def skip_layout?(layout)
      to_skip = @config['skip_layouts'] || []

      true if to_skip.empty?

      to_skip.include?(layout)
    end

    def skip_titles?(title)
      to_skip = @config['skip_titles'] || []
      to_skip.include?(title)
    end
  end
end
