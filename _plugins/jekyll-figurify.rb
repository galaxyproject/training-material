# frozen_string_literal: true

require 'jekyll'

module Jekyll
  # Our modifications to the markdown renderer to process images with figure captions
  class Figurify < Jekyll::Generator
    safe true

    def initialize(config)
      super
      @config = config['figurify'] ||= {}
    end

    def generate(site)
      site.pages
          .reject { |page| skip_layout? page.data['layout'] }
          .each { |page| figurify page, site }
      site.posts.docs
          .reject { |post| skip_layout? post.data['layout'] }
          .each { |post| figurify post, site }
    end

    private

    def insert_image(url, alt, style, dimensions, actual_path)
      # If it's a *local* SVG (we don't want to do this with remote SVGs, doesn't work right)
      if url =~ (/svg$/) && !actual_path.nil?
        fallback = ''
        if actual_path.nil?
          # External image, no fallback possible
          fallback = ''
        elsif File.exist?(actual_path.gsub(/svg$/, 'png'))
          fallback = "<img src=\"#{url.gsub(/svg$/, 'png')}\" alt=\"#{alt}\">"
        elsif File.exist?(actual_path.gsub(/svg$/, 'jpg'))
          fallback = "<img src=\"#{url.gsub(/svg$/, 'jpg')}\" alt=\"#{alt}\">"
        elsif File.exist?(actual_path.gsub(/svg$/, 'jpeg'))
          fallback = "<img src=\"#{url.gsub(/svg$/, 'jpeg')}\" alt=\"#{alt}\">"
        end

        %(
          <div style="overflow-x: auto">
          <object data="#{url}" #{style} type="image/svg+xml" alt="#{alt}">
            #{fallback}
            #{alt}
          </object>
          </div>
        )
      else
        %(
          <img src="#{url}"  alt="#{alt}" #{style} #{dimensions} loading="lazy">
        )
      end
    end

    def figurify(page, site)
      num_figure = 0
      return if page.content.nil?

      tuto_dir = File.dirname(page.path)
      page.content = page.content.gsub(/!\[([^\]]*)\]\((.+?)\s*(?:"(.*)")\)({:(.*)})?/) do
        alt = ::Regexp.last_match(1)
        url = ::Regexp.last_match(2)
        title = ::Regexp.last_match(3)
        style = ::Regexp.last_match(5)

        if skip_titles?(title) || (title.to_s.empty? && skip_empty?)
          Regexp.last_match
        else
          num_figure += 1

          alt.gsub!(/"/, '&quot;')
          if alt.strip.length.positive? && !(alt.end_with?('.') || alt.end_with?('!') || alt.end_with?('?'))
            alt = "#{alt}. "
          end

          dimensions, actual_path = Gtn::Images.html_image_dimensions(tuto_dir, url)
          prefix = figcaption_prefix(page, site)
          image = insert_image(url, alt, style, dimensions, actual_path)

          %(
            <figure id="figure-#{num_figure}" style="max-width: 90%;">
              #{image}
              <a target="_blank" href="#{url}" rel="noopener noreferrer"><small>Open image in new tab</small></a><br/><br/>
              <figcaption>
                <span class="figcaption-prefix"><strong>#{prefix}#{num_figure}</strong>:</span> #{title}
              </figcaption>
            </figure>
          ).split("\n").map(&:strip).join
        end
      end

      page.content = page.content.gsub(/!\[([^\]]*)\]\((.+?)?\)({:(.*)})?/) do
        alt = ::Regexp.last_match(1)
        url = ::Regexp.last_match(2)
        style = ::Regexp.last_match(4)

        dimensions, _actual_path = Gtn::Images.html_image_dimensions(tuto_dir, url)

        alt.gsub!(/"/, '&quot;')
        if alt.strip.length.positive? && !(alt.end_with?('.') || alt.end_with?('!') || alt.end_with?('?'))
          alt = "#{alt}. "
        end

        %(
        <a href="#{url}" rel="noopener noreferrer">
          <img src="#{url}"  alt="#{alt}" #{style} #{dimensions} loading="lazy">
        </a>
        ).split("\n").map(&:strip).join
      end
    end

    def figcaption_prefix(page, site)
      fig = 'Figure'
      if page['lang']
        lang = page['lang']
        fig = site.data['lang'][lang]['figure']
      end
      @config['prefix'] || "#{fig} "
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
