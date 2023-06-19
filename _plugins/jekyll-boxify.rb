# frozen_string_literal: true

require 'jekyll'
require './_plugins/gtn'

module Jekyll
  # The GTN Box generation process
  class Boxify < Jekyll::Generator
    def initialize(config) # :nodoc:
      super
      @config = config['boxify'] ||= {}
    end

    def generate(site) # :nodoc:
      puts '[GTN/Boxify]'
      site.pages.each { |page| boxify page, site }
      site.posts.docs.each { |post| boxify post, site }
    end

    ##
    # This function adds boxes to the page content.
    # Params:
    # +page+:: The page to add boxes to
    # +site+:: The +Jekyll::Site+ object
    def boxify(page, _site)
      return if page.content.nil?

      lang = page['lang'] || 'en'

      # Interim solution, fancier box titles
      # rubocop:disable Layout/LineLength
      page.content = page.content.gsub(%r{<(?<boxclass>#{Gtn::Boxify.box_classes})-title( ?(?<noprefix>noprefix))>(?<title>.*?)</\s*\k<boxclass>-title\s*>}) do
        # rubocop:enable Layout/LineLength
        m = ::Regexp.last_match
        box_type = m[:boxclass]
        title = m[:title]
        noprefix = m[:noprefix]
        if page.data['citation_target'] == 'jupyter'
          title = Gtn::Boxify.safe_title(title)
          title = Gtn::Boxify.format_box_title(title, box_type, lang, noprefix: noprefix)
          icon = Gtn::Boxify.get_icon(box_type, emoji: true)
          box = "<div class=\"box-title\" aria-description=\"#{box_type} box: " \
                "#{title}\" style=\"font-size: 150%\">#{icon} #{title}</div>"
          box.gsub!(/\\&quot/, '&quot')
          box.gsub!(/([^\\])"/, '\1\\"')
        else
          _, box = Gtn::Boxify.generate_title(box_type, title, lang, page.path, noprefix: noprefix)
        end

        box
      end

      # Long term solution, proper new boxes
      # BUT: does not work with <details></details> that are actual HTML elements, so we'll need to rename those.
      # page.content = page.content.gsub(/<(#{Gtn::Boxify.box_classes})>/) {
      # box_type = $1
      # box = Gtn::Boxify.generate_box(box_type, nil, lang, page.path)
      # box
      # }

      # page.content = page.content.gsub(/<(#{Gtn::Boxify.box_classes}) title="([^"]*)">/) {
      # box_type = $1
      # title = $2
      # box = Gtn::Boxify.generate_box(box_type, title, lang, page.path)
      # box
      # }

      # page.content = page.content.gsub(/<\/\s*(#{Gtn::Boxify::box_classes})\s*>/) {
      # box_type = $1
      # "\n</div></div><!--#{box_type}-->"
      # }
    end
  end
end
