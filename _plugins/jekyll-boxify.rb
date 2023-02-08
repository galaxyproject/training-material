require 'jekyll'
require './_plugins/gtn.rb'

module Jekyll
  class Boxify < Jekyll::Generator
    def initialize(config)
      @config = config['boxify'] ||= {}
    end

    def generate(site)
      puts "[GTN/Boxify]"
      site.pages.each { |page| boxify page,site }
      site.posts.docs.each { |post| boxify post, site }
    end

    def boxify(page, site)
      if page.content.nil?
        return
      end

      if page['lang']
        lang = page['lang']
      else
        lang = "en"
      end

      # Interim solution, fancier box titles
      page.content = page.content.gsub(/<(#{Gtn::Boxify.box_classes})-title>(.*?)<\/\s*\1-title\s*>/) {
        box_type = $1
        title = $2
        if page.data['citation_target'] == 'jupyter'
          title = Gtn::Boxify.safe_title(title)
          title = Gtn::Boxify.format_box_title(title, box_type, lang=lang)
          icon = Gtn::Boxify.get_icon(box_type, emoji: true)
          box = "<div class=\"box-title\" aria-label=\"#{box_type} box: #{title}\" style=\"font-size: 150%\">#{icon} #{title}</div>"
          box.gsub!(/\\&quot/, '&quot')
          box.gsub!(/([^\\])"/, '\1\\"')
        else
          _, box = Gtn::Boxify.generate_title(box_type, title, lang, page.path)
        end

        box
      }

      # Long term solution, proper new boxes
      # BUT: does not work with <details></details> that are actual HTML elements, so we'll need to rename those.
      #page.content = page.content.gsub(/<(#{Gtn::Boxify.box_classes})>/) {
        #box_type = $1
        #box = Gtn::Boxify.generate_box(box_type, nil, lang, page.path)
        #box
      #}

      #page.content = page.content.gsub(/<(#{Gtn::Boxify.box_classes}) title="([^"]*)">/) {
        #box_type = $1
        #title = $2
        #box = Gtn::Boxify.generate_box(box_type, title, lang, page.path)
        #box
      #}

      #page.content = page.content.gsub(/<\/\s*(#{Gtn::Boxify::box_classes})\s*>/) {
        #box_type = $1
        #"\n</div></div><!--#{box_type}-->"
      #}
    end
  end
end
