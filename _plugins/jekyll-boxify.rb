require 'jekyll'

module Jekyll
  class Figurify < Jekyll::Generator

    ICONS = {
      "tip" => "far fa-lightbulb",
      "comment" => "far fa-comment-dots",
    }

    BOX_TITLES = {
      "tip" => "Tip",
      "comment" => "Comment",
    }

    def initialize(config)
      @config = config['boxify'] ||= {}
    end

    def generate(site)
      site.pages
        .select { |page| not skip_layout? page.data['layout'] }
        .each { |page| boxify page,site }
      site.posts.docs
        .select { |post| not skip_layout? post.data['layout'] }
        .each { |post| boxify post, site }
    end

    def get_icon(icon)
       if icon.start_with?("fa")
        %Q(<i class="#{icon}" aria-hidden="true"></i><span class="visually-hidden">#{@text}</span>)
       elsif icon.start_with?("ai")
        %Q(<i class="ai #{icon}" aria-hidden="true"></i><span class="visually-hidden">#{@text}</span>)
       end
    end

    def boxify(page, site)
      if page.content.nil?
        return
      end
      count = 0

      page.content = page.content.gsub(/<(tip|comment) title="([^"]*)">/) {
        box_type = $1
        title = $2

        count += 1

        "<div class=\"box #{box_type}\">" +
          "<button class=\"box-title\" type=\"button\" aria-controls=\"box-no-#{count}\" aria-expanded=\"true\" aria-label=\"Toggle #{box_type} box: #{title}\">" + 
            "#{get_icon(ICONS[box_type])} #{BOX_TITLES[box_type]}: #{title}" + 
            "<span role='button' class='fold-unfold fa fa-plus-square'></span>" + 
          "</button>" +
          "<div id=\"box-no-#{count}\" class=\"box-content\" markdown=1>"
      }

      page.content = page.content.gsub(/<\/\s*(tip|comment)\s*>/) {
        "</div>" +  # box-content
        "</div>" # box itself
      }
    end
  end
end
