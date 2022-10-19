require 'jekyll'

module Gtn
  module Boxify
    @@ICONS = {
      "agenda" => "",
      "code-in" => "far fa-keyboard",
      "code-out" => "fas fa-laptop-code",
      "comment" => "far fa-comment-dots",
      "details" => "fas fa-info-circle",
      "feedback" => "far fa-comments",
      "hands-on" => "fas fa-pencil-alt",
      "hands_on" => "fas fa-pencil-alt",
      "hidden" => "",
      "matrix" => "",
      "overview" => "",
      "question" => "far fa-question-circle",
      "quote" => "",
      "solution" => "far fa-eye",
      "spoken" => "",
      "tip" => "far fa-lightbulb",
      "warning" => "fas fa-exclamation-triangle",
    }

    @@ICONS_EMOJI = {
      'tip' => '💡',
      'code-in' => '⌨️',
      'code-out' => '🖥',
      'question' => '❓',
      'solution' => '👁',
      'warning' => '⚠️',
      'comment' => '💬',
      'feedback' => '⁉️',
      'details' => '💬',
      'hands_on' => '✏️',
    }

    @@BOX_TITLES = {
      "en" => {
        "agenda" => "Agenda",
        "code-in" => "Input",
        "code-out" => "Output",
        "comment" => "Comment",
        "details" => "Details",
        "hands-on" => "Hands-on",
        "hands_on" => "Hands-on",
        "question" => "Question",
        "solution" => "Solution",
        "tip" => "Tip",
        "warning" => "Warning",
      },
      "es" => {
        # Just google translated these.
        "agenda" => "Agenda",
        "code-in" => "Entrada",
        "code-out" => "Salida",
        "comment" => "Comentario",
        "details" => "Detalles",
        "solution" => "Solución",
        "warning" => "Aviso",

        # The only ones we have translations for??
        "hands-on" => "Práctica",
        "hands_on" => "Práctica",
        "question" => "Preguntas",
        "tip" => "Consejo",
      }
    }

    @title_unique_offsets = Hash.new

    @@COLLAPSIBLE_BOXES = [
      "details", "solution", "tip",
    ]

    @@BOX_CLASSES = @@ICONS.keys.join "|"
    @@TITLE_CLASSES = @@ICONS.keys.map{|x| "#{x}-title" }.join "|"

    def self.box_classes
      @@BOX_CLASSES
    end

    def self.title_classes
      @@TITLE_CLASSES
    end

    def self.get_icon(icon_cls, emoji: false)
      if emoji
        return @@ICONS_EMOJI.fetch(icon_cls, '')
      end

      icon = @@ICONS[icon_cls]
      if !icon.nil?
       if icon.start_with?("fa")
        %Q(<i class="#{icon}" aria-hidden="true"></i><span class="visually-hidden">#{@text}</span>)
       elsif icon.start_with?("ai")
        %Q(<i class="ai #{icon}" aria-hidden="true"></i><span class="visually-hidden">#{@text}</span>)
       end
      else
        %Q(<span class="visually-hidden"></span>)
      end
    end

    def self.get_id(box_type, title, path)
      box_id = "#{box_type}-#{title}"
      box_safe = Jekyll::Utils::slugify(box_id)
      key = "#{path}|#{box_type}|#{box_safe}"

      if not @title_unique_offsets.has_key?(key)
        @title_unique_offsets[key] = 0
        box_safe_final = box_safe
      else
        box_safe_final = "#{box_safe}-#{@title_unique_offsets[key]}"
      end
      @title_unique_offsets[key] += 1

      box_safe_final
    end

    def self.format_box_title(title, box_type, lang="en")
      title_fmted = ((!title.nil?) && title.length > 0 ? ": #{title}" : "")
      "#{@@BOX_TITLES[lang][box_type]}#{title_fmted}"
    end

    def self.generate_collapsible_title(box_type, title, lang="en", key)
      box_id = self.get_id(box_type, title, key)
      box_title = self.format_box_title(title, box_type, lang=lang)
      return [box_id, %Q(
        <div id="#{box_id}" class="box-title">
        <button type="button" aria-controls="#{box_id}-contents" aria-expanded="true" aria-label="Toggle #{box_type} box: #{title}">
          #{self.get_icon(box_type)} #{box_title}
          <span role="button" class="fold-unfold fa fa-minus-square"></span>
        </button>
        </div>
      ).split(/\n/).map{|x| x.lstrip.rstrip}.join("").lstrip.rstrip]
    end

    def self.generate_static_title(box_type, title, lang="en", key)
      box_id = self.get_id(box_type, title, key)
      box_title = self.format_box_title(title, box_type, lang=lang)

      if title.nil?
        puts "Static | typ=#{box_type} | t=#{title} | l=#{lang} | k=#{key}"
      end

      return [box_id, %Q(
        <div id="#{box_id}" class="box-title" aria-label="#{box_type} box: #{title}">
          #{self.get_icon(box_type)} #{box_title}
        </div>
      ).split(/\n/).map{|x| x.lstrip.rstrip}.join("").lstrip.rstrip]
    end

    def self.safe_title(title)
      if title.nil?
        title = ""
      end

      title = title.gsub(/"/, '&quot;')
      title
    end

    def self.generate_title(box_type, title, lang="en", key)
      title = self.safe_title(title)
      if @@COLLAPSIBLE_BOXES.include?(box_type)
        self.generate_collapsible_title(box_type, title, lang, key)
      else
        self.generate_static_title(box_type, title, lang, key)
      end
    end

    def self.generate_box(box_type, title, lang="en", key)
      title = self.safe_title(title)
      box_id, box_title = generate_title(box_type, title, lang, key)
      return %Q(
        <div class="box #{box_type}" markdown=0>
        #{box_title}
        <div id=#{box_id}-contents class="box-content" markdown=1>
      ).split(/\n/).map{|x| x.lstrip.rstrip}.join("").lstrip.rstrip
    end


  end
end

if $0 == __FILE__
  require 'test/unit'
  class BoxIdTest < Test::Unit::TestCase
    def test_single_page
      assert_equal(Gtn::Boxify::get_id('hands-on', 'a box', 'index.md'), 'hands-on-a-box')
      assert_equal(Gtn::Boxify::get_id('hands-on', 'a box', 'index.md'), 'hands-on-a-box-1')
      assert_equal(Gtn::Boxify::get_id('hands-on', 'a box', 'index.md'), 'hands-on-a-box-2')
      assert_equal(Gtn::Boxify::get_id('hands-on', 'a box', 'index2.md'), 'hands-on-a-box')
      assert_equal(Gtn::Boxify::get_id('hands-on', 'a box', 'index2.md'), 'hands-on-a-box-1')
      assert_equal(Gtn::Boxify::get_id('hands-on', 'a box', 'index2.md'), 'hands-on-a-box-2')

      assert_equal(Gtn::Boxify::get_id('hands-on', 'z-w-e', 'index2.md'), 'hands-on-z-w-e')
      assert_equal(Gtn::Boxify::get_id('hands-on', 'z-w-e', 'index2.md'), 'hands-on-z-w-e-1')
      assert_equal(Gtn::Boxify::get_id('hands-on', 'z-w-e', 'index2.md'), 'hands-on-z-w-e-2')
    end
  end
end
