# frozen_string_literal: true

require 'jekyll'

module Gtn
  # Generate boxes
  module Boxify
    @@ICONS = {
      'agenda' => '',
      'code-in' => 'far fa-keyboard',
      'code-out' => 'fas fa-laptop-code',
      'comment' => 'far fa-comment-dots',
      'details' => 'fas fa-info-circle',
      'feedback' => 'far fa-comments',
      'hands-on' => 'fas fa-pencil-alt',
      'hands_on' => 'fas fa-pencil-alt',
      'hidden' => '',
      'matrix' => '',
      'overview' => '',
      'question' => 'far fa-question-circle',
      'quote' => '',
      'solution' => 'far fa-eye',
      'spoken' => '',
      'tip' => 'far fa-lightbulb',
      'warning' => 'fas fa-exclamation-triangle',
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
      'en' => {
        'agenda' => 'Agenda',
        'code-in' => 'Input',
        'code-out' => 'Output',
        'comment' => 'Comment',
        'details' => 'Details',
        'hands-on' => 'Hands-on',
        'hands_on' => 'Hands-on',
        'question' => 'Question',
        'solution' => 'Solution',
        'tip' => 'Tip',
        'warning' => 'Warning',
      },
      'es' => {
        # Just google translated these.
        'agenda' => 'Agenda',
        'code-in' => 'Entrada',
        'code-out' => 'Salida',
        'comment' => 'Comentario',
        'details' => 'Detalles',
        'solution' => 'Solución',
        'warning' => '¡Precaucion!',

        # The only ones we have translations for??
        'hands-on' => 'Práctica',
        'hands_on' => 'Práctica',
        'question' => 'Preguntas',
        'tip' => 'Consejo',
      },
      'fr' => {
        'agenda' => 'Agenda',
        'code-in' => 'Entrée',
        'code-out' => 'Sortie',
        'comment' => 'Commentaire',
        'details' => 'Détails',
        'hands-on' => 'En pratique',
        'hands_on' => 'En pratique',
        'question' => 'Question',
        'solution' => 'Solution',
        'tip' => 'Astuce',
        'warning' => 'Attention',
      }
    }

    @title_unique_offsets = {}

    @@COLLAPSIBLE_BOXES = %w[
      details solution tip
    ]

    @@BOX_CLASSES = @@ICONS.keys.join '|'
    @@TITLE_CLASSES = @@ICONS.keys.map { |x| "#{x}-title" }.join '|'

    def self.box_classes
      @@BOX_CLASSES
    end

    def self.title_classes
      @@TITLE_CLASSES
    end

    def self.get_icon(icon_cls, emoji: false, a11y: false)
      return @@ICONS_EMOJI.fetch(icon_cls, '') if emoji

      icon = @@ICONS[icon_cls]

      # We support announcing the proper label of the box, e.g. 'hands on box',
      # but default to hiding this, as the icons are *mostly* decorative.
      icon_a11y_title = icon_cls.gsub(/[-_]/, ' ')
      icon_aria_label = a11y ? "title=\"#{icon_a11y_title} box\"" : ''
      accessible_addition = a11y ? %(<span class="sr-only">#{icon_a11y_title} box</span>) : ''

      if icon.nil?
        %(<span class="visually-hidden"></span>)
      elsif icon.start_with?('fa')
        %(<i class="#{icon}" aria-hidden="true" #{icon_aria_label}></i>#{accessible_addition})
      elsif icon.start_with?('ai')
        %(<i class="ai #{icon}" aria-hidden="true" #{icon_aria_label}></i>#{accessible_addition})
      end
    end

    def self.get_id(box_type, title, path)
      box_id = "#{box_type}-#{title}"
      box_safe = Jekyll::Utils.slugify(box_id)
      key = "#{path}|#{box_type}|#{box_safe}"

      if @title_unique_offsets.key?(key)
        box_safe_final = "#{box_safe}-#{@title_unique_offsets[key]}"
      else
        @title_unique_offsets[key] = 0
        box_safe_final = box_safe
      end
      @title_unique_offsets[key] += 1

      box_safe_final
    end

    def self.format_box_title(title, box_type, lang = 'en', noprefix: false)
      lang = 'en' if (lang == '') || lang.nil?
      title_fmted = (!title.nil? && title.length.positive? ? ": #{title}" : '')
      if noprefix && !title.nil?
        title
      else
        "#{@@BOX_TITLES[lang][box_type]}#{title_fmted}"
      end
    end

    def self.generate_collapsible_title(box_type, title, lang = 'en', key, contents: false, noprefix: false)
      box_id = get_id(box_type, title, key)
      box_title = format_box_title(title, box_type, lang, noprefix: noprefix)
      refers_to_contents = contents ? '-contents' : ''
      # These are all collapsed by default, details, tip, and solution.
      # rubocop:disable Layout/LineLength
      [box_id, %(
        <div class="box-title #{box_type}-title" id="#{box_id}">
        <button class="gtn-boxify-button #{box_type}" type="button" aria-controls="#{box_id}#{refers_to_contents}" aria-expanded="true">
          #{get_icon(box_type)} <span>#{box_title}</span>
          <span class="fold-unfold fa fa-minus-square"></span>
        </button>
        </div>
      ).split(/\n/).map(&:strip).join.strip]
      # rubocop:enable Layout/LineLength
    end

    def self.generate_static_title(box_type, title, lang = 'en', key, noprefix: false)
      box_id = get_id(box_type, title, key)
      box_title = format_box_title(title, box_type, lang, noprefix: noprefix)

      puts "Static | typ=#{box_type} | t=#{title} | l=#{lang} | k=#{key}" if title.nil?

      [box_id, %(
        <div class="box-title #{box_type}-title" id="#{box_id}">
          #{get_icon(box_type)} #{box_title}
        </div>
      ).split(/\n/).map(&:strip).join.strip]
    end

    def self.safe_title(title)
      title = '' if title.nil?

      title.gsub(/"/, '&quot;')
    end

    def self.generate_title(box_type, title, lang = 'en', key, contents: false, noprefix: false)
      title = safe_title(title)
      if @@COLLAPSIBLE_BOXES.include?(box_type)
        generate_collapsible_title(box_type, title, lang, key, contents: contents, noprefix: noprefix)
      else
        generate_static_title(box_type, title, lang, key, noprefix: noprefix)
      end
    end

    def self.generate_box(box_type, title, lang = 'en', key)
      title = safe_title(title)
      box_id, box_title = generate_title(box_type, title, lang, key, contents: true)
      %(
        <div class="box #{box_type}" markdown=0>
        #{box_title}
        <div class="box-content" id="#{box_id}-contents" markdown=1>
      ).split(/\n/).map(&:strip).join.strip
    end

    def self.replace_elements(text, lang = 'en', key)
      # We want to replace any `<x-title>(.*)</x-title>` bits
      # And replace them one by one with "proper" boxes, based on generate_title.
      #
      # We're going to rely on never having two on one line
      text.split("\n").map do |line|
        line.gsub(%r{<(?<type>[a-z-]*)-title>(?<title>.*?)</[a-z-]*-title>}) do |_m|
          title = Regexp.last_match[:title]
          type = Regexp.last_match[:type]
          _, box = generate_title(type, title, lang, key)
          box
        end
      end.join("\n")
    end
  end
end

if $PROGRAM_NAME == __FILE__
  require 'test/unit'
  # Test the box ID algorithm
  class BoxIdTest < Test::Unit::TestCase
    def test_single_page
      assert_equal(Gtn::Boxify.get_id('hands-on', 'a box', 'index.md'), 'hands-on-a-box')
      assert_equal(Gtn::Boxify.get_id('hands-on', 'a box', 'index.md'), 'hands-on-a-box-1')
      assert_equal(Gtn::Boxify.get_id('hands-on', 'a box', 'index.md'), 'hands-on-a-box-2')
      assert_equal(Gtn::Boxify.get_id('hands-on', 'a box', 'index2.md'), 'hands-on-a-box')
      assert_equal(Gtn::Boxify.get_id('hands-on', 'a box', 'index2.md'), 'hands-on-a-box-1')
      assert_equal(Gtn::Boxify.get_id('hands-on', 'a box', 'index2.md'), 'hands-on-a-box-2')

      assert_equal(Gtn::Boxify.get_id('hands-on', 'z-w-e', 'index2.md'), 'hands-on-z-w-e')
      assert_equal(Gtn::Boxify.get_id('hands-on', 'z-w-e', 'index2.md'), 'hands-on-z-w-e-1')
      assert_equal(Gtn::Boxify.get_id('hands-on', 'z-w-e', 'index2.md'), 'hands-on-z-w-e-2')
    end
  end
end
