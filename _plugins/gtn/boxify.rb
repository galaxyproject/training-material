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

    @@BOX_TITLES = {
      "en" => {
        "agenda" => "Agenda",
        "code-in" => "Input",
        "code-out" => "Output",
        "comment" => "Comment",
        "details" => "Details",
        "hands-on" => "Hands-on",
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
        "question" => "Preguntas",
        "tip" => "Consejo",
      }
    }

    @@COLLAPSIBLE_BOXES = [
      "details", "solution", "tip",
    ]

    @@BOX_CLASSES = @@ICONS.keys.join "|"
    @@TITLE_CLASSES = @@ICONS.keys.map{|x| "#{x}-title" }.join "|"

    def self.box_classes
      @@BOX_CLASSES
    end

    def self.get_icon(icon)
       if icon.start_with?("fa")
        %Q(<i class="#{icon}" aria-hidden="true"></i><span class="visually-hidden">#{@text}</span>)
       elsif icon.start_with?("ai")
        %Q(<i class="ai #{icon}" aria-hidden="true"></i><span class="visually-hidden">#{@text}</span>)
       end
    end

    def self.generate_collapsible_title(box_type, count, title, lang="en")
      title_fmted = (title.length > 0 ? ": #{title}" : "")
      return %Q(
        <div class="box-title">
        <button type="button" aria-controls="box-#{box_type}-#{count}" aria-expanded="true" aria-label="Toggle #{box_type} box: #{title}">
          #{self.get_icon(@@ICONS[box_type])} #{@@BOX_TITLES[lang][box_type]}#{title_fmted}
          <span role="button" class="fold-unfold fa fa-plus-square"></span>
        </button>
        </div>
      ).split(/\n/).map{|x| x.lstrip.rstrip}.join("").lstrip.rstrip
    end

    def self.generate_static_title(box_type, count, title, lang="en")
      title_fmted = (title.length > 0 ? ": #{title}" : "")
      return %Q(
        <div class="box-title" aria-label="#{box_type} box: #{title}">
          #{self.get_icon(@@ICONS[box_type])} #{@@BOX_TITLES[lang][box_type]}#{title_fmted}
        </div>
      ).split(/\n/).map{|x| x.lstrip.rstrip}.join("").lstrip.rstrip
    end

    def self.generate_title(box_type, count, title, lang="en")
      if @@COLLAPSIBLE_BOXES.include?(box_type)
        self.generate_collapsible_title(box_type, count, title, lang)
      else
        self.generate_static_title(box_type, count, title, lang)
      end
    end

    def self.generate_box(box_type, count, title, lang="en")
      box_title = generate_title(box_type, count, title, lang)
      return %Q(
        <div class="box #{box_type}" markdown=0>
        #{box_title}
        <div id="box-#{box_type}-#{count}" class="box-content" markdown=1>
      ).split(/\n/).map{|x| x.lstrip.rstrip}.join("").lstrip.rstrip
    end


  end
end
