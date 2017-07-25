class Vimeo < Liquid::Tag
  Syntax = /^\s*([^\s]+)(\s+(\d+)\s+(\d+)\s*)?/

  def initialize(tagName, markup, tokens)
    super

    if markup =~ Syntax then
      @id = $1

      if $2.nil? then
          @width = 560
          @height = 420
      else
          @width = $2.to_i
          @height = $3.to_i
      end
    else
      raise "No Vimeo ID provided in the \"youtube\" tag"
    end
  end

  def render(context)
    "<div class=\"video-container\"><iframe width=\"#{@width}\" height=\"#{@height}\" frameborder=\"0\" webkitallowfullscreen mozallowfullscreen allowfullscreen src=\"https://player.vimeo.com/video/#{@id}?color=2980b9&title=0&byline=0&portrait=0\"></iframe></div>"
  end

  Liquid::Template.register_tag "vimeo", self
end
