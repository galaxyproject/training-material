module Jekyll
  class CiteTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(context)
      # Which page is rendering this tag?
      source_page = context.registers[:page]['path']
      # If the overall cache is nil, create it
      if context.registers[:site].config['citation_cache'].nil?
        context.registers[:site].config['citation_cache'] = Hash.new
      end
      # If the individual page in the chace is nil, create it.
      if context.registers[:site].config['citation_cache'][source_page].nil?
        context.registers[:site].config['citation_cache'][source_page] = Array.new
      end

      # Push it to our cache.
      context.registers[:site].config['citation_cache'][source_page].push(@text)
      # Mark this page as having citations
      context.registers[:page]['cited'] = true

      begin
        citation_text = context.registers[:site].config['cached_citeproc'].render(:citation, id: @text)
        res = %Q(<span class="citation"><a href="##{@text}">#{citation_text}</a></span>)
      rescue
        puts "[GTN/scholar] Could not render #{@text} from #{source_page}"
        res = %Q(<span>ERROR INVALID CITATION #{@text}</span>)
      end

      res
    end

    private
  end

  class BibTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(context)
      # Which page is rendering this tag?
      source_page = context.registers[:page]['path']
      global_bib = context.registers[:site].config['cached_global_bib']
      citeproc = context.registers[:site].config['cached_citeproc']
      # We have our page's citations
      citations = context.registers[:site].config['citation_cache'][source_page]
      # For each of these citation IDs, we need to get the formatted version + pull out
      # year, month for sorting.
      unique_citations = citations.reduce(Hash.new(0)) { |a, b| a[b] += 1; a }.keys
      # Remove nil citations
      unique_citations = unique_citations.select{|c| ! global_bib[c].nil? }
      # And now sort them by date + names
      sorted_citations = unique_citations.sort{|a, b|
        global_bib[a].date.to_s + global_bib[a].names.join(' ') <=>
        global_bib[b].date.to_s + global_bib[b].names.join(' ')
      }

      out = "<ol class=\"bibliography\">"
      out += sorted_citations.map{|c|
        r = citeproc.render(:bibliography, id: c)[0]
        entry = global_bib[c]
        if entry.note
          r += " #{entry.note}."
        end
        doi = entry.fetch('doi', nil)
        if doi
          r += " <a href=\"https://doi.org/#{doi}\">#{doi}</a>"
        end
        url = entry.fetch('url', nil)
        if url
          if ! (url.index('doi.org') and entry.doi)
            r += " <a href=\"#{url}\">#{url}</a>"
          end
        end
        %Q(<li id="#{c}">#{r}</li>)
      }.join("\n")
      out += "</ol>"
      out
    end

    private
  end

end

Liquid::Template.register_tag('cite', Jekyll::CiteTag)
Liquid::Template.register_tag('bibliography', Jekyll::BibTag)
