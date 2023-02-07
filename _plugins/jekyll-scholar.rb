require './_plugins/gtn/scholar'

module Jekyll
  class CiteTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(context)
      page = context.registers[:page]
      site = context.registers[:site]
      Gtn::Scholar.load_bib(site)

      # Mark this page as having citations
      page['cited'] = true

      if page['citation_target'] == 'R'
        return "@#{@text}"
      end

      # Which page is rendering this tag?
      source_page = page['path']

      # Citation Frequency
      if ! site.config.has_key?("citation_count")
        site.config["citation_count"] = Hash.new(0)
      end
      site.config["citation_count"][@text] += 1


      # If the overall cache is nil, create it
      if site.config['citation_cache'].nil?
        site.config['citation_cache'] = Hash.new
      end
      # If the individual page in the chace is nil, create it.
      if site.config['citation_cache'][source_page].nil?
        site.config['citation_cache'][source_page] = Array.new
      end

       # Push it to our cache.
      site.config['citation_cache'][source_page].push(@text)

      begin
        citation_text = site.config['cached_citeproc'].render(:citation, id: @text)
        layout = page.fetch('layout', nil)
        if ['tutorial_slides', 'base_slides', 'introduction_slides'].include? layout
          doi = site.config['cached_citeproc'].items[@text].doi
          url = site.config['cached_citeproc'].items[@text].url
          if ! doi.nil?
            furl = "https://doi.org/#{doi}"
          elsif ! url.nil?
            furl = url
          else
            furl = nil
          end

          if furl.nil?
            res = %Q(<span class="citation">#{citation_text}</span>)
          else
            res = %Q(<span class="citation"><a href="#{furl}">#{citation_text}</a></span>)
          end
        else
          res = %Q(<span class="citation"><a href="##{@text}">#{citation_text}</a></span>)
        end

      rescue => error
        puts "[GTN/scholar] Could not render #{@text} from #{source_page} (#{error})"
        res = %Q(<span>ERROR INVALID CITATION #{@text}</span>)
      end

      if page['citation_target'] == 'jupyter'
        res.gsub!(/"/, '\"')
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
      site = context.registers[:site]
      Gtn::Scholar.load_bib(site)
      # Which page is rendering this tag?
      source_page = context.registers[:page]['path']
      global_bib = site.config['cached_global_bib']
      citeproc = site.config['cached_citeproc']
      # We have our page's citations
      citations = site.config['citation_cache'][source_page] || []
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
        r = Gtn::Scholar.render_citation(c)
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
