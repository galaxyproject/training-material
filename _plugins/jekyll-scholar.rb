# frozen_string_literal: true

require './_plugins/gtn/scholar'

module Jekyll
  # {% cite X %} which generates the link to the bib + text
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

      return "@#{@text}" if page['citation_target'] == 'R'

      # Which page is rendering this tag?
      source_page = page['path']

      # Citation Frequency
      site.config['citation_count'] = Hash.new(0) if !site.config.key?('citation_count')
      site.config['citation_count'][@text] += 1

      # If the overall cache is nil, create it
      site.config['citation_cache'] = {} if site.config['citation_cache'].nil?
      # If the individual page in the chace is nil, create it.
      site.config['citation_cache'][source_page] = [] if site.config['citation_cache'][source_page].nil?

      # Push it to our cache.
      site.config['citation_cache'][source_page].push(@text)

      begin
        citation_text = site.config['cached_citeproc'].render(:citation, id: @text)
        layout = page.fetch('layout', nil)
        if %w[tutorial_slides base_slides introduction_slides].include? layout
          doi = site.config['cached_citeproc'].items[@text].doi
          url = site.config['cached_citeproc'].items[@text].url
          furl = if !doi.nil?
                   "https://doi.org/#{doi}"
                 elsif !url.nil?
                   url
                 end

          res = if furl.nil?
                  %(<span class="citation">#{citation_text}</span>)
                else
                  %(<span class="citation"><a href="#{furl}">#{citation_text}</a></span>)
                end
        else
          res = %(<span class="citation"><a href="##{@text}">#{citation_text}</a></span>)
        end
      rescue StandardError => e
        puts "[GTN/scholar] Could not render #{@text} from #{source_page} (#{e})"
        res = %(<span>ERROR INVALID CITATION #{@text}</span>)
      end

      res.gsub!(/"/, '\"') if page['citation_target'] == 'jupyter'

      res
    end
  end

  # {% cite_url X %} which generates URL for the article
  class CiteUrlTag < Liquid::Tag
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

      return "@#{@text}" if page['citation_target'] == 'R'

      # Which page is rendering this tag?
      source_page = page['path']

      # Citation Frequency
      site.config['citation_count'] = Hash.new(0) if !site.config.key?('citation_count')
      site.config['citation_count'][@text] += 1

      # If the overall cache is nil, create it
      site.config['citation_cache'] = {} if site.config['citation_cache'].nil?
      # If the individual page in the chace is nil, create it.
      site.config['citation_cache'][source_page] = [] if site.config['citation_cache'][source_page].nil?

      # Push it to our cache.
      site.config['citation_cache'][source_page].push(@text)

      begin
        doi = site.config['cached_citeproc'].items[@text].doi
        url = site.config['cached_citeproc'].items[@text].url
        if !doi.nil?
          "https://doi.org/#{doi}"
        elsif !url.nil?
          url
        end
        res = url
      rescue StandardError => e
        puts "[GTN/scholar] Could not render #{@text} from #{source_page} (#{e})"
        res = %(<span>https://example.com/ERROR+INVALID+CITATION+#{@text}</span>)
      end
      res
    end
  end

  # {% bibliography %} which generates the bibliography
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
      # citeproc = site.config['cached_citeproc']
      # We have our page's citations
      citations = site.config['citation_cache'][source_page] || []
      # For each of these citation IDs, we need to get the formatted version + pull out
      # year, month for sorting.
      unique_citations = citations.each_with_object(Hash.new(0)) do |b, a|
        a[b] += 1
      end.keys
      # Remove nil citations
      unique_citations = unique_citations.reject { |c| global_bib[c].nil? }
      # And now sort them by date + names
      sorted_citations = unique_citations.sort do |a, b|
        global_bib[a].date.to_s + global_bib[a].names.join(' ') <=>
          global_bib[b].date.to_s + global_bib[b].names.join(' ')
      end

      out = '<ol class="bibliography">'
      out += sorted_citations.map do |c|
        r = Gtn::Scholar.render_citation(c)
        %(<li id="#{c}">#{r}</li>)
      end.join("\n")
      out += '</ol>'
      out
    end
  end
end

Liquid::Template.register_tag('cite', Jekyll::CiteTag)
Liquid::Template.register_tag('cite_url', Jekyll::CiteUrlTag)
Liquid::Template.register_tag('bibliography', Jekyll::BibTag)
