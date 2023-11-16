# frozen_string_literal: true

module Gtn
  # GTN implementation of Jekyll::Scholar except faster.
  module Scholar
    def self.load_bib(site)
      return if site.config.key?('cached_global_bib')

      (global_bib, cp) = populate_cache
      site.config['cached_global_bib'] = global_bib
      site.config['cached_citeproc'] = cp
    end

    def self.populate_cache
      @@cache ||= discover_bib
    end

    def self.render_citation(key)
      (global_bib, citeproc) = populate_cache

      text = citeproc.render(:bibliography, id: key)[0]
      entry = global_bib[key]
      text += " #{entry.note}." if entry.note
      doi = entry.fetch('doi', nil)
      text += " <a href=\"https://doi.org/#{doi}\">#{doi}</a>" if doi
      url = entry.fetch('url', nil)
      text += " <a href=\"#{url}\">#{url}</a>" if url && !(url.index('doi.org') && entry.doi)

      text
    end

    def self.discover_bib
      puts '[GTN/scholar] Creating global bib cache'
      global_bib = BibTeX::Bibliography.new
      bib_paths = [Find.find('./topics'), Find.find('./faqs'), Find.find('./news')].lazy.flat_map(&:lazy).grep(/bib$/)
      bib_paths.each do |path|
        BibTeX.open(path).each do |x|
          x = x.convert_latex
          global_bib << x
        end
      end
      puts "[GTN/scholar] Done: #{global_bib.length}"
      style = CSL::Style.load('_layouts/g3.csl')
      cp = CiteProc::Processor.new style: style,
                                   format: 'html', locale: 'en'
      cp.import global_bib.to_citeproc

      [global_bib, cp]
    end
  end
end
