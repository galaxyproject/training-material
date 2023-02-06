module Gtn
  module Scholar
    def self.load_bib(site)
      if not site.config.has_key?("cached_global_bib")
        (global_bib, cp) = self.populate_cache()
        site.config['cached_global_bib'] = global_bib
        site.config['cached_citeproc'] = cp
      end
    end

    #def self.bib()
    #  @@cache ||= self.populate_cache()
    #end

    def self.populate_cache()
      puts "[GTN/scholar] Creating global bib cache"
      global_bib = BibTeX::Bibliography.new
      bib_paths = [Find.find('./topics'), Find.find('./faqs'), Find.find('./news')].lazy.flat_map(&:lazy).select{|x| x =~ /bib$/}
      bib_paths.each{|path|
        for x in BibTeX.open(path)
          x = x.convert_latex
          global_bib << x
        end
      }
      puts "[GTN/scholar] Done: #{global_bib.length}"
      style = CSL::Style.load("_layouts/g3.csl")
      cp = CiteProc::Processor.new style: style,
                                   format: 'html', locale: 'en'
      cp.import global_bib.to_citeproc

      [global_bib, cp]
    end
  end
end
