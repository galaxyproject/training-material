# frozen_string_literal: true

require 'date'
require './_plugins/gtn'

module Jekyll
  # Generate a sitemap like Jekyll::Sitemap
  class SitemapGenerator < Generator
    safe true

    ##
    # Generate a sitemap.xml file
    # We reimplement the default Jekyll sitemap generator, because we want to
    # leverage the GTN::ModificationTimes class to obtain the last modification
    # date of a page, in a more efficient way than the default Jekyll sitemap
    #
    # Params:
    # +site+:: The +Jekyll::Site+ object
    def generate(site)
      puts '[GTN/Sitemap] Generating'
      result = '<?xml version="1.0" encoding="UTF-8"?>'
      result += '<urlset xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ' \
                'xsi:schemaLocation="http://www.sitemaps.org/schemas/sitemap/0.9 ' \
                'http://www.sitemaps.org/schemas/sitemap/0.9/sitemap.xsd" ' \
                'xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'

      site.pages.reject { |t| t.path =~ /ipynb$/ }.each do |t|
        begin
          d = Gtn::ModificationTimes.obtain_time(t.path)
          d.format = '%FT%T%:z'
          formatted_date = d.to_s
        rescue StandardError
          d = Time.new
          formatted_date = d.strftime('%FT%T%:z')
        end

        result += "<url><loc>#{site.config['url'] + site.config['baseurl'] + t.url}</loc>" \
                  "<lastmod>#{formatted_date}</lastmod></url>"
      end
      result += '</urlset>'

      page2 = PageWithoutAFile.new(site, '', '.', 'sitemap.xml')
      page2.content = result
      site.pages << page2
    end
  end
end
