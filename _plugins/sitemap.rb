require 'date'

module Jekyll
  class SitemapGenerator < Generator
    safe true

    def generate(site)
      puts "[GTN/Sitemap] Generating"
      result = '<?xml version="1.0" encoding="UTF-8"?>'
      result += '<urlset xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.sitemaps.org/schemas/sitemap/0.9 http://www.sitemaps.org/schemas/sitemap/0.9/sitemap.xsd" xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'

      site.pages.select{|t| t.path !~ /ipynb$/ }.each {|t|
        begin
          d = t.data['last_modified_at']
          d.format = '%FT%T%:z'
          formatted_date = d.to_s
        rescue
          d = Time.new
          formatted_date = d.strftime('%FT%T%:z')
        end

        result += "<url><loc>#{site.config['url'] + site.config['baseurl'] + t.url}</loc><lastmod>#{formatted_date}</lastmod></url>"
      }
      result += "</urlset>"

      page2 = PageWithoutAFile.new(site, "", '.', "sitemap.xml")
      page2.content = result
      site.pages << page2
    end
  end
end
