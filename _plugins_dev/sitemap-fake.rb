module Jekyll
  class SitemapGenerator < Generator
    safe true

    def generate(site)
      result = '<?xml version="1.0" encoding="UTF-8"?>'
      result += '<!-- This sitemap has FAKE modification dates intentionally, that step is slow. -->'
      result += '<urlset xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.sitemaps.org/schemas/sitemap/0.9 http://www.sitemaps.org/schemas/sitemap/0.9/sitemap.xsd" xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'

      site.pages.each {|t|
        result += "<url><loc>#{site.config['url'] + site.config['baseurl'] + t.url}</loc><lastmod>2016-06-30T18:00:00-07:00</lastmod></url>"
      }
      result += "</urlset>"

      page2 = PageWithoutAFile.new(site, "", '.', "sitemap.xml")
      page2.content = result
      site.pages << page2
    end
  end
end

