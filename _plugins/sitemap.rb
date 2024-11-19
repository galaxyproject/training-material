# frozen_string_literal: true

module Jekyll
  # Generate a sitemap like Jekyll::Sitemap
  class SitemapGenerator2 < Generator
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
      if Jekyll.env == 'production'
        _build(site)
      else
        Jekyll.logger.info '[GTN/Sitemap] Skipping in development mode'
      end
    end

    def _build(site)
      # We import later in case we don't need to bother importing in the first place.
      require 'date'
      require './_plugins/gtn'

      Jekyll.logger.info '[GTN/Sitemap] Generating'
      result = '<?xml version="1.0" encoding="UTF-8"?>'
      result += '<urlset xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ' \
                'xsi:schemaLocation="http://www.sitemaps.org/schemas/sitemap/0.9 ' \
                'http://www.sitemaps.org/schemas/sitemap/0.9/sitemap.xsd" ' \
                'xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'

      subset_pages = site.pages
        .reject { |t| t.path =~ /ipynb$/ || t.path =~ /api\/ga4gh\/trs\/v2/}
        .reject { |t| t.data.fetch('layout', 'page') =~ /external/}
        .reject { |t| t.data.fetch('hands_on', '') == 'external'}

      subset_pages.each do |t|
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
