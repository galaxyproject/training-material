# frozen_string_literal: true

require './_plugins/jekyll-topic-filter'
require './_plugins/gtn'

def generate_topic_feeds(site)
  TopicFilter.list_topics(site).each do |topic|
    feed_path = File.join(site.dest, 'topics', topic, 'feed.xml')
    Jekyll.logger.debug "Generating feed for #{topic} => #{feed_path}"

    topic_pages = site.pages
                      .select { |x| x.path =~ %r{^\.?/?topics/#{topic}} }
                      .select { |x| x.path =~ %r{(tutorial.md|slides.html|faqs/.*.md)} }
                      .reject { |x| x.path =~ /index.md/ }
                      .reject { |x| x.data.fetch('draft', '').to_s == 'true' }
                      .reject { |x| x.url =~ /slides-plain.html/ }
                      .reject { |x| File.symlink?(x.path) } # Remove symlinks to other faqs/tutorials
                      .uniq(&:path)
                      .sort_by { |page| Gtn::PublicationTimes.obtain_time(page.path) }
                      .reverse

    if topic_pages.empty?
      Jekyll.logger.warn "No pages for #{topic}"
      next
    else
      Jekyll.logger.debug "Found #{topic_pages.length} pages for #{topic}"
    end

    builder = Nokogiri::XML::Builder.new(encoding: 'UTF-8') do |xml|
      # Set stylesheet
      xml.feed(xmlns: 'http://www.w3.org/2005/Atom') do
        # Set generator also needs a URI attribute
        xml.generator('Jekyll', uri: 'https://jekyllrb.com/')
        xml.link(href: "#{site.config['url']}#{site.baseurl}/topics/#{topic}/feed.xml", rel: 'self')
        xml.updated(Gtn::ModificationTimes.obtain_time(topic_pages.first.path).to_datetime.rfc3339)
        xml.id("#{site.config['url']}#{site.baseurl}/topics/#{topic}/feed.xml")
        topic_title = site.data[topic]['title']
        xml.title("Galaxy Training Network - #{topic_title}")
        xml.subtitle("Recently added tutorials, slides, and FAQs in the #{topic} topic")

        topic_pages.each do |page|
          page_type = if page.path =~ %r{faqs/.*.md}
                        'faq'
                      else
                        page.path.split('/').last.split('.').first
                      end

          xml.entry do
            xml.title(page.data['title'])
            link = "#{site.config['url']}#{site.baseurl}#{page.url}"
            xml.link(href: link)
            # Our links are stable
            xml.id(link)

            # This is a feed of only NEW tutorials, so we only include publication times.
            # xml.published(Gtn::PublicationTimes.obtain_time(page.path).to_datetime.rfc3339)
            xml.updated(Gtn::PublicationTimes.obtain_time(page.path).to_datetime.rfc3339)

            # xml.path(page.path)
            xml.category(term: "new #{page_type}")
            # xml.content(page.content, type: "html")
            xml.summary(page.content.strip.split("\n").first, type: 'html')

            Gtn::Contributors.get_authors(page.data).each do |c|
              xml.author do
                xml.name(Gtn::Contributors.fetch_name(site, c))
                xml.uri("#{site.config['url']}#{site.baseurl}/hall-of-fame/#{c}/")
              end
            end

            Gtn::Contributors.get_non_authors(page.data).each do |c|
              xml.contributor do
                xml.name(Gtn::Contributors.fetch_name(site, c))
                xml.uri("#{site.config['url']}#{site.baseurl}/hall-of-fame/#{c}/")
              end
            end
          end
        end
      end
    end

    # The builder won't let you add a processing instruction, so we have to
    # serialise it to a string and then parse it again. Ridiculous.
    finalised = Nokogiri::XML builder.to_xml
    pi = Nokogiri::XML::ProcessingInstruction.new(
      finalised, 'xml-stylesheet',
      %(type="text/xml" href="#{site.config['url']}#{site.baseurl}/feed.xslt.xml")
    )
    finalised.root.add_previous_sibling pi
    File.write(feed_path, finalised.to_xml)
  end

  nil
end

# Basically like `PageWithoutAFile`
Jekyll::Hooks.register :site, :post_write do |site|
  if Jekyll.env == 'production'
    generate_topic_feeds(site)
  end
end
