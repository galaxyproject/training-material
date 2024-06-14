# frozen_string_literal: true

require './_plugins/jekyll-topic-filter'
require './_plugins/gtn'

# TODO: move into a lib somewhere.
def collapse_date_pretty(event)
  s = event['date_start']
  e = if event['date_end'].nil?
        s
      else
        event['date_end']
      end
  # want dates like "Mar 22-25, 2024" or "Mar 22-May 1, 2024"
  if s.year == e.year
    if s.month == e.month
      if s.day == e.day
        "#{s.strftime('%B')} #{s.day}, #{s.year}"
      else
        "#{s.strftime('%B')} #{s.day}-#{e.day}, #{s.year}"
      end
    else
      "#{s.strftime('%B')} #{s.day}-#{e.strftime('%B')} #{e.day}, #{s.year}"
    end
  else
    "#{s.strftime('%B')} #{s.day}, #{s.year}-#{e.strftime('%B')} #{e.day}, #{e.year}"
  end
end

def generate_topic_feeds(site)
  TopicFilter.list_topics(site).each do |topic|
    feed_path = File.join(site.dest, 'topics', topic, 'feed.xml')
    Jekyll.logger.info "[GTN/Feeds] Generating feed for #{topic}"

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

def generate_tag_topic_feeds(_site)
  # Any new materials in a topic with the equivalent tag
  # Any new materials tagged with that tag
  # Any news by tag
  ''
end

def all_date_sorted_materials(site)
  events = site.pages.select { |x| x['layout'] == 'event' || x['layout'] == 'event-external' }
  materials = TopicFilter.list_all_materials(site).reject { |k, _v| k['draft'] }
  news = site.posts.select { |x| x['layout'] == 'news' }

  bucket = events.map do |e|
    [Gtn::PublicationTimes.obtain_time(e.path).to_datetime, 'events', e, ['event'] + e.data.fetch('tags', [])]
  end

  materials.each do |m|
    tags = [m['topic_name']] + (m['tags'] || [])
    bucket += m.fetch('ref_tutorials', []).map do |t|
      [Gtn::PublicationTimes.obtain_time(t.path).to_datetime, 'tutorials', t, tags]
    end

    bucket += m.fetch('ref_tutorials', []).reject { |s| s.url =~ /-plain.html/ }.map do |s|
      [Gtn::PublicationTimes.obtain_time(s.path).to_datetime, 'slides', s, tags]
    end
  end

  bucket += news.map do |n|
    [n.date.to_datetime, 'news', n, ['news'] + n.data.fetch('tags', [])]
  end

  bucket += site.data['contributors'].map do |k, v|
    [DateTime.parse("#{v['joined']}-01"), 'contributors', k, ['contributor']]
  end
  bucket += site.data['funders'].map do |k, v|
    # TODO: backdate funders, organisations
    if v['joined']
      [DateTime.parse("#{v['joined']}-01"), 'funders', k, ['funder']]
    end
  end.compact
  bucket += site.data['organisations'].map do |k, v|
    if v['joined']
      [DateTime.parse("#{v['joined']}-01"), 'organisations', k, ['organisation']]
    end
  end.compact

  bucket.sort_by { |x| x[0] }.reverse
end

def group_bucket_by(bucket, group_by: 'day')
  case group_by
  when 'day'
    bucket
      .group_by { |x| x[0].strftime('%Y-%m-%d') }
      .to_h { |_k, v| [v.map { |x| x[0] }.min, v] }
  when 'week'
    bucket
      .group_by { |x| x[0].strftime('%Y-%W') }
      .to_h { |_k, v| [v.map { |x| x[0] }.min, v] }
  when 'month'
    bucket
      .group_by { |x| x[0].strftime('%Y-%m') }
      .to_h { |_k, v| [v.map { |x| x[0] }.min, v] }
  else
    raise "Unknown group_by: #{group_by}"
  end
end

def format_contents(xml, site, parts, title, group_by: 'day')
  # output += '<div xmlns="http://www.w3.org/1999/xhtml">'
end

# Our old style matrix bot postsx
def generate_matrix_feed(site, mats, group_by: 'day', filter_by: nil)
  # new materials (tut + sli)
  # new funders/contributors/orgs
  # new news posts(?)
  filter_title = nil
  if !filter_by.nil?
    mats = mats.select { |x| x[3].include?(filter_by) }
    filter_title = filter_by.gsub('-', ' ').capitalize
  end

  case group_by
  when 'day'
    # Reject anything that is today
    mats = mats.reject { |x| x[0].strftime('%Y-%m-%d') == Date.today.strftime('%Y-%m-%d') }
  when 'week'
    mats = mats.reject { |x| x[0].strftime('%Y-%W') == Date.today.strftime('%Y-%W') }
  when 'month'
    mats = mats.reject { |x| x[0].strftime('%Y-%m') == Date.today.strftime('%Y-%m') }
  end

  bucket = group_bucket_by(mats, group_by: group_by)
  lookup = {
    'day' => 'Daily',
    'week' => 'Weekly',
    'month' => 'Monthly'
  }

  path = "feeds/matrix-#{group_by}.xml"
  if filter_by
    path = "feeds/#{filter_by}-#{group_by}.xml"
  end

  feed_path = File.join(site.dest, path)
  Jekyll.logger.info '[GTN/Feeds] Generating matrix feed'

  dir = File.dirname(feed_path)
  FileUtils.mkdir_p(dir) unless File.directory?(dir)

  # Group by days

  builder = Nokogiri::XML::Builder.new(encoding: 'UTF-8') do |xml|
    # Set stylesheet
    xml.feed(xmlns: 'http://www.w3.org/2005/Atom') do
      # Set generator also needs a URI attribute
      xml.generator('Jekyll', uri: 'https://jekyllrb.com/')
      xml.link(href: "#{site.config['url']}#{site.baseurl}/#{path}", rel: 'self')
      # convert '2024-01-01' to date
      xml.updated(DateTime.now.rfc3339)
      xml.id("#{site.config['url']}#{site.baseurl}/#{path}")
      title_parts = ['GTN', filter_title, "#{lookup[group_by]} Updates"].compact
      xml.title(title_parts.join(' — '))
      xml.subtitle('An RSS feed with the latest materials, events, news in the GTN.')

      bucket.each do |date, parts|
        xml.entry do
          case group_by
          when 'day'
            title = "GTN #{lookup[group_by]} updates for #{date.strftime('%B %d, %Y')}"
          when 'week'
            title = "GTN #{lookup[group_by]} updates for #{date.strftime('W%W, %Y')}"
          when 'month'
            title = "GTN #{lookup[group_by]} updates for #{date.strftime('%B %Y')}"
          end
          xml.title(title)
          # Our IDs should be stable
          xml.id("#{site.config['url']}#{site.baseurl}/#{group_by}/#{date.strftime('%Y-%m-%d')}")

          # This is a feed of only NEW tutorials, so we only include publication times.
          xml.published(parts.map { |x| x[0] }.min.to_datetime.rfc3339)
          xml.updated(parts.map { |x| x[0] }.max.to_datetime.rfc3339)

          # xml.category(term: "new #{type}")
          xml.content(type: 'xhtml') do
            xml.div(xmlns: 'http://www.w3.org/1999/xhtml') do
              # xml.h4 title

              icon_for = {
                'contributors' => '🧑‍🏫',
                'funders' => '💰',
                'organisations' => '🏢',
                'events' => '📅',
                'tutorials' => '📚',
                'slides' => '🖼️',
                'news' => '📰',
              }

              prio = {
                'news' => 0,
                'events' => 1,
                'contributors' => 4,
                'funders' => 5,
                'organisations' => 6,
                'tutorials' => 2,
                'slides' => 3,
              }

              parts.group_by { |x| x[1] }.sort_by { |x| prio[x[1]] }.each do |type, items|
                xml.h4 "#{icon_for[type]} #{type.capitalize}"
                if items.length.positive?
                  xml.ul do
                    items.each do |date, _type, page, _tags|
                      xml.li do
                        if page.is_a?(String)
                          href = "#{site.config['url']}#{site.config['baseurl']}/hall-of-fame/#{page}/"
                          text = "@#{page}"
                        else
                          text = page.data['title']
                          href = "#{site.config['url']}#{site.config['baseurl']}#{page.url}"
                        end
                        if group_by != 'day'
                          text += " (#{date.strftime('%B %d, %Y')})"
                        end

                        xml.a(text, href: href)
                      end
                    end
                  end
                end
              end

              if group_by != 'day'
                xml.small do
                  xml.span 'Powered by '
                  xml.a('GTN RSS Feeds', href: 'https://training.galaxyproject.org/training-material/news/2024/06/04/gtn-standards-rss.html')
                end
              end
            end
          end

          xml.author do
            xml.name('GTN')
            xml.uri("#{site.config['url']}#{site.baseurl}/hall-of-fame/")
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
    %(type="text/xml" href="#{site.config['url']}#{site.baseurl}/feed-html.xslt.xml")
  )
  finalised.root.add_previous_sibling pi
  File.write(feed_path, finalised.to_xml)
end

def generate_event_feeds(site)
  events = site.pages.select { |x| x['layout'] == 'event' || x['layout'] == 'event-external' }
  feed_path = File.join(site.dest, 'events', 'feed.xml')
  Jekyll.logger.info '[GTN/Feeds] Generating event feed'

  # Pre-filering.
  updated = events.map { |x| Gtn::PublicationTimes.obtain_time(x.path) }.max

  events = events
           .reject { |x| x.data.fetch('draft', '').to_s == 'true' }
           .reject { |x| x.data['event_over'] == true } # Remove past events, prunes our feed nicely.
           .sort_by { |page| Gtn::PublicationTimes.obtain_time(page.path) }
           .reverse

  if !events.empty?
    Jekyll.logger.debug "Found #{events.length} events"
  end

  builder = Nokogiri::XML::Builder.new(encoding: 'UTF-8') do |xml|
    # Set stylesheet
    xml.feed(xmlns: 'http://www.w3.org/2005/Atom') do
      # Set generator also needs a URI attribute
      xml.generator('Jekyll', uri: 'https://jekyllrb.com/')
      xml.link(href: "#{site.config['url']}#{site.baseurl}/events/feed.xml", rel: 'self')
      xml.updated(updated.to_datetime.rfc3339)
      xml.id("#{site.config['url']}#{site.baseurl}/events/feed.xml")
      xml.title('Galaxy Training Network - Events')
      xml.subtitle('Events in the Inter-Galactic Network')

      events.each do |page|
        xml.entry do
          pdate = collapse_date_pretty(page.data)
          xml.title("[#{pdate}] #{page.data['title']}")
          link = "#{site.config['url']}#{site.baseurl}#{page.url}"
          xml.link(href: link)
          # Our links are stable
          xml.id(link)

          # This is a feed of only NEW tutorials, so we only include publication times.
          # xml.published(Gtn::PublicationTimes.obtain_time(page.path).to_datetime.rfc3339)
          xml.published(Gtn::PublicationTimes.obtain_time(page.path).to_datetime.rfc3339)
          xml.updated(Gtn::PublicationTimes.obtain_time(page.path).to_datetime.rfc3339)

          # TODO: find a better solution maybe with namespaces?
          # xml.category(term: "starts:#{page.data['date_start'].to_datetime.rfc3339}")
          # xml.category(term: "ends:#{(page.data['date_end'] || page.data['date_start']).to_datetime.rfc3339}")
          # xml.category(term: "days:#{page.data['duration']}")

          # xml.path(page.path)
          xml.category(term: "new #{page['layout']}")
          # xml.content(page.content, type: "html")
          xml.summary(page.data['description'])

          if page.data['location'] && page.data['location']['geo']
            lat = page.data['location']['geo']['lat']
            lon = page.data['location']['geo']['lon']
            xml.georss('point', "#{lat} #{lon}")
          end

          Gtn::Contributors.get_organisers(page.data).each do |c|
            xml.author do
              xml.name(Gtn::Contributors.fetch_name(site, c))
              xml.uri("#{site.config['url']}#{site.baseurl}/hall-of-fame/#{c}/")
              if page.data['contact_email']
                xml.email(page.data['contact_email'])
              end
            end
          end

          Gtn::Contributors.get_instructors(page.data).each do |c|
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

# Basically like `PageWithoutAFile`
Jekyll::Hooks.register :site, :post_write do |site|
  if Jekyll.env == 'production'
    generate_topic_feeds(site)
    generate_event_feeds(site)

    bucket = all_date_sorted_materials(site)
    bucket.freeze
    generate_matrix_feed(site, bucket, group_by: 'day')
    generate_matrix_feed(site, bucket, group_by: 'week')
    generate_matrix_feed(site, bucket, group_by: 'month')
    generate_matrix_feed(site, bucket, group_by: 'month', filter_by: 'single-cell')
    generate_matrix_feed(site, bucket, group_by: 'month', filter_by: 'genome-annotation')
    generate_matrix_feed(site, bucket, group_by: 'month', filter_by: 'one-health')
  end
end
