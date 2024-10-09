require './_plugins/gtn'

Jekyll::Hooks.register :site, :post_write do |site|
  # No need to run this except in prod.
  if Jekyll.env != 'production'
    # Build our Images!
    count = 0
    social_start = Time.now
    templates = {}
    meta = {}

    # Pre-existing templates
    Dir.glob('assets/branding/social/*.svg').each do |file|
      templates[File.basename(file, '.svg')] = File.read(file)
    end

    Jekyll.logger.info "[GTN/Socials] Writing social media cards, #{templates.keys} discovered"
    # svgs = []

    (site.pages + site.posts).each do |page|
      if page.data['layout'].nil? || page.data['layout'] == ''
        next
      end

      if ['faq-page', 'slides-plain', 'workflow-list', 'base', 'event-track', 'embed', 'none', 'page'].include? page.data['layout']
        next
      end

      title_override = nil
      parts = page.path.split('/')
      if parts.include?('topics')
        topic_id = parts[parts.index('topics') + 1]
        topic_name = site.data[topic_id]['title'].upcase
      else
        topic_id = nil
        topic_name = ""
      end
      mod_time = Gtn::ModificationTimes.obtain_time(page.path)
      # Default is usually fine
      authors = Gtn::Contributors.get_authors(page.data).map { |c| Gtn::Contributors.fetch_name(site, c) }

      # Handle specific layouts.
      case page.data['layout']
      when 'event'
        tpl = templates['event'].dup
        tpl = tpl.gsub('TOPIC_NAME', 'EVENT')
        # Event specific
        authors = Gtn::Contributors.get_organisers(page.data).map { |c| Gtn::Contributors.fetch_name(site, c) }
        if page.data['mode'] == 'online'
          tpl = tpl.gsub('In-Person', 'Online')
        end
        if page.data['async'] == false
          tpl = tpl.gsub('Asynchronous', 'Synchronous')
        end
      when 'workflow'
        tpl = templates['workflow'].dup
        tpl = tpl.gsub('TOPIC_NAME', site.data[page.data['workflow']['topic_id']]['title'].upcase)
        authors = page.data['workflow']['creators'].map { |c| c['name'] }

        mod_time = page.data['workflow']['modified']
      when 'faq'
        tpl = templates['faq'].dup
        if page.url =~ /faqs\/gtn/
          tpl = tpl.gsub('TOPIC_NAME', 'GTN FAQ')
        elsif page.url =~ /faqs\/galaxy/
          tpl = tpl.gsub('TOPIC_NAME', 'Galaxy FAQ')
        end
      when 'recordings'
        tpl = templates['recording'].dup
        authors = page.data['material']['ref'].data.fetch('recordings', [])
          .reverse.map { |r| r.fetch('speakers', []) + r.fetch('captioners',[]) }.flatten.uniq
        authors = authors.map { |c| Gtn::Contributors.fetch_name(site, c) }

        mod_time = [
          Gtn::ModificationTimes.obtain_time(page.path),
          Gtn::ModificationTimes.obtain_time(page.data['material']['ref'].path)
        ].max
      when 'tutorial_hands_on'
        tpl = templates['tutorial'].dup
      when 'tutorial_slides'
        tpl = templates['slides'].dup
      when 'introduction_slides'
        tpl = templates['slides'].dup
      when 'learning-pathway'
        tpl = templates['learning-pathway'].dup
        tpl = tpl.gsub('TOPIC_NAME', 'LP')
        authors = page.data['editorial_board'].map { |c| Gtn::Contributors.fetch_name(site, c) }
      when 'topic'
        tpl = templates['topic'].dup
        tpl = tpl.gsub('TOPIC_NAME', 'Topic')
        topic = site.data[page.data['topic_name']]
        title_override = topic['title']
        authors = topic.fetch('editorial_board', []).map { |c| Gtn::Contributors.fetch_name(site, c) }

        mod_time = [
          Gtn::ModificationTimes.obtain_time(page.path),
          Gtn::ModificationTimes.obtain_time(page.path.gsub('index.md', 'metadata.yaml')),
          Gtn::ModificationTimes.obtain_time("metadata/#{page.data['topic_name']}.yaml"),
        ].max
      else
        Jekyll.logger.debug "[GTN/Socials] Skipping #{page.data['layout']} => #{page.path}"
        next
      end

      if ! topic_name.nil?
        tpl = tpl.gsub('TOPIC_NAME', topic_name)
      end

      tpl = tpl.gsub(/gxy.io\/GTN:....../, "gxy.io/GTN:#{page.data['short_id']}")

      # Title
      tpl = tpl
        .gsub('TITLE2TITLE2', '')
        .gsub('TITLE3TITLE3', '')
        .gsub('TITLE4TITLE4', '')
        .gsub('TITLE1TITLE1', (title_override || page.data['title'].gsub('_', ' ')).gsub('&', '&amp;'))

      # Authors
      if authors.length > 0
        tpl = tpl.gsub('AUTHOR1', authors[0])
      else
        tpl = tpl.gsub('AUTHOR1', '')
      end
      if authors.length > 1
        tpl = tpl.gsub('AUTHOR2', authors[1])
      else
        tpl = tpl.gsub('AUTHOR2', '')
      end
      if authors.length > 2
        tpl = tpl.gsub('AUTHOR3', authors[2])
      else
        tpl = tpl.gsub('AUTHOR3', '')
      end
      if authors.length > 4
        tpl = tpl.gsub('AUTHOR4', 'et al.')
      elsif authors.length > 3
        tpl = tpl.gsub('AUTHOR4', authors[3])
      else
        tpl = tpl.gsub('AUTHOR4', '')
      end

      if page.url.end_with?('.html')
        svg_path = page.url.sub('.html', '.svg')
      else
        svg_path = File.join(page.url, 'index.svg')
      end
      svg = File.join(site.dest, svg_path)
      if ! File.directory?(File.dirname(svg))
        FileUtils.mkdir_p(File.dirname(svg))
      end
      File.open(svg, 'w') do |f|
        f.write(tpl)
      end
      # svgs << svg

      meta[svg_path] = mod_time.to_i

      # Convert to PNG
      count += 1
    end
    Jekyll.logger.info "[GTN/Socials] Social media cards written in #{Time.now - social_start} seconds for #{count} pages"

    File.open(File.join(site.dest, 'api', 'social-meta.json'), 'w') do |f|
      f.write(meta.to_json)
    end

    # This takes on the order of 30 minutes.
    # png_start = Time.now
    # svgs.each do |svg|
    #   # I'm not sure about this.
    #   `magick -density 50 #{svg} #{svg.sub('.svg', '.png')}`
    # end
    # Jekyll.logger.info "[GTN/Socials] Social media cards convert to png in #{Time.now - png_start} seconds"
  end
end
