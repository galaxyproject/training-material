# frozen_string_literal: true

module Gtn
  # This module is responsible for generating shortlinks for tutorials and FAQs and any other pages we add.
  #
  # Every category gets its own prefix letter.
  module Shortlinks
    CATEGORY_TUTORIAL = 'T'
    CATEGORY_SLIDES = 'S'
    CATEGORY_FAQ = 'F'
    CATEGORY_NEWS = 'N'
    CATEGORY_PATHWAYS = 'P'
    CATEGORY_EVENTS = 'E'
    CATEGORY_WORKFLOW = 'W'

    REDIRECT_TEMPLATE = <<~REDIR
      <!DOCTYPE html>
      <html lang="en-US">
        <meta charset="utf-8">
        <title>Redirecting&hellip;</title>
        <link rel="canonical" href="REDIRECT_URL">
        <script>location="REDIRECT_URL"</script>
        <meta http-equiv="refresh" content="0; url=REDIRECT_URL">
        <meta name="robots" content="noindex">
        <h1>Redirecting&hellip;</h1>
        <a href="REDIRECT_URL">Click here if you are not redirected.</a>
      </html>
    REDIR

    def self.mapped?(tutorial, current_mapping)
      current_mapping['id'].values.include? tutorial
    end

    ##
    # Duplicate of the jekyll-redirect-from plugin template.
    # We can't use that for, reasons.
    def self.html_redirect(target)
      REDIRECT_TEMPLATE.gsub('REDIRECT_URL', target)
    end

    ##
    # Fix missing symlinks (usually exist because the target file has been
    # renamed and doesn't exist anymore.) However, a redirect *will* be present
    # for the original filename so we just fix the missing symlink.
    #
    # Params:
    # +site+:: The Jekyll site object
    def self.fix_missing_redirs(site)
      missing_redirs = site.data['shortlinks']['id'].select do |id, target|
        short_link = "short/#{id}.html"
        ! File.exist?(site.in_dest_dir(short_link))
      end

      missing_redirs.each do |id, target|
        short_link = "short/#{id}.html"
        Jekyll.logger.warn "[GTN/Shortlink]" "Shortlink target #{target} does not exist for shortlink #{short_link}, fixing."
        File.write(site.in_dest_dir(short_link), Gtn::Shortlinks.html_redirect(target))
      end
    end

    def self.update(current_mapping)
      current_mapping['id'] = {} if !current_mapping.key? 'id'

      current_mapping['name'] = {} if !current_mapping.key? 'name'

      # Discover tutorials
      num_tutorials = current_mapping['id'].select { |x| x[0] == CATEGORY_TUTORIAL }.length
      Dir.glob('topics/*/tutorials/*/tutorial.md').each do |tutorial|
        html_path = "/#{tutorial.gsub(/md$/, 'html')}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = num_tutorials.to_s.rjust(5, '0')
          short_code = CATEGORY_TUTORIAL + short_code_number
          puts "Discovered tutorial #{short_code} (#{html_path})"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
          num_tutorials += 1
        end

        # Also generate one from topic/tutorial name
        # These are idempotent and safe
        short_code2 = tutorial.split('/')[1..3].join('/').gsub(%r{/tutorials}, '')
        current_mapping['name'][short_code2] = "/#{tutorial.gsub(/md$/, 'html')}"
      end

      # Discover slides
      num_slides = current_mapping['id'].select { |x| x[0] == CATEGORY_SLIDES }.length
      Dir.glob('topics/*/tutorials/*/slides.html').each do |tutorial|
        html_path = "/#{tutorial}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = num_slides.to_s.rjust(5, '0')
          short_code = CATEGORY_SLIDES + short_code_number
          puts "Discovered slides #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
          num_slides += 1
        end

        # Also generate one from topic/tutorial name
        # These are idempotent and safe
        short_code2 = "#{tutorial.split('/')[1..3].join('/').gsub(%r{/tutorials}, '')}/slides"
        current_mapping['name'][short_code2] = "/#{tutorial.gsub(/md$/, 'html')}"
      end

      # Discover FAQs
      all_faqs = Dir.glob('faqs/**/*.md') + Dir.glob('topics/*/faqs/**/*.md') + \
                 Dir.glob('topics/*/tutorials/*/faqs/*.md')
      # Remove symlinked files
      all_faqs = all_faqs.reject { |x| File.symlink?(x) }
      # Reject indexes, readme, etc.
      all_faqs = all_faqs.grep_v(/index.md$/)
      all_faqs = all_faqs.grep_v(/README.md$/)

      num_faqs = current_mapping['id'].select { |x| x[0] == CATEGORY_FAQ }.length
      all_faqs.each do |tutorial|
        html_path = "/#{tutorial.gsub(/md$/, 'html')}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = num_faqs.to_s.rjust(5, '0')
          short_code = CATEGORY_FAQ + short_code_number
          puts "Discovered FAQ #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
          num_faqs += 1
        end
      end

      # Discover news
      num_news = current_mapping['id'].select { |x| x[0] == CATEGORY_NEWS }.length
      Dir.glob('news/_posts/*.md').each do |material|
        m = material.match(%r{news/_posts/(?<year>\d\d\d\d)-(?<month>\d\d)-(?<day>\d\d)-(?<title>.*)\.md})
        html_path = "/news/#{m[:year]}/#{m[:month]}/#{m[:day]}/#{m[:title]}.html"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = num_news.to_s.rjust(5, '0')
          short_code = CATEGORY_NEWS + short_code_number
          puts "Discovered news #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
          num_news += 1
        end
      end

      # Discover learning pathways
      lps = Dir.glob('learning-pathways/*.md')
      lps.reject! { |t| t =~ /index.md/ }
      lps.reject! { |t| t =~ /pathway-example.md/ }

      num_lps = current_mapping['id'].select { |x| x[0] == CATEGORY_PATHWAYS }.length
      lps.each do |tutorial|
        html_path = "/#{tutorial.gsub(/md$/, 'html')}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = num_lps.to_s.rjust(5, '0')
          short_code = CATEGORY_PATHWAYS + short_code_number
          puts "Discovered learning pathway #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
          num_lps += 1
        end
      end

      # Discover events
      events = Dir.glob('events/*.md')
      events.reject! { |t| t =~ /index.md/ }
      events.reject! { |t| t =~ /pathway-example.md/ }

      num_events = current_mapping['id'].select { |x| x[0] == CATEGORY_EVENTS }.length
      events.each do |event|
        html_path = "/#{event.gsub(/md$/, 'html')}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = num_events.to_s.rjust(5, '0')
          short_code = CATEGORY_EVENTS + short_code_number
          puts "Discovered event #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
          num_events += 1
        end
      end

      # Discover workflows
      workflows = Dir.glob('topics/**/workflows/*.ga')

      num_workflows = current_mapping['id'].select { |x| x[0] == CATEGORY_WORKFLOW }.length
      workflows.each do |workflow|
        html_path = "/#{workflow.gsub(/ga$/, 'html')}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = num_workflows.to_s.rjust(5, '0')
          short_code = CATEGORY_WORKFLOW + short_code_number
          puts "Discovered workflow #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
          num_workflows += 1
        end
      end

      current_mapping
    end
  end
end

if $PROGRAM_NAME == __FILE__
  require 'yaml'
  current_mapping = YAML.load_file('metadata/shortlinks.yaml')
  Gtn::Shortlinks.update(current_mapping)
  File.write('metadata/shortlinks.yaml', current_mapping.to_yaml)
end
