# frozen_string_literal: true

module Gtn
  # This module is responsible for generating shortlinks for tutorials and FAQs
  module Shortlinks
    CATEGORY_TUTORIAL = 'T'
    CATEGORY_SLIDES = 'S'
    CATEGORY_FAQ = 'F'
    CATEGORY_NEWS = 'N'
    CATEGORY_PATHWAYS = 'P'

    def self.mapped?(tutorial, current_mapping)
      current_mapping['id'].values.include? tutorial
    end

    def self.update(current_mapping)
      current_mapping['id'] = {} if !current_mapping.key? 'id'

      current_mapping['name'] = {} if !current_mapping.key? 'name'

      # Discover tutorials
      Dir.glob('topics/*/tutorials/*/tutorial.md').each do |tutorial|
        html_path = "/#{tutorial.gsub(/md$/, 'html')}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = current_mapping['id'].select { |x| x[0] == CATEGORY_TUTORIAL }.length.to_s.rjust(5, '0')
          short_code = CATEGORY_TUTORIAL + short_code_number
          puts "Discovered tutorial #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
        end

        # Also generate one from topic/tutorial name
        # These are idempotent and safe
        short_code2 = tutorial.split('/')[1..3].join('/').gsub(%r{/tutorials}, '')
        current_mapping['name'][short_code2] = "/#{tutorial.gsub(/md$/, 'html')}"
      end

      # Discover slides
      Dir.glob('topics/*/tutorials/*/slides.html').each do |tutorial|
        html_path = "/#{tutorial}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = current_mapping['id'].select { |x| x[0] == CATEGORY_SLIDES }.length.to_s.rjust(5, '0')
          short_code = CATEGORY_SLIDES + short_code_number
          puts "Discovered slides #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
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

      all_faqs.each do |tutorial|
        html_path = "/#{tutorial.gsub(/md$/, 'html')}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = current_mapping['id'].select { |x| x[0] == CATEGORY_FAQ }.length.to_s.rjust(5, '0')
          short_code = CATEGORY_FAQ + short_code_number
          puts "Discovered FAQ #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
        end
      end

      # Discover news
      Dir.glob('news/_posts/*.md').each do |material|
        m = material.match(%r{news/_posts/(?<year>\d\d\d\d)-(?<month>\d\d)-(?<day>\d\d)-(?<title>.*)\.md})
        html_path = "/news/#{m[:year]}/#{m[:month]}/#{m[:day]}/#{m[:title]}.html"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = current_mapping['id'].select { |x| x[0] == CATEGORY_NEWS }.length.to_s.rjust(5, '0')
          short_code = CATEGORY_NEWS + short_code_number
          puts "Discovered news #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
        end
      end

      # Discover learning pathways
      lps = Dir.glob('learning-pathways/*.md')
      lps.reject! { |t| t =~ /index.md/ }
      lps.reject! { |t| t =~ /pathway-example.md/ }

      lps.each do |tutorial|
        html_path = "/#{tutorial.gsub(/md$/, 'html')}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = current_mapping['id'].select { |x| x[0] == CATEGORY_PATHWAYS }.length.to_s.rjust(5, '0')
          short_code = CATEGORY_PATHWAYS + short_code_number
          puts "Discovered slides #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
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
