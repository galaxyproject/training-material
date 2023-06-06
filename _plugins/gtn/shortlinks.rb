# frozen_string_literal: true

module Gtn
  # This module is responsible for generating shortlinks for tutorials and FAQs
  module Shortlinks
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
          short_code_number = current_mapping['id'].select { |x| x[0] == 'T' }.length.to_s.rjust(5, '0')
          short_code = "T#{short_code_number}"
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
          short_code_number = current_mapping['id'].select { |x| x[0] == 'S' }.length.to_s.rjust(5, '0')
          short_code = "S#{short_code_number}"
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
      Dir.glob('faqs/**/*.md').each do |tutorial|
        html_path = "/#{tutorial.gsub(/md$/, 'html')}"
        # If it's not already mapped by a key, add it.
        if !mapped?(html_path, current_mapping)
          # Generate a short code
          short_code_number = current_mapping['id'].select { |x| x[0] == 'F' }.length.to_s.rjust(5, '0')
          short_code = "F#{short_code_number}"
          puts "Discovered FAQ #{short_code}"
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
