module Gtn
  module Shortlinks
    def self.is_mapped(tutorial, current_mapping)
      return current_mapping['id'].values.include? tutorial
    end

    def self.update(current_mapping)
      if not current_mapping.has_key? 'id'
        current_mapping['id'] = {}
      end

      if not current_mapping.has_key? 'name'
        current_mapping['name'] = {}
      end

      # Discover tutorials
      Dir.glob('topics/*/tutorials/*/tutorial.md').each do |tutorial|
        html_path = '/' + tutorial.gsub(/md$/, 'html')
        # If it's not already mapped by a key, add it.
        if not is_mapped(html_path, current_mapping)
          # Generate a short code
          short_code = 'T' + (current_mapping['id'].select{|x| x[0] == 'T'}.length).to_s.rjust(5, '0')
          puts "Discovered tutorial #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
        end

        # Also generate one from topic/tutorial name
        # These are idempotent and safe
        short_code2 = tutorial.split('/')[1..3].join('/').gsub(/\/tutorials/, '')
        current_mapping['name'][short_code2] = '/' + tutorial.gsub(/md$/, 'html')
      end

      # Discover slides
      Dir.glob('topics/*/tutorials/*/slides.html').each do |tutorial|
        html_path = '/' + tutorial
        # If it's not already mapped by a key, add it.
        if not is_mapped(html_path, current_mapping)
          # Generate a short code
          short_code = 'S' + (current_mapping['id'].select{|x| x[0] == 'S'}.length).to_s.rjust(5, '0')
          puts "Discovered slides #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
        end

        # Also generate one from topic/tutorial name
        # These are idempotent and safe
        short_code2 = tutorial.split('/')[1..3].join('/').gsub(/\/tutorials/, '') + '/slides'
        current_mapping['name'][short_code2] = '/' + tutorial.gsub(/md$/, 'html')
      end

      # Discover FAQs
      Dir.glob('faqs/**/*.md').each do |tutorial|
        html_path = '/' + tutorial.gsub(/md$/, 'html')
        # If it's not already mapped by a key, add it.
        if not is_mapped(html_path, current_mapping)
          # Generate a short code
          short_code = 'F' + (current_mapping['id'].select{|x| x[0] == 'F'}.length).to_s.rjust(5, '0')
          puts "Discovered FAQ #{short_code}"
          # If the target of this flavour of short code isn't already in here, then add it
          current_mapping['id'][short_code] = html_path
        end
      end
      current_mapping
    end
  end
end

if $0 == __FILE__
  require 'yaml'
  current_mapping = YAML.load_file('metadata/shortlinks.yaml')
  Gtn::Shortlinks.update(current_mapping)
  File.open('metadata/shortlinks.yaml', 'w') do |file|
    file.write(current_mapping.to_yaml)
  end
end
