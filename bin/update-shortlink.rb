#!/usr/bin/env ruby
require 'yaml'
require './_plugins/gtn/shortlinks.rb'

current_mapping = YAML.load_file('metadata/shortlinks.yaml')
Gtn::Shortlinks.update(current_mapping)
File.open('metadata/shortlinks.yaml', 'w') do |file|
  file.write(current_mapping.to_yaml)
end
