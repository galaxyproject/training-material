#!/usr/bin/env ruby
require 'yaml'
require './_plugins/gtn/shortlinks'

current_mapping = YAML.load_file('metadata/shortlinks.yaml')
Gtn::Shortlinks.update(current_mapping)
File.write('metadata/shortlinks.yaml', current_mapping.to_yaml)
