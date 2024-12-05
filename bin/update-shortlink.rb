#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require './_plugins/gtn/shortlinks'

current_mapping = YAML.load_file('metadata/shortlinks.yaml')
Gtn::Shortlinks.update(current_mapping)
File.write('metadata/shortlinks.yaml', current_mapping.to_yaml)
