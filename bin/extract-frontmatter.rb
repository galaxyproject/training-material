#!/usr/bin/env ruby
require 'yaml'
require 'json'
fn = ARGV[0]

puts JSON.pretty_generate(YAML.load_file(fn))
