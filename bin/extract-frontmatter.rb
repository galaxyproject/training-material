#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'json'
fn = ARGV[0]

puts JSON.pretty_generate(YAML.load_file(fn))
