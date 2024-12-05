#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'json'
# [sys.stdout.write(json.dumps(doc, indent=2)) for doc in yaml.safe_load_all(sys.stdin)]
YAML.load_stream($stdin.read).each { |doc| puts JSON.pretty_generate(doc) }
