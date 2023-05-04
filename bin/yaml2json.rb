#!/usr/bin/env ruby
require 'yaml'
require 'json'
# [sys.stdout.write(json.dumps(doc, indent=2)) for doc in yaml.safe_load_all(sys.stdin)]
YAML.load_stream(STDIN.read).each { |doc| puts JSON.pretty_generate(doc) }
