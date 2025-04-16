#!/usr/bin/env ruby
require 'json'
# Ensure that some mandatory files are here before we deploy. Mostly to catch
# people moving/renaming files that aren't as visibly used.

# _site/training-material/api/top-tools.json must exist + have correct structure.
raise 'top-tools.json missing' unless File.exist?('_site/training-material/api/top-tools.json')
top_tools = JSON.parse(File.read('_site/training-material/api/top-tools.json'))
raise 'top-tools.json is not an array' unless top_tools.is_a?(Hash)
raise 'Missing important tools' unless top_tools.key? 'Grep1'
raise 'Wrong structure: missing tool_id' unless top_tools['iuc/circos/circos'].key? 'tool_id'
raise 'Wrong structure: missing tutorials' unless top_tools['iuc/circos/circos'].key? 'tutorials' and top_tools['iuc/circos/circos']['tutorials'].length.positive?
