#!/usr/bin/env ruby
require 'json'
require 'net/http'
require 'uri'
require 'yaml'

# Get the list of toolcats
def fetch_toolcats(server)
  uri = URI.parse("#{server}")
  request = Net::HTTP::Get.new(uri)
  req_options = {
    use_ssl: uri.scheme == 'https',
  }
  response = Net::HTTP.start(uri.hostname, uri.port, req_options) do |http|
    http.request(request)
  end

  begin
    JSON.parse(response.body) do |w|
      w
    end
  rescue StandardError
    {}
  end
end

# Parse the response
toolcats_eu = fetch_toolcats('https://usegalaxy-eu.github.io/usegalaxy-eu-tools/api/labels.json')
# toolcats_eu = File.open('/tmp/tmp.ccFYsrbAa5/usegalaxy-eu-tools/api/labels.json') { |f| YAML.safe_load(f) }
toolcats_org = fetch_toolcats('https://galaxyproject.github.io/usegalaxy-tools/api/labels.json')
toolcats_aus = fetch_toolcats('https://usegalaxy-au.github.io/usegalaxy-au-tools/api/labels.json')
# toolcats_aus = File.open('/tmp/tmp.ccFYsrbAa5/usegalaxy-au-tools/api/labels.json') { |f| YAML.safe_load(f) }
tool_ids = (toolcats_org.keys + toolcats_eu.keys + toolcats_aus.keys).uniq
# tool_ids = toolcats_org.keys
tool_ids.sort!

toolcats = {}
# Cleanup the list
tool_ids.each do |k|
  eu = toolcats_eu[k] || nil
  org = toolcats_org[k] || nil
  aus = toolcats_aus[k] || nil

  # We get N 'votes' for the categories
  values = [eu, org, aus].compact
  # values = [org].compact


  # Majority answer wins
  # set that value to toolcats[k]
  # If there is no majority, pick one.
  # print("#{k} - #{values.length} => #{values.uniq.compact.length}\n")
  if values.length.positive?
    toolcats[k] = values.max_by { |v| v['count'] }
  else
    toolcats[k] = nil
  end
end

# Write the list to a file
File.write('metadata/toolcats.yml', toolcats.to_yaml)
