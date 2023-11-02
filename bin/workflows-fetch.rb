#!/usr/bin/env ruby
require 'json'
require 'net/http'
require 'uri'
require 'yaml'

# Get the list of workflows
def fetch_workflows(server)
  uri = URI.parse("#{server}/api/workflows/")
  request = Net::HTTP::Get.new(uri)
  req_options = {
    use_ssl: uri.scheme == 'https',
  }
  response = Net::HTTP.start(uri.hostname, uri.port, req_options) do |http|
    http.request(request)
  end

  begin
    JSON.parse(response.body).map do |w|
      w['server'] = server
      w
    end
  rescue StandardError
    []
  end
end

# Parse the response
workflows_eu = fetch_workflows('https://usegalaxy.eu')
workflows_org = fetch_workflows('https://usegalaxy.org')
workflows_aus = fetch_workflows('https://usegalaxy.org.au')
workflows = workflows_eu + workflows_org + workflows_aus

# Cleanup the list
workflows.filter! do |w|
  w['published'] == true && w['importable'] == true && w['deleted'] == false && w['hidden'] == false
end

# Group by name + owner
cleaned = workflows.group_by { |w| "#{w['name']}<WFID>#{w['owner']}" }
cleaned = cleaned.map do |_k, v|
  {
    'name' => v[0]['name'],
    'owner' => v[0]['owner'],
    'steps' => v[0]['number_of_steps'],
    'ids' => v.map { |w| [w['server'], w['id']] },
    'tags' => v.map { |w| w['tags'] }.flatten.uniq,
    'updated' => v.map { |w| w['update_time'] }.max,
  }
end

# Write the list to a file
File.write('metadata/workflows.yml', cleaned.to_yaml)
