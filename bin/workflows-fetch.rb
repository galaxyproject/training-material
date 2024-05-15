#!/usr/bin/env ruby
require 'json'
require 'net/http'
require 'uri'
require 'yaml'
require './_plugins/gtn/usegalaxy'

def request(url)
  uri = URI.parse(url)
  request = Net::HTTP::Get.new(uri)
  request['Accept'] = 'application/json'
  req_options = {
    use_ssl: uri.scheme == 'https',
  }
  Net::HTTP.start(uri.hostname, uri.port, req_options) do |http|
    http.request(request)
  end
end

# Get the list of workflows
def fetch_workflows(server)
  begin
    response = request("#{server}/api/workflows/")
  rescue StandardError
    puts "ERROR: Failed to fetch workflows from #{server}"
    return []
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

def fetch_workflowhub
  projects = JSON.parse(request('https://workflowhub.eu/projects').body)
  project_mapping = projects['data'].to_h { |p| [p['id'], p['attributes']['title']] }

  response = request('https://workflowhub.eu/workflows?filter[workflow_type]=galaxy')
  data = JSON.parse(response.body)
  if !data['links']['next'].nil?
    puts 'ERROR: Cannot yet handle multiple pages'
    exit 42
  end
  puts "INFO: Fetching #{data['data'].length} workflows from WorkflowHub"
  data['data'].map.with_index do |w, _i|
    # {"id"=>"14", "type"=>"workflows", "attributes"=>{"title"=>"Cheminformatics - Docking"}, "links"=>{"self"=>"/workflows/14"}}
    wf_info = JSON.parse(request("https://workflowhub.eu#{w['links']['self']}").body)
    creator_list = []

    creator0 = wf_info['data']['attributes']['creators'][0]
    if creator0.nil?
      # Other creators
      other = wf_info['data']['attributes']['other_creators']
      if !other.nil? && other.length.positive?
        creator_list.push(wf_info['data']['attributes']['other_creators'].split(',').map(&:strip))
      end
    else
      # Primary
      creator_list.push("#{creator0['given_name']} #{creator0['family_name']}")
    end
    # Projects
    wf_info['data']['relationships']['projects']['data'].each do |p|
      creator_list.push(project_mapping[p['id']])
    end

    creator_list = creator_list.flatten.compact.uniq

    begin
      r = {
        'name' => wf_info['data']['attributes']['title'],
        'owner' => creator_list.join(', '),
        'number_of_steps' => wf_info['data']['attributes']['internals']['steps'].length,
        'server' => 'https://workflowhub.eu',
        'id' => wf_info['data']['id'],
        'tags' => wf_info['data']['attributes']['tags'].map { |t| t.gsub(/^name:/, '') },
        'update_time' => wf_info['data']['attributes']['updated_at'],
      }
    rescue StandardError
      r = nil
    end
    r
  end.compact
end

# Parse the response
workflows = Gtn::Usegalaxy.servers.map do |server|
  workflows = fetch_workflows(server[:url])
  puts "INFO: Fetched #{workflows.length} workflows from #{server[:name]}"
  workflows
end.flatten

# Cleanup the list
workflows.filter! do |w|
  w['published'] == true && w['importable'] == true && w['deleted'] == false && w['hidden'] == false
end

# Add in WFHub workflows
workflows += fetch_workflowhub

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
