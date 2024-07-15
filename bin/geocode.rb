#!/usr/bin/env ruby
require 'yaml'
require 'net/http'
require 'json'
require 'uri'
require 'date'

def request(url)
  uri = URI.parse(url)
  request = Net::HTTP::Get.new(uri)
  request['Accept'] = 'application/json'
  request['User-Agent'] = 'GTN-geocode/1.0 (+https://github.com/galaxyproject/training-material)'
  req_options = {
    use_ssl: uri.scheme == 'https',
  }
  Net::HTTP.start(uri.hostname, uri.port, req_options) do |http|
    http.request(request)
  end
end

# loop across events
events = Dir.glob('events/*.md')
events.each do |event|
  # Load as yaml
  begin
    event_data = YAML.load_file(event, permitted_classes: [Date])
  rescue StandardError
    event_data = YAML.load_file(event)
  end
  # Check if it is already geocoded
  if !event_data.key?('location')
    next
  end

  if !event_data['location'].key?('geo')
    # Geocode
    loc = {
      'street' => event_data['location'].fetch('address', nil),
      'city' => event_data['location'].fetch('city', nil),
      'country' => event_data['location'].fetch('country', nil),
      'postcode' => event_data['location'].fetch('postcode', nil),
    }

    if loc.values.compact.empty?
      # Online event
      next
    end

    # Nominatim
    #
    # https://nominatim.openstreetmap.org/search.php?street=Dr+Molewaterplein+40&city=Rotterdam&state=Zuid+Holland&country=the+Netherlands&postalcode=3015GD&dedupe=0&limit=1&format=jsonv2
    kv = loc.map { |k, v| "#{k}=#{URI.encode_www_form_component(v)}" }.join('&')
    url = "https://nominatim.openstreetmap.org/search.php?#{kv}&format=jsonv2"
    data = JSON.parse(request(url).body)

    # Open the file, replace ^location:$ with the geo text appended
    if data.length.positive?
      event_data['location']['geo'] = {
        'lat' => data[0]['lat'],
        'lon' => data[0]['lon'],
      }
      contents = File.read(event)
      contents = contents.gsub(/^location:$/,
                               "location:\n  geo:\n    lat: #{data[0]['lat']}\n    lon: #{data[0]['lon']}")
      File.open(event, 'w') { |file| file.puts contents }
    end
  end
end
