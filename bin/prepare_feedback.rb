#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'csv'

TEST_PHRASES = [
  'TESTING',
  'test',
  'Helena testing'
]

TITLE_TOPIC = /^(?<tutorial>.*) \((?<topic>.*)\)$/
TUTO_ID = %r{^(?<topic>[A-Za-z0-9_-]*)/(?<tutorial>[A-Za-z0-9_-]*):?(?<lang>[a-zA-Z][a-zA-Z])?$}

data = CSV.parse($stdin.read, col_sep: "\t", headers: true, quote_char: '|')
data = data.map do |x|
  x[:timestamp] = x['Timestamp']
  x[:rating] = x['How much did you like this tutorial?']
  x[:pro] = x['What did you like?']
  x[:con] = x['What could be improved?']
  x[:id] = x['Tutorial']

  x.delete('Timestamp')
  x.delete('How much did you like this tutorial?')
  x.delete('What did you like?')
  x.delete('What could be improved?')
  x.delete('Tutorial')
  x.delete('Your feedback is always anonymous. Also make it confidential (only visible to admins)?')
  x
end

def lookup_tuto(topic_id, tuto_id)
  @cache ||= {}
  @cache.fetch("#{topic_id}/#{tuto_id}") do |key|
    @cache[key] = nil

    file = "topics/#{topic_id}/tutorials/#{tuto_id}/tutorial.md"
    if File.exist? file
      data = YAML.load_file(file)
      @cache[key] = data['title']
    else
      file = "topics/#{topic_id}/tutorials/#{tuto_id}/slides.html"
      if File.exist? file
        data = YAML.load_file(file)
        @cache[key] = data['title']
      else
        puts "No file for #{topic_id}/#{tuto_id}"
      end
    end
  end
end

def lookup_topic(topic_id)
  @cache ||= {}
  @cache.fetch(topic_id) do |key|
    file = "metadata/#{topic_id}.yaml"
    return nil unless File.exist? file

    data = YAML.load_file(file)
    @cache[key] = data['title']
  end
end

data = data.reject { |x| x[:timestamp].nil? }
           # remove rows with NaN on note, pro and con
           .reject { |x| (x[:rating].nil? && x[:pro].nil? && x[:con].nil?) }
           # Ensure that they properly have a tutorial name and title
           .reject { |x| x[:id].nil? or x[:id].empty? }
           .select { |x| TUTO_ID.match(x[:id]) }

# Mutate
data.each do |x| # Various cleanups/extractions
  # replace NaN in rating by 0, convert to int.
  x[:rating] = x[:rating].to_i

  # Extract dates
  x[:date] = DateTime.parse(x[:timestamp]).strftime('%Y-%m-%d')
  x[:month] = DateTime.parse(x[:timestamp]).strftime('%Y-%m')

  # Extract tutorial/topic
  mq = TUTO_ID.match(x[:id])
  x[:lang] = mq[:lang].downcase if mq[:lang]
  x[:topic_id] = mq[:topic]
  x[:tutorial_id] = mq[:tutorial]

  x[:topic] = lookup_topic(x[:topic_id])
  x[:tutorial] = lookup_tuto(x[:topic_id], x[:tutorial_id])
  # p "Lookup #{x[:topic_id]}/#{x[:tutorial_id]} => #{x[:tutorial]} (#{x[:topic]})"

  # Replace N/As with nil
  x[:pro] = x[:pro] == 'N/A' || x[:pro] == 'NA' ? nil : x[:pro]
  x[:con] = x[:con] == 'N/A' || x[:con] == 'NA' ? nil : x[:con]
end
    .select do |x|
  !TEST_PHRASES.include? x[:pro] and !TEST_PHRASES.include? x[:con] and !x[:pro].to_s.include? 'Helena Testing'
end

CSV.open('metadata/feedback.csv', 'wb') do |csv|
  csv << [nil, 'note', 'pro', 'con', 'anonymous', 'tutorial', 'topic', 'month', 'date']
  data.each.with_index do |row, index|
    csv << [index, row[:rating], row[:pro], row[:con], nil, row[:tutorial], row[:topic], row[:month], row[:date]]
  end
end

# Also save as yaml
yaml_data = data.group_by { |x| x[:topic_id] }.to_h do |topic, feedback|
  fg = feedback.group_by { |x| x[:tutorial_id] }
               .map do |tutorial, feedback2|
    [tutorial, feedback2.map do |x|
                 { 'rating' => x[:rating], 'pro' => x[:pro], 'con' => x[:con], 'date' => x[:date], 'lang' => x[:lang] }
               end]
  end
  [topic, fg.to_h]
end.to_yaml
File.write('metadata/feedback2.yaml', yaml_data)
