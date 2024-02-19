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
TUTO_ID = /^(?<tutorial>[A-Za-z0-9_-]*)\/(?<topic>[A-Za-z0-9_-]*):?(?<lang>[a-zA-Z][a-zA-Z])?$/

data = CSV.parse($stdin.read, col_sep: "\t", headers: true, quote_char: '|')
data = data.map do |x|
  x[:timestamp] = x['Timestamp']
  x[:rating] = x['How much did you like this tutorial?']
  x[:pro] = x['What did you like?']
  x[:con] = x['What could be improved?']
  x[:tutorial_topic] = x['Tutorial']
  x[:id] = x['Temporary Tutorial ID']

  x.delete('Timestamp')
  x.delete('How much did you like this tutorial?')
  x.delete('What did you like?')
  x.delete('What could be improved?')
  x.delete('Tutorial')
  x.delete('Your feedback is always anonymous. Also make it confidential (only visible to admins)?')
  x
end

data = data.reject { |x| x[:timestamp].nil? }
           # remove rows with NaN on note, pro and con
           .reject { |x| (x[:rating].nil? && x[:pro].nil? && x[:con].nil?) }
           # Ensure that they properly have a tutorial name and title
           # If this returns NIL, then we'll ignore those.
           # .select { |x| TITLE_TOPIC.match(x[:tutorial_topic]) }
           .reject { |x| x[:id].nil? or x[:id].empty? }
  # .reject{ |x| x[:id] =~ / / }
           .select { |x| TUTO_ID.match(x[:id]) }

# Mutate
data.each do |x| # Various cleanups/extractions
  # replace NaN in rating by 0, convert to int.
  x[:rating] = x[:rating].to_i

  # Extract dates
  x[:date] = DateTime.parse(x[:timestamp]).strftime('%Y-%m-%d')
  x[:month] = DateTime.parse(x[:timestamp]).strftime('%Y-%m')

  # Extract tutorial/topic
  m = TITLE_TOPIC.match(x[:tutorial_topic])
  x.delete('tutorial_topic')
  mq = TUTO_ID.match(x[:id])
  x[:lang] = mq[:lang].downcase if mq[:lang]
  x[:topic] = mq[:topic]
  x[:tutorial] = mq[:tutorial]

  # Replace N/As with nil
  x[:pro] = x[:pro] == 'N/A' || x[:pro] == 'NA' ? nil : x[:pro]
  x[:con] = x[:con] == 'N/A' || x[:con] == 'NA' ? nil : x[:con]

end
  .select do |x|
    !TEST_PHRASES.include? x[:pro] and !TEST_PHRASES.include? x[:con] and !x[:pro].to_s.include? 'Helena Testing'
  end

CSV.open('metadata/feedback.csv', 'wb') do |csv|
  csv << [nil, 'note', 'pro', 'con', 'anonymous', 'tutorial', 'topic', 'month', 'date', 'lang']
  data.each.with_index do |row, index|
    csv << [index, row[:rating], row[:pro], row[:con], nil, row[:tutorial], row[:topic], row[:month], row[:date], row[:lang]]
  end
end

# Also save as yaml
File.open('metadata/feedback2.yaml', 'w') do |file|
  file.write(
    data.group_by { |x| x[:tutorial] }.map do |tutorial, feedback|
      fg = feedback.group_by { |x| x[:topic] }
        .map do |topic, feedback|
          [topic, feedback.map { |x| {"rating" => x[:rating], "pro" => x[:pro], "con" => x[:con], "date" => x[:date], "lang" => x[:lang]} }]
        end
      [tutorial, fg.to_h]
    end.to_h.to_yaml
  )
end

