#!/usr/bin/env ruby
# frozen_string_literal: true

require 'json'
require 'kramdown'
require 'uri'
require 'net/http'
require 'yaml'
require 'optparse'

options = {}
OptionParser.new do |opt|
  opt.on('-p', '--previous-commit PREVIOUS_COMMIT') { |o| options[:previousCommit] = o }
  opt.on('--matrix-post') { |o| options[:postToMatrix] = o }
  opt.on('--matrix-use-test-room') { |o| options[:useTestRoom] = o }
end.parse!

if options[:previousCommit].nil?
  puts 'Missing a previous commit to compare against'
  exit 1
end

$rooms = {
  'test' => {
    server: 'https://matrix.org',
    room: '!LcoypptdBmTAvbxeJm:matrix.org',
  },
  'default' => {
    server: 'https://matrix.org',
    room: '!yShkfqncgjcRcRztBU:gitter.im',
  },
  'single-cell' => {
    server: 'https://matrix.org',
    room: '!yuLoaCWKpFHkWPmVEO:gitter.im',
  }
}

addedfiles = `git diff --cached --name-only --ignore-all-space --diff-filter=A #{options[:previousCommit]}`.split("\n")
modifiedfiles = `git diff --cached --name-only --ignore-all-space --diff-filter=M #{options[:previousCommit]}`.split("\n")
# modifiedfiles = `git diff --cached --name-only --ignore-all-space
# --diff-filter=M #{options[:previousCommit]}`.split("\n")

NOW = Time.now
CONTRIBUTORS = YAML.load_file('CONTRIBUTORS.yaml')

# new   news
# new   slidevideos
# new   contributors   Done
# new   tutorials      Done
# new   slides         Done

def filterTutorials(x)
  x =~ %r{topics/.*/tutorials/.*/tutorial.*\.md}
end

def filterSlides(x)
  x =~ %r{topics/.*/tutorials/.*/slides.*\.html}
end

def onlyEnabled(x)
  tutorial_meta = YAML.load_file(x)
  tutorial_enabled = tutorial_meta.fetch('enable', true)

  topic = x.split('/')[1]
  topic_meta = YAML.load_file("metadata/#{topic}.yaml")
  topic_enabled = topic_meta.fetch('enable', true)

  tutorial_enabled and topic_enabled
end

def linkify(text, path)
  "[#{text}](https://training.galaxyproject.org/training-material/#{path}?utm_source=matrix&utm_medium=newsbot&utm_campaign=matrix-news)"
end

def printableMaterial(path)
  d = YAML.load_file(path)
  {md: linkify(d['title'], path.gsub(/.md/, '.html')),
  path: path}
end

def fixNews(n)
  # news/_posts/2021-11-10-api.html => news/2021/11/10/api.html
  n[:md].gsub(%r{news/_posts/(....)-(..)-(..)-(.*.html)}, 'news/\1/\2/\3/\4')
end

data = {
  added: {
    slides: addedfiles
       .select { |x| filterSlides(x) }
       .select { |x| onlyEnabled(x) }
       .map { |x| printableMaterial(x) },
    tutorials: addedfiles
       .select { |x| filterTutorials(x) }
       .select { |x| onlyEnabled(x) }
       .map { |x| printableMaterial(x) },
    news: addedfiles.grep(%r{news/_posts/.*\.md}).map { |x| printableMaterial(x) }.map { |n| fixNews(n) }
  },
  modified: {
    slides: modifiedfiles
       .select { |x| filterSlides(x) }
       .select { |x| onlyEnabled(x) }
       .map { |x| printableMaterial(x) },
    tutorials: modifiedfiles
       .select { |x| filterTutorials(x) }
       .select { |x| onlyEnabled(x) }
       .map { |x| printableMaterial(x) },
  },
  contributors: `git diff --unified --ignore-all-space #{options[:previousCommit]} CONTRIBUTORS.yaml`
       .split("\n").grep(/^\+[^ ]+:\s*$/).map { |x| x.strip[1..-2] }
}

def titleize(t)
  t.gsub(/-/, ' ').gsub(/\w+/) do |word|
    word.capitalize
  end
end

def format_news(news)
  output = ""
  if news.length.positive?
    output += "\n\n## Big News!\n\n"
    output += news.join("\n").gsub(/^/, '- ')
  end
  output
end

def format_tutorials(added, modified, kind: "tutorials")
  output =""
  if added.length.positive? or modified.length.positive?
    output += "\n\n## #{added.length + modified.length} #{kind}!"
  end

  if added.length.positive?
    output += "\n\nNew Tutorials:\n\n"
    output += added.map{|n| n[:md]}.join("\n").gsub(/^/, '- ')
  end

  if modified.length.positive?
    output += "\n\nUpdated Tutorials:\n\n"
    output += modified.map{|n| n[:md]}.join("\n").gsub(/^/, '- ')
  end
  output
end

def build_news(data, filter: nil)
  infix = filter.nil? ? '' : titleize(filter)
  output = "# GTN #{infix} News for #{NOW.strftime('%b %d')}"
  newsworthy = false

  if filter.nil?
  output += format_news(data[:added][:news])
  newsworthy = newsworthy | format_news(data[:added][:news]).length
  end

  o = format_tutorials(
    data[:added][:tutorials].select{|n| filter.nil? || n[:path] =~ /topics\/#{filter}/ },
    data[:modified][:tutorials].select{|n| filter.nil? || n[:path] =~ /topics\/#{filter}/ }
  )

  output += o
  newsworthy = newsworthy | o.length

  o = format_tutorials(
    data[:added][:slides].select{|n| filter.nil? || n[:path] =~ /topics\/#{filter}/ },
    data[:modified][:slides].select{|n| filter.nil? || n[:path] =~ /topics\/#{filter}/ }
  )
  output += o
  newsworthy = newsworthy | o.length

  if filter.nil? && data[:contributors].length.positive?
    newsworthy = true
    output += "\n\n## #{data[:contributors].length} new contributors!\n\n"
    output += data[:contributors].map { |c| linkify("@#{c}", "hall-of-fame/#{c}") }.join("\n").gsub(/^/, '- ')
  end

  return output, newsworthy
end

def send_news(output, options, channel: 'default')
  if options[:postToMatrix]
    homeserver = $rooms[channel]
    pp homeserver

    data = {
      'msgtype' => 'm.notice',
      'body' => output,
      'format' => 'org.matrix.custom.html',
      'formatted_body' => Kramdown::Document.new(output).to_html,
    }

    headers = {
      'Authorization' => "Bearer #{ENV.fetch('MATRIX_ACCESS_TOKEN', nil)}",
      'Content-type' => 'application/json',
    }

    uri_send_message = URI("#{homeserver[:server]}/_matrix/client/r0/rooms/#{homeserver[:room]}/send/m.room.message")
    req = Net::HTTP.post(uri_send_message, JSON.generate(data), headers)
    # Parse response
    resp = JSON.parse(req.body)
    puts resp

    if resp['errcode'] == 'M_FORBIDDEN'
      if resp['error'] =~ /not in room/
        puts "Not in room, attempting to join"
        # Join room
        #  POST /_matrix/client/v3/join/{roomIdOrAlias}
        uri_join = URI("#{homeserver[:server]}/_matrix/client/v3/join/#{homeserver[:room]}")
        req = Net::HTTP.post(uri_join, JSON.generate({}), headers)
        # Parse response
        resp = JSON.parse(req.body)

        # Now we're safe to re-try
        if resp.key?('room_id')
          req = Net::HTTP.post(uri_send_message, JSON.generate(data), headers)
          # Parse response
          resp = JSON.parse(req.body)
          puts resp
        end
      end
    end
  else
    puts output
  end
end

output, newsworthy = build_news(data)
if newsworthy
  channel = options[:useTestRoom] ? 'test' : 'default'
  send_news(output, options, channel: channel)
end

output, newsworthy = build_news(data, filter: 'single-cell')
if newsworthy
  send_news(output, options, channel: 'single-cell')
end
