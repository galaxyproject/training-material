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

addedfiles = `git diff --cached --name-only --ignore-all-space --diff-filter=A #{options[:previousCommit]}`.split("\n")
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
  linkify(d['title'], path.gsub(/.md/, '.html'))
end

def fixNews(n)
  # news/_posts/2021-11-10-api.html => news/2021/11/10/api.html
  n.gsub(%r{news/_posts/(....)-(..)-(..)-(.*.html)}, 'news/\1/\2/\3/\4')
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
  # 'modified': {
  # 'slides': modifiedfiles.select{|x| filterSlides(x)},
  # 'tutorials': modifiedfiles.select{|x| filterTutorials(x)},
  # },
  contributors: `git diff --unified --ignore-all-space #{options[:previousCommit]} CONTRIBUTORS.yaml`
       .split("\n").grep(/^\+[^ ]+:\s*$/).map { |x| x.strip[1..-2] }
}

output = "# GTN News for #{NOW.strftime('%b %d')}"
newsworthy = false

if data[:added][:news].length.positive?
  newsworthy = true
  output += "\n\n## Big News!\n\n"
  output += data[:added][:news].join("\n").gsub(/^/, '- ')
end

if data[:added][:tutorials].length.positive?
  newsworthy = true
  output += "\n\n## #{data[:added][:tutorials].length} new tutorials!\n\n"
  output += data[:added][:tutorials].join("\n").gsub(/^/, '- ')
end

if data[:added][:slides].length.positive?
  newsworthy = true
  output += "\n\n## #{data[:added][:slides].length} new slides!\n\n"
  output += data[:added][:slides].join("\n").gsub(/^/, '- ')
end

if data[:contributors].length.positive?
  newsworthy = true
  output += "\n\n## #{data[:contributors].length} new contributors!\n\n"
  output += data[:contributors].map { |c| linkify("@#{c}", "hall-of-fame/#{c}") }.join("\n").gsub(/^/, '- ')
end

if newsworthy
  if options[:postToMatrix]
    homeserver = if options[:useTestRoom]
                   'https://matrix.org/_matrix/client/r0/rooms/!LcoypptdBmTAvbxeJm:matrix.org/send/m.room.message'
                 else
                   'https://matrix.org/_matrix/client/r0/rooms/!yShkfqncgjcRcRztBU:gitter.im/send/m.room.message'
                 end

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

    uri = URI(homeserver)
    req = Net::HTTP.post(uri, JSON.generate(data), headers)
    puts req
    puts req.body
  else
    puts output
  end
end
