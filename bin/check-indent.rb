#!/usr/bin/env ruby
require 'find'
require 'nokogiri'

def check_indent(file)
  doc = Nokogiri::HTML(File.open(file))
  # Find all <pre> tags
  # <div class="language-plaintext highlighter-rouge"><div class="highlight"><pre class="highlight"><code>
  ec = false
  doc.css('div.language-plaintext.highlighter-rouge div.highlight pre.highlight code').each do |pre|
    # Get the text content of the <pre> tag
    content = pre.text
    # Split the content by newlines
    lines = content.split("\n")

    # If all lines look like URLs:
    if lines.all? { |line| line =~ %r{://} }
      # If any are space indented
      lines.each do |line|
        if line =~ /^\s+/
          puts "#{file}: Indentation error: #{line}"
          ec = true
        end
      end
    end
  end
  ec
end

should_exit = false
Find.find('./_site/training-material/topics/') do |path|
  if path =~ (/tutorial.*\.html$/) && check_indent(path)
    should_exit = true
  end
end

if should_exit
  exit 1
end
#
# check_indent(ARGV[0])
