#!/usr/bin/env ruby
require 'find'
require 'nokogiri'

def check_indent(file)
  doc = Nokogiri::HTML(File.open(file))
  # Find all <pre> tags
  # <div class="language-plaintext highlighter-rouge"><div class="highlight"><pre class="highlight"><code>
  doc.css('div.language-plaintext.highlighter-rouge div.highlight pre.highlight code').each do |pre|
    # Get the text content of the <pre> tag
    content = pre.text
    # Split the content by newlines
    lines = content.split("\n")

    # If all lines look like URLs:
    if lines.all? { |line| line =~ /:\/\// }
      # If any are space indented
      lines.each do |line|
        if line =~ /^\s+/
          puts "#{file}: Indentation error: #{line}"
        end
      end
    end
  end
end

Find.find('./_site/training-material/topics/') do |path|
  if path =~ /tutorial.*\.html$/
    check_indent(path)
  end
end
#
# check_indent(ARGV[0])
