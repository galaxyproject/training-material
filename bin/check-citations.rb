require 'find'
require 'bibtex'
require 'citeproc/ruby'
require 'csl/styles'

# Load our citation library
a = BibTeX::Bibliography.new
Find.find('./topics') do |path|
  if FileTest.directory?(path)
    if File.basename(path).start_with?('.')
      Find.prune       # Don't look any further into this directory.
    else
      next
    end
  else
    if path =~ /bib$/ then
      b = BibTeX.open(path)
      for x in b
        # Record the bib path.
        x._path = path
        a << x
      end
    end
  end
end


count = 0
# Get actual used citations
results = `grep '{% cite' -R topics --line-number`
results.split("\n").each{ |x|
  (file, line, text) = x.split(":", 3)

  matches = []
  text.scan(/{% cite\s+([^%]*)\s*%}/) do |c|
    matches << [c[0].strip, $~.offset(0)[0]]
  end

  matches.select{|(id, _)| a[id].nil? }.each{|(id, offset)|
    count += 1
    puts %Q(#{file}:#{line}:#{offset} This citation was not found "#{id}")
  }
}

# TODO: Re-enable when #2788 is fixed.
#if count > 0
#  exit 1
#end
