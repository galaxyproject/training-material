require 'net/http'
require 'date'
require './_plugins/util'
require 'uri'

# Tutorials, FAQ, Video, News, Events are requested.
#
# tutorials, faqs, news, and events are all included in
# https://training.galaxyproject.org/training-material/feeds/single-cell-month.xml
# so we can just look through that for most of the data we need.
# The only thing missing is videos which we'll get separately.

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

# Read the feed from the url
feed = request('https://training.galaxyproject.org/training-material/feeds/single-cell-month.xml').body
# The body of the feed has structure but it's slightly inconvenient.
# The URLs however...
# <a href="https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing/workflows/scrna_mp_celseq.html?utm_source=matrix&amp;utm_medium=newsbot&amp;utm_campaign=matrix-news">CelSeq2: Multi Batch (mm10) (February 22, 2019)</a>
#
# That has everything and a date, can just figure out what 'type' of thing it was by the folder/name.

def classify(url)
  if url =~ /tutorial(_[A-Z_]*)?\.html/
    return 'tutorial'
  elsif url.include? '/events/'
    return 'event'
  elsif url =~ /slides(_[A-Z_]*)?\.html/
    return 'slide'
  elsif url.include? '/faqs/'
    return 'faq'
  elsif url.include? '/workflows/'
    return 'workflow'
  elsif url.include? '/news/'
    return 'news'
  end

  p url

  1/0
end

out = []

feed.scan(/<a href="([^"]+)">([^<]+)<\/a>/).each do |url, title|
  next if url.include? 'gtn-standards-rss.html'
  date = title.match(/\(([^)]+)\)$/)[1]
  date = Date.parse(date).to_s
  out << [date, classify(url)]
end

# Let's add in videos, they're stored in the tutorials
tutos = Dir.glob("topics/single-cell/tutorials/*/tutorial.md")
tutos.each do |tuto|
  meta = safe_load_yaml(tuto)
  next if meta['recordings'].nil? || meta['recordings'].empty?

  meta['recordings'].each do |rec|
    out << [rec['date'], 'video']
  end
end


File.open("single-cell-over-time.tsv", "w") do |f|
  f.puts "date\ttype"
  out.sort_by{|d, t| d}.each do |date, type|
    f.puts "#{date}\t#{type}"
  end
end
