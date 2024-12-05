require 'nokogiri'

data = File.read('_site/training-material/feeds/matrix-month.xml')
doc = Nokogiri::XML(data)
res = doc.css('entry content')[0].children[1].children.to_s
res.gsub!(/utm_source=matrix.*utm_campaign=matrix-news/, 'utm_source=github&utm_medium=release&utm_campaign=release-notes')
puts res
