require 'json'
require './_plugins/jekyll-topic-filter.rb'

module Jekyll
  class APIGenerator < Generator

    def generate(site)
      puts "[GTN/API] Disabled"
    end
  end
end

