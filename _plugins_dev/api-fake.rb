require 'json'

module Jekyll
  # API Generation Disabled
  class APIGenerator < Generator
    def generate(_site)
      puts '[GTN/API] Disabled'
    end
  end
end
