# frozen_string_literal: true

require 'json'

module Jekyll
  # A generator that creates author pages for each author in the site
  # normally. But this one just disables it.
  class AuthorPageGenerator < Generator
    def generate(_site)
      puts '[GTN/AuthorPages] Disabled'
    end
  end
end
