require 'jekyll'

module Jekyll
  class Abbreviate < Jekyll::Generator
    safe true

    def initialize(config)
      @config = config['abbreviate'] ||= {}
    end

    def generate(site)
      site.pages
        .select { |page| not skip_layout?(page.data['layout']) }
        .each { |page| abbreviate page }
      site.posts.docs
        .select { |post| not skip_layout?(post.data['layout']) }
        .each { |post| abbreviate post }
    end

    private

    def abbreviate(page)
      if page.data.key?('abbreviations') then
        seen = Hash.new
        page.data['abbreviations'].each{|abbr, definition|
          page.content = page.content.gsub(/\{(#{abbr})\}/) {
            if seen.key?(abbr) then
              firstdef = false
            else
              firstdef = true
              seen[abbr] = true
            end

            if firstdef then
              "#{definition} (#{abbr})"
            else
              "<abbr title=\"#{definition}\">#{abbr}</abbr>"
            end
          }
        }
      end
    end

    def skip_layout?(layout)
      to_skip = @config['skip_layouts'] || []

      if to_skip.empty?
        true
      end

      to_skip.include?(layout)
    end
  end
end
