module Jekyll
  class AuthorPageGenerator < Generator
    safe true

    def generate(site)
      if site.layouts.key? 'contributor_index'
        dir = 'hall-of-fame'
        site.data['contributors'].each_key do |contributor|
          site.pages << AuthorPage.new(site, site.source, File.join(dir, contributor), contributor)
        end
      end
    end
  end

  # A Page subclass used in the `AuthorPageGenerator`
  class AuthorPage < Page
    def initialize(site, base, dir, contributor)
      @site = site
      @base = base
      @dir  = dir
      @name = 'index.html'

      self.process(@name)
      self.read_yaml(File.join(base, '_layouts'), 'contributor_index.html')
      self.data['contributor'] = contributor

      if @site.data['contributors'][contributor] and @site.data['contributors'][contributor].has_key?('name') then
        name = @site.data['contributors'][contributor]['name']
      else
        name = contributor
      end

      self.data['personname'] = name

      self.data['title'] = "GTN Contributor: #{name}"
    end
  end
end
