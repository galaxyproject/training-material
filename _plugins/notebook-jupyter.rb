require 'json'
require 'mkmf'
require 'fileutils'
require './_plugins/notebook'

module MakeMakefile::Logging
  @logfile = File::NULL
end

module Jekyll
  class JupyterNotebookGenerator < Generator
    safe true

    def generate(site)
      if find_executable('notedown').nil?
        puts '[GTN/Notebooks] We could not find the notedown executable, so, notebooks will not be rendered.'
        return
      end

      # For every tutorial with the 'notebook' key in the page data
      site.pages.select { |page| GTNNotebooks.notebook_filter(page.data) }.each do |page|
        # We get the path to the tutorial source
        dir = File.dirname(File.join('.', page.url))
        fn = File.join('.', page.url).sub(/html$/, 'md')
        notebook_language = page.data['notebook'].fetch('language', 'python')

        # Tag our source page
        page.data['tags'] = [] unless page.data.has_key? 'tags'
        page.data['tags'].push('jupyter-notebook')

        puts "[GTN/Notebooks] Rendering #{notebook_language} #{fn}"
        last_modified = begin page.last_modified.to_s rescue Time.new.to_s end
        notebook = GTNNotebooks.render_jupyter_notebook(page.data, page.content, page.url, last_modified, notebook_language, site)

        # Write it out!
        page2 = PageWithoutAFile.new(site, '', dir, 'tutorial.md.ipynb')
        page2.content = notebook
        page2.data['layout'] = nil
        page2.data['citation_target'] = 'jupyter'
        site.pages << page2
      end
    end
  end
end
