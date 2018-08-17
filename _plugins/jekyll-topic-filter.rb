module Jekyll
  module TopicFilter
    def topic_filter(pages, topic_name)
      relevant_pages = []

      # Must provide a topic name.
      if topic_name.nil? then
        return relevant_pages
      end

      interesting = {}
      for page in pages do
        page_parts = page.url.split('/')
        # Skip anything outside of topics.
        if not page.url.include?('/topics/') then next end

        if page_parts.length > 3 and page_parts[2] == topic_name
          key = page.url.sub(/^\//, '')
          interesting[key] = page
        end
      end

      # Theory is as follows:
      #
      # for each folder we need to extract
      #   if any slides/workflows/tutorials/etc
      #
      # and we need to update this object with the tutorial's title in order to
      # make the object complete.

      intro_slides = Dir.glob("topics/#{topic_name}/slides/*")
      tutorial_folders = Dir.glob("topics/#{topic_name}/tutorials/*").sort
      puts "=> #{topic_name}"

      for intro_slide in intro_slides do
        page = interesting[intro_slide]
        #puts "int> #{page.data}"
        page_obj = page.data.dup
        page_obj['slides'] = true
        page_obj['type'] = 'introduction'
        relevant_pages.push(page_obj)
      end

      for folder in tutorial_folders do
        # Discover resources
        resources = Dir.glob("#{folder}/*").map{ |a| a.split('/')[-1] }
        #puts "fol> #{folder} #{resources}"

        # First pull out all of the things in 'interesting' that are in the
        # folder we're currently examining
        known_pages = interesting.select{|a| a.include?(folder + '/')}
        # Next we'll look for the specific page keys. Should not encounter multiple.
        tutorial_page_keys = known_pages.keys.select{|a| a.include?('tutorial.html')}
        slides_page_keys   = known_pages.keys.select{|a| a.include?('slides.html')}
        # And finally obtain the main place for metadata, slies if they exist,
        # or better a hands-on
        page = false
        # We'll handle slides first and have hands-on override.
        if slides_page_keys.length == 1 then
          page = interesting[slides_page_keys[0]]
        end

        if tutorial_page_keys.length == 1 then
          page = interesting[tutorial_page_keys[0]]
        end

        if page == false then
          puts "Error? No tutorial OR slides found in #{folder}. We saw #{known_pages.keys}"
          next
        end

        page_obj = page.data.dup
        page_obj['slides'] = resources.include?('slides.html')
        page_obj['hands-on'] = resources.include?('tutorial.md')
        page_obj['workflows'] = resources.include?('workflows')
        page_obj['tours'] = resources.include?('tours')
        page_obj['type'] = 'tutorial'

        relevant_pages.push(page_obj)
      end

      if relevant_pages.length == 0 then
        puts "Error? Could not find any relevant pages for #{topic_name}"
      end

      sorted = relevant_pages.sort_by{ |k| k["title"] }
      sorted
    end
  end
end

Liquid::Template.register_filter(Jekyll::TopicFilter)
