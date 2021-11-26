module TopicFilter
  def self.list_topics(site)
    site.data.select{|k, v| v.is_a?(Hash) && v.has_key?('type')}.map{|k, v| k}
  end

  def self.topic_filter(site, topic_name)
    if not site.data.has_key?('cache_topic_filter')
      site.data['cache_topic_filter'] = Hash.new

      # For each topic
      self.list_topics(site).each{|topic|
        site.data['cache_topic_filter'][topic] = self.run_topic_filter(site.pages, topic)
      }
    end

    site.data['cache_topic_filter'][topic_name]
  end

  def self.run_topic_filter(pages, topic_name)
    # Arrays that will store all introduction slides and tutorials we discover.
    resource_intro = []
    resource_pages = []

    # In order to speed up queries later, we'll store a set of "interesting"
    # pages (i.e. things that are under `topic_name`)
    interesting = {}
    for page in pages do
      page_parts = page.url.split('/')
      # Skip anything outside of topics.
      if not page.url.include?('/topics/') then next end
      page.data['url'] = page.url
      # If the path is long enough and it is under the topic name
      if page_parts.length > 3 and page_parts[2] == topic_name
        # Automate the tutorial-name thing. This writes back to the shared
        # data structure.
        page.data['topic_name'] = page_parts[2]

        # Slides are one directory level shorter and re-use name for their identity.
        if page_parts[3] == 'slides' then
          page.data['tutorial_name'] = page_parts[4].sub(/\.html$/, '')
        else
          page.data['tutorial_name'] = page_parts[4]
        end

        # And then store in our interesting stuff
        key = page.url.sub(/^\//, '') # strip leading slash since later queries don't have it.
        interesting[key] = page
      end
    end

    # Theory is as follows:
    #
    # for each folder we need to extract
    #   if any slides/workflows/tutorials/etc
    #
    # and we need to update this object with the tutorial's title in order to
    # make the object complete. It is called a page_obj or a resource within
    # the context of the script but basically it's the same as the old
    # entries from the metadata.yaml file. A thing that has slides:yes/no,
    # hands_on:yes/no, and all of the other associated data so we can
    # construct the single line in the topic page with the tutorial's
    # resources.

    # All of the files in the intro slides and tutorials directories
    intro_slides = Dir.glob("topics/#{topic_name}/slides/*").sort
    tutorial_folders = Dir.glob("topics/#{topic_name}/tutorials/*").sort

    # Turn these into an object. Override the slides to true (since they're
    # slides) and the type to introduction, to prevent people from having to
    # set those variables. All other variables will be copied directly with
    # the `page.data.dup` so people can use external hands ons and similar.
    for intro_slide in intro_slides do
      page = interesting[intro_slide]
      page_obj = page.data.dup
      page_obj['slides'] = true
      page_obj['type'] = 'introduction'
      resource_intro.push(page_obj)
    end

    # Look in every /topic/*/tutorials/* folder, and turn these disparate
    # resources into a page_obj as well. Most variables are copied directly,
    # either from a tutorial, or a slides (if no tutorial is available.) This
    # means we do not (cannot) support external_slides AND external_handson.
    # This is probably a sub-optimal situation we'll end up fixing someday.
    for folder in tutorial_folders do
      # Discover resources, defined as any file in a given tutorial folder
      resources = Dir.glob("#{folder}/*").map{ |a| a.split('/')[-1] }

      # We will compare this against the 'interesting' hash we built earlier,
      # let's find everything from interesting that is in this specific topic
      # AND tutorial (interesting only contained things in this topic).
      # We add the trailing slash because tutorial directories sharing a
      # common prefix would cause issues.
      known_pages = interesting.select{|a| a.include?("#{folder}/")}

      # We need to extract metadata from some of the pages to determine the
      # proper 'title' attribute for a tutorial. We can only pull this if we
      # have the page object in jekyll to work with. So we'll look for the
      # specific keys in `known_pages` (keyed on filesystem paths) and select
      # for the things we're interested in.
      tutorial_page_keys = known_pages.keys.select{|a| a.include?('tutorial.html')}
      slides_page_keys   = known_pages.keys.select{|a| a.include?('slides.html')}

      # We'll handle slides first and have hands-on override.
      page = false
      slide_has_video = false
      slide_translations = []
      if slides_page_keys.length == 1 then
        page = interesting[slides_page_keys[0]]
        slide_has_video = page.data.has_key?('video') and page['video']
        slide_translations = page.data.fetch('translations', [])
      end

      tutorial_translations = []
      if tutorial_page_keys.length == 1 then
        page = interesting[tutorial_page_keys[0]]
        tutorial_translations = page.data.fetch('translations', [])
      end

      # If no tutorial OR slides are found, then we have an issue
      if page == false then
        puts "Error? No tutorial OR slides found in #{folder}. We saw #{known_pages.keys}"
        next
      end

      # Otherwise clone the metadata from it which works well enough.
      page_obj = page.data.dup

      # Sometimes `hands_on` is set to something like `external`, in which
      # case it is important to not override it. So we only do that if the
      # key isn't already set. Then we choose to set it to a test for the
      # tutorial being present. We probably don't need to test both, but it
      # is hard to follow which keys are which and safer to test for both in
      # case someone edits the code later. If either of these exist, we can
      # automatically set `hands_on: true`
      if not page_obj.has_key?("hands_on") then
        page_obj['hands_on'] = resources.include?('tutorial.md') or resources.include?('tutorial.html')
      end

      # Same for slides, if there's a resource by that name, we can
      # automatically set `slides: true`
      if not page_obj.has_key?("slides") then
        page_obj['slides'] = resources.include?('slides.html')
      end

      # Similar as above.
      page_obj['workflows'] = resources.include?('workflows')
      page_obj['tours'] = resources.include?('tours')
      page_obj['video'] = slide_has_video
      page_obj['translations'] = Hash.new
      page_obj['translations']["tutorial"] = tutorial_translations
      page_obj['translations']["slides"] = slide_translations
      page_obj['translations']['video'] = slide_has_video # Just demand it?
      # I feel less certain about this override, but it works well enough in
      # practice, and I did not find any examples of `type: <anything other
      # than tutorial>` in topics/*/tutorials/*/tutorial.md but that doesn't
      # make it future proof.
      page_obj['type'] = 'tutorial'

      if page_obj.has_key?("enable") and !page_obj['enable'] then
        if ! page_obj.has_key? 'tags'
          page_obj['tags'] = Array.new
        end
        page_obj['tags'].push('work-in-progress')
      end

      # Push onto our stack.
      resource_pages.push(page_obj)
    end

    # The complete resources we'll return is the introduction slides first
    # (regardless of alphabetisation), and then the rest of the pages.
    resource_pages = resource_intro + resource_pages.sort_by{ |k| k["title"].downcase}

    if resource_pages.length == 0 then
      puts "Error? Could not find any relevant pages for #{topic_name}"
    end

    # Apparently return is optional?
    resource_pages
  end
end


module Jekyll
  module ImplTopicFilter
    def most_recent_contributors(contributors, count)
      # Remove non-hof
      hof = contributors.select{ |k, v| v.fetch("halloffame", "yes") != "no" }
      # Get keys + sort by joined date
      hof_k = hof.keys.sort{ |x, y|
        hof[y].fetch('joined', '2016-01') <=> hof[x].fetch('joined', '2016-01')
      }

      # Transform back into hash
      Hash[hof_k.slice(0, count).collect{|k| [k, hof[k]]}]
    end

    def filter_recent_modified(tutorials, count)
      latest = tutorials.sort{ |x, y|
        x.data['last_modified_at'].format = '%s' # Originally %d-%b-%y
        y.data['last_modified_at'].format = '%s' # Originally %d-%b-%y

        y.data['last_modified_at'].to_s <=> x.data['last_modified_at'].to_s
      }
      latest.slice(0, count)
    end

    def topic_count(resources)
      # Count lines in the table except introduction slides
      resources.select{ |a| a['type'] != 'introduction' }.length
    end

    def topic_filter(site, topic_name)
      TopicFilter.topic_filter(site, topic_name)
    end

  end
end

Liquid::Template.register_filter(Jekyll::ImplTopicFilter)
