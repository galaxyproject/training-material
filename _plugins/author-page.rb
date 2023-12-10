# frozen_string_literal: true

require './_plugins/gtn'

module Jekyll
  ##
  # This class generates the GTN's author pags
  class AuthorPageGenerator < Generator
    safe true

    ##
    # This extracts the contributions and pushes them on to an existing
    # datastructure, modifying it in the process. It's pretty gross.
    #
    # Params
    # +t+:: The tutorial, slide, or news item
    # +datastructure+:: The hash of contributors that the author information should be pushed onto
    # +flat+:: Whether the datastructure is a flat array or a nested array
    #
    # Returns
    # +datastructure+:: The modified datastructure
    def pusher(t, datastructure, flat)
      if t.data.key?('contributors')
        if flat
          t.data['contributors'].each { |c| datastructure[c].push(t) }
        else
          t.data['contributors'].each { |c| datastructure[c].push([t, nil]) }
        end
      elsif t.data.key?('contributions')
        t.data['contributions'].each do |contribution_type, contributor|
          contributor.each do |c|
            if flat
              datastructure[c].push(t)
            else
              datastructure[c].push([t, contribution_type])
            end
          end
        end
      end

      t.data['maintainers'].each { |c| datastructure[c].push([t, 'maintainer']) } if t.data.key?('maintainers')
      t.data['funding'].each { |c| datastructure[c].push([t, 'funding']) } if t.data.key?('funding')

      datastructure
    end

    ##
    # This generates the author pages
    # Params
    # +site+:: The site object
    def generate(site)
      return unless site.layouts.key? 'contributor_index'

      dir = 'hall-of-fame'

      # pre-calculating this hash saves about 4.9 seconds off the previous
      # build time of 5 seconds.
      tutorials_by_author = Hash.new { |hash, key| hash[key] = [] }
      learning_pathways_by_author = Hash.new { |hash, key| hash[key] = [] }
      slides_by_author = Hash.new { |hash, key| hash[key] = [] }
      news_by_author = Hash.new { |hash, key| hash[key] = [] }
      has_philosophy = Hash.new { false }

      site.pages.each do |t|
        # Tutorials
        pusher(t, tutorials_by_author, false) if t['layout'] == 'tutorial_hands_on'

        # Slides
        if !%w[base_slides introduction_slides tutorial_slides].index(t['layout']).nil?
          pusher(t, slides_by_author, false)
        end

        pusher(t, learning_pathways_by_author, false) if t['layout'] == 'learning-pathway'

        # Philosophies
        has_philosophy[t.data['username']] = true if t['layout'] == 'training_philosophy' && !t.data['username'].nil?
      end

      site.posts.docs.each do |t|
        # News
        pusher(t, news_by_author, true) if t['layout'] == 'news'
      end

      Gtn::Contributors.list(site).reject { |c| c['halloffame'] == 'no' }.each_key do |contributor|
        # Using PageWithoutAFile instead of a custom class which reads files
        # from disk each time, saves some time, but it is unclear how much
        # due to how the previous was accounted. But assuming 0.040s per page * 193 should be about 8 seconds.
        page2 = PageWithoutAFile.new(site, '', File.join(dir, contributor), 'index.html')
        page2.content = nil
        name = Gtn::Contributors.fetch_contributor(site, contributor).fetch('name', contributor)

        # Their tutorials
        page2.data['contributor'] = contributor
        page2.data['personname'] = name
        page2.data['title'] = "GTN Contributor: #{name}"
        page2.data['layout'] = 'contributor_index'

        page2.data['tutorials'] = tutorials_by_author[contributor]
        page2.data['slides'] = slides_by_author[contributor]
        page2.data['news'] = news_by_author[contributor]
        page2.data['learning_pathways'] = learning_pathways_by_author[contributor]

        page2.data['tutorials_count'] = tutorials_by_author[contributor].length
        page2.data['slides_count'] = slides_by_author[contributor].length
        page2.data['news_count'] = news_by_author[contributor].length
        page2.data['learning_pathways_count'] = learning_pathways_by_author[contributor].length

        page2.data['editors'] = TopicFilter.enumerate_topics(site).select do |t|
          t.fetch('editorial_board', []).include?(contributor)
        end
        # Also their learning pathways
        page2.data['editors'] += site.pages.select do |t|
          t['layout'] == 'learning-pathway' && t.data.fetch('editorial_board', []).include?(contributor)
        end
        page2.data['editor_count'] = page2.data['editors'].length

        page2.data['has_philosophy'] = has_philosophy[contributor]

        site.pages << page2
      end
    end
  end
end
