# frozen_string_literal: true

require 'json'
require 'liquid'
require './_plugins/colour-tags'
require './_plugins/jekyll-topic-filter'

module Jekyll
  # Export search data as JSON
  class DumpSearchDataTag < Liquid::Tag
    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def getlist(tutorial, attr)
      tutorial[attr] || []
    end

    def render(context)
      puts '[GTN/Search]'

      site = context.registers[:site]
      topics = TopicFilter.list_topics_h(site)

      results = {}
      topics.each do |k, topic|
        tutorials = site.data['cache_topic_filter'][k]
        tutorials.each do |tutorial|
          results[tutorial['url']] = {
            'type' => 'Tutorial',
            'topic' => topic['title'],
            'title' => tutorial['title'],
            'contributors' => getlist(tutorial, 'contributors').map do |c|
              site.data['contributors'].fetch(c, {}).fetch('name', c)
            end.join(', '),
            'tags' => getlist(tutorial, 'tags').map do |tag|
              href = "#{site.baseurl}/search?query=#{tag}"
              title = "Show all tutorials tagged #{tag}"
              style = ColourTag.colour_tag tag
              %(<a class="label label-default" title="#{title}" href="#{href}" style="#{style}">#{tag}</a>)
            end,
            'url' => site.baseurl + tutorial['url'],
          }
        end
      end

      faqs = site.pages.select { |p| p.data['layout'] == 'faq' }
      faqs.each do |resource|
        results[resource['url']] = {
          'type' => 'FAQ',
          'topic' => 'FAQ',
          'title' => resource['title'],
          'contributors' => getlist(resource.data, 'contributors').map do |c|
            site.data['contributors'].fetch(c, {}).fetch('name', c)
          end.join(', '),
          'tags' => [],
          'url' => site.baseurl + resource['url'],
        }
      end

      JSON.pretty_generate(results)
    end
  end
end

Liquid::Template.register_tag('dump_search_view', Jekyll::DumpSearchDataTag)
