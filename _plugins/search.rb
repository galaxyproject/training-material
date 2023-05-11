# frozen_string_literal: true

require 'json'
require 'liquid'
require './_plugins/colour-tags'

module Jekyll
  class DumpSearchDataTag < Liquid::Tag
    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def getlist(tutorial, attr)
      [] if tutorial[attr].nil?

      tutorial.fetch(attr, [])
    end

    def render(context)
      puts '[GTN/Search]'

      site = context.registers[:site]
      topics = site.data.select { |_k, v| v.is_a?(Hash) && v.key?('type') }
                   .select { |_k, v| %w[use admin-dev basics].include? v['type'] }

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
              %(<a class="label label-default" title="Show all tutorials tagged #{tag}" href="#{site.baseurl}/search?query=#{tag}" style="#{ColourTag.colour_tag tag}">#{tag}</a>)
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
