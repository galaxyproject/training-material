require './_plugins/jekyll-topic-filter'

module Jekyll
  ##
  # This class generates the GTN's author pags
  class TagPageGenerator < Generator
    safe true

    ##
    # This generates the author pages
    # Params
    # +site+:: The site object
    def generate(site)
      Jekyll.logger.info '[GTN/SyntheticTopics] Generating By-Tag Indexes'
      TopicFilter.list_all_tags(site).map do |tag|
        site.data["by_tag_#{tag}"] = {
          'name' => "by_tag_#{tag}",
          'type' => 'use',
          'title' => tag,
          'summary' => "Tutorials covering #{tag}",
          'tag_based' => true,
          'hidden' => true,
        }

        topic_index = PageWithoutAFile.new(site, '', "tags/#{tag}", 'index.md')
        topic_index.content = ''
        topic_index.data['layout'] = 'topic'
        topic_index.data['topic_name'] = "by_tag_#{tag}"
        topic_index.data['topic'] = site.data["by_tag_#{tag}"]

        site.pages << topic_index
      end

      Jekyll.logger.info '[GTN/SyntheticTopics] Generating By-Tag Embeds'
      TopicFilter.list_all_tags(site).map do |tag|
        topic_index = PageWithoutAFile.new(site, '', "tags/#{tag}", 'embed.html')
        topic_index.content = ''
        topic_index.data['layout'] = 'topic-embed'
        topic_index.data['topic_name'] = "by_tag_#{tag}"
        topic_index.data['topic'] = site.data["by_tag_#{tag}"]

        site.pages << topic_index
      end
    end
  end
end
