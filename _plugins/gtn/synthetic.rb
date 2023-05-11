require './_plugins/jekyll-topic-filter.rb'

module Jekyll
  class SyntheticTopicGenerator < Generator
    def generate(site)
      # Full Bibliography
      puts '[GTN/SyntheticTopics] Generating Indexes'

      TopicFilter.list_topics(site).select { |t| site.data[t]['tag_based'] }.each do |topic|
        puts "[GTN/SyntheticTopics] Creating #{topic} topic"

        topic_index = PageWithoutAFile.new(site, '', "topics/#{topic}", 'index.md')
        topic_index.content = ''
        topic_index.data['layout'] = 'topic'
        topic_index.data['topic_name'] = topic
        site.pages << topic_index

        # For now, intentionally no FAQ
        # faq_index = PageWithoutAFile.new(site, "", "topics/#{topic}/faqs", "index.md")
        # faq_index.content = ""
        # faq_index.data["layout"] = "faq-page"
        # site.pages << faq_index
      end
    end
  end
end
