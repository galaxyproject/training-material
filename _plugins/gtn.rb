require './_plugins/gtn/boxify'
require './_plugins/gtn/mod'
require './_plugins/gtn/synthetic'
require './_plugins/gtn/metrics'


module Jekyll
  module GtnFunctions

    def self.cache
      @@cache ||= Jekyll::Cache.new("GtnFunctions")
    end

    ELIXIR_NODES = {
      "au" => "Australia",
      "be" => "Belgium",
      "ch" => "Switzerland",
      "cz" => "Czechia",
      "de" => "Germany",
      "dk" => "Denmark",
      "ee" => "Estonia",
      "es" => "Spain",
      "fi" => "Finland",
      "fr" => "France",
      "gr" => "Greece",
      "hu" => "Hungary",
      "ie" => "Ireland",
      "il" => "Israel",
      "it" => "Italy",
      "lu" => "Luxembourg",
      "nl" => "the Netherlands",
      "no" => "Norway",
      "pt" => "Portugal",
      "se" => "Sweden",
      "si" => "Slovenia",
      "uk" => "United Kingdom",
    }

    def elixirnode2name(name)
      ELIXIR_NODES[name]
    end

    def slugify_unsafe(text)
      # Gets rid of *most* things without making it completely unusable?
      text.gsub(/["'\\\/-;:,.!@#$%^&*()-]/, '').gsub(/\s/, '-')
    end

    def humanize_types(type)
      data = {
        "seq" => "List of Items",
        "str" => "Free Text",
        "map" => "A dictionary/map",
        "float" => "Decimal Number",
        "int" => "Integer Number",
        "bool" => "Boolean"
      }
      data[type]
    end

    def replace_newline_doublespace(text)
      text.gsub(/\n/, "\n  ")
    end

    # These two could be unified tbh
    def last_modified_at(page)
      Gtn::ModificationTimes.obtain_time(page['path'])
    end

    def gtn_mod_date(path)
      # Automatically strips any leading slashes.
      Gtn::ModificationTimes.obtain_time(path.gsub(/^\//, ''))
    end

    def how_many_topic_feedbacks(feedback, name)
      feedback.select{|x| x["topic"] == name}.length
    end

    def how_many_tutorial_feedbacks(feedback, name)
      feedback.select{|x| x["tutorial"] == name}.length
    end

    def filter_authors(contributors, contributions)
      if not contributors.nil?
        return contributors
      else
        return contributions["authorship"]
      end
    end

    def get_default_link(material)
      url = nil

      if material['slides']
        url = "topics/#{material['topic_name']}/tutorials/#{material['tutorial_name']}/slides.html"
      end

      if material['hands_on']
        if material['hands_on'] != "external" && material['hands_on'] != ""
          url = "topics/#{material['topic_name']}/tutorials/#{material['tutorial_name']}/tutorial.html"
        end
      end

      url
    end

  end
end

Liquid::Template.register_filter(Jekyll::GtnFunctions)

if $0 == __FILE__
  result = Gtn::ModificationTimes.obtain_time(ARGV[0].gsub(/^\//, ''))
  puts "Modification time of #{ARGV[0].gsub(/^\//, '')} is #{result}"
end
