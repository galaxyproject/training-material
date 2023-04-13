require './_plugins/gtn/boxify'
require './_plugins/gtn/mod'
require './_plugins/gtn/synthetic'
require './_plugins/gtn/metrics'
require './_plugins/gtn/scholar'


puts "[GTN] You are running #{RUBY_VERSION} released on #{RUBY_RELEASE_DATE} for #{RUBY_PLATFORM}"
version_parts = RUBY_VERSION.split(".")
if version_parts[0].to_i < 3
  puts "[GTN] WARNING: This Ruby is pretty old, you might want to update."
end


## 
# This module contains functions that are used in the GTN, our internal functions that is.

##
# This function returns the authors of a material, if it has any. (via contributors or contribution)
# Params:
# +material+:: The material to get the authors of
# Returns:
# +Array+:: The authors of the material
def get_authors(material)
  if material.key?('contributors') then
    material['contributors']
  elsif material.key?('contributions') then
    material['contributions']['authorship']
  else
    []
  end
end

##
# This function returns the name of a user, if it is known. Otherwise, it returns the user name.
# Params:
# +user+:: The user to get the name of
# +site+:: The +Jekyll::Site+ object
# Returns:
# +String+:: The name of the user
def lookup_name(user, site)
  if site.data['contributors'].has_key?(user) then
    site.data['contributors'][user].fetch('name', user)
  else
    user
  end
end

module Jekyll
  module GtnFunctions

    def self.cache
      @@cache ||= Jekyll::Cache.new("GtnFunctions")
    end

    ##
    # List of elixir node country IDs (ISO 3166-1 alpha-2) and their names
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

    ##
    # Returns the name of an elixir node, given its country ID
    # Params:
    # +name+:: The country ID of the node (ISO 3166-1 alpha-2)
    # Returns:
    # +String+:: The name of the node
    def elixirnode2name(name)
      ELIXIR_NODES[name]
    end

    ##
    # Obtain the most cited paper in the GTN
    # Params:
    # +citations+:: The citations to search through
    #
    # Returns:
    # +Hash+:: The papers including their text citation and citation count
    def top_citations(citations)
      if citations.nil?
        {}
      else
        citations.sort_by{|k, v| v}.reverse.to_h.first(20).map{|k, v| 
          [k, {"count" => v, "text" => Gtn::Scholar.render_citation(k)}]
        }.to_h
      end
    end

    ##
    # A slightly more unsafe slugify function
    # Params:
    # +text+:: The text to slugify
    # Returns:
    # +String+:: The slugified text
    #
    # Example:
    #  slugify_unsafe("Hello, World!") # => "Hello-World"
    def slugify_unsafe(text)
      # Gets rid of *most* things without making it completely unusable?
      text.gsub(/["'\\\/-;:,.!@#$%^&*()-]/, '').gsub(/\s/, '-')
    end

    ##
    # Return human text for ruby types
    # Params:
    # +type+:: The type to humanize
    # Returns:
    # +String+:: The humanized type
    #
    # Example:
    #  humanize_types("seq") # => "List of Items"
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

    ##
    # Replaces newlines with newline + two spaces
    def replace_newline_doublespace(text)
      text.gsub(/\n/, "\n  ")
    end

    ##
    # Returns the last modified date of a page
    # Params:
    # +page+:: The page to get the last modified date of
    # Returns:
    # +String+:: The last modified date of the page
    #
    # TODO: These two could be unified tbh
    def last_modified_at(page)
      Gtn::ModificationTimes.obtain_time(page['path'])
    end

    ##
    # Returns the last modified date of a page
    # Params:
    # +page+:: The page to get the last modified date of
    # Returns:
    # +String+:: The last modified date of the page
    #
    def gtn_mod_date(path)
      # Automatically strips any leading slashes.
      Gtn::ModificationTimes.obtain_time(path.gsub(/^\//, ''))
    end

    ##
    # How many times has a topic been mentioned in feedback?
    # Params:
    # +feedback+:: The feedback to search through
    # +name+:: The name of the topic to search for
    # Returns:
    # +Integer+:: The number of times the topic has been mentioned
    def how_many_topic_feedbacks(feedback, name)
      feedback.select{|x| x["topic"] == name}.length
    end

    ## 
    # How many times has a tutorial been mentioned in feedback?
    # Params:
    # +feedback+:: The feedback to search through
    # +name+:: The name of the tutorial to search for
    # Returns:
    # +Integer+:: The number of times the tutorial has been mentioned
    def how_many_tutorial_feedbacks(feedback, name)
      feedback.select{|x| x["tutorial"] == name}.length
    end

    ##
    # Fix the titles of boxes in a page
    # Params:
    # +content+:: The content to fix
    # +lang+:: The language of the content
    # +key+:: The key of the content
    # Returns:
    # +String+:: The fixed content
    def fix_box_titles(content, lang, key)
      Gtn::Boxify.replace_elements(content, lang, key)
    end

    ##
    # Basically a dupe of 'get_authors'
    # Params:
    # +contributors+:: The contributors to the material
    # +contributions+:: The contributions to the material
    # Returns:
    # +Array+:: The "authors" of the material
    #
    # TODO(hexylena) de-duplicate
    #
    # Example:
    #  {% assign authors = page.contributors | filter_authors:page.contributions -%}
    def filter_authors(contributors, contributions)
      if not contributors.nil?
        return contributors
      else
        return contributions["authorship"]
      end
    end

    ##
    # Convert a fedi address to a link
    # Params:
    # +fedi_address+:: The fedi address to convert
    # Returns:
    # +String+:: The URL at which their profile is accessible
    #
    # Example:
    #  {{ contributors[page.contributor].fediverse | fedi2link }}
    #
    #  fedi2link("@hexylena@galaxians.garden") => "https://galaxians.garden/@hexylena"
    def fedi2link(fedi_address)
      fedi_address.gsub(/^@?(?<user>.*)@(?<host>.*)$/){|m| "https://#{$~[:host]}/@#{$~[:user]}" }
    end

    ##
    # Load an SVG file directly into the page
    # Params:
    # +path+:: The path of the SVG file (relative to GTN workspace root)
    # Returns:
    # +String+:: The SVG file contents
    #
    # Example:
    #  {{ "assets/images/mastodon.svg" | load_svg }}
    def load_svg(path)
      File.open(path).read.gsub(/\R+/, '')
    end

    def regex_replace(str, regex_search, value_replace)
      regex = /#{regex_search}/m
      return str.gsub(regex, value_replace)
    end

    def regex_replace_once(str, regex_search, value_replace)
      regex = /#{regex_search}/m
      return str.sub(regex, value_replace)
    end

    ##
    # Convert a workflow path to a TRS path
    # Params:
    # +str+:: The workflow path
    # Returns:
    # +String+:: The TRS path
    #
    # Example:
    #  {{ "topics/metagenomics/tutorials/mothur-miseq-sop-short/workflows/workflow1_quality_control.ga" | convert_workflow_path_to_trs }}
    #  => "/api/ga4gh/trs/v2/tools/metagenomics-mothur-miseq-sop-short/versions/workflow1_quality_control"
    def convert_workflow_path_to_trs(str)
      if str.nil?
        return "GTN_TRS_ERROR_NIL"
      end

      m = str.match(/topics\/(?<topic>.*)\/tutorials\/(?<tutorial>.*)\/workflows\/(?<workflow>.*)\.ga/)
      if m
        puts "str=#{str} m=#{m} #{m[:topic]} #{m[:tutorial]} #{m[:workflow]}"
        return "/api/ga4gh/trs/v2/tools/#{m[:topic]}-#{m[:tutorial]}/versions/#{m[:workflow]}"
      end
      return "GTN_TRS_ERROR"
    end

    ##
    # Get the topic of a page's path
    # Params:
    # +page+:: The page to get the topic of, it will inspect page['path']
    # Returns:
    # +String+:: The topic of the page
    #
    # Example:
    #  {{ page | get_topic }}
    def get_topic(page)
      # Arrays that will store all introduction slides and tutorials we discover.
      page['path'].split('/')[1]
    end

    ##
    # Gets the 'default' link for a material, hands on if it exists, otherwise slides.
    # Params:
    # +material+:: The material to get the link for
    # Returns:
    # +String+:: The URL of the default link
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


##
# This does post-modification to every page
# Mapping the authors to their human names, and copying the cover (when present) to 'image'
#
# This exists because the jekyll-feed plugin expects those fields to look like that.
Jekyll::Hooks.register :posts, :pre_render do |post, out|
  post.data['author'] = get_authors(post.data).map{|c| lookup_name(c, post.site)}.join(", ")
  post.data['image'] = post.data['cover']
end

if $0 == __FILE__
  result = Gtn::ModificationTimes.obtain_time(ARGV[0].gsub(/^\//, ''))
  puts "Modification time of #{ARGV[0].gsub(/^\//, '')} is #{result}"
end
