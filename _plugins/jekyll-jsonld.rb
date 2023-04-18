require 'json'
require './_plugins/gtn.rb'

module Jekyll
  module JsonldFilter
    GTN = {
      "@type": "Organization",
      "email": "galaxytrainingnetwork@gmail.com",
      "name": "Galaxy Training Network",
      "url": "https://galaxyproject.org/teach/gtn/"
    }

    A11Y = {
      "accessMode": ["textual", "visual"],
      "accessModeSufficient": ["textual", "visual"],
      #"accessibilityAPI": ,
      "accessibilityControl": ["fullKeyboardControl", "fullMouseControl"],
      "accessibilityFeature": ["alternativeText", "tableOfContents"],
      #"accessibilityHazard": [],
      "accessibilitySummary": "The text aims to be as accessible as possible. Image descriptions will vary per tutorial, from images being completely inaccessible, to images with good descriptions for non-visual users.",
    }

    ##
    # Generate the Dublin Core metadata for a material.
    # Parmaeters:
    # +material+:: The material to generate the metadata for.
    # +site+:: The site object.
    # Returns:
    # A string containing the metadata.
    #
    # Example:
    #  {{ material | generate_dublin_core: site }}
    #  => <meta name="DC.identifier" content="..." />
    def generate_dublin_core(material, site)
      if material.key?('data') && material['data'].fetch('type', 'none') != "tutorial_hands_on"
        return
      end

      attributes = [
        ["DC.identifier", site['github_repository']],
        ["DC.type", "text"],
        ["DC.title", material['title']],
        ["DC.publisher", "Galaxy Training Network"],
        ["DC.date", Gtn::ModificationTimes.obtain_time(material['path'])],
      ]

      attributes += get_authors(material).map{|user|
        if site['data']['contributors'].has_key?(user) then
          ['DC.creator', site['data']['contributors'][user].fetch('name', user)]
        else
          puts "[GTN/Meta] #{user} not found in CONTRIBUTORS.yaml"
          ['DC.creator', user]
        end
      }

      return attributes.map{|a, b| "<meta name=\"#{a}\" content=\"#{b}\">" }.join("\n")
    end

    ##
    # Get the authors of a material.
    # Time number 4 i've seen this
    # TODO(hexylena): make this a function in gtn.rb
    # Parmeters:
    # +material+:: The material to get the authors for.
    # Returns:
    # An array of authors.
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
    # Generate the JSON-LD metadata for a person
    # Parameters:
    # +id+:: The id of the person.
    # +contributor+:: The contributor object from CONTRIBUTORS.yaml.
    # +site+:: The site object.
    # Returns:
    # +Hash+:: The JSON-LD metadata.
    #
    # Example:
    #  generate_person_jsonld("hexylena", site['data']['contributors']['hexylena'], site)
    #  => {
    #    "@context": "https://schema.org",
    #    "@type": "Person",
    #    "http://purl.org/dc/terms/conformsTo": {
    #      # Bioschemas profile
    #      "@id": "https://bioschemas.org/profiles/Person/0.2-DRAFT-2019_07_19",
    #      "@type": "Person"
    #    },
    #    "url": "https://training.galaxyproject.org/hall-of-fame/hexylena/",
    #    "mainEntityOfPage": "https://training.galaxyproject.org/hall-of-fame/hexylena/",
    #    "name": "hexylena",
    #    "image": "https://avatars.githubusercontent.com/hexylena",
    #    "description": "A contributor to the GTN project.",
    #    "memberOf": [...],
    #    "identifier": "https://orcid.org/0000-0002-6601-2165",
    #    "orcid": "https://orcid.org/0000-0002-6601-2165"
    #  }
    #
    def generate_person_jsonld(id, contributor, site)
      person = {
        "@context": "https://schema.org",
        "@type": "Person",
        "http://purl.org/dc/terms/conformsTo": {
            "@id": "https://bioschemas.org/profiles/Person/0.2-DRAFT-2019_07_19",
            "@type": "Person"
        },
        # I guess these are identical?
        "url": "#{site['url']}#{site['baseurl']}/hall-of-fame/#{id}/",
        "mainEntityOfPage": "#{site['url']}#{site['baseurl']}/hall-of-fame/#{id}/",
        "name": contributor.nil? ? id : contributor.fetch('name', id),
        "image": "https://avatars.githubusercontent.com/#{id}",
        # No clue what to put here it's a person.
        "description": contributor.nil? ? "A contributor to the GTN project." : contributor.fetch("bio", "A contributor to the GTN project."),
        "memberOf": [GTN],
      }
      if ! contributor.nil? && contributor.has_key?('orcid')
        person['identifier'] = "https://orcid.org/" + contributor['orcid']
        person['orcid'] = "https://orcid.org/" + contributor['orcid']
      end

      person
    end

    ##
    # Generate the JSON-LD metadata for a person as JSON.
    # Parameters:
    # +id+:: The id of the person.
    # +contributor+:: The contributor object from CONTRIBUTORS.yaml.
    # +site+:: The site object.
    # Returns:
    # +String+:: The JSON-LD metadata.
    def to_person_jsonld(id, contributor, site)
      JSON.pretty_generate(generate_person_jsonld(id, contributor, site))
    end

    ##
    # Generate the JSON-LD metadata for a news article (blog)
    # Parameters:
    # +page+:: The page object.
    # +site+:: The +Jekyll::Site+ site object.
    # Returns:
    # +Hash+:: The JSON-LD metadata.
    def generate_news_jsonld(page, site)
      authors = get_authors(page.to_h).map{ |x| generate_person_jsonld(x, site['data']['contributors'][x], site) }

      data = {
        "@context": "https://schema.org",
        "@type": "BlogPosting",
        "url": "#{site['url']}#{site['baseurl']}#{page['url']}",
        "name": page['title'],
        "headline": page.excerpt[0..100].gsub(/\n/, ' '), # todo remove html tags.
        "keywords": page.fetch('tags', []),
        "description": page.excerpt[0..100].gsub(/\n/, ' '), # todo remove html tags
        "articleBody": page.content, # todo remove html tags
        "datePublished": page.date,
        "dateModified": Gtn::ModificationTimes.obtain_time(page.path),
        "author": authors,
        "publisher": GTN,
        "mainEntityOfPage": {
          "@type": "WebPage",
          "@id": "#{site['url']}#{page['url']}"
        },
        "image": {
          "@type": "ImageObject",
          "width": 60,
          "height": 60,
          "url": "#{site['baseurl']}/assets/images/GTN-60px.png"
        }
      }
      data.update(A11Y)

      JSON.pretty_generate(data)
    end

    ##
    # Convert a material to JSON-LD.
    # Parameters:
    # +material+:: The material object.
    # +topic+:: The topic object.
    # +site+:: The +Jekyll::Site+ site object.
    #
    # Returns:
    # +String+:: The JSON-LD metadata.
    def to_jsonld(material, topic, site)
      langCodeMap = {
        'en': "English",
        'es': "Español",
      }


      eduLevel = {
        "Introductory" => "Beginner",
        "Intermediate" => "Intermediate",
        "Advanced"     => "Advanced"
      }
      if not topic then
        return '{}'
      end

      topic_desc = {
        "@type": "CreativeWork",
        "name": "#{topic['title']}",
        "description": "#{topic['summary']}",
        "url": "#{site['url']}#{site['baseurl']}/topics/#{topic['name']}/"
      }

      # aggregate everything
      data = {
        # Properties from Course
        "@context": "http://schema.org",
        "@type": "LearningResource",

        # Required for BioSchemas
        "http://purl.org/dc/terms/conformsTo": {
            "@id": "https://bioschemas.org/profiles/TrainingMaterial/1.0-RELEASE",
            "@type": "CreativeWork"
        },

        # Properties from CreativeWork
        #"about" described below
        #
        #"accountablePerson":,
        #"aggregateRating":,
        #"alternativeHeadline":,
        #"associatedMedia":,
        "audience": {
          "@type": "EducationalAudience",
          "educationalRole": "Students"
        },
        #"audio":,
        #"award":,
        #"author" described below
        #"character":,
        "citation": {
          "@type": "CreativeWork",
          "name": "Community-Driven Data Analysis Training for Biology",
          "url": "https://doi.org/10.1016/j.cels.2018.05.012"
        },
        #"comment":,
        #"commentCount":,
        #"contentLocation":,
        #"contentRating":,
        #"contentReferenceTime":,
        #"contributor" described below
        "copyrightHolder": GTN,
        #"copyrightYear":,
        #"correction":,
        #"creator":,
        #"dateCreated":,
        "dateModified": Gtn::ModificationTimes.obtain_time(material['path']),
        #"datePublished":,
        "discussionUrl": site["gitter_url"],
        #"editor":,
        #"educationalAlignment":,
        #"educationalUse":,
        #"encoding":,
        #"encodingFormat":,
        #"exampleOfWork":,
        #"expires":,
        #"funder":,
        #"genre":,
        #"hasPart" described below
        "headline": "#{material['title']}",
        #"interactionStatistic":,
        "interactivityType": "mixed",
        "isAccessibleForFree": true,
        #"isBasedOn":,
        "isFamilyFriendly": true,
        #"isPartOf" described below
        #"keywords": described below
        #"learningResourceType" described below
        "license": "https://spdx.org/licenses/CC-BY-4.0.html",
        #"locationCreated":,
        #"mainEntity":,
        #"material":,
        #"mentions" described below
        #"offers":,
        #"position":,
        "producer": GTN,
        "provider": GTN,
        #"publication":,
        #"publisher":,
        #"publisherImprint":,
        #"publishingPrinciples":,
        #"recordedAt":,
        #"releasedEvent":,
        #"review":,
        #"schemaVersion":,
        #"sdDatePublished":,
        #"sdLicense":,
        #"sdPublisher":,
        "sourceOrganization": GTN,
        #"spatialCoverage":,
        #"sponsor":,
        #"temporalCoverage":,
        #"text":,
        #"thumbnailUrl":,
        #"timeRequired" described below
        #"translationOfWork":,
        #"translator": Google Translate???,
        #"typicalAgeRange":,
        #"version":,
        #"video":,
        #"workExample":,
        #"workTranslation":,

        # Properties from Thing
        #"additionalType":,
        #"alternateName":,
        #"description" described below
        #"disambiguatingDescription":,
        "identifier": site['github_repository'],
        #"image":,
        #"mainEntityOfPage":,
        #"name" described below
        #"potentialAction":,
        #"sameAs":,
        #"subjectOf":,
        # "url" described below
      }
      data.update(A11Y)

      #info depending if tutorial, hands-on or slide level
      parts = []
      #data['hasPart'] = parts

      mentions = []
      description = []

      data['isPartOf'] = topic_desc

      if material['name'] == 'tutorial.md' or material['name'] == 'slides.html' then
        if material['name'] == 'tutorial.md' then
          data['learningResourceType'] = "hands-on tutorial"
          data['name'] = "Hands-on for '#{material['title']}' tutorial"
        else
          data['learningResourceType'] = "slides"
          data['name'] = "Slides for '#{material['title']}' tutorial"
        end
        data['url'] = "#{site['url']}#{site['baseurl']}#{material['url']}"

        # Time required
        if material.key?('time_estimation') and not material['time_estimation'].nil? then
          data['timeRequired'] = "PT#{material['time_estimation'].upcase}"
        end

        # Description with questions, objectives and keypoints
        if material.key?('questions') and not material['questions'].nil? and material['questions'].length > 0 then
          questions = material['questions'].join("\n - ")
          description.push("The questions this #{material['type']} addresses are:\n - #{questions}\n\n")
        end
        if material.key?('objectives') and not material['objectives'].nil? and material['objectives'].length > 0 then
          objectives = material['objectives'].join("\n - ")
          description.push("The objectives are:\n - #{objectives}\n\n")
        end
        if material.key?('keypoints') and not material['keypoints'].nil? and material['keypoints'].length > 0 then
          keypoints = material['keypoints'].join("\n - ")
          description.push("The keypoints are:\n - #{keypoints}\n\n")
        end

        # Keywords
        data['keywords'] = [topic['name']] + material.fetch('tags', [])
        data['keywords'] = data['keywords'].join(', ')
        #Zenodo links
        if material.key?('zenodo_link') then
          mentions = mentions.push({
            "@type": "Thing",
            "url": "#{material['zenodo_link']}",
            "name": "Training data for #{material['title']} tutorial"
          })
        end
      end
      data['description'] = description.join("\n")

      if material.key?("lang") then
        data['inLanguage'] = {
          "@type": "Language",
          "name": langCodeMap[material['lang']],
          "alternateName": material['lang']
        }
      else
        data['inLanguage'] = {
          "@type": "Language",
          "name": "English",
          "alternateName": "en"
        }
      end

      # Course requirements (material + topic)
      reqs = []
      if topic.key?('requirements') then
        reqs.push(*topic['requirements'])
      end
      if material.key?('requirements') then
        reqs.push(*material['requirements'])
      end
      if !reqs.empty?
        coursePrerequisites = []
        for req in reqs do
          if req['type'] == "internal" then
            if req.key?('tutorials') then
              for tuto in req['tutorials'] do
                for page in site['pages'] do
                  if page['name'] == 'tutorial.md' or page['name'] == 'slides.html' then
                    if page['topic_name'] == req['topic_name'] and page['tutorial_name'] == tuto then
                      #slides
                      if page['name'] == 'slides.html' then
                        coursePrerequisites.push({
                          "@context": "http://schema.org",
                          "@type": "LearningResource",
                          "url": "#{site['url']}#{site['baseurl']}/topics/#{req['topic_name']}/tutorials/#{tuto}/slides.html",
                          "name": "#{page['title']}",
                          "description": "Slides for '#{page['title']}' tutorial",
                          "learningResourceType": "slides",
                          "interactivityType": "expositive",
                          "provider": GTN
                        })
                        if page['hands_on_url'] then
                          coursePrerequisites.push({
                            "@context": "http://schema.org",
                            "@type": "LearningResource",
                            "url": "#{page['hands_on_url']}",
                            "learningResourceType": "hands-on tutorial",
                            "interactivityType": "expositive",
                          })
                        end
                      end
                      #hands-on
                      if page['name'] == 'tutorial.md' then
                        coursePrerequisites.push({
                          "@context": "http://schema.org",
                          "@type": "LearningResource",
                          "url": "#{site['url']}#{site['baseurl']}/topics/#{req['topic_name']}/tutorials/#{tuto}/tutorial.html",
                          "name": "#{page['title']}",
                          "description": "Hands-on for '#{page['title']}' tutorial",
                          "learningResourceType": "hands-on tutorial",
                          "interactivityType": "expositive",
                          "provider": GTN
                        })
                      end
                    end
                  end
                end
              end
            else
              coursePrerequisites.push({
                "@context": "http://schema.org",
                "@type": "LearningResource",
                "url": "#{site['url']}#{site['baseurl']}/topics/#{req['topic_name']}/",
                "name": "#{site['data'][req['topic_name']]['title']}",
                "description": "#{site['data'][req['topic_name']]['title']}",
                "provider": GTN
              })
            end
          elsif req['type'] == "external" then
            coursePrerequisites.push({
              "@type": "CreativeWork",
              "url": "#{req['link']}",
              "name": "#{req['title']}"
            })
          else
            coursePrerequisites.push("#{req['title']}")
          end
        end
        data['competencyRequired'] = coursePrerequisites.uniq
      end

      # Add contributors/authors
      if material.key?('contributors') then
        contributors = material['contributors'].map{ |x| generate_person_jsonld(x, site['data']['contributors'][x], site) }
        data['author'] = contributors
        data['contributor'] = contributors
      end

      about = []
      about.push(topic_desc)
      if topic.key?('edam_ontology') then
        about.push({
          "@type": "DefinedTerm",
          "@id": "http://edamontology.org/#{topic['edam_ontology']}",
          "inDefinedTermSet": "http://edamontology.org",
          "termCode": "#{topic['edam_ontology']}",
          #"name": ,
          "url": "https://bioportal.bioontology.org/ontologies/EDAM/?p=classes&conceptid=http%3A%2F%2Fedamontology.org%2F#{topic['edam_ontology']}"
        })
      end
      data['about'] = about

      if material.key?('level') then
        data['educationalLevel'] = eduLevel[material['level']]
      end

      return JSON.pretty_generate(data)
    end
  end
end

Liquid::Template.register_filter(Jekyll::JsonldFilter)
