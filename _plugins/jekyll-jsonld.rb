# frozen_string_literal: true

require 'json'
require './_plugins/gtn'
require './_plugins/gtn/git'
require './_plugins/util'

module Jekyll
  # Generate JSON-LD metadata for the GTN.
  module JsonldFilter
    GTN = {
      '@type': 'Organization',
      'http://purl.org/dc/terms/conformsTo': {
        # Bioschemas profile
        '@id': 'https://bioschemas.org/profiles/Organization/0.2-DRAFT-2019_07_19',
        '@type': 'Organization'
      },
      id: 'https://training.galaxyproject.org',
      email: 'galaxytrainingnetwork@gmail.com',
      name: 'Galaxy Training Network',
      legalName: 'Galaxy Training Network',
      alternateName: 'GTN',
      url: 'https://training.galaxyproject.org',
      logo: 'https://training.galaxyproject.org/training-material/assets/images/GTNLogo1000.png',
      fundingModel: "The GTN's infrastructure relies on GitHub and the Galaxy Project for hosting costs. " \
                    'There are no full time paid staff members of the GTN. Individuals are occasionally funded on ' \
                    'GTN-adjacent projects.',
      keywords: %w[galaxy bioinformatics training fair accessible],
      status: 'active',
      foundingDate: Gtn::Git.discover['founding_date'].to_s,
      socialMedia: 'https://mstdn.science/@gtn',
      type: 'project',
    }.freeze

    A11Y = {
      accessMode: %w[textual visual],
      accessModeSufficient: %w[textual visual],
      # "accessibilityAPI": ,
      accessibilityControl: %w[fullKeyboardControl fullMouseControl],
      accessibilityFeature: %w[alternativeText tableOfContents],
      # "accessibilityHazard": [],
      accessibilitySummary: 'The text aims to be as accessible as possible. Image descriptions will vary per ' \
                            'tutorial, from images being completely inaccessible, to images with good descriptions ' \
                            'for non-visual users.',
    }.freeze

    EDU_ROLES = {
      'use' => 'Students',
      'admin-dev' => 'Galaxy Administrators',
      'basics' => 'Students',
      'data-science' => 'Data-Science Students',
      'instructors' => 'Instructors',
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
      return if material.key?('data') && material['data'].fetch('type', 'none') != 'tutorial_hands_on'

      attributes = [
        ['DC.identifier', site['github_repository']],
        ['DC.type', 'text'],
        ['DC.title', material['title']],
        ['DC.publisher', 'Galaxy Training Network'],
        ['DC.date', Gtn::ModificationTimes.obtain_time(material['path'])]
      ]

      attributes += Gtn::Contributors.get_authors(material).map do |user|
        ['DC.creator', Gtn::Contributors.fetch_name(site, user)]
      end

      attributes.map { |a, b| "<meta name=\"#{a}\" content=\"#{b}\">" }.join("\n")
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
      member_of = Gtn::Contributors.fetch_contributor(site, id)['affiliations'] || []
      member_of = member_of.map do |org_id|
        org = Gtn::Contributors.fetch_contributor(site, org_id)
        generate_org_jsonld(org_id, org, site)
      end

      person = {
        '@context': 'https://schema.org',
        '@type': 'Person',
        'http://purl.org/dc/terms/conformsTo': {
          '@id': 'https://bioschemas.org/profiles/Person/0.3-DRAFT',
          '@type': 'CreativeWork'
        },
        # I guess these are identical?
        url: "#{site['url']}#{site['baseurl']}/hall-of-fame/#{id}/",
        mainEntityOfPage: "#{site['url']}#{site['baseurl']}/hall-of-fame/#{id}/",
        name: Gtn::Contributors.fetch_name(site, id),
        image: "https://avatars.githubusercontent.com/#{id}",
        # No clue what to put here it's a person.
        description: if contributor.nil?
                       'A contributor to the GTN project.'
                     else
                       contributor.fetch('bio',
                                         'A contributor to the GTN project.')
                     end,
        memberOf: [GTN] + member_of,
      }
      if !contributor.nil? && contributor.key?('orcid') && contributor['orcid']
        person['identifier'] = "https://orcid.org/#{contributor['orcid']}"
        person['orcid'] = "https://orcid.org/#{contributor['orcid']}"
      end

      person
    end

    def generate_org_jsonld(id, contributor, site)
      organization = {
        '@context': 'https://schema.org',
        '@type': 'Organization',
        'http://purl.org/dc/terms/conformsTo': {
          '@id': 'https://bioschemas.org/profiles/Organization/0.3-DRAFT',
          '@type': 'CreativeWork'
        },
        id: "#{site['url']}#{site['baseurl']}/hall-of-fame/#{id}/",
        name: Gtn::Contributors.fetch_name(site, id),
        description: 'An organization supporting the Galaxy Training Network',
      }

      organization['url'] = contributor['url'] if contributor.key?('url') && contributor['url']

      organization
    end

    def generate_funder_jsonld(id, contributor, site)
      {
        '@context': 'https://schema.org',
        '@type': 'Organization',
        'http://purl.org/dc/terms/conformsTo': {
          '@id': 'https://bioschemas.org/profiles/Organization/0.3-DRAFT',
          '@type': 'CreativeWork'
        },
        name: Gtn::Contributors.fetch_name(site, id),
        description: contributor.fetch('funding_statement', 'An organization supporting the Galaxy Training Network'),
        url: contributor.fetch('url', "https://training.galaxyproject.org/training-material/hall-of-fame/#{id}/"),
        logo: contributor.fetch('avatar', "https://github.com/#{id}.png"),
      }
    end

    def generate_grant_jsonld(id, contributor, site)
      organization = {
        '@context': 'https://schema.org',
        '@type': 'Grant',
        identifier: contributor['funding_id'],
        url: Gtn::Contributors.fetch_funding_url(contributor) || contributor['url'],
        funder: generate_funder_jsonld(id, contributor, site)
      }

      organization['startDate'] = contributor['start_date'] if contributor.key?('start_date')
      organization['endDate'] = contributor['end_date'] if contributor.key?('end_date')

      organization
    end

    ##
    # Generate the JSON-LD metadata for a person, funder, or organisation as JSON.
    # Parameters:
    # +id+:: The id of the person.
    # +site+:: The site object.
    # +json+:: Should the output be rendered as JSON (only really used in contributor page.)
    # Returns:
    # +String+:: The JSON-LD metadata.
    def to_pfo_jsonld(id, site, json: true)
      contributor = Gtn::Contributors.fetch_contributor(site, id)
      d = if Gtn::Contributors.person?(site, id)
            generate_person_jsonld(id, contributor, site)
          elsif Gtn::Contributors.grant?(site, id)
            generate_grant_jsonld(id, contributor, site)
          else
            generate_org_jsonld(id, contributor, site)
          end

      if json
        JSON.pretty_generate(d)
      else
        d
      end
    end

    ##
    # Generate the JSON-LD metadata for a news article (blog)
    # Parameters:
    # +page+:: The page object.
    # +site+:: The +Jekyll::Site+ site object.
    # Returns:
    # +Hash+:: The JSON-LD metadata.
    def generate_news_jsonld(page, site)
      authors = Gtn::Contributors.get_authors(page.to_h).map do |x|
        to_pfo_jsonld(x, site, json: false)
      end

      data = {
        '@context': 'https://schema.org',
        '@type': 'BlogPosting',
        url: "#{site['url']}#{site['baseurl']}#{page['url']}",
        name: page['title'],
        headline: page.excerpt[0..100].gsub(/\n/, ' '), # TODO: remove html tags.
        keywords: page['tags'] || [],
        description: page.excerpt[0..100].gsub(/\n/, ' '), # TODO: remove html tags
        articleBody: page.content, # TODO: remove html tags
        datePublished: page.date,
        dateModified: Gtn::ModificationTimes.obtain_time(page.path),
        author: authors,
        publisher: GTN,
        mainEntityOfPage: {
          '@type': 'WebPage',
          '@id': "#{site['url']}#{page['url']}"
        },
        image: {
          '@type': 'ImageObject',
          width: 60,
          height: 60,
          url: "#{site['baseurl']}/assets/images/GTN-60px.png"
        }
      }
      data.update(A11Y)

      JSON.pretty_generate(data)
    end

    ##
    # Generate the JSON-LD metadata for an event
    # Parameters:
    # +page+:: The page object.
    # +site+:: The +Jekyll::Site+ site object.
    # Returns:
    # +Hash+:: The JSON-LD metadata.
    def generate_event_jsonld(page, site)
      organisers = Gtn::Contributors.get_organisers(page.to_h).map do |x|
        to_pfo_jsonld(x, site, json: false)
      end
      instructors = Gtn::Contributors.get_instructors(page.to_h).map do |x|
        to_pfo_jsonld(x, site, json: false)
      end
      funders = Gtn::Contributors.get_funders(site, page.to_h).map do |x|
        to_pfo_jsonld(x, site, json: false)
      end
      funding = Gtn::Contributors.get_grants(site, page.to_h).map do |x|
        to_pfo_jsonld(x, site, json: false)
      end

      materials = []
      if page['program']
        page['program'].each do |section|
          if !section.key? 'tutorials'
            next
          end

          section['tutorials'].each do |tutorial|
            if tutorial.key?('custom')
              next
            end

            material = TopicFilter.fetch_tutorial_material(site, tutorial['topic'], tutorial['name'])
            materials.push(material)
          end
        end
      end
      materials.compact!

      # Extract EDAM terms from all materials
      edam_terms = materials.map do |material|
        material.fetch('edam_ontology', []).map do |term|
          {
            '@type': 'DefinedTerm',
            '@id': "http://edamontology.org/#{term}",
            inDefinedTermSet: 'http://edamontology.org',
            termCode: term,
          }
        end
      end.flatten.uniq

      learning_objectives = materials.map do |material|
        material.fetch('objectives', [])
      end.flatten.compact

      # TODO: add topic edam terms too? Not sure.
      parts = []
      materials.each do |material|
        mat = generate_material_jsonld(material, site['data'][material['topic_name']], site)
        if !mat.nil? && !mat.empty?
          parts.push(mat)
        end
      end

      if page['program']
        syllab = page['program'].reject { |s| s['section'].nil? }.map do |section|
          {
            '@type': 'Syllabus',
            name: section['section'],
            description: section.fetch('description', nil),
          }
        end
      end

      data = {
        '@context': 'https://schema.org',
        '@type': 'Course',
        url: "#{site['url']}#{site['baseurl']}#{page['url']}",
        name: page['title'],
        keywords: page['tags'] || [],
        description: page['description'],

        about: edam_terms, # TeSS, "scientific topics".
        audience: page['audience'], # TeSS: target audience
        # If 'online' is present in the mode, the course is online.
        # Will fail on "this is NOT an online course"
        # Acceptable.
        courseMode: page['mode'],
        startDate: page['date_start'],
        endDate: page['date_end'],
        organizer: organisers, # TeSS only, US spelling, non-standard

        location: page['location'], # TODO, TeSS location
        teaches: learning_objectives, # TeSS, "learning objectives"
        # timeRequired: 'P1D', # TeSS, "duration", TODO: calculate from start/end date, not implemented in scraper currently.

        availableLanguage: ['en'], # TODO: support other languages
        inLanguage: ['en'], # TODO: support other languages
        # courseCode
        # coursePrerequisites
        # educationalCredentialAwarded
        # financialAidEligible
        # hasCourseInstance
        # numberOfCredits
        # occupationalCredentialAwarded
        # syllabusSections
        # totalHistoricalEnrollment

        # assesses
        # competencyRequired
        # educationalAlignment
        # educationalLevel
        # educationalUse
        # learningResourceType
        # teaches

        funder: funders, # Org or person
        funding: funding, # Grant
        publisher: GTN,
        provider: GTN,
        syllabusSections: syllab,
        # Session materials
        # TODO: not currently parsed by TeSS, google just complains about it, so we're leaving it out.
        # hasPart: parts,
      }

      begin
        data['dateModified'] = Gtn::ModificationTimes.obtain_time(page.path)
        data['datePublished'] = Gtn::PublicationTimes.obtain_time(page.path)
      rescue StandardError
        data['dateModified'] = Gtn::ModificationTimes.obtain_time(page['path'])
        data['datePublished'] = Gtn::PublicationTimes.obtain_time(page['path'])
      end

      if page['cover']
        data['image'] = if page['cover'] =~ /^http/
                          [page['cover']]
                        else
                          ["#{site['url']}#{site['baseurl']}#{page['cover']}"]
                        end
      end

      # We CANNOT guarantee A11Y
      # data.update(A11Y)
      if page['cost'] and page['cost'].downcase == 'free'
        data['isAccessibleForFree'] = true
        offer = {
          '@type': 'Offer',
          price: 0,
          priceCurrency: 'EUR',
          category: 'Free',
          isAccessibleForFree: true,
        }
      elsif page['cost']
        data['isAccessibleForFree'] = false
        offer = {
          '@type': 'Offer',
          price: page['cost'].split[0],
          priceCurrency: page['cost'].split[1],
          isAccessibleForFree: false,
          category: 'Paid',
          # TODO: this can be more advanced but we need to collect start/end times, and timezone.
        }
      end

      # TODO: this is wrong in a whole host of scenarios like incl weekends.
      course_days = (page.fetch('date_end', page['date_start']) - page['date_start']).to_i
      if course_days < 1
        course_days = 1
      end
      data['hasCourseInstance'] = [
        {
          '@type': 'CourseInstance',
          courseMode: page['mode'],
          # courseWorkload: "A daily course running from #{page['date_start']} to #{page['date_end']}",
          offers: offer,
          instructor: instructors,
          isAccessibleForFree: data['isAccessibleForFree'],
          courseSchedule: {
            '@type': 'Schedule',
            startDate: page['date_start'],
            endDate: page.fetch('date_end', page['date_start']),
            repeatCount: course_days,
            repeatFrequency: 'daily', # Contrary to schema.org spec, this is what Google wants.
          },
          courseWorkload: "P#{course_days}D",
        }
      ]

      data['offers'] = [offer]

      if page.key?('location') && page['location'].keys.length > 1
        data['location'] = {
          '@type': 'Place',
          name: page['location']['name'],
          address: {
            '@type': 'PostalAddress',
            streetAddress: page['location'].fetch('address', nil),
            addressLocality: page['location'].fetch('city', nil),
            addressRegion: page['location'].fetch('region', nil),
            postalCode: page['location'].fetch('postcode', nil),
            addressCountry: page['location'].fetch('country', nil)
          }
        }
      end

      JSON.pretty_generate(data)
    end

    ##
    # Generate the JSON-LD metadata for a learning pathway
    # Parameters:
    # +page+:: The page object.
    # +site+:: The +Jekyll::Site+ site object.
    # Returns:
    # +Hash+:: The JSON-LD metadata.
    def generate_learning_pathway_jsonld(page, site)
      materials = []
      page['pathway'].each do |section|
        if !section.key? 'tutorials'
          next
        end

        section['tutorials'].each do |tutorial|
          if tutorial.key?('custom')
            next
          end

          material = TopicFilter.fetch_tutorial_material(site, tutorial['topic'], tutorial['name'])
          materials.push(material)
        end
      end
      materials.compact!

      # Extract EDAM terms from all materials
      edam_terms = materials.map do |material|
        material.fetch('edam_ontology', []).map do |term|
          {
            '@type': 'DefinedTerm',
            '@id': "http://edamontology.org/#{term}",
            inDefinedTermSet: 'http://edamontology.org',
            termCode: term,
          }
        end
      end.flatten.uniq

      learning_objectives = materials.map do |material|
        material.fetch('objectives', [])
      end.flatten.compact

      funders = materials.map do |material|
        Gtn::Contributors.get_funders(site, material).map do |x|
          to_pfo_jsonld(x, site, json: false)
        end
      end.flatten.uniq.compact

      funding = materials.map do |material|
        Gtn::Contributors.get_grants(site, material).map do |x|
          to_pfo_jsonld(x, site, json: false)
        end
      end.flatten.uniq.compact

      # TODO: add topic edam terms too? Not sure.
      parts = []
      materials.each do |material|
        mat = generate_material_jsonld(material, site['data'][material['topic_name']], site)
        if !mat.nil? && !mat.empty?
          parts.push(mat)
        end
      end

      syllab = page['pathway'].reject { |s| s['section'].nil? }.map do |section|
        {
          '@type': 'Syllabus',
          name: section['section'],
          description: section.fetch('description', nil),
        }
      end

      data = {
        '@context': 'https://schema.org',
        '@type': 'Course',
        url: "#{site['url']}#{site['baseurl']}#{page['url']}",
        name: "Learning Pathway #{page['title']}",
        keywords: page['tags'] || [],
        description: page['description'],
        about: edam_terms, # TeSS, "scientific topics".
        audience: page['audience'], # TeSS: target audience
        teaches: learning_objectives, # TeSS, "learning objectives"
        availableLanguage: ['en'], # TODO: support other languages
        inLanguage: ['en'], # TODO: support other languages
        # courseCode
        # coursePrerequisites
        # educationalCredentialAwarded
        # financialAidEligible
        # hasCourseInstance
        # numberOfCredits
        # occupationalCredentialAwarded
        # syllabusSections
        # totalHistoricalEnrollment

        # assesses
        # competencyRequired
        # educationalAlignment
        # educationalLevel
        # educationalUse
        # learningResourceType
        # teaches

        funder: funders, # Org or person
        funding: funding, # Grant
        publisher: GTN,
        provider: GTN,
        syllabusSections: syllab,
        # Session materials
        # TODO: not currently parsed by TeSS, google just complains about it, so we're leaving it out.
        # hasPart: parts,
      }

      begin
        data['dateModified'] = Gtn::ModificationTimes.obtain_time(page.path)
        data['datePublished'] = Gtn::PublicationTimes.obtain_time(page.path)
      rescue StandardError
        data['dateModified'] = Gtn::ModificationTimes.obtain_time(page['path'])
        data['datePublished'] = Gtn::PublicationTimes.obtain_time(page['path'])
      end

      if page['cover']
        data['image'] = if page['cover'] =~ /^http/
                          [page['cover']]
                        else
                          ["#{site['url']}#{site['baseurl']}#{page['cover']}"]
                        end
      end

      # We CANNOT guarantee A11Y
      # data.update(A11Y)
      data['isAccessibleForFree'] = true
      offer = {
        '@type': 'Offer',
        price: 0,
        priceCurrency: 'EUR',
        category: 'Free',
        isAccessibleForFree: true,
      }
      data['offers'] = [offer]

      # TODO: this is basically just wrong.
      data['hasCourseInstance'] = [
        {
          '@type': 'CourseInstance',
          courseMode: 'online',
          offers: offer,
          isAccessibleForFree: data['isAccessibleForFree'],
        }
      ]

      JSON.pretty_generate(data)
    end

    ##
    # Convert a material to JSON-LD, intended to be used in Jekyll Liquid templates.
    # Parameters:
    # +material+:: The material object.
    # +topic+:: The topic object.
    # +site+:: The +Jekyll::Site+ site object.
    #
    # Returns:
    # +String+:: The JSON-LD metadata.
    def to_jsonld(material, topic, site)
      JSON.pretty_generate(generate_material_jsonld(material, topic, site))
    end

    ##
    # Convert a material to JSON-LD.
    # Parameters:
    # +material+:: The material object.
    # +topic+:: The topic object.
    # +site+:: The +Jekyll::Site+ site object.
    #
    # Returns:
    # +Hash+:: The JSON-LD metadata.
    def generate_material_jsonld(material, topic, site)
      langCodeMap = {
        "en" => 'English',
        "es" => 'Español',
        "fr" => 'Français',
      }

      eduLevel = {
        'Introductory' => 'Beginner',
        'Intermediate' => 'Intermediate',
        'Advanced' => 'Advanced'
      }
      return '{}' if !topic

      topic_desc = {
        '@type': 'CreativeWork',
        name: (topic['title']).to_s,
        description: (topic['summary']).to_s,
        url: "#{site['url']}#{site['baseurl']}/topics/#{topic['name']}/"
      }

      # aggregate everything
      data = {
        # Properties from Course
        '@context': 'http://schema.org',
        '@type': 'LearningResource',

        # Required for BioSchemas
        'http://purl.org/dc/terms/conformsTo': {
          '@id': 'https://bioschemas.org/profiles/TrainingMaterial/1.0-RELEASE',
          '@type': 'CreativeWork'
        },

        # Properties from CreativeWork
        # "about" described below
        #
        # "accountablePerson":,
        # "aggregateRating":,
        # "alternativeHeadline":,
        # "associatedMedia":,
        audience: {
          '@type': 'EducationalAudience',
          educationalRole: EDU_ROLES[topic['type']]
        },
        # "audio":,
        # "award":,
        # "author" described below
        # "character":,
        citation: {
          '@type': 'CreativeWork',
          name: 'Community-Driven Data Analysis Training for Biology',
          url: 'https://doi.org/10.1016/j.cels.2018.05.012'
        },
        # "comment":,
        # "commentCount":,
        # "contentLocation":,
        # "contentRating":,
        # "contentReferenceTime":,
        # "contributor" described below
        # copyrightHolder: GTN,
        # copyrightNotice: m
        # "copyrightYear":,
        # "correction":,
        # "creator":,
        # "dateCreated":,
        # "datePublished":,
        discussionUrl: site['gitter_url'],
        # "editor":,
        # "educationalAlignment":,
        # "educationalUse":,
        # "encoding":,
        # "encodingFormat":,
        # "exampleOfWork":,
        # "expires":,
        # "funder": funding,
        # "genre":,
        # "hasPart" described below
        headline: (material['title']).to_s,
        # "interactionStatistic":,
        interactivityType: 'mixed',
        isAccessibleForFree: true,
        # "isBasedOn":,
        isFamilyFriendly: true,
        # "isPartOf" described below
        # "keywords": described below
        # "learningResourceType" described below
        license: 'https://spdx.org/licenses/CC-BY-4.0.html',
        # "locationCreated":,
        # "mainEntity":,
        # "material":,
        # "mentions" described below
        # "offers":,
        # "position":,
        producer: GTN,
        provider: GTN,
        # "publication":,
        # "publisher":,
        # "publisherImprint":,
        # "publishingPrinciples":,
        # "recordedAt":,
        # "releasedEvent":,
        # "review":,
        # "schemaVersion":,
        # "sdDatePublished":,
        # "sdLicense":,
        # "sdPublisher":,
        sourceOrganization: GTN,
        # "spatialCoverage":,
        # "sponsor":,
        # "temporalCoverage":,
        # "text":,
        # "thumbnailUrl":,
        # "timeRequired" described below
        # "translationOfWork":,
        # "translator": Google Translate???,
        # "typicalAgeRange":,
        # "version":,
        # "video":,
        # "workExample":,
        # "workTranslation":,

        # Properties from Thing
        # "additionalType":,
        # "alternateName":,
        # "description" described below
        # "disambiguatingDescription":,
        # "image":,
        # "mainEntityOfPage":,
        # "name" described below
        # "potentialAction":,
        # "sameAs":,
        # "subjectOf":,
        # "url" described below
        workTranslation: [],
        creativeWorkStatus: material['draft'] ? 'Draft' : 'Active',
      }

      if material.key?('pub_date')
        data['dateModified'] = material['mod_date']
        data['datePublished'] = material['pub_date']
      else
        begin
          data['dateModified'] = Gtn::ModificationTimes.obtain_time(material.path)
          data['datePublished'] = Gtn::PublicationTimes.obtain_time(material.path)
        rescue StandardError
          data['dateModified'] = Gtn::ModificationTimes.obtain_time(material['path'])
          data['datePublished'] = Gtn::PublicationTimes.obtain_time(material['path'])
        end
      end

      if material.key?('copyright')
        # copyrightHolder: GTN,
        data['copyrightNotice'] = material['copyright']
      else
        # I'm not sure this is accurate.
        data['copyrightHolder'] = GTN
      end

      funders = Gtn::Contributors.get_funders(site, material).map do |x|
        to_pfo_jsonld(x, site, json: false)
      end
      grants = Gtn::Contributors.get_grants(site, material).map do |x|
        to_pfo_jsonld(x, site, json: false)
      end

      data['funder'] = funders
      data['funding'] = grants

      data['identifier'] = "https://gxy.io/GTN:#{material['short_id']}" if material.key?('short_id')

      data.update(A11Y)

      actual_material = TopicFilter.fetch_tutorial_material(site, material['topic_name'], material['tutorial_name'])

      # info depending if tutorial, hands-on or slide level
      # parts = []
      # data['hasPart'] = parts

      mentions = []
      description = []

      data['isPartOf'] = topic_desc

      data['abstract'] = material
        .fetch('content', '')
        .strip
        .split("\n")
        .first

      if ! data['abstract'].nil?
        data['abstract'] = data['abstract']
          .gsub(/\{\{\s*site.baseurl\s*\}\}/, url_prefix(site))
          .gsub(/\[{{\s*site.url\s*}}/, '[' + url_prefix(site))
          .gsub(/{% link (topics[^%]*).md %}/, url_prefix(site) + '\1.html')
          .gsub(/{% link (topics[^%]*).html %}/, url_prefix(site) + '\1.html')
          .gsub(/\s*\(?{%\s*cite [^}]+\s*%}\)?/, '')
          .gsub('{{ site.github_repository }}', safe_site_config(site, 'github_repository', 'https://example.com'))
          .gsub(/{% snippet ([^%]*) %}/, '')
          .gsub(/{% include ([^%]*) %}/, '')
      end

      description.push("## Abstract\n\n#{data['abstract']}\n\n")

      if (material['name'] == 'tutorial.md') || (material['name'] == 'slides.html')

        if material['name'] == 'tutorial.md'
          data['learningResourceType'] = 'e-learning'
          description.push("## About This Material\n\nThis is a Hands-on Tutorial from the GTN which is usable either for individual self-study, or as a teaching material in a classroom.\n\n")
        else
          data['learningResourceType'] = 'slides'
        end

        data['name'] = material['title']
        data['url'] = "#{site['url']}#{site['baseurl']}#{material['url']}"

        # Requires https://github.com/galaxyproject/training-material/pull/4271
        data['version'] = Gtn::ModificationTimes.obtain_modification_count(material['path'])

        # Time required
        if material.key?('time_estimation') && !material['time_estimation'].nil?
          data['timeRequired'] = "PT#{material['time_estimation'].upcase}"
        end

        # Description with questions, objectives and keypoints
        if material.key?('questions') && !material['questions'].nil? && material['questions'].length.positive?
          questions = material['questions'].join("\n - ")
          description.push("## Questions this #{material['type']} will address\n\n - #{questions}\n\n")
        end
        if material.key?('objectives') && !material['objectives'].nil? && material['objectives'].length.positive?
          objectives = material['objectives'].map{|x| "- #{x}"}.join("\n")
          description.push("## Learning Objectives\n\n#{objectives}\n\n")
          data['teaches'] = objectives
        end
        if material.key?('keypoints') && !material['keypoints'].nil? && material['keypoints'].length.positive?
          keypoints = material['keypoints'].join("\n - ")
          description.push("## Key Points\n\n - #{keypoints}\n\n")
        end

        # Keywords
        data['keywords'] = [topic['title']] + (material['tags'] || [])
        # Zenodo links
      end

      # Mentions are 'external resources' in TeSS.
      # This could be expanded with
      # - supported servers
      # - tools and resources used (e.g. Galaxy) or tools linked to the TS.
      # - slides (if tutorial) and tutorial (if slides)
      # - other materials in the same topic?
      if actual_material.key?('workflows')
        mentions.push({
                        '@type': 'Thing',
                        url: "#{site['url']}#{site['baseurl']}#{material['dir']}workflows/",
                        name: "Associated Workflows"
                      })
      end

      # Notebooks
      if actual_material.key?('notebook')
        if actual_material['notebook']['language'] != 'r'
          # Python, Bash, SQL (all via jupyter)
          url = "#{site['url']}#{site['baseurl']}#{material['dir']}#{material['topic_name']}-#{material['tutorial_name']}.ipynb"
          mentions.push({
                          '@type': 'Thing',
                          url: url,
                          name: "Jupyter Notebook (with Solutions)"
                        })
          mentions.push({
                          '@type': 'Thing',
                          url: url.gsub(/\.ipynb$/, '-course.ipynb'),
                          name: "Jupyter Notebook (without Solutions)"
                        })
        elsif actual_material['notebook']['language'] == 'r' # Actual R
          url = "#{site['url']}#{site['baseurl']}#{material['dir']}#{material['topic_name']}-#{material['tutorial_name']}.Rmd"
          mentions.push({
                          '@type': 'Thing',
                          url: url,
                          name: "Quarto/RMarkdown Notebook"
                        })
          mentions.push({
                          '@type': 'Thing',
                          url: "https://bio.tools/tool/rstudio",
                          name: "RStudio"
                        })
        end
      end

      # Zenodo link out
      if actual_material.key?('zenodo_link') && ! actual_material['zenodo_link'].nil?
        if actual_material['zenodo_link'].length.positive?
          mentions.push({
                          '@type': 'Thing',
                          url: (actual_material['zenodo_link']).to_s,
                          name: "Associated Training Datasets"
                        })
        end
      end

      if description.empty?
        description.push(material.fetch('content', '').strip.split("\n").first)
      end
      data['description'] = description.join("\n")

      data['inLanguage'] = if material.key?('lang')
                             {
                               '@type': 'Language',
                               name: langCodeMap[material['lang']],
                               alternateName: material['lang']
                             }
                           else
                             {
                               '@type': 'Language',
                               name: 'English',
                               alternateName: 'en'
                             }
                           end

      # Course requirements (material + topic)
      reqs = []
      reqs.push(*topic['requirements']) if topic.key?('requirements')
      reqs.push(*material['requirements']) if material.key?('requirements')
      if !reqs.empty?
        coursePrerequisites = []
        reqs.each do |req|
          if req['type'] == 'internal'
            if req.key?('tutorials')
              (req['tutorials']).each do |tuto|
                (site['pages']).each do |page|
                  if ((page['name'] == 'tutorial.md') || (page['name'] == 'slides.html')) &&
                     ((page['topic_name'] == req['topic_name']) && (page['tutorial_name'] == tuto))
                    # slides
                    if page['name'] == 'slides.html'
                      coursePrerequisites.push(
                        {
                          '@context': 'http://schema.org',
                          '@type': 'LearningResource',
                          url: "#{site['url']}#{site['baseurl']}/topics/#{req['topic_name']}/" \
                               "tutorials/#{tuto}/slides.html",
                          name: (page['title']).to_s,
                          description: "Slides for '#{page['title']}' tutorial",
                          learningResourceType: 'slides',
                          interactivityType: 'expositive',
                          provider: GTN
                        }
                      )
                      if page['hands_on_url']
                        coursePrerequisites.push(
                          {
                            '@context': 'http://schema.org',
                            '@type': 'LearningResource',
                            url: (page['hands_on_url']).to_s,
                            learningResourceType: 'e-learning',
                            interactivityType: 'expositive',
                          }
                        )
                      end
                    end
                    # hands-on
                    if page['name'] == 'tutorial.md'
                      coursePrerequisites.push(
                        {
                          '@context': 'http://schema.org',
                          '@type': 'LearningResource',
                          url: "#{site['url']}#{site['baseurl']}/topics/#{req['topic_name']}/tutorials" \
                               "/#{tuto}/tutorial.html",
                          name: (page['title']).to_s,
                          description: "Hands-on for '#{page['title']}' tutorial",
                          learningResourceType: 'e-learning',
                          interactivityType: 'expositive',
                          provider: GTN
                        }
                      )
                    end
                  end
                end
              end
            else
              coursePrerequisites.push(
                {
                  '@context': 'http://schema.org',
                  '@type': 'LearningResource',
                  url: "#{site['url']}#{site['baseurl']}/topics/#{req['topic_name']}/",
                  name: (site['data'][req['topic_name']]['title']).to_s,
                  description: (site['data'][req['topic_name']]['title']).to_s,
                  provider: GTN
                }
              )
            end
          elsif req['type'] == 'external'
            coursePrerequisites.push({
                                       '@type': 'CreativeWork',
                                       url: (req['link']).to_s,
                                       name: (req['title']).to_s
                                     })
          else
            coursePrerequisites.push((req['title']).to_s)
          end
        end
        data['competencyRequired'] = coursePrerequisites.uniq
      end

      # Add contributors/authors
      if material.key?('contributors') || material.key?('contributions')
        authors = Gtn::Contributors.get_authors(material).map do |x|
          generate_person_jsonld(x, Gtn::Contributors.fetch_contributor(site, x), site)
        end

        data['author'] = authors
      end

      # Add non-author contributors
      if material.key?('contributions')
        data['contributor'] = Gtn::Contributors.get_non_authors(material).map do |x|
          generate_person_jsonld(x, site['data']['contributors'][x], site)
        end
      end

      about = []
      about.push(topic_desc)
      edam_terms = topic.fetch('edam_ontology', []) | material.fetch('edam_ontology', [])

      about += edam_terms.map do |term|
        {
          '@type': 'DefinedTerm',
          '@id': "http://edamontology.org/#{term}",
          inDefinedTermSet: 'http://edamontology.org',
          termCode: term,
          # "name": ,
          url: 'https://bioportal.bioontology.org/ontologies/EDAM/?p=classes&conceptid=' \
               "http%3A%2F%2Fedamontology.org%2F#{term}"
        }
      end

      data['about'] = about

      data['educationalLevel'] = material.key?('level') ? eduLevel[material['level']] : 'Beginner'
      data['mentions'] = mentions

      data
    end
  end
end

Liquid::Template.register_filter(Jekyll::JsonldFilter)
