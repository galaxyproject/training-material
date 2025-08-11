#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'pathname'
require 'find'
require 'bibtex'
require 'json'
require 'kramdown'
require 'kramdown-parser-gfm'
require 'citeproc/ruby'
require 'csl/styles'
require './_plugins/util'

GTN_HOME = Pathname.new(__dir__).parent.to_s


module Gtn

  # A custom module to properly format reviewdog json output
  module ReviewDogEmitter
    @CODE_URL = 'https://training.galaxyproject.org/training-material/gtn_rdoc/Gtn/Linter.html'

    def self.delete_text(path: '', idx: 0, text: '', message: 'No message', code: 'GTN000', full_line: '', fn: '')
      error(
        path: path,
        idx: idx,
        match_start: 0,
        match_end: text.length,
        replacement: '',
        message: message,
        code: code,
        full_line: full_line,
        fn: fn,
      )
    end

    def self.file_error(path: '', message: 'None', code: 'GTN:000', fn: '')
      error(
        path: path,
        idx: 0,
        match_start: 0,
        match_end: 1,
        replacement: nil,
        message: message,
        code: code,
        full_line: '',
        fn: fn
      )
    end

    def self.warning(path: '', idx: 0, match_start: 0, match_end: 1,
                     replacement: nil, message: 'No message', code: 'GTN000', full_line: '', fn: '')
      self.message(
        path: path,
        idx: idx,
        match_start: match_start,
        match_end: match_end,
        replacement: replacement,
        message: message,
        level: 'WARNING',
        code: code,
        full_line: full_line,
        fn: fn,
      )
    end

    def self.error(path: '', idx: 0, match_start: 0, match_end: 1, replacement: nil, message: 'No message',
                   code: 'GTN000', full_line: '', fn: '')
      self.message(
        path: path,
        idx: idx,
        match_start: match_start,
        match_end: match_end,
        replacement: replacement,
        message: message,
        level: 'ERROR',
        code: code,
        full_line: full_line,
        fn: fn,
      )
    end

    def self.message(path: '', idx: 0, match_start: 0, match_end: 1, replacement: nil, message: 'No message',level: 'WARNING', code: 'GTN000', full_line: '', fn: '')
      end_area = { 'line' => idx + 1, 'column' => match_end }
      end_area = { 'line' => idx + 2, 'column' => 1 } if match_end == full_line.length

      res = {
        'message' => message,
        'location' => {
          'path' => path,
          'range' => {
            'start' => { 'line' => idx + 1, 'column' => match_start + 1 },
            'end' => end_area
          }
        },
        'severity' => level
      }
      if !code.nil?
        res['code'] = {
          'value' => code
        }
        if !fn.nil?
          res['code']['url'] = "#{@CODE_URL}#method-c-#{fn}"
        end
      end
      if !replacement.nil?
        res['suggestions'] = [{
          'text' => replacement,
          'range' => {
            'start' => { 'line' => idx + 1, 'column' => match_start + 1 },
            'end' => end_area
          }
        }]
      end
      res
    end
  end

  # This is our ONE central linting script that handles EVERYTHING.
  module Linter
    @BAD_TOOL_LINK = /{% tool (\[[^\]]*\])\(\s*https?.*tool_id=([^)]*)\)\s*%}/i
    @BAD_TOOL_LINK2 = %r{{% tool (\[[^\]]*\])\(\s*https://toolshed.g2([^)]*)\)\s*%}}i
    @MAYBE_OK_TOOL_LINK = /{% tool (\[[^\]]*\])\(([^)]*)\)\s*%}/i

    def self.find_matching_texts(contents, query)
      contents.map.with_index do |text, idx|
        [idx, text, text.match(query)]
      end.select { |_idx, _text, selected| selected }
    end

    ##
    # GTN:001 - Setting no_toc is discouraged as headers are useful for learners to link to and to jump to. Setting no_toc removes it from the table of contents which is generally inadvisable.
    #
    # Remediation: remove {: .no_toc}
    def self.fix_notoc(contents)
      find_matching_texts(contents, /{:\s*.no_toc\s*}/)
        .map do |idx, text, _selected|
        ReviewDogEmitter.delete_text(
          path: @path,
          idx: idx,
          text: text,
          message: 'Setting no_toc is discouraged, these headings provide useful places for readers to jump to.',
          code: 'GTN:001',
          full_line: text,
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:002 - YouTube links are discouraged. Please consider using our include for it:
    #
    # E.g, instead of
    #
    #   <iframe ... youtube.../>
    #
    # Consider:
    #
    #   {% include _includes/youtube.html id="e0vj-0imOLw" title="Difference between climate and weather" %}
    def self.youtube_bad(contents)
      find_matching_texts(contents, %r{<iframe.*youtu.?be.*</iframe>})
        .map do |idx, _text, selected|
        ReviewDogEmitter.warning(
          path: @path,
          idx: idx,
          match_start: selected.begin(0),
          match_end: selected.end(0) + 1,
          replacement: '',
          message: 'Instead of embedding IFrames to YouTube contents, consider adding this video to the ' \
                   'GTN tutorial "recordings" metadata where it will ' \
                   'be more visible for others.',
          code: 'GTN:002',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:003 - We discourage linking to training.galaxyproject.org or
    # galaxyproject.github.io/training-material as those are "external" links,
    # which are slower for us to validate. Every build we run tests to be sure
    # that every link is valid, but we cannot do that for every external site to
    # avoid putting unnecessary pressure on them.
    #
    # Instead of
    #
    #   [see this other tutorial(https://training.galaxyproject.org/training-material/topics/admin/tutorials/ansible/tutorial.html)
    #
    # Consider:
    #
    #   [see this other tutorial({% link topics/admin/tutorials/ansible/tutorial.md %})
    def self.link_gtn_tutorial_external(contents)
      find_matching_texts(
        contents,
        %r{\(https?://(training.galaxyproject.org|galaxyproject.github.io)/training-material/([^)]*)\)}
      )
        .map do |idx, _text, selected|
        # puts "#{idx} 0 #{selected[0]} 1 #{selected[1]} 2 #{selected[2]} 3 #{selected[3]}"
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          # We wrap the entire URL (inside the explicit () in a matching group to make it easy to select/replace)
          match_start: selected.begin(0) + 1,
          match_end: selected.end(0),
          replacement: "{% link #{selected[2].gsub('.html', '.md')} %}",
          message: 'Please use the link function to link to other pages within the GTN. ' \
                   'It helps us ensure that all links are correct',
          code: 'GTN:003',
          fn: __method__.to_s,
        )
      end
    end


    ##
    # GTN:003 - We discourage linking to training.galaxyproject.org or
    # galaxyproject.github.io/training-material as those are "external" links,
    # which are slower for us to validate. Every build we run tests to be sure
    # that every link is valid, but we cannot do that for every external site to
    # avoid putting unnecessary pressure on them.
    #
    # Instead of
    #
    #   [see this other tutorial(https://training.galaxyproject.org/training-material/topics/admin/tutorials/ansible/slides.html)
    #
    # Consider:
    #
    #   [see this other tutorial({% link topics/admin/tutorials/ansible/slides.html %})
    def self.link_gtn_slides_external(contents)
      find_matching_texts(
        contents,
        %r{\((https?://(training.galaxyproject.org|galaxyproject.github.io)/training-material/(.*slides.html))\)}
      )
        .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: "{% link #{selected[3]} %}",
          message: 'Please use the link function to link to other pages within the GTN. ' \
                   'It helps us ensure that all links are correct',
          code: 'GTN:003',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:004 - Instead of linking directly to a DOI and citing it yourself, consider obtaining the BibTeX formatted citation and adding it to a tutorial.bib (or slides.bib) file. Then we can generate a full set of references for the citations and give proper credit.
    #
    # Companion function to Gtn::Linter.check_pmids
    def self.check_dois(contents)
      find_matching_texts(contents, %r{(\[[^\]]*\]\(https?://doi.org/[^)]*\))})
        .reject { |_idx, _text, selected| selected[0].match(%r{10.5281/zenodo}) } # Ignoring zenodo
        .map do |idx, _text, selected|
        ReviewDogEmitter.warning(
          path: @path,
          idx: idx,
          match_start: selected.begin(0),
          match_end: selected.end(0) + 2,
          replacement: '{% cite ... %}',
          message: 'This looks like a DOI which could be better served by using the built-in Citations mechanism. ' \
                   'You can use https://doi2bib.org to convert your DOI into a .bib formatted entry, ' \
                   'and add to your tutorial.md',
          code: 'GTN:004',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:004 - Instead of linking directly to a PMID URL and citing it yourself, consider obtaining the BibTeX formatted citation and adding it to a tutorial.bib (or slides.bib) file. Then we can generate a full set of references for the citations and give proper credit.
    #
    # Companion function to Gtn::Linter.check_dois
    def self.check_pmids(contents)
      # https://www.ncbi.nlm.nih.gov/pubmed/24678044
      find_matching_texts(contents,
                          %r{(\[[^\]]*\]\(https?://www.ncbi.nlm.nih.gov/pubmed//[0-9]*\))}).map do |idx, _text, selected|
        ReviewDogEmitter.warning(
          path: @path,
          idx: idx,
          match_start: selected.begin(0),
          match_end: selected.end(0) + 2,
          replacement: '{% cite ... %}',
          message: 'This looks like a PMID which could be better served by using the built-in Citations mechanism. ' \
                   'You can use https://doi2bib.org to convert your PMID/PMCID into a .bib formatted entry, ' \
                   'and add to your tutorial.md',
          code: 'GTN:004',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:005 - Using link names like 'here' are unhelpful for learners who are progressing through the material with a screenreader. Please use a more descriptive text for your linke
    #
    # Instead of
    #
    #   see the documentation [here](https://example.com)
    #
    # Consider
    #
    #   see [edgeR's documentation](https://example.com)
    def self.check_bad_link_text(contents)
      find_matching_texts(contents, /\[\s*(here|link)\s*\]/i)
        .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(0),
          match_end: selected.end(0) + 1,
          replacement: '[Something better here]',
          message: "Please do not use 'here' as your link title, it is " \
                   '[bad for accessibility](https://usability.yale.edu/web-accessibility/articles/links#link-text). ' \
                   'Instead try restructuring your sentence to have useful descriptive text in the link.',
          code: 'GTN:005',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:006 - This is a potentially incorrect Jekyll/Liquid template function/variable access.
    #
    # Variables can be placed into your template like so:
    #
    #   {{ page.title }}
    #
    # And functions can be called like so:
    #
    #   {% if page.title %}
    #
    # So please be sure {{ }} and {% %} are matching.
    def self.incorrect_calls(contents)
      a = find_matching_texts(contents, /([^{]|^)(%\s*[^%]*%})/i)
          .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(2),
          match_end: selected.end(2) + 1,
          replacement: "{#{selected[2]}",
          message: 'It looks like you might be missing the opening { of a jekyll function',
          code: 'GTN:006',
          fn: __method__.to_s,
        )
      end
      b = find_matching_texts(contents, /{([^%]\s*[^%]* %})/i)
          .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: "%#{selected[1]}",
          message: 'It looks like you might be missing the opening % of a jekyll function',
          code: 'GTN:006',
          fn: __method__.to_s,
        )
      end

      c = find_matching_texts(contents, /({%\s*[^%]*%)([^}]|$)/i)
          .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 2,
          replacement: "#{selected[1]}}#{selected[2]}",
          message: 'It looks like you might be missing the closing } of a jekyll function',
          code: 'GTN:006',
          fn: __method__.to_s,
        )
      end

      d = find_matching_texts(contents, /({%\s*[^}]*[^%])}/i)
          .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: "#{selected[1]}%",
          message: 'It looks like you might be missing the closing % of a jekyll function',
          code: 'GTN:006',
          fn: __method__.to_s,
        )
      end
      a + b + c + d
    end

    @CITATION_LIBRARY = nil

    def self.citation_library
      if @CITATION_LIBRARY.nil?
        lib = BibTeX::Bibliography.new
        (enumerate_type(/bib$/) + enumerate_type(/bib$/, root_dir: 'faqs')).each do |path|
          b = BibTeX.open(path)
          b.each do |x|
            # Record the bib path.
            x._path = path
            lib << x
          end
        end
        @CITATION_LIBRARY = lib
      end

      @CITATION_LIBRARY
    end

    @JEKYLL_CONFIG = nil

    def self.jekyll_config
      if @JEKYLL_CONFIG.nil?
        # Load
        @JEKYLL_CONFIG = YAML.load_file('_config.yml')
      end
      @JEKYLL_CONFIG
    end

    ##
    # GTN:007 - We could not find a citation key, please be sure it is used in a bibliography somewhere.
    def self.check_bad_cite(contents)
      find_matching_texts(contents, /{%\s*cite\s+([^%]*)\s*%}/i)
        .map do |idx, _text, selected|
        citation_key = selected[1].strip
        if citation_library[citation_key].nil?
          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: selected.begin(0),
            match_end: selected.end(0),
            replacement: nil,
            message: "The citation (#{citation_key}) could not be found.",
            code: 'GTN:007',
            fn: __method__.to_s,
          )
        end
      end
    end

    ##
    # GTN:033 - This icon is not known to use. If it is new, please add it to {our configuration.}[https://training.galaxyproject.org/training-material/topics/contributing/tutorials/create-new-tutorial-content/faqs/icons_list.html]
    def self.check_bad_icon(contents)
      find_matching_texts(contents, /{%\s*icon\s+([^%]*)\s*%}/i)
        .map do |idx, _text, selected|
        icon_key = selected[1].strip.split[0]
        if jekyll_config['icon-tag'][icon_key].nil?
          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: selected.begin(0),
            match_end: selected.end(0),
            replacement: nil,
            message: "The icon (#{icon_key}) could not be found, please add it to _config.yml.",
            code: 'GTN:033',
            fn: __method__.to_s,
          )
        end
      end
    end

    ##
    # GTN:008 - This snippet is not known to us, please check that it exists somewhere in the snippets/ folder.
    def self.non_existent_snippet(contents)
      find_matching_texts(contents, /{%\s*snippet\s+([^ ]*)/i)
        .reject do |_idx, _text, selected|
        File.exist?(selected[1])
      end
        .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(0),
          match_end: selected.end(0),
          replacement: nil,
          message: "This snippet (`#{selected[1]}`) does not seem to exist",
          code: 'GTN:008',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:009 - This looks like an invalid tool link. There are several ways that tool links can be invalid, and only one correct way to reference a tool
    #
    # Correct
    #
    #   {% tool [JBrowse genome browser](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.4+galaxy3) %}
    #
    # Incorrect
    #
    #   {% tool [JBrowse genome browser](https://toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.4+galaxy3) %}
    #   {% tool [JBrowse genome browser](https://toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse) %}
    #   {% tool [JBrowse genome browser](jbrowse/1.16.4+galaxy3) %}
    #   {% tool [JBrowse genome browser](https://toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/deadbeefcafe) %}
    def self.bad_tool_links(contents)
      find_matching_texts(contents, @BAD_TOOL_LINK) + \
        find_matching_texts(contents, @BAD_TOOL_LINK2)
        .map do |idx, _text, selected|
          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: selected.begin(0),
            match_end: selected.end(0) + 1,
            replacement: "{% tool #{selected[1]}(#{selected[2]}) %}",
            message: 'You have used the full tool URL to a specific server, here we only need the tool ID portion.',
            code: 'GTN:009',
            fn: __method__.to_s,
          )
        end

      find_matching_texts(contents, @MAYBE_OK_TOOL_LINK)
        .map do |idx, _text, selected|

          if acceptable_tool?(selected[2])
            next
          end

          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: selected.begin(0),
            match_end: selected.end(0) + 1,
            replacement: "{% tool #{selected[1]}(#{selected[2]}) %}",
            message: 'You have used an invalid tool URL, it should be of the form "toolshed.g2.bx.psu.edu/repos/{owner}/{repo}/{tool}/{version}" (or an internal tool ID) so, please double check.',
            code: 'GTN:009'
          )
        end
    end

    ##
    # GTN:040 - zenodo.org/api links are invalid in the GTN, please use the zenodo.org/records/id/files/<filename> format instead. This ensures that when users download files from zenodo into Galaxy, they appear correctly, with a useful filename.
    def self.bad_zenodo_links(contents)
      find_matching_texts(contents, /https:\/\/zenodo.org\/api\//)
        .reject { |_idx, _text, selected| _text =~ /files-archive/ }
        .map do |idx, _text, selected|
          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: selected.begin(0),
            match_end: selected.end(0) + 1,
            replacement: nil,
            message: 'Please do not use zenodo.org/api/ links, instead it should look like zenodo.org/records/id/files/<filename>',
            code: 'GTN:040',
            fn: __method__.to_s,
          )
        end
    end

    ##
    # GTN:032 - Snippets are too close together which sometimes breaks snippet rendering. Please ensure snippets are separated by one line.
    def self.snippets_too_close_together(contents)
      prev_line = -2
      res = []
      find_matching_texts(contents, /^[> ]*{% snippet/)
        .each do |idx, _text, selected|
        if idx == prev_line + 1
          res.push(ReviewDogEmitter.error(
                     path: @path,
                     idx: idx,
                     match_start: selected.begin(0),
                     match_end: selected.end(0) + 1,
                     replacement: nil,
                     message: 'Snippets too close together',
                     code: 'GTN:032',
                     fn: __method__.to_s,
                   ))
        end
        prev_line = idx
      end
      res
    end

    ##
    # GTN:009 - See Gtn::Linter.bad_tool_links
    def self.check_tool_link(contents)
      find_matching_texts(contents, /{%\s*tool \[([^\]]*)\]\(([^)]*)\)\s*%}/)
        .map do |idx, _text, selected|
        # text = selected[1]
        link = selected[2]

        errs = []
        if link.match(%r{/})
          if link.count('/') < 5
            errs.push(ReviewDogEmitter.error(
                        path: @path,
                        idx: idx,
                        match_start: selected.begin(2),
                        match_end: selected.end(2) + 1,
                        replacement: nil,
                        message: "This tool identifier looks incorrect, it doesn't have the right number of segments.",
                        code: 'GTN:009'
                      ))
          end

          if link.match(/testtoolshed/)
            errs.push(ReviewDogEmitter.warning(
                        path: @path,
                        idx: idx,
                        match_start: selected.begin(2),
                        match_end: selected.end(2) + 1,
                        replacement: nil,
                        message: 'The GTN strongly avoids using testtoolshed tools in your tutorials or workflows',
                        code: 'GTN:009'
                      ))
          end
        else
          if link.match(/\+/)
            errs.push(ReviewDogEmitter.error(
                        path: @path,
                        idx: idx,
                        match_start: selected.begin(2),
                        match_end: selected.end(2) + 1,
                        replacement: nil,
                        message: 'Broken tool link, unnecessary +',
                        code: 'GTN:009'
                      ))
          end


          if !acceptable_tool?(link)
            errs.push(ReviewDogEmitter.error(
                        path: @path,
                        idx: idx,
                        match_start: selected.begin(2),
                        match_end: selected.end(2) + 1,
                        replacement: nil,
                        message: 'Unknown short tool ID. Please use the full tool ID, or check bin/lint.rb ' \
                                 'if you believe this is correct.',
                        code: 'GTN:009'
                      ))
          end
        end

        errs
      end
    end

    ##
    # GTN:010 - We have a new, more accessible syntax for box titles. Please use this instead:
    #
    #   > <box-title>Some Title</box-title>
    #   > ...
    #   {: .box}
    def self.new_more_accessible_boxes(contents)
      #  \#\#\#
      find_matching_texts(contents, /> (### {%\s*icon ([^%]*)\s*%}[^:]*:?(.*))/)
        .map do |idx, _text, selected|
        key = selected[2].strip.gsub(/_/, '-')
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: "<#{key}-title>#{selected[3].strip}</#{key}-title>",
          message: 'We have developed a new syntax for box titles, please consider using this instead.',
          code: 'GTN:010',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:010 - See Gtn::Linter.new_more_accessible_boxes_agenda
    def self.new_more_accessible_boxes_agenda(contents)
      #  \#\#\#
      find_matching_texts(contents, /> (###\s+Agenda\s*)/)
        .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: '<agenda-title></agenda-title>',
          message: 'We have developed a new syntax for box titles, please consider using this instead.',
          code: 'GTN:010',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:011 - Do not use target="_blank" it is bad for accessibility.
    def self.no_target_blank(contents)
      find_matching_texts(contents, /target=("_blank"|'_blank')/)
        .map do |idx, _text, selected|
        ReviewDogEmitter.warning(
          path: @path,
          idx: idx,
          match_start: selected.begin(0),
          match_end: selected.end(0),
          replacement: nil,
          message: 'Please do not use `target="_blank"`, [it is bad for accessibility.]' \
                   '(https://www.a11yproject.com/checklist/#identify-links-that-open-in-a-new-tab-or-window)',
          code: 'GTN:011',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:034 - Alternative text or alt-text is mandatory for every image in the GTN.
    def self.empty_alt_text(contents)
      find_matching_texts(contents, /!\[\]\(/i)
        .map do |idx, _text, selected|
        path = selected[1].to_s.strip
        if !File.exist?(path.gsub(%r{^/}, ''))
          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: selected.begin(0),
            match_end: selected.end(0),
            replacement: nil,
            message: 'The alt text for this image seems to be empty',
            code: 'GTN:034',
            fn: __method__.to_s,
          )
        end
      end
    end

    ##
    # GTN:018 - You have linked to a file but this file could not be found. Check your link to make sure the path exists.
    #
    # Note that we use a customised version of Jekyll's link function which will not work correctly for _posts/news items, which should be corrected at some point. We should remove our custom link function and go back to the official one.
    def self.check_bad_link(contents)
      find_matching_texts(contents, /{%\s*link\s+([^%]*)\s*%}/i)
        .map do |idx, _text, selected|
        path = selected[1].to_s.strip
        if !File.exist?(path.gsub(%r{^/}, ''))
          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: selected.begin(0),
            match_end: selected.end(0),
            replacement: nil,
            message: "The linked file (`#{selected[1].strip}`) could not be found.",
            code: 'GTN:018',
            fn: __method__.to_s,
          )
        end
      end

      find_matching_texts(contents, /\]\(\)/i)
        .map do |idx, _text, selected|
        path = selected[1].to_s.strip
        if !File.exist?(path.gsub(%r{^/}, ''))
          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: selected.begin(0),
            match_end: selected.end(0),
            replacement: nil,
            message: 'The link does not seem to have a target.',
            code: 'GTN:018',
            fn: __method__.to_s,
          )
        end
      end
    end

    ##
    # GTN:036 - You have used the TRS snippet to link to a TRS ID but the link does not seem to be correct.
    def self.check_bad_trs_link(contents)
      find_matching_texts(contents, %r{snippet faqs/galaxy/workflows_run_trs.md path="([^"]*)"}i)
        .map do |idx, _text, selected|
        path = selected[1].to_s.strip
        if !File.exist?(path)
          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: selected.begin(0),
            match_end: selected.end(0),
            replacement: nil,
            message: "The linked file (`#{path}`) could not be found.",
            code: 'GTN:036',
            fn: __method__.to_s,
          )
        end
      end
    end

    ##
    # GTN:020 - Please do not bold random lines, use a heading properly.
    def self.check_looks_like_heading(contents)
      # TODO: we should remove this someday, but, we need to have a good solution
      # and we're still a ways from that.
      #
      # There's no clear way to say "this subsection of the content has its own hierarchy"
      return if @path.match(/faq/)

      find_matching_texts(contents, /^\*\*(.*)\*\*$/)
        .map do |idx, _text, selected|
        ReviewDogEmitter.warning(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: "### #{selected[1]}",
          message: "This looks like a heading, but isn't. Please use proper semantic headings where possible. " \
                   'You should check the heading level of this suggestion, rather than accepting the change as-is.',
          code: 'GTN:020',
          fn: __method__.to_s,
        )
      end
    end

    @KNOWN_TAGS = [
      # GTN
      'cite',
      'snippet',
      'link',
      'icon',
      'tool',
      'color',

      'set', # This isn't strictly GTN, it's seen inside a raw in a tool tutorial.
      # Jekyll
      'if', 'else', 'elsif', 'endif',
      'capture', 'assign', 'include',
      'comment', 'endcomment',
      'for', 'endfor',
      'unless', 'endunless',
      'raw', 'endraw'
    ].freeze

    ##
    # GTN:021 - We are not sure this tag is correct, there is a very limited set of Jekyll/liquid tags that are used in GTN tutorials, and this checks for surprises.
    def self.check_bad_tag(contents)
      find_matching_texts(contents, /{%\s*(?<tag>[a-z]+)/)
        .reject { |_idx, _text, selected| @KNOWN_TAGS.include? selected[:tag] }
        .map do |idx, _text, selected|
        ReviewDogEmitter.warning(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: nil,
          message: "We're not sure this tag is correct (#{selected[:tag]}), it isn't one of the known tags.",
          code: 'GTN:021',
          fn: __method__.to_s,
        )
      end
    end

    @BOX_CLASSES = %w[
      agenda
      code-in
      code-out
      comment
      details
      feedback
      hands-on
      hands_on
      question
      solution
      tip
      warning
    ].freeze

    ##
    # GTN:022 - Please do not prefix your boxes of a type with the box name.
    #
    # Do not do:
    #
    #   > <question-title>Question: Some question!</question-title>
    #
    # Instead:
    #
    #   > <question-title>Some question!</question-title>
    #
    # As the Question: prefix will be added automatically when necessary. This goes also for tip/comment/etc.
    def self.check_useless_box_prefix(contents)
      find_matching_texts(contents, /<(?<tag>[a-z_-]+)-title>(?<fw>[a-zA-Z_-]+:?\s*)/)
        .select do |_idx, _text, selected|
        @BOX_CLASSES.include?(selected[:tag]) and selected[:tag] == selected[:fw].gsub(/:\s*$/, '').downcase
      end
        .map do |idx, _text, selected|
        ReviewDogEmitter.warning(
          path: @path,
          idx: idx,
          match_start: selected.begin(2),
          match_end: selected.end(2) + 1,
          replacement: '',
          message: "It is no longer necessary to prefix your #{selected[:tag]} box titles with " \
                   "#{selected[:tag].capitalize}, this is done automatically.",
          code: 'GTN:022',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:028 - Your headings are out of order. Please check it properly, that you do not skip levels.
    def self.check_bad_heading_order(contents)
      doc = Kramdown::Document.new(contents.join("\n"), input: 'GFM')
      headers = doc.root.children.select{|k| k.type == :header}

      bad_depth = headers
        .each_cons(2) # Two items at a time
        .select{|k1, k2| k2.options[:level] - k1.options[:level] > 1} # All that have a >1 shift in heading depth
        .map{|_,b | b} # Only the second, failing one.

      all_headings = headers
        .map{|k| "#" * k.options[:level] + " "+ k.options[:raw_text] }

      bad_depth.map{|k|
        ReviewDogEmitter.error(
          path: @path,
          idx: k.options[:location] - 1,
          match_start: 0,
          match_end: k.options[:raw_text].length + k.options[:level] + 1,
          replacement: '#' * (k.options[:level] - 1),
          message: "You have skipped a heading level, please correct this.\n<details>" \
                   "<summary>Listing of Heading Levels</summary>\n\n```\n#{all_headings.join("\n")}\n```\n</details>",
          code: 'GTN:028',
          fn: __method__.to_s,
        )
      }
    end

    ##
    # GTN:029 - Please do not bold headings
    def self.check_bolded_heading(contents)
      find_matching_texts(contents, /^#+ (?<title>\*\*.*\*\*)$/)
        .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: selected[:title][2..-3],
          message: 'Please do not bold headings, it is unncessary ' \
                   'and will potentially cause screen readers to shout them.',
          code: 'GTN:029',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:032 - zenodo.org/api links are invalid in the GTN, please use the zenodo.org/records/id/files/<filename> format instead. This ensures that when users download files from zenodo into Galaxy, they appear correctly, with a useful filename.
    #
    # Seems to be a duplicate of Gtn::Linter.bad_zenodo_links
    def self.zenodo_api(contents)
      find_matching_texts(contents, %r{(zenodo\.org/api/files/)})
        .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: nil,
          message: 'The Zenodo.org/api URLs are not stable, you must use a URL of the format zenodo.org/record/...',
          code: 'GTN:032'
        )
      end
    end

    ##
    # GTN:035 - This is a non-semantic list which is bad for accessibility and screenreaders.
    #
    # Do not do:
    #
    #   * Step 1. Some text
    #   * Step 2. some other thing
    #
    # Do not do:
    #
    #   Step 1. Some text
    #   Step 2. some other thing
    #
    # Instead:
    #
    #   1. some text
    #   2. some other
    #
    # That is a proper semantic list.
    def self.nonsemantic_list(contents)
      find_matching_texts(contents, />\s*(\*\*\s*[Ss]tep)/)
        .map do |idx, _text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: nil,
          message: 'This is a non-semantic list which is bad for accessibility and bad for screenreaders. ' \
                   'It results in poorly structured HTML and as a result is not allowed.',
          code: 'GTN:035',
          fn: __method__.to_s,
        )
      end
    end

    ##
    # GTN:041, GTN:042, GTN:043, GTN:044, GTN:045 - This checks for a myriad variety of CYOA issues. Please see the error message for help resolving them.
    def self.cyoa_branches(contents)
      joined_contents = contents.join("\n")
      cyoa_branches = joined_contents.scan(/_includes\/cyoa-choices[^%]*%}/m)
        .map{|cyoa_line|
          cyoa_line.gsub(/\n/, ' ') # Remove newlines, want it all one one line.
            .gsub(/\s+/, ' ') # Collapse multiple whitespace for simplicity
            .gsub(/_includes\/cyoa-choices.html/, '').gsub(/%}$/, '') # Strip start/end
            .strip
            .split('" ') # Split on the end of an option to get the individual option groups
            .map{|p| p.gsub(/="/, '=').split('=')}.to_h} # convert it into a convenient hash
      # NOTE: Errors on this line usually mean that folks have used ' instead of " in their CYOA.


      # cyoa_branches =
      # [{"option1"=>"Quick one tool method",
      #   "option2"=>"Convert to AnnData object compatible with Filter, Plot, Explore workflow",
      #   "default"=>"Quick one tool method",
      #   "text"=>"Choose below if you just want to convert your object quickly or see how it all happens behind the scenes!",
      #   "disambiguation"=>"seurat2anndata\""},

      # We use slugify_unsafe to convert it to a slug, now we should check:
      # 1. Is it unique in the file? No duplicate options?
      # 2. Is every branch used?

      # Uniqueness:
      options = cyoa_branches.map{|o| o.select{|k, v| k =~ /option/}.values}.flatten
      slugified = options.map{|o| [o, unsafe_slugify(o)]}
      slugified_grouped = slugified.group_by{|before, after| after}
        .map{|k, pairs| [k, pairs.map{|p| p[0]}]}.to_h

      errors = []
      if slugified_grouped.values.any?{|v| v.length > 1}
        dupes = slugified_grouped.select{|k, v| v.length > 1}
        msg = "We identified the following duplicate options in your CYOA: "
        msg += dupes.map do |slug, options|
          "Options #{options.join(', ')} became the key: #{slug}"
        end.join("; ")

        errors << ReviewDogEmitter.error(
          path: @path,
          idx: 0,
          match_start: 0,
          match_end: 1,
          replacement: nil,
          message: 'You have non-unique options in your Choose Your Own Adventure. Please ensure that each option is unique in its text. Unfortunately we do not currently support re-using the same option text across differently disambiguated CYOA branches, so, please inform us if this is a requirement for you.' + msg,
          code: 'GTN:041',
          fn: __method__.to_s,
        )
      end

      # Missing default
      cyoa_branches.each do |branch|
        if branch['default'].nil?
          errors << ReviewDogEmitter.error(
            path: @path,
            idx: 0,
            match_start: 0,
            match_end: 1,
            replacement: nil,
            message: 'We recommend specifying a default for every branch',
            code: 'GTN:042',
            fn: __method__.to_s,
          )
        end

        # Checking default/options correspondence.
        options = branch.select{|k, v| k =~ /option/}.values
        if branch.key?("default") && ! options.include?(branch['default'])
          if options.any?{|o| unsafe_slugify(o) == unsafe_slugify(branch['default'])}
            errors << ReviewDogEmitter.warning(
              path: @path,
              idx: 0,
              match_start: 0,
              match_end: 1,
              replacement: nil,
              message: "We did not see a corresponding option# for the default: «#{branch['default']}», but this could have been written before we automatically slugified the options. If you like, please consider making your default option match the option text exactly.",
              code: 'GTN:043',
              fn: __method__.to_s,
            )
          else
            errors << ReviewDogEmitter.warning(
              path: @path,
              idx: 0,
              match_start: 0,
              match_end: 1,
              replacement: nil,
              message: "We did not see a corresponding option# for the default: «#{branch['default']}», please ensure the text matches one of the branches.",
              code: 'GTN:044',
              fn: __method__.to_s,
            )
          end
        end
      end

      # Branch testing.
      cyoa_branches.each do |branch|
        options = branch
          .select{|k, v| k =~ /option/}
          .values

        # Check for matching lines in the file.
        options.each do |option|
          slug_option = unsafe_slugify(option)
          if !joined_contents.match(/#{slug_option}/)
            errors << ReviewDogEmitter.warning(
              path: @path,
              idx: 0,
              match_start: 0,
              match_end: 1,
              replacement: nil,
              message: "We did not see a branch for #{option} (#{slug_option}) in the file. Please consider ensuring that all options are used.",
              code: 'GTN:045',
              fn: __method__.to_s,
            )
          end
        end
      end

      # find_matching_texts(contents, />\s*(\*\*\s*[Ss]tep)/) .map do |idx, _text, selected|
      #   ReviewDogEmitter.error(
      #     path: @path,
      #     idx: idx,
      #     match_start: selected.begin(1),
      #     match_end: selected.end(1) + 1,
      #     replacement: nil,
      #     message: 'This is a non-semantic list which is bad for accessibility and bad for screenreaders. ' \
      #              'It results in poorly structured HTML and as a result is not allowed.',
      #     code: 'GTN:035'
      #   )
      # end
      errors
    end

    ##
    # GTN:046 - Please do not add an # Introduction section, as it is unnecessary, please start directly into an abstract or hook for your tutorial that will get the learner interested in the material.
    def self.useless_intro(contents)
      joined_contents = contents.join("\n")
      joined_contents.scan(/\n---\n+# Introduction/m)
        .map do |line|
        ReviewDogEmitter.error(
          path: @path,
          idx: 0,
          match_start: 0,
          match_end: 0,
          replacement: '',
          message: 'Please do not include an # Introduction section, it is unnecessary here, just start directly into your text. The first paragraph that is seen by our infrastructure will automatically be shown in a few places as an abstract.',
          code: 'GTN:046',
          fn: __method__.to_s,
        )
      end
    end

    def self.fix_md(contents)
      [
        *fix_notoc(contents),
        *youtube_bad(contents),
        *link_gtn_slides_external(contents),
        *link_gtn_tutorial_external(contents),
        *check_dois(contents),
        *check_pmids(contents),
        *check_bad_link_text(contents),
        *incorrect_calls(contents),
        *check_bad_cite(contents),
        *non_existent_snippet(contents),
        *bad_tool_links(contents),
        *check_tool_link(contents),
        *new_more_accessible_boxes(contents),
        *new_more_accessible_boxes_agenda(contents),
        *no_target_blank(contents),
        *check_bad_link(contents),
        *check_bad_icon(contents),
        *check_looks_like_heading(contents),
        *check_bad_tag(contents),
        *check_useless_box_prefix(contents),
        *check_bad_heading_order(contents),
        *check_bolded_heading(contents),
        *snippets_too_close_together(contents),
        *bad_zenodo_links(contents),
        *zenodo_api(contents),
        *empty_alt_text(contents),
        *check_bad_trs_link(contents),
        *nonsemantic_list(contents),
        *cyoa_branches(contents),
        *useless_intro(contents)
      ]
    end

    def self.bib_missing_mandatory_fields(bib)
      results = []
      bib.each do |x|
        begin
          doi = x.doi
        rescue StandardError
          doi = nil
        end

        begin
          url = x.url
        rescue StandardError
          url = nil
        end

        begin
          isbn = x.isbn
        rescue StandardError
          isbn = nil
        end

        results.push([x.key, 'Missing a DOI, URL or ISBN. Please add one of the three.']) if doi.nil? && url.nil? && isbn.nil?

        begin
          x.title
          results.push([x.key, 'This entry is missing a title attribute. Please add it.']) if !x.title
        rescue StandardError
          results.push([x.key, 'This entry is missing a title attribute. Please add it.'])
        end
      end
      results
    end

    ##
    # GTN:015, GTN:016, GTN:025, GTN:026, GTN:027 others.
    # These error messages indicate something is amiss with your workflow. Please consult the error message to correct it.
    def self.fix_ga_wf(contents)
      results = []
      if !contents.key?('tags') or contents['tags'].empty?
        path_parts = @path.split('/')
        topic = path_parts[path_parts.index('topics') + 1]

        results.push(ReviewDogEmitter.file_error(
                       path: @path, message: "This workflow is missing required tags. Please add `\"tags\": [\"#{topic}\"]`",
                       code: 'GTN:015',
                       fn: __method__.to_s,
                     ))
      end

      if !contents.key?('annotation')
        results.push(ReviewDogEmitter.file_error(
                       path: @path,
                       message: 'This workflow is missing an annotation. Please add `"annotation": "title of tutorial"`',
                       code: 'GTN:016',
                       fn: __method__.to_s,
                     ))
      end

      if !contents.key?('license')
        results.push(ReviewDogEmitter.file_error(
                       path: @path,
                       message: 'This workflow is missing a license. Please select a valid OSI license. ' \
                                'You can correct this in the Galaxy workflow editor.',
                       code: 'GTN:026',
                       fn: __method__.to_s,
                     ))
      end

      tool_ids = tool_id_extractor(contents)

      # Check if they use TS tools, we do this here because it's easier to look at the plain text.
      tool_ids.each do |step_id, id|
        if ! acceptable_tool?(id)
          results += [
            ReviewDogEmitter.error(
              path: @path,
              idx: 0,
              match_start: 0,
              match_end: 0,
              replacement: nil,
              message: "A step in your workflow (#{step_id}) uses an invalid tool ID (#{id}) or a tool ID from the testtoolshed. These are not permitted in GTN tutorials. If this is in error, you can add it to the top of _plugins/utils.rb",
              code: 'GTN:017',
              fn: __method__.to_s,
            )
          ]
        end
      end




      if contents.key?('creator')
        contents['creator']
          .select { |c| c['class'] == 'Person' }
          .each do |p|
            if !p.key?('identifier') || (p['identifier'] == '')
              results.push(ReviewDogEmitter.file_error(
                             path: @path,
                             message: 'This workflow has a creator but is missing an identifier for them. ' \
                                      'Please ensure all creators have valid ORCIDs.',
                             code: 'GTN:025'
                           ))
            end

            if !p.key?('name') || (p['name'] == '')
              results.push(ReviewDogEmitter.file_error(
                             path: @path, message: 'This workflow has a creator but is a name, please add it.',
                             code: 'GTN:025'
                           ))
            end
          end
      else
        results.push(ReviewDogEmitter.file_error(
                       path: @path,
                       message: 'This workflow is missing a Creator. Please edit this workflow in ' \
                                'Galaxy to add the correct creator entities',
                       code: 'GTN:024'
                     ))
      end
      results
    end

    ##
    # GTN:012 - Your bibliography is missing mandatory fields (either a URL or DOI).
    # GTN:031 - Your bibliography unnecessarily fills the DOI field with https://doi.org, you can just directly specify the DOI.
    def self.fix_bib(contents, bib)
      bad_keys = bib_missing_mandatory_fields(bib)
      results = []
      bad_keys.each do |key, reason|
        results += find_matching_texts(contents, /^\s*@.*{#{key},/)
                   .map do |idx, text, _selected|
          ReviewDogEmitter.error(
            path: @path,
            idx: idx,
            match_start: 0,
            match_end: text.length,
            replacement: nil,
            message: reason,
            code: 'GTN:012',
            fn: __method__.to_s,
          )
        end
      end

      # 13:  doi = {https://doi.org/10.1016/j.cmpbup.2021.100007},
      results += find_matching_texts(contents, %r{doi\s*=\s*\{(https?://doi.org/)})
                 .map do |idx, _text, selected|
        ReviewDogEmitter.warning(
          path: @path,
          idx: idx,
          match_start: selected.begin(1),
          match_end: selected.end(1) + 1,
          replacement: '',
          message: 'Unnecessary use of URL in DOI-only field, please just use the doi component itself',
          code: 'GTN:031'
        )
      end
      results
    end

    @PLAIN_OUTPUT = false

    def self.set_plain_output
      @PLAIN_OUTPUT = true
    end

    def self.set_rdjson_output
      @PLAIN_OUTPUT = false
    end

    @SHORT_PATH = false
    def self.set_short_path
      @SHORT_PATH = true
    end

    @LIMIT_EMITTED_CODES = nil
    def self.code_limits(codes)
      @LIMIT_EMITTED_CODES = codes
    end

    @AUTO_APPLY_FIXES = false
    def self.enable_auto_fix
      @AUTO_APPLY_FIXES = true
    end

    def self.format_reviewdog_output(message)
      return if message.nil? || message.empty?
      return if !@LIMIT_EMITTED_CODES.nil? && !@LIMIT_EMITTED_CODES.include?(message['code']['value'])


      if !message.nil? && (message != []) && message.is_a?(Hash)
        path = message['location']['path']
        if @SHORT_PATH && path.include?(GTN_HOME + '/')
          path = path.gsub(GTN_HOME + '/', '')
        end
        if @PLAIN_OUTPUT # $stdout.tty? or
          parts = [
            path,
            message['location']['range']['start']['line'],
            message['location']['range']['start']['column'],
            message['location']['range']['end']['line'],
            message['location']['range']['end']['column'],
            "#{message['code']['value'].gsub(/:/, '')} #{message['message'].split("\n")[0]}"
          ]
          puts parts.join(':')
        else
          puts JSON.generate(message)
        end
      end

      return unless @AUTO_APPLY_FIXES && message['suggestions'].length.positive?

      start_line = message['location']['range']['start']['line']
      start_coln = message['location']['range']['start']['column']
      end_line = message['location']['range']['end']['line']
      end_coln = message['location']['range']['end']['column']

      if start_line == end_line
        # We only really support single-line changes. This will probs fuck up
        lines = File.read(message['location']['path']).split("\n")
        original = lines[start_line - 1].dup

        repl = message['suggestions'][0]['text']

        # puts "orig #{original}"
        # puts "before #{original[0..start_coln - 2]}"
        # puts "selected '#{original[start_coln-1..end_coln-2]}'"
        # puts "after #{original[end_coln-2..-1]}"
        # puts "replace: #{repl}"

        # puts "#{original[0..start_coln - 2]} + #{repl} + #{original[end_coln-1..-1]}"
        fixed = original[0..start_coln - 2] + repl + original[end_coln - 1..]
        warn "DIFF\n-#{original}\n+#{fixed}"
        lines[start_line - 1] = fixed

        # Save our changes
        File.write(message['location']['path'], (lines + ['']).join("\n"))
      else
        warn 'Cannot apply this suggestion sorry'
      end
    end

    def self.emit_results(results)
      return unless !results.nil? && results.length.positive?

      results.compact.flatten
        .select{|r| r.is_a? Hash }
        .each { |r| format_reviewdog_output(r) }
    end

    def self.should_ignore(contents)
      contents.select { |x| x.match(/GTN:IGNORE:(\d\d\d)/) }.map { |x| "GTN:#{x.match(/GTN:IGNORE:(\d\d\d)/)[1]}" }.uniq
    end

    def self.filter_results(results, ignores)
      if !results.nil?
        # Remove any empty lists
        results = results.select { |x| !x.nil? && x.length.positive? }.flatten
        # Before ignoring anything matching GTN:IGNORE:###
        return results if ignores.nil? or ignores.empty?

        results = results.select { |x| ignores.index(x['code']['value']).nil? } if results.length.positive?
        return results
      end
      nil
    end

    def self.fix_file(path)
      @path = path

      if path.match(/\s/)
        emit_results([ReviewDogEmitter.file_error(path: path,
                                                  message: 'There are spaces in this filename, that is forbidden.',
                                                  code: 'GTN:014')])
      end

      if path.match(/\?/)
        emit_results([ReviewDogEmitter.file_error(path: path,
                                                  message: 'There ?s in this filename, that is forbidden.',
                                                  code: 'GTN:014')])
      end

      case path
      when /md$/
        handle = File.open(path, 'r')
        contents = handle.read.split("\n")
        ignores = should_ignore(contents)
        results = fix_md(contents)

        results = filter_results(results, ignores)
        emit_results(results)
      when /.bib$/
        handle = File.open(path, 'r')
        contents = handle.read.split("\n")

        bib = BibTeX.open(path)
        results = fix_bib(contents, bib)

        results = filter_results(results, ignores)
        emit_results(results)
      when /.ga$/
        handle = File.open(path, 'r')
        begin
          contents = handle.read
          data = JSON.parse(contents)
        rescue StandardError => e
          warn "Error parsing #{path}: #{e}"
          emit_results([ReviewDogEmitter.file_error(path: path, message: 'Unparseable JSON in this workflow file.',
                                                    code: 'GTN:019')])
        end

        results = []
        # Check if there's a missing workflow test
        folder = File.dirname(path)
        basename = File.basename(path).gsub(/.ga$/, '')
        possible_tests = Dir.glob("#{folder}/#{Regexp.escape(basename)}*ym*")
        possible_tests = possible_tests.grep(/#{Regexp.escape(basename)}[_-]tests?.ya?ml/)

        contains_interactive_tool = contents.match(/interactive_tool_/)

        if possible_tests.empty?
          if !contains_interactive_tool
            results += [
              ReviewDogEmitter.file_error(path: path,
                                          message: 'This workflow is missing a test, which is now mandatory. Please ' \
                                                   'see [the FAQ on how to add tests to your workflows](' \
                                                   'https://training.galaxyproject.org/training-material/faqs/' \
                                                   'gtn/gtn_workflow_testing.html).',
                                          code: 'GTN:027')
            ]
          end
        else
          # Load tests and run some quick checks:
          possible_tests.each do |test_file|
            if !test_file.match(/-tests.yml/)
              results += [
                ReviewDogEmitter.file_error(path: path,
                                            message: 'Please use the extension -tests.yml ' \
                                                     'for this test file.',
                                            code: 'GTN:032')
              ]
            end

            test = YAML.safe_load(File.open(test_file))
            test_plain = File.read(test_file)
            # check that for each test, the outputs is non-empty
            unless test.is_a?(Array)
              next
            end
            test.each do |test_job|
              if (test_job['outputs'].nil? || test_job['outputs'].empty?) && !test_plain.match(/GTN_RUN_SKIP_REASON/)
                results += [
                  ReviewDogEmitter.file_error(path: path,
                                              message: 'This workflow test does not test the contents of outputs, ' \
                                                       'which is now mandatory. Please see [the FAQ on how to add ' \
                                                       'tests to your workflows](' \
                                                       'https://training.galaxyproject.org/training-material/faqs/' \
                                                       'gtn/gtn_workflow_testing.html).',
                                              code: 'GTN:030')
                ]
              end
            end
          end

        end

        results += fix_ga_wf(data)

        results = filter_results(results, ignores)
        emit_results(results)
      end
    end

    def self.enumerate_type(filter, root_dir: 'topics')
      paths = []
      Find.find("./#{root_dir}") do |path|
        if FileTest.directory?(path)
          next unless File.basename(path).start_with?('.')

          Find.prune       # Don't look any further into this directory.

        elsif path.match(filter)
          paths.push(path)
        end
      end
      paths
    end

    def self.enumerate_symlinks
      paths = []
      Find.find('./topics') do |path|
        if FileTest.directory?(path)
          next unless File.basename(path).start_with?('.')

          Find.prune       # Don't look any further into this directory.

        elsif File.symlink?(path)
          paths.push(path)
        end
      end
      paths
    end

    def self.enumerate_lintable
      enumerate_type(/bib$/) + enumerate_type(/md$/) + enumerate_type(/md$/,
                                                                      root_dir: 'faqs') + enumerate_type(/md$/,
                                                                                                         root_dir: 'news')
    end

    def self.enumerate_all
      enumerate_type(/.*/)
    end

    ##
    # GTN:014 - please do not use : colon in your filename.
    # GTN:013 - Please fix this symlink
    # GTN:023 - data libraries must be named data-library.yaml
    def self.run_linter_global
      enumerate_type(/:/).each do |path|
        format_reviewdog_output(
          ReviewDogEmitter.file_error(path: path,
                                      message: 'There are colons in this filename, that is forbidden.', 
                                      code: 'GTN:014',
                                      fn: __method__.to_s,
                                     )
        )
      end

      enumerate_symlinks.each do |path|
        if !File.exist?(Pathname.new(path).realpath)
          format_reviewdog_output(
            ReviewDogEmitter.file_error(path: path, message: 'This is a BAD symlink', 
                                        code: 'GTN:013',
                                        fn: __method__.to_s)
          )
        end
      rescue StandardError
        format_reviewdog_output(
          ReviewDogEmitter.file_error(path: path, message: 'This is a BAD symlink', code: 'GTN:013',
                                      fn: __method__.to_s)
        )
      end
      enumerate_type(/data[_-]library.ya?ml/).each do |path|
        if path.split('/')[-1] != 'data-library.yaml'
          format_reviewdog_output(
            ReviewDogEmitter.file_error(path: path,
                                        message: 'This file must be named data-library.yaml. Please rename it.',
                                        code: 'GTN:023')
          )
        end
      end
      enumerate_type(/\.ga$/).each do |path|
        fix_file(path)
      end
      enumerate_lintable.each do |path|
        fix_file(path)
      end
    end
  end
end

if $PROGRAM_NAME == __FILE__
  linter = Gtn::Linter

  require 'optparse'
  require 'ostruct'

  options = {}
  OptionParser.new do |opt|
    # Mutually exclusive
    opt.on('-f', '--format [plain|rdjson]', 'Preferred output format, defaults to plain') { |o| options[:format] = o }
    opt.on('-p', '--path file.md', 'Specify a single file to check instead of the entire repository') do |o|
      options[:path] = o
    end
    opt.on('-l', '--limit GTN:001,...', 'Limit output to specific codes') { |o| options[:limit] = o }
    opt.on('-a', '--auto-fix', 'I am not sure this is really safe, be careful') { |_o| options[:apply] = true }
    opt.on('-s', '--short-path', 'Use short path in outputs') { |_o| options[:short] = true }
  end.parse!

  options[:format] = 'plain' if options[:format].nil?

  if options[:format] == 'plain'
    linter.set_plain_output
  else
    linter.set_rdjson_output
  end

  linter.set_short_path if options[:short]
  linter.code_limits(options[:limit].split(',')) if options[:limit]

  linter.enable_auto_fix if options[:apply]

  if options[:path].nil?
    linter.run_linter_global
  else
    linter.fix_file(options[:path])
  end
end
