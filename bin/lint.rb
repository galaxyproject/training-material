#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'pathname'
require 'find'
require 'bibtex'
require 'json'
require 'citeproc/ruby'
require 'csl/styles'

# This is our ONE central linting script that handles EVERYTHING.

# A custom module to properly format reviewdog json output
module ReviewDogEmitter
  @CODE_URL = 'https://github.com/galaxyproject/training-material/wiki/Error-Codes'
  def self.delete_text(path: '', idx: 0, text: '', message: 'No message', code: 'GTN000', full_line: '')
    error(
      path: path,
      idx: idx,
      match_start: 0,
      match_end: text.length,
      replacement: '',
      message: message,
      code: code,
      full_line: full_line
    )
  end

  def self.file_error(path: '', message: 'None', code: 'GTN:000')
    error(
      path: path,
      idx: 0,
      match_start: 0,
      match_end: 1,
      replacement: nil,
      message: message,
      code: code,
      full_line: ''
    )
  end

  def self.warning(path: '', idx: 0, match_start: 0, match_end: 1,
                   replacement: nil, message: 'No message', code: 'GTN000', full_line: '')
    self.message(
      path: path,
      idx: idx,
      match_start: match_start,
      match_end: match_end,
      replacement: replacement,
      message: message,
      level: 'WARNING',
      code: code,
      full_line: full_line
    )
  end

  def self.error(path: '', idx: 0, match_start: 0, match_end: 1, replacement: nil, message: 'No message',
                 code: 'GTN000', full_line: '')
    self.message(
      path: path,
      idx: idx,
      match_start: match_start,
      match_end: match_end,
      replacement: replacement,
      message: message,
      level: 'ERROR',
      code: code,
      full_line: full_line
    )
  end

  def self.message(path: '', idx: 0, match_start: 0, match_end: 1, replacement: nil, message: 'No message',
                   level: 'WARNING', code: 'GTN000', full_line: '')
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
        'value' => code,
        'url' => "#{@CODE_URL}##{code.gsub(/:/, '').downcase}",
      }
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

# Linting functions for the GTN
module GtnLinter
  @BAD_TOOL_LINK = /{% tool (\[[^\]]*\])\(https?.*tool_id=([^)]*)\)\s*%}/i
  @BAD_TOOL_LINK2 = %r{{% tool (\[[^\]]*\])\(https://toolshed.g2([^)]*)\)\s*%}}i

  def self.find_matching_texts(contents, query)
    contents.map.with_index do |text, idx|
      [idx, text, text.match(query)]
    end.select { |_idx, _text, selected| selected }
  end

  def self.fix_notoc(contents)
    find_matching_texts(contents, /{:\s*.no_toc\s*}/)
      .map do |idx, text, _selected|
      ReviewDogEmitter.delete_text(
        path: @path,
        idx: idx,
        text: text,
        message: 'Setting no_toc is discouraged, these headings provide useful places for readers to jump to.',
        code: 'GTN:001',
        full_line: text
      )
    end
  end

  # GTN:002 youtube discouraged
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
                 '[GTN Video Library](https://github.com/gallantries/video-library/issues/) where it will ' \
                 'be more visible for others.',
        code: 'GTN:002'
      )
    end
  end

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
        code: 'GTN:003'
      )
    end
  end

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
        code: 'GTN:003'
      )
    end
  end

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
        code: 'GTN:004'
      )
    end
  end

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
        code: 'GTN:004'
      )
    end
  end

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
        code: 'GTN:005'
      )
    end
  end

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
        code: 'GTN:006'
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
        code: 'GTN:006'
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
        code: 'GTN:006'
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
        code: 'GTN:006'
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
          code: 'GTN:007'
        )
      end
    end
  end

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
        code: 'GTN:008'
      )
    end
  end

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
          code: 'GTN:009'
        )
      end
  end

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
                   code: 'GTN:032'
                 ))
      end
      prev_line = idx
    end
    res
  end

  ALLOWED_SHORT_IDS = [
    'ChangeCase',
    'Convert characters1',
    'Count1',
    'Cut1',
    'Extract_features1',
    'Filter1',
    'Grep1',
    'Grouping1',
    'Paste1',
    'Remove beginning1',
    'Show beginning1',
    'Summary_Statistics1',
    'addValue',
    'cat1',
    'comp1',
    'gene2exon1',
    'gff2bed1',
    'intermine',
    'join1',
    'param_value_from_file',
    'random_lines1',
    'sort1',
    # 'ucsc_table_direct1', # This does not work, surprisingly.
    'upload1',
    'wc_gnu',
    'wig_to_bigWig'
  ].freeze

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

        if !ALLOWED_SHORT_IDS.include?(link) &&
           !link.match(/^interactive_tool_/) &&
           !link.match(/__[A-Z_]+__/) &&
           !link.match(/^{{.*}}$/) &&
           !link.match(/^CONVERTER_/)
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
        code: 'GTN:010'
      )
    end
  end

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
        code: 'GTN:010'
      )
    end
  end

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
        code: 'GTN:011'
      )
    end
  end

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
          code: 'GTN:018'
        )
      end
    end
  end

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
        code: 'GTN:020'
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
        code: 'GTN:021'
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
        code: 'GTN:022'
      )
    end
  end

  def self.check_bad_heading_order(contents)
    depth = 1
    headings = find_matching_texts(contents, /^(?<level>#+)\s?(?<title>.*)/)
               .map do |idx, text, selected|
      new_depth = selected[:level].length
      depth_change = new_depth - depth
      depth = new_depth
      [idx, text, selected, depth_change, new_depth]
    end

    all_headings = headings.map do |_idx, _text, selected, _depth_change, _new_depth|
      "#{selected[:level]} #{selected[:title]}"
    end

    headings.select do |_idx, _text, _selected, depth_change, _new_depth|
      depth_change > 1
    end.map do |idx, _text, selected, depth_change, new_depth|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(1),
        match_end: selected.end(1) + 1,
        replacement: '#' * (new_depth - depth_change + 1),
        message: "You have skipped a heading level, please correct this.\n<details>" \
                 "<summary>Listing of Heading Levels</summary>\n\n```\n#{all_headings.join("\n")}\n```\n</details>",
        code: 'GTN:028'
      )
    end
  end

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
        code: 'GTN:029'
      )
    end
  end

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
      *check_looks_like_heading(contents),
      *check_bad_tag(contents),
      *check_useless_box_prefix(contents),
      *check_bad_heading_order(contents),
      *check_bolded_heading(contents),
      *snippets_too_close_together(contents),
      *zenodo_api(contents)
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

      results.push([x.key, 'Missing both a DOI and a URL. Please add one of the two.']) if doi.nil? && url.nil?

      begin
        x.title
        results.push([x.key, 'This entry is missing a title attribute. Please add it.']) if !x.title
      rescue StandardError
        results.push([x.key, 'This entry is missing a title attribute. Please add it.'])
      end
    end
    results
  end

  def self.fix_ga_wf(contents)
    results = []
    if !contents.key?('tags')
      topic = @path.split('/')[1]
      results.push(ReviewDogEmitter.file_error(
                     path: @path, message: "This workflow is missing tags. Please add `\"tags\": [\"#{topic}\"]`",
                     code: 'GTN:015'
                   ))
    end

    if !contents.key?('annotation')
      results.push(ReviewDogEmitter.file_error(
                     path: @path,
                     message: 'This workflow is missing an annotation. Please add `"annotation": "title of tutorial"`',
                     code: 'GTN:016'
                   ))
    end

    if !contents.key?('license')
      results.push(ReviewDogEmitter.file_error(
                     path: @path,
                     message: 'This workflow is missing a license. Please select a valid OSI license. ' \
                              'You can correct this in the Galaxy workflow editor.',
                     code: 'GTN:026'
                   ))
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
          code: 'GTN:012'
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

  @LIMIT_EMITTED_CODES = nil
  def self.code_limits(codes)
    @LIMIT_EMITTED_CODES = codes
  end

  @AUTO_APPLY_FIXES = false
  def self.enable_auto_fix
    @AUTO_APPLY_FIXES = true
  end

  def self.format_reviewdog_output(message)
    return if !@LIMIT_EMITTED_CODES.nil? && !@LIMIT_EMITTED_CODES.include?(message['code']['value'])

    if !message.nil? && (message != [])
      if @PLAIN_OUTPUT # $stdout.tty? or
        parts = [
          message['location']['path'],
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

    results.compact.flatten.each { |r| format_reviewdog_output(r) }
  end

  def self.should_ignore(contents)
    contents.select { |x| x.match(/GTN:IGNORE:(\d\d\d)/) }.map { |x| "GTN:#{x.match(/GTN:IGNORE:(\d\d\d)/)[1]}" }.uniq
  end

  def self.filter_results(results, ignores)
    if !results.nil?
      # Remove any empty lists
      results = results.select { |x| !x.nil? && x.length.positive? }.flatten
      # Before ignoring anything matching GTN:IGNORE:###
      return results if ignores.nil?

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

      if possible_tests.empty?
        results += [
          ReviewDogEmitter.file_error(path: path,
                                      message: 'This workflow is missing a test, which is now mandatory. Please ' \
                                               'see [the FAQ on how to add tests to your workflows](' \
                                               'https://training.galaxyproject.org/training-material/faqs/' \
                                               'gtn/gtn_workflow_testing.html).',
                                      code: 'GTN:027')
        ]
      else
        # Load tests and run some quick checks:
        possible_tests.each do |test_file|
          if !test_file.match(/-tests?.yml/)
            results += [
              ReviewDogEmitter.file_error(path: path,
                                          message: 'Please use the extension -test.yml ' \
                                                   'or -tests.yml for this test file.',
                                          code: 'GTN:032')
            ]
          end

          test = YAML.safe_load(File.open(test_file))
          # check that for each test, the outputs is non-empty
          test.each do |test_job|
            if test_job['outputs'].nil? || test_job['outputs'].empty?
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

      # Check if they use TS tools, we do this here because it's easier to look at the plain text.
      contents.split("\n").each.with_index do |text, linenumber|
        if text.match(/testtoolshed/)
          results += [
            ReviewDogEmitter.error(
              path: @path,
              idx: linenumber,
              match_start: 0,
              match_end: text.length,
              replacement: nil,
              message: 'This step uses a tool from the testtoolshed. These are not permitted in GTN tutorials.',
              code: 'GTN:017'
            )
          ]
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
    enumerate_type(/bib$/) + enumerate_type(/md$/) + enumerate_type(/md$/, root_dir: 'faqs')
  end

  def self.enumerate_all
    enumerate_type(/.*/)
  end

  def self.run_linter_global
    enumerate_type(/:/).each do |path|
      format_reviewdog_output(
        ReviewDogEmitter.file_error(path: path,
                                    message: 'There are colons in this filename, that is forbidden.', code: 'GTN:014')
      )
    end

    enumerate_symlinks.each do |path|
      if !File.exist?(Pathname.new(path).realpath)
        format_reviewdog_output(
          ReviewDogEmitter.file_error(path: path, message: 'This is a BAD symlink', code: 'GTN:013')
        )
      end
    rescue StandardError
      format_reviewdog_output(
        ReviewDogEmitter.file_error(path: path, message: 'This is a BAD symlink', code: 'GTN:013')
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

if $PROGRAM_NAME == __FILE__
  linter = GtnLinter

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
  end.parse!

  options[:format] = 'plain' if options[:format].nil?

  if options[:format] == 'plain'
    linter.set_plain_output
  else
    linter.set_rdjson_output
  end

  linter.code_limits(options[:limit].split(',')) if options[:limit]

  linter.enable_auto_fix if options[:apply]

  if options[:path].nil?
    linter.run_linter_global
  else
    linter.fix_file(options[:path])
  end
end
