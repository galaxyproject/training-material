#!/usr/bin/env ruby
require 'pathname'
require 'find'
require 'bibtex'
require 'json'
require 'find'
require 'bibtex'
require 'citeproc/ruby'
require 'csl/styles'

# Load our citation library
CITATION_LIBRARY = BibTeX::Bibliography.new
Find.find('./') do |path|
  if FileTest.directory?(path)
    if File.basename(path).start_with?('.')
      Find.prune       # Don't look any further into this directory.
    else
      next
    end
  else
    if path =~ /bib$/ then
      b = BibTeX.open(path)
      for x in b
        # Record the bib path.
        x._path = path
        CITATION_LIBRARY << x
      end
    end
  end
end



# This is our ONE central linting script that handles EVERYTHING.

module ReviewDogEmitter

  @CODE_URL = "https://github.com/galaxyproject/training-material/wiki/Error-Codes"
  def self.delete_text(path: "", idx: 0, text: "", message: "No message", code: "GTN000", full_line: "")
    self.error(
      path: path,
      idx: idx,
      match_start: 0,
      match_end: text.length,
      replacement: "",
      message: message,
      code: code,
      full_line: full_line,
    )
  end

  def self.file_error(path: "", message: "None", code: "GTN:000")
    self.error(
      path: path,
      idx: 0,
      match_start: 0,
      match_end: 1,
      replacement: nil,
      message: message,
      code: code,
      full_line: "",
    )
  end

  def self.warning(path: "", idx: 0, match_start: 0, match_end: 1, replacement: nil, message: "No message", code: "GTN000", full_line: "")
    self.message(
      path: path,
      idx: idx,
      match_start: match_start,
      match_end: match_end,
      replacement: replacement,
      message: message,
      level:"WARNING",
      code: code,
      full_line: full_line,
    )
  end

  def self.error(path: "", idx: 0, match_start: 0, match_end: 1, replacement: nil, message: "No message", code: "GTN000", full_line: "")
    self.message(
      path: path,
      idx: idx,
      match_start: match_start,
      match_end: match_end,
      replacement: replacement,
      message: message,
      level:"ERROR",
      code: code,
      full_line: full_line,
    )
  end

  def self.message(path: "", idx: 0, match_start: 0, match_end: 1, replacement: nil, message: "No message", level: "WARNING", code: "GTN000", full_line: "")
    end_area = { "line" => idx + 1, "column" => match_end}
    if match_end == full_line.length
      end_area = { "line" => idx + 2, "column" => 1}
    end

    res = {
      "message" => message,
      'location' => {
        'path' => path,
        'range' => {
          'start' => { "line" => idx + 1, "column" => match_start + 1},
          'end' => end_area
        }
      },
      "severity" => level
    }
    if !code.nil?
      res["code"] = {
        "value" => code,
        "url" => @CODE_URL + "#" + code.gsub(/:/, '').downcase,
      }
    end
    if !replacement.nil?
      res['suggestions'] = [{
        'text' => replacement,
        'range' => {
          'start' => { "line" => idx + 1, "column" => match_start + 1 },
          'end' => end_area
        }
      }]
    end
    res
  end
end

module GtnLinter
  @BAD_TOOL_LINK = /{% tool (\[[^\]]*\])\(https?.*tool_id=([^)]*)\)\s*%}/i

  def self.find_matching_texts(contents, query)
    contents.map.with_index { |text, idx|
      [idx, text, text.match(query)]
    }.select { |idx, text, selected| selected }
  end

  def self.fix_notoc(contents)
    self.find_matching_texts(contents, /{:\s*.no_toc\s*}/)
        .map { |idx, text, selected|
      ReviewDogEmitter.delete_text(
        path: @path,
        idx: idx,
        text: text,
        message: "Setting  is discouraged, these headings provide useful places for readers to jump to.",
        code: "GTN:001",
        full_line: text,
      )
    }
  end

  # GTN:002 youtube discouraged
  def self.youtube_bad(contents)
    self.find_matching_texts(contents, /<iframe.*youtu.?be.*<\/iframe>/)
        .map { |idx, text, selected|
      ReviewDogEmitter.warning(
        path: @path,
        idx: idx,
        match_start: selected.begin(0),
        match_end: selected.end(0) + 1,
        replacement: "",
        message: "Instead of embedding IFrames to YouTube contents, consider adding this video to the [GTN Video Library](https://github.com/gallantries/video-library/issues/) where it will be more visible for others.",
        code: "GTN:002"
      )
    }
  end

  def self.link_gtn_tutorial_external(contents)
    self.find_matching_texts(
      contents,
      /\((https?:\/\/(training.galaxyproject.org|galaxyproject.github.io)\/training-material\/(.*tutorial).html)\)/
    )
    .map { |idx, text, selected|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        # We wrap the entire URL (inside the explicit () in a matching group to make it easy to select/replace)
        match_start: selected.begin(1),
        match_end: selected.end(1) + 1,
        replacement: "{% link #{selected[3]}.md %}",
        message: "Please use the link function to link to other pages within the GTN. It helps us ensure that all links are correct",
        code: "GTN:003",
      )
    }
  end

  def self.link_gtn_slides_external(contents)
    self.find_matching_texts(
      contents,
      /\((https?:\/\/(training.galaxyproject.org|galaxyproject.github.io)\/training-material\/(.*slides.html))\)/
    )
    .map { |idx, text, selected|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(1),
        match_end: selected.end(1) + 1,
        replacement:"{% link #{selected[3]} %}",
        message: "Please use the link function to link to other pages within the GTN. It helps us ensure that all links are correct",
        code: "GTN:003",
      )
    }
  end

  def self.check_dois(contents)
    self.find_matching_texts(contents, /(\[[^\]]*\]\(https?:\/\/doi.org\/[^)]*\))/)
      .select{|idx, text, selected| ! selected[0].match(/10.5281\/zenodo/) } # Ignoring zenodo
        .map { |idx, text, selected|
      ReviewDogEmitter.warning(
        path: @path,
        idx: idx,
        match_start: selected.begin(0),
        match_end: selected.end(0) + 2,
        replacement: "{% cite ... %}",
        message: "This looks like a DOI which could be better served by using the built-in Citations mechanism. You can use https://doi2bib.org to convert your DOI into a .bib formatted entry, and add to your tutorial.md",
        code: "GTN:004"
      )
    }
  end

  def self.check_bad_link_text(contents)
    self.find_matching_texts(contents, /\[\s*(here|link)\s*\]/i)
        .map { |idx, text, selected|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(0),
        match_end: selected.end(0) + 1,
        replacement: "[Something better here]",
        message: "Please do not use 'here' as your link title, it is " +
         "[bad for accessibility](https://usability.yale.edu/web-accessibility/articles/links#link-text). " +
         "Instead try restructuring your sentence to have useful descriptive text in the link.",
        code: "GTN:005",
      )
    }
  end

  def self.incorrect_calls(contents)
    a = self.find_matching_texts(contents, /([^{]|^)(%\s*[^%]*%})/i)
            .map { |idx, text, selected|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(2),
        match_end: selected.end(2) + 1,
        replacement: "{#{selected[2]}",
        message: "It looks like you might be missing the opening { of a jekyll function",
        code: "GTN:006",
      )
    }
    b = self.find_matching_texts(contents, /{([^%]\s*[^%]* %})/i)
            .map { |idx, text, selected|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(1),
        match_end: selected.end(1) + 1,
        replacement: "%#{selected[1]}",
        message: "It looks like you might be missing the opening % of a jekyll function",
        code: "GTN:006",
      )
    }

    c = self.find_matching_texts(contents, /({%\s*[^%]*%)([^}]|$)/i)
            .map { |idx, text, selected|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(1),
        match_end: selected.end(1) + 2,
        replacement: "#{selected[1]}}#{selected[2]}",
        message: "It looks like you might be missing the closing } of a jekyll function",
        code: "GTN:006",
      )
    }

    d = self.find_matching_texts(contents, /({%\s*[^}]*[^%])}/i)
            .map { |idx, text, selected|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(1),
        match_end: selected.end(1) + 1,
        replacement: "#{selected[1]}%",
        message: "It looks like you might be missing the closing % of a jekyll function",
        code: "GTN:006",
      )
    }
    a + b + c + d
  end

  def self.check_bad_cite(contents)
    self.find_matching_texts(contents, /{%\s*cite\s+([^%]*)\s*%}/i)
    .map { |idx, text, selected|
      if CITATION_LIBRARY[selected[1].strip].nil?
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(0),
          match_end: selected.end(0),
          replacement: nil,
          message: "The citation (`#{selected[1].strip}`) could not be found.",
          code: "GTN:007",
        )
      end
    }
  end

  def self.non_existent_snippet(contents)
    self.find_matching_texts(contents, /{%\s*snippet\s+([^ ]*)/i)
        .select { |idx, text, selected|
      !File.exists?(selected[1])
    }
        .map { |idx, text, selected|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(0),
        match_end: selected.end(0),
        replacement: nil,
        message: "This snippet (`#{selected[1]}`) does not seem to exist",
        code: "GTN:008",
      )
    }
  end

  def self.bad_tool_links(contents)
    self.find_matching_texts(contents, @BAD_TOOL_LINK)
        .map { |idx, text, selected|
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(0),
        match_end: selected.end(0) + 1,
        replacement: "{% tool #{selected[1]}(#{selected[2]}) %}",
        message: "You have used the full tool URL to a specific server, here we only need the tool ID portion.",
        code: "GTN:009",
      )
    }
  end

  ALLOWED_SHORT_IDS = [
      'ChangeCase',
      'Convert characters1',
      'Count1',
      'Cut1',
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
      'intermine',
      'join1',
      'param_value_from_file',
      'random_lines1',
      'sort1',
      'ucsc_table_direct1',
      'upload1',
      'wc_gnu',
      'wig_to_bigWig',
  ]

  def self.check_tool_link(contents)
    self.find_matching_texts(contents, /{%\s*tool \[([^\]]*)\]\(([^)]*)\)\s*%}/)
      .map { |idx, text, selected|
        text = selected[1]
        link = selected[2]

        errs = []

        if link.match(/\//)
          if link.count('/') < 5
            errs.push(ReviewDogEmitter.error(
              path: @path,
              idx: idx,
              match_start: selected.begin(2),
              match_end: selected.end(2) + 1,
              replacement: nil,
              message: "This tool identifier looks incorrect, it doesn't have the right number of segments.",
              code: "GTN:009",
            ))
          end

          if link.match(/testtoolshed/)
            errs.push(ReviewDogEmitter.warning(
              path: @path,
              idx: idx,
              match_start: selected.begin(2),
              match_end: selected.end(2) + 1,
              replacement: nil,
              message: "The GTN strongly avoids using testtoolshed tools in your tutorials or workflows",
              code: "GTN:009",
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
              message: "Broken tool link, unnecessary +",
              code: "GTN:009",
            ))
          end

          if not ALLOWED_SHORT_IDS.include?(link) and not link.match(/^interactive_tool_/) and not link.match(/__[A-Z_]+__/)
            if not link.match(/^{{.*}}$/)
              errs.push(ReviewDogEmitter.error(
                path: @path,
                idx: idx,
                match_start: selected.begin(2),
                match_end: selected.end(2) + 1,
                replacement: nil,
                message: "Unknown short tool ID. Please use the full tool ID, or check bin/lint.rb if you believe this is correct.",
                code: "GTN:009",
              ))
            end
          end
        end

        errs
    }
  end

  def self.new_more_accessible_boxes(contents)
    #  \#\#\#
    self.find_matching_texts(contents, /> (### {%\s*icon ([^%]*)\s*%}[^:]*:?(.*))/)
        .map { |idx, text, selected|
      key = selected[2].strip.gsub(/_/, '-')
      ReviewDogEmitter.error(
        path: @path,
        idx: idx,
        match_start: selected.begin(1),
        match_end: selected.end(1) + 1,
        replacement: "<#{key}-title>#{selected[3].strip}</#{key}-title>",
        message: "We have developed a new syntax for box titles, please consider using this instead.",
        code: "GTN:010",
      )
    }
  end

  def self.no_target_blank(contents)
    self.find_matching_texts(contents, /target=("_blank"|'_blank')/)
        .map { |idx, text, selected|
      ReviewDogEmitter.warning(
        path: @path,
        idx: idx,
        match_start: selected.begin(0),
        match_end: selected.end(0),
        replacement: nil,
        message: "Please do not use `target=\"_blank\"`, [it is bad for accessibility.](https://www.a11yproject.com/checklist/#identify-links-that-open-in-a-new-tab-or-window)",
        code: "GTN:011",
      )
    }
  end

  def self.check_bad_link(contents)
    self.find_matching_texts(contents, /{%\s*link\s+([^%]*)\s*%}/i)
    .map { |idx, text, selected|
      path = selected[1].to_s.strip
      if ! File.exist?(path.gsub(/^\//, ''))
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: selected.begin(0),
          match_end: selected.end(0),
          replacement: nil,
          message: "The linked file (`#{selected[1].strip}`) could not be found.",
          code: "GTN:018",
        )
      end
    }
  end

  def self.check_looks_like_heading(contents)
    # TODO: we should remove this someday, but, we need to have a good solution
    # and we're still a ways from that.
    #
    # There's no clear way to say "this subsection of the content has its own hierarchy"
    if @path.match(/faq/)
      return
    end
    self.find_matching_texts(contents, /^\*\*(.*)\*\*$/)
    .map { |idx, text, selected|
      ReviewDogEmitter.warning(
        path: @path,
        idx: idx,
        match_start: selected.begin(1),
        match_end: selected.end(1) + 1,
        replacement: "### #{selected[1]}",
        message: "This looks like a heading, but isn't. Please use proper semantic headings where possible. You should check the heading level of this suggestion, rather than accepting the change as-is.",
        code: "GTN:020",
      )
    }
  end

  def self.fix_md(contents)
    [
      *fix_notoc(contents),
      *youtube_bad(contents),
      *link_gtn_slides_external(contents),
      *link_gtn_tutorial_external(contents),
      *check_dois(contents),
      *check_bad_link_text(contents),
      *incorrect_calls(contents),
      #*check_bad_cite(contents),
      *non_existent_snippet(contents),
      *bad_tool_links(contents),
      *check_tool_link(contents),
      *new_more_accessible_boxes(contents),
      *no_target_blank(contents),
      *check_bad_link(contents),
      *check_looks_like_heading(contents),
    ]
  end

  def self.bib_missing_mandatory_fields(bib)
    results = []
    for x in bib
      begin
        doi = x.doi
      rescue
        doi = nil
      end

      begin
        url = x.url
      rescue
        url = nil
      end

      if doi.nil? && url.nil?
        results.push([x.key, "Missing both a DOI and a URL. Please add one of the two."])
      end

      begin
        x.title
        if !x.title
          results.push([x.key, "This entry is missing a title attribute. Please add it."])
        end
      rescue
        results.push([x.key, "This entry is missing a title attribute. Please add it."])
      end
    end
    return results
  end

  def self.fix_ga_wf(contents)
    results = []
    if ! contents.has_key?("tags")
      topic = @path.split('/')[1]
      results.push(ReviewDogEmitter.file_error(
        path: @path, message: "This workflow is missing tags. Please add `\"tags\": [\"#{topic}\"]`", code: "GTN:015"))
    end

    if ! contents.has_key?("annotation")
      topic = @path.split('/')[1]
      results.push(ReviewDogEmitter.file_error(
        path: @path, message: "This workflow is missing an annotation. Please add `\"annotation\": \"title of tutorial\"`", code: "GTN:016"))
    end

    results
  end

  def self.fix_bib(contents, bib)
    bad_keys = bib_missing_mandatory_fields(bib)
    results = []
    bad_keys.each { |key, reason|
      results += self.find_matching_texts(contents, /^\s*@.*{#{key},/)
                     .map { |idx, text, selected|
        ReviewDogEmitter.error(
          path: @path,
          idx: idx,
          match_start: 0,
          match_end: text.length,
          replacement:  nil,
          message: reason,
          code: "GTN:012",
        )
      }
    }
    results
  end

  @PLAIN_OUTPUT = false

  def self.set_plain_output
    @PLAIN_OUTPUT = true
  end

  def self.format_reviewdog_output(message)
    if !message.nil? and message != []
      if $stdout.tty? or @PLAIN_OUTPUT
        parts = [
          message['location']['path'],
          message['location']['range']['start']['line'],
          message['location']['range']['start']['column'],
          message['location']['range']['end']['line'],
          message['location']['range']['end']['column'],
          "#{message['code']['value'].gsub(/:/, '')} #{message['message']}"
        ]
        puts parts.join(":")
      else
        puts JSON.generate(message)
      end
    end
  end

  def self.emit_results(results)
    if ! results.nil? and results.length > 0
      results.select{|r| !r.nil? }.flatten.each { |r| format_reviewdog_output(r) }
    end
  end

  def self.should_ignore(contents)
    contents.select{|x| x.match(/GTN:IGNORE:(\d\d\d)/)}.map{|x| "GTN:" + x.match(/GTN:IGNORE:(\d\d\d)/)[1] }.uniq
  end

  def self.fix_file(path)
    @path = path

    if path.match(/\s/)
      emit_results([ReviewDogEmitter.file_error(path: path, message: "There are spaces in this filename, that is forbidden.", code: "GTN:014")])
    end

    if path.match(/md$/)
      handle = File.open(path)
      contents = handle.read.split("\n")
      ignores = should_ignore(contents)
      results = fix_md(contents)

      if ! results.nil?
        # Remove any empty lists
        results = results.select{|x| !x.nil? && x.length > 0 }.flatten
        # Before ignoring anything matching GTN:IGNORE:###
        if results.length > 0
            results = results.select{|x| ignores.index(x['code']['value']).nil?}
        end
      end
      emit_results(results)
    elsif path.match(/.bib$/)
      handle = File.open(path)
      contents = handle.read.split("\n")

      bib = BibTeX.open(path)
      results = fix_bib(contents, bib)
      emit_results(results)
    elsif path.match(/.ga$/)
      handle = File.open(path)
      begin
        contents = handle.read
        data = JSON.parse(contents)
        contents.split("\n").each.with_index{|text, linenumber|
          if text.match(/testtoolshed/)
            emit_results([
              ReviewDogEmitter.error(
                path: @path,
                idx: linenumber,
                match_start: 0,
                match_end: text.length,
                replacement: nil,
                message: "This step uses a tool from the testtoolshed. These are not permitted in GTN tutorials.",
                code: "GTN:017",
              )
            ])
          end
        }
        results = fix_ga_wf(data)
        emit_results(results)
      rescue
        emit_results([ReviewDogEmitter.file_error(path: path, message: "Unparseable JSON in this workflow file.", code: "GTN:019")])
      end
      bib = BibTeX.open(path)
    end
  end

  def self.enumerate_type(filter, root_dir: "topics")
    paths = []
    Find.find("./#{root_dir}") do |path|
      if FileTest.directory?(path)
        if File.basename(path).start_with?('.')
          Find.prune       # Don't look any further into this directory.
        else
          next
        end
      else
        if path.match(filter)
          paths.push(path)
        end
      end
    end
    paths
  end

  def self.enumerate_symlinks
    paths = []
    Find.find('./topics') do |path|
      if FileTest.directory?(path)
        if File.basename(path).start_with?('.')
          Find.prune       # Don't look any further into this directory.
        else
          next
        end
      else
        if File.symlink?(path) then
          paths.push(path)
        end
      end
    end
    paths
  end

  def self.enumerate_lintable
    self.enumerate_type(/bib$/) + self.enumerate_type(/md$/) + self.enumerate_type(/md$/, root_dir: "faqs")
  end

  def self.run_linter_global
    self.enumerate_symlinks.each{|path|
      begin
        if ! File.exists?(Pathname.new(path).realpath)
          self.format_reviewdog_output(
            ReviewDogEmitter.file_error(path: path, message: "This is a BAD symlink", code: "GTN:013")
          )
        end
      rescue
          self.format_reviewdog_output(
            ReviewDogEmitter.file_error(path: path, message: "This is a BAD symlink", code: "GTN:013")
          )
      end
    }
    self.enumerate_type(/\.ga$/).each{|path|
      self.fix_file(path)
    }
    self.enumerate_lintable.each{|path|
      self.fix_file(path)
    }
  end
end

if $0 == __FILE__
  linter = GtnLinter

  if ARGV.length > 0
    if ARGV[0] == '--plain'
      linter.set_plain_output()
      linter.run_linter_global()
    else
      linter.fix_file(ARGV[0])
    end
  else
    linter.run_linter_global()
  end
end
