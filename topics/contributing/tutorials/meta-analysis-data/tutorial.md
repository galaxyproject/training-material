---
layout: tutorial_hands_on
title: Single Cell Publication - Data Analysis
level: Advanced
requirements: []
follow_up_training:
  - type: "internal"
    topic_name: contributing
    tutorials:
      - meta-analysis-plot

questions:
objectives:
time_estimation:  1H
key_points:
contributions:
  authorship:
  - hexylena
subtopic: meta
#notebook:
#    language: ruby
priority: 3
tags:
- Ruby

---

Extracting data from the GTN's Git history isn't that difficult, but it requires some internal knowledge of how the GTN's Jekyll-based codebase works. Here we'll document what we've done!

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

Our imports and metadata (all merged github PRs), and the list of all historical names of single-cell tutorial folders.

The author recommends running this code in a 'Jekyll Console' context. Jekyll does not natively have support for a console, but [there is an open Pull Request to add it](https://github.com/jekyll/jekyll/pull/8164). We recommend you install this yourself to most easily run the following code. You can do that by:

> <hands-on-title>Installing Jekyll's Console</hands-on-title>
> 1. View [the open Pull Request to add it](https://github.com/jekyll/jekyll/pull/8164)
> 2. Download `lib/jekyll/commands/console.rb` to somewhere on your computer.
> 3. Find out where jekyll is installed:
>    
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > gem which jekyll
>    > ```
>    {: .code-in}
>    
>    > <code-out-title></code-out-title>
>    > ```bash
>    > /home/user/galaxy/training-material/.direnv/ruby/gems/jekyll-4.3.3/lib/jekyll.rb
>    > ```
>    {: .code-out}
>    
>    That means the commands should be in:
>    
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ls /home/user/galaxy/training-material/.direnv/ruby/gems/jekyll-4.3.3/lib/jekyll/commands/
>    > ```
>    {: .code-in}
>    
>    > <code-out-title></code-out-title>
>    > ```bash
>    > build.rb  clean.rb  doctor.rb  help.rb  new.rb  new_theme.rb  serve  serve.rb
>    > ```
>    {: .code-out}
> 
> 4. Copy the `console.rb` you downloaded to that folder
> 
>     ```bash
>     cp ~/Downloads/console.rb $(dirname $(gem which jekyll))/jekyll/commands/
>     ```
> 
>     This was written with commands in case folks want to copy paste it, to reduce error.
>
> 5. Launch the console with 
>    
>    ```bash
>    jekyll console
>    ```
> 
{: .hands_on}

## Which PRs are Single Cell RPs?

```ruby
require 'yaml'
data = YAML.load_file('metadata/github.yml')
```

Let's fetch data from within the GTN's infrastructure. You can see documentation for some of our APIs in the [RDoc](https://training.galaxyproject.org/training-material/gtn_rdoc/). E.g. here is how we document [`TopicFilter.list_materials_structured`](https://training.galaxyproject.org/training-material/gtn_rdoc/TopicFilter.html#method-c-list_materials_structured).

```ruby
# Obtain all single cell materials
mats = TopicFilter
 .list_materials_structured(site, 'single-cell')
 .map { |k, v| v['materials'] }
 .flatten
 .uniq { |x| x['id'] }

# And flatten them into a useful list of old folder names
@sc = mats.map{|x| x['ref_tutorials'].map{|t| t['redirect_from']} + x['ref_slides'].map{|t| t['redirect_from']}}
    .flatten.uniq
    .reject{|x| x=~ /\/short\//} # /short/ is a folder of redirects
    .map{|x| x.split('/')[1..-2].join('/')} # Remove the filename.

```

Let's go ahead and patch array to let us calculate a mean, because laziness is great actually

```ruby
class Array
  def mean
    return self.sum / (0.0 + self.length)
  end
end
```

Here's how we'll define what is or isn't a single cell tutorial, based on URL:

```ruby
def is_sc(path)
  p = path.gsub(/-ES/, '_ES').gsub(/-CAT/, '_CAT'); 

  if p =~ /\/single-cell\// then
    return 1
  end

  # If anything in @sc is a prefix for p, then it's a single-cell file.
  if @sc.any?{|sc| p.start_with?(sc)} then
    return 1
  end
  
  return 0
end
```

Now let's obtain everything that IS a single cell PR. Here we define everything with 50% or more of the files being "single cell" files, as a single cell PR. How did we arrive at 50%? We made some plots and spot-checked individual results to see what made sense.

```ruby
sc_prs = data
  .reject{|num, pr| pr['author']['login'] == 'github-actions'} # Remove all automation
  .map{|num, pr| 
    [
      num,
      pr['files']
        .reject{|f| f['path'].split('/')[2] == 'images' && f['path'] !~ /scrna/ }
        .reject{|f| f['path'] =~ /^assets/}
        .map{|f| is_sc(f['path'])}
      .mean
    ]
  }
  .reject{|num, sc| sc < 0.5}
  # Reject NaN
  .select{|num, sc| sc == sc}
```

Output number 1 done!

```ruby
File.open('dist.txt', 'w') do |f|
  f.puts "num\tdist"
  sc_prs.each do |num, sc|
    f.puts "#{num}\t#{sc}"
  end
end
```

## Classifying PR Content

Let's write a classifier for each file type, to enhance our statistics:

```ruby
def classify(path)
  if path =~ /tutorial[A-Z_]*?\.md/ then
    return 'tutorial'
  elsif path =~ /slides.*\.html/ then
    return 'slides'
  elsif path =~ /faqs.*md/ || path =~ /^snippets/ then
    return 'faq'
  elsif path =~ /metadata.yaml/ || path == 'CONTRIBUTORS.yaml' || path =~ /index.md$/ || path =~ /README.md$/ then
    return 'metadata'
  elsif path =~ /\/workflows\// then
    return 'workflows'
  elsif path =~ /data-(library|manager)/ then
    return 'data-library'
  elsif path =~ /.bib$/ then
    return 'bibliography'
  elsif path =~ /\/images\// then
    return 'image'
  elsif path =~ /tutorials\/.*md/ then
    return 'tutorial'
  elsif path =~ /_plugins/ || path =~ /^bin/ || path =~ /_layouts/ || path =~ /_include/ || path == '_config.yml' || path =~ /assets/ || path =~ /Gemfile/ || path =~/shared/ then
    return 'framework'
  elsif path =~ /metadata\/.*.yaml/ then
    return 'metadata'
  elsif path =~ /metadata\/.*.csv/ || path =~ /Dockerfile/ then
    return 'ignore'
  elsif path =~ /^news/ then
    return 'news'
  end

  # This will raise an exception which will ensure we catch the case where we
  # haven't defined a classification rule for a file yet.
  1/ 0
end
```

And now we'll classify each of the single cell pull requests by their file type:

```ruby
results = []
results << [
  "num", "path", "class", "additions", "deletions", "createdAt", "mergedAt"
]
sc_prs.each do |num, _|
  data[num]['files'].reject{|f| f['path'] =~ /test-data/}.each do |f|
    results << [
      num, f['path'], 
      classify(f['path']),
      f['additions'], f['deletions'],
      data[num]['createdAt'], data[num]['mergedAt']
    ]
  end
end
```

Output number 2!

```ruby
# save to file.csv
File.open('sc.tsv', 'w') do |f|
  results.each do |r|
    f.puts r.join("\t")
  end
end
```

Preview of that data:

num | path | class | additions | deletions | createdAt | mergedAt
--- | ---- | ----- | --------- | --------- | --------- | --------
5484 | topics/single-cell/faqs/single_cell_omics.md | faq | 3 | 3 | 2024-10-29T12:24:51Z | 2024-10-29T12:47:23Z
5473 | topics/single-cell/tutorials/alevin-commandline/tutorial.md | tutorial | 48 | 47 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case-jupyter_basic-pipeline/tutorial.md | tutorial | 3 | 2 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case_FilterPlotandExploreRStudio/tutorial.md | tutorial | 1 | 0 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case_FilterPlotandExplore_SeuratTools/tutorial.md | tutorial | 1 | 0 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case_JUPYTER-trajectories/tutorial.md | tutorial | 6 | 5 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case_alevin-combine-datasets/tutorial.md | tutorial | 4 | 4 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case_alevin/tutorial.md | tutorial | 3 | 2 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md | tutorial | 1 | 0 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case_monocle3-rstudio/tutorial.md | tutorial | 1 | 0 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case_monocle3-trajectories/tutorial.md | tutorial | 1 | 0 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5473 | topics/single-cell/tutorials/scrna-case_trajectories/tutorial.md | tutorial | 1 | 0 | 2024-10-25T11:24:40Z | 2024-10-28T12:59:16Z
5447 | topics/single-cell/tutorials/GO-enrichment/tutorial.md | tutorial | 7 | 3 | 2024-10-11T11:22:34Z | 2024-10-11T13:36:22Z
5445 | topics/single-cell/tutorials/scrna-case-cell-annotation/slides.html | slides | 1 | 1 | 2024-10-11T10:52:08Z | 2024-10-12T15:44:00Z
5443 | topics/single-cell/faqs/single_cell_omics.md | faq | 4 | 2 | 2024-10-11T10:16:18Z | 2024-10-12T15:43:24Z
5416 | topics/single-cell/metadata.yaml | metadata | 4 | 2 | 2024-10-07T22:26:36Z | 2024-10-12T15:35:11Z
5416 | topics/single-cell/tutorials/GO-enrichment/tutorial.md | tutorial | 1 | 0 | 2024-10-07T22:26:36Z | 2024-10-12T15:35:11Z

## Contributors over time

Additionally we want to figure out how our contributions and contributors changed over time.

```ruby
require 'yaml'
require 'date'
require 'pp'
```

Here are the current classifications of contributors:

```ruby
KEYS = %w[authorship editing testing ux infrastructure translation data]
```

Let's get all the single cell tutorials:

```ruby
tutorials = Dir.glob("topics/single-cell/tutorials/*/tutorial.md")
```

We'll want to get all the timepoints from 2019 (when single cell was added) to
2025, so we'll setup an empty data structure for this.

**NB** This is NOT the best solution, this is just a simple brute force solution because the runtime is fine actually.

```ruby
timepoints = []
(2019..2025).each do |year|
  (1..12).each do |month|
    s = "#{year}-#{month}-01T00:00:00Z"
    timepoints << [
      s,
      DateTime.parse(s).to_time.to_i
    ]
  end
end
contribs_over_time = timepoints.map{|n, t| [t, KEYS.map{|k| [k, []]}.to_h]}.to_h
```

Ok, let's get the history of each tutorial

```ruby
tutorials.each do |tutorial|
  # if tutorial !~ /bulk/
  #   next
  # end

  git=`git log --follow --name-only --format="GTN_GTN %H %at" #{tutorial}`
  commits = git.split("GTN_GTN ")
  commits.reject!{|c| c.empty?}
  commits.map!{|c| 
    res = c.gsub(/\n+/, "\t").split(/\t/)
    if res.size > 2
      puts "ERROR: #{res}"
    end

    hash = res[0].split(' ')[0]
    time = res[0].split(' ')[1].to_i

    f = res[1]
    contents_at_time = `git show #{hash}:#{f}`
    begin
      contents_meta = YAML.load(contents_at_time)
    rescue
      next
    end

    if contents_meta.nil?
      next
    end

    if contents_meta.key?("contributors")
      c = {
        'authorship' => contents_meta["contributors"],
      }
    else
      c = contents_meta["contributions"]
    end

    squashed_i = DateTime.parse(Time.at(time).strftime("%Y-%m-01T00:00:00Z")).to_time.to_i

    {
      :hash => hash,
      :time => time,
      :date => Time.at(time),
      :sqsh => squashed_i, # The time rounded to the month
      :path => res[1],
      :role => c
    }
  }

  # For every commit
  commits.reverse.compact.each do |c|
    KEYS.each do |k|
      # For every role
      if c[:role].key?(k)
        # add to contribs now and at every time point in the future
        now_and_future_keys = contribs_over_time.keys.select{|t| t >= c[:sqsh] }
        now_and_future_keys.each do |t|
          contribs_over_time[t][k] << c[:role][k]
          contribs_over_time[t][k].flatten!
          contribs_over_time[t][k].uniq!
        end
      end
    end
  end
end
```

See, terrible, but it works!

```ruby
contribs_over_time.reject!{|k, v| v.values.all?{|vv| vv.empty?}}
pp contribs_over_time

File.open("sc-roles.tsv", "w") do |f|
  f.write("date\tarea\tcount\tcontributors\n")

  KEYS.each do |k|
    f.write(contribs_over_time.map{|date, roles| [date, roles[k]]}.map{|date, contribs| "#{Time.at(date).strftime("%Y-%m-01")}\t#{k}\t#{contribs.count}\t#{contribs.join(',')}"}.join("\n"))
  end
end
```

With that, I think we're done! Ready to plot.
