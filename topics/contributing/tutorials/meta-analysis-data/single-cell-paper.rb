require 'yaml'

data = YAML.load_file('metadata/github.yml')

@sc = File.read('sc-historical.txt').split("\n").map{|x| x[1..]}

class Array
  def mean
    return self.sum / (0.0 + self.length)
  end
end


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

File.open('dist.txt', 'w') do |f|
  f.puts "num\tdist"
  sc_prs.each do |num, sc|
    f.puts "#{num}\t#{sc}"
  end
end

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

  1/ 0
end

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

# save to file.csv
File.open('sc.tsv', 'w') do |f|
  results.each do |r|
    f.puts r.join("\t")
  end
end
