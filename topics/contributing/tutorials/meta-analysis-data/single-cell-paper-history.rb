require 'yaml'
require 'date'
require 'pp'

KEYS = %w[authorship editing testing ux infrastructure translation data]

tutorials = Dir.glob("topics/single-cell/tutorials/*/tutorial.md")

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

contribs_over_time.reject!{|k, v| v.values.all?{|vv| vv.empty?}}
pp contribs_over_time

File.open("tmp/data.tsv", "w") do |f|
  f.write("date\tarea\tcount\tcontributors\n")

  KEYS.each do |k|
    f.write(contribs_over_time.map{|date, roles| [date, roles[k]]}.map{|date, contribs| "#{Time.at(date).strftime("%Y-%m-01")}\t#{k}\t#{contribs.count}\t#{contribs.join(',')}"}.join("\n"))
  end
end
