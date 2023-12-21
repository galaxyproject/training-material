# frozen_string_literal: true

BASE_REF = ARGV[0]
redirects = `git diff --name-status -C #{BASE_REF}`.split(/\n/).select { |x| x[0] == 'R' }

redirects.each do |x|
  _, redir_from, redir_to = x.split(/\t/)
  redir_from = redir_from.gsub(/.md/, '').gsub(/.html/, '')

  puts "Adding redirect to #{redir_to}"
  f = File.open(redir_to, 'r')
  contents = f.read.split(/\n/)
  contents = contents[0..0] + ['redirect_from:', "- /#{redir_from}"] + contents[1..]
  f.close

  File.write(redir_to, contents.join("\n"))
end
