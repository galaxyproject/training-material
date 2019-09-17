#!/usr/bin/env ruby
require 'yaml'
require 'pathname'

fn = ARGV[0]

# Any error messages
#errs = []
data = YAML.load_file(fn)
current_contributors = data['contributors']


# Full Contributors Data
CONTRIBUTORS = YAML.load_file('CONTRIBUTORS.yaml')
contributor_emails = CONTRIBUTORS.map{ |k, v|
  if v
    if v.has_key?("email")
      [v['email'], k]
    end
  end
}.compact.to_h

# Private map of emails to github usernames
if Pathname.new('private-contrib-map.yaml').exist?
  private_contrib_map = YAML.load_file('private-contrib-map.yaml')
  contributor_emails.merge!(private_contrib_map)
end

file_contributors = `git log --follow --pretty=%aE #{fn}`.lines.sort.uniq

fixed_contribs = file_contributors.map{ |email|
  email = email.strip
  if /users.noreply.github.com/.match(email)
    parts = /^([0-9]+\+)?(?<id>.*)@users.noreply.github.com/.match(email)
    # we just want their gh id
    parts[:id]
  else
    if contributor_emails.has_key?(email)
      contributor_emails[email]
    else
      email
    end
  end
}

# known contributors
known = fixed_contribs.select{ |x| ! /@/.match(x) }
unknown = fixed_contribs.select{ |x| /@/.match(x) }

# These contributors not yet recognised
puts "Missing contributors: #{(known - current_contributors).sort.uniq}"
if unknown.length > 0
  puts "Unknown emails: #{unknown.sort.uniq}"
end
