#!/usr/bin/env ruby
require 'yaml'


# If there are unknown contributors, write a .mailmap file like:
# <gh-username> <user@earlham.ac.uk>
# <gh-username> <user@gmail.com>

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

file_contributors = `git log --use-mailmap --follow --pretty=%aE #{fn}`.lines.sort.uniq

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

missing = (known - current_contributors).sort.uniq
# These contributors not yet recognised
if missing.length > 0
  puts "Missing contributors: #{missing}"
end

# These emails might map to known users, but we don't know yet.
if unknown.length > 0
  puts "Unknown emails: #{unknown.sort.uniq}"
end
