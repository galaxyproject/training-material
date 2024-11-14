#!/usr/bin/env ruby
require 'json'
require 'net/http'
require 'uri'
require 'yaml'

# Existing known PRs in our database
data = YAML.load_file('metadata/github.yml') || {}

def gh_cli_pr_list(count: 3)
  d = JSON.parse(`gh pr list --search 'sort:updated-desc is:merged' --limit #{count} --json number`)
  d.map { |pr| pr['number'] }
end

def gh_cli_pr_info(num, fields)
  JSON.parse(`gh pr view #{num} --json #{fields.join(',')}`)
end

# Fetch all PRs in the range in order to backfill anything missing. Shouldn't be needed anymore.
# todo = data.keys.min..data.keys.max
# todo = todo.reject { |pr| data.key?(pr) }

# The most recently touched 30 merged PRs
todo = gh_cli_pr_list(count: 30)

todo.each do |pull|
  puts "Fetching PR #{pull}"
  # We already have this PR.
  next if data.key?(pull)

  # Fetch data
  fields = %w[state title author createdAt updatedAt closedAt mergedAt mergedBy
              labels headRefName headRepository reactionGroups reviews files url]

  begin
    pr = gh_cli_pr_info(pull, fields)
  rescue
    puts "Failed to fetch PR #{pull}. Is it an issue?"
    next
  end

  pr['labels'] = pr['labels'].map { |l| l['name'] }
  pr['reviews'] = pr['reviews'].map { |r| r.slice('author', 'state', 'submittedAt', 'reactionGroups') }

  data[pull] = pr
end
File.write('metadata/github.yml', data.to_yaml)
