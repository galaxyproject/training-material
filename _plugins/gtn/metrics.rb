# frozen_string_literal: true

require './_plugins/gtn/mod'
require 'time'

# Monkey patch Hash to make it Prometheus compatible
class Hash
  def to_prometheus
    if keys.length.positive?
      inner = map { |k, v| "#{k}=\"#{v}\"" }.join(',')
      "{#{inner}}"
    else
      ''
    end
  end
end

module Gtn
  # Module for generating metrics for Prometheus
  module Metrics
    def self.iqr(array)
      # https://stackoverflow.com/questions/8856716/calculate-interquartile-mean-from-ruby-array/8856863#8856863
      arr = array.sort
      length = arr.size
      quart = (length / 4.0).floor
      fraction = 1 - ((length / 4.0) - quart)
      new_arr = arr[quart..-(quart + 1)]
      ((fraction * (new_arr[0] + new_arr[-1])) + new_arr[1..-2].inject(:+)) / (length / 2.0)
    end

    def self.bin_width(values)
      # https://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule
      2 * iqr(values) / (values.length**(1 / 3))
    end

    def self.histogram(values)
      values.map!(&:to_i)
      width = bin_width(values)
      bins = ((values.max - values.min) / width).ceil

      (0..bins).map do |bin_idx|
        left = values.min + (bin_idx * width)
        right = values.min + ((bin_idx + 1) * width)
        count = values.select { |x| left <= x and x < right }.length
        {
          :value => count,
          'le' => bin_idx == bins ? '+Inf' : right.to_s
        }
      end
    end

    def self.histogram_dates(values)
      day_bins = [1, 7, 28, 90, 365, 365 * 2, 365 * 3, 365 * 5]
      last_bin = 0
      day_bins.map do |bin, _idx|
        count = values.select { |x| last_bin <= x and x < bin }.length
        last_bin = bin
        {
          :value => count,
          'le' => bin == day_bins[-1] ? '+Inf' : bin
        }
      end
    end

    def self.segment(values, attr)
      [
        {
          :value => values.select { |v| v.key? attr }.length,
          attr.to_s => true,
        },
        {
          :value => values.reject { |v| v.key? attr }.length,
          attr.to_s => false,
        }
      ]
    end

    def self.segment_page_by_key(values, key)
      possible_keys = values.map { |v| v.data[key].to_s }.sort.uniq
      possible_keys.map do |k|
        {
          :value => values.select { |v| v.data[key] == k }.length,
          key.to_s => k,
        }
      end
    end

    def self.collect_metrics(site)
      tutorials = site.pages.select { |x| x.data['layout'] == 'tutorial_hands_on' }
      first_commit = Date.parse('2015-06-29')
      today = Date.today

      {
        'gtn_pages_total' => {
          value: segment_page_by_key(site.pages, 'layout'),
          help: 'Total number of Pages',
          type: 'counter'
        },
        'gtn_contributors_total' => {
          value: segment(site.data['contributors'].values.reject { |x| x['halloffame'] == 'no' }, 'orcid'),
          help: 'Total number of contributors',
          type: 'counter'
        },
        'gtn_organisations_total' => {
          value: segment(site.data['organisations'].values.reject { |x| x['halloffame'] == 'no' }, 'orcid'),
          help: 'Total number of organisations',
          type: 'counter'
        },
        'gtn_funders_total' => {
          value: segment(site.data['funders'].values.reject { |x| x['halloffame'] == 'no' }, 'orcid'),
          help: 'Total number of funders',
          type: 'counter'
        },
        'gtn_tutorials_total' => {
          value: tutorials.length,
          help: 'Total number of Hands-on Tutorials',
          type: 'counter'
        },
        'gtn_faqs_total' => {
          value: site.pages.select { |x| x.data['layout'] == 'faq' }.length,
          help: 'Total number of FAQs',
          type: 'counter'
        },
        'gtn_project_years_total' => {
          value: (today - first_commit).to_f / 365,
          help: 'Total years of project lifetime',
          type: 'counter'
        },
        'tutorial_update_age_days' => {
          type: 'histogram',
          help: 'How recently was every single Hands-on Tutorial touched within the GTN, ' \
                'grouped by days since last edited.',
          value: histogram_dates(
            tutorials
            .map do |page|
              Time.now - Gtn::ModificationTimes.obtain_time(page['path'].gsub(%r{^/}, ''))
            end
            .map { |seconds| seconds / 3600.0 / 24.0 }
          )
        },
        'contributor_join_age_days' => {
          type: 'histogram',
          help: 'When did contributors join? Buckets of their ages by days since joined.',
          value: histogram_dates(
            site.data['contributors']
            .reject { |x| x['halloffame'] == 'no' }
            .map do |_, contributor|
              (Date.today - Date.parse("#{contributor['joined']}-01")).to_i
            end
          )
        },
      }
    end

    def self.generate_metrics(site)
      data = collect_metrics(site)
      data.map do |k, v|
        out = "# HELP #{k} #{v[:help]}\n# TYPE #{k} #{v[:type]}\n"

        if v[:value].is_a?(Array)
          v[:value].each do |val|
            attrs = val.except(:value).to_h
            out += "#{k}#{attrs.to_prometheus} #{val[:value]}\n"
          end
        else
          attrs = v.select { |k2, _v2| k2 != :value and k2 != :help and k2 != :type }.to_h
          out += "#{k}#{attrs.to_prometheus} #{v[:value]}\n"
        end

        out
      end.join("\n")
    end
  end
end
