require './_plugins/gtn/mod'
require 'time'

class Hash
  def to_prometheus
    if self.keys.length > 0
      '{' + self.map{|k, v| "#{k}=\"#{v}\""}.join(',') + '}'
    else
      ''
    end
  end
end


module Gtn
  module Metrics
    def self.iqr(array)
      # https://stackoverflow.com/questions/8856716/calculate-interquartile-mean-from-ruby-array/8856863#8856863
      arr = array.sort
      length = arr.size
      quart = (length/4.0).floor
      fraction = 1-((length/4.0)-quart)
      new_arr = arr[quart..-(quart + 1)]
      (fraction*(new_arr[0]+new_arr[-1]) + new_arr[1..-2].inject(:+))/(length/2.0)
    end

    def self.bin_width(values)
      # https://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule
      2 * self.iqr(values) / (values.length ** (1/3))
    end

    def self.histogram(values)
      values.map!{|v| v.to_i}
      width = self.bin_width(values)
      bins = ((values.max - values.min) / width).ceil

      (0..bins).map{|bin_idx|
        left = values.min + (bin_idx * width)
        right = values.min + ((bin_idx + 1) * width)
        count = values.select{|x| left <= x and x < right }.length
        {
          :value => count,
          'le' => bin_idx == bins ? '+Inf' : right.to_s
        }
      }
    end

    def self.histogram_dates(values)
      day_bins = [1, 7, 28, 90, 365, 365 * 2, 365 * 3, 365 * 5]
      last_bin = 0
      day_bins.map{|bin, idx|
        count = values.select{|x| last_bin <= x and x < bin }.length
        last_bin = bin
        {
          :value => count,
          'le' => bin == day_bins[-1] ? '+Inf' : bin
        }
      }
    end

    def self.segment(values, attr)
      [
        {
          :value => values.select{|v| v.has_key? attr }.length,
          "#{attr}" => true,
        },
        {
          :value => values.select{|v| ! v.has_key? attr }.length,
          "#{attr}" => false,
        }
      ]
    end

    def self.segment_page_by_key(values, key)
      possible_keys = values.map{|v| v.data[key].to_s}.sort.uniq
      possible_keys.map{|k|
        {
          :value => values.select{|v| v.data[key] == k }.length,
          "#{key}" => k,
        }
      }
    end

    def self.collect_metrics(site)
      tutorials = site.pages.select{|x| x.data['layout'] == 'tutorial_hands_on'}
      first_commit = Date.parse('2015-06-29')
      today = Date.today()

      output = {
        'gtn_pages_total' => {
          :value => self.segment_page_by_key(site.pages, 'layout'),
          :help => 'Total number of Pages',
          :type => 'counter'
        }, 'gtn_contributors_total' => {
          :value => self.segment(site.data['contributors'].values.select{|x| x['halloffame'] != 'no'}, 'orcid'),
          :help => 'Total number of Contributors',
          :type => 'counter'
        },
        'gtn_tutorials_total' => {
          :value => tutorials.length,
          :help => 'Total number of Hands-on Tutorials',
          :type => 'counter'
        },
        'gtn_faqs_total' => {
          :value => site.pages.select{|x| x.data['layout'] == 'faq'}.length,
          :help => 'Total number of FAQs',
          :type => 'counter'
        },
        'gtn_project_years_total' => {
          :value => (today - first_commit).to_f / 365,
          :help => 'Total years of project lifetime',
          :type => 'counter'
        },
        'tutorial_update_age_days' => {
          :type => 'histogram',
          :help => 'How recently was every single Hands-on Tutorial touched within the GTN, grouped by days since last edited.',
          :value => self.histogram_dates(
            tutorials
            .map{|page|
              Time.now() - Gtn::ModificationTimes.obtain_time(page['path'].gsub(/^\//, ''))
            }
            .map{|seconds| seconds / 3600.0 / 24.0}
          )
        },
        'contributor_join_age_days' => {
          :type => 'histogram',
          :help => 'When did contributors join? Buckets of their ages by days since joined.',
          :value => self.histogram_dates(
            site.data['contributors']
            .select{|x| x['halloffame'] != 'no'}
            .map{|_, contributor|
              (Date.today() - Date.parse("#{contributor['joined']}-01")).to_i
            }
          )
        },
      }
      output
    end

    def self.generate_metrics(site)
      data =  self.collect_metrics(site)
      output = data.map{|k, v|
        out = "# HELP #{k} #{v[:help]}\n# TYPE #{k} #{v[:type]}\n"

        if v[:value].is_a?(Array)
          v[:value].each{|val|
            attrs = val.select{|k, v| k != :value}.to_h
            out += "#{k}#{attrs.to_prometheus} #{val[:value]}\n"
          }
        else
          attrs = v.select{|k, v| k != :value and k != :help and k != :type}.to_h
          out += "#{k}#{attrs.to_prometheus} #{v[:value]}\n"
        end

        out
      }.join("\n")
    end
  end
end
