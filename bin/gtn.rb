# frozen_string_literal: true

CONTRIBUTORS = YAML.load_file('CONTRIBUTORS.yaml')
ORGANISATIONS = YAML.load_file('ORGANISATIONS.yaml')
FUNDERS = YAML.load_file('FUNDERS.yaml')

def automagic_loading(f)
  # Remove our documentation
  f.reject! { |k, v| k == 'description' and v.is_a?(String) }
  f.reject! { |k| k == 'examples' }

  # Auto-replace CONTRIBUTORS in enums.
  f.each do |k, v|
    if v.is_a?(Hash)
      automagic_loading(v)
    elsif v.is_a?(Array)
      v.replace CONTRIBUTORS.keys if (k == 'enum') && (v[0] == 'CONTRIBUTORS')
      v.replace FUNDERS.keys if (k == 'enum') && (v[0] == 'FUNDERS')
      v.replace ORGANISATIONS.keys if (k == 'enum') && (v[0] == 'ORGANISATIONS')
      v.flatten.each { |x| automagic_loading(x) if x.is_a?(Hash) }
    end
  end
  f
end
