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
      if k == 'enum'
        repl = []
        # If one of the elements in this array is CONTRIBUTORS, replace it with the same named variable
        if v.find { |x| x == 'CONTRIBUTORS' }
          repl << CONTRIBUTORS.keys
        end
        if v.find { |x| x == 'FUNDERS' }
          repl << FUNDERS.keys
        end
        if v.find { |x| x == 'ORGANISATIONS' }
          repl << ORGANISATIONS.keys
        end
        if repl.length.positive?
          v.replace repl.flatten
        end
      end
      v.flatten.each { |x| automagic_loading(x) if x.is_a?(Hash) }
    end
  end
  f
end
