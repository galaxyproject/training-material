require 'yaml'

module Jekyll
  module Tags
    class SnippetIncludeTag < IncludeTag

      def markdownify(text)
        @site.find_converter_instance(
          Jekyll::Converters::Markdown
        ).convert(text.to_s)
      end

      def get_icon(icon)
         if icon.start_with?("fa")
          %Q(<i class="#{icon}" aria-hidden="true"></i><span class="visually-hidden">#{@text}</span>)
         elsif icon.start_with?("ai")
          %Q(<i class="ai #{icon}" aria-hidden="true"></i><span class="visually-hidden">#{@text}</span>)
         end
      end

       def get_config(context)
           context.registers[:site].config['icon-tag']
       end

      def render(context)

        @site ||= context.registers[:site]

        file = render_variable(context) || @file
        validate_file_name(file)

        @site.inclusions[file] ||= locate_include_file(file)
        inclusion = @site.inclusions[file]

        add_include_to_dependency(inclusion, context) if @site.config["incremental"]

        context.stack do
          context["include"] = parse_params(context) if @params
          x = "#{inclusion.render(context)}"
          p = context["include"]

          box_start=""
          box_end=""
          if x.slice(0, 3) == '---'
            metadata = YAML.load(x)

            # allow overriding box type with include parameter ("none" to render without a box)
            if not p.nil? and p["box_type"]
                box_type = p["box_type"]
            else
                box_type = metadata['box_type']
            end
            icons = get_config(context)

            if box_type == 'tip'
                icon_text = icons['tip']
                box_start = '> ### '+get_icon(icon_text)+' Tip: ' + metadata['title']
                box_end   = "\n{: .tip}"
            end
            if box_type == 'hands_on'
                icon_text = icons['hands_on']
                box_start = '> ### '+get_icon(icon_text)+' Hands-on: ' + metadata['title']
                box_end   = "\n{: .hands_on}"
            end
            if box_type == 'comment'
                icon_text = icons['comment']
                box_start = '> ### '+get_icon(icon_text)+' ' + metadata['title']
                box_end   = "\n{: .comment}"
            end
          end
          y = x.gsub(/\A---(.|\n)*?---/, '')
          if box_start != ""
             y = y.gsub(/\R+/,"\n> ")
             #puts box_start+y+box_end
          end

          '<!--SNIPPET-->' + markdownify(box_start+y+box_end).gsub(/\R+/, '').gsub('<h3','<h3 data-toc-skip')
        end
      end

      private

      def locate_include_file(file)
        @site.includes_load_paths.each do |dir|
          path = PathManager.join(dir, file)
          return Inclusion.new(@site, dir, file) if valid_include_file?(path, dir)
        end
        raise IOError, could_not_locate_message(file, @site.includes_load_paths, @site.safe)
      end

      def valid_include_file?(path, dir)
        File.file?(path) && !outside_scope?(path, dir)
      end

      def outside_scope?(path, dir)
        @site.safe && !realpath_prefixed_with?(path, dir)
      end

      def realpath_prefixed_with?(path, dir)
        File.realpath(path).start_with?(dir)
      rescue StandardError
        false
      end

      def add_include_to_dependency(inclusion, context)
        return unless context.registers[:page]&.key?("path")

        @site.regenerator.add_dependency(
          @site.in_source_dir(context.registers[:page]["path"]),
          inclusion.path
        )
      end
    end
  end

  module RegexReplace
    def regex_replace(str, regex_search, value_replace)
      regex = /#{regex_search}/m
      return str.gsub(regex, value_replace)
    end

    def regex_replace_once(str, regex_search, value_replace)
      regex = /#{regex_search}/m
      return str.sub(regex, value_replace)
    end
  end

end

Liquid::Template.register_tag("snippet", Jekyll::Tags::SnippetIncludeTag)
Liquid::Template.register_filter(Jekyll::RegexReplace)
