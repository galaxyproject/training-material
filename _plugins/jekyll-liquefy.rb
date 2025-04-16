module Jekyll
    module LiquefyFilter
        def liquefy(input)
            Liquid::Template.parse(input).render(@context)
        end
    end
end

Liquid::Template.register_filter(Jekyll::LiquefyFilter)
