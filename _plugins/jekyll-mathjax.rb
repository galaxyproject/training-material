
Jekyll::Hooks.register :pages, :post_init do |page|
  page.data['js_requirements'] = {
    'mathjax' => page.content =~ /\$\$/ || page.content =~ /\\\(/,
    'mermaid' => page.content =~ /```mermaid/ || page.content =~ /pre class="mermaid"/ || page.data['layout'] == 'workflow-list',
  }

  # some fixes for mathjax: escape underscores
  # both in inline mode \\( .. \\) and block mode $$ ..$$
  if page.content
    page.content = page.content.gsub(/\$\$(.*?)\$\$/) { |m| m.gsub('_', '\\\\_') }
    page.content = page.content.gsub(/\\\\\((.*?)\\\\\)/) { |m| m.gsub('_', '\\\\_') }

  end
end
