
Jekyll::Hooks.register :pages, :post_init do |page|
  page.data['js_requirements'] = {
    'mathjax' => page.content =~ /\$\$/,
    'mermaid' => page.content =~ /```mermaid/ || page.content =~ /pre class="mermaid"/ || page.data['layout'] == 'workflow-list',
  }

  # simplify math expressions; automatically replace
  # $$ with \\( and \\) to support inline math with $$
  if page.content
    page.content = page.content.gsub(/\$\$(.*?)\$\$/) {|m| m.gsub('_','\\\\_')}
    page.content = page.content.gsub(/\$\$(.*?)\$\$/,'\\\\\\\\( \1 \\\\\\\\)')
  end

end
