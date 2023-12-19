def jupyter_pre_render(_site)
  nil
end

def jupyter_post_write(_site)
  nil
end

module Jekyll
  # Notebook Generation Disabled
  class RmarkdownGenerator < Generator
    def generate(_site)
      Jekyll.logger.info 'Notebooks disabled'
    end
  end

  # Notebook Generation Disabled
  class JupyterNotebookGenerator < Generator
    def generate(_site)
      Jekyll.logger.info 'Notebooks disabled'
    end
  end
end
