>
>    > ### {% icon tip %} Tip: Creating a new file
>    >
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the file contents into the text field
>    > {% if include.format %}
>    > * Change **Type** from "Auto-detect" to `{{ include.format }}`
>    > {% endif %}
>    > {% if include.genome %}
>    > * Change **Genome** to `{{ include.genome }}`
>    > {% endif %}
>    > {% if include.convertspaces %} * From the Settings menu ({% icon galaxy-gear %}) select **Convert spaces to tabs** {% endif %}
>    > * Press **Start** and **Close** the window
>    {: .tip}
>