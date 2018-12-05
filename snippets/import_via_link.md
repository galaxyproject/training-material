>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager ({% icon galaxy-upload %} on the top-right of the tool panel)
>    > {% if include.collection %}
>    > * Click on **Collection** on the top
>    > {% endif %}
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > {% if include.format %}
>    > * Change **Type** from "Auto-detect" to `{{ include.format }}`
>    > {% endif %}
>    > {% if include.genome %}
>    > * Change **Genome** to `{{ include.genome }}`
>    > {% endif %}
>    > * Press **Start**
>    > {% if include.collection %}
>    > * Click on **Build** when available
>    > * Enter a name for the collection
>    > * Click on **Create list** (and wait a bit)
>    > {% else %}
>    > * **Close** the window
>    > {% endif %}
>    > By default, Galaxy uses the URL as the name, so rename the files with a more useful name.
>    {: .tip}
>