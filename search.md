---
layout: base
---

{% include _includes/default-header.html %}

<h2> Search Tutorials </h2>

<!-- Html Elements for Search -->

<div id="tutorial_list">

<div id="search-container">
<input type="text" id="search-input" placeholder="search...">
<table class="table table-responsive table-striped" >
  <thead>
    <tr>
      <th>Topic</th>
      <th>Lesson</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody class="list" id="results-container">
</tbody>
</table>
</div>

</div>

<!-- Script pointing to search-script.js -->
<script src="assets/js/search-script.js" type="text/javascript"></script>

<!-- Configuration -->
<script>
SimpleJekyllSearch({
  searchInput: document.getElementById('search-input'),
  resultsContainer: document.getElementById('results-container'),
  json: 'search.json',
  noResultsText: ("No result found!"),
  searchResultTemplate: '{entry}'

})

</script>

{% include _includes/default-footer.html %}
