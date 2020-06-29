---
layout: base
---

{% include _includes/default-header.html %}

<div class="container main-content">
<section>
<h2> Search Tutorials </h2>

<!-- Html Elements for Search -->


<div id="search-container">
 <input type="text" id="search-input" placeholder="Search..." style="margin-right:1em;">{% icon search %}
 <div class="search-results row" id="results-container">

 </div>
</div>
</section>
</div>

<!-- Script pointing to search-script.js -->
<script src="assets/js/search-script.js" type="text/javascript"></script>

<!-- Configuration -->
<script>

var data= [ {% for topic in site.data %}
    {% unless topic[0] == 'contributors' %}
      {% assign topic_material = site.pages | topic_filter:topic[0] %}
      {% assign topic_title = topic[1].title %}
      {% for tutorial in topic_material %}

       {% capture result_entry %}
        <div class='col-sm-6'>
        <div class='card'>
        <div class='card-body'>
          <h5 class='card-title'>{{ tutorial.title | escape }}</h5>
          <h6 class='card-subtitle text-muted'>{{ topic_title}}</h6>
          <p class='card-text'> {{tutorial.description}}</p>
          {% if tutorial.tags %}
            <p>
            {% for tag in tutorial.tags %}
              <span class='label label-default tutorial_tag' id='{{ tag }}' style='{{ tag | colour_tag }}' title='Click to show only tutorials with this tag'>{{ tag  }}</span>
            {% endfor %}
            </p>
          {% endif %}
          <p>{% include _includes/contributor-badge-list.html contributors=tutorial.contributors %}</p>
          <a class='btn btn-primary' href='{{ site.baseurl }}{{ tutorial.url }}'>View Tutorial</a>
          </div>
          </div>
          </div>
          {% endcapture %}
      {
        "topic"    : "{{ topic_title }}",
        "title"    : "{{ tutorial.title | escape }}",
        "description": "{{ tutorial.description }}",
        "question" : "{{ tutoral.questions | join: ', '}}",
        "objectives"  : "{{ tutorial.objectives | join: ', ' }}",
        "tags"     : "{{ tutorial.tags | join: ', ' }}",
        "level"     : "{{ tutorial.level }}",
        "time_estimation": "{{ tutorial.time_estimation }}",
        "url"      : "{{ site.baseurl }}{{ tutorial.url }}",
        "level"     : "{{ tutorial.level}}",
        "contributors": "{{ tutorial.contributors | join: ', '}}",
        "entry"      : "{{ result_entry | strip_newlines | replace: '"',"'" }}"
      }{% unless forloop.last %},{% endunless %}
    {% endfor %}
    {% unless forloop.last %},{% endunless %}
    {% endunless %}
  {% endfor %}
]

var sjs = SimpleJekyllSearch({
  searchInput: document.getElementById('search-input'),
  resultsContainer: document.getElementById('results-container'),
  json: data,
  limit: '50',
  noResultsText: ("No result found!"),
  success: function(){},
  searchResultTemplate: '{entry}'
});

window.addEventListener('DOMContentLoaded', (event) => {
    console.log('DOM fully loaded and parsed');
    params = (new URL(document.location)).searchParams;
    paramQuery = params.get('query');
    if(paramQuery){
      document.getElementById('search-input').value = paramQuery;
      sjs.search(paramQuery);
    }
});


</script>

{% include _includes/default-footer.html %}
