---
layout: page
title: Search Tutorials
---

<script src="https://cdnjs.cloudflare.com/ajax/libs/lunr.js/2.3.9/lunr.min.js" integrity="sha512-4xUl/d6D6THrAnXAwGajXkoWaeMNwEKK4iNfq5DotEbLPAfk6FSxSP3ydNxqDgCw1c/0Z1Jg6L8h2j+++9BZmg==" crossorigin="anonymous"></script>

<div id="search-container">
	<input type="text" id="search-input" placeholder=" search..." class="nicer">{% icon search %}
	<div class="search-results row" id="results-container"></div>
</div>


<!-- Configuration -->
<script>
var tutorials = { {% for topic in site.data %}
    {% if topic[1].type == 'use' or topic[1].type == 'admin-dev' or topic[1].type == 'basics' %}
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
              <form method="GET" action="{{site.baseurl}}/search" style="display:inline"><input type="hidden" name="query" value="{{tag}}"><button class='label label-default tutorial_tag' id='{{ tag }}' style='{{ tag | colour_tag }}' title='Click to show all tutorials with this tag'>{{ tag  }}</button></form>
            {% endfor %}
            </p>
          {% endif %}
          <p>{% include _includes/contributor-badge-list.html contributors=tutorial.contributors %}</p>
          <a class='btn btn-primary' href='{{ site.baseurl }}{{ tutorial.url }}'>View Tutorial</a>
          </div>
          </div>
          </div>
          {% endcapture %}
      "{{ tutorial.url }}": {
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
    {% endif %}
  {% endfor %}
};



function search(idx, q){
	if(q.length > 2){
        var results_partial = idx.search(`*${q}*`),
            results_exact = idx.search(`${q}`),
            results_fuzzy = idx.search(`${q}~3`);

        thereMap  = Object.assign({}, ...results_partial.map((x) => ({[x.ref]: x.score})));

        results_exact.forEach(x => {
            if(thereMap[x.ref] !== undefined){
                if(thereMap[x.ref] < x.score + 4){
                    thereMap[x.ref] = x.score + 4
                }
            } else {
                    thereMap[x.ref] = x.score + 4
            }
        })
        results_fuzzy.forEach(x => {
            if(thereMap[x.ref] !== undefined){
                if(thereMap[x.ref] < x.score - 2){
                    thereMap[x.ref] = x.score - 2
                }
            } else {
                    thereMap[x.ref] = x.score - 2
            }
        })

        combined_results = Object.getOwnPropertyNames(thereMap);
        combined_results.sort((a, b) => {
            if (thereMap[a] > thereMap[b]) {
                return -1;
            }
            if (thereMap[a] < thereMap[b]) {
                return 1;
            }
            return 0;
        });

		var results_final = combined_results.map(x => {
			return tutorials['/' + x.replaceAll(".md", ".html")];
		}).filter(x => x !== undefined);

		$("#results-container").html(results_final.map(x => x.entry));
	}
}

fetch('{{ site.baseurl }}/search.json')
	.then(response => response.json())
	.then(data => {
		var idx = lunr.Index.load(data);

		var  params = (new URL(document.location)).searchParams;
		paramQuery = params.get('query');
		if(paramQuery){
			document.getElementById('search-input').value = paramQuery;
			search(idx, paramQuery);
		}

		$("#search-input").on("change keyup paste", function(){
			search(idx, $("#search-input").val());
		})
});
</script>
