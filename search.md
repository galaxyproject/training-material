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

var tutorials = {% dump_search_view testing %};

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

        $("#results-container").html(results_final.map(x => `
        <div class='col-sm-6'>
        <div class='card'>
        <div class='card-body'>
          <h5 class='card-title'>${x.title}</h5>
          <h6 class='card-subtitle text-muted'>${x.topic}</h6>
          <p>${x.tags.join(' ')}</p>
          <p>${x.contributors}</p>
          <a class='btn btn-primary' href='${x.url}'>View Tutorial</a>
          </div>
          </div>
          </div>
                    `));
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
