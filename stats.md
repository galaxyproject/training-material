---
layout: base
---
{% include _includes/default-header.html %}


<!-- tutorial stats -->
{% assign tutorials = site.pages | where:"layout", "tutorial_hands_on" %}

<!-- topic stats -->
{% assign topics = site.data | where_exp: "item", "item.type" %}
{% assign topics_science = topics | where: "type","use" | sort: "name" %}
{% assign topics_technical = topics | where_exp: "item", "item.type != 'use'"%}

<!-- contributors stats -->
{% assign contributors = site.data['contributors'] | where_exp: "item", "item.halloffame != 'no'" | sort: "joined" %}
{% assign contributors_by_month = contributors | group_by: "joined" %}
{% assign contributors_over_time = "" | split: ',' %}
{% assign contributors_over_time_labels = "" | split: ',' %}
{% assign contributors_count = 0 %}
{% for m in contributors_by_month %}
      {% assign new_contributors = m.items | size %}
      {% assign contributors_count = contributors_count | plus: new_contributors%}
      {% assign contributors_over_time = contributors_over_time | push: contributors_count %}
      {% assign contributors_over_time_labels = contributors_over_time_labels | push: m.name %}
{% endfor %}


<!-- use chart.js for graphs -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.9.3/Chart.bundle.js"></script>
<!-- plugin for adding data labels to charts -->
<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels@0.7.0"></script>


<style type="text/css" media="all">
.card-title {
	font-size: 2rem;
}
</style>


<div class="container main-content">
<section>
<h1>GTN Statistics</h1>

<div class="stats">
<div class="row">

 <!-- stats cards -->

 <!-- number of topics -->
 <div class="card" style="width: 30%;margin:1em">
  <div class="card-body">
   <h5 class="card-title">{{ topics | size }} Topics</h5>
   </div>
 </div>

 <!-- number of tutorials -->
 <div class="card" style="width: 30%;margin:1em">
  <div class="card-body">
   <h5 class="card-title">{{ tutorials | size }} Tutorials</h5>
   <!--<p class="card-text">amazing!.</p>-->
  </div>
 </div>

<!-- number of contributors -->
 <div class="card" style="width: 30%;margin:1em">
  <div class="card-body">
   <h5 class="card-title">{{ contributors | size }} Contributors</h5>
  </div>
 </div>

<!-- tutorials per topic -->
<div class="card" style="width: 47%;margin:1em">
  <div class="card-body">
   <h5 class="card-title">{{ topics_science | size }} Scientific Topics</h5>
   <canvas id="tutorialsBar" width="400" height="400"></canvas>
   </div>
 </div>

 <!-- conttributors over time  -->
 <div class="card" style="width: 47%;margin:1em">
  <div class="card-body">
   <h5 class="card-title">Growing Community</h5>
   <canvas id="contributorsGraph" width="400" height="400"></canvas>
   </div>
 </div>

 <!-- tutorials per topic (technical topics) -->
<div class="card" style="width: 75%;margin:1em">
  <div class="card-body">
   <h5 class="card-title">{{ topics_technical | size }} Technical Topics</h5>
   <canvas id="tutorialsBarTechnical" width="400"></canvas>
   </div>
 </div>
 <!-- end stats cards -->

</div>
</div>

</section>
</div>


<!-- make the charts -->
<script type="text/javascript">
Chart.plugins.unregister(ChartDataLabels);

function genColors(size) {
	var o = [];
	for(i = 0; i < size; i++){
		o.push(`hsl(${ i * 360 / size }, 100%, 50%)`)
	}
	return o;
}


// Charts displaying number of tutorials per topic
// Scientific Topics
var tutoBar = document.getElementById('tutorialsBar');

var data_tutos = [{% for topic in topics_science %}{{site.pages | topic_filter: topic.name | size }}{%unless forloop.last%},{%endunless%}{% endfor %}];
var labels_topics = [{% for topic in topics_science %}"{{ topic.title }}"{%unless forloop.last%},{%endunless%}{% endfor %}];

var tutorialsBar = new Chart(tutoBar, {
  type: 'horizontalBar',

  data: {
    labels: labels_topics,
    datasets: [{
      backgroundColor: genColors(data_tutos.length),
      data: data_tutos
    }]
  },
  plugins: [ChartDataLabels],

  options: {
    legend: {
	  display: false
	},
    title: {
      display: true,
      text: 'Tutorials per Topic'
    },
    plugins: {
      datalabels: {
        anchor: 'end',
        align: 'end'
      }
    }
  }
});

// Technical Topics
// Chart displaying number of tutorials per topic
var tutoBarTechnical = document.getElementById('tutorialsBarTechnical');

var data_tutos = [{% for topic in topics_technical %}{{site.pages | topic_filter: topic.name | size }}{%unless forloop.last%},{%endunless%}{% endfor %}];
var labels_topics = [{% for topic in topics_technical %}"{{ topic.title }}"{%unless forloop.last%},{%endunless%}{% endfor %}];

var tutorialsBar = new Chart(tutoBarTechnical, {
  type: 'horizontalBar',

  data: {
    labels: labels_topics,
    datasets: [{
      backgroundColor: genColors(data_tutos.length),
      data: data_tutos
    }]
  },
  plugins: [ChartDataLabels],

  options: {
    legend: {
	  display: false
	},
    title: {
      display: true,
      text: 'Tutorials per Topic'
    },
    plugins: {
      datalabels: {
        anchor: 'end',
        align: 'end'
      }
    }
  }
});



// Contributors chart
var contributorsGraph = document.getElementById('contributorsGraph');

var data_contributors = [{%for c in contributors_over_time %}{x:"{{contributors_over_time_labels[forloop.index]}}" , y: {{c}} } {%unless forloop.last%},{%endunless%}{%endfor%}];

var labels_contributors = [{%for l in contributors_over_time_labels %}"{{l}}"{%unless forloop.last%},{%endunless%}{%endfor%}];

var tutorialsBar = new Chart(contributorsGraph, {
  type: 'line',
  data: {
    datasets: [{
      data: data_contributors,
    }]
  },

  options: {
    scales: {
      yAxes: [{
        ticks: {
          beginAtZero: true
        }
      }],
      xAxes: [{
        type: 'time',
        time: {
          displayFormats:{month:'YYYY-MM'},
          min:'2017-10',
          unit: 'month',
          distribution: 'linear'
        }
      }]
    },
    legend: {
	  display: false
	},
    title: {
      display: true,
      text: 'Contributors over time'
    }
  }
});

</script>

{% include _includes/default-footer.html %}
