---
layout: base
---

{% assign fdbk= site.data['feedback'] %}
{% assign start = fdbk | first  %}
{% assign topics = site.data | where_exp: "item", "item.title" %}
{% assign tutorials = site.pages | where:"layout", "tutorial_hands_on" %}

<!-- use chart.js for graphs -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.9.3/Chart.bundle.js"></script>
<!-- plugin for adding data labels to charts -->
<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels@0.7.0"></script>
<!-- plugin for palette -->
<script src="https://github.com/nagix/chartjs-plugin-colorschemes/releases/download/v0.4.0/chartjs-plugin-colorschemes.js"></script>
<!-- define the charts -->
<script type="text/javascript">
    Chart.plugins.unregister(ChartDataLabels);

    function genColors(size) {
        var o = [];
        for(i = 0; i < size; i++){
            o.push(`hsl(${ i * 360 / size }, 100%, 50%)`)
        }
        return o;
    }

    // Feedback over time
    function drawFdbkOverMonthGraph(id, datasets) {
        var fdbkOverMonthGraph = document.getElementById(id);
        var tutorialsBar = new Chart(fdbkOverMonthGraph, {
            type: 'line',
            data: {
                datasets: datasets
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
                            min: '2018-09',
                            unit: 'month',
                            distribution: 'linear'
                        }
                    }]
                },
                legend: {
                    display: false
                },
                title: {
                    display: false,
                },
                plugins: {
                    colorschemes: {
                        scheme: 'tableau.Classic20'
                    }
                }
            }
        });
    }

    function drawFdbkNoteGraph(id, data, labels) {
        var tutoBar = document.getElementById(id);
        var fdbkNoteGraph = new Chart(tutoBar, {
            type: 'horizontalBar',
            data: {
                labels: labels,
                datasets: [{
                    backgroundColor: genColors(data.length),
                    data: data
                }]
            },
            plugins: [ChartDataLabels],
            options: {
                legend: {
                    display: false
                },
                title: {
                    display: false,
                },
                plugins: {
                    datalabels: {
                        anchor: 'end',
                        align: 'end'
                    }
                }
            }
        });
    }
</script>

<!-- prepare data -->
{% assign fdbk_by_month = fdbk | group_by: 'month' %}
{% assign fdbk_over_time = "" | split: ',' %}
{% assign fdbk_over_time_labels = "" | split: ',' %}
{% assign fdbk_count = 0 %}
{% for f in fdbk_by_month %}
    {% assign new_fdbk = f.items | size %}
    {% assign fdbk_count = fdbk_count | plus: new_fdbk %}
    {% assign fdbk_over_time = fdbk_over_time | push: fdbk_count %}
    {% assign fdbk_over_time_labels = fdbk_over_time_labels | push: f.name %}
{% endfor %}
{% assign fdbk_by_note = fdbk | sort:'note' | group_by: 'note' %}
<script>
    var fdbk_over_month = {
        all_topics: [
            {
                label: "All topic",
                data: [{% for c in fdbk_over_time %}{x:"{{ fdbk_over_time_labels[forloop.index0] }}", y: {{ c }} } {% unless forloop.last %},{% endunless %}{% endfor %}],
                fill: false,
            }
        ]
    };
    var data_notes = {
        all_topics: [{% for n in (0..5) %}{% assign f_f_n = fdbk_by_note | where: "name",n | first %}{% if f_f_n.size == 0 %}0{% else %}{{ f_f_n.items | size }}{% endif %}{% unless forloop.last %},{% endunless %}{% endfor %}]
    };
    var labels_notes = [{% for n in (0..5) %}"{{ n }}"{% unless forloop.last %},{% endunless %}{% endfor %}];
</script>

{% assign fdbk_by_topic = fdbk | group_by: 'topic' | sort: 'name' %}
{% for tn in fdbk_by_topic %}
    {% assign t_metadata = topics | where: "title",tn.name | first %}
    {% assign topic_name = t_metadata['name'] %}
    {% assign fdbk_by_month = tn.items | group_by: 'month' %}
    {% assign fdbk_over_time = "" | split: ',' %}
    {% assign fdbk_over_time_labels = "" | split: ',' %}
    {% assign fdbk_count = 0 %}
    {% for f in fdbk_by_month %}
        {% assign new_fdbk = f.items | size %}
        {% assign fdbk_count = fdbk_count | plus: new_fdbk %}
        {% assign fdbk_over_time = fdbk_over_time | push: fdbk_count %}
        {% assign fdbk_over_time_labels = fdbk_over_time_labels | push: f.name %}
    {% endfor %}
    {% assign fdbk_by_note = tn.items | sort:'note' | group_by: 'note' %}
    <script>
        var topic_data =
        {
            label: "{{ tn.name }}",
            data: [{% for c in fdbk_over_time %}{x:"{{ fdbk_over_time_labels[forloop.index0] }}", y: {{ c }} } {% unless forloop.last %},{% endunless %}{% endfor %}],
            fill: false,
        };
        fdbk_over_month['topic-{{ topic_name }}'] = [topic_data];
        fdbk_over_month['all_topics'].push(topic_data);
        data_notes['topic-{{ topic_name }}'] = [{% for n in (0..5) %}{% assign f_f_n = fdbk_by_note | where: "name",n | first %}{% if f_f_n.size == 0 %}0{% else %}{{ f_f_n.items | size }}{% endif %}{% unless forloop.last %},{% endunless %}{% endfor %}];
    </script>

    {% assign fdbk_by_tutorial = tn.items | group_by: 'tutorial' | sort: 'name' %}
    {% for tun in fdbk_by_tutorial %}
        {% assign tuto_metadata = tutorials | where: "title",tun.name | first %}
        {% assign tuto_name = tuto_metadata['dir'] | split: "/" | last %}
        {% if tuto_metadata['name'] != 'tutorial.md' %}
            {% assign lang = tuto_metadata['name'] | split: "." | first | split: "-" | last %}
            {% assign tuto_name = tuto_name | append: lang %}
        {% endif %}
        {% if tuto_name != '' %}
            {% assign fdbk_by_month = tun.items | group_by: 'month' %}
            {% assign fdbk_over_time = "" | split: ',' %}
            {% assign fdbk_over_time_labels = "" | split: ',' %}
            {% assign fdbk_count = 0 %}
            {% for f in fdbk_by_month %}
                {% assign new_fdbk = f.items | size %}
                {% assign fdbk_count = fdbk_count | plus: new_fdbk %}
                {% assign fdbk_over_time = fdbk_over_time | push: fdbk_count %}
                {% assign fdbk_over_time_labels = fdbk_over_time_labels | push: f.name %}
            {% endfor %}
            {% assign fdbk_by_note = tun.items | sort:'note' | group_by: 'note' %}
            <script>
                var tuto_data =
                {
                    label: "{{ tun.name }}",
                    data: [{% for c in fdbk_over_time %}{x:"{{ fdbk_over_time_labels[forloop.index0] }}", y: {{ c }} } {% unless forloop.last %},{% endunless %}{% endfor %}],
                    fill: false,
                };
                fdbk_over_month['tuto-{{ tuto_name }}'] = [tuto_data];
                fdbk_over_month['topic-{{ topic_name }}'].push(tuto_data);
                data_notes['tuto-{{ tuto_name }}'] = [{% for n in (0..5) %}{% assign f_f_n = fdbk_by_note | where: "name",n | first %}{% if f_f_n.size == 0 %}0{% else %}{{ f_f_n.items | size }}{% endif %}{% unless forloop.last %},{% endunless %}{% endfor %}];
            </script>
        {% endif %}
    {% endfor %}
{% endfor %}


<div class="container main-content">
    <section>
        <h1>GTN Feedback</h1>

        <p>Aggregation of feedback submitted since {{ start.month }} using the embed feedback form
            at the bottom of tutorials. Thank you everyone who submitted feedback!
        </p>

        <div class="fdbk">
            <h2>Overall</h2>
            <div class="row">
                <!-- feedback over time  -->
                <div class="col-md-6">
                    <div class="card">
                        <div class="card-body">
                            <h5 class="card-title">{{ fdbk | size }} responses</h5>
                            {% assign chart_id = 'allFdbkOverMonthGraph' %}
                            <canvas id="{{ chart_id }}" width="400" height="400">
                                <script>
                                    drawFdbkOverMonthGraph("{{ chart_id }}", fdbk_over_month.all_topics);
                                </script>
                            </canvas>
                        </div>
                    </div>
                </div>
                <!-- feedback notes  -->
                <div class="col-md-6">
                    <div class="card">
                        <div class="card-body">
                            <h5 class="card-title">Rating Distribution</h5>
                            {% assign chart_id = 'allFdbkNoteGraph' %}
                            <canvas id="{{ chart_id }}" width="400" height="400">
                                <script>
                                    drawFdbkNoteGraph("{{ chart_id }}", data_notes.all_topics, labels_notes);
                                    var counter = 0;
                                </script>
                            </canvas>
                        </div>
                    </div>
                </div>
            </div>

            <h2>By topic</h2>
            {% for tn in fdbk_by_topic %}
                {% assign t_metadata = topics | where: "title",tn.name | first %}
                {% assign topic_name = t_metadata['name'] %}
                <h3 id="topic-{{ topic_name }}">{{ tn.name }}  <a href="{{ site.baseurl }}/topics/{{ topic_name }}">{% icon galaxy_instance %}</a></h3>
                <div class="row">
                    <!-- feedback over time  -->
                    <div class="col-md-6">
                        <div class="card">
                            <div class="card-body">
                                <h5 class="card-title">{{ tn.items | size }} responses</h5>
                                {% assign chart_id = topic_name | append: '-topic-fdbkOverMonthGraph' %}
                                <canvas id="{{ chart_id }}" width="400" height="400">
                                    <script>
                                        if ( counter > 4) {
                                            window.addEventListener('DOMContentLoaded', (event) => {
                                                drawFdbkOverMonthGraph("{{ chart_id }}", fdbk_over_month['topic-{{ topic_name }}']);
                                            });
                                        } else {
                                            drawFdbkOverMonthGraph("{{ chart_id }}", fdbk_over_month['topic-{{ topic_name }}']);
                                        }
                                    </script>
                                </canvas>
                            </div>
                        </div>
                    </div>
                    <!-- feedback notes  -->
                    <div class="col-md-6">
                        <div class="card">
                            <div class="card-body">
                                <h5 class="card-title">Rating Distribution</h5>
                                {% assign chart_id = topic_name | append: '-topic-fdbkNoteGraph' %}
                                <canvas id="{{ chart_id }}" width="400" height="400">
                                    <script>
                                        if ( counter > 4) {
                                            window.addEventListener('DOMContentLoaded', (event) => {
                                                drawFdbkNoteGraph("{{ chart_id }}", data_notes['topic-{{ topic_name }}'], labels_notes);
                                            });
                                        } else {
                                            drawFdbkNoteGraph("{{ chart_id }}", data_notes['topic-{{ topic_name }}'], labels_notes);
                                        }
                                        counter = counter + 1;
                                    </script>
                                </canvas>
                            </div>
                        </div>
                    </div>
                </div>

                {% assign fdbk_by_tutorial = tn.items | group_by: 'tutorial' %}
                <div class="accordion" id="accordion-{{ topic_name }}">
                {% for tun in fdbk_by_tutorial %}
                    {% assign tuto_metadata = tutorials | where: "title",tun.name | first %}
                    {% assign tuto_name = tuto_metadata['dir'] | split: "/" | last %}
                    {% if tuto_name == '' %}
                        <div class="alert alert-warning" role="alert">
                            <p>Tutorial "{{ tun.name }}" is not available anymore.</p>
                        </div>
                    {% else %}
                        {% if tuto_metadata['name'] != 'tutorial.md' %}
                            {% assign lang = tuto_metadata['name'] | split: "." | first %}
                            {% assign tuto_name = tuto_name | append: lang %}
                        {% endif %}
                        <div class="accordion-card card" id="tutorial-{{ tuto_name }}">
                            <div class="card-header" id="heading-{{ tuto_name }}">
                                <div class="mb-0">
                                    <button class="btn btn-link collapsed" type="button" data-toggle="collapse" data-target="#collapse-{{ tuto_name }}" aria-expanded="false" aria-controls="collapse-{{ tuto_name }}">
                                        {{ tun.name }}
                                    </button>
                                </div>
                            </div>

                            <div id="collapse-{{ tuto_name }}" class="collapse" aria-labelledby="heading-{{ tuto_name }}" data-parent="#accordion-{{ topic_name }}">
                                <div class="card-body">
                                    {% if tuto_metadata.enable == false %}
                                    <div class="alert alert-warning" role="alert">
                                        <p>This tutorial is not listed in the topic pages.</p>
                                    </div>
                                    {% endif %}
                                    <div class="row">
                                        <!-- feedback over time  -->
                                        <div class="col-md-6">
                                            <div class="card">
                                                <div class="card-body">
                                                    <h5 class="card-title">{{ tun.items | size }} responses</h5>
                                                    {% assign chart_id = tuto_name | append: '-tutorial-fdbkOverMonthGraph' %}
                                                    <canvas id="{{ chart_id }}" width="400" height="400">
                                                        <script>
                                                            window.addEventListener('DOMContentLoaded', (event) => {
                                                                drawFdbkOverMonthGraph("{{ chart_id }}", fdbk_over_month['tuto-{{ tuto_name }}']);
                                                            });
                                                        </script>
                                                    </canvas>
                                                </div>
                                            </div>
                                        </div>
                                        <!-- feedback notes  -->
                                        <div class="col-md-6">
                                            <div class="card">
                                                <div class="card-body">
                                                    <h5 class="card-title">Rating Distribution</h5>
                                                    {% assign chart_id = tuto_name | append: '-tutorial-fdbkNoteGraph' %}
                                                    <canvas id="{{ chart_id }}" width="400" height="400">
                                                        <script>
                                                            window.addEventListener('DOMContentLoaded', (event) => {
                                                                drawFdbkNoteGraph("{{ chart_id }}", data_notes['tuto-{{ tuto_name }}'], labels_notes);
                                                            });
                                                        </script>
                                                    </canvas>
                                                </div>
                                            </div>
                                        </div>
                                    </div>

                                    <div class="row">
                                        <div class="col">
                                            <div class="card">
                                                <div class="card-body">
                                                    <h5 class="card-title">Detailed feedback</h5>
                                                    <table class="table table-striped">
                                                        <thead>
                                                            <tr>
                                                                <th>Date</th>
                                                                <th>What did you like?</th>
                                                                <th>What could be improved?</th>
                                                            </tr>
                                                        </thead>
                                                        <tbody>
                                                            {% assign d_f = tun.items | reverse %}
                                                            {% for f in d_f %}
                                                                {% if f.anonymous != "yes" %}
                                                                    {% if f.pro != null or f.con != null %}
                                                                    <tr>
                                                                        <td>{{ f.date | date: "%b %-d, %Y"  }}</td>
                                                                        <td>{{ f.pro }}</td>
                                                                        <td>{{ f.con }}</td>
                                                                    </tr>
                                                                    {% endif %}
                                                                {% endif %}
                                                            {% endfor %}
                                                        </tbody>
                                                    </table>
                                                </div>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    {% endif %}
                {% endfor %}
                </div>
            {% endfor %}
        </div>

    </section>
</div>
