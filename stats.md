---
layout: base
---
{% include _includes/default-header.html %}




{% assign contributors = site.data['contributors'] %}

<style type="text/css" media="all">
.card-title {
	font-size: 2rem;
}

</style>
<div class="container main-content">
	<section>
		<h1>Statistics</h1>

		<div class="stats">
			<div class="row">
				<div class="card" style="width: 60%;margin:1em">
					<div class="card-body">
						<h5 class="card-title">{{ site.data | science_topics | size }} Scientific Topics</h5>
						<p class="card-text">Covering a wide range of scientific domains: {{ site.data | science_topic_links | join: ', ' }}</p>

						<!--<a href="#" class="card-link">Card link</a>-->
						<!--<a href="#" class="card-link">Another link</a>-->
					</div>
				</div>

				<div class="card" style="width: 30%;margin:1em">
					<div class="card-body">
						<h5 class="card-title">Tutorials</h5>
						<p class="card-text">amazing!.</p>
						<!--<a href="#" class="card-link">Card link</a>-->
						<!--<a href="#" class="card-link">Another link</a>-->
					</div>
				</div>



<script type="text/javascript">
	data = {{ site.data |jsonify }}
</script>
			</div>
		</div>
	</section>
</div>






{% include _includes/default-footer.html %}
