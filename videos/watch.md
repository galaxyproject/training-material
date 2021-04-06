---
layout: base
---
{% include _includes/default-header.html %}

<div class="container main-content">
	<div class="row">
		<div class="col-md-12">
			<video id="player" width="100%" height="610" controls preload="metadata" style="background: black">
			</video>
		</div>
	</div>
	<div class="row">
		<div class="col-sm-8">
			<div class="row">
				<h2 id="title"></h2>
			</div>
			<div class="row">
				<div class="col-sm-6">
					<h3>Transcript</h3>
				</div>
				<div class="col-sm-6" id="transcript-edit">
					Edit
				</div>
			</div>
			<div class="row" id="transcript">
			</div>
		</div>
		<div class="col-sm-4">
			<div class="row">
				<div class="col-sm-12">
					<h3>Other Videos</h3>
					<div><a href="{% link videos/index.md %}">See all GTN Videos</a></div>
					<div id="playlist" class="vertical">
					</div>
				</div>
			</div>
		</div>
	</div>
</div>


<script type="text/javascript">
var params = (new URL(document.location)).searchParams,
	videoid = params.get('v'),
	seekTo = params.get('t'),
	videohost = 'http://localhost:4002/training-material/',
	vtt = `${videohost}/videos/topics/${videoid}.en.vtt`,
	mp4 = `${videohost}/videos/topics/${videoid}.mp4`,
	png = `${videohost}/videos/topics/${videoid}.mp4.png`,
	player = document.getElementById("player");
	//videohost = 'https://training.galaxyproject.org';


player.setAttribute('poster', png);
player.innerHTML = `
	<source src="${mp4}" type="video/mp4">
	<track label="English" kind="captions" srclang="en" src="${vtt}" default>
`;

document.getElementById("transcript-edit").innerHTML = `<a href="https://github.com/galaxyproject/training-material/edit/main/topics/${videoid}.html">Found a typo? Edit</a>`

if(seekTo !== null){
	if(seekTo.indexOf(":") > -1){
		var seekToparts = seekTo.split(":");
		if(seekToparts.length == 2) {
			player.currentTime = (parseInt(seekToparts[0]) * 60) + parseInt(seekToparts[1]);
		} else if (seekToparts.length == 3){
			player.currentTime = (parseInt(seekToparts[0]) * 3600) + (parseInt(seekToparts[1]) * 60) + parseInt(seekToparts[2]);
		} else {
			console.error("Could not parse time")
		}
	} else {
		player.currentTime = parseInt(seekTo);
	}
}


fetch(vtt)
	.then(response => response.text())
	.then(data => {
		lines = data.split("\n").slice(4).filter((x, i) => { return i % 4 == 0 || i % 4 == 1});

		timestamps = lines.filter((x, i) => i % 2 == 0).map(x => x.split(' ')[0]);
		words = lines.filter((x, i) => i % 2 == 1);

		var zipped = timestamps.map(function(e, i) {
			return [e, words[i]];
		});
		lines = zipped.map(x => { return `<tr><td>${x[0]}</td><td>${x[1]}</td></tr>` }).join('');
		document.getElementById("transcript").innerHTML = '<table>' + lines + '</table>';
	});

fetch('{{ site.baseurl }}/api/videos.json')
	.then(response => response.json())
	.then(data => {
		// Remove empty
		data = data.filter(x => x !== null);
		// We've got a 'list' of video, we'll pretend this is a 'ring' buffer.

		var idx = data.findIndex(x => x.id == videoid);
		var videoSelf = data[idx];
		document.getElementById("title").innerHTML = videoSelf.title;


		var ring = [...data.slice(idx + 1), ...data.slice(0, idx)].slice(0, 8);
		var fmt = ring.map(x => {
			return `
			<div class="pl-item">
				<a href="?v=${x.id}">
					<div class="cover">
						<img src="{{site.baseurl}}/videos/topics/${x.id}.mp4.png" width="200px"/>
					</div>
					<div>
						<div class="title">${x.title}</div>
						<div class="topic">${x.topic}</div>
					</div>
				</a>
			</div>
			`;
		})
		document.getElementById("playlist").innerHTML = fmt;
	});

</script>
{% include _includes/default-footer.html %}
