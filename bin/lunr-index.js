var lunr = require("lunr"),
	metadataParser = require("markdown-yaml-metadata-parser"),
	path = require("path"),
	fs = require("fs"),
	stdout = process.stdout,
	properties = ["tags", "questions", "objectives", "key_points", "contributors"];

function findTutorials(startPath, filter) {
	if (!fs.existsSync(startPath)) {
		return;
	}
	var tutorials = [];

	var files = fs.readdirSync(startPath);
	for (var i = 0; i < files.length; i++) {
		var filename = path.join(startPath, files[i]);
		var stat = fs.lstatSync(filename);
		if (stat.isDirectory()) {
			findTutorials(filename, filter).forEach(x => {
				tutorials.push(x);
			});
		} else if (filename.indexOf(filter) >= 0) {
			tutorials.push(filename);
		}
	}

	return tutorials;
}

tutorials = findTutorials("topics", "tutorial.md").map(tuto => {
	const contents = fs.readFileSync(tuto, "utf8");
	const result = metadataParser(contents);
	result.metadata.id = tuto;
	result.metadata.content = result.content;
	filtered = result.content.split('\n')
		.filter(x => x[0] !== '>') // Remove boxes
		.filter(x => x[0] !== '{') // remove closing tag
		.slice(0, 25) // first 25 lines
		.join('\n')

	output = {
		id: tuto,
		title: result.metadata.title,
		content: filtered,
	};

	// Array style properties
	properties.forEach(prop => {
		if (Array.isArray(result.metadata[prop])) {
			output[prop] = result.metadata[prop].join(" ");
		}
	});
	return output;
});

var idx = lunr(function() {
	this.ref("id");
	this.field("title");
	this.field("tags");
	this.field("questions");
	this.field("objectives");
	this.field("key_points");
	this.field("contributors");
	this.field("content");

	tutorials.forEach(function(doc) {
		this.add(doc);
	}, this);
});

stdout.write(JSON.stringify(idx));
