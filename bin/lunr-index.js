var lunr = require("lunr"),
	metadataParser = require("markdown-yaml-metadata-parser"),
	path = require("path"),
	fs = require("fs"),
	stdout = process.stdout,
	properties = ["tags", "questions", "objectives", "key_points"];

function findResources(startPath) {
	if (!fs.existsSync(startPath)) {
		return;
	}
	var resources = [];

	var files = fs.readdirSync(startPath);
	for (var i = 0; i < files.length; i++) {
		var filename = path.join(startPath, files[i]);
		var stat = fs.lstatSync(filename);
		if (stat.isDirectory()) {
			findResources(filename).forEach(x => {
				resources.push(x);
			});
		} else {
			if(filename.indexOf('.') != 0 && filename.indexOf('_site') != 0 && !stat.isSymbolicLink()){
				resources.push(filename);
			}
		}
	}

	return resources;
}

var resources = findResources('.');
tutorials = resources.filter(x => x.indexOf('topics/') == 0).filter(x => x.endsWith('tutorial.md'))
faqs = resources.filter(x => x.indexOf('faqs/') > 0 ).filter(x => x.endsWith('.md'))


tutorials = tutorials.map(tuto => {
	const contents = fs.readFileSync(tuto, "utf8");
	const result = metadataParser(contents);
	result.metadata.id = tuto;
	result.metadata.content = result.content;
	filtered = result.content.split('\n')
		.filter(x => x[0] !== '>') // Remove boxes
		.filter(x => x[0] !== '{') // remove closing tag
		.filter(x => x != '') // don't index whitespace
		//.slice(0, 25) // first 25 lines
		.join('\n')

	output = {
		id: tuto,
		type: "tutorial",
		title: result.metadata.title,
		content: filtered,
		contributors: (result.metadata.contributors || []).join(" "),
		meta: "",
	};

	// Array style properties
	properties.forEach(prop => {
		if (Array.isArray(result.metadata[prop])) {
			output['meta'] += result.metadata[prop].join(" ") + "\n";
		}
	});
	return output;
});

faqs = faqs.map(faqPath => {
	const contents = fs.readFileSync(faqPath, "utf8");
	const result = metadataParser(contents);
	return [result, faqPath];
}).filter(x => {
	return x[0].metadata.layout === 'faq'
}).map(faqdata => {
	result = faqdata[0];

	output = {
		id: faqdata[1],
		type: "faq",
		title: result.metadata.title,
		content: result.content,
		contributors: (result.metadata.contributors || []).join(" "),
		meta: "",
	};

	// Array style properties
	properties.forEach(prop => {
		if (Array.isArray(result.metadata[prop])) {
			output['meta'] += result.metadata[prop].join(" ") + "\n";
		}
	});
	return output;
});

var idx = lunr(function() {
	this.ref("id");
	this.field("title");
	this.field("meta");
	this.field("contributors");
	this.field("content");

	tutorials.forEach(function(doc) {
		this.add(doc);
	}, this);

	faqs.forEach(function(doc) {
		this.add(doc);
	}, this);
});

stdout.write(JSON.stringify(idx));
