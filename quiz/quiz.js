/**
This code was originally licensed under MIT, and then heavily modified. The MIT license is retained below.

MIT License

Copyright (c) 2017 Jack McKernan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 * */

var lastPeerId = null,
	peer = null, // Own peer object
	peerId = null,
	conns = [],
	conn = null,
	recvId = document.getElementById("receiver-id"),
	status = document.getElementById("status"),
	message = document.getElementById("message"),
	advanceButton = document.getElementById("advance"),
	questionArea = document.getElementById("questionArea"),
	lobby = document.getElementById("lobby"),
	progressBar = document.getElementById("progress"),
	displayTime = 5000;

// Teacher Only
var players = {};
var slideTimer;
var slides, quiz_title;
var currentSlide = -1;

// Student only
var player = {};
var answersGiven = {};


function getUrlParam(name){
	return (new URL(document.location)).searchParams.get(name);
}

const roomIdBase = 'e8a2e4ea-galaxy-training-network-';
const mode = getUrlParam('mode');

if(mode === 'teacher' || mode === 'self') {
	advanceButton.style.visibility = 'visible';
} else {
	progressBar.style.visibility = 'hidden';
}
if(mode === 'self'){
	displayTime = 0;
}

function generateRoomId(){
	var roomNumber = (Math.random() * 10000000).toString().substring(0, 4),
		roomId = roomIdBase + roomNumber;
	return [roomNumber, roomId]
}

function updateQrCode(roomNumber, roomId){
	var docurl = new URL(document.location);
	var studentUrl = docurl.origin + docurl.pathname.replace('teacher', 'student') +  `?id=${roomNumber}`
	var qrcode = new QRious({
		element: document.getElementById("qrcode"),
		background: '#ffffff',
		backgroundAlpha: 1,
		foreground: '#5868bf',
		foregroundAlpha: 1,
		level: 'H',
		padding: 0,
		size: 300,
		value: studentUrl
	});
	document.getElementById("join-url").innerHTML = `<div><a href="${studentUrl}">Copy This link</a></div><div class='room-number'>${roomNumber}</div>`
}


function postInit(){
	if(mode === 'self'){
		// DO SOMETGIN
	} else if(mode === 'teacher'){
		status.innerHTML = "Awaiting connections...";
	} else {
		var urlOption = getUrlParam("id");
		if(urlOption !== null){
			//recvIdInput.value = urlOption;
			join(null, urlOption);
			//setTimeout(function() { join(null, urlOption) }, 500);
		} else {
			showRoomCode();
		}
	}
}

function peerOn(peer, id){
	// Workaround for peer.reconnect deleting previous id
	if (peer.id === null) {
		console.log('Received null id from peer open');
		peer.id = lastPeerId;
	} else {
		lastPeerId = peer.id;
	}

	console.log('ID: ' + peer.id, mode);
	postInit();
}

function peerClose(){
	conn = null;
	status.innerHTML = "Connection destroyed. Please refresh";
	console.log('Connection destroyed');
}
function peerError (err) {
	console.log(err);
	questionArea.innerHTML = err;
}

function peerConnect(c){
	console.log("Peer connected", c, mode)
	if(mode === 'teacher'){
		conns.push(c);
		players[c.peer] = {}
		updateStudentList();
		ready();

		setTimeout(function(){
			c.send({type: 'setup', title: quiz_title});
		}, 1000);
	} else {
		// Disallow incoming connections
		c.on('open', function() {
			c.send("Sender does not accept incoming connections");
			setTimeout(function() { c.close(); }, 500);
		});
	}
}

function peerDisconnect() {
	status.innerHTML = "Connection lost. Please reconnect";
	console.log('Connection lost. Please reconnect');

	// Workaround for peer.reconnect deleting previous id
	peer.id = lastPeerId;
	peer._lastServerId = lastPeerId;
	peer.reconnect();
}

function showWelcome(){
	questionArea.innerHTML = `
		<div class="display: flex; align-items: center; justify-content: center;">
			<h1>Enter your name</h1>
			<input type="text" id="name">
			<button id="submit-name" class="btn btn-primary">Join!</button>
		</div>
	`

	var name_input = document.getElementById(`submit-name`);
	name_input.addEventListener('click', () => {
		player.name = document.getElementById("name").value;
		safeSend({
			"event": "registerPlayer",
			"player": player,
		})
		showPostWelcome();
	})
}

function showPostWelcome(){
	questionArea.innerHTML = `
		<div class="display: flex; align-items: center; justify-content: center;">
			<h1>Get Ready ${player.name}</h1>
		</div>
	`
}

/**
 * Create the connection between the two Peers.
 *
 * Sets up callbacks that handle any events related to the
 * connection and data received on it.
 */
function join(evt, joinId) {
	// Close old connection
	if (conn) {
		conn.close();
	}

	// Create connection to destination peer specified in the input field
	var connectToId = joinId || recvIdInput.value
	console.log('connecting to e8a2e4ea-galaxy-training-network-' + connectToId)
	conn = peer.connect('e8a2e4ea-galaxy-training-network-' + connectToId, {
		reliable: true
	});

	conn.on('open', function () {
		status.innerHTML = "Connected to: gtn-" + conn.peer.substring(33);
		console.log("Connected to: " + conn.peer);

		// Check URL params for comamnds that should be sent immediately
		var command = getUrlParam("command");
		if (command)
			conn.send(command);
	});
	// Handle incoming data (messages only since this is the signal sender)
	conn.on('data', function (data) {
		if(data.type === "choose-1"){
			showQuestion(data);
		} else if(data.type === "choose-many"){
			showQuestion(data);
		} else if (data.type === "poll") {
			showQuestion(data);
		} else if (data.type === "free-text") {
			showQuestion(data);
		} else if (data.type === "setup") {
			document.getElementById("title").innerHTML = data.title;
			document.getElementsByTagName("title")[0].innerHTML = data.title;
		} else if (data.type === "clear") {
			questionArea.innerHTML = `<h2>Question ${data.question + 1}</h2>`;
		} else if (data.type === "answer") {
			if(Array.isArray(data.answer)){
				numCorrect = 0;
				data.answer.forEach(correctAnswer => {
					if(answersGiven[data.question].indexOf(correctAnswer) > -1){
						numCorrect = numCorrect + 1;
					}
				})

				questionArea.innerHTML = `
					<h2>You got ${numCorrect}/${data.answer.length} Correct</h2>
				`
			} else {
				if(data.answer === answersGiven[data.question]){
					// correct!
					questionArea.innerHTML = `
						<h2>Congratulations</h2>
					`
				} else {
					console.log(data)
					if(data.slide_type === "free-text") {
						questionArea.innerHTML = `
							<h2>Answer</h2>
							${data.answer}
						`
					}
					else {
						questionArea.innerHTML = `
							<h2>Too bad :(</h2>
						`
					}
				}
			}
		} else {
			console.log("Unknown message")
		}
	});
	conn.on('close', function () {
		status.innerHTML = "Connection closed";
	});

	showWelcome();
};

/**
 * Send a signal via the peer connection and add it to the log.
 * This will only occur if the connection is still alive.
 */
 function signal(sigName) {
	if (conn && conn.open) {
		conn.send(sigName);
		console.log(sigName + " signal sent");
	} else {
		console.log('Connection is closed');
	}
}

function safeSend(msg){
	if(mode === 'self'){
		console.log('SELF', msg);
		processStudentMessage('SELF', msg);
		return;
	}

	console.log("Sending", msg)
	if (conn && conn.open) {
		conn.send(msg);
	} else {
		console.log('Connection is closed');
	}
}

function showResult(result, score){
	message.innerHTML = `
		<div class="display: flex; align-items: center; justify-content: center;">
			<h1>${result}</h1>
			<div>
				Your score is: ${score}
			</div>
		</div>
	`
}

function showQuestion(data){
	if(data.type === "choose-many"){
		showQuestionMany(data)
	} else if(data.type === "free-text"){
		showQuestionText(data)
	} else {
		showQuestionOne(data);
	}
}

function showQuestionMany(data){
	var show = '';
	if(data.image){
		show += `<div style="display: flex;flex-direction: column"><h2>${data.title}</h2><img src="${data.image}" /></div>`
		questionArea.style['flex-wrap'] = 'unset';
	} else {
		show += `<h2>${data.title}</h2>`
		questionArea.style['flex-wrap'] = 'wrap';
	}

	show += '<p>Choose as many as you think are applicable</p>'
	show += `<div class="answer-group">`;
	show += data.answers.map((q, idx) => {
		return `
			<label for="answer-${data.id}-${idx}-r" class="btn answer-button">
			<div id="answer-${data.id}-${idx}">
			<input type="checkbox" id="answer-${data.id}-${idx}-r" name="radio-answer" value="${q}">
				${q}
			</div>
			</label>
		`
	}).join("");
	show += '</div>';
	show += `<div><button id="answer-${data.id}-submit" class="btn btn-primary">Submit</button></div>`
	questionArea.innerHTML = show;

	// Only the teacher doesn't need them to be clickable.
	if(mode !== 'teacher'){
		var e = document.getElementById(`answer-${data.id}-submit`);
		console.log(e)
		e.addEventListener('click', function () {
			var answers = Array.from(document.querySelectorAll("input[name='radio-answer']")).filter(x => x.checked).map(x => x.value);
			safeSend({
				"event": "answer",
				"question": data.id,
				"result": answers,
			})
			answersGiven[data.id] = answers
			console.log(data);
			if(data.type !== "poll" || ! data.live) {
				Array.from(document.getElementsByClassName("answer-button")).forEach(x => x.style.display = 'none')
			}
		})
	}
}

function showQuestionOne(data){
	var show = '';
	if(data.image){
		show += `<div style="display: flex;flex-direction: column"><h2>${data.title}</h2><img src="${data.image}" /></div>`
		questionArea.style['flex-wrap'] = 'unset';
	} else {
		show += `<h2>${data.title}</h2>`
		questionArea.style['flex-wrap'] = 'wrap';
	}

	show += `<div class="answer-group">`;
	if(data.answers !== undefined){
		show += data.answers.map((q, idx) => {
			return `
				<button id="answer-${data.id}-${idx}" value="${q}" class="btn answer-button">${q}</button>
			`
		}).join("");
	}
	show += '</div>';
	questionArea.innerHTML = show;

	// Only the teacher doesn't need them to be clickable.
	if(mode !== 'teacher'){
		data.answers.forEach((q, idx) => {
			var e = document.getElementById(`answer-${data.id}-${idx}`);
			e.addEventListener('click', function () {
				safeSend({
					"event": "answer",
					"question": data.id,
					"result": e.value,
				})
				answersGiven[data.id] = e.value
				console.log(data);
				if(data.type !== "poll" || ! data.live) {
					Array.from(document.getElementsByClassName("answer-button")).forEach(x => x.style.display = 'none')
				}
			})
		})
	}
}

function showQuestionText(data){
	var show = '';
	if(data.image){
		show += `<div style="display: flex;flex-direction: column"><h2>${data.title}</h2><img src="${data.image}" /></div>`
		questionArea.style['flex-wrap'] = 'unset';
	} else {
		show += `<h2>${data.title}</h2>`
		questionArea.style['flex-wrap'] = 'wrap';
	}

	if(mode !== 'teacher'){
		show += `<div class="answer-group"><textarea rows="8" cols="32" id="answer-${data.id}-text"></textarea></div>`;
		show += `<div class="answer-group">`;
		show += `<button id="answer-${data.id}-submit" value="answered" class="btn answer-button">Submit</button>`
		show += `</div>`;
	}
	questionArea.innerHTML = show;

	if(mode !== 'teacher'){
		var e = document.getElementById(`answer-${data.id}-submit`);
		console.log(e)
		e.addEventListener('click', function () {
			var answer = document.getElementById(`answer-${data.id}-submit`).value
			safeSend({
				"event": "answer",
				"question": data.id,
				"result": answer,
			})
			answersGiven[data.id] = answer
			document.getElementById(`answer-${data.id}-submit`).disabled = true
			document.getElementById(`answer-${data.id}-text`).disabled = true
			console.log(data);
		})
	}
}

function updateStudentList(){
	status.innerHTML = `${conns.length} connections`;
	if(currentSlide === -1) {
		lobby.innerHTML =
			Object.keys(players).map(playerId => {
				return `<div class="player-name">${playerName(playerId)}</div>`
			}).join(" ") + '</div>';
	}
}

function processStudentMessage(connId, message){
	console.log(connId, message)
	if(message.event == "registerPlayer"){
		players[connId] = {
			"name": message.player.name,
		}
		updateStudentList();
	} else if (message.event === "answer") {
		if(message.question !== currentSlide){
			console.log("Attempting to answer wrong question! Ignoring.")
			return
		}

		if(slides[currentSlide].type === "choose-many"){
			slides[currentSlide].results[connId] = message.result.map(ans => slides[currentSlide].answers.indexOf(ans))
		} else if(slides[currentSlide].type === "free-text"){
			slides[currentSlide].results[connId] = message.result // TODO: will people submit crappy answers here?
		} else {
			slides[currentSlide].results[connId] =
				slides[currentSlide].answers.indexOf(message.result)
		}

		if(slides[currentSlide].live){
			showResults(slides[currentSlide])
		}

		console.log(slides[currentSlide])
	} else {
		console.log("Unknown message", message);
	}
}

function chunkArray(array, chunkSize){
	var chunks = []
	for (let i = 0; i < array.length; i += chunkSize) {
		chunks.push(array.slice(i, i + chunkSize));
	// do whatever
	}
	return chunks;
}

function broadcast(msg){
	console.log("Broadcast", msg)
	conns.filter(conn => conn && conn.open).forEach(conn => {
		conn.send(msg);
	})
}

function showResults(){
	var slide = slides[currentSlide];
	var show = '<h1>Results</h1>';
	var counts = {}
	var final_count = 0;
	Object.keys(slide.results).forEach(connId => {
		var tmp = slide.results[connId];

		// I'm sorry.
		var asdf = [];
		if(Array.isArray(tmp)){
			asdf = tmp
		} else {
			asdf = [tmp]
		}

		asdf.forEach(theirAnswer => {
			var answerKey = "";
			if(slide.type !== 'free-text'){
				if(theirAnswer < 0){
					answerKey = "SOMETHING ODD";
				} else {
					answerKey = slide.answers[theirAnswer]
					final_count += 1;
				}
			} else {
				// We'll give everyone a pass.
				answerKey = slide.correct
				final_count += 1;
			}

			console.log(theirAnswer, answerKey, slide.correct)
			if(answerKey === slide.correct){
				console.log(players[connId])
				players[connId]['score'] = 1 + (players[connId]['score'] || 0);
			}

			counts[answerKey] = 1 + (counts[answerKey] || 0)
		})

	})
	slide.final_results = counts;
	slide.final_count = final_count;
	console.log(slide)
	console.log(players)

	show += '<table class="table table-striped">'
	if(slide.type === 'free-text'){
			show += `<tr class="correct-answer"><td>${slide.correct}</td> <td>n/a</td></tr>`
	} else {
		slide.answers.forEach(x => {
			show += `<tr ${isCorrectAnswer(x) ? 'class="correct-answer"' : ''}><td>${x}</td> <td><div class="bar-chart" style="width: ${25 * (counts[x] || 0) / final_count}em">${counts[x] || 0}</div></td></tr>`
		})
	}
	show += '</table>'
	questionArea.innerHTML = show;
}

function isCorrectAnswer(x){
	if( Array.isArray(slides[currentSlide].correct)){
		return slides[currentSlide].correct.indexOf(x) >= 0
	} else {
		return slides[currentSlide].correct === x
	}
}

function renderTable(arr, headers){
	show = '<table class="table table-striped">'
	if(headers !== undefined){
		show += '<thead>'
		show += "<tr>" + headers.map(col => `<th>${col}</th>`).join("") + "</tr>"
		show += '</thead>'
	}
	show += '<tbody>'
	show += arr.map(row => {
		return "<tr>" + row.map(col => `<td>${col}</td>`).join("") + "</tr>"
	}).join("")
	show += '</tbody>'
	show += '</table>'
	show += '</div>'
	return show
}

function playerName(playerId){
	return players[playerId].name || "Anonymous " + playerId.substr(0, 6)
}

function showFinalResults(){
	var show = '<h1>Final Results</h1>';
	var counts = {}
	// Ensure they have a score.
	var playerIds = Object.keys(players);
	playerIds.forEach(playerId => { players[playerId].score = players[playerId].score || 0 })

	// Sort them
	playerIds.sort(function(a, b){return players[b].score - players[a].score});

	// Display
	show += renderTable(playerIds.map(playerId => {
		return [
			playerName(playerId),
			players[playerId].score || 0
		]
	}).slice(0, 3), ['Name', 'Score'])
	show += '</table>'


	// Difficult questions
	show += '<h2>Difficult Questions</h2>';
	// What's our def here? Top 3 worst?
	worstQuestions = slides
		.filter(s => s.type !== "poll")
		.map(s => {
			var pc = 0;
			if (!s.final_results || !s.final_results[s.correct]) {
				pc = 0;
			} else {
				if(s.final_count == 0){
					pc = 0;
				} else {
					pc = s.final_results[s.correct] / s.final_count;
				}
			}

			return {
				"title": s.title,
				"correct": s.correct,
				"final_results": s.final_results,
				"pc": pc,
				"percent_correct": (pc * 100).toFixed(1) + '%',
			}
		})
		.filter(s => s.pc < 1)

	// Sort by % correct
	worstQuestions.sort((a, b) => b.percent_correct < a.percent_correct)

	// Render
	show += renderTable(worstQuestions.slice(0, 3).map(wq => {
		return [
			wq.title,
			wq.percent_correct,
			JSON.stringify(wq.final_results)
			//renderTable(Object.keys(wq.final_results).map(x => [x, wq.final_results[x]]))
		]
	}), ['Question', '% Correct', 'Answers'])

	show += '</div>'

	console.log(worstQuestions)


	questionArea.innerHTML = show;
}

function handleCurrentSlide(){
	var studentSlide = {
		id: currentSlide,
		type: slides[currentSlide].type,
		title: slides[currentSlide].title,
		answers: slides[currentSlide].answers,
		image: slides[currentSlide].image,
		started: new Date().getTime(),
		timeout: mode === 'self' ? 600 : slides[currentSlide].timeout,
		live: slides[currentSlide].live,
	}

	showQuestion(studentSlide);
	broadcast({type: 'clear', question: currentSlide})

	var haveBroadcast = false;

	// Update the count down every 1 second
	slideTimer = setInterval(function() {
		var now = new Date().getTime(),
			showQuestionLeft = studentSlide.started + displayTime - now
			timeLeft = studentSlide.started + (studentSlide.timeout * 1000) + displayTime - now;

		if(showQuestionLeft < 0){
			if(!haveBroadcast){
				haveBroadcast = true;
				broadcast(studentSlide)
				if(currentSlide.answers !== null) {
					document.getElementsByClassName("answer-group")[0].style.display = '';
				}
			}
		} else {
			progressBar.style.width = ((displayTime - showQuestionLeft) / 50) + "vw"
			return;
		}

		var doneCondition = timeLeft < 0;
		// How many students have answered?
		if(Object.keys(players).length === Object.keys(slides[currentSlide].results).length && ! slides[currentSlide].live){
			doneCondition = true
		}

		// Check the time left
		if(doneCondition){
			progressBar.style.width = "100vw"
			progressBar.innerHTML = "&nbsp;";
			clearInterval(slideTimer);
			showResults(studentSlide);
			if(studentSlide.type !== "poll"){
				broadcast({
					type: "answer",
					question: currentSlide,
					answer: slides[currentSlide].correct,
					slide_type: slides[currentSlide].type,
				});
			}
		} else {
			progressBar.style.width = timeLeft / studentSlide.timeout / 10 + "vw"
			progressBar.innerHTML = Math.round(timeLeft / 1000, 2);
		}
	}, 25);

}

function ready() {
	conns
		.filter(conn => {
			return Object.keys(conn._events).indexOf('data') == -1
		})
		.forEach(conn => {
			console.log(conn.peer, Object.keys(conn._events))
			conn.on('data', function (data) {
				console.log("Data recieved");
				console.log(data)
				processStudentMessage(conn.peer, data);
				var cueString = "<span class=\"cueMsg\">Cue: </span>";
			});
		})

	conns
		.filter(conn => {
			return Object.keys(conn._events).indexOf('data') == -1
		})
		.forEach(conn => {
			conn.on('close', function () {
				status.innerHTML = "Connection reset<br>Awaiting connection...";
				conns = conns.filter(c => c.peer != conn.peer);
			});
		});
}

function advanceListener(){
	// Remove any existing timer if there is one.
	clearInterval(slideTimer);
	advanceButton.innerHTML = 'Next Question'
	currentSlide += 1;
	if(currentSlide === slides.length - 1){
		advanceButton.innerHTML = 'Final Results'
	}
	if(currentSlide === slides.length){
		advanceButton.style.display = 'none';
		showFinalResults();
		return
	}
	handleCurrentSlide();
}

function hideJoin(){
	document.getElementById("connect-area").style.display = 'none';
}

function showDebug(){
	document.getElementsByClassName("debug-info")[0].style.display = '';
}

function loadQuiz(url){
	data = YAML.load(url);
	slides = data.questions.map(x => {
		x.results = {};
		return x
	});

	quiz_title = data.title;
	document.getElementById("title").innerHTML = data.title;
	document.getElementsByTagName("title")[0].innerHTML = data.title;
}

function showRoomCode(){
	questionArea.innerHTML = `
		<div class="display: flex; align-items: center; justify-content: center;">
			<h1>Enter Your Room Number</h1>
			<input type="text" id="roomnumber">
			<button id="submit-room" class="btn btn-primary">Join!</button>
		</div>
	`

	var name_input = document.getElementById(`submit-room`);
	name_input.addEventListener('click', () => {
		join(null, document.getElementById("roomnumber").value);
	})
}

(function () {
	 function initialize() {
		var docurl = new URL(document.location);
		var roomId;

		if(mode === 'self' || mode === 'teacher'){
			loadQuiz(docurl.searchParams.get('quiz'))
			var [roomNumber, roomId] = generateRoomId();
		}

		if(mode == 'self'){
			// Fake a 'self' connection
			players['SELF'] = {}
			safeSend({event: 'registerPlayer', player: {name: 'Self Study'}})
			return;
		}

		peer = new Peer(roomId, {
			debug: 2,
			//host: 'localhost',
			//port: 9000,
			//path: '/'
		});

		if(mode === 'teacher'){
			updateQrCode(roomNumber, roomId);
		}

		peer.on('open', (id) => peerOn(peer, id));
		peer.on('connection', (c) => peerConnect(c));
		peer.on('disconnected', peerDisconnect);
		peer.on('close', peerClose);
		peer.on('error', (err) => peerError(err));
	};

	advanceButton.addEventListener('click', advanceListener);

	initialize();
})();
