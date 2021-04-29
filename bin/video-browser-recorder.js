const { chromium } = require('playwright');
const { saveVideo } = require('playwright-video');
const fs = require('fs');
var actions;
var syncReport = [];
fs.readFile(process.argv[2], 'utf8', (err, data) => {
	actions = JSON.parse(data);
});
var video_output_name = process.argv[3];
var videoSpeed = 1000;
if(process.argv.length > 4){
	videoSpeed = 10;
}


function logtime(now, start, msg){
	var timestamp = now.getTime() - start.getTime();
	syncReport.push({
		'time': timestamp,
		'msg': msg,
	})
	//console.log(timestamp/1000, msg);
}

(async () => {
	var start = new Date();
	// this seems to be ignored? no way to set zoom?
	const browser = await chromium.launch({args: ["--force-device-scale-factor=1.5"]});
	const context = await browser.newContext({ignoreHTTPSErrors: true});
	const page = await context.newPage();
	await page.setViewportSize({
		width: 1920,
		height: 1080,
	});
	await saveVideo(page, video_output_name);

	for(var i = 0; i < actions.length; i++){
		var step = actions[i];
		//console.log(step);
		if(step.action == 'goto'){
			await page.goto(step.target);
			await page.waitForLoadState('networkidle');
			now = new Date();
			logtime(now, start, 'gone')
		} else if (step.action == 'scrollTo'){
			await page.evaluate((step) => document.getElementById(step.target.slice(1)).scrollIntoView({behavior: "smooth"}), step).catch((err) => console.log(err));
			now = new Date();
			logtime(now, start, 'scrolled')
		} else if (step.action == 'fill'){
			await page.fill(step.target, step.value)
			now = new Date();
			logtime(now, start, 'filled')
		} else if (step.action == 'click'){
			await page.click(step.target)
			now = new Date();
			logtime(now, start, 'clicked')
		} else {
			console.log("Unknown step type!", step)
		}

		if(step.sleep){
			await page.waitForTimeout(step.sleep * videoSpeed);
		}
	}
	// Sleep an extra 1.5s at the end.
	await page.waitForTimeout(1500);
	await browser.close();
	process.stdout.write(JSON.stringify(syncReport));
})();
