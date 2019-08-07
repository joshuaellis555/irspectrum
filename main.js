'use strict';

const express = require('express'),
	app = express(),
	bodyParser = require('body-parser'),
	formidable = require('formidable'),
	path = require('path'),
	fs = require('fs'),
	spawn = require('child_process').spawn;

app.use(bodyParser.json());
app.use(bodyParser.urlencoded({
	extended: true,
}));

app.use('/favicon.ico', express.static('./public/favicon.ico'));

// Setup our static resources folder
app.use(express.static('public'));

// Setup Pug
app.set('view engine', 'pug');
app.set('views', './views');

// Setup views
app.get('/', function (req, res) {
	res.render('home');
});

app.get('/aboutus', function (req, res) {
	res.render('aboutus');
});

app.post('/getmatch', function (req, res) {
	// set up recieving file
	let form = new formidable.IncomingForm();
	form.multiples = false;
	form.uploadDir = path.join(__dirname, '/temp');
	let filepath,
		oldFilePath,
		filename,
		timestamp;

	// parse the incoming file
	form.parse(req, (err, fields, files) => {
		timestamp = fields.timestamp;
	});

	// save the file as its name in the /temp folder
	form.on('file', (field, file) => {
		oldFilePath = file.path;
		filename = file.name;
	});

	// if something wonky happens
	form.on('error', function (err) {
		console.log('An error has occured: \n' + err);
	});

	// once its finished saving, do the business. This is where the magic comparing happens
	form.on('end', () => {
		// path to the file including the timestamp. etc: c:/something/somethingelse/90-40-3_a_123456789.pdf
		filepath = path.join(form.uploadDir, filename.substring(0, filename.length - 4) + '_' + timestamp + filename.substring(filename.length - 4));
		fs.rename(oldFilePath, filepath, err => {
			if (err) {
				console.log(err);
			}
			runJoshsPython(filepath, res, filename.substring(0, filename.length - 4) + '_' + timestamp + '.jpg');
		});
	});
});

// Start the server and have it listen on a port
const server = app.listen(3000, function () {
	console.log(`Server started on port ${server.address().port}!`);
});

// this function runs the python programs
let runJoshsPython = function (pathToPDF, res, filename) {
	// run 'Compare to Query.py'
	const comparePy = spawn('python', ['Query.py', pathToPDF, filename]);

	// if something wonky happens let me know
	comparePy.on('error', err => {
		console.log('error:', err);
	});

	// if the program has an exception let me know
	comparePy.on('uncaughtException', function (err) {
		console.log('Caught exception: ' + err);
	});

	// once the program gives any output
	comparePy.stdout.on('data', chunk => {
		// format output
		let textchunk = chunk.toString('utf8');

		// send results back to browser
		returnResults(textchunk, res);
	});
};

// this function returns the results back to the browser
let returnResults = function (winners, res) {
	// split string of winners into array
	winners = winners.split(' ');

	// send back JSON of top 20 results
	res.send(JSON.stringify({
		1: winners[0],
		2: winners[1],
		3: winners[2],
		4: winners[3],
		5: winners[4],
		6: winners[5],
		7: winners[6],
		8: winners[7],
		9: winners[8],
		10: winners[9],
		11: winners[10],
		12: winners[11],
		13: winners[12],
		14: winners[13],
		15: winners[14],
		16: winners[15],
		17: winners[16],
		18: winners[17],
		19: winners[18],
		20: winners[19],
	}));
};
