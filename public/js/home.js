'use strict';

// html css bootstrap updating
let tab = document.getElementById('home');
tab.classList.add('active');

// grabs a hold of the input for use later
let input = document.getElementById('fileInput');

// allows the input button to look pretty
document.getElementById('selectFile').addEventListener('click', function () {
	input.click();
});

// displays the image when you click select file
input.onchange = function () {
	// gets a hold of the image
	let imageInput = input.files[0];

	let ext = imageInput.name.substring(imageInput.name.length - 3);

	if (ext == 'pdf' || ext == 'jpg') {
		// create html stuff
		document.getElementById('filename').innerText = input.files[0].name;
	} else {
		document.getElementById('filename').innerText = 'Invalid File';
	}

	// enable the find button once file is uploaded
	document.getElementById('findButton').classList = 'btn btn-primary';
};

// This is where we send the picture to the server and the server will do the comparison
// and reply with the matches somehow.
document.getElementById('findButton').addEventListener('click', function () {
	$(' #bodyholder ').html('');

	let imageInput = input.files[0];
	let ext = imageInput.name.substring(imageInput.name.length - 3);
	if (ext == 'pdf' || ext == 'jpg') {
		// show spinny wheel until finished
		document.body.style.cursor = 'wait';

		// get the timestamp for the file
		let timeStampInMs = window.performance && window.performance.now && window.performance.timing && window.performance.timing.navigationStart ? window.performance.now() + window.performance.timing.navigationStart : Date.now();
		timeStampInMs = timeStampInMs.toString();

		// prepare file to be sent
		let myData = new FormData();
		myData.append('timestamp', timeStampInMs);
		myData.append('file', imageInput);

		// send file to server
		$.ajax({
			url: 'http://' + location.host + '/getmatch',
			data: myData,
			processData: false,
			contentType: false,
			type: 'POST',
			complete: function (data) {
				let response = JSON.parse(data.responseText);

				// create all the html stuff we need
				let winnersContainer = document.createElement('span');
				let winnersList = document.createElement('ul');
				winnersList.classList.add('list-group');

				let winnerULS = [];

				// create the winners nav list on left side
				for (let i = 0; i < 20; i++) {
					let li = document.createElement('li');
					li.classList.add('list-group-item');
					li.innerText = response[i + 1];
					li.addEventListener('click', updateResults);
					winnerULS.push(li);
				}

				// add those to the list
				for (let i = 0; i < 20; i++) {
					winnersList.appendChild(winnerULS[i]);
				}

				// add the list to the screen
				winnersContainer.appendChild(winnersList);
				winnersContainer.classList.add('floatleft');
				$(' #bodyholder ').append(winnersContainer);

				// add the graph image of the uploaded pic
				let uploaded = document.createElement('img');
				uploaded.src = '/uploads/' + imageInput.name.substring(0, imageInput.name.length - 4) + '_' + timeStampInMs + '.jpg';
				uploaded.classList.add('resize');
				$(' #bodyholder ').append(uploaded);

				// create holder of results
				let resultsSpan = document.createElement('span');
				resultsSpan.classList.add('card');

				// create the results picture
				let resultsImg = document.createElement('img');
				resultsImg.id = 'resultsImage';
				resultsImg.classList.add('resize');
				resultsImg.classList.add('card-img-top');
				resultsImg.src = '/images/' + response[1] + '.jpg';

				// create a holder for the results info
				let info = document.createElement('div');
				info.classList.add('card-body');

				// create the picture of the winner molecule
				let molecule = document.createElement('img');
				molecule.id = 'molecule';
				molecule.classList.add('smaller-resize');
				molecule.src = '/info/' + response[1] + '.png';
				info.appendChild(molecule);

				// create the name of the result
				let name = document.createElement('h5');
				name.innerText = response[1];
				name.id = 'name';
				info.appendChild(name);

				resultsSpan.appendChild(resultsImg);
				resultsSpan.appendChild(info);
				$(' #bodyholder ').append(resultsSpan);

				// get other info about the results
				$.ajax({
					url: 'http://' + location.host + '/info/' + response[1] + '.json',
					processData: false,
					contentType: false,
					type: 'GET',
					complete: function (json) {
						let moleculeinfo = json.responseJSON;

						// create html stuff
						let specID = document.createElement('p');
						let cas = document.createElement('p');
						let formula = document.createElement('p');
						let compound = document.createElement('p');

						// give the html stuff their words
						specID.innerText = 'Spectrum ID: ' + moleculeinfo.spectrumID;
						cas.innerText = 'CAS: ' + moleculeinfo.cas;
						formula.innerHTML = '<p>Formula: ' + getFormulaHTML(moleculeinfo.formula) + '</p>';
						compound.innerText = 'Name: ' + moleculeinfo.name;

						// give the html stuff id's
						specID.id = 'specid';
						cas.id = 'cas';
						formula.id = 'formula';
						compound.id = 'compound';

						// display the html stuff
						info.appendChild(specID);
						info.appendChild(cas);
						info.appendChild(formula);
						info.appendChild(compound);
					}
				});

				// set cursor back to normal
				document.body.style.cursor = 'default';
			}
		});
	}
});

let updateResults = evt => {
	$.ajax({
		url: 'http://' + location.host + '/info/' + evt.target.innerText + '.json',
		processData: false,
		contentType: false,
		type: 'GET',
		complete: function (json) {
			let moleculeinfo = json.responseJSON;

			$(' #resultsImage ').attr('src', '/images/' + evt.target.innerText + '.jpg');
			$(' #name ').text(evt.target.innerText);
			$(' #specid ').text('Spectrum ID: ' + moleculeinfo.spectrumID);
			$(' #cas ').text('CAS: ' + moleculeinfo.cas);
			$(' #formula ').html('<p>Formula: ' + getFormulaHTML(moleculeinfo.formula) + '</p>');
			$(' #compound ').text('Name: ' + moleculeinfo.name);
			$(' #molecule ').attr('src', '/info/' + evt.target.innerText + '.png');


		}
	});
};

let getFormulaHTML = formula => {
	let nums = false;
	let returnstring = '';

	for (let i = 0; i < formula.length; i++) {
		if (isNaN(formula.charAt(i))) {
			if (nums) {
				nums = false;
				returnstring += '</sub>' + formula.charAt(i);
			} else {
				returnstring += formula.charAt(i);
			}
		} else {
			if (!nums) {
				nums = true;
				returnstring += '<sub>' + formula.charAt(i);
			} else {
				returnstring += formula.charAt(i);
			}
		}
	}

	if (nums) {
		returnstring += '</sub>';
	}

	return returnstring;
};