echo off
pushd %~dp0
title IRSpectrum Init
echo Initializing...
pip install -r requirements.txt
npm install
python UpdateDB.py
node main.js
