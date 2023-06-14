# Introduction

The pgm tool is written in Python with the framework [flask](https://flask.palletsprojects.com/). You have two option to run this fool:

1. run it *with* flask
2. run it *without* flask

## Running the tool *with* flask

### Prerequisites

- python3 Installation
- pip3 (installation via: `apt-get install -y python3-pip`)
- flask (installation with pip: `pip3 install --user flask`)

To get all the functionalities of the provided application you should be running it with flask. You can start the flask server with the following commands:

```shell
export FLASK_APP=app.py
flask run
```

## Running the tool *without* flask

You can also use the plain `index.html`-file which is provided with this folder as well. This is a reduced version of the `index.html` which is located in the `template`-folder. However in this case you can't use the uploading / restoring function. The technical background is as follows: In order to handle an upload to a website, a webserver is required, which handles the `POST` request of the browser. The latter request can't be handled if you just open the `index.html` in your browser.

![](../../img/Flask_App.svg)