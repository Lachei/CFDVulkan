# Tool Description


## `convert_flask_template.sh`
This bash script converts the flask template to an index file, which can be used by a simple webserver. You have to run the file by `./convert_flask_template.sh` or `bash convert_flask_template.sh` in the misc folder itself, otherwise the directory levels do not fit for the converstion.

⚠️ Due to modularisation of the javascript files you can't open the file just by clicking on index.html. Just run in the directory of the index-html file:

```python
python3 -m http.server
```

Visit [localhost:8000](localhost:8000)