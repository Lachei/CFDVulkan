# Import relevant libraries and tools
import os
from flask import (
    Flask,
    redirect,
    url_for,
    render_template,
    send_from_directory,
    request,
    flash,
)

# Define Flask Application and Uploads folders
app = Flask(__name__)
app.secret_key = os.urandom(24)
app.config["UPLOADED_PATH"] = os.path.join(app.root_path, "uploads")

# Check if uploads folder exists
if os.path.isdir('uploads') == False:
    os.mkdir('uploads')

def prepare_array(array,pgm=False,csv=False):
    data = dict()
    if pgm:
        row_col_temp = array.split("\n")[1].split(" ")
        cols = row_col_temp[0]
        rows = row_col_temp[1]
        data_array = [i.strip() for i in array.split("\n")[2:-1]]
        for i in range(len(data_array)):
            data.update({str(i): [j for j in data_array[i].split(" ") if j.isdigit()]})
        return rows, cols, data
    elif csv:
        array_list = array.split("\n")
        for i in range(len(array_list)):
            data.update({str(i): [j for j in array_list[i].split(",")]})
        rows = len(data.keys())
        cols = len(data['0'])
        return rows, cols, data
    return None

@app.route("/", methods=["GET", "POST"])
def start():
    if request.method == "POST":
        for f in request.files.getlist("upload"):
            filetype = f.filename.split('.')[-1]
            array = f.read().decode("utf-8")
            if filetype == "pgm":
                rows, cols, data = prepare_array(array,pgm=True)
            elif filetype == "csv":
                rows, cols, data = prepare_array(array,csv=True)
            else:
                return render_template('index_flask.html')
        return render_template(
            "index_flask.html", array=data, rows=rows, cols=cols, filename=f.filename, filetype=filetype
        )
    return render_template("index_flask.html")


if __name__ == "__main__":
    app.run(debug=True, port=8000)
