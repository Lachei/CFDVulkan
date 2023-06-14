#!/usr/bin/env bash

# Get line number of starting point 
# cat -n prints the file with line numbers
# grep searchs for the respective start
# awk prints out the filenumber
start=$(cat -n ../templates/index_flask.html | grep 'Flask Template Start' | awk '{print $1}')
end=$(cat -n ../templates/index_flask.html | grep 'Flask Template End' | awk '{print $1}')

# sed deletes the line numbers from start till end
# the output is stored in the upper directory with the filename index.html
sed "${start},${end}d" ../templates/index_flask.html > ../index.html
