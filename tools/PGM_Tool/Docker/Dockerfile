# SRC: https://raw.githubusercontent.com/joyzoursky/docker-python-chromedriver/master/py3/py3.8-alpine3.10-selenium/Dockerfile
FROM python:3.8-alpine3.10

# update apk repo
RUN echo "https://dl-4.alpinelinux.org/alpine/v3.10/main" >> /etc/apk/repositories && \
    echo "https://dl-4.alpinelinux.org/alpine/v3.10/community" >> /etc/apk/repositories

# install chromedriver
RUN apk update
RUN apk add chromium chromium-chromedriver

# upgrade pip
RUN pip install --upgrade pip

# install selenium
RUN pip install selenium flask