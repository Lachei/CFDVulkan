# Documentation CFD Lab

The documentation is based on [mkdocs](https://www.mkdocs.org/). This should give you the possibility to write your own documentation and have an overview of the different tools that come along with the skeleton. Finally, you can have a look into the setup section which gives you information for the software setup of the CFD Lab. 

!!! warning
    This documentation should be seen as an extension to the worksheets and doesn't replace the worksheets.

## Tools 

Please see the provided Tools as **experimental**. 

## Folder Structure

The `docs` folder contains all the necessary Markdown Files which then are parsed by `mkdocs`. The `mkdocs.yml` file defines the main properties of the documentation.

## Theme  

The theme can be easily changed via the `theme` keyword. Currently the theme is selected to `materialize`, a CSS framework by google. One can also change other themes, an overview can be found with the following [link](https://github.com/mkdocs/mkdocs/wiki/MkDocs-Themes). Tipp: Make sure that you have installed the desired theme e.g. via `pip3 install mkdocs-material` otherwise you run will likely run into problems. 

## Gitlab Config

For serving the documentation via the [Gitlab LRZ Page](https://doku.lrz.de/display/PUBLIC/GitLab)-Feature, a CI pipeline is required. Furthermore the default folder which is used by gitlab for the `html`-files is the `public`. The setting in the `mkdocs.yml` is called `site_dir`.

## Mkdocs Instruction

The following commands should be executed in the top level of this repository.

- `mkdocs serve`: builds the documentation and serves it under [localhost:8000](localhost:8000)
- `mkdocs build`: create the whole documentation as html folder

## MkDocs Configs 

Code highlighting is enabled for the `materialize`-theme via the following setting:

```bash
markdown_extensions:
    - codehilite
```

in the `readthedocs`-theme one can simply extend the `theme`-settings:

```bash
theme:
    name: 'readthedocs'
    highlightjs: true
    hljs_languages:
        - cpp
```

## Other Documentation Frameworks


Tool Name | Readthedocs - Publication | Github & Gitlab Pages | Comment
---------|----------|---------|---------
[docusaurus](https://docusaurus.io/en/) | ❌ | ✅ | Development by Facebook, based on react, more complex folder structure
 [mkdocs](https://www.mkdocs.org/) | ✅ (limited to certain themes) | ✅ | Python based, installation via `pip`, simple folder structured
 [RBookdown](https://bookdown.org/yihui/bookdown/) | ❌ | ✅ | R based, documentation is written in Rmakdown books, typically written with Rstudio
 [Sphinx](https://www.sphinx-doc.org/en/master/) | ✅ | ✅ | Python based, documentation is written in `rst`-files, writing the rst-files requires some introduction

