# Build Documentaion with Gitlab CI

One advantage of using a documentation like `mkdocs` is that you can also create the docs in the Gitlab CI pipeline. The following snippets provides you with an example. Important to mention is the Gitlab Page functionality requires an html output with the folder `public`. Therefore the `mkdocs` output is directed to a folder called `public`. This is then collected as artifact of the pipeline which then can be stored. 

```yaml
pages:
    image: python:alpine
    before_script:
      - pip install mkdocs
      - pip install mkdocs-material
    script:
        - mkdocs build -d public
    artifacts:
        paths:
            - public
    only:
      changes:
        - docs/**/*
```

