### copyright 2017-2021 Regents of the University of California and the Broad Institute. All rights reserved.
FROM genepattern/docker-r-2-5:0.1

MAINTAINER Ted Liefeld <liefld@broadinstitute.org>

# creating CMS_HOME & copying src files there - the command line then puts the jars on the classpath
# the R file in src is needed to compute q-values and is included in the module zip so that it is in taskLib
RUN mkdir -p /build/audit
COPY src /build/audit/

# docker build --rm https://github.com/genepattern/ComparativeMarkerSelection.git#develop -f Dockerfile -t genepattern/comparativemarkerselection:v11
# make sure this repo and tag match the manifest & don't forget to docker push!
# docker push genepattern/genepattern/comparativemarkerselection:<tag>

# you can use this command to run Docker and iterate locally (update for your paths and module name, of course)
# docker run --rm -it -v /c/Users/myUSER/PathTo/ExampleModule/ComparativeMarkerSelection:/mnt/mydata:rw genepattern/comparativemarkerselection:v11 bash

#Note that this container is running as root - not ideal, but I'm leaving this for now so I don't have to refactor any Java

