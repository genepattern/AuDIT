FROM 	genepattern/docker-r-2-5:0.1
ENV R_HOME=/usr/local/lib64/R
RUN mkdir -p /build/audit

COPY src/run_audit.sh /build/audit/
RUN chmod a+x /build/audit/run_audit.sh
COPY src/*.R /build/audit/
