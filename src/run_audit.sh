#!/bin/bash

# to avoid mucking about with escaping args to java

java -Dr_flags="--no-save --quiet --slave --no-restore"  -DR_HOME=/usr/local/lib64/R -cp /build RunR /build/audit/AuDIT.R parse.cmdline /build/audit/ "$@"
java_exit=$?
echo "AuDIT finished" 
exit   $java_exit

