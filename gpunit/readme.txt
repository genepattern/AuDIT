Gp-Unit Test cases for the AuDIT module. 

Manually created yaml test cases, based on exported jobs. The jobs were run initially with AuDIT 5.2 installed 
on a GP 3.7.0 server running on the gdev CentOS VM at the Broad.

To run the gpunit tests from this directory
    ant -f ../../util/gp-unit/build.xml -Dgpunit.dir=../../util/gp-unit -Dbasedir=`pwd` -Dgpunit.properties=`pwd`/gpunit.AuDIT.properties gpunit 
To clean up
     delete the System.out file (It's a mystery to me why this shows up)
     delete the build directory and all of it's contents

This works best (it's only been tested) after you have run gpunit from the gpunit.dir, e.g.
    cd ../../util/gp-unit
    ant 
