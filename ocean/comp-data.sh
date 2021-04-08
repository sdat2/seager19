#!/bin/bash
# usage: bash ./comp-data.sh 
name='bk'
diff -rq RUN/DATA/ RUN/DATA-${name}/
diff -rq RUN/output/ RUN/output-${name}/
diff -rq SRC/output/ SRC/output-${name}/
diff -rq SRC/DATA/ SRC/DATA-${name}/
diff -rq DATA/ DATA-${name}/
