#!/bin/bash
# usage: sh ./backup-data.sh 
name='bk'
cp -avr RUN/DATA/ RUN/DATA-${name}/
cp -avr RUN/output/ RUN/output-${name}/
cp -avr SRC/output/ SRC/output-${name}/
cp -avr SRC/DATA/ SRC/DATA-${name}/
cp -avr DATA/ DATA-${name}/
