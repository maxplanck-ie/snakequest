#!/bin/bash

#Set up the mounting
python /usr/local/bin/mounts.py

#R --no-save --no-restore

R -e 'shiny::runApp("/root/snakequest",port=2525,host="0.0.0.0")' #> /data/prcessing/shiny/snakequest.Rout

