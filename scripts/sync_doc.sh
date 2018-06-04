#!/bin/bash

rsync -avu --delete doc/html/ sorbus.if.uj.edu.pl:Sites/ && \
ssh sorbus.if.uj.edu.pl 'sed -i \"\" -e s,src=\"http://,src=\"//,g Sites/*.html'
