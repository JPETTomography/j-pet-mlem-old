#!/bin/bash

exec rsync -avu --delete doc/html/ sorbus.if.uj.edu.pl:Sites/
