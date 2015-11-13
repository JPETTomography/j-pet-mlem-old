#!/bin/bash
if [ $# -ne 1 ]; then
	echo >&2 "usage: $(basename $0) host"
fi

cd >/dev/null 2>&1 $(dirname $0)/..

rsync --delete -avc src/ "$1":Projects/PET/src/
