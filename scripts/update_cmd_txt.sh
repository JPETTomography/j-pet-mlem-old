#!/bin/bash
cd $(dirname $0)/..

while read cmd; do
	cmd=${cmd:2}
	src=${cmd/_//}
	src=${src/_//}
	[ "$src" == "$cmd" ] && continue
	src=src/${src}_cmd.txt
	./$cmd 2> $src
done < <(find . -maxdepth 1 -type f -perm +0111)
