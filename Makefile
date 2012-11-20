# run make in all subdirectories
# with Makefile inside
subdirs := $(patsubst %/,%,$(dir $(wildcard **/Makefile)))
# use terminal colors if in terminal
bold := \033[1m
none := \033[0m

all:
%::
	@for dir in $(subdirs); do \
		[ -t 1 ] && printf "$(bold)## %s$(none)\n" "$$dir" || echo "## $$dir"; \
		make -C $$dir $*; \
	done
