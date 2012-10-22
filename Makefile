# run make in all subdirectories
# with Makefile inside
subdirs := $(patsubst %/,%,$(dir $(wildcard **/Makefile)))
# use terminal colors if in terminal
bold := \033[1m
none := \033[0m

all:
%::
	@for dir in $(subdirs); do \
		if [ -t 1 ]; \
			then echo "$(bold)## $$dir$(none)"; \
			else echo "## $$dir"; \
		fi; \
		make -C $$dir $*; \
	done
