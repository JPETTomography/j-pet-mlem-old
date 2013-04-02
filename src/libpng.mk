# libpng support

ifneq ($(LIBPNG),)
	# defined by modules
	LIBPNG_INCLUDE := $(LIBPNG)/include
	LIBPNG_LIB     := $(LIBPNG)/lib
else
	# look in well known places
	LIBPNG_INCLUDE := $(strip $(firstword \
		$(wildcard /usr/include/libpng) \
		$(wildcard /usr/local/include/libpng) \
		$(wildcard /opt/X11/include/libpng15) ) )
	ifneq ($(LIBPNG_INCLUDE),)
		LIBPNG_LIB := $(dir $(patsubst %/,%,$(dir $(LIBPNG))))lib
	endif
endif

ifneq ($(LIBPNG_INCLUDE),)
	CPPFLAGS += -I$(LIBPNG_INCLUDE) -DHAVE_LIBPNG
	LDFLAGS  += -L$(LIBPNG_LIB) -lpng
endif
