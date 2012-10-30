-include ../lib/common.mk

# libpng support & detection
LIBPNG := $(strip $(firstword \
		$(wildcard /usr/include/libpng) \
		$(wildcard /usr/local/include/libpng) \
		$(wildcard /opt/X11/include/libpng15) ) )
ifneq ($(LIBPNG),)
	CPPFLAGS += -I$(LIBPNG) -DHAVE_LIBPNG
	LDFLAGS  += -L$(dir $(patsubst %/,%,$(dir $(LIBPNG))))lib -lpng
endif
