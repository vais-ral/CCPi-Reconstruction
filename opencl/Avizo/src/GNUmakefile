#! gmake

include ../GNUmakefile.setroot

TARGET_DIRS = $(shell /bin/sh $(AVIZO_ROOT)/make/findtargets.sh \
			PACKAGE LIBS $(wildcard */GNUmakefile */*/GNUmakefile */*/*/GNUmakefile */*/*/*/GNUmakefile))

include $(AVIZO_ROOT)/make/multiple.rules
