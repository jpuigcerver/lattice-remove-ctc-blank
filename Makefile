.PHONY: all install

all:
ifndef KALDI_ROOT
  $(error KALDI_ROOT environment variable is undefined)
endif

EXTRA_CXXFLAGS = -Wno-sign-compare -Wno-unused-variable -I$(KALDI_ROOT)/src
include $(KALDI_ROOT)/src/kaldi.mk

BINFILES = lattice-remove-ctc-blank

OBJFILES =

TESTFILES =

ADDLIBS = \
	$(KALDI_ROOT)/src/lat/kaldi-lat.a \
	$(KALDI_ROOT)/src/fstext/kaldi-fstext.a \
        $(KALDI_ROOT)/src/hmm/kaldi-hmm.a \
        $(KALDI_ROOT)/src/util/kaldi-util.a \
	$(KALDI_ROOT)/src/matrix/kaldi-matrix.a \
	$(KALDI_ROOT)/src/thread/kaldi-thread.a \
	$(KALDI_ROOT)/src/base/kaldi-base.a

include $(KALDI_ROOT)/src/makefiles/default_rules.mk

PREFIX=/usr/local
install: lattice-remove-ctc-blank
	test -d $(PREFIX) || mkdir -p $(PREFIX)/bin
	install -m 0755 lattice-remove-ctc-blank $(PREFIX)/bin
