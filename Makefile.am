SUBDIRS = src

EXTRA_DIST = LICENSE \
             autogen.sh \
             README.md

DISTCLEANFILES = Makefile.in \
                 install-sh \
                 aclocal.m4 \
                 configure \
                 missing

distclean-local:
	rm -rf autom4te.cache \
	rm -f compile \
	rm -f depcomp \
	rm -rf src/.deps/ \
	rm -f config.guess{?,} \
	rm -f config.sub{?,}
