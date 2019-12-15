## Copyright (C) 2008 by Roberto Mosca.
##
## E-mail: info@librosa.org
##
## This file is part of Rosa.
##
## Rosa is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3, or (at your option)
## any later version.
##
## Rosa is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Rosa. If not, see <http://www.gnu.org/licenses/>.

MAINTAINERCLEANFILES	= Makefile.in

AM_CPPFLAGS=-I$(top_builddir)/src -I$(top_srcdir)/src

do_subst = sed -e 's,[@]datadir[@],$(datadir),g' \
	-e 's,[@]BASH[@],$(BASH),g' \
	-e 's,[@]PERL[@],$(PERL),g' \
	-e 's,[@]PACKAGE[@],$(PACKAGE),g' \
	-e 's,[@]VERSION[@],$(VERSION),g'
                 
check_PROGRAMS	= test_atom test_atomsel test_pdbcut test_transform

TESTS = test_atom.sh test_atomsel.sh test_pdbcut.sh test_transform.sh

dist_noinst_DATA = atoms1.pdb atoms2.pdb 1atp.pdb correct_atoms.pdb \
	incorrect_atoms.pdb

check_SCRIPTS = test_atom.sh test_atomsel.sh test_pdbcut.sh test_transform.sh

noinst_HEADERS  = test.h
test_atom_SOURCES	= test_atom.cpp test.cpp
test_atom_LDADD	= $(top_builddir)/src/rosa/librosa.a
test_atomsel_SOURCES	= test_atomsel.cpp test.cpp
test_atomsel_LDADD	= $(top_builddir)/src/rosa/librosa.a
test_pdbcut_SOURCES	= test_pdbcut.cpp test.cpp
test_pdbcut_LDADD	= $(top_builddir)/src/rosa/librosa.a
test_transform_SOURCES	= test_transform.cpp test.cpp
test_transform_LDADD	= $(top_builddir)/src/rosa/librosa.a

CLEANFILES = $(check_SCRIPTS)
EXTRA_DIST = test_atom.sh.in test_atomsel.sh.in test_pdbcut.sh.in \
	test_transform.sh.in

test_atom.sh: $(srcdir)/test_atom.sh.in
	$(do_subst) < $(srcdir)/test_atom.sh.in > test_atom.sh
	chmod +x test_atom.sh

test_atomsel.sh: $(srcdir)/test_atomsel.sh.in
	$(do_subst) < $(srcdir)/test_atomsel.sh.in > test_atomsel.sh
	chmod +x test_atomsel.sh

test_pdbcut.sh: $(srcdir)/test_pdbcut.sh.in
	$(do_subst) < $(srcdir)/test_pdbcut.sh.in > test_pdbcut.sh
	chmod +x test_pdbcut.sh

test_transform.sh: $(srcdir)/test_transform.sh.in
	$(do_subst) < $(srcdir)/test_transform.sh.in > test_transform.sh
	chmod +x test_transform.sh

##	rm -f test_atomsel.sh
##	echo "#! "$(BASH) > test_atomsel.sh
##	cat $(srcdir)/test_atomsel.sh.in >> test_atomsel.sh
##	chmod ugo+x test_atomsel.sh
