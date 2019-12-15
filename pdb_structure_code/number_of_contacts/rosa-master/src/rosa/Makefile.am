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

if INCLUDE_DEBUG
LEX_ARGS = -d
YACC_ARGS = --debug
endif

MAINTAINERCLEANFILES	= Makefile.in

VPATH = $(srcdir)

AM_CPPFLAGS=-I$(top_builddir)/src -I$(top_srcdir)/src 

pkginclude_HEADERS 	= atom.h common.h dbg.h elems.h mathutil.h \
  types.h util.h atomsel_driver.h atomsel_scanner.h pdb.h structure.h \
  error.h version.h geometry.h superpos.h seqalign.h mem_shared_ptr.h

nodist_pkginclude_HEADERS = atomsel_parser.h

lib_LIBRARIES 		= librosa.a
librosa_a_SOURCES = atom.cpp elems.cpp atomsel_driver.cpp structure.cpp \
	version.cpp error.cpp util.cpp geometry.cpp mathutil.cpp superpos.cpp\
	pdb.cpp seqalign.cpp

dist_noinst_DATA  =  atomsel_parser.yy atomsel_scanner.ll
nodist_librosa_a_SOURCES = atomsel_scanner.cpp atomsel_parser.cpp

atomsel_parser.h atomsel_parser.cpp: $(srcdir)/atomsel_parser.yy
	$(YACC) $(YACC_ARGS) --defines=atomsel_parser.h $(srcdir)/atomsel_parser.yy \
  -o atomsel_parser.cpp 
  
atomsel_scanner.cpp: atomsel_scanner.ll atomsel_parser.h
	$(LEX) $(LEX_ARGS) -f -oatomsel_scanner.cpp $(srcdir)/atomsel_scanner.ll

$(srcdir)/atomsel_driver.cpp: atomsel_parser.h

CLEANFILES = atomsel_parser.h atomsel_parser.cpp atomsel_scanner.cpp location.hh \
  position.hh stack.hh

		
