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

VPATH = $(srcdir)

AM_CPPFLAGS=-I$(top_builddir)/src -I$(top_srcdir)/src

bin_PROGRAMS = align bindinterf contacts distances extractmodel index2res \
  interact interpretscontacts mapmodels pdbcut pdbinfo pdbmodify pdbselect \
  rmsd superpos transform  
		
align_SOURCES = \
	align.cpp

bindinterf_SOURCES = \
	bindinterf.cpp

contacts_SOURCES = \
	contacts.cpp

distances_SOURCES = \
	distances.cpp

extractmodel_SOURCES = \
	extractmodel.cpp

index2res_SOURCES = \
	index2res.cpp

interact_SOURCES = \
	interact.cpp

interpretscontacts_SOURCES = \
	interpretscontacts.cpp

mapmodels_SOURCES = \
	mapmodels.cpp

pdbcut_SOURCES = \
	pdbcut.cpp

pdbinfo_SOURCES = \
	pdbinfo.cpp

pdbmodify_SOURCES = \
	pdbmodify.cpp

pdbselect_SOURCES = \
	pdbselect.cpp

rmsd_SOURCES = \
	rmsd.cpp

superpos_SOURCES = \
	superpos.cpp

transform_SOURCES = \
	transform.cpp

LDADD = \
	$(top_builddir)/src/rosa/librosa.a
