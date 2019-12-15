/***************************************************************************
 *                                                                         *
 *   Copyright (C) 2008 by Roberto Mosca.                                  *
 *                                                                         *
 *   E-mail: info@librosa.org                                              *
 *                                                                         *
 *   This file is part of Rosa.                                            *
 *                                                                         *
 *   Rosa is free software: you can redistribute it and/or modify          *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3, or (at your option)   *
 *   any later version.                                                    *
 *                                                                         *
 *   Rosa is distributed in the hope that it will be useful,               *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Rosa. If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                         *
 ***************************************************************************/

/*! \file memsharedptr.h
 *  \brief Compatibility stuff for the shared_ptr template class.
 */

#ifndef ROSA_MEM_SHARED_PTR_H_
#define ROSA_MEM_SHARED_PTR_H_

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef HAVE_LIBBOOST
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;
#else
#include <tr1/memory>
using std::tr1::shared_ptr;
#endif


#endif // ROSA_MEM_SHARED_PTR_H_
