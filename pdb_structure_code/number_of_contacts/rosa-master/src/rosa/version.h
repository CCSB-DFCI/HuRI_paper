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

/*! \file version.h
 *  \brief Declares functions used to print and refer to the version number
 *         of the library.
 */

#ifndef ROSA_VERSION_H_
#define ROSA_VERSION_H_

#include <string>
 
namespace rosa {
  //! Prints the version number of the rosa library. If an optional program name
  //! is passed as the argument, then the message is customized to include
  //! the program name. This is done in the assumptions that programs
  //! distributed with rosa are released with the same version number.
  void printRosaVersion( const std::string &programName = "" );
}

#endif // ROSA_VERSION_H_
