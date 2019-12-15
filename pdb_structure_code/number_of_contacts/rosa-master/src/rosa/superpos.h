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

/*! \file superpos.h
 *  \brief Declares functions used to superpose sets of positions.
 */

#ifndef ROSA_SUPERPOS_H_
#define ROSA_SUPERPOS_H_

#include <rosa/geometry.h>
#include <rosa/mathutil.h>

namespace rosa {
  //! Calculates the quaternion superposition for superposing the set v2 on the
  //! set v1. The RMSD of the superposition is returned by the function. The
  //! Transformation is saved in the three variables rotMatrix, cm1 and cm2. In
  //! order to superpose v2 on v1 you need to:
  //! 1) translate of -cm1
  //! 2) apply rotation matrix rotMatrix
  //! 3) translate of cm2
  double quaternionSuperposition( const std::vector<Point> &v1,
                                  const std::vector<Point> &v2,
                                  Mat2D<double> &rotMatrix,
                                  Point &cm1, Point &cm2 );

  //! Calculates the RMSD for the quaternion superposition for superposing the
  //! set v2 on the set v1.
  double getRmsdQuaternionSuperposition( const std::vector<Point> &v1,
                                         const std::vector<Point> &v2 );
}

#endif // ROSA_SUPERPOS_H_
