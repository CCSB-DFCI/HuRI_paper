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

#include <rosa/geometry.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

using namespace std;
using namespace rosa;

ostream &rosa::operator << ( ostream &os, const Point &p )
{
  os << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
  return os;
}

void rosa::getBaricenter( const vector<Point> &v, Point &b )
{
  typedef vector<Point>::const_iterator VPCit;
  
  VPCit vEnd = v.end();
  
  b = Point();
  for( VPCit it = v.begin(); it != vEnd; ++it )
    b += (*it);
  
  b /= (Coord)v.size();
}


double rosa::rmsd( const vector<Point> &v1, const vector<Point> &v2 )
{
  // Checks the two vectors have the same number of elements
  if( v1.size() != v2.size() ) {
    ostringstream s;
    s << "rmsd: sets of points with different sizes (" << v1.size()
      << "!=" << v2.size() << ")" << endl;
    throw logic_error( s.str() );
  }

  unsigned long vSize = v1.size();

  if( vSize == 0 ) return 0.0;
  
  double rmsd = 0.0;
  
  for( unsigned long i = 0UL; i < vSize; ++i )
    rmsd += pointSqrDistance( v1[i], v2[i] );
  
  return sqrt( rmsd / (double)vSize );
}
