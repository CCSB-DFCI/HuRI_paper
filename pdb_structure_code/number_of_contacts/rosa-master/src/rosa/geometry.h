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

/*! \file geometry.h
 *  \brief Declares data types and functions for handling geometrica entities
 *         like points, vectors, etc...
 */

#ifndef ROSA_GEOMETRY_H_
#define ROSA_GEOMETRY_H_

#include <rosa/types.h>
#include <iosfwd>
#include <vector>
#include <cmath>
 
namespace rosa {
  
  struct Point {
    Coord pos[3];
   
    Point( Coord aX, Coord aY, Coord aZ ) {
      pos[0] = aX; pos[1] = aY; pos[2] = aZ;
    }
    
    //! Initializes to a point in the origin of the cartesian system (0, 0, 0)
    Point() {
      pos[0] = 0.0; pos[1] = 0.0; pos[2] = 0.0;
    }
    
    Coord x() const { return pos[0]; }
    Coord y() const { return pos[1]; }
    Coord z() const { return pos[2]; }
        
    void x( Coord aX ) { pos[0] = aX; }
    void y( Coord aY ) { pos[1] = aY; }
    void z( Coord aZ ) { pos[2] = aZ; }
    
    Point &operator += ( const Point &p ) {
      pos[0] += p.pos[0];
      pos[1] += p.pos[1];
      pos[2] += p.pos[2];
      return (*this);
    }

    Point &operator -= ( const Point &p ) {
      pos[0] -= p.pos[0];
      pos[1] -= p.pos[1];
      pos[2] -= p.pos[2];
      return (*this);
    }

    Point &operator += ( Coord s ) {
      pos[0] += s;
      pos[1] += s;
      pos[2] += s;
      return (*this);
    }
    
    Point &operator /= ( double s ) {
      pos[0] /= s;
      pos[1] /= s;
      pos[2] /= s;
      return (*this);
    }

    double sqNorm() const {
      return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
    }

    double norm() const {
      return sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
    }

    Point unit() const {
      double pNorm = norm();
      Point r( pos[0] / pNorm, pos[1] / pNorm, pos[2] / pNorm );
      return r;
    }
  };
  
  inline Point operator - ( const Point &p ) {
    Point r( -p.x(), -p.y(), -p.z() );
    return r;
  }

  inline Point operator * ( const Point &p, double weight ) {
    Point r( p.x()*weight, p.y()*weight, p.z()*weight );
    return r;
  }

  inline Point operator * ( double weight, const Point &p ) {
    Point r( weight*p.x(), weight*p.y(), weight*p.z() );
    return r;
  }

  inline double operator * ( const Point &p1, const Point &p2 ) {
    return p1.x()*p2.x() + p1.y()*p2.y() + p1.z()*p2.z();
  }

  inline Point operator / ( const Point &p, double weight ) {
    Point r( p.x() / weight, p.y() / weight, p.z() / weight );
    return r;
  }

  inline Point operator + ( const Point &p1, const Point &p2 ) {
    Point r( p1.x()+p2.x(), p1.y()+p2.y(), p1.z()+p2.z() );
    return r;
  }

  inline Point operator - ( const Point &p1, const Point &p2 ) {
    Point r( p1.x()-p2.x(), p1.y()-p2.y(), p1.z()-p2.z() );
    return r;
  }

  std::ostream &operator << ( std::ostream &, const Point & );
  
  //! Returns the euclidean distance between two points
  inline double pointDistance( const Point &p, const Point &r )
  {
    double x = p.x() - r.x();
    double y = p.y() - r.y();
    double z = p.z() - r.z();
    return sqrt( x*x + y*y + z*z );
  }
  
  //! Returns the squared euclidean distance between two points
  inline double pointSqrDistance( const Point &p, const Point &r )
  {
    double x = p.x() - r.x();
    double y = p.y() - r.y();
    double z = p.z() - r.z();
    return ( x*x + y*y + z*z );
  }
  
  void getBaricenter( const std::vector<Point> &v, Point &b );

  double rmsd( const std::vector<Point> &v1, const std::vector<Point> &v2 );

} // namespace rosa

#endif // ROSA_GEOMETRY_H_
