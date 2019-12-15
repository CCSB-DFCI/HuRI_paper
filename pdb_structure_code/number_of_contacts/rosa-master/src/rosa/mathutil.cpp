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

#include <rosa/mathutil.h>
#include <cmath>

using namespace std;
using namespace rosa;

Mat2D<Coord> rosa::rotMatrix( double xAng, double yAng, double zAng,
                              RotRefSystem rotSystem )
{
  Mat2D<Coord> rotMat( 3, 3 );
  
  double cz = cos( zAng );
  double sz = sin( zAng );
  double cy = cos( yAng );
  double sy = sin( yAng );
  double cx = cos( xAng );
  double sx = sin( xAng );
  
  switch( rotSystem ) {
    case ROT_XYZ:
      rotMat[0][0] = cy * cz;
      rotMat[0][1] = sx * sy * cz - cx * sz;
      rotMat[0][2] = cx * sy * cz + sx * sz;
      
      rotMat[1][0] = cy * sz;
      rotMat[1][1] = sx * sy * sz + cx * cz;
      rotMat[1][2] = cx * sy * sz - sx * cz;
      
      rotMat[2][0] = -sy;
      rotMat[2][1] = sx * cy;
      rotMat[2][2] = cx * cy;
      break;
    case ROT_ZXZ:
      rotMat[0][0] = cz * cx - cy * sx * sz;
      rotMat[0][1] = cz * sx + cy * cx * sz;
      rotMat[0][2] = sz * sy;
      
      rotMat[1][0] = - sz * cx - cy * sx * cz;
      rotMat[1][1] = - sz * sx + cy * cx * cz;
      rotMat[1][2] = cz * sy;
      
      rotMat[2][0] = sy * sx;
      rotMat[2][1] = - sy * cx;
      rotMat[2][2] = cy;
      break;
  }
  
  return rotMat;
}


/** this conversion uses conventions as described on page:
*   http://www.euclideanspace.com/maths/geometry/rotations/euler/index.htm
*   Coordinate System: right hand
*   Positive angle: right hand
*   Order of euler angles: heading first, then attitude, then bank
*   matrix row column ordering:
*   [m00 m01 m02]
*   [m10 m11 m12]
*   [m20 m21 m22]*/
void rosa::rotAngles( const Mat2D<double> &m,
                      double &xAng, double &yAng, double &zAng )
{
  // Assuming the angles are in radians.
	if (m[1][0] > 0.998) { // singularity at north pole
		xAng = atan2(m[0][2],m[2][2]);
		yAng = M_PI/2.0;
		zAng = 0.0;
		return;
	}
	if (m[1][0] < -0.998) { // singularity at south pole
		xAng = atan2(m[0][2],m[2][2]);
		yAng = -M_PI/2.0;
		zAng = 0.0;
		return;
	}
	xAng = atan2(-m[2][0],m[0][0]);
	yAng = atan2(-m[1][2],m[1][1]);
	zAng = asin(m[1][0]);
}
void rosa::eigenValVectSort( vector<double> &d, Mat2D<double> &v)
{
  int i, j, k;
	double p;
 
	int n = d.size();
	for( i = 0; i < n-1; i++) {
		p = d[k=i];
    for( j = i; j < n; j++ )
			if( d[j] >= p )
        p = d[k=j];
    if (k != i) {
      d[k] = d[i];
			d[i] = p;
			for( j = 0; j < n; j++ ) {
        p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}

namespace rosa {
  inline void rot( Mat2D<double> &a, double s, double tau,
                         int i, int j, int k, int l )
  {
    double g, h;

    g = a[i][j];
    h = a[k][l];
    a[i][j] = g - s*( h + g*tau );
    a[k][l] = h + s*( g - h*tau );
  }
}

bool rosa::jacobiTransformation( Mat2D<double> &a, vector<double> &d,
                                 Mat2D<double> &v, int &nrot )
{
	int i, j, ip, iq;
	double tresh, theta, tau, t, sm, s, h, g, c;

	int n = d.size();
	vector<double> b(n), z(n);
  
	for( ip = 0; ip < n; ip++ ) {
		for( iq = 0; iq < n; iq++ ) v[ip][iq] = 0.0;
		v[ip][ip] = 1.0;
	}
	for(ip = 0; ip < n; ip++ ) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	nrot = 0;
	for( i = 1; i <= 50; i++ ) {
		sm = 0.0;
		for( ip = 0; ip < n-1; ip++ ) {
			for( iq = ip+1; iq < n; iq++)
				sm +=  fabs(a[ip][iq]);
		}
		if(sm  ==  0.0)
			return true;
		if(i < 4)
			tresh = 0.2*sm/(n*n);
		else
			tresh = 0.0;
		for( ip = 0; ip < n-1; ip++ ) {
			for( iq = ip+1; iq < n; iq++ ) {
				g = 100.0*fabs(a[ip][iq]);
				if(i > 4 && (fabs(d[ip])+g)  ==  fabs(d[ip])
					&& (fabs(d[iq])+g)  ==  fabs(d[iq]))
						a[ip][iq] = 0.0;
				else if(fabs(a[ip][iq]) > tresh) {
					h = d[iq]-d[ip];
					if((fabs(h)+g)  ==  fabs(h))
						t = (a[ip][iq])/h;
					else {
						theta = 0.5*h/(a[ip][iq]);
						t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if(theta < 0.0) t  =  -t;
					}
					c = 1.0/sqrt(1+t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*a[ip][iq];
					z[ip] -=  h;
					z[iq] +=  h;
					d[ip] -=  h;
					d[iq] +=  h;
					a[ip][iq] = 0.0;
					for( j = 0; j < ip; j++ )
						rot(a, s, tau, j, ip, j, iq);
					for( j = ip+1; j < iq; j++ )
						rot(a, s, tau, ip, j, j, iq);
					for( j = iq+1; j < n; j++ )
						rot(a, s, tau, ip, j, iq, j);
					for( j = 0; j < n; j++ )
						rot(v, s, tau, j, ip, j, iq);
					++nrot;
				}
			}
		}
		for(ip = 0;ip<n;ip++) {
			b[ip] +=  z[ip];
			d[ip]  = b[ip];
			z[ip]  = 0.0;
		}
	}
  return false;
}
