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

#include <rosa/superpos.h>
#include <rosa/util.h>
#include <stdexcept>

using namespace std;

double rosa::quaternionSuperposition( const vector<Point> &v1,
                                      const vector<Point> &v2,
                                      Mat2D<double> &rotMatrix,
                                      Point &cm1, Point &cm2 )
{
  // Checks the two vectors have the same number of elements
  if( v1.size() != v2.size() )
    throw logic_error( "quaternionSuperposition: sets of points with different "
                       "sizes ("+ltos(v1.size())+" != "+ltos(v2.size())+")" );

  unsigned long vSize = v1.size();

  if( vSize == 0 )
    return -1.0;

  getBaricenter( v1, cm1 );
  getBaricenter( v2, cm2 );
    
  // Creates coordinate differences delta x plus (dxp) and minus (dxm)
  // v1 -> x'
  // v2 -> x
  // x_m = x' - x
  // x_p = x' + x

  Mat2D<double> dxm(vSize,3), dxp(vSize,3);
  
  for( unsigned long k = 0; k < vSize; k++ )
    for( unsigned long i = 0; i < 3; i++ ) {
      dxm[k][i] = v1[k].pos[i] - cm1.pos[i] - (v2[k].pos[i] - cm2.pos[i]);
      dxp[k][i] = v1[k].pos[i] - cm1.pos[i] + (v2[k].pos[i] - cm2.pos[i]);
    }
  
  // Fills the quaternion matrix
  Mat2D<double> q( 4, 4, 0.0 );
  
  for( unsigned long k = 0; k < vSize; k++ ) {
    // Main diagonal
    q[0][0] += mSqr(dxm[k][0]) + mSqr(dxm[k][1]) + mSqr(dxm[k][2]);
    q[1][1] += mSqr(dxp[k][1]) + mSqr(dxp[k][2]) + mSqr(dxm[k][0]);
    q[2][2] += mSqr(dxp[k][0]) + mSqr(dxp[k][2]) + mSqr(dxm[k][1]);
    q[3][3] += mSqr(dxp[k][0]) + mSqr(dxp[k][1]) + mSqr(dxm[k][2]);
    
    // Cross differences
    q[0][1] += dxp[k][1]*dxm[k][2] - dxm[k][1]*dxp[k][2];
    q[0][2] += dxm[k][0]*dxp[k][2] - dxp[k][0]*dxm[k][2];
    q[0][3] += dxp[k][0]*dxm[k][1] - dxm[k][0]*dxp[k][1];
    q[1][2] += dxm[k][0]*dxm[k][1] - dxp[k][0]*dxp[k][1];
    q[1][3] += dxm[k][0]*dxm[k][2] - dxp[k][0]*dxp[k][2];
    q[2][3] += dxm[k][1]*dxm[k][2] - dxp[k][1]*dxp[k][2];
  }

  // Transposition
  q[1][0] = q[0][1];
  q[2][0] = q[0][2];
  q[3][0] = q[0][3];
  q[2][1] = q[1][2];
  q[3][1] = q[1][3];
  q[3][2] = q[2][3];

  // Diagonalize the matrix using the Jacobi method
  vector<double> dm(4);
  Mat2D<double>  vm(4,4);
  int            nrot;
  
  //if( C_MDBG_ON ) q.printf( cout );
  
  if( !jacobiTransformation( q, dm, vm, nrot ) )
    throw runtime_error( "quaternionSuperposition: jacobiTransformation has done too many iterations" );
  
  // Sort eigenvectors after eigenvalues, descending
  eigenValVectSort( dm, vm );
  
  vector<double> ev(4);
  
  for( int i = 0; i < 4; i++ )
    ev[i] = vm[i][3];

  // Fill the rotation matrix
  rotMatrix = Mat2D<double>(3,3);
  
  // ============================== FIRST ROW ==============================
  rotMatrix[0][0] = mSqr(ev[0]) + mSqr(ev[1]) - mSqr(ev[2]) - mSqr(ev[3]);
  rotMatrix[0][1] = 2.0 * (ev[1]*ev[2] + ev[0]*ev[3]);
  rotMatrix[0][2] = 2.0 * (ev[1]*ev[3] - ev[0]*ev[2]);
  
  // ============================= SECOND ROW ==============================
  // R[1][0]
  rotMatrix[1][0] = 2.0 * (ev[1]*ev[2] - ev[0]*ev[3]);
  // R[1][1]
  rotMatrix[1][1] = mSqr(ev[0]) + mSqr(ev[2]) - mSqr(ev[1]) - mSqr(ev[3]);
  // R[1][2]
  rotMatrix[1][2] = 2.0 * (ev[2]*ev[3] + ev[0]*ev[1]);
  
  // ============================== THIRD ROW ==============================
  // R[2][0]
  rotMatrix[2][0] = 2.0 * (ev[1]*ev[3] + ev[0]*ev[2]);
  // R[2][1]
  rotMatrix[2][1] = 2.0 * (ev[2]*ev[3] - ev[0]*ev[1]);
  // R[2][2]
  rotMatrix[2][2] = mSqr(ev[0]) + mSqr(ev[3]) - mSqr(ev[1]) - mSqr(ev[2]);
  
  // Calculate RMSD
  double rmsd = sqrt( dm[3] / (double)vSize);
  
  return rmsd;
}


double rosa::getRmsdQuaternionSuperposition( const vector<Point> &v1,
                                             const vector<Point> &v2 )
{
  // Checks the two vectors have the same number of elements
  if( v1.size() != v2.size() )
    throw logic_error( "quaternionSuperposition: sets of points with different "
                       "sizes ("+ltos(v1.size())+" != "+ltos(v2.size())+")" );

  unsigned long vSize = v1.size();

  if( vSize == 0 )
    return -1.0;

  Point cm1, cm2;
  
  getBaricenter( v1, cm1 );
  getBaricenter( v2, cm2 );
    
  // Creates coordinate differences delta x plus (dxp) and minus (dxm)
  // v1 -> x'
  // v2 -> x
  // x_m = x' - x
  // x_p = x' + x

  Mat2D<double> dxm(vSize,3), dxp(vSize,3);
  
  for( unsigned long k = 0; k < vSize; k++ )
    for( unsigned long i = 0; i < 3; i++ ) {
      dxm[k][i] = v1[k].pos[i] - cm1.pos[i] - (v2[k].pos[i] - cm2.pos[i]);
      dxp[k][i] = v1[k].pos[i] - cm1.pos[i] + (v2[k].pos[i] - cm2.pos[i]);
    }
  
  // Fills the quaternion matrix
  Mat2D<double> q( 4, 4, 0.0 );
  
  for( unsigned long k = 0; k < vSize; k++ ) {
    // Main diagonal
    q[0][0] += mSqr(dxm[k][0]) + mSqr(dxm[k][1]) + mSqr(dxm[k][2]);
    q[1][1] += mSqr(dxp[k][1]) + mSqr(dxp[k][2]) + mSqr(dxm[k][0]);
    q[2][2] += mSqr(dxp[k][0]) + mSqr(dxp[k][2]) + mSqr(dxm[k][1]);
    q[3][3] += mSqr(dxp[k][0]) + mSqr(dxp[k][1]) + mSqr(dxm[k][2]);
    
    // Cross differences
    q[0][1] += dxp[k][1]*dxm[k][2] - dxm[k][1]*dxp[k][2];
    q[0][2] += dxm[k][0]*dxp[k][2] - dxp[k][0]*dxm[k][2];
    q[0][3] += dxp[k][0]*dxm[k][1] - dxm[k][0]*dxp[k][1];
    q[1][2] += dxm[k][0]*dxm[k][1] - dxp[k][0]*dxp[k][1];
    q[1][3] += dxm[k][0]*dxm[k][2] - dxp[k][0]*dxp[k][2];
    q[2][3] += dxm[k][1]*dxm[k][2] - dxp[k][1]*dxp[k][2];
  }

  // Transposition
  q[1][0] = q[0][1];
  q[2][0] = q[0][2];
  q[3][0] = q[0][3];
  q[2][1] = q[1][2];
  q[3][1] = q[1][3];
  q[3][2] = q[2][3];

  // Diagonalize the matrix using the Jacobi method
  vector<double> dm(4);
  Mat2D<double>  vm(4,4);
  int            nrot;
  
  //if( C_MDBG_ON ) q.printf( cout );
  
  if( !jacobiTransformation( q, dm, vm, nrot ) )
    throw runtime_error( "quaternionSuperposition: jacobiTransformation has done too many iterations" );
  
  // Sort eigenvectors after eigenvalues, descending
  eigenValVectSort( dm, vm );
  
  // Calculate RMSD
  double rmsd = sqrt( dm[3] / (double)vSize);
  
  return rmsd;
}
