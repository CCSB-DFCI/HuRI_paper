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

/*! \file mathutil.h
 *  \brief Contains mathematical functions and data types.
 */

#ifndef ROSA_MATHUTIL_H_
#define ROSA_MATHUTIL_H_

#include <rosa/types.h>
#include <rosa/mem_shared_ptr.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <cstdlib>
#include <algorithm>

namespace rosa {
  
  //! Implements a 2 dimensional matrix
  template<typename T>
  class Mat2D {
  private:
    unsigned long rows, cols;
    
    shared_ptr<T> mm;
    
    void copy( const Mat2D<T> &rhs ) {
      if( rhs.mm != 0 ) {
        shared_ptr<T> temp( new T[rhs.rows*rhs.cols] );
        for( unsigned long i = 0; i < rhs.rows*rhs.cols; ++i )
          temp.get()[i] = rhs.mm.get()[i];
        mm = temp;
      } else
        mm.reset();
      
      rows = rhs.rows;
      cols = rhs.cols;
    }
    
  public:
    //! Creates an empty matrix
    Mat2D(): rows( 0 ), cols( 0 ) {}

    Mat2D( const Mat2D<T> &rhs ):
      rows( rhs.rows ), cols( rhs.cols ), mm( rhs.mm!=0?new T[rhs.rows*rhs.cols]:0 )
    {
      for( unsigned long i = 0; i < rows*cols; ++i )
        mm.get()[i] = rhs.mm.get()[i];
    }
    
    //! Creates a matrix of dimensions (aRows, aCols) without initializing
    //! the single elements.
    Mat2D( unsigned long aRows, unsigned long aCols ):
      rows( aRows ), cols( aCols ), mm( new T[aRows*aCols] )
    {}
    
    //! Creates a matrix of dimensions (aRows, aCols) and initialize every
    //! element to initValue.
    Mat2D( unsigned long aRows, unsigned long aCols, T initValue ):
      rows( aRows ), cols( aCols ), mm( new T[aRows*aCols] ) 
    {
      for( unsigned long i = 0; i < aRows*aCols; i++ )
        mm.get()[i] = initValue;
    }
    
    ~Mat2D() {
    }
    
    //! Prints the matrix to output. If printHeadings is true the indices for
    //! rows and columns are printed. By default the matrix is printed to cout
    //! but it can be printed to any other stream passed as the second argument.
    void print( bool printHeadings = true, std::ostream &os = std::cout ) const {
      std::ios::fmtflags bakFlags = os.flags();
      
      os << std::fixed;
      
      if( printHeadings ) {
        os << "   ";
        for( unsigned long j = 0; j < cols; ++j )
          os << std::setw(8) << j;
        os << std::endl;
      }
      
      for( unsigned long i = 0; i < rows; ++i ) {
        if( printHeadings ) os << std::setw(3) << i;
        for( unsigned long j = 0; j < cols; ++j )
          os << std::setw(8) << std::setprecision(2) << (mm.get())[i*cols+j];
        os << std::endl;
      }
      os.flags( bakFlags );
    }
    
    //! Returns a pointer to the first element of the i-th row.
    T *operator [] ( unsigned long i ) { return mm.get()+i*cols; }
    //! Returns a 'const' pointer to the first element of the i-th row.
    const T *operator [] ( unsigned long i ) const { return mm.get()+i*cols; }
    
    Mat2D<T> &operator = ( const Mat2D<T> &rhs ) {
      if( this != &rhs )
        copy( rhs );

      return (*this);
    }
    
    Mat2D<T> operator += ( const Mat2D<T> &rhs ) {
      if( rows != rhs.rows || cols != rhs.cols )
        throw std::logic_error( "Mat2D: cannot add two matrices with different sizes" );
      
      for( unsigned long i = 0; i < rows*cols; ++i )
        mm.get()[i] += rhs.mm.get()[i];
      
      return (*this);
    }
    
    //! Returns the number of rows.
    unsigned long nrows() const { return rows; }
    //! Returns the number of columns.
    unsigned long ncols() const { return cols; }
    
    //! Finds the position of the maximum element in the matrix
    T getMax( unsigned long &i, unsigned long &j ) const {
      unsigned long ij = std::max_element( mm.get(), mm.get()+rows*cols ) - mm.get();
      i = ij / cols;
      j = ij % cols;
      return mm.get()[ij];
    }
  };

  //! Returns an identity matrix of dimension (dim, dim)
  template<class T>
  Mat2D<T> identityMatrix( unsigned long dim ) {
    Mat2D<T> m( dim, dim, T(0.0) );
    for( unsigned long i = 0UL; i < dim; ++i ) m[i][i] = T(1.0);
    return m;
  }
  
  typedef enum RotRefSystem_tag{
    ROT_XYZ=1,
    ROT_ZXZ
  } RotRefSystem;
  
  //! Returns the rotation matrix corresponding to the transformation given by
  //! - one rotation around the first axis of aAng radiants
  //! - one rotation around the second axis of bAng radiants
  //! - one rotation around the third axis of cAng radiants
  //! The rotation are done around the axis depending on the rotation reference
  //! system (xyz or zxz).
  Mat2D<Coord> rotMatrix( double aAng, double bAng, double cAng,
                          RotRefSystem rotSystem = ROT_XYZ );
  
  //! Returns the three angles corresponding to the transformation given by the
  //! rotation matrix m.
  void rotAngles( const Mat2D<double> &m,
                  double &xAng, double &yAng, double &zAng );

  //! Returns a random number in the interval [0.0, 1.0[
  inline double rand01() { return ((double)std::rand() / (double)(RAND_MAX+1.0)); }

  //! Returns a random number in the interval [0.0, 1.0]
  inline double rand01c() { return ((double)std::rand() / (double)RAND_MAX); }

  //! Returns the square of the value
  template<class T>
  inline T mSqr( T a ) { return (a*a); }
  
  //! Sorts the eigenvalues stored in the vector d  in descending order. The
  //! corresponding eigenvectors contained in the matrix v are sorted accordingly.
  void eigenValVectSort( std::vector<double> &d, Mat2D<double> &v);
  
  //! Applies the Jacobi transformation to the matrix a returning the
  //! eigenvalues and eigenvectors in d and v respectively
  bool jacobiTransformation( Mat2D<double> &a, std::vector<double> &d,
                             Mat2D<double> &v, int &nrot );

} // namespace rosa

#endif // ROSA_MATHUTIL_H_
