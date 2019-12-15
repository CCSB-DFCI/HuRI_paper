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

#include "test.h"
#include <rosa/atom.h>
#include <rosa/structure.h>
#include <rosa/mathutil.h>
#include <rosa/atomsel_driver.h>
#include <sys/times.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <algorithm>

using namespace std;
using namespace rosa;
using namespace TestToolbox;

const double MY_EPS = 0.001;

class TestTransform: public Test {
private:
  Structure s;
  
  void loadAtomsFile( const string &atomsFilename ) {
    s.loadFromPDBFile( atomsFilename );
  }

public:
  TestTransform( const string &atomsFilename ) {
    loadAtomsFile( atomsFilename );
  }

  void run() {
    Point origBaricenter( s.getBaricenter() );
    
    Structure t( s );
    
    t.translate( -origBaricenter );
    Point newBaricenter1( t.getBaricenter() );
    Point origin;
    double dist1 = pointDistance( newBaricenter1, origin );
    SETUP_TEST_LEQ_LBL( dist1, MY_EPS, "Baricenter in origin" );

    Structure t1( t );
    t.rotate( rotMatrix( 2 * M_PI, 2 * M_PI, 2 * M_PI ) );
    unsigned long numDiff = t.compareAtoms( t1 );
    SETUP_TEST_EQ_LBL( numDiff, 0UL, "Rotation of (2pi, 2pi, 2pi)" );

    t.rotate( rotMatrix( (rand01() - 0.5) * M_PI,
                         (rand01() - 0.5) * M_PI,
                         (rand01() - 0.5) * M_PI ) );

    Point newBaricenter2( t.getBaricenter() );
    double dist2 = pointDistance( newBaricenter2, newBaricenter1 );
    SETUP_TEST_LEQ_LBL( dist2, MY_EPS, "Baricenter after rotation" );
    
    t.translate( Point( 1.0, 1.0, 1.0 ) );
    Point newBaricenter3( t.getBaricenter() );
    double dist3 = pointDistance( newBaricenter3, origin ) - sqrt(3.0);
    SETUP_TEST_LEQ_LBL( dist3, MY_EPS, "Baricenter in (1.0,1.0,1.0)" );

    t.translate( Point( 0.0, 0.0, 0.0 ) );
    Point newBaricenter4( t.getBaricenter() );
    double dist4 = pointDistance( newBaricenter4, newBaricenter3 );
    SETUP_TEST_LEQ_LBL( dist3, MY_EPS, "Translation of (0.0,0.0,0.0)" );

  }
              
};

int main( int argc, char **argv )
{
  if( argc != 2 ) {
    return EXIT_FAILURE;
  }
  
  string path( argv[1] );
  
  if( !path.empty() ) {
    if( path[path.size()-1] != '/' )
      path += '/';
  }
  
  string filename1( joinPath(path,"1atp.pdb") );
  
  TestGroup tgrp( "Transform tests" );
  shared_ptr<TestTransform> transTest1( new TestTransform( filename1 ) );
  
  tgrp.addTest( transTest1 );

  tgrp.run();
  
  long nFail = tgrp.outputReport();
  
  return nFail;
}
