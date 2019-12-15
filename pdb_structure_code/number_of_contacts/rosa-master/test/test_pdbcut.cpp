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
#include <rosa/structure.h>
#include <rosa/util.h>
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

class TestPdbcut: public Test {
  string filesPath;
public:
  TestPdbcut( const string &aPath ):
    filesPath(aPath)
  {}
  
  void run() {
    /* Reading correct lines */
    string inputPdbFilename = joinPath(filesPath,"1atp.pdb");
    Structure is( inputPdbFilename );
    StructurePtr os( is.crop("E",71,213) );

    SETUP_TEST_EQ_LBL((int)os->numChains(),1,"Cropped file has 1 chain");
    SETUP_TEST_EQ_LBL((int)os->numModels(),1,"Cropped file has 1 model");
    SETUP_TEST_EQ_LBL(os->getChain(0).getChainID(),string("E"),"Chain E");
    SETUP_TEST_EQ_LBL((int)os->getChain(0).numResidues(),142,"Cropped file has 142 residues");
    SETUP_TEST_EQ_LBL(os->getChain(0).getDescription(),string("E86-E227"),"Chain description");
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

  TestGroup tgrp( "Atom tests" );
  shared_ptr<TestPdbcut> pdbcutTest( new TestPdbcut(path) );
  
  tgrp.addTest( pdbcutTest );

  tgrp.run();
  
  long nFail = tgrp.outputReport();
  
  return nFail;
}
