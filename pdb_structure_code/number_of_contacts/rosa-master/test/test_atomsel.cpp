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

class TestAtomSelection: public Test {
  typedef pair<string,long> SelTestPair;
  typedef vector<SelTestPair> TestPairVect;
  typedef vector<SelTestPair>::iterator TestPairVectIt;
  typedef vector<SelTestPair>::const_iterator TestPairVectCIt;
  
  TestPairVect selections;
public:
  TestAtomSelection( const string &atomsFilename ) {
    loadAtomsFile( atomsFilename );
  }

  void addSelTestPair( const string &selString, long aCount ) {
    selections.push_back( SelTestPair( selString,aCount) );
  }
  
  void run() {
    TestPairVectCIt selStringsEnd = selections.end();
    for( TestPairVectCIt it = selections.begin(); it != selStringsEnd;
         ++it ) {
      long counter = countSelected( it->first );
      SETUP_TEST_EQ_LBL( counter, it->second, it->first );
    }
  }
              
private:
  AtomVect atomVector;
  
  long countSelected( const string &aSelStr ) {
    AtomSelectDriver selDriver( aSelStr );
    long counter = 0L;
    AtomVectCit atVectEnd = atomVector.end();
    for( AtomVectCit it = atomVector.begin(); it != atVectEnd; ++it )
      if( selDriver.evalOnAtom( *it ) )
        ++counter;
    
    return counter;
  }
  
  void loadAtomsFile( const string &atomsFilename ) {
    string atomLine;

    ifstream atomsFile( atomsFilename.c_str() );
    while( getline(atomsFile, atomLine) ) {
      Atom at;
      at.readPDBAtomRecord( atomLine );
      atomVector.push_back(at);
    }
    atomsFile.close();
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
  
  string filename1( joinPath(path,"atoms1.pdb") );
  
  TestGroup tgrp( "Atom selection tests" );
  shared_ptr<TestAtomSelection> selTest1( new TestAtomSelection( filename1 ) );
  
  selTest1->addSelTestPair( "resi 30:87", 2886 );
  selTest1->addSelTestPair( "resi -5:-1",  150 );
  selTest1->addSelTestPair( "resi -2:2",   165 );
  selTest1->addSelTestPair( "resi 23",      66 );
  selTest1->addSelTestPair( "resi -4",      30 );
  selTest1->addSelTestPair( "resi = 116",   33 );
  selTest1->addSelTestPair( "resi <> 157",9600 );
  selTest1->addSelTestPair( "resi < 152", 7350 );
  selTest1->addSelTestPair( "resi <= 62", 3212 );
  selTest1->addSelTestPair( "resi > 194",  237 );
  selTest1->addSelTestPair( "resi >= 194", 285 );
  selTest1->addSelTestPair( "Bf = 20.0",     5 );
  selTest1->addSelTestPair( "Bf <> 20.0", 9637 );
  selTest1->addSelTestPair( "Bf < 40.0",  6097 );
  selTest1->addSelTestPair( "Bf <= 50.0", 7465 );
  selTest1->addSelTestPair( "Bf > 30.0",  5836 );
  selTest1->addSelTestPair( "Bf >= 35.0", 4573 );
  selTest1->addSelTestPair( "name CA",    1187 );
  selTest1->addSelTestPair( "name = C",   1187 );
  selTest1->addSelTestPair( "name <> N",  8446 );
  selTest1->addSelTestPair( "resn CYS",    108 );
  selTest1->addSelTestPair( "resn = ALA",  187 );
  selTest1->addSelTestPair( "resn <> TYR",9102 );
  selTest1->addSelTestPair( "chain A",    1700 );
  selTest1->addSelTestPair( "chain = B",    92 );
  selTest1->addSelTestPair( "chain <> C", 7942 );
  selTest1->addSelTestPair( "chain A or resi 12",                 1736 );
  selTest1->addSelTestPair( "chain A or resi 12 or name CA",      2708 );
  selTest1->addSelTestPair( "resn ALA and name CA",                 37 );
  selTest1->addSelTestPair( "resn ALA and name CA and chain A",      8 );
  selTest1->addSelTestPair( "not resn ALA",                       9455 );
  selTest1->addSelTestPair( "not(not resn ALA)",                   187 );
  selTest1->addSelTestPair( "resn ALA and (name CA or name C)",     74 );
  selTest1->addSelTestPair( "resn ALA and not (name CA or name C)",113 );

  tgrp.addTest( selTest1 );

  string filename2( path+"atoms2.pdb" );
  shared_ptr<TestAtomSelection> selTest2( new TestAtomSelection( filename2 ) );

  selTest2->addSelTestPair( "occ = 1.0",  12380 );
  selTest2->addSelTestPair( "occ <> 1.0",    20 );
  selTest2->addSelTestPair( "occ < 0.5",     15 );
  selTest2->addSelTestPair( "occ <= 0.6",    16 );
  selTest2->addSelTestPair( "occ > 0.4",  12385 );
  selTest2->addSelTestPair( "occ >= 0.64",12384 );
  /*selTest2->addSelTestPair( "ins B",10 );
  selTest2->addSelTestPair( "ins = B",10 );
  selTest2->addSelTestPair( "ins <> A",10 );
  selTest2->addSelTestPair( "alt A",      15 );
  selTest2->addSelTestPair( "alt = B",    15 );
  selTest2->addSelTestPair( "alt <> A",12385 );
  */
  tgrp.addTest( selTest2 );

  tgrp.run();
  
  long nFail = tgrp.outputReport();
  
  return nFail;
}

/*
  cout << "Read " << atomVector.size() << " atoms." << endl << endl;

  //copy( atomVector.begin(), atomVector.end(), ostream_iterator<Atom>(cout, "\n"));

  cout << "Insert a selection string: ";
  
  string selString;
  getline( cin, selString );

  struct tms startTmsTime, stopTmsTime;
  const long SYS_CLK_TCK = sysconf(_SC_CLK_TCK);

  cout << "Selection: " << selString << endl;

  AtomSelectDriver selDriver( selString );
  
  if( !selDriver.good() ) {
    cerr << "Error while parsing the string!" << endl;
    return EXIT_FAILURE;
  }
  //selDriver.setDebug( true );

  long counter = 0;
  bool bSelected = false;
  times( &startTmsTime );
  for( int t = 0; t < 1000; ++t )
    for( long i = 0L; i < (long) atomVector.size(); ++i ) {
      if( selDriver.evalOnAtom( atomVector[i] ) ) ++counter;
      //  cout << atomVector[i] << endl;
      //cout << "Parser output: " << parseAtomSel( selString.c_str(), &atomVector[i], bSelected ) << endl;
      //cout << "Selected: " << (bSelected?"YES":"NO") << endl;
    }
  times( &stopTmsTime );
  float elapsedTime = stopTmsTime.tms_utime - startTmsTime.tms_utime +
                      stopTmsTime.tms_stime - startTmsTime.tms_stime +
                      stopTmsTime.tms_cutime - startTmsTime.tms_cutime +
                      stopTmsTime.tms_cstime - startTmsTime.tms_cstime;
  
  elapsedTime /= (double)SYS_CLK_TCK;
  
  cout << "Elapsed time: " << elapsedTime << endl;
  cout << "counter:      " << counter     << endl;
  
  return EXIT_SUCCESS;
}
*/
