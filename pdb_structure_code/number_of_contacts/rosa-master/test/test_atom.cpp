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
#include <rosa/elems.h>
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

class TestAtom: public Test {
  string filesPath;
public:
  TestAtom( const string &aPath ):
    filesPath(aPath)
  {}
  
  void run() {
    /* Reading correct lines */
    string correctFilename = joinPath(filesPath,"correct_atoms.pdb");
    ifstream correctFile( correctFilename.c_str() );
    if( !correctFile.good() )
      throw runtime_error( "error while opening file "+correctFilename );

    while( !correctFile.eof() ) {
      string line;
      if( getline( correctFile, line ) ) {
        Atom a;
        SETUP_TEST_EXEC( a.readPDBAtomRecord(line), string("Correct atom line \'")+line+"\'" );
      }
    }
    correctFile.close();

    /* Reading incorrect lines */
    string incorrectFilename = joinPath(filesPath,"incorrect_atoms.pdb");
    ifstream incorrectFile( incorrectFilename.c_str() );
    if( !incorrectFile.good() )
      throw runtime_error( "error while opening file "+incorrectFilename );

    while( !incorrectFile.eof() ) {
      string line;
      if( getline( incorrectFile, line ) ) {
        Atom a;
        SETUP_TEST_EXCPT( a.readPDBAtomRecord(line), string("Incorrect atom line \'")+line+"\'" );
      }
    }
    incorrectFile.close();

    // Parsing the fields
    string atStr = "ATOM    370  CB AARG A  45B     39.578   2.601  21.553  0.53 52.02           C  ";
    Atom a; a.readPDBAtomRecord(atStr);

    SETUP_TEST_EQ_LBL(a.getName(),string("CB"),"Name");
    SETUP_TEST_EQ_LBL(a.getOrigName(),string(" CB "),"OrigName");
    SETUP_TEST_EQ_LBL(a.getSerialNum(),370L,"SerialNum");
    SETUP_TEST_EQ_LBL(a.getAltLoc(),string("A"),"AltLoc");
    SETUP_TEST_EQ_LBL(a.getInsCode(),string("B"),"InsCode");
    SETUP_TEST_EQ_LBL(a.getResName(),string("ARG"),"ResName");
    SETUP_TEST_EQ_LBL(a.getResNum(),45L,"ResNum");
    SETUP_TEST_EQ_LBL(a.getChainID(),string("A"),"OrigName");
    SETUP_TEST_EQ_LBL(trim(a.getPdbRecordName()),string("ATOM"),"RecordName");
    SETUP_TEST_EQ_LBL(a.getPosX(),(Coord)39.578,"PosX");
    SETUP_TEST_EQ_LBL(a.getPosY(),(Coord)2.601,"PosY");
    SETUP_TEST_EQ_LBL(a.getPosZ(),(Coord)21.553,"PosZ");
    SETUP_TEST_EQ_LBL(a.getBfac(),(float)52.02,"Bfac");
    SETUP_TEST_EQ_LBL(a.getOcc(),(float)0.53,"Occ");
    SETUP_TEST_LBL(!a.isTer(),"isTer");
    SETUP_TEST_EQ_LBL(a.getElemSymbol(),string("C"),"ElemSymbol");
    SETUP_TEST_EQ_LBL(a.getElemNo(),(int)Elem_C,"ElemNo");
    SETUP_TEST_EQ_LBL(a.identifier(),string("_CB__A_ARG_45_B_A"),"identifier");
    SETUP_TEST_EQ_LBL(trim(Atom::getResLabel(a.getResNum(),a.getInsCode())),string("45B"),"ResLabel");

    std::ostringstream outStream;
    a.writePDBAtomRecord(outStream);
    SETUP_TEST_EQ_LBL(outStream.str(),string("ATOM    370  CB AARG A  45B     39.578   2.601  21.553  0.53 52.02           C  "),"writePDBAtomRecord");

    Atom b;
    b.readPDBAtomRecord(outStream.str());
    SETUP_TEST_LBL(!compareAtoms(a,b),"compareAtoms");
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
  shared_ptr<TestAtom> atomTest( new TestAtom(path) );
  
  tgrp.addTest( atomTest );

  tgrp.run();
  
  long nFail = tgrp.outputReport();
  
  return nFail;
}
