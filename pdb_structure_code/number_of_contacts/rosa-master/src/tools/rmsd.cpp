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

#include <rosa/structure.h>
#include <rosa/version.h>
#include <rosa/error.h>
#include <rosa/dbg.h>
#include <getopt.h>
#include <iostream>
#include <cstdlib>

using namespace rosa;
using namespace std;

// Options
static struct option long_options[] = {
  { "version",    no_argument,       0, 'v' },  // Outputs only the version number
  { "all",        no_argument,       0, 'a' },  // Rmsd on all the atoms
  { "help",       no_argument,       0, 'h' },  // Outputs a short guide  
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] pdbFile1 PdbFile2" << endl << endl;
  cout << "The RMSD is normally calculated only on the CAs. If you want to calculate it on" << endl;
  cout << "all the atoms use the '-a' option." << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -v       --version          prints the version number" << endl;
  cout << "   -a       --all              superpose on all the atoms" << endl << endl;
}

 
int main( int argc, char **argv )
{
  DBG_BLOCK_SUPPRESS( "rmsd::main" );
  
  string programName( argv[0] );
  bool bAll = false;
  
  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":vah", long_options, &option_index );

    if( ch == -1 ) // Finished with the options
      break;
    
    // Options parsing
    switch( ch ) {
      case 'v':
        printRosaVersion( programName );
        return EXIT_SUCCESS;
        break;
      case 'a':
        bAll = true;
        break;
      case 'h':
        printRosaVersion( programName );
        usage( programName );
        return EXIT_SUCCESS;
        break;
      case '?':
      default:
        printRosaVersion( programName );
        usage( programName );
        errorExit( EXIT_FAILURE, "Invalid command line argument" );
        break;
    }
  }

  if( bAll )
    cout << "Using all the atoms" << endl;
  else
    cout << "Using only CA atoms" << endl;
    
  argc -= optind;
  argv += optind;

  DBG_VDUMP( argc );
  
  if( argc != 2 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }
  
  DBG_MSG( "Mark1" );
  
  string pdbFilename1( argv[0] );
  string pdbFilename2( argv[1] );
  
  DBG_MSG( "Mark2" );

  Structure s1, s2;
  
  s1.loadFromPDBFile( pdbFilename1 );
  s2.loadFromPDBFile( pdbFilename2 );

  DBG_MSG( "Mark3" );

  double rmsd;
  long numAtoms;
  
  try {
    if( !bAll ) {
      DBG_VDUMP( s1.size() );
      DBG_VDUMP( s1.numModels() );
      
      StructurePtr selS1 = s1.select( "name CA" );
      StructurePtr selS2 = s2.select( "name CA" );

      reduceToCommonAtoms ( *selS1, *selS2 );
      
      numAtoms = selS1->size();
      
      DBG_VDUMP( selS1->size() );
      DBG_VDUMP( selS1->numModels() );

      rmsd = calcRmsd( *selS1, *selS2 );
      
      DBG_MSG( "Mark5" );
    } else {
        numAtoms = s1.size();
        rmsd = calcRmsd( s1, s2 );
    }
    
    cout << "Rmsd: " << fixed << setprecision(2) << rmsd << " Ang. (" << numAtoms << " atoms)" << endl;
  } catch( exception &e ) {
    errorExit( EXIT_FAILURE, e.what() );
  }
  
  return EXIT_SUCCESS;
}
