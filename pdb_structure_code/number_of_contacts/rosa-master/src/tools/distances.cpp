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
#include <rosa/mathutil.h>
#include <rosa/atomsel_driver.h>
#include <rosa/error.h>
#include <rosa/dbg.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace rosa;
using namespace std;

// Options
static struct option long_options[] = {
  { "version",    no_argument,         0, 'v' },  // Outputs only the version number
  { "first",      no_argument,         0, 'f' },  // Distance between the first atom and all the others
  { "surf",       no_argument,         0, 's' },  // Surface distance instead of euclidean distance
  { "help",       no_argument,         0, 'h' },  // Outputs a short guide
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] selection pdbFile" << endl << endl;
  cout << "This program selects a set of atoms in a structure and calculates the distances" << endl;
  cout << "between them. By default it calculates all pairwise distances." << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -v       --version          prints the version number" << endl;
  cout << "   -f       --first            only calculates distances between the first atom and all the others" << endl;
  cout << "   -s       --surf             calculates the minimum surface distances between the the atoms" << endl;
}


int main( int argc, char **argv )
{
  DBG_BLOCK_SUPPRESS( "distances::main" );

  string programName( argv[0] );

  bool bFirst    = false;
  bool bSurfDist = false;

  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, "vsfh", long_options, &option_index );

    if( ch == -1 ) // Finished with the options
      break;

    // Options parsing
    switch( ch ) {
      case 'v':
        printRosaVersion( programName );
        return EXIT_SUCCESS;
        break;
      case 'f':
        bFirst = true;
        break;
      case 's':
        bSurfDist = true;
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

  argc -= optind;
  argv += optind;

  DBG_VDUMP( argc );

  if( argc != 2 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }

  string selectionString( argv[0] );
  string aPdbFilename( argv[1] );

  AtomSelectDriver selDriver( selectionString );
  DBG_VDUMP( selDriver.good() );
  if( !selDriver.good() ) {
    errorExit( EXIT_FAILURE, "Invalid selection string. Please check the syntax..." );
  }

  Structure aStruct;

  aStruct.loadFromPDBFile( aPdbFilename );
  if( DBG_ON ) {
    DBG_MSG( "Loaded structure" )
    printChainDescription( aStruct );
  }

  aStruct.removeMultConformations();

  if( DBG_ON ) {
    DBG_MSG( "After remove mult conformations" )
    printChainDescription( aStruct );
  }

  try {
    StructurePtr selAtoms = aStruct.select( selectionString );
    if( DBG_ON )
      printChainDescription( *selAtoms );

    // First we check for disulfide bridges interactions
    cout << "NAME1\tRESI1\tRESN1\tCHAIN1\tNAME2\tRESI2\tRESN2\tCHAIN2\tDIST" << endl;
    for( long int i = 0; (i < (long int)(selAtoms->size())-1) && (!bFirst || i < 1); ++i ) {
      DBG_VDUMP(i);
      for( long int j = i+1; j < (long int)(selAtoms->size()); ++j ) {
        DBG_VDUMP(j);
        const Atom &iA = (*selAtoms)[i];
        const Atom &jA = (*selAtoms)[j];
        double lDist;
        if( bSurfDist ) {
          unsigned long currentChainIndex = 0;
          for( unsigned long ci = 0; ci < aStruct.numChains(); ++ci ) {
            if( aStruct.getChain(ci).getChainID() == iA.getChainID() ) {
              currentChainIndex = ci;
              break;
            }
          }

          const Structure::Chain &iChain = aStruct.getChain( currentChainIndex );
          lDist = iChain.getSurfaceDistance( iA.getPos(), jA.getPos() );
        } else {
          lDist = pointDistance( iA.getPos(), jA.getPos() );
        }
        cout << iA.getName() << "\t"<< trim(Atom::getResLabel(iA.getResNum(),iA.getInsCode())) << "\t" << iA.getResName() << "\t" << iA.getChainID() << "\t"
             << jA.getName() << "\t"<< trim(Atom::getResLabel(jA.getResNum(),jA.getInsCode())) << "\t" << jA.getResName() << "\t" << jA.getChainID() << "\t"
             << lDist << endl;
      }
    }

  } catch( exception &e ) {
    errorExit( EXIT_FAILURE, e.what() );
  }

  return EXIT_SUCCESS;
}
