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

#include <rosa/dbg.h>
#include <rosa/version.h>
#include <rosa/error.h>
#include <rosa/util.h>
#include <rosa/structure.h>
#include <getopt.h>
#include <vector>
#include <iostream>
#include <cstdlib>

using namespace rosa;
using namespace std;

// Options
static struct option long_options[] = {
  { "version",     no_argument,       0, 'v' },  // Outputs only the version number
  { "dist",        required_argument, 0, 'd' },  // Sets the max. distance of interface residues from the ligand
  { "help",        no_argument,       0, 'h' },  // Outputs a short guide
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] receptorFile ligandFile" << endl << endl;
  cout << "Calculates the interface to the ligand in the receptor" << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -d       --dist DIST        sets the max. distance of interface residues from the ligand (default 6.0 Ang)" << endl;
  cout << "   -v       --version          prints the version number" << endl;
}


string getSelStringChains( const string &chains )
{
  vector<string> tokens;
  split( chains, tokens, "," );
  string selString;

  for( vector<string>::iterator it = tokens.begin(); it != tokens.end(); ++it ) {
    if( it != tokens.begin() )
      selString += " or ";
    selString += "chain "+*it;
  }

  return selString;
}


int main( int argc, char *argv[] )
{
  DBG_BLOCK_SUPPRESS( "bindinterf::main" );

  double maxDist = 6.0;

  string programName( argv[0] );

  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":d:vh", long_options, &option_index );

    if( ch == -1 ) // Finished with the options
      break;

    // Options parsing
    switch( ch ) {
      case 'd':
        maxDist = stod( optarg );
        break;
      case 'v':
        printRosaVersion( programName );
        return EXIT_SUCCESS;
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

  if( argc != 2 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }

  string recFilename( argv[0] );
  string ligFilename( argv[1] );

  // Loads the relevant structures
  Structure recOrig( recFilename );
  //printChainDescription(recOrig);
  Structure ligOrig( ligFilename );
  //printChainDescription(ligOrig);

  StructurePtr interface( selectBindingSite( recOrig, ligOrig, maxDist ) );
  //printChainDescription(*interface);

  long recNumIntRes = 0L;
  for( unsigned long ic = 0UL; ic < interface->numChains(); ++ic ) {
    const Structure::Chain &c = interface->getChain(ic);
    recNumIntRes += c.numResidues();
    for( unsigned long ir = 0UL; ir < c.numResidues(); ++ir )
      cout << c.getChainID() << "\t" << c.getResLabel(ir) << endl;
  }
  return EXIT_SUCCESS;
}
