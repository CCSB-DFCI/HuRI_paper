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
#include <set>

using namespace rosa;
using namespace std;

// Options
static struct option long_options[] = {
  { "version",    no_argument,       0, 'v' },  // Outputs only the version number
  { "help",       no_argument,       0, 'h' },  // Outputs a short guide
  { "output",     required_argument, 0, 'o' },  // Output filename
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] pdbFile chainId start end" << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -v       --version          prints the version number" << endl << endl;
  cout << "   -o FILE  --output=FILE      set the output filename" << endl;
}


int main( int argc, char **argv )
{
  string programName( argv[0] );
  string outputFilename;

  DBG_BLOCK_SUPPRESS( programName+"::main" );

  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":o:vh", long_options, &option_index );

    if( ch == -1 ) // Finished with the options
      break;

    // Options parsing
    switch( ch ) {
      case 'v':
        printRosaVersion( programName );
        return EXIT_SUCCESS;
        break;
      case 'h':
        printRosaVersion( programName );
        usage( programName );
        return EXIT_SUCCESS;
        break;
      case 'o':
        outputFilename = string( optarg );
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

  if( argc < 4 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }

  string pdbFilename( argv[0] );
  string chainId( argv[1] );
  int startIndex = stol(argv[2]);
  int endIndex   = stol(argv[3]);

  if( startIndex <= 0 || endIndex <= 0 || endIndex < startIndex) {
    errorExit( EXIT_FAILURE, "Indices have to be such that 1 <= start <= end" );
  }

  if( outputFilename.empty() )
    outputFilename = baseFilename( pdbFilename ) + "-cut." +
      extFilename( pdbFilename );

  StructurePtr s( new Structure(pdbFilename, true) );

  StructurePtr cs( s->crop(chainId,startIndex-1,endIndex) );

  cs->saveToPDBFile( outputFilename );

  return EXIT_SUCCESS;
}
