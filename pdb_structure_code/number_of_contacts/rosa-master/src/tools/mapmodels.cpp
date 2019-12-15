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

#include <rosa/version.h>
#include <rosa/error.h>
#include <rosa/dbg.h>
#include <rosa/util.h>
#include <rosa/structure.h>
#include <getopt.h>
#include <cstdlib>

using namespace rosa;
using namespace std;

// Options
static struct option long_options[] = {
  { "version",    no_argument,       0, 'v' },  // Outputs only the version number
  { "help",       no_argument,       0, 'h' },  // Outputs a short guide  
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] pdbFile" << endl << endl;
  cout << "Prints a mapping between internal model indexes and the MODEL number" << endl;
  cout << "provided by the PDB file. If the PDB file contains no model the index 0" << endl;
  cout << "is indicating the default model and is mapped to \"-\"." << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -v       --version          prints the version number" << endl << endl;
}

 
int main( int argc, char **argv )
{
  string programName( argv[0] );
  
  DBG_BLOCK_SUPPRESS( programName+"::main" );
  
  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":vh", long_options, &option_index );

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
      case '?':
      default:
        printRosaVersion( programName );
        usage( programName );
        string opChar = string(1,(char)ch);
        errorExit( EXIT_FAILURE, "Invalid command line argument '"+opChar+"'" );
        break;
    }
  }

  argc -= optind;
  argv += optind;


  DBG_VDUMP( argc );
  
  if( argc != 1 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid number of arguments" );
  }
  
  DBG_MSG( "Mark1" );
  

  string pdbFilename( argv[0] );
  
  Structure s;
  s.loadFromPDBFile( pdbFilename, true );
  //cout << "INDEX    MODEL_ID" << endl;
  for( unsigned long i = 0; i < s.numModels(); ++i ) {
    int modelId = s.getModelId(i);
    cout << i << "\t";
    if( modelId < 0 )
      cout << "-" << endl;
    else
      cout << modelId << endl;
  }
  
  return EXIT_SUCCESS;
}
