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

#include <rosa/atomsel_driver.h>
#include <rosa/version.h>
#include <rosa/error.h>
#include <rosa/dbg.h>
#include <rosa/util.h>
#include <rosa/pdb.h>
#include <rosa/atom.h>
#include <getopt.h>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <set>

using namespace rosa;
using namespace std;

typedef set<string> StringSet;

typedef map<string,StringSet> AtomsInResMap;
typedef map<string,StringSet>::iterator AtomsInResMapIt;
typedef map<string,StringSet>::const_iterator AtomsInResMapCit;

// Options
static struct option long_options[] = {
  { "version",    no_argument,       0, 'v' },  // Outputs only the version number
  { "header",     no_argument,       0, 'H' },  // Include headers
  { "output",     required_argument, 0, 'o' },  // Output filename
  { "oneconf",    no_argument,       0, 'O' },  // Removes multiple conformations after loading
  { "stopatter",  no_argument,       0, 'T' },  // Do not include the contents of a chain after the first TER record
  { "model",      required_argument, 0, 'm' },  // Model
  { "help",       no_argument,       0, 'h' },  // Outputs a short guide
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] select_string pdbFile" << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -O       --oneconf          only select one conformation for every residue and clears the alt loc" << endl;
  cout << "   -H       --header           includes headers" << endl;
  cout << "   -T       --stopatter        do not include the contents of a chain after the first TER record" << endl;
  cout << "   -m MODEL --model            select only the specified model" << endl;
  cout << "   -o FILE  --output=FILE      set the output filename" << endl;
  cout << "   -v       --version          prints the version number" << endl << endl;
}


int main( int argc, char **argv )
{
  string programName( argv[0] );

  bool bHeader    = false;
  bool bOneConf   = false;
  bool bStopAtTer = false;
  long modelToSelect = -1;
  string outputFilename;

  DBG_BLOCK_SUPPRESS( programName+"::main" );

  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":o:m:HOTvh", long_options, &option_index );

    if( ch == -1 ) // Finished with the options
      break;

    // Options parsing
    switch( ch ) {
      case 'v':
        printRosaVersion( programName );
        return EXIT_SUCCESS;
        break;
      case 'O':
        bOneConf = true;
        break;
      case 'T':
        bStopAtTer = true;
        break;
      case 'H':
        bHeader = true;
        break;
      case 'o':
        outputFilename = string( optarg );
        break;
      case 'm':
        modelToSelect = stol( optarg );
        if( modelToSelect <= 0 )
          errorExit( EXIT_FAILURE, "Invalid model number, must be an integer "
                     "number greater than 0." );
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

  if( argc != 2 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid number of arguments" );
  }

  DBG_MSG( "Mark1" );

  string selectionString( argv[0] );
  string pdbFilename( argv[1] );

  if( outputFilename.empty() )
    outputFilename = baseFilename( pdbFilename ) + "-sel." +
                     extFilename( pdbFilename );

  DBG_MSG( "Mark2" );

  //AtomSelectDriver selDriver( selectionString );
  AtomSelectDriver selDriver;

  //if( DBG_ON )
  //  selDriver.setDebug( true );

  DBG_VDUMP( selectionString );
  selDriver.parseString( selectionString );


  if( !selDriver.good() )
    errorExit( EXIT_FAILURE, "Error while parsing selection string "+selectionString );

  ifstream inFile( pdbFilename.c_str() );

  if( !inFile.good() )
    errorExit( EXIT_FAILURE, "Error while opening file "+pdbFilename );

  ofstream outFile( outputFilename.c_str() );

  if( !outFile.good() )
    errorExit( EXIT_FAILURE, "Error while opening output file "+outputFilename );

  map<string,string > altLocInRes;
  AtomsInResMap       atomsSeenInRes;
  set<string>         terminatedChains;

  bool bOutput = true;
  while( !inFile.eof() ) {
    string line;
    DBG_MSG( "Reading next line" );
    if( getline( inFile, line ) ) {
      DBG_MSG( "Line read." );
      DBG_VDUMP( line.size() );

      // Pads the length of the line to reach the normal length of 80 characters
      line += string( max<int>( PDB_LINE_LENGTH - line.size(), 0 ), ' ' );

      DBG_MSG( "Mark1" );
      DBG_VDUMP( line );
      string recordName( line.substr( 0, 6 ) );


      if( recordName == TerRecordName ) {
        Atom at;
        at.readPDBAtomRecord( line );
        if( bOutput ) {
          if( selDriver.evalOnAtom( at ) ) {
            outFile << line << endl;
          }
        }
        terminatedChains.insert(at.getChainID());
      } else if( recordName == AtomRecordName or
                 recordName == HetatmRecordName ) {
        Atom at;
        at.readPDBAtomRecord( line );
        if( bStopAtTer && terminatedChains.find(at.getChainID()) != terminatedChains.end() )
          continue;
        if( bOutput ) {
          string currentRes = at.getChainID() + "/" + Atom::getResLabel( at.getResNum(), at.getInsCode() );

          if( bOneConf ) {
            AtomsInResMapIt cra = atomsSeenInRes.find(currentRes);
            if( cra != atomsSeenInRes.end() ) {
              if( cra->second.find(at.getName()) != cra->second.end() ) {
                continue;
              } else {
                if( altLocInRes.find(currentRes) != altLocInRes.end() ) {
                  if( at.getAltLoc() != altLocInRes[currentRes] )
                    continue;
                } else {
                  atomsSeenInRes[currentRes].insert(at.getName());
                  if( at.getAltLoc() != " " )
                    altLocInRes[currentRes] = at.getAltLoc();
                }
              }
            } else {
              atomsSeenInRes.insert(pair<string,StringSet >(currentRes,StringSet()));
              atomsSeenInRes[currentRes].insert(at.getName());
              if( at.getAltLoc() != " " )
                altLocInRes[currentRes] = at.getAltLoc();
            }
          }

          if( selDriver.evalOnAtom( at ) ) {
            if( bOneConf )
              line[16] = ' ';
            outFile << line << endl;
          }
        }
      } else if( recordName == ModelRecordName ) {
        long modelNumber = stol( trim(line.substr(6)) );
        if( modelNumber == modelToSelect or modelToSelect < 0 ) {
          bOutput = true;
          outFile << line << endl;
        } else
          bOutput = false;
        terminatedChains.clear();
      } else if( recordName == EndmdlRecordName ) {
        if( bOutput ) {
          outFile << line << endl;
          if( modelToSelect >= 0 )
            bOutput = false;
        }
        terminatedChains.clear();
      } else {
        if( bHeader )
          outFile << line << endl;
      }
    }
  }

  return EXIT_SUCCESS;
}
