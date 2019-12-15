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

using namespace rosa;
using namespace std;

const char rotChainIDs[] = {
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
  'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
  'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
  'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
  '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '!', '@', '#',
  '$', '%', '^', '&', '*', '(', ')', '-', '+', '=', '{', '}', '[',
  ']', ':', '"', ';', '\'','<', '>', ',', '.', '?', '/', '`', '~',
  };

const int numOfRotChainIDs = 62;

// Options
static struct option long_options[] = {
  { "version",    no_argument,       0, 'v' },  // Outputs only the version number
  { "flatten",    no_argument,       0, 'F' },  // flatten the models in the PDB file
  { "chain",      required_argument, 0, 'c' },  // Change the chain ID
  { "alt",        required_argument, 0, 'a' },  // Change the alt code
  { "offset",     required_argument, 0, 'f' },  // Sum an offset to the residue number
  { "renumber",   required_argument, 0, 'r' },  // Renumber residues
  { "output",     required_argument, 0, 'o' },  // Output filename
  { "help",       no_argument,       0, 'h' },  // Outputs a short guide
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] select_string pdbFile" << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h        --help             prints a short guide" << endl;
  cout << "   -F        --flatten          flatten the models" << endl;
  cout << "   -c CH_ID  --chain=CH_ID      set the chaind ID to CH_ID" << endl;
  cout << "   -a ALT    --alt=ALT          set the alt. loc. ID to ALT" << endl;
  cout << "   -f OFFSET --output=OFFSET    sum an offset to the residue number" << endl;
  cout << "   -r START  --renumber=START   renumber residues sequentially starting from START" << endl;
  cout << "   -o FILE   --output=FILE      set the output filename" << endl;
  cout << "   -v        --version          prints the version number" << endl << endl;
}


int main( int argc, char **argv )
{
  string programName( argv[0] );

  bool   bChain    = false;
  bool   bAlt      = false;
  bool   bOffset   = false;
  bool   bRenumber = false;
  bool   bFlatten  = false;

  long   offset = 0;
  long   start  = 0;

  string newAlt;
  string newChain;
  string outputFilename;

  DBG_BLOCK_SUPPRESS( programName+"::main" );

  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":o:c:f:r:a:vFh", long_options, &option_index );

    if( ch == -1 ) // Finished with the options
      break;

    // Options parsing
    switch( ch ) {
      case 'F':
        bFlatten = true;
        break;
      case 'v':
        printRosaVersion( programName );
        return EXIT_SUCCESS;
        break;
      case 'o':
        outputFilename = string( optarg );
        break;
      case 'c':
        bChain = true;
        newChain = string( optarg );
        break;
      case 'a':
        bAlt = true;
        newAlt = string( optarg );
        break;
      case 'f':
        bOffset = true;
        offset = stol( string(optarg) );
        break;
      case 'r':
        bRenumber = true;
        start = stol( string(optarg) );
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

  if( bFlatten && (bChain or bAlt or bOffset or bRenumber ) )
    errorExit( EXIT_FAILURE, "-F option can only go alone." );

  if( bOffset and bRenumber )
    errorExit( EXIT_FAILURE, "Cannot give -f and -r at the same time." );
  DBG_MSG( "Mark1" );

  string selectionString( argv[0] );
  string pdbFilename( argv[1] );

  if( outputFilename.empty() )
    outputFilename = baseFilename( pdbFilename ) + "-mod." +
                     extFilename( pdbFilename );

  DBG_MSG( "Mark2" );

  AtomSelectDriver selDriver( selectionString );

  if( !selDriver.good() )
    errorExit( EXIT_FAILURE, "Error while parsing selection string "+selectionString );

  ifstream inFile( pdbFilename.c_str() );

  if( !inFile.good() )
    errorExit( EXIT_FAILURE, "Error while opening file "+pdbFilename );

  ofstream outFile( outputFilename.c_str() );

  if( !outFile.good() )
    errorExit( EXIT_FAILURE, "Error while opening output file "+outputFilename );


  if( bFlatten ) {
    int chainIDindex = -1;
    map<string,string> chainIDremapping;
    bool bEndModelFound = false;
    while( !inFile.eof() ) {
      string line;
      DBG_MSG( "Reading next line" );
      if( getline( inFile, line ) ) {
        DBG_MSG( "Line read." );
        DBG_VDUMP( line.size() );

        // Pads the length of the line to reach the normal length of 80 characters
        line += string( max<int>( PDB_LINE_LENGTH - line.size(), 0 ), ' ' );

        DBG_VDUMP( line );
        string recordName( line.substr( 0, 6 ) );

        if( recordName == AtomRecordName or
            recordName == HetatmRecordName or
            recordName == TerRecordName ) {
          Atom at;
          at.readPDBAtomRecord( line );

          string curChainID = at.getChainID();
          string newChainID;
          map<string,string>::iterator cidIt = chainIDremapping.find(curChainID);
          if( cidIt != chainIDremapping.end() )
            newChainID = cidIt->second;
          else {
            if( chainIDindex < 0 )
              chainIDindex = 0;
            newChainID = rotChainIDs[chainIDindex];
            chainIDremapping.insert( pair<string,string>( curChainID, newChainID ) );
            cout << "MAP CHAIN '" << curChainID << "' TO CHAIN '" << newChainID << "'" << endl;
            chainIDindex++;
            if( chainIDindex >= numOfRotChainIDs )
              warning( string("Too many chains, restarting from ")+rotChainIDs[0] );
            chainIDindex %= numOfRotChainIDs;
          }
          DBG_VDUMP( chainIDindex )
          DBG_VDUMP( newChainID )
          DBG_VDUMP( curChainID )

          if( selDriver.evalOnAtom( at ) ) {
            at.setChainID( newChainID );
            outFile << at << endl;
          }
        } else if( recordName == ModelRecordName ) {
          if( chainIDindex != -1 and not bEndModelFound )
            errorExit( EXIT_FAILURE, "Wrong PDB format, MODEL record without "
                       "the corresponding ENDMDL." );
          //chainIDindex++;
          //if( chainIDindex >= numOfRotChainIDs )
          //  warning( string("Too many chains, restarting from ")+rotChainIDs[0] );
          //chainIDindex %= numOfRotChainIDs;
          chainIDremapping.clear();
          bEndModelFound = false;
        } else if( recordName == EndmdlRecordName ) {
          bEndModelFound = true;
        } else {
          outFile << line << endl;
        }
      }
    }
  } else {
    long currentResNum = start;
    string lastResLabel;
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

        if( recordName == AtomRecordName or
            recordName == HetatmRecordName or
            recordName == TerRecordName ) {
          Atom at;
          at.readPDBAtomRecord( line );

          if( selDriver.evalOnAtom( at ) ) {
            if( bChain )
              at.setChainID( newChain );
            if( bAlt )
              at.setAltLoc( newAlt );
            if( bOffset )
              at.setResNum( at.getResNum() + offset );
            if( bRenumber ) {
              string localResLabel = Atom::getResLabel( at.getResNum(), at.getInsCode() );
              if( localResLabel != lastResLabel ) {
                if( !lastResLabel.empty() ) ++currentResNum;
                lastResLabel = localResLabel;
              }
              at.setResNum( currentResNum );
              at.setInsCode( " " );
            }
          }
          outFile << at << endl;
        } else
          outFile << line << endl;
      }
    }
  }

  return EXIT_SUCCESS;
}
