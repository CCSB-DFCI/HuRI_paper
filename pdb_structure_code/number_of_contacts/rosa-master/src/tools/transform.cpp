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
#include <rosa/error.h>
#include <rosa/dbg.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace rosa;
using namespace std;

// Options
static struct option long_options[] = {
  { "version",    no_argument,       0, 'v' },  // Outputs only the version number
  { "verbose",    no_argument,       0, 'V' },  // Outputs only the version number
  { "output",     required_argument, 0, 'o' },  // Output filename
  { "rot-angles", required_argument, 0, 'R' },  // Eulerian angles for the rotation
  { "trans-vect", required_argument, 0, 'T' },  // Cartesian coords for the translation
  { "rot-sys",    required_argument, 0, 'f' },  // Rotation reference system
  { "rand-rot",   no_argument,       0, 'r' },  // Random rotation
  { "rand-trans", no_argument,       0, 't' },  // Random translation
  { "invert",     no_argument,       0, 'i' },  // Performs the translation and then the rotation 
  { "baricenter", no_argument,       0, 'b' },  // Transformations with respect to the baricenter
  { "origin",     no_argument,       0, 'O' },  // Brings the structure to the origin before applying transformations
  { "help",       no_argument,       0, 'h' },  // Outputs a short guide  
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] inputPdbFile" << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -v       --version          prints the version number" << endl;
  cout << "   -V       --verbose          prints informations about tranformations" << endl;
  cout << "   -o FILE  --output=FILE      set the output filename" << endl;
  cout << "   -b       --baricenter       performs rotation w.r.t. the baricenter" << endl;
  cout << "   -f zxz   --rot-ref zxz      rotation reference system (xyz, zxz)" << endl;
  cout << "   -r       --rand-rot         performs a random rotation" << endl;
  cout << "   -t       --rand-trans       performs a random translation" << endl;
  cout << "   -O       --origin           moves the struct. to origin before transf." << endl;
  cout << "   -i       --invert           performs the translation and then the rotation" << endl;
  cout << "   -R a,b,c --rot-angles=a,b,c angles for rotation (radiants)" << endl;
  cout << "   -T x,y,z --trans-vect=x,y,z vector for translation" << endl << endl;
}

 
int main( int argc, char **argv )
{
  DBG_BLOCK_SUPPRESS( "transform::main" );
  
  string programName( argv[0] );
  
  string outputFilename;
  bool   bRandomRotation = false;
  bool   bRandomTranslation = false;
  bool   bBaricenter = false;
  bool   bOrigin     = false;
  bool   bRotAngles  = false;
  bool   bInvert     = false;
  bool   bTransVect  = false;
  bool   bVerbose    = false;
  RotRefSystem rotRefSys = ROT_XYZ;
  
  double aAng = 0.0, bAng = 0.0, cAng = 0.0;
  double xVect = 0.0, yVect = 0.0, zVect = 0.0;
  
  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":o:R:T:f:vViObrth", long_options, &option_index );

    if( ch == -1 ) // Finished with the options
      break;
    
    // Options parsing
    switch( ch ) {
      case 'v':
        printRosaVersion( programName );
        return EXIT_SUCCESS;
        break;
      case 'V':
        bVerbose = true;
        break;          
      case 'o':
        outputFilename = string( optarg );
        break;
      case 'f':
        if( string( optarg ) == "zxz" ) {
          rotRefSys = ROT_ZXZ;
        } else if( string( optarg ) != "xyz" )
          errorExit( EXIT_FAILURE, "Invalid rotation system specification, use 'xyz' or 'zxz'" );
        break;
      case 'h':
        printRosaVersion( programName );
        usage( programName );
        return EXIT_SUCCESS;
        break;
      case 'b':
        bBaricenter = true;
        break;
      case 'r':
        bRandomRotation = true;
        break;
      case 't':
        bRandomTranslation = true;
        break;
      case 'O':
        bOrigin = true;
        break;
      case 'i':
        bInvert = true;
        break;
      case 'R': {
          bRotAngles = true;
          vector<string> angles;
          split( optarg, angles, "," );
          if( angles.size() != 3 )
            errorExit( EXIT_FAILURE, "Invalid angles specification, type '"
                       +programName+" -h' for help" );
          try {
            aAng = stod( angles[0] );
            bAng = stod( angles[1] );
            cAng = stod( angles[2] );
            DBG_VDUMP( aAng );
            DBG_VDUMP( bAng );
            DBG_VDUMP( cAng );
          } catch( exception &e ) {
            errorExit( EXIT_FAILURE, "Invalid angles specification:"+
                       string(e.what())+", type '"+programName+" -h' for help" );
          }
        }
        break;
      case 'T': {
          bTransVect = true;
          vector<string> comps;
          split( optarg, comps, "," );
          if( comps.size() != 3 )
            errorExit( EXIT_FAILURE, "Invalid components specification, type '"
                       +programName+" -h' for help" );
          try {
            xVect = stod( comps[0] );
            yVect = stod( comps[1] );
            zVect = stod( comps[2] );
            DBG_VDUMP( xVect );
            DBG_VDUMP( yVect );
            DBG_VDUMP( zVect );
          } catch( exception &e ) {
            errorExit( EXIT_FAILURE, "Invalid components specification:"+
                       string(e.what())+", type '"+programName+" -h' for help" );
          }
        }
        break;
      case '?':
      default:
        printRosaVersion( programName );
        usage( programName );
        errorExit( EXIT_FAILURE, "Invalid command line argument" );
        break;
    }
  }

  if( bRandomRotation && bRotAngles )
    errorExit( EXIT_FAILURE, "Options '-r' and '-A' are mutually exclusive." );
  if( bRandomTranslation && bTransVect )
    errorExit( EXIT_FAILURE, "Options '-t' and '-V' are mutually exclusive." );
  if( bOrigin && bBaricenter )
    warning( "With options '-O' the option '-b' is not taken into consideration." );

  argc -= optind;
  argv += optind;

  DBG_VDUMP( argc );
  
  if( argc != 1 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }
  
  string pdbFilename( argv[0] );
  
  if( outputFilename.empty() )
    outputFilename = woExtFilename ( pdbFilename ) + "-tr." + extFilename( pdbFilename );
  
  Structure s;
  
  DBG_MSG( "Mark1" );
  if( bVerbose )
    cout << "Loading file " << pdbFilename << endl;
  s.loadFromPDBFile( pdbFilename );
  
  DBG_MSG( "Mark2" );
  
  if( bRandomRotation || bRandomTranslation )
    srand( time(0) );
  
  if( bOrigin ) {
    Point baricenter( s.getBaricenter() );
    if( bVerbose )
      cout << "Translating of " << -baricenter << " to origin." << endl;
    s.translate( -baricenter );
  }
  
  if( (bRandomTranslation || bTransVect) and bInvert ) {
    Coord xt = 0.0, yt = 0.0, zt = 0.0;
    if( bRandomTranslation ) {
      double rGyr = s.calcRadiusOfGyration();
      xt = rand01() * 2.0 * rGyr;
      yt = rand01() * 2.0 * rGyr;
      zt = rand01() * 2.0 * rGyr;
    } else if( bTransVect ) {
      xt = xVect;
      yt = yVect;
      zt = zVect;
    }
    if( bVerbose )
      cout << "Translating of " << Point( xt, yt, zt ) << endl;
    s.translate( Point( xt, yt, zt ) );
  }

  if( bRandomRotation || bRotAngles ) {
    Point baricenter;
    if( bBaricenter && !bOrigin ) {
      baricenter = s.getBaricenter();
      if( bVerbose )
        cout << "Translating of " << -baricenter << " to origin." << endl;
      s.translate( -baricenter );
    }
    DBG_MSG( "Mark2" );
    if( bRandomRotation ) {
      aAng = (rand01() * 2.0 - 1.0) * M_PI;
      bAng = (rand01() * 2.0 - 1.0) * M_PI;
      cAng = (rand01() * 2.0 - 1.0) * M_PI;
    }
    Mat2D<Coord> rotMat( rotMatrix( aAng, bAng, cAng, rotRefSys ) );
    if( bVerbose )
      cout << "Rotating of " << Point( aAng, bAng, cAng )
           << " radians corresponding to "
           << Point( aAng/M_PI*180.0, bAng/M_PI*180.0, cAng/M_PI*180.0 )
           << " degrees" << endl;
    s.rotate( rotMat );
    DBG_MSG( "Mark3" );
    if( bBaricenter && !bOrigin ) {
      if( bVerbose )
        cout << "Translating back of " << baricenter << " to initial position." << endl;
      s.translate( baricenter );
    }
  }

  if( (bRandomTranslation || bTransVect) and !bInvert ) {
    Coord xt = 0.0, yt = 0.0, zt = 0.0;
    if( bRandomTranslation ) {
      double rGyr = s.calcRadiusOfGyration();
      xt = rand01() * 2.0 * rGyr;
      yt = rand01() * 2.0 * rGyr;
      zt = rand01() * 2.0 * rGyr;
    } else if( bTransVect ) {
      xt = xVect;
      yt = yVect;
      zt = zVect;
    }
    if( bVerbose )
      cout << "Translating of " << Point( xt, yt, zt ) << endl;
    s.translate( Point( xt, yt, zt ) );
  }
  
  if( bVerbose )
    cout << "Saving file " << outputFilename << endl;
  s.saveToPDBFile( outputFilename );
  
  return EXIT_SUCCESS;
}
