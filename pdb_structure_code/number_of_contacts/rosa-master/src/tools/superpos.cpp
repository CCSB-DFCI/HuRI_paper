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
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace rosa;
using namespace std;

// Options
static struct option long_options[] = {
  { "version",    no_argument,         0, 'v' },  // Outputs only the version number
  { "verbose",    no_argument,         0, 'V' },  // Outputs only the version number
  { "rmsd",       no_argument,         0, 'r' },  // Only RMSD of the superposition
  { "all",        no_argument,         0, 'a' },  // Superpose on all the atoms
  { "common",     no_argument,         0, 'c' },  // Superpose on the common residues/atoms
  { "preserve",   no_argument,         0, 'P' },  // preserve ligands (after TER)
  { "selref",     required_argument,   0, 'R' },  // Performs a selection on the first PDB before superposing
  { "selmov",     required_argument,   0, 'M' },  // Performs a selection on the second PDB before superposing
  { "output",     required_argument,   0, 'o' },  // Output filename
  { "tr-output",  required_argument,   0, 't' },  // Transformation output filename
  { "help",       no_argument,         0, 'h' },  // Outputs a short guide
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] refPdbFile PdbFileToMove" << endl << endl;
  cout << "The superposition is normally done only on the CAs. If you want to superpose on" << endl;
  cout << "all the atoms use the '-a' option." << endl << endl;
  cout << "* Please use PDB files with only one chain. *" << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -r       --rmsd             only calculates the rmsd of the best superposition" << endl;
  cout << "   -v       --version          prints the version number" << endl;
  cout << "   -V       --verbose          prints informations about tranformations" << endl;
  cout << "   -a       --all              superpose on all the atoms" << endl;
  cout << "   -P       --preserve         preserve ligands (after TER)" << endl;
  cout << "   -c       --common           superpose on the common CAs (does not work with -a)" << endl;
  cout << "   -R SEL   --selref=SEL       select on part of the first PDB for superposition" << endl;
  cout << "   -M SEL   --selmov=SEL       select on part of the second PDB for superposition" << endl;
  cout << "   -o FILE  --output=FILE      set the output filename" << endl;
  cout << "   -t FILE  --tr-output=FILE   saves the transformation to an output file" << endl << endl;
}


int main( int argc, char **argv )
{
  DBG_BLOCK_SUPPRESS( "superpos::main" );

  string programName( argv[0] );

  string outputFilename;
  string trOutputFilename;
  bool bVerbose  = false;
  bool bAll      = false;
  bool bRmsd     = false;
  bool bSelRef   = false;
  bool bSelMov   = false;
  bool bCommon   = false;
  bool bPreserve = false;

  RotRefSystem rotRefSys = ROT_XYZ;

  string selStringRef;
  string selStringMov;

  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":o:t:R:M:vacPVhr", long_options, &option_index );

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
      case 'a':
        bAll = true;
        break;
      case 'c':
        bCommon = true;
        break;
      case 'P':
        bPreserve = true;
        break;
      case 'r':
        bRmsd = true;
        break;
      case 'R':
        bSelRef = true;
        selStringRef = string( optarg );
        break;
      case 'M':
        bSelMov = true;
        selStringMov = string( optarg );
        break;
      case 'o':
        outputFilename = string( optarg );
        break;
      case 't':
        trOutputFilename = string( optarg );
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

  if( bRmsd && !trOutputFilename.empty() )
    errorExit( EXIT_FAILURE, "Options '-r' and '-t' are mutually exclusive." );

  if( bCommon && bAll )
    errorExit( EXIT_FAILURE, "Options '-c' and '-a' are mutually exclusive." );

  if( bVerbose ) {
    if( bAll )
      cout << "Using all the atoms" << endl;
    else
      cout << "Using only CA atoms" << endl;
  }

  argc -= optind;
  argv += optind;

  DBG_VDUMP( argc );

  if( argc != 2 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }

  DBG_MSG( "Mark1" );

  string refPdbFilename( argv[0] );
  string toMovePdbFilename( argv[1] );

  if( outputFilename.empty() )
    outputFilename = woExtFilename ( toMovePdbFilename ) + "-sup." +
                     extFilename( toMovePdbFilename );

  DBG_MSG( "Mark2" );

  Structure ref, toMove;

  if( bVerbose )
    cout << "Loading reference file: " << refPdbFilename << endl;
  ref.loadFromPDBFile( refPdbFilename );
  if( bVerbose )
    cout << "Loading file to move:   " << toMovePdbFilename << endl;
  toMove.loadFromPDBFile( toMovePdbFilename, false, !bPreserve );

  DBG_MSG( "Mark3" );

  if( DBG_ON ) {
    DBG_MSG( "Before remove mult conformations" )
    printChainDescription( ref );
    printChainDescription( toMove );
  }

  ref.removeMultConformations();
  toMove.removeMultConformations();

  if( DBG_ON ) {
    DBG_MSG( "After remove mult conformations" )
    printChainDescription( ref );
    printChainDescription( toMove );
  }

  Mat2D<double> rotMat( identityMatrix<double>( 3 ) );
  Point cm1, cm2;
  double rmsd;
  long numAtoms;

  try {
    if( !bAll ) {
      DBG_VDUMP( ref.size() );
      DBG_VDUMP( ref.numModels() );

      StructurePtr selRef    = ref.select( "name CA" );
      StructurePtr selToMove = toMove.select( "name CA" );

      if( bSelRef ) {
        if( bVerbose )
          cout << "Selecting portion of structure " << refPdbFilename << "...";
        StructurePtr tmpRefStruct = selRef->select( selStringRef );
        selRef.swap( tmpRefStruct );
        if( bVerbose )
          cout << selRef->size() << " atoms selected." << endl;
        if( DBG_ON )
          printChainDescription( *selRef );
      }

      if( bSelMov ) {
        if( bVerbose )
          cout << "Selecting portion of structure " << toMovePdbFilename << "...";
        StructurePtr tmpToMoveStruct = selToMove->select( selStringMov );
        selToMove.swap( tmpToMoveStruct );
        if( bVerbose )
          cout << selToMove->size() << " atoms selected." << endl;
        if( DBG_ON )
          printChainDescription( *selToMove );
      }

      DBG_VDUMP( selRef->size() );
      DBG_VDUMP( selRef->numModels() );

      if( bCommon ) {
        reduceToAlignedAtoms( *selRef, *selToMove );
      }

      numAtoms = selRef->size();

      if( bRmsd )
        rmsd = calcRmsdSuperposition( *selRef, *selToMove );
      else
        rmsd = calcSuperposition( *selRef, *selToMove, rotMat, cm1, cm2 );

      int loopCounter = 0;
      while( !isnormal(rmsd) && loopCounter++ < 10 ) {
        double aRndAng = 0.0, bRndAng = 0.0, cRndAng = 0.0;
        srand( time(0) );
        aRndAng = (rand01() * 2.0 - 1.0) * M_PI;
        bRndAng = (rand01() * 2.0 - 1.0) * M_PI;
        cRndAng = (rand01() * 2.0 - 1.0) * M_PI;
        Mat2D<Coord> rndRotMat( rotMatrix( aRndAng, bRndAng, cRndAng, rotRefSys ) );
        selToMove->rotate( rndRotMat );
        if( bRmsd )
          rmsd = calcRmsdSuperposition( *selRef, *selToMove );
        else
          rmsd = calcSuperposition( *selRef, *selToMove, rotMat, cm1, cm2 );
      }
      DBG_MSG( "Mark5" );
    } else {
      numAtoms = ref.size();
      if( bRmsd )
        rmsd = calcRmsdSuperposition( ref, toMove );
      else
        rmsd = calcSuperposition( ref, toMove, rotMat, cm1, cm2 );

      int loopCounter = 0;
      while( !isnormal(rmsd) && loopCounter++ < 10 ) {
        double aRndAng = 0.0, bRndAng = 0.0, cRndAng = 0.0;
        srand( time(0) );
        aRndAng = (rand01() * 2.0 - 1.0) * M_PI;
        bRndAng = (rand01() * 2.0 - 1.0) * M_PI;
        cRndAng = (rand01() * 2.0 - 1.0) * M_PI;
        Mat2D<Coord> rndRotMat( rotMatrix( aRndAng, bRndAng, cRndAng, rotRefSys ) );
        toMove.rotate( rndRotMat );
        if( bRmsd )
          rmsd = calcRmsdSuperposition( ref, toMove );
        else
          rmsd = calcSuperposition( ref, toMove, rotMat, cm1, cm2 );
      }
    }

    if( !trOutputFilename.empty() ) {
      ofstream trOutFile( trOutputFilename.c_str() );

      trOutFile << rotMat[0][0] << "\t" << rotMat[0][1] << "\t" << rotMat[0][2] << endl;
      trOutFile << rotMat[1][0] << "\t" << rotMat[1][1] << "\t" << rotMat[1][2] << endl;
      trOutFile << rotMat[2][0] << "\t" << rotMat[2][1] << "\t" << rotMat[2][2] << endl;
      trOutFile << endl;
      trOutFile << cm1.x() << "\t" << cm1.y() << "\t" << cm1.z() << endl;
      trOutFile << cm2.x() << "\t" << cm2.y() << "\t" << cm2.z() << endl;
      double xAng, yAng, zAng;
      rotAngles( rotMat, xAng, yAng, zAng );
      trOutFile << endl;
      trOutFile << xAng << "\t" << yAng << "\t" << zAng << endl;
    }

    if( bVerbose or bRmsd )
      cout << "Rmsd after superposition " << fixed << setprecision(2) << rmsd << " Ang. (" << numAtoms << " atoms)" << endl;
  } catch( exception &e ) {
    errorExit( EXIT_FAILURE, e.what() );
  }

  if( !bRmsd ) {
    DBG_MSG( "Mark6" );
    Structure toTransf( toMove );
    DBG_MSG( "Mark7" );

    toTransf.translate( -cm2 );
    toTransf.rotate( rotMat );
    toTransf.translate( cm1 );

    DBG_MSG( "Mark8" );
    if( bVerbose )
      cout << "Saving file " << outputFilename << endl;
    toTransf.saveToPDBFile( outputFilename );
  }

  return EXIT_SUCCESS;
}
