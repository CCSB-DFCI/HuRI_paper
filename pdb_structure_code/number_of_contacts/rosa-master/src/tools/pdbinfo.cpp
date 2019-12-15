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
  { "tab",        no_argument,       0, 't' },  // Output in tabular form
  { "fasta-seqs", no_argument,       0, 'F' },  // Output the sequences corresponding to the chains in fasta format
  { "non-std",    no_argument,       0, 'X' },  // Translates also non standard aminoacids
  { "seqs",       no_argument,       0, 'S' },  // Output the sequences with the corresponding residue numbers
  //{ "sstruct",    no_argument,       0, 'T' },  // Together with -S outputs the secondary structure elements
  { "sel",        required_argument, 0, 's' },  // Performs a selection before returning informations
  { "id",         required_argument, 0, 'i' },  // Sets the Pdb Id
  { "min-length", required_argument, 0, 'm' },  // Only outputs info for chains long at least min-length
  { "firstchain", no_argument,       0, 'f' },  // Output only the first chain when chainID repeats
  { "max-dist",   no_argument,       0, 'M' },  // Returns the maximum atom-atom distance in the PDB
  { "only-cas",   no_argument,       0, 'O' },  // Checks if the file only contains CA atoms and prints a message
  { "check",      no_argument,       0, 'C' },  // Performs different integrity checks on the PDB file
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] pdbFile [pdbfile [pdbfile [...]]" << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -t       --tab              prints output in tab form" << endl;
  cout << "   -s SEL   --sel=SEL          select only a portion of the PDB" << endl;
  cout << "   -i ID    --id=ID            sets the PDB id to ID" << endl;
  cout << "   -m LEN   --min-length=LEN   only outputs info for chains long at least LEN" << endl;
  cout << "   -f       --firstchain       outputs only the first of repeating chains" << endl;
  cout << "   -F       --fasta-seqs       outputs the sequences corresponding to the chains in fasta format" << endl;
  cout << "   -X       --non-std          works only with the -F option. Translates also non standard aminoacids" << endl;
  cout << "                               (normally they are tranlsated as 'X')" << endl;
  cout << "   -S       --seqs             outputs the residue numbers and the corresponding residue name" << endl;
  //cout << "   -T       --sstruct          with the -S option, outputs the secondary structure assignments" << endl;
  cout << "   -M       --max-dist         returns the maximum atom-atom distance in the PDB" << endl;
  cout << "   -O       --only-cas         checks if the file only contains CA atoms and prints a message" << endl;
  cout << "   -C       --check            performs different integrity checks on the PDB file" << endl;
  cout << "   -v       --version          prints the version number" << endl << endl;
}


int main( int argc, char **argv )
{
  string programName( argv[0] );

  DBG_BLOCK_SUPPRESS( programName+"::main" );

  bool bTab          = false;
  bool bFirstChain   = false;
  bool bSelect       = false;
  bool bFastaSeqs    = false;
  bool bTranslNonStd = false;
  bool bSeqs         = false;
  //bool bSecStruct    = false;
  bool bSetPdbId     = false;
  bool bMinLength    = false;
  bool bMaxDist      = false;
  bool bOnlyCAs      = false;
  bool bCheck        = false;

  string selString;
  string newPdbId;

  long minLength = 0;


  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":s:i:m:vSFXMOCfth", long_options, &option_index );

    if( ch == -1 ) // Finished with the options
      break;

    // Options parsing
    switch( ch ) {
      case 'v':
        printRosaVersion( programName );
        return EXIT_SUCCESS;
        break;
      case 't':
        bTab = true;
        break;
      case 'F':
        bFastaSeqs = true;
        break;
      case 'X':
        bTranslNonStd = true;
        break;
      //case 'T':
      //  bSecStruct = true;
      //  break;
      case 'S':
        bSeqs = true;
        break;
      case 'M':
        bMaxDist = true;
        break;
      case 'O':
        bOnlyCAs = true;
        break;
      case 'C':
        bCheck = true;
        break;
      case 's':
        bSelect = true;
        selString = string( optarg );
        break;
      case 'i':
        bSetPdbId = true;
        newPdbId = string( optarg );
        if( newPdbId.size() != 4 )
          errorExit( EXIT_FAILURE, "Please specify a valid PDB id with the '-i' option (4 characters code)" );
        break;
      case 'm':
        bMinLength = true;
        try {
          minLength = stol( string( optarg ) );
        } catch( exception &e ) {
          errorExit( EXIT_FAILURE, "Please specify a valid integer number with the '-m' option" );
        }
        if( minLength < 0 )
          errorExit( EXIT_FAILURE, "Please specify a positive integer number with the '-m' option" );
        break;
      case 'f':
        bFirstChain = true;
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

  if( argc < 1 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }

  if( bSetPdbId && argc > 1 ) {
    errorExit( EXIT_FAILURE, "Only one pdb file can be processed with the '-i' option" );
  }

  if( bTranslNonStd && !(bFastaSeqs || bSeqs) ) {
    warning( "'-X' option given without the '-F' option will be ignored." );
  }

  if( bCheck && bOnlyCAs ) {
    warning( "'-O' option given with the '-C' option will be ignored." );
  }

  //if( bSecStruct && !bSeqs ) ) {
  //  warning( "'-T' option given without the '-S' option will be ignored." );
  //}

  DBG_MSG( "Mark1" );

  for( int nCase = 0; nCase < argc; ++nCase ) {
    string pdbFilename = trim(string( argv[nCase] ));
    /*
    int pdbFilenameLen = pdbFilename.size();

    if( pdbFilenameLen <= 0 ) continue;
    if( pdbFilename[pdbFilenameLen-1] == '/' ) {
      cerr << "WARNING! Directory name skipped: " << pdbFilename << endl;
      continue;
    }
    if( pdbFilenameLen > 1 ) {
      if( pdbFilename.substr( pdbFilenameLen-2 ) == "/." ) {
        cerr << "WARNING! Directory name skipped: " << pdbFilename << endl;
        continue;
      }
      if( pdbFilenameLen > 2 ) {
        if( pdbFilename.substr( pdbFilenameLen-3 ) == "/.." ) {
          cerr << "WARNING! Directory name skipped: " << pdbFilename << endl;
          continue;
        }
      }
    }
    */
    if( !bTab && !bFastaSeqs ) {
      cout << "---------------------------------------------------------------" << endl;
      cout << "Information for file: " << pdbFilename << endl;
    }
    DBG_MSG( "Mark2" );

    StructurePtr s( new Structure(pdbFilename, true) );
    DBG_MSG( "Mark3" );

    if( bSetPdbId )
      s->setPdbID( newPdbId );

    if( bSelect ) {
      StructurePtr selStruct = s->select( selString );
      s.swap( selStruct );
    }

    if( bMaxDist ) {
      double maxDist = 0.0;
      for( std::size_t mi = (std::size_t)0; mi < s->numModels(); ++mi ) {
        std::vector<Point> v;
        s->getAtomPositions( v, mi );
        std::size_t vSize = v.size();
        for( std::size_t ai = 0UL; ai < vSize; ++ai )
          for( std::size_t aj = 0UL; aj < vSize; ++aj ) {
            double lDist = pointDistance( v[ai], v[aj] );
            if( lDist > maxDist )
              maxDist = lDist;
          }
      }
      cout << "PDB ID: " << s->getPdbID() << endl;
      cout << "Max atom-atom distance (Ang): " << fixed << setprecision(3) << maxDist << endl;
    } else {
      if( !(bFastaSeqs || bSeqs)) {
        if( !bTab ) {
          set<string> repChains;
          cout << "PDB ID: " << s->getPdbID() << endl;
          cout << "Number of models: " << s->numModels() << endl;
          for( std::size_t mi = (std::size_t)0; mi < s->numModels(); ++mi ) {
            cout << "**** MODEL " << s->getModelId(mi) << " ****" << endl;

            cout << "Number of chains: " << s->numChains(mi) << endl;
            for( unsigned long i = 0; i < s->numChains(mi); ++i ) {
              const Structure::Chain &c = s->getChain(i,mi);
              string chainID = c.getChainID() == " " ? "_" : c.getChainID();
              
              bool lChainHasOnlyCA = false; // Only CAs
              bool lChainHasOnlyMC = false; // Only main chain
              bool lChainHasOnlyAL = false; // Only alternative locations
              bool lChainHasOnlyUK = false; // Only unknown residues
              if( bOnlyCAs || bCheck ) {
                // Here I check if at least 20% of the residues are only CAs
                unsigned long numResWithOnlyCA = 0UL;
                unsigned long numResWithOnlyMC = 0UL;
                unsigned long numResWithOnlyAL = 0UL;
                unsigned long numUnknownRes    = 0UL;
                for( long r = 0L; r < c.numResidues(); ++r ) {
                  AtomVectCitPair avp = c.getResExtr(r);
                  int numAtomsInRes = avp.end - avp.beg;
                  if( c.hasResCA(r) && (numAtomsInRes == 1) )
                    ++numResWithOnlyCA;
                  if( bCheck ) {
                    if( c.hasOnlyMainChain(r) )
                      ++numResWithOnlyMC;
                    if( c.hasOnlyAltLoc(r) )
                      ++numResWithOnlyAL;
                    if( avp.beg->getResName() == "UNK" )
                      ++numUnknownRes;
                  }
                }
                if( float(numResWithOnlyCA) / float(c.numResidues()) >= 0.2 )
                  lChainHasOnlyCA = true;
                if( bCheck ) {
                  if( float(numResWithOnlyCA) / float(c.numResidues()) >= 0.5 )
                    lChainHasOnlyCA = true;
                  if( float(numResWithOnlyMC) / float(c.numResidues()) >= 0.5 )
                    lChainHasOnlyMC = true;
                  if( float(numResWithOnlyAL) / float(c.numResidues()) >= 0.5 )
                    lChainHasOnlyAL = true;
                  if( float(numUnknownRes) / float(c.numResidues()) >= 0.5 )
                    lChainHasOnlyUK = true;
                }
              }

              bool bOutput = true;

              if( bFirstChain )
                if( repChains.find( chainID ) != repChains.end() )
                  bOutput = false;

              if( bMinLength )
                if( c.numResidues() < minLength )
                  bOutput = false;

              if( bOutput ) {
                cout << "Chain " << chainID << " ("
                     << c.numResidues() << " res., "
                     << c.numAtoms() << " atoms, ";
                if( c.getChainType() == Structure::Chain::ProteinType )
                  cout << "protein): " << c.getDescription();
                else if( c.getChainType() == Structure::Chain::ProteinType )
                  cout << "nucleotides): " << c.getDescription();
                else
                  cout << "unknown)";
                cout << endl;
                if( (bOnlyCAs || bCheck) && lChainHasOnlyCA )
                  cout << "Chain " << chainID << " has only CA atoms" << endl;
                if( bCheck && lChainHasOnlyMC )
                  cout << "Chain " << chainID << " has only main chain atoms" << endl;
                if( bCheck && lChainHasOnlyAL )
                  cout << "Chain " << chainID << " has only atoms with alternative location set" << endl;
                if( bCheck && lChainHasOnlyUK )
                  cout << "Chain " << chainID << " has only unknown residues" << endl;
              }
            }
          }
        } else {
          set<string> repChains;
          for( unsigned long i = 0; i < s->numChains(); ++i ) {
            const Structure::Chain &c = s->getChain(i);
            const Structure::Chain::ChainType cType = c.getChainType();
            string chainID = c.getChainID() == " " ? "_" : c.getChainID();

            bool lChainHasOnlyCA = false; // Only CAs
            bool lChainHasOnlyMC = false; // Only main chain
            bool lChainHasOnlyAL = false; // Only alternative locations
            bool lChainHasOnlyUK = false; // Only unknown residues
            if( bOnlyCAs || bCheck ) {
              // Here I check if at least 20% of the residues are only CAs
              unsigned long numResWithOnlyCA = 0UL;
              unsigned long numResWithOnlyMC = 0UL;
              unsigned long numResWithOnlyAL = 0UL;
              unsigned long numUnknownRes    = 0UL;
              for( long r = 0L; r < c.numResidues(); ++r ) {
                AtomVectCitPair avp = c.getResExtr(r);
                int numAtomsInRes = avp.end - avp.beg;
                if( c.hasResCA(r) && (numAtomsInRes == 1) )
                  ++numResWithOnlyCA;
                if( bCheck ) {
                  if( c.hasOnlyMainChain(r) )
                    ++numResWithOnlyMC;
                  if( c.hasOnlyAltLoc(r) )
                    ++numResWithOnlyAL;
                  if( avp.beg->getResName() == "UNK" )
                    ++numUnknownRes;
                }
              }
              if( float(numResWithOnlyCA) / float(c.numResidues()) >= 0.2 )
                lChainHasOnlyCA = true;
              if( bCheck ) {
                if( float(numResWithOnlyCA) / float(c.numResidues()) >= 0.5 )
                  lChainHasOnlyCA = true;
                if( float(numResWithOnlyMC) / float(c.numResidues()) >= 0.5 )
                  lChainHasOnlyMC = true;
                if( float(numResWithOnlyAL) / float(c.numResidues()) >= 0.5 )
                  lChainHasOnlyAL = true;
                if( float(numUnknownRes) / float(c.numResidues()) >= 0.5 )
                  lChainHasOnlyUK = true;
              }
            }

            bool bOutput = true;

            if( bFirstChain )
              if( repChains.find( chainID ) != repChains.end() )
                bOutput = false;

            if( bMinLength )
              if( c.numResidues() < minLength )
                bOutput = false;

            if( bOutput ) {
              cout << s->getPdbID() << "\t" << chainID << "\t";
              if( cType == Structure::Chain::ProteinType )
                cout << "protein";
              else if( cType == Structure::Chain::ProteinType )
                cout << "nucleotides";
              else
                cout << "unknown";
              cout << "\t" << c.numResidues();
              cout << "\t" << c.numAtoms();
              cout << "\t";
              if( c.numResidues() > 0 )
                cout << c.getResLabel(0);
              else
                cout << "-";

              cout << "\t";
              if( c.numResidues() > 0 )
                cout << c.getResLabel(c.numResidues()-1);
              else
                cout << "-";
              if( bOnlyCAs || bCheck ) {
                cout << "\t";
                if( lChainHasOnlyCA )
                  cout << "ONLY_CAs";
              }
              if( bCheck ) {
                cout << "\t";
                if( lChainHasOnlyMC )
                  cout << "ONLY_MC";
                cout << "\t";
                if( lChainHasOnlyAL )
                  cout << "ONLY_AL";
                cout << "\t";
                if( lChainHasOnlyUK )
                  cout << "ONLY_UNK";
              }
              cout << endl;

              if( cType == Structure::Chain::ProteinType )
                repChains.insert( chainID );
            }
          }
        }
      }

      if( bFastaSeqs ) {
        set<string> repeatingChains;
        for( unsigned long i = 0; i < s->numChains(); ++i ) {
          const Structure::Chain &c = s->getChain(i);
          const Structure::Chain::ChainType cType = c.getChainType();
          string chainID = c.getChainID();
          if( chainID[0] < (int)(' ') )
            errorExit( EXIT_FAILURE, "Invalid character for chain "
                       "identifier 0x"+htos(chainID[0],2,'0') );
          //cout << (int)chainID[0] << "#" << chainID << "#" << endl;
          if( chainID == " " ) chainID = "_";
          bool bOutput = true;

          if( cType == Structure::Chain::ProteinType ) {
            if( bFirstChain )
              if( repeatingChains.find( chainID ) != repeatingChains.end() ) {
                cerr << "WARNING! Chain " << chainID << " is repeated." << endl;
                bOutput = false;
              }

            if( bMinLength )
              if( c.numResidues() < minLength )
                bOutput = false;

            if( bOutput ) {
              cout << ">" << s->getPdbID() << "_" << chainID << endl;
              DBG_VDUMP(bTranslNonStd);
              printInLines( c.getSequence(bTranslNonStd), 60 );
              cout << endl;
              repeatingChains.insert( chainID );
            }
          }
        }
      }

      if( bSeqs ) {
        set<string> repeatingChains;
        for( unsigned long i = 0; i < s->numChains(); ++i ) {
          const Structure::Chain &c = s->getChain(i);
          const Structure::Chain::ChainType cType = c.getChainType();
          string chainID = c.getChainID();
          if( chainID[0] < (int)(' ') )
            errorExit( EXIT_FAILURE, "Invalid character for chain "
                       "identifier 0x"+htos(chainID[0],2,'0') );
          //cout << (int)chainID[0] << "#" << chainID << "#" << endl;
          if( chainID == " " ) chainID = "_";
          bool bOutput = true;

          if( cType == Structure::Chain::ProteinType ) {
            if( bFirstChain )
              if( repeatingChains.find( chainID ) != repeatingChains.end() ) {
                cerr << "WARNING! Chain " << chainID << " is repeated." << endl;
                bOutput = false;
              }

            if( bMinLength )
              if( c.numResidues() < minLength )
                bOutput = false;

            if( bOutput ) {
              cout << ">" << s->getPdbID() << "_" << chainID << endl;
              //printInLines( c.getSequence(), 60 );
              for( long r = 0L; r < c.numResidues(); ++r ) {
                AtomVectCitPair avp = c.getResExtr(r);
                cout << avp.beg->getResNum() << "\t" << avp.beg->getInsCode()
                     << "\t" << avp.beg->getResName() << "\t" << res3to1(avp.beg->getResName(),bTranslNonStd);
                //if( bSecStruct )
                //  cout << avp.beg->getResNum()
                cout << endl;
              }
              repeatingChains.insert( chainID );
            }
          }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}
