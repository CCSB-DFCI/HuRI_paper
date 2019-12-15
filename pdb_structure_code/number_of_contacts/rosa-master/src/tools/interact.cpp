#include <rosa/structure.h>
#include <rosa/version.h>
#include <rosa/error.h>
#include <rosa/util.h>
#include <rosa/elems.h>
#include <rosa/dbg.h>
#include <getopt.h>
#include <iostream>
#include <cstdlib>

using namespace rosa;
using namespace std;

enum InteractionType {
  None = 0,
  Disulfide,
  HydrogenBond,
  VanDerWaals,
  Clash
};

// Options
static struct option long_options[] = {
  { "version",    no_argument,       0, 'v' },  // Outputs only the version number
  { "verbose",    no_argument,       0, 'V' },  // Verbose mode
  { "table",      no_argument,       0, 't' },  // Tabular form
  { "res-index",  no_argument,       0, 'i' },  // Uses residue index
  { "expand",     no_argument,       0, 'E' },  // Expand models
  { "help",       no_argument,       0, 'h' },  // Outputs a short guide  
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] pdbFile chain_id1 chain_id2" << endl << endl;
  cout << "Returns the number of interactions between the two chains." << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -t       --table            prints interactions in a tabular form" << endl;
  cout << "   -E       --expand           expands models" << endl;
  cout << "   -i       --res-index        uses residue index instead of label" << endl;
  cout << "   -v       --version          prints the version number" << endl;
  cout << "   -V       --verbose          activates verbose mode" << endl << endl;
}


int main( int argc, char *argv[] )
{
  DBG_BLOCK_SUPPRESS( "interact::main" );
  string programName( argv[0] );
  bool bVerbose  = false;
  bool bTable    = false;
  bool bExpand   = false;
  bool bResIndex = false;
  
  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":vtiEVh", long_options, &option_index );

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
      case 'E':
        bExpand = true;
        break;
      case 't':
        bTable = true;
        break;
      case 'i':
        bResIndex = true;
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

  if( argc != 3 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }

  string pdbFilename( argv[0] );
  string chainID1( argv[1] );
  string chainID2( argv[2] );
  
  if( chainID1 == "_" ) chainID1 = string(" ");
  if( chainID2 == "_" ) chainID2 = string(" ");
  
  if( bVerbose )
    cout << "Loading file " << pdbFilename << "..." << flush;
  Structure pdbFile;
  pdbFile.loadFromPDBFile( pdbFilename, true );
  if( DBG_ON ) {
    DBG_VDUMP(pdbFile.numModels());
    for( int imod = 0; imod < (int)pdbFile.numModels(); ++imod ) {
      DBG_VDUMP(imod);
      DBG_VDUMP(pdbFile.numChains(imod));
    }
  }
  if( bVerbose )
    cout << "done." << endl;
  
  if( bVerbose )
    cout << "Selecting chain \"" << chainID1 << "\"..." << flush;
  StructurePtr chain1( pdbFile.select( "chain \""+chainID1+"\"" ) );
  Structure &s1 = *(chain1.get());
  if( bVerbose ) {
    cout << "done." << endl;
    cout << "Selecting chain \"" << chainID2 << "\"..." << flush;
  }
  StructurePtr chain2( pdbFile.select( "chain \""+chainID2+"\"" ) );
  Structure &s2 = *(chain2.get());
  if( bVerbose )
    cout << "done." << endl;
  
  
  unsigned long numModels          = pdbFile.numModels();
  unsigned long reprModel          = pdbFile.getReprModel();

  unsigned long startModel         = (bExpand?0:reprModel);
  unsigned long endModel           = (bExpand?numModels-1:reprModel);

  set< pair<string,string> > alreadyDone;
  
  DBG_VDUMP(startModel);
  DBG_VDUMP(endModel);
  for( int imod = startModel; imod <= (int)endModel; ++imod ) {
    DBG_VDUMP(imod);
    
    for( int jmod = startModel; jmod <= (int)endModel; ++jmod ) {
      DBG_VDUMP(jmod);
      DBG_VDUMP(s1.numChains(imod));
      
      // We process the two structures residue by residue
      for( unsigned long ci = 0; ci < s1.numChains(imod); ++ci ) {
        const Structure::Chain &c1 = s1.getChain(ci,imod);
        DBG_VDUMP(s2.numChains(jmod));

        for( unsigned long cj = 0; cj < s2.numChains(jmod); ++cj ) {
          const Structure::Chain &c2 = s2.getChain(cj,jmod);
          DBG_MSG( "Mark3" );
          
          string doing1( c1.getChainID()+":"+ltos(imod) );
          string doing2( c2.getChainID()+":"+ltos(jmod) );

          if( alreadyDone.find( pair<string,string>(doing1, doing2) ) != alreadyDone.end() )
            continue;
          if( alreadyDone.find( pair<string,string>(doing2, doing1) ) != alreadyDone.end() )
            continue;
          
          DBG_MSG( "Mark4" );
          alreadyDone.insert( pair<string,string>(doing1, doing2) );
          
          if( (jmod == imod) && (c1.getChainID() == c2.getChainID()) )
            continue;
          
          if( bExpand ) {
            cout << endl << "Comparing " << c1.getChainID() << ":" << imod
                 << " and " << c2.getChainID() << ":" << jmod << endl;
          }
          
          unsigned long numClashes         = 0;
          unsigned long numDisulfides      = 0;
          unsigned long numVdWInteractions = 0;
          unsigned long numHbonds          = 0;

          if( bVerbose )
            cout << endl << "List of interactions:" << endl;
            
          if( bVerbose or bTable )
            cout << "    Residue     Atom   Residue     Atom  Type     Dist" << endl;
          // For every residue of the first chain
          for( long ri = 0; ri < c1.numResidues(); ++ri ) {
            
            // First we want to be sure that the two residues have a CA
            if( !(c1.hasResCA( ri )) ) continue;
            
            AtomVectCitPair vp1 = c1.getResExtr( ri );

            // For every residue of the second chain
            for( long rj = 0; rj < c2.numResidues(); ++rj ) {
              
              // First we want to be sure that the two residues have a CA
              if( !(c2.hasResCA( rj )) ) continue;

              AtomVectCitPair vp2 = c2.getResExtr( rj );
          
              bool bSkipToNextResPair = false;
              InteractionType intType = None;
              AtomVectCit selIt, selJt;
              double selDist = 0.0;
              DBG_MSG( "Residue pair: ("+ltos(ri)+c1.getChainID()+", "+ltos(rj)+c2.getChainID()+")" );
              
              // First we check for disulfide bridges interactions
              for( AtomVectCit it = vp1.beg; it != vp1.end and !bSkipToNextResPair and intType == None; ++it ) {
                for( AtomVectCit jt = vp2.beg; jt != vp2.end and intType == None; ++jt ) {
                  double lDist = pointDistance( it->getPos(), jt->getPos() );
                  if( lDist > 30.0 ) {
                    bSkipToNextResPair = true;
                    break;
                  }
                  //DBG_MSG( "distance("+it->identifier()+", "+jt->identifier()+") = "+dtos(lDist,-1,2) );

                  int iElem = it->getElemNo();
                  int jElem = jt->getElemNo();
                  double covMaxDist = Elements[iElem].crad +
                                      Elements[jElem].crad + 0.5;
                
                  
                  if( lDist <= covMaxDist ) {
                    // The two atoms are covalently bonded
                    DBG_MSG( "Residues "+Atom::getResLabel(it->getResNum(),it->getInsCode())+\
                             " and "+Atom::getResLabel(jt->getResNum(),jt->getInsCode())+\
                             " are potentially covalently bonded via atoms pair ("+\
                             it->identifier()+", "+jt->identifier()+")" );
                    if( it->getResName() == "CYS" and jt->getResName() == "CYS" and
                        it->getName() == "SG" and jt->getName() == "SG" ) {
                      numDisulfides++;
                      intType = Disulfide;
                    } else {
                      numClashes++;
                      intType = Clash;
                    }
                    selIt = it; selJt = jt;
                    selDist = lDist;
                    break;
                  }
                }
              }
              
              if( intType == None and !bSkipToNextResPair ) {
                for( AtomVectCit it = vp1.beg; it != vp1.end and intType == None; ++it ) {
                  for( AtomVectCit jt = vp2.beg; jt != vp2.end and intType == None; ++jt ) {
                    double lDist = pointDistance( it->getPos(), jt->getPos() );
                    //DBG_MSG( "distance("+it->identifier()+", "+jt->identifier()+") = "+dtos(lDist,-1,2) );

                    int iElem = it->getElemNo();
                    int jElem = jt->getElemNo();

                    // Is there any hydrogen bond?
                    double maxHBDist = -1.0;
                    
                    if( iElem == Elem_O ) {
                      if( jElem == Elem_O ) maxHBDist = 2.8;
                      else if( jElem == Elem_N ) maxHBDist = 3.1;
                    } else if( iElem == Elem_N ) {
                      if( jElem == Elem_O ) maxHBDist = 3.1;
                      else if( jElem == Elem_N )  maxHBDist = 3.2;
                    }
                    
                    if( maxHBDist > 0.0 ) {
                      if( lDist <= maxHBDist ) {
                        // The two atoms are hydrogen bonded
                        DBG_MSG( "Residues "+Atom::getResLabel(it->getResNum(),it->getInsCode())+\
                                 " and "+Atom::getResLabel(jt->getResNum(),jt->getInsCode())+\
                                 " are potentially hidrogen bonded via atoms pair ("+\
                                 it->identifier()+", "+jt->identifier()+")" );
                        numHbonds++;
                        intType = HydrogenBond;
                        selIt = it; selJt = jt;
                        selDist = lDist;
                        break;
                      }
                    }
                  }
                }
              }
                
              if( intType == None and !bSkipToNextResPair ) {
                for( AtomVectCit it = vp1.beg; it != vp1.end and intType == None; ++it ) {
                  for( AtomVectCit jt = vp2.beg; jt != vp2.end and intType == None; ++jt ) {
                    double lDist = pointDistance( it->getPos(), jt->getPos() );
                    //DBG_MSG( "distance("+it->identifier()+", "+jt->identifier()+") = "+dtos(lDist,-1,2) );

                    int iElem = it->getElemNo();
                    int jElem = jt->getElemNo();

                    // Check for Van Der Waals interaction
                    double vdw1 = Elements[iElem].vdw;
                    double vdw2 = Elements[jElem].vdw;

                    if( vdw1 > 0.0 && vdw2 > 0.0 ) {
                      if( lDist <= vdw1+vdw2 ) {
                        // The two atoms have VDW interactions
                        DBG_MSG( "Residues "+Atom::getResLabel(it->getResNum(),it->getInsCode())+\
                                 " and "+Atom::getResLabel(jt->getResNum(),jt->getInsCode())+\
                                 " have a potential VdW contact via atoms pair ("+\
                                 it->identifier()+", "+jt->identifier()+")" );
                        numVdWInteractions++;
                        intType = VanDerWaals;
                        selIt = it; selJt = jt;
                        selDist = lDist;
                        bSkipToNextResPair = true;
                        break;
                      }
                    }
                  }
                }
              }
              
              if( intType != None and (bVerbose or bTable)  ) {
                cout << setw(5) << selIt->getResName();
                if( bResIndex )
                  cout << setw(8) << ri;
                else
                  cout << setw(8) << trim(Atom::getResLabel(selIt->getResNum(), selIt->getInsCode()));
                cout << setw(6) << selIt->getName();
                cout << setw(5) << selJt->getResName();
                if( bResIndex )
                  cout << setw(8) << rj;
                else
                  cout << setw(8) << trim(Atom::getResLabel(selJt->getResNum(), selJt->getInsCode()));
                cout << setw(6) << selJt->getName();
                switch( intType ) {
                  case Clash:        cout << "   clash"; break;
                  case Disulfide:    cout << "   disul"; break;
                  case HydrogenBond: cout << "   hb   "; break;
                  case VanDerWaals:  cout << "   vdw  "; break;
                }
                cout << setw(8) << fixed << setprecision(3) << selDist;
                cout << endl;
              }
            }
          }
          unsigned long numInteractions = numDisulfides + numVdWInteractions + numHbonds;

          cout << endl;
          cout << "Number of potential clashes:    " << setw(8) << numClashes << endl;
          cout << "Number of disulfide bridges:    " << setw(8) << numDisulfides << endl;
          cout << "Number of hydrogen bonds:       " << setw(8) << numHbonds << endl;
          cout << "Number of vdw interactions:     " << setw(8) << numVdWInteractions << endl;
          cout << "Total number of interactions:   " << setw(8) << numInteractions << endl;
        }
      }
    }
  }
  
  return EXIT_SUCCESS;
}
