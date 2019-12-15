#include <rosa/structure.h>
#include <rosa/version.h>
#include <rosa/error.h>
#include <rosa/util.h>
#include <rosa/elems.h>
#include <rosa/dbg.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

using namespace rosa;
using namespace std;

// Minimum absolute solvent accessible surface for a residue to be considered on
// the surface, in square Angstrom
const double MIN_ACC_SURF_REL = 5.0;

// Maximum distance between the centers of residues for them to be processed
const double MAX_CENTER_DIST = 20.0;

const double MAX_HB_DIST  = 3.5;
const double MAX_SB_DIST  = 5.5;
const double MAX_VDW_DIST = 5.0;

enum InteractionType {
  None = 0,
  Disulfide,
  HydrogenBond,
  SaltBridge,
  VanDerWaals,
  Clash
};

class Contact {
private:
  unsigned long res1Idx;  // Index of residue 1 in the sequence (1 based)
  string        res1Lbl;  // Label of residue 1
  string        res1Name; // Name of residue 1

  unsigned long res2Idx;  // Index of residue 2 in the sequence (1 based)
  string        res2Lbl;  // Label of residue 2
  string        res2Name; // Name of residue 2

  string        type;     // Mainchain-Sidechain, etc. [mm ms sm ss]
public:
  Contact( unsigned long aRes1Idx, const string &aRes1Lbl, const string &aRes1Name,
           unsigned long aRes2Idx, const string &aRes2Lbl, const string &aRes2Name,
           const string &aType ):
    res1Idx(aRes1Idx), res1Lbl(aRes1Lbl), res1Name(aRes1Name),
    res2Idx(aRes2Idx), res2Lbl(aRes2Lbl), res2Name(aRes2Name), type(aType)
  {}

  string getType() const { return type; }

  /*! \brief Writes the contact in the 3did format to the stream passed as an
   *         argument.
   *
   *  \param outStream the stream where the output is sent to */
  void writeContact( ostream& outStream ) const
  {
    outStream << res1Name << "\t" << res2Name << "\t"
              << res1Idx  << "\t" << res2Idx  << "\t"
              << res1Lbl  << "\t" << res2Lbl  << "\t"
              << type;
  }

};


//! Calls the writePDBAtomRecord function on the stream "os"
ostream& operator<<( ostream& os, const Contact& ct )
{
  ct.writeContact(os);
  return os;
}


typedef vector<Contact>                 ContactVect;
typedef vector<Contact>::iterator       ContactVectIt;
typedef vector<Contact>::const_iterator ContactVectCit;


inline bool compareContactsByBonality( const Contact &a, const Contact &b )
{
  return (a.getType() > b.getType());
}


void sortVectorByBonality( ContactVect &cv )
{
  std::sort( cv.begin(), cv.end(), compareContactsByBonality );
}


inline int originalInterpretsBonality( const Atom &a )
{
  // The name has to contain no whitespaces for this to work
  string atomName( a.getName() );
  string resName( a.getResName() );

  if( atomName == "C"  ||
      atomName == "CA" ||
      atomName == "N"  ||
      atomName == "O" ) {
    return 1; // 1 = main chain
  }

  if( (atomName[0]  == 'C' && atomName.length() > 2) ||
      (atomName == "CB" && resName == "ALA") ||
      (atomName == "CE" && resName == "MET") ||
      (atomName[0] == 'N') ||
      (atomName[0] == 'O') ) {
    return 2; // 2 = side chain
  }

  return 0; // 0 = ignore
}


/*
0         1         2         3         4         5         6         7
01234567890123456789012345678901234567890123456789012345678901234567890123456789

REM  Relative accessibilites read from external file "/aloy/data/programs/pac2/naccess/standard.data"
REM  File of summed (Sum) and % (per.) accessibilities for
REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
RES VAL A  18    93.09  61.5  70.54  61.7  22.55  60.7  71.02  61.5  22.07  61.4
RES THR A  19    73.97  53.1  45.12  44.4  28.85  76.8  27.24  36.0  46.73  73.5
...
RES ILE A 102   125.76  71.8  76.57  55.5  49.19 132.4  90.78  65.2  34.98  97.2
END  Absolute sums over single chains surface
CHAIN  1 A     5221.1       4505.8        715.3       2831.1       2390.0
END  Absolute sums over all chains
TOTAL          5221.1       4505.8        715.3       2831.1       2390.0
*/
typedef pair<string,string> ResSpec;
typedef map<ResSpec,float>  RsaMap;


inline void getAtomsBaricenter( AtomVectCit &s, AtomVectCit &e, Point &center ) {
  vector<Point> v;
  for( AtomVectCit it = s; it != e; ++it ) {
    v.push_back(it->getPos());
  }
  getBaricenter( v, center );
}


inline void loadRsaFile( const string &aFilename, RsaMap &resRsa ) {
  DBG_BLOCK_SUPPRESS("loadRsaFile");
  ifstream inFile( aFilename.c_str() );

  if( !inFile.good() )
    throw runtime_error( "error while opening file "+aFilename );

  while( !inFile.eof() ) {
    string line;
    if( getline( inFile, line ) ) {
      DBG_VDUMP(line);
      if( line.substr(0, 3) == "RES" ) {
        string chainId = line.substr( 8, 1 );
        string resLbl  = trim( line.substr( 9, 5 ) );
        double rsa     = stod( line.substr( 22, 6 ) );
        resRsa[ResSpec(chainId,resLbl)] = rsa;
        DBG_MSG(chainId+resLbl+" -> "+dtos(rsa,-1,2));
      }
    }
  }
  inFile.close();
}


// Options
static struct option long_options[] = {
  { "version",      no_argument,       0, 'v' },  // Outputs only the version number
  { "help",         no_argument,       0, 'h' },  // Outputs a short guide
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] pdbFile1 rsaFile1 pdbFile2 rsaFile2" << endl << endl;
  cout << "Returns the number of contacts between the two chains." << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help             prints a short guide" << endl;
  cout << "   -v       --version          prints the version number" << endl;
}


int main( int argc, char *argv[] )
{
  DBG_BLOCK_SUPPRESS( "contact::main" );
  string programName( argv[0] );

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
        errorExit( EXIT_FAILURE, "Invalid command line argument" );
        break;
    }
  }

  argc -= optind;
  argv += optind;

  if( argc != 4 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }

  string struct1Filename( argv[0] );
  string rsa1Filename( argv[1] );
  string struct2Filename( argv[2] );
  string rsa2Filename( argv[3] );

  RsaMap rsaMap1, rsaMap2;

  Structure struct1( struct1Filename );
  loadRsaFile(rsa1Filename,rsaMap1);

  Structure struct2( struct2Filename );
  loadRsaFile(rsa2Filename,rsaMap2);

  if( struct1.numChains() > 1 || struct2.numChains() > 1 ) {
    errorExit( EXIT_FAILURE, "Works only on structures with one chain." );
  }

  if( struct1.numChains() < 1 || struct2.numChains() < 1 )
    return EXIT_SUCCESS;

  // We process the two structures residue by residue
  const Structure::Chain &c1 = struct1.getChain(0);
  const Structure::Chain &c2 = struct2.getChain(0);

  ContactVect contacts;

  for( long ri = 0; ri < c1.numResidues(); ++ri ) {

    AtomVectCitPair vp1 = c1.getResExtr( ri );
    ResSpec rs1(c1.getChainID(),trim(c1.getResLabel(ri)));


    if( rsaMap1.find(rs1) != rsaMap1.end() && rsaMap1[rs1] < MIN_ACC_SURF_REL )
      continue;

    Point rc1; getAtomsBaricenter( vp1.beg, vp1.end, rc1 );

    // For every residue of the second chain
    for( long rj = 0; rj < c2.numResidues(); ++rj ) {

      AtomVectCitPair vp2 = c2.getResExtr( rj );
      ResSpec rs2(c2.getChainID(),trim(c2.getResLabel(rj)));

      if( DBG_ON ) {
        if( trim(c2.getResLabel(rj)) == "138") {
          DBG_MSG("RSA res 138:"+dtos(rsaMap2[rs2],-1,2));
        }
      }

      if( rsaMap2.find(rs2) != rsaMap2.end() && rsaMap2[rs2] < MIN_ACC_SURF_REL )
        continue;

      Point rc2; getAtomsBaricenter( vp2.beg, vp2.end, rc2 );

      double rcDist = pointDistance(rc1,rc2);

      if( rcDist > MAX_CENTER_DIST )
        continue;


      /* This function will return the number of contacts found between two
       * residues as defined in 'Aloy et al. PNAS 2002':
       * 1. hydrogen bonds: N-O distances <= 3.5 A
       * 2. Salt bridges: N-O distances <= 5.5 A
       * 3. van der Waals interactions: C-C distances <= 5A
       *
       * hydrogen bonds seem to be included in salt bridges definition... ?
       */

      ContactVect rcv;

      for( AtomVectCit A = vp1.beg; A != vp1.end; ++A ) {
        for( AtomVectCit B = vp2.beg; B != vp2.end; ++B ) {
          int bA = originalInterpretsBonality( *A );
          int bB = originalInterpretsBonality( *B );

          if( bA > 0 && bB > 0 ) {

            string type = ((bA==1)?((bB==1)?"mm":"ms"):((bB==1)?"sm":"ss"));

            // van der Waals
            if( (A->getElemNo() == Elem_C) && (B->getElemNo() == Elem_C) ) {
              double atomDist = pointDistance( A->getPos(), B->getPos() );
              if( atomDist <= MAX_VDW_DIST ) {
                rcv.push_back( Contact( ri+1, trim(Atom::getResLabel(A->getResNum(),A->getInsCode())), A->getResName(),
                                        rj+1, trim(Atom::getResLabel(B->getResNum(),B->getInsCode())), B->getResName(),
                                        type ) );
              }
              continue;
            }
            //Hydrogen bond/salt bridge
            if( ((A->getElemNo() == Elem_N) && (B->getElemNo() == Elem_O)) ||
                ((A->getElemNo() == Elem_O) && (B->getElemNo() == Elem_N)) ) {
              double atomDist = pointDistance( A->getPos(), B->getPos() );
              // Electrostatic
              if( atomDist > MAX_SB_DIST )
                continue;
              /*
              if( atomDist <= MAX_HB_DIST ) {
                // Hydrogen bond
              }
              */
              rcv.push_back( Contact( ri+1, trim(Atom::getResLabel(A->getResNum(),A->getInsCode())), A->getResName(),
                                      rj+1, trim(Atom::getResLabel(B->getResNum(),B->getInsCode())), B->getResName(),
                                      type ) );
            }
          }
        }
      }

      if( rcv.size() > 0 ) {
        sortVectorByBonality(rcv);
        contacts.push_back(rcv.front());
      }
    }
  }

  for( ContactVectCit ci = contacts.begin(); ci != contacts.end(); ++ci ) {
    cout << *ci << endl;
  }

  return EXIT_SUCCESS;
}
