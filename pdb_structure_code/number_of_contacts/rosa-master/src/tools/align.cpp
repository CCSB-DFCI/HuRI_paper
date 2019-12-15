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
#include <rosa/seqalign.h>
#include <rosa/version.h>
#include <rosa/error.h>
#include <rosa/dbg.h>
#include <getopt.h>
#include <fstream>
#include <set>


using namespace rosa;
using namespace std;

// Options
static struct option long_options[] = {
  { "version",    no_argument,       0, 'v' },  // Outputs only the version number
  { "rmsd",       no_argument,       0, 'r' },  // RMSD of the seq. alignment based superposition
  { "pedro",      no_argument,       0, 'p' },  // RMSD of the seq. alignment based superposition
  { "outaligned", required_argument, 0, 'O' },  // saves the common atoms of the aligned structures
  { "noaltloc",   no_argument,       0, 'a' },  // removes multiple conformations of residues
  //{ "text",       no_argument,       0, 't' },  // Outputs the alignment in text format
  { "help",       no_argument,       0, 'h' },  // Outputs a short guide  
  { 0, 0, 0, 0 }
};


static void usage( const string &programName )
{
  cout << endl;
  cout << "Usage: " << programName << " [OPTIONS] pdbFile1 chain1 pdbFile2 chain2" << endl << endl;
  cout << "OPTIONS" << endl;
  cout << "   -h       --help                   prints a short guide" << endl;
  cout << "   -r       --rmsd                   returns the RMSD of the seq. alignment based superposition" << endl;
  cout << "   -O       --outaligned=FILE1,FILE2 saves the common atoms of the aligned structures" << endl;
  cout << "   -a       --noaltloc               removes multiple conformations of residues" << endl;
  //cout << "   -t       --text             outputs the alignment in text format" << endl;
  cout << "   -v       --version                prints the version number" << endl << endl;
  cout << "If the chain id is the string \"FIRST\" the first chain is selected." << endl << endl;

}

 
int main( int argc, char **argv )
{
  string programName( argv[0] );
  
  DBG_BLOCK_SUPPRESS( programName+"::main" );
  
  bool bRmsd       = false;
  bool bNoAltLoc   = false;
  bool bOutAligned = false;
  bool bPedro      = false;
  vector<string> outFiles;

  while( true ) {
    int option_index = 0;
    opterr = 0;

    int ch = getopt_long( argc, argv, ":vrapO:h", long_options, &option_index );

    if( ch == -1 ) // Finished with the options
      break;
    
    // Options parsing
    switch( ch ) {
      case 'v':
        printRosaVersion( programName );
        return EXIT_SUCCESS;
        break;
      case 'r':
        bRmsd = true;
        break;
      case 'a':
        bNoAltLoc = true;
        break;
      case 'p':
        bPedro = true;
        break;
      case 'O':
        bOutAligned = true;
        split(string(optarg),outFiles,",");
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
  
  if( argc != 4 ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid command line argument" );
  }

  if( bOutAligned && (outFiles.size() != 2) ) {
    printRosaVersion( programName );
    usage( programName );
    errorExit( EXIT_FAILURE, "Invalid -O argument. The correct format is \"FILE1,FILE2\"" );
  }

  DBG_MSG( "Mark1" );
  
  string pdbFilename1( argv[0] );
  string selChain1( argv[1] );

  string pdbFilename2( argv[2] );
  string selChain2( argv[3] );

  Structure s1;
  s1.loadFromPDBFile( pdbFilename1, false, true, 1 );
  if( selChain1 == "_" ) selChain1 = " ";
  DBG_VDUMP(s1.numChains());
  if( selChain1 == "FIRST" && s1.numChains() > 0 ) selChain1 = s1.getChain(0).getChainID();
  DBG_VDUMP(selChain1);
  if( DBG_ON ) {
    s1.saveToPDBFile("s1.pdb");
  }
  StructurePtr chain1( s1.select( "chain \""+selChain1+"\"" ) );
  DBG_VDUMP(chain1->size());
  if( bNoAltLoc ) {
    DBG_MSG( "s1.removeMultConformations()" );
    chain1->removeMultConformations();
  }
  
  Structure s2;
  s2.loadFromPDBFile( pdbFilename2, false, true, 1 );
  if( selChain2 == "_" ) selChain2 = " ";
  if( selChain2 == "FIRST" && s2.numChains() > 0 ) selChain2 = s2.getChain(0).getChainID();
  StructurePtr chain2( s2.select( "chain \""+selChain2+"\"" ) );
  if( bNoAltLoc ) {
    chain2->removeMultConformations();
  }

  DBG_VDUMP(chain1->numChains());
  DBG_VDUMP(chain2->numChains());

  string seq1 = chain1->getChain(0).getSequence();
  string seq2 = chain2->getChain(0).getSequence();

  SeqAlignment seqAl;
        
  if( chain1->getChain(0).getChainType() == Structure::Chain::NucleicAcidType)
    alignNucAcidSeq( seq1, seq2, seqAl );
  else
    alignProteinSeq( seq1, seq2, seqAl );

  DBG_MSG( "Mark3" );

  cout << "Length 1st struct:              " << setw(8) << seq1.size() << endl;
  cout << "Length 2nd struct:              " << setw(8) << seq2.size() << endl;
  cout << "# aligned residues:             " << setw(8) << seqAl.aln.size() << endl;
  cout << "# identical residues:           " << setw(8) << seqAl.idRes << endl;
  cout << "% coverage:                     " << fixed << setprecision(1) << setw(8) << float(seqAl.aln.size())/float(min(seq1.size(),seq2.size()))*100.0 << "%" << endl;
  cout << "% identity:                     " << fixed << setprecision(1) << setw(8) << float(seqAl.idRes)/float(seqAl.aln.size())*100.0 << "%" << endl;

  StructurePtr chain1CA( chain1->select( "name CA" ) );
  StructurePtr chain2CA( chain2->select( "name CA" ) );
  reduceToAlignedAtoms( *chain1CA, *chain2CA );
  cout << "Gaps in the aligned 1st struct: " << setw(8) << (chain1CA->isContinuous()?"FALSE":"TRUE") << endl;
  cout << "Gaps in the aligned 2nd struct: " << setw(8) << (chain2CA->isContinuous()?"FALSE":"TRUE") << endl;

  DBG_MSG( "Mark4" );

  if( bRmsd ) {
    cout << "RMSD superpos. on aligned res.: " << fixed << setprecision(3) << setw(8) << calcRmsdSuperposition( *chain1CA, *chain2CA ) << " Ang." << endl;
  }
  
  if( bOutAligned ) {
    if( !bPedro ) {
      reduceToAlignedAtoms( *chain1, *chain2 );
    }
    chain1->saveToPDBFile(outFiles[0]+".pdb",0);
    chain2->saveToPDBFile(outFiles[1]+".pdb",0);
    if( bPedro ) {
      ofstream al1File(string(outFiles[0]+".aln").c_str());
      for( int n = 0; n < chain1CA->size(); ++n )
        al1File << trim(Atom::getResLabel((*chain1CA)[n].getResNum(), (*chain1CA)[n].getInsCode())) << endl;
      al1File.close();
      ofstream al2File(string(outFiles[1]+".aln").c_str());
      for( int n = 0; n < chain2CA->size(); ++n )
        al2File << trim(Atom::getResLabel((*chain2CA)[n].getResNum(), (*chain2CA)[n].getInsCode())) << endl;
      al2File.close();
    }
  }
  return EXIT_SUCCESS;
}
