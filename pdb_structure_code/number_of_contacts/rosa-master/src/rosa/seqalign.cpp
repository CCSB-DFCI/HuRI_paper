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

#include <rosa/seqalign.h>
#include <rosa/mathutil.h>
#include <rosa/dbg.h>
#include <rosa/util.h>
#include <fstream>
#include <algorithm>

using namespace std;

namespace rosa {
  
const int protSubMatRowLength = 23;

//! PAM 250 substitution matrix
const int pam250[] = {
   3, -3, -1,  0, -3, -1,  0,  1, -3, -1, -3, -2, -2, -4,  1,  1,  1, -7, -4,  0,  0, -1, -1,
  -3,  6, -1, -3, -4,  1, -3, -4,  1, -2, -4,  2, -1, -5, -1, -1, -2,  1, -5, -3, -2, -1, -2,
  -1, -1,  4,  2, -5,  0,  1,  0,  2, -2, -4,  1, -3, -4, -2,  1,  0, -4, -2, -3,  3,  0, -1,
   0, -3,  2,  5, -7,  1,  3,  0,  0, -3, -5, -1, -4, -7, -3,  0, -1, -8, -5, -3,  4,  3, -2,
  -3, -4, -5, -7,  9, -7, -7, -4, -4, -3, -7, -7, -6, -6, -4,  0, -3, -8, -1, -3, -6, -7, -4,
  -1,  1,  0,  1, -7,  6,  2, -3,  3, -3, -2,  0, -1, -6,  0, -2, -2, -6, -5, -3,  0,  4, -1,
   0, -3,  1,  3, -7,  2,  5, -1, -1, -3, -4, -1, -3, -7, -2, -1, -2, -8, -5, -3,  3,  4, -1,
   1, -4,  0,  0, -4, -3, -1,  5, -4, -4, -5, -3, -4, -5, -2,  1, -1, -8, -6, -2,  0, -2, -2,
  -3,  1,  2,  0, -4,  3, -1, -4,  7, -4, -3, -2, -4, -3, -1, -2, -3, -3, -1, -3,  1,  1, -2,
  -1, -2, -2, -3, -3, -3, -3, -4, -4,  6,  1, -3,  1,  0, -3, -2,  0, -6, -2,  3, -3, -3, -1,
  -3, -4, -4, -5, -7, -2, -4, -5, -3,  1,  5, -4,  3,  0, -3, -4, -3, -3, -2,  1, -4, -3, -2,
  -2,  2,  1, -1, -7,  0, -1, -3, -2, -3, -4,  5,  0, -7, -2, -1, -1, -5, -5, -4,  0, -1, -2,
  -2, -1, -3, -4, -6, -1, -3, -4, -4,  1,  3,  0,  8, -1, -3, -2, -1, -6, -4,  1, -4, -2, -2,
  -4, -5, -4, -7, -6, -6, -7, -5, -3,  0,  0, -7, -1,  8, -5, -3, -4, -1,  4, -3, -5, -6, -3,
   1, -1, -2, -3, -4,  0, -2, -2, -1, -3, -3, -2, -3, -5,  6,  1, -1, -7, -6, -2, -2, -1, -2,
   1, -1,  1,  0,  0, -2, -1,  1, -2, -2, -4, -1, -2, -3,  1,  3,  2, -2, -3, -2,  0, -1, -1,
   1, -2,  0, -1, -3, -2, -2, -1, -3,  0, -3, -1, -1, -4, -1,  2,  4, -6, -3,  0,  0, -2, -1,
  -7,  1, -4, -8, -8, -6, -8, -8, -3, -6, -3, -5, -6, -1, -7, -2, -6, 12, -2, -8, -6, -7, -5,
  -4, -5, -2, -5, -1, -5, -5, -6, -1, -2, -2, -5, -4,  4, -6, -3, -3, -2,  8, -3, -3, -5, -3,
   0, -3, -3, -3, -3, -3, -3, -2, -3,  3,  1, -4,  1, -3, -2, -2,  0, -8, -3,  5, -3, -3, -1,
   0, -2,  3,  4, -6,  0,  3,  0,  1, -3, -4,  0, -4, -5, -2,  0,  0, -6, -3, -3,  4,  2, -1,
  -1, -1,  0,  3, -7,  4,  4, -2,  1, -3, -3, -1, -2, -6, -1, -1, -2, -7, -5, -3,  2,  4, -1,
  -1, -2, -1, -2, -4, -1, -1, -2, -2, -1, -2, -2, -2, -3, -2, -1, -1, -5, -3, -1, -1, -1, -2
};

//! PAM 120 substitution matrix
const int pam120[] = {
   2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0,  0,  0,
  -2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1,  0, -1,
   0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2,  1,  0,
   0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3,  3, -1,
  -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -3,
   0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1,  3, -1,
   0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3,  3, -1,
   1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0,  0, -1,
  -1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1,  2, -1,
  -1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2, -2, -1,
  -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3, -3, -1,
  -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1,  0, -1,
  -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2, -2, -1,
  -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4, -5, -2,
   1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1,  0, -1,
   1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0,  0,  0,
   1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1,  0,
  -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -6, -4,
  -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -4, -2,
   0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2, -2, -1,
   0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3,  2, -1,
   0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2,  3, -1,
   0, -1,  0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1,  0,  0, -4, -2, -1, -1, -1, -1
};

const int dnaSubMatRowLength = 7;

const int dnaRna[] = {
   5, -4, -4, -4, -4, -4, -1,
  -4,  5, -4, -4, -4, -4, -1,
  -4, -4,  5, -4, -4, -4, -1,
  -4, -4, -4,  5, -4, -4, -1,
  -4, -4, -4, -4,  5,  5, -1,
  -4, -4, -4, -4,  5,  5, -1,
  -1, -1, -1, -1, -1, -1, -1 
};


static void res2idx( const string &aStr, vector<int> &aInd )
{
  aInd.resize(aStr.size());
  
  for( unsigned long i = 0; i < aStr.size(); ++i )
    switch( aStr[i] ) {
      case 'A': aInd[i] =  0; break;
      case 'R': aInd[i] =  1; break;
      case 'N': aInd[i] =  2; break;
      case 'D': aInd[i] =  3; break;
      case 'C': aInd[i] =  4; break;
      case 'Q': aInd[i] =  5; break;
      case 'E': aInd[i] =  6; break;
      case 'G': aInd[i] =  7; break;
      case 'H': aInd[i] =  8; break;
      case 'I': aInd[i] =  9; break;
      case 'L': aInd[i] = 10; break;
      case 'K': aInd[i] = 11; break;
      case 'M': aInd[i] = 12; break;
      case 'F': aInd[i] = 13; break;
      case 'P': aInd[i] = 14; break;
      case 'S': aInd[i] = 15; break;
      case 'T': aInd[i] = 16; break;
      case 'W': aInd[i] = 17; break;
      case 'Y': aInd[i] = 18; break;
      case 'V': aInd[i] = 19; break;
      case 'B': aInd[i] = 20; break;
      case 'Z': aInd[i] = 21; break;
      case 'X': aInd[i] = 22; break;
      default:
        throw logic_error( string("res2idx: unknown character ") + aStr[i] );
   }
}

static void na2idx( const string &aStr, vector<int> &aInd )
{
  aInd.resize(aStr.size());
  
  for( unsigned long i = 0; i < aStr.size(); i++ )
    switch( aStr[i] ) {
      case 'A': aInd[i] = 0; break;
      case 'C': aInd[i] = 1; break;
      case 'G': aInd[i] = 2; break;
      case 'I': aInd[i] = 3; break;
      case 'T': aInd[i] = 4; break;
      case 'U': aInd[i] = 5; break;
      case 'X': aInd[i] = 6; break;
      default:
        throw logic_error( string("na2idx: unknown character ") + aStr[i] );
   }
}


long longestCommonSubstring( const string& str1, const string& str2,
                             long &index1, long &index2 )
{
  if(str1.empty() || str2.empty())
    return 0;

  long str1Length = str1.size();
  long str2Length = str2.size();
  
  Mat2D<long> table( str1Length, str2Length );
  long maxSubstr = 0;
  
  for(long i = 0; i < str1Length; ++i)
    for(long j = 0; j < str2Length; ++j) {
      if(str1[i] != str2[j])
        table[i][j] = 0;
      else {
        if(i == 0 || j == 0)
          table[i][j] = 1;
        else
          table[i][j] = 1 + table[i-1][j-1];
        if(maxSubstr < table[i][j]) {
          maxSubstr = table[i][j];
          index1 = i - maxSubstr + 1;
          index2 = j - maxSubstr + 1;
        }
      }
    }
 
  return maxSubstr;
}


//! Actually performs the calculation of the alignment using a maximization
//! function instead of a minimization
//! Given this example
//!
//!              A|A|C|---|A|ATTCCG|A|C|T
//!              A|C|T|ACC|T|------|C|G|C
//!               (a)  (b)     (c)
//! we have three kinds of blocks:
//! - the blocks of type (a) of matches and mismatches
//! - the blocks of type (b) of characters of the second sequence aligned
//!   with gaps
//! - the blocks of type (c) of characters of the first sequence aligned
//!   with gaps
//!
//! The matrices used in the algorithm are the following:
//! btMat: the matrix used for back tracing
//! scMat: the matrix holding the scores of each match and the final score of
//!        the alignments
//! eMat:  the matrix holding the scores of the alignments for blocks of type (b)
//! fMat:  the matrix holding the scores of the alignments for blocks of type (c)
//!
static long dynAlign( const vector<int> &Aseq, const vector<int> &Bseq,
                      SeqAlignment &aSeqAln, const int *substMat,
                      int matRowLength, int gop, int gex )
{
  DBG_BLOCK_SUPPRESS( "dynAlign" );

  int sc1a, sc1b, sc2a, sc2b, sc3a, sc3b, sc3c;
  
  long Asize = Aseq.size();
  long Bsize = Bseq.size();

  DBG_VDUMP( Asize );
  DBG_VDUMP( Bsize );

  DBG_VDUMP( matRowLength );
  DBG_VDUMP( gop );
  DBG_VDUMP( gex );

  /*********************************************************************
   *                          MATRIX INITIALIZATION                    *
   *********************************************************************/

  // Score matrices
  Mat2D<long> e( Asize+1, Bsize+1, 0 );
  Mat2D<long> f( Asize+1, Bsize+1, 0 );
  Mat2D<long> s( Asize+1, Bsize+1, 0 );
  
  // Backtracing matrices
  Mat2D<int> be( Asize+1, Bsize+1, -1 );
  Mat2D<int> bf( Asize+1, Bsize+1, -1 );
  Mat2D<int> bs( Asize+1, Bsize+1, -1 );
  

  /*********************************************************************
   *                            DEBUG STUFF                            *
   *********************************************************************/
  if( DBG_ON ) {
    ofstream dbgMatFile( "dbgMatrixInit.txt" );
    dbgMatFile << "E" << endl;
    e.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "F" << endl;
    f.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "S" << endl;
    s.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "BE" << endl;
    be.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "BF" << endl;
    bf.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "BS" << endl;
    bs.print( true, dbgMatFile );
    dbgMatFile.close();
  }
  
  
  /*********************************************************************
   *                           DP ALGORITHM                            *
   *********************************************************************/
  
  for( long i = 1L; i < Asize+1; ++i )
    for( long j = 1L; j < Bsize+1; ++j ) {
      
      // ======================== From the top ========================
      // Align first sequence chars with gaps in the second sequence
      
      sc1a = e[i-1][j] + gex;
      sc1b = s[i-1][j] + (gex + gop);
      
      if( sc1a >= sc1b ) {
        e[i][j]  = sc1a;
        be[i][j] = 0;
      } else {
        e[i][j]  = sc1b;
        be[i][j] = 3;
      }
      
      // ======================== From the left =======================
      // Align second sequence chars with gaps in the first sequence
      
      sc2a = f[i][j-1] + gex;
      sc2b = s[i][j-1] + (gex + gop);
      
      if( sc2a >= sc2b ) {
        f[i][j]  = sc2a;
        bf[i][j] = 2;
      } else {
        f[i][j]  = sc2b;
        bf[i][j] = 4;
      }
      
      // ========================== Diagonal ==========================
      sc3a = s[i-1][j-1] + substMat[Aseq[i-1]*matRowLength+Bseq[j-1]];
      sc3b = e[i][j];
      sc3c = f[i][j];
      
      
      //  ======================== Best one? =======================
      if( sc3a < 0 && sc3b < 0 && sc3c < 0 ) {
        s[i][j]  =  0;
        bs[i][j] = -1;
      } else if( sc3a >= sc3b && sc3a >= sc3c ) {
        s[i][j]  = sc3a;
        bs[i][j] =  1;
      } else if( sc3b >= sc3a && sc3b >= sc3c ) {
        s[i][j]  = sc3b;
        bs[i][j] =  5;
      } else {
        s[i][j]  = sc3c;
        bs[i][j] =  6;        
      }
    }
      
  if( DBG_ON ) {
    ofstream dbgMatFile( "dbgMatrixFinal.txt" );
    dbgMatFile << "E" << endl;
    e.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "F" << endl;
    f.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "S" << endl;
    s.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "BE" << endl;
    be.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "BF" << endl;
    bf.print( true, dbgMatFile );
    dbgMatFile << endl << endl << "BS" << endl;
    bs.print( true, dbgMatFile );
    dbgMatFile.close();
  }
      
  // Reconstruction of the alignment
  aSeqAln.aln.clear();

  unsigned long i = Asize;
  unsigned long j = Bsize;

  // First we have to find the maximum element of the score matrix
  aSeqAln.score = s.getMax( i, j );

  int op = bs[i][j];
  
  DBG_VDUMP( aSeqAln.score );
  DBG_VDUMP( i );
  DBG_VDUMP( j );
  
  while( op != -1 ) {
    if( DBG_ON ) {
      cout << "op: " << (int) op
      << "; i: " << setw(6) << i
      << "; j: " << setw(6) << j << endl; cout.flush();
    }
    
    switch( op ) {
      case 0:
        --i;
        op = be[i][j];
        break;
      case 1:
        aSeqAln.aln.push_back( IdxPair( i, j ) );
        if( Aseq[i-1] == Bseq[j-1] ) ++(aSeqAln.idRes);
        --i; --j;
        op = bs[i][j];
        break;
      case 2:
        --j;
        op = bf[i][j];
        break;
      case 3:
        --i;
        op = bs[i][j];
        break;
      case 4:
        --j;
        op = bs[i][j];
        break;
      case 5:
        op = be[i][j];
        break;
      case 6:
        op = bf[i][j];
        break;
      default:
        throw logic_error( "dynAlign: unexpected case "+ltos(op) );
        break;
    }
  }
  
  if( !aSeqAln.aln.empty() )
    reverse( aSeqAln.aln.begin(), aSeqAln.aln.end() );
  
  return aSeqAln.aln.size();
}


long alignProteinSeq( const std::string &aQuery, const std::string &aText,
                      SeqAlignment &aSeqAln, const int *substMat,
                      int gapOp, int gapEx )
{
  DBG_BLOCK_SUPPRESS( "alignProteinSeq" );
  
  vector<int> Aseq, Bseq;
  
  DBG_VDUMP( aQuery );
  DBG_VDUMP( aText );
  
  res2idx( aQuery, Aseq );
  res2idx( aText,  Bseq );

  DBG_VDUMP( Aseq.size() );
  DBG_VDUMP( Bseq.size() );
  
  long toReturn = dynAlign( Aseq, Bseq, aSeqAln, substMat,
                            protSubMatRowLength, gapOp, gapEx );
  
  if( DBG_ON ) {
    for( long k = 0UL; k < toReturn; ++k ) {
      cout << setw(3) << k << ": " << setw(3) << aSeqAln.aln[k].i1 << "  "
           << aQuery[aSeqAln.aln[k].i1-1] << " - " << aText[aSeqAln.aln[k].i2-1]
           << "  " << setw(3) << aSeqAln.aln[k].i2 << endl;
    }
    cout << "Similarity: " << aSeqAln.aln.size() / (float)max( Aseq.size(), Bseq.size() ) * 100.0 << endl;
    cout << "Identity:   " << aSeqAln.idRes      / (float)min( Aseq.size(), Bseq.size() ) * 100.0 << endl;
  }

  return toReturn;
}


long alignNucAcidSeq( const std::string &aQuery, const std::string &aText,
                      SeqAlignment &aSeqAln, const int *substMat,
                      int gapOp, int gapEx )
{
  DBG_BLOCK_SUPPRESS( "alignNucAcidSeq" );
  
  vector<int> Aseq, Bseq;
  
  DBG_VDUMP( aQuery );
  DBG_VDUMP( aText );
  
  na2idx( aQuery, Aseq );
  na2idx( aText,  Bseq );
  
  DBG_VDUMP( Aseq.size() );
  DBG_VDUMP( Bseq.size() );
  
  long toReturn = dynAlign( Aseq, Bseq, aSeqAln, substMat,
                            dnaSubMatRowLength, gapOp, gapEx );
  
  if( DBG_ON ) {
    for( long k = 0UL; k < toReturn; ++k ) {
      cout << setw(3) << k << ": " << setw(3) << aSeqAln.aln[k].i1 << "  "
           << aQuery[aSeqAln.aln[k].i1-1] << " - " << aText[aSeqAln.aln[k].i2-1]
           << "  " << setw(3) << aSeqAln.aln[k].i2 << endl;
    }
    cout << "Similarity: " << aSeqAln.aln.size() / (float)max( Aseq.size(), Bseq.size() ) * 100.0 << endl;
    cout << "Identity:   " << aSeqAln.idRes      / (float)min( Aseq.size(), Bseq.size() ) * 100.0 << endl;
  }
  
  return toReturn;
}


} // namespace rosa
