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

/*! \file seqalign.h
 *  \brief Routines for sequence alignment (including dynamic programming seq
 *         alignment 
 */

#ifndef ROSA_SEQALIGN_H_
#define ROSA_SEQALIGN_H_

#include <string>
#include <vector>


namespace rosa {
  //! Stores a pair of indexes
  struct IdxPair {
    long i1, i2;
    
    IdxPair( long aI1, long aI2 ):
      i1( aI1 ), i2( aI2 )
    {}
  };


  //! Stores an alignment as a vector of index pairs
  struct SeqAlignment {
    int score;          //!< Score of the alignment
    long idRes;         //!< Number of identical residues
    std::vector<IdxPair> aln;  //!< The alignment stored as a vector of pairs of indices
    SeqAlignment():
      score( 0 ), idRes( 0 )
    {}
  };


  //! Substitution matrix for sequence alignment with the dyn. progr. algorithm
  //! PAM 250 substitution matrix
  extern const int pam250[];
  //! PAM 120 substitution matrix
  extern const int pam120[];
  //! DNA/RNA substitution matrix
  extern const int dnaRna[];

  //! Finds the longest common substring between the two strings passed as the
  //! arguments. Returns the length of the longest common substring and, if it
  //! is longer than 0, puts in 'index1' and 'index2' the indices in the two
  //! strings of the beginning point.
  long longestCommonSubstring( const std::string& str1, const std::string& str2,
                               long &index1, long &index2 );

  //! Performs a dynamic programming alignment of the two sequences following
  //! the Smith-Waterman algorithm. To be used for protein sequences.
  long alignProteinSeq( const std::string &aQuery, const std::string &aText,
                        SeqAlignment &aSeqAln, const int *substMat = pam250,
                        int gapOp=-12, int gapEx=-2 );

  //! Performs a dynamic programming alignment of the two sequences following
  //! the Smith-Waterman algorithm. To be used for nucleic acids sequences.
  long alignNucAcidSeq( const std::string &aQuery, const std::string &aText,
                        SeqAlignment &aSeqAln, const int *substMat = dnaRna,
                        int gapOp=-12, int gapEx=-4 );

} // namespace rosa

#endif // ROSA_SEQALIGN_H_
