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

/*! \file pdb.h
 *  \brief Declares constants related to the PDB format
 */

#ifndef ROSA_PDB_H_
#define ROSA_PDB_H_

#include <string>
#include <set>

namespace rosa {
  
  const int PDB_LINE_LENGTH = 80;
  
  const std::string HeaderRecordName( "HEADER" );
  const std::string AtomRecordName(   "ATOM  " );
  const std::string TerRecordName(    "TER   " );
  const std::string HetatmRecordName( "HETATM" );
  const std::string HelixRecordName(  "HELIX " );
  const std::string SheetRecordName(  "SHEET " );
  const std::string TurnRecordName(   "TURN  " );
  const std::string SeqresRecordName( "SEQRES" );
  const std::string RemarkRecordName( "REMARK" );
  const std::string CompndRecordName( "COMPND" );
  const std::string JrnlRecordName(   "JRNL  " );
  const std::string TitleRecordName(  "TITLE " );
  const std::string RevdatRecordName( "REVDAT" );
  const std::string Cryst1RecordName( "CRYST1" );
  const std::string ModelRecordName(  "MODEL " );
  const std::string EndmdlRecordName( "ENDMDL" );
  
  //! Returns true if the three letter code that is passed as the second
  //! argument corresponds to an amino acid (see in the following the
  //! allowed three letter codes).
  bool isAminoacid( const std::string &resName );
  //! Returns true if the three letter code that is passed as the second
  //! argument corresponds to a nucleic acid (see in the following the
  //! allowed three letter codes).
  bool isNucleicacid( const std::string &resName );
  //! Returns true if the three letter code that is passed as the second
  //! argument corresponds to a solvent molecule
  bool isSolvent( const std::string &resName );
  
  //! Three letter code for every amino acid
  const char * const aminoAcidThreeLetterCode[] =
  {
    "ALA",  //  0
    "ARG",  //  1
    "ASN",  //  2
    "ASP",  //  3
    "CYS",  //  4
    "GLN",  //  5
    "GLU",  //  6
    "GLY",  //  7
    "HIS",  //  8
    "ILE",  //  9
    "LEU",  // 10
    "LYS",  // 11
    "MET",  // 12
    "PHE",  // 13
    "PRO",  // 14
    "SER",  // 15
    "THR",  // 16
    "TRP",  // 17
    "TYR",  // 18
    "VAL",  // 19
    "ASX",  // 20
    "GLX",  // 21
  };
  
  //! One letter code for every amino acid
  const char aminoAcidOneLetterCode[] =
  {
    'A', 'R', 'N', 'D', 'C',
    'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V',
    'B', 'Z'
  };

  //! Length of the aminoAcidThreeLetterCode array
  const int numAminoacids = 22;

  //! Three letter code for non standard amino acids
  const char * const nonStdAminoAcidThreeLetterCode[] =
  {
    "SEP",  //  0  PHOSPHOSERINE
    "SET",  //  1  AMINOSERINE
    "TPO",  //  2  PHOSPHOTHREONINE
    "THC",  //  3  N-METHYLCARBONYLTHREONINE
    "BMT",  //  4  4-METHYL-4-[(E)-2-BUTENYL]-4,N-METHYL-THREONINE
    "PTR",  //  5  PHOSPHOTYROSINE
    "TYS",  //  6  O-SULFO-L-TYROSINE
    "FLT",  //  7  FLUOROMALONYL TYROSINE
    "ALB",  //  8  DELTA-2-ALBOMYCIN A1
    "ABA",  //  9  ALPHA-AMINOBUTYRIC ACID
    "DAL",  // 10  D-ALANINE
    "MVA",  // 11  N-METHYLVALINE
    "SMC",  // 12  S-METHYLCYSTEINE
    "ASA",  // 13  ASPARTIC ALDEHYDE
    "ASK",  // 14  DEHYDROXYMETHYLASPARTIC ACID
    "NLE",  // 15  NORLEUCINE
    "MLE",  // 16  N-METHYLLEUCINE
    "CLE",  // 17  LEUCINE AMIDE
    "DLY",  // 18  D-LYSINE
    "ALY",  // 19  N(6)-ACETYLLYSINE
    "MLY",  // 20  N-DIMETHYL-LYSINE
    "M3L",  // 21  N-TRIMETHYLLYSINE
    "MLZ",  // 22  N-METHYL-LYSINE
    "HMR",  // 23  BETA-HOMOARGININE
    "HSE",  // 24  HOMOSERINE
    "HYP",  // 25  HYDROXYPROLINE
    "DPR",  // 26  D-PROLINE
    "HYL"   // 27  HYDROXYLYSINE
  };
  
  //! One letter code for every amino acid
  const char nonStdAminoAcidOneLetterCode[] =
  {
    'S', 'S', 'T', 'T', 'T',
    'Y', 'Y', 'Y', 'A', 'A',
    'A', 'V', 'C', 'D', 'D',
    'L', 'L', 'L', 'K', 'K',
    'K', 'K', 'K', 'R', 'S',
    'P', 'P', 'K'
  };

  //! Length of the aminoAcidThreeLetterCode array
  const int numNonStdAminoacids = 28;
  
  //! Three letter code for every nucleic acid
  const char * const nucleicAcidThreeLetterCode[] =
  {
    "  A",  //  0
    "  C",  //  1
    "  G",  //  2
    "  I",  //  3
    "  T",  //  4
    "  U",  //  5
    " +A",  //  6
    " +C",  //  7
    " +G",  //  8
    " +I",  //  9
    " +T",  // 10
    " +U"   // 11
  };

  //! One letter code for every nucleic acid
  static const char nucleicAcidOneLetterCode[] =
  {
    'A', 'C', 'G', 'I', 'T', 'U', 'A', 'C', 'G', 'I', 'T', 'U', 'X'
  };

  //! Length of the nucleicAcidThreeLetterCode array
  const int numNucleicacids = 12;

  //! The default tolerance for PDB coordinates (0.001)
  const double DEFAULT_PDB_TOLERANCE = 0.001;

  //! Convert the three letter code to one letter code
  char res3to1( const std::string &astr, bool bTranslNonStd = false );

  //! Check if the set of atoms is complete for the specified residue
  //! Return 1 if it is, 0 if it is not and -1 if the residue is unknown
  //! The input set of atoms must not include hydrogens
  int isResComplete( const std::string &aResName, const std::set<std::string> &aAtoms );

  //! Check if the set of atoms only contains main chain atoms
  //! Return true if it is, false otherwise
  bool hasResOnlyMainChain( const std::set<std::string> &aAtoms );

} // namespace rosa

#endif // ROSA_PDB_H_
