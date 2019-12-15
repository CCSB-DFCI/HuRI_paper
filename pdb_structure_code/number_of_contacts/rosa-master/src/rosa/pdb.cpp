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

#include <rosa/pdb.h>
#include <rosa/dbg.h>
#include <cstring>
#include <set>
#include <map>
#include <algorithm>
 
using namespace rosa;
using namespace std;
 
class AAcodes {
private:
  set<string> aa;
  
  AAcodes()
  {
    for( int i = 0; i < numAminoacids; i++ )
      aa.insert( aminoAcidThreeLetterCode[i] );
  }
  friend bool rosa::isAminoacid( const string &resName );
};

class NAcodes {
private:
  set<string> na;
  
  NAcodes()
  {
    for( int i = 0; i < numAminoacids; ++i )
      na.insert( aminoAcidThreeLetterCode[i] );
  }
      
  friend bool rosa::isNucleicacid( const string &resName );
};

class SolventCodes {
private:
  set<string> solvent;
  
  SolventCodes() {
    solvent.insert( "HOH" );
    solvent.insert( " CL" );
    solvent.insert( " KR" );
  }
      
  friend bool rosa::isSolvent( const string &resName );
};


typedef set<string> AtomSet;


string ALAatoms[] = {"N","CA","O","C","CB"};
string ARGatoms[] = {"N","CA","O","C","CB","CG","CD","NE","CZ","NH1","NH2"};
string ASNatoms[] = {"N","CA","O","C","CB","CG","OD1","ND2"};
string ASPatoms[] = {"N","CA","O","C","CB","CG","OD1","OD2"};
string CYSatoms[] = {"N","CA","O","C","CB","SG"};
string GLNatoms[] = {"N","CA","O","C","CB","CG","CD","OE1","NE2"};
string GLUatoms[] = {"N","CA","O","C","CB","CG","CD","OE1","OE2"};
string GLYatoms[] = {"N","CA","O","C"};
string HISatoms[] = {"N","CA","O","C","CB","CG","ND1","CD2","CE1","NE2"};
string ILEatoms[] = {"N","CA","O","C","CB","CG1","CG2","CD1"};
string LEUatoms[] = {"N","CA","O","C","CB","CG","CD1","CD2"};
string LYSatoms[] = {"N","CA","O","C","CB","CG","CD","CE","NZ"};
string METatoms[] = {"N","CA","O","C","CB","CG","SD","CE"};
string PHEatoms[] = {"N","CA","O","C","CB","CD1","CD2","CE1","CE2","CZ"};
string PROatoms[] = {"N","CA","O","C","CB","CG","CD"};
string SERatoms[] = {"N","CA","O","C","CB","OG"};
string THRatoms[] = {"N","CA","O","C","CB","OG1","CG2"};
string TRPatoms[] = {"N","CA","O","C","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"};
string TYRatoms[] = {"N","CA","O","C","CB","CG","CD1","CD2","CE1","CE2","CZ","OH"};
string VALatoms[] = {"N","CA","O","C","CB","CG1","CG2"};

//! Contains sets of atoms for every residue
class ResAtoms {
private:
  map<string,AtomSet> resAtomSets;

  ResAtoms()
  {
    resAtomSets["ALA"] = AtomSet(ALAatoms, ALAatoms + sizeof(ALAatoms) / sizeof(ALAatoms[0]));
    resAtomSets["ARG"] = AtomSet(ARGatoms, ARGatoms + sizeof(ARGatoms) / sizeof(ARGatoms[0]));
    resAtomSets["ASN"] = AtomSet(ASNatoms, ASNatoms + sizeof(ASNatoms) / sizeof(ASNatoms[0]));
    resAtomSets["ASP"] = AtomSet(ASPatoms, ASPatoms + sizeof(ASPatoms) / sizeof(ASPatoms[0]));
    resAtomSets["CYS"] = AtomSet(CYSatoms, CYSatoms + sizeof(CYSatoms) / sizeof(CYSatoms[0]));
    resAtomSets["GLN"] = AtomSet(GLNatoms, GLNatoms + sizeof(GLNatoms) / sizeof(GLNatoms[0]));
    resAtomSets["GLU"] = AtomSet(GLUatoms, GLUatoms + sizeof(GLUatoms) / sizeof(GLUatoms[0]));
    resAtomSets["GLY"] = AtomSet(GLYatoms, GLYatoms + sizeof(GLYatoms) / sizeof(GLYatoms[0]));
    resAtomSets["HIS"] = AtomSet(HISatoms, HISatoms + sizeof(HISatoms) / sizeof(HISatoms[0]));
    resAtomSets["ILE"] = AtomSet(ILEatoms, ILEatoms + sizeof(ILEatoms) / sizeof(ILEatoms[0]));
    resAtomSets["LEU"] = AtomSet(LEUatoms, LEUatoms + sizeof(LEUatoms) / sizeof(LEUatoms[0]));
    resAtomSets["LYS"] = AtomSet(LYSatoms, LYSatoms + sizeof(LYSatoms) / sizeof(LYSatoms[0]));
    resAtomSets["MET"] = AtomSet(METatoms, METatoms + sizeof(METatoms) / sizeof(METatoms[0]));
    resAtomSets["PHE"] = AtomSet(PHEatoms, PHEatoms + sizeof(PHEatoms) / sizeof(PHEatoms[0]));
    resAtomSets["PRO"] = AtomSet(PROatoms, PROatoms + sizeof(PROatoms) / sizeof(PROatoms[0]));
    resAtomSets["SER"] = AtomSet(SERatoms, SERatoms + sizeof(SERatoms) / sizeof(SERatoms[0]));
    resAtomSets["THR"] = AtomSet(THRatoms, THRatoms + sizeof(THRatoms) / sizeof(THRatoms[0]));
    resAtomSets["TRP"] = AtomSet(TRPatoms, TRPatoms + sizeof(TRPatoms) / sizeof(TRPatoms[0]));
    resAtomSets["TYR"] = AtomSet(TYRatoms, TYRatoms + sizeof(TYRatoms) / sizeof(TYRatoms[0]));
    resAtomSets["VAL"] = AtomSet(VALatoms, VALatoms + sizeof(VALatoms) / sizeof(VALatoms[0]));
  }

  friend int rosa::isResComplete( const string &aResName, const set<string> &aAtoms );
};


int rosa::isResComplete( const string &aResName, const set<string> &aAtoms )
{
  static ResAtoms resAtoms;

  if( isAminoacid(aResName) ) {
    const AtomSet &resAtomSet = resAtoms.resAtomSets[aResName];
    return (includes(resAtomSet.begin(),resAtomSet.end(),aAtoms.begin(),aAtoms.end())?1:0);
  }
  return -1;
}


string MCatoms[]  = {"N","CA","O","C"};


bool rosa::hasResOnlyMainChain( const set<string> &aAtoms )
{
  static AtomSet mca = AtomSet(MCatoms, MCatoms + sizeof(MCatoms) / sizeof(MCatoms[0]));

  return includes(mca.begin(),mca.end(),aAtoms.begin(),aAtoms.end());
}


bool rosa::isAminoacid( const std::string &resName )
{
  static AAcodes aaCodes;

  return (aaCodes.aa.find( resName ) != aaCodes.aa.end());
}


bool rosa::isNucleicacid( const std::string &resName )
{
  static NAcodes naCodes;

  return (naCodes.na.find( resName ) != naCodes.na.end());
}


bool rosa::isSolvent( const std::string &resName )
{
  static SolventCodes solventCodes;

  return (solventCodes.solvent.find( resName ) != solventCodes.solvent.end());
}


char rosa::res3to1( const string &res, bool bTranslNonStd )
{
  DBG_BLOCK_SUPPRESS("rosa::res3to1");

  DBG_VDUMP(res);
  DBG_VDUMP(bTranslNonStd);
  
  for( int i = 0; i < numAminoacids; ++i )
    if( strncasecmp( res.c_str(), aminoAcidThreeLetterCode[i], 3 ) == 0 )
      return aminoAcidOneLetterCode[i];
  for( int i = 0; i < numNucleicacids; ++i )
    if( strncasecmp( res.c_str(), nucleicAcidThreeLetterCode[i], 3 ) == 0 )
      return nucleicAcidOneLetterCode[i];
  if( bTranslNonStd ) {
    for( int i = 0; i < numNonStdAminoacids; ++i )
      if( strncasecmp( res.c_str(), nonStdAminoAcidThreeLetterCode[i], 3 ) == 0 )
        return nonStdAminoAcidOneLetterCode[i];
  }  
  return 'X';
}
