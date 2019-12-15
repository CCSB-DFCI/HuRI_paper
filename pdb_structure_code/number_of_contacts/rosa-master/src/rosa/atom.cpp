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

#include <rosa/common.h>
#include <rosa/atom.h>
#include <rosa/dbg.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>

using namespace std;
using namespace rosa;

void Atom::setElemNo()
{
  DBG_BLOCK_SUPPRESS( "Atom::setElemNo" );
  
  DBG_VDUMP( elemSymbol );
  DBG_VDUMP( origName );
  
  if( elemSymbol != "" ) {
    elemNo = elementIndex( elemSymbol.c_str() );
    
    if( elemNo != Elem_Unknown )
      return;
  }
  
  DBG_VDUMP( origName );
  
  if( origName[0] == 'H' ) {
    elemNo = Elem_H;
    return;
  }
  
  string lElemSymbol;
  if( origName[0] != ' ' )
    lElemSymbol = origName.substr( 0, 2 );
  else if( origName[1] != ' ' )
    lElemSymbol = origName.substr( 1, 1 );
  else
    lElemSymbol = "";
  
  if( lElemSymbol != "" ) {
    elemNo = elementIndex( lElemSymbol.c_str() );
    
    if( elemNo != Elem_Unknown )
      return;
  }
  
  elemNo = Elem_Unknown;
}


Atom::Atom():
  serialNum( -1 ),
  resSeq( -1 ),
  pos( 0.0, 0.0, 0.0 ),
  occ( -1.0 ),
  temp( -1.0 ),
  ter( false ),
  elemNo( -1 )
{}


string Atom::identifier() const
{
  return (replace( origName, " ", "_" ) + "_" + altLoc + "_" + resName + "_" +
          ltos(resSeq) + "_" + iCode + "_" + chainID);
}


void Atom::readPDBAtomRecord( const std::string &aStr )
{
  DBG_BLOCK_SUPPRESS( "Atom::readPDBAtomRecord( aStr )" );
  
  DBG_VDUMP( aStr );

  if( aStr.size() < 66 )
    throw runtime_error( "ATOM record not in standard PDB format!" );
  
  recordName = aStr.substr( 0, 6 );
  
  istringstream inStr( aStr.substr( 6, 5 ).c_str() );
  inStr >> serialNum;
  
  origName   = aStr.substr( 12, 4 );
  altLoc     = aStr.substr( 16, 1 );
  resName    = aStr.substr( 17, 3 );
  chainID    = aStr.substr( 21, 1 );
  
  inStr.str( aStr.substr( 22, 4 ).c_str() ); inStr.clear();
  inStr >> resSeq;
  iCode      = aStr.substr( 26, 1 );
  
  Coord ttmmpp;
  inStr.str( aStr.substr( 30, 8 ).c_str() ); inStr.clear();
  inStr >> ttmmpp; pos.x(ttmmpp);
  inStr.str( aStr.substr( 38, 8 ).c_str() ); inStr.clear();
  inStr >> ttmmpp; pos.y(ttmmpp);
  inStr.str( aStr.substr( 46, 8 ).c_str() ); inStr.clear();
  inStr >> ttmmpp; pos.z(ttmmpp);
  
  inStr.str( aStr.substr( 54, 6 ).c_str() ); inStr.clear();
  inStr >> occ;
  inStr.str( aStr.substr( 60, 6 ).c_str() ); inStr.clear();
  inStr >> temp;

  if( aStr.size() >= 78 )
    elemSymbol = trim(aStr.substr( 76, 2 ));
  
  if( aStr.size() >= 80 )
    charge     = aStr.substr( 78, 2 );

  trimmedName = trim( origName );
  
  setElemNo();
}


void Atom::writePDBAtomRecord( ostream& outStream ) const
{
  outStream << recordName;
  outStream << right << setw(5) << serialNum << " ";
  outStream << origName;
  outStream << altLoc[0];
  outStream << resName << " ";
  outStream << chainID[0];
  outStream << right << setw(4) << resSeq;
  outStream << iCode[0] << "   ";
  outStream << right << fixed << setw(8)   << setprecision(3) << pos.x();
  outStream << right << fixed << setw(8)   << setprecision(3) << pos.y();
  outStream << right << fixed << setw(8)   << setprecision(3) << pos.z();
  outStream << right << fixed << setw(6)   << setprecision(2) << occ;
  outStream << right << fixed << setw(6)   << setprecision(2) << temp;
  outStream << "          ";
  outStream << right << setw(2) << elemSymbol;
  outStream << right << setw(2) << charge;
}


void Atom::writePDBTerRecord( ostream& outStream ) const
{
  outStream << TerRecordName;
  outStream << right << setw(5) << serialNum+1 << " ";
  outStream << "    ";
  outStream << " ";
  outStream << resName << " ";
  outStream << chainID[0];
  outStream << right << setw(4) << resSeq;
  outStream << iCode[0];
  outStream << "                                                     ";
}


void Atom::rotate( const Mat2D<Coord> &rotMat )
{
  setPos( rotMat[0][0]*pos.x() + rotMat[0][1]*pos.y() + rotMat[0][2]*pos.z(),
          rotMat[1][0]*pos.x() + rotMat[1][1]*pos.y() + rotMat[1][2]*pos.z(),
          rotMat[2][0]*pos.x() + rotMat[2][1]*pos.y() + rotMat[2][2]*pos.z() );
}

        
void Atom::translate( const Point &v )
{
  pos += v;
}


namespace rosa {
  
  ostream& operator<<( ostream& os, const Atom& at ) {
    at.writePDBAtomRecord( os );
    return os;
  }

  unsigned int compareAtoms( const Atom &a1, const Atom &a2,
                             int compFlags, double tolerance )
  {
    unsigned int numDiff = 0U;
    
    if( a1.recordName != a2.recordName ) numDiff++; 
    if( a1.serialNum != a2.serialNum ) numDiff++;  

    if( (compFlags | CHECK_CHAIN    ) and (a1.origName != a2.origName) ) numDiff++;
    if( (compFlags | CHECK_RESNAME  ) and (a1.altLoc != a2.altLoc)     ) numDiff++;     
    if( (compFlags | CHECK_RESNUM   ) and (a1.resName != a2.resName)   ) numDiff++;    
    if( (compFlags | CHECK_INSCODE  ) and (a1.chainID != a2.chainID)   ) numDiff++;    
    if( (compFlags | CHECK_ATOMNAME ) and (a1.resSeq != a2.resSeq)     ) numDiff++;     
    if( (compFlags | CHECK_ALTLOC   ) and (a1.iCode != a2.iCode)       ) numDiff++;      

    if( pointDistance( a1.pos, a2.pos ) > tolerance ) numDiff++;
    if( a1.occ != a2.occ ) numDiff++;        
    if( a1.temp != a2.temp ) numDiff++;       
    
    if( a1.elemSymbol != a2.elemSymbol ) numDiff++; 
    if( a1.charge != a2.charge ) numDiff++;     
    
    return numDiff;
  }
} // namespace rosa
