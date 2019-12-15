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

/*! \file atom.h
 *  \brief Declares the class Atom implementing a Atom object. A atom is
 *         usually described by its 3-D coordinates plus all the informations
 *         like the element name, the atom name and all the other informations
 *         contained in a PDB file.
 */

#ifndef ROSA_ATOM_H_
#define ROSA_ATOM_H_

#include <rosa/types.h>
#include <rosa/geometry.h>
#include <rosa/mathutil.h>
#include <rosa/elems.h>
#include <rosa/util.h>
#include <rosa/pdb.h>

namespace rosa {

  //! The Atom class contains all the attributes of an atom.
  /*! Atom's attributes include the attributes of an ATOM or HETATM record
   *  in a PDB file, as well as other attributes like the uncertainty in the
   *  atom's coordinates, various flags and other useful informations.
   */
  class Atom {

  private:

    std::string  recordName;  //!< PDB record name
    long         serialNum;   //!< Serial number in the PDB entry
    std::string  origName;    //!< Name of the atom as it is saved in the PDB file
    std::string  altLoc;      //!< Alternative location code
    std::string  resName;     //!< Residue name
    std::string  chainID;     //!< Chain ID
    long         resSeq;      //!< Residue sequence number
    std::string  iCode;       //!< Insertion code

    Point        pos;         //!< Atom coordinates
    float        occ;         //!< Occupancy
    float        temp;        //!< Temperature factor

    bool         ter;         //!< Atom is followed by a TER record

    std::string  elemSymbol;  //!< Element symbol
    std::string  charge;      //!< Charge on the atom

    std::string  trimmedName; //!< Name of the atom without leading and trailing
    //!< spaces

    int          elemNo;      //!< Element number, e.g. index of the element in
    //!< the element table contained in the file elems.h

    float        radius;      //!< VdW radius of the atom (depends on the element)

    //! Sets the elemNo member variable to the correct value
    void setElemNo();

  public:
    Atom();

    //! Returns the atom's name without leading and trailing blanks
    std::string getName() const { return trimmedName; }

    //! Returns the name of the atom stored in the PDB record
    std::string getOrigName() const { return origName; }

    //! Changes the name of the atom stored in the PDB record
    void setOrigName( const std::string &aName ) {
      origName = aName;
      trimmedName = trim(aName);
    }

    //! Returns the atom serial number
    long getSerialNum() const { return serialNum; }

    //! Returns the atom alternate location code
    std::string getAltLoc() const { return altLoc; }

    //! Returns the atom alternate location code
    void setAltLoc( const std::string &newAltLoc ) { altLoc = newAltLoc; }

    //! Returns the atom's insertion code
    std::string getInsCode() const { return iCode; }

    //! Sets the atom insertion code
    void setInsCode( const std::string &newICode ) { iCode = newICode; }

    //! Returns the name of the residue
    std::string getResName() const { return resName; }

    //! Returns the residue sequence number
    long getResNum() const { return resSeq; }

    //! Sets the residue sequence number
    void  setResNum( long aNewResNum ) { resSeq = aNewResNum; }

    //! Returns the ChainID
    std::string getChainID() const { return chainID; }

    //! Sets the ChainID
    void setChainID( const std::string &newChainID ) { chainID = newChainID; }

    //! Returns the PDB record name
    std::string getPdbRecordName() const { return recordName; }

    //! Sets the PDB record type
    void setPdbRecordName( const std::string &aRecordName ) {
      recordName = aRecordName;
    }

    //! Returns the position of the atom as a vector of three REALs
    const Point &getPos() const { return pos; }

    //! Returns the X coordinate of the atom
    Coord getPosX() const { return pos.x(); }

    //! Returns the Y coordinate of the atom
    Coord getPosY() const { return pos.y(); }

    //! Returns the Z coordinate of the atom
    Coord getPosZ() const { return pos.z(); }

    //! Sets the position of the atom as a vector of three REALs
    void setPos( Coord x, Coord y, Coord z ) { pos.x(x); pos.y(y); pos.z(z); }

    //! Sets the position of the atom as a vector of three REALs
    void setPos( const Point &aPos ) {
      pos = aPos;
    }

    //! Returns the temperatur factor
    float getBfac() const { return temp; }

    //! Sets the temperatur factor
    void setBfac( float aTemp ) { temp = aTemp; }

    //! Returns the Occupancy
    float getOcc() const { return occ; }

    //! Sets the Occupancy
    void setOcc( float aOcc ) { occ = aOcc; }

    //! Returns true if the atom is followed by a TER record
    float isTer() const { return ter; }

    //! Sets the TER flag to true
    void setTer() { ter = true; }

    //! Sets the TER flag to false
    void cleanTer() { ter = false; }

    //! Returns the element symbol
    std::string getElemSymbol() const { return elemSymbol; }

    //! Returns the element type index as described in elems.h
    int getElemNo() const { return elemNo; }

    //! Returns the atom's radius
    float getRadius() const {
      if( elemNo != Elem_Unknown && elemNo < ElementsSize )
        return Elements[elemNo].vdw;
      else
        return 0.0;
    }

    //! Returns the atom's mass
    float getMass() const {
      if( elemNo != Elem_Unknown && elemNo < ElementsSize )
        return Elements[elemNo].mass;
      else
        return 0.0;
    }

    //! Creates an identifier that is unique for the atom. The identifier is
    //! composed by the atom name as originally saved in the PDB where the spaces
    //! are substituted by underscore characters followed by the alt loc, the
    //! residue name, the residue number, the insertion code, the chain ID. All
    //! the items are separated by underscore characters
    std::string identifier() const;

    //! Parse the atom record contained in the string aStr
    /*!
     *  \param aStr the string containing the atom PDB record
     *  Throws an exception if the atom is not in standard PDB format
     */
    void readPDBAtomRecord( const std::string &aStr );

    /*! \brief Writes the atom record in the standard PDB file format to the
     *         stream passed as the argument.
     *
     *  \param outStream the stream where the output is sent to */
    void writePDBAtomRecord( std::ostream& outStream ) const;

    /*! \brief Writes the TER record in the standard PDB file format to the
     *         stream passed as the argument.
     *
     *  \param outStream the stream where the output is sent to */
    void writePDBTerRecord( std::ostream& outStream ) const;

    //! Applies a rotation of a given rotation matrix
    void rotate( const Mat2D<Coord> &rotMat );

    //! Applies a translation of a given vector
    void translate( const Point &v );

    //! Builds a residue label by concatenating the residue number with the
    //! residue label
    static std::string getResLabel( long resNum, const std::string &insCode ) {
      return ltos( resNum, 5 ) + insCode;
    }

    friend unsigned int compareAtoms( const Atom &a1, const Atom &a2,
                                      int compFlags, double tolerance );

  };

  //! Calls the writePDBAtomRecord function on the stream "os"
  std::ostream& operator<<( std::ostream& os, const Atom& at );

  //! Constants used for atoms comparison
  enum {
    CHECK_CHAIN    = 0x0001,
    CHECK_RESNAME  = 0x0002,
    CHECK_RESNUM   = 0x0004,
    CHECK_INSCODE  = 0x0008,
    CHECK_ATOMNAME = 0x0010,
    CHECK_ALTLOC   = 0x0020,
    CHECK_ALL      = 0x003F
  };

  //! Compares the two atoms. Returns the number of fields that are different.
  //! The position is compared by calculating the distance between the two
  //! atoms. If this distance is less than the tolerance the two atoms are
  //! considered to be in the same position. If the function returns 0 than
  //! the two atoms are equal.
  unsigned int compareAtoms( const Atom &a1, const Atom &a2,
                             int compFlags = CHECK_ALL,
                             double tolerance = DEFAULT_PDB_TOLERANCE );

  //! The class residue simply contains the description of a residue including
  //! the residue name, the chainID, the residue number and the insertion code.
  struct Residue {
    std::string  chainID;     //!< Chain ID
    long         resNum;      //!< Residue sequence number
    std::string  iCode;       //!< Insertion code
    std::string  resName;     //!< Residue name

    Residue( const std::string &aChainID, long aResNum,
              const std::string &aInsCode, const std::string &aResName = "" ):
      chainID( aChainID ), resNum( aResNum ),
      iCode( aInsCode ), resName( aResName )
    {}
  };

  inline bool operator == ( const Residue &r1, const Residue &r2 ) {
    if( r1.chainID != r2.chainID ) return false;
    if( r1.resNum  != r2.resNum  ) return false;
    if( r1.iCode   != r2.iCode   ) return false;
    if( r1.resName != r2.resName ) return false;
    return true;
  }

  inline bool operator != ( const Residue &r1, const Residue &r2 ) {
    return not rosa::operator == ( r1, r2 );
  }

  inline bool operator < ( const Residue &r1, const Residue &r2 ) {
    if( r1.chainID < r2.chainID ) return true;
    else if( r1.chainID == r2.chainID ) {
      if( r1.resNum < r2.resNum ) return true;
      else if( r1.resNum == r2.resNum ) {
        if( r1.iCode < r2.iCode )
          return true;
      }
    }
    return false;
  }

} // namespace rosa

#endif // _ATOM_H__
