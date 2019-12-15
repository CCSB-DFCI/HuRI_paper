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

/*! \file structure.h
 *  \brief Declares the class Structure implementing a 'structure'. A structure
 *         is in an abstract sense a collection of atoms. They should be usually
 *         related by the fact that they together represent a molecule or an
 *         ensemble of molecules derived from an experiments of crystallography,
 *         NMR, ...
 */

#ifndef ROSA_STRUCTURE_H_
#define ROSA_STRUCTURE_H_

#include <rosa/dbg.h>
#include <rosa/common.h>
#include <rosa/atom.h>
#include <rosa/geometry.h>
#include <rosa/mathutil.h>
#include <rosa/mem_shared_ptr.h>
#include <vector>
#include <map>

namespace rosa {

  //! The distance threshold between two adjacient Ca to be considered continuous
  const double CA_CA_DIST_THRESHOLD = 4.4;

  //! An STL vector of "Atom" objects
  typedef std::vector<Atom> AtomVect;

  //! Iterator to an AtomVect
  typedef std::vector<Atom>::iterator AtomVectIt;

  //! Constant iterator to an AtomVect
  typedef std::vector<Atom>::const_iterator AtomVectCit;

  struct AminoAcidModel {
    Residue res;
    Point   center;
    float   radius;

    AminoAcidModel( const std::string &aChainID, long aResNum,
                    const std::string &aInsCode, const std::string &aResName,
                    Point aCenter, float aRadius ):
      res( aChainID, aResNum, aInsCode, aResName ),
      center( aCenter ),
      radius( aRadius )
    {}


    //! Initializes the amino acid model with the weighted center of mass
    //! and the equivalent radius of the set of atoms between the two
    //! iterators
    AminoAcidModel( const AtomVectCit &beg, const AtomVectCit &end );
  };

  //! An STL vector of "AminoAcidModel" objects
  typedef std::vector<AminoAcidModel> AAModelVect;

  //! Iterator to an AAModelVect
  typedef std::vector<AminoAcidModel>::iterator AAModelVectIt;

  //! Constant iterator to an AAModelVect
  typedef std::vector<AminoAcidModel>::const_iterator AAModelVectCit;

  //! A pair of begin,end iterators to atoms
  struct AtomVectCitPair {
    AtomVectCit beg;
    AtomVectCit end;
    AtomVectCitPair( const AtomVectCit &aBeg, const AtomVectCit &aEnd ):
      beg( aBeg ), end( aEnd )
     {}
  };

  //! An STL vector of "AtomVect" used to store different models
  typedef std::vector<AtomVect> ModelVect;

  //! Iterator to a ModelVect
  typedef std::vector<AtomVect>::iterator ModelVectIt;

  //! Constant iterator to a ModelVect
  typedef std::vector<AtomVect>::const_iterator ModelVectCit;

  //! Forward declaration
  class Structure;

  //! Smart pointer to Structure
  typedef shared_ptr<Structure> StructurePtr;

  //! This class represent a molecule as a collection or atoms.
  class Structure {
  public:
    //! Describes a single chain in the PDB
    class Chain {
    public:
      typedef enum  {
        UnknownType = 0,
        ProteinType,
        NucleicAcidType
      } ChainType;

    private:
      //! Iterator to the first atom in the atom vector
      AtomVectCit begIt;
      //! Iterator to one after the last atom in the atom vector
      AtomVectCit endIt;
      //! Chain ID
      std::string chainID;
      //! Chain type
      ChainType cType;
      //! Sequence of residues (without non-standard residues)
      std::string sequence;
      //! Sequence of residues (with non-standard residues)
      std::string sequenceNonStd;
      //! Vector of residue labels
      std::vector<std::string> resLbls;
      //! Vector of pairs of iterators pointing at the first and last+1 atoms of
      //! the residue
      std::vector<AtomVectCitPair> resExtr;
      //! Vector of bools tracking if the residue has a CA
      std::vector<bool> resHasCA;

      //! Parses the chain to determine the number of residues
      void parseChain();

      //! Initializes extra helper structures
      void buildSequence();

      //! Creates a new chain with the given indices
      Chain( const AtomVectCit &aBeg, const AtomVectCit &aEnd );

      //! Only Structure can create/copy/destroy a chain
      friend class Structure;

    public:
      //! Returns a constant iterator to the beginning of the chain
      AtomVectCit begin() const { return begIt; }
      //! Returns a constant iterator to the end of the chain
      AtomVectCit end() const { return endIt; }
      //! Returns the number of residues/nucleotides in the chain
      long numResidues() const { return sequence.size(); }
      //! Returns the chain ID
      const std::string &getChainID() const { return chainID; }
      //! Returns the number of atoms in the chain
      long numAtoms() const { return endIt - begIt; }
      //! Returns the chain type
      ChainType getChainType() const { return cType; }
      //! Returns a string describing the chain
      std::string getDescription() const;
      //! Returns the sequence of residues/amino acids in the chain
      std::string getSequence( bool bTranslNonStd = false ) const { if( !bTranslNonStd ) return sequence; else return sequenceNonStd; }
      //! Returns the label of the i-th residue
      std::string getResLabel( unsigned long resIdx ) const {
        if( resIdx < resLbls.size() ) return resLbls[resIdx];
        else return std::string();
      }
      //! Returns the extremes of the i-th residue
      AtomVectCitPair getResExtr( unsigned long resIdx ) const {
        if( resIdx < resExtr.size() ) return resExtr[resIdx];
        else return AtomVectCitPair( endIt, endIt );
      }
      //! Returns true if i-th residue has a CA atom
      bool hasResCA( unsigned long resIdx ) const {
        if( resIdx < resExtr.size() ) return resHasCA[resIdx];
        else return false;
      }
      //! Returns 1 if the i-th residue is complete, 0 if it is not and -1
      //! if the residue is unknown. Does not take into consideration hydrogen
      //! atoms
      int isComplete( unsigned long resIdx ) const;
      //! Returns true if the i-th residue only contains main chain atoms.
      //! Does not take into consideration hydrogen atoms
      bool hasOnlyMainChain( unsigned long resIdx ) const;
      //! Returns true if the i-th residue only contains atoms with the
      //! alternative location set
      bool hasOnlyAltLoc( unsigned long resIdx ) const;

      //! Returns true if the chain of CAs respect the constraint "foreach i: dist(Ca(i),Ca(i+1)) < CA_CA_DIST_THRESHOLD"
      bool isContinuous() const;

      //! Gets the baricenter of the chain
      Point getBaricenter() const;

      //! Calculates the average radius of the chain
      double calcAverageRadius() const;

      //! Returns a vector that contains an the 'amino acid model' of the
      //! corresponding chain in the original structure. This model is created
      //! such that for every residue in the original structure we have here one
      //! atom whose center is the weighted center of mass of the residue and
      //! whose radius is the equivalent radius of the residue.
      //! The amino acid model is calculated for the requested MODEL in the
      //! structure (by default the representative model, the modelIndex is 0
      //! based, so that MODEL 1 has index 0, MODEL 2 has index 1, ...)
      AAModelVect buildAminoAcidModel() const;

      //! Returns the surface distance between the projection of the two points
      //! to the surface of the chain, considering a distance radT from the
      //! surface (we imagine a rolling probe on the surface with radius radT
      double getSurfaceDistance( const Point &p1, const Point &p2, double radT = 1.4 ) const;
    };

  private:
    //! Different models are stored as vectors of atoms.
    ModelVect models;

    //! Ids for models are stored as a vector of strings
    std::vector<int> modelIds;

    //! Reference model. The index of the model to be used as a representative
    //! model for the structure.
    long reprModel;

    //! True if the MODEL record is present in the PDB file
    bool bModels;

    //! The model for which chains were built
    mutable int chainsModelIndex;

    //! Vector of chain objects describing the different chains in the file
    mutable std::vector<Chain> chains;

    //! A list containing one record for every missing residue (parsed from the PDB file)
    std::vector<Residue> missingRes;

    //! Classification string from the PDB file
    std::string classification;

    //! Deposition date from the PDB file
    std::string depDate;

    //! PDB ID
    std::string pdbID;

    //! Title from the PDB file
    std::string title;

    //! Reinitialize the structure to empty
    void reinit();

    //! Parse chains
    void buildChains( std::vector<Chain> &chains, int modelIndex = -1 ) const;

    //! Rebuild chains information
    void rebuildChains( int modelIndex = -1 ) const;

    //! Invalidate the chains
    void invalidateChains() const {
      chains.clear();
      chainsModelIndex = -1;
    }

    //! Need to build chains?
    bool needToBuildChains( int modelIndex ) const {
      return (((modelIndex == -1)?(chainsModelIndex!=reprModel):(chainsModelIndex!=modelIndex)) or
              (chains.empty()));
    }

    // Parses remark 465 (missing residues) in a PDB file
    void readRemark465( const std::string &line );

  public:
    Structure();

    //! Creates a structure by loading it from a PDB file
    Structure( const std::string &aFilename,
               bool bAllowDiffModels = false,
               bool bStopAtTer = true );

    //! Loads the content of a PDB file
    //! Checks that all the models have the same atoms. If not it aborts. If you
    //! want to read PDB files with models of different sequence of atoms put
    //! the argument 'bAllowDiffModels' to true. Normally the parser does not
    //! read the content of a chain that lies after a TER record. If you want to
    //! read it put the argument 'bStopAtTer' to false.
    void loadFromPDBFile( const std::string &aFilename,
                          bool bAllowDiffModels = false,
                          bool bStopAtTer = true,
                          int  loadOnlyModel = 0 );

    //! Saves the structure to a PDB file. If modelIndex > 0 only the
    //! corresponding model will be saved
    void saveToPDBFile( const std::string &aFilename, const int modelIndexToSave = -1 );

    //! Returns the number of atoms
    std::size_t size(const int modelIndexToSize = -1) const {
      int i = (modelIndexToSize<0)?reprModel:modelIndexToSize;
      return static_cast<std::size_t>(i<(int)models.size()?models[i].size():0);
    }

    //! Returns the number of models
    std::size_t numModels() const { return static_cast<std::size_t>(models.size()); }

    //! Returns the MODEL id of the model indexed by i. If there is no such model
    //! or if the PDB file contains no MODEL record -1 is returned
    int getModelId( int i ) const {
      if( i < 0 or i > int(modelIds.size())-1 )
        return -1;
      else
        return modelIds[i];
    }

    //! Returns a reference to the i-th atom
    const Atom &operator [] ( unsigned long i ) const { return models[reprModel][i]; }

    //! Returns the classification
    std::string getClassification() const { return classification; }

    //! Returns the deposition date
    std::string getDepDate() const { return depDate; }

    //! Returns the PDB ID
    std::string getPdbID() const { return pdbID; }

    //! Returns the title contained in the PDB
    std::string getTitle() const { return title; }

    //! Sets the PDB ID
    void setPdbID( const std::string &aPdbId ) { pdbID = aPdbId; }

    //! Applies a rotation of a given rotation matrix
    void rotate( const Mat2D<Coord> &rotMat );

    //! Applies a translation of a given vector
    void translate( const Point &v );

    //! Gets the baricenter of the structure
    Point getBaricenter() const;

    //! Compares the two structures. Returns the number of atoms beeing
    //! different positions (0 if the two structures are the same). Two atoms
    //! are considered to be in the same position if their distance is less
    //! than the tolerance passed as the argument. The two atoms also need to
    //! have the same name, residue, occupancy, ...
    unsigned long compareAtoms( const Structure &s2, int compFlags = CHECK_ALL,
                                double tolerance = DEFAULT_PDB_TOLERANCE ) const;

    //! Calculates and returns the radius of gyration of the model indexed by
    //! modelIndex (which is 0 based, so that MODEL 1 has index 0, MODEL 2 has
    //! index 1, ...)
    double calcRadiusOfGyration( int modelIndex = -1 ) const;

    //! Fills the vector of points with the positions of the atoms in the
    //! structure
    long getAtomPositions( std::vector<Point> &v, int modelIndex = -1 ) const;

    //! Returns a dynamically allocated Structure containing only the atoms
    //! selected by the selection string. In order to maintain consistency
    //! in the models the selection is calculated based on a particular model
    //! and applied to all the models.
    StructurePtr select( const std::string &selection ) const;

    //! Returns a dynamically allocated Structure containing only the atoms
    //! of residues in the range [startRes,endRes) for chain aChain. Indices
    //! are 0 based. Only the representative model is returned.
    StructurePtr crop( const std::string &aChain,
                       unsigned long startRes,
                       unsigned long endRes ) const;


    //! Returns the number of chains
    unsigned long numChains( int modelIndex = -1 ) const {
      if( needToBuildChains(modelIndex) ) rebuildChains( modelIndex );
      return chains.size();
    }

    //! Returns the index of the representative model
    unsigned long getReprModel() const { return reprModel; }

    //! Returns the i-th chain of the PDB file
    const Chain &getChain( unsigned long chainIdx, int modelIndex = -1 ) const;

    //! Add the atoms from another structure
    void addAtoms( const Structure &s );

    //! Remove multiple conformations
    void removeMultConformations();

    //! Returns true if all the chains are continuous
    bool isContinuous() const;

    friend long reduceToCommonAtoms( Structure &s1, Structure &s2 );

    friend long reduceToAlignedAtoms( Structure &s1, Structure &s2 );

    friend StructurePtr selectInterface( const Structure &s1,
                                         const Structure &s2,
                                         float dist );

    friend StructurePtr selectBindingSite( const Structure &rec,
                                           const Structure &lig,
                                           float dist );

    friend StructurePtr mapAtoms( const Structure &subsetEq1,
                                  const Structure &eq1,
                                  const Structure &eq2 );

    friend unsigned long getContacts( const Structure &s1, const Structure &s2,
                           float dist,
                           std::vector<std::pair<Residue,Residue> > &resPairs );

    friend unsigned long getResidueEquiv( const Structure &eq1,
                           const Structure &eq2,
                           std::map<Residue,Residue> &resEquivMap );

    friend unsigned long countAtomContacts( const Structure &s1,
                                            const Structure &s2,
                                            float dist );
  };

  //! Calculates the superposition between the two structures
  double calcSuperposition( const Structure &ref, const Structure &toMove,
                            Mat2D<double> &rotMat, Point &cm1, Point &cm2 );

  //! Superpose the second structure on the first
  double superpose( const Structure &ref, Structure &toMove );

  //! Superpose the second structure on the first only on the selected atoms
  double superpose( const Structure &ref, const std::string &refSelString,
                    Structure &toMove, const std::string &toMoveSelString );

  //! Returns the RMSD of the best superposition between the two structures
  double calcRmsdSuperposition( const Structure &ref, const Structure &toMove );

  //! Returns the RMSD between the two structures without superposition
  double calcRmsd( const Structure &s1, const Structure &s2 );

  //! Eliminates from the two structures the atoms that are not in common. This
  //! is done by identifying the atoms having the same residue number, residue
  //! name, atom name, insertion code and alt. loc. code. The two structures
  //! must contain the same number of chains. Even if the chainIDs do not
  //! correspond they will be matched in the order they appear in the two
  //! structures. Alternative location codes need to correspond for an atom to
  //! be preserved. Common atoms listed in different orders are reordered to
  //! match.
  //! At the end the two structures will contain only the common atoms listed in
  //! the same order. Returns the number of atoms preserved in a single
  //! structure. The two structures need to have the same number of models.
  //! Every model will be checked separately.
  long reduceToCommonAtoms( Structure &s1, Structure &s2 );

  //! Eliminates from the two structures the atoms that are not in common. This
  //! is done by identifying throw protein sequence alignment the common
  //! residues for every chain and then by selecting for every common residue
  //! all the atoms that are present in both the structures. The two structures
  //! must contain the same number of chains. Even if the chainIDs do not
  //! correspond they will be matched in the order they appear in the two
  //! structures. Alternative location codes need to correspond for an atom to
  //! be preserved. Common atoms listed in different orders are reordered to
  //! match.
  //! At the end the two structures will contain only the common atoms liste in
  //! the same order. Returns the number of atoms preserved in a single
  //! structure. The two structures need to have the same number of models.
  //! Every model will be checked separately.
  long reduceToAlignedAtoms( Structure &s1, Structure &s2 );

  //! Prints a description of the chains composing the structure
  void printChainDescription( const Structure &s, std::ostream &os = std::cout );

  //! Selects residues in the two strutures that have at least one atom that
  //! is within a certain distance from any atom in the other structure. Returns
  //! a Structure containing the union of the atoms in the two structures
  //! belonging to the interface. The distances are determined on the different
  //! models. At the end the different models can have different number of atoms.
  StructurePtr selectInterface( const Structure &s1, const Structure &s2,
                                float dist );

  //! Selects residues in the first struture that have at least one atom that
  //! is within a certain distance from any atom in the other structure. Returns
  //! a Structure containing the atoms in the first structure belonging to the
  //! binding site. The distances are determined on the different
  //! models. At the end the different models can have different number of atoms.
  StructurePtr selectBindingSite( const Structure &rec, const Structure &lig,
                                  float dist );

  //! Maps a set of atoms in one structure to a set of atoms in another
  //! structure. Only the atoms that are present both in the first set and
  //! in the equivalence sets are mapped. The secon and third structure must
  //! have the same number of models, chains and atoms to establish an
  //! equivalence. Returns the subset of the third structure corresponding to
  //! the mapping of the first structure through the equivalence map.
  StructurePtr mapAtoms( const Structure &subsetEq1, const Structure &eq1,
                         const Structure &eq2 );

  //! Maps between residues
  typedef std::map<Residue,Residue> ResResMap;

  //! Maps residues in the first structure to residues in the second structure
  //! by considering that the i-th atom in the first structure correspond
  //! to the i-th atom in the second structure. The new equivalences are added
  //! to the residue map without clearing it. Returns the number of equivalences
  //! added.
  unsigned long getResidueEquiv( const Structure &eq1, const Structure &eq2,
                                 ResResMap &resEquivMap );

  //! Vector of residue pairs
  typedef std::vector<std::pair<Residue,Residue> > ResPairVect;

  //! Fills the vector resPairs a list of residue pairs corresponding to the
  //! residue-residue contacts between the two structures given the distance
  //! threshold dist. Returns the number of new contacts found. The contacts
  //! are pushed on the back of the vector without clearing it.
  unsigned long getResContacts( const Structure &s1, const Structure &s2,
                                float dist,
                                ResPairVect &resPairs );

  //! Returns the number of contacts between atoms. Two atoms are considered
  //! to be in contact if the distance is less than dist
  unsigned long countAtomContacts( const Structure &s1, const Structure &s2,
                                   float dist );
}

#endif // ROSA_STRUCTURE_H_
