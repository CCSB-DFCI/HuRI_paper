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
#include <rosa/structure.h>
#include <rosa/superpos.h>
#include <rosa/seqalign.h>
#include <rosa/atomsel_driver.h>
#include <rosa/elems.h>
#include <set>
#include <string>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <list>
#include <queue>

using namespace std;
using namespace rosa;

const int NUM_SAMPLING_POINTS_SURF_DIST = 7;

Structure::Structure()
{
  // All the initialization and clean up is done in the reinit() function
  reinit();
}


Structure::Structure( const string &aFilename,
               bool bAllowDiffModels,
               bool bStopAtTer )
{
  loadFromPDBFile ( aFilename, bAllowDiffModels, bStopAtTer );
}


void Structure::reinit()
{
  models.clear();
  modelIds.clear();
  invalidateChains();
  reprModel = -1;
  bModels = false;
  classification = "";
  depDate = "";
  pdbID = "";
  title = "";
}


void Structure::loadFromPDBFile( const string &aFilename,
                                 bool bAllowDiffModels,
                                 bool bStopAtTer,
                                 int loadOnlyModel )
{
  DBG_BLOCK_SUPPRESS( "Structure::loadFromPDBFile" );

  DBG_MSG( aFilename );
  ifstream inFile( aFilename.c_str() );

  if( !inFile.good() )
    throw runtime_error( "error while opening file "+aFilename );

  // Re-initialize the structure to empty
  reinit();

  DBG_MSG( "After re-init" );

  // See the end of the function
  models.push_back( AtomVect() );

  // True if a MODEL record has been already encountered
  bModels = false;

  bool bAfterTer = true;
  set<string> terminatedChains;
  string oldChainID = "";
  while( !inFile.eof() ) {
    string line;
    if( getline( inFile, line ) ) {

      // Pads the length of the line to reach the normal length of 80 characters
      line += string( max<int>( PDB_LINE_LENGTH - line.size(), 0 ), ' ' );

      DBG_VDUMP( line );

      string recordName( line.substr( 0, 6 ) );

      if( recordName == HeaderRecordName ) {
        classification = line.substr( 10, 40 );
        depDate        = line.substr( 50, 9  );
        pdbID          = replace( line.substr( 62, 4  ), " ", "_" );
      } else if( recordName == TitleRecordName ) {
        title += line.substr( 10, 60 );
        title = trim( title );
      } else if( recordName == AtomRecordName ||
                 recordName == HetatmRecordName ) {
        Atom at;
        at.readPDBAtomRecord( line );
        if( at.getChainID() != oldChainID &&
            terminatedChains.find(at.getChainID()) == terminatedChains.end() &&
            bAfterTer )
          bAfterTer = false;
        if( !bAfterTer ) {
          models.back().push_back( at );
          DBG_MSG( "Atom added." );
        }
        oldChainID = at.getChainID();
      } else if( recordName == TerRecordName ) {
        if( line.length() > 21 )
          oldChainID = line.substr(21,1);
        if( oldChainID != "" )
          terminatedChains.insert(oldChainID);
        if( !bAfterTer ) {
          if( !models.back().empty() )
            models.back().back().setTer();
          if( bStopAtTer )
            bAfterTer = true;
        }
      } else if( recordName == ModelRecordName ) {
        if( bModels ) {
          models.push_back( AtomVect() );
        } else
          bModels = true;
        bAfterTer = false;
        terminatedChains.clear();
        modelIds.push_back((int)(stol(trim(line.substr(10,4)))));
      } else if( recordName == EndmdlRecordName ) {
        if( models.back().empty() )
          throw runtime_error( "empty model" );
        bAfterTer = false;
        terminatedChains.clear();
      } else if( recordName == RemarkRecordName ) {
          int remarkNum = -1;

          try {
            stol( line.substr( 7, 3 ) );
          } catch (runtime_error &e) { remarkNum = -1; }

          switch( remarkNum ) {
            case 465: readRemark465(line); break;
          }
      } else if( recordName == HelixRecordName ) {
      } else if( recordName == SheetRecordName ) {
      } else if( recordName == TurnRecordName ) {
      } else if( recordName == RemarkRecordName ) {
      }
    }
  }
  // If no atom was added to the models then we can cancel the only empty model
  // that was created at the beginning
  if( models.back().empty() ) {
    models.pop_back();
    if( !modelIds.empty() )
      modelIds.pop_back();
  }

  if( loadOnlyModel > 0 && models.size() > 1 ) {
    ModelVect newModels;
    ModelVectIt mvi = models.begin();
    for( vector<int>::const_iterator mi = modelIds.begin(); mi != modelIds.end(); ++mi, ++mvi ) {
      if( *mi == loadOnlyModel ) {
        newModels.push_back( AtomVect() );
        newModels.front().swap(*mvi);
      }
    }
    if( newModels.size() < 1 )
      throw runtime_error( "no model with index "+ltos(loadOnlyModel) );
    models.swap(newModels);
  }
  reprModel = 0;

  // Enforce the policy of having all models with the same number of atoms
  if( models.size() > 0 and !bAllowDiffModels ) {
    unsigned long numAtoms = models[0].size();
    ModelVectIt endMvit = models.end();
    for( ModelVectIt mvit = models.begin(); mvit != endMvit; ++mvit )
      if( mvit->size() != numAtoms )
        throw logic_error( "multiple models with different number of atoms" );
  }
}


void Structure::readRemark465( const string &line )
{
  DBG_BLOCK_SUPPRESS( "Structure::readRemark465" );
  // If we find the residue number in the right position
  if( isInteger( line.substr( 21, 5 ) ) ) {

    missingRes.push_back( Residue( line.substr(19,1),
                                   stol( line.substr( 21, 5 ) ),
                                   line.substr( 26, 1 ),
                                   line.substr(15,3) ) );
  }
}


void Structure::saveToPDBFile( const std::string &aFilename,
                               const int modelIndexToSave )
{
  ofstream outFile( aFilename.c_str() );

  if( !outFile.good() )
    throw runtime_error( "error while opening file " + aFilename );

  if( !classification.empty() or !depDate.empty() or !pdbID.empty() ) {
    outFile << HeaderRecordName << "    " << setw(40)  << classification;
    outFile << setw(9)   << depDate << "   " << setw(4)   << pdbID;
    outFile << "              " << endl;
  }

  if( !missingRes.empty() ) {
    outFile << "REMARK 465                                                                      " << endl;
    outFile << "REMARK 465 MISSING RESIDUES                                                     " << endl;
    outFile << "REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE                       " << endl;
    outFile << "REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN               " << endl;
    outFile << "REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)                " << endl;
    outFile << "REMARK 465                                                                      " << endl;
    outFile << "REMARK 465   M RES C SSSEQI                                                     " << endl;

    for( vector<Residue>::iterator rit = missingRes.begin();
         rit != missingRes.end(); ++rit ) {
      outFile << "REMARK 465     " << (*rit).resName << " "
              << (*rit).chainID << setw(6) << (*rit).resNum << (*rit).iCode
              << "                                                     " << endl;
    }
  }

  bool bWriteModelRecord = (((models.size() > 1) || bModels) && (modelIndexToSave < 0));
  string s66( 66, ' ' );
  string s74( 74, ' ' );

  ModelVectIt endMvit = models.end();
  int modelCounter = 1;
  for( ModelVectIt mvit = models.begin(); mvit != endMvit; ++mvit ) {
    if( (modelIndexToSave < 0) || (modelIndexToSave == modelCounter-1) ) {
      if( bWriteModelRecord )
        outFile << ModelRecordName << setw(8) << modelCounter << s66 << endl;
      AtomVectIt endAvit = mvit->end();
      for( AtomVectIt avit = mvit->begin(); avit != endAvit; ++avit ) {
        outFile << (*avit) << endl;
        if( avit->isTer() ) {
          avit->writePDBTerRecord(outFile);
          outFile << endl;
        }
      }
      if( bWriteModelRecord )
        outFile << EndmdlRecordName << s74 << endl;
    }
    modelCounter++;
  }
  outFile << "END                                                                             " << endl;
}


void Structure::rotate( const Mat2D<Coord> &rotMat )
{
  ModelVectIt endMvit = models.end();
  for( ModelVectIt mvit = models.begin(); mvit != endMvit; ++mvit ) {
    AtomVectIt endAvit = mvit->end();
    for( AtomVectIt avit = mvit->begin(); avit != endAvit; ++avit ) {
      avit->rotate( rotMat );
    }
  }
}


void Structure::translate( const Point &v )
{
  ModelVectIt endMvit = models.end();
  for( ModelVectIt mvit = models.begin(); mvit != endMvit; ++mvit ) {
    AtomVectIt endAvit = mvit->end();
    for( AtomVectIt avit = mvit->begin(); avit != endAvit; ++avit ) {
      avit->translate( v );
    }
  }
}

Point Structure::getBaricenter() const
{
  Point baricenter;
  unsigned long numAtoms = 0UL;

  ModelVectCit endMvit = models.end();
  for( ModelVectCit mvit = models.begin(); mvit != endMvit; ++mvit ) {
    AtomVectCit endAvit = mvit->end();
    for( AtomVectCit avit = mvit->begin(); avit != endAvit; ++avit )
      baricenter += avit->getPos();
    numAtoms += mvit->size();
  }
  baricenter /= (double)numAtoms;
  return baricenter;
}

unsigned long Structure::compareAtoms( const Structure &s2, int compFlags,
                                       double tolerance ) const
{
  DBG_BLOCK_SUPPRESS( "Structure::compareAtoms" );
  unsigned long numDiffAtoms = 0UL;

  ModelVectCit endMvit1 = models.end();
  ModelVectCit endMvit2 = s2.models.end();

  ModelVectCit mvit1 = models.begin();
  ModelVectCit mvit2 = s2.models.begin();

  while( mvit1 != endMvit1 and mvit2 != endMvit2 ) {
    unsigned long minNumberOfAtoms = min( mvit1->size(), mvit2->size() );
    numDiffAtoms += max( mvit1->size(), mvit2->size() ) - minNumberOfAtoms;

    for( unsigned long i = 0UL; i < minNumberOfAtoms; i++ ) {
      DBG_MSG( "Comparing "+(*mvit1)[i].identifier()+" with "+(*mvit2)[i].identifier() )
      if( rosa::compareAtoms( (*mvit1)[i], (*mvit2)[i], compFlags, tolerance ) )
        numDiffAtoms++;
    }
    ++mvit1; ++mvit2;
  }

  while( mvit1 != endMvit1 ) {
    numDiffAtoms += mvit1->size();
    ++mvit1;
  }

  while( mvit2 != endMvit2 ) {
    numDiffAtoms += mvit2->size();
    ++mvit2;
  }

  return numDiffAtoms;
}


double Structure::calcRadiusOfGyration( int modelIndex ) const
// Radius of Gyration:
//        Rgyr = sqrt( < ( r_i - < r_i > )^2 > )
{
  if( modelIndex < 0 ) modelIndex = reprModel;
  if( (int)models.size() <= modelIndex ) return -1.0;

  double radiusOfGyration = 0.0;

  Point baricenter( getBaricenter() );

  AtomVectCit endAvit = models[modelIndex].end();
    for( AtomVectCit avit = models[modelIndex].begin(); avit != endAvit; ++avit )
      radiusOfGyration += pointSqrDistance( avit->getPos(), baricenter );

  radiusOfGyration = sqrt(radiusOfGyration/(double)models[modelIndex].size());

  return radiusOfGyration;

}


long Structure::getAtomPositions( std::vector<Point> &v, int modelIndex ) const
{
  if( modelIndex < 0 ) modelIndex = reprModel;
  if( (int)models.size() <= modelIndex ) return 0;
  v.clear();
  v.reserve( models[modelIndex].size() );
  AtomVectCit endAvit = models[modelIndex].end();
  for( AtomVectCit avit = models[modelIndex].begin(); avit != endAvit; ++avit )
    v.push_back( avit->getPos() );
  return v.size();
}


StructurePtr Structure::select( const std::string &selection ) const
{
  DBG_BLOCK_SUPPRESS( "Structure::select" );

  DBG_VDUMP(selection);

  StructurePtr newStruct( new Structure() );

  AtomSelectDriver selDriver( selection );
  DBG_VDUMP( selDriver.good() );
  if( !selDriver.good() )
    return newStruct;

  DBG_VDUMP( reprModel );
  DBG_VDUMP( selection );

  DBG_MSG( "Mark1" );

  newStruct->reprModel      = reprModel;
  newStruct->bModels        = bModels;
  newStruct->classification = classification;
  newStruct->depDate        = depDate;
  newStruct->pdbID          = pdbID;
  newStruct->title          = title;

  DBG_MSG( "Mark2" );
  ModelVectCit endMvit = models.end();
  DBG_MSG( "Mark2a" );
  for( ModelVectCit mvit = models.begin(); mvit != endMvit; ++mvit ) {

    const AtomVect &lMod = *mvit;
    unsigned long numAtoms = lMod.size();
    DBG_VDUMP( numAtoms );

    newStruct->models.push_back( vector<Atom>() );
    AtomVect &av = newStruct->models.back();
    for( unsigned long i = 0UL; i < numAtoms; ++i )
      if( selDriver.evalOnAtom( lMod[i] ) )
        av.push_back( lMod[i] );
  }

  DBG_MSG( "Mark3" );

  return newStruct;
}

StructurePtr Structure::crop( const string &aChain,
                              unsigned long startRes,
                              unsigned long endRes ) const
{
  DBG_BLOCK_SUPPRESS( "Structure::crop" );

  StructurePtr newStruct( new Structure() );

  int chainIdx = -1;

  for( unsigned long ci = 0; ci < numChains(); ++ci ) {
    const Chain &c = getChain(ci);
    if( c.getChainID() == aChain) {
      chainIdx = ci;
      break;
    }
  }

  if( chainIdx < 0 || startRes >= endRes )
    return newStruct;

  const Chain &c = getChain(chainIdx);
  long numResidues = c.numResidues();

  if( (long)startRes >= numResidues )
    return newStruct;

  endRes = (long)endRes > numResidues ? numResidues : endRes;

  DBG_VDUMP( reprModel );

  DBG_MSG( "Mark1" );

  newStruct->reprModel      = 0;
  newStruct->bModels        = false;
  newStruct->classification = classification;
  newStruct->depDate        = depDate;
  newStruct->pdbID          = pdbID;
  newStruct->title          = title;

  AtomVectCit aStart = c.getResExtr(startRes).beg;
  AtomVectCit aEnd   = c.getResExtr(endRes-1).end;

  newStruct->models.push_back( vector<Atom>(aStart,aEnd) );

  DBG_MSG( "Mark3" );

  return newStruct;
}

/*
static void resId( const Atom &a, string &s )
{
  string chID  = a.getChainID();
  long resNum  = a.getChainID();
  string iCode = a.getInsCode();
  s = trim( chID + itos( resNum ) + iCode );
}
*/
void Structure::rebuildChains( int modelIndex ) const
{
  if( modelIndex < 0 )
    modelIndex = reprModel;
  buildChains( chains, modelIndex );
  chainsModelIndex = modelIndex;
}

void Structure::buildChains( std::vector<Chain> &aChains, int modelIndex ) const
{
  if( modelIndex < 0 )
    modelIndex = reprModel;
  aChains.clear();

  const AtomVect &av = models[modelIndex];
  if( av.size() == 0 ) return;
  AtomVectCit cBeg = av.begin();
  AtomVectCit cEnd = cBeg;
  string lChainID = cBeg->getChainID();
  while( cEnd != av.end() ) {
    if( cEnd->getChainID() != lChainID ) {
      aChains.push_back( Chain( cBeg, cEnd ) );
      lChainID = cEnd->getChainID();
      cBeg = cEnd;
    }
    ++cEnd;
  }
  aChains.push_back( Chain( cBeg, av.end() ) );
}


const Structure::Chain &Structure::getChain( unsigned long chainIdx,
                                             int modelIndex ) const
{
  if( needToBuildChains(modelIndex) ) rebuildChains( modelIndex );
  if( chainIdx >= chains.size() )
    throw logic_error( "Structure::getChain: invalid chain index" );

  return chains[chainIdx];
}

void Structure::addAtoms( const Structure &s )
{
  if( s.size() > 0 ) {
    if( s.models.size() != models.size() )
      throw runtime_error( "Structure::addAtoms: structures with different number of models" );

    ModelVectCit endMvit = s.models.end();
    ModelVectIt  it = models.begin();
    for( ModelVectCit mvit = s.models.begin(); mvit != endMvit; ++mvit, ++it )
      it->insert(it->end(), mvit->begin(), mvit->end() );
  }
}


void Structure::removeMultConformations()
{
  DBG_BLOCK_SUPPRESS( "Structure::removeMultConformations" );

  ModelVectIt endMvit = models.end();
  for( ModelVectIt mvit = models.begin(); mvit != endMvit; ++mvit ) {
    set<pair<string,string> > atomsSeen;
    vector<bool> toRemove;
    toRemove.reserve( mvit->size());
    long numberToRemove = 0;
    AtomVectCit endAvit = mvit->end();
    for( AtomVectCit avit = mvit->begin(); avit != endAvit; ++avit ) {
      string resLbl = Atom::getResLabel( avit->getResNum(), avit->getInsCode() );
      string chainID = avit->getChainID();
      pair<string,string> atSeen( chainID+resLbl, avit->getOrigName() );
      if( atomsSeen.find( atSeen ) != atomsSeen.end() ) {
        toRemove.push_back( true );
        ++numberToRemove;
      } else {
        toRemove.push_back( false );
        atomsSeen.insert( atSeen );
      }
    }
    if( numberToRemove > 0 ) {
      vector<Atom> newAtoms;
      newAtoms.reserve( mvit->size()-numberToRemove );
      AtomVectCit avit = mvit->begin();
      for( vector<bool>::const_iterator trit = toRemove.begin(); avit != endAvit; ++avit, ++trit ) {
        if( !*trit )
          newAtoms.push_back( *avit );
      }
      swap( *mvit, newAtoms );
    }
  }
  DBG_VDUMP(size());
  invalidateChains();
}


bool Structure::isContinuous() const
{
  long lNumChains = numChains();
  for( long i = 0; i < lNumChains; ++i )
    if( !(getChain(i).isContinuous()) ) return false;
  return true;
}


double rosa::calcSuperposition( const Structure &ref, const Structure &toMove,
                                Mat2D<double> &rotMat, Point &cm1, Point &cm2 )
{
  vector<Point> v1, v2;

  ref.getAtomPositions( v1 );
  toMove.getAtomPositions( v2 );

  return quaternionSuperposition( v1, v2, rotMat, cm1, cm2 );
}


double rosa::calcRmsdSuperposition( const Structure &ref,
                                    const Structure &toMove )
{
  vector<Point> v1, v2;

  ref.getAtomPositions( v1 );
  toMove.getAtomPositions( v2 );

  return getRmsdQuaternionSuperposition( v1, v2 );
}


double rosa::calcRmsd( const Structure &s1, const Structure &s2 )
{
  vector<Point> v1, v2;

  s1.getAtomPositions( v1 );
  s2.getAtomPositions( v2 );

  return rmsd( v1, v2 );
}


Structure::Chain::Chain( const AtomVectCit &aBeg, const AtomVectCit &aEnd ):
  begIt( aBeg ), endIt( aEnd ), cType( UnknownType )
{
  parseChain();
}


void Structure::Chain::parseChain()
{
  DBG_BLOCK_SUPPRESS( "Structure::Chain::parseChain" );

  if( numAtoms() <= 0 ) {
    cType = UnknownType;
    return;
  }

  chainID = begIt->getChainID();

  unsigned long aaCount     = 0UL; // Number of amino acids
  unsigned long naCount     = 0UL; // Number of nucleic acids
  unsigned long unCount     = 0UL; // Number of unknown res

  unsigned long caAtoms     = 0UL; // Number of CA atoms
  unsigned long pAtoms      = 0UL; // Number of P atoms
  unsigned long totResidues = 0UL; // Total number of different atoms

  long oldResNum = begIt->getResNum()-1;
  string oldInsCode = "";
  bool caPresent = false;
  bool pPresent  = false;

  for( AtomVectCit it = begIt; it != endIt; ++it ) {
    if( it->getResNum() != oldResNum ||
        it->getInsCode() != oldInsCode ) {
      // At the beginning of a new residue
      string resName = it->getResName();
      if( isAminoacid( resName ) )
        ++aaCount;
      else if( isNucleicacid( resName ) )
        ++naCount;
      else
        ++unCount;
      if( !isSolvent(resName) )
        totResidues++;
      oldResNum = it->getResNum();
      oldInsCode = it->getInsCode();
      caPresent = false;
      pPresent = false;
    }
    if( it->getOrigName() == " CA " and !caPresent ) {
      ++caAtoms;
      caPresent = true;
    }
    if( it->getOrigName() == " P  " and !pPresent ) {
      ++pAtoms;
      pPresent = true;
    }
  }

  DBG_VDUMP(aaCount);
  DBG_VDUMP(totResidues);
  DBG_VDUMP(caAtoms);
  DBG_VDUMP(naCount);

  // if the number of aminoacids is higher than 80% then it is a protein
  // if there are CA atoms and there is any residue and the number of CAs is more than 80% the number of
  // residues then it is a protein
  if( (float) (aaCount) / (float) totResidues > 0.8 ) {
    cType = ProteinType;
  } else if( caAtoms > 0 and totResidues > 0 and (float) (caAtoms) / (float) totResidues > 0.8 ) {
    cType = ProteinType;
  } else if( ((float) (naCount) / (float) totResidues > 0.8) or
           (naCount > 0 and totResidues > 0 and (float) (naCount) / (float) totResidues > 0.8) ) {
    // similar reasoning is applied to check if it is RNA...
    cType = NucleicAcidType;
  } else {
    // otherwise it is unknown
    cType = UnknownType;
  }

  buildSequence();
}


string Structure::Chain::getDescription() const
{
  string result;

  if( numAtoms() > 0 ) {
    long   oldResNum = begIt->getResNum();
    string oldICode  = begIt->getInsCode();

    result += trim( chainID + ltos( oldResNum ) + oldICode );
    long   writtenResNum = oldResNum;
    string writtenICode  = oldICode;
    for( AtomVectCit avit = begIt; avit != endIt; ++avit ) {
      long   newResNum = avit->getResNum();
      string newICode  = avit->getInsCode();

      if( abs(newResNum - oldResNum) > 1 or
          newICode != oldICode ) {
        if( abs(writtenResNum-oldResNum) != 0 or writtenICode != oldICode )
          result += "-" + trim( chainID + ltos( oldResNum ) + oldICode );
        result += ".." + trim( chainID + ltos( newResNum ) + newICode );
        writtenResNum = newResNum;
        writtenICode  = newICode;
      } else if( avit+1 == endIt ) {
        if( abs(writtenResNum-newResNum) != 0 or writtenICode != newICode )
          result += "-" + trim( chainID + ltos( newResNum ) + newICode );
      }
      oldResNum = newResNum;
      oldICode  = newICode;
    }
  }
  return result;
}


void Structure::Chain::buildSequence()
{
  DBG_BLOCK_SUPPRESS( "Structure::Chain::buildSequence" );

  long   oldResNum   = begIt->getResNum();
  string oldICode    = begIt->getInsCode();
  AtomVectCit oldBeg = begIt;
  string resLabel    = Atom::getResLabel( oldResNum, oldICode );

  sequence.clear();
  sequenceNonStd.clear();

  char nextChar = res3to1(begIt->getResName(),false);
  sequence       += nextChar;
  if( nextChar == 'X' )
    sequenceNonStd += res3to1(begIt->getResName(),true);
  else
    sequenceNonStd += nextChar;

  bool bHasCA = false;
  for( AtomVectCit avit = begIt; avit != endIt; ++avit ) {
    long   newResNum = avit->getResNum();
    string newICode  = avit->getInsCode();
    if( !bHasCA )
      if( avit->getOrigName() == " CA " )
        bHasCA = true;
    if( newResNum != oldResNum || newICode != oldICode ) {
      resExtr.push_back( AtomVectCitPair( oldBeg, avit ) );
      resHasCA.push_back( bHasCA );
      resLbls.push_back( resLabel );
      nextChar = res3to1(avit->getResName(),false);
      sequence += nextChar;
      if( nextChar == 'X' )
        sequenceNonStd += res3to1(avit->getResName(),true);
      else
        sequenceNonStd += nextChar;
      oldResNum = newResNum;
      oldICode  = newICode;
      oldBeg    = avit;
      bHasCA    = false;
      resLabel  = Atom::getResLabel( oldResNum, oldICode );
    }
  }
  resExtr.push_back( AtomVectCitPair( oldBeg, endIt ) );
  resHasCA.push_back( bHasCA );
  resLbls.push_back( resLabel );

  DBG_VDUMP( sequence.size() );
  DBG_VDUMP( resExtr.size() );
  DBG_VDUMP( resLbls.size() );
  DBG_VDUMP( resLbls.back() );
}


bool Structure::Chain::isContinuous() const
{
  long chainNumResidues = numResidues();

  if( chainNumResidues <= 0 )
    return true;

  AtomVectCitPair vp0 = getResExtr( 0 );
  Point prevPos;
  bool bOk = false;
  for( AtomVectCit avit = vp0.beg; avit != vp0.end; ++avit )
    if( avit->getOrigName() == " CA " ) {
      prevPos = avit->getPos();
      bOk = true;
      break;
    }
  if( !bOk ) return false;

  for( long rj = 1; rj < chainNumResidues; ++rj ) {
    AtomVectCitPair vp = getResExtr( rj );

    Point pos;
    bOk = false;
    for( AtomVectCit avit = vp.beg; avit != vp.end; ++avit )
      if( avit->getOrigName() == " CA " ) {
        pos = avit->getPos();
        bOk = true;
        break;
      }
    if( !bOk ) return false;

    if( pointDistance( prevPos, pos ) > CA_CA_DIST_THRESHOLD )
      return false;

    prevPos = pos;
  }

  return true;
}


int Structure::Chain::isComplete( unsigned long resIdx ) const
{
  set<string> resAtoms;
  AtomVectCitPair avp = getResExtr(resIdx);
  for( AtomVectCit ai = avp.beg; ai != avp.end; ++ai ) {
    if( ai->getElemNo() != Elem_H )
      resAtoms.insert(ai->getName());
  }
  return isResComplete(avp.beg->getResName(),resAtoms);
}


bool Structure::Chain::hasOnlyMainChain( unsigned long resIdx ) const
{
  set<string> resAtoms;
  AtomVectCitPair avp = getResExtr(resIdx);
  for( AtomVectCit ai = avp.beg; ai != avp.end; ++ai ) {
    if( ai->getElemNo() != Elem_H )
      resAtoms.insert(ai->getName());
  }
  return hasResOnlyMainChain(resAtoms);
}


bool Structure::Chain::hasOnlyAltLoc( unsigned long resIdx ) const
{
  AtomVectCitPair avp = getResExtr(resIdx);
  if( (avp.end-avp.beg) <= 0 )
    return false;

  for( AtomVectCit ai = avp.beg; ai != avp.end; ++ai ) {
    if( ai->getAltLoc() == " " )
      return false;
  }

  return true;
}


Point Structure::Chain::getBaricenter() const
{
  Point baricenter;
  unsigned long numAtoms = 0UL;

  for( AtomVectCit avit = begin(); avit != end(); ++avit ) {
    baricenter += avit->getPos();
    ++numAtoms;
  }
  baricenter /= (double)numAtoms;

  return baricenter;
}


double Structure::Chain::calcAverageRadius() const
{
  double avrgRadius = 0.0;
  unsigned long numAtoms = 0UL;

  Point baricenter( getBaricenter() );

  for( AtomVectCit avit = begin(); avit != end(); ++avit ) {
    avrgRadius += pointDistance( avit->getPos(), baricenter );
    ++numAtoms;
  }

  avrgRadius /= (double)numAtoms;

  return avrgRadius;
}


AAModelVect Structure::Chain::buildAminoAcidModel() const
{
  DBG_BLOCK_SUPPRESS( "Structure::buildAminoAcidModel" );

  AAModelVect aaModelVect;

  DBG_MSG( "Mark1" );
  long chainResSize = numResidues();

  DBG_MSG( "Mark2" );
  for( long ri = 0; ri < chainResSize; ++ri ) {
    AtomVectCitPair resExtr = getResExtr(ri);
    aaModelVect.push_back( AminoAcidModel( resExtr.beg, resExtr.end ) );
  }

  DBG_MSG( "Mark3" );
  return aaModelVect;

}


static double findSurfIntersectingDistance( const AAModelVect &m, const Point &unitP )
{
  double maxD = 0.0;

  // First I need to modify the coordinates of the aaModel to the new center
  for( AAModelVectCit it = m.begin(); it != m.end(); ++it ) {
    double t1    = (unitP * (it->center));
    double delta = t1 * t1 - it->center.sqNorm() + (it->radius * it->radius);
    if( delta > 0 ) {
      //is.push_back(*it);
      double d1 = t1 + sqrt(delta);
      if( d1 > maxD )
        maxD = d1;
    }
  }
  return maxD;
}


static Point projectPointOnSurface( const AAModelVect &m, const Point &p,
                                    double radT )
{
  Point unitP( p.unit() );

  return unitP * (findSurfIntersectingDistance(m,unitP)+radT);
}


typedef list<Point> PointList;
typedef list<Point>::iterator PointListIt;
typedef list<Point>::const_iterator PointListCit;


static double pathLength( const PointList &pList )
{
  double lPathLength = 0.0;

  if( pList.empty() ) return lPathLength;

  PointListCit p1it = pList.begin();
  PointListCit p2it = pList.begin();
  ++p2it;
  for( ; p2it != pList.end(); ++p1it, ++p2it ) {
    lPathLength += pointDistance(*p1it,*p2it);
  }
  return lPathLength;
}


static void createPath( PointList &pList, const AAModelVect &m, double radT )
{
  if( pList.size() != 2 )
    throw logic_error( "list with less than 2 elements given to createPath" );

  static const double DIST_EPS = 1.0;
  static const long   MAX_ITERATIONS = 8;

  for( long iteration = 0L; iteration < MAX_ITERATIONS; ++iteration ) {
    bool bAdded = false;

    PointListIt p1it = pList.begin();
    PointListIt p2it = pList.begin();
    ++p2it;
    while( p2it != pList.end() ) {
      Point p1( *p1it );
      Point p2( *p2it );

      if( pointDistance( p1, p2 ) > DIST_EPS ) {
        bAdded = true;

        Point halfwayPoint( projectPointOnSurface( m, p1+p2, radT ) );

        pList.insert( p2it, halfwayPoint );
      }
      p1it = p2it++;
    }
    if( !bAdded )
      break;
  }
}


double Structure::Chain::getSurfaceDistance( const Point &p1,
                                             const Point &p2,
                                             double radT ) const
{
  DBG_BLOCK_SUPPRESS( "Structure::getSurfaceDistance" );
  AAModelVect aaModelVect;

  // First we need the amino acid model of the chain
  buildAminoAcidModel().swap(aaModelVect);

  // Then we need a set of points around the baricenter for the sampling
  Point baricenter = getBaricenter();
  double avrgRadius = calcAverageRadius();
  double oneThirdAvrgRadius = avrgRadius / 3.0;

  vector<Point> samplingPoints;

  samplingPoints.push_back( baricenter );
  samplingPoints.push_back( baricenter + Point(oneThirdAvrgRadius,0.0,0.0) );
  samplingPoints.push_back( baricenter - Point(oneThirdAvrgRadius,0.0,0.0) );
  samplingPoints.push_back( baricenter + Point(0.0,oneThirdAvrgRadius,0.0) );
  samplingPoints.push_back( baricenter - Point(0.0,oneThirdAvrgRadius,0.0) );
  samplingPoints.push_back( baricenter + Point(0.0,0.0,oneThirdAvrgRadius) );
  samplingPoints.push_back( baricenter - Point(0.0,0.0,oneThirdAvrgRadius) );

  for( int i = 0; i < NUM_SAMPLING_POINTS_SURF_DIST; ++i ) {
    samplingPoints.push_back( baricenter + Point( (2.0 * rand01() - 1.0) * oneThirdAvrgRadius,
                                                  (2.0 * rand01() - 1.0) * oneThirdAvrgRadius,
                                                  (2.0 * rand01() - 1.0) * oneThirdAvrgRadius ) );
  }

  /*****************************************************************************
   *                    BEGIN DEBUGGING PYMOL SCRIPT                           *
   *****************************************************************************
  ofstream outputScript( "surf_distance.pym" );

  outputScript << "from pymol.cgo import *" << endl;
  outputScript << "from pymol import cmd" << endl << endl;

  outputScript << "aaModels = [" << endl;
  outputScript << "  COLOR,    0.800,    0.800,    0.800";

  for( AAModelVectCit it = aaModelVect.begin(); it != aaModelVect.end(); ++it ) {
    const Point &cc = it->center;

    outputScript << "," << endl << "  SPHERE, " << right << fixed << setw(8) << setprecision(3) << cc.x()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << cc.y()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << cc.z()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << it->radius;
  }

  outputScript << endl << "]" << endl;
  outputScript << "cmd.load_cgo(aaModels, 'aa_models',   1)" << endl << endl;

  outputScript << "avrgRadius = [" << endl;
  outputScript << "  COLOR,    0.750,    1.000,    0.750";
  outputScript << "," << endl << "  SPHERE, " << right << fixed << setw(8) << setprecision(3) << baricenter.x()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << baricenter.y()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << baricenter.z()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << avrgRadius;
  outputScript << endl << "]" << endl;
  outputScript << "cmd.load_cgo(avrgRadius, 'avrg_radius',   1)" << endl << endl;

  outputScript << "endPoints = [" << endl;
  outputScript << "  COLOR,    0.100,    1.000,    0.100";
  outputScript << "," << endl << "  SPHERE, " << right << fixed << setw(8) << setprecision(3) << p1.x()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << p1.y()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << p1.z()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << 0.5;
  outputScript << "," << endl << "  SPHERE, " << right << fixed << setw(8) << setprecision(3) << p2.x()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << p2.y()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << p2.z()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << 0.5;
  outputScript << "," << endl << "  LINEWIDTH, 1.0," << endl;
  outputScript << "  BEGIN, LINES," << endl;
  outputScript << "  VERTEX, " << right << fixed << setw(8) << setprecision(3) << baricenter.x()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << baricenter.y()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << baricenter.z()     << ", " << endl;
  outputScript << "  VERTEX, " << right << fixed << setw(8) << setprecision(3) << p1.x()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << p1.y()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << p1.z()     << ", " << endl;
  outputScript << "  VERTEX, " << right << fixed << setw(8) << setprecision(3) << baricenter.x()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << baricenter.y()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << baricenter.z()     << ", " << endl;
  outputScript << "  VERTEX, " << right << fixed << setw(8) << setprecision(3) << p2.x()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << p2.y()     << ", "
                                              << right << fixed << setw(8) << setprecision(3) << p2.z()     << ", " << endl;
  outputScript << "  END" << endl;
  outputScript << "]" << endl;
  outputScript << "cmd.load_cgo(endPoints, 'end_points',   1)" << endl << endl;

  outputScript << "samplingPoints = [" << endl;
  outputScript << "  COLOR,    1.000,    0.100,    0.100";
  for( vector<Point>::const_iterator it = samplingPoints.begin(); it != samplingPoints.end(); ++it ) {
    outputScript << "," << endl << "  SPHERE, " << right << fixed << setw(8) << setprecision(3) << it->x()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << it->y()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << it->z()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << 0.5;
  }
  outputScript << endl << "]" << endl;
  outputScript << "cmd.load_cgo(samplingPoints, 'centers',   1)" << endl << endl;

   *****************************************************************************
   *                      END DEBUGGING PYMOL SCRIPT                           *
   *****************************************************************************/

  vector<double> distances;
  long distCounter = 0;
  for( vector<Point>::const_iterator cpi  = samplingPoints.begin(); cpi != samplingPoints.end(); ++cpi, ++distCounter ) {

    AAModelVect modAAModel(aaModelVect);

    // First I need to modify the coordinates of the aaModel to the new center
    for( AAModelVectIt it = modAAModel.begin(); it != modAAModel.end(); ++it ) {
      it->center -= *cpi;
    }

    // Then I need to identify the position of the two points on the surface of
    // the modified model
    Point mp1( p1 - *cpi );
    Point mp2( p2 - *cpi );

    //AAModelVect intersectingSpheres;
    Point pp1( projectPointOnSurface(modAAModel,mp1,radT) );
    Point pp2( projectPointOnSurface(modAAModel,mp2,radT) );

    list<Point> path;

    path.push_back(pp1);
    path.push_back(pp2);

    createPath(path,modAAModel,radT);

    distances.push_back( pathLength( path ) );

    /*****************************************************************************
     *                    BEGIN DEBUGGING PYMOL SCRIPT                           *
     *****************************************************************************

    outputScript << "projEndPoints" << distCounter << " = [" << endl;
    outputScript << "  COLOR,    0.100,    0.100,    1.000";
    outputScript << "," << endl << "  SPHERE, " << right << fixed << setw(8) << setprecision(3) << (pp1+*cpi).x()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << (pp1+*cpi).y()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << (pp1+*cpi).z()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << 0.5;
    outputScript << "," << endl << "  SPHERE, " << right << fixed << setw(8) << setprecision(3) << (pp2+*cpi).x()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << (pp2+*cpi).y()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << (pp2+*cpi).z()     << ", "
                                                << right << fixed << setw(8) << setprecision(3) << 0.5;
    outputScript << endl << "]" << endl;
    outputScript << "cmd.load_cgo(projEndPoints" << distCounter << ", 'surf_end_points" << distCounter << "',   1)" << endl << endl;

    outputScript << "path" << distCounter << " = [" << endl;
    outputScript << "  COLOR,  " << right << fixed << setw(8) << setprecision(3) << rand01()     << ", "
                                 << right << fixed << setw(8) << setprecision(3) << rand01()     << ", "
                                 << right << fixed << setw(8) << setprecision(3) << rand01()     << ", " << endl;
    outputScript << "  LINEWIDTH, 1.0," << endl;
    outputScript << "  BEGIN, LINE_STRIP," << endl;

    if( path.size() > 1 ) {

      for(PointListCit p1it = path.begin(); p1it != path.end(); ++p1it ) {
        Point pathP(*cpi+*p1it);
        outputScript << "  VERTEX, " << right << fixed << setw(8) << setprecision(3) << pathP.x()     << ", "
                                     << right << fixed << setw(8) << setprecision(3) << pathP.y()     << ", "
                                     << right << fixed << setw(8) << setprecision(3) << pathP.z()     << ", " << endl;
      }
    }
    outputScript << "  END" << endl;
    outputScript << "]" << endl;
    outputScript << "cmd.load_cgo(path" << distCounter << ", 'path_" << distCounter << "',   1)" << endl << endl;

    outputScript << "intSpheres = [" << endl;
    outputScript << "  COLOR,    1.000,  0.700,  0.200";

    for( AAModelVectCit it = intersectingSpheres.begin(); it != intersectingSpheres.end(); ++it ) {
      const Point &cc = it->center + *cpi;

      outputScript << "," << endl << "  SPHERE, " << right << fixed << setw(8) << setprecision(3) << cc.x()     << ", "
                                                  << right << fixed << setw(8) << setprecision(3) << cc.y()     << ", "
                                                  << right << fixed << setw(8) << setprecision(3) << cc.z()     << ", "
                                                  << right << fixed << setw(8) << setprecision(3) << it->radius;
    }

    outputScript << endl << "]" << endl;
    outputScript << "cmd.load_cgo(intSpheres, 'int_spheres',   1)" << endl << endl;

     *****************************************************************************
     *                      END DEBUGGING PYMOL SCRIPT                           *
     *****************************************************************************/

  }

  /*****************************************************************************
   *                    BEGIN DEBUGGING PYMOL SCRIPT                           *
   *****************************************************************************
  outputScript.close();
   *****************************************************************************
   *                      END DEBUGGING PYMOL SCRIPT                           *
   *****************************************************************************/

  return *min_element(distances.begin(),distances.end());
}


double rosa::superpose( const Structure &ref, Structure &toMove )
{
  Mat2D<double> rotMat;
  Point cm1, cm2;

  double rmsd = calcSuperposition( ref, toMove, rotMat, cm1, cm2 );

  toMove.translate( -cm2 );
  toMove.rotate( rotMat );
  toMove.translate( cm1 );

  return rmsd;
}


double rosa::superpose( const Structure &ref, const std::string &refSelString,
                        Structure &toMove, const std::string &toMoveSelString )
{
  StructurePtr selRef    = ref.select( refSelString );
  StructurePtr selToMove = toMove.select( toMoveSelString );

  Mat2D<double> rotMat;
  Point cm1, cm2;

  double rmsd = calcSuperposition( *selRef, *selToMove, rotMat, cm1, cm2 );

  toMove.translate( -cm2 );
  toMove.rotate( rotMat );
  toMove.translate( cm1 );

  return rmsd;
}


long rosa::reduceToCommonAtoms( Structure &s1, Structure &s2 )
{
  DBG_BLOCK_SUPPRESS( "reduceToCommonAtoms" );

  if( s1.numModels() != s2.numModels() )
    throw logic_error( "reduceToCommonAtoms: structures with different number "
                       "of models ("+ltos(s1.numModels())+" != "
                       +ltos(s2.numModels())+")" );

  for( unsigned long i = 0; i < s1.numModels(); ++i ) {
    DBG_MSG( "Processing model "+ltos(i) );
    AtomVect newV1, newV2;

    std::vector<Structure::Chain> chains1, chains2;
    s1.buildChains( chains1, i );
    s2.buildChains( chains2, i );

    if( chains1.size() != chains2.size() )
      throw logic_error( "reduceToCommonAtoms: structures with different number "
                         "of chains ("+ltos(chains1.size())+" != "+
                         ltos(chains2.size())+") in model "+ltos(i) );

    unsigned long chainSize = chains1.size();
    for( unsigned long j = 0; j < chainSize; ++j ) {
      DBG_MSG( "Processing chain "+chains1[j].getChainID() );
      Structure::Chain &c1 = chains1[j];
      Structure::Chain &c2 = chains2[j];

      // First for every residue in the second structure we insert its index in
      // a map
      map<string,long> resMap;

      for( long rt = 0; rt < c2.numResidues(); ++rt )
        resMap[c2.getResLabel(rt)] = rt;

      DBG_MSG( "Now rebuilding vector of atoms" );

      // Now we rebuild the vector of atoms following the order of atoms
      // in the first of the two structures
      for( long k = 0; k < c1.numResidues(); ++k ) {
        // If the residue is present in the second structure
        map<string,long>::const_iterator rit = resMap.find( c1.getResLabel(k) );
        if( rit != resMap.end() ) {
          unsigned long rk = rit->second;
          AtomVectCitPair vp1 = c1.getResExtr( k );
          AtomVectCitPair vp2 = c2.getResExtr( rk );

          for( AtomVectCit it = vp1.beg; it != vp1.end; ++it )
            for( AtomVectCit jt = vp2.beg; jt != vp2.end; ++jt )
              if( it->getOrigName() == jt->getOrigName() &&
                  it->getAltLoc() == jt->getAltLoc() ) {
                newV1.push_back( *it );
                newV2.push_back( *jt );
                break;
              }
        }
      }
    }
    swap( s1.models[i], newV1 );
    swap( s2.models[i], newV2 );
  }

  s1.invalidateChains();
  s2.invalidateChains();

  return 0;
}


long rosa::reduceToAlignedAtoms( Structure &s1, Structure &s2 )
{
  DBG_BLOCK_SUPPRESS( "reduceToAlignedAtoms" );

  if( s1.numModels() != s2.numModels() )
    throw logic_error( "reduceToAlignedAtoms: structures with different number "
                       "of models ("+ltos(s1.numModels())+" != "
                       +ltos(s2.numModels())+")" );

  for( unsigned long i = 0; i < s1.numModels(); ++i ) {
    DBG_MSG( "Processing model "+ltos(i) );
    AtomVect newV1, newV2;

    std::vector<Structure::Chain> chains1, chains2;
    s1.buildChains( chains1, i );
    s2.buildChains( chains2, i );

    if( chains1.size() != chains2.size() )
      throw logic_error( "reduceToAlignedAtoms: structures with different number "
                         "of chains ("+ltos(chains1.size())+" != "+
                         ltos(chains2.size())+") in model "+ltos(i) );

    unsigned long chainSize = chains1.size();
    for( unsigned long j = 0; j < chainSize; ++j ) {
      DBG_MSG( "Processing chain "+chains1[j].getChainID() );
      Structure::Chain &c1 = chains1[j];
      Structure::Chain &c2 = chains2[j];

      string seq1 = c1.getSequence();
      string seq2 = c2.getSequence();

      DBG_VDUMP(seq1);
      DBG_VDUMP(seq2);

      SeqAlignment seqAl;

      if( c1.getChainType() == Structure::Chain::NucleicAcidType) {
        DBG_MSG( "Aligning nucleic acid sequences")
        alignNucAcidSeq( seq1, seq2, seqAl );
      } else {
        DBG_MSG( "Aligning residue sequences")
        alignProteinSeq( seq1, seq2, seqAl );
      }

      DBG_MSG( "Now rebuilding vector of atoms" );

      DBG_VDUMP(seqAl.aln.size());

      // Now we rebuild the vector of atoms following the order of atoms
      // in the first of the two structures
      for( unsigned long k = 0; k < seqAl.aln.size(); ++k ) {
        AtomVectCitPair vp1 = c1.getResExtr( seqAl.aln[k].i1-1 );
        AtomVectCitPair vp2 = c2.getResExtr( seqAl.aln[k].i2-1 );

        for( AtomVectCit it = vp1.beg; it != vp1.end; ++it )
          for( AtomVectCit jt = vp2.beg; jt != vp2.end; ++jt ) {
            if( it->getOrigName() == jt->getOrigName() &&
                it->getAltLoc() == jt->getAltLoc() ) {
              DBG_VDUMP(it->getOrigName());
              DBG_VDUMP(jt->getOrigName());
              DBG_VDUMP(it->getAltLoc());
              DBG_VDUMP(jt->getAltLoc());
              newV1.push_back( *it );
              newV2.push_back( *jt );
              break;
            }
          }
      }
    }
    swap( s1.models[i], newV1 );
    swap( s2.models[i], newV2 );
  }

  s1.invalidateChains();
  s2.invalidateChains();

  return 0;
}


void rosa::printChainDescription( const Structure &s, ostream &os )
{
  DBG_BLOCK_SUPPRESS( "printChainDescription" );
  DBG_VDUMP( s.numChains() );

  for( unsigned long i = 0; i < s.numChains(); ++i ) {
    const Structure::Chain &c = s.getChain( i );
    os << "Chain " << c.getChainID() << " [" << c.numResidues() << " res.]: ";
    os << c.getDescription() << endl;
  }
}


StructurePtr rosa::selectInterface( const Structure &s1, const Structure &s2,
                                    float dist )
{
  DBG_BLOCK_SUPPRESS( "selectInterface" );

  StructurePtr result( new Structure() );

  if( s1.numModels() != s2.numModels() )
    throw logic_error( "selectInterface: structures with different number "
                       "of models ("+ltos(s1.numModels())+" != "
                       +ltos(s2.numModels())+")" );

  for( unsigned long i = 0; i < s1.numModels(); ++i ) {

    DBG_MSG( "Processing model "+ltos(i) );
    AtomVect newV;

    std::vector<Structure::Chain> chains1;
    std::vector<Structure::Chain> chains2;
    s1.buildChains( chains1, i );
    s2.buildChains( chains2, i );

    for( unsigned long ci = 0; ci < chains1.size(); ++ci ) {
      Structure::Chain &c1 = chains1[ci];
      DBG_MSG( "Processing chain "+c1.getChainID() );
      DBG_VDUMP( c1.numResidues() );
      // For every residue of the first chain
      for( long ri = 0; ri < c1.numResidues(); ++ri ) {
        AtomVectCitPair vp1 = c1.getResExtr( ri );

        bool toInclude = false;

        for( unsigned long cj = 0; cj < chains2.size() and !toInclude; ++cj ) {
          const Structure::Chain &c2 = chains2[cj];

          // For every residue of the second chain
          for( long rj = 0; rj < c2.numResidues() and !toInclude; ++rj ) {
            AtomVectCitPair vp2 = c2.getResExtr( rj );

            bool abort     = false;
            for( AtomVectCit it = vp1.beg; it != vp1.end and !toInclude and !abort; ++it ) {
              for( AtomVectCit jt = vp2.beg; jt != vp2.end; ++jt ) {
                double lDist = pointDistance( it->getPos(), jt->getPos() );
                if( lDist > dist + 25.0 ) {
                  abort = true;
                  break;
                }
                DBG_MSG( "distance("+it->identifier()+", "+jt->identifier()+") = "+dtos(lDist,-1,2) );
                if( lDist <= dist ) {
                  toInclude = true;
                  break;
                }
              }
            }
          }
        }
        if( toInclude )
          for( AtomVectCit it = vp1.beg; it != vp1.end; ++it )
            newV.push_back( *it );
      }
    }


    // For every residue of the second chain
    for( unsigned long cj = 0; cj < chains2.size(); ++cj ) {
      Structure::Chain &c2 = chains2[cj];
      DBG_MSG( "Processing chain "+c2.getChainID() );

      // For every residue of the first chain
      for( long rj = 0; rj < c2.numResidues(); ++rj ) {
        AtomVectCitPair vp2 = c2.getResExtr( rj );

        bool toInclude = false;

        for( unsigned long ci = 0; ci < chains1.size() and !toInclude; ++ci ) {
          const Structure::Chain &c1 = chains1[ci];

          // For every residue of the second chain
          for( long ri = 0; ri < c1.numResidues() and !toInclude; ++ri ) {
            AtomVectCitPair vp1 = c1.getResExtr( ri );

            bool abort     = false;
            for( AtomVectCit jt = vp2.beg; jt != vp2.end and !toInclude and !abort; ++jt ) {
              for( AtomVectCit it = vp1.beg; it != vp1.end; ++it ) {
                double lDist = pointDistance( it->getPos(), jt->getPos() );
                if( lDist > dist + 25.0 ) {
                  abort = true;
                  break;
                }
                if( lDist <= dist ) {
                  toInclude = true;
                  break;
                }
              }
            }
          }
        }
        if( toInclude )
          for( AtomVectCit it = vp2.beg; it != vp2.end; ++it )
            newV.push_back( *it );
      }
    }
    result->models.push_back(newV);
  }
  result->reprModel = s1.reprModel;

  return result;
}

StructurePtr rosa::selectBindingSite( const Structure &rec, const Structure &lig,
                                      float dist )
{
  DBG_BLOCK_SUPPRESS( "selectBindingSite" );

  StructurePtr result( new Structure() );

  if( rec.numModels() != lig.numModels() )
    throw logic_error( "selectInterface: structures with different number "
                       "of models ("+ltos(rec.numModels())+" != "
                       +ltos(lig.numModels())+")" );

  for( unsigned long i = 0; i < rec.numModels(); ++i ) {

    DBG_MSG( "Processing model "+ltos(i) );
    AtomVect newV;

    std::vector<Structure::Chain> chains1;
    std::vector<Structure::Chain> chains2;
    rec.buildChains( chains1, i );
    lig.buildChains( chains2, i );

    for( unsigned long ci = 0; ci < chains1.size(); ++ci ) {
      Structure::Chain &c1 = chains1[ci];
      DBG_MSG( "Processing chain "+c1.getChainID() );
      DBG_VDUMP( c1.numResidues() );
      // For every residue of the first chain
      for( long ri = 0; ri < c1.numResidues(); ++ri ) {
        AtomVectCitPair vp1 = c1.getResExtr( ri );

        bool toInclude = false;

        for( unsigned long cj = 0; cj < chains2.size() and !toInclude; ++cj ) {
          const Structure::Chain &c2 = chains2[cj];

          // For every residue of the second chain
          for( long rj = 0; rj < c2.numResidues() and !toInclude; ++rj ) {
            AtomVectCitPair vp2 = c2.getResExtr( rj );

            bool abort     = false;
            for( AtomVectCit it = vp1.beg; it != vp1.end and !toInclude and !abort; ++it ) {
              for( AtomVectCit jt = vp2.beg; jt != vp2.end; ++jt ) {
                double lDist = pointDistance( it->getPos(), jt->getPos() );
                if( lDist > dist + 25.0 ) {
                  abort = true;
                  break;
                }
                DBG_MSG( "distance("+it->identifier()+", "+jt->identifier()+") = "+dtos(lDist,-1,2) );
                if( lDist <= dist ) {
                  toInclude = true;
                  break;
                }
              }
            }
          }
        }
        if( toInclude )
          for( AtomVectCit it = vp1.beg; it != vp1.end; ++it )
            newV.push_back( *it );
      }
    }

    result->models.push_back(newV);
  }
  result->reprModel = rec.reprModel;

  return result;
}


StructurePtr rosa::mapAtoms( const Structure &subsetEq1, const Structure &eq1,
                             const Structure &eq2 )
{
  DBG_BLOCK_SUPPRESS( "mapAtoms" );

  StructurePtr result( new Structure( eq2 ) );

  result->models.clear();
  result->invalidateChains();

  if( subsetEq1.numModels() != eq1.numModels() or
      subsetEq1.numModels() != eq2.numModels() )
    throw logic_error( "mapAtoms: structures with different number "
                       "of models ("+ltos(subsetEq1.numModels())+", "
                       +ltos(eq1.numModels())+", "
                       +ltos(eq2.numModels())+")" );

  if( eq1.size() != eq2.size() )
    throw logic_error( "mapAtoms: structures with different number "
                       "of atoms ("+ltos(eq1.size())+" != "
                       +ltos(eq2.size())+")" );

  for( unsigned long i = 0; i < subsetEq1.numModels(); ++i ) {

    const AtomVect &v1 = subsetEq1.models[i];
    const AtomVect &v2 = eq1.models[i];
    const AtomVect &v3 = eq2.models[i];

    DBG_MSG( "Processing model "+ltos(i) );
    AtomVect newV;

    AtomVectCit it = v1.begin();
    AtomVectCit v1End = v1.end();
    AtomVectCit v2Beg = v2.begin();
    AtomVectCit v2End = v2.end();
    AtomVectCit jtStart = v2Beg;
    while( it != v1End ) {
      DBG_MSG( "Now searching for atom "+it->identifier() );
      DBG_MSG( "Starting from atom "+jtStart->identifier() );
      AtomVectCit jt = jtStart;
      bool found = false;
      while( jt != v2End )
        if( it->getSerialNum() != jt->getSerialNum() ) ++jt;
        else {
          found = true;
          break;
        }
      if( !found ) { // Let us try from the beginning
        jt = v2Beg;
        while( jt != jtStart )
          if( it->getSerialNum() != jt->getSerialNum() ) ++jt;
          else {
          found = true;
          break;
        }
      }
      if( found ) {
        DBG_MSG( "Atom found!" );
        newV.push_back( v3[jt-v2Beg] );
        jtStart = jt+1;
      }
      //else
      //  throw logic_error( "mapAtoms: first Structure is not an exact subset "
      //                     "of the second" );
      ++it;
    }
    result->models.push_back( newV );
  }

  return result;
}

unsigned long rosa::getResContacts( const Structure &s1, const Structure &s2,
                                    float dist,
                                    ResPairVect &resPairs )
{
  DBG_BLOCK_SUPPRESS( "getContacts" );

  unsigned long contactsCounter = 0;
  for( unsigned long ci = 0; ci < s1.numChains(); ++ci ) {
    const Structure::Chain &c1 = s1.getChain(ci);

    // For every residue of the first chain
    for( long ri = 0; ri < c1.numResidues(); ++ri ) {
      AtomVectCitPair vp1 = c1.getResExtr( ri );

      for( unsigned long cj = 0; cj < s2.numChains(); ++cj ) {
        const Structure::Chain &c2 = s2.getChain(cj);

        // For every residue of the second chain
        for( long rj = 0; rj < c2.numResidues(); ++rj ) {
          AtomVectCitPair vp2 = c2.getResExtr( rj );

          bool toInclude = false;
          bool abort     = false;
          DBG_MSG( "Residue pair: ("+ltos(ri)+c1.getChainID()+", "+ltos(rj)+c2.getChainID()+")" );
          for( AtomVectCit it = vp1.beg; it != vp1.end && !toInclude && !abort; ++it ) {
            for( AtomVectCit jt = vp2.beg; jt != vp2.end; ++jt ) {
              double lDist = pointDistance( it->getPos(), jt->getPos() );
              if( lDist > dist + 25.0 ) {
                abort = true;
                break;
              }
              DBG_MSG( "distance("+it->identifier()+", "+jt->identifier()+") = "+dtos(lDist,-1,2) );
              if( lDist <= dist ) {
                toInclude = true;
                break;
              }
            }
          }

          if( toInclude ) {
            ++contactsCounter;
            resPairs.push_back( pair<Residue,Residue>( Residue( vp1.beg->getChainID(),
                                                                vp1.beg->getResNum(),
                                                                vp1.beg->getInsCode(),
                                                                vp1.beg->getResName() ),
                                                       Residue( vp2.beg->getChainID(),
                                                                vp2.beg->getResNum(),
                                                                vp2.beg->getInsCode(),
                                                                vp2.beg->getResName() ) ) );
          }
        }
      }
    }
  }

  return contactsCounter;

}

unsigned long rosa::getResidueEquiv( const Structure &eq1, const Structure &eq2,
                                     map<Residue,Residue> &resEquivMap )
{
  DBG_BLOCK_SUPPRESS( "getResidueEquiv" );

  if( eq1.size() != eq2.size() )
    throw logic_error( "getResidueEquiv: structures with different number "
                       "of atoms ("+ltos(eq1.size())+" != "
                       +ltos(eq2.size())+")" );

  const AtomVect &av1 = eq1.models[eq1.reprModel];
  const AtomVect &av2 = eq2.models[eq2.reprModel];

  unsigned long mappingsCounter = 0;
  for( unsigned long ci = 0; ci < eq1.numChains(); ++ci ) {
    const Structure::Chain &c1 = eq1.getChain(ci);

    // For every residue of the first chain
    for( long ri = 0; ri < c1.numResidues(); ++ri ) {
      AtomVectCitPair vp1 = c1.getResExtr( ri );
      const Atom &a = *(vp1.beg);
      const Atom &b = av2[vp1.beg-av1.begin()];

      Residue r1(a.getChainID(), a.getResNum(), a.getInsCode(), a.getResName());
      Residue r2(b.getChainID(), b.getResNum(), b.getInsCode(), b.getResName());

      resEquivMap.insert( map<Residue,Residue>::value_type(r1,r2) );
      ++mappingsCounter;
    }
  }

  return mappingsCounter;
}

unsigned long rosa::countAtomContacts( const Structure &s1, const Structure &s2,
                                       float dist )
{
  DBG_BLOCK_SUPPRESS( "countAtomContacts" );
  DBG_VDUMP( dist );
  const AtomVect &av1 = s1.models[s1.reprModel];
  const AtomVect &av2 = s2.models[s2.reprModel];

  unsigned atomContactsCounter = 0;
  AtomVectCit av1End = av1.end();
  AtomVectCit av2End = av2.end();
  for( AtomVectCit it = av1.begin(); it != av1End; ++it ) {
    for( AtomVectCit jt = av2.begin(); jt != av2End; ++jt ) {
      if( pointDistance ( it->getPos(), jt->getPos() ) <= dist ) {
        DBG_MSG( it->identifier()+" - "+jt->identifier()+" = "+dtos(pointDistance ( it->getPos(), jt->getPos() ),-1,2) );
        ++atomContactsCounter;
      }
    }
  }

  return atomContactsCounter;
}

rosa::AminoAcidModel::AminoAcidModel( const AtomVectCit &beg,
                                      const AtomVectCit &end ):
  res( "", -1, "", "" ), center( 0.0, 0.0, 0.0 ), radius( 0.0 )
{
  if( beg != end ) {
    res.chainID = beg->getChainID();
    res.resNum  = beg->getResNum();
    res.iCode   = beg->getInsCode();
    res.resName = beg->getResName();

    double totalMass = 0.0;

    for( AtomVectCit it = beg; it != end; ++it ) {
      double atomMass = it->getMass();
      center += it->getPos() * atomMass;
      totalMass  += atomMass;
    }
    if( totalMass > 0.0 )
      center /= (double)totalMass;

    vector<double> distances;

    for( AtomVectCit it = beg; it != end; ++it ) {
      double atomRadius = it->getRadius();
      distances.push_back(pointDistance( it->getPos(), center ) + atomRadius);
    }

    if( ! distances.empty() )
      radius = accumulate(distances.begin(),distances.end(),0.0) / (double)(distances.size());
  }
}
