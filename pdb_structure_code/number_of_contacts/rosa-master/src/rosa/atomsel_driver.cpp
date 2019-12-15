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

#include <rosa/atom.h>
#include <rosa/dbg.h>
#include <rosa/atomsel_driver.h>
#include <rosa/atomsel_scanner.h>
#include <rosa/atomsel_parser.h>
#include <sstream>
#include <stack>
#include <vector>

using namespace rosa;

AtomSelectDriver::AtomSelectDriver():
errorFlag( false )
#ifdef ACTIVATE_DEBUG
              , traceScanning( false ), traceParsing( false )
#endif
{
}

AtomSelectDriver::AtomSelectDriver( const std::string &input ):
errorFlag( false )
#ifdef ACTIVATE_DEBUG
              , traceScanning( false ), traceParsing( false )
#endif
{
  parseString( input );
}


bool AtomSelectDriver::parseString( const std::string& input )
{
  std::istringstream iss(input);
  
  //! Resets the error flag
  errorFlag = false;
  
  AtomSelectScanner lScanner(&iss);
#ifdef ACTIVATE_DEBUG
  lScanner.set_debug(traceScanning);
#endif

  ipnExpr.clear();
  
  AtomSelectParser parser( *this, lScanner );
#ifdef ACTIVATE_DEBUG
  parser.set_debug_level(traceParsing);
#endif
  errorFlag = (parser.parse() != 0);
  
  return !errorFlag;
}

bool AtomSelectDriver::evalOnAtom( const Atom &atom ) const
{
  DBG_BLOCK_SUPPRESS( "AtomSelectDriver::evalOnAtom" );
  std::stack<bool,std::vector<bool> > lStack;

  DBG_VDUMP( ipnExpr.size() );
  DBG_VDUMP( ipnExpr.front().op );
  DBG_VDUMP( ipnExpr.front().stringArg );
  
  if( ipnExpr.empty() )
    return true;
  
  for( IPNExprList::const_iterator it = ipnExpr.begin(); it != ipnExpr.end(); ++it ) {
    const IPNExprElem &el = (*it);
    switch( el.op ) {
      case IPNExprElem::TRUE_VAL:
        lStack.push( true );
        break;
      case IPNExprElem::NAME_EQ:
        lStack.push( atom.getName() == el.stringArg );
        break;
      case IPNExprElem::NAME_DIFF:
        lStack.push( atom.getName() != el.stringArg );
        break;
      case IPNExprElem::ELEM_EQ:
        lStack.push( atom.getElemSymbol() == el.stringArg );
        break;
      case IPNExprElem::ELEM_DIFF:
        lStack.push( atom.getElemSymbol() != el.stringArg );
        break;
      case IPNExprElem::RESI_RANGE: {
          long rNum = atom.getResNum();
          lStack.push( (rNum >= el.intArg1) && (rNum <= el.intArg2) );
        }
        break;
      case IPNExprElem::RESI_EQ:
        lStack.push( atom.getResNum() == el.intArg1 );
        break;
      case IPNExprElem::RESI_DIFF:
        lStack.push( atom.getResNum() != el.intArg1 );
        break;
      case IPNExprElem::RESI_LS:
        lStack.push( atom.getResNum() < el.intArg1 );
        break;
      case IPNExprElem::RESI_LEQ:
        lStack.push( atom.getResNum() <= el.intArg1 );
        break;
      case IPNExprElem::RESI_GR:
        lStack.push( atom.getResNum() > el.intArg1 );
        break;
      case IPNExprElem::RESI_GREQ:
        lStack.push( atom.getResNum() >= el.intArg1 );
        break;
      case IPNExprElem::OCC_EQ:
        lStack.push( atom.getOcc() == el.realArg );
        break;
      case IPNExprElem::OCC_DIFF:
        lStack.push( atom.getOcc() != el.realArg );
        break;
      case IPNExprElem::OCC_LS:
        lStack.push( atom.getOcc() < el.realArg );
        break;
      case IPNExprElem::OCC_LEQ:
        lStack.push( atom.getOcc() <= el.realArg );
        break;
      case IPNExprElem::OCC_GR:
        lStack.push( atom.getOcc() > el.realArg );
        break;
      case IPNExprElem::OCC_GREQ:
        lStack.push( atom.getOcc() >= el.realArg );
        break;
      case IPNExprElem::BFACT_EQ:
        lStack.push( atom.getBfac() == el.realArg );
        break;
      case IPNExprElem::BFACT_DIFF:
        lStack.push( atom.getBfac() != el.realArg );
        break;
      case IPNExprElem::BFACT_LS:
        lStack.push( atom.getBfac() < el.realArg );
        break;
      case IPNExprElem::BFACT_LEQ:
        lStack.push( atom.getBfac() <= el.realArg );
        break;
      case IPNExprElem::BFACT_GR:
        lStack.push( atom.getBfac() > el.realArg );
        break; 
      case IPNExprElem::BFACT_GREQ:
        lStack.push( atom.getBfac() >= el.realArg );
        break;
      case IPNExprElem::RESN_EQ:
        lStack.push( atom.getResName() == el.stringArg );
        break;
      case IPNExprElem::RESN_DIFF:
        lStack.push( atom.getResName() != el.stringArg );
        break;
      case IPNExprElem::INS_EQ:
        lStack.push( atom.getInsCode() == el.stringArg );
        break;
      case IPNExprElem::INS_DIFF:
        lStack.push( atom.getInsCode() != el.stringArg );
        break;
      case IPNExprElem::ALT_EQ:
        lStack.push( atom.getAltLoc() == el.stringArg );
        break;
      case IPNExprElem::ALT_DIFF:
        lStack.push( atom.getAltLoc() != el.stringArg );
        break;
      case IPNExprElem::CHAIN_EQ_STR:
        DBG_VDUMP( atom.getChainID() );
        DBG_VDUMP( el.stringArg );
        lStack.push( atom.getChainID() == el.stringArg );
        DBG_VDUMP( lStack.top() );
        break;
      case IPNExprElem::CHAIN_DIFF_STR:
        lStack.push( atom.getChainID() != el.stringArg );
        break;
      case IPNExprElem::CHAIN_EQ_NUM:
        DBG_VDUMP( atom.getChainID() );
        DBG_VDUMP( el.stringArg );
        lStack.push( atom.getChainID() == ltos(el.intArg1) );
        DBG_VDUMP( lStack.top() );
        break;
      case IPNExprElem::CHAIN_DIFF_NUM:
        lStack.push( atom.getChainID() != ltos(el.intArg1) );
        break;
      case IPNExprElem::OR: {
          bool lEl1 = lStack.top();
          lStack.pop();
          bool lEl2 = lStack.top();
          lStack.pop();
          lStack.push( lEl1 || lEl2 );
          //bool &lTop = lStack.top();
          //lTop = lTop || lEl1;
        }
        break;
      case IPNExprElem::AND: {
          bool lEl1 = lStack.top();
          lStack.pop();
          bool lEl2 = lStack.top();
          lStack.pop();
          lStack.push( lEl1 && lEl2 );
          //bool &lTop = lStack.top();
          //lTop = lTop && lEl1;
        }
        break;
      case IPNExprElem::NOT: {
          bool lEl = lStack.top();
          lStack.pop();
          lStack.push( !lEl );
          //bool &lTop = lStack.top();
          //lTop = !lTop;
        }
        break;
    }
  }
  
  return lStack.top();
}
