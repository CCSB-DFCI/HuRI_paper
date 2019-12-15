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

/*! \file atomsel_driver.h
 *  \brief Declares the class AtomSelectDriver which combines the lexer and the
 *         parser for atom selection strings.
 */

#ifndef ROSA_ATOMSEL_DRIVER_H_
#define ROSA_ATOMSEL_DRIVER_H_

#include <iostream>
#include <string>
#include <list>
#include "common.h"

namespace rosa {
  
  //! Description of an element of the inverse polish notation of a selection
  //! expression
  struct IPNExprElem {
    //! The operator can be one of the following
    enum OpType {
      TRUE_VAL       =  1,
      NAME_EQ        =  2,
      NAME_DIFF      =  3,
      ELEM_EQ        =  4,
      ELEM_DIFF      =  5,
      RESI_RANGE     =  6,
      RESI_EQ        =  7,
      RESI_DIFF      =  8,
      RESI_LS        =  9,
      RESI_LEQ       = 10,
      RESI_GR        = 11,
      RESI_GREQ      = 12,
      OCC_EQ         = 13,
      OCC_DIFF       = 14,
      OCC_LS         = 15,
      OCC_LEQ        = 16,
      OCC_GR         = 17,
      OCC_GREQ       = 18,
      BFACT_EQ       = 19,
      BFACT_DIFF     = 20,
      BFACT_LS       = 21,
      BFACT_LEQ      = 22,
      BFACT_GR       = 23,
      BFACT_GREQ     = 24,
      RESN_EQ        = 25,
      RESN_DIFF      = 26,
      INS_EQ         = 27,
      INS_DIFF       = 28,
      ALT_EQ         = 29,
      ALT_DIFF       = 30,
      CHAIN_EQ_STR   = 31,
      CHAIN_DIFF_STR = 32,
      CHAIN_EQ_NUM   = 33,
      CHAIN_DIFF_NUM = 34,
      OR             = 35,
      AND            = 36,
      NOT            = 37
    };
    
    //! Operator to be applied
    OpType      op;
    //! First integral argument
    int         intArg1;
    //! Second integral argument
    int         intArg2;
    //! Real argument
    double      realArg;
    //! Argument of type string
    std::string stringArg;
    
    //! For operations without arguments
    IPNExprElem( OpType aOp ):
    op( aOp ) {}
    //! For operations with one integral argument
    IPNExprElem( OpType aOp, int aIntArg1 ):
    op( aOp ), intArg1( aIntArg1 ) {}
    //! For operations with two integral arguments
    IPNExprElem( OpType aOp, int aIntArg1, int aIntArg2 ):
    op( aOp ), intArg1( aIntArg1 ), intArg2( aIntArg2 ) {}
    //! For operations with one real argument
    IPNExprElem( OpType aOp, double aRealArg ):
    op( aOp ), realArg( aRealArg ) {}
    //! For operations with one string argument
    IPNExprElem( OpType aOp, std::string aStringArg ):
    op( aOp ), stringArg( aStringArg ) {}
  };
  
  //! The list of operations in inverse polish notation is stored as an STL
  //! list of elements of type IPNExprElem
  typedef std::list<IPNExprElem> IPNExprList;
  
  //! Driver for the scanner/parser of atom selection strings
  class AtomSelectDriver {
  private:
    //! Pointer to the exaluated list of instructions
    IPNExprList ipnExpr;
    
    //! True if an error has occured while parsing the string
    bool errorFlag;
    
#ifdef ACTIVATE_DEBUG
    //! enable debug output in the flex scanner
    bool traceScanning;
    //! enable debug output in the bison parser
    bool traceParsing;
#endif
    
    //! Necessary for the parser to access the pointer to the scanner
    friend class AtomSelectParser;
    
  protected:
    
    /*! Error handling with associated line number. This can be modified to
     *  output the error e.g. to a dialog box. */
    virtual void error( std::string const& m )
    { errorFlag = true; }
    
  public:
    AtomSelectDriver();
    
    //! Constructor. Takes the string to parse as input.
    AtomSelectDriver( const std::string &input );
    
    //! Default destructor
    virtual ~AtomSelectDriver()
    {}
    
    //! Parse the string and evaluate it on the atom passed as the argument. The
    //! result of the evaluation is stored in the result variable.
    //! \param  input string to be evaluated
    //! \return false if the string is syntactically wrong, true otherwise
    bool parseString( const std::string &input );
    
    //! Evaluate the expression that has been parsed for the atom passed as the argument
    bool evalOnAtom( const class Atom &atom ) const;
    
    //! Return true if no error has occured, false otherwise
    bool good() const { return !errorFlag; }
    
#ifdef ACTIVATE_DEBUG
    //! Sets debug information on and off
    void setDebug( bool bFlag ) {
      traceScanning = bFlag;
      traceParsing = bFlag;
    }
#endif
    
  };
  
} // namespace rosa

#endif // ROSA_ATOMSEL_DRIVER_H_

