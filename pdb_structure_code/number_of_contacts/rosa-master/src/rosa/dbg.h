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

/*! \file dbg.h
 *  \brief Declares macros for inserting debugging code that can be later easily
 *         eliminated.
 */

#ifndef ROSA_DBG_H_
#define ROSA_DBG_H_

#include "common.h"
#include <iostream>

namespace rosa {
  
  class DebugObject {
  private:
    std::string scopeMsg;
    bool        bActivated;
    
  public:
    DebugObject( const std::string &aScopeMsg, bool bActivate ):
    scopeMsg( aScopeMsg ), bActivated( bActivate )
    {
      if( bActivated )
        std::cout << "DEB-ENTER " << scopeMsg << "..." << std::endl;
    }
    
    ~DebugObject()
    {
      if( bActivated )
        std::cout << "DEB-EXIT " << scopeMsg << "..." << std::endl;
    }
    
    template<typename T>
    void dumpVar( const std::string &varName, const T &val ) const
    {
      if( bActivated )
        std::cout << "DEB " << varName << ":" << val << std::endl; std::cout.flush();
    }
    
    void msg( const std::string &aMsg ) const
    {
      if( bActivated )
        std::cout << "DEB " << aMsg << std::endl;
    }
    
    const bool activated() const { return bActivated; }
  };
  
#ifdef ACTIVATE_DEBUG
  
  //! Activate debug code for the current block
#define DBG_BLOCK_ACTIVATE( name )   DebugObject _dbgObj__( (name), true );
  //! Activate debug code for the current block
#define DBG_BLOCK_SUPPRESS( name )   DebugObject _dbgObj__( (name), false );
  //! Prints a conditional debugging message
#define DBG_MSG( aMsg )              _dbgObj__.msg( (aMsg) );
  //! Dumps the content of a variable
#define DBG_VDUMP(a)                 _dbgObj__.dumpVar( (#a), (a) );
  //! Boolean variable indicating the debug code on
#define DBG_ON                       (_dbgObj__.activated())
  //! Boolean variable indicating the debug code off
#define DBG_OFF                      (!_dbgObj__.activated())

  // For the parser
#define ROSADEBUG 1
  
#else // ACTIVATE_DEBUG
  
#define DBG_BLOCK_ACTIVATE
#define DBG_BLOCK_SUPPRESS
#define DBG_MSG(a)
#define DBG_VDUMP(a)
#define DBG_ON                       (false)
#define DBG_OFF                      (false)
  
#define ROSADEBUG 0

#endif // ACTIVATE_DEBUG
  
} // namespace rosa

#endif // ROSA_DBG_H_
