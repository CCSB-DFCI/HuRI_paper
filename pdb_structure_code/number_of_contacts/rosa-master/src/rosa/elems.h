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

/*! \file elems.h
 *  \brief Contains structures and data related to the treatment of
 *         the chemical elements.
 */

#ifndef ROSA_ELEMS_H_
#define ROSA_ELEMS_H_

#include <cstring>

namespace rosa {
  
  //! Structure containing informations about one element
  struct Element {
    char  name[3];                  //!< short name of element
    float crad;                     //!< covalent radius in A
    float mass;                     //!< mass in Da
    float a1,b1,a2,b2,a3,b3,a4,b4;  //!< scattering factors
    float c;                        //!< ditto
    float vdw;                      //!< van der Waals radius in A
    //!< negative if no vdw radius available
  };
  
  //! Atom element types
  typedef enum {
    Elem_Unknown =  0,
    Elem_H,
    Elem_HE,
    Elem_LI,
    Elem_BE,
    Elem_B,
    Elem_C,
    Elem_N,
    Elem_O,
    Elem_F,
    Elem_NE,
    Elem_NA,
    Elem_MG,
    Elem_AL,
    Elem_SI,
    Elem_P,
    Elem_S,
    Elem_CL,
    Elem_AR,
    Elem_K,
    Elem_CA,
    Elem_SC,
    Elem_TI,
    Elem_V,
    Elem_CR,
    Elem_MN,
    Elem_FE,
    Elem_CO,
    Elem_NI,
    Elem_CU,
    Elem_ZN,
    Elem_GA,
    Elem_GE,
    Elem_AS,
    Elem_SE,
    Elem_BR,
    Elem_KR,
    Elem_RB,
    Elem_SR,
    Elem_Y,
    Elem_ZR,
    Elem_NB,
    Elem_MO,
    Elem_TC,
    Elem_RU,
    Elem_RH,
    Elem_PD,
    Elem_AG,
    Elem_CD,
    Elem_IN,
    Elem_SN,
    Elem_SB,
    Elem_TE,
    Elem_I,
    Elem_XE,
    Elem_CS,
    Elem_BA,
    Elem_LA,
    Elem_CE,
    Elem_PR,
    Elem_ND,
    Elem_PM,
    Elem_SM,
    Elem_EU,
    Elem_GD,
    Elem_TB,
    Elem_DY,
    Elem_HO,
    Elem_ER,
    Elem_TM,
    Elem_YB,
    Elem_LU,
    Elem_HF,
    Elem_TA,
    Elem_W,
    Elem_RE,
    Elem_OS,
    Elem_IR,
    Elem_PT,
    Elem_AU,
    Elem_HG,
    Elem_TL,
    Elem_PB,
    Elem_BI,
    Elem_PO,
    Elem_AT,
    Elem_RN,
    Elem_FR,
    Elem_RA,
    Elem_AC,
    Elem_TH,
    Elem_PA,
    Elem_U,
    Elem_NP,
    Elem_PU,
    ElementsSize
  } ElemType;
  
  //! Structure for keeping information about chemical elements 
  extern const Element Elements[ElementsSize];
  
  
  inline int elementIndex( const char *elemName ) {
    for( int i = 0; i < ElementsSize; ++i )
      if( std::strncmp( elemName, Elements[i].name, 2 ) == 0 )
        return i;
    return Elem_Unknown;
  }
  
  
  inline const char *elementName( ElemType elIndex ) {
    return Elements[elIndex].name;
  }
  
} // namespace rosa

#endif // ROSA_ELEMS_H_
