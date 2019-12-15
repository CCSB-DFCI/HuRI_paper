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

%{ /*** C/C++ Declarations ***/
#include <rosa/atomsel_driver.h>
#include <rosa/dbg.h>
#include <iostream>
#include <string>
%}

/*** yacc/bison Declarations ***/

/* Require bison 2.3 or later */
%require "2.3"

/* start symbol is named "start" */
%start input

/* use newer C++ skeleton file */
%skeleton "lalr1.cc"

/* namespace to enclose parser in */
/* %name-prefix="rosa" */
%define api.prefix {rosa}

/* set the parser's class identifier */
%define parser_class_name {AtomSelectParser}

/* The driver is passed by reference to the parser and to the scanner. This
 * provides a simple but effective pure interface, not relying on global
 * variables. */
%parse-param { class AtomSelectDriver&  driver }
%parse-param { class AtomSelectScanner& scanner }

/* verbose error messages */
%error-verbose

%union {
long        int_val;   //!< For returning integer numbers.
double      real_val;  //!< For returning real numbers.
std::string *str;      //!< For returning strings.
rosa::IPNExprList *ipn_expr; //!< For returning the inverse polish notation of an expression.
}

%token ALL ELEM NAME RESI RESN OCC BFACT INS ALT CHAIN LESS LEQ GREATER GREQ EQUAL DIFF OR AND NOT RANGE OPEN_PAR CLOSE_PAR
%token END        0        "end of file"
%token <int_val>  INT_NUM
%token <real_val> REAL_NUM
%token <str>      STRING
%type  <real_val> num
%type  <ipn_expr> exp
%type  <ipn_expr> input

%destructor { delete ($$); } STRING
%destructor { if( $$ ) delete ($$); } exp input

%left OR
%left AND
%left NOT     /* negation */


%{
  
#include <rosa/atomsel_driver.h>
#include <rosa/atomsel_scanner.h>

using namespace rosa;

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex scanner.lex
  
%}



%%

input     : /* empty */               { $$ = 0; }
          | exp                       { $$ = new IPNExprList( *$1 ); driver.ipnExpr = *$$; }
;

exp       : ALL												{ $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::TRUE_VAL) );                          }
          | RESI INT_NUM ':' INT_NUM  { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESI_RANGE, (int)($2), (int)($4)) );  }
          | RESI INT_NUM							{ $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESI_EQ, (int)($2)) );                }
          | RESI EQUAL INT_NUM  			{ $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESI_EQ, (int)($3)) );                }
          | RESI DIFF INT_NUM					{ $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESI_DIFF, (int)($3)) );              }
          | RESI LESS INT_NUM         { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESI_LS, (int)($3)) );                }
          | RESI LEQ  INT_NUM         { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESI_LEQ, (int)($3)) );               }
          | RESI GREATER INT_NUM      { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESI_GR, (int)($3)) );                }
          | RESI GREQ INT_NUM         { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESI_GREQ, (int)($3)) );              }
          | OCC EQUAL num             { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::OCC_EQ, (double)($3)) );              }
          | OCC DIFF num              { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::OCC_DIFF, (double)($3)) );            }
          | OCC LESS num              { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::OCC_LS, (double)($3)) );              }
          | OCC LEQ num               { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::OCC_LEQ, (double)($3)) );             }
          | OCC GREATER num           { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::OCC_GR, (double)($3)) );              }
          | OCC GREQ num              { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::OCC_GREQ, (double)($3)) );            }
          | BFACT EQUAL num           { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::BFACT_EQ, (double)($3)) );            }
          | BFACT DIFF num            { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::BFACT_DIFF, (double)($3)) );          }
          | BFACT LESS num            { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::BFACT_LS, (double)($3)) );            }
          | BFACT LEQ num             { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::BFACT_LEQ, (double)($3)) );           }
          | BFACT GREATER num         { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::BFACT_GR, (double)($3)) );            }
          | BFACT GREQ num            { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::BFACT_GREQ, (double)($3)) );          }
          | NAME STRING								{ $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::NAME_EQ, *($2)) );                    }
          | NAME EQUAL STRING					{ $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::NAME_EQ, *($3)) );                    }
          | NAME DIFF STRING					{ $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::NAME_DIFF, *($3)) );                  }
          | ELEM STRING               { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::ELEM_EQ, *($2)) );                    }
          | ELEM EQUAL STRING         { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::ELEM_EQ, *($3)) );                    }
          | ELEM DIFF STRING          { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::ELEM_DIFF, *($3)) );                  }
          | RESN STRING               { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESN_EQ, *($2)) );                    }
          | RESN EQUAL STRING         { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESN_EQ, *($3)) );                    }
          | RESN DIFF STRING          { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::RESN_DIFF, *($3)) );                  }
          | INS STRING                { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::INS_EQ, *($2)) );                     }
          | INS EQUAL STRING          { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::INS_EQ, *($3)) );                     }
          | INS DIFF STRING           { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::INS_DIFF, *($3)) );                   }
          | ALT STRING                { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::ALT_EQ, *($2)) );                     }
          | ALT EQUAL STRING          { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::ALT_EQ, *($3)) );                     }
          | ALT DIFF STRING           { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::ALT_DIFF, *($3)) );                   }
          | CHAIN STRING              { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::CHAIN_EQ_STR, *($2)) );               }
          | CHAIN EQUAL STRING        { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::CHAIN_EQ_STR, *($3)) );               }
          | CHAIN DIFF STRING         { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::CHAIN_DIFF_STR, *($3)) );             }
          | CHAIN INT_NUM             { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::CHAIN_EQ_NUM, (int)($2)) );               }
          | CHAIN EQUAL INT_NUM       { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::CHAIN_EQ_NUM, (int)($3)) );               }
          | CHAIN DIFF INT_NUM        { $$ = new IPNExprList( 1, IPNExprElem(IPNExprElem::CHAIN_DIFF_NUM, (int)($3)) );             }
          | exp OR exp                { $$ = new IPNExprList( *$1 ); $$->insert( $$->end(), $3->begin(), $3->end() ); $$->push_back( IPNExprElem(IPNExprElem::OR) );  }
          | exp AND exp               { $$ = new IPNExprList( *$1 ); $$->insert( $$->end(), $3->begin(), $3->end() ); $$->push_back( IPNExprElem(IPNExprElem::AND) );  }
          | NOT exp                   { $$ = new IPNExprList( *$2 ); $$->push_back( IPNExprElem(IPNExprElem::NOT) );            }
          | '(' exp ')'               { $$ = new IPNExprList( *$2 ); }
;

num       : INT_NUM                   { $$ = (double)$1; }
          | REAL_NUM                  { $$ = $1; }
;

%% /*** Additional Code ***/

void AtomSelectParser::error(std::string const& m)
{
  driver.error(m);
}
