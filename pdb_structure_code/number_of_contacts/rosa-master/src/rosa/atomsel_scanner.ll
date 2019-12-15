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
%{ 
  
#include <rosa/atomsel_driver.h>
#include <rosa/atomsel_parser.h>
#include <rosa/atomsel_scanner.h>
#include <string>
#include <cstdio>

/* At a certain point it was needed by an old version of flex (2.5.4). Solved
   by installing newer version (actually there was an older version installed
   by fink that hid the new version. I simply removed the fink version.
#include <iostream>
using namespace std;
*/
  
/* import the parser's token type into a local typedef */
typedef rosa::AtomSelectParser::token token;
typedef rosa::AtomSelectParser::token_type token_type;
  
/* Work around an incompatibility in flex (at least versions 2.5.31 through
 * 2.5.33): it generates code that does not conform to C89.  See Debian bug
 * 333231 <http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=333231>.  */
#undef yywrap
#define yywrap()	1
  
/* This disables inclusion of unistd.h, which is not available under Visual C++
 * on Win32. The C++ scanner uses STL streams instead. */
#ifndef HAVE_UNISTD_H
#define YY_NO_UNISTD_H
#endif

/* By default yylex returns int, we use token_type. Unfortunately yyterminate
 * by default returns 0, which is not of token_type. */
#define yyterminate() return token::END

%}

/*** Flex Declarations and Options ***/

/* enable c++ scanner class generation */
%option c++

/* change the name of the scanner class. results in "ExampleFlexLexer" */
%option prefix="AtomSelect"

/* the manual says "somewhat more optimized" */
%option batch

/* no support for include files is planned */
%option noyywrap nounput 

/* define a few variables */

string        \"[^\n"]+\"
string2       [A-Za-z0-9$#_]+
ws            [ \t]+
endline       \n
alpha         [A-Za-z]
dig           [0-9]
int_num       [+-]?{dig}+
exponant      [eE]{int_num}
real_num1     {int_num}\.?({exponant})?
real_num2     [-+]?{dig}*\.{dig}+({exponant})?
real_num      {real_num1}|{real_num2}

%%

or						{ return token::OR; }
and           { return token::AND; }
not           { return token::NOT; }
all           { return token::ALL; }
name          { return token::NAME; }
elem          { return token::ELEM; }
resi          { return token::RESI; }
resn          { return token::RESN; }
occ           { return token::OCC; }
Bf            { return token::BFACT; }
ins           { return token::INS; }
alt           { return token::ALT; }
chain         { return token::CHAIN; }
{int_num}     { yylval->int_val  = atoi(YYText()); return token::INT_NUM; }
{real_num}    { yylval->real_val = atof(YYText()); return token::REAL_NUM; }
\<            { return token::LESS; }
\<=           { return token::LEQ; }
\>            { return token::GREATER; }
\>=           { return token::GREQ; }
=             { return token::EQUAL; }
\<\>          { return token::DIFF; }
{string}      { yylval->str = new std::string(YYText()+1,YYLeng()-2); return token::STRING; }
{string2}     { yylval->str = new std::string(YYText()); return token::STRING; }
{endline}     /* ignore end of line */
{ws}          /* ignore whitespace */
.             { return static_cast<token_type>(*YYText()); }

%% /*** Additional Code ***/

namespace rosa {

AtomSelectScanner::AtomSelectScanner( std::istream* in,
                                      std::ostream* out )
: AtomSelectFlexLexer(in, out)
{
}


AtomSelectScanner::~AtomSelectScanner()
{
}

#ifdef ACTIVATE_DEBUG
void AtomSelectScanner::set_debug( bool b )
{
  yy_flex_debug = b;
}
#endif

} // namespace rosa

/* This implementation of ExampleFlexLexer::yylex() is required to fill the
 * vtable of the class ExampleFlexLexer. We define the scanner's main yylex
 * function via YY_DECL to reside in the Scanner class instead. */

#ifdef yylex
#undef yylex
#endif


int AtomSelectFlexLexer::yylex()
{
  std::cerr << "...in AtomSelectFlexLexer::yylex() !" << std::endl;
  return 0;
}

