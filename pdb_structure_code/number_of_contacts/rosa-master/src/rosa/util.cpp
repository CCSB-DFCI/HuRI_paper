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

#include <rosa/mem_shared_ptr.h>
#include <rosa/util.h>
#include <string>
#include <stdexcept>
#include <fstream>

using namespace std;
using namespace rosa;

namespace rosa {

string baseFilename( const string &filename )
{
  string::size_type slashPos = filename.rfind( '/' );
  string::size_type dotPos   = filename.rfind( '.' );
  string result;
  
  if( slashPos != string::npos ) {
  
    if( dotPos != string::npos && dotPos > slashPos )
      result = filename.substr( slashPos+1, dotPos - slashPos - 1 );
    else
      result = filename.substr( slashPos+1 );
  
  } else {
    
    if( dotPos != string::npos )
      result = filename.substr( 0, dotPos );
    else
      result = filename; 
  
  }
  
  return result;
}


string extFilename( const string &path )
{
  string::size_type slashPos = path.rfind( '/' );
  string::size_type dotPos   = path.rfind( '.' );
  string result;
  
  if( dotPos != string::npos &&
      (slashPos == string::npos || dotPos > slashPos) )
    result = path.substr( dotPos+1 );
  
  return result;
}


string pathFilename( const string &filename )
{
  string::size_type slashPos = filename.rfind( '/' );
  string result;
  
 if( slashPos != string::npos && slashPos != 0 )
   result = filename.substr( 0, slashPos-1 );
  
  return result;
}


string joinPath( const string &path, const string &filename )
{
  string result;
  
  if( !path.empty() ) {
    if( path[path.size()-1] != '/' )
      result = path + '/' + filename;
    else
      result = path + filename;
  } else
      result = filename;

  return result;
}

string woExtFilename( const string &path )
{
  string::size_type slashPos = path.rfind( '/' );
  string::size_type dotPos   = path.rfind( '.' );
  string result;
  
  if( dotPos == string::npos ||
      (slashPos != string::npos && dotPos < slashPos) )
    result = path;
  else
    result = path.substr( 0, dotPos );
  
  return result;
}


void split( const string& str, vector<string>& tokens, const string& delim )
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of( delim, 0 );
  // Find first "delimiter".
  string::size_type pos     = str.find_first_of( delim, lastPos );
  
  while( string::npos != pos || string::npos != lastPos ) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of( delim, pos );
    // Find next "delimiter"
    pos = str.find_first_of( delim, lastPos );
  }
}


double stod( const string &str )
{
  double result;
  istringstream inStr( str );
  inStr >> result;
  if( inStr.fail() )
    throw runtime_error( "not a valid floating point value" );
  return result;
}


long stol( const string &str )
{
  long result;
  istringstream inStr( str );
  inStr >> result;
  if( inStr.fail() )
    throw runtime_error( "not a valid long value" );
  return result;
}


bool readTable( const string &aFilename,
                vector< vector< string > > &outTable,
                const string &delimiters )
{
  ifstream inFile( aFilename.c_str() );
  
  if( inFile.fail() )
    throw runtime_error( string("error opening file ")+aFilename );
  
  outTable.clear();
  
  while( !inFile.eof() ) {
    string line;
    
    getline( inFile, line );
    
    if( !line.empty() ) {
        vector<string> newLine;
        split( line, newLine, delimiters );
        outTable.push_back( newLine );
    }
  }
  
  inFile.close();
  
  return true;
}


const string &IniFile::Section::operator [] ( const string &aName ) const
{
  map<string,string>::const_iterator it;
  if( ( it = entries.find( aName )) != entries.end() )
    return it->second;
  else
    throw runtime_error( title+": entry "+aName+" not found" );
}


void IniFile::load( const string &aIniFilename )
{
  ifstream iniFile( aIniFilename.c_str() );
  
  if( !iniFile.good() )
    throw runtime_error( "error while opening file "+aIniFilename );
  
  shared_ptr<Section> s;
  while( !iniFile.eof() ) {
    string line;
    if( getline( iniFile, line ) ) {
      line = trim( line );
      if( !line.empty() ) {
        if(line[0] == '[') {
          // New section
          if( line[line.size()-1] != ']' )
            throw runtime_error( aIniFilename+": error while parsing the file, "
                                 "wrong section title\n"+line );
          string title( line.substr( 1, line.size()-2 ) );
          if( sections.find( title ) != sections.end() )
            throw runtime_error( aIniFilename+": error while parsing the file, "
                                 "repeated section title\n"+title );
          if( s )
            sections.insert( SectionsMap::value_type( s->getTitle(), *s ) );
                            
          s.reset( new Section( title ) );
        } else {
          string::size_type eqPos = line.find( '=' );
          if( eqPos == string::npos || eqPos == 0 )
            throw runtime_error( aIniFilename+": error while parsing the file\n"+
                                 line );
          else { // Entry name=value
            string name = trim(line.substr( 0, eqPos ));
            string value = trim(line.substr( eqPos+1 ));
            if( !value.empty() )
              if( value[0] == '\"' && value[value.size()-1] == '\"' ) {
                /* Don't you think you're too picky????
                  if( value.find( "\"" ) != string::npos )
                  throw runtime_error( aIniFilename+": error while parsing the file, unescaped '\"' character\n"+
                                       line );
                */
                value = replace( value.substr( 1, value.size()-2 ), "\\\"", "\"" );
              }
            
            if( !s )
              throw runtime_error( aIniFilename+": error while parsing the file, "
                                   "entry outside section\n"+line );
            
            if( s->entries.find( name ) != s->entries.end() )
              throw runtime_error( aIniFilename+": error while parsing the file, "
                                   "repeated entry '"+name+"' in section"+s->getTitle()+"\n"+line );
            
            s->entries.insert( Section::EntriesMap::value_type( name, value ) );
          }
        }
      }
    }
  }
  
  if( s )
    sections.insert( SectionsMap::value_type( s->getTitle(), *s ) );

  filename = aIniFilename;
}


void IniFile::save( const string &aIniFilename ) const
{
  ofstream iniFile( aIniFilename.c_str() );
  
  if( !iniFile.good() )
    throw runtime_error( "error while opening file "+aIniFilename );
  
  SectionsMap::const_iterator sEnd = sections.end();
  for( SectionsMap::const_iterator it = sections.begin(); it != sEnd; ++it ) {
    const Section &sec = it->second;
    
    iniFile << "[" << sec.getTitle() << "]" << endl;
    Section::EntriesMap::const_iterator eEnd = sec.entries.end();
    for( Section::EntriesMap::const_iterator it = sec.entries.begin(); it != eEnd; ++it ) {
      string value( it->second );

      iniFile << it->first << "=";
      
      if( value.find( "\"" ) != string::npos or
          value.find( " " )  != string::npos or
          value.find( "\t" ) != string::npos or
          value.find( "\n" ) != string::npos or
          trim(value) != value  ) {
        value = replace( value, "\"", "\\\"" );
        iniFile << "\"" << value << "\"" << endl;
      } else {
        iniFile << value << endl;
      }
    }
    iniFile << endl;
  }
}


const IniFile::Section &IniFile::section( const string &aTitle ) const
{
  map<string,Section>::const_iterator it;
  
  if( (it = sections.find( aTitle )) != sections.end() )
    return it->second;
  else
    throw runtime_error( filename+": section "+aTitle+" not found" );
}

void printInLines( const string &pText, int lineWidth, ostream &os )
{
  int pTextLength = pText.size();
    
  for( int i = 0; i < pTextLength; i += lineWidth ) {
    if( i + lineWidth < pTextLength ) {
      os << pText.substr( i, lineWidth ) << endl;
    } else {
      os << pText.substr( i );
    }
  }
}

}
