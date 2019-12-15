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

/*! \file util.h
 *  \brief Contains general purpose routines to perform basic tasks.
 */

#ifndef ROSA_UTIL_H_
#define ROSA_UTIL_H_

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>

namespace rosa {
  
  //! Gets the basename of a file without the extension and the path information
  std::string baseFilename( const std::string &path );

  //! Gets the path from a filename omitting the base filename and extension
  std::string pathFilename( const std::string &path );

  //! Returns the filename without the extension
  std::string woExtFilename( const std::string &path );

  //! Returns the extension of the filename
  std::string extFilename( const std::string &filename );
  
  //! Joins the path and the filename
  std::string joinPath( const std::string &path, const std::string &filename );

  //! Eliminates leading and trailing white spaces in a string and returns the
  //! resulting string
  inline std::string trim( const std::string &s,
                           const std::string &charset = " \a\b\f\n\r\t\v" ) {
    if(s.length() == 0) return s;
    
    std::size_t beg = s.find_first_not_of( charset );
    std::size_t end = s.find_last_not_of( charset );
    
    if(beg == std::string::npos) return std::string();
    
    return std::string(s, beg, end - beg + 1);
  }  
  
  //! Replaces every occurrence of string 's' with string 'r' inside string src
  //! and returns the result
  inline std::string replace( const std::string &src, const std::string &s,
                             const std::string &r )
  {
    std::string result( src );
    std::string::size_type pos = result.find( s, 0 );
    
    while( pos != std::string::npos ) {
      result.replace( pos, s.size(), r );
      pos = result.find( s, pos+r.size() );
    }
    
    return result;
  }
  
  //! Splits a string into components using the separators specified in the
  //! argument. Fills the vector passed as the second argument with the
  //! components. The vector is not oerwritten, the components are just
  //! appended at the end.
  void split( const std::string& str, std::vector<std::string>& tokens,
              const std::string& delimiters = " \t" );
  
  //! Safely convert a string to a double. If the string doesn't contain a
  //! valid double the function throws an exception.
  double stod( const std::string &str );
  
  //! Safely convert a string to a long. If the string doesn't contain a
  //! valid long the function throws an exception.
  long stol( const std::string &str );

  //! Converts a long to a C++ string object
  inline std::string ltos( long ii, int width = -1, char fillChar = ' ' ) {
    std::ostringstream s;
    s.setf( std::ios::fixed, std::ios::floatfield );
    if( width > 0 ) s.width( width );
    if( fillChar != ' ' ) s.fill( fillChar );
    s << ii;
    return s.str();
  }

  //! Converts a long to a C++ string object in hexadecimal form
  inline std::string htos( long ii, int width = -1, char fillChar = ' ' ) {
    std::ostringstream s;
    if( width > 0 ) s.width( width );
    if( fillChar != ' ' ) s.fill( fillChar );
    s << std::hex << ii;
    return s.str();
  }
  
  
  //! Converts a double to a C++ string object
  inline std::string dtos( double rr, int width, int precision ) {
    std::ostringstream s;
    s.setf( std::ios::fixed, std::ios::floatfield );
    if( width >= 0 )     s.width(width);
    if( precision >= 0 ) s.precision(precision);
    s << rr;
    return s.str();
  }

  //! Returns true if the string passed as an argument contains a valid integer
  //! value, false otherwise
  inline bool isInteger( const std::string &aStr ) {
    std::string tmp( trim(aStr) );
    if( tmp.empty() ) return false;
    return (bool)(tmp.find_first_not_of("0123456789")==std::string::npos);
  }

  //! Prints a string splitting it in lines of fixed length
  void printInLines( const std::string &pText, int lineWidth,
                     std::ostream &os = std::cout );

  //! Reads a table in a CSV format. The delimiter used can be specified as the
  //! third argument and default to the TAB (\t) character. The tale is read in
  //! a vector of vector of strings.
  bool readTable( const std::string &aFilename,
                  std::vector< std::vector< std::string > > &outTable,
                  const std::string& delimiters = "\t" );

  //! Class used for parsing a file .ini with the usual .INI syntax. Values
  //! surrounded by " will be taken literally for the part sorrounded. Escaped
  //! \" will be unescaped if surrounded by ". Multiple sections
  //! with the same title are not allowed. Multiple entries with the same name
  //! are not allowed.
  class IniFile {
  public:
    //! Stores one section of the INI file
    class Section {
    private:
      typedef std::map<std::string,std::string> EntriesMap;
      typedef std::map<std::string,std::string>::iterator EntriesMapIt;
      //! Title of the section
      std::string title;
      //! Map of entries as couples (name,value)
      EntriesMap entries;
      //! Creates an empy section with a title
      Section( const std::string &aTitle ):
        title( aTitle )
      {}

      friend class IniFile;

    public:
      const std::string &getTitle() const { return title; }
      const std::string &operator [] ( const std::string &aName ) const;
    };
  private:
    typedef std::map<std::string,Section> SectionsMap;
        
    //! Map of sections. Two sections cannot share the same title
    SectionsMap sections;
    //! Filename
    std::string filename;
    
  public:
    //! Creates an empty IniFile
    IniFile() {}
    
    //! Opens and parses a file and loads the content in the current file
    void load( const std::string &aIniFilename );

    //! Saves the file. Now only for debugging purposes.
    void save( const std::string &aIniFilename ) const;

    //! Returns a reference to an object of type Section representing the
    //! found section. If the section is not present an exception id thrown
    const Section &section( const std::string &aTitle ) const;
    
  };
  
} // namespace rosa

#endif // ROSA_UTIL_H_
