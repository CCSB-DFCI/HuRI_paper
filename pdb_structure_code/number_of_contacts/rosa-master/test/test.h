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

#ifndef TEST_H__
#define TEST_H__

#include <rosa/mem_shared_ptr.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#define SETUP_TEST(cond) \
testCondition((cond), #cond, __FILE__, __LINE__)

#define SETUP_TEST_LBL(cond,label) \
testCondition((cond), (label), __FILE__, __LINE__)

#define SETUP_TEST_EQ_LBL(expr1,expr2,label) \
testEqualityCondition((expr1), (expr2), (label), __FILE__, __LINE__)

#define SETUP_TEST_LEQ_LBL(expr1,expr2,label) \
testLeqCondition((expr1), (expr2), (label), __FILE__, __LINE__)

#define SETUP_TEST_GEQ_LBL(expr1,expr2,label) \
testGeqCondition((expr1), (expr2), (label), __FILE__, __LINE__)

#define SETUP_TEST_EXEC(expr,label) \
{ try { expr; } catch(...) { reportFailure( label + ": " + #expr, __FILE__, __LINE__ ); return; }; reportSuccess(); }

#define SETUP_TEST_EXCPT(expr,label) \
{ try { expr; reportFailure( label + ": " + #expr, __FILE__, __LINE__ ); } catch(...) { reportSuccess(); } }

namespace TestToolbox {
  
  //! Tests a series of related conditions.
  class Test {
    //! Output stream
    std::ostream* outStreamPtr;
    //! Number of test conditions that succeded
    long numSucceeded;
    //! Number of test conditions that failed
    long numFailed;
    //! Avoid copy, assignment and pass by value
    Test( const Test& );
    Test& operator=( const Test& );
  protected:
    void testCondition( bool cond, const std::string& condLabel,
                        const char* fname, long lineno );
    template<typename T>
    void testEqualityCondition( const T &expr1, const T &expr2,
                                const std::string& condLabel,
                                const char* fname, long lineno ) {
      if( expr1 != expr2 ) {
        std::ostringstream outStream;
        outStream << condLabel << " [" << expr1 << "!=" << expr2 << "]";
        reportFailure( outStream.str(), fname, lineno );
      } else
        reportSuccess();
    }

    template<typename T>
    void testLeqCondition( const T &expr1, const T &expr2,
                           const std::string& condLabel,
                           const char* fname, long lineno ) {
      if( expr1 > expr2 ) {
        std::ostringstream outStream;
        outStream << condLabel << " [" << expr1 << ">" << expr2 << "]";
        reportFailure( outStream.str(), fname, lineno );
      } else
        reportSuccess();
    }

    template<typename T>
    void testGeqCondition( const T &expr1, const T &expr2,
                           const std::string& condLabel,
                           const char* fname, long lineno ) {
      if( expr1 < expr2 ) {
        std::ostringstream outStream;
        outStream << condLabel << " [" << expr1 << ">" << expr2 << "]";
        reportFailure( outStream.str(), fname, lineno );
      } else
        reportSuccess();
    }

    void reportFailure( const std::string& condLabel,
                        const char* fname, long lineno );
    void reportSuccess() { ++numSucceeded; }
  public:
    Test( std::ostream* aOsPtr = &std::cout ):
      outStreamPtr( aOsPtr ), numSucceeded( 0 ), numFailed( 0 )
    {}

    virtual ~Test() {}
    
    //! Run the tests. Should be overridden by derived classes to instantiate
    //! the actual tests to be performed
    virtual void run() = 0;
    
    //! Returns the number of successes
    long getNumSucceeded() const { return numSucceeded; }
    //! Returns the number of failures
    long getNumFailed()    const { return numFailed; }
    
    //! Returns the internal pointer to the output stream
    const std::ostream* getOutStream() const  { return outStreamPtr; }
    //! Sets the internal pointer to the output stream
    void setOutStream( std::ostream* aOsPtr ) { outStreamPtr = aOsPtr; }
    
    //! Writes a short report of the results
    long outputReport() const;

    //! Resets the tests
    virtual void reset() { numSucceeded = numFailed = 0; }
  };
  
  //! Vector of test pointers
  typedef shared_ptr<Test> SharedTestPtr;
  typedef std::vector<SharedTestPtr> TestPtrVect;
  
  //! Tests a group of related tests.
  class TestGroup {
    //! Title for the test group
    std::string title;
    //! Output stream
    std::ostream* outStreamPtr;
    //! Vector of tests
    TestPtrVect tests;
    //! Resets the tests
    void reset();
    // Avoid copy, assignment and pass by value
    TestGroup( const TestGroup& );
    TestGroup& operator=( const TestGroup& );
  public:
    TestGroup( const std::string& aTitle, std::ostream* aOsPtr = &std::cout ):
      title( aTitle ), outStreamPtr( aOsPtr )
    {}
    
    //! Returns the title
    std::string getTitle() const { return title; }
    
    //! Returns the number of successes
    long getNumSucceeded() const;
    //! Returns the number of failures
    long getNumFailed()    const;

    //! Returns the internal pointer to the output stream
    const std::ostream* getOutStream() const  { return outStreamPtr; }
    //! Sets the internal pointer to the output stream
    void setOutStream( std::ostream* aOsPtr ) { outStreamPtr = aOsPtr; }
    
    //! Writes a short report of the results
    long outputReport() const;

    //! Adds a test to the group
    void addTest( Test *t );
    //! Adds a test to the group
    void addTest( SharedTestPtr t );
    //! Adds a group of tests to the group
    void addTestGroup( const TestGroup &tg);
    
    //! Runs all the tests in the group by calling the Test::run() function
    void run();
    
    //! Frees all the tests
    void clear();
  };
  
} // namespace TestToolbox

#endif // TEST_H__
