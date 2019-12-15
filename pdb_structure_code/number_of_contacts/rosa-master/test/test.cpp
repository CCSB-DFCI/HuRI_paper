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

#include "test.h"
#include <iostream>
#include <typeinfo>
#include <cassert>
#include <stdexcept>

using namespace std;
using namespace TestToolbox;

// ******************* Test *******************
void Test::testCondition( bool cond, const string& condLabel,
                          const char* fname, long lineno ) {
  if( !cond )
    reportFailure( condLabel, fname, lineno );
  else
    reportSuccess();
}

void Test::reportFailure( const std::string& condLabel,
                          const char* fname, long lineno ) {
  ++numFailed;
  if( outStreamPtr ) {
    *outStreamPtr << typeid(*this).name() << "\t"
      << "Test failed: (" << condLabel << ") , " << fname
      << " (line " << lineno << ")" << endl;
  }
}

long Test::outputReport() const {
  if( outStreamPtr ) {
    *outStreamPtr << "Test \"" << typeid(*this).name()
      << "\":\n\tPassed: " << numSucceeded
      << "\tFailed: " << numFailed
      << endl;
  }
  return numFailed;
}


// ******************* TestGroup *******************
void TestGroup::addTest( Test *t )
{
  addTest( SharedTestPtr(t) );
}

void TestGroup::addTest( SharedTestPtr t )
{
  if(t.get() == 0)
    throw runtime_error( "Null test in TestGroup::addTest" );
  else if( outStreamPtr && !t->getOutStream())
    t->setOutStream( outStreamPtr );
  tests.push_back( t );
  t->reset();
}

void TestGroup::addTestGroup( const TestGroup& s ) {
  for(size_t i = 0; i < s.tests.size(); ++i) {
    assert( tests[i].get() !=0 );
    addTest( s.tests[i] );
  }
}

void TestGroup::clear() {
  tests.clear();
}

void TestGroup::run() {
  reset();
  for(size_t i = 0; i < tests.size(); ++i) {
    assert( tests[i].get() );
    tests[i]->run();
  }
}

long TestGroup::outputReport() const {
  if(outStreamPtr) {
    long totFail = 0;
    *outStreamPtr << "\nTest group \"" << title
    << "\"\n=======";
    size_t i;
    for(i = 0; i < title.size(); ++i)
      *outStreamPtr << '=';
    *outStreamPtr << "=" << endl;
    for(i = 0; i < tests.size(); ++i) {
      assert(tests[i].get());
      totFail += tests[i]->outputReport();
    }
    *outStreamPtr << "=======";
    for(i = 0; i < title.size(); ++i)
      *outStreamPtr << '=';
    *outStreamPtr << "=" << endl;
    return totFail;
  }
  else
    return getNumFailed();
}

long TestGroup::getNumSucceeded() const {
  long totPass = 0;
  for(size_t i = 0; i < tests.size(); ++i) {
    assert(tests[i].get());
    totPass += tests[i]->getNumSucceeded();
  }
  return totPass;
}

long TestGroup::getNumFailed() const {
  long totFail = 0;
  for(size_t i = 0; i < tests.size(); ++i) {
    assert(tests[i].get());
    totFail += tests[i]->getNumFailed();
  }
  return totFail;
}

void TestGroup::reset() {
  for(size_t i = 0; i < tests.size(); ++i) {
    assert(tests[i].get());
    tests[i]->reset();
  }
}
