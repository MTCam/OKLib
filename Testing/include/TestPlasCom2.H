///
/// @file
/// @ingroup plascom2_group
/// @brief Testing utilities for PlasCom2.
///
#ifndef __TEST_PC2_H__
#define __TEST_PC2_H__

#include "PlasCom2Test.H"

namespace plascom2 {
  
  ///
  /// ComLineObject for PlasCom2 testing command-line interface.
  ///
  class TestComLine : public ix::util::ComLineObject
  {
  public:
    TestComLine()
      : ix::util::ComLineObject()
    {};
    TestComLine(const char *args[])
      : ix::util::ComLineObject(args)
    {};
    void Initialize(){
      AddOption('h',"help");
      AddHelp("help","Print out long version of help and exit.");
      AddOption('v',"verblevel",1,"level");
      AddHelp("verblevel","Set the verbosity level. (default = 0)");
      AddOption('o',"output",2,"filename");
      AddHelp("output","Set the output file to <filename>. (default = stdout)");
      AddOption('l',"list",2,"filename");
      
      AddOption('L',"List_Suite_or_Test", 1, "Testname or SuiteName");
      //AddArgument("List_Suite_or_Test",0);
      AddHelp("List Suite or Test", "Print out the SuiteName for Test or all tests for SuiteName");
      
      std::ostringstream Ostr;
      Ostr << "Set the list file name to <filename>. (no default). "
           << "The list file should be a text file with one test name per line.";
      AddHelp("list",Ostr.str());
      AddOption('n',"name",2,"TestName");
      AddHelp("name","Run test by name. (no default)");
      Ostr.str("");
      Ostr << "Command-line interface for the test of the PlasCom2.";
      _description.assign(Ostr.str());
    };
  };
};
#endif
