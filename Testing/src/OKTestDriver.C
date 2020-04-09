///
/// @file 
/// @ingroup plascom2_group
/// @brief Implements a command-line interface for the PlasCom2 tests.
///
/// Note that in this file, the "main" function is separated off from 
/// the program implementation.  This is done to make the documentation
/// usable and should be done for every new program in %PlasCom2. 
///
#include "ComLine.H"
#include "TestPlasCom2.H"

namespace plascom2 {
  
  
  ///
  /// Drives the plascom2::TestObject. 
  /// 
  /// @param argc number of string command line tokens
  /// @param argv string command line tokens
  /// @returns 0 if successful, 1 otherwise
  ///
  /// Drives the PlasCom2::TestObject, which should encapsulate
  /// all the tests for the PlasCom2 namespace (and thus the project).
  /// 
  /// Command line documentation:
  ///
  ///           plascom2_test [-h] [-v [level] -o <filename> -l <filename> -n <TestName> ] 
  ///
  ///           -h,--help
  ///              Print out long version of help and exit.
  ///
  ///           -v,--verblevel [level]
  ///              Set the verbosity level. (default = 0)
  ///
  ///           -o,--output <filename>
  ///              Set the output file to <filename>. (default = stdout)
  ///
  ///           -l,--list <filename>
  ///              Set the list file name to <filename>. (no default). The list file should be a text file with one test name per line.
  ///
  ///           -n,--name <TestName>
  ///              Run test by name. (no default)
  ///
  ///           -L,--List_Test_or_Suite <TestName or SuiteName>.
  ///              print the suite-name of a specific test, or print all the test for the given suitename.(no default)

  int Test(int argc,char *argv[])
  {
    
    // The default verbosity is 0
    int verblevel = 0;

    // This line creates the PlasCom2::TestComLine object and passes in 
    // the command line arguments as a (const char **).
    TestComLine comline((const char **)(argv));
    // The call to comline.Initialize() reads the command line arguments
    // from the array passed in the previous line.
    comline.Initialize();
    // The ProcessOptions() call does detailed examination of the command
    // line arguments to check for user errors or other problems. This call
    // will return non-zero if there were errors on the commandline.
    int clerr = comline.ProcessOptions();
    // Check if the user just wanted to get the help and exit
    if(!comline.GetOption("help").empty()){
      // Print out the "long usage" (i.e. help) message to stdout
      std::cout << comline.LongUsage() << std::endl;
      if(verblevel > 2)
        std::cout << "PlasCom2::Test: Exiting test function (success)" << std::endl;
      return(0);
    }
    if(clerr){
      std::cout << comline.ErrorReport() << std::endl
                << std::endl << comline.ShortUsage() << std::endl;
      if(verblevel > 2)
        std::cout << "PlasCom2::Test: Exiting test function (fail)" << std::endl;
      return(1);
    }
    // These outstreams allow the output to file to be set up and separated
    // from the stdout.
    std::ofstream Ouf;
    std::ostream *Out = &std::cout;
    std::ostringstream outStream; // if a file is specified, buffer output

    // The next few lines populate some strings based on the 
    // users input from the commandline.
    std::string OutFileName(comline.GetOption("output"));
    std::string TestName(comline.GetOption("name"));
    std::string ListName(comline.GetOption("list"));
    std::string sverb(comline.GetOption("verblevel"));
    std::string printlist(comline.GetOption("List_Suite_or_Test"));
    // The following block parses and sets the verbosity level
    if(!sverb.empty()){
      verblevel = 1;
      if(sverb != ".true."){
        std::istringstream Istr(sverb);
        Istr >> verblevel;
        if(verblevel < 0)
          verblevel = 1;
      }
    }
    
    // This block sets up the output file if the user specified one
    if(!OutFileName.empty()){
      Out = &outStream;
    }
    
    if(verblevel > 2)
      std::cout << "PlasCom2::Test: Entering test function" << std::endl;
    
    // Make an instance of the PlasCom2 testing object, PlasCom2::TestingObject
    // plascom2::TestingObject<plascom2::TestResults> test;
    // Make an instance of the PlasCom2 results object, PlasCom2::TestResults
    plascom2::TestingObject<plascom2::TestResults> test;
    plascom2::TestResults results;

//     if (!printlist.empty()){
//       test.PrintTest(printlist);
//       return 0;
//     }
    // If the user specified a name, then run only the named test
    if(!TestName.empty()){
      // This call runs a test by name
      test.RunTest(TestName,results);
    }
    // Otherwise, if the user specified a list, then read the list and
    // run the listed tests.
    else if(!ListName.empty()){
      std::ifstream ListInf;
      ListInf.open(ListName.c_str());
      if(!ListInf){
        std::cout << "PlasCom2::Test> Error: Could not open list of tests in file " 
                  << ListName << ". Exiting (fail)." << std::endl;
        return(1);
      }
      std::string testname;
      while(std::getline(ListInf,testname))
        test.RunTest(testname,results);
      ListInf.close();
    }
    else {
      // This call runs all the tests for the PlasCom2 namespace.
      test.Process(results);
    }

    *Out << results << std::endl;

    if(verblevel > 2)
      *Out << "PlasCom2::Test: Exiting test function (success)" << std::endl;    
    
    // Don't open until the end - race condition resolution
    if(!OutFileName.empty()){
      Ouf.open(OutFileName.c_str(),std::ios::app);
      if(!Ouf){
        std::cout << "PlasCom2::Test> Error: Could not open output file, " 
                  <<  OutFileName << " for test output. Exiting (fail)." << std::endl;
        return(1);
      }
      Ouf << outStream.str();
      Ouf.flush();
      Ouf.close();
    } else {
      //       int numTests = results.size();
      //       Parameter::iterator pmit = results.begin();
      //       while(pmit != results.end()){
      //         if(pmit->second == 0){ // fail
      //         } else{ // pass
      //         }
      //           pmit++;
      //       }
      //       // Count up and report the percent pass/fail.
      //       *Out << "Results report:" << std::endl;
      //       // *Out << THE REST << std::endl;
    }
    // if(Ouf)
    //   Ouf.close();

    return(0);
  }
};

int main(int argc,char *argv[])
{
  return(plascom2::Test(argc,argv));
}
