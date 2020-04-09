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
#include "COMM.H"


namespace plascom2 {
  
  ///
  /// Convenience typedef for CommunicatorObject 
  ///
  typedef ix::comm::CommunicatorObject CommType;
  ///
  /// Drives the PlasCom2::TestObject. 
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

  ///   -L,--List_Test_or_Suite <TestName or SuiteName>.
  ///        print the suite-name of a specific test, or print all the test for the given suitename.(no default)
  int ParallelTest(int argc,char *argv[])
  {
    
    // The default verbosity is 0
    int verblevel = 0;

    // This sets everything up with MPI
    MPI_Comm myComm = MPI_COMM_WORLD;
    plascom2::CommType communicator(myComm);
    int rank  = communicator.Rank();
    int nproc = communicator.Size();
    bool do_stdout = !rank;
    
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
      if(do_stdout){
        std::cout << comline.LongUsage() << std::endl;
        if(verblevel > 1)
          std::cout << "PlasCom2::ParallelTest: Exiting test function (success)" 
                    << std::endl;
      }
      communicator.SetExit(1);
    }
    if(communicator.Check())
      return(1);
    if(clerr){
      if(do_stdout){
        std::cout << comline.ErrorReport() << std::endl
                  << std::endl << comline.ShortUsage() << std::endl;
        if(verblevel > 2)
          std::cout << "PlasCom2::ParallelTest: Exiting test function (fail)" << std::endl;
      }
      communicator.SetExit(1);
    }
    if(communicator.Check())
      return(2);

    // These outstreams allow the output to file to be set up and separated
    // from the stdout.
    std::ofstream Ouf;
    std::ostream *Out = NULL;
    if(do_stdout)
      Out = &std::cout;

    // The next few lines populate some strings based on the 
    // users input from the commandline.
    std::string OutFileName(comline.GetOption("output"));
    std::string TestName(comline.GetOption("name"));
    std::string ListName(comline.GetOption("list"));
    std::string sverb(comline.GetOption("verblevel"));
    std::string printtest(comline.GetOption("List_Suite_or_Test"));
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
    std::ostringstream outStream;
    if(!OutFileName.empty()){
      if(Out)
  Out = &outStream;
    }
    
    if(verblevel > 1 && Out)
      *Out << "PlasCom2::ParallelTest: Entering test function" << std::endl;
    
    // Make an instance of the PlasCom2 testing object, PlasCom2::ParallelTestingObject
    plascom2::ParallelTestingObject<plascom2::CommType,plascom2::TestResults> test(communicator);
    // Make an instance of the PlasCom2 results object, PlasCom2::TestResults
    plascom2::TestResults results;
    
    // If the user specified a name, then run only the named test
    //     if (!printtest.empty()){
    //       test.PrintTest(printtest);
    //       return 0;
    //     }
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
        if(Out)
          *Out << "PlasCom2::ParallelTest> Error: Could not open list of tests in file " 
               << ListName << ". Exiting (fail)." << std::endl;
        communicator.SetExit(1);
      }
      if(communicator.Check())
        return(3);
      std::string testname;
      while(std::getline(ListInf,testname))
        test.RunTest(testname,results);
      ListInf.close();
    }
    else {
      // This call runs all the tests for the PlasCom2 namespace.
      test.Process(results);
    }
    if(Out)
      *Out << results << std::endl;
    
    // if(Out && Ouf)
    //   Ouf.close();

    if((verblevel > 1) && Out)
      *Out << "PlasCom2::ParallelTest: Exiting test function (success)" << std::endl;
    
    if(!OutFileName.empty()){
      if(Out){
        Ouf.open(OutFileName.c_str(),std::ios::app);
        if(!Ouf){
          std::cout << "PlasCom2::ParallelTest> Error: Could not open output file, " 
                    <<  OutFileName << " for test output. Exiting (fail)." << std::endl;
          communicator.SetExit(1);
        } else {
    Ouf << outStream.str();
    Ouf.close();
  }
      }
    }
    if(communicator.Check())
      return(4);
    communicator.Barrier();

    return(0);
  }
};

int main(int argc,char *argv[])
{
  int provided;

  // Make the MPI runtime aware that we might want to run OpenMP threads,
  // but will do MPI calls only from the master thread.
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#ifdef USE_OMP
  assert(provided == MPI_THREAD_FUNNELED);
#endif

  // As HDF5 might run MPI functions as part of its shutdown, we can not call
  // MPI_Finalize from main directly.
  atexit((void (*)())MPI_Finalize);

  int returnCode = plascom2::ParallelTest(argc,argv);
  return(returnCode);
}
