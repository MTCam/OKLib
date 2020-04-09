///
/// @file
/// @ingroup plascom2_group
/// @brief Testing object for PlasCom2.
///
/// This file contains the implementation of the
/// overarching testing object for the project.
/// Every project should have one of these as
/// it allows external entities unified methodology
/// for access project or package-specific tests.
///
#ifndef __PC2_TEST_H__
#define __PC2_TEST_H__
#ifdef _PC2_PARALLEL_
#include "COMM.H"
#endif
#include "Testing.H"
//#include "ExampleHeader.H"
#include "Profiler.H"
#include "Simulation.H"
#include "ExternalTestingFunctions.H"

#include <vector>
#include <cmath>

namespace plascom2 {

  ///
  /// Namespace for storing simple test utility fixtures
  ///
  namespace TestFixture {
    ///
    /// Simple test function fixture returns \f$2x\f$.
    /// 
    /// F1 is a 1st order function in x that has unity integral over [0,1].
    ///
    double F1(double x) { return (2.0*x); };
    ///
    /// Simple test function fixture returns \f$3x^2\f$.
    ///
    /// F2 is quadratic function in x and integrates to unity over [0,1].
    ///
    double F2(double x) { return (3.0*x*x); };
  };
  ///
  /// Project-specific test results type.
  ///
  /// Test results types can be custom implemented
  /// by the project developer.  The requirement
  /// is that they produce the proper test results
  /// output when used with outstreams.
  ///
  typedef ix::test::results TestResults;

  /// 
  /// Project-specific testing object.
  ///
  /// The project-specific testing object encapsulates
  /// all tests for the project. This is the interface
  /// to all testing available for the given project.
  /// It implements the interface given by 
  /// ix::util::TestingObject.
  ///
  /// Every project should define a "TestingObject"
  /// which implements (or interfaces) all of the tests
  /// for the project.
  /// If necessary (e.g. for parallel tests through
  /// batch systems), this is the object responsible
  /// for spawning off children processes to conduct
  /// parallel tests and for gathering the results
  /// of those tests.
  ///
  template<typename ResultsType>
  class TestingObject : public ix::test::manager<ResultsType>
  {
    typedef ix::test::manager<ResultsType> TestingObjectBaseType;
  protected:
    // Example test fixtures for PlasCom2 testing.
    //
    // In a real project, there would be many test fixtures
    // in the TestingObject.
    ///
    /// A sample string test fixture. 
    ///
    std::string ExampleTestFixture;
    ///
    /// A set of values for testing quadrature methods.
    /// 
    std::vector<int> N;
  public:
    /// 
    /// Default constructor.
    ///
    TestingObject() : TestingObjectBaseType() {};

    ///
    /// Tears down the testing fixtures if needed.
    ///
    virtual void Epilogue() {};

    ///
    /// Sets up the data fixtures for the tests.
    /// 
    /// This function gives the object a chance to 
    /// allocate and set up any "fixtures" before the
    /// tests are run.
    ///
    virtual void  Prologue(){
      ExampleTestFixture.assign("ExampleTestData");
      for(int i = 10;i < 10000000;i*=10) N.push_back(i);
    }

    ///
    /// Simple test function fixture returns \f$2x\f$.
    /// 
    /// F1 is a 1st order function in x that has unity integral over [0,1].
    ///
    double F1(double x) { return (2.0*x); };

    ///
    /// Simple test function fixture returns \f$3x^2\f$.
    ///
    /// F2 is quadratic function in x and integrates to unity over [0,1].
    ///
    double F2(double x) { return (3.0*x*x); };
    
    virtual void Test__CodeParser(ix::test::results &serialUnitResults){TestCodeParser(serialUnitResults); }; 
    virtual void Test__Roe__Matrix(ix::test::results &serialUnitResults){TestRoe_Matrix(serialUnitResults); }; 
    virtual void Test__Various__FunctionCalls(ResultsType &serialUnitResults){TestFunctionCalls(serialUnitResults); };
    virtual void Test__NumUtil__MatMat(ResultsType &serialUnitResults){TestMatMat(serialUnitResults);};
    virtual void Test__NumUtil__MatVec(ResultsType &serialUnitResults){TestMatVec(serialUnitResults);};
    virtual void Test__PCPP__Various(ResultsType &serialUnitResults){TestVarious(serialUnitResults);};
    virtual void Test__PCPP__Interval__Basic(ResultsType &serialUnitResults){TestIntervalBasic(serialUnitResults);};
    virtual void Test__State__Basic(ResultsType &serialUnitResults) {TestStateBasic(serialUnitResults);};
    virtual void Test__State__Handles(ResultsType &serialUnitResults) {TestStateHandles(serialUnitResults);};
    virtual void Test__Halo__Basic(ResultsType &serialUnitResults){TestHaloBasic(serialUnitResults);};
    virtual void Test__Halo__Threaded(ResultsType &serialUnitResults){TestHaloThreaded(serialUnitResults);};
    virtual void Test__State__Operations(ResultsType &serialUnitResults){TestStateOperations(serialUnitResults);};
//     virtual void Test__FD1(ResultsType &serialUnitResults){TestFD1(serialUnitResults);};
//     virtual void Test__FD1_2(ResultsType &serialUnitResults){TestFD1_2(serialUnitResults);};
//     virtual void Test__FD1_3(ResultsType &serialUnitResults){TestFD1_3(serialUnitResults);};
//     virtual void Test__FD1_Div3All(ResultsType &serialUnitResults){TestFD1_Div3All(serialUnitResults);};
    //    virtual void Test__1DAdvectionAdvancer__Basic(ResultsType &serialUnitResults){Test1DAdvectionAdvancer(serialUnitResults);}
    virtual void Test__Advancer__RK4(ResultsType &serialUnitResults){TestRK4Advancer(serialUnitResults);}
    virtual void Test__Advancer__RK4_2(ResultsType &serialUnitResults){TestRK4Advancer2(serialUnitResults);}
    //    virtual void Test__Advancer__RK4Euler(ResultsType &serialUnitResults){TestRK4EulerAdvancer(serialUnitResults);}
//     virtual void Test__Euler__RHS(ResultsType &serialUnitResults){TestEulerRHS(serialUnitResults);}
    virtual void Test__Euler__Kernels(ResultsType &serialUnitResults){TestEulerKernels(serialUnitResults);};
    virtual void Test__Viscid__Kernels(ResultsType &serialUnitResults){TestViscidKernels(serialUnitResults);};
    virtual void Test__Viscid__Kernels__Metrics(ResultsType &serialUnitResults){TestViscidKernelsMetrics(serialUnitResults);};
    virtual void Test__Viscid__Kernels__Curvilinear(ResultsType &serialUnitResults){TestViscidKernelsCurvilinear(serialUnitResults);};
    virtual void Test__Probe(ResultsType &serialUnitResults){TestProbe(serialUnitResults);};
    virtual void Test__Probe__BufferProbe(ResultsType &serialUnitResults){TestBufferProbe(serialUnitResults);};
    //    virtual void Test__Grid__PUBS(ResultsType &serialUnitResults){TestGrid_PUBS(serialUnitResults);};
    virtual void Test__Grid__SubRegion(ix::test::results &serialUnitResults){TestGrid_SubRegion(serialUnitResults);};
    virtual void Test__Operators__SBPInitialize(ResultsType &serialUnitResults){TestSBPInitialize(serialUnitResults);};
    virtual void Test__Operators__SBP12(ResultsType &serialUnitResults){TestOperatorSBP12(serialUnitResults);};
    virtual void Test__Operators__SBP24(ResultsType &serialUnitResults){TestOperatorSBP24(serialUnitResults);};
    virtual void Test__Operators__SBP36(ResultsType &serialUnitResults){TestOperatorSBP36(serialUnitResults);};
    virtual void Test__Operators__ApplyBlobs(ResultsType &serialUnitResults){TestApplyOperatorBlobs(serialUnitResults);};
    virtual void Test__Operators__Arithmetic(ResultsType &serialUnitResults){TestOperators(serialUnitResults);};
    virtual void Test__Operators__HoleDetection(ResultsType &serialUnitResults){TestHoleDetection(serialUnitResults);};
    //    virtual void Test__ApplyOperator1(ResultsType &serialUnitResults){TestApplyOperator1(serialUnitResults);};
    ///
    /// Runs all test associated with the ICE CS-Tool
    ///
    /// @param serialUnitResults PlasCom2::TestResults object to store test results.
    ///
#ifdef HAVE_ICE
    virtual void Test__ICE(ResultsType &serialUnitResults){TestICE(serialUnitResults);};
#ifdef HAVE_ICE_ROSE
    virtual void Test__ICE__Rose(ResultsType &serialUnitResults){TestICE_Rose(serialUnitResults);};
#endif
#endif
    virtual void Test__WENO__Stencils(ResultsType &serialUnitResults){TestWENO_Stencils(serialUnitResults);};
    virtual void Test__WENO__ReconstPointVal(ResultsType &serialUnitResults){TestWENO_ReconstPointVal(serialUnitResults);};
    virtual void Test__WENO__ReconstPointValSten(ResultsType &serialUnitResults){TestWENO_ReconstPointValSten(serialUnitResults);};
    virtual void Test__WENO__SmoothInd(ResultsType &serialUnitResults){TestWENO_SmoothInd(serialUnitResults);};
    virtual void Test__WENO__Project(ResultsType &serialUnitResults){TestWENO_Project(serialUnitResults);};
    virtual void Test__WENO__EntropyFix(ResultsType &serialUnitResults){TestWENO_EntropyFix(serialUnitResults);};
    virtual void Test__Pooya__WENO1(ResultsType &serialUnitResults){TestPooyaWENO1(serialUnitResults); };
    virtual void Test__MetricOps__PointMetric(ResultsType &serialUnitResults){TestMetricOps_PointMetric(serialUnitResults); };

    // EOS tests
    virtual void Test__EOS__ComputeBuffers(ResultsType &serialUnitResults) 
    {TestEOS_ComputeBuffers(serialUnitResults);};
    
    // Commented until testing data examination [mtc]
    virtual void Test__EOS__MultiIdealCaloricallyPerfect(ResultsType &serialUnitResults) 
    {TestEOS_MultiIdealCaloricallyPerfect(serialUnitResults);};
#ifdef HAVE_CANTERA
    virtual void Test__Cantera__Agent(ResultsType &serialUnitResults)
    { TestCantera_Agent(serialUnitResults); }
#ifdef USE_CANTERA
    virtual void Test__Cantera__CanteraMultiIdeal(ResultsType &serialUnitResults) 
    {TestCantera_CanteraMultiIdeal(serialUnitResults);};
    virtual void Test__Cantera__CanteraMultiIdealND(ResultsType &serialUnitResults) 
    {TestCantera_CanteraMultiIdealND(serialUnitResults);};
    virtual void Test__Cantera__PlasmaEOS(ResultsType &serialUnitResults)
    {TestCantera_PlasmaEOS(serialUnitResults); };
    virtual void Test__Chemistry__CanteraSource(ResultsType &serialUnitResults) 
      {TestChemistry_CanteraSource(serialUnitResults);};
    virtual void Test__Chemistry__CanteraPlasmaSource(ResultsType &serialUnitResults) 
      {TestChemistry_CanteraPlasmaSource(serialUnitResults);};
#endif

    virtual void Test__Cantera__PrometheusMultiIdeal(ResultsType &serialUnitResults) 
      {TestCantera_PrometheusMultiIdeal(serialUnitResults);};
    virtual void Test__Cantera__PrometheusMultiIdealNew(ResultsType &serialUnitResults) 
      {TestCantera_PrometheusMultiIdealNew(serialUnitResults);};
    virtual void Test__Cantera__PrometheusExtended(ResultsType &serialUnitResults) 
      {TestCantera_PrometheusExtended(serialUnitResults);};
    virtual void Test__Cantera__PrometheusMultiIdealND(ResultsType &serialUnitResults) 
      {TestCantera_PrometheusMultiIdealND(serialUnitResults);};
    virtual void Test__Cantera__PrometheusGetTemperature(ResultsType &serialUnitResults) 
      {TestCantera_PrometheusGetTemperature(serialUnitResults);};
    virtual void Test__Cantera__PrometheusGetTemperatureExtended(ResultsType &serialUnitResults) 
      {TestCantera_PrometheusGetTemperatureExtended(serialUnitResults);};
    virtual void Test__Chemistry__PrometheusSource(ResultsType &serialUnitResults) 
      {TestChemistry_PrometheusSource(serialUnitResults);};
    virtual void Test__Chemistry__PrometheusPlasmaSource(ResultsType &serialUnitResults) 
      {TestChemistry_PrometheusPlasmaSource(serialUnitResults);};
    virtual void Test__Chemistry__PrometheusSourceExtended(ResultsType &serialUnitResults) 
      {TestChemistry_PrometheusSourceExtended(serialUnitResults);};
    virtual void Test__Cantera__PrometheusPlasmaEOS(ResultsType &serialUnitResults)
    {TestCantera_PrometheusPlasmaEOS(serialUnitResults); };
#endif
#ifdef HAVE_CUDA
    virtual void Test__CUDA__Basic(ResultsType &result)
    { TestCUDA_Basic(result); };
#endif

    virtual void Test__Examples__Serial(ResultsType &serialUnitResults)
    {ExamplesSerialTestDriver(serialUnitResults); };
    
    ///
    /// Runs all tests implemented by the PlasCom2::TestingObject.
    ///
    /// @param result PlasCom2::TestResults object to store test results.
    ///
    virtual void Process(ResultsType &result){
      Prologue();
      Test__Examples__Serial(result);
      //      for(int i = 0;i < 1000;i++){
      Test__PCPP__Interval__Basic(result);
      Test__State__Basic(result);
      Test__State__Handles(result);
      Test__State__Operations(result);
      Test__Halo__Basic(result);
      Test__Halo__Threaded(result);
//       Test__FD1(result);
//       Test__FD1_2(result);
//       Test__FD1_3(result);
//       Test__FD1_Div3All(result);
      //      Test__1DAdvectionAdvancer__Basic(result);
      Test__Advancer__RK4(result);
      Test__Advancer__RK4_2(result);
      //      Test__Advancer__RK4Euler(result);
//       Test__Euler__RHS(result);
      Test__Euler__Kernels(result);
      Test__Viscid__Kernels(result);
      Test__Viscid__Kernels__Metrics(result);
      Test__Viscid__Kernels__Curvilinear(result);
      Test__Probe(result);
      Test__Probe__BufferProbe(result);
      //      Test__Grid__PUBS(result);
      Test__Grid__SubRegion(result);
      Test__Operators__SBPInitialize(result);
      Test__Operators__SBP12(result);
      Test__Operators__SBP24(result);
      Test__Operators__SBP36(result);
      Test__Operators__ApplyBlobs(result);
      Test__Operators__Arithmetic(result);
      Test__Operators__HoleDetection(result);
      //      Test__ApplyOperator1(result);
      Test__PCPP__Various(result);
      Test__NumUtil__MatMat(result);
      Test__NumUtil__MatVec(result);
#ifdef HAVE_ICE
      Test__ICE(result);
#ifdef HAVE_ICE_ROSE
      Test__ICE__Rose(result);
#endif
#endif 
      Test__WENO__Stencils(result);
      Test__WENO__ReconstPointVal(result);
      Test__WENO__ReconstPointValSten(result);
      Test__WENO__SmoothInd(result);
      Test__WENO__Project(result);
      Test__WENO__EntropyFix(result);
      Test__Pooya__WENO1(result);
      Test__EOS__ComputeBuffers(result);
      Test__Various__FunctionCalls(result);
      Test__Roe__Matrix(result);
      Test__CodeParser(result);
      Test__MetricOps__PointMetric(result);
      Test__EOS__MultiIdealCaloricallyPerfect(result);
#ifdef HAVE_CANTERA
      //Test__Cantera__InitializeCantera(result);
      //Test__Chemistry__Cantera__UniformBoxPeriodic__Regression(result);
#ifdef USE_CANTERA
      Test__Cantera__CanteraMultiIdeal(result);
      Test__Cantera__CanteraMultiIdealND(result);
      Test__Chemistry__CanteraSource(result);
#endif
      Test__Cantera__Agent(result);

      Test__Cantera__PrometheusMultiIdeal(result);
      Test__Cantera__PrometheusMultiIdealND(result);
      Test__Cantera__PrometheusGetTemperature(result);
      Test__Cantera__PrometheusGetTemperatureExtended(result);
      Test__Chemistry__PrometheusSource(result);
      Test__Chemistry__PrometheusSourceExtended(result);
#endif
#ifdef HAVE_CUDA
      Test__CUDA__Basic(result);
#endif
      //      }
      Epilogue();
    }
    
    ///
    /// Runs a test specified by name
    ///
    /// @param name String name of the test to run.
    /// @param result PlasCom2::TestResults object to store test results.
    ///
    virtual void RunTest(const std::string &name,ResultsType &result)
    {
      Prologue();
      std::cout << "************************************************" << std::endl
                << "\t Running Tests from serial suite " << name << std::endl
                << "************************************************" << std::endl;
      if(name == "PCPP"){
        Test__PCPP__Various(result);
        Test__PCPP__Interval__Basic(result);
        Test__PCPP__Various(result);
      } else if(name == "Examples") {
        Test__Examples__Serial(result);
      } else if (name == "Cantera"){
#ifdef HAVE_CANTERA 
        Test__Cantera__Agent(result);
#endif
      } else if (name == "MetricOps") {
        Test__MetricOps__PointMetric(result);
      } else if (name == "Roe"){
        Test__Roe__Matrix(result);
      } else if(name == "Various"){
        Test__Various__FunctionCalls(result);
        Test__CodeParser(result);
      } else if(name == "NumUtil") {
        Test__NumUtil__MatMat(result);
        Test__NumUtil__MatVec(result);
      } else if(name == "State"){
        Test__State__Basic(result);
        Test__State__Handles(result);
        Test__State__Operations(result);
      } else if(name == "Operators"){
        Test__Operators__SBPInitialize(result);
        Test__Operators__SBP12(result);
        Test__Operators__SBP24(result);
        Test__Operators__SBP36(result);
        Test__Operators__ApplyBlobs(result);
        Test__Operators__Arithmetic(result);
        Test__Operators__HoleDetection(result);
        Test__MetricOps__PointMetric(result);
        //	Test__ApplyOperator1(result);
      } else if(name == "Halo"){
        Test__Halo__Basic(result);
        Test__Halo__Threaded(result);
      } else if(name == "HaloBasic") {
        Test__Halo__Basic(result);
      } else if(name == "HaloThreaded"){
        Test__Halo__Threaded(result);
//       else if(name == "Advect1D")
//         Test__1DAdvectionAdvancer__Basic(result);
      } else if(name == "Advancer"){
        Test__Advancer__RK4(result);
        Test__Advancer__RK4_2(result);
        //        Test__Advancer__RK4Euler(result);
      } else if(name == "Euler"){
        //        Test__Euler__RHS(result);
        Test__Euler__Kernels(result);
        //Test__RK4EulerAdvancer(result);
      } else if(name == "Viscid"){
        Test__Viscid__Kernels(result);
        Test__Viscid__Kernels__Metrics(result);
        Test__Viscid__Kernels__Curvilinear(result);
      }
      //       else if(name == "FD1"){
      //         Test__FD1(result);
      //         Test__FD1_2(result);
      //         Test__FD1_3(result);
      //         Test__FD1_Div3All(result);
      //       }
      else if(name == "Probe"){
        Test__Probe(result);
        Test__Probe__BufferProbe(result);
      } else if(name == "Grid"){
        //        Test__Grid__PUBS(result);
        Test__Grid__SubRegion(result);
#ifdef HAVE_ICE
      } else if(name == "ICE"){
        Test__ICE(result);
#ifdef HAVE_ICE_ROSE
      } else if(name == "ICE-Rose"){
        Test__ICE__Rose(result);
#endif
#endif
      } else if(name == "WENO"){
        Test__WENO__Stencils(result);
        Test__WENO__ReconstPointVal(result);
        Test__WENO__ReconstPointValSten(result);
        Test__WENO__SmoothInd(result);
        Test__WENO__Project(result);
        Test__WENO__EntropyFix(result);
        //        Test__WENO__Pooya1(result);
      } else if(name == "Pooya"){
        Test__Pooya__WENO1(result);
      } else if (name == "EOS"){
        Test__EOS__ComputeBuffers(result);
        Test__EOS__MultiIdealCaloricallyPerfect(result);
#ifdef HAVE_CANTERA
#ifdef USE_CANTERA
        Test__Cantera__CanteraMultiIdeal(result);
        Test__Cantera__CanteraMultiIdealND(result);
        Test__Cantera__PlasmaEOS(result);
#endif
        Test__Cantera__PrometheusPlasmaEOS(result);

        Test__Cantera__PrometheusMultiIdeal(result);
        Test__Cantera__PrometheusMultiIdealNew(result);
        Test__Cantera__PrometheusMultiIdealND(result);
        Test__Cantera__PrometheusGetTemperature(result);
        Test__Cantera__PrometheusExtended(result);
        Test__Cantera__PrometheusGetTemperatureExtended(result);
#endif
      } else if (name == "Chemistry"){
#ifdef HAVE_CANTERA
#ifdef USE_CANTERA
        Test__Chemistry__CanteraSource(result);
        Test__Chemistry__CanteraPlasmaSource(result);
#endif
        Test__Chemistry__PrometheusSource(result);
        Test__Chemistry__PrometheusPlasmaSource(result);
        Test__Chemistry__PrometheusSourceExtended(result);
#endif
      } else if (name == "CUDA") {
#ifdef HAVE_CUDA
      Test__CUDA__Basic(result);
#endif
      } else {
        //No such test
        std::cerr << "Error: No such suite '" << name <<"'. Exiting." << std::endl;
        std::exit(1);
      } 
      std::cout << "************************************************" << std::endl
                << "\t Done Running Tests from serial suite " << name << std::endl
                << "************************************************" << std::endl;
      
      Epilogue();
    }
    
    ///
    /// Runs a list of tests specified by name
    /// 
    /// @param test_names list of string test names
    /// @param result PlasCom2::TestResults object to store test results.
    ///
    virtual void ProcessTests(std::list<std::string> &test_names,ResultsType &result){
      Prologue();
      std::list<std::string>::iterator tni = test_names.begin();
      while(tni != test_names.end())
        RunTest(*tni++,result);
      Epilogue();
    }


  };

#ifdef _PC2_PARALLEL_
  /// 
  /// Project-specific parallel testing object.
  ///
  /// The project-specific parallel testing object encapsulates
  /// all parallel tests for the project. This is the interface
  /// to all "in-core" parallel testing available for the given project.
  /// It implements the interface given by 
  /// PlasCom2::TestingObject and has an 
  /// ix::comm::CommunicatorObject for parallel capabilities.
  ///
  template<typename CommType,typename ResultsType>
  class ParallelTestingObject : TestingObject<ResultsType>
  {
  protected:
    // Example test fixtures for parallel PlasCom2 testing.
    //
    /// 
    /// The CommunicatorObject for the parallel testing object.
    ///
    CommType _communicator;
    ///
    /// A set of values for testing quadrature methods.
    /// 
    std::vector<int> N;
  public:
    ///
    /// Constructor with incoming communicator.
    ///
    ParallelTestingObject(CommType &incomm) :
      TestingObject<ResultsType>(), _communicator(incomm) {};
    CommType &GetCommunicator() { return(_communicator); };

    ///
    /// Tears down the testing fixtures if needed.
    ///
    virtual void Epilogue() {};

    ///
    /// Sets up the data fixtures for the tests.
    /// 
    /// This function gives the object a chance to 
    /// allocate and set up any "fixtures" before the
    /// tests are run.
    ///
    virtual void  Prologue(){
      for(int i = _communicator.Size();i <= 1000000;i*=2) N.push_back(i);
    };



    ///
    /// Runs all tests implemented by the PlasCom2::TestingObject.
    ///
    /// @param result PlasCom2::TestResults object to store test results.
    ///
    virtual void Process(ResultsType &result){
      Prologue();
      Test__Examples__Parallel(result);
      //      Test__Interval__Basic(result);
      Test__Halo(result);
      Test__Halo__Periodic(result);
      
#ifdef HAVE_CVODE
      Test__CVODE__SimpleInit(result);
#endif

#ifdef ENABLE_HDF5
      Test__HDF5__LegacyFileInfo(result);
      Test__HDF5__Legacy2DBSGRead(result);
      Test__HDF5__PCPPRead(result);
      Test__HDF5__PCPPWrite(result);
      Test__HDF5__PCPPHyperSlab(result);
      Test__HDF5__GridIO(result);
      Test__HDF5__PBSGridIO(result);
      Test__HDF5__Output(result);
      Test__HDF5__Input(result);
#endif
      Test__Grid__PBSIH(result);
      Test__Grid__CartesianMetric(result);
      Test__Grid__RectilinearMetric(result);
      Test__Grid__CurvilinearMetric(result);
      Test__Grid__CurvilinearVGWavy(result);
      Test__Grid__UniformGridMetrics(result);
      Test__Grid__UniformPeriodicMetrics(result);

#ifdef USE_OVERKIT
      Test__Overkit__CreateDomain(result);
      Test__Overkit__Interface2D(result);
      Test__Overkit__Interface3D(result);
      Test__Overkit__AcousticPulse2D__BoxInBoxInjection(result);
#endif

      Test__Euler__RHS(result);
      Test__Euler__RHS2(result);
//       Test__WENO__ApplyWENO(result);
//       Test__WENO__RHS(result);
//       Test__WENO__StencilConnectivity(result);
//       Test__WENO__ApplyWENOBounded(result);
      Test__PlasCom2__AcousticPulse2D__Periodic(result);
      Test__PlasCom2__AcousticPulse2D__Farfield(result);
      Test__PlasCom2__AcousticPulse2D__Slipwall(result);
      Test__PlasCom2__AcousticPulse2D__Sponge(result);
      Test__PlasCom2__AcousticPulse3D__Periodic(result);
      Test__PlasCom2__AcousticPulse3D__Farfield(result);
      Test__PlasCom2__AcousticPulse3D__Slipwall(result);
      Test__PlasCom2__AcousticPulse3D__Sponge(result);
      Test__PlasCom2__Rectilinear2D__Periodic(result);
      Test__PlasCom2__Curvilinear2D__Periodic(result);
      Test__PlasCom2__Rectilinear2D__Farfield(result);
      Test__PlasCom2__Curvilinear2D__Farfield(result);
      Test__PlasCom2__Rectilinear2D__Slipwall(result);
      Test__PlasCom2__Curvilinear2D__Slipwall(result);
      Test__PlasCom2__Rectilinear2D__Sponge(result);
      Test__PlasCom2__Curvilinear2D__Sponge(result);
      Test__PlasCom2__Curvilinear3D__Sponge(result);
      Test__PlasCom2__Curvilinear3D__Slipwall(result);
      Test__PlasCom2__Curvilinear3D__Periodic(result);
      Test__PlasCom2__Curvilinear3D__Farfield(result);
      Test__PlasCom2__Rectilinear3D__Sponge(result);
      Test__PlasCom2__Rectilinear3D__Slipwall(result);
      Test__PlasCom2__Rectilinear3D__Periodic(result);
      Test__PlasCom2__Rectilinear3D__Farfield(result);
      Test__PlasCom2__Restart(result);

      Test__WENO__Convergence(result);
      Test__Integrated__ViscousShock2DX(result);
      Test__Integrated__VSRectilinear2DX(result);
      Test__Integrated__VSCurvilinear2DX(result);
      Test__Integrated__ViscousShock2DY(result);
      Test__Integrated__VSRectilinear2DY(result);
      Test__Integrated__VSCurvilinear2DY(result);
      Test__Integrated__ViscousShock3DZ(result);
      Test__Integrated__VSRectilinear3DZ(result);
      Test__Integrated__VSCurvilinear3DZ(result);
      Test__PlasCom2__ViscousShock2DX__Regression(result);
      Test__PlasCom2__ViscousShock2DY__Regression(result);
      Test__PlasCom2__ViscousShock3DZ__Regression(result);
      Test__PlasCom2__ViscousShockAD3DZ__Regression(result);
      Test__Integrated__Poiseuille2DX(result);
      Test__Integrated__PFRectilinear2DX(result);
      Test__Integrated__PFCurvilinear2DX(result);
      Test__Integrated__Poiseuille2DY(result);
      Test__Integrated__PFRectilinear2DY(result);
      Test__Integrated__PFCurvilinear2DY(result);
      Test__Integrated__Poiseuille3DZ(result);
      Test__Integrated__PFRectilinear3DZ(result);
      Test__Integrated__PFCurvilinear3DZ(result);
      Test__PlasCom2__Poiseuille2DX__Regression(result);
      Test__PlasCom2__Poiseuille2DY__Regression(result);
      Test__PlasCom2__Poiseuille3DZ__Regression(result);
      
      Test__PlasCom2__VSCurvilinear2DX__Regression(result);
      Test__PlasCom2__VSRectilinear2DX__Regression(result);
      Test__PlasCom2__VSCurvilinear2DY__Regression(result);
      Test__PlasCom2__VSRectilinear2DY__Regression(result);
      Test__PlasCom2__VSCurvilinear3DZ__Regression(result);
      Test__PlasCom2__VSRectilinear3DZ__Regression(result);

      Test__PlasCom2__PFCurvilinear2DX__Regression(result);
      Test__PlasCom2__PFRectilinear2DX__Regression(result);
      Test__PlasCom2__PFCurvilinear2DY__Regression(result);
      Test__PlasCom2__PFRectilinear2DY__Regression(result);
      Test__PlasCom2__PFCurvilinear3DZ__Regression(result);
      Test__PlasCom2__PFRectilinear3DZ__Regression(result);
      
      Test__PlasCom2__Scalar__Advection2DX(result);
      Test__PlasCom2__Scalar__Advection2DY(result);
      Test__PlasCom2__Scalar__AdvectionDiffusion2DX(result);
      Test__PlasCom2__Scalar__AdvectionDiffusion2DY(result);
      
      Test__PlasCom2__PooyaWENO(result);
      Test__PlasCom2__PavelAxisymmetryDriver(result);
#ifndef DISABLE_Y6      
      Test__Y6(result);
#endif
// #ifdef HAVE_CANTERA
//       Test__PlasCom2__PavelWENO(result);
//       Test__PlasCom2__PavelSBP(result);
// #endif
#ifdef USE_STELLA
#ifdef USE_PETSC
      Test__Efield__Single2D(result);
      Test__Efield__Single3D(result);
#ifdef USE_OVERKIT
      Test__EfieldOverset__Multiple2D(result);
      Test__EfieldOverset__Multiple3D(result);
      Test__EfieldOverset__FourGrid(result);
#endif
#endif
#endif

      Test__Comm__Split(result);
      Test__PCPP__SetupCartNeighborhood(result);
      Test__Various__PartitionIntervals(result);
      Test__PAPIX(result);
      Epilogue();
    }
    
    ///
    /// Runs a test specified by name
    ///
    /// @param name String name of the test to run.
    /// @param result PlasCom2::TestResults object to store test results.
    ///
    virtual void RunTest(const std::string &name,ResultsType &result)
    {

      if(_communicator.Rank() == 0) {
        std::cout << "************************************************" << std::endl
                  << "\t Running Tests from parallel suite " << name << std::endl
                  << "************************************************" << std::endl;
      }
      
      if(name == "Halo"){
        Test__Halo(result);
        Test__Halo__Periodic(result);
      } else if(name == "Examples") {
        Test__Examples__Parallel(result);
      } else if (name == "COMM") {
	Test__Comm__Split(result);
        Test__PCPP__SetupCartNeighborhood(result);
      } else if(name == "HDF5") {
#ifdef ENABLE_HDF5
        Test__HDF5__LegacyFileInfo(result);
        Test__HDF5__Legacy2DBSGRead(result);
        Test__HDF5__PCPPRead(result);
        Test__HDF5__PCPPWrite(result);
        Test__HDF5__PCPPHyperSlab(result);
        Test__HDF5__GridIO(result);
        Test__HDF5__PBSGridIO(result);
        Test__HDF5__Output(result);
        Test__HDF5__Input(result);
#else
        std::cout << "PlasCom2::ParallelTestingObject: HDF5 not enabled. Skipping HDF5 tests." << std::endl;
#endif
      } else if (name == "Grid") {
        Test__Grid__PBSIH(result);
        Test__Grid__CartesianMetric(result);
        Test__Grid__RectilinearMetric(result);
        Test__Grid__CurvilinearMetric(result);
        Test__Grid__CurvilinearVGWavy(result);
        Test__Grid__UniformGridMetrics(result);
        Test__Grid__UniformPeriodicMetrics(result);
#ifdef USE_OVERKIT
      } else if (name == "Overkit"){
        Test__Overkit__CreateDomain(result);
        Test__Overkit__Interface2D(result);
        Test__Overkit__Interface3D(result);
        Test__Overkit__AcousticPulse2D__BoxInBoxInjection(result);
#endif
      } else if (name == "Euler"){
        Test__Euler__RHS(result);
        Test__Euler__RHS2(result);
      } else if (name == "Viscid"){
        Test__Velocity__Gradient(result);
        Test__Velocity__Gradient__Periodic(result);
        Test__Viscid__RHS(result);
      } else if (name == "WENO"){
        //        Test__WENO__ApplyWENO(result);
        //        Test__WENO__RHS(result);
        //         Test__WENO__StencilConnectivity(result);
        //         Test__WENO__ApplyWENOBounded(result);
      } else if (name == "Chemistry"){
#ifdef HAVE_CANTERA
        Test__Chemistry__InitializeCantera(result);
        Test__Chemistry__Cantera__UniformBoxPeriodic__Regression(result);
        Test__Chemistry__Prometheus__UniformBoxPeriodic__Regression(result);
        Test__Chemistry__PrometheusPlasma__UniformBoxPeriodic(result);
//       Test__PlasCom2__PavelSBP(result);
//       Test__PlasCom2__PavelWENO(result);
#endif
      } else if (name == "MultiSpecies"){
        // These tests roughly duplicate the PlasCom2 tests, but exercise
        // the multi-species implementations. In most cases, results
        // are compared against the single-species equivalent
        //
        // Here we draw a distiction between passive scalars which have
        // no effect on the flow and mass fractions which directly affect the 
        // mixture material properties
        // 
        Test__MultiSpecies__Advection2DX(result);
        Test__MultiSpecies__Advection2DY(result);
        Test__MultiSpecies__AcousticPulse2D__Periodic(result);
        Test__MultiSpecies__AcousticPulse2D__SlipAdiabatic(result);
        Test__MultiSpecies__AcousticPulse2D__SlipAdiabaticCurvilinear(result);
        Test__MultiSpecies__AcousticPulse2D__Farfield(result);
        Test__MultiSpecies__AcousticPulse2D__NoSlipIso(result);
        Test__MultiSpecies__AcousticPulse3D__SlipAdiabatic_XY(result);
        Test__MultiSpecies__AcousticPulse3D__SlipAdiabaticCurvilinear_XY(result);
        Test__MultiSpecies__AcousticPulse3D__SlipAdiabatic_YZ(result);
        Test__MultiSpecies__AcousticPulse3D__SlipAdiabaticCurvilinear_YZ(result);

        //Test__MultiSpecies__VGWavy(result);
        //Test__MultiSpecies__VGWavy__Overset(result);
      } else if (name == "MultiSpeciesWENO"){
        // These tests roughly duplicate the MultiSpecies tests, but exercise
        // the WENO and NSCBC implementations. In most cases, results
        // are compared against the single-species equivalent
        //
        // Here we draw a distiction between passive scalars which have
        // no effect on the flow and mass fractions which directly affect the 
        // mixture material properties
        // 
        Test__MultiSpeciesWENO__AcousticPulse2D__Periodic(result);
        Test__MultiSpeciesWENO__AcousticPulse2D__Outflow(result);
        Test__MultiSpeciesWENO__AcousticPulse2D__SlipAdiabatic(result);
        Test__MultiSpeciesWENO__AcousticPulse2D__NoSlipIso(result);
#ifdef HAVE_CANTERA
        Test__MultiSpeciesWENO__Prometheus__AcousticPulse2D__NoSlipIso(result);
        Test__MultiSpeciesWENO__Prometheus__Shock1D(result);
#endif
#ifdef USE_STELLA
#ifdef USE_PETSC
      } else if (name == "Efield"){
        Test__Efield__Single2D(result);
        Test__Efield__Single3D(result);
#ifdef USE_OVERKIT
      } else if (name == "EfieldOverset"){
        Test__EfieldOverset__Multiple2D(result);
        Test__EfieldOverset__Multiple3D(result);
        Test__EfieldOverset__FourGrid(result);
#endif
#endif
#endif
      } else if (name == "PlasCom2"){
        Test__PlasCom2__AcousticPulse2D__Periodic(result);
        Test__PlasCom2__AcousticPulse2D__Farfield(result);
        Test__PlasCom2__AcousticPulse2D__Slipwall(result);
        Test__PlasCom2__AcousticPulse2D__Sponge(result);
        Test__PlasCom2__AcousticPulse3D__Periodic(result);
        Test__PlasCom2__AcousticPulse3D__Farfield(result);
        Test__PlasCom2__AcousticPulse3D__Slipwall(result);
        Test__PlasCom2__AcousticPulse3D__Sponge(result);
        Test__PlasCom2__Rectilinear2D__Periodic(result);
        Test__PlasCom2__Curvilinear2D__Periodic(result);
        Test__PlasCom2__Rectilinear2D__Farfield(result);
        Test__PlasCom2__Curvilinear2D__Farfield(result);
        Test__PlasCom2__Rectilinear2D__Slipwall(result);
        Test__PlasCom2__Curvilinear2D__Slipwall(result);
        Test__PlasCom2__Rectilinear2D__Sponge(result);
        Test__PlasCom2__Curvilinear2D__Sponge(result);

        Test__PlasCom2__Rectilinear3D__Slipwall(result);
        Test__PlasCom2__Rectilinear3D__Periodic(result);
        Test__PlasCom2__Rectilinear3D__Farfield(result);
        Test__PlasCom2__Rectilinear3D__Sponge(result);

        Test__PlasCom2__Curvilinear3D__Slipwall(result);
        Test__PlasCom2__Curvilinear3D__Periodic(result);
        Test__PlasCom2__Curvilinear3D__Farfield(result);
        Test__PlasCom2__Curvilinear3D__Sponge(result);
        Test__PlasCom2__Restart(result);

        Test__WENO__Convergence(result);
        Test__PlasCom2__ViscousShock2DX__Regression(result);
        Test__PlasCom2__ViscousShock2DY__Regression(result);
        Test__PlasCom2__ViscousShock3DZ__Regression(result);
        Test__PlasCom2__ViscousShockAD3DZ__Regression(result);
        Test__PlasCom2__Poiseuille2DX__Regression(result);
        Test__PlasCom2__Poiseuille2DY__Regression(result);
        Test__PlasCom2__Poiseuille3DZ__Regression(result);
        Test__PlasCom2__Scalar__Advection2DX(result);
        Test__PlasCom2__Scalar__Advection2DY(result);
        Test__PlasCom2__Scalar__AdvectionDiffusion2DX(result);
        Test__PlasCom2__Scalar__AdvectionDiffusion2DY(result);
        
        Test__PlasCom2__VSCurvilinear2DX__Regression(result);
        Test__PlasCom2__VSRectilinear2DX__Regression(result);
        Test__PlasCom2__VSCurvilinear2DY__Regression(result);
        Test__PlasCom2__VSRectilinear2DY__Regression(result);
        Test__PlasCom2__VSCurvilinear3DZ__Regression(result);
        Test__PlasCom2__VSRectilinear3DZ__Regression(result);

        Test__PlasCom2__PFCurvilinear2DX__Regression(result);
        Test__PlasCom2__PFRectilinear2DX__Regression(result);
        Test__PlasCom2__PFCurvilinear2DY__Regression(result);
        Test__PlasCom2__PFRectilinear2DY__Regression(result);
        Test__PlasCom2__PFCurvilinear3DZ__Regression(result);
        Test__PlasCom2__PFRectilinear3DZ__Regression(result);

      } else if (name == "PlasCom2-2D"){

        Test__PlasCom2__AcousticPulse2D__Periodic(result);
        Test__PlasCom2__AcousticPulse2D__Farfield(result);
        Test__PlasCom2__AcousticPulse2D__Slipwall(result);
        Test__PlasCom2__AcousticPulse2D__Sponge(result);
        Test__PlasCom2__Rectilinear2D__Periodic(result);
        Test__PlasCom2__Curvilinear2D__Periodic(result);
        Test__PlasCom2__Rectilinear2D__Farfield(result);
        Test__PlasCom2__Curvilinear2D__Farfield(result);
        Test__PlasCom2__Rectilinear2D__Slipwall(result);
        Test__PlasCom2__Curvilinear2D__Slipwall(result);
        Test__PlasCom2__Rectilinear2D__Sponge(result);
        Test__PlasCom2__Curvilinear2D__Sponge(result);

        Test__WENO__Convergence(result);
        Test__PlasCom2__ViscousShock2DX__Regression(result);
        Test__PlasCom2__ViscousShock2DY__Regression(result);
        Test__PlasCom2__Poiseuille2DX__Regression(result);
        Test__PlasCom2__Poiseuille2DY__Regression(result);
        Test__PlasCom2__Scalar__Advection2DX(result);
        Test__PlasCom2__Scalar__Advection2DY(result);
        Test__PlasCom2__Scalar__AdvectionDiffusion2DX(result);
        Test__PlasCom2__Scalar__AdvectionDiffusion2DY(result);
        
        Test__PlasCom2__VSCurvilinear2DX__Regression(result);
        Test__PlasCom2__VSRectilinear2DX__Regression(result);
        Test__PlasCom2__VSCurvilinear2DY__Regression(result);
        Test__PlasCom2__VSRectilinear2DY__Regression(result);

        Test__PlasCom2__PFCurvilinear2DX__Regression(result);
        Test__PlasCom2__PFRectilinear2DX__Regression(result);
        Test__PlasCom2__PFCurvilinear2DY__Regression(result);
        Test__PlasCom2__PFRectilinear2DY__Regression(result);

      } else if(name == "Axisymmetric"){

        Test__PlasCom2__PavelAxisymmetryDriver(result);
//         Test__PlasCom2__PavelSBP(result);
//         Test__PlasCom2__PavelWENO(result);

      } else if (name == "Pooya"){
        Test__PlasCom2__PooyaWENO(result);

#ifdef HAVE_CANTERA
        //        Test__PlasCom2__PooyaWENOOx(result);
#endif
        Test__PlasCom2__PooyaWENO1(result);
        Test__PlasCom2__PooyaWENO1x(result);
        Test__PlasCom2__PooyaWENO1y(result);
        Test__PlasCom2__PooyaWENO1z(result);
        Test__PlasCom2__PooyaWENO2(result);
        Test__PlasCom2__PooyaWENO3(result);
        Test__PlasCom2__PooyaWENO4(result);
        Test__PlasCom2__PooyaWENO5(result);
        Test__PlasCom2__PooyaWENOAxisymmetric(result);
        Test__PlasCom2__Pooyanoslipisop1(result);
        Test__PlasCom2__Pooyanoslipisop2(result);
        Test__PlasCom2__Pooyanoslipisom1(result);
        Test__PlasCom2__Pooyanoslipisom2(result);
        Test__PlasCom2__Pooyanoslipisop3__Regression(result);
        Test__PlasCom2__Pooyanoslipisom3__Regression(result);

        Test__PlasCom2__PooyaEnthalpyp1(result);
        Test__PlasCom2__PooyaEnthalpyp2(result);

        Test__PlasCom2__PooyaWavy3(result);

      } else if(name == "CVODE") {
#ifdef HAVE_CVODE
        Test__CVODE__SimpleInit(result);
#endif
      } else if (name == "Y6") {
#ifndef DISABLE_Y6
        Test__Y6(result);
#endif
      } else if (name == "PAPIX"){
        
        Test__PAPIX(result);

      } else if (name == "Various") {

	Test__Various__PartitionIntervals(result);

      } else if (name == "Integrated"){

        Test__Integrated__Poiseuille2DX(result);
        Test__Integrated__PFRectilinear2DX(result);
        Test__Integrated__PFCurvilinear2DX(result);
        Test__Integrated__Poiseuille2DY(result);
        Test__Integrated__PFRectilinear2DY(result);
        Test__Integrated__PFCurvilinear2DY(result);
        Test__Integrated__Poiseuille3DZ(result);
        Test__Integrated__PFRectilinear3DZ(result);
        Test__Integrated__PFCurvilinear3DZ(result);

        Test__Integrated__ViscousShock2DX(result);
        Test__Integrated__VSRectilinear2DX(result);
        Test__Integrated__VSCurvilinear2DX(result);
        Test__Integrated__ViscousShock2DY(result);
        Test__Integrated__VSRectilinear2DY(result);
        Test__Integrated__VSCurvilinear2DY(result);
        Test__Integrated__ViscousShock3DZ(result);
        Test__Integrated__VSRectilinear3DZ(result);
        Test__Integrated__VSCurvilinear3DZ(result);

        Test__Integrated__Scalar__Advection3DZ(result);
#ifdef HAVE_CANTERA
        Test__Integrated__ChemistryCantera__UniformBoxPeriodic(result);
        //Test__Integrated__ChemistryPrometheus__UniformBoxPeriodic(result);
#endif

/*
        Test__PlasCom2__Pooyanoslipisop3(result);
        Test__PlasCom2__Pooyanoslipisom3(result);
*/
      } else if (name == "Cantera") {
        // No parallel tests for just Cantera (yet)
      } else {
        //No such test
        if(_communicator.Rank() == 0) {
          std::cerr << "Error: No such suite '" << name <<"'. Exiting." << std::endl;
        }

        std::exit(1);

      }

    }

#ifndef DISABLE_Y6
    virtual void Test__Y6(ResultsType &result)
    { Y6TestDriver(result,_communicator); };
#endif

#ifdef HAVE_CVODE
    virtual void Test__CVODE__SimpleInit(ResultsType &result)
    { TestCVODE_SimpleInit(result,_communicator); };
#endif
   virtual void Test__PAPIX(ResultsType &result)
    { PAPIXTestDriver(result,_communicator); };
    virtual void Test__Halo(ResultsType &result) 
    { TestHaloParallel(result,_communicator); };
    virtual void Test__Halo__Periodic(ResultsType &result)
    {TestHaloPeriodic(result,_communicator);};
#ifdef ENABLE_HDF5
    virtual void Test__HDF5__LegacyFileInfo(ResultsType &result)
    {TestHDF5LegacyFileInfo(result,_communicator);};
    virtual void Test__HDF5__Legacy2DBSGRead(ResultsType &result)
    {TestHDF5Read2DBlockStructuredGrid(result,_communicator);};
    virtual void Test__HDF5__PCPPRead(ResultsType &result)
    {TestPCPPHDF5Read(result,_communicator);};
    virtual void Test__HDF5__PCPPWrite(ResultsType &result)
    {TestPCPPHDF5Write(result,_communicator);};
    virtual void Test__HDF5__PCPPHyperSlab(ResultsType &result)
    {TestHDF5HyperSlab(result,_communicator);};
    virtual void Test__HDF5__GridIO(ResultsType &result)
    {TestHDF5GridIO(result,_communicator);};
    virtual void Test__HDF5__PBSGridIO(ResultsType &result)
    {TestHDF5PBSGridIO(result,_communicator);};
    virtual void Test__HDF5__Output(ResultsType &result)
    { TestPlasCom2Output(result,_communicator); };
    virtual void Test__HDF5__Input(ResultsType &result)
    { TestPlasCom2Input(result,_communicator); };
#endif
    virtual void Test__Grid__PBSIH(ResultsType &result)
    {TestGrid_PBS_IntegratedHalo(result,_communicator); };
    virtual void Test__Grid__CartesianMetric(ResultsType &result)
    {TestGrid_CartesianMetric(result,_communicator); };
    virtual void Test__Grid__RectilinearMetric(ResultsType &result)
    {TestGrid_RectilinearMetric(result,_communicator); };
    virtual void Test__Grid__CurvilinearMetric(ResultsType &result)
    {TestGrid_CurvilinearMetric(result,_communicator); };
    virtual void Test__Grid__CurvilinearVGWavy(ResultsType &result)
    {TestGrid_CurvilinearVGWavy(result,_communicator); };
    virtual void Test__Grid__UniformGridMetrics(ResultsType &result)
    {TestGrid_UniformGridMetrics(result,_communicator); };
    virtual void Test__Grid__UniformPeriodicMetrics(ResultsType &result)
    {TestGrid_UniformPeriodicMetrics(result,_communicator); };
#ifdef USE_OVERKIT
    virtual void Test__Overkit__CreateDomain(ResultsType &result)
    {TestOverkit_CreateDomain(result,_communicator);};
    virtual void Test__Overkit__Interface2D(ResultsType &result)
    {TestOverkit_Interface2D(result,_communicator);};
    virtual void Test__Overkit__Interface3D(ResultsType &result)
    {TestOverkit_Interface3D(result,_communicator);};
#endif
    virtual void Test__Euler__RHS(ResultsType &result)
    {TestEulerRHS(result,_communicator); };
    virtual void Test__Euler__RHS2(ResultsType &result)
    {TestEulerRHS2(result,_communicator); };
    virtual void Test__Velocity__Gradient(ResultsType &result)
    {TestVelocityGradient(result,_communicator); };
    virtual void Test__Velocity__Gradient__Periodic(ResultsType &result)
    {TestVelocityGradientPeriodic(result,_communicator); };
    virtual void Test__Viscid__RHS(ResultsType &result)
    {TestViscidRHS(result,_communicator); };

//     virtual void Test__WENO__ApplyWENO(ResultsType &result)
//     {TestWENO_ApplyWENO(result,_communicator);};
//     virtual void Test__WENO__RHS(ResultsType &result)
//     {TestWENO_RHS(result,_communicator);};
//     virtual void Test__WENO__StencilConnectivity(ResultsType &result)
//     {TestWENO_StencilConnectivity(result,_communicator);};
//     virtual void Test__WENO__ApplyWENOBounded(ResultsType &result)
//     {TestWENO_ApplyWENOBounded(result,_communicator);};

    virtual void Test__Integrated__ViscousShock2DX(ResultsType &result)
    {TestIntegrated_ViscousShock2DX(result,_communicator);};
    virtual void Test__Integrated__VSRectilinear2DX(ResultsType &result)
    {TestIntegrated_VSRectilinear2DX(result,_communicator);};
    virtual void Test__Integrated__VSCurvilinear2DX(ResultsType &result)
    {TestIntegrated_VSCurvilinear2DX(result,_communicator);};
    virtual void Test__Integrated__ViscousShock2DY(ResultsType &result)
    {TestIntegrated_ViscousShock2DY(result,_communicator);};
    virtual void Test__Integrated__VSRectilinear2DY(ResultsType &result)
    {TestIntegrated_VSRectilinear2DY(result,_communicator);};
    virtual void Test__Integrated__VSCurvilinear2DY(ResultsType &result)
    {TestIntegrated_VSCurvilinear2DY(result,_communicator);};
    virtual void Test__Integrated__ViscousShock3DZ(ResultsType &result)
    {TestIntegrated_ViscousShock3DZ(result,_communicator);};
    virtual void Test__Integrated__VSRectilinear3DZ(ResultsType &result)
    {TestIntegrated_VSRectilinear3DZ(result,_communicator);};
    virtual void Test__Integrated__VSCurvilinear3DZ(ResultsType &result)
    {TestIntegrated_VSCurvilinear3DZ(result,_communicator);};

    virtual void Test__Integrated__Poiseuille2DX(ResultsType &result)
    {TestIntegrated_Poiseuille2DX(result,_communicator);};
    virtual void Test__Integrated__PFRectilinear2DX(ResultsType &result)
    {TestIntegrated_PFRectilinear2DX(result,_communicator);};
    virtual void Test__Integrated__PFCurvilinear2DX(ResultsType &result)
    {TestIntegrated_PFCurvilinear2DX(result,_communicator);};
    virtual void Test__Integrated__Poiseuille2DY(ResultsType &result)
    {TestIntegrated_Poiseuille2DY(result,_communicator);};
    virtual void Test__Integrated__PFRectilinear2DY(ResultsType &result)
    {TestIntegrated_PFRectilinear2DY(result,_communicator);};
    virtual void Test__Integrated__PFCurvilinear2DY(ResultsType &result)
    {TestIntegrated_PFCurvilinear2DY(result,_communicator);};
    virtual void Test__Integrated__Poiseuille3DZ(ResultsType &result)
    {TestIntegrated_Poiseuille3DZ(result,_communicator);};
    virtual void Test__Integrated__PFRectilinear3DZ(ResultsType &result)
    {TestIntegrated_PFRectilinear3DZ(result,_communicator);};
    virtual void Test__Integrated__PFCurvilinear3DZ(ResultsType &result)
    {TestIntegrated_PFCurvilinear3DZ(result,_communicator);};
#ifdef HAVE_CANTERA
    //virtual void Test__Integrated__MultiSpecies__Advection3DZ(ResultsType &result)
    //{TestIntegrated_MultiSpecies_Advection3DZ(result,_communicator);};
    virtual void Test__Integrated__ChemistryCantera__UniformBoxPeriodic(ResultsType &result)
    {TestIntegrated_ChemistryCantera_UniformBoxPeriodic(result,_communicator);};

    virtual void Test__Chemistry__InitializeCantera(ResultsType &result)
    {TestChemistry_InitializeCantera(result,_communicator);};

    virtual void Test__Chemistry__Cantera__UniformBoxPeriodic__Regression(ResultsType &result)
    {TestChemistry_Cantera_UniformBoxPeriodic_Regression(result,_communicator);};

    virtual void Test__Chemistry__Prometheus__UniformBoxPeriodic__Regression(ResultsType &result)
    {TestChemistry_Prometheus_UniformBoxPeriodic_Regression(result,_communicator);};
    virtual void Test__Chemistry__PrometheusPlasma__UniformBoxPeriodic(ResultsType &result)
    {TestChemistry_PrometheusPlasma_UniformBoxPeriodic(result,_communicator);};

    virtual void Test__Integrated__ChemistryPrometheus__UniformBoxPeriodic(ResultsType &result)
    {TestIntegrated_ChemistryPrometheus_UniformBoxPeriodic(result,_communicator);};
#endif



    virtual void Test__PlasCom2__Scalar__Advection2DX(ResultsType &result)
    {TestPlasCom2_Scalar_Advection2DX(result,_communicator);};
    virtual void Test__PlasCom2__Scalar__Advection2DY(ResultsType &result)
    {TestPlasCom2_Scalar_Advection2DY(result,_communicator);};
    virtual void Test__PlasCom2__Scalar__AdvectionDiffusion2DX(ResultsType &result)
    {TestPlasCom2_Scalar_AdvectionDiffusion2DX(result,_communicator);};
    virtual void Test__PlasCom2__Scalar__AdvectionDiffusion2DY(ResultsType &result)
    {TestPlasCom2_Scalar_AdvectionDiffusion2DY(result,_communicator);};
    virtual void Test__Integrated__Scalar__Advection3DZ(ResultsType &result)
    {TestIntegrated_Scalar_Advection3DZ(result,_communicator);};

    virtual void Test__MultiSpecies__Advection2DX(ResultsType &result)
    {TestMultiSpecies_Advection2DX(result,_communicator);};
    virtual void Test__MultiSpecies__Advection2DY(ResultsType &result)
    {TestMultiSpecies_Advection2DY(result,_communicator);};
    virtual void Test__MultiSpecies__AcousticPulse2D__Periodic(ResultsType &result)
    {TestMultiSpecies_AcousticPulsePeriodic2D(result,_communicator);};
    virtual void Test__MultiSpecies__AcousticPulse2D__Farfield(ResultsType &result)
    {TestMultiSpecies_AcousticPulseFarfield2D(result,_communicator);};
    virtual void Test__MultiSpecies__AcousticPulse2D__SlipAdiabatic(ResultsType &result)
    {TestMultiSpecies_AcousticPulseSlipAdiabatic2D(result,_communicator);};
    virtual void Test__MultiSpecies__AcousticPulse2D__SlipAdiabaticCurvilinear(ResultsType &result)
    {TestMultiSpecies_AcousticPulseSlipAdiabatic2DCurvilinear(result,_communicator);};
    virtual void Test__MultiSpecies__AcousticPulse2D__NoSlipIso(ResultsType &result)
    {TestMultiSpecies_AcousticPulseNoSlipIso2D(result,_communicator);};
   
    virtual void Test__MultiSpecies__AcousticPulse3D__SlipAdiabatic_XY(ResultsType &result)
    {TestMultiSpecies_AcousticPulseSlipAdiabatic3D_XY(result,_communicator);};
    virtual void Test__MultiSpecies__AcousticPulse3D__SlipAdiabaticCurvilinear_XY(ResultsType &result)
    {TestMultiSpecies_AcousticPulseSlipAdiabatic3DCurvilinear_XY(result,_communicator);};
    virtual void Test__MultiSpecies__AcousticPulse3D__SlipAdiabatic_YZ(ResultsType &result)
    {TestMultiSpecies_AcousticPulseSlipAdiabatic3D_YZ(result,_communicator);};
    virtual void Test__MultiSpecies__AcousticPulse3D__SlipAdiabaticCurvilinear_YZ(ResultsType &result)
    {TestMultiSpecies_AcousticPulseSlipAdiabatic3DCurvilinear_YZ(result,_communicator);};

    virtual void Test__MultiSpeciesWENO__AcousticPulse2D__Periodic(ResultsType &result)
    {TestMultiSpeciesWENO_AcousticPulsePeriodic2D(result,_communicator);};
    virtual void Test__MultiSpeciesWENO__AcousticPulse2D__Outflow(ResultsType &result)
    {TestMultiSpeciesWENO_AcousticPulseOutflow2D(result,_communicator);};
    virtual void Test__MultiSpeciesWENO__AcousticPulse2D__SlipAdiabatic(ResultsType &result)
    {TestMultiSpeciesWENO_AcousticPulseSlipAdiabatic2D(result,_communicator);};
    virtual void Test__MultiSpeciesWENO__AcousticPulse2D__NoSlipIso(ResultsType &result)
    {TestMultiSpeciesWENO_AcousticPulseNoSlipIso2D(result,_communicator);};
    virtual void Test__MultiSpecies__VGWavy(ResultsType &result)
    {TestMultiSpecies_VGWavy(result,_communicator);};  


#ifdef HAVE_CANTERA
    virtual void Test__MultiSpeciesWENO__Prometheus__AcousticPulse2D__NoSlipIso(ResultsType &result)
    {TestMultiSpeciesWENO_Prometheus_AcousticPulseNoSlipIso2D(result,_communicator);};
    virtual void Test__MultiSpeciesWENO__Prometheus__Shock1D(ResultsType &result)
    {TestMultiSpeciesWENO_Prometheus_Shock1D(result,_communicator);};
#endif
    //virtual void Test__MultiSpecies__VGWavy__Overset(ResultsType &result)
    //{TestMultiSpecies_VGWavy_Overset(result,_communicator);};
 
/*
    virtual void Test__MultiSpeciesWENO__AcousticPulse3D__SlipAdiabatic_XY(ResultsType &result)
    {TestMultiSpeciesWENO_AcousticPulseSlipAdiabatic3D_XY(result,_communicator);};
    virtual void Test__MultiSpeciesWENO__AcousticPulse3D__SlipAdiabatic_YZ(ResultsType &result)
    {TestMultiSpeciesWENO_AcousticPulseSlipAdiabatic3D_YZ(result,_communicator);};
*/

#ifdef USE_STELLA
#ifdef USE_PETSC
    virtual void Test__Efield__Single2D(ResultsType &result)
    {TestEfield_Single2D(result,_communicator);};
    virtual void Test__Efield__Single3D(ResultsType &result)
    {TestEfield_Single3D(result,_communicator);};
#ifdef USE_OVERKIT
    virtual void Test__EfieldOverset__Multiple2D(ResultsType &result)
    {TestEfieldOverset_Multiple2D(result,_communicator);};
    virtual void Test__EfieldOverset__Multiple3D(ResultsType &result)
    {TestEfieldOverset_Multiple3D(result,_communicator);};
    virtual void Test__EfieldOverset__FourGrid(ResultsType &result)
    {TestEfieldOverset_FourGrid(result,_communicator);};
#endif
#endif
#endif
    
    virtual void Test__PlasCom2__Poiseuille2DX__Regression(ResultsType &result)
    {TestPlasCom2_Poiseuille2DX_Regression(result,_communicator);};
    virtual void Test__PlasCom2__Poiseuille2DY__Regression(ResultsType &result)
    {TestPlasCom2_Poiseuille2DY_Regression(result,_communicator);};
    virtual void Test__PlasCom2__Poiseuille3DZ__Regression(ResultsType &result)
    {TestPlasCom2_Poiseuille3DZ_Regression(result,_communicator);};

    virtual void Test__PlasCom2__ViscousShock2DX__Regression(ResultsType &result)
    {TestPlasCom2_ViscousShock2DX_Regression(result,_communicator);};
    virtual void Test__PlasCom2__ViscousShock2DY__Regression(ResultsType &result)
    {TestPlasCom2_ViscousShock2DY_Regression(result,_communicator);};
    virtual void Test__PlasCom2__ViscousShock3DZ__Regression(ResultsType &result)
    {TestPlasCom2_ViscousShock3DZ_Regression(result,_communicator);};
    virtual void Test__PlasCom2__ViscousShockAD3DZ__Regression(ResultsType &result)
    {TestPlasCom2_VSCurvilinearAD3DZ_Regression(result,_communicator);};

    virtual void Test__PlasCom2__VSCurvilinear2DX__Regression(ResultsType &result)
    {TestPlasCom2_VSCurvilinear2DX_Regression(result,_communicator); };
    virtual void Test__PlasCom2__VSRectilinear2DX__Regression(ResultsType &result)
    {TestPlasCom2_VSRectilinear2DX_Regression(result,_communicator); };
    virtual void Test__PlasCom2__VSCurvilinear2DY__Regression(ResultsType &result)
    {TestPlasCom2_VSCurvilinear2DY_Regression(result,_communicator); };
    virtual void Test__PlasCom2__VSRectilinear2DY__Regression(ResultsType &result)
    {TestPlasCom2_VSRectilinear2DY_Regression(result,_communicator); };
    virtual void Test__PlasCom2__VSCurvilinear3DZ__Regression(ResultsType &result)
    {TestPlasCom2_VSCurvilinear3DZ_Regression(result,_communicator); };
    virtual void Test__PlasCom2__VSRectilinear3DZ__Regression(ResultsType &result)
    {TestPlasCom2_VSRectilinear3DZ_Regression(result,_communicator); };

    virtual void Test__PlasCom2__PFCurvilinear2DX__Regression(ResultsType &result)
    {TestPlasCom2_PFCurvilinear2DX_Regression(result,_communicator); };
    virtual void Test__PlasCom2__PFRectilinear2DX__Regression(ResultsType &result)
    {TestPlasCom2_PFRectilinear2DX_Regression(result,_communicator); };
    virtual void Test__PlasCom2__PFCurvilinear2DY__Regression(ResultsType &result)
    {TestPlasCom2_PFCurvilinear2DY_Regression(result,_communicator); };
    virtual void Test__PlasCom2__PFRectilinear2DY__Regression(ResultsType &result)
    {TestPlasCom2_PFRectilinear2DY_Regression(result,_communicator); };
    virtual void Test__PlasCom2__PFCurvilinear3DZ__Regression(ResultsType &result)
    {TestPlasCom2_PFCurvilinear3DZ_Regression(result,_communicator); };
    virtual void Test__PlasCom2__PFRectilinear3DZ__Regression(ResultsType &result)
    {TestPlasCom2_PFRectilinear3DZ_Regression(result,_communicator); };
    
    
    virtual void Test__PlasCom2__AcousticPulse2D__Periodic(ResultsType &result)
    {TestPlasCom2_AcousticPulsePeriodic2D(result,_communicator);};
    virtual void Test__PlasCom2__AcousticPulse2D__Farfield(ResultsType &result)
    {TestPlasCom2_AcousticPulseFarfield2D(result,_communicator);};
    virtual void Test__PlasCom2__AcousticPulse2D__Slipwall(ResultsType &result)
    {TestPlasCom2_AcousticPulseSlipwall2D(result,_communicator);};
    virtual void Test__PlasCom2__AcousticPulse2D__Sponge(ResultsType &result)
    {TestPlasCom2_AcousticPulseSponge2D(result,_communicator);};
    virtual void Test__PlasCom2__AcousticPulse3D__Periodic(ResultsType &result)
    {TestPlasCom2_AcousticPulsePeriodic3D(result,_communicator);};
    virtual void Test__PlasCom2__AcousticPulse3D__Farfield(ResultsType &result)
    {TestPlasCom2_AcousticPulseFarfield3D(result,_communicator);};
    virtual void Test__PlasCom2__AcousticPulse3D__Slipwall(ResultsType &result)
    {TestPlasCom2_AcousticPulseSlipwall3D(result,_communicator);};
    virtual void Test__PlasCom2__AcousticPulse3D__Sponge(ResultsType &result)
    {TestPlasCom2_AcousticPulseSponge3D(result,_communicator);};
    virtual void Test__PlasCom2__Rectilinear2D__Periodic(ResultsType &result)
    {TestPlasCom2_APRectilinearPeriodic2D(result,_communicator);};
    virtual void Test__PlasCom2__Curvilinear2D__Periodic(ResultsType &result)
    {TestPlasCom2_APCurvilinearPeriodic2D(result,_communicator);};
    virtual void Test__PlasCom2__Rectilinear3D__Periodic(ResultsType &result)
    {TestPlasCom2_APRectilinearPeriodic3D(result,_communicator);};
    virtual void Test__PlasCom2__Curvilinear3D__Periodic(ResultsType &result)
    {TestPlasCom2_APCurvilinearPeriodic3D(result,_communicator);};
    virtual void Test__PlasCom2__Rectilinear2D__Farfield(ResultsType &result)
    {TestPlasCom2_APRectilinearFarfield2D(result,_communicator);};
    virtual void Test__PlasCom2__Curvilinear2D__Farfield(ResultsType &result)
    {TestPlasCom2_APCurvilinearFarfield2D(result,_communicator);};
    virtual void Test__PlasCom2__Rectilinear3D__Farfield(ResultsType &result)
    {TestPlasCom2_APRectilinearFarfield3D(result,_communicator);};
    virtual void Test__PlasCom2__Curvilinear3D__Farfield(ResultsType &result)
    {TestPlasCom2_APCurvilinearFarfield3D(result,_communicator);};
    virtual void Test__PlasCom2__Rectilinear2D__Slipwall(ResultsType &result)
    {TestPlasCom2_APRectilinearSlipwall2D(result,_communicator);};
    virtual void Test__PlasCom2__Curvilinear2D__Slipwall(ResultsType &result)
    {TestPlasCom2_APCurvilinearSlipwall2D(result,_communicator);};
    virtual void Test__PlasCom2__Rectilinear3D__Slipwall(ResultsType &result)
    {TestPlasCom2_APRectilinearSlipwall3D(result,_communicator);};
    virtual void Test__PlasCom2__Curvilinear3D__Slipwall(ResultsType &result)
    {TestPlasCom2_APCurvilinearSlipwall3D(result,_communicator);};
    virtual void Test__PlasCom2__Rectilinear2D__Sponge(ResultsType &result)
    {TestPlasCom2_APRectilinearSponge2D(result,_communicator);};
    virtual void Test__PlasCom2__Curvilinear2D__Sponge(ResultsType &result)
    {TestPlasCom2_APCurvilinearSponge2D(result,_communicator);};
    virtual void Test__PlasCom2__Rectilinear3D__Sponge(ResultsType &result)
    {TestPlasCom2_APRectilinearSponge3D(result,_communicator);};
    virtual void Test__PlasCom2__Curvilinear3D__Sponge(ResultsType &result)
    {TestPlasCom2_APCurvilinearSponge3D(result,_communicator);};
    virtual void Test__PlasCom2__Restart(ResultsType &result)
    {TestPlasCom2_Restart(result,_communicator);};
#ifdef USE_OVERKIT
    virtual void Test__Overkit__AcousticPulse2D__BoxInBoxInjection(ResultsType &result)
    {TestOverkit_AcousticPulseBoxInBoxInjection2D(result,_communicator);};
#endif
    virtual void Test__WENO__Convergence(ResultsType &result)
    {TestWENO_Convergence(result,_communicator);};
    virtual void Test__PlasCom2__PooyaWENO(ResultsType &result)
    {TestPlasCom2_PooyaWENO(result,_communicator); };
    virtual void Test__PlasCom2__PavelAxisymmetryDriver(ResultsType &result)
    {TestPlasCom2_PavelAxisymmetryDriver(result,_communicator); };
    virtual void Test__PlasCom2__PavelWENO(ResultsType &result)
    {TestPlasCom2_PavelWENO(result,_communicator); };
    virtual void Test__PlasCom2__PavelSBP(ResultsType &result)
    {TestPlasCom2_PavelSBP(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWavy3(ResultsType &result)
    {TestPlasCom2_PooyaWavy3(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENO1(ResultsType &result)
    {TestPlasCom2_PooyaWENO1(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENO1x(ResultsType &result)
    {TestPlasCom2_PooyaWENO1x(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENOOx(ResultsType &result)
    {TestPlasCom2_PooyaWENOOx(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENO1y(ResultsType &result)
    {TestPlasCom2_PooyaWENO1y(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENO1z(ResultsType &result)
    {TestPlasCom2_PooyaWENO1z(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENO2(ResultsType &result)
    {TestPlasCom2_PooyaWENO2(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENO3(ResultsType &result)
    {TestPlasCom2_PooyaWENO3(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENO4(ResultsType &result)
    {TestPlasCom2_PooyaWENO4(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENO5(ResultsType &result)
    {TestPlasCom2_PooyaWENO5(result,_communicator); };
    virtual void Test__PlasCom2__PooyaWENOAxisymmetric(ResultsType &result)
    {TestPlasCom2_PooyaWENOAxisymmetric(result,_communicator); };
    virtual void Test__PlasCom2__Pooyanoslipisop1(ResultsType &result)
    {TestPlasCom2_Pooyanoslipisop1(result,_communicator); };
    virtual void Test__PlasCom2__Pooyanoslipisop2(ResultsType &result)
    {TestPlasCom2_Pooyanoslipisop2(result,_communicator); };
    virtual void Test__PlasCom2__Pooyanoslipisop3(ResultsType &result)
    {TestPlasCom2_Pooyanoslipisop3(result,_communicator); };
    virtual void Test__PlasCom2__Pooyanoslipisom1(ResultsType &result)
    {TestPlasCom2_Pooyanoslipisom1(result,_communicator); };
    virtual void Test__PlasCom2__Pooyanoslipisom2(ResultsType &result)
    {TestPlasCom2_Pooyanoslipisom2(result,_communicator); };
    virtual void Test__PlasCom2__Pooyanoslipisom3(ResultsType &result)
    {TestPlasCom2_Pooyanoslipisom3(result,_communicator); };
    virtual void Test__PlasCom2__PooyaEnthalpyp1(ResultsType &result)
    {TestPlasCom2_PooyaEnthalpyp1(result,_communicator); };
    virtual void Test__PlasCom2__PooyaEnthalpyp2(ResultsType &result)
    {TestPlasCom2_PooyaEnthalpyp2(result,_communicator); };

    virtual void Test__PlasCom2__Pooyanoslipisop3__Regression(ResultsType &result)
    {TestPlasCom2_Pooyanoslipisop3_Regression(result,_communicator); };
    virtual void Test__PlasCom2__Pooyanoslipisom3__Regression(ResultsType &result)
    {TestPlasCom2_Pooyanoslipisom3_Regression(result,_communicator); };

    virtual void Test__Comm__Split(ResultsType &result)
    {TestComm_Split(result,_communicator); };

    virtual void Test__PCPP__SetupCartNeighborhood(ResultsType &result)
    {TestCartNeighborhood(result,_communicator); };

    virtual void Test__Various__PartitionIntervals(ResultsType &result)
    {TestPartitionIntervals(result,_communicator); };

    virtual void Test__Examples__Parallel(ResultsType &result)
    {ExamplesParallelTestDriver(result,_communicator); };

    ///
    /// Runs a list of tests specified by name
    /// 
    /// @param test_names list of string test names
    /// @param result PlasCom2::TestResults object to store test results.
    ///
    virtual void ProcessTests(std::list<std::string> &test_names,ResultsType &result){
      Prologue();
      std::list<std::string>::iterator tni = test_names.begin();
      while(tni != test_names.end())
        RunTest(*tni++,result);
      Epilogue();
    }

  };
#endif
};
#endif





