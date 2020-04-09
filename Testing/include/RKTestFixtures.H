#ifndef __RK_TEST_FIXTURES_H__
#define __RK_TEST_FIXTURES_H__

#include "Simulation.H"

namespace rktestfixtures {

typedef simulation::grid::halo                                   halo_t;
  
template<typename GridType,typename StateType,typename OperatorType>
class testrhs { // : public DomainT::RHSType {
public:
  //   typedef typename DomainT::GridType GridType;
  //   typedef typename DomainT::StateType StateType;
  //   typedef typename DomainT::OperatorType OperatorType;
  typedef int BCType;
public:
  static double TestF(double time,int order){
    return(std::pow(time,static_cast<double>(order)));
  };
  
  static double RHSF(double inTime,int inOrder)
  {
    double coeff   = static_cast<double>(inOrder);
    double timePow = std::pow(inTime,static_cast<double>(inOrder-1));
    return(coeff*timePow);
  };

  virtual int SetupRHS(GridType &inGrid,StateType &inState,StateType &rhsState,int myThreadId = 0)
  {
    myGrid     = &inGrid;
    myState    = &inState;
    myRHSState = &rhsState;
    myDT = 0;
    return(0); 
  };
  void SetTimestep(double inDT) { myDT = inDT; };
  virtual int Initialize(GridType &inGrid,StateType &inState)
  { return(0);};
  virtual int RHS(int threadId = 0)
  {
    return(ProcessInterior(*myGrid,*myState,*myRHSState)); 
  };
  virtual int RHSFinalize(StateType &inState,int threadId) { return(0); };
  virtual int ProcessInterior2(GridType &inGrid,StateType &inState,StateType &rhsState)
  { return(0); };
  virtual int ProcessInterior(GridType &inGrid,StateType &inState,StateType &rhsState,int threadId = 1)
  {
    
    pcpp::field::dataset::DataBufferType &timeData(inState.Field("simTime"));
    double *inTime = timeData.Data<double>();

    pcpp::field::dataset::DataBufferType &aData(inState.Field("A"));
    double *inA    = aData.Data<double>();

    std::cout << "testrhs::ProcessInterior: Time = " << *inTime << std::endl;
    std::cout << "testrhs::ProcessInterior: A    = " << *inA << std::endl;
    std::cout << "testrhs::Address of time: " << inTime << std::endl;

    double *outF0  = rhsState.GetStateFieldData("F0");
    double *outF1  = rhsState.GetStateFieldData("F1");
    double *outF2  = rhsState.GetStateFieldData("F2");
    double *outF3  = rhsState.GetStateFieldData("F3");
    double *outF4  = rhsState.GetStateFieldData("F4");
    double *outF5  = rhsState.GetStateFieldData("F5");
    double *outF6  = rhsState.GetStateFieldData("F6");

    *outF0 = *inA*RHSF(*inTime,0);
    *outF1 = *inA*RHSF(*inTime,1);
    *outF2 = *inA*RHSF(*inTime,2);
    *outF3 = *inA*RHSF(*inTime,3);
    *outF4 = *inA*RHSF(*inTime,4);
    *outF5 = *inA*RHSF(*inTime,5);
    *outF6 = *inA*RHSF(*inTime,6);

    std::cout << "testrhs::Out0 = " << *outF0 << std::endl
              << "testrhs::Out1 = " << *outF1 << std::endl
              << "testrhs::Out2 = " << *outF2 << std::endl
              << "testrhs::Out3 = " << *outF3 << std::endl
              << "testrhs::Out4 = " << *outF4 << std::endl
              << "testrhs::Out5 = " << *outF5 << std::endl
              << "testrhs::Out6 = " << *outF6 << std::endl;

    return(0);
  };
  virtual double TimeStep(double *inCFL = NULL) { return(myDT); };
  virtual int ProcessHalo(halo_t &inHalo,GridType &inGrid,StateType &inState,StateType &rhsState,int threadId = 1)
  { return(0); };
  virtual int ProcessBoundary(GridType &inGrid,StateType &inState,StateType &rhsState,int threadId = 1)
  { return(0); };
protected:
  StateType *myState;
  StateType *myRHSState;
  GridType  *myGrid;
  double myDT;
};

//typedef testrhs RHSHandlerType;

}

#endif
