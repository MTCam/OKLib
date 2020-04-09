#ifndef __EFIELD_TEST_FIXTURES_H__
#define __EFIELD_TEST_FIXTURES_H__

#include "Simulation.H"

namespace efieldtestfixtures {

  typedef simulation::grid::halo halo_t;

  template <typename GridT, typename StateT> class testbc {

  public:

    typedef GridT GridType;
    typedef simulation::grid::subregion SubRegionType;
    typedef StateT StateType;
    typedef typename StateType::MetaDataSetType DataDictionaryType;

  protected:

    int         bcType;
    std::string bcName;

    StateType  *paramPtr;
    StateType  *statePtr;

    const double     *paramData;
    const int        *flagData;

    DataDictionaryType     stateDictionary;
    DataDictionaryType     paramDictionary;
    
    bool myParams;
    bool myState;

    GridType      *gridPtr;
    SubRegionType *subRegionPtr;

    const double  *bcParams;
    const int     *bcFlags;

    std::vector<double *> myRHSBuffers;
    std::vector<double *> myStateBuffers;
    std::vector<double *> myTargetBuffers;

    int numDim;
    int normalDir;
    int dirNormal;
    int signNormal;
    int dirIndex;
    size_t numPointsBuffer;
    size_t numPointsRegion;

  public:

    testbc() : bcType(0), paramPtr(NULL), statePtr(NULL), 
               myParams(false), myState(false),
               paramData(NULL), flagData(NULL),
               gridPtr(NULL), subRegionPtr(NULL),
               bcParams(NULL), bcFlags(NULL)
    {
      Create();
    };
    
    testbc(const testbc &inBC) : bcType(0), paramPtr(NULL), statePtr(NULL), 
                               myParams(false), myState(false),
                               paramData(NULL), flagData(NULL),
                               gridPtr(NULL), subRegionPtr(NULL),
                               bcParams(NULL), bcFlags(NULL)
    {
      Copy(inBC); 
    };

    testbc &operator=(const testbc &inBC)
    {
      Copy(inBC);
      return(*this);
    }

    void Copy(const testbc &inBC)
    {

      Destroy();

      bcType          = inBC.bcType;
      paramPtr        = inBC.paramPtr;
      statePtr        = inBC.statePtr;
      paramData       = inBC.paramData;
      flagData        = inBC.flagData;
      bcParams        = inBC.bcParams;
      bcFlags         = inBC.bcFlags;
      numDim          = inBC.numDim;

      gridPtr            = inBC.gridPtr;
      subRegionPtr       = inBC.subRegionPtr;
      numPointsBuffer    = inBC.numPointsBuffer;
      numPointsRegion    = inBC.numPointsRegion;

      normalDir       = inBC.normalDir;
      dirNormal       = inBC.dirNormal;
      signNormal      = inBC.signNormal;
      dirIndex        = inBC.dirIndex;

      myState  = false;
      myParams = false;

      int numBuffers = inBC.myRHSBuffers.size();
      myRHSBuffers.resize(numBuffers);
      for(int iBuffer = 0;iBuffer < numBuffers;iBuffer++)
        myRHSBuffers[iBuffer] = inBC.myRHSBuffers[iBuffer];
     
      numBuffers = myStateBuffers.size();
      myStateBuffers.resize(numBuffers);
      for(int iBuffer = 0;iBuffer < numBuffers;iBuffer++)
        myStateBuffers[iBuffer] = inBC.myStateBuffers[iBuffer];
      
      numBuffers = myTargetBuffers.size();
      myTargetBuffers.resize(numBuffers);
      for(int iBuffer = 0;iBuffer < numBuffers;iBuffer++)
        myTargetBuffers[iBuffer] = inBC.myTargetBuffers[iBuffer];

    };

    void SetParams(StateType &inParams)
    {
      DestroyParams();
      paramPtr = &inParams;
    };

    void SetState(StateType &inParams)
    {
      DestroyState();
      statePtr = &inParams;
    };
    
    void SetFlagData(const int *inFlags)
    { flagData = inFlags; };

    void SetParamData(const double *inParams)
    { paramData = inParams; };

    void DestroyParams()
    {
      if(myParams && (paramPtr != NULL)){
        delete paramPtr;
      }
      paramPtr = NULL;
      myParams = false;
    };

    void DestroyState()
    {
      if(myState && (statePtr != NULL)){
        delete statePtr;
      }
      statePtr = NULL;
      myState = false;
    };
    
    void Destroy()
    {
      DestroyParams();
      DestroyState();
    };
    
    void Create()
    {
      Destroy();
      paramPtr = new StateType;
      statePtr = new StateType;
      myParams = true;
      myState  = true;
    };
    
    const double *ParameterData() const
    { return(paramData); };
    
    double *ParameterData() { return(const_cast<double *>(paramData)); };

    const int *FlagData() const
    { return(flagData); };
    
    int *FlagData() { return(const_cast<int *>(flagData)); };
    
    StateType &Params() { return(*paramPtr); };
    const StateType &Params() const { return(*paramPtr); };
    
    int BCType() const { return(bcType); };
    void SetBCType(int inType) { bcType = inType; };

    std::string &BCName() { return(bcName); };
    const std::string &BCName() const { return(bcName); };
    void SetBCName(const std::string &inName){ bcName = inName; };
    
    DataDictionaryType &ParamDictionary() { return(paramDictionary); };
    const DataDictionaryType &ParamDictionary() const { return(paramDictionary); };
    DataDictionaryType &StateDictionary() { return(stateDictionary); };
    const DataDictionaryType &StateDictionary() const { return(stateDictionary); };
    
    int InitInternalData(std::ostream &messageStream) 
    {
      return(0);
  //     if(paramPtr){
  //       int dataIndex = paramPtr->GetDataIndex("Numbers");
  //       if(dataIndex >= 0){
  //         paramData = paramPtr->GetFieldBuffer<double>(dataIndex);
  //       }
  //       dataIndex = paramPtr->GetDataIndex("Flags");
  //       if(dataIndex >= 0){
  //         flagData = paramPtr->GetFieldBuffer<int>(dataIndex);
  //       }
  //     }
    };

    int SetParams(const double *inParams,
                  const int    *inFlags)
    {
      bcParams = inParams;
      bcFlags  = inFlags;
      return(0);
    }

  };
    
  template<typename GridType,typename StateType,typename OperatorType>
  class testrhs { // : public DomainT::RHSType {
  public:
    //   typedef typename DomainT::GridType GridType;
    //   typedef typename DomainT::StateType StateType;
    //   typedef typename DomainT::OperatorType OperatorType;
    typedef testbc<GridType, StateType> BCType;
  public:

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
    { return(0); };
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

}

#endif
