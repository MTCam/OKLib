#ifndef __TEST_DOMAIN_H__
#define __TEST_DOMAIN_H__

#include "PCPPTypes.H"
#include "PCPPUtil.H"
#include "Grid.H"
#include "TestBoundary.H"
#include "TestBoundaryCondition.H"
#include "OperatorTypes.H"
#include "ConfigUtil.H"
#include "NavierStokesBC.H"
#ifdef HAVE_CANTERA
#include "cantera/IdealGasMix.h"
#endif
#include "Geometry.H"
#include "NavierStokesInitialization.H"
#include "NavierStokesOptions.H"
#include "Initializer.H"
#include "Info.H"
#include "PC2Util.H"
#include "InterGrid.H"
#include "EfieldSolver.H"
#include "GridTimers.H"

#ifdef USE_OVERKIT
#include "overkit.h"
#endif

#ifdef USE_OVERKIT_FORTRAN
#include "OverkitFortranInterface.H"
#endif

using namespace pcpp::util;
using namespace fixtures;

namespace simulation {

namespace testdomain {

#define CONSTANT_DT  0
#define CONSTANT_CFL 1
  
  struct info {
    info() : numDomains(0) {};
    int numDomains;
    std::vector<std::string> domainNames;
    int numLocalDomains;
    std::vector<int> globalDomainIndex;
  };
  
  struct partitioninfo {
    std::vector<int>    numProcsGrid;
    std::vector<size_t> numPointsProc;
    std::vector<int>    gridColor;
    /// The global grid index of each local grid (has size of numLocalGrids)
    std::vector<int>                        gridIndexOfLocalGrids;
  };
  

  template<typename GridT,typename StateT,typename RhsT>
  class base {

    template<typename DomainT>
    friend class initializer::base; // DomainInitializer;
    
  public:

    typedef GridT                            GridType;
    typedef StateT                           StateType;
    typedef RhsT                             RHSType;
    typedef typename GridType::HaloType      HaloType;
    typedef intergrid::intergrid<GridType>   InterGridType;

    // some fixed types, smells bad
    typedef typename RHSType::BCType            BCType;
    typedef testboundary::base<BCType>              BoundaryType;
    typedef typename StateType::MetaDataSetType DataDictionaryType;
    typedef operators::sbp_operator_t           OperatorType;
    typedef geometry::info                      GeometryInfoType;
    typedef geometry::gridinfo                  GridInfoType;
    typedef simulation::grid::subregion         GridRegionType;
    typedef std::vector<GridType*>     gridvector;
    typedef std::vector<StateType*>    statevector;
    typedef std::vector<RHSType*>      rhsvector;

#ifdef HAVE_CANTERA
    typedef Cantera::IdealGasMix CanteraType;
#endif

    typedef efield::efield_solver<GridType, StateType, BCType> EfieldSolverType;
    typedef efield::efield<GridType, StateType> EfieldType;

  public:

    base() : numGrids(0), rhsPtr(&dummyRHS), domainTime(0), stepNumber(0), domainRank(0),
             geometryInfoPtr(NULL), interGridPtr(NULL), usesOverkit(false), virtualNodePtr(NULL),
             efieldSolverPtr(NULL), efieldSolveInterval(0), runMode(pcpp::runmode::SPMD)
    { };

    StateType   &State(int iGrid)  { return (*gridStates[iGrid]);  };
    StateType   &Target(int iGrid) { return (*gridTargets[iGrid]); };
    StateType   &Param(int iGrid)  { return (*gridParams[iGrid]);  };
    GridType    &Grid(int iGrid)   { return (*domainGrids[iGrid]); };
    RHSType     &RHS(int iGrid)    { return (*gridRHSs[iGrid]);    };

    gridvector  &Grids()          { return (domainGrids);         };
    statevector &States()         { return (gridStates);          };
    statevector &Targets()        { return (gridTargets);         };
    statevector &Params()         { return (gridParams);          };
    rhsvector   &RHSs()           { return (gridRHSs);            };

    std::vector<BoundaryType>  &GridBoundaries(int iGrid) { return (domainBoundaries[iGrid]); };
    std::vector<BoundaryType>  &DomainBoundary(int iGrid) { return (domainBoundaries[iGrid]); };
    std::vector<std::vector<BoundaryType> > &DomainBoundaries() { return (domainBoundaries); };
    
    const std::vector<std::string> &BCNames() const {return(bcNames); };
    const std::vector<std::string> &BoundaryNames() const {return(boundaryNames); };

    void SetVirtualNodes(pcpp::VirtualNodeType *inVNode)
    { virtualNodePtr = inVNode; };

#ifdef HAVE_CANTERA
    CanteraType *GetCanteraGasPtr() {return (canteraGasPtr);};
#endif
    //std::vector<std::string> GetSpeciesNames() {return speciesNames;};

    //    const std::vector<std::string> &BoundaryNames() const {return(boundaryNames); };

    bool HasInterGrid() { return interGridPtr != NULL; }
    InterGridType &InterGrid() { return *interGridPtr; }

    void SetGrid(int gridIndex,GridType &inGrid)
    {
      assert(gridIndex >= 0 && gridIndex < domainGrids.size());
      domainGrids[gridIndex] = &inGrid;
    };

    void SetGridState(int gridIndex,StateType &inState)
    {
      assert(gridIndex >= 0 && gridIndex < domainGrids.size());
      gridStates[gridIndex] = &inState;
    };

    void SetGridTarget(int gridIndex,StateType &inState)
    {
      assert(gridIndex >= 0 && gridIndex < domainGrids.size());
      gridTargets[gridIndex] = &inState;
    };

    void SetGridParams(int gridIndex,StateType &paramState)
    {
      assert(gridIndex >= 0 && gridIndex < domainGrids.size());
      gridParams[gridIndex] = &paramState;
    };
    
    void SetName(const std::string &name) { domainName = name; }
    std::string Name() const { return(domainName); };

    void SetGlobal(pcpp::ParallelGlobalType &inGlobal)
    {

      globalPtr = &inGlobal;

      if (interGridPtr) {
        interGridPtr->SetGlobal(inGlobal);
      }

    };

    int SetNumGrids(int numberOfGrids)
    {
      if(numberOfGrids <= 0)
        return(1);
      numGrids = numberOfGrids;
      ownGrid.resize(numGrids,false);
      ownState.resize(numGrids,false);
      ownRHS.resize(numGrids,false);
      domainGrids.resize(numGrids,NULL);
      gridStates.resize(numGrids,NULL);
      gridTargets.resize(numGrids,NULL);
      gridParams.resize(numGrids,NULL);
      gridRHSs.resize(numGrids,NULL);
      domainBoundaries.resize(numGrids);
      //      dummyComms.resize(numGrids);
      gridCommunicators.resize(numGrids,domainCommunicator);
      gridIOCommunicators.resize(numGrids,domainCommunicator);
      gridWorldCommunicators.resize(numGrids,domainCommunicator);
      std::vector<double> gridTimes(::timers::grid::NUMTIMERS,0.0);
      std::vector<size_t> gridCalls(::timers::grid::NUMTIMERS,0);
      gridTimers.resize(numGrids,gridTimes);
      gridCallCounts.resize(numGrids,gridCalls);

      // for(int iGrid = 0;iGrid < numGrids;iGrid++){
      // 	dummyComms[iGrid] = MPI_COMM_NULL;
      // 	gridCommunicators[iGrid].SetCommunicator(dummyComms[iGrid]);
      // }
      return(0);
    };

    int ConfigureDomain(fixtures::ConfigurationType &inConfig,std::ostream &inMessageStream)
    {

      SetMessageStream(inMessageStream);
      std::ostream &messageStream(*messageStreamPtr);

      domainConfig = inConfig;
      domainName = domainConfig.GetValue("ConfigName");

      if(domainName.empty()){
        messageStream << "Name of domain unspecified." << std::endl;
        return(1);
      }

      spatialOrder = domainConfig.GetValue<int>(ConfigKey(domainName,"SpatialOrder"));
      messageStream << "Setting up operators for spatial order " << spatialOrder << "."
                    << std::endl;
      int interiorOrder  = (spatialOrder - 1)*2;
      messageStream << "Setting up SBP"
                    << interiorOrder/2 <<  interiorOrder << " operator."
                    << std::endl;

      operators::sbp::Initialize(domainOperator,interiorOrder);
   
      // fortranification
      int *stencilStarts = domainOperator.stencilStarts;
      int numStencils    = domainOperator.numStencils;
      for(int iStencil = 0;iStencil < numStencils;iStencil++)
        stencilStarts[iStencil]++;
      
      

      if(domainConfig.IsSet(ConfigKey(domainName,"ScalarNames")))
        scalarNames = domainConfig.GetValueVector<std::string>(ConfigKey(domainName,"ScalarNames"));
      
      int numGridsConfigured = 0;
      gridNames = domainConfig.GetValueVector<std::string>(ConfigKey(domainName,"GridNames"));
      if(!domainConfig.IsSet(ConfigKey(domainName,"NumGrids"))){
        if(gridNames.empty())
          messageStream << "Number of grids not specified" << std::endl;
        else numGridsConfigured = gridNames.size();
      } else {
        numGridsConfigured = domainConfig.GetValue<int>(ConfigKey(domainName,"NumGrids"));
      }
      messageStream << "Domain found " << numGridsConfigured << " grids." << std::endl;

      if(numGridsConfigured < 0){
        messageStream << "Number of grids cannot be negative." << std::endl;
        return(1);
      } else if (numGridsConfigured == 0){
        messageStream << "Gridless domains are not supported (yet)." << std::endl;
        return(1);
      }

      SetNumGrids(numGridsConfigured);

      bool nativeGeometry = true;
      std::vector<std::string>::iterator gridNameIt = gridNames.begin();
      while(gridNameIt != gridNames.end()){
        int gridIndex = gridNameIt - gridNames.begin();
        std::string &gridName(*gridNameIt++);
        std::string geometryName;
        std::string::size_type geomPos = gridName.find(":");
        if(geomPos != std::string::npos){
          geometryName = gridName.substr(0,geomPos);
        }
        geometryNames.push_back(geometryName);
      }
      
      if(ConfigureDataDictionary(domainConfig,"State",stateDataDictionary,messageStream)){
        messageStream << "Domain state meta failed to configure."
                      << std::endl;
        return(1);
      }
      //       std::cout << "Configured state dictionary:" << std::endl
      //                 << stateDataDictionary << std::endl;
      if(ConfigureDataDictionary(domainConfig,"Param",paramDataDictionary,messageStream)){
        messageStream << "Domain param meta failed to configure."
                      << std::endl;
        return(1);
      }
      
      // Allowing use of named options to set the bitvector on the inside
      // This is internal-only data that is set up based on user options
      paramDataDictionary.AddField("optionFlags",'g',1,4,"");
      //       std::cout << "Configured state dictionary:" << std::endl
      //                 << paramDataDictionary << std::endl;
      
      if(domainConfig.IsSet(ConfigKey(domainName,"IO:Fields")))
        ioFieldNames = domainConfig.GetValueVector<std::string>(ConfigKey(domainName,"IO:Fields"));
      if(domainConfig.IsSet(ConfigKey(domainName,"Restart:Fields")))
        restartFieldNames = domainConfig.GetValueVector<std::string>(ConfigKey(domainName,"Restart:Fields"));      

      // Check which external tools are needed
      usesOverkit = false;
      usesPETSc = false;
      usesStella = false;
    
      // Check if Overkit is required
      if (inConfig.IsSet(ConfigKey(domainName, "InterGrid"))) {
        std::string interGridMode = inConfig.GetValue(ConfigKey(domainName, "InterGrid"));
        if (interGridMode == "Overkit") {
          usesOverkit = true;
          // set up dictionary entry if it does not exist, otherwise accept user's dictionary
          stateDataDictionary.ResetField("ovkState",'n',1,4,"OverkitData",false);
        }
      }

      if (domainConfig.IsSet(ConfigKey(domainName, "Efield:GridNames"))) {
        std::vector<std::string> efieldGridNames = domainConfig.GetValueVector(ConfigKey(domainName,
          "Efield:GridNames"));
        if (efieldGridNames.size() > 0) {
          usesPETSc = true;
          usesStella = true;
        }
      }

      PopulateDomainStateMeta();
      PopulateDomainParamMeta();

      if(ConfigureDomainBCs(domainConfig, messageStream)){
        messageStream << "Failed to configure BCs for domain " << domainName << "." << std::endl;
        return(1);
      }
      // initialize domain specific third party libraries (i.e. Cantera)
#ifdef HAVE_CANTERA
      if(domainConfig.IsSet(ConfigKey(domainName,"CanteraInputFile"))){
        
        std::string canteraFileName = 
          domainConfig.GetValue<std::string>(ConfigKey(domainName,"CanteraInputFile"));

        if(ix::sys::FILEEXISTS(canteraFileName)){
          //std::string canteraInputString;
          if(ReadCanteraInput(canteraFileName,canteraInputString)){
            messageStream << "Failed to successfully read Cantera input file " 
                          << canteraFileName << std::endl;
            return(1);
          }
          
          //std::string canteraId;
          if(domainConfig.IsSet(ConfigKey(domainName,"CanteraId"))){
            canteraId = domainConfig.GetValue<std::string>(ConfigKey(domainName,"CanteraId"));
          } else {
            messageStream << "canteraId not specified in input file. Defaulting to gas." << std::endl;
            canteraId = "gas";
          }
          
          /// @todo find a way to check that this happens for gasModel=1
          /// @todo add a check for number of scalars and number of cantera species
          // Create a Cantera Object, sets the canteraGasPtr memeber data
          // Note this is a specific type (IdealGasMix), for use with Chemistry
          // If a different type of object is needed, this will need to be more general
          if(CreateCanteraGasPhase(canteraInputString,canteraId)){
            messageStream << "Failed to successfully create the domain Cantera object " 
                          << canteraFileName << std::endl;
          }
        } else {
          messageStream << "Warning: Cantera input file (" << canteraFileName
                        << ") not found. Cantera setup skipped." << std::endl;
        }
      }
#endif

      return(0);
    }

    int ConfigureDomainBCs(fixtures::ConfigurationType &inConfig, std::ostream &messageStream) {

      bcNames = inConfig.GetValueVector<std::string>(ConfigKey(domainName,"BCNames"));
      int numBoundaries = bcNames.size();
      domainBCs.reserve(numBoundaries);
      domainBCs.resize(numBoundaries);

      //std::cout << "ConfigureDomainBCs" << std::endl;

      //      std::cout << "Number of boundary condition types: " << numBoundaries << std::endl;
      for(int iBoundary = 0;iBoundary < numBoundaries;iBoundary++){


        std::string &bcName(bcNames[iBoundary]);
        //        std::cout << "BCname = " << bcName << std::endl;
        std::string boundaryConfigName(ConfigKey(domainName,bcName));
        //         std::cout << "boundaryType (" << bcName << ") has bc config name: ("
        //                   << boundaryConfigName << ")" << std::endl;

        //std::cout << "bcName " << bcName << std::endl;

        fixtures::ConfigurationType 
          boundaryConfig(pcpp::util::ExtractConfigParams(inConfig,boundaryConfigName));
        
        std::string bcTypeName(boundaryConfig.GetValue(ConfigKey(boundaryConfigName,"BCType")));
        int bcType = testboundary::bc::ResolveBCName(bcTypeName);
        if(bcType == testboundary::bc::INVALIDBCTYPE){
          std::cerr << "INVALID BC!!" << std::endl;
          return(1);
        }
        BCType tempBC;
        BCType &myBC(domainBCs[iBoundary]);
        myBC = tempBC;
        myBC.SetBCName(bcTypeName);
        myBC.SetBCType(bcType);
        //         std::cout << "Setting BC Type Name (" 
        //                   << &myBC << ") = " << bcTypeName << std::endl;

        if(bcType != testboundary::bc::HOLE){ 
          
          std::ostringstream localMessageStream;
          //if(ConfigureDataDictionary(boundaryConfig,"State",myBC.StateDictionary(),messageStream)){
          if(ConfigureDataDictionary(boundaryConfig,"State",myBC.StateDictionary(),localMessageStream)){
            messageStream << "BC state dictionary failed to configure for BC " << bcName 
                          << std::endl;
            messageStream << localMessageStream.str();
            return(1);
          }
          pcpp::io::RenewStream(localMessageStream);

          if(ConfigureDataDictionary(boundaryConfig,"Param",myBC.ParamDictionary(),localMessageStream)){
            messageStream << "BC param dictionary failed to configure for BC " << bcName
                          << std::endl;
            messageStream << localMessageStream.str();
            return(1);
          }
          pcpp::io::RenewStream(localMessageStream);

          // inspect the dicitony for a specific boundry type
          

          if(myBC.ParamDictionary().empty() && myBC.StateDictionary().empty()){
            messageStream << "BC params and fields not configured for " << bcName
                          << std::endl;
            return(1);
          }
          
          // Populate any BC parameter metadata from config, and allocate
          StateType &bcParams(myBC.Params());
          bcParams.SetMetaData(myBC.ParamDictionary());
          bcParams.Create();
          
          // Populate parameter values from config file - any grid-connected state data done later
          std::string bcParamKey(ConfigKey(boundaryConfigName,"Param"));
          if(ConfigStateValues<StateType>(boundaryConfig,bcParamKey,bcParams,localMessageStream)){
            messageStream << "BC param values failed to populate from configuration." << std::endl;
            messageStream << localMessageStream.str();
            return(1);
          }
          pcpp::io::RenewStream(localMessageStream);

          //if(myBC.InitInternalData(localMessageStream)){
          int returnValue = myBC.InitInternalData(localMessageStream);

          //if(myBC.CheckInternalData(localMessageStream)){
          if (returnValue > 0){
            messageStream << "Initializing BC " << bcName << std::endl;
            messageStream << localMessageStream.str();
            if(returnValue > 1) 
              return returnValue;
          }
          pcpp::io::RenewStream(localMessageStream);

          // sanity check, each BC has a required number of data items that must be set
          // these are read in from the Numbers and Flags Param Fields
          
          
        } else {
          messageStream << "Deferring setup of hole region." << std::endl;
        }
        //        domainBCs.push_back(myBC);
      }

      // initialize domain specific third party libraries (i.e. Cantera)
#ifdef HAVE_CANTERA
      if(domainConfig.IsSet(ConfigKey(domainName,"CanteraInputFile"))){
        std::string canteraFileName = 
          domainConfig.GetValue<std::string>(ConfigKey(domainName,"CanteraInputFile"));
        if(ix::sys::FILEEXISTS(canteraFileName)){
         // std::string canteraInputString;
          if(ReadCanteraInput(canteraFileName,canteraInputString)){
            messageStream << "Failed to successfully read Cantera input file " 
                          << canteraFileName << std::endl;
            return(1);
          }
          
          //std::string canteraId;
          if(domainConfig.IsSet(ConfigKey(domainName,"CanteraId"))){
            canteraId = domainConfig.GetValue<std::string>(ConfigKey(domainName,"CanteraId"));
          } else {
            messageStream << "WARNING: canteraId not specified in input file. Defaulting to gas." << std::endl;
            canteraId = "gas";
          }
          
          /// @todo find a way to check that this happens for gasModel=1
          /// @todo add a check for number of scalars and number of cantera species
          // Create a Cantera Object, sets the canteraGasPtr memeber data
          // Note this is a specific type (IdealGasMix), for use with Chemistry
          // If a different type of object is needed, this will need to be more general
          if(CreateCanteraGasPhase(canteraInputString,canteraId)){
            messageStream << "Failed to successfully create the domain Cantera object " 
                          << canteraFileName << std::endl;
          }
        } else {
          messageStream << "Warning: Cantera input file (" << canteraFileName
                        << ") not found. Skipping Cantera setup." << std::endl;
        }
      }
#endif

      return(0);

    }
    
    int ConfigureGridDomainBoundaries(fixtures::ConfigurationType &inConfig,int iGrid,
                                      std::ostream &inMessageStream)
    {
      
      int numBoundaries = bcNames.size();
      std::ostream &messageStream(inMessageStream);
      
      const std::string &currentGridName(gridNames[iGrid]);
      //       std::cout <<  "Processing Grid(" << iGrid << "," 
      //                 << currentGridName << ") " << std::endl;
      for(int iBoundary = 0;iBoundary < numBoundaries;iBoundary++){
        
        std::string &bcName(bcNames[iBoundary]);
        //         std::cout << "BCName(" << iBoundary << ") = " << bcName << std::endl;
        BCType &boundaryBC(domainBCs[iBoundary]);
        //         std::cout << "BC Address (" << &boundaryBC << ")" << std::endl;
        std::string boundaryConfigKey(ConfigKey(domainName,bcName));
        //         std::cout << "Boundary config key:  " << boundaryConfigKey << std::endl;
        std::vector<int> boundaryGridIDs;
        
        // This is a vector of grid regions to which this boundary condition applies
        std::vector<std::string> boundaryGridRegionNames = 
          inConfig.GetValueVector<std::string>(ConfigKey(boundaryConfigKey,"RegionNames"));
        int numBoundaryGridRegions = boundaryGridRegionNames.size();
        //         std::cout << "NUmber of boundary grid regions: " << numBoundaryGridRegions << std::endl;
        std::vector<std::string>::iterator rnIt = boundaryGridRegionNames.begin();
        //         while(rnIt != boundaryGridRegionNames.end()){
        //           std::cout << "Applies to grid region: " << *rnIt++ << std::endl;
        //         }
        
        boundaryNames = 
          inConfig.GetValueVector<std::string>(ConfigKey(boundaryConfigKey,"BoundaryNames"));
        rnIt = boundaryNames.begin();
        //         while(rnIt != boundaryNames.end()){
        //           std::cout << "Applies to domain boundary name: " << *rnIt++ << std::endl;
        //         }
        if(boundaryNames.empty() || boundaryNames.size() != boundaryGridRegionNames.size()){
          messageStream << "WARNING: Setting default boundary names." << std::endl;
          boundaryNames.resize(numBoundaryGridRegions);
          for(int iRegion = 0;iRegion < numBoundaryGridRegions;iRegion++)
            boundaryNames[iRegion] = bcName + std::string(":") + boundaryGridRegionNames[iRegion]; 
        }
        // this boundary condition applies to (numBoundaryRegions) grid boundaries
        for(int iRegion = 0;iRegion < numBoundaryGridRegions;iRegion++){
          std::string &boundaryName(boundaryNames[iRegion]);
          std::string gridRegionName(boundaryGridRegionNames[iRegion]);
          //           std::cout << "IRegion(" << iRegion << ") = (" << boundaryName << "," << gridRegionName << ")" 
          //                     << std::endl;
          std::string::size_type x = gridRegionName.find_last_of(":");
          // std::string geometryName;
          std::string regionGridName;
          if(x != std::string::npos){
            regionGridName = gridRegionName.substr(0,x);
            gridRegionName = gridRegionName.substr(x+1);
          }
          if(x == std::string::npos || regionGridName.empty())
            regionGridName = currentGridName;
          std::string regionName = gridRegionName;
          int gridIndex = 0;
          //           std::cout << "RegionName: " << regionName << std::endl;
          //           std::cout << "regionGridName: " << regionGridName << std::endl;
          if(regionGridName != currentGridName){
            if(!regionGridName.empty()){ // make sure *some* grid region will match
              gridIndex = GetGridIndex(regionGridName);
              if(gridIndex < 0){
                messageStream << "Failed to find grid id for \""
                              << regionGridName << "\" in BC resolve."
                              << std::endl;
                return(1);
              }
            }
            continue;
          }
            
          gridIndex = iGrid;
            
          std::vector<BoundaryType> &domainGridBoundaries(GridBoundaries(gridIndex));
          std::string &bcName(boundaryBC.BCName());
          //           std::cout << "BCTypeName: " << bcName << std::endl;
          //  domainConfig.GetValue<std::string>(ConfigKey(boundaryConfigKey,"Type")));
          //  myBoundary.SetBCName(domainConfig.GetValue<std::string>(ConfigKey(boundaryConfigKey,"Type")));
          GridType &boundaryGrid(Grid(gridIndex));
          //           std::cout << "gridIndex = " << gridIndex << std::endl;
          //           std::cout << "Finding grid region for regionName = " << regionName
          //                     << std::endl; 
          int gridRegionIndex = boundaryGrid.GetRegionIndex(regionName);
	  //          std::cout << "gridRegionIndex = " << gridRegionIndex << std::endl;
          if(gridRegionIndex < 0){
            messageStream << "Failed to find grid region id for \""
                          << gridNames[gridIndex] << ":" << regionName 
                          << "\" in BC resolve." << std::endl;
            //            std::cout << "Did not find region index, aborting." << std::endl;
            return(1);
          }
            
          /// @todo Revisit this, probably we shouldn't be "pushing back" here [mtc]
          BoundaryType tempBoundary;
          domainGridBoundaries.push_back(tempBoundary);
          int boundaryIndex = domainGridBoundaries.size() - 1;

          BoundaryType &thisBoundary(domainGridBoundaries[boundaryIndex]);
          thisBoundary.SetGridRegionID(gridRegionIndex);
          thisBoundary.SetGrid(boundaryGrid);
          thisBoundary.SetBC(boundaryBC);
          thisBoundary.SetName(boundaryName);
          thisBoundary.CreateData();

          messageStream << "Found boundary " << boundaryName << " with BC " 
                        << bcName << " for (grid:region) = ("
                        << gridNames[gridIndex] << ":" << regionName
                        << ")" << std::endl;
          //          std::cout << "Adding boundary!" << std::endl;
          //          domainGridBoundaries.push_back(thisBoundary);
        }
      }      
      return(0);
    };

    // If this function returns non-zero, the grid(s) is *not* OK
    int ValidateGrids()
    {
      int returnCode = 0;

      std::ostringstream messageStream;

      std::vector<int> &localGridIndices(LocalGridIndices());
      std::vector<int>::iterator localGridIndexIt = localGridIndices.begin();
      while(localGridIndexIt != localGridIndices.end()){
        
        int localGridId = (localGridIndexIt - localGridIndices.begin()) + 1;
        int gridIndex = *localGridIndexIt++;

        if(gridIndex < 0){
          messageStream << "Domain::ValidateGrids:Error: Invalid index ("
                        << gridIndex << ") for local grid (id=" << localGridId
                        << ")." << std::endl;

          returnCode = 1;
          continue;
        }

        GridType &currentGrid(*domainGrids[gridIndex]);
        int gridDim = currentGrid.Dimension();
        if(gridDim <= 0){
          messageStream << "Domain::ValidateGrids:Error: Invalid dimension (" << gridDim
                        << ") for grid(id=" << gridIndex+1 << ")." << std::endl;
          returnCode = 1;
          continue;
        }

        std::vector<BoundaryType>   &gridDomainBoundaries(DomainBoundary(gridIndex));
        std::vector<GridRegionType> &gridRegions(currentGrid.SubRegions());

        int holeStencil          = domainOperator.numStencils;

        /// @note Need to validate the comm size against partitioning [mtc]
        pcpp::CommunicatorType &gridComm(currentGrid.Communicator());
        
      
        std::vector<size_t> &bufferSizes(currentGrid.BufferSizes());
        if(bufferSizes.size() != gridDim){
          messageStream << "Domain::ValidateGrids:Error: Invalid buffer dimension for grid(id="
                        << gridIndex + 1 << ")." << std::endl;
          returnCode = 1;
          continue;
        }

        size_t numPointsGridBuffer = currentGrid.BufferSize();
        size_t numPointsCheck = 1;
        for(int iDim = 0;iDim < gridDim;iDim++)
          numPointsCheck *= bufferSizes[iDim];
        if(numPointsGridBuffer != numPointsCheck){
          messageStream << "Domain::ValidateGrids:Error: Buffer sizes/numpoints mismatch grid(id="
                        << gridIndex + 1 << ")." << std::endl;
          returnCode = 1;
          continue;
        }

        pcpp::IndexIntervalType &partitionBufferInterval(currentGrid.PartitionBufferInterval());
        pcpp::IndexIntervalType &partitionInterval(currentGrid.PartitionInterval());

        if((partitionInterval.NNodes() == 0) || (partitionBufferInterval.NNodes() == 0)){
          messageStream << "Domain::ValidateGrids:Error: Empty partition of grid(id="
                        << gridIndex + 1 << ")." << std::endl;
          returnCode = 1;
          continue;
        }
      }

      if(returnCode){
        if(globalPtr){
          globalPtr->ErrOut(messageStream.str());
        } else {
          std::cerr << messageStream.str();
        }
      }
      
      return(returnCode);
      
    };

    // Allocates the grid mask, considers any iblank, and flags holes and boundaries
    int CreateGridMasks()
    {
      
      int returnCode = 0;

      returnCode = ValidateGrids();

      if(returnCode){
        domainCommunicator.SetErr(1);
      }
      if(domainCommunicator.Check()){
        if(globalPtr){
          globalPtr->ErrOut("CreateGridMasks:Error: Grid validation failed.\n");
        } else {
          std::cerr << "CreateGridMasks:Error: Grid validation failed." << std::endl;
        }
        return(1);
      }
      
      std::ostringstream messageStream;

      std::vector<int> &localGridIndices(LocalGridIndices());
      std::vector<int>::iterator localGridIndexIt = localGridIndices.begin();
      while(localGridIndexIt != localGridIndices.end()){

        int gridIndex = *localGridIndexIt++;
        GridType &currentGrid(*domainGrids[gridIndex]);
        int gridDim = currentGrid.Dimension();

        std::vector<BoundaryType>   &gridDomainBoundaries(DomainBoundary(gridIndex));
        std::vector<GridRegionType> &gridRegions(currentGrid.SubRegions());

        int holeStencil          = domainOperator.numStencils;

        pcpp::CommunicatorType &gridComm(currentGrid.Communicator());
      
      
        std::vector<size_t> &bufferSizes(currentGrid.BufferSizes());

        size_t numPointsGridBuffer = currentGrid.BufferSize();

        pcpp::IndexIntervalType &partitionBufferInterval(currentGrid.PartitionBufferInterval());
        pcpp::IndexIntervalType &partitionInterval(currentGrid.PartitionInterval());

        pcpp::IndexIntervalType bufferInterval(gridDim,&bufferSizes[0]);
        std::vector<size_t> partitionBufferExtent;
        partitionBufferInterval.Flatten(partitionBufferExtent);

        returnCode    = currentGrid.CreateMask();
        int *gridMask = currentGrid.Mask();

        // Set boundary bits for points that sit in prescribed holes and boundaries
        bool haveHoles = false;
        typename std::vector<BoundaryType>::iterator dbIt = gridDomainBoundaries.begin();
        while(dbIt != gridDomainBoundaries.end()){
          BoundaryType &domainBoundary(*dbIt++);
          GridRegionType &gridRegion(gridRegions[domainBoundary.GridRegionID()]);
          
          int boundaryNormalDirection = gridRegion.normalDirection;
          pcpp::IndexIntervalType &boundaryPartitionInterval(gridRegion.regionPartitionInterval);
          
          int orientationMask = mask::OrientedBoundaryMask(boundaryNormalDirection);

          BCType &boundaryCondition(domainBoundary.BC());
          
          // unique int id for bc
          int bcType = boundaryCondition.BCType();
          
          pcpp::IndexIntervalType boundaryBufferInterval;
          boundaryPartitionInterval.RelativeTranslation(partitionInterval,
                                                        partitionBufferInterval,
                                                        boundaryBufferInterval);
          if(boundaryBufferInterval.NNodes() > 0) {
            int maskValue = 0;
            if(bcType == testboundary::bc::HOLE) { 
              maskValue = mask::holeMask;
            } else {
              maskValue = (mask::domainBoundaryMask|orientationMask);
            }        
            mask::SetMask(bufferSizes,boundaryBufferInterval,maskValue,gridMask);
          }
        }        
        size_t numHoles = currentGrid.FlagMask(gridMask,0,mask::holeMask);
      }


      if(returnCode > 0){
        domainCommunicator.SetErr(1);
      }

      if(domainCommunicator.Check())
        returnCode = 1;
      if(returnCode > 0){
        messageStream << "NSDomain::SetGridMasks:Error: Failed to finalize grid setup." << std::endl;
        if(globalPtr){
          globalPtr->ErrOut(messageStream.str());
        } else {
          std::cerr << messageStream.str(); 
        }
      }
      return(returnCode);      
    }

    // Assumes mask is already created and/or populated. Creates the stencil connectivity
    // and populates it based on the mask values
    int SetupGridOperators()
    {
      int returnCode = 0;

      returnCode = ValidateGrids();

      if(returnCode){
        domainCommunicator.SetErr(1);
      }
      if(domainCommunicator.Check()){
        if(globalPtr){
          globalPtr->ErrOut("SetupGridOperators:Error: Grid validation failed.\n");
        } else {
          std::cerr << "SetupGridOperators:Error: Grid validation failed." << std::endl;
        }
        return(1);
      }
      
      std::ostringstream messageStream;

      std::vector<int> &localGridIndices(LocalGridIndices());
      std::vector<int>::iterator localGridIndexIt = localGridIndices.begin();
      while(localGridIndexIt != localGridIndices.end()){

        int gridIndex = *localGridIndexIt++;
        GridType &currentGrid(*domainGrids[gridIndex]);
        int gridDim = currentGrid.Dimension();
        int holeStencil          = domainOperator.numStencils;
        pcpp::CommunicatorType &gridComm(currentGrid.Communicator());
      
        if(currentGrid.AxiSymmetric()){
          operators::sbp::Initialize(domainOperator,domainOperator.interiorOrder,true);
          // fortranification
          int *stencilStarts = domainOperator.stencilStarts;
          int numStencils    = domainOperator.numStencils;
          for(int iStencil = 0;iStencil < numStencils;iStencil++)
            stencilStarts[iStencil]++;
        }
        currentGrid.SetDifferentialOperator(&domainOperator);
        if(currentGrid.InitializeStencilConnectivity()){
          gridComm.SetErr(1);
        }
        if(gridComm.Check()){
          messageStream << "NSDomain::SetupGridOperators:Error: Failed to initialize stencil connectivity"
                        << " for grid(id=" << gridIndex + 1 << ","
                        << gridNames[gridIndex] << ")." << std::endl;
          returnCode = 1;
          continue;
        }        

        int myReturnCode = currentGrid.DetectHoles();

        if(myReturnCode)
          gridComm.SetErr(1);
        if(gridComm.Check()){
          std::ostringstream failMessage;
          failMessage << "NSDomain::SetupGridOperators::Error: Placing stencils in grid("
                      << gridNames[gridIndex] << ")failed."
                      << std::endl;
          if(globalPtr){
            globalPtr->ErrOut(failMessage.str());
          } else {
            std::cerr << failMessage.str();
          }
          returnCode = 1;
        }
      }

      if(returnCode > 0){
        domainCommunicator.SetErr(1);
      }
      if(domainCommunicator.Check())
        returnCode = 1;
      
      return(returnCode);      
    }

    int ComputeGridMetrics()
    {
      
      int returnCode = 0;

      returnCode = ValidateGrids();

      if(returnCode){
        domainCommunicator.SetErr(1);
      }
      if(domainCommunicator.Check()){
        if(globalPtr){
          globalPtr->ErrOut("ComputeGridMetrics:Error: Grid validation failed.\n");
        } else {
          std::cerr << "ComputeGridMetrics:Error: Grid validation failed." << std::endl;
        }
        return(1);
      }
      

      std::ostringstream messageStream;

      std::vector<int> &localGridIndices(LocalGridIndices());
      std::vector<int>::iterator localGridIndexIt = localGridIndices.begin();
      while(localGridIndexIt != localGridIndices.end()){

        int gridIndex = *localGridIndexIt++;
        GridType &currentGrid(*domainGrids[gridIndex]);
        int gridDim = currentGrid.Dimension();
        pcpp::CommunicatorType &gridComm(currentGrid.Communicator());
        
        std::ostringstream metricMessages;
        int errorCode = currentGrid.ComputeMetrics(metricMessages);
        if(errorCode)
          gridComm.SetErr(1);
        if(gridComm.Check()){
          
          std::ostringstream failMessage;
          failMessage << "NSDomain::ComputeGridMetrics::Error: grid::ComputeMetrics failed "
                      << "for grid(id=" << gridIndex+1 << ") with messages:" 
                      << std::endl << metricMessages.str() << std::endl;
          if(globalPtr){
            globalPtr->ErrOut(failMessage.str());
          } else {
            std::cout << failMessage.str();
          }
        } else {
          if(globalPtr){
            globalPtr->StdOut("Grid metric messages:\n",3);
            globalPtr->StdOut(metricMessages.str(),3);
          } else {
            std::cout << "Grid metric messages:\n" << metricMessages.str() << std::endl;
          }
        }
        
        pcpp::io::RenewStream(metricMessages);
        int gridType = currentGrid.Type();

        if(gridType == simulation::grid::CURVILINEAR){
          int metricErrorCode = plascom2::util::CheckMetricIdentities(currentGrid,metricMessages);
          if(metricErrorCode){
            std::ostringstream failMessage;
            failMessage << "NSDomain::ComputeGridMetrics:WARNING: GCL test failed for grid(id="
                        << gridIndex+1 << ")" << " with messages:" << std::endl
                        << metricMessages.str() << std::endl;
            if(globalPtr){
              globalPtr->StdOut(failMessage.str());
            } else {
              std::cout << failMessage.str();
            }
          } else {
            if(globalPtr){
              globalPtr->StdOut(metricMessages.str(),2);
            } else {
              std::cout << metricMessages.str();
            }
          }
        }
        
      }
      if(returnCode > 0){
        domainCommunicator.SetErr(1);
      }
      if(domainCommunicator.Check())
        returnCode = 1;
    
      return(returnCode);      
    }
    
#ifdef USE_OVERKIT
    int InitializeOverkit(ovk_shared_context *inOverkitContextPtr, std::ostream &messageStream)
    {

      overkitContextPtr = inOverkitContextPtr;

      std::vector<int> &gridIndexOfLocalGrids(partitionInfo.gridIndexOfLocalGrids);

      // Overkit currently requires all grids to have the same number of dimensions
      int numDim = domainGrids[gridIndexOfLocalGrids[0]]->Dimension();

      int numLocalGrids = gridIndexOfLocalGrids.size();

      // Create the domain
      ovk_domain_params *domainParams;
      ovkCreateDomainParams(&domainParams);
      ovkSetDomainParamName(domainParams, domainName.c_str());
      ovkSetDomainParamDimension(domainParams, numDim);
      ovkSetDomainParamComm(domainParams, domainCommunicator.GetCommunicator());
      ovkCreateDomain(&overkitDomainPtr, overkitContextPtr, &domainParams);

      std::vector<int> gridIDs(numGrids);
      std::vector<ovk_grid_params *> maybeGridParams(numGrids, nullptr);

      // Fill the grid params
      for (int iGrid = 0; iGrid < numGrids; ++iGrid) {

        bool gridIsLocal(GridIsLocal(iGrid));

        gridIDs[iGrid] = iGrid+1;

        if (gridIsLocal) {

          const GridType &domainGrid(*domainGrids[iGrid]);

          const pcpp::CommunicatorType &gridComm(domainGrid.Communicator());
          const pcpp::IndexIntervalType &partitionInterval(domainGrid.PartitionInterval());
          const std::vector<size_t> &gridSizes(domainGrid.GridSizes());
          const std::vector<bool> &isPeriodic(domainGrid.PeriodicDirs());
          const std::vector<double> &gridPeriodicLengths(domainGrid.PeriodicLengths());
          int gridType = domainGrid.Type();

          // Global index range
          int globalBegin[3] = {0, 0, 0};
          int globalEnd[3] = {1, 1, 1};
          for (int d = 0; d < numDim; ++d) {
            globalEnd[d] = (int)gridSizes[d];
          }

          // Local index range
          int localBegin[3] = {0, 0, 0};
          int localEnd[3] = {1, 1, 1};
          for (int d = 0; d < numDim; ++d) {
            localBegin[d] = partitionInterval[d].first;
            localEnd[d] = partitionInterval[d].second + 1;
          }

          // Can't pass as &isPeriodic[0] because std::vector<bool> is special
          bool periodic[3] = {false, false, false};
          if (!isPeriodic.empty()) {
            for (int d = 0; d < numDim; ++d) {
              periodic[d] = isPeriodic[d];
            }
          }

          // Create grid
          ovkCreateGridParams(&maybeGridParams[iGrid]);
          ovk_grid_params *gridParams = maybeGridParams[iGrid];
          ovkSetGridParamName(gridParams, gridNames[iGrid].c_str());
          ovkSetGridParamDimension(gridParams, numDim);
          ovkSetGridParamComm(gridParams, gridComm.GetCommunicator());
          ovkSetGridParamGlobalRange(gridParams, globalBegin, globalEnd);
          ovkSetGridParamLocalRange(gridParams, localBegin, localEnd);
          ovkSetGridParamPeriodic(gridParams, periodic);

        }

      }

      // Create the grids
      ovkCreateGrids(overkitDomainPtr, numGrids, &gridIDs[0], &maybeGridParams[0]);

      // Create connectivity component (just use ID=1 for now)
      ovkCreateComponent(overkitDomainPtr, 1, OVK_COMPONENT_TYPE_CONNECTIVITY, NULL);

      return(0);
    }
#endif

    int ImportXINTOUT(const std::string &HOFileName, const std::string &XFileName,
      int readGranularityAdjust, std::ostream &messageStream) {

#ifdef USE_OVERKIT
      ovk_error error;
      ovkImportXINTOUT(overkitDomainPtr, 1, HOFileName.c_str(), XFileName.c_str(),
        readGranularityAdjust, MPI_INFO_NULL, &error);
      if (error != OVK_ERROR_NONE) {
        messageStream << "Reading XINTOUT failed." << std::endl;
        return(1);
      }
#endif

      return(0);

    }

    int CreateOverkitInterGrid(std::ostream &messageStream) {

#ifdef USE_OVERKIT
      interGridPtr = new intergrid::overkit_intergrid<GridType>(overkitDomainPtr, 1);
#endif

      int error = interGridPtr->Initialize(domainCommunicator, domainGrids,
        partitionInfo.gridIndexOfLocalGrids, messageStream);
      if (error) {
        messageStream << "Overkit inter-grid initialization failed." << std::endl;
        return(1);
      }

      return(0);

    }

    int ConfigureInterGrid(std::ostream &messageStream) {

      if (domainConfig.IsSet(ConfigKey(domainName, "InterGrid"))) {

        std::string interGridMode = domainConfig.GetValue(ConfigKey(domainName, "InterGrid"));

        if (interGridMode != "None") {

          if (interGridMode == "Overkit") {

#ifndef USE_OVERKIT
            messageStream << "Must build with Overkit to use Overkit inter-grid mode."
              << std::endl;
            return(1);
#endif

            std::string connectivitySource = domainConfig.GetValue(ConfigKey(domainName,
              "Overkit:ConnectivitySource"));

            if (connectivitySource == "Assemble") {

              messageStream << "Test domain doesn't support runtime overset assembly." << std::endl;
              return(1);

            } else if (connectivitySource == "Import") {

              fixtures::ConfigurationType importConfig = pcpp::util::GetSubConfig(
                ConfigKey(domainName, "Overkit:Import"), domainConfig);

              fixtures::ConfigurationType readConfig = pcpp::util::GetSubConfig(
                ConfigKey(domainName, "Overkit:Read"), domainConfig);

              std::string readFormat = readConfig.GetValue("Format");

              if (readFormat == "XINTOUT") {
                std::string HOFileName = readConfig.GetValue("XINTOUT:HOFileName");
                std::string XFileName = readConfig.GetValue("XINTOUT:XFileName");
                int readGranularityAdjust = 0;
                if (readConfig.IsSet("XINTOUT:ReadGranularityAdjust")) {
                  readGranularityAdjust = readConfig.GetValue<int>("XINTOUT:ReadGranularityAdjust");
                }
                int error = ImportXINTOUT(HOFileName, XFileName, readGranularityAdjust,
                  messageStream);
                if (error) {
                  messageStream << "Overkit import failed." << std::endl;
                  return(1);
                }
              } else {
                messageStream << "Unknown connectivity file format '" << readFormat << "'."
                  << std::endl;
                return(1);
              }

            } else {

              messageStream << "Unknown connectivity source '" << connectivitySource << "'."
                << std::endl;
              return(1);

            }

            int error = CreateOverkitInterGrid(messageStream);
            if (error) {
              messageStream << "Failed to created Overkit inter-grid object." << std::endl;
              return(1);
            }

          } else {

            messageStream << "Unknown inter-grid mode '" << interGridMode << "'." << std::endl;
            return(1);

          }

        }

      } else if (domainConfig.IsSet(ConfigKey(domainName, "Connectivity:Source"))) {

        // Don't silently fail when old connectivity options are specified
        std::string connectivitySource = domainConfig.GetValue(ConfigKey(domainName,
          "Connectivity:Source"));
        if (connectivitySource != "None") {
          messageStream << "Configuration uses old inter-grid connectivity options, which are "
            << "no longer supported." << std::endl;
          return(1);
        }

      }

      return(0);
    }

    int InitializeEfield(const std::vector<bool> &gridHasEfield, std::ostream &messageStream) {

      efieldSolverPtr = new EfieldSolverType();

      int error = efieldSolverPtr->Initialize(domainCommunicator, domainGrids, gridHasEfield,
        partitionInfo.gridIndexOfLocalGrids, messageStream);
      if (error) {
        messageStream << "Failed to initialize Efield solver." << std::endl;
        return(1);
      }

      return(0);

    }

    bool HasEfield(int iGrid = -1) {

      if (efieldSolverPtr) {
        if (iGrid >= 0) {
          return efieldSolverPtr->HasEfield(iGrid);
        } else {
          return true;
        }
      } else {
        return false;
      }

    }

    EfieldType &Efield(int iGrid) { return efieldSolverPtr->Efield(iGrid); }

    int SetEfieldState(int iGrid, std::ostream &messageStream) {

      int error = efieldSolverPtr->SetState(iGrid, *gridStates[iGrid]);
      if (error) {
        messageStream << "Failed to set Efield state." << std::endl;
        return(1);
      }

      return(0);

    }

    int SetupEfieldSolver(std::ostream &messageStream) {

      int error;

      error = efieldSolverPtr->SetDomainBCs(domainBoundaries);
      if (error) {
        messageStream << "Failed to set Efield solver boundary conditions." << std::endl;
        return(1);
      }

      if (interGridPtr) {
        error = efieldSolverPtr->SetInterGrid(*interGridPtr);
        if (error) {
          messageStream << "Failed to set Efield solver inter-grid communication." << std::endl;
          return(1);
        }
      }

      error = efieldSolverPtr->SetupSolver(messageStream);
      if (error) {
        messageStream << "Failed to set up Efield solver." << std::endl;
        return(1);
      }

      return(0);

    }

    void SetEfieldNumSchwarz(int numSchwarz) {
      efieldSolverPtr->SetNumSchwarz(numSchwarz);
    }

    int EfieldSolveInterval() { return efieldSolveInterval; }

    void SetEfieldSolveInterval(int solveInterval) {
      efieldSolveInterval = solveInterval;
    }

    int SolveEfield(int threadId) {

      int error = efieldSolverPtr->Solve(threadId);
      if (error) {
        return(1);
      }

      return(0);

    }

    int ConfigureEfield(std::ostream &messageStream) {

      fixtures::ConfigurationType efieldConfig = pcpp::util::GetSubConfig(
        ConfigKey(domainName, "Efield"), domainConfig);

      if (efieldConfig.IsSet("GridNames")) {

        std::vector<std::string> efieldGridNames = efieldConfig.GetValueVector("GridNames");

        if (efieldGridNames.size() > 0) {

          std::vector<bool> gridHasEfield(numGrids, false);

          for (int iGridName = 0; iGridName < efieldGridNames.size(); ++iGridName) {
            std::vector<std::string>::iterator Iter = std::find(gridNames.begin(),
              gridNames.end(), efieldGridNames[iGridName]);
            if (Iter != gridNames.end()) {
              int iGrid = std::distance(gridNames.begin(), Iter);
              gridHasEfield[iGrid] = true;
            } else {
              messageStream << "Invalid grid name " << efieldGridNames[iGridName] << "." << std::endl;
              return(1);
            }
          }

          int error = InitializeEfield(gridHasEfield, messageStream);
          if (error) {
            messageStream << "Efield solver initialization failed." << std::endl;
            return(1);
          }

          int numSchwarz = 0;
          if (efieldConfig.IsSet("NumSchwarz")) {
            numSchwarz = efieldConfig.GetValue<int>("NumSchwarz");
          }
          SetEfieldNumSchwarz(numSchwarz);

          int solveInterval = 0;
          if (efieldConfig.IsSet("SolveInterval")) {
            solveInterval = efieldConfig.GetValue<int>("SolveInterval");
          }
          SetEfieldSolveInterval(solveInterval);

        }

      }

      return(0);

    }

    int CreateData(std::ostream &messageStream)
    {

      int numGrids      = domainGrids.size();
      int numGridStates = gridStates.size();
      int numGridParams = gridParams.size();

      if(numGrids != numGridStates)
        return(1);
      if(numGridStates != numGridParams)
        return(1);

      std::string stateFields;
      if(domainConfig.IsSet("Advancer:Fields")){
        stateFields = domainConfig.GetValue(ConfigKey(domainName,"Advancer:Fields"));
        messageStream << "Setting state fields to (" << stateFields << ")" << std::endl;
      } else {
        messageStream << "WARNING: No fields set up as State Fields." << std::endl;
      }
      typename gridvector::iterator gridIt   = domainGrids.begin();
      typename statevector::iterator stateIt = gridStates.begin();
      typename statevector::iterator targIt  = gridTargets.begin();
      typename statevector::iterator paramIt = gridParams.begin();
      std::vector<std::string>::iterator gridNameIt = gridNames.begin();

      messageStream << "Creating state and parameter data for " << domainGrids.size() << " grids."
                    << std::endl;

      while(gridIt != domainGrids.end()){

        GridType  &thisGrid(**gridIt++);
        StateType &gridState(**stateIt++);
        //        StateType &gridTarget(**targIt++);
        StateType &gridParams(**paramIt++);
        std::string &gridName(*gridNameIt++);

        const std::vector<size_t> &bufferSizes(thisGrid.BufferSizes());
        const pcpp::IndexIntervalType &partitionInterval(thisGrid.PartitionInterval());
        const pcpp::IndexIntervalType &partitionBufferInterval(thisGrid.PartitionBufferInterval());  
        size_t numPointsBuffer = thisGrid.BufferSize();


        messageStream << "Creating buffers for " << gridName
                      << " with bufferSizes(";
        pcpp::io::DumpContents(messageStream,bufferSizes,",");
        messageStream << ")" << std::endl;

        gridState.SetMetaData(stateDataDictionary);
        gridParams.SetMetaData(paramDataDictionary);

        gridState.Create(numPointsBuffer,0);
        if(!stateFields.empty()) 
          gridState.SetStateFields(stateFields);
        if(*targIt){
          (*targIt)->CopyStateFields(gridState);
        }
        gridParams.Create(numPointsBuffer,0);

        std::string paramKeyRoot(ConfigKey(domainName,"Param"));
        messageStream << "Configuring parameters with KeyRoot(" 
                      << paramKeyRoot << ")." << std::endl;
        if(ConfigStateValues(domainConfig,paramKeyRoot,gridParams,messageStream)){
          messageStream << "ERROR! Domain parameters failed to populate." << std::endl;
          return(1);
        }

        int *optionsFlagPtr = gridParams.template GetFieldBuffer<int>("optionFlags");
        if(optionsFlagPtr == NULL){
          messageStream << "ERROR! Domain parameters missing internal data item for optionFlags!"
                        << std::endl;
        }
        int legacyOptions  = GetParameterValue<int>("Flags",gridParams,0);
        *optionsFlagPtr = legacyOptions;
        std::string optionsKey(ConfigKey(domainName,"Option"));
        bool optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"Pooya"));
        if(optionVal) *optionsFlagPtr |= navierstokes::USEPOOYA;
        optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"ExchangeFlux"));
        if(optionVal) *optionsFlagPtr |= navierstokes::EXCHANGEFLUX;
        optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"ExchangeDV"));
        if(optionVal) *optionsFlagPtr |= navierstokes::EXCHANGEDV;
        optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"WENO"));
        if(optionVal) *optionsFlagPtr |= navierstokes::USEWENO;
        optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"Axisymmetry"));
        if(optionVal) *optionsFlagPtr |= navierstokes::AXISYMM;
        optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"MetricIdentityCorrection"));
        if(optionVal) *optionsFlagPtr |= navierstokes::METRICID;
        
      }
      return(0);
    };

    void RemoveLegacyDictionary()
    {
      stateDataDictionary.RemoveMatchingFields("cv0");
      stateDataDictionary.RemoveMatchingFields("aux");
    }

    int CreateComputeData(int gridIndex,std::ostream &messageStream)
    {
      
      if(runMode != pcpp::runmode::SPMD){
        if(virtualNodePtr->performCompute < 1){
          messageStream << "Warning: CreateComputeData called by non-compute rank ("
                        << domainRank << " for grid index = " << gridIndex  
                        << std::endl;
          return(0);
        }
      }
      
      std::string stateFields(domainConfig.GetValue(ConfigKey(domainName,"Advancer:Fields")));
      messageStream << "Setting state fields to (" << stateFields << ")" << std::endl;      
      messageStream << "Creating state and parameter data for grid(" << gridIndex << ","
                    << gridNames[gridIndex] << ")" << std::endl;
      
      
      GridType  &thisGrid(Grid(gridIndex));
      StateType &gridState(State(gridIndex));
      //      StateType &gridTarget(Target(gridIndex));
      StateType &gridParams(Param(gridIndex));
      std::string &gridName(gridNames[gridIndex]);
      
      const std::vector<size_t> &bufferSizes(thisGrid.BufferSizes());
      const pcpp::IndexIntervalType &partitionInterval(thisGrid.PartitionInterval());
      const pcpp::IndexIntervalType &partitionBufferInterval(thisGrid.PartitionBufferInterval());  
      size_t numPointsBuffer = thisGrid.BufferSize();
      
      
      messageStream << "Creating buffers for " << gridName
                    << " with bufferSizes(";
      pcpp::io::DumpContents(messageStream,bufferSizes,",");
      messageStream << ")" << std::endl;
      
      gridState.SetMetaData(stateDataDictionary);
      gridParams.SetMetaData(paramDataDictionary);
      
      gridState.Create(numPointsBuffer,0);
      gridState.SetStateFields(stateFields);
      if(gridTargets[gridIndex] == NULL){
        gridTargets[gridIndex] = new StateType;
      }
      gridTargets[gridIndex]->CopyStateData(gridState);
      
      gridParams.Create(numPointsBuffer,0);
      
      std::string paramKeyRoot(ConfigKey(domainName,"Param"));
      if(ConfigStateValues(domainConfig,paramKeyRoot,gridParams,messageStream)){
        messageStream << "ERROR! Domain parameters failed to populate." << std::endl;
        return(1);
      }

      int *optionsFlagPtr = gridParams.template GetFieldBuffer<int>("optionFlags");
      if(optionsFlagPtr == NULL){
        messageStream << "ERROR! Domain parameters missing internal data item for optionFlags!"
                      << std::endl;
      }
      int legacyOptions  = GetParameterValue<int>("Flags",gridParams,0);
      *optionsFlagPtr = legacyOptions;
      std::string optionsKey(ConfigKey(domainName,"Option"));
      bool optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"Pooya"));
      if(optionVal) *optionsFlagPtr |= navierstokes::USEPOOYA;
      optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"ExchangeFlux"));
      if(optionVal) *optionsFlagPtr |= navierstokes::EXCHANGEFLUX;
      optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"ExchangeDV"));
      if(optionVal) *optionsFlagPtr |= navierstokes::EXCHANGEDV;
      optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"WENO"));
      if(optionVal) *optionsFlagPtr |= navierstokes::USEWENO;
      optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"Axisymmetry"));
      if(optionVal) *optionsFlagPtr |= navierstokes::AXISYMM;
      optionVal = domainConfig.GetFlagValue(ConfigKey(optionsKey,"MetricIdentityCorrection"));
      if(optionVal) *optionsFlagPtr |= navierstokes::METRICID;
      
      return(0);
    };
    
    int CreateIOData(int gridIndex,std::ostream &messageStream)
    {
      
      if(runMode == pcpp::runmode::SPMD || runMode == pcpp::runmode::SPAMD)
        return(0);

      /// @note This check will prevent mixed-duty ranks outside of SPAMD [mtc]
      if((virtualNodePtr->performCompute > 0) &&
         (virtualNodePtr->performIO      > 0))
        return(0);

      std::string stateFields(domainConfig.GetValue(ConfigKey(domainName,"Advancer:Fields")));
      messageStream << "Setting state fields to (" << stateFields << ")" << std::endl;      
      messageStream << "Creating state and parameter data for grid(" << gridIndex << ","
                    << gridNames[gridIndex] << ")" << std::endl;
      
      
      GridType  &thisGrid(Grid(gridIndex));
      StateType &gridState(State(gridIndex));
      //      StateType &gridTarget(Target(gridIndex));
      StateType &gridParams(Param(gridIndex));
      std::string &gridName(gridNames[gridIndex]);
      
      const std::vector<size_t> &bufferSizes(thisGrid.BufferSizes());
      const pcpp::IndexIntervalType &partitionInterval(thisGrid.PartitionInterval());
      const pcpp::IndexIntervalType &partitionBufferInterval(thisGrid.PartitionBufferInterval());  
      size_t numPointsBuffer = thisGrid.BufferSize();
      
      
      messageStream << "Creating buffers for " << gridName
                    << " with bufferSizes(";
      pcpp::io::DumpContents(messageStream,bufferSizes,",");
      messageStream << ")" << std::endl;
      
      gridState.SetMetaData(stateDataDictionary);
      gridParams.SetMetaData(paramDataDictionary);
      
      gridState.Create(numPointsBuffer,0);
      gridState.SetStateFields(stateFields);
      if(gridTargets[gridIndex] == NULL){
        gridTargets[gridIndex] = new StateType;
      }
      gridTargets[gridIndex]->CopyStateData(gridState);
      
      gridParams.Create(numPointsBuffer,0);
      
      std::string paramKeyRoot(ConfigKey(domainName,"Param"));
      if(ConfigStateValues(domainConfig,paramKeyRoot,gridParams,messageStream)){
        messageStream << "ERROR! Domain parameters failed to populate." << std::endl;
        return(1);
      }
      
      return(0);
    };

    /// CreateData configures and allocates main storage for state and parameter data
    int CreateData(int gridIndex,std::ostream &messageStream)
    {
      if(gridIndex > (domainGrids.size()-1))
        return(1);

      if(runMode == pcpp::runmode::SPMD){

        if(CreateComputeData(gridIndex,messageStream)){
          messageStream << "Domain::CreateData:Error: CreateComputeData failed." << std::endl;
          return(1);
        }

      } else {

        if(!virtualNodePtr){
          std::cout << "Domain::CreateData:Error: missing virtual node information for MPMD mode."
                    << std::endl;
          return(1);
        }

        if(virtualNodePtr->performIO > 0){
          if(CreateIOData(gridIndex,messageStream)){
            messageStream << "Domain::CreateData:Error: CreateIOData failed." << std::endl;
            return(1);
          }
        }

        if(virtualNodePtr->performCompute > 0){
          if(CreateComputeData(gridIndex,messageStream)){
            messageStream << "Domain::CreateData:Error: CreateComputeData failed." << std::endl;
            return(1);
          }
        }
        
      }
      
      return(0);

    };

    int GetGridIndex(const std::string &gridName) const
    {
      int gridIndex = 0;
      std::vector<std::string>::const_iterator nameIt = gridNames.begin();
      while(nameIt != gridNames.end()){
        if(*nameIt++ == gridName){
          return(gridIndex);
        }
        gridIndex++;
      }
      return(-1);
    };


    int NumberOfGrids(){ return(numGrids); };
    int NumberOfLocalGrids(){ return(partitionInfo.gridIndexOfLocalGrids.size()); };
    std::vector<int> &LocalGridIndices()
    { 
      return(partitionInfo.gridIndexOfLocalGrids); 
    };
    int LocalGridIndex(int localGridIndex)
    { 
      return(partitionInfo.gridIndexOfLocalGrids[localGridIndex]); 
    };
    bool GridIsLocal(int globalGridIndex) const
    {
      const std::vector<int> &globalGridIndices(partitionInfo.gridIndexOfLocalGrids);
      std::vector<int>::const_iterator indexIt = std::find(globalGridIndices.begin(),
        globalGridIndices.end(), globalGridIndex);
      return(indexIt != globalGridIndices.end());
    };
    int GetGridRootRank(int globalGridIndex)
    {
      int rootRank;
      if (GridIsLocal(globalGridIndex)) {
        GridType &domainGrid(*domainGrids[globalGridIndex]);
        if (domainGrid.Communicator().Rank() == 0) {
          rootRank = domainCommunicator.Rank();
        }
      }
      MPI_Bcast(&rootRank, 1, MPI_INT, 0, domainCommunicator.GetCommunicator());
      return(rootRank);
    };

    std::string Report()
    {
      std::ostringstream Ostr;
      Ostr << "+++Domain data information+++ " << std::endl
           << "State data dictionary:        " << std::endl
           << stateDataDictionary.Report()     << std::endl
           << "Parameter data dictionary:    " << std::endl
           << paramDataDictionary.Report()     << std::endl;
      if(!domainBCs.empty()){
        Ostr << "Boundary Conditions:          " << std::endl;
        typename std::vector<BCType>::iterator bcIt = domainBCs.begin();
        int bcIndex = 0;
        while(bcIt != domainBCs.end()){
          BCType &myBC(*bcIt++);
          Ostr << "BC(" << bcNames[bcIndex++] << ") = " 
               << myBC.BCName() << std::endl;
          if(!myBC.StateDictionary().empty()){
            Ostr << "  Boundary state data: " << std::endl
                 << myBC.StateDictionary().Report() << std::endl;
          } else {
            Ostr << "  No boundary state data. " << std::endl;
          }
          if(!myBC.ParamDictionary().empty()){
            Ostr << "  Boundary params: " << std::endl
                 << myBC.ParamDictionary().Report() << std::endl;
          } else {
            Ostr << "  No boundary parameters. " << std::endl;
          }
        }
      } else {
        Ostr << "No boundary conditions configured."      << std::endl;
      }
      Ostr << "Number of grids: " << numGrids  << std::endl
           << "... Grids ..." << std::endl;
      for(int iGrid = 0;iGrid < numGrids;iGrid++){
        Ostr << "Grid(" << iGrid+1 << "):" << std::endl
             << "Grid Name: " << gridNames[iGrid] << std::endl
             << "Grid Geometry: " << geometryNames[iGrid] << std::endl;
        if(domainGrids[iGrid] != NULL){
          Ostr << "Grid Configuration: " << std::endl
               << domainGrids[iGrid]->ReportConfiguration() << std::endl;
          std::vector<size_t> &gridSizes(domainGrids[iGrid]->GridSizes());
          std::vector<size_t> &bufferSizes(domainGrids[iGrid]->BufferSizes());
          pcpp::IndexIntervalType &partitionInterval(domainGrids[iGrid]->PartitionInterval());
          pcpp::IndexIntervalType &partitionBufferInterval(domainGrids[iGrid]->PartitionBufferInterval());
          size_t numPointsBuffer = domainGrids[iGrid]->BufferSize();
          Ostr << "Grid Sizes: (";
          pcpp::io::DumpContents(Ostr,gridSizes,",");
          Ostr << ")" << std::endl
               << "Partition Interval: ";
          partitionInterval.PrettyPrint(Ostr);
          Ostr << std::endl
               << "Buffer Sizes: (";
          pcpp::io::DumpContents(Ostr,bufferSizes,",");
          Ostr << ")" << std::endl
               << "Partition Buffer Interval: ";
          partitionBufferInterval.PrettyPrint(Ostr);
          Ostr << std::endl
               << "Buffer size: " << numPointsBuffer << std::endl;
        } 
      }
      return(Ostr.str());
    };
    DataDictionaryType &StateDictionary() { return(stateDataDictionary); };
    DataDictionaryType &ParamDictionary() { return(paramDataDictionary); };

    void PopulateGridMeta(DataDictionaryType &gridMeta)
    {     
      typename DataDictionaryType::iterator metaIt = stateDataDictionary.begin();
      while(metaIt != stateDataDictionary.end()){
        typename StateType::MetaDataType &metaData(*metaIt++);
        if(metaData.loc == 'g' || metaData.loc == 'm' ||
           metaData.loc == 'c' || metaData.loc == 'n')
          gridMeta.AddField(metaData.name,metaData);
      }
    };

    void PopulateGridParamMeta(DataDictionaryType &gridMeta)
    {     
      typename DataDictionaryType::iterator metaIt = paramDataDictionary.begin();
      while(metaIt != paramDataDictionary.end()){
        typename StateType::MetaDataType &metaData(*metaIt++);
        if(metaData.loc == 'g' || metaData.loc == 'm' ||
           metaData.loc == 'c' || metaData.loc == 'n')
          gridMeta.AddField(metaData.name,metaData);
      }
    };

    void PopulateDomainStateMeta()
    {     
      typename DataDictionaryType::iterator metaIt = stateDataDictionary.begin();
      while(metaIt != stateDataDictionary.end()){
        typename StateType::MetaDataType &metaData(*metaIt++);
        if(metaData.loc == 'g' || metaData.loc == 'd' ||
           metaData.loc == 's')
          domainState.Meta().AddField(metaData.name,metaData);
        //          gridMeta.AddField(metaData.name,metaData);
      }
    };
    
    void PopulateDomainParamMeta()
    {     
      typename DataDictionaryType::iterator metaIt = paramDataDictionary.begin();
      while(metaIt != paramDataDictionary.end()){
        typename StateType::MetaDataType &metaData(*metaIt++);
        if(metaData.loc == 'g' || metaData.loc == 'd' ||
           metaData.loc == 's')
          domainParams.Meta().AddField(metaData.name,metaData);
        //          gridMeta.AddField(metaData.name,metaData);
      }
    };

    fixtures::ConfigurationType &Config(){return(domainConfig);};
    std::vector<std::string> &GeometryNames() {return(geometryNames);};
    std::vector<std::string> &GridNames() {return(gridNames); };
    std::string GridName(int iGrid) {return(gridNames[iGrid]); };
    std::vector<std::string> IOFieldNames() { return(ioFieldNames); };
    std::vector<std::string> RestartFieldNames() { return(restartFieldNames); };
    void SetSpatialOrder(int inOrder){ spatialOrder = inOrder;};
    int SpatialOrder(){return(spatialOrder);};
    void SetTime(double inTime){ domainTime = inTime; };
    double Time(){return(domainTime); };
    OperatorType &Operator(){return(domainOperator); };
    const OperatorType &Operator() const {return(domainOperator); };
    void SetMessageStream(std::ostream &inStream){messageStreamPtr = &inStream;};
    std::vector<BCType> &BCs() { return(domainBCs); };
    const std::vector<BCType> &BCs() const { return(domainBCs); };
    RHSType &RHS() { return(*rhsPtr); };
    int SetRHS(RHSType &inRHS) 
    { 
      if(numGrids <= 0){
        return(1);
      }
      if(SetRHS(0,inRHS))
        return(1);
      rhsPtr = gridRHSs[0]; 
      return(0);
    };
    int SetRHS(int gridIndex,RHSType &inRHS) 
    {
      /// @note potential bug here!  revisit in the context of decomp maps [mtc]
      if(partitionInfo.gridIndexOfLocalGrids.empty()){
        if(PartitionDomain(0)){
          return(1);
        }
      }

      if(gridRHSs.size() != numGrids)
        gridRHSs.resize(numGrids,NULL);
      if(ownRHS.size() != numGrids)
        ownRHS.resize(numGrids,false);
      if(ownRHS[gridIndex] && (gridRHSs[gridIndex] != NULL)){
        delete gridRHSs[gridIndex];
        ownRHS[gridIndex] = false;
      }

      gridRHSs[gridIndex] = &inRHS;
      return(0);
    };
    fixtures::CommunicatorType &Communicator() { return(domainCommunicator); };
    fixtures::CommunicatorType &GridCommunicator(int gridIndex) 
    {
      return(gridCommunicators[gridIndex]); 
    };
    fixtures::CommunicatorType &GridComputeCommunicator(int gridIndex) 
    {
      return(gridCommunicators[gridIndex]); 
    };
    fixtures::CommunicatorType &GridWorldCommunicator(int gridIndex) 
    {
      return(gridWorldCommunicators[gridIndex]); 
    };
    fixtures::CommunicatorType &GridIOCommunicator(int gridIndex) 
    {
      return(gridIOCommunicators[gridIndex]); 
    };
    std::vector<fixtures::CommunicatorType> &GridCommunicators() {return(gridCommunicators); };
    fixtures::CommunicatorType &SetCommunicator(fixtures::CommunicatorType &inCommunicator) 
    { return(domainCommunicator = inCommunicator); };

    // This routine only determines which domain processors will calculate on which 
    // grids.  It *does not* decompose the grids.
    int PartitionDomain(int inOptions = 0)
    { 
      int returnCode = 0;
      if(inOptions == 0){ // Trivial "partitioning" where all grids are on all procs
        int numLocalGrids = domainGrids.size();
        partitionInfo.gridIndexOfLocalGrids.resize(numLocalGrids);
        for(int iGrid = 0;iGrid < numLocalGrids;iGrid++){
          partitionInfo.gridIndexOfLocalGrids[iGrid] = iGrid;
	  //	  gridCommunicators[iGrid] = domainCommunicator;
        }
	//	gridCommunicators.resize(numLocalGrids,domainCommunicator);
      } else if(inOptions == 1){ // Auto-split grids according to size
        /// @todo finish automatic domain splitting [mtc]
        int numDomainProcs = 1;
        if(domainCommunicator.Good()){
          numDomainProcs = domainCommunicator.NProc();
          domainRank = domainCommunicator.Rank();
        }
        if(numDomainProcs == 1 || numGrids > numDomainProcs)
          return(PartitionDomain(0));
        size_t numPointsTotal = 0;
        std::vector<size_t> numPointsGrid(numGrids,0);
        std::vector<int>  numProcsGrid(numGrids,0);
        std::vector<size_t> numPointsPerProc(numGrids,0);
        for(int iGrid = 0;iGrid < numGrids;iGrid++){
          GridType &currentGrid(*domainGrids[iGrid]);
          const std::vector<size_t> &gridSizes(currentGrid.GridSizes());
          int numDim = gridSizes.size();
          for(int iDim = 0;iDim < numDim;iDim++){
            numPointsGrid[iGrid] += gridSizes[iDim];
          }
          numPointsTotal += numPointsGrid[iGrid];
        }
        size_t targetPointsPerProc = numPointsTotal/numDomainProcs;
        int totalGridProcs = 0;
        for(int iGrid = 0;iGrid < numGrids;iGrid++){
          numProcsGrid[iGrid] = int(numPointsGrid[iGrid]*numDomainProcs/float(numPointsTotal));
          numProcsGrid[iGrid] = std::max(1,numProcsGrid[iGrid]);
          totalGridProcs += numProcsGrid[iGrid];
          numPointsPerProc[iGrid] = numPointsGrid[iGrid]/numProcsGrid[iGrid];
        }
        if(totalGridProcs > numDomainProcs){
          int numToScrape = totalGridProcs - numDomainProcs;
          std::vector<unsigned int> numProcsOrder;
          ix::util::SortPermutation(numProcsGrid,numProcsOrder);
          std::reverse(numProcsOrder.begin(),numProcsOrder.end());
          int orderIndex = 0;
          for(int iScrape = numToScrape;iScrape > 0;iScrape--){
            int gridIndexToScrape = numProcsOrder[orderIndex];
            numProcsGrid[gridIndexToScrape]--;
            numPointsPerProc[gridIndexToScrape] = 
              numPointsGrid[gridIndexToScrape]/numProcsGrid[gridIndexToScrape];
          }
        } else {
          int numToDoleOut = numDomainProcs - totalGridProcs;
        } 
        //        if(totalGridProcs == numDomainProcs)
          //          return(PartitionDomain(numProcsGrid));

      } else if(inOptions == 3){  // Fully user prescribed domain splitting

        if(!virtualNodePtr){
          std::cout << "Domain::PartitionDomain:Fatal error: missing virtual node information."
                    << std::endl;
          return(1);
        }

        pcpp::VirtualNodeType &virtualNode(*virtualNodePtr);

        int numVirtualNodes   = virtualNode.numNodes;
        int numVirtualCompute = virtualNode.computePerNode;
        int numVirtualIO      = virtualNode.ioPerNode;
        int numDomainCompute  = numVirtualNodes * numVirtualCompute;
        int numDomainIO       = numVirtualNodes * numVirtualIO;
        int virtualNodeId     = virtualNode.nodeId;
        int performIO         = virtualNode.performIO;
        int performCompute    = virtualNode.performCompute;
        int computeRank       = virtualNode.computeWorldRank;
        int ioRank            = virtualNode.ioWorldRank;
        int ppn               = virtualNode.procsPerNode;
        int mode              = virtualNode.mode;

        //         if(mode == 0) // Virtual nodes *OFF*
        //           return(PartitionDomain(3));

	//	std::cout << "splitting grids to procs" << std::endl;
        if(!geometryInfoPtr){
	  std::cout << "Domain::PartitionDomain:Fatal error: missing geometry info data." << std::endl;
          return(1);
	} 
	//        if((splitMap.size() != numGrids) || slitMap.empty())
        //           return(1);
	const simulation::geometry::info &geomInfo(*geometryInfoPtr); 
        int numDomainProcs  = 1;
	
        if(domainCommunicator.Good()){
          numDomainProcs = domainCommunicator.NProc();
          domainRank     = domainCommunicator.Rank();
        } else {
	  std::cout << "Domain::PartitionDomain:Fatal error: communicator no good." << std::endl;	  
	}
	
	bool alternate_processing = false;
	for(int iGrid = 0;iGrid < numGrids;iGrid++){
	  const std::string &gridName(gridNames[iGrid]);
	  std::pair<int,int> gridInfoIndex(geomInfo.GetGridInfoIndex(gridName));
	  
	  if(gridInfoIndex.first < 0 || gridInfoIndex.second < 0){
	    std::cout << "Domain::PartitionDomain:Fatal error: couldn't find grid index" << std::endl;
	    return(1);
	  }
	  int geomIndex = gridInfoIndex.first;
	  int gridIndex = gridInfoIndex.second;
	  
	  const simulation::geometry::gridinfo &gridInfo(geomInfo.gridInfos[geomIndex][gridIndex]);
	  if(gridInfo.gridSizes.empty()){
	    std::cout << "Domain::PartitionDomain:Fatal error: gridsizes = 0" << std::endl;
	    return(1);
	  }
	  if(gridInfo.decompSizes.empty())
	    alternate_processing = true;
	  
	}
	
	if(alternate_processing)
	  return(PartitionDomain(0));

	// Ask the geometry info data structure for the
	// cartesian topology guidance 
	int numDecomp = 0;
	std::vector<int> numProcsGrid;
        std::vector<std::vector<int> > numVirtualNodesDir;

	for(int iGrid = 0;iGrid < numGrids;iGrid++){
	  const std::string &gridName(gridNames[iGrid]);
	  std::pair<int,int> gridInfoIndex(geomInfo.GetGridInfoIndex(gridName));
	  if(gridInfoIndex.first < 0 || gridInfoIndex.second < 0){
	    std::cout << "Domain::PartitionDomain:Fatal error: couldn't find grid index" << std::endl;
	    return(1);
	  }
	  int geomIndex = gridInfoIndex.first;
	  int gridIndex = gridInfoIndex.second;
	  const simulation::geometry::gridinfo &gridInfo(geomInfo.gridInfos[geomIndex][gridIndex]);
	  if(gridInfo.gridSizes.empty() || gridInfo.decompSizes.empty()){
	    std::cout << "Domain::PartitionDomain:Fatal error: gridsizes or decompsizes = 0" << std::endl;
	    return(1);
	  }

	  const std::vector<size_t> &gridSizes(gridInfo.gridSizes);
	  const std::vector<int>    &decompSizes(gridInfo.decompSizes);	  
	  int numDim = gridSizes.size();
	  int numGridProcs = 1;
	  for(int iDim = 0;iDim < numDim;iDim++)
	    numGridProcs *= decompSizes[iDim];
	  if(numGridProcs <= 0){
	    std::cout << "Domain::PartitionDomain:Fatal error: numgridprocs <= 0" << std::endl;
	    return(1);
	  }
	  numDecomp += numGridProcs;
	  numProcsGrid.push_back(numGridProcs);
	}
	if(numDecomp != numDomainProcs){
	  std::cout << "Domain::PartitionDomain:Fatal error: numdecomp != numdomainprocs" << std::endl;
	  return(1);
	}

	// Now the numprocs for each grid is known - 
	// 2. Split the domain communicator into grid-specific ones
	//	bool gridAssigned = false;
	int procsUsed = 0;
	for(int iGrid = 0;iGrid < numGrids;iGrid++){
	  // std::cout << "Domain Rank (" << domainRank << ") has domain Communicator("
	  // 	    << domainCommunicator.GetCommunicator() << ")" << std::endl;
	  domainCommunicator.Barrier();
	  int gridColor = 0;
	  if(((domainRank >= procsUsed) && ((domainRank-procsUsed) < numProcsGrid[iGrid])))
	    gridColor = 1;
	  domainCommunicator.Split(gridColor,domainRank,gridCommunicators[iGrid]);
	  procsUsed += numProcsGrid[iGrid];
	  if(gridColor == 1){
	    //	    gridAssigned == true;
	    partitionInfo.gridIndexOfLocalGrids.resize(1,iGrid);
	    // std::cout << "domainRank(" << domainRank << ") taking Grid(" << iGrid
	    // 	      << ") with communicator " << gridCommunicators[iGrid].GetCommunicator()
	    // 	      << " and grid Rank(" << gridCommunicators[iGrid].Rank() << ")" << std::endl;
	    gridCommunicators[iGrid].Barrier();
	  }
	  domainCommunicator.Barrier();
        }
	returnCode = 0;
      } else if(inOptions == pcpp::runmode::SPAMD){
        return(PartitionDomain(3));
      } else if(inOptions == 10){  // Fully user prescribed domain splitting
        // ORIGINAL
        if(!virtualNodePtr){
          std::cout << "Domain::PartitionDomain:Fatal error: missing virtual node information."
                    << std::endl;
          return(1);
        }

        pcpp::VirtualNodeType &virtualNode(*virtualNodePtr);

        int numVirtualNodes   = virtualNode.numNodes;
        int numVirtualCompute = virtualNode.computePerNode;
        int numVirtualIO      = virtualNode.ioPerNode;
        int numDomainCompute  = numVirtualNodes * numVirtualCompute;
        int numDomainIO       = numVirtualNodes * numVirtualIO;
        int virtualNodeId     = virtualNode.nodeId;
        int performIO         = virtualNode.performIO;
        int performCompute    = virtualNode.performCompute;
        int computeRank       = virtualNode.computeRank;
        int ioRank            = virtualNode.ioRank;
        int ppn               = virtualNode.procsPerNode;
        int mode              = virtualNode.mode;

        //         if(mode == 0) // Virtual nodes *OFF*
        //           return(PartitionDomain(3));

	//	std::cout << "splitting grids to procs" << std::endl;
        if(!geometryInfoPtr){
	  std::cout << "Domain::PartitionDomain:Fatal error: missing geometry info data." << std::endl;
          return(1);
	} 
	//        if((splitMap.size() != numGrids) || slitMap.empty())
        //           return(1);
	const simulation::geometry::info &geomInfo(*geometryInfoPtr); 
        int numDomainProcs  = 1;
	
        if(domainCommunicator.Good()){
          numDomainProcs = domainCommunicator.NProc();
          domainRank     = domainCommunicator.Rank();
        } else {
	  std::cout << "Domain::PartitionDomain:Fatal error: communicator no good." << std::endl;	  
	}
	
	bool alternate_processing = false;
	for(int iGrid = 0;iGrid < numGrids;iGrid++){
	  const std::string &gridName(gridNames[iGrid]);
	  std::pair<int,int> gridInfoIndex(geomInfo.GetGridInfoIndex(gridName));
	  
	  if(gridInfoIndex.first < 0 || gridInfoIndex.second < 0){
	    std::cout << "Domain::PartitionDomain:Fatal error: couldn't find grid index" << std::endl;
	    return(1);
	  }
	  int geomIndex = gridInfoIndex.first;
	  int gridIndex = gridInfoIndex.second;
	  
	  const simulation::geometry::gridinfo &gridInfo(geomInfo.gridInfos[geomIndex][gridIndex]);
	  if(gridInfo.gridSizes.empty()){
	    std::cout << "Domain::PartitionDomain:Fatal error: gridsizes = 0" << std::endl;
	    return(1);
	  }
	  if(gridInfo.decompSizes.empty())
	    alternate_processing = true;
	  
	}
	
	if(alternate_processing)
	  return(PartitionDomain(0));

	// Ask the geometry info data structure for the
	// cartesian topology guidance 
	int numDecomp = 0;
	std::vector<int> numProcsGrid;
	for(int iGrid = 0;iGrid < numGrids;iGrid++){
	  const std::string &gridName(gridNames[iGrid]);
	  std::pair<int,int> gridInfoIndex(geomInfo.GetGridInfoIndex(gridName));
	  if(gridInfoIndex.first < 0 || gridInfoIndex.second < 0){
	    std::cout << "Domain::PartitionDomain:Fatal error: couldn't find grid index" << std::endl;
	    return(1);
	  }
	  int geomIndex = gridInfoIndex.first;
	  int gridIndex = gridInfoIndex.second;
	  const simulation::geometry::gridinfo &gridInfo(geomInfo.gridInfos[geomIndex][gridIndex]);
	  if(gridInfo.gridSizes.empty() || gridInfo.decompSizes.empty()){
	    std::cout << "Domain::PartitionDomain:Fatal error: gridsizes or decompsizes = 0" << std::endl;
	    return(1);
	  }
	  const std::vector<size_t> &gridSizes(gridInfo.gridSizes);
	  const std::vector<int>    &decompSizes(gridInfo.decompSizes);	  
	  int numDim = gridSizes.size();
	  int numGridProcs = 1;
	  for(int iDim = 0;iDim < numDim;iDim++)
	    numGridProcs *= decompSizes[iDim];
	  if(numGridProcs <= 0){
	    std::cout << "Domain::PartitionDomain:Fatal error: numgridprocs <= 0" << std::endl;
	    return(1);
	  }
	  numDecomp += numGridProcs;
	  numProcsGrid.push_back(numGridProcs);
	}
	if(numDecomp != numDomainProcs){
	  std::cout << "Domain::PartitionDomain:Fatal error: numdecomp != numdomainprocs" << std::endl;
	  return(1);
	}

	// Now the numprocs for each grid is known - 
	// 2. Split the domain communicator into grid-specific ones
	//	bool gridAssigned = false;
	int procsUsed = 0;
	for(int iGrid = 0;iGrid < numGrids;iGrid++){
	  // std::cout << "Domain Rank (" << domainRank << ") has domain Communicator("
	  // 	    << domainCommunicator.GetCommunicator() << ")" << std::endl;
	  domainCommunicator.Barrier();
	  int gridColor = 0;
	  if(((domainRank >= procsUsed) && ((domainRank-procsUsed) < numProcsGrid[iGrid])))
	    gridColor = 1;
	  domainCommunicator.Split(gridColor,domainRank,gridCommunicators[iGrid]);
	  procsUsed += numProcsGrid[iGrid];
	  if(gridColor == 1){
	    //	    gridAssigned == true;
	    partitionInfo.gridIndexOfLocalGrids.resize(1,iGrid);
	    // std::cout << "domainRank(" << domainRank << ") taking Grid(" << iGrid
	    // 	      << ") with communicator " << gridCommunicators[iGrid].GetCommunicator()
	    // 	      << " and grid Rank(" << gridCommunicators[iGrid].Rank() << ")" << std::endl;
	    gridCommunicators[iGrid].Barrier();
	  }
	  domainCommunicator.Barrier();
        }
	returnCode = 0;
      } else if(inOptions == 4){  // Split among virtual nodes

        if(!virtualNodePtr){
          std::cout << "Domain::PartitionDomain:Fatal error: missing virtual node information."
                    << std::endl;
          return(1);
        }

        pcpp::VirtualNodeType &virtualNode(*virtualNodePtr);

        int numVirtualNodes   = virtualNode.numNodes;
        int numVirtualCompute = virtualNode.computePerNode;
        int numVirtualIO      = virtualNode.ioPerNode;
        int numDomainCompute  = numVirtualNodes * numVirtualCompute;
        int numDomainIO       = numVirtualNodes * numVirtualIO;
        int virtualNodeId     = virtualNode.nodeId;
        int performIO         = virtualNode.performIO;
        int performCompute    = virtualNode.performCompute;
        int computeRank       = virtualNode.computeRank;
        int ioRank            = virtualNode.ioRank;
        int ppn               = virtualNode.procsPerNode;
        int mode              = virtualNode.mode;

        if(mode == 0) // Virtual nodes *OFF*
          return(PartitionDomain(3));

        fixtures::CommunicatorType &computeCommunicator(virtualNode.computeCommunicator);
        fixtures::CommunicatorType &ioCommunicator(virtualNode.ioCommunicator);
        
        if((numVirtualNodes  <= 0) || (numVirtualCompute <= 0) ||
           (numVirtualIO     <= 0) || (numDomainCompute  <= 0) ||
           (numDomainIO      <= 0) || (virtualNodeId     <  0)){
          std::cout << "Domain::PartitionDomain:Fatal error: misconfigured virtual node." << std::endl;
          return(1);
        }
        
        if(!geometryInfoPtr){
	  std::cout << "Domain::PartitionDomain:Fatal error: missing geometry info data." << std::endl;
          return(1);
	} 
	const simulation::geometry::info &geomInfo(*geometryInfoPtr); 
        int numDomainProcs  = 1;
	
        if(domainCommunicator.Good()){
          numDomainProcs = domainCommunicator.NProc();
          domainRank     = domainCommunicator.Rank();
        } else {
	  std::cout << "Domain::PartitionDomain:Fatal error: communicator no good." << std::endl;
          return(1);
	}

	// VN
	bool autoPartition(false);

	for(int iGrid = 0;iGrid < numGrids;iGrid++){

	  const std::string &gridName(gridNames[iGrid]);
	  std::pair<int,int> gridInfoIndex(geomInfo.GetGridInfoIndex(gridName));
	  
	  if(gridInfoIndex.first < 0 || gridInfoIndex.second < 0){
	    std::cout << "Domain::PartitionDomain:Fatal error: couldn't find grid index" << std::endl;
	    return(1);
	  }
          
	  int geomIndex = gridInfoIndex.first;
	  int gridIndex = gridInfoIndex.second;
	  
	  const simulation::geometry::gridinfo &gridInfo(geomInfo.gridInfos[geomIndex][gridIndex]);
	  if(gridInfo.gridSizes.empty()){
	    std::cout << "Domain::PartitionDomain:Fatal error: gridsizes = 0" << std::endl;
	    return(1);
	  }
          
          if(gridInfo.decompSizes.empty())
            autoPartition = true; 
	}
	
        if(autoPartition){
          // messageStream something
          //           std::cout << "Domain::PartitionDomain:Status: "
          //                     << "Automatic partitioning of grids will be required." 
          //                     << std::endl;
          if(numGrids > numDomainProcs){
            std::cout << "Switching to SPMD mode." << std::endl;
            virtualNode.mode = pcpp::runmode::SPMD;
            runMode = pcpp::runmode::SPMD;
            return(PartitionDomain(0));
          }
        }
        
        int numDecompProcs = 0;
        int numDecompNodes = 0;
        std::vector<int> numProcsGrid(numGrids,0);
        std::vector<int> numNodesGrid(numGrids,0);
        runMode = virtualNode.mode;

        if(!autoPartition){ // User-prescribed partitioning

          for(int iGrid = 0;iGrid < numGrids;iGrid++){
            const std::string &gridName(gridNames[iGrid]);
            std::pair<int,int> gridInfoIndex(geomInfo.GetGridInfoIndex(gridName));
            int geomIndex = gridInfoIndex.first;
            int gridIndex = gridInfoIndex.second;
            const simulation::geometry::gridinfo &gridInfo(geomInfo.gridInfos[geomIndex][gridIndex]);
            const std::vector<size_t> &gridSizes(gridInfo.gridSizes);
            const std::vector<int>    &decompSizes(gridInfo.decompSizes);	  
            int numDim = gridSizes.size();
            int numGridProcs = 1;
            int numGridNodes = 1;
            for(int iDim = 0;iDim < numDim;iDim++){
              //              numGridProcs *= decompSizes[iDim];
              numGridNodes *= decompSizes[iDim];
            }
            //             if(numGridProcs <= 0){
            //               std::cout << "Domain::PartitionDomain:Fatal error: numgridprocs <= 0" << std::endl;
            //               return(1);
            //             }
            if(numGridNodes <= 0){
              std::cout << "Domain::PartitionDomain:Fatal error: numgridnodes <= 0" << std::endl;
              return(1);
            }
            numGridProcs        = numGridNodes*numVirtualCompute;
            numProcsGrid[iGrid] = numGridProcs;
            numNodesGrid[iGrid] = numGridNodes;
            numDecompProcs += numGridProcs;
            numDecompNodes += numGridNodes;
          }
          if(numDecompNodes != numVirtualNodes){
            std::cout << "Domain::PartitionDomain:Fatal error: numdecomp != numvirtualnodes" << std::endl;
            return(1);
          }
        } else { // Auto-partition
          
          // first count some stuff
          double totalCost = 0.0;
          int totalNumNodes = 0;
          std::vector<double> gridCost(numGrids,0.0);
          //          std::vector<int> numNodesGrid(numGrids,0);
          for(int iGrid = 0;iGrid < numGrids;iGrid++){
            const std::string &gridName(gridNames[iGrid]);
            std::pair<int,int> gridInfoIndex(geomInfo.GetGridInfoIndex(gridName));
            int geomIndex = gridInfoIndex.first;
            int gridIndex = gridInfoIndex.second;
            const simulation::geometry::gridinfo &gridInfo(geomInfo.gridInfos[geomIndex][gridIndex]);
            gridCost[iGrid] = gridInfo.costFactor;
            const std::vector<size_t> &gridSizes(gridInfo.gridSizes);
            size_t numGridPoints = 1;            
            int numDim = gridSizes.size();
            for(int iDim = 0;iDim < numDim;iDim++)
              numGridPoints *= gridSizes[iDim];
            gridCost[iGrid] *= numGridPoints;
            totalCost += gridCost[iGrid];
          }
          std::cout << "Absolute cost of domain grids = " << totalCost << std::endl;
          for(int iGrid = 0;iGrid < numGrids;iGrid++){
            gridCost[iGrid] *= 1.0/totalCost;
            std::cout << "Cost of Grid[" << iGrid << "] = " 
                      << gridCost[iGrid] << std::endl;
            numNodesGrid[iGrid] = gridCost[iGrid]*numVirtualNodes;
            if(numNodesGrid[iGrid] == 0)
              numNodesGrid[iGrid] = 1;
            totalNumNodes += numNodesGrid[iGrid];
            std::cout << "Initial Virtual Node Partition(s): Grid[" << iGrid << "] / " << numNodesGrid[iGrid]
                      << std::endl;
          }
          std::cout << "Allocated " << totalNumNodes << "/" << numVirtualNodes 
                    << " virtual nodes." << std::endl;
          while(totalNumNodes != numVirtualNodes){
            int minCostIndex = 0;
            int maxCostIndex = 0;
            double maxCost = 0;
            double minCost = std::numeric_limits<double>::max();            
            for(int iGrid = 0;iGrid < numGrids;iGrid++){
              double gCost = gridCost[iGrid]/double(numNodesGrid[iGrid]);
              if(gCost > maxCost) {
                maxCost = gCost;
                maxCostIndex = iGrid;
              }
              if(gCost < minCost) {
                minCost = gCost;
                minCostIndex = iGrid;
              }
            }
            if(totalNumNodes < numVirtualNodes){
              std::cout << "Refining grid[" << maxCostIndex 
                        << "]" << std::endl; 
              numNodesGrid[maxCostIndex]++;
              totalNumNodes++;
            } else {
              if(numNodesGrid[minCostIndex] == 1){
                std::cout << "Domain::PartitionDomain:Fatal error: Not enough virtual nodes to distribute grids."
                          << std::endl;
                return(1);
              }
              std::cout << "Coarsening grid[" << maxCostIndex 
                        << "]" << std::endl; 
              numNodesGrid[minCostIndex]--;
              totalNumNodes--;
            }
          }
          for(int iGrid = 0;iGrid < numGrids;iGrid++){
            numProcsGrid[iGrid] = numNodesGrid[iGrid]*numVirtualCompute;
          }
        }
        
        // Now the numComputeProcs for each grid is known  
	// 2. Split the domain communicator into grid-specific ones
	//	bool gridAssigned = false;
	int procsUsed = 0;
        int computeProcsUsed = 0;
        int ioProcsUsed  = 0;
        int nodesUsed = 0;
	for(int iGrid = 0;iGrid < numGrids;iGrid++){
	  // std::cout << "Domain Rank (" << domainRank << ") has domain Communicator("
	  // 	    << domainCommunicator.GetCommunicator() << ")" << std::endl;
	  domainCommunicator.Barrier();

          int numGridCompute = numProcsGrid[iGrid];
          int numGridNodes   = numNodesGrid[iGrid];
          int numGridIO      = numGridNodes*numVirtualIO;

          if(numGridNodes == 0){
            std::cout << "PartitionDomain:Error: numGridNodes = 0 for grid "
                      << iGrid << std::endl;
            return(1);
          }

	  int gridColor        = 0;
          int gridIOColor      = 0;
          int gridComputeColor = 0;

          bool gridNode(false);
          if(((virtualNodeId >= nodesUsed) && ((virtualNodeId - nodesUsed) < numGridNodes))){
            gridColor = 1;
            gridIOColor = performIO;
            gridComputeColor = performCompute;
          }
          if(performCompute){
            computeCommunicator.Split(gridComputeColor,computeRank,gridCommunicators[iGrid]);
          }
          if(performIO){
            ioCommunicator.Split(gridIOColor,ioRank,gridIOCommunicators[iGrid]);
          }
	  domainCommunicator.Split(gridColor,domainRank,gridWorldCommunicators[iGrid]);
	  computeProcsUsed += numGridCompute;
          nodesUsed        += numGridNodes;
          ioProcsUsed      += numGridIO;
          procsUsed        += (numGridNodes*ppn);
	  if(gridColor == 1){
	    //	    gridAssigned == true;
	    partitionInfo.gridIndexOfLocalGrids.resize(1,iGrid);
	    // std::cout << "domainRank(" << domainRank << ") taking Grid(" << iGrid
	    // 	      << ") with communicator " << gridCommunicators[iGrid].GetCommunicator()
	    // 	      << " and grid Rank(" << gridCommunicators[iGrid].Rank() << ")" << std::endl;
	    gridWorldCommunicators[iGrid].Barrier();
	  }
	  domainCommunicator.Barrier();
        }
	returnCode = 0;
      }

      // std::cout << "returncode = " << returnCode << std::endl;
      return(returnCode);
    };

    int InitializeRHS(int inOptions)
    {
      if(gridRHSs.size() != numGrids){
        gridRHSs.resize(numGrids,NULL);
      }
      if(ownRHS.size() != numGrids){
        ownRHS.resize(numGrids,false);
      }
      std::vector<int>::iterator gridIdIt = partitionInfo.gridIndexOfLocalGrids.begin();
      while(gridIdIt != partitionInfo.gridIndexOfLocalGrids.end()){
        int gridIndex = *gridIdIt++;
        gridRHSs[gridIndex] = new RHSType;
        ownRHS[gridIndex] = true;
        RHSType &gridRHS(*gridRHSs[gridIndex]);
        
      }
      return(0);
    };
    
    int Step() { return(stepNumber); };
    void SetStep(int inStep) { stepNumber = inStep; };
    int PrepareToStep(std::ostream &messageStream) { return(0); };

    // Must be called on all threads
    double GetDT(double *outCFL = NULL)
    {
      double timeStep = std::numeric_limits<double>::max();
      double cfl = 1.0;
      //      if(timeSteppingMode == CONSTANT_DT){
      if(true){
        std::vector<int>::iterator gridIdIt = partitionInfo.gridIndexOfLocalGrids.begin();
        while(gridIdIt != partitionInfo.gridIndexOfLocalGrids.end()){
          int gridIndex = *gridIdIt++;
          double gridDT = gridRHSs[gridIndex]->TimeStep();
          if(gridDT < timeStep) timeStep = gridDT;
        }
      }
      return(timeStep);
    };

    // Must be called on all threads
    void ComputeDV(int threadId = 0)
    {
      std::vector<int> &gridIndexOfLocalGrids(partitionInfo.gridIndexOfLocalGrids);
      std::vector<int>::iterator gridIdIt = gridIndexOfLocalGrids.begin();
      while(gridIdIt != gridIndexOfLocalGrids.end()){
        int gridIndex = *gridIdIt++;
        gridRHSs[gridIndex]->ComputeDV(threadId);
      }
    };

    // Must be called on all threads
    // std::string StatusMessage(int threadId = 0)
    // {
    //   std::vector<int> &gridIndexOfLocalGrids(partitionInfo.gridIndexOfLocalGrids);
    //   std::vector<int>::iterator gridIdIt = gridIndexOfLocalGrids.begin();
    //   while(gridIdIt != gridIndexOfLocalGrids.end()){
    //     int gridIndex = *gridIdIt++;
    //     gridRHSs[gridIndex]->ComputeDV(threadId);
    //   }
    // };

    int SetNumThreads(int numThreads)
    {
      std::vector<int> &gridIndexOfLocalGrids(partitionInfo.gridIndexOfLocalGrids);
      std::vector<int>::iterator gridIdIt = gridIndexOfLocalGrids.begin();
      while(gridIdIt != gridIndexOfLocalGrids.end()){
        int gridIndex = *gridIdIt++;
        domainGrids[gridIndex]->SetNumThreads(numThreads);
	if(domainGrids[gridIndex]->PartitionThreadIntervals()){
	  std::cout << "Domain::SetNumThreads:Fatal error: Thread interval partitioning failed on grid "
		    << gridIndex+1 << "." << std::endl;
	  return(1);
	}
        if(gridRHSs[gridIndex]->SetupThreadIntervals()){
	  std::cout << "Domain::SetNumThreads:Fatal error: Thread interval setup failed on rhs "
		    << gridIndex+1 << "." << std::endl;
	  return(1);
	}
      }
      return(0);
    };

    int SetDomainBCs()
    {

      //    int SetDomainBCs(std::vector<BoundaryType> &domainBoundaries,std::vector<BCType> &domainBCs)
      std::vector<int> &gridIndexOfLocalGrids(partitionInfo.gridIndexOfLocalGrids);
      std::vector<int>::iterator gridIdIt = gridIndexOfLocalGrids.begin();
      while(gridIdIt != gridIndexOfLocalGrids.end()){
        int gridIndex = *gridIdIt++;
        if(gridRHSs[gridIndex]->SetDomainBCs(DomainBoundary(gridIndex),BCs())){
          std::cerr << "Domain::SetDomainBCs::Fatal error: Something went wrong setting up "
                    << "boundary stencils for grid(" << gridIndex+1 << ")." << std::endl;
          return(1);
        }
        gridCommunicators[gridIndex].Barrier();
      }
      return(0);
    };

    int InitializeSolution(std::ostream &messageStream)
    {
      std::string modeKey(domainName+":InitMode");
      if(!domainConfig.IsSet(modeKey)){
	messageStream << "simulation::domain: No initialization mode specified.\n";
	return(0);
      }

      std::string initializationMode = domainConfig.GetValue(modeKey);
      std::vector<int> &gridIndexOfLocalGrids(partitionInfo.gridIndexOfLocalGrids);

      if(initializationMode == "HotSpot"){
	std::string hotSpotKey(domainName+":Init:HotSpot");
	std::vector<double> hotSpotParams;
	if(domainConfig.IsSet(ConfigKey(hotSpotKey,"Parameters")))
	  hotSpotParams = domainConfig.GetValueVector<double>(ConfigKey(hotSpotKey,"Parameters"));
	std::vector<int> hotSpotFlags;
	if(domainConfig.IsSet(ConfigKey(hotSpotKey,"Flags")))
	  hotSpotFlags  = domainConfig.GetValueVector<int>(ConfigKey(hotSpotKey,"Flags"));
	
	std::vector<int>::iterator gridIdIt = gridIndexOfLocalGrids.begin();
	while(gridIdIt != gridIndexOfLocalGrids.end()){
	  int gridIndex = *gridIdIt++;
	  
	  if(initflow::DepositBlob(*this,gridIndex,hotSpotParams,
				   hotSpotFlags,messageStream)){
	    messageStream << "simulation::domain: HotSpot initialization failed."
			  << std::endl;
	    return(1);
	  }
	}
	return(0);
      }
      return(0);
    };
    
    int InitThreads()
    {
      std::vector<int> &gridIndexOfLocalGrids(partitionInfo.gridIndexOfLocalGrids);
      std::vector<int>::iterator gridIdIt = gridIndexOfLocalGrids.begin();
      while(gridIdIt != gridIndexOfLocalGrids.end()){
        int gridIndex = *gridIdIt++;
        gridRHSs[gridIndex]->InitThreadIntervals();
      }
      return(0);
    };

    bool UsesOverkit() { return usesOverkit; }
    bool UsesPETSc() { return usesPETSc; }
    bool UsesStella() { return usesStella; }

#ifdef USE_OVERKIT
    ovk_shared_context *OverkitContext() {
      return overkitContextPtr;
    }
    ovk_domain *OverkitDomain() {
      return overkitDomainPtr;
    }
#endif
    int FinalizeGridSetup(int gridIndex,std::ostream &messageStream)
    {

      if(CreateData(gridIndex,messageStream)){
        messageStream << "Domain data creation failed." << std::endl;
        return(1);
      } else {
        messageStream << "Domain data created." << std::endl;
      }

      return(0);
    }

    int FinalizeAdvance(double newTime){
      int advancementMode = 0;
      if(advancementMode == 0){ // just barrier until all grids are done
        // DOH!  this caused failure because MPI not initialized - added a
        // check to see if the domain communicator is "good" before 
        // invoking barrier.
        if(domainCommunicator.Good())
          domainCommunicator.Barrier();
        // could/should validate here that all grids have the same
        // time - at least for this mode.
	//	std::cout << "NewTime = " <<  newTime << std::endl;
        SetTime(newTime);
        return(0);
      } else {
        return(1);
      }
      return(0);
    };
 
    
#ifdef HAVE_CANTERA

   int BuildCanteraObjects(int numThreads, RHSType *inRHS){

     std::vector<CanteraType *> &ptrs(inRHS->GetCanteraGasPtrs());
     ptrs.resize(numThreads,NULL);
     if (!canteraInputString.empty()){
      Cantera::XML_Node *canteraXMLNode = Cantera::get_XML_from_string(canteraInputString);
      for (int i=0; i < numThreads ; i++){
       ptrs[i] = new CanteraType(*canteraXMLNode,canteraId);
      }
      //std::cout << "size of cantera ptr vector " << ptrs.size() << std::endl;
     }

     return(0);

    }

    int CreateCanteraGasPhase(std::string &inputFileBuffer,std::string &id){

      Cantera::XML_Node *canteraXMLNode = Cantera::get_XML_from_string(inputFileBuffer);
      canteraGasPtr = new CanteraType(*canteraXMLNode,id);

      return 0;
    }

    // read the contents of inputFileName and stick them into inputFileBuffer
    // overloaded by parallel_base
    int ReadCanteraInput(std::string &inputFileName, std::string &inputFileBuffer){

      int myRank = 0;
      if(domainCommunicator.Good())
        myRank = domainCommunicator.Rank();

      if(myRank == 0) {
        if(!ix::sys::FILEEXISTS(inputFileName))
          return(1);
        std::string fileBufferString;
        std::ifstream inputFile;
        inputFile.open(inputFileName);
  
        std::ostringstream fileStream;
        fileStream << inputFile.rdbuf();
        inputFileBuffer=fileStream.str();
    
        inputFile.close();
      }

      if(domainCommunicator.Good())
        domainCommunicator.BroadCast(inputFileBuffer);
   
      // MJA
      //std::cout << inputFileBuffer << std::endl;

      return 0;
    }
#endif

    void SetGeometryInfo(geometry::info &geometryInfo)
    { geometryInfoPtr = &geometryInfo; };

    geometry::info GeometryInfo(){
      if(geometryInfoPtr == NULL)
	geometryInfoPtr = new geometry::info;
      return(*geometryInfoPtr);
    };
    
    void SetScalarNames(const std::string &inNames){ scalarNames = inNames; };
    std::vector<std::string> &ScalarNames() { return(scalarNames); };
    const std::vector<std::string> &ScalarNames() const { return(scalarNames); };
    const fixtures::ConfigurationType &GetConfig() const { return(domainConfig); };

    void StatusMessage(const std::vector<std::string> &fieldNames,std::ostream &outStream,std::ofstream *outFile=NULL)
    {
      if(fieldNames.empty())
        return;
      
      double initVal = std::numeric_limits<double>::max();
      

      int numStatusFields = fieldNames.size();

      // First, determine the total number of status components
      int numStatusComponents = 0;
      typename DataDictionaryType::iterator metaIt = stateDataDictionary.begin();
      for(int iField = 0;iField < numStatusFields;iField++){
        const std::string &fieldName(fieldNames[iField]);
        typename StateType::MetaDataType fieldMeta(stateDataDictionary.MetaData(fieldName));
        numStatusComponents += fieldMeta.ncomp;
      }
     
      domainCommunicator.Barrier();


      domainCommunicator.Barrier();

      std::vector<double> fieldMins(numGrids*numStatusComponents,initVal);
      std::vector<double> fieldMaxs(numGrids*numStatusComponents,-initVal);
      std::vector<double> localFieldMins(numGrids*numStatusComponents,initVal);
      std::vector<double> localFieldMaxs(numGrids*numStatusComponents,-initVal);
      
      // loop through the grids
      std::vector<int> &gridIndexOfLocalGrids(partitionInfo.gridIndexOfLocalGrids);
      std::vector<int>::iterator gridIdIt = gridIndexOfLocalGrids.begin();
      while(gridIdIt != gridIndexOfLocalGrids.end()){
        
        int gridIndex = *gridIdIt++;
        int gridID = gridIndex + 1;
        

        GridT &domainGrid(*domainGrids[gridIndex]);
        StateT &gridState(*gridStates[gridIndex]);
        pcpp::CommunicatorType &gridComm(domainGrid.Communicator());
        const pcpp::IndexIntervalType &partitionBufferInterval(domainGrid.PartitionBufferInterval());
        const std::vector<size_t> &bufferSizes(domainGrid.BufferSizes());
        
        pcpp::IndexIntervalType bufferInterval;
        bufferInterval.InitSimple(bufferSizes);
        std::vector<size_t> localIndices;
        bufferInterval.GetFlatIndices(partitionBufferInterval,localIndices);
        size_t numPointsStatus = localIndices.size();
        
        std::vector<double> fieldStats;
        
        // This gets statistics for local partition of current grid
        if(gridState.FieldStats(fieldNames,localIndices,fieldStats)){
          outStream << "Domain::StatusMessage:Error: Could not get field stats for (";
          pcpp::io::DumpContents(outStream,fieldNames," ");
          outStream << ") on Grid(" << gridID << ")" << std::endl;
          domainCommunicator.SetErr(1);
        } else {

          // Parallelize stats
          int numStatMetrics = 4;
          int numStatii = fieldStats.size()/numStatMetrics; // this should be = numFields (for all scalar fields)

          std::vector<double> localMins(numStatusComponents,initVal);
          std::vector<double> localMaxs(numStatusComponents,-initVal);
          std::vector<double> gridMins(numStatusComponents,initVal);
          std::vector<double> gridMaxs(numStatusComponents,-initVal);

          for(int iStat = 0;iStat < numStatusComponents;iStat++){
            localMins[iStat]  = fieldStats[iStat*numStatMetrics];
            localMaxs[iStat]  = fieldStats[iStat*numStatMetrics+1];
          }
          // Reduce over the grid to get grid-specific min/max for current grid
          gridComm.Reduce(localMins,gridMins,ix::comm::DTDOUBLE,ix::comm::MINOP,0);
          gridComm.Reduce(localMaxs,gridMaxs,ix::comm::DTDOUBLE,ix::comm::MAXOP,0);

          // Pack the grid-specific min/max in array with outer scope
          for(int iStat = 0;iStat < numStatusComponents;iStat++){
            localFieldMins[gridIndex*numStatusComponents+iStat] = gridMins[iStat];
            localFieldMaxs[gridIndex*numStatusComponents+iStat] = gridMaxs[iStat];
          }
        }
      }

      if(domainCommunicator.Check()){
        outStream << "Domain::StatusMessage:Error: Could not get field stats on some grid, skipping status."
                  << std::endl;
        return;
      }
      // At this point, all local grid numbers have been collected - use the global 
      // communicator to reduce across all ranks to global rank 0. After this, rank 0 will
      // have the min/max for each of all the grids.
      if(numGrids > 0){
        domainCommunicator.Reduce(localFieldMins,fieldMins,ix::comm::DTDOUBLE,ix::comm::MINOP,0);
        domainCommunicator.Reduce(localFieldMaxs,fieldMaxs,ix::comm::DTDOUBLE,ix::comm::MAXOP,0);
      } else {
        localFieldMins = fieldMins;
        localFieldMaxs = fieldMaxs;
      }

      // fieldMins, fieldMaxs are field-and-grid specific

      //       outStream << "Domain(" << Name() << ") status for fields [ ";
      //       pcpp::io::DumpContents(outStream,fieldNames,"\t");
      //       outStream << " ]" << std::endl;
      std::vector<double> overallMins(numStatusComponents,initVal);
      std::vector<double> overallMaxs(numStatusComponents,-initVal);

      for(int iGrid = 0;iGrid < numGrids;iGrid++){

        if(numGrids > 1)
          outStream << "Grid(" << iGrid+1 << ") : ";
        else
          outStream << "Domain Wide : "; 
        int iStat = 0;

        std::vector<std::string>::const_iterator fieldNameIt = fieldNames.begin();
        while(fieldNameIt != fieldNames.end()){
          size_t iField = fieldNameIt - fieldNames.begin();
          const std::string &fieldName(*fieldNameIt++);
          typename StateType::MetaDataType fieldMeta(stateDataDictionary.MetaData(fieldName));
          int numFieldComponents = fieldMeta.ncomp;
          //          size_t statIndex = iGrid*numStatusFields+iField;
          for(int iComponent = 0;iComponent < numFieldComponents;iComponent++){
            double fieldMin = fieldMins[iGrid*numStatusComponents+iStat];
            double fieldMax = fieldMaxs[iGrid*numStatusComponents+iStat];
            std::string fieldComponentName;
            if(numFieldComponents == 1)
              fieldComponentName = fieldName;
            else {
              std::ostringstream Ostr;
              Ostr << fieldName << "-" << iComponent+1;
              fieldComponentName = Ostr.str();
            }
            if(fieldMin < overallMins[iStat]) overallMins[iStat] = fieldMin;
            if(fieldMax > overallMaxs[iStat]) overallMaxs[iStat] = fieldMax;
            outStream << fieldComponentName << "(" << fieldMin << "," << fieldMax << ")\t";
            iStat++;
          }
        }
        outStream << std::endl;
      }
      if(numGrids > 1){
        outStream << "DomainWide: ";
        std::vector<std::string>::const_iterator fieldNameIt = fieldNames.begin();
        int iStat = 0;
        while(fieldNameIt != fieldNames.end()){
          size_t iField = fieldNameIt - fieldNames.begin();
          const std::string &fieldName(*fieldNameIt++);
          typename StateType::MetaDataType fieldMeta(stateDataDictionary.MetaData(fieldName));
          int numFieldComponents = fieldMeta.ncomp;
          for(int iComp = 0;iComp < numFieldComponents;iComp++){
            std::string fieldComponentName(fieldName);
            if(numFieldComponents > 1){
              std::ostringstream Ostr;
              Ostr << fieldName << "-" << iComp+1;
              fieldComponentName = Ostr.str();
            }
            outStream << fieldComponentName << "(" << overallMins[iStat] << "," 
                      << overallMaxs[iStat] << ")\t";
            iStat++;
          }
        }
      }
      if(outFile){
        for(size_t iField = 0;iField < numStatusComponents;iField++)
          *outFile << overallMins[iField] << " " 
                   << overallMaxs[iField] << " ";
        *outFile << std::endl;
      }
      outStream << std::endl;
    };

    std::vector<std::vector<double> > &GridTimers() { return(gridTimers); };
    std::vector<std::vector<size_t> > &GridCallCounts() { return(gridCallCounts); };

  protected:
    
    int                                     runMode;
    geometry::info                          *geometryInfoPtr;
    pcpp::ParallelGlobalType                *globalPtr;
    std::ostream                            *messageStreamPtr;
    std::string                             domainName;
    OperatorType                            domainOperator;    
    StateType                               domainState;
    StateType                               domainParams;
    gridvector                              domainGrids;
    statevector                             gridStates;
    statevector                             gridTargets;
    statevector                             gridParams;
    rhsvector                               gridRHSs;
    int                                     spatialOrder;
    int                                     stepNumber;
    double                                  domainTime;
    std::vector<bool>                       ownRHS;
    std::vector<bool>                       ownGrid;
    std::vector<bool>                       ownState;
    fixtures::ConfigurationType             domainConfig;
    DataDictionaryType                      stateDataDictionary;
    DataDictionaryType                      paramDataDictionary;
    std::vector<std::string>                gridNames;
    std::vector<std::string>                geometryNames;
    std::vector<std::string>                ioFieldNames;
    std::vector<std::string>                restartFieldNames;
    std::vector<std::vector<BoundaryType> > domainBoundaries;
    std::vector<std::string>                bcNames;
    std::vector<std::string>                boundaryNames;
    std::vector<std::string>                scalarNames;
    std::vector<BCType>                     domainBCs;
    RHSType                                 *rhsPtr;
    RHSType                                 dummyRHS;
    fixtures::CommunicatorType              domainCommunicator;
    std::vector<fixtures::CommunicatorType> gridCommunicators;
    std::vector<fixtures::CommunicatorType> gridIOCommunicators;
    std::vector<fixtures::CommunicatorType> gridWorldCommunicators;    
    int                                     domainRank;
    int                                     numGrids;
    partitioninfo                           partitionInfo;
    pcpp::VirtualNodeType                   *virtualNodePtr;

    bool                                    usesOverkit;
    bool                                    usesPETSc;
    bool                                    usesStella;

#ifdef USE_OVERKIT
    ovk_shared_context                      *overkitContextPtr;
    ovk_domain                              *overkitDomainPtr;
#endif

    InterGridType                           *interGridPtr;

#ifdef HAVE_CANTERA
    std::string                             canteraInputString;
    std::string                             canteraId; 
    CanteraType                             *canteraGasPtr;
#ifdef USE_OMP
    std::vector<CanteraType* >              threadCanteraGasPtr;
#endif
#endif

    EfieldSolverType                        *efieldSolverPtr;
    int                                     efieldSolveInterval;

    std::vector<std::vector<double> >       gridTimers;
    std::vector<std::vector<size_t> >       gridCallCounts;
  };
  
  
  // namespace initializer {
  //   class base {
  //   protected:
  //     ConfigurationType initConfig;
  //   public:
  //     virtual int Configure(ConfigurationType &inConfig)  { return (0); };
  //     template<typename GridT,typename StateT>
  //     virtual int InitializeDomain(domain::base<GridT,StateT> &inDomain){ return (0); };
  //   };
    
  // };
};
};
#endif

