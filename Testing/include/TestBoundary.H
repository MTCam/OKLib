#ifndef __TEST_BOUNDARY_H__
#define __TEST_BOUNDARY_H__

#include "State.H"
#include "Grid.H"

namespace simulation {
  namespace testdomain {
    namespace testboundary {  

      // class boundary::base
      template <typename BCType> class base {
        
      public:
        
        typedef simulation::state::base                    StateType;
        typedef simulation::grid::parallel_blockstructured GridType;

        base() : statePtr(NULL), paramPtr(NULL), bcPtr(NULL), 
                 myState(false), myParams(false), myBC(false)
        {
          Create();
        };
        
        
//         base(const base &inBase) : statePtr(NULL), paramPtr(NULL), bcPtr(NULL), 
//                                    myState(false), myParams(false), myBC(false)
//         {
//           Copy(inBase);
//         };
        
//         void Copy(const base &inBase)
//         {
//           Destroy();
//           statePtr = inBase.statePtr;
//           paramPtr = inBase.paramPtr;
//           bcPtr    = inBase.bcPtr;
//           gridPtr  = inBase.gridPtr;
//         };

        void SetBC(BCType &inBC)
        {
          DestroyBC();
          bcPtr = &inBC;
        }

        void SetState(StateType &inState)
        {
          DestroyState();
          statePtr = &inState;
        };
        
        void SetParams(StateType &inParams)
        {
          DestroyParams();
          paramPtr = &inParams;
        };
        
        void DestroyState()
        {
          if(myState && (statePtr != NULL))
            delete statePtr;
          statePtr = NULL;
          myState = false;
        };
        
        void DestroyParams()
        {
          if(myParams && (paramPtr != NULL))
            delete paramPtr;
          paramPtr = NULL;
          myParams = false;
        };

        void DestroyBC()
        {
          if(myBC && (bcPtr != NULL))
            delete bcPtr;
          bcPtr = NULL;
          myBC = false;
        };
        
        void Destroy()
        {
          DestroyState();
          DestroyParams();
          DestroyBC();
        };
        
        void Create()
        {
          Destroy();
          statePtr = new StateType;
          myState = true;
          paramPtr = new StateType;
          myParams = true;
          bcPtr    = new BCType;
          myBC     = true;
        };
        

        void CreateData()
        {
          if(bcPtr){
            SetParams(bcPtr->Params()); // Set the actual bc params object to global one
            if(statePtr){
              statePtr->SetMetaData(bcPtr->StateDictionary());// Copy state metadata from BC prototype
              if(!(bcPtr->StateDictionary().empty())){ 
                if(gridPtr){
                  if(gridRegionID >= 0) { // create boundary-specific state for this BC
                    std::vector<grid::subregion> &boundaryRegions(gridPtr->SubRegions());
                    grid::subregion &boundaryRegion(boundaryRegions[gridRegionID]);
                    size_t numPointsBoundary = boundaryRegion.regionInterval.NNodes();
                    if(numPointsBoundary > 0)
                      statePtr->Create(numPointsBoundary,0);
                  }
                }
              }
            }
          }
        };

        StateType &State() { return(*statePtr); };
        const StateType &State() const {return(*statePtr); };
        StateType &Params() { return(*paramPtr); };
        const StateType &Params() const { return(*paramPtr); };
        
        GridType &Grid() { return(*gridPtr); };
        const GridType &Grid() const {return(*gridPtr); };
        void SetGrid(GridType &inGrid) { gridPtr = &inGrid; };

        BCType &BC(){ return(*bcPtr); };
        const BCType &BC() const {return(*bcPtr); };

        int GridRegionID() const {return(gridRegionID); };
        void SetGridRegionID(int inID) { gridRegionID = inID; };

        const std::vector<pcpp::IndexIntervalType> &ThreadPartitionIntervals() const
        { return(threadPartitionIntervals); };
        std::vector<pcpp::IndexIntervalType> &ThreadPartitionIntervals()
        { return(threadPartitionIntervals); };
        const std::vector<pcpp::IndexIntervalType> &ThreadBufferIntervals() const
        { return(threadBufferIntervals); };
        std::vector<pcpp::IndexIntervalType> &ThreadBufferIntervals()
        { return(threadBufferIntervals); };
        const pcpp::IndexIntervalType &BoundaryPartitionInterval() const
        { return(boundaryPartitionInterval);};
        pcpp::IndexIntervalType &BoundaryPartitionInterval()
        { return(boundaryPartitionInterval);};
        const pcpp::IndexIntervalType &PartitionBufferInterval() const
        { return(partitionBufferInterval); };
        pcpp::IndexIntervalType &PartitionBufferInterval()
        { return(partitionBufferInterval); };
        
        void SetName(const std::string &inName) { boundaryName = inName; };
        const std::string &Name() const {return(boundaryName); };

      protected:
        
        int         gridRegionID;
        BCType     *bcPtr;
        StateType  *statePtr;
        StateType  *paramPtr;
        GridType   *gridPtr;
        std::string boundaryName;

        pcpp::IndexIntervalType boundaryPartitionInterval;
        pcpp::IndexIntervalType partitionBufferInterval;
        std::vector<pcpp::IndexIntervalType> threadPartitionIntervals;
        std::vector<pcpp::IndexIntervalType> threadBufferIntervals;
        
        //         StateType::MetaDataSetType     stateDataDictionary;
        //         StateType::MetaDataSetType     paramDataDictionary;
        
        bool myState;
        bool myParams;
        bool myBC;
        
      };

//       template<typename BcT>
//       class bcbase : public base
//       {
//       public:
//         typedef typename BcT BCType; 
//         bcbase() : base() {};
//         BCType &BC(){ return(*bcPtr); };
//         const BCType &BC() const {return(*bcPtr); };        
//         void SetBC(BCType &inBC)
//         {
//           DestroyBC();
//           bcPtr = &inBC;
//         }
//         void Create()
//         {
//           DestroyBC();
//           base::Create();
//           bcPtr    = new BCType;
//           myBC     = true;
//         };
//         void Destroy()
//         {
//           DestroyBC();
//           base::Destroy();
//         };
//         void DestroyBC()
//         {
//           if(myBC && (bcPtr != NULL))
//             delete bcPtr;
//           bcPtr = NULL;
//           myBC = false;
//         };
//       protected:
//         BCType     *bcPtr;
//         bool myBC;
//       };

    }
  }
}
#endif
