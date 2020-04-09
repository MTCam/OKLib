#ifndef __TEST_BOUNDARY_CONDITION_H__
#define __TEST_BOUNDARY_CONDITION_H__

#include "State.H"
#include "Grid.H"

namespace simulation {
  namespace testdomain {
    namespace testboundary {  

      namespace nsbc {
      // Each new domain type needs a namespace for its
      // own boundary conditions.  The first one *must*
      // be HOLE=0, named "HOLE".
      //
      // Caution, NSCBC boundaries are hard-coded to their index (starting with 6)
      // inside NavierStokesRHSWENO.f90. If new boundaries are added before any NSCBC 
      // boundaries, this must be updated to match the new mapping!  
      //
        enum BCType { HOLE=0,SAT_FARFIELD,SAT_NOSLIP_ISOTHERMAL,SAT_SLIP_ADIABATIC, 
                      SPONGE,SATBC,NSCBC_NOSLIP_ISOTHERMAL,NSCBC_SLIP_ADIABATIC, 
                      NSCBC_INFLOW_SUPERSONIC,NSCBC_INFLOW_SUBSONIC,NSCBC_INFLOW_VELTEMP, 
                      NSCBC_OUTFLOW_NOREFLECT,AXISYMMETRY,PERIODIC,EF_DIRICHLET,EF_NEUMANN,
                      NUMBCTYPE};

/*
        // support for legacy boundary condition names, 
        // when adding a new type make sure to increment NUMBCTYPELEGACY
        enum BCTypeLegacy {NOSLIP_ISOTHERMAL=NSCBC_NOSLIP_ISOTHERMAL,
                           SLIP_ADIABATIC=NSCBC_SLIP_ADIABATIC, 
                           INFLOW_SUPERSONIC=NSCBC_INFLOW_SUPERSONIC,
                           INFLOW_SUBSONIC=NSCBC_INFLOW_SUBSONIC,
                           INFLOW_VELTEMP=NSCBC_INFLOW_VELTEMP, 
                           OUTFLOW_NOREFLECT=NSCBC_OUTFLOW_NOREFLECT,
                           NUMBCTYPELEGACY=6};
*/
        static const int numLegacyBCNames = 14;

        // boundary conditions names
        static const char *bcNames[] = 
        {
          "HOLE",
          "SAT_FARFIELD",
          "SAT_NOSLIP_ISOTHERMAL",
          "SAT_SLIP_ADIABATIC", 
          "SPONGE",
          "DONOTUSE",
          "NSCBC_NOSLIP_ISOTHERMAL",
          "NSCBC_SLIP_ADIABATIC", 
          "NSCBC_INFLOW_SUPERSONIC",
          "NSCBC_INFLOW_SUBSONIC",
          "NSCBC_INFLOW_VELTEMP", 
          "NSCBC_OUTFLOW_NOREFLECT",
          "AXISYMMETRY",
          "PERIODIC",
          "EF_DIRICHLET",
          "EF_NEUMANN",
          "UNKNOWN"
        };

        // legacy names for backwards compatability
        static const char *bcNamesLegacy[] = 
        {
          "HOLE",
          "SAT_FARFIELD",
          "SAT_NOSLIP_ISOTHERMAL",
          "SAT_SLIP_ADIABATIC", 
          "SPONGE",
          "DONOTUSE",
          "NOSLIP_ISOTHERMAL",
          "SLIP_ADIABATIC", 
          "INFLOW_SUPERSONIC",
          "INFLOW_SUBSONIC",
          "INFLOW_VELTEMP", 
          "OUTFLOW_NOREFLECT",
          "AXISYMMETRY",
          "PERIODIC",
          "UNKNOWN"
        };
        
        // number of expected parameter numbers
        static const int bcNumParamNumbers[] = 
        {
          0, //"HOLE",
          2, //"SAT_FARFIELD",
          3, //"SAT_NOSLIP_ISOTHERMAL",
          2, //"SAT_SLIP_ADIABATIC", 
          2, //"SPONGE",
          0, //"DONOTUSE",
          1, //"NSCBC_NOSLIP_ISOTHERMAL",
          0, //"NSCBC_SLIP_ADIABATIC", 
          0, //"NSCBC_INFLOW_SUPERSONIC",
          1, //"NSCBC_INFLOW_SUBSONIC",
          1, //"NSCBC_INFLOW_VELTEMP", 
          0, //"NSCBC_OUTFLOW_NOREFLECT",
          2, //"AXISYMMETRY",
          0, //"PERIODIC",
          0, //"EF_DIRICHLET",
          0, //"EF_NEUMANN",
          0 //"UNKNOWN"
        };

        // number of expected parameter flags
        static const int bcNumParamFlags[] = 
        {
          0, //"HOLE",
          1, //"SAT_FARFIELD",
          1, //"SAT_NOSLIP_ISOTHERMAL",
          1, //"SAT_SLIP_ADIABATIC", 
          1, //"SPONGE",
          0, //"DONOTUSE",
          1, //"NSCBC_NOSLIP_ISOTHERMAL",
          1, //"NSCBC_SLIP_ADIABATIC", 
          1, //"NSCBC_INFLOW_SUPERSONIC",
          1, //"NSCBC_INFLOW_SUBSONIC",
          1, //"NSCBC_INFLOW_VELTEMP", 
          1, //"NSCBC_OUTFLOW_NOREFLECT",
          1, //"AXISYMMETRY",
          1, //"PERIODIC",
          1, //"EF_DIRICHLET",
          1, //"EF_NEUMANN",
          0 //"UNKNOWN"
        };
      }

      namespace bc {
        
        const static int INVALIDBCTYPE = -1;
        const static int HOLE = 0;

        int ResolveBCName(const std::string &inName);
        std::string ResolveBCType(int inType);
        
        class base {
          
        public:
          
          typedef state::base StateType;
          typedef StateType::MetaDataSetType DataDictionaryType;

          base() : bcType(0), paramPtr(NULL), 
                   myParams(false), myState(false),
                   paramData(NULL), flagData(NULL)
          {
            Create();
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
            if(myParams && (paramPtr != NULL))
              delete paramPtr;
            paramPtr = NULL;
            myParams = false;
          };
          void DestroyState()
          {
            if(myState && (statePtr != NULL))
              delete statePtr;
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
          int NumParamNumbers() const {return(numParamNumbers);};
          void SetBCType(int inType) { bcType = inType; };

          std::string &BCName() { return(bcName); };
          const std::string &BCName() const { return(bcName); };
          void SetBCName(const std::string &inName){ bcName = inName; };
          
          DataDictionaryType &ParamDictionary() { return(paramDictionary); };
          const DataDictionaryType &ParamDictionary() const { return(paramDictionary); };
          DataDictionaryType &StateDictionary() { return(stateDictionary); };
          const DataDictionaryType &StateDictionary() const { return(stateDictionary); };
          
          void InitInternalData() 
          {
            if(paramPtr){
              int dataIndex = paramPtr->GetDataIndex("Numbers");
              numParamNumbers = 0;
              if(dataIndex >= 0){
                paramData = paramPtr->GetFieldBuffer<double>(dataIndex);
                numParamNumbers = paramPtr->Meta()[dataIndex].ncomp;
              }
              numParamFlags = 0;
              dataIndex = paramPtr->GetDataIndex("Flags");
              if(dataIndex >= 0){
                flagData = paramPtr->GetFieldBuffer<int>(dataIndex);
                numParamFlags = paramPtr->Meta()[dataIndex].ncomp;
              }
            }
          };

        protected:
          
          int         bcType;
          std::string bcName;

          StateType  *paramPtr;
          StateType  *statePtr;

          const double     *paramData;
          const int        *flagData;
          int numParamNumbers;
          int numParamFlags;

          DataDictionaryType     stateDictionary;
          DataDictionaryType     paramDictionary;
          
          bool myParams;
          bool myState;

        }; // class bc::base

      } // namespace bc
    }
  }
}
#endif
