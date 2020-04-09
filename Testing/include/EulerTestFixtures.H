#ifndef __EULER_TEST_FIXTURES_H__
#define __EULER_TEST_FIXTURES_H__

namespace testfixtures {
  namespace euler {

    inline double Gauss(const double &amplitude,const std::vector<double> &width,
                        const std::vector<double> &center, const std::vector<double> &xPos)
    {
      std::vector<double>::const_iterator wIt = width.begin();
      std::vector<double>::const_iterator cIt = center.begin();
      std::vector<double>::const_iterator xIt = xPos.begin();
      double expTerm = 0;
      while(wIt != width.end()){
        double xLoc = (*xIt-*cIt)/(*wIt);
        xIt++;
        wIt++;
        cIt++;
        expTerm += xLoc*xLoc;
      }
      expTerm *= -1.0;
      return(amplitude*std::exp(expTerm));
    };
    
    template<typename GridType,typename StateType>
    int SetupEulerState(const GridType &inGrid,StateType &inState,StateType &inParams,
                        int numScalars,bool withFluxes = false)
    {
      size_t numPointsBuffer = inGrid.BufferSize();
      const std::vector<size_t> &bufferSizes(inGrid.BufferSizes());
      const std::vector<size_t> &gridSizes(inGrid.GridSizes());
      const pcpp::IndexIntervalType &partitionInterval(inGrid.PartitionInterval());
      const pcpp::IndexIntervalType &partitionBufferInterval(inGrid.PartitionBufferInterval());
      int numDim = bufferSizes.size();

      inState.AddField("simTime",'s',1,8,"s");
      
      // Conserved quantities, mass, momentum, and energy densities
      inState.AddField("rho",'n',1,8,"mass");
      inState.AddField("rhoV",'n',numDim,8,"momentum");
      inState.AddField("rhoE",'n',1,8,"energy");
      if(numScalars>0)
        inState.AddField("scalarVars",'n',numScalars,8,"m/M");
      
      // Dependent quantities, pressure, temperature, velocity
      inState.AddField("pressure",'n',1,8,"pressure");
      inState.AddField("temperature",'n',1,8,"temperature");
      inState.AddField("rhom1",'n',1,8,"volume");
      inState.AddField("velocity",'n',numDim,8,"velocity");

      //       if(withFluxes){
      //         // Flux quanitites (if tracking desired)
      //         inState.AddField("massFlux",'n',3,8,"");
      //         inState.AddField("mom1Flux",'n',3,8,"");
      //         inState.AddField("mom2Flux",'n',3,8,"");
      //         inState.AddField("mom3Flux",'n',3,8,"");
      //         inState.AddField("energyFlux",'n',3,8,"");
      //       }
      
      inState.Create(numPointsBuffer,0);
      if(numScalars == 0)
        inState.SetStateFields("rho rhoV rhoE");
      else
        inState.SetStateFields("rho rhoV rhoE scalarVars");
      inState.InitializeFieldHandles();
      
      inParams.AddField("gamma",'s',1,8,"");
      inParams.AddField("inputDT",'s',1,8,"s");
      inParams.AddField("inputCFL",'s',1,8,"");
      inParams.AddField("refRe",'s',1,8,"");
      inParams.AddField("Numbers",'s',3,8,"");
      inParams.AddField("Flags",'s',1,4,"");
      inParams.AddField("optionFlags",'s',1,4,"");
   
      inParams.Create(numPointsBuffer,0);
      
      return(0);
    };


    template<typename GridType,typename StateType>
    int InitializeQuiescentFlow(const GridType &inGrid,StateType &inState,bool everyWhere = false)
    {
      int rhoHandle = inState.GetDataIndex("rho");
      if(rhoHandle < 0)
        return(1);
      pcpp::field::dataset::DataBufferType &rhoData(inState.Field(rhoHandle));
      double *rhoPtr    = rhoData.Data<double>();
      if(!rhoPtr)
        return(1);
      
      int rhoVHandle = inState.GetDataIndex("rhoV");
      if(rhoVHandle < 0)
        return(1);
      pcpp::field::dataset::DataBufferType &rhoVData(inState.Field(rhoVHandle));
      double *rhoVPtr    = rhoVData.Data<double>();
      if(!rhoVPtr)
        return(1);
      
      int rhoEHandle = inState.GetDataIndex("rhoE");
      if(rhoEHandle < 0)
        return(1);
      pcpp::field::dataset::DataBufferType &rhoEData(inState.Field(rhoEHandle));
      double *rhoEPtr    = rhoEData.Data<double>();
      if(!rhoPtr)
        return(1);
      
      size_t numPointsBuffer = inGrid.BufferSize();
      const std::vector<size_t> &bufferSizes(inGrid.BufferSizes());
      const std::vector<size_t> &gridSizes(inGrid.GridSizes());
      const pcpp::IndexIntervalType &partitionInterval(inGrid.PartitionInterval());
      const pcpp::IndexIntervalType &partitionBufferInterval(inGrid.PartitionBufferInterval());
      int numDim = gridSizes.size();
      double rhoE = 1.0/.4;
      if(everyWhere){
        for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
          rhoPtr[iPoint]  = 1.0;
          rhoEPtr[iPoint] = rhoE;
        }
        for(size_t iPoint = 0;iPoint < numDim*numPointsBuffer;iPoint++)
          rhoVPtr[iPoint] = 0.0;
      } else {
        if(numDim == 3){
          size_t iStart = partitionBufferInterval[0].first;
          size_t iEnd   = partitionBufferInterval[0].second;
          size_t jStart = partitionBufferInterval[1].first;
          size_t jEnd   = partitionBufferInterval[1].second;
          size_t kStart = partitionBufferInterval[2].first;
          size_t kEnd   = partitionBufferInterval[2].second;
          size_t nPlane = bufferSizes[0]*bufferSizes[1];
          for(size_t kIndex = kStart;kIndex <= kEnd;kIndex++){
            size_t kBufferIndex = kIndex*nPlane;
            for(size_t jIndex = jStart;jIndex <= jEnd;jIndex++){
              size_t jkBufferIndex = kBufferIndex + jIndex*bufferSizes[0];
              for(size_t iIndex = iStart;iIndex <= iEnd;iIndex++){
                size_t bufferIndex = jkBufferIndex + iIndex;
                rhoPtr[bufferIndex]  = 1.0;
                rhoEPtr[bufferIndex] = rhoE;
                rhoVPtr[bufferIndex] = 0.0;
                rhoVPtr[bufferIndex+numPointsBuffer] = 0.0;
                rhoVPtr[bufferIndex+2*numPointsBuffer] = 0.0;
              }
            }
          }
        } else if (numDim == 2) {
          size_t iStart = partitionBufferInterval[0].first;
          size_t iEnd   = partitionBufferInterval[0].second;
          size_t jStart = partitionBufferInterval[1].first;
          size_t jEnd   = partitionBufferInterval[1].second;
          for(size_t jIndex = jStart;jIndex <= jEnd;jIndex++){
            size_t jBufferIndex = jIndex*bufferSizes[0];
            for(size_t iIndex = iStart;iIndex <= iEnd;iIndex++){
              size_t bufferIndex = jBufferIndex + iIndex;
              rhoPtr[bufferIndex]  = 1.0;
              rhoEPtr[bufferIndex] = rhoE;
              rhoVPtr[bufferIndex] = 0.0;
              rhoVPtr[bufferIndex+numPointsBuffer] = 0.0;
            }
          }
        } else if (numDim == 1) {
          size_t iStart = partitionBufferInterval[0].first;
          size_t iEnd   = partitionBufferInterval[0].second;
          for(size_t iIndex = iStart;iIndex <= iEnd;iIndex++){
            rhoPtr[iIndex]  = 1.0;
            rhoEPtr[iIndex] = rhoE;
            rhoVPtr[iIndex] = 0.0;
          }
        }
      }
      return(0);
    };

  }
}
#endif
