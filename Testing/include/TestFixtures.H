#ifndef __TEST_FIXTURES_H__
#define __TEST_FIXTURES_H__

#include "OperatorTypes.H"
#include "WENO.H"

namespace testfixtures {

  template<typename GridType>
  int CreateSimpleGrid(GridType &inGrid,
                       const std::vector<size_t> &gridSizes)
  {
    int numDim = gridSizes.size();
    std::vector<double> dX(numDim,0);
    size_t numGlobalNodes = 1;
    size_t numGlobalCells = 1;
  
    for(int iDim = 0;iDim < numDim;iDim++){
      numGlobalNodes *= gridSizes[iDim];
      numGlobalCells *= (gridSizes[iDim]-1);
      dX[iDim] = 1.0/static_cast<double>(gridSizes[iDim]-1);
    }

    inGrid.SetGridSizes(gridSizes);
    inGrid.SetGridSpacings(dX);
    inGrid.Finalize();
    return(0);
  };

  template<typename GridType,typename HaloType>
  int CreateSimulationFixtures(GridType &inGrid,HaloType &inHalo,
                               operators::sbp_operator_t &inOperator,
                               const std::vector<size_t> &gridSizes,
                               int interiorOrder,
                               pcpp::CommunicatorType &inComm,
                               std::ostream *messageStreamPtr = NULL)
  {
  
    int numProc = inComm.Size();
    int myRank  = inComm.Rank();
    inGrid.SetType(simulation::grid::UNIRECT);

    if(!messageStreamPtr)
      messageStreamPtr = &std::cout;

    pcpp::CommunicatorType &gridComm(inGrid.Communicator());

    int numDim = gridSizes.size();
    std::vector<double> dX(numDim,0);
    std::vector<double> physicalExtent(2*numDim,1.0);
    size_t numGlobalNodes = 1;
    size_t numGlobalCells = 1;
  
    for(int iDim = 0;iDim < numDim;iDim++){
      physicalExtent[2*iDim] = 0.0;
      numGlobalNodes *= gridSizes[iDim];
      numGlobalCells *= (gridSizes[iDim]-1);
      dX[iDim] = 1.0/static_cast<double>(gridSizes[iDim]-1);
    }
    
    inGrid.SetGridSizes(gridSizes);
    inGrid.SetGridSpacings(dX);
    inGrid.SetPhysicalExtent(physicalExtent);


    operators::sbp::Initialize(inOperator,interiorOrder);
    int boundaryDepth = (inOperator.numStencils-1)/2;
    std::vector<int> haloDepths(2*numDim,boundaryDepth);
    int *stencilStarts = inOperator.stencilStarts;
    int numStencils = inOperator.numStencils;
    for(int iStencil = 0;iStencil < numStencils;iStencil++)
      stencilStarts[iStencil]++;
    inGrid.SetDifferentialOperator(&inOperator);
 
    pcpp::ParallelTopologyInfoType cartInfo;
    cartInfo.numDimensions = numDim;
    pcpp::ParallelTopologyInfoType pbsCartInfo;
    pbsCartInfo.numDimensions = numDim;
    pbsCartInfo.isPeriodic.resize(numDim,1);
    pcpp::comm::SetupCartesianTopology(inComm,pbsCartInfo,gridComm,*messageStreamPtr);

    // Grab the local coordinates and global dimension of the Cartesian topo
    std::vector<int> &cartCoords(gridComm.CartCoordinates());
    std::vector<int> &cartDims(gridComm.CartDimensions());

    // Generate a report of the Cartesian setup and put it the messageStream
    pcpp::report::CartesianSetup(*messageStreamPtr,pbsCartInfo);

    pcpp::IndexIntervalType &partitionInterval(inGrid.PartitionInterval());
    pcpp::IndexIntervalType globalInterval;
    globalInterval.InitSimple(gridSizes);

    // === Partition the grid ===
    if(pcpp::util::PartitionCartesianInterval(globalInterval,cartDims,cartCoords,
                                              partitionInterval,*messageStreamPtr)){
      *messageStreamPtr << "TestFixtures:SetupSimulationFixtures:Error: Partitioning failed." << std::endl;
      return(1);
    }

    std::vector<size_t> extendGrid(2*numDim,0);
    std::vector<bool> isPeriodic(numDim,true);
    size_t numPointsPart = partitionInterval.NNodes();
    std::vector<bool> haveNeighbors(2*numDim,true);

    for(int iDim = 0;iDim < numDim; iDim++){
      if(isPeriodic[iDim]){
        extendGrid[2*iDim]   = haloDepths[2*iDim];
        extendGrid[2*iDim+1] = haloDepths[2*iDim+1];
      } else {
        if(partitionInterval[iDim].first == 0) {
          haveNeighbors[2*iDim] = false;
        } else { 
          extendGrid[2*iDim] = haloDepths[2*iDim];
        }
        if(partitionInterval[iDim].second == (gridSizes[iDim]-1)){
          haveNeighbors[2*iDim+1] = false;
        } else {
          extendGrid[2*iDim + 1] = haloDepths[2*iDim+1];
        }
      }
    }
    
    inGrid.SetDimensionExtensions(extendGrid);

    if(inGrid.Finalize()){
      *messageStreamPtr << "TestFixtures:SetupSimulationFixtures:Error: Grid finalization failed." << std::endl;
      return(1);
    }
    inGrid.PartitionThreadIntervals();
    inGrid.GenerateCoordinates(*messageStreamPtr);
    
    const std::vector<size_t> &bufferSizes(inGrid.BufferSizes());
    const pcpp::IndexIntervalType &partitionBufferInterval(inGrid.PartitionBufferInterval());  
    size_t numPointsBuffer = inGrid.BufferSize();

    // Set up the HALO data structure 
    //  halo_t testHalo(globalInterval,pbsPartInterval);
    inHalo.SetGridInterval(globalInterval);
    inHalo.SetPartitionInterval(partitionInterval);

    // This forces the halo to have neighbors even on 
    // partitions with domain boundaries.
    inHalo.SetNeighbors(haveNeighbors);
    inHalo.SetPeriodicDirs(pbsCartInfo.isPeriodic);

    std::vector<pcpp::IndexIntervalType> remoteHaloExtents(inHalo.CreateRemoteHaloExtents(extendGrid));
    std::vector<pcpp::IndexIntervalType> localHaloExtents(inHalo.CreateLocalHaloExtents(extendGrid));
  
    inHalo.SetLocalBufferSizes(bufferSizes);
    inHalo.SetLocalPartitionExtent(partitionBufferInterval);
  
    inHalo.SetRemoteHaloExtents(remoteHaloExtents);
    inHalo.SetLocalHaloExtents(localHaloExtents);
  
    const std::vector<pcpp::IndexIntervalType> &remoteHaloBufferExtents(inHalo.RemoteHaloBufferExtents());
    const std::vector<pcpp::IndexIntervalType> &localHaloBufferExtents(inHalo.LocalHaloBufferExtents());
  
    inHalo.CreateSimpleSendIndices();
    inHalo.CreateSimpleRecvIndices();

    //     inGrid.CreateStencilConnectivity();
    //     inGrid.CreateMask();
    inGrid.InitializeStencilConnectivity();
    std::ostringstream outStream;
    inGrid.ComputeMetrics(outStream);
    
    return(0);
    
  };

  template<typename GridType,typename HaloType>
  int CreateSimulationFixtures(GridType &inGrid,HaloType &inHalo,
                               operators::sbp_operator_t &inOperator,
                               WENO::CoeffsWENO &coeffsWENO,
                               const std::vector<size_t> &gridSizes,
                               int interiorOrder,
                               pcpp::CommunicatorType &inComm,
                               std::ostream *messageStreamPtr = NULL,
                               bool periodic = true)
  {
  
    int numProc = inComm.Size();
    int myRank  = inComm.Rank();
  
    if(!messageStreamPtr)
      messageStreamPtr = &std::cout;

    pcpp::CommunicatorType &gridComm(inGrid.Communicator());
    inGrid.SetType(simulation::grid::UNIRECT);
    int numDim = gridSizes.size();
    std::vector<double> dX(numDim,0);
    size_t numGlobalNodes = 1;
    size_t numGlobalCells = 1;
  
    std::vector<double> &physicalExtent(inGrid.PhysicalExtent());
    if(physicalExtent.empty()){
      physicalExtent.resize(numDim*2,0);
      for(int iDim = 0;iDim < numDim;iDim++)
        physicalExtent[2*iDim+1] = 20.0;
    }

    for(int iDim = 0;iDim < numDim;iDim++){
      numGlobalNodes *= gridSizes[iDim];
      numGlobalCells *= (gridSizes[iDim]-1);
      dX[iDim] = (physicalExtent[2*iDim+1] - physicalExtent[2*iDim])/static_cast<double>(gridSizes[iDim]-1);
    }

    inGrid.SetGridSizes(gridSizes);
    inGrid.SetGridSpacings(dX);

    operators::sbp::Initialize(inOperator,interiorOrder);
    int boundaryDepth = coeffsWENO.HaloWidth();
    int boundaryWidth = inOperator.boundaryWidth;
    boundaryDepth = std::max(boundaryWidth,boundaryDepth);
    int numStencils = inOperator.numStencils;
    int *stencilStarts = inOperator.stencilStarts;
    for(int iStencil = 0;iStencil < numStencils;iStencil++)
      stencilStarts[iStencil]++;
    
    std::vector<int> haloDepths(2*numDim,boundaryDepth);

    pcpp::ParallelTopologyInfoType cartInfo;
    cartInfo.numDimensions = numDim;
    pcpp::ParallelTopologyInfoType pbsCartInfo;
    pbsCartInfo.numDimensions = numDim;
    if (periodic) {
      pbsCartInfo.isPeriodic.resize(numDim,1);
    } else {
      pbsCartInfo.isPeriodic.resize(numDim,0);
    }
    pcpp::comm::SetupCartesianTopology(inComm,pbsCartInfo,gridComm,*messageStreamPtr);

    // Grab the local coordinates and global dimension of the Cartesian topo
    std::vector<int> &cartCoords(gridComm.CartCoordinates());
    std::vector<int> &cartDims(gridComm.CartDimensions());

    // Generate a report of the Cartesian setup and put it the messageStream
    pcpp::report::CartesianSetup(*messageStreamPtr,pbsCartInfo);

    pcpp::IndexIntervalType &partitionInterval(inGrid.PartitionInterval());
    pcpp::IndexIntervalType globalInterval;
    globalInterval.InitSimple(gridSizes);

    // === Partition the grid ===
    if(pcpp::util::PartitionCartesianInterval(globalInterval,cartDims,cartCoords,
                                              partitionInterval,*messageStreamPtr)){
      *messageStreamPtr << "TestFixtures:SetupSimulationFixtures:Error: Partitioning failed." << std::endl;
      return(1);
    }

    std::vector<size_t> extendGrid(2*numDim,0);
    std::vector<bool> isPeriodic(numDim, periodic);
    size_t numPointsPart = partitionInterval.NNodes();
  
    for(int iDim = 0;iDim < numDim; iDim++){
      if(isPeriodic[iDim]){
        extendGrid[2*iDim]   = haloDepths[2*iDim];
        extendGrid[2*iDim+1] = haloDepths[2*iDim+1];
      } else {
        if(partitionInterval[iDim].first != 0)
          extendGrid[2*iDim] = haloDepths[2*iDim];
        if(partitionInterval[iDim].second != (gridSizes[iDim]-1))
          extendGrid[2*iDim + 1] = haloDepths[2*iDim+1];
      }
    }

    inGrid.SetDimensionExtensions(extendGrid);
    if(inGrid.Finalize()){
      *messageStreamPtr << "TestFixtures:SetupSimulationFixtures:Error: Grid finalization failed." << std::endl;
      return(1);
    }
    inGrid.PartitionThreadIntervals();
    inGrid.SetDifferentialOperator(&inOperator);
    //     inGrid.CreateStencilConnectivity();
    //     inGrid.CreateMask();
    inGrid.InitializeStencilConnectivity();
    std::ostringstream outStream;
    inGrid.ComputeMetrics(outStream);

    const std::vector<size_t> &bufferSizes(inGrid.BufferSizes());
    const pcpp::IndexIntervalType &partitionBufferInterval(inGrid.PartitionBufferInterval());  
    size_t numPointsBuffer = inGrid.BufferSize();

    // Set up the HALO data structure 
    //  halo_t testHalo(globalInterval,pbsPartInterval);
    inHalo.SetGridInterval(globalInterval);
    inHalo.SetPartitionInterval(partitionInterval);

    // This forces the halo to have neighbors even on 
    // partitions with domain boundaries.
    std::vector<bool> haveNeighbors(2*numDim,true);
    inHalo.SetNeighbors(haveNeighbors);
    inHalo.SetPeriodicDirs(pbsCartInfo.isPeriodic);

    std::vector<pcpp::IndexIntervalType> remoteHaloExtents(inHalo.CreateRemoteHaloExtents(extendGrid));
    std::vector<pcpp::IndexIntervalType> localHaloExtents(inHalo.CreateLocalHaloExtents(extendGrid));
  
    inHalo.SetLocalBufferSizes(bufferSizes);
    inHalo.SetLocalPartitionExtent(partitionBufferInterval);
  
    inHalo.SetRemoteHaloExtents(remoteHaloExtents);
    inHalo.SetLocalHaloExtents(localHaloExtents);
  
    const std::vector<pcpp::IndexIntervalType> &remoteHaloBufferExtents(inHalo.RemoteHaloBufferExtents());
    const std::vector<pcpp::IndexIntervalType> &localHaloBufferExtents(inHalo.LocalHaloBufferExtents());
  
    inHalo.CreateSimpleSendIndices();
    inHalo.CreateSimpleRecvIndices();
  
    return(0);

  };

  template<typename GridType>
  int CreateSerialSimulationFixtures(GridType &inGrid,
                                     operators::sbp_operator_t &inOperator,
                                     const std::vector<size_t> &gridSizes,
                                     int interiorOrder,
                                     std::ostream *messageStreamPtr = NULL)
  {
  
    
    if(!messageStreamPtr)
      messageStreamPtr = &std::cout;
    
    int numDim = gridSizes.size();
    std::vector<double> dX(numDim,0);
    size_t numGlobalNodes = 1;
    size_t numGlobalCells = 1;
    inGrid.SetType(simulation::grid::UNIRECT);

    for(int iDim = 0;iDim < numDim;iDim++){
      numGlobalNodes *= gridSizes[iDim];
      numGlobalCells *= (gridSizes[iDim]-1);
      dX[iDim] = 1.0/static_cast<double>(gridSizes[iDim]-1);
    }
    
    inGrid.SetGridSizes(gridSizes);

    operators::sbp::Initialize(inOperator,interiorOrder);
    int boundaryDepth = (inOperator.numStencils-1)/2;
    std::vector<int> haloDepths(2*numDim,boundaryDepth);

    pcpp::IndexIntervalType &partitionInterval(inGrid.PartitionInterval());
    pcpp::IndexIntervalType globalInterval;
    globalInterval.InitSimple(gridSizes);

    std::vector<size_t> extendGrid(2*numDim,0);
    std::vector<bool> isPeriodic(numDim,true);
    size_t numPointsPart = partitionInterval.NNodes();
  
    for(int iDim = 0;iDim < numDim; iDim++){
      if(isPeriodic[iDim]){
        extendGrid[2*iDim]   = haloDepths[2*iDim];
        extendGrid[2*iDim+1] = haloDepths[2*iDim+1];
      } else {
        if(partitionInterval[iDim].first == 0)
          extendGrid[2*iDim] = haloDepths[2*iDim];
        if(partitionInterval[iDim].second == (gridSizes[iDim]-1))
          extendGrid[2*iDim + 1] = haloDepths[2*iDim+1];
      }
    }

    inGrid.SetDimensionExtensions(extendGrid);
    if(inGrid.Finalize()){
      *messageStreamPtr << "TestFixtures:SetupSimulationFixtures:Error: Grid finalization failed." << std::endl;
      return(1);
    }
    inGrid.PartitionThreadIntervals();
    const std::vector<size_t> &bufferSizes(inGrid.BufferSizes());
    const pcpp::IndexIntervalType &partitionBufferInterval(inGrid.PartitionBufferInterval());  
    size_t numPointsBuffer = inGrid.BufferSize();

    return(0);

  };

#ifdef USE_OVERKIT
  inline int CreateOverkitContext(pcpp::CommunicatorType &contextComm, ovk_shared_context
    *&overkitContext) {

    ovk_context_params *contextParams;
    ovkCreateContextParams(&contextParams);
    ovkSetContextParamComm(contextParams, contextComm.GetCommunicator());
    ovkSetContextParamStatusLoggingThreshold(contextParams, 4);

    ovk_context *tempContextPtr;
    ovk_error error;
    ovkCreateContext(&tempContextPtr, &contextParams, &error);
    if (error != OVK_ERROR_NONE) return(1);

    ovkShareContext(&tempContextPtr, &overkitContext);

    return(0);

  }

  inline void DestroyOverkitContext(ovk_shared_context *&overkitContext) {

    ovkResetSharedContext(&overkitContext);

  }
#endif

}
#endif
