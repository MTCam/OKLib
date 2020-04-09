#include "Testing.H"
#include "PlasCom2.H"
#include "PlasComCM.H"
#include "PC2Util.H"
#include "Report.H"
#include "MetricKernels.H"


typedef plascom2::grid_t               grid_t;
typedef plascom2::halo_t               halo_t;
typedef pcpp::ParallelGlobalType       global_t;
typedef pcpp::CommunicatorType         comm_t;
typedef pcpp::IndexIntervalType        interval_t;
typedef simulation::state::base        state_t;
typedef operators::sbp_operator_t      operator_t;


void TestMetricOps_PointMetric(ix::test::results &serialUnitResults)
{
  std::string testFunctionName("TestMetricOps_PointMetric");
  std::cout << testFunctionName << std::endl;

  global_t myGlobal;

  std::ostringstream messageStream;

  myGlobal.Init(testFunctionName);
  myGlobal.SetVerbLevel(3);
  myGlobal.Profiling(true);
 
  size_t numPointsBuffer = 10;
  
  int uniformResult[]     = {1,1,1,1,1,1,1,1};
  int rectilinearResult[] = {1,1,1,1,1,1,1,1};
  int curvilinearResult[] = {1,1,1,1,1,1,1,1};

  
  double TOL = 5e-16;

  for(int numDim = 2;numDim < 4;numDim++){
    size_t numUniformComponents = 0;
    size_t numRectilinearComponents = numDim;
    size_t numCurvilinearComponents = numDim*numDim;

    size_t uniformMetricSize = numDim;
    size_t rectilinearMetricSize = numRectilinearComponents*numPointsBuffer;
    size_t curvilinearMetricSize = numCurvilinearComponents*numPointsBuffer;

    std::vector<double> uniformGridMetric(numDim,1.0);
    std::vector<double> rectilinearGridMetric(rectilinearMetricSize,0.0);
    std::vector<double> curvilinearGridMetric(curvilinearMetricSize,0.0);
    std::vector<double> curvilinearMags(numPointsBuffer,0);
   
    for(int iDim = 0;iDim < numDim;iDim++){
      uniformGridMetric[iDim] = double(iDim+1);
      size_t iOffset  = iDim*numPointsBuffer;
      size_t iOffset2 = iDim*numDim*numPointsBuffer;
      for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
        size_t rectilinearOffset = iOffset + iPoint;
        rectilinearGridMetric[rectilinearOffset] = double(iDim*numPointsBuffer+iPoint+1);
        for(int jDim = 0;jDim < numDim;jDim++){
          size_t jOffset = jDim*numPointsBuffer;
          size_t curvilinearOffset = iOffset2 + jOffset + iPoint;
          curvilinearGridMetric[curvilinearOffset] = double((iDim*numDim + jDim)*numPointsBuffer + iPoint + 1);
          curvilinearMags[iPoint] += curvilinearGridMetric[curvilinearOffset]*curvilinearGridMetric[curvilinearOffset];
        }
        curvilinearMags[iPoint] = std::sqrt(curvilinearMags[iPoint]);
      }
    }

    std::vector<double> uniformPointMetric(numDim*numDim,0);
    std::vector<double> rectilinearPointMetric(numDim*numDim,0);
    std::vector<double> curvilinearPointMetric(numDim*numDim,0);
    std::vector<double> uniformBoundaryMetric(numDim*numDim,0);
    std::vector<double> rectilinearBoundaryMetric(numDim*numDim,0);
    std::vector<double> curvilinearBoundaryMetric(numDim*numDim,0);
    std::vector<double> uniformEigenVectors(numDim*numDim,0);
    std::vector<double> rectilinearEigenVectors(numDim*numDim,0);
    std::vector<double> curvilinearEigenVectors(numDim*numDim,0);
    std::vector<double> uniformMetricMags(numDim,0);
    std::vector<double> rectilinearMetricMags(numDim,0);
    std::vector<double> curvilinearMetricMags(numDim,0);
  
    for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
      size_t pointID = iPoint + 1;
      int gridType = 1;

      FC_MODULE(metricops,getpointmetric,METRICOPS,GETPOINTMETRIC)
        (&numDim,&numPointsBuffer,&pointID,&gridType,&uniformGridMetric[0],&uniformPointMetric[0]);
      gridType = 2;
      FC_MODULE(metricops,getpointmetric,METRICOPS,GETPOINTMETRIC)
        (&numDim,&numPointsBuffer,&pointID,&gridType,&rectilinearGridMetric[0],&rectilinearPointMetric[0]);
      gridType = 3;
      FC_MODULE(metricops,getpointmetric,METRICOPS,GETPOINTMETRIC)
        (&numDim,&numPointsBuffer,&pointID,&gridType,&curvilinearGridMetric[0],&curvilinearPointMetric[0]);
      FC_MODULE(metricops,pointeigenvectors,METRICOPS,POINTEIGENVECTORS)
        (&numDim,&uniformPointMetric[0],&uniformEigenVectors[0],&uniformMetricMags[0]);
      FC_MODULE(metricops,pointeigenvectors,METRICOPS,POINTEIGENVECTORS)
        (&numDim,&rectilinearPointMetric[0],&rectilinearEigenVectors[0],&rectilinearMetricMags[0]);
      FC_MODULE(metricops,pointeigenvectors,METRICOPS,POINTEIGENVECTORS)
        (&numDim,&curvilinearPointMetric[0],&curvilinearEigenVectors[0],&curvilinearMetricMags[0]);
     
      for(int iDim = 0;iDim < numDim;iDim++){
        size_t iOffset = iDim*numPointsBuffer;
        size_t iOffset2 = iDim*numDim*numPointsBuffer;
        size_t uniformOffset = iDim;
        size_t rectilinearOffset = iOffset + iPoint;
        int normDir = iDim+1;
        gridType = 1;

        FC_MODULE(metricops,boundarypointmetric,METRICOPS,BOUNDARYPOINTMETRIC)
          (&numDim,&normDir,&numPointsBuffer,&pointID,&gridType,&uniformGridMetric[0],&uniformBoundaryMetric[0]);
        gridType = 2;
        FC_MODULE(metricops,boundarypointmetric,METRICOPS,BOUNDARYPOINTMETRIC)
          (&numDim,&normDir,&numPointsBuffer,&pointID,&gridType,&rectilinearGridMetric[0],&rectilinearBoundaryMetric[0]);
        gridType = 3;
        FC_MODULE(metricops,boundarypointmetric,METRICOPS,BOUNDARYPOINTMETRIC)
          (&numDim,&normDir,&numPointsBuffer,&pointID,&gridType,&curvilinearGridMetric[0],&curvilinearBoundaryMetric[0]);
        
        double myError = std::abs(uniformMetricMags[iDim] - uniformGridMetric[uniformOffset]);
        if(myError > TOL){
          std::cout << "Uniform Metric Magnitude (exp,act) = (" << uniformGridMetric[uniformOffset]
                    << "," << uniformMetricMags[iDim] << "), Error("  << myError << ")" << std::endl;
          uniformResult[(numDim-2)+4] = 0;
        }
        myError = std::abs(rectilinearMetricMags[iDim]-rectilinearGridMetric[rectilinearOffset]);
        if(myError > TOL){
          std::cout << "Rectilinear Metric Magnitude (exp,act) = (" << rectilinearGridMetric[rectilinearOffset]
                    << "," << rectilinearMetricMags[iDim] << "), Error(" << myError << ")" << std::endl;
          rectilinearResult[(numDim-2)+4] = 0;
        }
        
        double curvilinearPointMagnitude = 0.0;
        for(int jDim = 0;jDim < numDim;jDim++){
          size_t jOffset = jDim*numPointsBuffer;
          size_t pointMetricOffset = iDim*numDim + jDim;
          size_t curvilinearOffset = iOffset2 + jOffset + iPoint;

          if(iDim == jDim){
            if(uniformPointMetric[pointMetricOffset] != uniformGridMetric[uniformOffset])
              uniformResult[numDim-2] = 0;
            if(std::abs(uniformEigenVectors[pointMetricOffset] - 1.0) > TOL)
              rectilinearResult[(numDim-2)+2] = 0;
            if(std::abs(rectilinearEigenVectors[pointMetricOffset] - 1.0) > TOL)
              uniformResult[(numDim-2)+2] = 0;
            if(rectilinearPointMetric[pointMetricOffset] != rectilinearGridMetric[rectilinearOffset]){
              //              rectilinearDiagError++;
              std::cout << "rectidiag(" << rectilinearPointMetric[pointMetricOffset]
                        << "," << rectilinearGridMetric[rectilinearOffset] << ")" 
                        << std::endl;
              rectilinearResult[numDim-2] = 0;
            }
          } else {
            if(std::abs(uniformEigenVectors[pointMetricOffset]) > TOL)
              uniformResult[(numDim-2)+2] = 0;
            if(std::abs(rectilinearEigenVectors[pointMetricOffset]) > TOL)
              rectilinearResult[(numDim-2)+2] = 0;
            if(uniformPointMetric[pointMetricOffset] > 0)
              uniformResult[numDim-2] = 0;
            if(rectilinearPointMetric[pointMetricOffset] > 0){
              //              rectilinearOffDiagError++;
              std::cout << "Rect off diag ";
              rectilinearResult[numDim-2] = 0;
            }
          }

          curvilinearPointMagnitude += curvilinearGridMetric[curvilinearOffset]*curvilinearGridMetric[curvilinearOffset];
          if(curvilinearPointMetric[pointMetricOffset] != curvilinearGridMetric[curvilinearOffset])
            curvilinearResult[numDim-2] = 0;

          // Check the boundary metric for the "idim" direction
          size_t boundaryRow = (iDim + jDim)%numDim;
          size_t boundaryOffset = boundaryRow*numDim;
          double *uniformStart = &uniformPointMetric[boundaryOffset];
          double *rectilinearStart = &rectilinearPointMetric[boundaryOffset];
          double *curvilinearStart = &curvilinearPointMetric[boundaryOffset];
          std::vector<double> uniformPointVector(uniformStart,uniformStart+numDim);
          std::vector<double> rectilinearPointVector(rectilinearStart,rectilinearStart+numDim);
          std::vector<double> curvilinearPointVector(curvilinearStart,curvilinearStart+numDim);
          std::vector<double> uniformBoundaryVector(&uniformBoundaryMetric[jDim*numDim],
                                                    &uniformBoundaryMetric[jDim*numDim]+numDim);
          std::vector<double> rectilinearBoundaryVector(&rectilinearBoundaryMetric[jDim*numDim],
                                                        &rectilinearBoundaryMetric[jDim*numDim]+numDim);
          std::vector<double> curvilinearBoundaryVector(&curvilinearBoundaryMetric[jDim*numDim],
                                                        &curvilinearBoundaryMetric[jDim*numDim]+numDim);
          if(uniformPointVector != uniformBoundaryVector)
            uniformResult[(numDim-2)+6] = 0;
          if(rectilinearPointVector != rectilinearBoundaryVector)
            rectilinearResult[(numDim+2)+6] = 0;
          if(curvilinearPointVector != curvilinearBoundaryVector)
            curvilinearResult[(numDim+2)+6] = 0;


        }

        curvilinearPointMagnitude = std::sqrt(curvilinearPointMagnitude);
        double myErr = std::abs(curvilinearMetricMags[iDim] - curvilinearPointMagnitude);
        if(myErr > TOL){
          std::cout << "Curvilinear Metric Magnitude (exp,act) = (" << curvilinearPointMagnitude
                    << "," << curvilinearMetricMags[iDim] << "), Error(" << myErr << ")" 
                    << std::endl;
          curvilinearResult[(numDim-2)+4] = 0;
        }

        for(int jDim = 0;jDim < numDim;jDim++){
          size_t jOffset = jDim*numPointsBuffer;
          size_t pointMetricOffset = iDim*numDim + jDim;
          size_t curvilinearOffset = iOffset2 + jOffset + iPoint;
          double expectedValue = curvilinearGridMetric[curvilinearOffset]/curvilinearPointMagnitude;
          myErr = std::abs(curvilinearEigenVectors[pointMetricOffset]-expectedValue);
          if(myErr > TOL){
            curvilinearResult[(numDim-2)+2] = 0;
            std::cout << "Curvilinear EigenVector (exp,act) = (" 
                      << curvilinearEigenVectors[pointMetricOffset]
                      << "," << expectedValue << "), Error(" 
                      << ")" << std::endl;
          }
        }        
      }
    } 
  }
  
  serialUnitResults.UpdateResult("Operators:MetricOps:UniformPointMetric2D",uniformResult[0]);
  serialUnitResults.UpdateResult("Operators:MetricOps:UniformPointMetric3D",uniformResult[1]);
  serialUnitResults.UpdateResult("Operators:MetricOps:RectilinearPointMetric2D",rectilinearResult[0]);
  serialUnitResults.UpdateResult("Operators:MetricOps:RectilinearPointMetric3D",rectilinearResult[1]);
  serialUnitResults.UpdateResult("Operators:MetricOps:CurvilinearPointMetric2D",curvilinearResult[0]);
  serialUnitResults.UpdateResult("Operators:MetricOps:CurvilinearPointMetric3D",curvilinearResult[1]);

  serialUnitResults.UpdateResult("Operators:MetricOps:UniformEigenVector2D",uniformResult[2]);
  serialUnitResults.UpdateResult("Operators:MetricOps:UniformEigenVector3D",uniformResult[3]);
  serialUnitResults.UpdateResult("Operators:MetricOps:RectilinearEigenVector2D",rectilinearResult[2]);
  serialUnitResults.UpdateResult("Operators:MetricOps:RectilinearEigenVector3D",rectilinearResult[3]);
  serialUnitResults.UpdateResult("Operators:MetricOps:CurvilinearEigenVector2D",curvilinearResult[2]);
  serialUnitResults.UpdateResult("Operators:MetricOps:CurvilinearEigenVector3D",curvilinearResult[3]);

  serialUnitResults.UpdateResult("Operators:MetricOps:UniformMetricMags2D",uniformResult[4]);
  serialUnitResults.UpdateResult("Operators:MetricOps:UniformMetricMags3D",uniformResult[5]);
  serialUnitResults.UpdateResult("Operators:MetricOps:RectilinearMetricMags2D",rectilinearResult[4]);
  serialUnitResults.UpdateResult("Operators:MetricOps:RectilinearMetricMags3D",rectilinearResult[5]);
  serialUnitResults.UpdateResult("Operators:MetricOps:CurvilinearMetricMags2D",curvilinearResult[4]);
  serialUnitResults.UpdateResult("Operators:MetricOps:CurvilinearMetricMags3D",curvilinearResult[5]);

  serialUnitResults.UpdateResult("Operators:MetricOps:UniformBoundaryMetric2D",uniformResult[6]);
  serialUnitResults.UpdateResult("Operators:MetricOps:UniformBoundaryMetric3D",uniformResult[7]);
  serialUnitResults.UpdateResult("Operators:MetricOps:RectilinearBoundaryMetric2D",rectilinearResult[6]);
  serialUnitResults.UpdateResult("Operators:MetricOps:RectilinearBoundaryMetric3D",rectilinearResult[7]);
  serialUnitResults.UpdateResult("Operators:MetricOps:CurvilinearBoundaryMetric2D",curvilinearResult[6]);
  serialUnitResults.UpdateResult("Operators:MetricOps:CurvilinearBoundaryMetric3D",curvilinearResult[7]);

  return;
}
