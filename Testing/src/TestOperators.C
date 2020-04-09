#include "Testing.H"
#include "Simulation.H"
#include "OperatorKernels.H"
#include "MetricKernels.H"
#include "OperatorTypes.H"




void TestOperators(ix::test::results &serialUnitResults)
{
  int numDim = 3;

  std::vector<size_t> numX(numDim,10);
  size_t numPoints = 1;
  for(int iDim = 0;iDim  < numDim;iDim++){
    numPoints *= numX[iDim];
  }
  
  pcpp::IndexIntervalType opInterval;
  opInterval.InitSimple(numX);
  std::vector<size_t> flatInterval;
  opInterval.Flatten(flatInterval);
  // Forticate the interval
  for(int iDim = 0;iDim < numDim;iDim++){
    flatInterval[2*iDim]++;
    flatInterval[2*iDim+1]++;
  }
  size_t iStart   = opInterval[0].first;
  size_t iEnd     = opInterval[0].second;
  size_t jStart   = opInterval[1].first;
  size_t jEnd     = opInterval[1].second;
  size_t kStart   = opInterval[2].first;
  size_t kEnd     = opInterval[2].second;
  size_t numPlane = numX[0]*numX[1];
  
  double a = 1.5;
  double x = 2.0;
  double y = 1.0;
  double z = 0.0;
  double b = 2.0;
  double w = 3.0;
  double v = 4.0;

  std::vector<double> W0(numPoints,w);
  std::vector<double> X0(numPoints,x);
  std::vector<double> Y0(numPoints,y);
  std::vector<double> Z0(numPoints,z);
  std::vector<double> X3(3*numPoints,0);
  std::vector<double> Y3(3*numPoints,0);
  std::vector<double> X3Len(numPoints,0);
  std::vector<double> Y3LenNeg(numPoints,0);
  std::vector<double> Y3Len(numPoints,0);
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    X0[iPoint] = x*iPoint;
    Y0[iPoint] = y*iPoint;
    W0[iPoint] = w*iPoint;
    X3[iPoint] = x*iPoint;
    X3[iPoint+numPoints]   = 2*x*iPoint;
    X3[iPoint+2*numPoints] = 3*x*iPoint;
    Y3[iPoint] = y*iPoint;
    Y3[iPoint+numPoints]   = 2*y*iPoint;
    Y3[iPoint+2*numPoints] = 3*y*iPoint;
    Y3LenNeg[iPoint] = Y3[iPoint]*Y3[iPoint] +
      Y3[iPoint+numPoints]*Y3[iPoint+numPoints] +
      Y3[iPoint+2*numPoints]*Y3[iPoint+2*numPoints];
    Y3LenNeg[iPoint] = -1.0*std::sqrt(Y3LenNeg[iPoint]);
    Y3Len[iPoint] = -1.0*Y3LenNeg[iPoint];

    X3Len[iPoint] = X3[iPoint]*X3[iPoint] +
      X3[iPoint+numPoints]*X3[iPoint+numPoints] +
      X3[iPoint+2*numPoints]*X3[iPoint+2*numPoints];
    X3Len[iPoint] = std::sqrt(X3Len[iPoint]);
  }
  
  std::vector<double> Yprime(Y0);
  std::vector<double> V(numPoints,v);
  std::vector<double> W(numPoints,w);
  std::vector<double> X(numPoints,x);
  std::vector<double> Y(numPoints,y);
  std::vector<double> Z(numPoints,z);
  std::vector<double> A(numPoints,a);
  std::vector<double> B(numPoints,b);
  
  Y = Y0;
  X = X0;
  
  size_t *bufferSize = &numX[0];
  size_t *bufferInterval = &flatInterval[0];
  int *numComponents = &numDim;
  
  // Test Z = X (dot) Y 
  FC_MODULE(operators,zxdoty,OPERATORS,ZXDOTY)
    (&numDim,&numPoints,bufferSize,bufferInterval,
     numComponents,&X3[0],&Y3[0],&Z[0]);
  bool zxdoty = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Z[xyzIndex] != (X3[xyzIndex]*Y3[xyzIndex]) +
           (X3[xyzIndex+numPoints]*Y3[xyzIndex+numPoints]) +
           (X3[xyzIndex+2*numPoints]*Y3[xyzIndex+2*numPoints]))
          zxdoty = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ZXDOTY",zxdoty);
  
  
  Z = Z0;
  
  // veclen = SQRT(V[0]**2 + .... V[n]**2) 
  bool veclen = true;
  std::vector<double> testLen(numPoints,0);
  FC_MODULE(operators,veclen,OPERATORS,VECLEN)
    (&numDim,&numPoints,bufferSize,bufferInterval,
     numComponents,&X3[0],&testLen[0]);
  if(testLen != X3Len){
    veclen = false;
    for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
      if(testLen[iPoint] != X3Len[iPoint])
        std::cout << "VecLenTest (Value,Expected)[" << iPoint << "] = ("
                  << testLen[iPoint] << "," << X3Len[iPoint] 
                  << ")" << std::endl;
    }
  }
  serialUnitResults.UpdateResult("Operators:VECLEN",veclen);

  // Y = ABS(X)
  bool assignyabsx = true;
  FC_MODULE(operators,assignmentyabsx,OPERATORS,ASSIGNMENTYABSX)
    (&numDim,&numPoints,bufferSize,bufferInterval,
     &Y3LenNeg[0],&testLen[0]);
  if(testLen != Y3Len)
    assignyabsx = false;
  serialUnitResults.UpdateResult("Operators:ASSIGNMENTYABSX",assignyabsx);

  // Y = XY
  FC_MODULE(operators,yxy,OPERATORS,YXY)
    (&numDim,&numPoints,bufferSize,bufferInterval,
     &X[0], &Y[0]);
  bool yxy = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Y[xyzIndex] != (X[xyzIndex]*Y0[xyzIndex]))
          yxy = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:YXY",yxy);
  Y = Y0;

  // Z = aW+XY
  FC_MODULE(operators,zawpxy,OPERATORS,ZAWPXY)
    (&numDim,&numPoints,bufferSize,bufferInterval,
     &a,&W[0],&X[0],&Y[0],&Z[0]);
  bool zawpxy = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Z[xyzIndex] != (a*W[xyzIndex] + X[xyzIndex]*Y[xyzIndex]))
          zawpxy = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ZAWPXY",zawpxy);
  Z = Z0;

  // Z = VW + XY
  FC_MODULE(operators,zvwpxy,OPERATORS,ZVWPXY)
    (&numDim,&numPoints,bufferSize,bufferInterval,
     &V[0],&W[0],&X[0],&Y[0],&Z[0]);
  bool zvwpxy = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Z[xyzIndex] != (V[xyzIndex]*W[xyzIndex] + X[xyzIndex]*Y[xyzIndex]))
          zvwpxy = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ZVWPXY",zvwpxy);
  Z = Z0;
  
  std::vector<double> inMatrix(9*numPoints,0);
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    for(int iComp = 0;iComp < 9;iComp++){
      inMatrix[iPoint+iComp*numPoints] = (iComp+1)*9*(iPoint+1);
    }
  }
  std::vector<double> matrixDeterminant(numPoints,0);

  bool det3 = true;
  // Test 3x3 determinant (i.e. at each point)
  FC_MODULE(operators,determinant3x3,OPERATORS,DETERMINANT3X3)
    (&numPoints,bufferSize,bufferInterval,
     &inMatrix[0],&matrixDeterminant[0]);
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    if(matrixDeterminant[iPoint] != inMatrix[iPoint]*
       (inMatrix[iPoint+4*numPoints]*inMatrix[iPoint+8*numPoints] -
        inMatrix[iPoint+5*numPoints]*inMatrix[iPoint+7*numPoints]) -
       inMatrix[iPoint+3*numPoints]*
       (inMatrix[iPoint+numPoints]*inMatrix[iPoint+8*numPoints] -
        inMatrix[iPoint+2*numPoints]*inMatrix[iPoint+7*numPoints]) +
       inMatrix[iPoint+6*numPoints]*
       (inMatrix[iPoint+numPoints]*inMatrix[iPoint+5*numPoints] -
        inMatrix[iPoint+2*numPoints]*inMatrix[iPoint+4*numPoints]))
      det3 = false;
  }
  serialUnitResults.UpdateResult("Operators:Determinant3x3",det3);
  
  size_t numPoints2 = bufferSize[0]*bufferSize[1];
  std::vector<double> inMatrix2(4*numPoints2,0);
  for(size_t iPoint = 0;iPoint < numPoints2;iPoint++){
    for(int iComp = 0;iComp < 4;iComp++){
      inMatrix2[iPoint+iComp*numPoints2] = (iComp+1)*4*(iPoint+1);
    }
  }
  bool det2 = true;
  // Test 2x2 determinant (i.e. at each point)
  FC_MODULE(operators,determinant2x2,OPERATORS,DETERMINANT2X2)
    (&numPoints2,bufferSize,bufferInterval,
     &inMatrix2[0],&matrixDeterminant[0]);
  for(size_t iPoint = 0;iPoint < numPoints2;iPoint++){
    double expected = (inMatrix2[iPoint]*inMatrix2[iPoint+3*numPoints2]) -
      (inMatrix2[iPoint+numPoints2]*inMatrix2[iPoint+2*numPoints2]);
    if(matrixDeterminant[iPoint] != expected){
      det2 = false;
    }
  }
  serialUnitResults.UpdateResult("Operators:Determinant2x2",det2);


  //   Test buf1*buf4 - buf2*buf3 + buf7*(buf5 - buf6)
  bool metricsum = true;
  FC_MODULE(operators,metricsum4,OPERATORS,METRICSUM4)
    (&numDim,&numPoints,bufferSize,bufferInterval,
     &A[0],&B[0],&V[0],&W[0],&X[0],&Y[0],&Yprime[0],&matrixDeterminant[0]);
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    if(matrixDeterminant[iPoint] != A[iPoint]*W[iPoint] -
       B[iPoint]*V[iPoint] + Yprime[iPoint]*(X[iPoint]-Y[iPoint]))
      metricsum = false;
  }
  serialUnitResults.UpdateResult("Operators:MetricSum4",metricsum);

  // Test YAXPY (Y = aX + Y) [scalar a]
  FC_MODULE(operators,yaxpy,OPERATORS,YAXPY)
    (&numDim,&numPoints,&numX[0],&flatInterval[0],&a,&X[0],&Y[0]);
  
  bool yaxpy = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Y[xyzIndex] != a*X[xyzIndex]+Y0[xyzIndex])
          yaxpy = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:YAXPY",yaxpy);

  Y = Y0;
  Yprime = Y0;
  
  // Test Y = W*X + Y  
  FC_MODULE(operators,ywxpy,OPERATORS,YWXPY)
    (&numDim,&numPoints,&numX[0],&flatInterval[0],&W[0],&X[0],&Yprime[0]);
  bool ywxpy = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Yprime[xyzIndex] != W[xyzIndex]*X[xyzIndex]+Y0[xyzIndex])
          ywxpy = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:YWXPY",ywxpy);

  Yprime = Y0;

  // Test Z = aX + bY [scalar a,b]
  FC_MODULE(operators,zaxpby,OPERATORS,ZAXPBY)
    (&numDim,&numPoints,&numX[0],&flatInterval[0],&a,&b,&X[0],&Y[0],&Z[0]);
  bool zaxpby = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Z[xyzIndex] != a*X[xyzIndex]+b*Y[xyzIndex])
          zaxpby = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ZAXPBY",zaxpby);
  Z = Z0;
  Yprime = Y0;

  // Test Y = aX + bY [scalar a,b]
  FC_MODULE(operators,yaxpby,OPERATORS,YAXPBY)
    (&numDim,&numPoints,&numX[0],&flatInterval[0],&a,&b,&X[0],&Yprime[0]);
  bool yaxpby = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Yprime[xyzIndex] != a*X[xyzIndex]+b*Y0[xyzIndex])
          yaxpby = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:YAXPBY",yaxpby);
  Yprime = Y0;

  // Test ZAXPY (Z = aX + Y) [scalar a]
  FC_MODULE(operators,zaxpy,OPERATORS,ZAXPY)(&numDim,&numPoints,&numX[0],&flatInterval[0],&a,&X[0],&Y[0],&Z[0]);
  bool zaxpy = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Z[xyzIndex] != a*X[xyzIndex]+Y[xyzIndex])
          zaxpy = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ZAXPY",zaxpy);
  
  // Test YAX (Y = a*X) [scalar a]
  FC_MODULE(operators,yax,OPERATORS,YAX)(&numDim,&numPoints,&numX[0],&flatInterval[0],&a,&X[0],&Y[0]);
  bool yax = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Y[xyzIndex] != a*X[xyzIndex])
          yax = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:YAX",yax);
  Y = Y0;
  
  // Test XAX (X = a*X) [scalar a]
  FC_MODULE(operators,xax,OPERATORS,XAX)(&numDim,&numPoints,&numX[0],&flatInterval[0],&a,&X[0]);
  bool xax = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(X[xyzIndex] != a*X0[xyzIndex])
          xax = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:XAX",xax);  
  X = X0;

  // Test ZXY (Z = X*Y)
  FC_MODULE(operators,zxy,OPERATORS,ZXY)(&numDim,&numPoints,&numX[0],&flatInterval[0],&X[0],&Y[0],&Z[0]);
  bool zxy = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Z[xyzIndex] != X[xyzIndex]*Y[xyzIndex])
          zxy = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ZXY",zxy);

  // Test ZWXPY (Z = W*X+Y)
  FC_MODULE(operators,zwxpy,OPERATORS,ZWXPY)(&numDim,&numPoints,&numX[0],&flatInterval[0],&W[0],&X[0],&Y[0],&Z[0]);
  bool zwxpy = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Z[xyzIndex] != W[xyzIndex]*X[xyzIndex]+Y[xyzIndex])
          zwxpy = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ZWXPY",zwxpy);


  // Test ZWMXPY (Z = W*(X+Y))
  FC_MODULE(operators,zwmxpy,OPERATORS,ZWMXPY)(&numDim,&numPoints,&numX[0],&flatInterval[0],&W[0],&X[0],&Y[0],&Z[0]);
  bool zwmxpy = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Z[xyzIndex] != W[xyzIndex]*(X[xyzIndex]+Y[xyzIndex])){
          zwmxpy = false;
          std::cout << "Z[" << xyzIndex << "] = " << Z[xyzIndex] << std::endl;
        }
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ZWMXPY",zwmxpy);
  
  
  // Test ASSIGNMENTYX (Y = X) [note the order of arguments to this operator, the output var is *always* last]
  FC_MODULE(operators,assignmentyx,OPERATORS,ASSIGNMENTYX)(&numDim,&numPoints,&numX[0],&flatInterval[0],&X[0],&Y[0]);
  bool yx = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Y[xyzIndex] != X[xyzIndex])
          yx = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ASSIGNMENTYX",yx);


  // Test ASSIGNMENTXA (X = a) [scalar assignment]
  X = X0;
  FC_MODULE(operators,assignmentxa,OPERATORS,ASSIGNMENTXA)(&numDim,&numPoints,&numX[0],&flatInterval[0],&a,&X[0]);
  bool xa = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(X[xyzIndex] != a)
          xa = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ASSIGNMENTXA",xa);

  X = X0;

  // Test YAXM1 (Y = a/X) [vec recip]
  FC_MODULE(operators,yaxm1,OPERATORS,YAXM1)(&numDim,&numPoints,&numX[0],&flatInterval[0],&a,&X[0],&Y[0]);
  bool xaym1 = true;
  for(size_t iZ = kStart;iZ <= kEnd;iZ++){
    size_t zIndex = iZ*numPlane;
    for(size_t iY = jStart;iY <= jEnd;iY++){
      size_t yzIndex = zIndex + iY*numX[0];
      for(size_t iX = iStart;iX <= iEnd;iX++){
        size_t xyzIndex = yzIndex + iX;
        if(Y[xyzIndex] != a/X[xyzIndex])
          xaym1 = false;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:XAYM1",xaym1);


  // Test Dilatation for all the different sorts of grid metrics
  bool graddiv = true;
  int gridType = simulation::grid::UNIRECT; // Uniform Rectangular
  std::vector<double> gridMetric0(3,1.0);
  gridMetric0[1] = 2.0;
  gridMetric0[2] = 4.0;
  std::vector<double> gridJacobian0(1,1.0/8.0);
  std::vector<double> gradV(9*numPoints,1.0);
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    for(int iComp = 0;iComp < 9;iComp++){
      gradV[iPoint + iComp*numPoints] = std::pow(2.0,static_cast<double>(iComp))*(iPoint+1);
    }
  }
  std::vector<double> vDiv(numPoints,0.0);

  FC_MODULE(metricops,ijkgradtoxyzdiv,METRICOPS,IJKGRADTOXYZDIV)
    (&numDim,&numPoints,&numX[0],&flatInterval[0],&gridType,&gridJacobian0[0], 
     &gridMetric0[0],&gradV[0],&vDiv[0]);
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    double expected = gridMetric0[0]*gradV[iPoint] + 
      gridMetric0[1]*gradV[iPoint+4*numPoints] + gridMetric0[2]*gradV[iPoint+8*numPoints];
    expected *= gridJacobian0[0];
    if(vDiv[iPoint] != expected){
      graddiv = false;
    }
  }
  serialUnitResults.UpdateResult("Operators:MetricOps:Dilatation:Uniform",graddiv);
  vDiv.resize(0);
  vDiv.resize(numPoints,0.0);

  graddiv = true;
  gridType++; // Rectilinear
  std::vector<double> gridMetric1(3*numPoints,1.0);
  std::vector<double> gridJacobian1(numPoints,1.0);
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    for(int iComp = 0;iComp < 3;iComp++){
      gridMetric1[iPoint+iComp*numPoints] = std::pow(2.0,static_cast<double>(iComp))*(iPoint+1);
      gridJacobian1[iPoint] *= 1.0/gridMetric1[iPoint+iComp*numPoints];
    }
  }

  FC_MODULE(metricops,ijkgradtoxyzdiv,METRICOPS,IJKGRADTOXYZDIV)
    (&numDim,&numPoints,&numX[0],&flatInterval[0],&gridType,&gridJacobian1[0], 
     &gridMetric1[0],&gradV[0],&vDiv[0]);

  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    double expected = gridMetric1[iPoint]*gradV[iPoint] + 
      gridMetric1[iPoint+numPoints]*gradV[iPoint+4*numPoints] + 
      gridMetric1[iPoint+2*numPoints]*gradV[iPoint+8*numPoints];
    expected *= gridJacobian1[iPoint];
    if(vDiv[iPoint] != expected){
      graddiv = false;
    }
  }
  serialUnitResults.UpdateResult("Operators:MetricOps:Dilatation:Rectilinear",graddiv);
  vDiv.resize(0);
  vDiv.resize(numPoints,0.0);

  gridType++; // Curvilinear
  graddiv = true;
  std::vector<double> gridMetric2(9*numPoints,1.0);
  std::vector<double> gridJacobian2(numPoints,1.0);
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    for(int iComp = 0;iComp < 9;iComp++){
      gridMetric2[iPoint+iComp*numPoints] = std::pow(2.0,static_cast<double>(iComp))*(iPoint+1);
      gridJacobian2[iPoint] *= 1.0/gridMetric2[iPoint+iComp*numPoints];
    }
  }

  FC_MODULE(metricops,ijkgradtoxyzdiv,METRICOPS,IJKGRADTOXYZDIV)
    (&numDim,&numPoints,&numX[0],&flatInterval[0],&gridType,&gridJacobian2[0], 
     &gridMetric2[0],&gradV[0],&vDiv[0]);

  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    double expecteddudx = gridMetric2[iPoint]*gradV[iPoint] + 
      gridMetric2[iPoint+3*numPoints]*gradV[iPoint+numPoints] +
      gridMetric2[iPoint+6*numPoints]*gradV[iPoint+2*numPoints];
    double expecteddvdy = gridMetric2[iPoint+numPoints]*gradV[iPoint+3*numPoints] +
      gridMetric2[iPoint+4*numPoints]*gradV[iPoint+4*numPoints] +
      gridMetric2[iPoint+7*numPoints]*gradV[iPoint+5*numPoints];
    double expecteddwdz = gridMetric2[iPoint+2*numPoints]*gradV[iPoint+6*numPoints] +
      gridMetric2[iPoint+5*numPoints]*gradV[iPoint+7*numPoints] +
      gridMetric2[iPoint+8*numPoints]*gradV[iPoint+8*numPoints];
    double expected = expecteddudx + expecteddvdy + expecteddwdz;
    expected *= gridJacobian2[iPoint];
    if(vDiv[iPoint] != expected){
      graddiv = false;
    }
  }
  serialUnitResults.UpdateResult("Operators:MetricOps:Dilatation:Curvilinear",graddiv);
  
  gridType = simulation::grid::UNIRECT;
  //   std::vector<double> gridMetric0(3,1.0);
  //   gridMetric0[1] = 2.0;
  //   gridMetric0[2] = 4.0;
  //   std::vector<double> gridJacobian0(1,1.0/8.0);
  bool alphawt = true;
  std::vector<double> alphaWeight(numPoints,0.0);
  for(int alphaDir = 1;alphaDir <= 3;alphaDir++){
    FC_MODULE(metricops,alphaweight,METRICOPS,ALPHAWEIGHT)
      (&numDim,&numPoints,&numX[0],&flatInterval[0],&gridType,&gridMetric0[0],
       &alphaDir,&alphaWeight[0]);
    for(size_t iPoint = 0;iPoint < numPoints;iPoint++)
      if(alphaWeight[iPoint] != gridMetric0[alphaDir-1])
        alphawt = false;
  }
  serialUnitResults.UpdateResult("Operators:MetricOps:AlphaWeight:Uniform",alphawt);
 
  gridType++;
  alphawt = true;
  for(int alphaDir = 1;alphaDir <= 3;alphaDir++){
    FC_MODULE(metricops,alphaweight,METRICOPS,ALPHAWEIGHT)
      (&numDim,&numPoints,&numX[0],&flatInterval[0],&gridType,&gridMetric1[0],
       &alphaDir,&alphaWeight[0]);
    for(size_t iPoint = 0;iPoint < numPoints;iPoint++)
      if(alphaWeight[iPoint] != gridMetric1[(alphaDir-1)*numPoints+iPoint])
        alphawt = false;
  }
  serialUnitResults.UpdateResult("Operators:MetricOps:AlphaWeight:Rectilinear",alphawt);
  
  gridType++;
  alphawt = true;
  for(int alphaDir = 1;alphaDir <= 3;alphaDir++){
    std::vector<double> metricLen(numPoints,0.0);
    FC_MODULE(operators,veclen,OPERATORS,VECLEN)
      (&numDim,&numPoints,bufferSize,bufferInterval,
       &numDim,&gridMetric2[(alphaDir-1)*numDim*numPoints],&metricLen[0]);
    FC_MODULE(metricops,alphaweight,METRICOPS,ALPHAWEIGHT)
      (&numDim,&numPoints,&numX[0],&flatInterval[0],&gridType,&gridMetric2[0],
       &alphaDir,&alphaWeight[0]);
    if(alphaWeight != metricLen)
      alphawt = false;
  }
  serialUnitResults.UpdateResult("Operators:MetricOps:AlphaWeight:Curvilinear",alphawt);

  
  

}

size_t factorial(size_t n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
//void TestApplySingleStencil(ix::test::results &serialUnitResults);

// int SBPOperatorTest1(plascom2::stencilset &inOperator,int interiorOrder,int numDim,
//                      int numComponents,int numTrials,std::vector<bool> &testResults);
/// Tests boundary stencil setting around holes
void TestHoleDetection(ix::test::results &serialUnitResults)
{
  operators::sbp_operator_t order12Operator;
  operators::sbp::Initialize(order12Operator,2);


  // Create a test geometry
  std::vector<bool> test1Pass(3,true);
  bool test1Result = true;
  std::vector<bool> test2Pass(3,true);
  bool test2Result = true;
  std::vector<bool> test3Pass(3,true);
  bool test3Result = true;
  std::vector<bool> test4Pass(3,true);
  bool test4Result = true;
  
  for(int numDim = 1;numDim <= 3;numDim++){

    std::cout << "Testing hole detection in " << numDim << " dimensions." 
              << std::endl;
    int dimId = numDim-1;

    size_t numPoints       = 10;
    size_t numPointsBuffer = 1;

    std::vector<size_t> bufferSizes(numDim,0);
    size_t *numPointsDim = &bufferSizes[0];
    size_t *opInterval = new size_t [2*numDim];
    std::vector<int> periodicDirs(numDim,0);
    
    for(int iDim = 0;iDim < numDim;iDim++){

      numPointsDim[iDim]   = numPoints;
      numPointsBuffer     *= numPointsDim[iDim];
      opInterval[2*iDim]   = 0;
      opInterval[2*iDim+1] = numPointsDim[iDim]-1;

    }

    std::vector<size_t> gridSizes(bufferSizes);
    std::vector<int> gridMask;

    pcpp::IndexIntervalType partitionInterval(opInterval,numDim);
    pcpp::IndexIntervalType partitionBufferInterval(opInterval,numDim);
    pcpp::IndexIntervalType bufferInterval;
    bufferInterval.InitSimple(bufferSizes);

    int holeMask      = 1;

    for(int boundaryDepth = 1;boundaryDepth <= 5;boundaryDepth++){

      // Reset the mask for each test
      gridMask.resize(0);
      gridMask.resize(numPointsBuffer,0);
      int *iMask        = &gridMask[0];
      
      std::cout << "Testing with boundaryDepth = " << boundaryDepth << std::endl;

      int numStencils = 2*(boundaryDepth+1);

      // Create stencil connectivity
      std::vector<int> stencilConnectivity(numDim*numPointsBuffer,0);
      int *stencilID    = &stencilConnectivity[0];

      operators::sbp::CreateStencilConnectivity(numDim,numPointsDim,opInterval,
                                                          boundaryDepth,&periodicDirs[0],stencilID,false);
      
      std::vector<int> expectedStencilConn(stencilConnectivity);
      
      // Output to make sure stencil connectivity is as it should be
      if(numDim == 1 || numDim == 2){
        std::cout << "Default/initial stencil connectivity:" << std::endl;
        for(int iDim = 0;iDim < numDim;iDim++){
          std::cout << "Dimension-" << iDim+1;
          size_t pointOffset = iDim*numPointsBuffer;
          for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
            if(!(iPoint%numPoints)) std::cout << std::endl;
            std::cout << stencilID[pointOffset+iPoint] << " ";
          }
          std::cout << std::endl;
        }
      }
      
      std::cout << "Calling DetectHoles on unholy grid." << std::endl;

      if(operators::sbp::DetectHoles(bufferSizes,partitionInterval,partitionBufferInterval,
                                               partitionInterval,boundaryDepth,holeMask,iMask,stencilID)){
        test1Pass[dimId] = false;
        std::cout << "DetectHoles failed to run." << std::endl;
      }    
      
      if(stencilConnectivity != expectedStencilConn){
        test1Pass[dimId] = false;
        std::cout << "DetectHoles produced wrong result." << std::endl;
        if(numDim == 1 || numDim == 2 || numDim == 3){
          std::cout << "Conn after hole detection: " << std::endl;
          for(int iDim = 0;iDim < numDim;iDim++){
            int *expectedConn = &stencilConnectivity[iDim*numPointsBuffer];
            std::cout << "Dimension-" << iDim+1 << ":";
            for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
              if(!(iPoint%(numPoints*numPoints))) {
                std::cout << std::endl << " --------------- " << std::endl;
              } else if(!(iPoint%numPoints)) std::cout << std::endl;
              std::cout << expectedConn[iPoint] << " ";
            }
            std::cout << std::endl;
          }
        }        
      } else {
        std::cout << "Unholy test passed." << std::endl;
      }
      
      test1Result = test1Result && test1Pass[dimId];
      
      // Now add a "hole" on the edge of the domain 
      // [X000000X]  in each dimension (X=hole)
      // Select each domain surface and add a hole
      // Left boundary i = 0, for all of the rest
      expectedStencilConn.resize(0);
      expectedStencilConn.resize(numDim*numPointsBuffer,1);
      
      for(int iDim = 0;iDim < numDim;iDim++){
        int *expectedConn = &expectedStencilConn[iDim*numPointsBuffer];
        pcpp::IndexIntervalType leftBoundaryInterval(partitionInterval);
        if(leftBoundaryInterval[iDim].first == 0) { // then I have a left boundary
          leftBoundaryInterval[iDim].second = 0;
          // set hole bits in mask
          pcpp::IndexIntervalType boundaryOpInterval;
          leftBoundaryInterval.RelativeTranslation(partitionInterval,
                                                   partitionBufferInterval,
                                                   boundaryOpInterval);
          if(boundaryOpInterval.NNodes() > 0) {
            
            // Set up left boundary
            mask::SetMask(bufferSizes,boundaryOpInterval,holeMask,iMask);
            
            // Set up expected stencilConn
            std::vector<size_t> nodeList;
            bufferInterval.GetFlatIndices(boundaryOpInterval,nodeList);
            for(int iPoint = 0;iPoint < nodeList.size();iPoint++)
              expectedConn[nodeList[iPoint]] = numStencils;
            for(int iBoundary = 0;iBoundary < boundaryDepth;iBoundary++){
              nodeList.resize(0);
              boundaryOpInterval[iDim].first++;
              boundaryOpInterval[iDim].second++;
              bufferInterval.GetFlatIndices(boundaryOpInterval,nodeList);
              for(int iPoint = 0;iPoint < nodeList.size();iPoint++)
                expectedConn[nodeList[iPoint]] = 2 + iBoundary;
            }
          }
        }
        pcpp::IndexIntervalType rightBoundaryInterval(partitionInterval);
        if(rightBoundaryInterval[iDim].second == (gridSizes[iDim]-1)){
          rightBoundaryInterval[iDim].first = gridSizes[iDim]-1;
          pcpp::IndexIntervalType boundaryOpInterval;
          rightBoundaryInterval.RelativeTranslation(partitionInterval,
                                                    partitionBufferInterval,
                                                    boundaryOpInterval);
          if(boundaryOpInterval.NNodes() > 0) {
            
            mask::SetMask(bufferSizes,boundaryOpInterval,holeMask,iMask);
            
            // Set up expected stencilConn
            std::vector<size_t> nodeList;
            bufferInterval.GetFlatIndices(boundaryOpInterval,nodeList);
            for(int iPoint = 0;iPoint < nodeList.size();iPoint++)
              expectedConn[nodeList[iPoint]] = numStencils;
            for(int iBoundary = 0;iBoundary < boundaryDepth;iBoundary++){
              nodeList.resize(0);
              boundaryOpInterval[iDim].second--;
              boundaryOpInterval[iDim].first--;
              bufferInterval.GetFlatIndices(boundaryOpInterval,nodeList);
              for(int iPoint = 0;iPoint < nodeList.size();iPoint++)
                expectedConn[nodeList[iPoint]] = 2 + boundaryDepth + iBoundary;
            } 
          }        
        }
      }
      if(numDim == 1 || numDim == 2 || numDim == 3){
        std::cout << "Setting holy border Mask: " << std::endl;
        for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
          if(!(iPoint%(numPoints*numPoints))) {
            std::cout << std::endl << "-----------" << std::endl;
          } else  if(!(iPoint%numPoints)) std::cout << std::endl;
          std::cout << iMask[iPoint] << " ";
          if(iMask[iPoint]&holeMask){
            for(int iDim = 0;iDim < numDim;iDim++){
              expectedStencilConn[iDim*numPointsBuffer+iPoint] = numStencils;
            }
          }
        }
        std::cout << std::endl;
        std::cout << "Expected Conn: " << std::endl;
        for(int iDim = 0;iDim < numDim;iDim++){
          int *expectedConn = &expectedStencilConn[iDim*numPointsBuffer];
          std::cout << "Dimension-" << iDim+1 << ":";
          for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
            if(!(iPoint%(numPoints*numPoints))) {
              std::cout << std::endl << "-----------" << std::endl;
            } else if(!(iPoint%numPoints)) std::cout << std::endl;
            std::cout << expectedConn[iPoint] << " ";
          }
          std::cout << std::endl;
        }
      }
      
      int returnCode = operators::sbp::DetectHoles(bufferSizes,partitionInterval,partitionBufferInterval,
                                                             partitionInterval,boundaryDepth,holeMask,iMask,stencilID);

      if(numStencils > numPoints){ // detect holes should fail
        if(returnCode == 0){
          test3Pass[dimId] = false;
        } else {
          std::cout << "DetectHoles detected too much holiness." << std::endl;
        }
      } else {
        if(returnCode > 0){
          test2Pass[dimId] = false;
          std::cout << "DetectHoles found spurious holiness." << std::endl;
        } else {
          if(stencilConnectivity != expectedStencilConn){
            test2Pass[dimId] = false;
            std::cout << "Test failed." << std::endl;
          } else {
            std::cout << "Holy border test passed." << std::endl;
          }
          if(numDim == 1 || numDim == 2 || numDim == 3){
            std::cout << "Conn after hole detection: " << std::endl;
            for(int iDim = 0;iDim < numDim;iDim++){
              int *expectedConn = &stencilConnectivity[iDim*numPointsBuffer];
              std::cout << "Dimension-" << iDim+1 << ":";
              for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
                if(!(iPoint%(numPoints*numPoints))) {
                  std::cout << std::endl << " --------------- " << std::endl;
                } else if(!(iPoint%numPoints)) std::cout << std::endl;
                std::cout << expectedConn[iPoint] << " ";
              }
              std::cout << std::endl;
            }
          }
        }
      }
      
      test2Result = test2Result && test2Pass[dimId];
      test3Result = test3Result && test3Pass[dimId];

      // Do a convex hole in the middle of the mesh
      // Reset the mask for each test
      gridMask.resize(0);
      gridMask.resize(numPointsBuffer,0);
      iMask        = &gridMask[0];

      // Create stencil connectivity
      stencilConnectivity.resize(0);
      stencilConnectivity.resize(numDim*numPointsBuffer,0);
      operators::sbp::CreateStencilConnectivity(numDim,numPointsDim,opInterval,
                                                          boundaryDepth,&periodicDirs[0],stencilID,false);
      expectedStencilConn = stencilConnectivity;
      
      // Now add a "hole" in the middle of the domain 
      // [000XXX000]  in each dimension (X=hole)
      expectedStencilConn.resize(0);
      expectedStencilConn.resize(numDim*numPointsBuffer,1);

      pcpp::IndexIntervalType holeInterval;
      holeInterval.InitSimple(gridSizes);
      pcpp::IndexIntervalType holeBoundaryZone(holeInterval);

      for(int iDim = 0;iDim < numDim;iDim++){
        size_t middleIndex = gridSizes[iDim]/2;
        holeInterval[iDim].first = middleIndex-1;
        holeInterval[iDim].second = middleIndex;
        if(boundaryDepth > holeInterval[iDim].first)
          holeBoundaryZone[iDim].first = 0;
        else
          holeBoundaryZone[iDim].first = holeInterval[iDim].first - boundaryDepth;
      }

      pcpp::IndexIntervalType myHoleInterval;
      partitionInterval.Overlap(holeInterval,myHoleInterval);

      if(myHoleInterval.NNodes() > 0){
          pcpp::IndexIntervalType holeOpInterval;
          myHoleInterval.RelativeTranslation(partitionInterval,
                                             partitionBufferInterval,
                                             holeOpInterval);
          mask::SetMask(bufferSizes,holeOpInterval,holeMask,iMask);
      }
      
      if(numDim == 1 || numDim == 2 || numDim == 3){
        std::cout << "Setting central holy border Mask: " << std::endl;
        for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
          if(!(iPoint%(numPoints*numPoints))) {
            std::cout << std::endl << "-----------" << std::endl;
          } else  if(!(iPoint%numPoints)) std::cout << std::endl;
          std::cout << iMask[iPoint] << " ";
        }
        std::cout << std::endl;
      }

     returnCode = operators::sbp::DetectHoles(bufferSizes,partitionInterval,partitionBufferInterval,
                                                        partitionInterval,boundaryDepth,holeMask,iMask,stencilID);
      if(returnCode != 0){
        std::cout << "Too holy." << std::endl;
        if(boundaryDepth < 3)
          test4Pass[dimId] = false;
      } else {
        if(numDim == 1 || numDim == 2 || numDim == 3){
          std::cout << "Conn after hole detection: " << std::endl;
          for(int iDim = 0;iDim < numDim;iDim++){
            int *expectedConn = &stencilConnectivity[iDim*numPointsBuffer];
            std::cout << "Dimension-" << iDim+1 << ":";
            for(size_t iPoint = 0;iPoint < numPointsBuffer;iPoint++){
              if(!(iPoint%(numPoints*numPoints))) {
                std::cout << std::endl << " --------------- " << std::endl;
              } else if(!(iPoint%numPoints)) std::cout << std::endl;
              std::cout << expectedConn[iPoint] << " ";
            }
            std::cout << std::endl;
          }
        }
      }
      test4Result = test4Result && test4Pass[dimId];
      // ======= middle hole ====      
    }

    delete[] opInterval;
  }
  serialUnitResults.UpdateResult("Operators:DetectHoliness:UnHoly",test1Result);
  serialUnitResults.UpdateResult("Operators:DetectHoliness:HolyBoundaries",test2Result);
  serialUnitResults.UpdateResult("Operators:DetectHoliness:TooHoly",test3Result);
  serialUnitResults.UpdateResult("Operators:DetectHoliness:HolyCenter",test4Result);
}

void TestSBPInitialize(ix::test::results &serialUnitResults)
{
  // Operator
  operators::sbp_operator_t order12Operator;
  operators::sbp::Initialize(order12Operator,2);
  
  bool operatorInitialization12 = true;

  if(order12Operator.numStencils != 4) {
    operatorInitialization12 = false;
  }

  if(order12Operator.numValues != 7) {
    operatorInitialization12 = false;
  }
  
  if(order12Operator.stencilSizes[0] != 2 ||
     order12Operator.stencilSizes[1] != 2 ||
     order12Operator.stencilSizes[2] != 2 ||
     order12Operator.stencilSizes[3] != 1) {
    operatorInitialization12 = false;
  }
  
  int start = 0;
  for(int iStencil = 0;iStencil < order12Operator.numStencils;iStencil++){
    if(order12Operator.stencilStarts[iStencil] != start) operatorInitialization12 = false;
    start += order12Operator.stencilSizes[iStencil];
  }
  
  if(order12Operator.stencilWeights[0] != -0.5 ||
     order12Operator.stencilWeights[1] !=  0.5 ||
     order12Operator.stencilWeights[2] != -1.0 ||
     order12Operator.stencilWeights[3] !=  1.0 ||
     order12Operator.stencilWeights[4] != -1.0 ||
     order12Operator.stencilWeights[5] !=  1.0 ||
     order12Operator.stencilWeights[6] != 0.0) {
    operatorInitialization12 = false;
  }
  
  if(order12Operator.stencilOffsets[0] != -1 ||
     order12Operator.stencilOffsets[1] !=  1 ||
     order12Operator.stencilOffsets[2] !=  0 ||
     order12Operator.stencilOffsets[3] !=  1 ||
     order12Operator.stencilOffsets[4] != -1 ||
     order12Operator.stencilOffsets[5] !=  0 ||
     order12Operator.stencilOffsets[6] !=  0) {
    operatorInitialization12 = false;
  }
  
  serialUnitResults.UpdateResult("Operators:SBP:Initialize12",operatorInitialization12);
  
  operators::sbp_operator_t op24;
  operators::sbp::Initialize(op24,4);
  
  bool op24InitWorks = true;
  if(op24.numStencils != 10)
    op24InitWorks = false;
  if(op24.numValues != 33)
    op24InitWorks = false;
  
  if(op24.stencilSizes[0] != 4 ||
     op24.stencilSizes[1] != 4 ||
     op24.stencilSizes[2] != 2 ||
     op24.stencilSizes[3] != 4 ||
     op24.stencilSizes[4] != 4 ||
     op24.stencilSizes[5] != 4 ||
     op24.stencilSizes[6] != 2 ||
     op24.stencilSizes[7] != 4 ||
     op24.stencilSizes[8] != 4 ||
     op24.stencilSizes[9] != 1)
    op24InitWorks = false;
  
  serialUnitResults.UpdateResult("Operators:SBP:Initialize24",op24InitWorks);

  operators::sbp_operator_t op36;
  operators::sbp::Initialize(op36,6);

  bool op36InitWorks = true;

  op36InitWorks = op36InitWorks && (op36.numStencils == 2 + 2*op36.boundaryDepth);

  for (int i=0; i<op36.boundaryDepth; i++) {
    op36InitWorks = op36InitWorks && (op36.stencilSizes[i+1] == op36.stencilSizes[i+1+op36.boundaryDepth]);
  }

  // tests if the right and left boundary conditions are mirrored
  int offsetSum = 0;
  double weightSum = 0.0;
  for (int i=op36.stencilStarts[1]; i<op36.stencilStarts[op36.numStencils-1]; i++) {
    offsetSum += op36.stencilOffsets[i];
    weightSum += op36.stencilWeights[i];
  }
  op36InitWorks = op36InitWorks && (offsetSum == 0);
  if (!op36InitWorks) std::cout << "sbp36 offsetSum " << offsetSum << std::endl;
  op36InitWorks = op36InitWorks && (std::abs(weightSum) < 1e-14);
  if (!op36InitWorks) std::cout << "sbp36 weightSum " << weightSum << std::endl;

  serialUnitResults.UpdateResult("Operators:SBP:Initialize36",op36InitWorks);
}

void TestOperatorSBP12(ix::test::results &serialUnitResults)
{
  operators::sbp_operator_t sbp12;
  operators::sbp::Initialize(sbp12,2);
  // Forticate the stencil starts
  for(int iStencil = 0;iStencil < sbp12.numStencils;iStencil++)
    sbp12.stencilStarts[iStencil]++;
  
  std::vector<bool> testResults;
  for(int iDim = 0;iDim < 3;iDim++){

    std::cout << "Brute force ApplyOperator SBP12 in " << iDim+1 
              << " dimensions." << std::endl;
    operators::sbp::BruteTest1(sbp12,2,iDim+1,1,4,testResults);
    std::vector<bool>::iterator trIt = testResults.begin();
    bool testResult = true;
    while(trIt != testResults.end()) {
      testResult = testResult&&(*trIt);
      trIt++;
    }
    switch(iDim){
    case 0:
      serialUnitResults.UpdateResult("Operators:SBP12:ApplyOperator:1D",testResult);
      break;
    case 1:
      serialUnitResults.UpdateResult("Operators:SBP12:ApplyOperator:2D",testResult);
      break;
    case 2:
      serialUnitResults.UpdateResult("Operators:SBP12:ApplyOperator:3D",testResult);
      break;
    default: // never get here
      break;
    }
  }

  
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:InteriorXAccurate",testResults[0]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:InteriorXOrder2",testResults[1]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:InteriorYAccurate",testResults[2]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:InteriorYOrder2",testResults[3]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:InteriorZAccurate",testResults[4]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:InteriorZOrder2",testResults[5]);

  serialUnitResults.UpdateResult("Operators:SBP12:Brute:LeftBoundaryXAccurate",testResults[6]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:LeftBoundaryXOrder1",testResults[7]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:LeftBoundaryYAccurate",testResults[8]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:LeftBoundaryYOrder1",testResults[9]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:LeftBoundaryZAccurate",testResults[10]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:LeftBoundaryZOrder1",testResults[11]);

  serialUnitResults.UpdateResult("Operators:SBP12:Brute:RightBoundaryXAccurate",testResults[12]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:RightBoundaryXOrder1",testResults[13]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:RightBoundaryYAccurate",testResults[14]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:RightBoundaryYOrder1",testResults[15]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:RightBoundaryZAccurate",testResults[16]);
  serialUnitResults.UpdateResult("Operators:SBP12:Brute:RightBoundaryZOrder1",testResults[17]);
  
//   std::cout << "New test results:" << std::endl;
//   pcpp::io::DumpContents(std::cout,testResults);
//   std::cout << std::endl;
  
}

void TestOperatorSBP24(ix::test::results &serialUnitResults)
{
  operators::sbp_operator_t sbp24;
  operators::sbp::Initialize(sbp24,4);

  // Forticate the stencil starts
  for(int iStencil = 0;iStencil < sbp24.numStencils;iStencil++)
    sbp24.stencilStarts[iStencil]++;
  
  // std::vector<bool> testResults;
  // operators::sbp::BruteTest1(sbp24,4,3,1,4,testResults);

  std::vector<bool> testResults;
  for(int iDim = 0;iDim < 3;iDim++){
    std::cout << "Brute force ApplyOperator SBP24 in " << iDim+1 
              << " dimensions." << std::endl;
    operators::sbp::BruteTest1(sbp24,4,iDim+1,1,4,testResults);
    std::vector<bool>::iterator trIt = testResults.begin();
    bool testResult = true;
    while(trIt != testResults.end()) {
      testResult = testResult&&(*trIt);
      trIt++;
    }
    switch(iDim){
    case 0:
      serialUnitResults.UpdateResult("Operators:SBP24:ApplyOperator:1D",testResult);
      break;
    case 1:
      serialUnitResults.UpdateResult("Operators:SBP24:ApplyOperator:2D",testResult);
      break;
    case 2:
      serialUnitResults.UpdateResult("Operators:SBP24:ApplyOperator:3D",testResult);
      break;
    default: // never get here
      break;
    }
  }

  serialUnitResults.UpdateResult("Operators:SBP24:Brute:InteriorXAccurate",testResults[0]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:InteriorXOrder4",testResults[1]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:InteriorYAccurate",testResults[2]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:InteriorYOrder4",testResults[3]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:InteriorZAccurate",testResults[4]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:InteriorZOrder4",testResults[5]);

  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBoundaryXAccurate",testResults[6]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBoundaryXOrder2",testResults[7]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBoundaryYAccurate",testResults[8]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBoundaryYOrder2",testResults[9]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBoundaryZAccurate",testResults[10]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBoundaryZOrder2",testResults[11]);

  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias1XAccurate",testResults[12]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias1XOrder2",testResults[13]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias1YAccurate",testResults[14]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias1YOrder2",testResults[15]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias1ZAccurate",testResults[16]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias1ZOrder2",testResults[17]);

  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias2XAccurate",testResults[18]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias2XOrder2",testResults[19]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias2YAccurate",testResults[20]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias2YOrder2",testResults[21]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias2ZAccurate",testResults[22]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias2ZOrder2",testResults[23]);

  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias3XAccurate",testResults[24]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias3XOrder2",testResults[25]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias3YAccurate",testResults[26]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias3YOrder2",testResults[27]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias3ZAccurate",testResults[28]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:LeftBias3ZOrder2",testResults[29]);

  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBoundaryXAccurate",testResults[30]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBoundaryXOrder2",testResults[31]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBoundaryYAccurate",testResults[32]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBoundaryYOrder2",testResults[33]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBoundaryZAccurate",testResults[34]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBoundaryZOrder2",testResults[35]);

  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias1XAccurate",testResults[36]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias1XOrder2",testResults[37]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias1YAccurate",testResults[38]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias1YOrder2",testResults[39]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias1ZAccurate",testResults[40]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias1ZOrder2",testResults[41]);

  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias2XAccurate",testResults[42]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias2XOrder2",testResults[43]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias2YAccurate",testResults[44]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias2YOrder2",testResults[45]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias2ZAccurate",testResults[46]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias2ZOrder2",testResults[47]);

  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias3XAccurate",testResults[48]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias3XOrder2",testResults[49]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias3YAccurate",testResults[50]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias3YOrder2",testResults[51]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias3ZAccurate",testResults[52]);
  serialUnitResults.UpdateResult("Operators:SBP24:Brute:RightBias3ZOrder2",testResults[53]);
  
//   std::cout << "New test results:" << std::endl;
//   pcpp::io::DumpContents(std::cout,testResults);
//   std::cout << std::endl;
}

void TestOperatorSBP36(ix::test::results &serialUnitResults)
{
  operators::sbp_operator_t sbp36;
  operators::sbp::Initialize(sbp36,6);

  // Forticate the stencil starts
  for(int iStencil = 0;iStencil < sbp36.numStencils;iStencil++)
    sbp36.stencilStarts[iStencil]++;

  std::vector<bool> testResults;
  for(int iDim = 0;iDim < 3;iDim++){
    std::cout << "Brute force ApplyOperator SBP36 in " << iDim+1
              << " dimensions." << std::endl;
    operators::sbp::BruteTest1(sbp36,6,iDim+1,1,4,testResults);

    std::vector<bool>::iterator trIt = testResults.begin();
    bool testResult = true;
    while(trIt != testResults.end()) {
      testResult = testResult&&(*trIt);
      trIt++;
    }

    switch(iDim){
    case 0:
      serialUnitResults.UpdateResult("Operators:SBP36:ApplyOperator:1D",testResult);
      break;
    case 1:
      serialUnitResults.UpdateResult("Operators:SBP36:ApplyOperator:2D",testResult);
      break;
    case 2:
      serialUnitResults.UpdateResult("Operators:SBP36:ApplyOperator:3D",testResult);
      break;
    default: // never get here
      break;
    }
  }

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:InteriorXAccurate",testResults[0]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:InteriorXOrder6",testResults[1]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:InteriorYAccurate",testResults[2]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:InteriorYOrder6",testResults[3]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:InteriorZAccurate",testResults[4]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:InteriorZOrder6",testResults[5]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBoundaryXAccurate",testResults[6]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBoundaryXOrder3",testResults[7]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBoundaryYAccurate",testResults[8]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBoundaryYOrder3",testResults[9]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBoundaryZAccurate",testResults[10]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBoundaryZOrder3",testResults[11]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias1XAccurate",testResults[12]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias1XOrder3",testResults[13]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias1YAccurate",testResults[14]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias1YOrder3",testResults[15]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias1ZAccurate",testResults[16]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias1ZOrder3",testResults[17]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias2XAccurate",testResults[18]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias2XOrder3",testResults[19]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias2YAccurate",testResults[20]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias2YOrder3",testResults[21]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias2ZAccurate",testResults[22]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias2ZOrder3",testResults[23]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias3XAccurate",testResults[24]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias3XOrder3",testResults[25]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias3YAccurate",testResults[26]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias3YOrder3",testResults[27]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias3ZAccurate",testResults[28]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias3ZOrder3",testResults[29]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias4XAccurate",testResults[30]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias4XOrder3",testResults[31]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias4YAccurate",testResults[32]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias4YOrder3",testResults[33]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias4ZAccurate",testResults[34]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias4ZOrder3",testResults[35]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias5XAccurate",testResults[36]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias5XOrder3",testResults[37]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias5YAccurate",testResults[38]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias5YOrder3",testResults[39]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias5ZAccurate",testResults[40]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:LeftBias5ZOrder3",testResults[41]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBoundaryXAccurate",testResults[42]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBoundaryXOrder3",testResults[43]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBoundaryYAccurate",testResults[44]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBoundaryYOrder3",testResults[45]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBoundaryZAccurate",testResults[46]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBoundaryZOrder3",testResults[47]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias1XAccurate",testResults[48]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias1XOrder3",testResults[49]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias1YAccurate",testResults[50]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias1YOrder3",testResults[51]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias1ZAccurate",testResults[52]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias1ZOrder3",testResults[53]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias2XAccurate",testResults[54]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias2XOrder3",testResults[55]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias2YAccurate",testResults[56]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias2YOrder3",testResults[57]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias2ZAccurate",testResults[58]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias2ZOrder3",testResults[59]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias3XAccurate",testResults[60]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias3XOrder3",testResults[61]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias3YAccurate",testResults[62]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias3YOrder3",testResults[63]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias3ZAccurate",testResults[64]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias3ZOrder3",testResults[65]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias4XAccurate",testResults[66]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias4XOrder3",testResults[67]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias4YAccurate",testResults[68]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias4YOrder3",testResults[69]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias4ZAccurate",testResults[70]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias4ZOrder3",testResults[71]);

  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias5XAccurate",testResults[72]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias5XOrder3",testResults[73]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias5YAccurate",testResults[74]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias5YOrder3",testResults[75]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias5ZAccurate",testResults[76]);
  serialUnitResults.UpdateResult("Operators:SBP36:Brute:RightBias5ZOrder3",testResults[77]);

}


void TestApplyOperatorBlobs(ix::test::results &serialUnitResults)
{
  
  // Create a test geometry
  int numDim = 3;
  size_t numPoints = 1;
  size_t *numX = new size_t [numDim];
  double *dX   = new double [numDim];
  size_t *opInterval = new size_t [2*numDim];
  int *periodicDirs = new int [numDim];
  for(int iDim = 0;iDim < numDim;iDim++){
    periodicDirs[iDim] = 0;
    //    numX[iDim] = 10 * (std::pow(2,static_cast<double>(numDim-iDim)));
    numX[iDim] = 3;
    dX[iDim] = 1.0/(numX[iDim]-1);
    numPoints *= numX[iDim];
    opInterval[2*iDim]   = 1;
    opInterval[2*iDim+1] = numX[iDim];
  }
  double *U       = new double [numPoints];
  double *dU      = new double [numDim*numPoints];
  double *exactDU = new double [numDim*numPoints];
  int *stencilID  = new int    [numDim*numPoints];

  // Initialize U with some field
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    size_t xPos = iPoint%numX[0];
    U[iPoint] = xPos*dX[0];
    exactDU[iPoint] = dX[0];
    dU[iPoint] = 0.0;
    stencilID[iPoint] = 0;
    for(int iDim = 1;iDim < numDim;iDim++){
      exactDU[iDim*numPoints+iPoint] = 0.0;
      dU[iDim*numPoints+iPoint] = 0.0;
      stencilID[iPoint+iDim*numPoints] = 0;
    }
  }
    
  // Create an SBP operator/stencilset
  operators::sbp_operator_t sbp12;
  operators::sbp::Initialize(sbp12,2);
  
  int numStencils = sbp12.numStencils;
  // Forticate the stencil starts
  for(int iStencil = 0;iStencil < numStencils;iStencil++)
    sbp12.stencilStarts[iStencil]++;
  int boundaryDepth = 1;
  
  // Create stencil connectivity
  operators::sbp::CreateStencilConnectivity(numDim,numX,opInterval,boundaryDepth,
                                                      periodicDirs,stencilID,true);
  bool createStencilConn = true;
  //   for(int iDim = 0;iDim < numDim;iDim++){
  //     std::cout << "Stencil Connectivity for dimension " << iDim+1 << "(" << numX[iDim] << "):" << std::endl;
  //     for(int iPoint = 0;iPoint < numPoints;iPoint++)
  //       std::cout << stencilID[iDim*numPoints+iPoint] << " ";
  //     std::cout << std::endl;
  //   }

  // Invert the connectivity
  //  operators::InvertStencilConnectivity(numDim,opInterval,stencilID,dualStencilConn,true);
  size_t *dualStencilConn = new size_t [numDim*numPoints];    // [DIM1[stencil1Points,stencil2Points....stencil(numStencils)Points]...]
  size_t *numPointsStencil = new size_t [numDim*numStencils]; // [DIM1[numPoints1,numPoints2....numPoints(numStencils)]...]

  operators::sbp::InvertStencilConnectivity(numDim,numX,opInterval,numStencils,stencilID,
                                                      dualStencilConn,numPointsStencil,true);
  size_t stencilPoint = 0;
  bool invertStencilConn = true;
  for(int iDim = 0;iDim < numDim;iDim++){
    std::cout << "Dual stencil conn for dimension " << iDim+1 << "(" << numX[iDim] << ")" << std::endl;
    for(int iStencil = 0;iStencil < numStencils;iStencil++){
      size_t numPointsThisStencil = numPointsStencil[iDim*numStencils+iStencil];
      std::cout << "Stencil[" << iStencil+1 << "] Points(" << numPointsThisStencil << "): ";
      for(size_t iPoint = 0;iPoint < numPointsThisStencil;iPoint++)
        std::cout << dualStencilConn[stencilPoint++] << " ";
    }
    std::cout << std::endl;
  }
  
  int numComponents = 1;
  size_t blobOffset = 0;
  bool applySingleStencilTest = true;
  for(int iDim = 0;iDim < numDim;iDim++){
    int opDir = iDim+1;
    double *dUPtr = &dU[iDim*numPoints];
    for(int iStencil = 0;iStencil < numStencils;iStencil++){
      size_t numPointsApply = numPointsStencil[iDim*numStencils+iStencil];
      int stencilSize = sbp12.stencilSizes[iStencil];
      int stencilStart = sbp12.stencilStarts[iStencil] - 1;
      double *stencilWeights = &sbp12.stencilWeights[stencilStart];
      int *stencilOffsets = &sbp12.stencilOffsets[stencilStart];
      size_t *blobPoints = &dualStencilConn[blobOffset];
      blobOffset += numPointsApply;
      FC_MODULE(operators,applysinglestencil,OPERATORS,APPLYSINGLESTENCIL)(&numDim,numX,&numComponents,&numPoints,&opDir,
                                                                           &numPointsApply,blobPoints,&stencilSize,
                                                                           stencilWeights,stencilOffsets,U,dUPtr);
      
      
    }
    double *exactdU = &exactDU[iDim*numPoints];
    for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
      if(std::abs(dUPtr[iPoint]-exactdU[iPoint]) > 1e-15) {
        applySingleStencilTest = false;
        std::cout << "DU[" << iDim << "][" << iPoint << "] = (" << dUPtr[iPoint] 
                  << "," << exactdU[iPoint] << ")" << std::endl;
      }
    }
  }

  serialUnitResults.UpdateResult("Operators:ApplySingleStencil",applySingleStencilTest);
  
  
  // Initialize U with some field
  for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
    size_t xPos = iPoint%numX[0];
    U[iPoint] = xPos*dX[0];
    exactDU[iPoint] = dX[0];
    dU[iPoint] = 0.0;
    for(int iDim = 1;iDim < numDim;iDim++){
      exactDU[iDim*numPoints+iPoint] = 0.0;
      dU[iDim*numPoints+iPoint] = 0.0;
    }
  }

  int *stencilSizes = sbp12.stencilSizes;
  int *stencilStarts = sbp12.stencilStarts;
  int numStencilValues = sbp12.numValues;
  double *stencilWeights = sbp12.stencilWeights;
  int *stencilOffsets = sbp12.stencilOffsets;

  bool applyOperatorBlobs = true;
  for(int iDim = 0;iDim < numDim;iDim++){
    int opDir = iDim+1;
    size_t *numPointsDimStencil = &numPointsStencil[iDim*numStencils];
    size_t *stencilPoints = &dualStencilConn[iDim*numPoints];
    double *dUPtr = &dU[iDim*numPoints];
    FC_MODULE(operators,applyoperatorblobs,OPERATORS,APPLYOPERATORBLOBS)(&numDim,numX,&numComponents,&numPoints,&opDir,&numStencils,
                                                                         stencilSizes,stencilStarts,&numStencilValues,stencilWeights,
                                                                         stencilOffsets,numPointsDimStencil,&numPoints,
                                                                         stencilPoints,U,dUPtr);
    double *exactdU = &exactDU[iDim*numPoints];
    for(size_t iPoint = 0;iPoint < numPoints;iPoint++){
      if(std::abs(dUPtr[iPoint]-exactdU[iPoint]) > 1e-12) {
        applyOperatorBlobs = false;
        std::cout << "DU[" << iDim << "][" << iPoint << "] = (" << dUPtr[iPoint] 
                  << "," << exactdU[iPoint] << ")" << std::endl;
      }
    }
  }
  serialUnitResults.UpdateResult("Operators:ApplyOperatorBlobs",applyOperatorBlobs);
  serialUnitResults.UpdateResult("Operators:CreateStencilConnectivity",createStencilConn);
  serialUnitResults.UpdateResult("Operators:InvertStencilConnectivity",invertStencilConn);
 
  delete [] dualStencilConn;
  delete [] numPointsStencil;
  delete [] stencilID;
  delete [] numX;
  delete [] dX;
  delete [] U;
  delete [] dU;
  delete [] exactDU;
  delete [] opInterval;
  delete [] periodicDirs;
                           
}
