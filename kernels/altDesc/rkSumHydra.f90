!hash:c68f7fc0b0fdf1b75f769e9812486a3f5c7719d3
SUBROUTINE RK4SUM(numDim,numPoints,bufferSize,bufferInterval,h,K1,K2,K3,K4,stateData)
    
    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: h
    REAL(KIND=8),    INTENT(IN)      :: K1(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: k2(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: K3(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: k4(numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: stateData(numPoints)
    
    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd
    REAL(KIND=8)    :: fac1,fac2

    INTEGER(KIND=4) :: hyb_tid_to_device_f, omp_get_thread_num, mydev
    mydev = hyb_tid_to_device_f(omp_get_thread_num())

    fac1 = h/6.0_8
    fac2 = h/3.0_8
    
    iStart = bufferInterval(1)
    iEnd = bufferInterval(2)
    
    xSize = bufferSize(1)
    
    IF(numDim == 3) THEN
       
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       
       nPlane = xSize*bufferSize(2)

       !$omp target teams distribute parallel do device(mydev) if(mydev>=0) is_device_ptr(K1,K2,K3,K4,stateData,bufferInterval,bufferSize,h,numPoints,numDim) shared(stateData,kStart,kEnd,jStart,jEnd,iStart,iEnd,fac1,fac2,xSize,nPlane) private(yzIndex,zIndex,bufferIndex) default(none) collapse(3)
       DO K = kStart, kEnd
          DO J = jStart, jEnd
             DO I = iStart, iEnd
          zIndex = (K-1)*nPlane
             yzIndex = (J-1)*xSize + zIndex
                bufferIndex = yzIndex + I
                stateData(bufferIndex) = fac1*(K1(bufferIndex)+K4(bufferIndex)) + &
                     fac2*(K2(bufferIndex) + K3(bufferIndex)) + stateData(bufferIndex)
             END DO
          END DO
       END DO
       
    ELSE IF (numDim == 2) THEN
       
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)

       !$omp target teams distribute parallel do device(mydev) if(mydev>=0) is_device_ptr(K1,K2,K3,K4,stateData,bufferInterval,bufferSize,h,numPoints,numDim) shared(jStart,jEnd,iStart,iEnd,fac1,fac2,xSize) private(yIndex,bufferIndex) default(none) collapse(2)
       DO J = jStart, jEnd
          DO I = iStart, iEnd
          yIndex = (J-1)*xSize
             bufferIndex = yIndex + I
             stateData(bufferIndex) = fac1*(K1(bufferIndex)+K4(bufferIndex)) + &
                  fac2*(K2(bufferIndex) + K3(bufferIndex)) + stateData(bufferIndex)
          END DO
       END DO
       
    ELSE IF (numDim == 1) THEN
       
       !$omp target teams distribute parallel do device(mydev) if(mydev>=0) is_device_ptr(K1,K2,K3,K4,stateData,bufferInterval,bufferSize,h,numPoints,numDim) shared(iStart,iEnd,fac1,fac2) default(none)
       DO bufferIndex = iStart, iEnd
          stateData(bufferIndex) = fac1*(K1(bufferIndex)+K4(bufferIndex)) + &
               fac2*(K2(bufferIndex) + K3(bufferIndex)) + stateData(bufferIndex)
       END DO
       
    ENDIF
    
    
  END SUBROUTINE RK4SUM
