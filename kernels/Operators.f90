MODULE OPERATORS

  IMPLICIT NONE

CONTAINS

  !> @file Operators.f90
  !> @brief Implementation of routines that apply differential (or other) operators to all points of the state data in the indicated interval

! 
  !> @brief applyoperator applies an operator specified as a stencil set to the provided state data
  !>
  !> applyoperator is a brute-force method of applying a set of stencils to a given data buffer. The
  !> The stencilset and operand are given as flat, primitive arrays, with supporting data to indicate
  !> the size of the data structures. An additional <em>stencil connectivity</em> is given which 
  !> indicates which stencil (i.e. which stencil from the stencil set) to apply to each point.
  !> This brute-force method simply loops through all points, and applies the stencil indicated by
  !> the resolved stencil connectivity to each point.
  !> @image html images/ApplyOperatorBrute.png "applyoperator brute-force method cartoon" width=256
  !> @param[in] numDim indicates the number of dimensions for the input data (dimSizes,opInterval)
  !> @param[in] dimSizes indicates the number of points in each dimension [Nx Ny Nz]
  !> @param[in] numComponents indicates the number of components in the input/output data (U,dU)
  !> @param[in] numPoints indicates the total number of points [Nx*Ny*Nz] (needed for C/Fortran interface)
  !> @param[in] opDir indicates in which direction the operator will be applied [X=1 | Y=2 | Z=3]
  !> @param[in] opInterval array of size (2*numDim) which indicates the index interval on which to operate [xStart,xEnd,yStart,yEnd...]
  !> @param[in] numStencils indicates the number of stencils in the input stencilset
  !> @param[in] stencilSizes array of size (numStencils) which indicates the number of weight values for each stencil
  !> @param[in] stencilStarts array of size (numStencils) which indicates the starting index into (stencilWeights and stencilOffsets) for each stencil
  !> @param[in] numValues total number of stencil weight values (numStencils[1]*numStencils[2]*...numStecils[numStencils]) (req'd for C/Fort interface)
  !> @param[in] stencilWeights array of size (numValues) which contains the weights for all the stencils in the stencilset
  !> @param[in] stencilOffsets array of size (numValues) which indicates the offsets from the @e considered point to the point where each weight is applied
  !> @param[in] stencilID array of size (numPoints) which indicates the stencil ID for each point 
  !> @param[in] U the data on which to operate
  !> @param[out] dU where to stuff the result
  SUBROUTINE applyoperator(numDim,dimSizes,numComponents,numPoints,opDir,opInterval,numStencils, &
       stencilSizes,stencilStarts,numValues,stencilWeights,stencilOffsets,stencilID,U,dU)

    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)  :: numDim, opDir,numStencils, numValues, numComponents
    INTEGER(KIND=8), INTENT(IN)  :: dimSizes(numDim),  numPoints
    INTEGER(KIND=8), INTENT(IN)  :: opInterval(2*numDim) 
    INTEGER(KIND=4), INTENT(IN)  :: stencilSizes(numStencils),stencilStarts(numStencils)
    INTEGER(KIND=4), INTENT(IN)  :: stencilOffsets(numValues)
    REAL(KIND=8),    INTENT(IN)  :: stencilWeights(numValues)
    INTEGER(KIND=4), INTENT(IN)  :: stencilID(numPoints)

    REAL(KIND=8),    INTENT(IN),  TARGET :: U(numPoints*numComponents)
    REAL(KIND=8),    INTENT(OUT), TARGET :: dU(numPoints*numComponents)

    REAL(KIND=8)    :: fac
    INTEGER(KIND=4) :: iStencil, iWeight, iComp
    INTEGER(KIND=8) :: plane, pointOffset, compOffset
    INTEGER(KIND=8) :: I, J, K, jIndex, jkIndex, kIndex, iPoint

    REAL(KIND=8),    DIMENSION(:), POINTER     :: uCompPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: duCompPtr

    !    DO iComp = 1, numComponents

    !       compOffset = (iComp-1)*numPoints
    !       uCompPtr   => U(compOffset+1:compOffset+numPoints)
    !       duCompPtr  => dU(compOffset+1:compOffset+numPoints)

    IF(numDim == 1) THEN
       DO I = opInterval(1), opInterval(2)
          iStencil = stencilID(I)
          dU(I) = 0.0_8
          DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
             dU(I) = dU(I) + stencilWeights(iWeight)*U(I+stencilOffsets(iWeight))
          END DO
       END DO
    ELSE IF(numDim == 2) THEN
       IF(opDir == 1) THEN
          DO J = opInterval(3), opInterval(4)
             jIndex = (J-1)*dimSizes(1)
             DO I = opInterval(1), opInterval(2)
                iPoint = jIndex + I
                iStencil = stencilID(iPoint)
                dU(iPoint) = 0.0_8
                DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
                   dU(iPoint) = dU(iPoint) + &
                        stencilWeights(iWeight)*U(iPoint+stencilOffsets(iWeight))
                END DO
             END DO
          END DO
       ELSE IF(opDir == 2) THEN
          plane = dimSizes(1)
          DO J = opInterval(3), opInterval(4)
             jIndex = (J-1)*plane
             DO I = opInterval(1), opInterval(2)
                iPoint = jIndex + I
                iStencil = stencilID(iPoint)
                dU(iPoint) = 0.0_8
                DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
                   dU(iPoint) = dU(iPoint) + &
                        stencilWeights(iWeight)*U(iPoint+stencilOffsets(iWeight)*plane)
                END DO
             END DO
          END DO
       ENDIF
    ELSE IF(numDim == 3) THEN
       plane = dimSizes(1) * dimSizes(2)
       IF(opDir == 1) THEN
          DO K = opInterval(5), opInterval(6)
             kIndex = (K-1)*plane
             DO J = opInterval(3), opInterval(4)
                jkIndex = kIndex + (J-1)*dimSizes(1)
                DO I = opInterval(1), opInterval(2)
                   iPoint   = jkIndex + I
                   iStencil = stencilID(iPoint) 
                   dU(iPoint) = 0.0_8
                   DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
                      dU(iPoint) = dU(iPoint) + &
                           stencilWeights(iWeight)*U(iPoint+stencilOffsets(iWeight))
                   END DO
                END DO
             END DO
          END DO
       ELSE 
          pointOffset = 1
          IF(opDir == 2) THEN
             pointOffset = dimSizes(1)
          ELSE
             pointOffset = plane
          ENDIF
          DO K = opInterval(5), opInterval(6)
             kIndex = (K-1)*plane
             DO J = opInterval(3), opInterval(4)
                jkIndex = kIndex + (J-1)*dimSizes(1)
                DO I = opInterval(1), opInterval(2)
                   iPoint   = jkIndex + I
                   iStencil = stencilID(iPoint)
                   dU(iPoint) = 0.0_8
                   DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
                      dU(iPoint) = dU(iPoint) + &
                           stencilWeights(iWeight)*U(iPoint+stencilOffsets(iWeight)*pointOffset)
                   END DO ! weight
                END DO ! I
             END DO ! J 
          END DO ! K
       ENDIF ! (opDir)
    ENDIF ! (numDim)
  END SUBROUTINE applyoperator
  
  !> @file ApplyOperator.f90
  !> @brief Implementation of routines that apply differential (or other) operators to all points of the state data in the indicated interval

  !> @brief applyoperator applies an operator specified as a stencil set to the provided state data
  !>
  !> applyoperator is a brute-force method of applying a set of stencils to a given data buffer. The
  !> The stencilset and operand are given as flat, primitive arrays, with supporting data to indicate
  !> the size of the data structures. An additional <em>stencil connectivity</em> is given which 
  !> indicates which stencil (i.e. which stencil from the stencil set) to apply to each point.
  !> This brute-force method simply loops through all points, and applies the stencil indicated by
  !> the resolved stencil connectivity to each point.
  !> @image html images/ApplyOperatorBrute.png "applyoperator brute-force method cartoon" width=256
  !> @param[in] numDim indicates the number of dimensions for the input data (dimSizes,opInterval)
  !> @param[in] dimSizes indicates the number of points in each dimension [Nx Ny Nz]
  !> @param[in] numComponents indicates the number of components in the input/output data (U,dU)
  !> @param[in] numPoints indicates the total number of points [Nx*Ny*Nz] (needed for C/Fortran interface)
  !> @param[in] opDir indicates in which direction the operator will be applied [X=1 | Y=2 | Z=3]
  !> @param[in] opInterval array of size (2*numDim) which indicates the index interval on which to operate [xStart,xEnd,yStart,yEnd...]
  !> @param[in] numStencils indicates the number of stencils in the input stencilset
  !> @param[in] stencilSizes array of size (numStencils) which indicates the number of weight values for each stencil
  !> @param[in] stencilStarts array of size (numStencils) which indicates the starting index into (stencilWeights and stencilOffsets) for each stencil
  !> @param[in] numValues total number of stencil weight values (numStencils[1]*numStencils[2]*...numStecils[numStencils]) (req'd for C/Fort interface)
  !> @param[in] stencilWeights array of size (numValues) which contains the weights for all the stencils in the stencilset
  !> @param[in] stencilOffsets array of size (numValues) which indicates the offsets from the @e considered point to the point where each weight is applied
  !> @param[in] stencilID array of size (numPoints) which indicates the stencil ID for each point 
  !> @param[in] U the data on which to operate
  !> @param[out] dU where to stuff the result
  SUBROUTINE applyoperatorv(numDim,dimSizes,numComponents,numPoints,opDir,opInterval,numStencils, &
       stencilSizes,stencilStarts,numValues,stencilWeights,stencilOffsets,stencilID,U,dU)

    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)  :: numDim, opDir,numStencils, numValues, numComponents
    INTEGER(KIND=8), INTENT(IN)  :: dimSizes(numDim),  numPoints
    INTEGER(KIND=8), INTENT(IN)  :: opInterval(2*numDim) 
    INTEGER(KIND=4), INTENT(IN)  :: stencilSizes(numStencils),stencilStarts(numStencils)
    INTEGER(KIND=4), INTENT(IN)  :: stencilOffsets(numValues)
    REAL(KIND=8),    INTENT(IN)  :: stencilWeights(numValues)
    INTEGER(KIND=4), INTENT(IN)  :: stencilID(numPoints)

    REAL(KIND=8),    INTENT(IN),  TARGET :: U(numPoints*numComponents)
    REAL(KIND=8),    INTENT(OUT), TARGET :: dU(numPoints*numComponents)

    REAL(KIND=8)    :: fac
    INTEGER(KIND=4) :: iStencil, iWeight, iComp
    INTEGER(KIND=8) :: plane, pointOffset, compOffset
    INTEGER(KIND=8) :: I, J, K, jIndex, jkIndex, kIndex, iPoint

    REAL(KIND=8),    DIMENSION(:), POINTER     :: uCompPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: duCompPtr

    DO iComp = 1, numComponents

       compOffset = (iComp-1)*numPoints
       uCompPtr   => U(compOffset+1:compOffset+numPoints)
       duCompPtr  => dU(compOffset+1:compOffset+numPoints)

       IF(numDim == 1) THEN
          DO I = opInterval(1), opInterval(2)
             iStencil = stencilID(I)
             duCompPtr(I) = 0.0_8
             DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
                duCompPtr(I) = duCompPtr(I) + stencilWeights(iWeight)*uCompPtr(I+stencilOffsets(iWeight))
             END DO
          END DO
       ELSE IF(numDim == 2) THEN
          IF(opDir == 1) THEN
             DO J = opInterval(3), opInterval(4)
                jIndex = (J-1)*dimSizes(1)
                DO I = opInterval(1), opInterval(2)
                   iPoint = jIndex + I
                   iStencil = stencilID(iPoint)
                   duCompPtr(iPoint) = 0.0_8
                   DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
                      duCompPtr(iPoint) = duCompPtr(iPoint) + &
                           stencilWeights(iWeight)*uCompPtr(iPoint+stencilOffsets(iWeight))
                   END DO
                END DO
             END DO
          ELSE IF(opDir == 2) THEN
             plane = dimSizes(1)
             DO J = opInterval(3), opInterval(4)
                jIndex = (J-1)*plane
                DO I = opInterval(1), opInterval(2)
                   iPoint = jIndex + I
                   iStencil = stencilID(iPoint)
                   duCompPtr(iPoint) = 0.0_8
                   DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
                      duCompPtr(iPoint) = duCompPtr(iPoint) + &
                           stencilWeights(iWeight)*uCompPtr(iPoint+stencilOffsets(iWeight)*plane)
                   END DO
                END DO
             END DO
          ENDIF
       ELSE IF(numDim == 3) THEN
          plane = dimSizes(1) * dimSizes(2)
          IF(opDir == 1) THEN
             DO K = opInterval(5), opInterval(6)
                kIndex = (K-1)*plane
                DO J = opInterval(3), opInterval(4)
                   jkIndex = kIndex + (J-1)*dimSizes(1)
                   DO I = opInterval(1), opInterval(2)
                      iPoint   = jkIndex + I
                      iStencil = stencilID(iPoint) 
                      duCompPtr(iPoint) = 0.0_8
                      DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
                         duCompPtr(iPoint) = duCompPtr(iPoint) + &
                              stencilWeights(iWeight)*uCompPtr(iPoint+stencilOffsets(iWeight))
                      END DO
                   END DO
                END DO
             END DO
          ELSE 
             pointOffset = 1
             IF(opDir == 2) THEN
                pointOffset = dimSizes(1)
             ELSE
                pointOffset = plane
             ENDIF
             DO K = opInterval(5), opInterval(6)
                kIndex = (K-1)*plane
                DO J = opInterval(3), opInterval(4)
                   jkIndex = kIndex + (J-1)*dimSizes(1)
                   DO I = opInterval(1), opInterval(2)
                      iPoint   = jkIndex + I
                      iStencil = stencilID(iPoint)
                      duCompPtr(iPoint) = 0.0_8
                      DO  iWeight = stencilStarts(iStencil),stencilStarts(iStencil)+stencilSizes(iStencil) - 1
                         duCompPtr(iPoint) = duCompPtr(iPoint) + &
                              stencilWeights(iWeight)*uCompPtr(iPoint+stencilOffsets(iWeight)*pointOffset)
                      END DO ! weight
                   END DO ! I
                END DO ! J 
             END DO ! K
          ENDIF ! (opDir)
       ENDIF ! (numDim)
    END DO ! (loop over components)
  END SUBROUTINE applyoperatorv

  !> @brief applyoperatorblobs applies an operator by applying each stencil in turn to all the points to which it applies
  !>
  !> applyoperatorblobs is a @e blobbed method of applying a set of stencils to a given data buffer. 
  !> The stencilset and operand are given as flat, primitive arrays, with supporting data to indicate
  !> the size of the data structures. An additional <em>dual stencil connectivity</em> is given which 
  !> indicates which points to apply to a given stencil.
  !> This blobbed method loops through all the stencils, and applies each stencil to the set of pionts indicated by
  !> the resolved dual stencil connectivity for each stencil.
  !> @image html images/ApplyOperatorBlobs.png "applyoperatorBLOBS - blobbed method cartoon" width=256
  !> @param[in] numDim indicates the number of dimensions for the input data (dimSizes,opInterval)
  !> @param[in] dimSizes indicates the number of points in each dimension [Nx Ny Nz]
  !> @param[in] numComponents indicates the number of components in the input/output data (U,dU)
  !> @param[in] numPointsBuffer indicates the total number of points [Nx*Ny*Nz] (needed for C/Fortran interface)
  !> @param[in] opDir indicates in which direction the operator will be applied [X=1 | Y=2 | Z=3]
  !> @param[in] numStencils indicates the number of stencils in the input stencilset
  !> @param[in] stencilSizes array of size (numStencils) which indicates the number of weight values for each stencil
  !> @param[in] stencilStarts array of size (numStencils) which indicates the starting index into (stencilWeights and stencilOffsets) for each stencil
  !> @param[in] numStencilValues total number of stencil weight values (numStencils[1]*numStencils[2]*...numStecils[numStencils]) (req'd for C/Fort interface)
  !> @param[in] stencilWeights array of size (numValues) which contains the weights for all the stencils in the stencilset
  !> @param[in] stencilOffsets array of size (numValues) which indicates the offsets from the @e considered point to the point where each weight is applied
  !> @param[in] numPointsStencil array of size (numStencils) which indicates how many points to apply each stencil
  !> @param[in] numPointsApply total number of points in the stencilPoints array (needed for C/Fortran interface)
  !> @param[in] stencilPoints array of size (numPointsStencil(1)*numPointsStencil(2)*...numPointsStencil(numStencils)) indicating the points to which each stencil applies
  !> @param[in] U the data on which to operate
  !> @param[out] dU where to stuff the result
  SUBROUTINE applyoperatorblobs(numDim,dimSizes,numComponents,numPointsBuffer,opDir,numStencils, &
       stencilSizes,stencilStarts,numStencilValues,stencilWeights,stencilOffsets,numPointsStencil,&
       numPointsApply,stencilPoints,U,dU)

    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)  :: numDim, opDir,numStencils, numStencilValues, numComponents
    INTEGER(KIND=8), INTENT(IN)  :: numPointsApply, numPointsBuffer
    INTEGER(KIND=8), INTENT(IN)  :: dimSizes(numDim),numPointsStencil(numStencils)
    INTEGER(KIND=8), INTENT(IN)  :: stencilPoints(numPointsApply)
    INTEGER(KIND=4), INTENT(IN)  :: stencilSizes(numStencils),stencilStarts(numStencils)
    INTEGER(KIND=4), INTENT(IN)  :: stencilOffsets(numStencilValues)
    REAL(KIND=8),    INTENT(IN)  :: stencilWeights(numStencilValues)
    REAL(KIND=8),    INTENT(IN)  :: U(numPointsBuffer)
    REAL(KIND=8),    INTENT(OUT) :: dU(numPointsBuffer)

    INTEGER(KIND=8) :: stencilPointsOffset, iPoint, pointsStart,pointsEnd,numPointsThisStencil
    INTEGER(KIND=4) :: iStencil, iWeight, stencilSize, stencilStart, stencilEnd

    !  dU(:) = 0.0_8
    stencilPointsOffset = 1
    !  WRITE(*,*) 'NumDim: ',numDim 
    !  WRITE(*,*) 'dimSizes:',dimSizes
    !  WRITE(*,*) 'numComponents:',numComponents
    !  WRITE(*,*) 'numPointsBuffer:',numPointsBuffer
    !  WRITE(*,*) 'opDir:',opDir
    !  WRITE(*,*) 'numStencils:',numStencils
    !  WRITE(*,*) 'stencilSizes:',stencilSizes 
    !  WRITE(*,*) 'stencilStarts:',stencilStarts
    !  WRITE(*,*) 'numStencilValues:',numStencilValues
    !  WRITE(*,*) 'stencilWeights:',stencilWeights
    !  WRITE(*,*) 'stencilOffsets:',stencilOffsets
    !  WRITE(*,*) 'numPointsStencil:',numPointsStencil
    !    WRITE(*,*) 'numPointsApply:',numPointsApply
    DO iStencil = 1, numStencils
       numPointsThisStencil = numPointsStencil(iStencil)
       stencilSize          = stencilSizes(iStencil)
       stencilStart         = stencilStarts(iStencil)
       stencilEnd           = stencilStart+stencilSize-1
       pointsStart          = stencilPointsOffset
       pointsEnd            = pointsStart + numPointsThisStencil-1
       CALL applysinglestencil(numDim,dimSizes,numComponents,numPointsBuffer,opDir,numPointsThisStencil, &
            stencilPoints(pointsStart:pointsEnd),stencilSize,stencilWeights(stencilStart:stencilEnd),&
            stencilOffsets(stencilStart:stencilEnd),U,dU)
       stencilPointsOffset = stencilPointsOffset+numPointsThisStencil
    END DO

  END SUBROUTINE applyoperatorblobs



  !> @brief applysinglestencil applies an operator by applying a given stencil to the specified points
  !>
  !> applysinglestencil is a single-stencil method which operates on the given points
  !> The stencil and operand are given as flat, primitive arrays, with supporting data to indicate
  !> the size of the data structures. An additional array of points is given which 
  !> indicates the points on which to operate. This single-stencil method loops through all 
  !> the specified points and applies the stencil to each.
  !> @image html images/ApplyOperatorBlobs.png "applyoperatorBLOBS - blobbed method cartoon" width=256
  !> @param[in] numDim indicates the number of dimensions for the input data (dimSizes,opInterval)
  !> @param[in] dimSizes indicates the number of points in each dimension [Nx Ny Nz]
  !> @param[in] numComponents indicates the number of components in the input/output data (U,dU)
  !> @param[in] numPointsBuffer indicates the total number of points [Nx*Ny*Nz] (needed for C/Fortran interface)
  !> @param[in] opDir indicates in which direction the operator will be applied [X=1 | Y=2 | Z=3]
  !> @param[in] numPointsApply total number of points in the stencilPoints array (needed for C/Fortran interface)
  !> @param[in] applyPoints array of size (numPointsApply) indicating the points to which to apply the stencil
  !> @param[in] stencilSize number of stencil weights
  !> @param[in] stencilWeights array of size (stencilSize) which contains the weights for all the stencils in the stencilset
  !> @param[in] stencilOffsets array of size (stencilSize) which indicates the offsets from the @e considered point to the point where each weight is applied
  !> @param[in] U the data on which to operate
  !> @param[out] dU where to stuff the result
  SUBROUTINE applysinglestencil(numDim,dimSizes,numComponents,numPointsBuffer,opDir,numPointsApply, &
       applyPoints,stencilSize,stencilWeights,stencilOffsets,U,dU)

    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)  :: numDim, opDir, stencilSize, numComponents
    INTEGER(KIND=8), INTENT(IN)  :: dimSizes(numDim),  numPointsBuffer, numPointsApply
    INTEGER(KIND=8), INTENT(IN)  :: applyPoints(numPointsApply) 
    INTEGER(KIND=4), INTENT(IN)  :: stencilOffsets(stencilSize)
    REAL(KIND=8),    INTENT(IN)  :: stencilWeights(stencilSize)
    REAL(KIND=8),    INTENT(IN)  :: U(numPointsBuffer)
    REAL(KIND=8),    INTENT(OUT) :: dU(numPointsBuffer)

    INTEGER(KIND=8) :: I, iPoint, plane
    INTEGER(KIND=4) :: iStencil, iWeight

    IF(opDir == 1) THEN
       DO I = 1,numPointsApply
          iPoint = applyPoints(I)
          dU(iPoint) = 0.0_8
          DO  iWeight = 1, stencilSize
             dU(iPoint) = dU(iPoint) + stencilWeights(iWeight)*U(iPoint+stencilOffsets(iWeight))
          END DO
       END DO
    ELSE
       IF(opDir == 2) THEN
          plane = dimSizes(1)
       ELSE IF(opDir == 3) THEN
          plane = dimSizes(1)*dimSizes(2)
       ENDIF
       DO I = 1, numPointsApply
          iPoint = applyPoints(I) 
          dU(iPoint) = 0.0_8
          DO  iWeight = 1,stencilSize
             dU(iPoint) = dU(iPoint) + stencilWeights(iWeight)*U(iPoint+stencilOffsets(iWeight)*plane)
          END DO
       END DO
    END IF

  END SUBROUTINE applysinglestencil


  !> @brief YAXPY point-wise operator performing Y = aX + Y (scalar a)
  SUBROUTINE YAXPY(numDim,numPoints,bufferSize,bufferInterval,a,X,Y)
    
    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a,X(numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: Y(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = Y(I) + a*X(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Y(bufferIndex) = Y(bufferIndex) + a*X(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Y(bufferIndex) = Y(bufferIndex) + a*X(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE YAXPY

  !> @brief ZAXPY point-wise operator performing Z = aX + Y (scalar a)
  SUBROUTINE ZAXPY(numDim,numPoints,bufferSize,bufferInterval,a,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a,X(numPoints),Y(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Z(I) = Y(I) + a*X(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Z(bufferIndex) = Y(bufferIndex) + a*X(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Z(bufferIndex) = Y(bufferIndex) + a*X(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE ZAXPY
  
  !> @brief ZAXPBY point-wise operator performing Z = aX + bY (scalar a,b)
  SUBROUTINE ZAXPBY(numDim,numPoints,bufferSize,bufferInterval,a,b,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a,b,X(numPoints),Y(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Z(I) = b*Y(I) + a*X(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Z(bufferIndex) = b*Y(bufferIndex) + a*X(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Z(bufferIndex) = b*Y(bufferIndex) + a*X(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE ZAXPBY

  !> @brief XAXM1 point-wise operator performing X = a/X (scalar a)
  SUBROUTINE XAXM1(numDim,numPoints,bufferSize,bufferInterval,a,X)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(INOUT)   :: X(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          X(I) = a/X(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             X(bufferIndex) = a/X(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                X(bufferIndex) = a/X(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE XAXM1
  

  !> @brief XAM1X point-wise operator performing X = X/a (scalar a)
  SUBROUTINE XAM1X(numDim,numPoints,bufferSize,bufferInterval,a,X)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(INOUT)   :: X(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          X(I) = X(I)/a
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             X(bufferIndex) = X(bufferIndex)/a
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                X(bufferIndex) = X(bufferIndex)/a
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE XAM1X
  


  !> @brief YAXPBY point-wise operator performing Y = aX + bY (scalar a,b)
  SUBROUTINE YAXPBY(numDim,numPoints,bufferSize,bufferInterval,a,b,X,Y)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a,b,X(numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: Y(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = b*Y(I) + a*X(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Y(bufferIndex) = b*Y(bufferIndex) + a*X(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Y(bufferIndex) = b*Y(bufferIndex) + a*X(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE YAXPBY

  !> @brief YAX point-wise operator performing Y = aX  (scalar a)
  SUBROUTINE YAX(numDim,numPoints,bufferSize,bufferInterval,a,X,Y)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a,X(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Y(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = a*X(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Y(bufferIndex) = a*X(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Y(bufferIndex) = a*X(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE YAX

  

  !> @brief ZXYM1 point-wise operator performing Z = X/Y (all vectors)
  SUBROUTINE ZXYM1(numDim,numPoints,bufferSize,bufferInterval,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints),Y(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Z(I) = X(I)/Y(I)
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Z(bufferIndex) = X(bufferIndex)/Y(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Z(bufferIndex) = X(bufferIndex)/Y(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE ZXYM1

  !> @brief ZXY point-wise operator performing Z = XY (all vectors)
  !>
  !> ZXY performs Z = XY where @em X, @em Y are each contiguous arrays.
  !> Operand arrays are contiguous of size @em numPoints.
  !> The kernel operates on the rectangular interval specified by
  !> @em bufferInterval.
  !> The shape of the input and output arrays are
  !> specified by @em bufferSize, which is an @em numDim - dimensional
  !> array that specifies the size in each of @em numDim dimesions.
  !>
  !> @param numDim - const integer input specifies the number of dimensions of the input and output arrays
  !> @param numPoints - const 64-bit integer input specifies the total size of the input and output arrays
  !> @param bufferSize - const 64-bit integer array of size @em numDim specifies the size of the input and output arrays in each of @em numDim dimensions 
  !> @param bufferInterval - const 64-bit integer array of size 2 x @em numDim indicating the rectangular interval in which the kernel should operate; e.g. [ @em iStart @em iEnd @em jStart @em jEnd ]
  !> @param X - const double precision input array
  !> @param Y - const double precision input array
  !> @param Z - double precision output array
  SUBROUTINE ZXY(numDim,numPoints,bufferSize,bufferInterval,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: Y(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Z(I) = X(I)*Y(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Z(bufferIndex) = X(bufferIndex)*Y(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Z(bufferIndex) = X(bufferIndex)*Y(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE ZXY

  !> @brief ZAXY point-wise operator performing Z = aXY (scalar a, vectors X,Y)
  SUBROUTINE ZAXY(numDim,numPoints,bufferSize,bufferInterval,a,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: Y(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Z(I) = a*X(I)*Y(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Z(bufferIndex) = a*X(bufferIndex)*Y(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Z(bufferIndex) = a*X(bufferIndex)*Y(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE ZAXY

  !> @brief YXY point-wise operator performing Y = XY (all vectors)
  SUBROUTINE YXY(numDim,numPoints,bufferSize,bufferInterval,X,Y)

    INTEGER(KIND=4), INTENT(IN)    :: numDim
    INTEGER(KIND=8), INTENT(IN)    :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)    :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)    :: X(numPoints)
    REAL(KIND=8),    INTENT(INOUT) :: Y(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = X(I)*Y(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Y(bufferIndex) = X(bufferIndex)*Y(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Y(bufferIndex) = X(bufferIndex)*Y(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE YXY

  !> @brief ZWXPY point-wise operator performing Z = WX + Y, where all are vectors
  SUBROUTINE ZWXPY(numDim,numPoints,bufferSize,bufferInterval,W,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: W(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: Y(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Z(I) = X(I)*W(I) + Y(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Z(bufferIndex) = X(bufferIndex)*W(bufferIndex) + Y(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Z(bufferIndex) = X(bufferIndex)*W(bufferIndex) + Y(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE ZWXPY

  !> @brief YWXPY point-wise operator performing Y = WX + Y, where all are vectors
  SUBROUTINE YWXPY(numDim,numPoints,bufferSize,bufferInterval,W,X,Y)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: W(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: Y(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = X(I)*W(I) + Y(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Y(bufferIndex) = X(bufferIndex)*W(bufferIndex) + Y(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Y(bufferIndex) = X(bufferIndex)*W(bufferIndex) + Y(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE YWXPY

  !> @brief ZAWPXY point-wise operator performing Z = aW + XY 
  !>
  !> ZAWPXY kernel performs Z = aW + XY with scalar @em a, contiguous input vectors
  !> @em W, @em X, and @em Y, and output vector @em Z.   
  !> Operand arrays are contiguous of size @em numPoints.
  !> The kernel operates on the rectangular interval specified by
  !> @em bufferInterval.
  !> The shape of the input and output arrays are
  !> specified by @em bufferSize, which is an @em numDim - dimensional
  !> array that specifies the size in each of @em numDim dimesions.
  !>
  !> @param numDim - const integer input specifies the number of dimensions of the input and output arrays
  !> @param numPoints - const 64-bit integer input specifies the total size of the input and output arrays
  !> @param bufferSize - const 64-bit integer array of size @em numDim specifies the size of the input and output arrays in each of @em numDim dimensions 
  !> @param bufferInterval - const 64-bit integer array of size 2 x @em numDim indicating the rectangular interval in which the kernel should operate; e.g. [ @em iStart @em iEnd @em jStart @em jEnd ]
  !> @param a - const input double precision scalar
  !> @param W - const double precision contiguous array of size @em numPoints
  !> @param X - const double precision contiguous array of size @em numPoints
  !> @param Y - const double precision contiguous array of size @em numPoints
  !> @param Z - output double precision contiguous array of size @em numpoints
  SUBROUTINE ZAWPXY(numDim,numPoints,bufferSize,bufferInterval,a,W,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(IN)      :: W(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: Y(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Z(I) = a*W(I) + X(I)*Y(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Z(bufferIndex) = a*W(bufferIndex) + X(bufferIndex)*Y(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Z(bufferIndex) = a*W(bufferIndex) + X(bufferIndex)*Y(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE ZAWPXY

  !> @brief ZVWPXY point-wise operator performing Z = VW + XY 
  !>
  !> ZVWPXY kernel performs Z = VW + XY with contiguous double precision
  !> input arrays  @em V, @em W, @em X, and @em Y, and output array @em Z.
  !> Operand arrays are contiguous of size @em numPoints.
  !> The kernel operates on the rectangular interval specified by
  !> @em bufferInterval.
  !> The shape of the input and output arrays are
  !> specified by @em bufferSize, which is an @em numDim - dimensional
  !> array that specifies the size in each of @em numDim dimesions.
  !>
  !> @param numDim - const integer input specifies the number of dimensions of the input and output arrays
  !> @param numPoints - const 64-bit integer input specifies the total size of the input and output arrays
  !> @param bufferSize - const 64-bit integer array of size @em numDim specifies the size of the input and output arrays in each of @em numDim dimensions 
  !> @param bufferInterval - const 64-bit integer array of size 2 x @em numDim indicating the rectangular interval in which the kernel should operate; e.g. [ @em iStart @em iEnd @em jStart @em jEnd ]
  !> @param V - input const double precision contiguous array of size @em numPoints
  !> @param W - input const double precision contiguous array of size @em numPoints
  !> @param X - input const double precision contiguous array of size @em numPoints
  !> @param Y - input const double precision contiguous array of size @em numPoints
  !> @param Z - output double precision contiguous array of size @em numpoints
  SUBROUTINE ZVWPXY(numDim,numPoints,bufferSize,bufferInterval,V,W,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: V(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: W(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: Y(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Z(I) = V(I)*W(I) + X(I)*Y(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Z(bufferIndex) = V(bufferIndex)*W(bufferIndex) + &
                  X(bufferIndex)*Y(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Z(bufferIndex) = V(bufferIndex)*W(bufferIndex) + &
                     X(bufferIndex)*Y(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE ZVWPXY
  
  !> @brief ZWMXPY point-wise operator performing Z = W(X+Y) where all are vectors
  !
  !> ZVWPXY kernel performs Z = W(X+Y) with contiguous double precision
  !> input arrays  @em W, @em X, and @em Y, and output array @em Z.
  !> Operand arrays are contiguous of size @em numPoints.
  !> The kernel operates on the rectangular interval specified by
  !> @em bufferInterval.
  !> The shape of the input and output arrays are
  !> specified by @em bufferSize, which is an @em numDim - dimensional
  !> array that specifies the size in each of @em numDim dimesions.
  !>
  !> @param numDim - const integer input specifies the number of dimensions of the input and output arrays
  !> @param numPoints - const 64-bit integer input specifies the total size of the input and output arrays
  !> @param bufferSize - const 64-bit integer array of size @em numDim specifies the size of the input and output arrays in each of @em numDim dimensions 
  !> @param bufferInterval - const 64-bit integer array of size 2 x @em numDim indicating the rectangular interval in which the kernel should operate; e.g. [ @em iStart @em iEnd @em jStart @em jEnd ]
  !> @param W - input const double precision contiguous array of size @em numPoints
  !> @param X - input const double precision contiguous array of size @em numPoints
  !> @param Y - input const double precision contiguous array of size @em numPoints
  !> @param Z - output double precision contiguous array of size @em numpoints
  SUBROUTINE ZWMXPY(numDim,numPoints,bufferSize,bufferInterval,W,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: W(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: Y(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Z(I) = W(I)*(X(I)+Y(I))
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Z(bufferIndex) = W(bufferIndex)*(X(bufferIndex)+Y(bufferIndex))
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Z(bufferIndex) = W(bufferIndex)*(X(bufferIndex) + Y(bufferIndex))
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE ZWMXPY

  !> @brief ZXDOTY numComponents-vector inner product Z =  X * Y 
  SUBROUTINE ZXDOTY(numDim,numPoints,bufferSize,bufferInterval,numComponents,X,Y,Z)

    INTEGER(KIND=4), INTENT(IN)      :: numDim, numComponents
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints*numComponents)
    REAL(KIND=8),    INTENT(IN)      :: Y(numPoints*numComponents)
    REAL(KIND=8),    INTENT(OUT)     :: Z(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=4) :: iComp
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd, compOffset

    REAL(KIND=8)    :: zero = 0.0_8

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)
    
    CALL ASSIGNMENTXA(numDim,numPoints,bufferSize,bufferInterval, &
         zero,Z)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          DO iComp = 0, numComponents-1
             compOffset = iComp*numPoints
             Z(I) = Z(I) + X(I+compOffset)*Y(I+compOffset)
          END DO
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             DO iComp = 0,numComponents-1
                compOffset = iComp*numPoints
                
                bufferIndex = yIndex + I
                Z(bufferIndex) = Z(bufferIndex) + &
                     X(bufferIndex+compOffset)*Y(bufferIndex+compOffset)
             END DO
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                DO iComp = 0,numComponents-1
                   compOffset = iComp*numPoints
                   
                   bufferIndex = yzIndex + I
                   Z(bufferIndex) = Z(bufferIndex) + &
                        X(bufferIndex+compOffset)*Y(bufferIndex+compOffset)
                END DO
             END DO
          END DO
       END DO
    ENDIF
  END SUBROUTINE ZXDOTY

  !> @brief XAX point-wise operator performing X = aX (scalar a)
  SUBROUTINE XAX(numDim,numPoints,bufferSize,bufferInterval,a,X)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(INOUT)   :: X(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          X(I) = a*X(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             X(bufferIndex) = a*X(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                X(bufferIndex) = a*X(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE XAX

  !> @brief ASSIGNMENTYX point-wise operator performing Y = X
  SUBROUTINE ASSIGNMENTYX(numDim,numPoints,bufferSize,bufferInterval,X,Y)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints
    INTEGER(KIND=8), INTENT(IN)      :: bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: Y(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = X(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Y(bufferIndex) = X(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Y(bufferIndex) = X(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE ASSIGNMENTYX

  !> @brief ASSIGNMENTXA point-wise operator performing X = scalar a
  SUBROUTINE ASSIGNMENTXA(numDim,numPoints,bufferSize,bufferInterval,a,X)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints
    INTEGER(KIND=8), INTENT(IN)      :: bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(INOUT)   :: X(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd


    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          X(I) = a
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             X(bufferIndex) = a
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                X(bufferIndex) = a
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE ASSIGNMENTXA

  !> @brief ASSIGNMENTYABSX point-wise operator performing X = scalar a
  SUBROUTINE ASSIGNMENTYABSX(numDim,numPoints,bufferSize,bufferInterval,X,Y)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints
    INTEGER(KIND=8), INTENT(IN)      :: bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: X(numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: Y(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd


    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = ABS(X(I))
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Y(bufferIndex) = ABS(X(bufferIndex))
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Y(bufferIndex) = ABS(X(bufferIndex))
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE ASSIGNMENTYABSX

  !> @brief VECLEN point-wise operator returning the length of a numComp-dimensional vector
  SUBROUTINE VECLEN(numDim,numPoints,bufferSize,bufferInterval,numComp,V,lenV)

    INTEGER(KIND=4), INTENT(IN)      :: numDim, numComp
    INTEGER(KIND=8), INTENT(IN)      :: numPoints
    INTEGER(KIND=8), INTENT(IN)      :: bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: V(numComp*numPoints)
    REAL(KIND=8),    INTENT(OUT)     :: lenV(numPoints)

    INTEGER(KIND=8) :: I, J, K, L
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd


    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          lenV(I) = 0.0_8
          DO L = 1, numComp
             lenV(I) = lenV(I) + V((L-1)*numPoints+I)*V((L-1)*numPoints+I)
          END DO
          lenV(I) = SQRT(lenV(I))
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             lenV(bufferIndex) = 0.0_8
             DO L = 1,numComp
                lenV(bufferIndex) = lenV(bufferIndex) +  &
                     V((L-1)*numPoints+bufferIndex)*V((L-1)*numPoints+bufferIndex)
             END DO
             lenV(bufferIndex) = SQRT(lenV(bufferIndex))
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                lenV(bufferIndex) = 0.0_8
                DO L = 1,numComp
                   lenV(bufferIndex) = lenV(bufferIndex) +  &
                        V((L-1)*numPoints+bufferIndex)*V((L-1)*numPoints+bufferIndex)
                END DO
                lenV(bufferIndex) = SQRT(lenV(bufferIndex))
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE VECLEN
  

  !> @brief YAXM1 point-wise operator performing Y = a/X (scalar a)
  SUBROUTINE YAXM1(numDim,numPoints,bufferSize,bufferInterval,a,X,Y)

    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(INOUT)   :: X(numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: Y(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = a/X(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             Y(bufferIndex) = a/X(bufferIndex)
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                Y(bufferIndex) = a/X(bufferIndex)
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE YAXM1


  !> @brief VASUPV point-wise operator performing V = V + aW where
  !>        a is a constant scalar, S is a scalar field, and V and W 
  !>        are vectors with numDim components.
  SUBROUTINE VASUPV(numDim,numPoints,bufferSize,bufferInterval,a,S,U,V)
    
    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(IN)      :: S(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: U(numDim*numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: V(numDim*numPoints)

    INTEGER(KIND=4) :: iDim
    INTEGER(KIND=8) :: I, J, K,vectorIndex
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd
    REAL(KIND=8)    :: scalFac

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          V(I) = V(I) + a*S(I)*U(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             scalFac = a*S(bufferIndex)
             DO iDim = 1, numDim
                vectorIndex = (iDim-1)*numPoints+bufferIndex
                V(vectorIndex) = V(vectorIndex) + scalFac*U(vectorIndex)
             END DO
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                scalFac = a*S(bufferIndex)
                DO iDim = 1, numDim
                   vectorIndex = (iDim-1)*numPoints+bufferIndex
                   V(vectorIndex) = V(vectorIndex) + scalFac*U(vectorIndex)
                END DO
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE VASUPV
  
  !> @brief YASSMWDOTUPY point-wise operator performing Y = Y + a*(S1+S2)*(W dot U) where
  !>        W, and U are vector fields with numDim components, and a is
  !>        a constant scalar, and S1, S2, and Y are scalar fields
  SUBROUTINE YASSMWDOTUPY(numDim,numPoints,bufferSize,bufferInterval,a,S1,S2,W,U,Y)
    
    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(IN)      :: S1(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: S2(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: W(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)      :: U(numDim*numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: Y(numPoints)

    INTEGER(KIND=4) :: iDim
    INTEGER(KIND=8) :: I, J, K,vectorIndex
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd
    REAL(KIND=8)    :: scalFac

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = Y(I) + a*(S1(I)+S2(I))*W(I)*U(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I  
             scalFac = a*(S1(bufferIndex)+S2(bufferIndex))
             DO iDim = 1, numDim
                vectorIndex = (iDim-1)*numPoints+bufferIndex
                Y(bufferIndex) = Y(bufferIndex) + scalFac*W(vectorIndex)*U(vectorIndex)
             END DO
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                scalFac = a*(S1(bufferIndex)+S2(bufferIndex))
                DO iDim = 1, numDim
                   vectorIndex = (iDim-1)*numPoints+bufferIndex
                   Y(bufferIndex) = Y(bufferIndex) + scalFac*W(vectorIndex)*U(vectorIndex)
                END DO
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE YASSMWDOTUPY
  

  !> @brief YASMWDOTUPY point-wise operator performing Y = Y + a*S*(W dot U) where
  !>        W, and U are vector fields with numDim components, and a is
  !>        a constant scalar, and S,Y are scalar fields
  SUBROUTINE YASMWDOTUPY(numDim,numPoints,bufferSize,bufferInterval,a,S,W,U,Y)
    
    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(IN)      :: S(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: W(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)      :: U(numDim*numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: Y(numPoints)

    INTEGER(KIND=4) :: iDim
    INTEGER(KIND=8) :: I, J, K,vectorIndex
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd
    REAL(KIND=8)    :: scalFac

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = Y(I) + a*S(I)*W(I)*U(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I  
             scalFac = a*S(bufferIndex)
             DO iDim = 1, numDim
                vectorIndex = (iDim-1)*numPoints+bufferIndex
                Y(bufferIndex) = Y(bufferIndex) + scalFac*W(vectorIndex)*U(vectorIndex)
             END DO
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                scalFac = a*S(bufferIndex)
                DO iDim = 1, numDim
                   vectorIndex = (iDim-1)*numPoints+bufferIndex
                   Y(bufferIndex) = Y(bufferIndex) + scalFac*W(vectorIndex)*U(vectorIndex)
                END DO
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE YASMWDOTUPY

  !> @brief YAVDOTWPY point-wise operator performing Y = Y + a(V dot W) where
  !>        Y is a scalar field, a is a constant scalar, V and W are vector fields
  !>        with numDim components
  SUBROUTINE YAVDOTWPY(numDim,numPoints,bufferSize,bufferInterval,a,V,W,Y)
    
    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: a
    REAL(KIND=8),    INTENT(IN)      :: V(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)      :: W(numDim*numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: Y(numPoints)

    INTEGER(KIND=4) :: iDim
    INTEGER(KIND=8) :: I, J, K,vectorIndex
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)

    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          Y(I) = Y(I) + a*V(I)*W(I)  
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             DO iDim = 1, numDim
                vectorIndex = (iDim-1)*numPoints+bufferIndex
                Y(bufferIndex) = Y(bufferIndex) + a*V(vectorIndex)*W(vectorIndex)
             END DO
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                DO iDim = 1, numDim
                   vectorIndex = (iDim-1)*numPoints+bufferIndex
                   Y(bufferIndex) = Y(bufferIndex) + a*V(vectorIndex)*W(vectorIndex)
                END DO
             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE YAVDOTWPY

  !> @brief Computes determinant of 3x3 matrix
  SUBROUTINE DETERMINANT3X3(numPoints,bufferSize,bufferInterval,inMatrix,matrixDeterminant)
    
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(3)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(6)
    REAL(KIND=8),    INTENT(IN)      :: inMatrix(9*numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: matrixDeterminant(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)
    jStart = bufferInterval(3)
    jEnd   = bufferInterval(4)
    kStart = bufferInterval(5)
    kEnd   = bufferInterval(6)
    nPlane = xSize * bufferSize(2)

    DO K = kStart, kEnd
       zIndex = (K-1)*nPlane
       DO J = jStart, jEnd
          yzIndex = zIndex + (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yzIndex + I

             matrixDeterminant(bufferIndex) =                &
                  inMatrix(bufferIndex)*                     &
                  (inMatrix(bufferIndex+8*numPoints) * &
                  inMatrix(bufferIndex+4*numPoints)-   &
                  inMatrix(bufferIndex+7*numPoints)*   &
                  inMatrix(bufferIndex+5*numPoints))-  &

                  inMatrix(bufferIndex+3*numPoints)*   &
                  (inMatrix(bufferIndex+8*numPoints)*  &
                  inMatrix(bufferIndex+numPoints) -    &
                  inMatrix(bufferIndex+7*numPoints) *  &
                  inMatrix(bufferIndex+2*numPoints)) + &

                  inMatrix(bufferIndex+6*numPoints)*   &
                  (inMatrix(bufferIndex+5*numPoints)*  &
                  inMatrix(bufferIndex+numPoints)-     &
                  inMatrix(bufferIndex+4*numPoints)*   &
                  inMatrix(bufferIndex+2*numPoints))

          END DO
       END DO
    END DO
  END SUBROUTINE DETERMINANT3X3

  !> @brief Computes determinant of 2x2 matrix
  SUBROUTINE DETERMINANT2X2(numPoints,bufferSize,bufferInterval,inMatrix,matrixDeterminant)
    
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(2)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(4)
    REAL(KIND=8),    INTENT(IN)      :: inMatrix(4*numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: matrixDeterminant(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)
    jStart = bufferInterval(3)
    jEnd   = bufferInterval(4)
    nPlane = xSize * bufferSize(2)

    DO J = jStart, jEnd
       yIndex = (J-1)*xSize
       DO I = iStart, iEnd
          bufferIndex = yIndex + I
          
          matrixDeterminant(bufferIndex) =                &
               (inMatrix(bufferIndex)*                    &
               inMatrix(bufferIndex+3*numPoints)) -       &
               (inMatrix(bufferIndex+numPoints)*          &
               inMatrix(bufferIndex+2*numPoints))
       END DO
    END DO

  END SUBROUTINE DETERMINANT2X2

  !> @brief Computes buf1*buf4 - buf2*buf3 + buf7*(buf5 - buf6)
  SUBROUTINE METRICSUM4(numDim,numPoints,bufferSize,bufferInterval,&
       buf1,buf2,buf3,buf4,buf5,buf6,buf7,metricSum)
    
    INTEGER(KIND=4), INTENT(IN)      :: numDim
    INTEGER(KIND=8), INTENT(IN)      :: numPoints,bufferSize(numDim)
    INTEGER(KIND=8), INTENT(IN)      :: bufferInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)      :: buf1(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: buf2(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: buf3(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: buf4(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: buf5(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: buf6(numPoints)
    REAL(KIND=8),    INTENT(IN)      :: buf7(numPoints)
    REAL(KIND=8),    INTENT(INOUT)   :: metricSum(numPoints)

    INTEGER(KIND=8) :: I, J, K
    INTEGER(KIND=8) :: nPlane, zIndex, yIndex, yzIndex, bufferIndex, xSize
    INTEGER(KIND=8) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

    iStart = bufferInterval(1)
    iEnd   = bufferInterval(2)
    xSize  = bufferSize(1)
    
    IF(numDim == 1) THEN
       DO I = iStart, iEnd
          metricSum(I) = buf1(I)*buf4(I) - buf2(I)*buf3(I) + buf7(I)*(buf5(I)-buf6(I))
       END DO
    ELSE IF(numDim == 2) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       DO J = jStart, jEnd
          yIndex = (J-1)*xSize
          DO I = iStart, iEnd
             bufferIndex = yIndex + I
             metricSum(bufferIndex) = buf1(bufferIndex)*buf4(bufferIndex) - &
                  buf2(bufferIndex)*buf3(bufferIndex) + &
                  buf7(bufferIndex)*(buf5(bufferIndex)-buf6(bufferIndex))
          END DO
       END DO
    ELSE IF(numDim == 3) THEN
       jStart = bufferInterval(3)
       jEnd   = bufferInterval(4)
       kStart = bufferInterval(5)
       kEnd   = bufferInterval(6)
       nPlane = xSize * bufferSize(2)
       DO K = kStart, kEnd
          zIndex = (K-1)*nPlane
          DO J = jStart, jEnd
             yzIndex = zIndex + (J-1)*xSize
             DO I = iStart, iEnd
                bufferIndex = yzIndex + I
                metricSum(bufferIndex) = buf1(bufferIndex)*buf4(bufferIndex) - &
                     buf2(bufferIndex)*buf3(bufferIndex) + &
                     buf7(bufferIndex)*(buf5(bufferIndex)-buf6(bufferIndex))
             END DO
          END DO
       END DO
    ENDIF
  END SUBROUTINE METRICSUM4

  SUBROUTINE VECTORCROSSPRODUCT(v1,v2,y)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN)  :: v1(3), v2(3)
    REAL(KIND=8), INTENT(OUT) :: y(3)

    y(1) = v1(2)*v2(3) - v1(3)*v2(2)
    y(2) = v1(3)*v2(1) - v1(1)*v2(3)
    y(3) = v1(1)*v2(2) - v1(2)*v2(1)
    
  END SUBROUTINE VECTORCROSSPRODUCT

END MODULE OPERATORS
