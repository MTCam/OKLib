MODULE Euler

  USE OPERATORS
  USE GRID

  IMPLICIT NONE

CONTAINS

  SUBROUTINE UniformRHS(                                           &
       numDim, gridSizes, numPoints,                               &
       fullInterval, opInterval, gridMetric,                       & 
       numStencils, numStencilValues, stencilSizes, stencilStarts, &
       stencilOffsets, stencilWeights, stencilID,                  &
       rhoBuffer,rhoVBuffer,rhoEBuffer,velHat,                     &
       pressureBuffer, rhoRHS, rhoVRHS, rhoERHS)
    
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)         :: numDim, numStencils, numStencilValues
    INTEGER(KIND=8), INTENT(IN)         :: gridSizes(numDim),  numPoints
    INTEGER(KIND=8), INTENT(IN)         :: fullInterval(2*numDim), opInterval(2*numDim)
    INTEGER(KIND=4), INTENT(IN)         :: stencilSizes(numStencils),stencilStarts(numStencils)
    INTEGER(KIND=4), INTENT(IN)         :: stencilOffsets(numStencilValues)
    INTEGER(KIND=4), INTENT(IN), TARGET :: stencilID(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: stencilWeights(numStencilValues)
    REAL(KIND=8),    INTENT(IN)         :: gridMetric(numDim)
    REAL(KIND=8),    INTENT(IN), TARGET :: velHat(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: rhoBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: rhoVBuffer(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: rhoEBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN)         :: pressureBuffer(numPoints)
    REAL(KIND=8),    INTENT(OUT)        :: rhoRHS(numPoints)
    REAL(KIND=8),    INTENT(OUT),TARGET :: rhoVRHS(numDim*numPoints)
    REAL(KIND=8),    INTENT(OUT)        :: rhoERHS(numPoints)


    INTEGER         :: iDim, numComponents, iVelDim
    INTEGER(KIND=8) :: iPoint,iX,iY,iZ,iStart,iEnd,jStart,jEnd,kStart,kEnd
    INTEGER(KIND=8) :: xIndex,zIndex,yIndex,yzIndex,xSize,ySize,zSize
    INTEGER(KIND=8) :: iPoint2, iPoint3, pointOffset, pointIndex,vectorPointIndex
    INTEGER(KIND=8) :: velIndex, fluxIndex, velOffset
    REAL(KIND=8)    :: gridScale

    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: fluxBuffer, dFluxBuffer, scaledPressure
    REAL(KIND=8),    DIMENSION(:), POINTER     :: rhoVPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: rhoVRHSPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: bufPtr
    INTEGER(KIND=4), DIMENSION(:), POINTER     :: stencilConnectivity

    ALLOCATE(fluxBuffer(numPoints))
    ALLOCATE(dFluxBuffer(numPoints))
    ALLOCATE(scaledPressure(numPoints))

    numComponents = 1

    ! Continuity 
    DO iDim = 1,numDim

       pointOffset = (iDim-1)*numPoints

       stencilConnectivity => stencilID((pointOffset+1):(pointOffset+numPoints))

       ! Calculates fluxBuffer = rho * vHat(iDim) [full Interval]
       bufPtr => velHat((pointOffset+1):(pointOffset+numPoints))
       CALL ZXY(numDim,numPoints,gridSizes,fullInterval,rhoBuffer,bufPtr,fluxBuffer)
       !       WRITE(*,*) 'FFLUXBUFFERRHO',iDim,fluxBuffer
       ! # need OMP barrier here to fix OMP issue?
       ! Deriv in X(n) direction of the flux:   d/dX(iDim) * [fluxBuffer]
       CALL APPLYOPERATOR(numDim, gridSizes, numComponents, numPoints, iDim, opInterval, &
            numStencils, stencilSizes, stencilStarts, numStencilValues,                  &
            stencilWeights, stencilOffsets, stencilConnectivity, fluxBuffer, dFluxBuffer)
       
       ! Sums dFlux into RHS [only on local interval]
       ! rhoRHS = rhoRHS - dFlux
       CALL YAXPY(numDim,numPoints,gridSizes,opInterval,-1.0_8,dFluxBuffer,rhoRHS)
        !       Call Save_Patch_Deriv(region, ng, 1, i, dflux)

    END DO ! iDim (Continuity)

    ! ... momentum
    DO iDim = 1, numDim
       
       pointOffset =  (iDim-1)*numPoints
       gridScale   =  gridMetric(iDim)
       rhoVPtr     => rhoVBuffer((pointOffset+1):(pointOffset+numPoints))
       rhoVRHSPtr  => rhoVRHS((pointOffset+1):(pointOffset+numPoints))
       
       ! Gets gridscaled pressure for diagonal (scaledPressure = gridScale*pressure)
       CALL YAX(numDim,numPoints,gridSizes,fullInterval,gridScale,pressureBuffer,scaledPressure)
       
       DO iVelDim = 1, numDim
          
          velOffset = (iVelDim-1)*numPoints
          bufPtr => velHat((velOffset+1):(velOffset+numPoints))  
          stencilConnectivity => stencilID((velOffset+1):(velOffset+numPoints))
          
          IF(iDim == iVelDim) THEN
             ! Diagonal term: flux = rhoV(iDim)*vHat(iDim) + scaledPressure
             CALL ZWXPY(numDim,numPoints,gridSizes,fullInterval,rhoVPtr,bufPtr,scaledPressure,fluxBuffer)
          ELSE
             ! Cross terms: flux = rhoV(iDim)*vHat(iVelDim)
             CALL ZXY(numDim,numPoints,gridSizes,fullInterval,rhoVPtr,bufPtr,fluxBuffer)
          ENDIF
          !          WRITE(*,*) 'FFLUXBUFFERRHOV',iDim,iVelDim,fluxBuffer
          ! # need OMP barrier here 
          
          CALL APPLYOPERATOR(numDim, gridSizes, numComponents, numPoints, iVelDim, opInterval, &
               numStencils, stencilSizes, stencilStarts, numStencilValues,                     &
               stencilWeights, stencilOffsets, stencilConnectivity, fluxBuffer, dFluxBuffer)
          
          ! Sum dFlux to RHS (using Y = a*X + Y)
          CALL YAXPY(numDim,numPoints,gridSizes,opInterval,-1.0_8,dFluxBuffer,rhoVRHSPtr)
          !        Call Save_Patch_Deriv(region, ng, i+1, j, dflux)
       END DO ! iVelDim
    END DO ! iDim (momentum)

    
    ! ... energy
    DO iDim = 1, numDim
       
       velOffset =  (iDim-1)*numPoints
       bufPtr => velHat((velOffset+1):(velOffset+numPoints))
       
       stencilConnectivity => stencilID((velOffset+1):(velOffset+numPoints))
       
       ! Calculate flux = vHat(iDim)*(rhoE + pressure)
       CALL ZWMXPY(numDim,numPoints,gridSizes,fullInterval,bufPtr,rhoEBuffer,pressureBuffer,fluxBuffer)
!       WRITE(*,*) 'FFLUXBUFFERRHOE',iDim,fluxBuffer
       
       ! # Need OMP barrier
       CALL APPLYOPERATOR(numDim, gridSizes, numComponents, numPoints, iDim, opInterval, &
            numStencils, stencilSizes, stencilStarts, numStencilValues,                  &
            stencilWeights, stencilOffsets, stencilConnectivity, fluxBuffer, dFluxBuffer)

       CALL YAXPY(numDim,numPoints,gridSizes,opInterval,-1.0_8,dFluxBuffer,rhoERHS)
          !        Call Save_Patch_Deriv(region, ng, i+1, j, dflux)
    END DO ! iDim (energy)


    DEALLOCATE(scaledPressure)
    DEALLOCATE(fluxBuffer)
    DEALLOCATE(dFluxBuffer)

  END SUBROUTINE UniformRHS



  SUBROUTINE UniformFlux(                                          &
       numDim, fluxDim, gridSizes, numPoints,                      &
       opInterval, gridMetric,                                     & 
       rhoBuffer,rhoVBuffer,rhoEBuffer,velHat,                     &
       pressureBuffer, scaledPressure, fluxBuffer)
    
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)         :: numDim, fluxDim
    INTEGER(KIND=8), INTENT(IN)         :: gridSizes(numDim),  numPoints
    INTEGER(KIND=8), INTENT(IN)         :: opInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)         :: gridMetric(numDim)
    REAL(KIND=8),    INTENT(IN)         :: rhoBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: rhoVBuffer(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: rhoEBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: velHat(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: pressureBuffer(numPoints)
    REAL(KIND=8),    INTENT(OUT)        :: scaledPressure(numPoints)
    REAL(KIND=8),    INTENT(OUT),TARGET :: fluxBuffer(numPoints*(numDim+2))

    INTEGER         :: iDim, numComponents, iVelDim
    INTEGER(KIND=8) :: iPoint,iX,iY,iZ,iStart,iEnd,jStart,jEnd,kStart,kEnd
    INTEGER(KIND=8) :: xIndex,zIndex,yIndex,yzIndex,xSize,ySize,zSize,fluxOffset
    INTEGER(KIND=8) :: iPoint2, iPoint3, pointOffset, pointIndex,vectorPointIndex
    INTEGER(KIND=8) :: velIndex, fluxIndex, velOffset, dimOffset
    REAL(KIND=8)    :: gridScale

    REAL(KIND=8),    DIMENSION(:), POINTER     :: rhoVPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: rhoVRHSPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: bufPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: velHatPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: fluxPtr
    INTEGER(KIND=4), DIMENSION(:), POINTER     :: stencilConnectivity

    numComponents = 1

    iDim = fluxDim
    dimOffset  =  (iDim-1)*numPoints
    fluxOffset =  0
    velHatPtr  => velHat(dimOffset+1:dimOffset+numPoints)
    fluxPtr    => fluxBuffer(1:numPoints)
    gridScale  =  gridMetric(iDim)


    ! Continuity 
    ! Calculates fluxBuffer = rho * vHat(iDim) [full Interval]
    CALL ZXY(numDim,numPoints,gridSizes,opInterval,rhoBuffer,velHatPtr,fluxPtr)
    !    WRITE(*,*) 'FFLUXBUFFERRHO',iDim,fluxPtr
    ! numDim components of Momentum
    ! Gets gridscaled pressure for diagonal (scaledPressure = gridScale*pressure)
    CALL YAX(numDim,numPoints,gridSizes,opInterval,gridScale,pressureBuffer,scaledPressure)       

    DO iVelDim = 1, numDim

       fluxOffset =  fluxOffset + numPoints
       velOffset  =  (iVelDim-1)*numPoints
       fluxPtr    => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)
       rhoVPtr    => rhoVBuffer((velOffset+1):(velOffset+numPoints))    

       IF(iDim == iVelDim) THEN
          ! Diagonal term: flux = rhoV(iDim)*vHat(iDim) + scaledPressure
          CALL ZWXPY(numDim,numPoints,gridSizes,opInterval,rhoVPtr,velHatPtr,scaledPressure,fluxPtr)
       ELSE
          ! Cross terms: flux = rhoV(iVelDim)*vHat(iDim)
          CALL ZXY(numDim,numPoints,gridSizes,opInterval,rhoVPtr,velHatPtr,fluxPtr)
       ENDIF
       !       WRITE(*,*) 'FFLUXBUFFERRHOV',iVelDim,iDim,fluxPtr
    END DO

    ! Energy
    fluxOffset = fluxOffset + numPoints
    fluxPtr => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)       

    ! Calculate flux = vHat(iDim)*(rhoE + pressure)
    CALL ZWMXPY(numDim,numPoints,gridSizes,opInterval,velHatPtr,rhoEBuffer,pressureBuffer,fluxPtr)
           !WRITE(*,*) 'FFLUXBUFFERRHOE Euler',iDim,fluxPtr
    
  END SUBROUTINE UniformFlux

  SUBROUTINE Flux1D(                                &
       numDim, numPoints, gridSizes, opInterval,    &
       fluxDir, gridType, gridMetric,               & 
       rhoBuffer,rhoVBuffer,rhoEBuffer,velHat,   &
       pressureBuffer, fluxBuffer)
    
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)         :: numDim, fluxDir, gridType
    INTEGER(KIND=8), INTENT(IN)         :: gridSizes(numDim),  numPoints
    INTEGER(KIND=8), INTENT(IN)         :: opInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN), TARGET :: gridMetric(numDim*numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: rhoBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: rhoVBuffer(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: rhoEBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN)         :: velHat(numPoints)
    REAL(KIND=8),    INTENT(IN)         :: pressureBuffer(numPoints)
    REAL(KIND=8),    INTENT(OUT),TARGET :: fluxBuffer(numPoints*(numDim+2))

    INTEGER         :: iDim, numComponents, iVelDim
    INTEGER(KIND=8) :: iPoint,iX,iY,iZ,iStart,iEnd,jStart,jEnd,kStart,kEnd
    INTEGER(KIND=8) :: xIndex,zIndex,yIndex,yzIndex,xSize,ySize,zSize,fluxOffset
    INTEGER(KIND=8) :: iPoint2, iPoint3, pointOffset, pointIndex,vectorPointIndex
    INTEGER(KIND=8) :: velIndex, fluxIndex, velOffset, dimOffset, metricOffset
    REAL(KIND=8)    :: gridScale

    REAL(KIND=8),    DIMENSION(:), POINTER     :: rhoVPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: rhoVRHSPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: fluxPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: metricPtr
    INTEGER(KIND=4), DIMENSION(:), POINTER     :: stencilConnectivity


    fluxOffset =  0
    fluxPtr    => fluxBuffer(1:numPoints)

    ! Continuity 
    ! Calculates fluxBuffer = rho * vHat(iDim)
    CALL ZXY(numDim,numPoints,gridSizes,opInterval,rhoBuffer,velHat,fluxPtr)

    ! numDim components of Momentum
    IF(gridType < RECTILINEAR) THEN

       gridScale = gridMetric(fluxDir)

       DO iVelDim = 1, numDim
          
          fluxOffset =  fluxOffset + numPoints
          velOffset  =  (iVelDim-1)*numPoints
          fluxPtr    => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)
          rhoVPtr    => rhoVBuffer((velOffset+1):(velOffset+numPoints))    
          
          IF(fluxDir == iVelDim) THEN
             ! gridMetric(iDim) is single constant over all points (0 for other dimensions)
             ! Diagonal term: flux = rhoV(iDim)*vHat(iDim) + gridMetric(iDim)*pressure
             CALL ZAWPXY(numDim,numPoints,gridSizes,opInterval,gridScale,pressureBuffer,&
                  rhoVPtr,velHat,fluxPtr)
          ELSE
             ! Cross terms: flux = rhoV(iVelDim)*vHat(iDim)
             CALL ZXY(numDim,numPoints,gridSizes,opInterval,rhoVPtr,velHat,fluxPtr)
          ENDIF
       END DO

    ELSE IF(gridType < CURVILINEAR) THEN

       metricOffset = (fluxDir-1)*numPoints
       metricPtr => gridMetric(metricOffset+1:metricOffset+numPoints)

       DO iVelDim = 1, numDim
          
          fluxOffset =  fluxOffset + numPoints
          velOffset  =  (iVelDim-1)*numPoints
          fluxPtr    => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)
          rhoVPtr    => rhoVBuffer((velOffset+1):(velOffset+numPoints))    
          
          IF(fluxDir == iVelDim) THEN
             ! gridMetric(iDim) is nodal scalar (0 for other dimensions)
             ! Diagonal term: flux = rhoV(iDim)*vHat(iDim) + gridMetric(iDim)*pressure
             CALL ZVWPXY(numDim,numPoints,gridSizes,opInterval,metricPtr,pressureBuffer,rhoVPtr,velHat,fluxPtr)
          ELSE
             ! Cross terms: flux = rhoV(iVelDim)*vHat(iDim)
             CALL ZXY(numDim,numPoints,gridSizes,opInterval,rhoVPtr,velHat,fluxPtr)
          ENDIF

       END DO

    ELSE ! gridMetric(iDim) is a nodal numDim-vector

       metricOffset = (fluxDir-1)*numDim*numPoints
       
       DO iVelDim = 1, numDim
          
          fluxOffset =  fluxOffset + numPoints
          velOffset  =  (iVelDim-1)*numPoints
          fluxPtr    => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)
          rhoVPtr    => rhoVBuffer((velOffset+1):(velOffset+numPoints))
          metricPtr  => gridMetric(metricOffset+velOffset+1:metricOffset+velOffset+numPoints)

          ! flux = rhoV(iDim)*vHat(iDim) + p*gridMetric(iDim)  
          CALL ZVWPXY(numDim,numPoints,gridSizes,opInterval,metricPtr,pressureBuffer,rhoVPtr,velHat,fluxPtr)

       END DO

    ENDIF
    
    ! Energy
    fluxOffset = fluxOffset + numPoints
    fluxPtr => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)       

    ! Calculate flux = vHat(iDim)*(rhoE + pressure)
    CALL ZWMXPY(numDim,numPoints,gridSizes,opInterval,velHat,rhoEBuffer,pressureBuffer,fluxPtr)
    
  END SUBROUTINE Flux1D

  SUBROUTINE UniformScalarRHS(                                           &
       numDim, gridSizes, numPoints, opInterval,                         &
       numStencils, numStencilValues, stencilSizes, stencilStarts,       &
       stencilOffsets, stencilWeights, stencilID,                        &
       numPointsApply,applyPoints,numScalar,                          &
       scalarBuffer,velHat,scalarRHS)

    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)         :: numDim, numStencils, numStencilValues, numScalar
    INTEGER(KIND=8), INTENT(IN)         :: gridSizes(numDim),  numPoints, numPointsApply
    INTEGER(KIND=8), INTENT(IN)         :: opInterval(2*numDim), applyPoints(numPointsApply)
    INTEGER(KIND=4), INTENT(IN)         :: stencilSizes(numStencils),stencilStarts(numStencils)
    INTEGER(KIND=4), INTENT(IN)         :: stencilOffsets(numStencilValues)
    INTEGER(KIND=4), INTENT(IN), TARGET :: stencilID(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: stencilWeights(numStencilValues)
    REAL(KIND=8),    INTENT(IN)         :: velHat(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: scalarBuffer(numScalar*numPoints)
    REAL(KIND=8),    INTENT(OUT)        :: scalarRHS(numScalar*numPoints)


    INTEGER         :: iDim, numComponents, iVelDim, iScalar
    INTEGER(KIND=8) :: iPoint,iX,iY,iZ,iStart,iEnd,jStart,jEnd,kStart,kEnd
    INTEGER(KIND=8) :: xIndex,zIndex,yIndex,yzIndex,xSize,ySize,zSize
    INTEGER(KIND=8) :: iPoint2, iPoint3, pointOffset, pointIndex,scalarIndex
    INTEGER(KIND=8) :: velIndex, fluxIndex, velOffset, scalarOffset


    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: fluxBuffer, dFluxBuffer
    INTEGER(KIND=4), DIMENSION(:), POINTER     :: stencilConnectivity

    ALLOCATE(fluxBuffer(numPoints))
    ALLOCATE(dFluxBuffer(numPoints))

    numComponents = 1

    ! Scalar scalar

    DO iScalar = 1, numScalar

       scalarOffset = (iScalar-1)*numPoints

       DO iDim = 1,numDim

          pointOffset = (iDim-1)*numPoints

          stencilConnectivity => stencilID((pointOffset+1):(pointOffset+numPoints))

          DO iPoint = 1, numPointsApply
             fluxIndex = applyPoints(iPoint)
             scalarIndex = fluxIndex + scalarOffset
             fluxBuffer(fluxIndex) = scalarBuffer(scalarIndex) * velHat(fluxIndex+pointOffset)
          END DO

          CALL APPLYOPERATOR(numDim, gridSizes, numComponents, numPoints, iDim, opInterval, &
               numStencils, stencilSizes, stencilStarts, numStencilValues,                  &
               stencilWeights, stencilOffsets, stencilConnectivity, fluxBuffer, dFluxBuffer)

          DO iPoint = 1, numPointsApply
             fluxIndex = applyPoints(iPoint)
             scalarIndex = fluxIndex + scalarOffset
             scalarRHS(scalarIndex) = scalarRHS(scalarIndex) - dFluxBuffer(fluxIndex)
          END DO
          !       Call Save_Patch_Deriv(region, ng, 1, i, dflux)
       END DO ! iDim
    END DO ! iScalar

    DEALLOCATE(fluxBuffer)
    DEALLOCATE(dFluxBuffer)

  END SUBROUTINE UniformScalarRHS


  SUBROUTINE UniformScalarFlux(                                    & 
       numDim, fluxDim, gridSizes, numPoints, opInterval,          &
       numScalar,scalarBuffer,velHat,scalarFlux)

    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)         :: numDim, fluxDim, numScalar
    INTEGER(KIND=8), INTENT(IN)         :: gridSizes(numDim),  numPoints
    INTEGER(KIND=8), INTENT(IN)         :: opInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN), TARGET :: velHat(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: scalarBuffer(numScalar*numPoints)
    REAL(KIND=8),    INTENT(OUT),TARGET :: scalarFlux(numScalar*numPoints)


    INTEGER         :: iDim, numComponents, iVelDim, iScalar
    INTEGER(KIND=8) :: iPoint,iX,iY,iZ,iStart,iEnd,jStart,jEnd,kStart,kEnd
    INTEGER(KIND=8) :: xIndex,zIndex,yIndex,yzIndex,xSize,ySize,zSize
    INTEGER(KIND=8) :: iPoint2, iPoint3, pointOffset, pointIndex,scalarIndex
    INTEGER(KIND=8) :: velIndex, fluxIndex, velOffset, scalarOffset

    REAL(KIND=8), DIMENSION(:), POINTER     :: fluxPtr,scalarPtr,velHatPtr


    numComponents = 1

    ! Scalar scalar flux

    iDim = fluxDim
    pointOffset =  (iDim-1)*numPoints
    velHatPtr   => velHat(pointOffset+1:pointOffset+numPoints)
       
    DO iScalar = 1, numScalar
       
       scalarOffset = (iScalar-1)*numPoints
       
       fluxPtr       => scalarFlux(scalarOffset+1:scalarOffset+numPoints)
       scalarPtr  => scalarBuffer(scalarOffset+1:scalarOffset+numPoints)
       
       CALL ZXY(numDim,numPoints,gridSizes,opInterval,scalarPtr,velHatPtr,fluxPtr)

    END DO ! iScalar

  END SUBROUTINE UniformScalarFlux

  !> @brief Flux for scalar transport
  SUBROUTINE ScalarFlux1D(                          &
       numDim, numPoints, gridSizes, opInterval,    &
       numScalars,scalarBuffer,velHat, fluxBuffer)
    
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)         :: numDim, numScalars
    INTEGER(KIND=8), INTENT(IN)         :: gridSizes(numDim),  numPoints
    INTEGER(KIND=8), INTENT(IN)         :: opInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN), TARGET :: scalarBuffer(numScalars*numPoints)
    REAL(KIND=8),    INTENT(IN)         :: velHat(numPoints)
    REAL(KIND=8),    INTENT(OUT),TARGET :: fluxBuffer(numScalars*numPoints)

    INTEGER         :: iScalar
    INTEGER(KIND=8) :: scalarOffset

    REAL(KIND=8),    DIMENSION(:), POINTER     :: fluxPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: scalarPtr


    scalarOffset =  0
    fluxPtr    => fluxBuffer(1:numPoints)

    ! Scalar transport
    ! Calculates fluxBuffer = scalar * vHat
    DO iScalar = 1, numScalars
       scalarPtr => scalarBuffer(scalarOffset+1:scalarOffset+numPoints)
       fluxPtr   => fluxBuffer(scalarOffset+1:scalarOffset+numPoints)
       CALL ZXY(numDim,numPoints,gridSizes,opInterval,scalarPtr,velHat,fluxPtr)
       scalarOffset = scalarOffset + numPoints
    END DO
    
  END SUBROUTINE ScalarFlux1D

END MODULE Euler
