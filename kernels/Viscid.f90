MODULE Viscid

  USE OPERATORS
  USE GRID

  IMPLICIT NONE

CONTAINS

  !
  ! Compute the viscous fluxes in 1 dimension is a uniform cartesian
  ! coordinate system
  !
  ! obsolete: see StrongFlux1D, which is generalized for curvilinear coordiantes systems
  SUBROUTINE ViscidStrongUniformFlux(                                     &
       numDim, fluxDim, gridSizes, numPoints, opInterval, velocity,       &
       gridMetric, tauOneBuffer, tauTwoBuffer, tauThreeBuffer, scaledTau, &
       heatFluxBuffer, fluxBuffer)
    
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)         :: numDim, fluxDim
    INTEGER(KIND=8), INTENT(IN)         :: gridSizes(numDim),  numPoints
    INTEGER(KIND=8), INTENT(IN)         :: opInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN)         :: gridMetric(numDim)
    REAL(KIND=8),    INTENT(IN), TARGET :: velocity(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: tauOneBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: tauTwoBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: tauThreeBuffer(numPoints)
    REAL(KIND=8),    INTENT(OUT)        :: scaledTau(numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: heatFluxBuffer(numPoints*numDim)
    REAL(KIND=8),    INTENT(OUT),TARGET :: fluxBuffer(numPoints*(numDim+2))

    INTEGER         :: iDim, numComponents, iVelDim
    INTEGER(KIND=8) :: iPoint,iX,iY,iZ,iStart,iEnd,jStart,jEnd,kStart,kEnd
    INTEGER(KIND=8) :: xIndex,zIndex,yIndex,yzIndex,xSize,ySize,zSize,fluxOffset
    INTEGER(KIND=8) :: iPoint2, iPoint3, pointOffset, pointIndex,vectorPointIndex
    INTEGER(KIND=8) :: velIndex, fluxIndex, velOffset, heatOffset, dimOffset
    REAL(KIND=8)    :: gridScale
    REAL(KIND=8)    :: minusOne

    REAL(KIND=8),    DIMENSION(:), POINTER     :: bufPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: velocityPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: fluxPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: tauOnePtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: tauTwoPtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: tauThreePtr
    REAL(KIND=8),    DIMENSION(:), POINTER     :: heatFluxPtr

    numComponents = 1
    minusOne = 1.0_8

    iDim = fluxDim
    fluxOffset =  0
    velocityPtr  => velocity(1:numPoints)
    heatFluxPtr => heatFluxBuffer((fluxDim-1)*numPoints+1:numPoints+(fluxDim-1)*numPoints)
    gridScale  =  gridMetric(iDim)
    fluxPtr    => fluxBuffer(1:numPoints)
    tauOnePtr  => tauOneBuffer(1:numPoints)
    tauTwoPtr  => tauTwoBuffer(1:numPoints)
    tauThreePtr=> tauThreeBuffer(1:numPoints)

    DO iVelDim = 1, numDim

       !velOffset  =  (iVelDim-1)*numPoints
       ! skip density
       fluxOffset =  fluxOffset + numPoints
       fluxPtr    => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)
       ! flux = tau(iDim)
       ! assuming tau(1,2,3) are the j components in tau_ij
       !gridScale = minusOne*gridMetric(iDim)*gridMetric(iVelDim)
       if(iVelDim == 1) THEN
         ! Gets grid (metric) scaled stress tensor (scaledPressure = gridScale*pressure)
         CALL YAX(numDim,numPoints,gridSizes,opInterval,gridScale,tauOnePtr,fluxPtr)       
         !CALL YAX(numDim,numPoints,gridSizes,opInterval,minusOne,tauOnePtr,fluxPtr)       
       ELSEIF(iVelDim == 2) THEN
         CALL YAX(numDim,numPoints,gridSizes,opInterval,gridScale,tauTwoPtr,fluxPtr)       
         !CALL YAX(numDim,numPoints,gridSizes,opInterval,minusOne,tauTwoPtr,fluxPtr)       
       ELSEIF(iVelDim == 3) THEN
         CALL YAX(numDim,numPoints,gridSizes,opInterval,gridScale,tauThreePtr,fluxPtr)       
         !CALL YAX(numDim,numPoints,gridSizes,opInterval,minusOne,tauThreePtr,fluxPtr)       
       ENDIF
              !WRITE(*,*) 'FFLUXBUFFERRHOV',iVelDim,iDim,fluxPtr
    END DO

    ! Energy
    fluxOffset =  fluxOffset + numPoints
    fluxPtr => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)  
    ! Calculate flux = vHat(iVelDim)*(tau(iVelDim,iDim))
    !WRITE(*,*) 'Before ZXY scaledTau',iDim,scaledTau
    !WRITE(*,*) 'Before ZXY velocityPtr',iDim,velocityPtr
    !WRITE(*,*) 'Before ZXY tauOneBuffer',iDim,tauOneBuffer

    ! viscous dissipation
    CALL ZXY(numDim,numPoints,gridSizes,opInterval,velocityPtr,tauOneBuffer,scaledTau)
    !WRITE(*,*) 'After ZXY scaledTau',iDim,scaledTau
    velocityPtr  => velocity(numPoints+1:2*numPoints)
    CALL ZWXPY(numDim,numPoints,gridSizes,opInterval,velocityPtr,tauTwoBuffer,scaledTau,scaledTau)
    IF(numDim == 3) THEN
      velocityPtr  => velocity(2*numPoints+1:3*numPoints)
      CALL ZWXPY(numDim,numPoints,gridSizes,opInterval,velocityPtr,tauThreeBuffer,scaledTau,scaledTau)
    ENDIF

    ! apply the grid metric
    !gridScale  =  minusOne*gridMetric(iDim)
             !WRITE(*,*) 'Before FFLUXBUFFERRHOE',iDim,fluxPtr(555)
    CALL YAX(numDim,numPoints,gridSizes,opInterval,gridScale,scaledTau,fluxPtr)       
             !WRITE(*,*) 'After FFLUXBUFFERRHOE Viscid',iDim,fluxPtr(555)

    ! heat conduction
    CALL YAXPY(numDim,numPoints,gridSizes,opInterval,gridScale,heatFluxPtr,fluxPtr)       
             !WRITE(*,*) 'After FFLUXBUFFERRHOE Conduction',iDim,fluxPtr(555)

  END SUBROUTINE ViscidStrongUniformFlux

  !>
  !! @brief Compute the curvilinear cartesian viscous fluxes in 1 dimension
  !!
  !! retains logic to simplify calculations on uniform or stretched cartesian grids
  !!
  !! See /ref conserve for theory
  !!
  SUBROUTINE StrongFlux1D                                            &
       (numDim, fluxDir, gridSizes, numPoints, opInterval, gridType, & 
       gridMetric, tauBuffer, energyBuffer, fluxBuffer)
    
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)         :: numDim, fluxDir, gridType
    INTEGER(KIND=8), INTENT(IN)         :: gridSizes(numDim),numPoints
    INTEGER(KIND=8), INTENT(IN)         :: opInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN), TARGET :: gridMetric(numDim*numDim*numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: tauBuffer(numPoints*numDim*(numDim+1)/2)
    REAL(KIND=8),    INTENT(IN), TARGET :: energyBuffer(numPoints*numDim)
    REAL(KIND=8),    INTENT(OUT),TARGET :: fluxBuffer(numPoints*(numDim+2))

    INTEGER         :: iDim, tensorIndex, iVel
    INTEGER(KIND=8) :: fluxOffset, tensorOffset, dirOffset, metricOffset
    INTEGER(KIND=8) :: velIndex, fluxIndex, velOffset, heatOffset, dimOffset
    REAL(KIND=8)    :: gridScale
    REAL(KIND=8)    :: minusOne

    REAL(KIND=8),    DIMENSION(:), POINTER :: bufPtr
    REAL(KIND=8),    DIMENSION(:), POINTER :: metricPtr
    REAL(KIND=8),    DIMENSION(:), POINTER :: fluxPtr
    REAL(KIND=8),    DIMENSION(:), POINTER :: tauPtr
    REAL(KIND=8),    DIMENSION(:), POINTER :: energyPtr

    INTEGER, DIMENSION(2,2) :: index2D
    INTEGER, DIMENSION(3,3) :: index3D

    ! Map 0-based index for (2,3)-dimensional symmetric tensor
    index2D = RESHAPE((/ 0, 1, 1, 2 /), SHAPE(index2D))
    index3D = RESHAPE((/ 0, 1, 2, 1, 3, 4, 2, 4, 5 /), SHAPE(index3D))

    dirOffset = (fluxDir-1)*numPoints
    fluxOffset = 0

    IF(gridType < RECTILINEAR) THEN

       gridScale  =  gridMetric(fluxDir)
       
       ! Momentum terms 
       DO iDim = 1, numDim
          
          IF(numDim == 2) THEN
             tensorIndex = index2D(fluxDir,iDim)
          ELSE
             tensorIndex = index3D(fluxDir,iDim)
          ENDIF
          
          tensorOffset =  tensorIndex*numPoints
          fluxOffset   =  fluxOffset + numPoints
          fluxPtr      => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)
          tauPtr       => tauBuffer(tensorOffset+1:tensorOffset+numPoints)

          ! Grid metric is scalar and constant, this gets flux = metric * tau_(iDim)(
          CALL YAX(numDim,numPoints,gridSizes,opInterval,gridScale,tauPtr,fluxPtr)       

       END DO
       
       ! Energy
       fluxOffset  =  fluxOffset + numPoints
       fluxPtr     => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)  
       energyPtr   => energyBuffer(dirOffset+1:dirOffset+numPoints)

       CALL YAX(numDim,numPoints,gridSizes,opInterval,gridScale,energyPtr,fluxPtr)

    ELSE IF (gridType == RECTILINEAR) THEN
       
       metricPtr    => gridMetric(dirOffset+1:dirOffset+numPoints)

       ! Momentum terms 
       DO iDim = 1, numDim
          
          IF(numDim == 2) THEN
             tensorIndex = index2D(fluxDir,iDim)
          ELSE
             tensorIndex = index3D(fluxDir,iDim)
          ENDIF
          
          tensorOffset =  tensorIndex*numPoints
          fluxOffset   =  fluxOffset + numPoints
          fluxPtr      => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)
          tauPtr       => tauBuffer(tensorOffset+1:tensorOffset+numPoints)

          ! Grid metric is scalar, this gets flux = metric * tau_(iDim)(
          CALL ZXY(numDim,numPoints,gridSizes,opInterval,metricPtr,tauPtr,fluxPtr)       

       END DO
       
       ! Energy
       fluxOffset  =  fluxOffset + numPoints
       fluxPtr     => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)  
       energyPtr   => energyBuffer(dirOffset+1:dirOffset+numPoints)

       CALL ZXY(numDim,numPoints,gridSizes,opInterval,metricPtr,energyPtr,fluxPtr)

    ELSE ! must be curvilinear
       
       ! need row 'fluxDir' of the metric tensor
       metricOffset = (fluxDir-1)*numDim*numPoints

       ! For curvilinear
       ! Momentum terms flux_i  = (metric_fluxDir_j * tau_j_i) i=[1:numDim]
       DO iDim = 1, numDim
          
          fluxOffset   =  fluxOffset + numPoints
          fluxPtr      => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)
          
          CALL ASSIGNMENTXA(numDim,numPoints,gridSizes,opInterval,0.0_8,fluxPtr)
          
          DO iVel = 1, numDim
             
             velOffset = (iVel-1)*numPoints+metricOffset
             metricPtr => gridMetric(velOffset+1:velOffset+numPoints)

             IF(numDim == 2) THEN
                tensorIndex = index2D(iVel,iDim)
             ELSE
                tensorIndex = index3D(iVel,iDim)
             ENDIF
          
             tensorOffset =  tensorIndex*numPoints
             tauPtr       => tauBuffer(tensorOffset+1:tensorOffset+numPoints)

             ! Get fluxPtr += (metric_fluxDir_iVel * tau_iVel_iDim) 
             CALL YWXPY(numDim,numPoints,gridSizes,opInterval,metricPtr,tauPtr,fluxPtr)
             
          END DO

       END DO
       
       ! Energy flux is just (metric_fluxDir_i * Q_i), since 
       ! both have proper vector storage the DOT operator 
       ! can be used
       fluxOffset  =  fluxOffset + numPoints
       fluxPtr     => fluxBuffer(fluxOffset+1:fluxOffset+numPoints)  
       metricPtr   => gridMetric(metricOffset+1:metricOffset+numDim*numPoints)

       CALL ZXDOTY(numDim,numPoints,gridSizes,opInterval,numDim,metricPtr,energyBuffer,fluxPtr)

    ENDIF

  END SUBROUTINE StrongFlux1D

  !>
  !! @brief Compute the curvilinear cartesian viscous fluxes in 1 dimension
  !!
  !! retains logic to simplify calculations on uniform or stretched cartesian grids
  !!
  !! See /ref conserve for theory
  !!
  SUBROUTINE ScalarFlux1D                                            &
       (numDim, fluxDir, gridSizes, numPoints, opInterval, gridType, & 
       gridMetric, gradScalar, fluxBuffer)
    
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)         :: numDim, fluxDir, gridType
    INTEGER(KIND=8), INTENT(IN)         :: gridSizes(numDim),numPoints
    INTEGER(KIND=8), INTENT(IN)         :: opInterval(2*numDim)
    REAL(KIND=8),    INTENT(IN), TARGET :: gridMetric(numDim*numDim*numPoints)
    REAL(KIND=8),    INTENT(IN), TARGET :: gradScalar(numPoints*numDim)
    REAL(KIND=8),    INTENT(OUT)        :: fluxBuffer(numPoints)

    INTEGER(KIND=8) :: dirOffset, metricOffset
    INTEGER(KIND=8) :: velIndex, fluxIndex, velOffset, heatOffset, dimOffset
    REAL(KIND=8)    :: gridScale
    REAL(KIND=8)    :: minusOne

    REAL(KIND=8),    DIMENSION(:), POINTER :: gradPtr
    REAL(KIND=8),    DIMENSION(:), POINTER :: metricPtr

    dirOffset = (fluxDir-1)*numPoints

    IF(gridType < RECTILINEAR) THEN

       gridScale    =  gridMetric(fluxDir)
       gradPtr      => gradScalar(dirOffset+1:dirOffset+numPoints)

       ! Grid metric is scalar and constant, this gets flux = metric * grad(Y)
       CALL YAX(numDim,numPoints,gridSizes,opInterval,gridScale,gradPtr,fluxBuffer)       

    ELSE IF (gridType == RECTILINEAR) THEN
       
       metricPtr    => gridMetric(dirOffset+1:dirOffset+numPoints)
       gradPtr      => gradScalar(dirOffset+1:dirOffset+numPoints)

       ! Grid metric is scalar, this gets flux = metric * tau_(iDim)(
       CALL ZXY(numDim,numPoints,gridSizes,opInterval,metricPtr,gradPtr,fluxBuffer)       

    ELSE ! must be curvilinear
       
       metricOffset = dirOffset*numDim
       metricPtr => gridMetric(metricOffset+1:metricOffset+numDim*numPoints)
       
       CALL ZXDOTY(numDim,numPoints,gridSizes,opInterval,numDim,metricPtr,gradScalar,fluxBuffer)
       
    ENDIF

  END SUBROUTINE ScalarFlux1D

END MODULE Viscid
