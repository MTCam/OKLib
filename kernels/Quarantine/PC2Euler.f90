MODULE Euler
  
  USE OPERATORS

  IMPLICIT NONE

CONTAINS

  SUBROUTINE UniformRHS( numDim, gridSizes, numPoints, opInterval, gridMetric,       &
                         numStencils, numStencilValues, stencilSizes, stencilStarts, &
                         stencilOffsets, stencilWeights, stencilID, numTransport,    &
                         numPointsApply,applyPoints,                                 &
                         rhoBuffer,rhoVBuffer,rhoEBuffer,transportBuffer,velHat,     &
                         pressureBuffer, rhsBuffer)

    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN)  :: numDim, numStencils, numStencilValues, numComponents, numTransport
    INTEGER(KIND=8), INTENT(IN)  :: gridSizes(numDim),  numPoints, numPointsApply
    INTEGER(KIND=8), INTENT(IN)  :: opInterval(2*numDim), applyPoints(numPointsApply)
    INTEGER(KIND=4), INTENT(IN)  :: stencilSizes(numStencils),stencilStarts(numStencils)
    INTEGER(KIND=4), INTENT(IN)  :: stencilOffsets(numValues)
    INTEGER(KIND=4), INTENT(IN)  :: stencilID(numPoints)
    REAL(KIND=8),    INTENT(IN)  :: stencilWeights(numValues)
    REAL(KIND=8),    INTENT(IN)  :: velHat(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)  :: rhoBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN)  :: rhoVBuffer(numDim*numPoints)
    REAL(KIND=8),    INTENT(IN)  :: rhoEBuffer(numPoints)
    REAL(KIND=8),    INTENT(IN)  :: transportBuffer(numTransport*numPoints)
    REAL(KIND=8),    INTENT(IN)  :: pressureBuffer(numPoints)
    REAL(KIND=8),    INTENT(OUT) :: rhsBuffer(numDim*numPoints)


    INTEGER         :: iDim, numComponents
    INTEGER(KIND=8) :: iPoint,iX,iY,iZ,iStart,iEnd,jStart,jEnd,kStart,kEnd
    INTEGER(KIND=8) :: xIndex,zIndex,yIndex,yzIndex,xSize,ySize,zSize
    INTEGER(KIND=8) :: iPoint2, iPoint3, pointOffset



    numComponents = 1

    ! Continuity 
    DO iDim = 1,numDim

       pointOffset = (iDim-1)*numPoints

       DO iPoint = 1, numPointsApply
          pointIndex = applyPoints(iPoint)
          vectorPointIndex = pointIndex + pointOffset
          flux(pointIndex) = rhoBuffer(pointIndex) * velHat(vectorPointIndex)
       END DO

       CALL APPLYOPERATOR(region, ng, iFirstDeriv, i, flux, dflux, .FALSE., 0)
       do ii = 1, Nc
#ifdef AXISYMMETRIC
          if(ND == 2 .and. i == 1  .and. grid%XYZ(ii,1) < grid%DRmin) dflux(ii)=0d0
#endif
          rhs(ii,1) = rhs(ii,1) - dflux(ii)
       end do
       Call Save_Patch_Deriv(region, ng, 1, i, dflux)
    end do ! i


#ifndef AXISYMMETRIC
    ! ... momentum
    do i = 1, ND
       ip1 = i+1
       do j = 1, ND
          t2Mapji = t2Map(j,i)
#ifdef ENABLE_SIMD
          !DIR$ SIMD
#endif
          do ii = 1, Nc
             flux(ii) = cv(ii,ip1) * UVWhat(ii,j) + MT1(ii,t2Mapji) * dv(ii,1)
          end do
          call APPLY_OPERATOR_box(region, ng, iFirstDeriv, j, flux, dflux, .FALSE., 0)
#ifdef ENABLE_SIMD
          !DIR$ SIMD
#endif
          do ii = 1, Nc
             rhs(ii,ip1) = rhs(ii,ip1) - dflux(ii)
          end do
          Call Save_Patch_Deriv(region, ng, i+1, j, dflux)
       end do ! j
    end do ! i
#else

    ! ... momentum
    do i = 1, ND
       ip1 = i+1
       do j = 1, ND
          t2Mapji = t2Map(j,i)
#ifdef ENABLE_SIMD
          !DIR$ SIMD
#endif
          do ii = 1, Nc
             flux(ii) = cv(ii,ip1) * UVWhat(ii,j)
          end do
          call APPLY_OPERATOR_box(region, ng, iFirstDeriv, j, flux, dflux, .FALSE., 0)
          do ii = 1, Nc
             if(ND == 2 .and. j == 1  .and. grid%XYZ(ii,1) < grid%DRmin) dflux(ii)=0d0
             rhs(ii,ip1) = rhs(ii,ip1) - dflux(ii)
          end do
          Call Save_Patch_Deriv(region, ng, i+1, j, dflux)
          do ii = 1, Nc
             flux(ii) =  MT0(ii,t2Mapji) * (dv(ii,1) - 1d0/input%GamRef)
          end do
          call APPLY_OPERATOR_box(region, ng, iFirstDeriv, j, flux, dflux, .FALSE., 0)
          do ii = 1, Nc
             dflux(ii) =  dflux(ii)*JAC0(ii)*INVJAC(ii)
          end do
          do ii = 1, Nc
             rhs(ii,ip1) = rhs(ii,ip1) - dflux(ii)
          end do
       end do ! j

    end do ! i
#endif

    ! ... energy
    do i = 1, ND
#ifdef ENABLE_SIMD
       !DIR$ SIMD
#endif
       do ii = 1, Nc
          flux(ii) = (cv(ii,NDp2) + dv(ii,1)) * UVWhat(ii,i) - XI_TAU(ii,i) * dv(ii,1)
       end do
       call APPLY_OPERATOR_box(region, ng, iFirstDeriv, i, flux, dflux, .FALSE., 0)
       do ii = 1, Nc
#ifdef AXISYMMETRIC
          if(ND == 2 .and. i == 1  .and. grid%XYZ(ii,1) < grid%DRmin) dflux(ii)=0d0
#endif
          rhs(ii,NDp2) = rhs(ii,NDp2) - dflux(ii)
       end do
       Call Save_Patch_Deriv(region, ng, NDp2, i, dflux)
    end do ! i

    ! ... scalar advection
    do k = 1, nAuxVars
       do i = 1, ND
#ifdef ENABLE_SIMD
          !DIR$ SIMD
#endif
          do ii = 1, Nc
             flux(ii) = auxVars(ii,k) * UVWhat(ii,i)
          end do
          call APPLY_OPERATOR_box(region, ng, iFirstDeriv, i, flux, dflux, .FALSE., 0)
          do ii = 1, Nc
#ifdef AXISYMMETRIC
             if(ND == 2 .and. i == 1  .and. grid%XYZ(ii,1) < grid%DRmin) dflux(ii)=0d0
#endif
             rhs_auxVars(ii,k) = rhs_auxVars(ii,k) - dflux(ii)
          end do
       end do ! i

    end do ! k

#ifdef AXISYMMETRIC
    Call NS_RHS_Axisymmetric_Hopital(ng, Nc, nAuxVars, ND, iFirstDeriv, region, grid, input, flux, dflux, MT0, INVJAC_TAU, JAC0, INVJAC, rhs, rhs_auxVars, cv, dv, auxVars, 0)
#endif

    ! ... deallocate
    deallocate( UVWhat )

    ! CAA Benchmark 2, Category 1, Problem 1
    if (.false.) then
#ifdef ENABLE_SIMD
       !DIR$ SIMD
#endif
       do i = 1, Nc
          source = 0.1_dp * exp(-(log(2.0_dp))*( ((XYZ(i,1)-4.0_dp)**2+XYZ(i,2)**2)/(0.4_dp)**2 )) &
               * dsin(4.0_dp * TWOPI * time(i)) * INVJAC(i)
          rhs(i,4) = rhs(i,4) + source / (gv(i,1)-1.0_dp)
       end do
    end if

    ! ... subtract off (violation of) metric identity

    if (useMetricIdentities == TRUE) then

       if( ND == 2 ) then

#ifdef ENABLE_SIMD
          !DIR$ SIMD
#endif
          do ii = 1, Nc
             rhs(ii,1)  = rhs(ii,1)    + cv(ii,2) * ident(ii,1) &
                  + cv(ii,3) * ident(ii,2)
             rhs(ii,2)  = rhs(ii,2)    + dv(ii,1) * ident(ii,1)
             rhs(ii,3)  = rhs(ii,3)    + dv(ii,1) * ident(ii,2)
             rhs(ii,4)  = rhs(ii,4)    + (cv(ii,4) + dv(ii,1)) * cv(ii,2)/cv(ii,1) * ident(ii,1) &
                  + (cv(ii,4) + dv(ii,1)) * cv(ii,3)/cv(ii,1) * ident(ii,2) 
          end do

       else if( ND == 3 ) then

#ifdef ENABLE_SIMD
          !DIR$ SIMD
#endif
          do ii = 1, Nc
             rhs(ii,1)  = rhs(ii,1)    + cv(ii,2) * ident(ii,1) &
                  + cv(ii,3) * ident(ii,2) &
                  + cv(ii,4) * ident(ii,3)
             rhs(ii,2)  = rhs(ii,2)    + dv(ii,1) * ident(ii,1)
             rhs(ii,3)  = rhs(ii,3)    + dv(ii,1) * ident(ii,2)
             rhs(ii,4)  = rhs(ii,4)    + dv(ii,1) * ident(ii,3)
             rhs(ii,5)  = rhs(ii,5)    + (cv(ii,5) + dv(ii,1)) * cv(ii,2)/cv(ii,1) * ident(ii,1) &
                  + (cv(ii,5) + dv(ii,1)) * cv(ii,3)/cv(ii,1) * ident(ii,2) &
                  + (cv(ii,5) + dv(ii,1)) * cv(ii,4)/cv(ii,1) * ident(ii,3)
          end do

       end if ! ND

    end if ! useMetricIdentities

    ! ... gravity

    if (Froude > 0.0_DP) then
       norm(1) = sin(gravity_angle_phi) * cos(gravity_angle_theta) * invFroude
       norm(2) = sin(gravity_angle_phi) * sin(gravity_angle_theta) * invFroude

       if( ND == 2 ) then

#ifdef ENABLE_SIMD
          !DIR$ SIMD
#endif
          do i = 1, Nc
             rhs(i,2) = rhs(i,2) + cv(i,1) * norm(1) * INVJAC(i)
             rhs(i,3) = rhs(i,3) + cv(i,1) * norm(2) * INVJAC(i)
          end do

       else if( ND == 3 ) then 

          norm(3) = cos(gravity_angle_phi) * invFroude
#ifdef ENABLE_SIMD
          !DIR$ SIMD
#endif
          do i = 1, Nc
             rhs(i,2) = rhs(i,2) + cv(i,1) * norm(1) * INVJAC(i)
             rhs(i,3) = rhs(i,3) + cv(i,1) * norm(2) * INVJAC(i)
             rhs(i,4) = rhs(i,4) + cv(i,1) * norm(3) * INVJAC(i)
          end do

       end if ! ND

    end if ! Froude

    return


  END SUBROUTINE EulerRHS

END MODULE Euler
