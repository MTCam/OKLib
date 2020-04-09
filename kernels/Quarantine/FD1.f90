SUBROUTINE FD1LOAD
  WRITE(*,*) 'FD1 Kernel library loaded'

END SUBROUTINE FD1LOAD

SUBROUTINE FD1INTERIOR(N1,N2,N3,S1,E1,S2,E2,S3,E3,DX, U, KN)

  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: DX
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)
  
  INTEGER(KIND=8) :: I, J, K, LI, LIK, LIJ, IS1, IE1
  REAL(KIND=8)    :: fac
  INTEGER(KIND=8) :: plane
  
  IS1 = S1
  IE1 = E1
  IF(S1 == 1)  IS1 = 2
  IF(E1 == N1) IE1 = N1 - 1

  fac = 0.5_8 / DX
  plane = (N1*N2)

  DO K = S3, E3
     LIK = (K-1)*plane
    DO J = S2, E2
       LIJ = (J-1)*N1 + LIK
      DO I = IS1, IE1
        LI = LIJ + I
        KN(LI) = ( U(LI+1) - U(LI-1) ) * FAC
      END DO
    END DO
  END DO


END SUBROUTINE FD1INTERIOR

SUBROUTINE FD1INTERIORX(N1,N2,N3,S1,E1,S2,E2,S3,E3,DX, U, KN)

  IMPLICIT NONE
  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: DX
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)
  
  INTEGER(KIND=8) :: I, J, K, LI, LIK, LIJ, IS1, IE1
  REAL(KIND=8)    :: fac
  INTEGER(KIND=8) :: plane
  
  IS1 = S1
  IE1 = E1
  IF(S1 == 1)  IS1 = 2
  IF(E1 == N1) IE1 = N1 - 1

  fac = 0.5_8 / DX
  plane = (N1*N2)

  DO K = S3, E3
     LIK = (K-1)*plane
    DO J = S2, E2
       LIJ = (J-1)*N1 + LIK
      DO I = IS1, IE1
        LI = LIJ + I
        KN(LI) = ( U(LI+1) - U(LI-1) ) * FAC
      END DO
    END DO
  END DO


END SUBROUTINE FD1INTERIORX

SUBROUTINE FD1INTERIORY(N1,N2,N3,S1,E1,S2,E2,S3,E3,DY,U,KN)

  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: DY
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)
  
  INTEGER(KIND=8) :: I, J, K, LI, LIK, LIJ
  INTEGER(KIND=8) :: IS1, IS2, IS3, IE1, IE2, IE3
  REAL(KIND=8)    :: fac
  INTEGER(KIND=8) :: plane
  
  IS1 = S1
  IE1 = E1
  IS2 = S2
  IE2 = E2
  IS3 = S3
  IE3 = E3

  IF(S2 == 1)  IS2 = 2
  IF(E2 == N2) IE2 = N2 - 1

  fac = 0.5_8 / DY
  plane = (N1*N2)

  DO K = IS3, IE3
     LIK = (K-1)*plane
    DO J = IS2, IE2
       LIJ = (J-1)*N1 + LIK
      DO I = IS1, IE1
        LI = LIJ + I
        KN(LI) = ( U(LI+N1) - U(LI-N1) ) * FAC
      END DO
    END DO
  END DO


END SUBROUTINE FD1INTERIORY

SUBROUTINE FD1INTERIORZ(N1,N2,N3,S1,E1,S2,E2,S3,E3,DZ,U,KN)

  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: DZ
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)
  
  INTEGER(KIND=8) :: I, J, K, LI, LIK, LIJ
  INTEGER(KIND=8) :: IS1, IS2, IS3, IE1, IE2, IE3
  REAL(KIND=8)    :: fac
  INTEGER(KIND=8) :: plane
  
  IS1 = S1
  IE1 = E1
  IS2 = S2
  IE2 = E2
  IS3 = S3
  IE3 = E3

  IF(S3 == 1)  IS3 = 2
  IF(E3 == N3) IE3 = N3 - 1

  fac = 0.5_8 / DZ
  plane = (N1*N2)

  DO K = IS3, IE3
     LIK = (K-1)*plane
    DO J = IS2, IE2
       LIJ = (J-1)*N1 + LIK
      DO I = IS1, IE1
        LI = LIJ + I
        KN(LI) = ( U(LI+plane) - U(LI-plane) ) * FAC
      END DO
    END DO
  END DO


END SUBROUTINE FD1INTERIORZ

SUBROUTINE FD1INTERIOR2(N1,N2,N3,S1,E1,S2,E2,S3,E3,DX,DY,U,KX,KY)

  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: DX,DY
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: KX(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: KY(N1*N2*N3)
  
  INTEGER(KIND=8) :: I, J, K
  INTEGER(KIND=8) :: LI, LIK, LIJ, IS1, IE1
  INTEGER(KIND=8) :: LII, IS2, IE2
  REAL(KIND=8)    :: fac1,fac2
  INTEGER(KIND=8) :: plane
  
  IS1 = S1
  IE1 = E1
  IS2 = S2
  IE2 = E2

  IF(S1 == 1)  IS1 = 2
  IF(E1 == N1) IE1 = N1 - 1
  IF(S2 == 1)  IS2 = 2
  IF(E2 == N2) IE2 = N2 - 1
  


  fac1 = 0.5_8 / DX
  fac2 = 0.5_8 / DY

  plane = (N1*N2)

  DO K = S3, E3
     LIK = (K-1)*plane

    DO J = S2, E2
       LIJ = (J-1)*N1 + LIK
      DO I = IS1, IE1
        LI = LIJ + I
        KX(LI) = ( U(LI+1) - U(LI-1) ) * FAC1
      END DO
    END DO

    DO I = S1, E1
       LII = LIK + I
       DO J = IS2, IE2
          LIJ = (J-1)*N1 + LII
          KY(LIJ) = ( U(LIJ+N1) - U(LIJ-N1) ) * FAC2
       END DO
    END DO

  END DO

END SUBROUTINE FD1INTERIOR2

SUBROUTINE FD1INTERIOR3(N1,N2,N3,S1,E1,S2,E2,S3,E3,DX,DY,DZ,U,KX,KY,KZ)

  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: DX,DY,DZ
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: KX(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: KY(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: KZ(N1*N2*N3)
  
  INTEGER(KIND=8) :: I, J, K
  INTEGER(KIND=8) :: LI, LIK, LIJ, IS1, IE1
  INTEGER(KIND=8) :: LII, IS2, IE2, IS3, IE3
  REAL(KIND=8)    :: fac1,fac2,fac3
  INTEGER(KIND=8) :: plane
  
  IS1 = S1
  IE1 = E1

  IS2 = S2
  IE2 = E2

  IS3 = S3
  IE3 = E3

  ! Operate on the interior only
  IF(S1 == 1)  IS1 = 2
  IF(E1 == N1) IE1 = N1 - 1

  IF(S2 == 1)  IS2 = 2
  IF(E2 == N2) IE2 = N2 - 1

  IF(S3 == 1)  IS3 = 2
  IF(E3 == N3) IE3 = N3 - 1


  fac1 = 0.5_8 / DX
  fac2 = 0.5_8 / DY
  fac3 = 0.5_8 / DZ

  plane = (N1*N2)

  DO K = S3, E3

     LIK = (K-1)*plane

    DO J = S2, E2
       LIJ = (J-1)*N1 + LIK
      DO I = IS1, IE1
        LI = LIJ + I
        KX(LI) = ( U(LI+1) - U(LI-1) ) * FAC1
      END DO
    END DO

    DO I = S1, E1
       LII = LIK + I
       DO J = IS2, IE2
          LIJ = (J-1)*N1 + LII
          KY(LIJ) = ( U(LIJ+N1) - U(LIJ-N1) ) * FAC2
       END DO
    END DO

  END DO
  
  DO J = S2, E2
     LIJ = (J-1)*N1 
     DO I = S1, E1
        LI = LIJ + I
        DO K = IS3, IE3
           LIK = (K-1)*plane + LI
           KZ(LIK) = ( U(LIK+plane) - U(LIK-plane) ) * FAC3
        END DO
     END DO
  END DO

END SUBROUTINE FD1INTERIOR3

SUBROUTINE FD1GRAD3(N1,N2,N3,S1,E1,S2,E2,S3,E3,DX,DY,DZ,U,GRADU)

  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: DX,DY,DZ
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: GRADU(3*N1*N2*N3)
!  REAL(KIND=8),    INTENT(OUT) :: KY(N1*N2*N3)
!  REAL(KIND=8),    INTENT(OUT) :: KZ(N1*N2*N3)
  
  INTEGER(KIND=8) :: I, J, K
  INTEGER(KIND=8) :: LI, LIK, LIJ, IS1, IE1
  INTEGER(KIND=8) :: LII, IS2, IE2, IS3, IE3
  REAL(KIND=8)    :: fac1,fac2,fac3
  INTEGER(KIND=8) :: plane
  INTEGER(KIND=8) :: numPoints,numPoints2

  numPoints = N1*N2*N3
  numPoints2 = 2*N1*N2*N3

  IS1 = S1
  IE1 = E1

  IS2 = S2
  IE2 = E2

  IS3 = S3
  IE3 = E3

  ! Operate on the interior only
  IF(S1 == 1)  IS1 = 2
  IF(E1 == N1) IE1 = N1 - 1

  IF(S2 == 1)  IS2 = 2
  IF(E2 == N2) IE2 = N2 - 1

  IF(S3 == 1)  IS3 = 2
  IF(E3 == N3) IE3 = N3 - 1


  fac1 = 0.5_8 / DX
  fac2 = 0.5_8 / DY
  fac3 = 0.5_8 / DZ

  plane = (N1*N2)

  DO K = S3, E3

     LIK = (K-1)*plane

    DO J = S2, E2
       LIJ = (J-1)*N1 + LIK
      DO I = IS1, IE1
        LI = LIJ + I
        GRADU(LI) = ( U(LI+1) - U(LI-1) ) * FAC1
      END DO
    END DO

    DO I = S1, E1
       LII = LIK + I
       DO J = IS2, IE2
          LIJ = (J-1)*N1 + LII
          GRADU(LIJ+numPoints) = ( U(LIJ+N1) - U(LIJ-N1) ) * FAC2
       END DO
    END DO

  END DO
  
  DO J = S2, E2
     LIJ = (J-1)*N1 
     DO I = S1, E1
        LI = LIJ + I
        DO K = IS3, IE3
           LIK = (K-1)*plane + LI
           GRADU(LIK+numPoints2) = ( U(LIK+plane) - U(LIK-plane) ) * FAC3
        END DO
     END DO
  END DO

END SUBROUTINE FD1GRAD3

SUBROUTINE FD1DIV3(N1,N2,N3,S1,E1,S2,E2,S3,E3,DX,DY,DZ,UVEC,DIVU)

  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: DX,DY,DZ
  REAL(KIND=8),    INTENT(IN)  :: UVEC(3*N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: DIVU(N1*N2*N3)
!  REAL(KIND=8),    INTENT(OUT) :: KY(N1*N2*N3)
!  REAL(KIND=8),    INTENT(OUT) :: KZ(N1*N2*N3)
  
  INTEGER(KIND=8) :: I, J, K
  INTEGER(KIND=8) :: LI, LIK, LIJ, IS1, IE1
  INTEGER(KIND=8) :: LII, IS2, IE2, IS3, IE3
  REAL(KIND=8)    :: fac1,fac2,fac3
  INTEGER(KIND=8) :: plane
  INTEGER(KIND=8) :: numPoints,numPoints2

  numPoints = N1*N2*N3
  numPoints2 = 2*N1*N2*N3

  IS1 = S1
  IE1 = E1

  IS2 = S2
  IE2 = E2

  IS3 = S3
  IE3 = E3

  ! Operate on the interior only
  IF(S1 == 1)  IS1 = 2
  IF(E1 == N1) IE1 = N1 - 1

  IF(S2 == 1)  IS2 = 2
  IF(E2 == N2) IE2 = N2 - 1

  IF(S3 == 1)  IS3 = 2
  IF(E3 == N3) IE3 = N3 - 1


  fac1 = 0.5_8 / DX
  fac2 = 0.5_8 / DY
  fac3 = 0.5_8 / DZ

  plane = (N1*N2)

  DO K = S3, E3

     LIK = (K-1)*plane

    DO J = S2, E2
       LIJ = (J-1)*N1 + LIK
      DO I = IS1, IE1
        LI = LIJ + I
        DIVU(LI) = ( UVEC(LI+1) - UVEC(LI-1) ) * FAC1
      END DO
    END DO

    DO I = S1, E1
       LII = LIK + I
       DO J = IS2, IE2
          LIJ = (J-1)*N1 + LII
          DIVU(LIJ) = DIVU(LIJ) + ( UVEC((LIJ+N1)+numPoints) - UVEC((LIJ-N1)+numPoints) ) * FAC2
       END DO
    END DO

  END DO
  
  DO J = S2, E2
     LIJ = (J-1)*N1 
     DO I = S1, E1
        LI = LIJ + I
        DO K = IS3, IE3
           LIK = (K-1)*plane + LI
           DIVU(LIK) = DIVU(LIK) + ( UVEC((LIK+plane)+numPoints2) - UVEC((LIK-plane)+numPoints2) ) * FAC3
        END DO
     END DO
  END DO

END SUBROUTINE FD1DIV3

! @ICE block=FD1DIV3ALL
SUBROUTINE FD1DIV3ALL(N1,N2,N3,OPINTERVAL,DX,DY,DZ,&
     UVHLEFT,UVHRIGHT,UVHFRONT,UVHBACK,UVHBOTTOM,UVHTOP, &
     UVEC,DIVU)
  
  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3
  INTEGER(KIND=8), INTENT(IN)  :: OPINTERVAL(6) ! S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: DX,DY,DZ
  REAL(KIND=8),    INTENT(IN)  :: UVHLEFT(3*N2*N3)
  REAL(KIND=8),    INTENT(IN)  :: UVHRIGHT(3*N2*N3)
  REAL(KIND=8),    INTENT(IN)  :: UVHFRONT(3*N1*N3)
  REAL(KIND=8),    INTENT(IN)  :: UVHBACK(3*N1*N3)
  REAL(KIND=8),    INTENT(IN)  :: UVHBOTTOM(3*N1*N2)
  REAL(KIND=8),    INTENT(IN)  :: UVHTOP(3*N1*N2)
  REAL(KIND=8),    INTENT(IN)  :: UVEC(3*N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: DIVU(N1*N2*N3)
  !  REAL(KIND=8),    INTENT(OUT) :: KY(N1*N2*N3)
  !  REAL(KIND=8),    INTENT(OUT) :: KZ(N1*N2*N3)
  
  INTEGER(KIND=8) :: I, J, K
  INTEGER(KIND=8) :: IH,IU,JH,JU,KH,KU
  INTEGER(KIND=8) :: LI, LIK, LIJ, IS1, IE1
  INTEGER(KIND=8) :: LII, IS2, IE2, IS3, IE3
  INTEGER(KIND=8) :: S1,E1,S2,E2,S3,E3
  REAL(KIND=8)    :: fac1,fac2,fac3
  INTEGER(KIND=8) :: plane,plane2,numHalo
  INTEGER(KIND=8) :: numPoints,numPoints2
  
  DIVU(:) = 0.0_8
  
  numPoints = N1*N2*N3
  numPoints2 = 2*N1*N2*N3
  
  fac1 = 0.5_8 / DX
  fac2 = 0.5_8 / DY
  fac3 = 0.5_8 / DZ
  
  plane  = (N1*N2)
  plane2 = 2*N1*N2
  
  S1 = OPINTERVAL(1)
  E1 = OPINTERVAL(2)
  
  S2 = OPINTERVAL(3)
  E2 = OPINTERVAL(4)
  
  S3 = OPINTERVAL(5)
  E3 = OPINTERVAL(6)
  
  IS1 = S1
  IE1 = E1
  
  IS2 = S2
  IE2 = E2
  
  IS3 = S3
  IE3 = E3
!  WRITE(*,*) 'DIV INTERVAL: ',IS1,IE1,IS2,IE2,IS3,IE3
  ! ============== X DIRECTION =================
  
  ! PROCESS LEFT BOUNDARY
  IF(IS1 == 1) THEN ! select only threads with intervals which include the left boundary
     
     DO K = S3, E3
        KH = (K-1)*N2
        KU = KH*N1 + 1
        DO J = S2, E2
           JH = KH + J
           JU = KU + (J-1)*N1
           DIVU(JU) = ( UVEC(JU+1) - UVHLEFT(JH) ) * FAC1
!           WRITE(*,*) 'LEFT BOUNDARY: DIV(',JU,') = ',DIVU(JU)
        END DO
     END DO
     
  ENDIF
  
  ! PROCESS RIGHT BOUNDARY
  IF(IE1 == N1) THEN ! select only threads touching right boundary
     
     DO K = S3, E3
        KH = (K-1)*N2
        KU = KH*N1 + N1
        DO J = S2, E2
           JH = KH + J
           JU = KU + (J-1)*N1
           DIVU(JU) = ( UVHRIGHT(JH) - UVEC(JU-1) ) * FAC1
!           WRITE(*,*) 'RIGHT BOUNDARY: DIV(',JU,') = ',DIVU(JU)
        END DO
     END DO
     
  ENDIF
  
  ! PROCESS INTERIOR
  IF((E1 > S1) .OR. ((E1 < N1) .AND. (S1 > 1))) THEN ! Necessary thread selection!!! (awful bug resolved)
     
     IF(S1 == 1)  IS1 = 2
     IF(E1 == N1) IE1 = N1 - 1
     
     DO K = S3, E3
        KU = (K-1)*plane
        DO J = S2, E2
           JU = (J-1)*N1 + KU
           DO I = IS1, IE1
              IU = JU + I
              DIVU(IU) = ( UVEC(IU+1) - UVEC(IU-1) ) * FAC1
           END DO
        END DO
     END DO
  ENDIF
  
  IS1 = S1
  IE1 = E1
  
!  DO K = S3, E3
!     KU = (K-1)*plane
!     DO J = S2, E2
!        JU = (J-1)*N1 + KU
!        DO I = IS1, IE1
!           IU = JU + I
!           WRITE(*,*) 'XDIVU(',IU,') = ', DIVU(IU)
!        END DO
!     END DO
!  END DO
  
  ! ========== Y DIRECTION ==============
  
  numHalo = N3*N1
  ! PROCESS FRONT BOUNDARY
  IF(S2 == 1) THEN
     DO K = S3, E3
        KH = (K-1)*N1
        KU = KH*N2
        DO I = S1, E1
           IH = KH + I
           IU = KU + I
           DIVU(IU) = DIVU(IU) + ( UVEC(IU+N1+numPoints) - UVHFRONT(IH+numHalo) ) * FAC2
        END DO
     END DO
  ENDIF
!  DO K = S3, E3
!     KU = (K-1)*plane
!     DO J = S2, E2
!        JU = (J-1)*N1 + KU
!        DO I = IS1, IE1
!           IU = JU + I
!           WRITE(*,*) 'after front YDIVU(',IU,') = ', DIVU(IU)
!        END DO
!     END DO
!  END DO
  
  ! PROCESS BACK BOUNDARY
  IF(E2 == N2) THEN
     J = (N2-1)*N1
     DO K = S3, E3
        KH = (K-1)*N1
        KU = KH*N2 + J 
        DO I = S1, E1
           IH = KH + I
           IU = KU + I
           DIVU(IU) = DIVU(IU) + ( UVHBACK(IH+numHalo) - UVEC((IU-N1)+numPoints) ) * FAC2
!           WRITE(*,*) 'DIVU(',IU,') = uvhyBACK(',IH,',',UVHBACK(IH+numHalo),') - UVECY(',IU-N1,',',UVEC((IU-N1)+numPoints),')'
        END DO
     END DO
  ENDIF

!  DO K = S3, E3
!     KU = (K-1)*plane
!     DO J = S2, E2
!        JU = (J-1)*N1 + KU
!        DO I = IS1, IE1
!           IU = JU + I
!           WRITE(*,*) 'after back YDIVU(',IU,') = ', DIVU(IU)
!        END DO
!     END DO
!  END DO
  
  ! PROCESS INTERIOR in Y
  IF((E2 > S2) .OR. ((E2 < N2) .AND. (S2 > 1))) THEN 
     
     IF(S2 == 1)  IS2 = 2
     IF(E2 == N2) IE2 = N2 - 1
     
     DO K = S3, E3
        KU = (K-1)*plane
        DO J = IS2, IE2
           JU = (J-1)*N1 + KU
           DO I = S1, E1
              IU = JU + I
              DIVU(IU) = DIVU(IU) + ( UVEC(IU+N1+numPoints) - UVEC(IU-N1+numPoints) ) * FAC2
           END DO
        END DO
     END DO
     
     IS2 = S2
     IE2 = E2
     
  ENDIF
  
  ! ========== Z DIRECTION ==============
  
  ! PROCESS BOTTOM BOUNDARY
  IF(S3 == 1) THEN
     K = 1
     DO J = S2, E2
        JU = (J-1)*N1
        DO I = S1, E1
           IU = JU + I
           DIVU(IU) = DIVU(IU) + ( UVEC(IU+plane+numPoints2) - UVHBOTTOM(IU+plane2) ) * FAC3
!           WRITE(*,*) 'BOTTOM BOUNDARY: DIV(',IU,') = ',DIVU(IU)
        END DO
     END DO
  ENDIF
  
  ! PROCESS TOP BOUNDARY
  IF(E3 == N3) THEN
     K = (N3-1)*plane
     DO J = S2, E2
        JH = (J-1)*N1
        JU = JH + K
        DO I = S1, E1
           IH = JH + I
           IU = JU + I
           DIVU(IU) = DIVU(IU) + ( UVHTOP(IH+plane2) - UVEC(IU-plane+numPoints2) ) * FAC3
 !          WRITE(*,*) 'TOP BOUNDARY: DIV(',IU,') = ',DIVU(IU)
        END DO
     END DO
  ENDIF
  
  ! PROCESS INTERIOR in Z
  IF((E3 > S3) .OR. ((E3 < N3) .AND. (S3 > 1))) THEN 
     
     IF(S3 == 1)  IS3 = 2
     IF(E3 == N3) IE3 = N3 - 1
     
     DO K = IS3, IE3
        KU = (K-1)*plane
        DO J = S2, E2
           JU = (J-1)*N1 + KU
           DO I = S1, E1
              IU = JU + I
              DIVU(IU) = DIVU(IU) + ( UVEC(IU+plane+numPoints2) - UVEC(IU-plane+numPoints2) ) * FAC3
           END DO
        END DO
     END DO
     
  ENDIF
  
END SUBROUTINE FD1DIV3ALL



! Use non-zero tau to specify boundary data
SUBROUTINE FD1LEFTBOUNDARY(N1,N2,N3,S1,E1,S2,E2,S3,E3, DX, u0, Tau, U, KN)

  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3), DX, Tau, u0
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)

  INTEGER(KIND=8) :: I, J, K, LI, plane, KI
  REAL(KIND=8)    :: fac,bcfac,bcterm

  fac = 1.0_8 / DX
  plane = (N1*N2)
  bcfac = Tau*fac
  bcterm = bcfac*u0

  DO K = S3, E3
    KI = (K-1)*plane + 1
    DO J = S2, E2
       LI = KI + (J-1)*N1
       KN(LI) = ( U(LI+1) - U(LI) ) * FAC + bcfac*U(LI) - bcterm
    END DO
 END DO

END SUBROUTINE FD1LEFTBOUNDARY


! Use non-zero tau to specify boundary data
SUBROUTINE FD1RIGHTBOUNDARY(N1,N2,N3,S1,E1,S2,E2,S3,E3,DX, u0, Tau, U, KN)

  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3), DX, Tau, u0
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)
  INTEGER(KIND=8)              :: J, K, LI, plane, KI
  REAL(KIND=8)                 :: fac,bcfac,bcterm


  fac = 1.0_8 / DX
  plane = (N1*N2)
  bcfac = Tau*fac
  bcterm = bcfac*u0

  DO K = S3, E3
     KI = (K-1)*plane
    DO J = S2, E2
        LI = KI + J*N1
        KN(LI) = ( U(LI) - U(LI-1) ) * FAC + bcfac*U(LI) - bcterm
    END DO
  END DO
  
END SUBROUTINE FD1RIGHTBOUNDARY


SUBROUTINE FD1HALOLEFT(N1,N2,N3,S1,E1,S2,E2,S3,E3, DX, U, KN, HL)


  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3), HL(N2*N3), DX
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)

  INTEGER(KIND=8) :: J, K, LI, plane, LII, KI, KII
  REAL(KIND=8)    :: FAC

  FAC = 0.5_8/DX
  plane = (N1*N2)

  DO K = S3, E3
     KII = (K-1)*N2
     KI  = (K-1)*plane + 1
     DO J = S2, E2
        LII = KII + J
        LI  = KI  + (J-1)*N1
        KN(LI) = ( U(LI+1) - HL(LII) ) * FAC
     END DO
  END DO
  
END SUBROUTINE FD1HALOLEFT


SUBROUTINE FD1HALORIGHT(N1, N2, N3, S1, E1, S2, E2, S3, E3, DX, U, KN, HR)
  
  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8), INTENT(IN) :: U(N1*N2*N3), HR(N2*N3), DX
  REAL(KIND=8), INTENT(OUT) :: KN(N1*N2*N3)
  
  INTEGER(KIND=8) :: J, K, LI, plane, LII, KI, KII
  REAL(KIND=8)    :: fac
  

  fac = 0.5_8/DX
  plane = (N1*N2)
  
  DO K = S3, E3
     KII = (K-1)*N2
     KI  = (K-1)*plane
     DO J = S2, E2
        LII = KII + J
        LI = KI + J*N1
        KN(LI) = ( HR(LII) - U(LI-1) ) * FAC
     END DO
  END DO
  
END SUBROUTINE FD1HALORIGHT

SUBROUTINE FD1HALOFRONT(N1,N2,N3,S1,E1,S2,E2,S3,E3, DY, U, KN, HL)


  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3), HL(N1*N3), DY
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)

  INTEGER(KIND=8) :: I, J, K, LI, plane, LII, KI, KII
  REAL(KIND=8)    :: FAC

  FAC = 0.5_8/DY
  plane = N1*N2

  DO K = S3, E3
     KII = (K-1)*N1
     KI  = (K-1)*plane
     DO I = S1, E1
        LII = KII + I
        LI  = KI  + I
        KN(LI) = ( U(LI+N1) - HL(LII) ) * FAC
     END DO
  END DO
  
END SUBROUTINE FD1HALOFRONT

SUBROUTINE FD1HALOBACK(N1,N2,N3,S1,E1,S2,E2,S3,E3, DY, U, KN, HL)


  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3), HL(N1*N3), DY
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)

  INTEGER(KIND=8) :: I, J, K, LI, plane, LII, KI, KII
  REAL(KIND=8)    :: FAC

  FAC = 0.5_8/DY
  plane = N1*N2
  J = (N2 - 1)*N1

  DO K = S3, E3
     KII = (K-1)*N1 
     KI  = (K-1)*plane + J 
     DO I = S1, E1
        LII = KII + I
        LI  = KI  + I
        KN(LI) = ( HL(LII) - U(LI-N1)) * FAC
     END DO
  END DO
  
END SUBROUTINE FD1HALOBACK

SUBROUTINE FD1HALOTOP(N1,N2,N3,S1,E1,S2,E2,S3,E3, DZ, U, KN, HL)


  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3), HL(N1*N2), DZ
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)

  INTEGER(KIND=8) :: I, J, K, LI, plane, LII, KI, KII
  REAL(KIND=8)    :: FAC

  FAC = 0.5_8/DZ
  plane = N1*N2
  K = (N3-1)*plane
  
  DO J = S2, E2
     KII = (J-1)*N1 
     KI  = K + KII
     DO I = S1, E1
        LII = KII + I
        LI  = KI  + I
        KN(LI) = ( HL(LII) - U(LI-plane) ) * FAC
     END DO
  END DO
  
END SUBROUTINE FD1HALOTOP

SUBROUTINE FD1HALOBOTTOM(N1,N2,N3,S1,E1,S2,E2,S3,E3, DZ, U, KN, HL)


  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, S1, E1, S2, E2, S3, E3
  REAL(KIND=8),    INTENT(IN)  :: U(N1*N2*N3), HL(N1*N2), DZ
  REAL(KIND=8),    INTENT(OUT) :: KN(N1*N2*N3)

  INTEGER(KIND=8) :: I, J, K, LI, plane, LII, KI, KII
  REAL(KIND=8)    :: FAC

  FAC = 0.5_8/DZ
  plane = N1*N2

  
  DO J = S2, E2
     KII = (J-1)*N1 
     DO I = S1, E1
        LII = KII + I
        KN(LII) = ( U(LI+plane) - HL(LII) ) * FAC
     END DO
  END DO
  
END SUBROUTINE FD1HALOBOTTOM

