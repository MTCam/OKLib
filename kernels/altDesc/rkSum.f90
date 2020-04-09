SUBROUTINE RKSUM(N1,N2,N3,localInterval,h,K1,K2,K3,K4,stateData)

  INTEGER(KIND=8), INTENT(IN)    :: N1, N2, N3, localInterval(6)
  REAL(KIND=8),    INTENT(IN)    :: h
  REAL(KIND=8),    INTENT(IN)    :: K1(N1*N2*N3)
  REAL(KIND=8),    INTENT(IN)    :: K2(N1*N2*N3)
  REAL(KIND=8),    INTENT(IN)    :: K3(N1*N2*N3)
  REAL(KIND=8),    INTENT(IN)    :: K4(N1*N2*N3)
  REAL(KIND=8),    INTENT(INOUT) :: stateData(N1*N2*N3)

  INTEGER(KIND=8) :: I, J, K, LI, LIK, LIJ, IS1, IE1
  INTEGER(KIND=8) :: xyPlane
  INTEGER(KIND=8) :: S1,E1,S2,E2,S3,E3
  REAL(KIND=8)    :: fac1,fac2

  fac1 = h/6.0_8
  fac2 = h/3.0_8

  S1 = localInterval(1)
  E1 = localInterval(2)
  S2 = localInterval(3)
  E2 = localInterval(4)
  S3 = localInterval(5)
  E3 = localInterval(6)

  xyPlane = (N1*N2)
  
  DO K = S3, E3
     LIK = (K-1)*xyPlane
     DO J = S2, E2
        LIJ = (J-1)*N1 + LIK
        DO I = S1, E1
           LI = LIJ + I
           stateData(LI) = fac1*(K1(LI)+K4(LI)) + fac2*(K2(LI) + K3(LI)) + stateData(LI)
        END DO
     END DO
  END DO
  

END SUBROUTINE RKSUM
