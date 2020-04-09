SUBROUTINE AXPY(N1,N2,N3,localInterval,a,X,Y,Z)

  INTEGER(KIND=8), INTENT(IN)  :: N1, N2, N3, localInterval(6)
  REAL(KIND=8),    INTENT(IN)  :: a,X(N1*N2*N3)
  REAL(KIND=8),    INTENT(IN)  :: Y(N1*N2*N3)
  REAL(KIND=8),    INTENT(OUT) :: Z(N1*N2*N3)

  INTEGER(KIND=8) :: I, J, K, LI, LIK, LIJ, IS1, IE1
  INTEGER(KIND=8) :: xyPlane
  INTEGER(KIND=8) :: S1,E1,S2,E2,S3,E3

  S1 = localInterval(1)
  E1 = localInterval(2)
  S2 = localInterval(3)
  E2 = localInterval(4)
  S3 = localInterval(5)
  E3 = localInterval(6)

  yzPlane = (N2*N3)
  
  DO I = S1, E1
     !LI = (I-1)
     DO J = S2, E2
        LIJ = (J-1)*N1 + I
        DO K = S3, E3
           LIJK = (K-1)*yzPlane + LIJ
           Z(LIJK) = a*X(LIJK) + Y(LIJK)
        END DO
     END DO
  END DO
  

END SUBROUTINE AXPY
