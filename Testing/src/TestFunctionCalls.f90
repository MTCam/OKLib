MODULE TestModule

  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER :: TEST0   = 0
  INTEGER(KIND=4), PARAMETER :: TEST1   = 1

  INTERFACE

     SUBROUTINE CFunc(a)
       IMPLICIT NONE
       INTEGER(8) :: a
     END SUBROUTINE CFunc

     SUBROUTINE CFunc2(a)
       IMPLICIT NONE
       REAL(8) :: a
     END SUBROUTINE CFunc2

  END INTERFACE

CONTAINS

  SUBROUTINE Function1(cppAddress)
    
    INTEGER(8), INTENT(IN) :: cppAddress
    
    CALL CFunc(cppAddress)
    
  END SUBROUTINE Function1

  SUBROUTINE Function2(a)
    
    REAL(8), INTENT(INOUT) :: a(2)
    
    REAL(8) :: b(2)

    a(1) = 1.0_8
    a(2) = 2.0_8
    b(1) = 1.0_8
    b(2) = 2.0_8

    CALL CFunc2(b(1))
    
  END SUBROUTINE Function2

END MODULE TestModule
