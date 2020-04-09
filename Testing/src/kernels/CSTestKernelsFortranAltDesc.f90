!
! @file Kernels for testing
!
MODULE TESTKERNELS

  IMPLICIT NONE

CONTAINS

  ! @ICE block=simpleReplaceFortran
  SUBROUTINE TESTKERNEL(result)

    INTEGER :: result

    result = 0
    WRITE(*,*) 'Hello PlasCom2 World.'

  END SUBROUTINE TESTKERNEL
  ! @ICE endblock

END MODULE TESTKERNELS
