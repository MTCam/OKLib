!
! @file Kernels for testing
! 
SUBROUTINE ICEUNROLLFORTRANTESTKERNEL(result)

  INTEGER :: result
  INTEGER :: iBegin, iEnd, jBegin, jEnd, kBegin, kEnd

  iBegin = 0
  iEnd = 20
  jBegin = 0
  jEnd = 10 
  kBegin = 0
  kEnd = 5 

  result = 0
  
  ! @ICE loop=unrollTestFortran
  do i = iBegin, iEnd
    do j = jBegin, jEnd
      do k = kBegin, kEnd
        result = result + 1
      enddo
    enddo
  enddo
  ! @ICE endloop

END SUBROUTINE ICEUNROLLFORTRANTESTKERNEL
! 
! 
SUBROUTINE ICEINTERCHANGEFORTRANTESTKERNEL(result)

  INTEGER :: result
  INTEGER :: iBegin, iEnd, jBegin, jEnd, kBegin, kEnd

  iBegin = 0
  iEnd = 20
  jBegin = 0
  jEnd = 10 
  kBegin = 0
  kEnd = 5 

  result = 0
  
  ! @ICE loop=interchangeTestFortran
  do i = iBegin, iEnd
    do j = jBegin, jEnd
      do k = kBegin, kEnd
        result = result + 1
      enddo
    enddo
  enddo
  ! @ICE endloop

END SUBROUTINE ICEINTERCHANGEFORTRANTESTKERNEL
! 
! 
SUBROUTINE ICETILEFORTRANTESTKERNEL(result)

  INTEGER :: result
  INTEGER :: iBegin, iEnd, jBegin, jEnd, kBegin, kEnd

  iBegin = 0
  iEnd = 20
  jBegin = 0
  jEnd = 10 
  kBegin = 0
  kEnd = 5 

  result = 0
  
  ! @ICE loop=tileTestFortran
  do i = iBegin, iEnd
    do j = jBegin, jEnd
      do k = kBegin, kEnd
        result = result + 1
      enddo
    enddo
  enddo
  ! @ICE endloop

END SUBROUTINE ICETILEFORTRANTESTKERNEL
! 
! 
SUBROUTINE ICESTRIPMINEFORTRANTESTKERNEL(result)

  INTEGER :: result
  INTEGER :: iBegin, iEnd, jBegin, jEnd, kBegin, kEnd

  iBegin = 0
  iEnd = 20
  jBegin = 0
  jEnd = 10 
  kBegin = 0
  kEnd = 5 

  result = 0
  
  ! @ICE loop=stripMineTestFortran
  do i = iBegin, iEnd
    do j = jBegin, jEnd
      do k = kBegin, kEnd
        result = result + 1
      enddo
    enddo
  enddo
  ! @ICE endloop

END SUBROUTINE ICESTRIPMINEFORTRANTESTKERNEL
