MODULE vars
  IMPLICIT NONE

  PUBLIC

! Public parameters
  INTEGER,PARAMETER :: q1=SELECTED_REAL_KIND(6,30)         ! single precision
  INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)       ! double precision
  REAL(q2),PARAMETER :: pi=3.141592653589793238462643_q2    

! Public, allocatable variables
  INTEGER,ALLOCATABLE,DIMENSION(:) :: num_atom,nel,addup
    
! Public, static variables
  REAL(q2),DIMENSION(3) :: corner,steps
  CHARACTER(LEN=20) :: chargefile

  INTEGER :: ndim,bdim,nrho,wdim

END MODULE vars

