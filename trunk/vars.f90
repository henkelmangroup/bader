  MODULE varsM
    IMPLICIT NONE

    PUBLIC

! Public parameters
    INTEGER,PARAMETER :: q1=SELECTED_REAL_KIND(6,30)         ! single precision
    INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)       ! double precision
    REAL(q2),PARAMETER :: pi=3.141592653589793238462643_q2    
    REAL(q2),PARAMETER :: bader_tol=1.0e-4_q2

! Public, allocatable variables
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: max_rho  
    REAL(q2),ALLOCATABLE,DIMENSION(:,:,:) :: rho
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: bader_charge,voronoi_charge,Rcar,Rdir,dipole,min_dist
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: bader_dist,bader_achg
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: path
    INTEGER,ALLOCATABLE,DIMENSION(:) :: bader_atom,num_atom,nel,addup
    
! Public, static variables
    REAL(q2),DIMENSION(3,3) :: Lattice
    REAL(q2),DIMENSION(3) :: corner,steps
    CHARACTER(LEN=20) :: chargefile


    LOGICAL :: vasp
    INTEGER :: ndim,ngxf,ngyf,ngzf,bdim,nrho,wdim,natypes,na

  END MODULE varsM


