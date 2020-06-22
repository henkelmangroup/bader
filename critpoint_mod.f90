  MODULE critpoints_mod
    USE kind_mod
    USE matrix_mod
    USE bader_mod
!    USE linsolve_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    USE ions_mod
    USE weight_mod
    USE dsyevj3_mod
    IMPLICIT NONE

    PRIVATE 
    PUBLIC :: critpoint_find

    TYPE hessian

      ! du dv dw are derivatives of the three original lattice vectors read from
      ! CHGCAR
      REAL(q2),DIMENSION(3) ::  du, dv, dw
      REAL(q2) :: dudu, dvdv, dwdw, dudv, dudw, dvdw
      ! eigval and eigvec are eigenvalues and eigvectors of hessian matrix
    END TYPE

    TYPE cpc ! stands for critical point candidate
      INTEGER,DIMENSION(3) :: ind  ! these are the indices of the cp
      REAL(q2),DIMENSION(3) :: trueind, truer
      REAL(q2),DIMENSION(3) :: grad
      REAL(q2), DIMENSION(3,3) :: hessianMatrix
      REAL(q2),DIMENSION(3) :: tempcart, tempind
      REAL(q2),DIMENSION(3,3,3) :: dx, dy, dz ! first derivatives of neighbors
!      REAL(q2),DIMENSION(3,3,3) :: du, dv, dw
      REAL(q2),DIMENSION(3) :: du, dv, dw ! 1 2 3 are backward, present, forward
      REAL(q2),DIMENSION(3) :: eigvals, r, cocart, colat, tempr
      INTEGER,DIMENSION(8,3) :: nnind ! indices of neighbors 
      ! indices from 0 to 8 are zyx 000 001 010 011 100 101 110 111
!      REAL(q2),DIMENSION(8,3) :: nngrad ! gradients of nn mentioned above.
!      REAL(q2),DIMENSION(6,3) :: intnngrad ! gradients of interpolated neighbors
      ! used to find interpolated hessians.
      REAL(q2),DIMENSION(3,3) :: eigvecs
      INTEGER :: negcount
      LOGICAL :: hasProxy, isunique
    END TYPE
    CONTAINS

!-----------------------------------------------------------------------------------!
!critpoint_find: find critical points on the edge of the Bader volumes
!NOTE: this subroutine should be called after refine_edge
!      in order to restrict the calculation to edge points
!-----------------------------------------------------------------------------------!
  SUBROUTINE critpoint_find(bdr,chg,opts,ions,stat)
! These are for screening CP due to numerical error. 

    
    TYPE(hessian) :: hes
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(weight_obj), ALLOCATABLE, DIMENSION(:) :: gradlist
    TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpcl, cpl, cpclt ! critical point list and a temporary
    ! The above three are CP candidate list, CP list and CP list temp
! copy
! for points, 1 and 2 are +1, -1
    INTEGER :: stat, setcount, counter
    INTEGER,DIMENSION(3) :: p, pt, ptt, ptx1, ptx2, pty1, pty2, ptz1, ptz2
    INTEGER,DIMENSION(3) :: tempind
    INTEGER :: n1, n2, n3, d1, d2, d3, cptnum, ucptnum, i, j, k, debugnum
    INTEGER :: negcount, bondcount, ringcount, maxcount, cagecount
    INTEGER :: ubondcount, uringcount, umaxcount, ucagecount
    REAL(q2),DIMENSION(3) :: eigvec1, eigvec2, eigvec3, tempVec
    REAL(q2),DIMENSION(3) :: truer, interpolGrad
    REAL(q2),DIMENSION(3) :: grad, prevgrad
    INTEGER, DIMENSION(3) :: tempr
    REAL(q2), DIMENSION(3) :: temprealr, tempreal3d
    ! to be used in newton method in finding unique critical points.
    REAL(q2),DIMENSION(3,3) ::  hessianMatrix, bkhessianMatrix
    ! these are vectors orthogonal to eigenvectors
    REAL(q2),DIMENSION(3) :: tem, tem2a,tem2b,tem2c, eigvals,carts
    REAL(q2) :: threshhold, tempreal, rndr
    REAL(q2),DIMENSION(3,3) :: eigvecs, inverseHessian
    ! linearized approximated derivatives for proxy critical screening
    REAL(q2) :: dx0,dx1,dy0,dy1,dz0,dz1 ! outputs from interpolated gradients
    REAL(q2),DIMENSION(6,3) :: intcarts ! positions in cart of 6 interpolated points
    ! row 1 2 are + and 1 x, then + and - y, then + and - z
    REAL(q2),DIMENSION(6,3) :: intgrads ! gradients of interpolated points
    REAL(q2),DIMENSION(6) :: intrhos ! rhos of interpolated points 
    REAL(q2),DIMENSION(6,3) :: intinds ! fraction indicies for interpolated
    REAL(q2) :: rhocur ! rho of current point
    REAL(q2) :: stepsize, temnormcap
    REAL(q2),DIMENSION(3) :: distance ! vector to 000 in trilinear
    REAL(q2),DIMENSION(3) :: preal
    INTEGER,DIMENSION(8,3) :: nn ! alternative trilinear approx.
    REAL(q2),DIMENSION(8) :: vals
    ! points
    LOGICAL,DIMENSION(3) :: cartcoor ! check if axis are alone cartesian.
    LOGICAL :: invac ! this point is in vacuum
    LOGICAL :: proxy, isReduced, phmrCompliant
    ! The followings are for finding unique critical points
    REAL(q2),DIMENSION(8,3) :: nngrad
    REAL(q2),DIMENSION(8,3,3) :: nnhes !hessian of 8 nn
    INTEGER, DIMENSION(:,:),ALLOCATABLE :: nnind
    REAL(q2),DIMENSION(6,3) :: intnngrad
    REAL(q2),DIMENSION(3,3) :: interpolHessian
    REAL(q2),DIMENSION(3) :: nexttem, previoustem, averager, temcap, temscale
    INTEGER :: averagecount, repeatcount
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: nucleiInd
    REAL(q2), DIMENSION(:,:), ALLOCATABLE :: cpRoster
    INTEGER :: stepcount, nnlayers
    ! The following are for least square calculations
    INTEGER, DIMENSION(3,26) :: vi, matw
    INTEGER, DIMENSION(26,3) :: vit
    REAL(q2), DIMENSION(3,3) :: ggrid
    REAL(q2), DIMENSION(26) :: wi
    REAL(q2), DIMENSION(3,13) :: matwprime
    REAL(q2), DIMENSION(3,3) :: matm, outerproduct
    REAL(q2),DIMENSION(10,3) :: rList,temList
    CHARACTER(128) :: string
    ! below are variables for least sqaures gradient
    stat = 0 ! 0 means nothing
    PRINT *, ''//achar(27)//'[31m Finding Critical points'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[91m Interrogation of the soul:'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[33m Did I turn on vacuum ?'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[92m Did I tell if this is a crystall or molecule?'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[36m Did I use the CHGCAR_sum ?'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[34m Did I use the least square flag -ls ? '//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[95m Critical point is like a box of chocolates. &
              You never know what you are gonna get.'//achar(27)//'[0m'
    !WRITE(*,'(A)')  'FINDING CRITICAL POINTS'
    IF (opts%leastsquare_flag .EQV. .TRUE. ) THEN
      PRINT *, 'Using least square gradient'
    END IF
    ! get expected nucleus indices
!    WRITE (97,*) , 'expecting nuclei at:'
    ALLOCATE(nnInd(8,3))
!    ALLOCATE(nnInd((nnlayers*2+1)**3,3))
    ALLOCATE(nucleiInd(ions%nions,3))
    DO d1 = 1, ions%nions
      nucleiInd(d1,1) = NINT(ions%r_lat(d1,1))
      nucleiInd(d1,2) = NINT(ions%r_lat(d1,2))
      nucleiInd(d1,3) = NINT(ions%r_lat(d1,3))
!      WRITE (97,*) nucleiInd(d1,:)
    END DO
    setcount = 0
    bondcount = 0
    ubondcount = 0
    ringcount = 0
    uringcount = 0
    maxcount = 0
    ucagecount = 0
    cagecount = 0
    ucptnum = 0
    cptnum = 0

!    PRINT *, 'lat2car is'
!    PRINT *, chg%lat2car(1,:)
!    PRINT *, chg%lat2car(2,:)
!    PRINT *, chg%lat2car(3,:)
!    PRINT *, 'lattice is'
!    PRINT *, ions%lattice(1,:)
!    PRINT *, ions%lattice(2,:)
!    PRINT *, ions%lattice(3,:)
!    PRINT *, 'ions locations are at' 
!    DO i = 1, ions%nions
!      PRINT *, ions%r_lat(i,:)
!    END DO

!    PRINT *, 'Printing the matrix '
!    PRINT *, '(/ (/1,2,3/) '
!    PRINT *, '   (/4,5,6/) '
!    PRINT *, '   (/7,8,9/) /)'
    hessianMatrix(1,:)=(/1.,2.,3./)
    hessianMatrix(2,:) = (/4.,5.,6./)
    hessianMatrix(3,:) = (/7.,8.,9./)
!    PRINT *, hessianMatrix
!    PRINT *, 'multiplying the matrix with (/1.,1.,1./)'
!    PRINT *, MATMUL(hessianMatrix,(/1.,1.,1./))

!    PRINT * , "These code requires -vac auto or -vac #"
!    PRINT *, '-----------                  WARNING             -----------'
!    PRINT *, 'Using valence charge may yield useless and confusing results'
!    PRINT *, '    It is recommended to use total charge for finding CPs   '
!    PRINT *, '____________________________________________________________'
!    OPEN(97,FILE='CPF.dat',STATUS='REPLACE',ACTION='WRITE')
!    OPEN(98,FILE='CPFU.dat',STATUS='REPLACE',ACTION='WRITE')
    OPEN(1,FILE='GRAD.dat',STATUS='REPLACE',ACTION='WRITE')
    WRITE(1,*) 'norm      ','| del a      | ','    del b     |','    del c'
    OPEN(2,FILE='HES.dat',STATUS='REPLACE',ACTION='WRITE')
!    OPEN(3,FILE='CPFU.dat',STATUS='REPLACE',ACTION='WRITE')
!    OPEN(4,FILE='CPFU.dat',STATUS='REPLACE',ACTION='WRITE')
    debugnum = 0

    ! check if axis are cartesian
!    cartcoor = coorcheck(ions%lattice)
    IF ( opts%leastsquare_flag .EQV. .TRUE. )THEN
      vi = makevi()
      vit = TRANSPOSE(vi)
      ggrid = makeggrid(chg,ions)
!      PRINT *, 'ggrid is'
!      PRINT *, ggrid(1,:)
!      PRINT *, ggrid(2,:)
!      PRINT *, ggrid(3,:)
!      PRINT *, 'wi is'
!      DO i = 1,26
!        wi(i) = DOT_PRODUCT(MATMUL(vit(i,:),ggrid),vi(:,i))
!        wi(i) = 1./wi(i)
!      END DO
      ! The following method checks wi directly. The above code is working
      ! properly.
      !wi = makewi(chg,nnlayers,vi)
      !PRINT *, 'new wi are'
      !DO i = 1,26
      !  PRINT *, wi(i)
      !END DO

      matm = 0.0_q2
!      PRINT *, 'matwprime is'
!      DO i = 1,13
!        matwprime(:,i) = vi(:,i) * wi(i)
!        PRINT *, matwprime(:,i)
!      END DO
      DO i = 1, 26
        outerproduct(1,1)= vi(1,i) * vit(i,1)
        outerproduct(1,2)= vi(1,i) * vit(i,2)
        outerproduct(1,3)= vi(1,i) * vit(i,3)
        outerproduct(2,1)= vi(2,i) * vit(i,1)
        outerproduct(2,2)= vi(2,i) * vit(i,2)
        outerproduct(2,3)= vi(2,i) * vit(i,3)
        outerproduct(3,1)= vi(3,i) * vit(i,1)
        outerproduct(3,2)= vi(3,i) * vit(i,2)
        outerproduct(3,3)= vi(3,i) * vit(i,3)
!        PRINT *, 'i is', i
!        PRINT *, 'vi is'
!        PRINT *, vi(:,i)
!        PRINT *, 'outer product is'
!        PRINT *, outerproduct(1,:)
!        PRINT *, outerproduct(2,:)
!        PRINT *, outerproduct(3,:)
!        PRINT *, 'wi vi x vi is'
!        PRINT *, outerproduct * wi(i)
        matm = matm + wi(i) * outerproduct
!        PRINT *, 'matm is now'
!        PRINT *, matm(1,:)
!        PRINT *, matm(2,:)
!        PRINT *, matm(3,:)
!        PRINT *, 'vi is' 
!        PRINT *, vi(:,i)
!        PRINT *, 'vit is'
!        PRINT *, vit(i,:)
!        PRINT *, 'outer product is'
!        PRINT *, outerproduct

      END DO
    END IF
    nnlayers = findnnlayers(ions)
    nnind = findnn(truer,nnlayers,chg,ions)



!    PRINT *, 'density run around the center on y=x direction '
!    ! This is the on grid part
!    p = (/1,7,7/)
!    !p = (/4,4,10/)
!    DO n1 = 1,66
!      CALL pbc(p,chg%npts)
!      PRINT *, rho_val(chg,p(1),p(2),p(3))
!!      PRINT *, 'getting density at'
!!      PRINT *, p
!!      PRINT *, 'gradient is'
!      p = p + (/1,1,0/)
!    END DO
!    PRINT *, 'printing debug output for hex he case'
!    PRINT *, 'Printing on grid gradient from 07 43 63 to 43 07 08'
!    truer = (/17.,33.,53./)
!    p = (/17,33,53/)
!    DO n1 = 1,200
!      CALL pbc(p,chg%npts)
!      !PRINT *, rho_val(chg,p(1),p(2),p(3))
!      PRINT *, CDGrad(p,chg)
!      truer = truer + (/0.08,-0.08,-0.175/)
!      p(1) = NINT(truer(1))
!      p(2) = NINT(truer(2))
!      p(3) = NINT(truer(3))
!    END DO
!
!    PRINT *, 'Printing on grid gradient from 07 43 63 to 43 07 08'
!    truer = (/17.,33.,53./)
!    p = (/17,33,53/)
!    DO n1 = 1,200
!      CALL pbc(p,chg%npts)
!      !PRINT *, rho_val(chg,p(1),p(2),p(3))
!      PRINT *, CDGrad(p,chg)
!      truer = truer + (/0.08,-0.08,-0.175/)
!      p(1) = NINT(truer(1))
!      p(2) = NINT(truer(2))
!      p(3) = NINT(truer(3))
!    END DO


!    PRINT *, 'Printing on grid gradient from 1 28 1 to 28 1 1'
!    p = (/1,28,1/)
!    DO n1 = 1,27
!      CALL pbc(p,chg%npts)
!      !PRINT *, rho_val(chg,p(1),p(2),p(3))
!      PRINT *, CDGrad(p,chg)
!      p = p + (/1,-1,0/)
!    END DO

!    PRINT *, 'These are interpolated values'
!!    ! These will give you interpolated results
!    truer = (/2.,16.,30./)
!    DO n1 = 1, 280
!      CALL pbc_r_lat(truer,chg%npts)
!      nnind = FindNN(truer,nnLayers,chg,ions)
!      !nnind = SimpleNN(truer,chg)
!      PRINT *, R2RhoInterpol(nnind,truer,chg,nnLayers)
!      truer = truer + (/0.2,0.1,0./)
!      IF (n1==2) STOP
!    END DO
!    

!    PRINT *, 'printing interpolated gradient from 17 33 53 to 33 17 18'
!    truer = (/17.,33.,53./)
!    DO n1 = 1, 200
!      CALL pbc_r_lat(truer,chg%npts)
!      !nnind = FindNN(truer,nnLayers,chg,ions)
!      !PRINT *, trilinear_interpol_rho(chg,truer)
!      nnind = SimpleNN(trueR,chg)
!      DO i = 1, 8
!        nnGrad(i,:) = CDGrad(nnInd(i,:),chg)
!      END DO 
!      distance = trueR - nnInd(1,:)
!      PRINT *, trilinear_interpol_grad(nnGrad,distance)
!      truer = truer + (/0.08,-0.08,-0.175/)
!    END DO

!    PRINT *, 'printing tem from 17 33 53 to 33 17 18'
!    truer = (/17.,33.,53./)
!    DO n1 = 1, 200
!      CALL pbc_r_lat(truer,chg%npts)
!      !nnind = FindNN(truer,nnLayers,chg,ions)
!      !PRINT *, trilinear_interpol_rho(chg,truer)
!      nnind = SimpleNN(trueR,chg)
!      DO i = 1, 8
!        nnGrad(i,:) = CDGrad(nnInd(i,:),chg)
!        nnHes(i,:,:) = CDHessian(nnInd(i,:),chg)
!      END DO 
!      distance = trueR - nnInd(1,:)
!      grad = trilinear_interpol_grad(nnGrad,distance)
!      hessianMatrix = trilinear_interpol_hes(nnHes,distance)
!      PRINT *, -MATMUL(INVERSE(hessianMatrix),grad)
!      !PRINT *, trilinear_interpol_grad(nnGrad,distance)
!      truer = truer + (/0.08,-0.08,-0.175/)
!    END DO

!    PRINT *, 'printing Hessian from 17 33 53 to 33 17 18'
!    truer = (/17.,33.,53./)
!    DO n1 = 1, 200
!      CALL pbc_r_lat(truer,chg%npts)
!      !nnind = FindNN(truer,nnLayers,chg,ions)
!      !PRINT *, trilinear_interpol_rho(chg,truer)
!      nnind = SimpleNN(trueR,chg)
!      DO i = 1, 8
!        nnHes(i,:,:) = CDHessian(nnInd(i,:),chg)
!      END DO 
!      distance = trueR - nnInd(1,:)
!      hessianMatrix = trilinear_interpol_hes(nnHes,distance)
!      PRINT *, hessianMatrix(1,:)
!      PRINT *, hessianMatrix(2,:)
!      PRINT *, hessianMatrix(3,:)
!      truer = truer + (/0.08,-0.08,-0.175/)
!    END DO

!    PRINT *, 'printing interpolated charge from 17 33 53 to 33 17 18'
!    truer = (/17.,33.,53./)
!    DO n1 = 1, 200
!      CALL pbc_r_lat(truer,chg%npts)
!      !nnind = FindNN(truer,nnLayers,chg,ions)
!      PRINT *, trilinear_interpol_rho(chg,truer)
!      !nnind = SimpleNN(trueR,chg)
!      !DO i = 1, 8
!      !  nnGrad(i,:) = CDGrad(nnInd(i,:),chg)
!      !END DO
!      !distance = trueR - nnInd(1,:)
!      !PRINT *, trilinear_interpol_grad(nnGrad,distance)
!      truer = truer + (/0.08,-0.08,-0.175/)
!    END DO


!
!    PRINT *, 'tem run'
!
!    truer = (/1.,7.,7./)
!    DO n1 = 1, 660
!      CALL pbc_r_lat(truer,chg%npts)
!      nnind = simpleNN(truer,chg)
!      DO n2 = 1, 8
!        nngrad(n2,:) = CDGrad(nnind(n2,:),chg)
!        nnhes(n2,:,:) = CDHessian(nnind(n2,:),chg)
!      END DO
!      distance = truer - nnind(1,:)
!      grad = trilinear_interpol_grad(nngrad,distance) 
!      hessianMatrix = trilinear_interpol_hes(nnhes,distance)
!      PRINT *, - MATMUL(INVERSE(hessianMatrix),grad)
!      truer = truer + (/0.1,0.1,0./)
!    END DO

!
!    PRINT *, 'central difference gradient is'
!    !p = (/10,4,10/)
!    p = (/chg%npts(1)*0.2,chg%npts(2)*0.2,chg%npts(3)*0.5/)
!    !p = (/chg%npts(1)*0.5,chg%npts(2)*0.2,chg%npts(3)*0.5/)
!    DO n1 = 1, chg%npts(1)*0.6
!      PRINT *, CDGrad(p,chg)
!      !p = p + (/0,1,0/)
!      p = p + (/1,1,0/)
!    END DO
!    
!    PRINT *, 'central difference hessian is'
!    !p = (/10,4,10/)
!    p = (/chg%npts(1)*0.2,chg%npts(2)*0.2,chg%npts(3)*0.5/)
!    !p = (/chg%npts(1)*0.5,chg%npts(2)*0.2,chg%npts(3)*0.5/)
!    DO n1 = 0, chg%npts(1)*0.6
!      hessianMatrix = CDHessian(p,chg)
!      PRINT *, hessianMatrix(1,:)
!      PRINT *, hessianMatrix(2,:)
!      PRINT *, hessianMatrix(3,:)
!      !p = p + (/0,1,0/)
!      p = p + (/1,1,0/)
!    END DO
!    STOP
    
    PRINT *, '******************************************************' 

    ALLOCATE (cpl(10000)) ! start with 10000 capacity
    ALLOCATE (cpcl(10000))
      ! Maxima should always be found first with gradient ascension
    PRINT *, 'Finding maxima'
    DO n1 = 1, ions%nions
      tempr(1) = NINT(ions%r_dir(n1,1) * chg%npts(1)) 
      tempr(2) = NINT(ions%r_dir(n1,2) * chg%npts(2)) 
      tempr(3) = NINT(ions%r_dir(n1,3) * chg%npts(3)) 
      temprealr = ascension(tempr,chg,matm,matwprime, &
                 wi,vi,vit,ggrid,outerproduct,opts,nnLayers,ions)
      ucptnum = ucptnum + 1
      cpl(ucptnum)%trueind = temprealr
      ! these points are to stay in the list so cpl is used instead of cpcl
      p(1) = NINT(cpl(ucptnum)%trueind(1))
      p(2) = NINT(cpl(ucptnum)%trueind(2))
      p(3) = NINT(cpl(ucptnum)%trueind(3))
!      PRINT *, 'ascension has point ', p
!      CALL RecordCP(p,chg,matm,matwprime,wi,vi,vit,ggrid &
!        ,outerproduct,cpl,ucptnum,eigvals,eigvecs, maxcount, uringcount, &
!        ubondcount, ucagecount, opts)
      CALL RecordCPR(temprealr,chg,cpl,ucptnum,eigvals,eigvecs, maxcount, uringcount, &
        ubondcount, ucagecount,opts,grad,hessianMatrix)
      ! without interpolation, atomic maxima may not be correctly characterized. 
      ! The solution is to manually correct it
!      PRINT *, 'ucptnum is', ucptnum
      IF (opts%noInterpolation_flag) THEN
        maxcount = ucptnum
        uringcount = 0
        ubondcount = 0
        ucagecount = 0
        cpl(ucptnum)%grad = 0
        cpl(ucptnum)%hessianMatrix = 0
        cpl(ucptnum)%negcount = 3
      END IF
    END DO
!    PRINT *, 'CPs after ascension is, max, bond, ring, min'
!    PRINT *, maxcount, ubondcount, uringcount, ucagecount
    ! ascension results may not give the best hessian. so the types are set
    ! manually.
    IF (maxcount /= ucptnum) THEN
      PRINT *, 'WARNING: It was detected that the number of maxima found &
        does not equal to trials started. The found critical points are &
        being manually set as nuclear critical points'
      DO n1 = 1, ucptnum
        cpl(n1)%negcount = 3
      END DO
      ubondcount = 0
      ucagecount = 0
      uringcount = 0
      maxcount = ucptnum
    END IF
    PRINT *, 'after adjustments, CP counts are'
    PRINT *, maxcount, ubondcount, uringcount, ucagecount

    !This part is on hold for the time being 
    !IF (opts%noInterpolation_flag) THEN
    !  ! It has been shown that without interpolation and Newton's method, we
    !  ! can't really find the CPs so these code should never be run. 
    !  PRINT *, 'find critical points without interpolations'
    !  PRINT *, '--------------------------------------------------'
    !  PRINT *, 'WITHOUT INTERPOLATION AND NEWTONS METHOD'
    !  PRINT *, 'CRITICAL POINTS CANNOT ALL BE FOUND WITHOUT ERRORS'
    !  PRINT *, '--------------------------------------------------'
    !  CALL SLEEP(10)
    !  PRINT *, 'YOU HAVE BEEN WARNED. PROCCEED AT YOUR OWN COST'
    !  CALL SLEEP(5)
    !  counter = 0
    !  ! first find out how many none vacuum points there are
    !  DO n1 = 1, chg%npts(1) 
    !    DO n2 = 1, chg%npts(2)
    !      DO n3 = 1, chg%npts(3)
    !        ! ignore points in vacuum
    !        IF (bdr%volnum(n1,n2,n3) == bdr%bnum + 1) THEN
    !          CYCLE
    !        END IF
    !        counter = counter + 1
    !      END DO
    !    END DO 
    !  END DO  
    !  PRINT *, 'found number of points to analyze: ', counter
    !  ALLOCATE(gradlist(counter))
    !  counter = 0
    !  ! now get the gradient of the points
    !  DO n1 = 1, chg%npts(1) 
    !    DO n2 = 1, chg%npts(2)
    !      DO n3 = 1, chg%npts(3)
    !        ! ignore points in vacuum
    !        IF (bdr%volnum(n1,n2,n3) == bdr%bnum + 1) THEN
    !          CYCLE
    !        END IF
    !        counter = counter + 1
!   !         PRINT *, lsg( &
!   !               (/n1,n2,n3/),chg,matm,matwprime, &
!   !               wi,vi,vit,ggrid,outerproduct)
!   !         PRINT *, mag( lsg( &
!   !               (/n1,n2,n3/),chg,matm,matwprime, &
!   !               wi,vi,vit,ggrid,outerproduct) &
!   !               )
    !        IF (opts%leastsquare_flag) THEN
    !          gradlist(counter)%rho = &
    !            mag( lsg( &
    !              (/n1,n2,n3/),chg,matm,matwprime, &
    !              wi,vi,vit,ggrid,outerproduct) &
    !              )
    !        ELSE
    !          gradlist(counter)%rho = &
    !            mag(CDGrad((/n1,n2,n3/),chg))
    !        END IF
    !        gradlist(counter)%pos=(/n1,n2,n3/)
    !      END DO
    !    END DO 
    !  END DO  
    !  PRINT *, 'calculated all gradients'
    !  CALL sort_weight(counter,gradlist)
    !  PRINT *, 'sorted'
    !  ! starting from the smallest gradient, characterize each point
    !  DO i = 0, counter - 1
    !    p = (/gradlist(counter - i)%pos(1), &
    !          gradlist(counter - i)%pos(2), &
    !          gradlist(counter - i)%pos(3)/)
    !    ucptnum = ucptnum + 1
    !    CALL RecordCP(p,chg,matm,matwprime,wi,vi,vit,ggrid &
    !        ,outerproduct,cpl,ucptnum,eigvals,eigvecs, maxcount, &
    !        uringcount,ubondcount,ucagecount, opts)
    !    PRINT *, 'ucptnum is', ucptnum
    !    PRINT *, gradlist(counter - i)%rho
    !    PRINT *, cpl(ucptnum)%negcount
    !    proxy = .FALSE.
!   !     PRINT *, 'point ',p, 'has gradient'
!   !     PRINT *, gradlist(counter - i)%rho
    !    ! start with already known CPs
    !    DO j = 1, ucptnum
    !      pt = cpl(j)%ind
    !      IF (GetPointDistance(p,pt,chg,nnlayers)<= 0.2) THEN
    !        proxy = .TRUE.
    !        EXIT
    !      END IF
    !    END DO
    !    ! IF this point is close to one already insepcted, skip this point
    !    DO j = 1, i-1
    !      pt = (/gradlist(counter - j)%pos(1), &
    !          gradlist(counter - j)%pos(2), &
    !          gradlist(counter - j)%pos(3)/)
    !      IF (GetPointDistance(p, pt, chg, nnlayers) <= 0.2 ) THEN
    !        proxy = .TRUE.
    !        EXIT
    !      END IF
    !    END DO
    !    IF (.NOT. proxy) THEN
    !      ! a new CP
    !      ucptnum = ucptnum + 1
    !      ! characterize and record the CP
    !      CALL RecordCP(p,chg,matm,matwprime,wi,vi,vit,ggrid &
    !        ,outerproduct,cpl,ucptnum,eigvals,eigvecs, maxcount, &
    !        uringcount,ubondcount,ucagecount, opts)
    !      PRINT *, 'ucptnum is', ucptnum
    !      IF (phRuleChecker(maxcount, ubondcount, uringcount, ucagecount,&
    !        opts, ions )) THEN
    !        setcount = setcount + 1
    !        PRINT *, 'set count is ', setcount
    !        PRINT *, 'counts are, max, bond, ring, cage:'
    !        PRINT *, maxcount, ubondcount, uringcount, ucagecount
    !        CALL outputCP(cpl,opts,ucptnum,chg,setcount,ubondcount, &
    !          uringcount, ucagecount, maxcount)
    !      END IF
    !    END IF
    !  END DO 
    !  
    !  DEALLOCATE(gradlist)
    !END IF
    IF (.NOT. opts%noInterpolation_flag) THEN
      CALL minimasearch(chg,cptnum,cpl,bdr,matm,matwprime, &
                            wi,vi,vit,ggrid,outerproduct,opts,ucagecount,&
                            nnLayers,ions)
      !DO n3 = 1, chg%npts(3)
      !  DO n2 = 1, chg%npts(2)
      !    DO n1 = 1, chg%npts(1)
      !          p = (/n1,n2,n3/)
      !          IF (opts%leastSquare_flag) THEN
      !            grad = lsg(p,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      !            hessianMatrix = &
      !              lsh(p,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      !          ELSE 
      !            grad = CDGrad(p,chg)
      !            hessianMatrix = CDHessian(p,chg)
      !          END IF
      !          WRITE(1,*) p
      !          WRITE(1,*) SQRT(SUM(grad*grad))
      !          WRITE(1,*) grad
      !          WRITE(2,*) p
      !          WRITE(2,*) hessianMatrix(1,:)
      !          WRITE(2,*) hessianMatrix(2,:)
      !          WRITE(2,*) hessianMatrix(3,:)
      !    END DO
      !  END DO
      !END DO
      DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
              ! check to see if this point is in the vacuum
            IF (bdr%volnum(n1,n2,n3) == bdr%bnum + 1) THEN
              debugnum = debugnum + 1
              CYCLE
            END IF
            ! check to see if this point is on the boundary
!            IF ( .NOT. is_vol_edge(bdr,chg,(/n1,n2,n3/))) THEN
!              CYCLE
!            END IF
            p = (/n1,n2,n3/)
            ! check to see if a critical point is already near by
            ! this should only speed things up and have no effect on accuracy
            preal = p
            IF (.FALSE.) THEN
            ! ! This determins if the initial step uses newton wrapping to determine
            ! ! candidency or other
            !   IF ( gradientfilter(p,chg,matm, &
            !        matwprime,wi,vi,vit,ggrid,outerproduct) ) THEN
            !     ! get the gradient and hessian to do initial characterization of
            !     ! the critical point
            !     IF (opts%leastSquare_flag) THEN
            !       grad = lsg(p,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
            !       hessianMatrix = &
            !         lsh(p,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
            !     ELSE
            !       grad = CDGrad(p,chg)
            !       hessianMatrix = CDHessian(p,chg)
            !     END IF
            !     CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
            !     negcount = 0
            !     DO i = 1,3
            !       IF (eigvals(i) < 0) THEN
            !         negcount = negcount + 1
            !       END IF
            !     END DO
            !     IF (negcount == 3) THEN
            !       cagecount = cagecount + 1
            !       PRINT *, 'Found a unique cage critical point'
            !        WRITE(98,*) 'This is a cage critical point'
            !        WRITE(98,*) ' '
            !     END IF
            !     IF (negcount == 2) THEN
            !       bondcount = bondcount + 1
            !       PRINT *, 'Found a unique bond critical point'
            !        WRITE(98,*) 'This is a bond critical point'
            !        WRITE(98,*) ' '
            !     ELSEIF(negcount == 1) THEN
            !       ringcount = ringcount + 1
            !       PRINT *, 'Found a unique ring critical point'
            !        WRITE(98,*) 'This is a ring critical point'
            !        WRITE(98,*) ' '
            !     ELSEIF(negcount == 0) THEN
            !       maxcount = maxcount + 1
            !       PRINT *, 'Found a unique nuclear critical point'
            !        WRITE(98,*) 'This is a nuclear critical point'
            !        WRITE(98,*) ' ' 
            !       cptnum = cptnum + 1
            !       PRINT *, cptnum
            !     END IF
            !   END IF
            ELSE  
            ! Below starts the defualt screening of candidency
  !-  ----------------------------------------------------------------------------------!
  !   now that this subroutine can find the correct amount of edge points, lets have
  !   it find the hessian
  !-  ----------------------------------------------------------------------------------!
              IF (opts%leastsquare_flag .EQV. .TRUE.) THEN
                grad = lsg(p,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
                hessianMatrix = &
                  lsh(p,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
                tem = - MATMUL(INVERSE(hessianMatrix),grad)
                ! tem is now in cartesian. convert it back to lattice
                tem = MATMUL(chg%car2lat,tem)
              ELSE 
                ! use central difference
                grad = CDGrad(p,chg)
                hessianMatrix = CDHessian(p,chg)
                tem = - MATMUL(INVERSE(hessianMatrix),grad)
                ! this tem is in cartesian. Transform it to lattice 
                tem = MATMUL(chg%car2lat,tem)
!                IF (p(1)==66.AND.p(2)==66.AND.p(3)==10) THEN
!                  PRINT *,'66 66 10'
!                  PRINT *, 'grad is'
!                  PRINT *, grad
!                  PRINT *, 'hessian is'
!                  PRINT *, hessianMatrix
!                  PRINT *, 'tem is'
!                  PRINT *, tem
!                END IF
              END IF
!              IF (mag(grad) >= 50 ) THEN
                ! This threshhold should be more flexible.
!                CYCLE
!              END IF
              IF ( (ABS(tem(1)) <= 1 + opts%knob_tem .AND. &
                   ABS(tem(2)) <= 1 + opts%knob_tem .AND. &
                   ABS(tem(3)) <= 1 + opts%knob_tem) .OR. &
                  (SUM(grad*grad) <= (0.1*opts%knob_gradfloor)**2 )) THEN              
                ! finding proximity could potentially be costly
                IF (ProxyToCPCandidate(p,opts,cpcl,cptnum,chg,nnLayers)) THEN
!                  PRINT *, 'for proximity to a candidate, skipped '
!                  PRINT *, n1,n2,n3
                  CYCLE
                END IF
!                PRINT *, 'candidate found at'
!                PRINT *, p
!                PRINT *, 'gradient is'
!                PRINT *, grad
!                PRINT *, 'tem is'
!                PRINT *, tem
                cptnum = cptnum + 1
!                PRINT *, 'cptnum updated to'
!                PRINT *, cptnum
                bkhessianMatrix = hessianMatrix
                ! Check if the candidate list needs to be expanded.
                IF (cptnum < SIZE(cpcl) - 1 ) THEN
                  cpcl(cptnum)%du = hes%du
                  cpcl(cptnum)%dv = hes%dv  
                  cpcl(cptnum)%dw = hes%dw
                  cpcl(cptnum)%ind(1) = n1
                  cpcl(cptnum)%ind(2) = n2
                  cpcl(cptnum)%ind(3) = n3
                  cpcl(cptnum)%grad = grad
                  cpcl(cptnum)%hasProxy = .FALSE.
                  cpcl(cptnum)%r = tem
                  cpcl(cptnum)%tempcart = MATMUL(chg%car2lat,tem + p)
                ELSE 
                  PRINT *, 'expanding cpcl size'
                  ALLOCATE(cpclt(cptnum + 1))
                  DO i = 1, cptnum - 1
                    cpclt(i) = cpcl(i)
                  END DO
                  DEALLOCATE(cpcl)
!                  PRINT *, 'copied'
                  ALLOCATE(cpcl(cptnum*2))
                  PRINT *, 'cpcl size now is', SIZE(cpcl)
                  DO i = 1, cptnum - 1
                    cpcl(i)=cpclt(i)
                  END DO
!                  PRINT *, 'copied back'
                  DEALLOCATE(cpclt)
                  cpcl(cptnum)%du = hes%du 
                  cpcl(cptnum)%dv = hes%dv
                  cpcl(cptnum)%dw = hes%dw
                  cpcl(cptnum)%ind(1) = n1
                  cpcl(cptnum)%ind(2) = n2
                  cpcl(cptnum)%ind(3) = n3
                  cpcl(cptnum)%grad = grad
                  cpcl(cptnum)%hasproxy = .FALSE.
                  cpcl(cptnum)%r = tem
                  cpcl(cptnum)%tempcart = MATMUL( chg%car2lat,tem + p)
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO
      PRINT *, "Number of Newton Rhapson trajectory needed: ", cptnum 
!!**  *****************************************************************
      ! To find critical points (unique), start with a cell that contains a
      ! critical point and its hessian and force. Use Newton's method to make a
      ! move. Interpolate the force inside the voxel. 
      ! Once moved, get the new force through trilinear interpolation, and
      ! get the new hessian which will be a matrix of constants, make moves until
      ! r is zero. get the coordinates of the new true critical point. If this
      ! point is within half lattice to another, do not record this new point.
!      WRITE(98,*) '***************************************'
!      WRITE(98,*) 'Below are points found through Newton Raphson'
!      WRITE(98,*) '***************************************'
      IF (.TRUE.) THEN
        ALLOCATE(cpRoster(cptnum,3))
      ! This turns critical point candidency validation on and off
      ! ??? exo me?
        DO i = 1, cptnum
          cpcl(i)%isunique = .FALSE.
          temcap = (/1.,1.,1./)
          temscale = (/1.,1.,1./)
          temnormcap = 1.
!          PRINT *, 'this is candidate number ', i
!          PRINT *, 'looking at candidate of indices'
!          PRINT *, cpcl(i)%ind
          IF (.FALSE.) THEN
          ! This determins if validation is done with gradient descend
          !    PRINT *, 'looking at critical point candidate # ', i
          !    PRINT *, 'indices are '
          !   p = cpcl(i)%ind
          !    PRINT *, p
          !   cpcl(i)%trueind = 0.
          !    PRINT *, 'starting sqgradientdescend'
          !   CALL sqgradientdescend(p,chg,matm,matwprime,wi,vi,vit, &
          !                      ggrid,outerproduct,opts,cpcl(i)%trueind,invac,bdr,nnLayers,ions)
          !    PRINT *, 'ending sqgradientdescend'
          !    PRINT *, 'after descend trueind is'
          !    PRINT *, cpcl(i)%trueind
          !   IF (cpcl(i)%trueind(1) == -1.) THEN
          !     cpcl(i)%isunique = .FALSE.
          !     PRINT *, 'set to be not unique'
          !     CYCLE
          !   END IF
          !   truer = cpcl(i)%trueind
          !   cpcl(i)%isunique = .TRUE.
          !   nnind = simpleNN(truer,chg)
          !   distance = truer - nnind(1,:)
          !   DO j = 1,8
          !     IF (opts%leastsquare_flag .EQV. .true.) THEN
          !       nngrad(j,:) = lsg(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
          !       hessianmatrix = lsh(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
          !       inversehessian = inverse(hessianmatrix)
          !       nnhes(j,1,:) = hessianmatrix(1,:)
          !       nnhes(j,2,:) = hessianmatrix(2,:)
          !       nnhes(j,3,:) = hessianmatrix(3,:)
          !     ELSE 
          !       nngrad(j,:) = CDGrad(p,chg)
          !       hessianMatrix = CDHessian(p,chg)
          !       nnhes(j,:,:) = hessianMatrix(:,:)
          !     END IF
          !   END DO
          !   interpolgrad = trilinear_interpol_grad(nngrad,distance) ! val r interpol
          !   interpolHessian = trilinear_interpol_hes(nnhes,distance)
          !   !nnind = FindNN(truer,nnlayers,chg,ions)
          !   !interpolGrad = R2GradInterpol(nnind,truer,chg,nnlayers)
          !   !interpolHessian = R2HesInterpol(nnind,truer,chg,nnlayers)
          ELSE
          ! Begins newton raphson validation process
!            PRINT *, 'looking at candidate ', i
!            PRINT *, 'starting position is'
!            PRINT *, cpcl(i)%ind
            stepcount = 1
            averagecount = 0
            averageR = (/-1.,-1.,-1./)
            n1 = cpcl(i)%ind(1)
            n2 = cpcl(i)%ind(2)
            n3 = cpcl(i)%ind(3)
            cpcl(i)%nnind = simpleNN(cpcl(i)%r,chg)
            ! get gradients " of all nearest neighbors
            DO j = 1, 8
              IF (opts%leastsquare_flag .EQV. .TRUE.) THEN
                nngrad(j,:) = lsg(cpcl(i)%nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
              ELSE 
                nngrad(j,:) = CDGrad(cpcl(i)%nnind(j,:),chg)
              END IF
              ! this force is not adjusted with lattice
              !cpl(i)%nngrad(j,:) = MATMUL(chg%car2lat,cpl(i)%nngrad(j,:))
              !PRINT *, cpl(i)%nngrad(j,:)
              ! now things are cartesian. forces are checked to be correct.
              ! correct: if I make things cartesian, forces are checked out to be
              ! correct. However I still need to find hessian in lattice first. 
              ! So I do not convert force to cartesian yet.
              ! force is checked to be correct.
            END DO
            ! Now start newton method iterations
            truer =  cpcl(i)%ind + cpcl(i)%r
            CALL pbc_r_lat(truer,chg%npts)
            previoustem = cpcl(i)%r
            DO stepcount = 1,1000
              !IF (truer(1) - cpcl(i)%ind(1) >= 2 + 2 * opts%knob_tem .OR. &
              !    truer(2) - cpcl(i)%ind(2) >= 2 + 2 * opts%knob_tem .OR. &
              !    truer(3) - cpcl(i)%ind(3) >= 2 + 2 * opts%knob_tem) THEN
              !  PRINT *,'too far'
              !  EXIT
              !END IF
!              IF (n1==23.AND.n2==25.AND.n3==35) THEN
!                PRINT *, 'truer is '
!                PRINT *, truer
!                PRINT *, 'grad is'
!                PRINT *, grad
!                PRINT *, 'nexttem is'
!                PRINT *, nexttem
!              END IF 
              IF (stepcount >= 1) THEN
                prevgrad = grad
              END IF
              CALL pbc_r_lat(truer,chg%npts)
              nnind = simpleNN(truer,chg)
              distance = truer - nnind(1,:)
              DO j = 1,8
                IF (opts%leastsquare_flag .EQV. .true.) THEN
                  nngrad(j,:) = lsg(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
                  hessianmatrix = lsh(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
                  inversehessian = inverse(hessianmatrix)
                  nnhes(j,1,:) = hessianmatrix(1,:)
                  nnhes(j,2,:) = hessianmatrix(2,:)
                  nnhes(j,3,:) = hessianmatrix(3,:)
                ELSE 
                  nngrad(j,:) = CDGrad(nnind(j,:),chg)
                  nnhes(j,:,:) = CDHessian(nnind(j,:),chg)
                END IF
              END DO
               !row find the nearest neighbors at this new locaiton
               !first update critical point location
               !the next big step is to interpolate the force at predicted critical
               !point.
              grad = trilinear_interpol_grad(nnGrad,distance) ! val r interpol
              interpolHessian = trilinear_interpol_hes(nnHes,distance)
              !grad = R2GradInterpol(nnind,truer,chg,nnLayers)
              IF (ABS(grad(1)) <= 0.1*opts%knob_gradfloor .AND. &
                  ABS(grad(2)) <= 0.1*opts%knob_gradfloor .AND. &
                  ABS(grad(3)) <= 0.1*opts%knob_gradfloor) THEN

                cpcl(i)%trueind = truer 
                cpcl(i)%isunique = .TRUE.
                EXIT
              END IF
              nexttem = - MATMUL(inverse(interpolHessian),grad)
              nexttem = MATMUL(chg%car2lat,nexttem)
              ! see if movement capping is necessary
              nexttem = TemMods(nexttem,temscale,temnormcap)
              temscale = scaleinspector( nexttem, previoustem, temscale)
              CALL DetectCircling(stepCount,rList,temList,trueR,nextTem,averageR)
              IF (ALL(averageR /= -1.,1)) THEN
                !PRINT *, 'averageR is'
                !PRINT *, averageR
                cpcl(i)%isUnique = .TRUE.
                cpcl(i)%trueInd = averageR
                trueR = averageR
                !PRINT *, 'averaged location at'
                !PRINT *, averageR
                ! the following code is temporary. it disables averaging.
                cpcl(i)%isUnique = .FALSE.
                EXIT
              END IF
              !nexttem = unstuck( nexttem, previoustem, temscale, temnormcap)
              ! There is no reason why the following should work
              !IF (prevgrad(1)*grad(1)<=0) THEN
              !  temcap(1) = temcap(1) * 0.5
              !END IF
              !IF (prevgrad(2)*grad(2)<=0) THEN
              !  temcap(2) = temcap(2) * 0.5
              !END IF
              !IF (prevgrad(3)*grad(3)<=0) THEN
              !  temcap(3) = temcap(3) * 0.5
              !END IF
              previoustem = nexttem
              tempr(1) = NINT(truer(1))
              tempr(2) = NINT(truer(2))
              tempr(3) = NINT(truer(3))
              CALL pbc(tempr,chg%npts)
!              IF (stepcount >= 50) THEN
!                ! take the average over 100 steps as the critical point position.
!                averager = averager + truer
!                averagecount = averagecount + 1
!              END IF
!              IF (stepcount == 100) THEN
!                averager = averager / averagecount
!                PRINT *, 'averager averaged to be', averager
!                truer = averager
!                nexttem = 0.
!              END IF
              IF (bdr%volnum(tempr(1), &
                  tempr(2),tempr(3)) == bdr%bnum + 1) THEN
                ! We are heading into the vacuum space, cosmonaughts! 
                cpcl(i)%isunique = .FALSE.
!                PRINT *, 'not unique because it wondered off to vacuum'
                EXIT
              END IF
              IF ( ABS(nexttem(1)) <= 0.1*opts%knob_newtonr .AND. &
                   ABS(nexttem(2)) <= 0.1*opts%knob_newtonr .AND. &
                   ABS(nexttem(3)) <= 0.1*opts%knob_newtonr ) THEN
                cpcl(i)%trueind = truer 
                cpcl(i)%isUnique = .TRUE.
!                PRINT *, 'vector small'
!                PRINT *, nexttem
!                IF (mag(grad)>=0.1) THEN 
!                  PRINT *, 'gradient still over 0.1'
!                  cpcl(i)%isunique = .FALSE.
!                END IF
                EXIT
              END IF
              truer = truer + nexttem
            END DO
          END IF
!          PRINT *, 'stepcount : ', stepcount
!          PRINT *, 'residue gradient is'
!          PRINT *,  grad
!          PRINT *, 'last nexttem is'
!          PRINT *,  nexttem
!          IF (mag(grad) >= 0.1) THEN
!            PRINT *, 'flagged not unique for residue gradient'
!            cpcl(i)%isunique = .FALSE.
!          END IF
          IF (cpcl(i)%isunique ) THEN
            CALL MakeCPRoster(cpRoster,i,truer)
            ! check distance to existing critical points.
!            DO j = 1, i - 1 
!              IF (GetPointDistanceR(truer, cpcl(j)%trueind, &
!                  chg, nnlayers) <= 0.2) THEN
!                PRINT *, 'not unique because a duplicate of'
!                PRINT *, cpcl(j)%trueind
!                cpcl(i)%isunique = .FALSE.
!                EXIT
!              END IF
!            END DO
            cpcl(i)%trueind = truer
!            DO j = 1, ucptnum
!              IF (GetPointDistanceR(truer,cpl(j)%truer, &
!                  chg, nnlayers) <= 0.2 ) THEN
!                cpcl(i)%isunique = .FALSE.
!                EXIT
!              END IF
!            END DO
          END IF
          IF (cpcl(i)%isunique ) THEN
            ucptnum = ucptnum + 1
            CALL RecordCPR(truer,chg,cpl,ucptnum,eigvals,eigvecs, maxcount, &
              uringcount,ubondcount,ucagecount, opts, grad, interpolHessian)
!            PRINT *, 'ucptnum is', ucptnum
            CYCLE
          END IF
          truer = truer + nexttem ! this keeps track the total movement
          CALL pbc_r_lat(truer,chg%npts)
        ! moving on to the next critical pint candidate
        END DO
      END IF
      PRINT *, 'Number of critical point count: ', ucptnum
      PRINT *, 'Number of nucl critical point count: ', maxcount
      PRINT *, 'Number of bond critical point count: ', ubondcount
      PRINT *, 'Number of ring critical point count: ', uringcount
      PRINT *, 'Number of cage critical point count: ', ucagecount

      ! remove duplicate CPs
      isReduced = .FALSE.
      CALL ReduceCP(cpl,opts,ucptnum,chg,uBondCount, &
        uringcount,uCageCount,maxCount,isReduced)
      PRINT *, 'After a round of reduction'
      PRINT *, 'Numbver of atoms: ', ions%nions
      PRINT *, 'Number of critical point count: ', ucptnum
      PRINT *, 'Number of nucl critical point count: ', maxcount
      PRINT *, 'Number of bond critical point count: ', ubondcount
      PRINT *, 'Number of ring critical point count: ', uringcount
      PRINT *, 'Number of cage critical point count: ', ucagecount

      CALL  PHRuleExam(maxCount,ubondCount,uringCount,ucageCount,opts,&
        ions,phmrCompliant)
      ! output the cp to files
      CALL outputCP(cpl,opts,ucptnum,chg,setcount, ubondcount, &
        uringcount, ucagecount, maxcount)
!      IF (opts%ismolecule) THEN
!        IF (maxcount - bondcount + ringcount - cagecount == 1) THEN
!          PRINT *, 'Satisfies The Poincare-Hopf rule for a molecule'
!!          WRITE(98,*) 'Satisfies The Poincare-Hopf rule for a molecule'
!        ELSE
!          PRINT *, 'WARNING 3: Fails The poincare-Hopf rule for a molecule'
!!          WRITE(98,*)  'WARNING 3: Fails The poincare-Hopf rule for a molecule'
!        END IF
!      ELSE IF (opts%iscrystal) THEN
!        IF (maxcount - bondcount + ringcount - cagecount == 0) THEN
!          PRINT *, 'Satisfies The Poincare-Hopf rule for a crystal'
!!          WRITE(98,*) 'Satisfies The Poincare-Hopf rule for a crystal'
!        ELSE
!          PRINT *, 'WARNING 4: Fails The poincare-Hopf rule for a crystal'
!!          WRITE(98,*)  'WARNING 4: Fails The poincare-Hopf rule for a crystal'
!
!        END IF
!      END IF
      ! ending noInterpolation if
!      DO i = 1, cptnum
!        PRINT *, cpRoster(i,:)
!      END DO
      DEALLOCATE(cpRoster)
    END IF
    PRINT *, 'outputting debugging information to allcpPOSCAR'
!    CALL VisAllCP(cpcl,cptnum,chg,ions,opts)
    CALL VisAllCP(cpl,ucptnum,chg,ions,opts,uringCount,ubondCount,ucageCount)
    DEALLOCATE (cpl)
    DEALLOCATE (cpcl)
!    CLOSE(97)
!    CLOSE(98)
    CLOSE(1)
    CLOSE(2)

    END SUBROUTINE critpoint_find

    ! get cartesian coordinates of a point
    FUNCTION getcart(ind,lat2car)
      !ind is indicies of the current point
      !cart is the cartisian coordinates of the current point
      INTEGER,DIMENSION(3),INTENT(IN) :: ind
      REAL(q2),DIMENSION(3) :: getcart
      REAL(q2),DIMENSION(3,3) :: lat2car
      getcart(1) = ind(1) * lat2car(1,1) + & 
        ind(2) * lat2car(1,2) + & 
        ind(3) *  lat2car(1,3)
      getcart(2) = ind(1) * lat2car(2,1) + &
        ind(2) * lat2car(2,2) + &
        ind(3) *  lat2car(2,3)
      getcart(3) = ind(1) * lat2car(3,1) + &
        ind(2) * lat2car(3,2) + &
        ind(3) *  lat2car(3,3)

      RETURN
    END FUNCTION
    
    ! get indices from cartesian coordinates.
    FUNCTION  getinds(car2lat,intcarts) 
      REAL(q2),DIMENSION(6,3) :: getinds, intcarts
      REAL(q2),DIMENSION(3,3) :: car2lat,W
      REAL(q2),DIMENSION(3) :: Z
      INTEGER :: i,j,k
      DO i = 1, 6 
        getinds(i,:) = MATMUL(car2lat,intcarts(i,:))
      END DO
      RETURN
    END FUNCTION

    ! this function determins when looking for nn, how many layers to search
    ! within. It looks for the smallest vector sum of lattice vectors, and the
    ! largest vector
    FUNCTION findnnlayers(ions)
      INTEGER :: findnnlayers
      TYPE(ions_obj) :: ions
      REAL(q2), DIMENSION(4,3) :: latsums!lat12, lat13, lat23, lat123
      REAL(q2) :: latmag, minmag, maxmag
      INTEGER :: i
      latsums(1,:) = ions%lattice(1,:) + ions%lattice(2,:)
      latsums(2,:) = ions%lattice(1,:) + ions%lattice(3,:)
      latsums(3,:) = ions%lattice(2,:) + ions%lattice(3,:)
      latsums(4,:) = ions%lattice(1,:) + ions%lattice(2,:) + ions%lattice(3,:)
      minmag = mag(latsums(1,:))
      maxmag = mag(latsums(1,:))
      DO i = 2, 4
        latmag = mag(latsums(i,:))
        IF (latmag >= maxmag) THEN
          maxmag = latmag
        END IF
        IF (latmag <= minmag) THEN
          minmag = latmag
        END IF
      END DO
      DO i = 1, 3
        latmag = mag(ions%lattice(i,:))
        IF (latmag >= maxmag) THEN
          maxmag = latmag
        END IF
        IF (latmag <= minmag) THEN
          minmag = latmag
        END IF
      END DO
      findnnlayers = CEILING(maxmag/minmag)

    END FUNCTION findnnlayers
  
    ! This function gives the simple box for trilinear interpolation.
    FUNCTION SimpleNN(p,chg)
      REAL(q2), DIMENSION(3) :: p
      TYPE(charge_obj) :: chg
      INTEGER, DIMENSION(8,3) :: SimpleNN
      INTEGER :: i
      SimpleNN(1,:) = (/FLOOR(p(1)),FLOOR(p(2)),FLOOR(p(3))/)
      SimpleNN(2,:) = (/CEILING(p(1)),FLOOR(p(2)),FLOOR(p(3))/)
      SimpleNN(3,:) = (/FLOOR(p(1)),CEILING(p(2)),FLOOR(p(3))/)
      SimpleNN(4,:) = (/CEILING(p(1)),CEILING(p(2)),FLOOR(p(3))/)
      SimpleNN(5,:) = (/floor(p(1)),FLOOR(p(2)),CEILING(p(3))/)
      SimpleNN(6,:) = (/ceiling(p(1)),FLOOR(p(2)),CEILING(p(3))/)
      SimpleNN(7,:) = (/floor(p(1)),CEILING(p(2)),CEILING(p(3))/)
      SimpleNN(8,:) = (/ceiling(p(1)),CEILING(p(2)),CEILING(p(3))/)
!      DO i = 1, 8
!        CALL pbc(SimpleNN(i,:),chg%npts)
!      END DO
    END FUNCTION


    ! find neares on grid points for a interpolated point
    ! to do trilinear interpolation
    ! p is the off grid point we want
    FUNCTION FindNN(truer,nnlayers,chg,ions) ! THIS FUNCTION IS NOT STABLE
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nnind
      !grid points
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      TYPE(weight_obj), ALLOCATABLE, DIMENSION(:) :: nndist 
      ! instead of storing weight, this is used to
      ! store distances so that they can be sorted
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: FindNN,tfindnn
      REAL(q2),DIMENSION(3) :: intcart,pcart
      INTEGER,DIMENSION(3) :: p,maxs
      REAL(q2), DIMENSION(3) :: truer, truecart
      INTEGER,DIMENSION(3) :: baseind
      INTEGER, ALLOCATABLE, DIMENSION(:) :: rank
      REAL(q2),DIMENSION(27) :: dist
      INTEGER,DIMENSION(27) :: scores
      INTEGER,DIMENSION(27,3) :: p2
      INTEGER :: i,j,k,counter, nnlayers
      INTEGER, DIMENSION(8,3) :: corners, tempcorners
      ! since only closest neigherbos of current grid point
      ! will be used for interpolation
      counter = 0
      nnlayers = FindNNLayers(ions)
      ALLOCATE(nndist((nnlayers*2 + 1)**3 ))
      ALLOCATE(nnind((nnlayers*2 + 1)**3 ,3))
      ALLOCATE(rank((nnlayers*2 + 1)**3 ))
      ALLOCATE(FindNN((nnlayers*2+1)**3,3))
      !to calculate distance it is not necessary to run pbc
      !infact pbc should be avoided at this stage
      baseind = (/FLOOR(truer(1)),FLOOR(truer(2)),FLOOR(truer(3))/)
      truecart = MATMUL(chg%lat2car,truer)
      DO i = -nnlayers, nnlayers
        DO j = -nnlayers, nnlayers
          DO k = -nnlayers, nnlayers
            !IF (i == 0 .AND. j == 0 .AND. k == 0) THEN
            !  CYCLE
            !END IF
            counter = counter + 1
            nnind(counter,:) = (/i,j,k/) + baseind
            pcart = MATMUL(chg%lat2car,nnind(counter,:))
            nndist(counter)%rho = mag(truecart - pcart)
            nndist(counter)%pos(:) = (/i,j,k/) + baseind
            ! remember, rho here is really the distance
            FindNN(counter,:) = (/i,j,k/) + baseind
          END DO
        END DO
      END DO
      !PRINT *, 'nnind are'
      !DO i = 1,  ((nnlayers * 2 + 1)**3 - 1)
      !  PRINT *, 'i = ', i
      !  PRINT *, nnind(i,:)
      !END DO
      ! now sort the list with distance high to low
!      CALL sort_weight((nnlayers*2+1)**3,nndist)
!      ! Keep the closest 8 points
!      DO i = 1,8
!        FindNN(i,:) = nndist((nnlayers*2+1)**3 +1 - i)%pos
!        CALL pbc(FindNN(i,:),chg%npts)
!      END DO
      ! Mission accomplished.
      ! Trying to organize these nn into a box is job for another function.

      !PRINT *, 'these are the places recorded'
      !DO i = 1,8
      !  PRINT *, findnn(i,:)
      !END DO 
      !PRINT *, 'these are all places looked at'
      !DO i = 1, (nnlayers*2+1)**3
      !  PRINT *, nndist(i)%pos
      !END DO
      !PRINT *, 'temcorners are'
      !DO i = 1,8
      !  tempcorners(i,:) = nndist( (nnlayers*2 + 1 )**3 -i)%pos
      !  PRINT *, tempcorners(i,:)
      !END DO
      !
      ! then these needs to be rearranged 
      ! so that is follows this format:
      !  1   2   3   4   5   6   7   8
      ! 000 001 010 100 011 101 110 111
      ! indices should have only 2 values
      ! on each axis. min gives 0 on that 
      ! axis and max gives 1. 
      !maxs = (/0,0,0/)
      !DO i = 1, 8
      !  DO j = 1,3
      !    IF (tfindnn(i,j)>= maxs(j) ) THEN
      !      maxs(j) = tfindnn(i,j)
      !    END IF
      !  END DO
      !END DO
      !PRINT *, 'maxs is', maxs
      !DO i = 1, 8
      !  PRINT *, 'tfindnn here is'
      !  PRINT * ,tfindnn(i,:)
      !  IF (ALL(tfindnn(i,:) == maxs)) THEN
      !    PRINT *, '8'
      !    findnn(8,:) = tfindnn(i,:)
      !    CYCLE
      !  END IF
      !  IF (ALL(tfindnn(i,:) == maxs - (/1,1,1/))) THEN
      !    findnn(1,:) = tfindnn(i,:)
      !    PRINT *, '1'
      !    CYCLE
      !  END IF
      !  IF (ALL(tfindnn(i,:) == maxs - (/1,1,0/))) THEN
      !    findnn(2,:) = tfindnn(i,:)
      !    PRINT *, '2'
      !    CYCLE
      !  END IF
      !  IF (ALL(tfindnn(i,:) == maxs - (/1,0,1/))) THEN
      !    findnn(3,:) = tfindnn(i,:)
      !    PRINT *, '3'
      !    CYCLE
      !  END IF
      !  IF (ALL(tfindnn(i,:) == maxs - (/0,1,1/))) THEN
      !    findnn(4,:) = tfindnn(i,:)
      !    PRINT *, '4'
      !    CYCLE
      !  END IF
      !  IF (ALL(tfindnn(i,:) == maxs - (/1,0,0/))) THEN
      !    findnn(5,:) = tfindnn(i,:)
      !    PRINT *, '5'
      !    CYCLE
      !  END IF
      !  IF (ALL(tfindnn(i,:) == maxs - (/0,1,0/))) THEN
      !    findnn(6,:) = tfindnn(i,:)
      !    PRINT *, '6'
      !    CYCLE
      !  END IF
      !  IF (ALL(tfindnn(i,:) == maxs - (/0,0,1/))) THEN
      !    findnn(7,:) = tfindnn(i,:)
      !    PRINT *, '7'
      !    CYCLE
      !  END IF
      !END DO 
      DEALLOCATE(nnind)
      DEALLOCATE(nndist)
      RETURN
    END FUNCTION FindNN

    ! This function should take in a list of nearest neighbors predetermined.
    FUNCTION nn_grad(chg,r,rho,nn)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),INTENT(OUT) :: rho
      REAL(q2),DIMENSION(3) :: nn_grad
      INTEGER :: p1, p2, p3
      REAL(q2),DIMENSION(3) :: rho_grad_lat
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2) :: rho00_, rho01_, rho10_, rho11_
      REAL(q2) :: rho0__, rho1__, rho_0_, rho_1_, rho__0, rho__1
      REAL(q2) :: rho_00, rho_01, rho_10, rho_11
      INTEGER,DIMENSION(8,3) :: nn
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      rho000 = rho_val(chg,nn(1,1),nn(1,2),nn(1,3))
      rho001 = rho_val(chg,nn(2,1),nn(2,2),nn(2,3))
      rho010 = rho_val(chg,nn(3,1),nn(3,2),nn(3,3))
      rho100 = rho_val(chg,nn(4,1),nn(4,2),nn(4,3))
      rho011 = rho_val(chg,nn(5,1),nn(5,2),nn(5,3))
      rho101 = rho_val(chg,nn(6,1),nn(6,2),nn(6,3))
      rho110 = rho_val(chg,nn(7,1),nn(7,2),nn(7,3))
      rho111 = rho_val(chg,nn(8,1),nn(8,2),nn(8,3))
      rho00_ = rho000*g3 + rho001*f3
      rho01_ = rho010*g3 + rho011*f3
      rho10_ = rho100*g3 + rho101*f3
      rho11_ = rho110*g3 + rho111*f3
      rho0__ = rho00_*g2 + rho01_*f2
      rho1__ = rho10_*g2 + rho11_*f2
      rho = rho0__*g1 + rho1__*f1
  ! More work for gradients
      rho_0_ = rho00_*g1 + rho10_*f1
      rho_1_ = rho01_*g1 + rho11_*f1
      rho_00 = rho000*g1 + rho100*f1
      rho_01 = rho001*g1 + rho101*f1
      rho_10 = rho010*g1 + rho110*f1
      rho_11 = rho011*g1 + rho111*f1
      rho__0 = rho_00*g2 + rho_10*f2
      rho__1 = rho_01*g2 + rho_11*f2
      rho_grad_lat(1) = rho1__ - rho0__
      rho_grad_lat(2) = rho_1_ - rho_0_
      rho_grad_lat(3) = rho__1 - rho__0
  !   CALL vector_matrix(rho_grad_lat, chg%car2lat, rho_grad)
      nn_grad = MATMUL(chg%car2lat,rho_grad_lat)
    RETURN
    END FUNCTION nn_grad

    ! this funciton takes in 8 values, return a
    ! trilinear interpolated gradient of the values.
    ! the 8 value list value order is 
    ! 000 001 010 100 011 101 110 111
    ! Note 02042019: the above order is what I wrote previously
    ! Note 02042019: I believe the actuall order is the following
    ! 000 100 010 110 001 101 011 111
    !  1   2   3   4   5   6   7   8
    ! r is the indice of the predicted critical point
    ! The interpolation result is checked to be OK by mathematica
    FUNCTION trilinear_interpol_grad(vals,r)
      ! varls come nngrad
      ! could r be the problem?
      ! r needs to be a vector where each component 
      ! number is between 0 and 1
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),DIMENSION(3) :: trilinear_interpol_grad
      INTEGER :: p1, p2, p3
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2),DIMENSION(8,3) :: vals
      !REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2),DIMENSION(3) :: val00_, val01_, val10_, val11_
      REAL(q2),DIMENSION(3) :: val0__, val1__, val_0_, val_1_, val__0, val__1
      REAL(q2),DIMENSION(3) :: val_00, val_01, val_10, val_11
 
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      ! f1 f2 f3 are checked to be correct. 
      ! they should equal to tem for the first step
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      val00_ = vals(1,:)*g3 + vals(5,:)*f3
      val01_ = vals(3,:)*g3 + vals(7,:)*f3
      val10_ = vals(2,:)*g3 + vals(6,:)*f3
      val11_ = vals(4,:)*g3 + vals(8,:)*f3
      val0__ = val00_*g2 + val01_*f2
      val1__ = val10_*g2 + val11_*f2
      trilinear_interpol_grad = val0__*g1 + val1__*f1
    RETURN
    END FUNCTION trilinear_interpol_grad
 
    FUNCTION trilinear_interpol_rho(chg,r)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2) :: trilinear_interpol_rho
      INTEGER :: p1, p2, p3
      INTEGER, DIMENSION(3) :: tempp
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2),DIMENSION(8) :: vals
      !REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2) :: val00_, val01_, val10_, val11_
      REAL(q2) :: val0__, val1__, val_0_, val_1_, val__0, val__1
      REAL(q2) :: val_00, val_01, val_10, val_11
 
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      tempp = (/p1,p2,p3/)
      CALL pbc(tempp,chg%npts)
      vals(1) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1+1,p2,p3/)
      CALL pbc(tempp,chg%npts)
      vals(2) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1,p2+1,p3/)
      CALL pbc(tempp,chg%npts)
      vals(3) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1+1,p2+1,p3/)
      CALL pbc(tempp,chg%npts)
      vals(4) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1,p2,p3+1/)
      CALL pbc(tempp,chg%npts)
      vals(5) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1+1,p2,p3+1/)
      CALL pbc(tempp,chg%npts)
      vals(6) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1,p2+1,p3+1/)
      CALL pbc(tempp,chg%npts)
      vals(7) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      tempp = (/p1+1,p2+1,p3+1/)
      CALL pbc(tempp,chg%npts)
      vals(8) = rho_val(chg,tempp(1),tempp(2),tempp(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      ! f1 f2 f3 are checked to be correct. 
      ! they should equal to tem for the first step
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      val00_ = vals(1)*g3 + vals(5)*f3
      val01_ = vals(3)*g3 + vals(7)*f3
      val10_ = vals(2)*g3 + vals(6)*f3
      val11_ = vals(4)*g3 + vals(8)*f3
      val0__ = val00_*g2 + val01_*f2
      val1__ = val10_*g2 + val11_*f2
      trilinear_interpol_rho = val0__*g1 + val1__*f1
    RETURN
    END FUNCTION

    FUNCTION trilinear_interpol_hes(vals,r)
      ! varls come nnhes
      ! could r be the problem?
      ! r needs to be a vector where each component
      ! number is between 0 and 1
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),DIMENSION(3,3) :: trilinear_interpol_hes
      INTEGER :: p1, p2, p3
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2),DIMENSION(8,3,3) :: vals
      !REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110,
      !rho111
      REAL(q2),DIMENSION(3,3) :: val00_, val01_, val10_, val11_
      REAL(q2),DIMENSION(3,3) :: val0__, val1__, val_0_, val_1_, val__0, val__1
      REAL(q2),DIMENSION(3,3) :: val_00, val_01, val_10, val_11

      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      ! f1 f2 f3 are checked to be correct.
      ! they should equal to tem for the first step
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      val00_(1,:) = vals(1,1,:)*g3 + vals(5,1,:)*f3
      val00_(2,:) = vals(1,2,:)*g3 + vals(5,2,:)*f3
      val00_(3,:) = vals(1,3,:)*g3 + vals(5,3,:)*f3
      val01_(1,:) = vals(3,1,:)*g3 + vals(7,1,:)*f3
      val01_(2,:) = vals(3,2,:)*g3 + vals(7,2,:)*f3
      val01_(3,:) = vals(3,3,:)*g3 + vals(7,3,:)*f3
      val10_(1,:) = vals(2,1,:)*g3 + vals(6,1,:)*f3
      val10_(2,:) = vals(2,2,:)*g3 + vals(6,2,:)*f3
      val10_(3,:) = vals(2,3,:)*g3 + vals(6,3,:)*f3
      val11_(1,:) = vals(4,1,:)*g3 + vals(8,1,:)*f3
      val11_(2,:) = vals(4,2,:)*g3 + vals(8,2,:)*f3
      val11_(3,:) = vals(4,3,:)*g3 + vals(8,3,:)*f3
      val0__ = val00_*g2 + val01_*f2
      val1__ = val10_*g2 + val11_*f2
      trilinear_interpol_hes = val0__*g1 + val1__*f1
      
    END FUNCTION


    FUNCTION coorcheck(lattice)
      LOGICAL,DIMENSION(3) :: coorcheck
      REAL(q2) :: prod1,prod2,prod3
      REAL(q2),DIMENSION(3,3) :: lattice
      INTEGER :: i
      ! elements 1, 2, 3 are for whether coordinates align with x, y, z axis.
      coorcheck = (/.TRUE.,.TRUE.,.TRUE./)
      IF (lattice(1,2) .NE. 0 .OR. lattice(1,3) .NE. 0) THEN
        coorcheck(1)=.FALSE.
      END IF 
      IF (lattice(2,1) .NE. 0 .OR. lattice(2,3) .NE. 0) THEN
        coorcheck(2)=.FALSE.
      END IF
      IF (lattice(3,2) .NE. 0 .OR. lattice(3,1) .NE. 0) THEN
        coorcheck(3)=.FALSE.
      END IF
    END FUNCTION
    
    ! The following subroutine gets gradient using central difference
    FUNCTION CDGrad(p,chg)
      TYPE(charge_obj) :: chg
      INTEGER, DIMENSION(3) :: p
      INTEGER, DIMENSION(3) :: pzm,pzp,pxm,pxp,pym,pyp
      REAL(q2), DIMENSION(3) :: CDGrad
      pzm = p + (/0,0,-1/)
      pzp = p + (/0,0,1/)
      pxm = p + (/-1,0,0/)
      pxp = p + (/1,0,0/)
      pym = p + (/0,-1,0/)
      pyp = p + (/0,1,0/)
      CALL pbc(pxm,chg%npts)
      CALL pbc(pym,chg%npts)
      CALL pbc(pzm,chg%npts)
      CALL pbc(pxp,chg%npts)
      CALL pbc(pyp,chg%npts)
      CALL pbc(pzp,chg%npts)
      CDGrad(3) = 0.5*(rho_val(chg,pzp(1),pzp(2),pzp(3)) - &
                  rho_val(chg,pzm(1),pzm(2),pzm(3)))
      CDGrad(2) = 0.5*(rho_val(chg,pyp(1),pyp(2),pyp(3)) - &
                  rho_val(chg,pym(1),pym(2),pym(3)))
      CDGrad(1) = 0.5*(rho_val(chg,pxp(1),pxp(2),pxp(3)) - &
                  rho_val(chg,pxm(1),pxm(2),pxm(3)))
      CDGrad = MATMUL(chg%car2lat,CDGrad)
      RETURN
      ! now the gradient should be in cartesian
    END FUNCTION
    
    ! the following subroutine gets hes and force in lattice units
    FUNCTION CDHessian(p,chg)
      REAL(q2),DIMENSION(3,3) :: CDHessian
      TYPE(charge_obj) :: chg
      INTEGER,DIMENSION(3) :: ptxy1, ptxy2, ptxy3, ptxy4, ptxz1, ptxz2, ptxz3 
      INTEGER,DIMENSION(3) :: ptxz4, ptyz1, ptyz2, ptyz3, ptyz4
      INTEGER,DIMENSION(3) :: p, ptx1, ptx2, pty1, pty2, ptz1, ptz2
      CALL pbc(p,chg%npts)
      ptx1 = p + (/1,0,0/)
      ptx2 = p + (/-1,0,0/)
      pty1 = p + (/0,1,0/)
      pty2 = p + (/0,-1,0/)
      ptz1 = p + (/0,0,1/)
      ptz2 = p + (/0,0,-1/)
      CALL pbc(ptx1,chg%npts)
      CALL pbc(ptx2,chg%npts)
      CALL pbc(pty1,chg%npts)
      CALL pbc(pty2,chg%npts)
      CALL pbc(ptz1,chg%npts)
      CALL pbc(ptz2,chg%npts)
      CDHessian(1,1) = &
         (rho_val(chg,ptx1(1),ptx1(2),ptx1(3)) - &
         rho_val(chg,p(1),p(2),p(3))) - &
         (rho_val(chg,p(1),p(2),p(3)) - &
         rho_val(chg,ptx2(1),ptx2(2),ptx2(3)))
      CDHessian(2,2) = &
        (rho_val(chg,pty1(1),pty1(2),pty1(3)) - &
         rho_val(chg,p(1),p(2),p(3))) - &
        (rho_val(chg,p(1),p(2),p(3)) - &
         rho_val(chg,pty2(1),pty2(2),pty2(3)))
      CDHessian(3,3) = &
        (rho_val(chg,ptz1(1),ptz1(2),ptz1(3)) - &
         rho_val(chg,p(1),p(2),p(3))) - &
        (rho_val(chg,p(1),p(2),p(3)) - &
         rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))
      ptxy1 = p + (/-1,-1,0/)
      ptxy2 = p + (/-1,+1,0/)
      ptxy3 = p + (/+1,+1,0/)
      ptxy4 = p + (/+1,-1,0/)
      ptxz1 = p + (/-1,0,-1/)
      ptxz2 = p + (/-1,0,+1/)
      ptxz3 = p + (/+1,0,+1/)
      ptxz4 = p + (/+1,0,-1/)
      ptyz1 = p + (/0,-1,-1/)
      ptyz2 = p + (/0,-1,+1/)
      ptyz3 = p + (/0,+1,+1/)
      ptyz4 = p + (/0,+1,-1/)
      CALL pbc(ptxy1,chg%npts)
      CALL pbc(ptxy2,chg%npts)
      CALL pbc(ptxy3,chg%npts)
      CALL pbc(ptxy4,chg%npts)
      CALL pbc(ptxz1,chg%npts)
      CALL pbc(ptxz2,chg%npts)
      CALL pbc(ptxz3,chg%npts)
      CALL pbc(ptxz4,chg%npts)
      CALL pbc(ptyz1,chg%npts)
      CALL pbc(ptyz2,chg%npts)
      CALL pbc(ptyz3,chg%npts)
      CALL pbc(ptyz4,chg%npts)
      CDHessian(1,2) = &
        ( &
        ! this is the backward dv
!        - 0.5_q2 / vmag * &
        - 0.25_q2 * & 
        ((rho_val(chg,ptxy2(1),ptxy2(2),ptxy2(3)) + &
            rho_val(chg,pty1(1),pty1(2),pty1(3)))  - &
          (rho_val(chg,ptxy1(1),ptxy1(2),ptxy1(3)) + &
            rho_val(chg,pty2(1),pty2(2),pty2(3)))  ) & 
        ! this is the forward dv
        + 0.25_q2 * & 
        ((rho_val(chg,ptxy3(1),ptxy3(2),ptxy3(3)) + &
            rho_val(chg,pty1(1),pty1(2),pty1(3)))  - &
          (rho_val(chg,ptxy4(1),ptxy4(2),ptxy4(3)) + &
            rho_val(chg,pty2(1),pty2(2),pty2(3)))  ) &
        )
      CDHessian(2,1) = CDHessian(1,2)
      CDHessian(1,3) = &
        ( &
        ! this is the bacward dw
        - 0.25_q2 * &
        ((rho_val(chg,ptxz2(1),ptxz2(2),ptxz2(3)) + &
            rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))   - &
          (rho_val(chg,ptxz1(1),ptxz1(2),ptxz1(3)) + &
            rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))   )  &
        ! this is the forward dw
        + 0.25_q2 * & 
        ((rho_val(chg,ptxz3(1),ptxz3(2),ptxz3(3)) + &
            rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))   - &
          (rho_val(chg,ptxz4(1),ptxz4(2),ptxz4(3)) + &
            rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))   ) &
        ) 
      CDHessian(3,1) = CDHessian(1,3)
      CDHessian(2,3) = &
        ( &
        ! this is the bacward dw
        - 0.25_q2 * &
        ((rho_val(chg,ptyz2(1),ptyz2(2),ptyz2(3)) + &
            rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))  - &
          (rho_val(chg,ptyz1(1),ptyz1(2),ptyz1(3)) + &
            rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))  ) & 
        ! this is the forward dw
        + 0.25_q2 * & 
        ((rho_val(chg,ptyz3(1),ptyz3(2),ptyz3(3)) + &
            rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))  - &
          (rho_val(chg,ptyz4(1),ptyz4(2),ptyz4(3)) + &
            rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))  ) &
        )
      CDHessian(3,2) = CDHessian(2,3)
      ! Convert the hessian, which is now in lattice coordinates, to cartesian
      !CDHessian = MATMUL(MATMUL(chg%car2lat,CDHessian),TRANSPOSE(chg%car2lat))
      CDHessian = MATMUL(TRANSPOSE(chg%car2lat),MATMUL(CDHessian,chg%car2lat))
      RETURN
    END FUNCTION
    

    ! This function takes in current position and grid point position in
    ! lattice, finds the distance this point is to the nearest cell, take half
    ! the distance as step size. 
    FUNCTION findstepsize(r)
      REAL(q2) :: findstepsize    
      REAL(q2),DIMENSION(3) :: r
      REAL(q2) :: f1,f2,f3,c1,c2,c3
      
      f1 = MIN(ABS(r(1)),ABS(0.5-ABS(r(1)))) 
      f2 = MIN(ABS(r(2)),ABS(0.5-ABS(r(2))))
      f3 = MIN(ABS(r(3)),ABS(0.5-ABS(r(3))))
      findstepsize = MIN(MIN(f1,f2),MIN(f1,f3))/2
    RETURN
    END FUNCTION findstepsize

    ! this funciton finds hessian of a interpolated point using interpolated
    ! nearest neighbor gradiants. Also gradiants taken in here should be in
    ! lattice.
    FUNCTION inthessian(grad,stepsize)
      REAL(q2),DIMENSION(6,3) :: grad
      REAL(q2),DIMENSION(3,3) :: inthessian
      REAL(q2) :: stepsize
      ! again, intnngrad is following this order:
      ! +x -x +y -y +z -z
      inthessian(1,1) = (grad(1,1)-grad(2,1))*0.5_q2/stepsize
      inthessian(2,2) = (grad(3,2)-grad(4,2))*0.5_q2/stepsize
      inthessian(3,3) = (grad(5,3)-grad(6,3))*0.5_q2/stepsize
      inthessian(1,2) = (grad(3,1)-grad(4,1))*0.5_q2/stepsize
      inthessian(2,1) = inthessian(1,2)
      inthessian(1,3) = (grad(5,1)-grad(6,1))*0.5_q2/stepsize
      inthessian(3,1) = inthessian(1,3)
      inthessian(2,3) = (grad(5,2)-grad(6,2))*0.5_q2/stepsize
      inthessian(3,2) = inthessian(2,3)
    ! assuming that this function is fine
    RETURN
    END FUNCTION inthessian

   
    ! below is the function for least sqare gradient
    FUNCTION lsg(r0,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct) 
      TYPE(charge_obj) :: chg
      REAL(q2), DIMENSION(3) :: lsg
      !INTEGER, DIMENSION(3) :: p
      ! input may either be integers or not. May need to change later on.
      INTEGER, DIMENSION(3) :: r0
      ! r0 is the position of the current grid point
      !INTEGER, DIMENSION(26,3) :: vi , nbp
      INTEGER,DIMENSION(3,26) :: vi, nbp
      ! v is a column of coordinates. 
      ! nbp is indecies of neighbors, to be used for pbc
      INTEGER, DIMENSION(26,3) :: vit
      ! vit is the transpose of vi
      ! vi are the vectors from point r0 to ri, where ri are all neighbors of r0
      REAL(q2), DIMENSION(3) :: frakturavi
      REAL(q2), DIMENSION(26) :: w, di 
      REAL(q2), DIMENSION(3,26) :: matw
      !REAL(q2), DIMENSION(26,3) :: matw
      REAL(q2), DIMENSION(3,13) :: matwprime
      !REAL(q2), DIMENSION(13,3) :: matwprime
      REAL(q2), DIMENSION(3,3) ::  matm
      REAL(q2), DIMENSION(26) :: wi
      ! wi are the weights of each neighbor
      REAL(q2), DIMENSION(3,26) :: bi
      ! bi is weight times direction times difference towards each neighbor
      REAL(q2), DIMENSION(3) :: vecb
      REAL(q2), DIMENSION(13) :: deltarho
      ! deltarho is the difference between +ri and -ri
      REAL(q2), DIMENSION(3,3) :: ggrid
      ! ggrid is the metric tensor, from grid basis to dual basis
      INTEGER :: i, j, k
      REAL(q2), DIMENSION(3,3) :: testermatrix
      REAL(q2), DIMENSION(3,3) :: outerproduct
      ! get the value differences to all neighbors
      PRINT *, 'position is'
      PRINT *, r0
      PRINT *, 'in cartesian is'
      PRINT *, MATMUL(chg%lat2car,r0)
      CALL pbc(r0,chg%npts)
      DO i = 1, 26
        nbp(:,i) = vi(:,i) + r0
        CALL pbc(nbp(:,i),chg%npts)
!        PRINT *, nbp(:,i)
!        PRINT *, MATMUL(chg%lat2car,nbp(:,i))
      END DO
      PRINT *, 'deltarho is'
      DO i = 1 , 13
        PRINT *, 'nbp i is'
        PRINT *, nbp(:,i)
        PRINT *, MATMUL(chg%lat2car,nbp(:,i))
        PRINT *, 'nbp i + 13 is'
        PRINT *, nbp(:,i + 13)
        PRINT *, MATMUL(chg%lat2car,nbp(:,i+13))
        deltarho(i) = rho_val(chg,nbp(1,i),nbp(2,i),nbp(3,i)) - &
          rho_val(chg,nbp(1,i+13),nbp(2,i+13),nbp(3,i+13))
        PRINT *, rho_val(chg,nbp(1,i),nbp(2,i),nbp(3,i)), &
          rho_val(chg,nbp(1,i+13),nbp(2,i+13),nbp(3,i+13))
        PRINT *, deltarho(i)
      END DO
      lsg = 0._q2
      PRINT *, 'ggrid is'
      PRINT *, ggrid(1,:)
      PRINT *, ggrid(2,:)
      PRINT *, ggrid(3,:)
      PRINT *, 'matm is'
      DO i = 1, 3
        PRINT *, matm(i,:)
      END DO
      PRINT *, 'matwprime is'
      DO i = 1,13
        PRINT *, matwprime(:,i)
      END DO
      PRINT *, 'product without delta rho is'
      PRINT *, MATMUL( &
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime)
      PRINT *, 'inverse ggrid mul inverse matm is'
      PRINT *, MATMUL(INVERSE(ggrid),INVERSE(matm))
      ! The following code is for debugging, which calculates gradient in direct
      ! coordinates.
      lsg = MATMUL(MATMUL(INVERSE(matm),matwprime),deltarho)
      PRINT *, 'gradient in direct is'
      PRINT *, lsg 
      lsg = &
        MATMUL( &
          MATMUL( & 
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime) & 
        , deltarho)
      PRINT *, 'calculated lsg is'
      PRINT *, lsg
      ! lsg up till now seems to be larger than force in cartesian units by car2lat
      lsg = MATMUL(chg%lat2car,lsg) ! now force should be in cartesian units
      PRINT *, 'multiplied with lat2car again is'
      PRINT *, lsg
      RETURN
    END FUNCTION lsg

    ! This function finds hessian by 
    ! finding lsg of gradients found with lsg
    FUNCTION lsh(r0,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      TYPE(charge_obj) :: chg
      TYPE(hessian) :: hes
      REAL(q2), DIMENSION(3,3) :: lsh
      INTEGER, DIMENSION(3) :: r0
      ! to store neighbor lsg's
      REAL(q2), DIMENSION(3,26) :: lsg_val
      ! vi, as usual
      INTEGER, DIMENSION(3,26) :: vi
      INTEGER, DIMENSION(26,3) :: vit
      ! ggrid, as usual
      REAL(q2), DIMENSION(3,3) :: ggrid
      ! one time use for neighbor positions
      INTEGER, DIMENSION(3) :: nbp
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
 
      ! differences in gradient components
      REAL(q2), DIMENSION(13) :: deltadx, deltady, deltadz
      REAL(q2), DIMENSION(3,3) :: matm
      REAL(q2), DIMENSION(3,3) :: outerproduct
      REAL(q2) :: average
      INTEGER :: i, j, k
      DO i = 1, 26
        nbp = vi(:,i) + r0
        CALL pbc(nbp,chg%npts) ! nbps are fine
        lsg_val(:,i) = lsg(nbp,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      END DO
      DO i = 1, 13
        ! matwprime matches with lsg
        deltadx(i) = lsg_val(1,i) - lsg_val(1,i+13)
        deltady(i) = lsg_val(2,i) - lsg_val(2,i+13)
        deltadz(i) = lsg_val(3,i) - lsg_val(3,i+13)
      END DO
      lsh(1,:) = &
        MATMUL( &
          MATMUL( &
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime &  
          ), &        
        deltadx &
        )
      lsh(2,:) = &
        MATMUL( &
          MATMUL( &
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime &  
          ), &        
        deltady &
        )
      lsh(3,:) = &
        MATMUL( &
          MATMUL( &
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime &  
          ), &        
        deltadz &
        )

      ! the hessian up till now seems larger than cartesian by car2lat
      lsh(1,:) = MATMUL(chg%lat2car,lsh(1,:))
      lsh(2,:) = MATMUL(chg%lat2car,lsh(2,:))
      lsh(3,:) = MATMUL(chg%lat2car,lsh(3,:))
      average = (lsh(1,2) + lsh(2,1))/2
      lsh(1,2) = average
      lsh(2,1) = average
      average = (lsh(1,3) + lsh(3,1))/2
      lsh(1,3) = average
      lsh(3,1) = average
      average = (lsh(2,3) + lsh(3,2)) /2
      lsh(2,3) = average
      lsh(3,2) = average
      RETURN
    END FUNCTION lsh

    FUNCTION makevi()
      INTEGER, DIMENSION(3,26) :: makevi
      makevi(:,1)=(/-1,-1,-1/)
      makevi(:,2)=(/-1,-1,0/)
      makevi(:,3)=(/-1,-1,1/)
      makevi(:,4)=(/-1,0,-1/)
      makevi(:,5)=(/-1,0,0/)
      makevi(:,6)=(/-1,0,1/)
      makevi(:,7)=(/-1,1,-1/)
      makevi(:,8)=(/-1,1,0/)
      makevi(:,9)=(/-1,1,1/)
      makevi(:,10)=(/0,-1,-1/)
      makevi(:,11)=(/0,-1,0/)
      makevi(:,12)=(/0,-1,1/)
      makevi(:,13)=(/0,0,-1/)
      makevi(:,14)=(/1,1,1/)
      makevi(:,15)=(/1,1,0/)
      makevi(:,16)=(/1,1,-1/)
      makevi(:,17)=(/1,0,1/)
      makevi(:,18)=(/1,0,0/)
      makevi(:,19)=(/1,0,-1/)
      makevi(:,20)=(/1,-1,1/)
      makevi(:,21)=(/1,-1,0/)
      makevi(:,22)=(/1,-1,-1/)
      makevi(:,23)=(/0,1,1/)
      makevi(:,24)=(/0,1,0/)
      makevi(:,25)=(/0,1,-1/)
      makevi(:,26)=(/0,0,1/)
      RETURN
    END FUNCTION

    FUNCTION makeggrid(chg,ions)
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      REAL(q2),DIMENSION(3,3) :: makeggrid
      makeggrid(1,1) = DOT_PRODUCT(ions%lattice(1,:),ions%lattice(1,:))/ &
        (chg%npts(1)*chg%npts(1))
      makeggrid(1,2) = DOT_PRODUCT(ions%lattice(1,:),ions%lattice(2,:))/ &
        (chg%npts(1)*chg%npts(2))
      makeggrid(1,3) = DOT_PRODUCT(ions%lattice(1,:),ions%lattice(3,:))/ &
        (chg%npts(1)*chg%npts(3))
      makeggrid(2,1) = DOT_PRODUCT(ions%lattice(2,:),ions%lattice(1,:))/ &
        (chg%npts(2)*chg%npts(1))
      makeggrid(2,2) = DOT_PRODUCT(ions%lattice(2,:),ions%lattice(2,:))/ &
        (chg%npts(2)*chg%npts(2))
      makeggrid(2,3) = DOT_PRODUCT(ions%lattice(2,:),ions%lattice(3,:))/ &
        (chg%npts(2)*chg%npts(3))
      makeggrid(3,1) = DOT_PRODUCT(ions%lattice(3,:),ions%lattice(1,:))/ &
        (chg%npts(3)*chg%npts(1))
      makeggrid(3,2) = DOT_PRODUCT(ions%lattice(3,:),ions%lattice(2,:))/ &
        (chg%npts(3)*chg%npts(2))
      makeggrid(3,3) = DOT_PRODUCT(ions%lattice(3,:),ions%lattice(3,:))/ &
        (chg%npts(3)*chg%npts(3))
      RETURN
    END FUNCTION

    FUNCTION ascension(ind,chg,matm,matwprime,wi,vi,vit, & 
                       ggrid,outerproduct,opts,nnLayers,ions)
      ! this function finds nucleus critical points. 
      INTEGER, DIMENSION(3) :: ind
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      TYPE(ions_obj) :: ions
      REAL(q2),DIMENSION(3) :: ascension, distance
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: nnind
      REAL(q2), DIMENSION(8,3) :: nngrad
      INTEGER :: j, stepcount, nnLayers
      REAL(q2), DIMENSION(3) :: tempr, rn, rnm1 ! rn minus 1
      REAL(q2), DIMENSION(3) :: grad, stepsize, gradnm1
      INTEGER, DIMENSION(3,26) :: vi, matw
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(3,3) :: matm, outerproduct
      stepcount = 0
      ALLOCATE(nnInd(8,3))
!      ALLOCATE(nnInd((nnlayers*2+1)**3,3))
      ! initialize the process
      CALL pbc(ind,chg%npts)
      rn(1) = REAL(ind(1),q2)
      rn(2) = REAL(ind(2),q2)
      rn(3) = REAL(ind(3),q2)
      stepsize(1) = 0.5
      stepsize(2) = 0.5
      stepsize(3) = 0.5
      IF (opts%leastSquare_flag) THEN
        grad = lsg(ind,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      ELSE
        grad = cdgrad(ind,chg)
      END IF
      ! this gradient is in cartesian. convert it to lattice
      grad = MATMUL(chg%lat2car,grad)
      DO WHILE (stepsize(1) >= 0.01 .AND. &
                stepsize(2) >= 0.01 .AND. &
                stepsize(3) >= 0.01 )
        ! the gradient is in cartesian. 
        gradnm1 = grad
        rnm1 = rn
        ! determine where to go, a unit vector
        tempr = grad / SQRT(SUM(grad*grad))
        tempr(1) = stepsize(1) * tempr(1)
        tempr(2) = stepsize(2) * tempr(2)
        tempr(3) = stepsize(3) * tempr(3)
        rn = rn + tempr
        CALL pbc_r_lat(rn,chg%npts)
        nnind = simpleNN(rn,chg)
        DO j = 1,8
          CALL pbc(nnind(j,:),chg%npts)
          !IF (opts%leastsquare_flag .EQV. .TRUE.) THEN
          IF ( opts%leastSquare_flag ) THEN
            nngrad(j,:) = lsg(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
          ELSE 
            nngrad(j,:) = cdgrad(nnind(j,:),chg)
          END IF
        END DO
        distance = rn - nnind(1,:)
         !Row find the nearest neighbors at this new locaiton
         !First update critical point location
         !The next big step is to interpolate the force at predicted critical
         !point.
        grad = trilinear_interpol_grad(nngrad,distance) ! val r interpol
        !nnind = FindNN(rn,nnLayers,chg,ions)
        !grad = R2GradInterpol(nnind,rn,chg,nnLayers)
        ! this grad is in cartesian. convert it to lattice
        grad = MATMUL(chg%lat2car,grad)
        IF (ABS(grad(1)) < 0.001 &
            .AND. ABS(grad(2)) <= 0.001 &
            .AND. ABS(grad(3)) <= 0.001)  THEN
          ! we are at a critical point!
!          PRINT *, 'ascention gradient sufficiently small'
!          PRINT *, 'ascension rn is'
!          PRINT *, rn
          EXIT
        END IF
        ! if grad points backwards, reduce stepsize
        IF (DOT_PRODUCT(grad,gradnm1) <= 0) THEN
          stepsize = 0.5*stepsize
        END IF
        !IF (grad(1) * gradnm1(1) < 0) THEN
        !  stepsize(1) = stepsize(1)/2
        !END IF
        !IF (grad(2) * gradnm1(2) < 0) THEN
        !  stepsize(2) = stepsize(2)/2
        !END IF
        !IF (grad(3) * gradnm1(3) < 0) THEN
        !  stepsize(3) = stepsize(3)/2
        !END IF
        stepcount = stepcount + 1
!        PRINT *, rn
      END DO 
      ascension = rn
!      PRINT *, 'ascension step count is ', stepcount
      RETURN 
    END FUNCTION ascension
    
    FUNCTION descension(ind,chg,matm,matwprime,wi,vi,vit, & 
                       ggrid,outerproduct,opts,nnLayers,ions)
      ! this function finds nucleus critical points. 
      INTEGER, DIMENSION(3) :: ind
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      TYPE(ions_obj) :: ions
      REAL(q2),DIMENSION(3) :: descension, distance
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: nnind
      REAL(q2), DIMENSION(8,3) :: nngrad
      INTEGER :: j, stepcount, nnLayers
      REAL(q2), DIMENSION(3) :: tempr, rn, rnm1 ! rn minus 1
      REAL(q2), DIMENSION(3) :: grad, stepsize, gradnm1
      INTEGER, DIMENSION(3,26) :: vi, matw
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(3,3) :: matm, outerproduct
      stepcount = 0
      ALLOCATE(nnInd(8,3))
!      ALLOCATE(nnInd((nnlayers*2+1)**3,3))
      ! initialize the process
      CALL pbc(ind,chg%npts)
      rn(1) = REAL(ind(1),q2)
      rn(2) = REAL(ind(2),q2)
      rn(3) = REAL(ind(3),q2)
      stepsize(1) = 0.5
      stepsize(2) = 0.5
      stepsize(3) = 0.5
      IF (opts%leastsquare_flag) THEN
        grad = lsg(ind,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      ELSE
        grad = cdgrad(ind,chg)
      END IF
!      PRINT *, 'initial grad is'
!      PRINT *, grad
      ! this gradient is in cartesian. convert it to lattice
      grad = MATMUL(chg%lat2car,grad)
      DO WHILE (stepsize(1) >= 0.1 .AND. &
                stepsize(2) >= 0.1 .AND. &
                stepsize(3) >= 0.1 )
        ! the gradient is in cartesian. 
!        PRINT *, 'grad is'
!        PRINT *, grad
!        PRINT *, 'rn is'
!        PRINT *, rn
        IF (SUM(grad*grad) <= (0.1*opts%knob_gradfloor)**2) THEN
          PRINT *, 'descension grad small enough'
          EXIT
        END IF
        gradnm1 = grad
        rnm1 = rn
        ! determine where to go, a unit vector
        tempr = grad / SQRT(SUM(grad*grad))
        
        tempr(1) = stepsize(1) * tempr(1)
        tempr(2) = stepsize(2) * tempr(2)
        tempr(3) = stepsize(3) * tempr(3)
        rn = rn - tempr
        CALL pbc_r_lat(rn,chg%npts)
        nnind = simpleNN(rn,chg)
        DO j = 1,8
          CALL pbc(nnind(j,:),chg%npts)
          !IF (opts%leastsquare_flag .EQV. .TRUE.) THEN
          IF ( .TRUE. ) THEN
            !nngrad(j,:) = lsg(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
            nngrad(j,:) = cdgrad(nnind(j,:),chg)
          ELSE 
          END IF
        END DO
        distance = rn - nnind(1,:)
        grad = trilinear_interpol_grad(nngrad,distance) ! val r interpol
        !nnInd = FindNN(rn,nnLayers,chg,ions)
        !grad = R2GradInterpol(nnInd,rn,chg,nnLayers)
        ! this grad is in cartesian. convert it to lattice
        grad = MATMUL(chg%lat2car,grad)
        IF (SUM(grad*grad)<=(0.1*opts%knob_gradfloor)**2) THEN
          ! we are at a critical point!
!          PRINT *, 'descention gradient sufficiently small'
!          PRINT *, 'descension rn is'
!          PRINT *, rn
          EXIT
        END IF
        ! if grad points backwards, reduce stepsize
        IF (grad(1) * gradnm1(1) < 0) THEN
          stepsize(1) = stepsize(1)/2
        END IF
        IF (grad(2) * gradnm1(2) < 0) THEN
          stepsize(2) = stepsize(2)/2
        END IF
        IF (grad(3) * gradnm1(3) < 0) THEN
          stepsize(3) = stepsize(3)/2
        END IF
        stepcount = stepcount + 1
      END DO 
      descension = rn
      RETURN 
    END FUNCTION 


    SUBROUTINE minimasearch(chg,cptnum,cpl,bdr,matm,matwprime, &
                          wi,vi,vit,ggrid,outerproduct,opts,ucagecount,&
                          nnLayers, ions)
      TYPE(charge_obj) :: chg
      TYPE(bader_obj) :: bdr
      TYPE(options_obj) :: opts
      TYPE(ions_obj) :: ions
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      INTEGER :: cptnum, ucagecount, nnLayers
      INTEGER :: n1,n2,n3,i
      INTEGER, DIMENSION(26,3) :: nn
      LOGICAL :: isminimum
      INTEGER, DIMENSION(3) :: minpos
      REAL(q2), DIMENSION(3) :: minr
      INTEGER, DIMENSION(3,26) :: vi, matw
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(3,3) :: matm, outerproduct
      PRINT *, 'Performing initial search for minima'
      DO n1 = 1 , chg%npts(1)  
        DO n2 = 1 , chg%npts(2)
          DO n3 = 1 , chg%npts(3)
            ! do not search in vacuum
            IF (bdr%volnum(n1,n2,n3) == bdr%bnum + 1) THEN
              CYCLE
            END IF
!            IF ( .NOT. is_vol_edge(bdr,chg,(/n1,n2,n3/))) THEN
!              CYCLE
!            END IF
            IF ( ProxyToCPCandidate((/n1,n2,n3/),opts,cpl,cptnum,chg,nnLayers)) THEN
              CYCLE
            END IF
            isminimum = .TRUE.
            DO i = 1,26
              nn(i,:) = (/n1,n2,n3/) + vi(:,i)
            END DO
            DO i = 1, 26
              CALL pbc(nn(i,:),chg%npts)
              IF (rho_val(chg,n1,n2,n3) >= rho_val(chg,nn(i,1),nn(i,2),nn(i,3))) THEN
                isminimum = .FALSE.
              END IF
            END DO
            IF ( isminimum .EQV. .TRUE.  ) THEN
              cptnum = cptnum + 1
              ucagecount = ucagecount + 1
              minpos = (/n1,n2,n3/)
!              WRITE(98,*) '_________________________________________'
!              WRITE(98,*) 'Though minimasearch, critical point found at'
!              WRITE(98,*) n1,n2,n3
!              WRITE(98,*) 'Coordinates in cartesian are'
!              WRITE(98,*) MATMUL((/n1,n2,n3/),chg%lat2car)
!              WRITE(98,*) 'Direct coordinates are'
!              WRITE(98,*)  n1/chg%npts(1),n2/chg%npts(2), &
!                           n3/chg%npts(3)
!              WRITE(98,*) 'Charge density here is'
!              WRITE(98,*) rho_val(chg,n1,n2,n3)
              minr = descension(minpos,chg,matm,matwprime, &
                          wi,vi,vit,ggrid,outerproduct,opts,nnLayers,ions)
!              WRITE(98,*) 'descends to'
!              WRITE(98,*) minr
!              WRITE(98,*) '_________________________________________'
!              WRITE(98,*) ' ' 
              
              cpl(cptnum)%ind = (/n1,n2,n3/)
              

            END IF
          END DO
        END DO
      END DO
      PRINT *, 'minima search completed'
    END SUBROUTINE

    ! follow the gradient down to a minimum of the squared gradient of the
    ! charge density
    SUBROUTINE sqgradientdescend(r0,chg,matm,matwprime,wi,vi,vit, &
                       ggrid,outerproduct,opts,rn,invac,bdr,nnLayers,ions)
      ! this function should find all critical points
      TYPE(bader_obj) :: bdr
      TYPE(ions_obj) :: ions
      INTEGER, DIMENSION(3) :: r0, rt
      INTEGER, DIMENSION(3) :: crossings ! count how many times pbc crossed
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      REAL(q2),DIMENSION(3) ::  distance
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: nnInd
      REAL(q2), DIMENSION(8,3) :: nngrad
      INTEGER :: i, j, stepcount, loopcount, nnLayers
      REAL(q2), DIMENSION(3) :: tempr, rn, trn, rnm1 ! rn minus 1
      REAL(q2), DIMENSION(3) :: grad,  gradnm1, avgrn
      INTEGER, DIMENSION(3,26) :: vi, matw
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(3,3) :: matm, outerproduct
      REAL(q2), DIMENSION(3) :: p, pbccorrectionr ! keep track of how much is 
                                ! added to the r to have it fall within the pbc
      REAL(q2), DIMENSION(3) :: stepsize
      LOGICAL :: invac, isaveraging
      ALLOCATE(nnInd(8,3))
!      ALLOCATE(nnInd((nnlayers*2+1)**3,3))
      loopcount = 0
      pbccorrectionr = 0
      PRINT *, 'r0 is'
      PRINT *, r0
      invac = .FALSE.
      isaveraging = .FALSE.
      crossings = 0
      CALL pbc(r0,chg%npts)
      ! this is the initial grad
      grad = lsgsqlsg(r0,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      ! this gradient is in cartesian. 
      stepcount = 0
      ! initialize the process
      rn(1) = REAL(r0(1),q2)
      rn(2) = REAL(r0(2),q2)
      rn(3) = REAL(r0(3),q2)
      stepsize = (/0.5,0.5,0.5/)
      DO WHILE (stepsize(1) >= 0.01 .OR. &
                stepsize(2) >= 0.01 .OR. &
                stepsize(3) >= 0.01)
        loopcount = loopcount + 1
        PRINT *, 'loop ', loopcount
        PRINT *, 'stepsize is'
        PRINT *, stepsize
        ! the gradient is in cartesian. 
        gradnm1 = grad
        rnm1 = rn
        ! determine where to go, a unit vector
        IF (ABS(grad(1)) < 0.001 &
            .AND. ABS(grad(2)) <= 0.001 &
            .AND. ABS(grad(3)) <= 0.001)  THEN
          ! we are at a critical point!
          PRINT *, 'gradient sufficiently small'
          PRINT *, 'grad is'
          PRINT *, grad
          PRINT *, 'rn is'
          PRINT *, rn
          isaveraging = .FALSE.
          EXIT
        END IF
        tempr = grad / SQRT(SUM(grad*grad))
        tempr(1) = MIN(stepsize(1), ABS(grad(1))) *  tempr(1)
        tempr(2) = MIN(stepsize(2), ABS(grad(2))) *  tempr(2)
        tempr(3) = MIN(stepsize(3), ABS(grad(3))) *  tempr(3)
        PRINT *, 'tempr is'
        PRINT *, tempr
        rn = rn - tempr
        PRINT *, 'rn is'
        PRINT *, rn
        trn = rn
        CALL pbc_r_lat(rn,chg%npts)
        pbccorrectionr = pbccorrectionr + rn - trn
        IF (isaveraging) THEN
          PRINT *, 'accumulating average'
          stepcount = stepcount + 1
          PRINT *, 'step ', stepcount
          ! if rn is hopping back and forth at the boundry
          IF (trn(1) > rn(1) ) THEN
            crossings(1) = crossings(1) + 1
          ELSE IF (trn(1) < rn(1)) THEN
            crossings(1) = crossings(1) - 1
          END IF
          IF (trn(2) > rn(2) ) THEN
            crossings(2) = crossings(2) + 1
          ELSE IF (trn(2) < rn(2)) THEN
            crossings(2) = crossings(2) - 1
          END IF
          IF (trn(3) > rn(3) ) THEN
            crossings(3) = crossings(3) + 1
          ELSE IF (trn(3) < rn(3)) THEN
            crossings(3) = crossings(3) - 1
          END IF
          IF (crossings(1) > 0) THEN
            avgrn(1) = avgrn(1) + rn(1) + chg%npts(1)
            PRINT *, 'compensating for leftward PBC'
          ELSE IF (crossings(1) < 0) THEN
            avgrn(1) = avgrn(1) + rn(1) - chg%npts(1)
            PRINT *, 'compensating for rightward PBC'
          ELSE 
            avgrn(1) = avgrn(1) + rn(1)
          END IF
          IF (crossings(2) > 0) THEN
            avgrn(2) = avgrn(2) + rn(2) + chg%npts(2)
            PRINT *, 'compensating for leftward PBC'
          ELSE IF (crossings(2) < 0) THEN
            avgrn(2) = avgrn(2) + rn(2) - chg%npts(2)
            PRINT *, 'compensating for rightward PBC'
          ELSE 
            avgrn(2) = avgrn(2) + rn(2)
          END IF
          IF (crossings(3) > 0) THEN
            avgrn(3) = avgrn(3) + rn(3) + chg%npts(3)
            PRINT *, 'compensating for leftward PBC'
          ELSE IF (crossings(3) < 0) THEN
            avgrn(3) = avgrn(3) + rn(3) - chg%npts(3)
            PRINT *, 'compensating for rightward PBC'
          ELSE 
            PRINT *, 'new sum'
            avgrn = avgrn + rn
            PRINT *, avgrn
          END IF
        END IF
        ! if the point is in vacuum, stop the descend
        rt(1) = NINT(rn(1))
        rt(2) = NINT(rn(2))
        rt(3) = NINT(rn(3))
        CALL pbc(rt,chg%npts)
        IF (bdr%volnum(rt(1),rt(2),rt(3)) == bdr%bnum + 1) THEN
          rn = (/-1.,-1.,-1./)
          invac = .TRUE.
          PRINT *, 'legacy of the void'
          EXIT 
        END IF
        ! if traveled too far, stop the descend
        PRINT *, 'total distance travelled'
        PRINT *, ABS(rn(1) - REAL(r0(1),q2) - pbccorrectionr(1)), & 
                 ABS(rn(2) - REAL(r0(2),q2) - pbccorrectionr(2)), &
                 ABS(rn(3) - REAL(r0(3),q2) - pbccorrectionr(3))
        IF (ABS(rn(1) - REAL(r0(1),q2) - pbccorrectionr(1)) >= 3 + opts%knob_tem * 3 .OR. &
            ABS(rn(2) - REAL(r0(2),q2) - pbccorrectionr(2)) >= 3 + opts%knob_tem * 3 .OR. &
            ABS(rn(3) - REAL(r0(3),q2) - pbccorrectionr(3)) >= 3 + opts%knob_tem * 3 ) THEN
            rn = (/-1.,-1.,-1./)
            PRINT *, 'youve gone too far'
            EXIT
        END IF
        nnind = simpleNN(rn,chg)
        DO j = 1,8
          CALL pbc(nnind(j,:),chg%npts)
          !IF (opts%leastsquare_flag .EQV. .TRUE.) THEN
          IF ( opts%leastSquare_flag ) THEN
            !nngrad(j,:) = lsg(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
            nngrad(j,:) =  lsgsqlsg(nnind(j,:),chg,matm,matwprime,wi &
                           ,vi,vit,ggrid,outerproduct)
          ELSE 
            nngrad(j,:) = CDGrad(nnind(j,:),chg)
          END IF
        END DO
        distance = rn - nnind(1,:)
        grad = trilinear_interpol_grad(nngrad,distance) ! val r interpol
        !nnind = FindNN(rn,nnLayers,chg,ions)
        !grad = R2GradInterpol(nnind,rn,chg,nnLayers)
                ! this grad is in cartesian. convert it to lattice
        ! if grad points backwards, reduce stepsize
!        IF (grad(1) * gradnm1(1) <= 0 .OR. &
!            grad(2) * gradnm1(2) <= 0 .OR. & 
!            grad(3) * gradnm1(3) <= 0 ) THEN
!            stepsize = stepsize * 0.5
!            isaveraging = .TRUE.
!        END IF        
        IF (grad(1) * gradnm1(1) <= 0 ) THEN
          stepsize(1) = stepsize(1) * 0.5
        END IF
        IF (grad(2) * gradnm1(2) <= 0 ) THEN
          stepsize(2) = stepsize(2) * 0.5
        END IF
        IF (grad(3) * gradnm1(3) <= 0 ) THEN
          stepsize(3) = stepsize(3) * 0.5
        END IF
      END DO 
      PRINT *, 'finished looping, rn is'
      PRINT *, rn
      PRINT *, 'final stepsize is'
      PRINT *, stepsize
      IF (isaveraging) THEN
        rn(1) = avgrn(1) / stepcount
        rn(2) = avgrn(2) / stepcount
        rn(3) = avgrn(3) / stepcount
        CALL pbc_r_lat(rn,chg%npts)
        PRINT *, 'rn is averaged to be'
        PRINT *, rn
      END IF
    END SUBROUTINE sqgradientdescend



    ! this function gives the ls gradient of the squared ls gradient
    FUNCTION lsgsqlsg(r0,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      TYPE(charge_obj) :: chg
      REAL(q2), DIMENSION(3,3) :: lsh
      INTEGER, DIMENSION(3) :: r0
      REAL(q2), DIMENSION(3,26) :: lsg_val
      INTEGER, DIMENSION(3,26) :: vi
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      INTEGER, DIMENSION(3) :: nbp
      REAL(q2), DIMENSION(26) :: wi, totlsg_val
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(13) :: deltadx, deltady, deltadz, delta
      REAL(q2), DIMENSION(3,3) :: matm
      REAL(q2), DIMENSION(3,3) :: outerproduct
      REAL(q2) :: average
      REAL(q2), DIMENSION(3) :: lsgsqlsg
      INTEGER :: i, j, k
      ! first need the gradient of the local neighbors, like the parts of lsh.
      DO i = 1, 26
        nbp = vi(:,i) + r0
        CALL pbc(nbp,chg%npts) 
        lsg_val(:,i) = lsg(nbp,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
        ! square the gradient
        lsg_val(1,i) = lsg_val(1,i) ** 2 
        lsg_val(2,i) = lsg_val(2,i) ** 2
        lsg_val(3,i) = lsg_val(3,i) ** 2
        totlsg_val(i) = lsg_val(1,i) + lsg_val(2,i) + lsg_val(3,i)
      END DO
      ! get the gradient of the squared gradient
      deltadx = 0.
      deltady = 0.
      deltadz = 0.
      DO i = 1, 13
!        matwprime(:,i) = vi(:,i) * wi(i)
        ! matwprime matches with lsg
!        deltadx(i) = lsg_val(1,i) - lsg_val(1,i+13)
!        deltady(i) = lsg_val(2,i) - lsg_val(2,i+13)
!        deltadz(i) = lsg_val(3,i) - lsg_val(3,i+13)
        delta(i) = totlsg_val(i) - totlsg_val(i+13)
      END DO
      lsgsqlsg= 0._q2
      lsgsqlsg = &
        MATMUL( &
          MATMUL( & 
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime) & 
        , delta)
      ! lsg up till now seems to be larger than force in cartesian units by car2lat
      lsgsqlsg = MATMUL(chg%lat2car,lsgsqlsg) ! now force should be in cartesian units
      END FUNCTION lsgsqlsg



      FUNCTION gradientfilter(p,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      TYPE(charge_obj) :: chg
      REAL(q2), DIMENSION(3,3) :: lsh
      INTEGER, DIMENSION(3) :: p
      REAL(q2), DIMENSION(3,26) :: lsg_val
      INTEGER, DIMENSION(3,26) :: vi, ps ! indexes of 26 neighbors
      REAL(q2), DIMENSION(3,26) ::  grads ! gradient of 26 neighbors
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      INTEGER, DIMENSION(3) :: nbp
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(13) :: deltadx, deltady, deltadz
      REAL(q2), DIMENSION(3,3) :: matm
      REAL(q2), DIMENSION(3,3) :: outerproduct
      REAL(q2) :: average
      REAL(q2), DIMENSION(3) :: grad
      INTEGER :: i, j, k
      LOGICAL :: gradientfilter
      gradientfilter = .TRUE.
      DO i = 1, 26
        ps(:,i) = p + vi(:,i)
        CALL pbc(ps(:,i),chg%npts)
        grads(:,i) = lsg(ps(:,i),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      END DO
!      DO i = 1, 13
!        IF (grads(1,i)*grads(1,i+13) >= 0 .OR. &
!            grads(2,i)*grads(2,i+13) >= 0 .OR. & 
!            grads(3,i)*grads(3,i+13) >= 0 ) THEN
!           gradientfilter = .FALSE.
!           EXIT
!        END IF
!      END DO
      IF (grads(1,5)*grads(1,18) >= 0 &
          .OR. grads(2,11)*grads(2,24)>=0 &
          .OR. grads(3,13)*grads(3,26)>=0 ) THEN
        gradientfilter = .FALSE.
      END IF

      RETURN
    END FUNCTION
  
    ! this function takes in the current point, see if a near by point has
    ! already been marked as a critical point. if yes, skip it. The criteria is
    ! set to be based on halfo the search radius, which is 1 + knob_tem
    ! THIS FUNCTION AS OF 20191221 WONT WORK WELL AROUND PBC
    ! As a loose first round checking, it is OK that a few points near the PBC
    ! are permitted into the candidacy. 
    ! Just like the Democratic 2020 primary, not everyone has to be super
    ! qualified to enter. 
    FUNCTION ProxyToCPCandidate(p,opts,cpl,cptnum,chg,nnLayers)
      LOGICAL :: ProxyToCPCandidate
      INTEGER :: i
      INTEGER, DIMENSION(3) :: p
      TYPE(options_obj) :: opts
      TYPE(cpc), ALLOCATABLE, DIMENSION(:) :: cpl
      TYPE(charge_obj) :: chg
      INTEGER :: cptnum, nnLayers
      ProxyToCPCandidate = .FALSE.
      ! These are lattice based and do not take PBC into consideration
      DO i = 1, cptnum
        IF (ABS(cpl(i)%ind(1) - p(1)) <= opts%par_sr .AND. &
            ABS(cpl(i)%ind(2) - p(2)) <= opts%par_sr .AND. &
            ABS(cpl(i)%ind(3) - p(3)) <= opts%par_sr ) THEN
          ProxyToCPCandidate = .TRUE.
        END IF
      END DO  
      ! These are very costly
      !DO i = 1, cptnum
      !  IF ( GetPointDistance(p,cpl(i)%ind,chg,nnLayers) <= 0.2) THEN
      !    ProxyToCPCandidate = .TRUE.
      !    EXIT
      !  END IF
      !END DO
      RETURN
    END FUNCTION ProxyToCPCandidate

    ! this function should always run. it first checks if things are stuck
    ! this function is for when newton raphson hops between two points
    FUNCTION unstuck( nexttem, previoustem, temscale, temnormcap)
      REAL(q2), DIMENSION(3) :: nexttem, previoustem, temscale, unstuck
      REAL(q2) :: temnormcap
      LOGICAL :: stuck
      LOGICAL :: capprev,capnext ! see if tem was or will be capped
      stuck = .FALSE. ! we check for stuckness
      capprev = .FALSE.
      capnext = .FALSE.
      !temnormcap is all rounded cap designed to prevent overstepping.
      !temnormcap should not change at all
      !temscale is direction sensitive. It gradually narrows in if tem changes
      !sign a lot.


      IF (mag(previoustem) == temnormcap) THEN
        capprev = .TRUE.
      END IF
      IF (mag(nexttem) == temnormcap ) THEN
        capnext = .TRUE.
      END IF

      ! if nexttem and prev tem are complete opposites, it's stuck for sure
      IF (nexttem(1) == -previoustem(1) .AND. &
          nexttem(2) == -previoustem(2) .AND. &
          nexttem(3) == -previoustem(3) ) THEN
        stuck = .TRUE.
      END IF

      IF ( stuck ) THEN
        nexttem = rngkick(nexttem)
      END IF
 
      unstuck = nexttem
      ! stuck is defined as nexttem and previous tem being equal but complete
      ! opposite

      ! if hopping back and forth without gradient changing sign
!      IF ( (nexttem(1)*previoustem(1)<=0 .OR.     &
!            nexttem(2)*previoustem(2)<=0 .OR.     &
!            nexttem(3)*previoustem(3)<=0 ) .AND. (&
!            (nexttem(1) + previoustem(1))**2 +    &
!            (nexttem(2) + previoustem(2))**2 +    &
!            (nexttem(3) + previoustem(3))**2      &
!            <= 0.5) .AND. .NOT. (                 &
!            ABS(nexttem(1)) == temcap(1) .AND.    &
!            ABS(nexttem(2)) == temcap(2) .AND.    &
!            ABS(nexttem(3)) == temcap(3)          &
!            ))  THEN
!!        PRINT *, '' 
!        CALL RANDOM_NUMBER(rndr)
!        nexttem(1) = nexttem(1) + rndr - 0.5
!        CALL RANDOM_NUMBER(rndr)
!        nexttem(2) = nexttem(2) + rndr - 0.5
!        CALL RANDOM_NUMBER(rndr)
!        nexttem(3) = nexttem(3) + rndr - 0.5
!      END IF
!      PRINT *, 'after comparing with previoustem, tem became'
!      PRINT *, nexttem
!      IF (nexttem(1) == previoustem(1) .AND. &
!          nexttem(2) == previoustem(2) .AND. &
!          nexttem(3) == previoustem(3) .AND. & 
!          .NOT. (                            &
!            ABS(nexttem(1)) == temcap(1) .AND.    &
!            ABS(nexttem(2)) == temcap(2) .AND.    &
!            ABS(nexttem(3)) == temcap(3)          &
!         )) THEN
!        CALL RANDOM_NUMBER(rndr)
!        nexttem(1) = nexttem(1) + rndr - 0.5
!        CALL RANDOM_NUMBER(rndr)
!        nexttem(2) = nexttem(2) + rndr - 0.5
!        CALL RANDOM_NUMBER(rndr)
!        nexttem(3) = nexttem(3) + rndr - 0.5
!      END IF
!      PRINT *, 'after another check, nexttem eventually came out as'
!      PRINT *, nexttem
      RETURN
    END FUNCTION

    ! this function takes in a vector, give it a random scale within 10%
    FUNCTION rngkick(nexttem)
      REAL(q2), DIMENSION(3) :: nexttem, rngkick
      REAL(q2) :: ran
      INTEGER :: i
      DO i = 1,3
        CALL RANDOM_NUMBER(ran)
        nexttem(1) = ((0.5-ran)*0.2 + 1)*nexttem(1)
      END DO
      rngkick = nexttem
      RETURN
    END FUNCTION

    ! give the magnitude of a 3d vector
    FUNCTION mag(vec3d)
      REAL(q2), DIMENSION(3) :: vec3d
      REAL(q2) :: mag
      mag = SQRT(SUM(vec3d*vec3d))
      RETURN
    END FUNCTION

    ! this function takes in the movement factor, modifies it as necessary
    FUNCTION TemMods(nexttem,temscale,temnormcap)
      REAL(q2), DIMENSION(3) :: TemMods, nexttem, temscale
      REAL(q2) :: temnormcap
      INTEGER :: i
      TemMods = nexttem
      DO i = 1 , 3
        TemMods(i) = TemMods(i) * temscale(i)
      END DO
      IF (mag(TemMods) > temNormCap) THEN
        TemMods = TemMods/mag(TemMods) * temScale
      END IF
      RETURN
    END FUNCTION
   
    ! this function inspects if it is beneficial to reduce temscale
    ! right now it does nothing because I'm not sure if limiting it helps in any
    ! way at all
    FUNCTION scaleinspector( nexttem, previoustem, temScale)
      REAL(q2), DIMENSION(3) :: nexttem, previoustem, temScale, scaleinspector
      INTEGER :: i
      !DO i = 1, 3
      !  IF (nexttem(i)*previoustem(i) <= 0) THEN
      !    temScale(i) = temScale(i) * 0.5
      !  END IF
      !END DO
      !IF (nexttem(1)*previoustem(1)<=0 .AND. &
      !    nexttem(1)*previoustem(1)<=0 .AND. &
      !    nexttem(1)*previoustem(1)<=0 ) THEN
      !    temScale = temScale * 0.5
      !END IF
      scaleinspector = temScale
      RETURN
    END FUNCTION

    ! if it is detected that the hessian 
    FUNCTION unZeroHessian(truer,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      TYPE(charge_obj) :: chg
      TYPE(hessian) :: hes
      REAL(q2), DIMENSION(3,3) :: unZeroHessian, lsd
      INTEGER, DIMENSION(3) :: truer
      INTEGER, DIMENSION(3,26) :: vi
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      INTEGER, DIMENSION(3) :: nbp
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(3,3) :: matm
      REAL(q2), DIMENSION(3,3) :: outerproduct
      INTEGER, DIMENSION(8,3) :: nnind
      REAL(q2) :: average, ran
      INTEGER :: i, j, k
      
      ! first step is to move a little bit
      DO i = 1 , 3
        CALL RANDOM_NUMBER(ran)
        truer(i) = truer(i) + (ran - 0.5)*0.000001
      END DO
      ! second step is to find the new nearest neighbors
      PRINT *, truer
      lsd = lsh(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      unZeroHessian = 0
      RETURN
    END FUNCTION 

    ! check if the number of each type of critical point satisfies the
    ! Poincare-Hopf rule
    FUNCTION PHRuleChecker(maxcount, ubondcount, uringcount, ucagecount, opts, &
      ions )
      LOGICAL :: PHRuleChecker
      INTEGER :: maxcount, ubondcount, uringcount, ucagecount
      TYPE(options_obj) :: opts
      TYPE(ions_obj) :: ions
      phrulechecker = .FALSE.
      IF ( maxcount < ions%nions) THEN
        PHRULECHECKER = .FALSE.
      ELSEIF (maxcount > ions%nions) THEN
        PRINT *, 'WARNING : MORE MAXIMA FOUND THAN NUMBER OF NUCLEI!'
      ELSE 
        IF (opts%iscrystal ) THEN
          IF (maxcount - ubondcount + uringcount - ucagecount == 0) THEN
            phrulechecker = .TRUE.
          END IF
        ELSE 
          IF (maxcount - ubondcount + uringcount - ucagecount == 1) THEN
            phrulechecker = .TRUE.
          END IF 
        END IF
      END IF
      RETURN
    END FUNCTION PHRuleChecker
    
    ! This subroutine checks if PH rule is satisfied given crystal/molecule
    ! inport or not
    SUBROUTINE PHRuleExam(maxCount,bondCount,ringCount,cageCount,opts,ions,&
      phmrCompliant)
      TYPE(options_obj) :: opts
      TYPE(ions_obj) :: ions
      INTEGER :: maxCount,bondCount,ringCount,cageCount,phSum,iphsum
      LOGICAL :: phmrCompliant
      phSum = maxCount - bondCount + ringCount - cageCount
      iphSum = ions%nions - bondCount + ringCount - cageCount
      phmrCompliant = .FALSE.
      IF (opts%isCrystal) THEN
        PRINT *, 'The system is assigned as a Crystal'
        IF (phSum == 0) THEN
          PRINT *, ''//achar(27)//'[32m Satisfies the Morse Relationship' &
            //achar(27)//'[0m'
          phmrCompliant = .TRUE.
        ELSE IF (phSum == 1) THEN
          PRINT *, ''//achar(27)//'[31m ERROR: The result satisfies the  & 
            Poincare Hopf Rule for a & molecule, not a crystal.' &
            //achar(27)//'[0m'
        ELSE IF (iphSum == 0 .AND. iphSum /= phSum) THEN
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[32m Satisfies the Morse Relationship' &
            //achar(27)//'[0m'
          phmrCompliant = .TRUE.
        ELSE IF (iphSum == 1 .AND. iphSum /= phSum) THEN
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[31m ERROR: The result satisfies the  & 
            Poincare Hopf Rule for a & molecule, not a crystal.' &
            //achar(27)//'[0m'
        ELSE 
          PRINT *, ''//achar(27)//'[31m ERROR: FAILED Morse relationship' &
            //achar(27)//'[0m'
        END IF
      ELSE IF (opts%isMolecule) THEN
        PRINT *, 'The system is assigned as a Molecule'
        IF (phSum == 0) THEN
          PRINT *, ''//achar(27)//'[31m ERROR: The result satisfies the Morse Relationship &
            for a crystal, not a molecule.'//char(27)//'[0m'
        ELSE IF (phSum == 1) THEN
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[32m Satisfies the Poincare Hopf Rule' &
            //achar(27)//'[0m'
        ELSE IF (iphSum == 0 .AND. iphSum /= phSum) THEN
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[31m ERROR: The result satisfies the Morse Relationship &
            for a crystal, not a molecule.'//char(27)//'[0m'
        ELSE IF (iphSum == 1 .AND. iphSum /= phSum) THEN
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[32m Satisfies the Poincare Hopf Rule' &
            //achar(27)//'[0m'
        ELSE 
          PRINT *, ''//achar(27)//'[31m ERROR: FAILED Poincare Hopf Rule' &
            //achar(27)//'[0m'
        END IF
      ELSE 
        IF (phSum == 0) THEN
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[32m This system has not be designated & 
            as a molecule or crystal but the Morse relationship for & 
            crystals are satisfied.' //achar(27)//'[0m'
        ELSE IF (phSum == 1) THEN
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[32m This system has not be designated & 
            as a molecule or crystal but the Poincare-Hopf rule for &
            moleculess are satisfied.' //achar(27)//'[0m'
        END IF
        IF (iphSum == 0 .AND. iphSum /= phSum) THEN
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[32m This system has not be designated & 
            as a molecule or crystal but the Morse relationship for & 
            crystals are satisfied.' //achar(27)//'[0m'
        ELSE IF (iphSum == 1 .AND. iphSum /= phSum) THEN
          phmrCompliant = .TRUE.
          PRINT *, ''//achar(27)//'[33m Using the number of atoms  &
            instead of the number of the number of maxima found' &
            //achar(27)//'[0m'
          PRINT *, ''//achar(27)//'[32m This system has not be designated & 
            as a molecule or crystal but the Poincare-Hopf rule for &
            moleculess are satisfied.' //achar(27)//'[0m'
        END IF
        IF (.NOT. phmrCompliant) THEN
          PRINT *, ''//achar(27)//'[31m ERROR: FAILED Poincare Hopf Rule & 
            and Morse Relationship' //achar(27)//'[0m'
        END IF
      END IF
    END SUBROUTINE 
    ! count the number of negative eigenvalues to characterize a critical point
    SUBROUTINE RecordCP(p,chg,matm,matwprime,wi,vi,vit,ggrid &
      ,outerproduct,cpl,ucptnum,eigvals,eigvecs, maxcount, uringcount, &
      ubondcount, ucagecount,opts)
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      REAL(q2), DIMENSION(3) :: eigvals, grad
      REAL(q2), DIMENSION(3,3) :: eigvecs, hessianMatrix
      INTEGER, DIMENSION(3) :: p
      INTEGER, DIMENSION(3,26) :: vi
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(3,3) :: matm
      REAL(q2), DIMENSION(3,3) :: outerproduct
      INTEGER :: i, j, k, ucptnum, negcount
      INTEGER :: maxcount, uringcount, ubondcount, ucagecount
      IF (opts%leastSquare_flag) THEN
        grad = lsg( &
          p,chg,matm,matwprime, &
          wi,vi,vit,ggrid,outerproduct)
        hessianMatrix = lsh( &
          p,chg,matm,matwprime, &
          wi,vi,vit,ggrid,outerproduct)
      ELSE
        grad = CDGrad(p,chg)
        hessianMatrix = CDHessian(p,chg)
      END IF
      CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
      negCount = CountNegModes(eigvals)
      CALL UpDateCounts(negCount,maxCount,uBondCount,uRingCount,uCageCount)
      cpl(ucptnum)%ind = p
      cpl(ucptnum)%truer(1) = REAL(p(1),q2)
      cpl(ucptnum)%truer(2) = REAL(p(2),q2)
      cpl(ucptnum)%truer(3) = REAL(p(3),q2)
      cpl(ucptnum)%grad = grad
      cpl(ucptnum)%hessianMatrix = hessianMatrix
      cpl(ucptnum)%eigvecs = eigvecs
      cpl(ucptnum)%eigvals = eigvals
      cpl(ucptnum)%negCount = negCount
    END SUBROUTINE RecordCP
    ! the version of the above subroutine where p is real not integer
    SUBROUTINE RecordCPR(p,chg,cpl,ucptnum,eigvals,eigvecs, maxcount, uringcount, &
      ubondcount, ucagecount,opts,grad,hessianMatrix)
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      REAL(q2), DIMENSION(3) :: eigvals, grad
      REAL(q2), DIMENSION(3,3) :: eigvecs, hessianMatrix
      REAL(q2), DIMENSION(3) :: p
      INTEGER, DIMENSION(3,26) :: vi
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(3,3) :: matm
      REAL(q2), DIMENSION(3,3) :: outerproduct
      INTEGER :: i, j, k, ucptnum, negCount
      INTEGER :: maxcount, uringcount, ubondcount, ucagecount
      cpl(ucptnum)%hessianMatrix = hessianMatrix
      CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
      negCount = CountNegModes(eigvals)
      CALL UpDateCounts(negCount,maxCount,uBondCount,uRingCount,uCageCount)
!      PRINT *, 'in recprd CP, updated counts are'
!      PRINT *, 'bond, ring, cage, max'
!      PRINT *, ubondcount, uringcount, ucagecount, maxcount
      cpl(ucptnum)%truer = p
      cpl(ucptnum)%grad = grad
      cpl(ucptnum)%eigvecs = eigvecs
      cpl(ucptnum)%eigvals = eigvals
      cpl(ucptnum)%negcount = negcount
    END SUBROUTINE RecordCPR
    
    SUBROUTINE RecordCPRLight(p,chg,cpl,ucptnum, maxcount, uringcount, &
      ubondcount, ucagecount)
      TYPE(charge_obj) :: chg
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      REAL(q2), DIMENSION(3) :: eigvals, grad
      REAL(q2), DIMENSION(3,3) :: eigvecs, hessianMatrix
      REAL(q2), DIMENSION(3) :: p
      INTEGER, DIMENSION(3,26) :: vi
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(3,3) :: matm
      REAL(q2), DIMENSION(3,3) :: outerproduct
      INTEGER :: i, j, k, ucptnum, negCount
      INTEGER :: maxcount, uringcount, ubondcount, ucagecount
      hessianMatrix = CDHessianR(p,chg)
      grad = CDGradR(p,chg)
      cpl(ucptnum)%hessianMatrix = hessianMatrix
      CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
      negCount = CountNegModes(eigvals)
      CALL UpDateCounts(negCount,maxCount,uBondCount,uRingCount,uCageCount)
!      PRINT *, 'in recprd CP, updated counts are'
!      PRINT *, 'bond, ring, cage, max'
!      PRINT *, ubondcount, uringcount, ucagecount, maxcount
      cpl(ucptnum)%truer = p
      cpl(ucptnum)%grad = grad
      cpl(ucptnum)%eigvecs = eigvecs
      cpl(ucptnum)%eigvals = eigvals
      cpl(ucptnum)%negcount = negcount
    END SUBROUTINE RecordCPRLight

    SUBROUTINE  outputCP(cpl,opts,ucptnum,chg,setcount,ubondcount, &
      uringcount,ucagecount, maxcount)
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl
      TYPE(options_obj) :: opts
      TYPE(charge_obj) :: chg
      INTEGER :: ucptnum, i, setcount
      INTEGER :: ubondcount, uringcount, ucagecount, maxcount
      CHARACTER(10) :: fileName
      PRINT *, 'Writting critical point output files'
      PRINT *, 'setcount is', setcount
      WRITE(fileName,fmt='(a,i2.2,a)') TRIM('CPFU'), setcount,TRIM('.dat')
      PRINT *, 'filename is ', filename
      OPEN(98,FILE=filename,STATUS='REPLACE',ACTION='WRITE')
      WRITE(98,*) 'Number of critical points written: ', ucptnum
      WRITE(98,*) 'Number of bond critical point written: ', ubondcount
      WRITE(98,*) 'Number of ring critical point written: ', uringcount
      WRITE(98,*) 'Number of cage critical point written: ', ucagecount
      WRITE(98,*) 'Number of nucl critical point written: ', maxcount
      DO i = 1, ucptnum
        WRITE (98,*) '_______________________________________'
        WRITE (98,*) 'Unique CP # : ', i
        WRITE (98,*) 'Critical point is found at indices'
        IF (opts%noInterpolation_flag) THEN
          WRITE (98,*) cpl(i)%ind
        ELSE 
          WRITE (98,*) cpl(i)%truer
        END IF
        WRITE (98,*) 'Coordinates in cartesian are'
        IF (opts%noInterpolation_flag) THEN
          WRITE (98,*) MATMUL(chg%lat2car,cpl(i)%ind)
        ELSE 
          WRITE (98,*) MATMUL(chg%lat2car,cpl(i)%truer)
        END IF
        WRITE (98,*) 'Direct coordinates are'
        IF (opts%noInterpolation_flag) THEN
          WRITE (98,*) cpl(i)%ind(1)/chg%npts(1), &
            cpl(i)%ind(2)/chg%npts(2), &
            cpl(i)%ind(3)/chg%npts(3)
        ELSE 
          WRITE (98,*) cpl(i)%truer(1)/chg%npts(1), &
            cpl(i)%truer(2)/chg%npts(2), &
            cpl(i)%truer(3)/chg%npts(3)
        END IF
        WRITE (98,*) 'Gradiant is'
        WRITE (98,*) cpl(i)%grad
        WRITE (98,*) 'Hessian is'
        WRITE (98,*) cpl(i)%hessianMatrix(1,:)
        WRITE (98,*) cpl(i)%hessianMatrix(2,:)
        WRITE (98,*) cpl(i)%hessianMatrix(3,:)
        WRITE (98,*) 'Eigenvalues are'
        WRITE (98,*) cpl(i)%eigvals
        WRITE (98,*) 'Eigenvectors are'
        WRITE (98,*) cpl(i)%eigvecs
        IF (cpl(i)%negcount == 0) THEN
          WRITE(98,*) 'This is a cage critical point'
          WRITE(98,*) ' '
        END IF
        IF (cpl(i)%negcount == 2) THEN
          WRITE(98,*) 'This is a bond critical point'
          WRITE(98,*) ' '
        ELSEIF(cpl(i)%negcount == 1) THEN
          WRITE(98,*) 'This is a ring critical point'
          WRITE(98,*) ' '
        ELSEIF(cpl(i)%negcount == 3) THEN
          WRITE(98,*) 'This is a nuclear critical point'
          WRITE(98,*) ' ' 
        END IF
        WRITE(98,*) '________________________________________'
      END DO
    WRITE(98,*) ''
    CLOSE(98)
    END SUBROUTINE outputCP


    ! get the distance between two points, considering that these two points can
    ! be neighbors across the periodic boundary
    ! The output distance is in cartesian
    FUNCTION GetPointDistance(p, pt, chg, nnlayer)
      REAL(q2) :: GetPointDistance, distance
      ! these are indices
      INTEGER, DIMENSION(3) :: p, pt
      INTEGER, DIMENSION( 3) :: nnp
      ! these are cartesian coordinates
      REAL(q2), DIMENSION(3) :: r, rt
      REAL(q2), DIMENSION(3) :: nnr
      TYPE(charge_obj) :: chg
      INTEGER :: nnlayer, i, j, k
      ! look for all equivalents of p within nnlayer, find the smallest distance
      rt = MATMUL(chg%lat2car,pt)
      r = MATMUL(chg%lat2car,p)
      GetPointDistance = mag( r - rt)
      DO i = -nnlayer, nnlayer
        DO j = -nnlayer, nnlayer
          DO k = -nnlayer, nnlayer
            IF (i == 0.AND.j == 0.AND.k == 0) CYCLE
            nnp = p + (/i * chg%npts(1), j * chg%npts(2), k * chg%npts(3)/)
            nnr = MATMUL(chg%lat2car,nnp)
            distance = mag( nnr - rt)
            GetPointDistance = MIN(GetPointDistance, distance)
          END DO
        END DO
      END DO
      RETURN
    END FUNCTION GetPointDistance

    ! This function takes two fractional lattice coordinate sets and produces
    ! the distance in cartesian
    FUNCTION GetPointDistanceR(p, pt, chg, nnlayer)
      REAL(q2) :: GetPointDistanceR, distance
      ! these are indices
      REAL(q2), DIMENSION(3) :: p, pt
      INTEGER, DIMENSION( 3) :: nnp
      ! these are cartesian coordinates
      REAL(q2), DIMENSION(3) :: r, rt
      REAL(q2), DIMENSION(3) :: nnr
      TYPE(charge_obj) :: chg
      INTEGER :: nnlayer, i, j, k
      ! look for all equivalents of p within nnlayer, find the smallest distance
      rt = MATMUL(chg%lat2car,pt)
      r = MATMUL(chg%lat2car,p)
      GetPointDistanceR = mag( r - rt)
      DO i = -nnlayer, nnlayer
        DO j = -nnlayer, nnlayer
          DO k = -nnlayer, nnlayer
            IF (i == 0.AND.j == 0.AND.k == 0) CYCLE
            nnp = p + (/REAL(i * chg%npts(1),q2), REAL(j * chg%npts(2),q2), &
                  REAL(k * chg%npts(3),q2)/)
            nnr = MATMUL(chg%lat2car,nnp)
            distance = mag( nnr - rt)
            GetPointDistanceR = MIN(GetPointDistanceR, distance)
          END DO
        END DO
      END DO
      RETURN
    END FUNCTION GetPointDistanceR

    FUNCTION makewi(chg,nnlayers,vi)
      REAL,DIMENSION(26) :: makewi
      INTEGER, DIMENSION(3,26) :: vi
      INTEGER, DIMENSION(3) :: p1, p2
      INTEGER :: i, nnlayers
      TYPE(charge_obj) :: chg
      p1 = (/chg%npts(1)/2,chg%npts(2)/2,chg%npts(3)/2/)
      DO i  = 1, 26
        p2 = p1 + vi(:,i)
        makewi(i) = 1/(GetPointDistance(p1,p2,chg,nnlayers)**2)
      END DO 
    END FUNCTION

    SUBROUTINE MakeCPRoster(cpr,cptnum,r)
      REAL(q2), DIMENSION(:,:), ALLOCATABLE :: cpr
      INTEGER :: cptnum
      REAL(q2), DIMENSION(3) :: r
      cpr(cptnum,:) = r
    END SUBROUTINE MakeCPRoster


    ! this function interpolates the gradient of a point. weight towards each nn
    ! is determined by its distance to that neighbor. 
    FUNCTION R2GradInterpol(nnInd,r,chg,nnlayers)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3) :: r, R2GradInterpol
      REAL(q2),DIMENSION(:),ALLOCATABLE :: weight
      REAL(q2) :: normalizer,distance
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: nnInd
      INTEGER :: nnlayers, i
      LOGICAL :: onGrid
      normalizer = 0
      onGrid = .FALSE.
      ALLOCATE(weight(SIZE(nnInd)/3))
      DO i = 1,SIZE(nnInd)/3
        distance = GetPointDistanceR( &
          MATMUL(chg%lat2car,r),MATMUL(chg%lat2car,nnInd(i,:)) &
          ,chg,nnlayers)
        IF (distance == 0) THEN
          onGrid = .TRUE.
          EXIT
        END IF
        weight(i) = (1/distance )**2
        normalizer = normalizer + weight(i)
      END DO   
      IF (ongrid) THEN
        R2GradInterpol = CDGrad(nnind(i,:),chg)
      ELSE
        weight = weight / normalizer
        R2GradInterpol = 0
        DO i = 1,SIZE(nnInd)/3
          R2GradInterpol = R2GradInterpol + CDGrad(nnInd(i,:),chg) * weight(i) 
        END DO
      END IF
      DEALLOCATE(weight)
      RETURN
    END FUNCTION R2GradInterpol
   
    FUNCTION R2HesInterpol(nnInd,r,chg,nnLayers)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3,3) :: R2HesInterpol
      REAL(q2),DIMENSION(:),ALLOCATABLE :: weight
      REAL(q2),DIMENSION(3) :: r
      REAL :: normalizer, distance
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: nnInd
      INTEGER :: nnlayers, i
      LOGICAL :: onGrid
      normalizer = 0
      onGrid = .FALSE.
      ALLOCATE(weight(SIZE(nnInd)/3))
      DO i = 1,SIZE(nnInd)/3
        distance = GetPointDistanceR( &
          MATMUL(chg%lat2car,r),MATMUL(chg%lat2car,nnInd(i,:)) &
          ,chg,nnlayers)
        IF (distance == 0) THEN
          onGrid = .TRUE.
          EXIT
        END IF
        weight(i) =  (1/distance)**2
        normalizer = normalizer + weight(i)
      END DO
      IF (onGrid) THEN
        R2HesInterpol = CDHessian(nnInd(i,:),chg)
      ELSE
        weight = weight / normalizer
        R2HesInterpol = 0
        DO i = 1,SIZE(nnInd)/3
          R2HesInterpol = R2HesInterpol + CDHessian(nnInd(i,:),chg) * weight(i)
        END DO
      END IF
      DEALLOCATE(weight)
      RETURN
    END FUNCTION R2HesInterpol

    FUNCTION R2RhoInterpol(nnInd,r,chg,nnLayers)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(:),ALLOCATABLE :: weight
      REAL(q2),DIMENSION(3) :: r
      REAL :: normalizer, R2RhoInterpol, distance
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: nnInd
      INTEGER :: nnlayers, i
      LOGICAL :: onGrid
      normalizer = 0
      onGrid = .FALSE.
      PRINT *, '30 30 30 in cartesian is'
      PRINT *, MATMUL(chg%lat2car,(/30.,30.,30./))
      PRINT *, 'This point is in cartesian'
      PRINT *, MATMUL(chg%lat2car,r)
      PRINT *, 'the nns are'
      ALLOCATE(weight(SIZE(nnInd)/3))
      DO i = 1,SIZE(nnInd)/3
        PRINT *, nnInd(i,:)
        distance = GetPointDistanceR( &
          MATMUL(chg%lat2car,r),MATMUL(chg%lat2car,nnInd(i,:)) &
          ,chg,nnlayers)
        IF (distance == 0) THEN
          onGrid = .TRUE.
          EXIT
        END IF
        weight(i) = (1/distance)**2
        normalizer = normalizer + weight(i)
      END DO
      PRINT *, 'the nn in cartesian are'
      DO i = 1, SIZE(nnInd)/3
        PRINT *, MATMUL(chg%lat2car,nnInd(i,:))
      END DO
      PRINT *, 'the distance are'
      DO i = 1, SIZE(nnInd)/3
        PRINT *, GetPointDistanceR( &
          MATMUL(chg%lat2car,r),MATMUL(chg%lat2car,nnInd(i,:)) &
          ,chg,nnlayers)
      END DO
      PRINT *, 'the weight are'
      DO i = 1, SIZE(nnInd)/3
        PRINT *, weight(i)
      END DO
      IF (onGrid) THEN
        R2RhoInterpol = rho_val(chg,nnind(i,1),nnind(i,2),nnind(i,3))
      ELSE
        weight = weight / normalizer
        R2RhoInterpol = 0
        DO i = 1,SIZE(nnInd)/3
          R2RhoInterpol = R2RhoInterpol + rho_val(chg,nnind(i,1),nnind(i,2),&
            nnind(i,3)) * weight(i)
        END DO
      END IF
      DEALLOCATE(weight)
      RETURN
    END FUNCTION R2RhoInterpol
  
    ! Counts the amount of negative modes in eigenvalues
    FUNCTION CountNegModes(eigvals)
      REAL(q2),DIMENSION(3) :: eigvals
      INTEGER :: CountNegModes,i
      CountNegModes = 0
      DO i = 1,3
        IF (eigvals(i) < 0) CountNegModes = CountNegModes + 1
      END DO
      RETURN
    END FUNCTION CountNegModes
 
    ! updates the count on all types of CPs
    SUBROUTINE UpDateCounts(negCount,maxCount,bondCount,ringCount,cageCount)
      INTEGER :: negCount,maxCount,bondCount,ringCount,cageCount
      IF (negCount == 3) THEN
        maxCount = maxCount + 1
      ELSE IF (negCount == 2) THEN
        bondCount = bondCount + 1
      ELSE IF (negCount == 1) THEN
        ringCount = ringCount + 1
      ELSE 
        cageCount = cageCount + 1
      END IF
    END SUBROUTINE UpDateCounts

    ! This is a subroutine for debugging purpose. It takes all converged
    ! locations of all CP candidates, and outputs them. Bond CP are marked as He,
    ! Ring CP are marked as Ne, cage CP are marked as Ar. Nucleus are written as
    ! normal
    SUBROUTINE VisAllCP (cpcl,cptnum,chg,ions,opts,&
      ringCount,bondCount,cageCount)
      TYPE(charge_obj) :: chg
      TYPE(ions_obj) :: ions
      TYPE(options_obj) :: opts
      TYPE(cpc),DIMENSION(:),ALLOCATABLE :: cpcl
      REAL(q2),DIMENSION(8,3,3) :: nnHes
      REAL(q2),DIMENSION(8,3) :: nnGrad
      REAL(q2),DIMENSION(3,3) :: hessianMatrix,eigvecs
      REAL(q2),DIMENSION(3) :: grad,eigvals,distance, dir
      INTEGER,DIMENSION(8,3) :: nnInd
      INTEGER :: cptnum,i,n1,negCount, j
      INTEGER :: maxCount,bondCount,ringCount,cageCount
      INTEGER :: im,ib,ir,ic
      CHARACTER(LEN=128) :: atoms, natoms
      LOGICAL :: isCart
      OPEN(100,FILE=opts%chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
      IF(opts%in_opt == opts%in_chgcar5) THEN
        DO n1 = 1, 5
          READ(100,'(a)')
        END DO
        READ(100,'(a)') atoms
        READ(100,'(a)') natoms
      ELSE
        ! have to hope that somewhere in the CHGCAR elements are specified.
        READ(100,'(a)') atoms
        DO n1 = 1,4
          READ(100,'(a)')
          ! skipping 1 line of scaling factor and 3 lines of lattice coordinates
        END DO
        READ(100,'(a)') natoms
      END IF
!      PRINT *, 'atoms are'
!      PRINT *, atoms
!      PRINT *, 'natoms are'
!      PRINT *, natoms
      CLOSE(100)

!      DO n1 = 1, cptnum
!        CALL  UpDatecOUNTS(cpcl(n1)%negCount,maxCount,bondCount,ringCount, &
!          cageCount)
!      END DO
!      PRINT *, 'The number of each type CP is'
!      PRINT *, maxCount, bondCount,ringCount,cageCount
!      PRINT *, 'in VisAllCP, cpcl 1 trueind is'
!      PRINT *, cpcl(1)%trueInd
      ! Bond, ring , cage CP are represented with He Ne Ar

      !First get the number of each type of CP

!      DO n1 = 1, cptnum
!        IF (.NOT.cpcl(n1)%isunique) CYCLE
!!        PRINT *, n1
!!        PRINT *, 'trueind is'
!!        PRINT *, cpcl(n1)%trueInd
!        nnInd = SimpleNN(cpcl(n1)%trueInd,chg)
!!        PRINT *, 'nnInd is'
!!        PRINT *, nnInd
!        DO i = 1, 8
!          nnGrad(i,:) = CDGrad(nnInd(i,:),chg)
!          nnHes(i,:,:) = CDHessian(nnInd(i,:),chg)
!        END DO
!        distance = cpcl(n1)%trueind - nnInd(1,:)
!!        PRINT *, 'distance is'
!!        PRINT *, distance
!        grad = trilinear_interpol_grad(nnGrad,distance)
!!        PRINT *, 'grad is'
!!        PRINT *, grad
!        hessianMatrix = trilinear_interpol_hes(nnHes,distance)
!!        PRINT *, 'hessian is'
!!        PRINT *, hessianMatrix
!        CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
!!        PRINT *, 'eigvals are'
!!        PRINT *, eigvals
!        negCount = CountNegModes(eigvals)
!        cpcl(n1)%negCount = negCount
!!        PRINT *, 'negcount is'
!!        PRINT *, negCount
!        CALL UpDateCounts(cpcl(n1)%negCount,maxCount,bondCount,ringCount, & 
!          cageCount)
!      END DO
!      PRINT *, 'updated counts are'
!      PRINT *, maxCount,bondCount,ringCount,cageCount
      ! Write header of allcpPOSCAR
      OPEN(11,FILE='allcpPOSCAR',STATUS='REPLACE',ACTION='WRITE')
      IF (cageCount>0) WRITE(11,'(a)',ADVANCE='NO') 'Ar '
      IF (ringCount>0) WRITE(11,'(a)',ADVANCE='NO') 'Ne '
      IF (bondCount>0) WRITE(11,'(a)',ADVANCE='NO') 'He '
      !IF (maxCount>0) WRITE(11,'(a)',ADVANCE='NO') 'Kr'
!      WRITE(11,'(a)',ADVANCE='NO') ' He  Ne  Ar'
      WRITE(11,'(a)') TRIM(atoms)
      WRITE(11,*) ions%scalefactor
      WRITE(11,*) '     ',ions%lattice(1,:)
      WRITE(11,*) '     ',ions%lattice(2,:)
      WRITE(11,*) '     ',ions%lattice(3,:)
      IF (cageCount>0) WRITE(11,'(a)',ADVANCE='NO') 'Ar '
      IF (ringCount>0) WRITE(11,'(a)',ADVANCE='NO') 'Ne '
      IF (bondCount>0) WRITE(11,'(a)',ADVANCE='NO') 'He '
      !IF (maxCount>0) WRITE(11,'(a)',ADVANCE='NO') ' Kr'
      WRITE(11,'(a)') TRIM(atoms)
      IF (cageCount>0) WRITE(11,'(I)',ADVANCE='NO') cageCount
      IF (ringCount>0) WRITE(11,'(I)',ADVANCE='NO') ringCount
      IF (bondCount>0) WRITE(11,'(I)',ADVANCE='NO') bondCount
      !IF (maxCount>0) WRITE(11,'(I)',ADVANCE='NO') maxCount
      WRITE(11,'(a)') TRIM(natoms)
      WRITE(11,'(a)') 'Cartesian'
      ! Loop three times, each time writes out one type of CP
      DO j = 0, 2
!      PRINT *, 'looking for negcount ' , j 
        DO n1 = 1, cptnum
!          PRINT *, n1, j, cpcl(i)%negCount
          IF (.NOT.cpcl(n1)%isunique) CYCLE
          IF (cpcl(n1)%negCount == j) THEN
!            PRINT *, 'negcount ', j, ' found at '
!            PRINT *, cpcl(n1)%trueInd
!            PRINT *, 'cptnum is', n1
            WRITE(11,*) MATMUL(chg%lat2car,cpcl(n1)%trueR)
            dir(1) = cpcl(n1)%trueR(1)/chg%npts(1)
            dir(2) = cpcl(n1)%trueR(2)/chg%npts(2)
            dir(3) = cpcl(n1)%trueR(3)/chg%npts(3)
          END IF
        END DO
      END DO
      PRINT *, 'writting atomic locations'
      DO j = 1, ions%nions
        WRITE(11,*) ions%r_car(j,:)
      END DO
      CLOSE(11)
    END SUBROUTINE VisAllCP
  
    ! This subroutine tracks steps in up to the past 10 steps. If the next step is
    ! identical to one taken before, it gives the location when repeat is
    ! detected 
    ! it could also give averaged location of the past
    ! 10 steps in ther future, given treatments to PBC. 
    SUBROUTINE DetectCircling(stepCount,rList,temList,trueR,nextTem,averageR)
      REAL(q2),DIMENSION(10,3) :: rList,temList
      REAL(q2),DIMENSION(3) :: trueR,nextTem,averageR
      INTEGER :: stepCount,i,j
      LOGICAL :: isRunningCircles
      ! establish lists if stepCount is low
      isRunningCircles = .FALSE.
      IF ( stepCount <= 10 ) THEN
        rList(stepCount,:) = trueR
        temList(stepCount,:) = nextTem
      ELSE 
        ! update lists
        DO i = 1, 9
          rList(i,:) = rList(i+1,:)
          temList(i,:) = temList(i+1,:)
        END DO
        rList(10,:) = trueR
        temList(10,:) = nextTem
        ! Detect if tem has repeated after 10 steps
        outer: DO i = 1, 10
          DO j = 1, 10
            IF ( j == i ) CYCLE
            IF (ALL(temList(i,:) == temList(j,:),1) ) THEN 
              IF (mag(temList(j,:))== 1.) CYCLE
              isRunningCircles = .TRUE.
!              PRINT *, 'a tem found to be familiar!'
!              PRINT *, 'step ', i, ' has tem'
!              PRINT *, temList(i,:)
!              PRINT *, 'step ', j, ' has tem'
!              PRINT *, temList(j,:)
              EXIT outer
            END IF
          END DO
        END DO outer
        IF (isRunningCircles) THEN
          averageR = 0.
          averageR = trueR
! need to figure out how to deal with PBC first
!          DO i = 1,10
!            averageR = averageR + rList(i,:)
!          END DO
!          averageR = averageR/10
        END IF
      END IF


    END SUBROUTINE DetectCircling

    ! This subroutine looks for critical points in the list that is too close to
    ! another, and averages the same types into one to remove duplicate critical
    ! points. 
    SUBROUTINE ReduceCP(cpl,opts,ucptnum,chg,uBondCount, &
        uRingCount,uCageCount,maxCount,isReduced)
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl,rcpl
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      REAL(q2),DIMENSION(3) :: avgR
      INTEGER :: uCPTNum, uBondCount,uRingCount,uCageCount,maxCount
      INTEGER :: rCPTNum, rBondCount,rRingCount,rCageCount,rmaxCount
      INTEGER :: i,j, nUCPTNum, weight, dupCount
      LOGICAL :: isReduced
      PRINT *, 'Checking for duplicate CP'
      isReduced = .TRUE. 
      dupCount = 0
      rmaxCount = 0
      rRingCount = 0
      rBondCount = 0
      rCageCount = 0
      ! Give the reduced list same length as before, it's ok if a little goes to
      ! waste
      ALLOCATE(rcpl(ucptnum))
      ! start from a critical point. loop through the entire list, look for
      ! another entry 
      nUCPTNum = 0
      DO i = 1, ucptnum
        ! check if this point is already determined as a proxy to some other
        ! point before
        weight = 1
        avgR = cpl(i)%truer
        IF (cpl(i)%hasProxy ) THEN 
          CYCLE
        END IF
        DO j = i, ucptnum
          ! Periodic boundary condition will come to haunt me, periodically. 
          IF (j == i) CYCLE
          IF ( mag(cpl(i)%truer - cpl(j)%truer) <= opts%knob_distance ) THEN
            cpl(j)%hasProxy = .TRUE.
            avgR = (avgR * weight + cpl(j)%truer )/(weight + 1)
            weight = weight + 1
            dupCount = dupCount + 1
            ! The two CP should be the same type!
            IF (cpl(i)%negCount /= cpl(j)%negCount) THEN
              PRINT *,'ERROR: TWO TYPES OF CP ARE TOO CLOSE TO EACH OTHER'
            END IF 
            isReduced = .FALSE.
          END IF
        END DO
        nUCPTNum = nUCPTNum + 1
        ! record the reduced CP
        CALL RecordCPRLight(avgR,chg,rcpl,nUCPTnum, rMaxCount, rRingCount, &
          rBondCount, rCageCount)
      END DO
      CALL ReplaceCPL(cpl,rcpl)
      DO i = 1, SIZE(cpl)
        cpl(i)%isUnique = .TRUE.
      END DO
      maxCount = rMaxCount
      uRingCount = rRingCount
      uBondCount = rBondCount
      uCageCount = rCageCount
      ucptnum = nUCPTnum
      PRINT *, 'The number of duplicate CP found is', dupcount
      DEALLOCATE(rcpl)
    END SUBROUTINE ReduceCP
   
    SUBROUTINE ReplaceCPL(replacee,replacer)
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: replacee,replacer
      INTEGER :: i
      DEALLOCATE (replacee)
      ALLOCATE (replacee(SIZE(replacer)))
      DO i = 1, SIZE(replacer)
        replacee(i) = replacer(i)
      END DO
    END SUBROUTINE ReplaceCPL
 
    SUBROUTINE ResizeCPL(cpl,newSize)
      TYPE(cpc),ALLOCATABLE,DIMENSION(:) ::cpl,newcpl
      INTEGER :: newSize, i, oldSize
      ALLOCATE (newcpl(newSize))
      oldSize = SIZE(cpl)
      DO i = 1, oldSize
        newcpl(i) = cpl(i)
      END DO
      DEALLOCATE (cpl)
      ALLOCATE (cpl(newSize))
      DO i = 1, oldSize
        cpl(i) = newcpl(i)
      END DO
      DEALLOCATE (newCPL)
    END SUBROUTINE ResizeCPL
 
    ! takes in coordinates, gives out interpolated Hessian using central
    ! difference
    FUNCTION CDHessianR(r,chg)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(8,3,3) :: nnhes
      REAL(q2),DIMENSION(3,3) :: CDHessianR
      REAL(q2),DIMENSION(3) :: r, distance
      INTEGER,DIMENSION(8,3) :: nnind
      INTEGER :: i
      nnind = simpleNN(r,chg)
      distance = r - nnind(1,:)
      DO i = 1,8
        nnHes(i,:,:) = CDHessian(nnind(i,:),chg)
      END DO
       !row find the nearest neighbors at this new locaiton
       !first update critical point location
       !the next big step is to interpolate the force at predicted critical
       !point.
      CDHessianR = trilinear_interpol_hes(nnHes,distance)
      RETURN
    END FUNCTION CDHessianR

    FUNCTION CDGradR(r,chg)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(8,3) :: nnGrad
      REAL(q2),DIMENSION(3) :: CDGradR
      REAL(q2),DIMENSION(3) :: r, distance
      INTEGER,DIMENSION(8,3) :: nnind
      INTEGER :: i
      nnind = simpleNN(r,chg)
      distance = r - nnind(1,:)
      DO i = 1,8
        nnGrad(i,:) = CDGrad(nnind(i,:),chg)
      END DO
       !row find the nearest neighbors at this new locaiton
       !first update critical point location
       !the next big step is to interpolate the force at predicted critical
       !point.
      CDGradR = trilinear_interpol_grad(nnGrad,distance)
      RETURN
    END FUNCTION CDGradR
  END MODULE


