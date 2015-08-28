  !-----------------------------------------------------------------------------------!
  ! Bader charge density analysis program
  !  Module implementing the weight method by Yu and Trinkle [JCP 134, 064111 (2011)]
  !-----------------------------------------------------------------------------------!

  MODULE weight_mod

    USE kind_mod
    USE matrix_mod
    USE bader_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    USE chgcar_mod
    IMPLICIT NONE
    PRIVATE

    TYPE weight_obj
      REAL(q2) :: rho
      INTEGER,DIMENSION(3) :: pos
      INTEGER :: volnum !  similar to bader volnum
      REAL(q2),DIMENSION(:),ALLOCATABLE :: volwgt ! weight to the bader volume
      LOGICAL :: isInterior
    END TYPE
    
    TYPE rvert_obj
      REAL(q2),DIMENSION(3) :: ot3
      REAL(q2) :: angle
    END TYPE
    PUBLIC :: weight_obj
    PUBLIC :: bader_weight_calc,quick_sort

  CONTAINS 

  SUBROUTINE bader_weight_calc(bdr,ions,chgval,opts)

    TYPE(weight_obj),ALLOCATABlE,DIMENSION(:) :: chgList, sortedList
    TYPE(weight_obj) :: tempwobj
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgval,chgtemp
    TYPE(weight_obj) :: wt
    TYPE(ions_obj) :: ions,ionstemp
    TYPE(options_obj) :: opts
    INTEGER(KIND=8) :: totalLength,temp
    INTEGER :: i, j, k, l, n1, n2, n3, walker,nnvect
    REAL(q2) :: tempflx,totalE,tempw,vol
    INTEGER :: t1,t2,cr,cm,Nneighvect
    LOGICAL :: nbb ! never boundry before
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: indList
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: atomList ! stores pos of sig. max. 
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: vect
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: alpha
    nbb=.TRUE.
    tempw=0
    totalE=0
    bdr%nvols=0
    CALL ws_voronoi(ions,nnvect,vect,alpha)
    CALL SYSTEM_CLOCK(t1,cr,cm)
    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
    ELSE
      chgtemp = chgval
    END IF
    totalLength=chgtemp%npts(1)*chgtemp%npts(2)*chgtemp%npts(3)
    ALLOCATE (sortedList(totalLength)) ! sortedList should have the correct size
!    temp=totalLength
    ALLOCATE (chgList(totalLength))
    ALLOCATE (bdr%volnum(chgtemp%npts(1),chgtemp%npts(2),chgtemp%npts(3)))
    ALLOCATE (indlist(chgtemp%npts(1),chgtemp%npts(2),chgtemp%npts(3)))
    bdr%bdim=1
    ALLOCATE(bdr%volpos_lat(bdr%bdim,3)) ! will be expanded as needed
    walker=1  
    DO n1=1, chgtemp%npts(1)
      DO n2=1, chgtemp%npts(2)
        DO n3=1,chgtemp%npts(3)
          chgList(walker)%rho=chgtemp%rho(n1,n2,n3)
          chgList(walker)%pos=(/n1,n2,n3/)
          chgList(walker)%isInterior=.TRUE.
          chgList(walker)%volnum=0
          walker=walker+1
        END DO
      END DO
    END DO
    PRINT *, 'sorting...'
    CALL quick_sort(chgList)
    ! this is a ascending loop. so future loop starts from the end
    ! or flipping it should be faster
    PRINT *,'DONE, recording indicies'
    DO walker=1,totalLength/2
      tempwobj=chgList(totalLength+1-walker)
      chgList(totalLength+1-walker)=chgList(walker)
      chgList(walker)=tempwobj
    END DO
    PRINT *,'DONE.'
    ! firsst loop,deal with all interior points.
    PRINT *,'looking through for interior points'
    i=0
    DO walker=1,totalLength
      indList(chgList(walker)%pos(1), &
              chgList(walker)%pos(2), &
              chgList(walker)%pos(3) &
             )=walker
      CALL interiors(bdr,chgtemp,chgList(walker),walker,indList,nnvect,vect)
      
    END DO
    ALLOCATE (bdr%volchg(bdr%nvols))
    ALLOCATE (bdr%ionvol(bdr%nvols))
    DO walker=1,bdr%nvols
      bdr%volchg(walker)=0
    END DO
    !second loop
    PRINT *,'DONE'
    PRINT *,'Running through all boundry points'
    DO walker=1,totalLength
      CALL boundries(bdr,chgtemp,chgList,walker,ions,indList,nnvect,vect,alpha)
    END DO
    PRINT *,'DONE'
    bdr%volchg = bdr%volchg/REAL(chgval%nrho,q2)
    vol=matrix_volume(ions%lattice)
    vol=vol/chgtemp%nrho
    bdr%ionvol = bdr%ionvol*vol
!    totalE=0
!    DO n1=1,bdr%nvols
!      totalE=totalE+bdr%volchg(n1)
!    END DO
    CALL SYSTEM_CLOCK(t2,cr,cm)
    PRINT *,'Time elapsed:',(t2-t1)/REAL(cr,q2),'seconds'
    ALLOCATE(bdr%nnion(bdr%nvols), bdr%iondist(bdr%nvols), bdr%ionchg(ions%nions))
    ALLOCATE(bdr%volpos_dir(bdr%nvols,3))
    ALLOCATE(bdr%volpos_car(bdr%nvols,3))
    DO i=1,bdr%nvols
      bdr%volpos_dir(i,:) = lat2dir(chgtemp, bdr%volpos_lat(i,:))
      bdr%volpos_car(i,:) = lat2car(chgtemp, bdr%volpos_lat(i,:))
    END DO
    CALL assign_chg2atom(bdr,ions,chgval)
    CALL bader_mindist(bdr,ions,chgtemp)

    ! output part
    
  END SUBROUTINE bader_weight_calc

  !-----------------------------------------------------------------------------------!
  ! interiors: this subroutine loops through all grid points, finding out which
  ! are interiors and which are boundries, and number of basins.
  !-----------------------------------------------------------------------------------!
  SUBROUTINE interiors(bdr,chgtemp,cLW,walker,indList,nnvect,vect) ! chgListWalker
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: indList
    INTEGER,DIMENSION(:,:) :: vect
    INTEGER :: nnvect
    TYPE(weight_obj) :: cLW
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgtemp,chgval
    INTEGER :: length,i,n1,n2,n3,temp,tempnvol,tempvolnum,temp2,j,walker
    REAL(q2) :: temprho,denom,rho,areasum,tempw
    INTEGER :: neivolnum,max_count
    INTEGER,DIMENSION(3) :: xyz,p
    LOGICAL :: ismax
    ismax=.TRUE.
    neivolnum=0
    DO n1=1,nnvect
      p=cLW%pos+(/vect(n1,1),vect(n1,2),vect(n1,3)/)
      CALL pbc(p,chgtemp%npts)
      IF (cLW%rho < chgtemp%rho(p(1),p(2),p(3))) THEN
        ismax=.FALSE.
        IF (bdr%volnum(p(1),p(2),p(3))==0) THEN 
          ! that neighbor must be a boundary point
          cLW%isInterior=.FALSE.
          cLW%volnum=0
          EXIT
          ! this subroutine can now be stopped
        ELSEIF (neivolnum==0) THEN 
          neivolnum=bdr%volnum(p(1),p(2),p(3))
          cLW%volnum=bdr%volnum(p(1),p(2),p(3))
          bdr%volnum(cLW%pos(1),cLW%pos(2),cLW%pos(3))=bdr%volnum(p(1),p(2),p(3))
          ! if this point turnout to be a boundary point, this will be
          ! erased. 
          CYCLE
        ELSEIF (neivolnum/=bdr%volnum(p(1),p(2),p(3))) THEN
          cLW%isInterior=.FALSE.
          cLW%volnum=0
          EXIT
        END IF
      END IF        
    END DO
    ! by the end if it gets no basin assignment, it is boundary.
    IF (ismax .AND. cLW%isInterior) THEN
        ! a maximum. need to test if it is significant
        IF (clW%rho<bdr%tol ) THEN
          ! probably just flucuations in vacuum
          ! assign it to closest real max
          max_count=0
        ELSE
          bdr%nvols=bdr%nvols+1
          bdr%bnum=bdr%nvols
          cLW%volnum=bdr%nvols
          bdr%volnum(cLW%pos(1),cLW%pos(2),cLW%pos(3))=bdr%nvols
          bdr%volpos_lat(bdr%nvols,1) = REAL(cLW%pos(1),q2)
          bdr%volpos_lat(bdr%nvols,2) = REAL(cLW%pos(2),q2)
          bdr%volpos_lat(bdr%nvols,3) = REAL(cLW%pos(3),q2)
          CALL reallocate_volpos(bdr, bdr%nvols+1)
        END IF
    ELSEIF (neivolnum==0) THEN
      cLW%isInterior=.FALSE.
    END IF
    IF (.NOT. cLW%isInterior) THEN
      cLW%volnum=0
      bdr%volnum(cLW%pos(1),cLW%pos(2),cLW%pos(3))=0
    END IF      
  END SUBROUTINE interiors

  ! This subroutine goes through all grid points, adding up rhos of interior
  ! points, calculating flx and weight for boundry points and then adding their
  ! rhos up
  SUBROUTINE boundries(bdr,chgtemp,chgList,walker,ions,indList,nnvect,vect,alpha)
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj):: chgtemp
    TYPE(ions_obj):: ions
    INTEGER         :: walker,i,nnvect
    INTEGER,DIMENSION(:,:) :: vect
    INTEGER,DIMENSION(:,:,:) :: indList
    TYPE(weight_obj),DIMENSION(:) :: chgList
    REAL(q2),DIMENSION(:) :: alpha 
    IF (chgList(walker)%isInterior) THEN
      bdr%volchg(chgList(walker)%volnum)=bdr%volchg(chgList(walker)%volnum)+chgList(walker)%rho
      bdr%ionvol(chgList(walker)%volnum)=bdr%ionvol(chgList(walker)%volnum)+1
    ELSE
      CALL flux_weight(bdr,chgtemp,chgList,walker,ions,indList,nnvect,vect,alpha)
    END IF
  END SUBROUTINE boundries
!
!-------------------------------------------------------------------!
! flux_weight: It takes the current position which is a boundry pooint, looks at
! its neighbors to find flx and weight to all basins.
!-------------------------------------------------------------------!

  SUBROUTINE flux_weight(bdr,chgtemp,chgList,walker,ions,indList,nnvect,vect,alpha)
    TYPE(weight_obj),DIMENSION(:) :: chgList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) ::chgtemp
    TYPE(ions_obj) :: ions
    INTEGER :: n1,n2,n3,n4,walker,nein,counter
    INTEGER,DIMENSION(:,:,:) :: indList
    INTEGER,DIMENSION(:,:) :: vect
    INTEGER :: nnvect
    REAL(q2),DIMENSION(:) :: alpha
    INTEGER,DIMENSION(3) :: p
    REAL(q2) :: denom,area,length,temp
    REAL(q2),DIMENSION(nnvect) :: flx,nom
    denom=0
    counter=0
    temp=0
    ALLOCATE(chgList(walker)%volwgt(bdr%nvols))
    DO n1=1,bdr%nvols
      chgList(walker)%volwgt(n1)=0
    END DO
    n4=0   
    DO n1=1,nnvect
      nom(n1)=0
      flx(n1)=0
      p=chgList(walker)%pos+(/vect(n1,1),vect(n1,2),vect(n1,3)/)
      n4=n4+1
      CALL pbc(p,chgtemp%npts)
      IF (chgtemp%rho(p(1),p(2),p(3))>chgList(walker)%rho) THEN
        nom(n1)=alpha(n1)*(chgtemp%rho(p(1),p(2),p(3))-chgList(walker)%rho) !*length
        denom=denom+nom(n1)
      END IF 
    END DO
    counter=0
    DO n1=1,nnvect
      p=chgList(walker)%pos+(/vect(n1,1),vect(n1,2),vect(n1,3)/)
      CALL pbc(p,chgtemp%npts)
      IF (chgtemp%rho(p(1),p(2),p(3))>chgList(walker)%rho) THEN
        flx(n1)=nom(n1)/denom
        ! interior neighbor and boundry neighbor needs to be treated
        ! separately
        IF (bdr%volnum(p(1),p(2),p(3))/=0) THEN !interior neighbor
          chgList(walker)%volwgt(bdr%volnum(p(1),p(2),p(3)))= & 
            chgList(walker)%volwgt(bdr%volnum(p(1),p(2),p(3)))+ flx(n1)
        ELSE
          ! need to find out that point's position in chgList
          nein=indList(p(1),p(2),p(3))
          DO n4=1,bdr%nvols
            chgList(walker)%volwgt(n4)=chgList(walker)%volwgt(n4)+ &
              chgList(nein)%volwgt(n4)*flx(n1)
          END DO
        END IF
      END IF          
    END DO
    DO n1=1,bdr%nvols
      bdr%volchg(n1)=bdr%volchg(n1)+chgList(walker)%rho*chgList(walker)%volwgt(n1)
      bdr%ionvol(n1)=bdr%ionvol(n1)+chgList(walker)%volwgt(n1)
      temp=temp+chgList(walker)%volwgt(n1)
    END DO
    ! once its all finished, rewrite bader%volnum to enable output generating
    temp=0
    DO n1=1,bdr%nvols
      IF (chgList(walker)%volwgt(n1)>temp) THEN
        temp=chgList(walker)%volwgt(n1)
        bdr%volnum(chgList(walker)%pos(1),chgList(walker)%pos(2),chgList(walker)%pos(3)) &
            = n1
      END IF
    END DO

  END SUBROUTINE flux_weight
 
  ! modified version to sort rvert
  RECURSIVE SUBROUTINE rvert_sort(list,length)
  
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    
    IMPLICIT NONE
    TYPE(rvert_obj), DIMENSION (:), INTENT(IN OUT)  :: list
    ! Local variable
    INTEGER :: i,length
    
    CALL rvert_sort_1(1, length)
    
    CONTAINS
    
      RECURSIVE SUBROUTINE rvert_sort_1(left_end, right_end)
      
        INTEGER, INTENT(IN) :: left_end, right_end
        
        !     Local variables
        INTEGER             :: i, j, itemp
        TYPE(rvert_obj) :: temp
        REAL(q2)                :: reference
        INTEGER, PARAMETER  :: max_simple_sort_size = 6
        
        IF (right_end < left_end + max_simple_sort_size) THEN
          ! Use interchange sort for small lists
          CALL rinterchange_sort(left_end, right_end)
        
        ELSE
          ! Use partition ("quick") sort
          reference = list((left_end + right_end)/2)%angle
          i = left_end - 1; j = right_end + 1
        
          DO
            ! Scan list from left end until element >= reference is found
            DO
              i = i + 1
              IF (list(i)%angle >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
              j = j - 1
              IF (list(j)%angle <= reference) EXIT
            END DO
        
        
            IF (i < j) THEN
              ! Swap two out-of-order elements
              temp = list(i); list(i) = list(j); list(j) = temp
            ELSE IF (i == j) THEN
              i = i + 1
              EXIT
            ELSE
              EXIT
            END IF
          END DO
        
          IF (left_end < j) CALL rvert_sort_1(left_end, j)
          IF (i < right_end) CALL rvert_sort_1(i, right_end)
        END IF
      
      END SUBROUTINE rvert_sort_1
    
    
      SUBROUTINE rinterchange_sort(left_end, right_end)
      
        INTEGER, INTENT(IN) :: left_end, right_end
        
        !     Local variables
        INTEGER             :: i, j, itemp
        TYPE(rvert_obj)    :: temp
        
        DO i = left_end, right_end - 1
          DO j = i+1, right_end
            IF (list(i)%angle > list(j)%angle) THEN
              temp = list(i); list(i) = list(j); list(j) = temp
            END IF
          END DO
        END DO
      
      END SUBROUTINE rinterchange_sort
    
  END SUBROUTINE rvert_sort

  !-----------------------------------------------------------------------------------!
  !Quick_sort: Modified from the original program (info given below)
  !-----------------------------------------------------------------------------------!
  RECURSIVE SUBROUTINE quick_sort(list)
  
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    
    IMPLICIT NONE
    TYPE(weight_obj), DIMENSION (:), INTENT(IN OUT)  :: list
    ! Local variable
    INTEGER(KIND=8) :: i,length
    
    length=SIZE(list)
    CALL quick_sort_1(1, length)
    
    CONTAINS
    
      RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
      
        INTEGER(KIND=8), INTENT(IN) :: left_end, right_end
        
        !     Local variables
        INTEGER(KIND=8)             :: i, j, itemp
        TYPE(weight_obj) :: temp
        REAL(q2)                :: reference
        INTEGER, PARAMETER  :: max_simple_sort_size = 6
        
        IF (right_end < left_end + max_simple_sort_size) THEN
          ! Use interchange sort for small lists
          CALL interchange_sort(left_end, right_end)
        
        ELSE
          ! Use partition ("quick") sort
          reference = list((left_end + right_end)/2)%rho
          i = left_end - 1; j = right_end + 1
        
          DO
            ! Scan list from left end until element >= reference is found
            DO
              i = i + 1
              IF (list(i)%rho >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
              j = j - 1
              IF (list(j)%rho <= reference) EXIT
            END DO
        
        
            IF (i < j) THEN
              ! Swap two out-of-order elements
              temp = list(i); list(i) = list(j); list(j) = temp
            ELSE IF (i == j) THEN
              i = i + 1
              EXIT
            ELSE
              EXIT
            END IF
          END DO
        
          IF (left_end < j) CALL quick_sort_1(left_end, j)
          IF (i < right_end) CALL quick_sort_1(i, right_end)
        END IF
      
      END SUBROUTINE quick_sort_1
    
    
      SUBROUTINE interchange_sort(left_end, right_end)
      
        INTEGER(KIND=8), INTENT(IN) :: left_end, right_end
        
        !     Local variables
        INTEGER(KIND=8)             :: i, j, itemp
        TYPE(weight_obj)    :: temp
        
        DO i = left_end, right_end - 1
          DO j = i+1, right_end
            IF (list(i)%rho > list(j)%rho) THEN
              temp = list(i); list(i) = list(j); list(j) = temp
            END IF
          END DO
        END DO
      
      END SUBROUTINE interchange_sort
    
  END SUBROUTINE quick_sort

!---------------------------
! code from Dallas Trinkle
!---------------------------
  SUBROUTINE ws_voronoi(ions,nnvect,vect,alpha)
    TYPE(ions_obj) :: ions
    TYPE(rvert_obj),ALLOCATABLE,DIMENSION(:) :: rvert
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: R,tR
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: nvect,tnvect  
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: vect
    REAL(q2),DIMENSION(3) :: nv,tempR,R2,rx,ry,tempR2,tempR3
    REAL(q2),DIMENSION(3,3) :: Rdot,Rinv,tempM
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: alph,alpha
    REAL(q2) :: detR,tol,temp,rdRn
    INTEGER :: nv1,nv2,nv3,i,n,j
    INTEGER :: Nrange,Nneigh,maxvert,nnvect,nvert
    INTEGER :: prunN,zeroarea,nvi
    LOGICAL :: ic
    tol=1E-8
    Nrange=3
    Nneigh=(2*Nrange+1)**3-1
    ALLOCATE(nvect(Nneigh,3))
    ALLOCATE(R(Nneigh,5))
    Nneigh=0
    DO nv1 = -Nrange, Nrange
      DO nv2 = -Nrange,Nrange
        DO nv3 = -Nrange,Nrange
          IF (nv1==0 .AND. nv2==0 .AND. nv3==0) THEN
            CYCLE
          ELSE
            Nneigh=Nneigh+1
            nvect(Nneigh,1)=nv1
            nvect(Nneigh,2)=nv2
            nvect(Nneigh,3)=nv3
            nv=(/nv1,nv2,nv3/)
            CALL mult_vect(ions%lattice,nv,tempR)
            DO i=1,3
              R(Nneigh,i)=tempR(i)
            END DO
!            R(Nneigh,1)=temp
!            R(Nneigh,1)=R(Nneigh,2)
!            R(Nneigh,2)=temp
            R(Nneigh,4)=0.5*(tempR(1)**2+tempR(2)**2+tempR(3)**2)
            R(Nneigh,5)=0
          END IF
        END DO
      END DO
    END DO
    DO nv1=1,Nneigh
      DO nv2=1,3
        tempR(nv2)=R(Nneigh-nv1+1,nv2)*0.5
      END DO
      CALL INCELL(tempR,Nneigh,R,ic)
      IF (.NOT. ic) THEN
        R(nv1,5)=1
      END IF
    END DO
    CALL TRIMARRAY(R,nvect,Nneigh)
    ! next step is to find vortex points
    maxvert=(Nneigh-2)*(Nneigh-4)
    ALLOCATE (rvert(maxvert))
    ALLOCATE (alph(Nneigh))
    nnvect=Nneigh
    DO n=1,Nneigh
      nvert=1
      DO nv1=1,3
        Rdot(nv1,1)=R(n,nv1)
      END DO
      R2(1)=R(N,4)
      DO nv1=1,Nneigh
        DO i=1,3
          Rdot(i,2)=R(nv1,i)
        END DO
        R2(2)=R(nv1,4)
        DO nv2=nv1+1,Nneigh
          DO i=1,3
            Rdot(i,3)=R(nv2,i)
          END DO
          R2(3)=R(nv2,4)
          CALL inverse(Rdot,Rinv)
          CALL det(Rdot,detR)
          IF (ABS(detR)>=tol) THEN
            CALL matrix_transpose(Rinv,tempM)
            Rinv=tempM
            tempR=(/0,0,0/)
            CALL mult_vect2(Rinv,R2,tempR)
            DO i=1,3
              rvert(nvert)%ot3(i)=tempR(i)/detR
            END DO
            CALL INCELL(tempR/detR,Nneigh,R,ic)
            IF (ic) nvert=nvert+1
          END IF
        END DO  
      END DO
      zeroarea=0
      DO nvi=1,nvert
        IF (ABS(rvert(nvi)%ot3(1)**2+rvert(nvi)%ot3(2)**2+rvert(nvi)%ot3(3)**2-0.5*R(n,4))<tol) & 
          zeroarea=1
      END DO
      IF (zeroarea ==1 .OR. nvert==0) THEN
        alph(n)=0
        nnvect=nnvect-1
        CYCLE
      END IF
      DO i=1,3
        rx(i)=rvert(1)%ot3(i)
      END DO
      rdRn= (rx(1)*R(n,1)+rx(2)*R(n,2)+rx(3)*R(n,3) ) /  & 
            (R(n,1)**2+R(n,2)**2+R(n,3)**2)
      DO i=1,3
        rx(i)=rx(i)-rdRn*R(n,i)
      END DO
      rdRn=SQRT(rx(1)**2+rx(2)**2+rx(3)**2)
      DO i=1,3
        rx(i)=rx(i)*1/rdRn
        tempR(i)=R(n,i)
      END DO
      CALL CROSS_PRODUCT(tempR,rx,ry)    
      rdRn=SQRT(ry(1)**2+ry(2)**2+ry(3)**2)
      DO i=1,3
        ry(i)=ry(i)*1/rdRn
      END DO
      DO nvi=1,nvert-1
        rvert(nvi)%angle = ATAN2( (rvert(nvi)%ot3(1)*ry(1)+rvert(nvi)%ot3(2)*ry(2)+rvert(nvi)%ot3(3)*ry(3) ), &
                       (rvert(nvi)%ot3(1)*rx(1)+rvert(nvi)%ot3(2)*rx(2)+rvert(nvi)%ot3(3)*rx(3)))
      END DO
      CALL rvert_sort(rvert,nvert-1)
      alph(n)=0
      DO nvi=1,nvert-1
        DO i=1,3
          tempR(i)=rvert(nvi)%ot3(i)
          tempR2(i)=rvert(MOD(nvi,(nvert-1))+1)%ot3(i)
          tempR3(i)=R(n,i)
        END DO
        CALL triple_product(tempR,tempR2,tempR3,temp)
        alph(n)=alph(n)+temp
      END DO
      alph(n)=alph(n)*0.25/R(n,4)
      IF (ABS(alph(n))<tol) THEN
        alph(n)=0
        nnvect=nnvect-1
      END IF
    END DO 
    DO n=1,Nneigh
    END DO
    ALLOCATE (vect(nnvect,3))
    ALLOCATE (alpha(nnvect))
    nvi=1
    DO n=1,Nneigh
      IF (alph(n) /=0 ) THEN
        DO i=1,3
          vect(nvi,i)=nvect(n,i)
        END DO
        alpha(nvi)=alph(n)
        nvi=nvi+1
      END IF
    END DO 
     
  END SUBROUTINE ws_voronoi
 

  ! this subroutine takes a given array and trim its size, delte unwanted
  ! entries.  
  SUBROUTINE trimarray(R,nvect,Nneigh)
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: R,tR
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: nvect,tnvect
    INTEGER :: i,Nneigh,j,nNneigh
    INTEGER :: n1
    nNneigh=0
    DO i=1,Nneigh
      IF ( R(i,5)/=1 )THEN
        nNneigh=nNneigh+1
      END IF
    END DO
    ALLOCATE (tR(nNneigh,5))
    ALLOCATE (tnvect(nNneigh,3))
    j=0
    DO i=1,Nneigh
      IF (R(i,5)==1) THEN
        CYCLE
      ELSE
        j=j+1
        DO n1=1,3
          tR(j,n1)=R(i,n1)
          tnvect(j,n1)=nvect(i,n1)
        END DO 
        tR(j,4)=R(i,4)
      END IF
    END DO
    DEALLOCATE(R)
    DEALLOCATE(nvect)
    ALLOCATE(R(nNneigh,5))
    ALLOCATE(nvect(nNneigh,3))
    DO i=1,nNneigh
      DO j=1,3
        R(i,j)=tR(i,j)
        nvect(i,j)=tnvect(i,j)
      END DO
      R(i,4)=tR(i,4)
      R(i,5)=0

    END DO
    DEALLOCATE(tR)
    DEALLOCATE(tnvect)
    Nneigh=nNneigh
  END SUBROUTINE TRIMARRAY

  SUBROUTINE  INCELL(r,Nneigh,Ra,ic)
    REAL(q2),DIMENSION(:,:) :: Ra
    REAL(q2),DIMENSION(3) :: r
    INTEGER :: Nneigh,i
    LOGICAL :: ic
    REAL(q2) :: tol,dot
    tol=1E-8
    ic=.TRUE.
    DO i=1,Nneigh 
      dot=Ra(i,1)*r(1)+Ra(i,2)*r(2)+Ra(i,3)*r(3)
      IF (dot>Ra(i,4)+tol) THEN 
        ic= .FALSE. 
        RETURN
      END IF
    END DO
     
  END SUBROUTINE INCELL 

  END MODULE
