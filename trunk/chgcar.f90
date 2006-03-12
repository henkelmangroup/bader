!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!    Module for reading and writing VASP CHGCAR files
!
! By Andri Arnaldson and Graeme Henkelman
! Last modified by 
!-----------------------------------------------------------------------------------!

MODULE chgcar
  USE vars , ONLY : q1,q2
  USE matrix , ONLY : transpose_matrix,matrix_vector
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_charge_chgcar,write_max_rho_chgcar,write_all_bader_chgcar,         &
  &         write_sel_bader_chgcar,write_all_atom_chgcar,write_sel_atom_chgcar
  CONTAINS
!-----------------------------------------------------------------------------------!
! read_charge: Reads the charge density from a file in vasp or Gaussian cube format,
!    by first reading the header, and then charges
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge_chgcar()

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: box,v
    REAL(q2) :: side,vol,t
    INTEGER :: i
    INTEGER,DIMENSION(110) :: elements
    CHARACTER(LEN=7) :: text
    INTEGER :: cr,count_max,t1,t2,nx,ny,nz

    CALL system_clock(t1,cr,count_max)

    vasp=.false.

    OPEN(100,FILE=chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
    READ(100,'(6/,1A7)') text
    REWIND(100)
    IF (text == 'Direct') THEN
      elements=0
      vasp=.true.
      WRITE(*,'(2x,A)') 'VASP-STYLE INPUT FILE'
      READ(100,'(/,1F20.16)') side
      READ(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
      READ(100,'(110I4)') elements
      READ(100,*)
      DO i=1,110
        if(elements(i).eq.0) exit
      ENDDO
      natypes=i-1
      ALLOCATE(num_atom(natypes))
      DO i=1,natypes
        num_atom(i)=elements(i)
      END DO
      ndim=SUM(elements)
      lattice=side*lattice
      CALL transpose_matrix(lattice,B,3,3)
      wdim=ndim
!      ALLOCATE(Rcar(ndim,3),Rdir(ndim,3),voronoi_charge(wdim,4))
! Do not allocate voronoi_charge here
      ALLOCATE(Rcar(ndim,3),Rdir(ndim,3))
      DO i=1,ndim
!   Shouldn't Rdir be multiplied by side?
        READ(100,'(3(2X,1F8.6))') Rdir(i,:)
        CALL matrix_vector(B,Rdir(i,:),v,3,3)
        Rcar(i,:)=v
      END DO
      READ(100,*) 
      READ(100,*) ngxf,ngyf,ngzf
      ALLOCATE(max_rho(ngxf,ngyf,ngzf))
    ELSE
! temp. readin
      WRITE(*,'(1A27)') 'GAUSSIAN-STYLE INPUT FILE'
! Skip the first two lines
      READ(100,*) text
      READ(100,*) text 
      READ(100,*) ndim,box
      corner=box
      wdim=ndim
!      ALLOCATE(Rcar(ndim,3),Rdir(ndim,3),voronoi_charge(wdim,4),nel(ndim))
! Do not allocate voronoi_charge here
      ALLOCATE(Rcar(ndim,3),Rdir(ndim,3),nel(ndim))
      READ(100,*) ngxf,steps(1),t,t
      READ(100,*) ngyf,t,steps(2),t
      READ(100,*) ngzf,t,t,steps(3)
      ALLOCATE(max_rho(ngxf,ngyf,ngzf))
      IF(ngxf<0) ngxf=(-1)*ngxf  ! This should really indicate the units (Bohr/Ang)
      box(1)=REAL((ngxf),q2)*steps(1)
      box(2)=REAL((ngyf),q2)*steps(2)
      box(3)=REAL((ngzf),q2)*steps(3)
      lattice=0.0_q2
      DO i=1,3
        lattice(i,i)=box(i)
      END DO
      vol=volume(lattice)
      DO i=1,3
        lattice(i,i)=lattice(i,i)-steps(i)
      END DO
      CALL transpose_matrix(lattice,B,3,3)
      DO i=1,ndim
        READ(100,*) nel(i),t,Rdir(i,:)
        Rdir(i,:)=(Rdir(i,:)-corner)/(box-steps)
        CALL matrix_vector(B,Rdir(i,:),v,3,3)
        Rcar(i,:)=v
      END DO
    END IF
    nrho=ngxf*ngyf*ngzf
    ALLOCATE(rho(ngxf,ngyf,ngzf))
    IF (vasp) THEN
      READ(100,*) (((rho(nx,ny,nz),nx=1,ngxf),ny=1,ngyf),nz=1,ngzf)
    ELSE
      READ(100,*) (((rho(nx,ny,nz),nz=1,ngzf),ny=1,ngyf),nx=1,ngxf)
      rho=rho*vol 
    END IF
    WRITE(*,'(1A12,1I5,1A2,1I4,1A2,1I4)') 'FFT-grid: ',ngxf,'x',ngyf,'x',ngzf
    WRITE(*,'(2x,A,1A20)') 'CLOSE ... ', chargefile
    CLOSE(100)
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE read_charge
!
!------------------------------------------------------------------------------------!
! output: Write out a summary of the bader analysis.
!         AtomVolumes.dat: Stores the 'significant' Bader volumes associated with
!                          each atom.
!         ACF.dat        : Stores the main output to the screen.
!         BCF.dat        : Stores 'significant' Bader volumes, their coordinates and
!                          charge, atom associated and distance to it. 
!------------------------------------------------------------------------------------!

  SUBROUTINE output()

    REAL(q2) :: sum_achg
    INTEGER :: i,bdimsig,mib,mab,cc,j,nmax
    INTEGER,DIMENSION(bdim) :: rck

    mab=MAXVAL(bader_atom)
    mib=MINVAL(bader_atom)
    OPEN(100,FILE='AtomVolumes.dat',STATUS='replace',ACTION='write')
    WRITE(100,'(A)') '   Atom                     Volume(s)   '
    WRITE(100,'(A,A)') '-----------------------------------------------------------',&
  &                    '-------------'

    DO i=mib,mab
      cc=0
      rck=0
      nmax=0
      DO j=1,bdim
        IF (bader_charge(j,4) > bader_tol) THEN 
          nmax=nmax+1
          IF(bader_atom(j) == i) THEN
            cc=cc+1
            rck(cc)=nmax
          END IF
        END IF
      END DO
      IF (cc == 0) CYCLE
      WRITE(100,'(2X,1I4,2X,A,2X,10000I5)') i,' ... ',rck(1:cc)
    END DO    
    CLOSE(100)

    OPEN(100,FILE='ACF.dat',STATUS='replace',ACTION='write')
    WRITE(*,555) '#','X','Y','Z','VORONOI','BADER','%','MIN DIST'
    WRITE(100,555) '#','X','Y','Z','VORONOI','BADER','%','MIN DIST'
    555 FORMAT(/,4X,1A,9X,1A1,2(11X,1A1),8X,1A7,5X,1A5,9X,1A1,6X,1A10)
    WRITE(*,666) '----------------------------------------------------------------', &
  &              '------------------------------'
    WRITE(100,667) '---------------------------------------------------------------',&
  &              '------------------------------'
    666 FORMAT(1A66,1A26)
    667 FORMAT(1A65,1A27)

    sum_achg=SUM(bader_achg)
    DO i=1,wdim
      WRITE(*,'(1I5,7F12.4)') i,voronoi_charge(i,1:4),bader_achg(i),                 &
  &                           100.*bader_achg(i)/sum_achg,min_dist(i,4)
      WRITE(100,'(1I5,7F12.4)') i,voronoi_charge(i,1:4),bader_achg(i),               &
  &                           100.*bader_achg(i)/sum_achg,min_dist(i,4)
    END DO
    CLOSE(100)

    bdimsig=0
    OPEN(200,FILE='BCF.dat',STATUS='replace',ACTION='write')

    WRITE(200,556) '#','X','Y','Z','CHARGE','NEAREST ATOM','DISTANCE'
    556 FORMAT(/,4X,1A1,13X,1A1,2(16X,1A1),13X,1A7,2X,1A14,2X,1A10)   

    WRITE(200,668) '---------------------------------------------------------------',&
  &                '----------------------------------------'
    668 FORMAT(1A65,1A37)

    DO i=1,bdim
        IF(bader_charge(i,4) > bader_tol) THEN
           bdimsig=bdimsig+1
           WRITE(200,777) bdimsig,bader_charge(i,:),bader_atom(i),bader_dist(i)
           777 FORMAT(1I5,4(5X,F12.4),5X,1I5,5X,1F12.4)
        END IF
    END DO
    CLOSE(200)

    OPEN(300,FILE='dipole.dat',STATUS='replace',ACTION='write')
    WRITE(300,557) '#','X','Y','Z','MAGNITUDE'
    557 FORMAT(/,4X,1A1,10X,1A1,2(15X,1A1),10X,1A10)
    WRITE(300,*) '--------------------------------------------------------------------'
    DO i=1,ndim
      WRITE(300,888) i,dipole(i,:)*4.803_q2,                                         &
  &                  sqrt(DOT_PRODUCT(dipole(i,:),dipole(i,:)))*4.803_q2
!      888 FORMAT(1I5,4ES16.5)
      888 FORMAT(1I5,4F16.6)
    END DO
    CLOSE(300)

    WRITE(*,'(/,2x,A,6X,1I8)')     'NUMBER OF BADER MAXIMA FOUND: ',bdim
    WRITE(*,'(2x,A,6X,1I8)')       '    SIGNIFICANT MAXIMA FOUND: ',bdimsig
    WRITE(*,'(2x,A,2X,1F12.5,/)')  '         NUMBER OF ELECTRONS: ',                 &
  &                                          SUM(bader_charge(1:bdim,4))

  RETURN
  END SUBROUTINE output

!------------------------------------------------------------------------------------!
! write_max_rho: Write out a CHGCAR type file with each entry containing an integer
!    indicating the associated Bader maximum.
!------------------------------------------------------------------------------------!

  SUBROUTINE write_max_rho()

    INTEGER nx,ny,nz

    OPEN(100,FILE='bader_rho.dat',STATUS='replace',ACTION='write')        
    WRITE(*,'(2x,A)') 'WRITING BADER VOLUMES TO BADER_RHO.DAT'
    WRITE(100,'(5E18.11)')                                                           &
  &              (((REAL(max_rho(nx,ny,nz),q2),nx=1,ngxf),ny=1,ngyf),nz=1,ngzf)
    CLOSE(100)

  RETURN
  END SUBROUTINE write_max_rho

!------------------------------------------------------------------------------------!
! Vallvolume: Write out a CHGCAR type file for each of the Bader volumes found.
!------------------------------------------------------------------------------------!

  SUBROUTINE Vallvolume()

    INTEGER :: nx,ny,nz,i,AtomNum,BaderCur,tenths_done,t1,t2,cr,count_max
    CHARACTER(15) :: AtomFileName,AtomNumText
    REAL(q1),dimension(ngxf,ngyf,ngzf) :: rho_tmp

    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    CALL system_clock(t1,cr,count_max)
    AtomNum=0
    tenths_done=0
    DO BaderCur=1,bdim
      DO WHILE ((BaderCur*10/bdim) > tenths_done)
        tenths_done=tenths_done+1
        WRITE(*,'(A,$)') '**'
      ENDDO
      IF(Bader_charge(BaderCur,4)>bader_tol) THEN
        AtomNum=AtomNum+1
        IF(AtomNum<10) THEN
          WRITE(AtomNumText,'(1A3,I1)') '000',AtomNum
        ELSE IF(AtomNum<100) THEN
          WRITE(AtomNumText,'(1A2,I2)') '00',AtomNum
        ELSE IF(AtomNum<1000) THEN
          WRITE(AtomNumText,'(1A1,I3)') '0',AtomNum
        ELSE
          WRITE(AtomNumText,'(I4)') AtomNum
        END IF
        AtomFileName = "Bvol"//Trim(AtomNumText(1:))//".dat"
        OPEN(100,FILE=AtomFileName)
        WRITE(100,*)'Bader volume of Atom',AtomNum
        WRITE(100,*)'1.00'
        WRITE(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
        WRITE(100,'(110I4)') num_atom
        WRITE(100,*)'DIRECT'
        WRITE(100,'(3(2X,1F8.6))') (Rdir(i,:) , i=1,ndim)
        WRITE(100,*)
        WRITE(100,*) ngxf,ngyf,ngzf
        rho_tmp=0.0_q1
        WHERE(max_rho == BaderCur) rho_tmp=rho        
        WRITE(100,'(5E18.11)') (((rho_tmp(nx,ny,nz),nx=1,ngxf),ny=1,ngyf),nz=1,ngzf)
        CLOSE(100)
      END IF
    END DO
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'
    
  RETURN
  END SUBROUTINE Vallvolume

!------------------------------------------------------------------------------------!
! Vatomvolume: Write out a CHGCAR type file for each atom where all Bader volumes
!              assigned to that atom are added together. This is only done if the 
!              atoms has any 'significant' bader volumes associated with it.
!------------------------------------------------------------------------------------!

  SUBROUTINE Vatomvolume()

    INTEGER :: nx,ny,nz,i,j,b,mab,mib,ik,sc,cc,tenths_done,t1,t2,cr,count_max
    INTEGER,DIMENSION(bdim) :: rck
    CHARACTER(15) :: AtomFileName,AtomNumText
    REAL(q1),dimension(ngxf,ngyf,ngzf) :: rho_tmp
 
    CALL system_clock(t1,cr,count_max)
    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES '
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    tenths_done=0
    mab=MAXVAL(bader_atom)
    mib=MINVAL(bader_atom)
    sc=0
    DO ik=mib,mab
      cc=0
      rck=0
      DO j=1,bdim
        IF (bader_charge(j,4) < bader_tol) CYCLE
        IF (bader_atom(j) == ik) THEN
          cc=cc+1
          rck(cc)=j
        END IF                   
      END DO
      sc=sc+cc
      DO WHILE ((ik*10/(mab-mib+1)) > tenths_done)
        tenths_done=tenths_done+1
        WRITE(*,'(A,$)') '**'
      ENDDO
      IF(cc == 0) CYCLE
      IF(ik < 10) THEN
        WRITE(AtomNumText,'(1A3,I1)') '000',ik
      ELSE IF(ik<100) THEN
        WRITE(AtomNumText,'(1A2,I2)') '00',ik
      ELSE IF(ik<1000) THEN
        WRITE(AtomNumText,'(1A1,I3)') '0',ik
      ELSE
        WRITE(AtomNumText,'(I4)') ik
      END IF
      AtomFileName = "BvAt"//Trim(AtomNumText(1:))//".dat"
      OPEN(100,FILE=AtomFileName) 
      WRITE(100,*)'Sum of Bader Volumes'
      WRITE(100,*)'1.00'
      WRITE(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
      WRITE(100,'(110I4)') num_atom
      WRITE(100,*)'DIRECT'
      WRITE(100,'(3(2X,1F8.6))') (Rdir(i,:) , i=1,ndim)
      WRITE(100,*)
      WRITE(100,*) ngxf,ngyf,ngzf
      rho_tmp=0.0_q1
      DO b=1,cc
        WHERE(max_rho == rck(b)) rho_tmp=rho
      END DO
      WRITE(100,'(5E18.11)') (((rho_tmp(nx,ny,nz),nx=1,ngxf),ny=1,ngyf),nz=1,ngzf)
      CLOSE(100)
    END DO
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE Vatomvolume

!------------------------------------------------------------------------------------!
! Vspecvolume: Write out a CHGCAR type file for specified Bader volumes by the user.
!              Volumes associated with a atom can be read from AtomVolumes.dat
!------------------------------------------------------------------------------------!

  SUBROUTINE Vspecvolume()

    INTEGER :: nx,ny,nz,i,b,bdimsig,t1,t2,cr,count_max
    CHARACTER(15) :: AtomFileName
    INTEGER,DIMENSION(bdim,2) :: volsig
    INTEGER,DIMENSION(na) :: vols
    REAL(q1),dimension(ngxf,ngyf,ngzf) :: rho_tmp

    CALL system_clock(t1,cr,count_max)
! Correlate the number for each 'significant' bader volumeto its real number
    bdimsig=0
    volsig=0
    DO i=1,bdim  
      IF (bader_charge(i,4) > bader_tol) THEN
        bdimsig=bdimsig+1
        volsig(bdimsig,1)=bdimsig
        volsig(bdimsig,2)=i
      END IF
    END DO
    vols=volsig(addup,2)
    WRITE(*,'(/,2x,A)') 'WRITING SPECIFIED BADER VOLUMES '
    AtomFileName = "Bvsm.dat"
    OPEN(100,FILE=AtomFileName)
    WRITE(100,*)'Sum of Bader Volumes'
    WRITE(100,*)'1.00'
    WRITE(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
    WRITE(100,'(110I4)') num_atom
    WRITE(100,*)'DIRECT'
    WRITE(100,'(3(2X,1F8.6))') (Rdir(i,:) , i=1,ndim)
    WRITE(100,*)
    WRITE(100,*) ngxf,ngyf,ngzf
    rho_tmp=0.0_q2   
    DO b=1,na
      WHERE(max_rho == vols(b)) rho_tmp=rho
    END DO
    WRITE(100,'(5E18.11)') (((rho_tmp(nx,ny,nz),nx=1,ngxf),ny=1,ngyf),nz=1,ngzf)
    CLOSE(100)
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE Vspecvolume
  
!------------------------------------------------------------------------------------!
! Gallvolume: Write out a Gaussian cube  type file for each of the Bader volumes 
!             found.
!------------------------------------------------------------------------------------!

  SUBROUTINE Gallvolume()

    INTEGER :: nx,ny,nz,i,AtomNum,BaderCur,tenths_done,t1,t2,cr,count_max
    CHARACTER(15) :: AtomFileName,AtomNumText
    INTEGER,DIMENSION(3) :: nxyz
 
    REAL(q1),dimension(ngxf,ngyf,ngzf) :: rho_tmp
    REAL(q2),DIMENSION(3) :: box
    REAL(q2) :: vol

! GH: I added this so that we could write the actual charge
    box(1)=REAL((ngxf),q2)*steps(1)    
    box(2)=REAL((ngyf),q2)*steps(2)
    box(3)=REAL((ngzf),q2)*steps(3)
    lattice=0.0_q2
    DO i=1,3
      lattice(i,i)=box(i)
    END DO
    vol=volume(lattice)
    CALL system_clock(t1,cr,count_max)
    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    nxyz=(/ngxf,ngyf,ngzf/)
            AtomNum=0
    tenths_done=0
    DO BaderCur=1,bdim
      DO WHILE ((BaderCur*10/bdim) > tenths_done)
        tenths_done=tenths_done+1
        WRITE(*,'(A,$)') '**'
      ENDDO          
!      IF ((BaderCur*10/bdim) > tenths_done) THEN
!        tenths_done=(BaderCur*10/bdim)
!        WRITE(*,'(A,$)') '**'
!      END IF          
      IF(Bader_charge(BaderCur,4)>bader_tol) THEN
        AtomNum=AtomNum+1
        IF(AtomNum<10) THEN
          WRITE(AtomNumText,'(1A3,I1)') '000',AtomNum
        ELSE IF(AtomNum<100) THEN
          WRITE(AtomNumText,'(1A2,I2)') '00',AtomNum
        ELSE IF(AtomNum<1000) THEN
          WRITE(AtomNumText,'(1A1,I3)') '0',AtomNum
        ELSE  
          WRITE(AtomNumText,'(I4)') AtomNum
        END IF
        AtomFileName = "Bvol"//Trim(AtomNumText(1:))//".dat"
        OPEN(100,FILE=AtomFileName)       
        WRITE(100,*) 'Gaussian cube file'
        WRITE(100,*) 'Bader volume of Atom',AtomNum
        WRITE(100,'(1I5,3(3X,1F9.6))') ndim,corner
        WRITE(100,'(1I5,3(3X,1F9.6))') ngxf,steps(1),0.000000,0.000000
        WRITE(100,'(1I5,3(3X,1F9.6))') ngyf,0.000000,steps(2),0.000000
        WRITE(100,'(1I5,3(3X,1F9.6))') ngzf,0.000000,0.000000,steps(3)
        DO i=1,ndim
          WRITE(100,'(1I5,3X,1F9.6,3(3X,1F9.6))')                                    & 
  &                  nel(i),REAL(nel(i),q2),steps*REAL((nxyz-1),q2)*Rdir(i,:)+corner
        END DO
        rho_tmp=0.0_q1
        WHERE(max_rho == BaderCur) rho_tmp=rho/vol
! Why do I need to interchange nz and nx in the loops?
        WRITE(100,'(6E13.5)') (((rho_tmp(nx,ny,nz),nz=1,ngzf),ny=1,ngyf),nx=1,ngxf)
        CLOSE(100)
      END IF
    END DO
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'
                
  RETURN
  END SUBROUTINE Gallvolume

!------------------------------------------------------------------------------------!
! Gatomvolume: Write out a Gaussian cube type file for each atom where all Bader 
!              volumes assigned to that atom are added together. This is only done 
!              if the atoms has any 'significant' bader volumes associated with it.
!------------------------------------------------------------------------------------!

  SUBROUTINE Gatomvolume()

    INTEGER :: nx,ny,nz,i,j,ik,mib,mab,tenths_done,sc,cc,b,t1,t2,cr,count_max
    CHARACTER(15) :: volfilename,volnumtext
    INTEGER,DIMENSION(3) :: nxyz
    INTEGER,DIMENSION(bdim) :: rck    
    REAL(q2),DIMENSION(3) :: box
    REAL(q2) :: vol
    REAL(q1),dimension(ngxf,ngyf,ngzf) :: rho_tmp
        
! GH: I added this so that we could write the actual charge
    box(1)=REAL((ngxf),q2)*steps(1)
    box(2)=REAL((ngyf),q2)*steps(2)
    box(3)=REAL((ngzf),q2)*steps(3)
    lattice=0.0_q2
    DO i=1,3   
      lattice(i,i)=box(i)
    END DO
    vol=volume(lattice)  
    CALL system_clock(t1,cr,count_max)
    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES '
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'        
    nxyz=(/ngxf,ngyf,ngzf/)        
    tenths_done=0
    mab=MAXVAL(bader_atom)
    mib=MINVAL(bader_atom)
    sc=0
    DO ik=mib,mab
      DO WHILE ((ik*10/(mab-mib+1)) > tenths_done) 
        tenths_done=tenths_done+1
        WRITE(*,'(A,$)') '**'
      ENDDO
      cc=0
      rck=0
      DO j=1,bdim
        IF (bader_charge(j,4) < bader_tol) CYCLE
        IF (bader_atom(j) == ik) THEN
          cc=cc+1
          rck(cc)=j
        END IF
      END DO
      sc=sc+cc
      IF(cc == 0) CYCLE   
      IF(ik < 10) THEN
        WRITE(volnumtext,'(1A3,I1)') '000',ik
      ELSE IF(ik < 100) THEN
        WRITE(volnumtext,'(1A2,I2)') '00',ik
      ELSE IF(ik < 1000) THEN
        WRITE(volnumtext,'(1A1,I3)') '0',ik
      ELSE   
        WRITE(volnumtext,'(I4)') ik
      END IF
      volfilename = "BvAt"//Trim(volnumtext(1:))//".dat"
      OPEN(100,FILE=volfilename)   
      WRITE(100,*) 'Gaussian cube file'
      WRITE(100,*) 'Bader volume number',ik
      WRITE(100,'(1I5,3(3X,1F9.6))') ndim,corner
      WRITE(100,'(1I5,3(3X,1F9.6))') ngxf,steps(1),0.000000,0.000000
      WRITE(100,'(1I5,3(3X,1F9.6))') ngyf,0.000000,steps(2),0.000000
      WRITE(100,'(1I5,3(3X,1F9.6))') ngzf,0.000000,0.000000,steps(3)
      DO i=1,ndim
        WRITE(100,'(1I5,3X,1F9.6,3(3X,1F9.6))')                                    &
&             nel(i),REAL(nel(i),q2),steps*REAL((nxyz-1),q2)*Rdir(i,:)+corner
      END DO  
      rho_tmp=0.0_q1
      DO b=1,cc
        WHERE(max_rho == rck(b)) rho_tmp=rho/vol  
      END DO
! Why do I need to interchange nz and nx in the loops?
      WRITE(100,'(6E13.5)') (((rho_tmp(nx,ny,nz),nz=1,ngzf),ny=1,ngyf),nx=1,ngxf)
      CLOSE(100)
    END DO
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'
    
  RETURN
  END SUBROUTINE Gatomvolume

!------------------------------------------------------------------------------------!
! Gspecvolume: Write out a Gaussian cube  type file for specified Bader volumes by 
!              the user. Volumes associated with a atom can be read from 
!              AtomVolumes.dat
!------------------------------------------------------------------------------------!

  SUBROUTINE Gspecvolume()

    INTEGER :: nx,ny,nz,i,b,bdimsig,t1,t2,cr,count_max
    CHARACTER(15) :: volfilename
    INTEGER,DIMENSION(3) :: nxyz
    REAL(q2),DIMENSION(3) :: box
    REAL(q2) :: vol
    INTEGER,DIMENSION(bdim,2) :: volsig
    INTEGER,DIMENSION(na) :: vols
    REAL(q1),DIMENSION(ngxf,ngyf,ngzf) :: rho_tmp

    CALL system_clock(t1,cr,count_max)
! Correlate the number for each 'significant' bader volumeto its real number 
    bdimsig=0
    volsig=0
    DO i=1,bdim
      IF (bader_charge(i,4) > bader_tol) THEN
        bdimsig=bdimsig+1
        volsig(bdimsig,1)=bdimsig
        volsig(bdimsig,2)=i
      END IF     
    END DO
    vols=volsig(addup,2)

! GH: I added this so that we could write the actual charge
    box(1)=REAL((ngxf),q2)*steps(1)    
    box(2)=REAL((ngyf),q2)*steps(2)
    box(3)=REAL((ngzf),q2)*steps(3)
    lattice=0.0_q2
    DO i=1,3
      lattice(i,i)=box(i)
    END DO
    vol=volume(lattice)
    WRITE(*,'(/,2x,A)') 'WRITING SPECIFIED BADER VOLUMES '
    nxyz=(/ngxf,ngyf,ngzf/)
!    tenths_done=0
    volfilename='Bvsm.dat'
    OPEN(100,FILE=volfilename)
    WRITE(100,*) 'Gaussian cube file'
    WRITE(100,*) 'Bader volume number',1
    WRITE(100,'(1I5,3(3X,1F9.6))') ndim,corner
    WRITE(100,'(1I5,3(3X,1F9.6))') ngxf,steps(1),0.000000,0.000000
    WRITE(100,'(1I5,3(3X,1F9.6))') ngyf,0.000000,steps(2),0.000000
    WRITE(100,'(1I5,3(3X,1F9.6))') ngzf,0.000000,0.000000,steps(3)
    DO i=1,ndim
        WRITE(100,'(1I5,3X,1F9.6,3(3X,1F9.6))')                                    &
&             nel(i),REAL(nel(i),q2),steps*REAL((nxyz-1),q2)*Rdir(i,:)+corner
    END DO
    rho_tmp=0.0_q2
    DO b=1,na
      WHERE(max_rho == vols(b)) rho_tmp=rho/vol
    END DO
! Why do I need to interchange nz and nx in the loops?
    WRITE(100,'(6E13.5)') (((rho_tmp(nx,ny,nz),nz=1,ngzf),ny=1,ngyf),nx=1,ngxf)
    CLOSE(100)
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'
    
  RETURN
  END SUBROUTINE Gspecvolume

!------------------------------------------------------------------------------------!
! volume: Function returning the triple product of the lattice vectors.
!------------------------------------------------------------------------------------!

  FUNCTION volume(h)
    REAL(q2),INTENT(IN),DIMENSION(3,3) :: h
    REAL(q2) :: volume

    volume = h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2))                                    &
  &       -h(1,2)*(h(2,1)*h(3,3)-h(3,1)*h(2,3))                                      &
  &       +h(1,3)*(h(2,1)*h(3,2)-h(3,1)*h(2,2))

  RETURN
  END FUNCTION volume

!------------------------------------------------------------------------------------!

END MODULE ChargeIOM
