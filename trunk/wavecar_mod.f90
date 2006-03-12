
	program readwf
	implicit none
! --------------------------------------------------------------------
! Finds the overlap between adiabatic basis states at adjacent time steps
! FINAL OVERLAP IS OFF BY FACTOR OF 2*hbar
! BUT WE TAKE CARE OF THAT IN A LATER SCRIPT
! --------------------------------------------------------------------

! Declare variables
!---------------------------------------------------------------------
	integer     nkpt, nband, nbandmin, nbandmax, ndiff, ic
	real(8)     emax, A(3,3), overlap
	character*5 code
	integer     npw
	real(8)     kpt(3)

	complex(8), allocatable :: coef(:)
	complex(8)  eval
	real(8)     fweight, gweight

	integer i,j,k, iband, ikpt,nxn

        complex(8), allocatable :: ac(:,:)
        complex(8), allocatable :: dc(:,:)

! Set bands we're interested in
!---------------------------------------------------------------------

	nbandmin = 261
	nbandmax = 325

	ndiff=nbandmax-nbandmin

	!????????????????????????????????????????????????????????
	!	if(ndiff.gt.20.or.ndiff.lt.1) stop
	!????????????????????????????????????????????????????????

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HERE'S MY TESTING FILE
!
!	open(unit=19,file="wattest",status="new",form="formatted")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Read in first WAVECAR
!---------------------------------------------------------------------

	open(unit=12,file="WAVECAROLD",status="old",form="unformatted")
	open(unit=13,file="guessit",status="new",form="formatted")

! Read in from WAVECAR:
! the number of k-points
! the number of bands
! the energy maximum
! the cell dimensions
! the type of VASP code

	read(12) nkpt,nband,emax,((A(i,j),i=1,3),j=1,3)

	!????????????????????????????????????????????????????????
	!      write(19,*) 
	!      write(19,*) 'nkpt  =',nkpt
	!      write(19,*) 'nband =',nband
	!      write(19,*) 'emax  =',emax
	!      write(19,*) 'A='
	!      write(19,'(3X,3(1X,f8.3))') (A(i,1),i=1,3)
	!      write(19,'(3X,3(1X,f8.3))') (A(i,2),i=1,3)
	!      write(19,'(3X,3(1X,f8.3))') (A(i,3),i=1,3)
	!????????????????????????????????????????????????????????

	read(12) code


	do ikpt=1, nkpt

! Read in from WAVECAR:
! the number of pw wave fns?
! the k-point? 
      
            read(12) npw, kpt(1:3)


	!????????????????????????????????????????????????????????
	! write(19,*)
 	!write(19,'("k-point #",I3,":  (",3f7.4,")    npw=",I6)') 
	!		& ikpt, (kpt(i),i=1,3),npw 
       	! write(19,*) 'kpt =',kpt
	! write(19,*) 'npw =',npw
	! write(19,*) "  band       energy        weight"
	!
	!????????????????????????????????????????????????????????

            allocate(coef(npw))
            allocate(ac(npw,200))
	    ic = 1 

            do iband = 1, nband
            
! Read in from WAVECAR:
! the energy (eval)
! the occupation (fweight)

                read(12) eval, fweight, (coef(i),i=1,npw)

	!????????????????????????????????????????????????????????
	! write(19,'(5X,I3,5X,f8.4,5x,f8.4)') iband, dreal(eval), fweight
	!????????????????????????????????????????????????????????

! Create matrix (ac):
! row: coefficients for plane waves
! column: coeff of 1st plane wave for bands interested in 

		if(iband.ge.nbandmin.and.iband.le.nbandmax) then
	            do i = 1,npw
	    	        ac(i,ic) = coef(i)
	            enddo
	            ic = ic+1
	        end if

            enddo
            deallocate(coef)

	enddo

	close(unit=12)


! Read second WAVECAR into memory
!---------------------------------------------------------------------

	open(unit=12,file="WAVECARNEW",status="old",form="unformatted")

! Read in from WAVECAR:
! the number of k-points
! the number of bands
! the energy maximum
! the cell dimensions
! the type of VASP code

        read(12) nkpt,nband,emax,((A(i,j),i=1,3),j=1,3)

        !????????????????????????????????????????????????????????
        !      write(19,*)
        !      write(19,*) 'nkpt  =',nkpt
        !      write(19,*) 'nband =',nband
        !      write(19,*) 'emax  =',emax
        !      write(19,*) 'A='
        !      write(19,'(3X,3(1X,f8.3))') (A(i,1),i=1,3)
        !      write(19,'(3X,3(1X,f8.3))') (A(i,2),i=1,3)
        !      write(19,'(3X,3(1X,f8.3))') (A(i,3),i=1,3)
        !????????????????????????????????????????????????????????

        read(12) code


        do ikpt=1, nkpt

! Read in from WAVECAR:
! the number of pw wave fns?
! the k-point?

            read(12) npw, kpt(1:3)

        !????????????????????????????????????????????????????????
        ! write(19,*)
        !write(19,'("k-point #",I3,":  (",3f7.4,")    npw=",I6)')
        !               & ikpt, (kpt(i),i=1,3),npw
        ! write(19,*) 'kpt =',kpt
        ! write(19,*) 'npw =',npw
        ! write(19,*) "  band       energy        weight"
	!
        !????????????????????????????????????????????????????????

            allocate(coef(npw))
            allocate(dc(npw,200))
            ic = 1

            do iband = 1, nband

! Read in from WAVECAR:
! the energy (eval)
! the occupation

                read(12) eval, fweight, (coef(i),i=1,npw)

        !????????????????????????????????????????????????????????
        ! write(19,'(5X,I3,5X,f8.4,5x,f8.4)') iband, dreal(eval), fweight
        !????????????????????????????????????????????????????????

! Create matrix (dc):
! row: coefficients for plane waves
! column: coeff of 1st plane wave for bands interested in

		if(iband.ge.nbandmin.and.iband.le.nbandmax) then
		    do i = 1,npw
			dc(i,ic) = coef(i)
		    enddo
		    ic = ic+1
		end if

	    enddo
	    deallocate(coef)

	enddo

	close(unit=12)

! THESE OVERLAPS SHOULD BE DIVIDED BY 2

!---------------------------------------------------------------------

! Overlap of "old wavecar" bands with other "old wavecar" bands
! matrix nxn by nxn:
! row1: band1 overlap with bands bandbegin through bandend

	    write(13,*)"old wavecar"
	nxn=ndiff+1
	do i = 1,nxn
	  do j = 1,nxn
	    overlap = 0.
	    do k = 1,npw
	      overlap=overlap+conjg(ac(k,i))*ac(k,j)
	    end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
	  end do
	    write(13,*)" "
	end do
!--------------------------------------------------------

! Overlap of "new wavecar" bands with other "old wavecar" bands
! matrix nxn by nxn:
! row1: band1 overlap with bands bandbegin through bandend

	    write(13,*)"new wavecar"
	nxn=ndiff+1
	 do i = 1,nxn
	  do j = 1,nxn
	    overlap = 0.
	    do k = 1,npw
	      overlap=overlap+conjg(dc(k,i))*dc(k,j)
	    end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
	  end do
	    write(13,*)" "
	end do
!--------------------------------------------------------

! Overlap of "old wavecar" bands with other "old wavecar" bands
! matrix nxn by nxn:
! row1: band1 overlap with bands bandbegin through bandend

	    write(13,*)"mix wavecars"
	nxn=ndiff+1
	do i = 1,nxn
	  do j = 1,nxn
	    overlap = 0.
	    do k = 1,npw
	      overlap=overlap+conjg(ac(k,i))*dc(k,j)
	    end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
	  end do
	    write(13,*)" "
	end do
!--------------------------------------------------------

! Calculate overlap we actually use (real### files)
! matrix nxn by nxn:

	    write(13,*)"first term"
	nxn=ndiff+1
	do i = 1,nxn
	  do j = 1,nxn
	    overlap = 0.
	    do k = 1,npw
	      overlap=overlap+conjg(ac(k,i))*dc(k,j) 
	      overlap=overlap-conjg(dc(k,i))*ac(k,j) 
	    end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
	  end do
	    write(13,*)" "
	end do
!--------------------------------------------------------

! Yet more overlap
! matrix nxn by nxn:

	    write(13,*)"second term"
	nxn=ndiff+1
	do i = 1,nxn
	  do j = 1,nxn
	    overlap = 0.
	    do k = 1,npw
	      overlap=overlap+conjg(dc(k,i))*dc(k,j) 
	      overlap=overlap-conjg(ac(k,i))*ac(k,j) 
	    end do
            write(13,'(E25.18)', ADVANCE='NO') overlap
	  end do
	    write(13,*)" "
	end do
!--------------------------------------------------------
      end


