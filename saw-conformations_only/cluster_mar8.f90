	module cluster
	use vrbls; use mumbers
	implicit none

	public 

! a link to a conformation
	integer, allocatable :: whichsite(:)     ! whichsite(i) is the site # of the link i 

! 'physical' vrbls
	real*8  :: JT    ! Ising coupling constant
	real*8 :: magn   ! \sum of all the spins

! cluster update arrays
	integer, allocatable :: cluster_(:),pocket(:)
	integer :: n_cluster, n_pocket   ! # of site in the cluster and the pocket

!------------ blocking statistics
!	integer, parameter :: b_n_max=1000    ! Max # of blocks
	integer  :: spin_Z_b         ! block size
	integer  :: spin_b_n, spin_i_b    ! filled block nbr. & current nbr. of measurements 
	                                          ! in the (b_n+1)-th block
	real*8 :: spin_Z
	real*8 :: esti_magn, esti_magn2, esti_absmagn
	real*8 :: av_n_cl, i_cl_upd    ! av. cluster size & # of cluster updates

	real*8, allocatable :: magn_stat(:),magn2_stat(:), absmagn_stat(:)



! just a couple of shorthands
	real*8 :: L1           ! L+1
	real*8 :: probb        ! Wolff's probability to occupy a bond

	integer, parameter :: marksite=32         ! add to the charge to mark the sites added to the cluster
!	integer, parameter :: maxcharge=dd*dd+dd  ! max possible charge for an unmarked site
	


! flow control
	real*8 :: cl_step



	contains

!---------------------------------------------------------------------
!----------- main routine for cluster updates given the conformation
!---------------------------------------------------------------------
	subroutine cluster_stuff
	implicit none
	integer :: j, yesorno
	logical :: enough


! init
	call cluster_init
	print*,'cl init ok'

! thermalize
	do j=1, int(1d6)
			!print*,'j = ', j
		!call local_update
		call wolff_update
	enddo

	print*,'therm ok'

!------------- main MC loop  :) --------
 111    continue
	do j=1, int(1d6)

		cl_step = cl_step+1.d0
		!call local_update
		call wolff_update
		call cluster_measure

		call cluster_check
	enddo
	call cluster_prnt(enough)        ! print out & check convergence

	if(.not.enough)goto 111
!--------------------------------------


! cleanup
	call cluster_cleanup

	end subroutine cluster_stuff


!--------------------------------------------
!--- A cluster update [Wolff's]
!--------------------------------------------
	subroutine wolff_update
	implicit none
	integer :: name, site, site0, site1, j

! which one to start from
	name = L1*rndm()+1; if(name>L1)name=L1
	site0=whichsite(name)

! add to the cluster and the pocket 
	n_pocket =1; pocket(n_pocket)=site0 ; 
	n_cluster=1; cluster_(n_cluster)=site0 ; 
	call mark_site(site0)
	

	do while(n_pocket>0)
		site=pocket(n_pocket); n_pocket=n_pocket-1  ! out of the pocket

		do j=1,dd				    ! check the nn
			site1 = ass(j,site)
			if(charge(site1)==0)cycle           
			if(.not.incluster(site1))then
			   if(spin(site1)==spin(site))then  ! can try adding site1 to the cluster
				if( rndm() < probb )then
					n_cluster = n_cluster+1      ! put site1 into the cluster....
					cluster_(n_cluster)=site1
					call mark_site(site1)
					n_pocket=n_pocket+1          ! ...and into the pocket
					pocket(n_pocket)=site1
				endif
			   endif
			endif
		enddo


	enddo

	!print*,'n_cl = ', n_cluster; pause

! update
	do j=1,n_cluster
	   site = cluster_(j) ;  charge(site)=-charge(site)  ! flip the spin
	   call unmark_site(site) ;  ! restore the charges to normal
	enddo
	magn = magn + 2.*spin(site)*n_cluster

	i_cl_upd = i_cl_upd +1.d0
	av_n_cl = av_n_cl + 1.d0*n_cluster


      end subroutine wolff_update




!--------------------------------------------
!--- local update
!--------------------------------------------
      subroutine local_update
	implicit none
	integer :: site0, site1, j, name
	real*8 :: dE, ratio

! which one to update
	name = L1*rndm()+1; if(name>L1)name=L1
	site0=whichsite(name)

! en. change
	dE = 0.d0
	do j=1,dd
	    site1 = ass(j,site0) ;  if(charge(site1)==0)cycle
	    dE = dE + spin(site1)
	enddo

	dE = dE*2.d0*spin(site0)
	ratio = exp(-dE*JT)


! metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif	

	charge(site0)=-charge(site0)
	magn = magn + 2.d0* spin(site0)

      end subroutine local_update




!----------------- return spin(site) given the charge(site) -------
	integer function spin(site)
	implicit none
	integer, intent(in) :: site

	
	if( charge(site)>0 )then; spin=  1 ;
			    else; spin= -1 ;
	endif

	end function spin


!------------  mark the site as a member of the cluster     	
	subroutine mark_site(site)
	implicit none
	integer :: site, s

	s = spin(site)	
	charge(site)=charge(site)+s*marksite   

	end subroutine mark_site


!------------ unmark the site
	subroutine unmark_site(site)
	implicit none
	integer :: site, s

	s = spin(site)	
	charge(site)=charge(site)-s*marksite

	end subroutine unmark_site


!------------ check whether the site is marked
	logical function incluster(site)
	implicit none
	integer :: site
	
	incluster=.false.
	if( abs(charge(site)) > marksite ) incluster=.true.

	end function incluster





!--------------------------------------------
!--- measurements
!--------------------------------------------
      subroutine cluster_measure

!--------- 'naive' statistics
	esti_magn = esti_magn + 1.d0*magn/(L1)
	esti_magn2 = esti_magn2 + (1.d0*magn/(L1))**2
	esti_absmagn = esti_absmagn + abs(1.d0*magn/(L1))	
	spin_Z = spin_Z+1.d0
!----------------------------


!------- blocking stuff

	   spin_i_b = spin_i_b+1

!-------------- is a current block full?
	if( spin_i_b == spin_Z_b )then    ! wrap it up
	    spin_b_n = spin_b_n + 1; spin_i_b=0
	    magn_stat(spin_b_n) = esti_magn; esti_magn = 0.d0
	    magn2_stat(spin_b_n) = esti_magn2; esti_magn2 = 0.d0
	    absmagn_stat(spin_b_n)=esti_absmagn; esti_absmagn=0.d0

		if(spin_b_n == b_n_max)then     ! collate blocks: 1000 -> 500 twice bigger blocks

		    ! check for overflow at Z_b
		  if(spin_Z_b > HUGE(spin_Z_b)/2 )then
	            call prnt
	            call wrt
	            print*; print*; print*
	            print*,'***********************************************'
			  print*,'***********  BLOCK SIZE IS TOO LARGE, exiting'
	            print*,'***********************************************'
	            stop
	          endif

	          call collate( magn_stat(:), spin_b_n )
	          call collate( magn2_stat(:), spin_b_n )
	          call collate( absmagn_stat(:), spin_b_n )
		  spin_b_n = spin_b_n/2; spin_Z_b = spin_Z_b * 2.d0
		endif
	endif


      end subroutine cluster_measure
	

!--------------------------------------------
!--- printout
!--------------------------------------------
      subroutine cluster_prnt(conv_all)
      use ctrls
      implicit none
	logical :: conv, conv_all
	real*8 :: av_o, err_o
	real*8 :: av_m, err_m, av_absm, err_absm, av_m2, err_m2	
	character*99 :: fname


	print*
	print*,'------------------- #sweeps = ',cl_step/L1
	print*,' J/T   = ', JT, ' L1 =  ', L1
!	print*,' <m>   = ', esti_magn/spin_Z
!	print*,' <m^2> = ', esti_magn2/spin_Z
!	print*,' <|m|> = ', esti_absmagn/spin_Z
!	print*
	print*,' av cl size = ', 1.*av_n_cl/i_cl_upd 
	print*,' spin_b_n, _i_b = ', spin_b_n, spin_i_b, spin_Z_b

	if(spin_b_n<4)then
		print*,' too few statistics : spin_b_n =', spin_b_n
		return
	endif

! magn
	print*
	!call mrg_conv(magn_stat(1:spin_b_n),spin_b_n,spin_Z_b,magn_av,magn_err,dummy)
	call mrg_conv(magn_stat(1:spin_b_n),spin_b_n,spin_Z_b,av_o,err_o,conv)

	if(conv)then 
        	print 778, av_o, err_o 
	else
		print 779, av_o, err_o
	endif
 778    format(4x,'<M>   = ',g12.5,2x,' +/- ',g12.5,8x,'   /conv/ ')
 779    format(4x,'<M>   = ',g12.5,2x,' +/- ',g12.5,8x,'   /unconv/ ')

	av_m=av_o ; err_m = err_o
	conv_all = conv	

	  !print*; 
	  !print*,'            doing < M > : '
	  !call mrg(magn_stat(1:spin_b_n),spin_b_n,spin_Z_b)



! abs(magn)
	call mrg_conv(absmagn_stat(1:spin_b_n),spin_b_n,spin_Z_b,av_o, err_o,conv)

	if(conv)then 
        	print 780, av_o, err_o 
	else
		print 781, av_o, err_o
	endif
 780    format(4x,'<|M|> = ',g12.5,2x,' +/- ',g12.5,8x,'   /conv/ ')
 781    format(4x,'<|M|> = ',g12.5,2x,' +/- ',g12.5,8x,'   /unconv/ ')

	av_absm=av_o ; err_absm = err_o
	conv_all = conv.and.conv_all

!	  print*; 
!	  print*,'            doing < |M| > : '
!	  call mrg(absmagn_stat(1:spin_b_n),spin_b_n,spin_Z_b)




! |magn|^2
	call mrg_conv(magn2_stat(1:spin_b_n),spin_b_n,spin_Z_b,av_o,err_o,conv)

	if(conv)then 
        	print 782, av_o, err_o 
	else
		print 783, av_o, err_o
	endif
 782    format(4x,'<M^2> = ',g12.5,2x,' +/- ',g12.5,8x,'   /conv/ ')
 783    format(4x,'<M^2> = ',g12.5,2x,' +/- ',g12.5,8x,'   /unconv/ ')

	av_m2=av_o ; err_m2 = err_o
	conv_all = conv.and.conv_all

!	  print*; 
!	  print*,'            doing < M^2 > : '
!	  call mrg(magn2_stat(1:spin_b_n),spin_b_n,spin_Z_b)


!------------- check convergence
	if(conv_all)then


	! save the results for the current conformation 
	   fname='replicas'//trim(suffix)//'.dat'
	   open(1,file=trim(fname), position='append')
		write(1,*)step
		write(1,*)nn_tot, 1.d0*sum(vect**2)
		write(1,*)av_m, err_m
		write(1,*)av_absm, err_absm
		write(1,*)av_m2, err_m2
		write(1,*)'-----'
	   close(1)


	endif


      end subroutine cluster_prnt



!----------------- cluster_init --------------------
	subroutine cluster_init
	implicit none
	integer :: site, i, dir

! allocate memory
	allocate( whichsite(1:L+1) )
	allocate( pocket(1:L+1), cluster_(1:L+1) )


! give names 
	site=ira; 
	i=1; whichsite(i)=ira ;
	do
	  dir = ch2o(charge(site)) ; site=ass(dir,site)
	  i=i+1; whichsite(i)=site
	  if(site==masha)exit
	enddo

! seed spins : spin(site) = sign( charge(site) )
	magn=0
	do i=1,L+1
	  site=whichsite(i)
	  if(rndm()>0.5)then; charge(site)=-charge(site)
		              magn = magn-1
			else; magn = magn+1
	  endif
	enddo

! init spin statistics
	spin_Z=0.d0
	esti_magn=0.d0; 
	esti_magn2=0.d0 
	esti_absmagn=0.d0

	av_n_cl = 0
	i_cl_upd = 0


	allocate( magn_stat(1:b_n_max) )
	allocate( magn2_stat(1:b_n_max) )
	allocate( absmagn_stat(1:b_n_max) )

	spin_b_n=0;   spin_i_b=0;   spin_Z_b=L+1;  ! something to start with


! shorthands
	L1 = L+1 ! # of sites 
	probb = 1.d0 - exp(-2.d0*JT)

! flow control
	cl_step=0.d0

! update arrays
	cluster_=0; pocket=0


!==========================
!	do i=1,L1
!	   site=whichsite(i)
!	   print*,i,spin(site)
!	enddo 
!	print*,'magn=', magn
!	pause
!============================

	end subroutine cluster_init



!----------------- cluster_cleanup --------------------------
	subroutine cluster_cleanup
	use vrbls
	implicit none
	integer :: i,site

! set all the charges back to >0
	do i=1,L+1
	  site=whichsite(i)
	  charge(site)=abs(charge(site))
	enddo

! deallocate memory
	deallocate( whichsite )
	deallocate( pocket, cluster_ )

        deallocate( magn_stat, magn2_stat, absmagn_stat )

	end subroutine cluster_cleanup


!----------------- cluster_cleanup --------------------------
	subroutine cluster_check
	use vrbls
	implicit none
	integer :: i,site, mmm

! check magn vs. sum over spinz
	mmm=0
	do i=1, int(L1)
	  site=whichsite(i)
	  mmm = mmm+spin(site)
	enddo

	if(mmm/=magn)then
		print*,'cluster_check : mmm, magn = ',mmm,magn
		stop
	endif


	end subroutine cluster_check


	end module cluster
