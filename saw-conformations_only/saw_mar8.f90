!---------------
!--- Common for global variables
!---------------
      module vrbls 
      implicit none

      REAL*8, PARAMETER  :: un1=1.0d0, nul=0.0d0 

! seems like N~400 is still OK on pallas ???

      integer            :: N  !130        ! Number of sites per dimension 
      integer            :: d, dd          ! Dimension 
      integer            :: Nsite          ! Total number of sites
      real*8             :: U, U_in,U_fin  ! interaction & initial for thermalization
      integer            :: L              ! SAW length  

	real*8, allocatable :: eu(:)    ! exp(n*U),  -2d<n<2d
	real*8              :: L23      ! L**(2/3), for the scaling

!--- Sites and Associations
      integer, allocatable :: ass(:,:)
	integer, allocatable :: back(:)
!     Array ass specifies associations between sites:
!     For a site site0, site1=ass(i,site0) yields the nearest
!     site in the direction i. Direction i is understood as the direction
!     of the axis i, if i is less than/equal to d, and as the direction
!     opposite to that of the axis (i-d), if i > d.  

!--- configuration
	integer*1, allocatable :: charge(:)
!     charge = 0 -- no worm on this site, 
!     charge > dd**2 -- masha or ira (direction = charge-dd**2) 
!     charge > 0 -- an inner part of the worm, see ch2io array for directions


	integer :: masha, ira   ! Worm masha and ira sites

! stat. 'keep-track'ers
	integer :: nn_tot       ! total number of n. neighbors
	integer, allocatable  ::  vect(:)  ! masha-to-ira vector
	integer, allocatable  ::  ndir(:)  ! # of bonds in a given direction


!--- Translation tables: charge <--> direction
	integer, allocatable   :: io2ch(:,:)
	integer, allocatable   :: ch2i(:),ch2o(:)

!--- Statistics
	real*8 :: Z, Z_e	 	       ! "Partition function" (normal & 'expensive')  

!------------ blocking statistics
	integer, parameter :: b_n_max=1000    ! Max # of blocks
	integer  :: Z_b         ! block size
	integer  :: b_n, i_b    ! filled block nbr. & current nbr. of measurements 
	                                          ! in the (b_n+1)-th block
	real*8   :: nn_,nn2_    ! # of self contacts
	real*8   :: dist2_      ! end-to-end distance
	real*8, allocatable :: nn_stat(:),nn2_stat(:), dist2_stat(:)

! measurements, 'expensive'
	integer :: Z_b_e      
	integer :: b_n_e, i_b_e


! distribution of the contact # & the end-to-end distance
	real*8, allocatable  ::  dist2_distr(:),nn_distr(:) 
	integer :: dist2_max,nn_distr_max


!------------------ autocorr stat
	integer, parameter :: n_acorr = 10000
	integer            :: i_acorr
	real*8 :: ene_acorr(n_acorr), R2_acorr(n_acorr), grsb_acorr(n_acorr)
	integer :: ac_freq    ! frequency to measure autocorr


	end module vrbls


!----------------------
!---- Common for control parameters
!----------------------
    MODULE ctrls
    use vrbls
    IMPLICIT NONE; SAVE

    CHARACTER*12 :: PAR_VERSION = "conf_only_v0"
    CHARACTER*12 :: version_str

    integer :: n_print
    real*8 :: step_p, step_w, step_m, step_t, step
    real*8 :: i_p, i_w, i_m, i_t

    real*8, dimension(1:5) :: prob    ! Update Manager probabilities

    character*99 :: suffix


! Address/Accept counters      
    REAL*8  :: c_r=nul, a_r=nul       ! for reconnect
    real*8  :: c_m_h, a_m_h, c_m_t, a_m_t ! draw/erase separate counters

    real*8  :: i_rec   ! recalculation counter

    CHARACTER*10 :: status='--------'


!---------- "how much on your watch?" -- "Six o'cloch"
	character*8  :: date
	character*10 :: time
	character*5  :: zone

	integer :: tvalues(8)
      
	real*8 :: time_prev, time_curr, time_elapsed, time_limit
	integer :: hr_curr, hr_prev


!----------- random number seeds
	integer :: r_ij, r_kl

	
    END MODULE ctrls

!------------
!---- Main block
!------------	
	PROGRAM main
	USE vrbls;  use cluster !use mumbers
	IMPLICIT NONE
	   
	CALL INIT  ! initializations
	CALL MC    ! MC looop 

	END PROGRAM main

!---------------------------
!--- Initializations
!---------------------------
	SUBROUTINE Init
	USE vrbls; USE ctrls;  use mumbers; use cluster 
	IMPLICIT NONE 
	integer :: ans1, ans2, i, o,   site, dir
	real*8  :: dp, n_sw
	character*99 :: fname, timestr


!--------------- Set the clock 
	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)
	time_prev=0.d0; time_elapsed = 0.d0; hr_prev = 0
!----------------------------------------------


!--------------- parse the command line : 1nd arg SUFFIX 
!
!                         ! IARGC() doesn't seem to work under Dev.Studio :(
!    if(iargc()/=2)then; 
!        print*,'Expecting TWO args: time(hrs) & SUFFIX'; stop
!    endif

    ! get the suffix
    call getarg(1,suffix)
    fname = 'par'//trim(adjustl(suffix))
    print*,'reading ',trim(adjustl(fname))
!---------------------------------------


! read parameter file
	open(1,file=trim(fname))

    read(1, *) version_str
    if(trim(adjustl(version_str)) /= PAR_VERSION)then
        print*, "par version: expected ", PAR_VERSION, ", got ", version_str
        STOP
    endif
	read(1,*)d;    dd = d+d;    
	read(1,*)N;    Nsite=N**d
	read(1,*)L;    L23 = L**(2.d0/3.d0); 
	            dist2_max=10*L23; nn_distr_max=dd*L
	            ac_freq = floor(0.1* L**2)
	read(1,*)U
	read(1,*)JT
	read(1,*)ans1
	read(1,*)ans2
	read(1,*)U_in ;  U_fin=U  
	read(1,*)n_sw;   step_t = n_sw * L**2   ! thermalization length
	read(1,*)step_p
	read(1,*)step_w
	read(1,*)step_m

	read(1,*)dp; prob(1) = dp            ! move_head (masha)
	             prob(2) = prob(1)+dp    ! move_tail (ira)
	read(1,*)dp; prob(3) = prob(2) + dp  ! reconnect
	if(abs(prob(3)-un1)>1.0d-10)then; print*, 'probabilities..'; stop
	endif

	read(1,*)time_limit;    time_limit=time_limit*3600  ! seconds
	read(1,*)r_ij, r_kl


	close(1)


	allocate( ass(dd,Nsite), back(dd))	
	allocate( charge(Nsite), vect(d), eu(-dd:dd), ndir(d) )
	allocate( io2ch(dd,dd), ch2i(1:dd**2+dd), ch2o(1:dd**2+dd) )

!--- Arrange associations
	call assa


!--- Init the transl. tables ch <--> io

!****************************************************************
! Consider d=2 => dd=4. The table of correspondence ch<-->io then reads:
!
!    \ o! 1 ! 2 ! 3 ! 4 !        Directions:      
!     i\!   !   !   !   !                         2 
!    --------------------                         !
!     1 ! 1 ! 2 ! 3 ! 4 !                     3 --!-- 1   
!    --------------------                         ! 
!     2 ! 5 ! 6 ! 7 ! 8 !                         4  
!    -------------------- 
!      ............. 
!
! ira& masha: their charges are dd**2 + dir !!! 
!
!
!****************************************************************
	do i=1,dd; do o=1,dd; io2ch(i,o)=(i-1)*dd+o
	enddo; enddo

	do i=1,dd**2 
	   ch2i(i)=i/dd; if(mod(i,dd)/=0)ch2i(i)=ch2i(i)+1 
	   ch2o(i)=mod(i,dd); if(ch2o(i)==0)ch2o(i)=dd
	enddo

	do i=1,dd; ch2i(dd**2+i)=i; ch2o(dd**2+i)=i
	enddo

!--- exp(j*U), -2d<j<2d
	do i=-dd,dd; 
	   eu(i)=exp(i*U)
	enddo


!--- Init the random number generator  
	call init_rndm(r_ij, r_kl) !(4836,2748)


!--- Nullifying flow contol counters
      i_w=0; i_p=0; i_m=0; i_t=0; step=0.d0

!--- Nullify address/accept counters
	a_r=nul; c_r=nul
	a_m_h=0.d0; c_m_h=0.d0; a_m_t=0.d0; c_m_t=0.d0
	i_rec=0.d0
	

!--- configuration
	if (ans1 == 1) then
         call read_conf
      else
         call init_conf
      end if


!--- statistics
      if (ans2 == 1) then
         call read_stat
      else
         call init_stat
      end if

      PRINT*,''
      PRINT "(' Init: d = ',I5)",d
      PRINT "('       N = ',I5)",N
      PRINT "('       L = ',I5)",L
      PRINT "('     U/T = ',G11.5)",U
      PRINT*,'#########################################################'


	END SUBROUTINE Init



!---------------------------
!--- MC loop
!---------------------------
	SUBROUTINE MC
	USE ctrls; use vrbls ; use mumbers; use cluster
	IMPLICIT NONE 

	real*8, external :: checktime
	real*8 :: r, pex
	integer :: i, site, tr

! ------------------------------------ Thermalization stage

	if(step_t>0)then

	print*,'thermalization will be done by',step_t/1d6,' mln. steps'


	U=U_in

	DO;

	if(i_t == step_t) exit; i_t = i_t + 1
    i_p=i_p+1; i_w=i_w+1; i_m=i_m+1

	if(mod(floor(i_t),L)==0)then   ! anneal once per L


	      pex=step_t/i_t-1.d0

	      if (pex < 0.05) then
	         U=U_fin
	         do i=-dd,dd;  eu(i) = exp(i*U); enddo

!	       print*,'pex  = ', pex; pause

	      else
	        U = U_fin + (U_in-U_fin)*exp(-1.d0/pex)
	        do i=-dd,dd; eu(i) = exp(i*U); enddo
	      end if

	endif



! update
    r=rndm() 	 
    IF     (r<=prob(1))THEN; CALL move_masha
	else if(r<=prob(2))THEN; CALL move_ira
    else if(r<=prob(3))THEN; CALL reconnect
    ELSE; PRINT*,'Update-> Something wrong: r=',r; STOP
    ENDIF

	if(i_p==step_p)then; i_p=0; CALL prnt;  endif
    if(i_w==step_w)then; i_w=0; call wrt; endif

    ENDDO


	call init_stat
	print*, 'thermalization done by ', i_t/1.d6,' mln. steps..'
	print*, '  U_final = ',log(eu(1))

	endif

    i_p=0; i_w=0; i_m=0



! -------------------------------- Main MC loop :)

	DO;

        step=step+1.d0; i_p=i_p+1; 	i_w=i_w+1;  i_m=i_m+1

        r=rndm() 	 
        IF     (r<=prob(1))THEN; CALL move_masha
	    else if(r<=prob(2))THEN; CALL move_ira
        else if(r<=prob(3))THEN; CALL reconnect
        ELSE; PRINT*,'Update-> Something wrong: r=',r; STOP
        ENDIF

        ! if the size is > N, then recalc.
	    if( any( abs(vect)>=N ) ) call recalc
	
	    call measure_cheap

        !    if(i_m==step_m)then; i_m=0; call measure_expensive; endif
        if(i_p==step_p)then; i_p=0; call prnt; endif
        if(i_w==step_w)then; i_w=0; call wrt; endif

!----------------------   spins

!    ENDDO
!		do;

!	step = step + 1.d0      ! replica #

!! do 10^6 conformational updates
!	do i_m=1, int(1d6)
!	      r=rndm() 	 
!	      if     (r<=prob(1))THEN; CALL move_masha
!	      else if(r<=prob(2))THEN; CALL move_ira
!	      else if(r<=prob(3))THEN; CALL reconnect
!	      else; PRINT*,'Update-> Something wrong: r=',r; STOP
!	      endif
!	enddo

!	print*
!	print*,'====================================================', step
!	print*,' current conformation : nn = ',nn_tot,'  R2 = ',sum(vect**2)
!!	if(nn_tot==0)then
!	  call cluster_stuff
!	  call check ; print*,' post-cluster check ok !'
!!	  pause
!!	endif
!-------------------- spins

			if( checktime() > time_limit )then

	call wrt
	print*
	print*,' *****************The End*******************'
	print*,''
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     V '
	stop
	endif

		enddo    !  main MC loop


	END SUBROUTINE MC


!*******************************************************
!
! Movement:        
!
!      for move_head 1 == tail, 4 == head, and i/o directions are:
!                                           
!           in     out                   For move_tail everyth. is vice-versa
!           <--- --->      
!
!       O--....--3---4           O--....--3---4---NY
!       !                 ===>   !
!       !                        !   
!  1----LA                       LA
! 
!  if $n_i$ is the number of filled nearest neighbours of the site i, then
!
!  $log(old)/U = (n_1 - 1) + (n_2 - 2) + (n_4 - 1)$
!  $log(new)/U = (n_2 - 1) + (n_4 - 2) + (n_N - 1)$ =>
!
!  $ new / old = exp[-U (n_N - n_1) ]$ 
!
!
! ========================== notes : ===============================
!
!    [1]  if 1 is a neighbor of NY, it's spurious as 1 is gonna be erased
!               by the time for the reverse move
!
!*******************************************************
	SUBROUTINE move_masha
	use vrbls; use ctrls ; use mumbers
      implicit none
	logical, external :: bondocc
	logical :: yes
	integer :: j,dir,NY,n_1,n_NY,LA, dir_1,site, site1,site2
	integer :: dir_,site_
	real*8 :: ratio,ff,ww

	status='move_h'; c_m_h=c_m_h+un1; 

! Play a direction of a new link
	dir=dd*rndm()+1; if(dir>dd)dir=dd; 
	NY=ass(dir,masha)
	if(charge(NY)/=0)return              ! SAW 

! Number of neighbours 
	n_1=0; n_NY=0
	do j=1,dd
	    if(charge(ass(j,ira))/=0)n_1=n_1+1
		site1=ass(j,NY)                        ! see note [1]
		if((site1/=ira).and.(charge(site1)/=0))n_NY=n_NY+1
	enddo

	dir_1 = ch2o(charge(ira))
	LA=ass(dir_1,ira)


! ratio
	ratio = eu(n_NY-n_1) 

      if (ratio < 1.d0) then; if (rndm() > ratio) return; end if

! Accepted...changing
	charge(LA)=dd**2+ch2o(charge(LA)); charge(ira)=0 ! LA becomes a new ira

	charge(NY)=dd**2+back(dir)
	charge(masha)=io2ch(charge(masha)-dd**2,dir)          ! NY becomes a new masha


!------- keep-track'ers 
	nn_tot = nn_tot + N_NY - n_1         ! update the nn_tot

! update the ira-masha vector
	if(dir>d)then; vect(back(dir)) = vect(back(dir))-1        ! masha
	               ndir(back(dir)) = ndir(back(dir))+1
	else;   	   vect(dir) = vect(dir)+1
	               ndir(dir) = ndir(dir)+1
	endif

	if(dir_1>d)then; vect(back(dir_1)) = vect(back(dir_1))+1  ! ira
	                 ndir(back(dir_1)) = ndir(back(dir_1))-1
	else;         	 vect(dir_1) = vect(dir_1)-1
	                 ndir(dir_1) = ndir(dir_1)-1
	endif

	masha=NY; ira=LA

	a_m_h=a_m_h+un1


      END SUBROUTINE move_masha





!-------------------
!--- Move ira rather than masha (see above)
!-------------------
	SUBROUTINE move_ira
	use vrbls; use ctrls ; use mumbers
	implicit none 
	logical, external :: bondocc
	logical :: yes
	integer :: j,dir,NY,n_1,n_NY,LA, dir_1,site !, nb, nb1, site1,site2
	integer :: site_, dir_ 
	real*8 :: ratio,ff,ww


	status='move_t'; c_m_t=c_m_t+un1; !print*,status

! Play a direction of a new link
	dir= dd*rndm()+1; if(dir>dd)dir=dd; ! print*,dir
	NY=ass(dir,ira)
	if(charge(NY)/=0)return            ! SAW

! Number of neighbours
	n_1=0; n_NY=0
	do j=1,dd
	    if(charge(ass(j,masha))/=0)n_1=n_1+1
	    site=ass(j,NY)
	    if((site/=masha).and.(charge(site)/=0))n_NY=n_NY+1
	enddo

	dir_1=ch2i(charge(masha))
	LA=ass(dir_1,masha)


! ratio
	ratio = eu(n_NY-n_1) 

      if (ratio < 1.d0) then; if (rndm() > ratio) return; end if

! Accepted...changing
	charge(LA)=dd**2+ch2i(charge(LA)); charge(masha)=0 ! LA becomes a new masha

	charge(NY)=dd**2+back(dir)
	charge(ira)=io2ch(dir,charge(ira)-dd**2)          ! NY becomes a new ira

! update the self contact #
	nn_tot = nn_tot + n_NY - n_1

! update the ira-masha vector
	if(dir>d)then; vect(back(dir)) = vect(back(dir))+1       ! masha
	               ndir(back(dir)) = ndir(back(dir))+1
	else; 	       vect(dir) = vect(dir)-1
	               ndir(dir) = ndir(dir)+1
	endif

	if(dir_1>d)then; vect(back(dir_1)) = vect(back(dir_1))-1 ! ira
			         ndir(back(dir_1)) = ndir(back(dir_1))-1
	else; 	         vect(dir_1) = vect(dir_1)+1
	                 ndir(dir_1) = ndir(dir_1)-1
	endif

! update ira&masha	
	ira=NY; masha=LA

	a_m_t=a_m_t+un1


      END SUBROUTINE move_ira




!*******************************************************
!
! Reconnection:                  
!
!        ...--NY---LA---...---ira            ...--NY   LA--...---ira  
!        .                                    .         !
!        .                                    .         !
!  in /\ .        masha                out /\ .        masha    
!     !! O         !           ===>        !! O         ! 
! out \/ .         !                   in  \/ .         !
!        .         O                          .         O
!        .         !                          .         !
!        ..........!                          ..........!   
!
! 1) ratio is unity
! 2) directions between old masha and NY are to be exchanged
!*******************************************************
      subroutine reconnect
	use vrbls; use ctrls ; use mumbers
      implicit none                
	logical, external :: bondocc
	integer :: site,LA,NY,dir,dir1,d1, nb,nb1,site_,dir_,site1,site2,j
	real*8  :: ratio
	logical :: yes

	status='reco'; c_r=c_r+un1; 

! Check the possibility 
	dir=dd*rndm()+1; if(dir>dd)dir=dd          ! masha --> LA direction
 	LA=ass(dir,masha); if(charge(LA)==0)return  ! no worm at LA ...
	if(LA==ass(ch2o(charge(masha)),masha))return ! LA == O :( 

	dir1=ch2o(charge(LA));   NY=ass(dir1,LA)
	

! Accepted... change the topology

	charge(masha)=io2ch(dir,ch2i(charge(masha)))

	site=ass(ch2o(charge(masha)),masha)  ! site = O in the sketch above

	do; ! exchange directions until NY is reached
	   charge(site)=io2ch(ch2o(charge(site)),ch2i(charge(site))) 
         if(site==NY)exit;
         site=ass(ch2o(charge(site)),site)
	enddo 

	charge(NY)=dd**2+ch2i(charge(NY)); masha=NY           ! NY becomes a new masha

      if(dir<=d)then; d1=d+dir; else; d1=dir-d; endif ! a direction LA --> masha

	if(LA/=ira)then; charge(LA)=io2ch(ch2i(charge(LA)),d1)
	else; charge(LA)=dd**2+d1      ! a special care is required if LA coincides with ira
	endif

! ira-to-masha vector
	if(dir >d)then; vect(dir-d) = vect(dir-d) -1
	                ndir(dir-d) = ndir(dir-d) +1
	else; vect(dir) = vect(dir)+1
	      ndir(dir) = ndir(dir)+1
	endif

	if(dir1>d)then; vect(dir1-d) = vect(dir1-d) -1
	                ndir(dir1-d) = ndir(dir1-d)-1
	else;           vect(dir1) = vect(dir1)+1
	                ndir(dir1) = ndir(dir1)-1
	endif

	a_r=a_r+un1

      end subroutine reconnect





!-------------------
!---- Measurements
!-------------------
      SUBROUTINE measure_cheap
	use vrbls; use ctrls
      implicit none 
	real*8, external :: measure_grassb
	integer :: tmp,j
	real*8  :: tt1
	character*99 :: fname


	tmp =sum(vect**2);

	tt1= 1.d0 - 1.d0*minval(ndir)/maxval(ndir)

	   Z = Z + 1.d0
	   i_b = i_b+1

	   nn_ = nn_ + nn_tot                   ! # of self contacts
	   nn2_ = nn2_ + nn_tot**2
	   dist2_ = dist2_ + 1.d0*tmp           ! masha-to-ira distance squared

! distributions
	   if(tmp<dist2_max)dist2_distr(tmp) = dist2_distr(tmp) + 1.d0
	   if(nn_tot<nn_distr_max)nn_distr(nn_tot)=nn_distr(nn_tot)+1.d0


!-------------- is a current block full?
	if( i_b == Z_b )then    ! wrap it up
	    b_n = b_n + 1; i_b=0
		nn_stat(b_n) = nn_; nn_ = 0.d0
	    nn2_stat(b_n) = nn2_; nn2_ = 0.d0
		dist2_stat(b_n) = dist2_ ; dist2_ = 0.d0  

		if(b_n == b_n_max)then     ! collate blocks: 1000 -> 500 twice bigger blocks

		    ! check for overflow at Z_b
		    if(Z_b > HUGE(Z_b)/2 )then
	          call prnt
	          call wrt
	          print*; print*; print*
	          print*,'***********************************************'
			  print*,'***********  BLOCK SIZE IS TOO LARGE, exiting'
	          print*,'***********************************************'
	          stop
	        endif

	        call collate( nn_stat(:), b_n )
	        call collate( nn2_stat(:), b_n )
	        call collate( dist2_stat(:), b_n )
			b_n = b_n/2; Z_b = Z_b * 2.d0
		endif


	endif


! ----------------------- a/corr stat
	
	

!	if( mod(floor(step),ac_freq)==0  ) then           ! measure once per ac_freq MC steps
!
!	  i_acorr=i_acorr+1
!	  ene_acorr(i_acorr)  = U*nn_tot + W*mw_tot
!	  R2_acorr(i_acorr)   = 1.d0*tmp
!	  grsb_acorr(i_acorr) = tt1
!
!
!	  if(i_acorr==n_acorr)then     ! dump & clear
!
!	     fname='acorr_ene'//trim(suffix)//'.dat'
!	     open(1,file=trim(fname),position='append')
!	      do j=1,n_acorr
!	        write(1,*)ene_acorr(j)
!	      enddo
!	     close(1)
!	     
!	     fname='acorr_R2'//trim(suffix)//'.dat'
!	     open(1,file=trim(fname),position='append')
!	      do j=1,n_acorr
!	        write(1,*)R2_acorr(j)
!	      enddo
!	     close(1)
!           
!	     fname='acorr_grsb'//trim(suffix)//'.dat'
!	     open(1,file=trim(fname),position='append')
!	      do j=1,n_acorr
!	        write(1,*)grsb_acorr(j)
!	      enddo
!	     close(1)
!
!	     i_acorr=0.d0
!	     ene_acorr=0.d0; R2_acorr=0.d0; grsb_acorr=0.d0
!
!	  endif
!
!
!	endif


      END SUBROUTINE measure_cheap


!-------------------
!---- Measurements
!-------------------
      SUBROUTINE measure_expensive
	use vrbls; use ctrls
      implicit none 
	real*8 :: xxx
	integer :: j



      END SUBROUTINE measure_expensive



!------------------------
!--- Init configuration
!------------------------	
	subroutine init_conf
	use vrbls ; use mumbers
	implicit none
	logical, external :: bondocc
!	logical :: yes 
	integer :: i,j,o,oo,dprev, site, next,nb(1:dd),r, k, site1
	integer :: x(1:d), dir, nn,dir_,site_,site2


!--- Nullify charges over the lattice
	charge=0

	if(L<4)then
	  print*,'What''s the point of simulating L = ',L,'  ???? '
	  stop
	endif

!--- Specify a starting configuration at random 
	print*,'set up configuration...'

	ira=3; site=ass(2,ira); charge(ira)=dd**2+2; dprev=d+2  ! set up a first unit

	i=0  !,L-1    ! L is the number of _bonds_ hence the number of _sites_ is L+1
	do; 
	  o=0; nb(:)=0
	  do j=1,dd   ! search for unoccupied site among the neighbours
	     next=ass(j,site) 
		 if(charge(next)==0)then; o=o+1; nb(o)=j; endif
 	  enddo

	  if(o>0)then    ! next link

	     oo=o*rndm()+1; if(oo>o)oo=o; next=ass(nb(oo),site)
	     charge(site)=io2ch(dprev,nb(oo))
	  
    ! set up dprev: a direction next -> site 
          if(nb(oo)<=d)then; dprev=nb(oo)+d; else; dprev=nb(oo)-d; endif 
	
	     site=next
		 i=i+1; if(i>=L-1)exit	

	  else         ! no place to put the next link; need to reconnect
	    
	    print*,'reco ', i
		masha=site; charge(masha)=dd**2+dprev; !	    call dump_struct; 
	    call reco_at_init
	    site=masha; dprev = charge(masha)-dd**2  ! keep drawing from the new masha's position
	    

	  endif
	
	enddo

	masha=site; charge(masha)=dd**2+dprev

	call dump_struct
	print*,'...done'


!-------------- figure the starting vect(:) & nn_tot
	site=ira; x(:)=0; nn=2; ndir=0;
	do;

	dir=ch2o(charge(site));	site=ass(dir,site) 

	! vect(:) & ndir(:)
	if(dir<=d)then;
	   x(dir)=x(dir)+1
	   ndir(dir)=ndir(dir)+1
	else; 
	   x(dir-d)=x(dir-d)-1
	   ndir(dir-d)=ndir(dir-d)+1
	endif

	! nn
	do j=1,dd; if(charge(ass(j,site))/=0)nn=nn+1 ! number of neighbours
	enddo; nn=nn-2;                              ! two trivial neighbours


	  	if(site==masha)exit

	enddo;

	vect=x
	nn_tot = nn/2

	call check
	print*,'init_cnf done'

 
	end subroutine init_conf


!---------------------
! 'Reconnect' @ init
!---------------------
      subroutine reco_at_init
	use vrbls; use ctrls ; use mumbers
      implicit none       
	integer :: ddir         
	integer :: site,LA,NY,dir,dir1,d1, nb,nb1,site_,dir_,site1,site2,j


! choose direction
1111  continue
	dir=dd*rndm()+1; if(dir>dd)dir=dd          ! masha --> LA direction
 	LA=ass(dir,masha); if(charge(LA)==0)goto 1111  ! no worm at LA ...
	if(LA==ass(ch2o(charge(masha)),masha))goto 1111 ! LA == O :( 


! Accepted... change the topology

	charge(masha)=io2ch(dir,ch2i(charge(masha)))

	dir1=ch2o(charge(LA));    NY=ass(dir1,LA)    ! LA -> NY direction
	site=ass(ch2o(charge(masha)),masha)  ! site = O in the sketch above

	do; ! exchange directions until NY is reached
	   charge(site)=io2ch(ch2o(charge(site)),ch2i(charge(site))) 
         if(site==NY)exit;
         site=ass(ch2o(charge(site)),site)
	enddo 

	charge(NY)=dd**2+ch2i(charge(NY)); masha=NY           ! NY becomes a new masha

      if(dir<=d)then; d1=d+dir; else; d1=dir-d; endif ! a direction LA --> masha

	if(LA/=ira)then; charge(LA)=io2ch(ch2i(charge(LA)),d1)
	else; charge(LA)=dd**2+d1      ! a special care is required if LA coincides with ira
	endif

	call dump_struct


      end subroutine reco_at_init



!----------------------
!---- Read configuration 
!----------------------
      SUBROUTINE read_conf
	USE vrbls; USE ctrls
      IMPLICIT NONE
	integer :: site, j
	character*99 :: fname

	fname = 'conf'//trim(adjustl(suffix))//'.dat'

      open(4,file=trim(fname))
	   read(4,*)masha,ira
	   read(4,*)nn_tot
	   read(4,*)vect
	   read(4,*)ndir
	   read(4,*)L
	   charge=0
	   do j=1,L+1
 	     read(4,*)site,charge(site)
	   enddo
      close(4)

	print*,'rd_conf done'

      END SUBROUTINE read_conf


!----------------------
!---- Read statistics	   
!----------------------
      SUBROUTINE Read_Stat
	USE vrbls; USE ctrls
      IMPLICIT NONE
      CHARACTER*5 :: cvoid
	real*8 ::  uu
	integer :: ll
	character*99 :: fname

	fname = 'stat'//trim(adjustl(suffix))//'.dat'

	open(1,file=trim(fname))
	read(1,*)uu, ll
	read(1,*)Z,Z_e

! foolproof
	if(uu/=U)then
	  print*,'rd_stat: U = ',U, ' uu = ',uu
	  stop
	endif

	if(ll/=L)then
	  print*,'rd_stat: L = ',L, ' ll = ',ll
	  stop
	endif



! allocate
	if(allocated(nn_stat))deallocate(nn_stat)
	if(allocated(nn2_stat))deallocate(nn2_stat)
	if(allocated(dist2_stat))deallocate(dist2_stat)

	if(allocated(dist2_distr))deallocate(dist2_distr)
	if(allocated(nn_distr))deallocate(nn_distr)

	allocate( nn_stat(1:b_n_max) )
	allocate( nn2_stat(1:b_n_max) )
	allocate( dist2_stat(1:b_n_max) )


	allocate( dist2_distr(1:dist2_max) )
	allocate( nn_distr(0:nn_distr_max))

	read(1,*)Z_b,   b_n,   i_b 
	read(1,*)Z_b_e, b_n_e, i_b_e 
	read(1,*)nn_,nn_stat
	read(1,*)nn2_,nn2_stat
	read(1,*)dist2_, dist2_stat
	read(1,*)dist2_distr
	read(1,*)nn_distr
	read(1,*)i_acorr     ! a/corr stat
	read(1,*)ene_acorr
	read(1,*)R2_acorr
	read(1,*)grsb_acorr

	close(1)



	print*,'rd_stat done.'

      END SUBROUTINE Read_Stat


!----------------------
!---- initialize statistics	   
!----------------------
      SUBROUTINE init_stat
	USE vrbls; USE ctrls; use cluster
      IMPLICIT NONE
	character*99 :: fname

	Z=0.d0; Z_e = 0.d0

!----------- blocking statistics
	if(allocated(nn_stat))deallocate(nn_stat)
	if(allocated(nn2_stat))deallocate(nn2_stat)
	if(allocated(dist2_stat))deallocate(dist2_stat)

	if(allocated(dist2_distr))deallocate(dist2_distr)
	if(allocated(nn_distr))deallocate(nn_distr)

	allocate( nn_stat(1:b_n_max) )
	allocate( nn2_stat(1:b_n_max) )
	allocate( dist2_stat(1:b_n_max) )

	allocate( dist2_distr(1:dist2_max) )
	allocate( nn_distr(0:nn_distr_max) )

	b_n=0;   i_b=0;   Z_b=L;  ! something to start with
	b_n_e=0; i_b_e=0; Z_b_e=L;  

	nn_=0.d0; nn_stat=0.d0
	nn2_=0.d0; nn2_stat=0.d0
	dist2_=0.d0; dist2_stat=0.d0

	dist2_distr=0.d0
	nn_distr=0.d0

!------------- a/corr stat
	i_acorr=0
	ene_acorr=0.d0; R2_acorr=0.d0; grsb_acorr=0.d0 

!	fname='acorr_ene'//trim(suffix)//'.dat'      ! clear the files if they exist
!	open(1,file=trim(fname))
!	write(1,*)' a/c energy '
!	close(1)
	     
!	fname='acorr_R2'//trim(suffix)//'.dat'
!	open(1,file=trim(fname))
!	write(1,*)' a/c R**2 '
!	close(1)
           
!	fname='acorr_grsb'//trim(suffix)//'.dat'
!	open(1,file=trim(fname))
!	write(1,*)' a/c grsb '
 !     	close(1)

	fname='replicas'//trim(suffix)//'.dat'
	open(1,file=trim(fname))
	write(1,*)U, JT, L, N
	write(1,*)'-----'
	close(1)


      end subroutine init_stat



!----------------------
!---- Write 
!----------------------
	subroutine wrt
	use vrbls; use ctrls
      IMPLICIT NONE
	integer, dimension(1:d) :: x
	integer :: site,dir, j, ccc
	character*99 :: fname
	real*8 :: nrm,thr

	print*,'.............  Writing!!!!...........'

!--- statistics
	fname = 'stat'//trim(adjustl(suffix))//'.dat'
	open(1,file=trim(fname))
	write(1,*)U, L
	write(1,*)Z,Z_e
	write(1,*)Z_b, b_n, i_b
	write(1,*)Z_b_e, b_n_e, i_b_e
	write(1,*)nn_,nn_stat
	write(1,*)nn2_,nn2_stat
	write(1,*)dist2_, dist2_stat
	write(1,*)dist2_distr
	write(1,*)nn_distr
	write(1,*)i_acorr     ! a/corr stat
	write(1,*)ene_acorr
	write(1,*)R2_acorr
	close(1)


!------------------------ distributions
!	if(allocated(dist2_distr))then
!	fname='dist2_distr'//trim(adjustl(suffix))//'.dat'
!	nrm=sum(dist2_distr)
!	thr = maxval(dist2_distr)*1d-10
!	open(1,file=trim(fname))
!	do j=1,dist2_max
!	  if(dist2_distr(j)>thr) write(1,*)j,dist2_distr(j)/nrm
!	enddo
!	close(1)
!	endif


!	fname='grassb_distr'//trim(adjustl(suffix))//'.dat'
!	nrm=sum(grsb_distr)/grassb_n
!	open(1,file=trim(fname))
!	do j=1,grassb_n
!	  write(1,*)1.d0*j/grassb_n,grsb_distr(j)/nrm
!	enddo
!	close(1)

!	fname='nn_distr'//trim(adjustl(suffix))//'.dat'
!	nrm=sum(nn_distr)
!	open(1,file=trim(fname))
!	do j=0,nn_distr_max
!	  write(1,*)j,nn_distr(j)/nrm
!	enddo
!	close(1)

!	fname='mw_distr'//trim(adjustl(suffix))//'.dat'
!	nrm=sum(mw_distr)
!	open(1,file=trim(fname))
!	do j=0,nn_distr_max
!	  write(1,*)j,mw_distr(j)/nrm
!	enddo
!	close(1)
!-----------------------------------------------------


!--- configuration
	fname = 'conf'//trim(adjustl(suffix))//'.dat'
      open(4,file=trim(fname))
	 write(4,*)masha,ira
	 write(4,*)nn_tot
	 write(4,*)vect
	 write(4,*)ndir
	 write(4,*)L
       ccc=0
	 do site=1,Nsite
	   if(charge(site)/=0)then; write(4,*)site,charge(site)
	       ccc = ccc + 1
	   endif
	 enddo
      close(4)
	
!--- structure 
	call dump_struct
!	pause

	print*,'............ Writing done.................'


	end subroutine wrt

!----------------------
!---- Dump the structure
!----------------------
	subroutine dump_struct
	use vrbls; use ctrls;
	implicit none
	integer :: site,dir, j, x(d)
	character*99 :: fname

	
	fname = 'struct'//trim(adjustl(suffix))//'.dat'
	open(2,file=trim(fname))
	write(2,*)'ira = ',ira,'  masha = ',masha
	site=ira; x(:)=0
	do;

	write(2,*)x  , site; 
	  	if(site==masha)exit
	dir=ch2o(charge(site)) !;	print*,x,"  ",site, charge(site)
	site=ass(dir,site) 
	if(dir<=d)then;x(dir)=x(dir)+1
	else; x(dir-d)=x(dir-d)-1
	endif

	enddo;

	close(2)


	end subroutine dump_struct

!----------------------
!---- check 
!----------------------
	subroutine check
	use vrbls; use ctrls
      IMPLICIT NONE
	logical, external :: bondocc
	integer, dimension(1:d) :: x, nd
	integer :: site,dir, nn, j,mmw, zprev
	integer :: dir_,site1,site2

! negative charges
	if( any(charge<0) )then
	        print*
		print*,' OOOPS : negative charges found !!!!'
		stop
	endif


!-------------- walk the walk [from ira to masha]
	site=ira; x(:)=0;  nn=2; zprev=x(d); nd=0
	do;

	do j=1,dd; if(charge(ass(j,site))/=0)nn=nn+1 ! number of neighbours
	enddo; nn=nn-2;                              ! two trivial neighbours

	if(site==masha)exit
	dir=ch2o(charge(site));	site=ass(dir,site) 
!	if(dir<=d)then;x(dir)=x(dir)+1
!	else; x(dir-d)=x(dir-d)-1
!	endif

	if(dir<=d)then;
	   x(dir)=x(dir)+1
	   nd(dir)=nd(dir)+1
	else; 
	   x(dir-d)=x(dir-d)-1
	   nd(dir-d)=nd(dir-d)+1
	endif

	zprev = x(d)

	enddo;
!---------------------------


! double-check nn vs nn_tot
	nn = nn/2
	if(nn/=nn_tot)then; print*,'check: nn /= nn_tot ', status,step
	print*,'nn = ',nn,'   nn_tot =', nn_tot; 
	print*
	call wrt
	
	stop
	endif
 
! nd vs ndir
	if(any(nd/=ndir))then; print*,'check: nd /= ndir ', status,step
	print*,' nd   = ',nd
	print*,' ndir = ',ndir
	print*
	call wrt	
	stop
	endif


	end subroutine check



!----------------------
!---- Print intermediate results 
!----------------------
      SUBROUTINE prnt
      USE vrbls; USE ctrls
      IMPLICIT NONE
	real*8, external :: checktime
	real*8 :: corr
	integer :: j

	integer :: n_ar
	real*8, allocatable :: nn_av(:),nn_err(:), nn2_av(:),nn2_err(:)
	real*8  :: var,err

	if(i_t<step_t)then
	  print*,'-----------------------', checktime()/3600,' hrs'
	  print*,' thermalization: ', 100.d0*i_t/step_t,' %'
	  print*,' U = ', log(eu(1)), '   U_fin = ', U_fin
	  print*,'-----------------------'
	  return
	endif


      if (Z < 1.d0) return

      n_print = n_print+1.d0


!--- Printing out
	print*
      print*,'-------------------------', checktime()/3600,' hrs'
	print*
	print "(' MC step =  ',G9.2,' Z(mln) = ',G9.2,'Z_b/L = ',G9.2)",step/1d6, Z/1.d6, Z_b/L
	print "(' N_lat = ',I4,'  L  = ',I6,'  U = ',G11.5)",N,L,U
	print "(' exp(U) = ',G9.3 )",eu(1)

	if(b_n>4)then

! nn
	  print*; 
	  print*,'            doing nn_: '
	  call mrg(nn_stat(1:b_n),b_n,Z_b)

! nn2
!	  print*; 
!	  print*,'            doing nn2_: '
!	  call mrg(nn2_stat(1:b_n),b_n,Z_b)



!-------------------------------------------------- var_nn
	  print*; 
	  print*,'            doing var_nn_: '

	  n_ar = int( log(1.*b_n)/log(2.) )+1
	  allocate( nn_av(1:n_ar),nn_err(1:n_ar),nn2_av(n_ar),nn2_err(n_ar) )
	  call mrg_save(nn_stat(1:b_n),b_n,Z_b,n_ar,nn_av,nn_err)
	  call mrg_save(nn2_stat(1:b_n),b_n,Z_b,n_ar,nn2_av,nn2_err)

	  do j=1,n_ar
	     var = nn2_av(j)-nn_av(j)**2 
	     err = 2.d0*nn_err(j)*nn_av(j); err = err**2
	     err = sqrt( nn2_err(j)**2 + err )
	     print 888, var/L,err/L ,j
	  enddo
 888    format(4x,g12.5,4x,' +/- ',g12.5,8x,I3)

	  deallocate( nn_av,nn_err,nn2_av,nn2_err )
!----------------------------------------------------------


! dist2
	  print*; 
	  print*,'            doing <R**2>_: '
	  call mrg(dist2_stat(1:b_n),b_n,Z_b)


	else
	  print*
	  print*,'             b_n = ',b_n,i_b 
	  print*
	endif

! address / accept
	print*
	print*,'--- addr / accpt '
	print "(' move masha / ira : ', G9.3,' / ',G9.3 )", a_m_h/(c_m_h+1d-6),a_m_t/(c_m_t+1d-6) 
	print "(' reconnect      : ',G9.3)",a_r/( c_r + 1d-6)

	print "(' recalculate    : ',G9.3)",i_rec/( step + 1d-6)


! overall check	
	call check

!	pause


! check time
	if( checktime() > time_limit )then

	call wrt
	print*
	print*,' *****************The End*******************'
	print*,''
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     | '
	print*,'                     V '
	stop
	endif


      
      end subroutine PRNT



!-----------------
!--- Arranges associations between sites
!-----------------
      subroutine ASSA 
	use vrbls 
      integer site, site1, i, i1, i2 
      integer :: ic(d) 
      
	ic(2:d)=1; ic(1)=0

      DO site=1, Nsite

!------------ Coordinates for site ----------
         i1=1 
         DO
         if (ic(i1) < N) then
            ic(i1)=ic(i1)+1
            DO i2=1,i1-1; ic(i2)=1; END DO
            EXIT
         else; i1=i1+1;  end if 

         END DO
!-------------------------------------------------------


         DO i=1, d
         back(i)=i+d; back(i+d)=i
            if (ic(i) < N) then
               site1=site+N**(i-1)
            else
               site1=site-N**(i-1)*(N-1)
            end if

            ass(i,site)=site1
            ass(d+i,site1)=site

         END DO

      END DO
      
      
      
      end subroutine ASSA



!--------------------------
!  .true. if the bond's occupied
!--------------------------
	logical function bondocc(dir,site)
	use vrbls
	implicit none
	integer :: dir,site
	integer :: c,i,o

	bondocc=.false.

	print*,'bondocc : charge = ', charge(site)

	c = charge(site)
	i=ch2i(c); o = ch2o(c)

	if((i==dir).or.(o==dir))bondocc=.true.

	end function bondocc



!--------------------
!--- recalculate end-to-end vector
!--------------------
	subroutine recalc
	use vrbls; use ctrls
	implicit none
	integer :: site,dir,x(d)

	i_rec = i_rec + 1.d0

!-------------- walk the walk [from ira to masha]
	site=ira; x(:)=0;  
	do;

	if(site==masha)exit
	dir=ch2o(charge(site));	site=ass(dir,site) 
	if(dir<=d)then;x(dir)=x(dir)+1
	else; x(dir-d)=x(dir-d)-1
	endif

	enddo;

	vect=x

	end subroutine recalc





!====================== Statistics management routines =======================

!----------- Collate the array ------------------
      subroutine collate(arr,n)
      real*8, dimension(*) :: arr
      integer :: n, Zb        ! array length & # of measurements per array entry

      integer :: i

      do i=1,n/2
          arr(i)=arr(2*i)+arr(2*i-1)
      enddo

      end subroutine collate

!-------------------------------
!--- Analyze block statistics
!-------------------------------
      subroutine bSTAT(arr,n,Zb,av,err)
      integer              :: n
	real*8               :: Zb
      real*8, dimension(1:n) :: arr
      real*8               :: av, err
      real*8 :: av2, dif
!
!  2 Sept 2007 :   Zb argument is now real*8 to avoid overflow
!
      av  = sum( arr(1:n) )/Zb/n    
      av2 = sum( arr(1:n)**2 )/Zb/Zb/n

	dif=av2-av*av; if(dif<0)dif=0.d0
      err = sqrt( dif ) / sqrt(1.d0*n)                 

      end subroutine bSTAT


!-------------------------------
!--- Merge blocks & emit av +/- err
!-------------------------------
      subroutine mrg(arr,n,Zb)
      integer, intent(in)              :: n, Zb
      real*8, dimension(1:n), intent(in) :: arr

      real*8  :: av, err, arr1(1:n), zb1
      integer :: i, n1

!
!  2 Sept 2007 :   Calls bSTAT(..) w/ real*8 Zb1
!

      zb1 = 1.d0*zb;       arr1(1:n) = arr(1:n); n1=n            

      print*,'-----------' !, int(log(1.d0*n)/log(2.d0))+1

      do;

! emit
        call bSTAT(arr1,n1,zb1,av,err)
        print 777, av, err,n1               
 777    format(4x,g12.5,4x,' +/- ',g12.5,8x,I3)

! enough?
        if(n1<3)exit

! merge
        n1=INT(n1/2); zb1=zb1*2.d0
        do i=1,n1
            arr1(i) =  arr1(2*i-1) + arr1(2*i)
        enddo

      enddo

      print*,'------------'; print*;

      end subroutine mrg



!-------------------------------
!--- Merge blocks & RETURN av +/- err 
!-------------------------------
      subroutine mrg_save(arr,n,Zb,n_o,av_o,err_o)
      integer, intent(in)              :: n, Zb
      real*8, dimension(1:n), intent(in) :: arr

	integer :: n_o
	real*8, dimension(*), intent(out) :: av_o,err_o

      real*8  :: av, err, arr1(1:n), zb1
      integer :: i, n1, jj
!
!  2 Sept 2007 :   Calls bSTAT(..) w/ real*8 Zb1
!
      zb1 = 1.d0*zb;       arr1(1:n) = arr(1:n); n1=n;  jj=1

      do;

! emit
        call bSTAT(arr1,n1,zb1,av,err)
	   av_o(jj) = av; err_o(jj)=err
	   jj=jj+1

! enough?
        if(n1<3)then; 
	   exit
	  endif

! merge
        n1=INT(n1/2); zb1=zb1*2.d0
        do i=1,n1
            arr1(i) =  arr1(2*i-1) + arr1(2*i)
        enddo

      enddo

      end subroutine mrg_save




!-------------------------------
!--- Merge blocks, check convergence & RETURN av +/- err [if conv]
!-------------------------------
	subroutine mrg_conv(arr,n,Zb,av_o,err_o, conv)
	implicit none
	integer, intent(in)              :: n, Zb
	real*8, dimension(1:n), intent(in) :: arr
	real*8, dimension(1:n)  :: av_t,err_t
	real*8 :: av_o,err_o
	logical :: conv

    real*8, parameter :: LIMIT = 0.05    ! convergence limit: if an errorbar fluctuates within 5%, then it has probably converged :)

    real*8  :: av, err, arr1(1:n), zb1
    integer :: i, n1, jj
    real*8 :: prev, dummy
    real*8, allocatable :: delta(:)

      zb1 = 1.d0*zb;       arr1(1:n) = arr(1:n); n1=n;  jj=1

      do;

! emit
        call bSTAT(arr1,n1,zb1,av,err)
	   av_t(jj) = av; err_t(jj)=err
	   jj=jj+1

! enough?
        if(n1<3)then; 
	   exit
	  endif

! merge
        n1=INT(n1/2); zb1=zb1*2.d0
        do i=1,n1
            arr1(i) =  arr1(2*i-1) + arr1(2*i)
        enddo

      enddo


! av and err: av is 2nd last
	jj = jj-1              ! jj is now the # of entry in the blocking array
	av_o = av_t(jj-1)
	err_o = maxval( err_t(1:jj) )

! check convergence
	conv=.false.

	if( jj>4 )then          ! it makes sense checking
		allocate(delta(1:jj-1))
 		prev=err_t(jj)	
		do i=jj-1,1,-1
		    dummy = 0.5*( prev+err_t(i) )
		    delta(i) = (prev-err_t(i))/dummy  
		enddo
		
		do i=jj-3,1,-1
			if( all( (/delta(i),delta(i+1),delta(i+2)/) <LIMIT ) )then
				conv=.true.
				exit
			endif
		enddo

		deallocate(delta)
	endif

      end subroutine mrg_conv




!--------------------------
!  how much on your watch -- six o'cloch
!--------------------------
	real*8 function checktime()
	use ctrls
	real*8 :: dt
			
	time_prev = time_curr; hr_prev = hr_curr

	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	dt = time_curr - time_prev
	if( hr_curr < hr_prev )dt = dt + 24*3600.d0   ! across the midnight

	time_elapsed = time_elapsed + dt

	checktime=time_elapsed	

	end function checktime

