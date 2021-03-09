!-----------------------------
!     Random mumbers
!-----------------------------
	module mumbers
	implicit none

	private

! variables
      DOUBLE PRECISION :: UGEN(97), C, CD, CM
      integer :: ij,kl,I97, J97

! export names
	public :: sranmar, rndm, rn, init_rndm, save_rndm, read_rndm



	contains



!************************************************************************           
       SUBROUTINE SRANMAR(IJ,KL)
C***********************************************************************
C This is the initialization routine for the random number generator RANMAR()
C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
C                                                      0 <= KL <= 30081
C The random number sequences created by these two seeds are of sufficient
C length to complete an entire calculation with. For example, if sveral
C different groups are working on different parts of the same calculation,
C each group could be assigned its own IJ seed. This would leave each group
C with 30000 choices for the second seed. That is to say, this random
C number generator can create 900 million different subsequences -- with
C each subsequence having a length of approximately 10^30.
C
C Use IJ = 1802 & KL = 9373 to test the random number generator. The
C subroutine RANMAR should be used to generate 20000 random numbers.
C Then display the next six random numbers generated multiplied by 4096*4096
C If the random number generator is working properly, the random numbers
C should be:
C             6533892.0  14220222.0  7275067.0
C             6172232.0  8354498.0   10633180.0
C***********************************************************************
c     INTEGER IRAND

      INTEGER :: i7, k7, l7, j7,ii7, jj7, m7 , ij, kl
      DOUBLE PRECISION :: S, T

      IF ( IJ<0 .OR. IJ>31328 .OR. KL<0 .OR. KL>30081) THEN
      PRINT '(A)',' The first seed must be between 0 and 31328'
      PRINT '(A)',' The second seed must be between 0 and 30081'
      STOP;  ENDIF

      I7 = MOD(IJ/177, 177) + 2; J7 = MOD(IJ    , 177) + 2
      K7 = MOD(KL/169, 178) + 1; L7 = MOD(KL,     169)

      DO  ii7 = 1, 97
      S = 0.D0;  T = 0.5D0
         DO  jj7 = 1, 24
            M7 = MOD(MOD(I7*J7, 179)*K7, 179)
            I7 = J7; J7 = K7; K7 = M7; L7 = MOD(53*L7+1, 169)
            IF (MOD(L7*M7, 64) > 32) S = S + T
            T = 0.5D0 * T
         ENDDO
      UGEN(ii7) = S
      ENDDO

      C = 362436.D0 / 16777216.D0
      CD = 7654321.D0 / 16777216.D0
      CM = 16777213.D0 /16777216.D0

      I97 = 97
      J97 = 33; RETURN ;  END subroutine SRANMAR

!**************************************************************************** 
      DOUBLE PRECISION FUNCTION RNDM()
      DOUBLE PRECISION, parameter :: r1=5.d-15, r2=1.d-14     
C***********************************************************************
      DOUBLE PRECISION :: RVAL

      RVAL = UGEN(I97) - UGEN(J97)
      IF ( RVAL < 0.D0 ) RVAL = RVAL + 1.0
      UGEN(I97) = RVAL; I97=I97-1; IF (I97 == 0) I97 =97
      J97=J97-1; IF (J97 == 0) J97=97
      C=C-CD; IF (C < 0.D0 ) C=C+CM
      RVAL = RVAL - C
      IF ( RVAL .LT. 0.D0 ) RVAL = RVAL + 1.0

      RNDM = max(RVAL-r1,r2)
      END FUNCTION RNDM


!-------------------------------------------------------------------------
!     Selecting at random a number RN out of the set {1,2,3, ... ,nmax}.
!-------------------------------------------------------------------------
      integer function RN(nmax)  
      integer :: nmax
      RN=nmax*rndm()+1.d0; if(RN>nmax) RN=nmax
      end function RN

!-----------------------
!    Initializations
!-----------------------
      subroutine init_rndm(a,b)
      integer :: a,b

	ij=a; kl=b; call SRANMAR(ij,kl)

      end subroutine init_rndm


!-----------------------
!    Save the state of the generator
!-----------------------
      subroutine save_rndm(fname)
      character*50 :: fname

	open(23,file=trim(fname))
	 write(23,*)ij,kl,i97,j97
	 write(23,*)c,cd,cm
	 write(23,*)ugen
	close(23)

      end subroutine save_rndm


!-----------------------
!    Read the state of the generator
!-----------------------
      subroutine read_rndm(fname)
      character*50 :: fname

	open(23,file=trim(fname))
	 read(23,*)ij,kl,i97,j97
	 read(23,*)c,cd,cm
	 read(23,*)ugen
	close(23)

      end subroutine read_rndm


	end module
