!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----Utility subroutines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine assert(condition, message)
      logical, parameter :: DEBUG=.true. 
      logical, intent(in) :: condition
      character(*), intent(in) :: message
      if (.not. DEBUG) then
         return
      else if (.not. condition) then
      	write(*,*) "Assertion failed: ", message
      	stop
      endif
      return
      end subroutine assert

!-----------------------------------------------------------------------
      subroutine invariantJ3(s11, s22, s33, s12, s23, s31, J3)
      implicit none
      real*8 s11, s22, s33, s12, s23, s31, J3
      J3 = 2 * s11**3 / 27
     +   - s11**2 * s22 / 9
     +   - s11**2 * s33 / 9
     +   + s11 * s12**2 / 3
     +   - s11 * s22**2 / 9
     +   + 4 * s11 * s22 * s33 / 9
     +   - 2 * s11 * s23**2 / 3
     +   + s11 * s31**2 / 3
     +   - s11 * s33**2 / 9
     +   + s12**2 * s22 / 3
     +   - 2 * s12**2 * s33 / 3
     +   + 2 * s12 * s23 * s31
     +   + 2 * s22**3 / 27
     +   - s22**2 * s33 / 9
     +   + s22 * s23**2 / 3
     +   - 2 * s22 * s31**2 / 3
     +   - s22 * s33**2 / 9
     +   + s23**2 * s33 / 3
     +   + s31**2 * s33 / 3
     +   + 2 * s33**3 / 27
      return
      end subroutine invariantJ3

!-----------------------------------------------------------------------
      subroutine dInvariantJ3_dSigma(s11, s22, s33, s12, s23, s31, 
     +                               dJ3ds)
      implicit none
      real*8 s11, s22, s33, s12, s23, s31, dJ3ds(6)
      dJ3ds(1) = 2 * s11**2 / 9
     +         - 2 * s11 * s22 / 9
     +         - 2 * s11 * s33 / 9
     +         + s12**2 / 3
     +         - s22**2 / 9
     +         + 4 * s22 * s33 / 9
     +         - 2 * s23**2 / 3
     +         + s31**2 / 3
     +         - s33**2 / 9
    
      dJ3ds(2) = -(s11**2) / 9
     +         - 2 * s11 * s22 / 9
     +         + 4 * s11 * s33 / 9
     +         + s12**2 / 3
     +         + 2 * s22**2 / 9
     +         - 2 * s22 * s33 / 9
     +         + s23**2 / 3
     +         - 2 * s31**2 / 3
     +         - s33**2 / 9
    
      dJ3ds(3) = -(s11**2) / 9
     +         + 4 * s11 * s22 / 9
     +         - 2 * s11 * s33 / 9
     +         - 2 * s12**2 / 3
     +         - s22**2 / 9
     +         - 2 * s22 * s33 / 9
     +         + s23**2 / 3
     +         + s31**2 / 3
     +         + 2 * s33**2 / 9
!-----------------------------------------------------------------------
      dJ3ds(4) = s11*s12/3 + s12*s22/3 - 2*s12*s33/3 + s23*s31

      dJ3ds(5) = -2*s11*s23/3 + s12*s31 + s22*s23 / 3 + s23*s33/3

      dJ3ds(6) = s11*s31/3 + s12*s23 - 2*s31*s22/3 + s31*s33/3
      
      return
      end subroutine dInvariantJ3_dSigma



		subroutine phi_mises(s11, s22, s33, s12, s23, s31, phi)
		implicit none
      real*8 s11, s22, s33, s12, s23, s31, phi
		real*8 J2
		real*8 sd11, sd22, sd33, sd12, sd23, sd31
		real*8 sH

		sH = (s11 + s22 + s33)/3.

		sd11 = s11 - sH
		sd22 = s22 - sH
		sd33 = s33 - sH
		sd12 = s12
		sd23 = s23
		sd31 = s31

      J2 = 0.5*(sd11**2 + sd22**2 + sd33**2)
     <        + sd12**2 + sd23**2 + sd31**2
	  	phi = sqrt(3*J2) 
		return
		end subroutine phi_mises

		subroutine dphi_mises_ds(s11, s22, s33, s12, s23, s31, dphids)
		implicit none
      real*8 s11, s22, s33, s12, s23, s31, phi
		real*8 dphids(6)
		real*8 J2
		real*8 sd(6)
		real*8 sH

		sH = (s11 + s22 + s33)/3.

		sd(1) = s11 - sH
		sd(2) = s22 - sH
		sd(3) = s33 - sH
		sd(4) = s12
		sd(5) = s23
		sd(6) = s31

      J2 = 0.5*(sd(1)**2 + sd(2)**2 + sd(3)**2)
     <        + sd(4)**2 + sd(5)**2 + sd(6)**2
	  	
		dphids = 3.0/(2.0*sqrt(3*J2))*sd
		return
		end subroutine dphi_mises_ds