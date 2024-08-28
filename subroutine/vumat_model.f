!-----------------------------------------------------------------------
!-----MJC Hershey implementation
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      subroutine vumat_model(sigma, deps, statev, props, ntens,
     <                        nstatev, nprops, rho, dt, time)
      implicit none
!-----------------------------------------------------------------------
!-----Declaration variables
!-----------------------------------------------------------------------
      real*8 sigma(ntens), deps(ntens), statev(nstatev), props(nprops)
      real*8 rho,dt,time
      integer ntens, nstatev, nprops
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      integer i
!-----Max number of iterations when plasticity occurs
      integer, parameter :: iter_max = 1000
!-----tolerance for update scheme
      real*8, parameter :: err_tol = 1e-8
!-----Used for special case when finding dfds
      real*8, parameter :: sing_tol = 1e-5
!-----PI
      real*8, parameter :: PI = 3.14159265358979323846264338327
!-----Elasticity constants
      real*8 E, nu
!-----Lamè parameters (lame1 = lambda, lame2 = mu = G)
      real*8 lame1, lame2
!-----Unique components isotropic elasticity matrix
      real*8 C11, C12, C44
!-----Yield function
      real*8 f
!-----Sigma equivalent
      real*8 phi
!-----Initial yield stress
      real*8 sigma0
!-----Hardening
      real*8 R
!-----Equivalent plastic strain 
      real*8 p, pold
!-----Temperature
      real*8 T
!-----Voce hardening constants
      real*8 Q1, C1, Q2, C2, Q3, C3 
!-----Hershey exponent
      real*8 n
!-----Thermal parameters
      real*8 cp, betaTQ, T0, Tm, m
!-----Viscoplastic parameters
      real*8 pdot0, c_visc
!-----Old stress and strain increment
      real*8 sold(6), de(6)
!-----Trial stress, deviatoric stress and stress components 
      real*8 tr(6), sdev(6), s11, s22, s33, s12, s23, s31
!-----hydrostatic stress
      real*8 sH
!-----Plastic multiplier increments
      real*8 ddlambda, dlambda
!-----gradient of yield function with respect to stresses
      real*8 dfds(6)
!-----Product of df_ds, C and df_ds
      real*8 dfds_C_dfds
!-----Product of C and dfds
      real*8 C_dfds(6)
!----- dzeta:h contribution from inner variables in cutting plane method
      real*8 dfdzeta_h
!-----Yield stress (sigma0 - R)
      real*8 sigmaY
!-----Gamma (temperature softening term)
      real*8 Gamma
!-----vp - Plastic strain rate term: (1+pdot)^c
      real*8 vp
!-----Hardening law change rate (dR/dp)
      real*8 hR
!-----Temperature softening rate (dGamma/dp)
      real*8 hGamma
!-----dvp/dp
      real*8 hv
!-----Lode angle
      real*8 Lode
!-----Invariants
      real*8 J2, J3
!-----Invariant gradients with respect to stress components
      real*8 dJ3ds(6)
!-----Principal stresses
      real*8 s1, s2, s3
!-----A = (s1-s2), B = (s2-s3), C=(s1-s3)
      real*8 A, B, C
!-----gradient of phi with respect to principal stresses
      real*8 dfds1, dfds2, dfds3
!-----Gradients of principal stresses with respect to stress components
      real*8 ds1ds(6), ds2ds(6), ds3ds(6)
!-----Gradient of Lode with respect to stress components
      real*8 dLodeds(6)
!-----Kronecker delta
      real*8 kronecker(6)
      data kronecker / 1., 1., 1., 0., 0., 0. /
!-----temporary real
      real*8 tmp
!-----------------------------------------------------------------------
!-----Read parameters and define constants
!-----------------------------------------------------------------------

      call assert((ntens.ne.4),
     <     "Plane stress implementation is not finished")
      call assert((nprops.eq.17),"nprops==17")
      call assert((nstatev.eq.2), "statev==2")
      E = props(1)
      nu = props(2)
      sigma0 = props(3)
      Q1 = props(4)
      C1 = props(5)
      Q2 = props(6)
      C2 = props(7)
      Q3 = props(8)
      C3 = props(9)
      n = props(10)
      cp = props(11)
      betaTQ = props(12)
      T0 = props(13)
      Tm = props(14)
      m = props(15)
      pdot0 = props(16)
      c_visc = props(17)
      call assert((n.ge.1).and.(n.le.100), "1 <= n (Hershey) <= 100")
      pold = statev(1)
      T = statev(2)
      p = pold
      dlambda = 0.0
      lame1 = nu*E/((1+nu)*(1-2*nu))
      lame2 = E/(2*(1+nu))

      if (ntens.eq.6) then
         C11 = lame1 + 2.*lame2
         C12 = lame1 
         C44 = 2.*lame2
      else
!-----Plane stress
         C11 = E/(1.-nu**2)
         C12 = nu*E/(1.-nu**2)
         C44 = E/(1.+nu)
      endif
!-----------------------------------------------------------------------
!-----Unpack old stresses
!-----------------------------------------------------------------------
      if (ntens.eq.6) then
         de = deps
         sold = sigma
      else
!-----Plane stress formulation
         call assert((ntens.eq.4),
     <   "ntens should be equal to 4 if shells are modelled")
         de(1) = deps(1)
         de(2) = deps(2)
         de(4) = deps(4)
         de(5) = 0.0
         de(6) = 0.0
         sold(1) = sigma(1)
         sold(2) = sigma(2)
         sold(3) = 0.0
         sold(4) = sigma(4)
         sold(5) = 0.0
         sold(6) = 0.0
      endif
!-----------------------------------------------------------------------
!-----Trial stress
!----------------------------------------------------------------------- 
      if (ntens.eq.6) then
         tr(1) = sold(1) + C11*de(1) + C12*de(2) + C12*de(3)
         tr(2) = sold(2) + C12*de(1) + C11*de(2) + C12*de(3)
         tr(3) = sold(3) + C12*de(1) + C12*de(2) + C11*de(3)
         tr(4) = sold(4) + C44*de(4)
         tr(5) = sold(5) + C44*de(5)
         tr(6) = sold(6) + C44*de(6)
      else
         tr(1) = sold(1) + C11*de(1) + C12*de(2) 
         tr(2) = sold(2) + C12*de(1) + C11*de(2) 
         tr(3) = 0.0
         tr(4) = sold(4) + C44*de(4)
         tr(5) = 0.0
         tr(6) = 0.0
      endif
!-----Setting the viscoplastic parameter to 1.0 initally, which means 
!-----no rate dependence in the trial step
      vp = 1.0
!-----------------------------------------------------------------------
      !Gather stress components
      s11 = tr(1) 
      s22 = tr(2) 
      s33 = tr(3) 
      s12 = tr(4)
      s23 = tr(5)
      s31 = tr(6)
!-----------------------------------------------------------------------
!-----Start inner loop (loop breaks at i=1 if the increment is elastic)
!-----------------------------------------------------------------------
      do i=0,iter_max 
!-----J2 
         J2 = 0.5 * (s11**2 + s22**2 + s33**2 + 
     <        2 * (s12**2 + s23**2 + s31**2)) -
     <        1.0 / 6 * (s11 + s22 + s33) ** 2
!-----J3  
         call invariantJ3(s11,s22,s33,s12,s23,s31,J3)
!-----------------------------------------------------------------------
!-----Calulating Lode angle: (This includes a fix in case the argument 
!-----is slightly out of the allowed range for arccos (-1 to 1) 
!-----------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////
CODE REMOVED - COMPLETE VERSION IS AVAILABLE AT CAE-ASSISTANT.COM 
////////////////////////////////////////////////////////////////////////

!-----------------------------------------------------------------------
!-----Computing yield function f
!-----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////
CODE REMOVED - COMPLETE VERSION IS AVAILABLE AT CAE-ASSISTANT.COM 
////////////////////////////////////////////////////////////////////////

!-----------------------------------------------------------------------
!-----CHECK IF YIELDING OCCURS
!-----------------------------------------------------------------------
         f = phi - sigmaY*Gamma*vp
         if (i.eq.0) then
!-----Yield function
            if (f.le.0) then
!-----------------------------------------------------------------------
!-----f<=0: Elastic increment
!-----------------------------------------------------------------------
               exit
            endif
!-----------------------------------------------------------------------
         else !----f>0 - Plastic increment
!-----------------------------------------------------------------------
!-----Convergence check
!-----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////
CODE REMOVED - COMPLETE VERSION IS AVAILABLE AT CAE-ASSISTANT.COM 
////////////////////////////////////////////////////////////////////////

!-----------------------------------------------------------------------
!-----Gradient of f with respect to stress components
!-----------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////
CODE REMOVED - COMPLETE VERSION IS AVAILABLE AT CAE-ASSISTANT.COM 
////////////////////////////////////////////////////////////////////////

!-----------------------------------------------------------------------
!-----Computing dfds:C:dfds
!-----------------------------------------------------------------------

            C_dfds(1) = C11*dfds(1)+ C12*dfds(2)+ C12*dfds(3)
            C_dfds(2) = C12*dfds(1)+ C11*dfds(2)+ C12*dfds(3)
            C_dfds(3) = C12*dfds(1)+ C12*dfds(2)+ C11*dfds(3)
            C_dfds(4) = C44*dfds(4)
            C_dfds(5) = C44*dfds(5)
            C_dfds(6) = C44*dfds(6)
!
            dfds_C_dfds = dfds(1)*C_dfds(1)+ 
     <                    dfds(2)*C_dfds(2)+ 
     <                    dfds(3)*C_dfds(3)+ 
     <                2.*(dfds(4)*C_dfds(4)+
     <                    dfds(5)*C_dfds(5)+ 
     <                    dfds(6)*C_dfds(6))     
!-----------------------------------------------------------------------
!-----Cumputing various contributions from inner variables
!-----------------------------------------------------------------------
         hR = Q1*C1*exp(-C1*p) + Q2*C2*exp(-C2*p) + C3*Q3*exp(-C3*p)
         vp = (1.0 + dlambda/(pdot0*dt))**c_visc
         hv = c_visc/(pdot0*dt)*(1.0 + dlambda/(pdot0*dt))**(c_visc-1.0) 
!-----Adding all inner variable contributions
         dfdzeta_h = -hR*Gamma*vp - sigmaY*Gamma*hv
!-----------------------------------------------------------------------
!-----Computing dlambda increment
!-----------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////
CODE REMOVED - COMPLETE VERSION IS AVAILABLE AT CAE-ASSISTANT.COM 
////////////////////////////////////////////////////////////////////////

!-----------------------------------------------------------------------
!-----Updating stresses and internal variables 
!-----------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////
CODE REMOVED - COMPLETE VERSION IS AVAILABLE AT CAE-ASSISTANT.COM 
////////////////////////////////////////////////////////////////////////

!-----------------------------------------------------------------------
!-----Pack internal variables
!-----------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////
CODE REMOVED - COMPLETE VERSION IS AVAILABLE AT CAE-ASSISTANT.COM 
////////////////////////////////////////////////////////////////////////

      end