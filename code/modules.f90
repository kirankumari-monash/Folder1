!!! Time-stamp: <modules.f90 17:41, 31 Jan 2004 by P Sunthar>

!________________________________________________________________
!   Global modules and interface headers
!________________________________________________________________

!!! $Log: modules.f90,v $
!!! Revision 1.1  2004/01/29 22:13:11  pxs565
!!! Initial revision
!!!


Module Bspglocons ! Bead Spring simulation, Global constant (parameters)
  Save

  Integer, Parameter :: Ndim = 3 ! Dimension of simulation
  Integer, Parameter ::  MaxXdist = 101
  Integer, Parameter :: k4b = Selected_int_kind(9)
  Integer, Parameter :: DBprec = Selected_real_kind(8)
  Integer, Parameter :: NProps = 17
  Integer, Parameter :: MAXCHEB = 500
  Real, Parameter :: PI = 3.14159265358979323846
  Real, Parameter :: TINI = 1e-25
  Real, Parameter :: MYEPS = 1e-6
  Integer, Parameter :: HOOK = 1, FENE = 2, ILC = 3, WLC = 4, Fraenkel = 5
End Module Bspglocons


Module Bspglovars ! Bead Spring simulation, Global variables
  Use Bspglocons


!!$  Real Cubsoln_lu(0:1000), Gama_inc, Gama_max
!!$  Real, Allocatable ::Cubsoln_lu_2d(:,:), Gama_inc_2d(:), Gama_max_2d(:)

  Character(10) :: Prop_names(NProps) 

  Character(10) , Parameter :: Correl_names(1) = (/"Sxy"/)

  Real Imploop_tol, Tstep_conv_tol
End Module Bspglovars

Module Flowvars
  Integer, Parameter :: EQ = 0, SH = 1, UA = 2, PL = 3,  UR = 4, PU = 5, PP = 6
  Integer FlowType
       ! EQ equilibrium, no flow
       ! SH Planar shear
       ! UA Uniaxial Elongational
       ! PL Planar Elongation
       ! UR Uniaxial extension followed by relaxation
       ! PU Periodic uniaxial extension 
  Real gdots
end Module Flowvars


Module Bspintfacs
  Use Bspglocons

  !---------------------------------------------------------------------c	
  !     Driver subroutines                                              c
  !---------------------------------------------------------------------c

  Interface
     Subroutine Initial_position(Stype,N,L0s,R,seed)
       Use bspglocons
       Integer, Intent (in) :: Stype
       Integer(k4b), Intent (in) :: N
  Real, intent (in) :: L0s
       Integer(k4b), Intent (inout) :: seed
       Real, Intent (out) :: R(:,:)
       !Real, intent (out) :: R(Ndim,N)
     End Subroutine Initial_position
  End Interface

  Interface
     Subroutine GaussRand(N,GX,seed)
       Use bspglocons
       Integer(k4b), Intent (in) :: N
       Integer(k4b), Intent (inout) :: seed
       Real, Intent (out), Dimension(N) :: GX
     End Subroutine GaussRand
  End Interface


  Interface
     Subroutine Time_Integrate_Chain(NBeads, R_Bead, spring_type, &
          tcur, tmax, Delts, &
          Hstar, Zstar, Dstar, L0s, Q0s, &
          seed1, Nsamples, times, samples, phi)

       Use bspglocons
       Integer, Intent(in) :: NBeads
       Real, Intent (inout), Dimension(:,:) :: R_Bead ! Pos. vector of Beads
       
       Integer, Intent (in) :: spring_type ! Spring force law type
       ! Hookean = 1, FENE = 2, ILC = 3, WLC = 4
       
       
       Real, Intent (in) :: tcur,tmax,Delts  ! Integration time interval specs
       
       Real,  Intent (in) :: Hstar         ! HI parameter (non-dim)
       Real,  Intent (in) :: Zstar,Dstar   ! EV parameters (non-dim)
       Real,  Intent (in) :: L0s           ! finite ext., param sqrt(b)
       Real,  Intent (in) :: Q0s
       
       Integer (k4b), Intent (inout) :: seed1 ! seed for rnd number
       
       Integer, Intent(in) :: Nsamples     ! Number of sampling points
       Real,  Intent (in), Dimension(:) :: times    ! Sampling instances
       Real, Intent (inout), Dimension(:,:) :: samples
       Real, Intent (in), Dimension(:,:) :: phi
     End Subroutine Time_Integrate_Chain
  End Interface
  

  Interface
     Subroutine polish_poly_root(c,xin,atol)
       Implicit None
       ! use newton raphson to polish the root of a polynomial
       Integer, Parameter :: n=4
       Real, Intent (in) :: c(n)
       Real, Intent (inout) :: xin
       Real, Intent (in) :: atol
     End Subroutine polish_poly_root
  End Interface

  !---------------------------------------------------------------------c	
  !     Utility subroutines                                             c
  !---------------------------------------------------------------------c

  Interface
     Subroutine ran_1(n, R, idum)
       Use bspglocons
       Integer(k4b), Intent(inout) ::idum
       Integer(k4b), Intent(in) :: n
       Real, Intent(inout) :: R(n)
     End Subroutine ran_1
  End Interface

  Interface
     Subroutine maxminev_fi(n, A, maxev, minev)
       Integer n
       Real A(:,:,:,:), maxev, minev
     End Subroutine maxminev_fi
  End Interface

  Interface
     Subroutine chbyshv (L, d_a, d_b, a)
       Use bspglocons
       Integer L
       Real d_a, d_b, a(0:MAXCHEB)
     End Subroutine chbyshv
  End Interface

  Interface
     Subroutine cublu
       Use bspglocons
     End Subroutine cublu
  End Interface
  !---------------------------------------------------------------------c
  !     Properties estimators                                           c
  !---------------------------------------------------------------------c	

  Interface
     Subroutine chain_props(N, Rbead, F, props)
       Use bspglocons
       Implicit None
       Integer, intent (in) :: N
       Real, Intent (in), Dimension(:,:)  :: Rbead, F
       Real, Intent (out),Dimension(:) ::props
     End Subroutine chain_props
  End Interface

  Interface
     Subroutine time_correl(N, R, F, t0, t, correl)
       Use bspglocons
       Use bspglovars
       Implicit None
       Integer, intent (in) :: N
       Real, Intent (in), Dimension(:,:)  :: R, F
       Real, Intent (out)  :: correl
       Real, Intent (in) :: t0, t
     End Subroutine time_correl
  End Interface

  Interface
     Subroutine numint(f,dx,nmin,nmax,nord,sumf)
       Implicit None

       Real, Intent(in), Dimension(:) :: f
       Real, Intent(in) :: dx
       Real, Intent(out) :: sumf
       Integer, Intent (in) :: nmin, nmax,nord
     End Subroutine numint
  End Interface

  Interface
     Subroutine meanerr(vec,mean,err)
       Real, Intent(in) :: vec(:)
       Real, Intent(out) :: mean
       Real, Intent(out) :: err
     End Subroutine meanerr
  End Interface

  Interface
     Function normeqpr(Stype, r, b)
       Implicit None
       Integer, Intent (in) :: Stype
       Real, Intent (in) :: b, r
       Real  normeqpr
     End Function normeqpr
  End Interface

  Interface
     Function  lam1_th(hs, N)
       Real lam1_th
       Real,  Intent (in) :: hs
       Integer,  Intent (in) :: N
     End Function lam1_th
  End Interface

  Interface
     subroutine get_kappa(t,K)
       use bspglocons
       Implicit none
       Real, intent (in) :: t
       Real, intent (in), dimension(:,:) :: K
     end subroutine get_kappa
  end Interface
  
  
    
  End Module Bspintfacs
