! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Global Data Module File
! 
! Generated by KPP-2.2.1_rs5 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : gckpp_Global.f90
! Time                 : Mon Jun  2 12:32:46 2014
! Working directory    : /home/mpayer/KPP
! Equation file        : gckpp.kpp
! Output root filename : gckpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE gckpp_Global

  USE gckpp_Parameters, ONLY: dp, NSPEC, NVAR, NFIX, NREACT
  PUBLIC
  SAVE


! Declaration of global variables

! C - Concentration of all species
  REAL(kind=dp) :: C(NSPEC)
! VAR - Concentrations of variable species (global)
  REAL(kind=dp) :: VAR(NVAR)
! FIX - Concentrations of fixed species (global)
  REAL(kind=dp) :: FIX(NFIX)
! VAR, FIX are chunks of array C
!      EQUIVALENCE( C(1),VAR(1) )
!      EQUIVALENCE( C(140),FIX(1) )
! RCONST - Rate constants (global)
  REAL(kind=dp) :: RCONST(NREACT)
! TIME - Current integration time
  REAL(kind=dp) :: TIME
! SUN - Sunlight intensity between [0,1]
  REAL(kind=dp) :: SUN
! TEMP - Temperature
  REAL(kind=dp) :: TEMP
! TSTART - Integration start time
  REAL(kind=dp) :: TSTART
! TEND - Integration end time
  REAL(kind=dp) :: TEND
! DT - Integration step
  REAL(kind=dp) :: DT
! ATOL - Absolute tolerance
  REAL(kind=dp) :: ATOL(NVAR)
! RTOL - Relative tolerance
  REAL(kind=dp) :: RTOL(NVAR)
! STEPMIN - Lower bound for integration step
  REAL(kind=dp) :: STEPMIN
! STEPMAX - Upper bound for integration step
  REAL(kind=dp) :: STEPMAX
! CFACTOR - Conversion factor for concentration units
  REAL(kind=dp) :: CFACTOR
! DDMTYPE - DDM sensitivity w.r.t.: 0=init.val., 1=params
  INTEGER :: DDMTYPE

! VAR_ADJ - Concentrations of variable species (global) [**]
  REAL(kind=dp) :: VAR_ADJ(NVAR)
! V_CSPEC - Concentrations of variable species (global)
  REAL(kind=dp) :: V_CSPEC(NVAR)
! V_CSPEC_ADJ - Concentrations of variable species (global) [**]
  REAL(kind=dp) :: V_CSPEC_ADJ(NVAR)
   
! NJ - Number of cost function being evaluated [**]
  INTEGER, PARAMETER :: NJ = 1 
! NTT - Total number of tropospheric grid cells
  INTEGER :: NTT
! JLOOP - Total number of tropospheric grid cells
  INTEGER :: JLOOP
! SMAL2 - Parameter for insuring positive tracer values, same as in reader.f
  REAL(kind=dp), PARAMETER :: SMAL2 = 1.0d-99 
! NCOEFF - Number of reaction rate coeff adjoints [**]
  INTEGER, PARAMETER :: NCOEFF = 24
! VAR_R_ADJ - Concentrations of reaction rate adjoint (global) [**]
  REAL(kind=dp) :: VAR_R_ADJ(NCOEFF)
! JCOEFF - Reaction numbers for each (define in INIT_KPP) [**]
  INTEGER :: JCOEFF(NCOEFF)
! IND - Reaction numbers for each (define in INIT_KPP)
  INTEGER :: IND(NREACT)

! INLINED global variable declarations

! INLINED global variable declarations

!$OMP THREADPRIVATE(VAR,VAR_ADJ,VAR_R_ADJ,V_CSPEC,V_CSPEC_ADJ, C )
!$OMP THREADPRIVATE(FIX,JLOOP,RCONST,TIME)

! Move stack_ptr here and make THREADPRIVATE for OMP parallelization (dkh, 07/28/09) 
  INTEGER :: stack_ptr = 0 ! last written entry
!$OMP THREADPRIVATE( stack_ptr )

END MODULE gckpp_Global

