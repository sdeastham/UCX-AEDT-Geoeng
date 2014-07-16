!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: longlived_mod
!     
! !DESCRIPTION: Module LONGLIVED\_MOD contains routines and variables
!  used to apply surface mixing ratio boundary conditions on
!  long-lived species.
!\\   
!\\   
! !INTERFACE: 
!
      MODULE LONGLIVED_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN
      USE ERROR_MOD,  ONLY : DEBUG_MSG
      USE CMN_SIZE_MOD

      ! NcdfUtil modules for netCDF I/O
      USE m_netcdf_io_open                    ! netCDF open
      USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
      USE m_netcdf_io_read                    ! netCDF data reads
      USE m_netcdf_io_close                   ! netCDF close

      IMPLICIT NONE
#     include "netcdf.inc"
      PRIVATE

!
! !PUBLIC DATA MEMBERS:
!
      INTEGER,            PUBLIC :: LLSYEAR ! Year for CFC emissions
      INTEGER,            PUBLIC :: OFFSETYEAR
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: EMISS_BASIC
      PUBLIC  :: LLS_SCALEFAC
      PUBLIC  :: INIT_LLS
      PUBLIC  :: CLEANUP_LLS
      PUBLIC  :: GET_SCALED_CH4
      PUBLIC  :: GET_SCALED_CH3Br
!
! PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: READ_SFC
      PRIVATE :: GET_RCP_NCDF
      PRIVATE :: GET_WMO_NCDF
      PRIVATE :: GET_LLS_MR
      PRIVATE :: CHECK_LLSYEAR
!
! !REVISION HISTORY: 
!  27 Mar 2014 - S. D. Eastham - Extracted routines from UCX_MOD and
!                                added RCP functionality
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! MODULE PARAMETERS
      !
      ! MAX_EMIT           : Limit on tracer emission count
      !
      !=================================================================

      INTEGER, PARAMETER                :: MAX_EMIT=50

!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Scalars
      !
      ! N_EMIT             : Number of tracers to emit
      ! SFC_MONTH          : Current month
      ! OFFSETYEAR         : Number of year-ends since simulation start
      !
      ! Numerical arrays
      ! IDT_EMIT           : Tracer ID to be emitted to
      ! GRID_EMIT          : Vector to hold output surface MR (by lat)
      !
      ! Strings
      !
      ! WMO_FILE_NCDF      : Location of file containing WMO data
      ! RCP_FILE_NCDF      : Location of file containing RCP data
      ! 
      !=================================================================

      ! Scalars
      INTEGER                           :: N_EMIT
      INTEGER                           :: SFC_MONTH

      ! Numerical arrays
      INTEGER,DIMENSION(MAX_EMIT)       :: IDT_EMIT
      REAL*8,DIMENSION(:,:),ALLOCATABLE :: GRID_EMIT
 
      ! Strings
      CHARACTER(LEN=255)                :: WMO_FILE_NCDF
      CHARACTER(LEN=255)                :: RCP_FILE_NCDF

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_sfc
!
! !DESCRIPTION: Subroutine READ\_SFC fills out the surface emissions for some
!  species introduced as part of the unified chemistry upgrade.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_SFC( am_I_Root, Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TRACERID_MOD,       ONLY : IDTOCS, IDTH2, IDTCH3Br
      USE TRACERID_MOD,       ONLY : IDTHCFCX, IDTCFCX, IDTCH4
      USE GRID_MOD,           ONLY : GET_YMID
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TRACER_MOD,         ONLY : TRACER_NAME
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT VARIABLES:
!
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt object via the arg list
!  27 Mar 2014 - S. D. Eastham - Updated for RCP and made NetCDF only
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8,PARAMETER :: OCS_MR = 5.00d-10 ! OCS is fixed to 500 pptv
      REAL*8,PARAMETER :: H2_MR = 5.00d-7  ! H2 is fixed to 500 ppbv
      REAL*8           :: YLAT
      REAL*8           :: MR90S,MR30S,MR30N,MR90N,MR55S,MR55N
      INTEGER          :: N, TARG_LINE, J
      CHARACTER(LEN=20):: LOCAL_NAME
      INTEGER          :: NCYR, NCMO, TARG_IDT
      CHARACTER(LEN=14):: LUCNAME, UCNAME
      LOGICAL          :: LRCPLLS

      !=================================================================
      ! READ_SFC begins here!
      !=================================================================

      ! Re-zero GRID_EMIT
      GRID_EMIT(:,:) = 0d0
         
      NCYR = LLSYEAR+OFFSETYEAR
      NCMO = GET_MONTH()
      DO N=1,N_EMIT
         TARG_IDT = IDT_EMIT(N)
         
         IF (TARG_IDT == IDTCFCX) THEN
            ! Add 3 CFCs
            GRID_EMIT(:,N) = 
     &         GET_LLS_MR('CFC-113',NCYR,NCMO,Input_Opt,am_I_Root) + 
     &         GET_LLS_MR('CFC-114',NCYR,NCMO,Input_Opt,am_I_Root) + 
     &         GET_LLS_MR('CFC-115',NCYR,NCMO,Input_Opt,am_I_Root)
         ELSEIF (TARG_IDT == IDTHCFCX) THEN 
            ! Add 2 HCFCs
            GRID_EMIT(:,N) = 
     &         GET_LLS_MR('HCFC-141b',NCYR,NCMO,Input_Opt,am_I_Root) + 
     &         GET_LLS_MR('HCFC-142b',NCYR,NCMO,Input_Opt,am_I_Root)
         ELSEIF (TARG_IDT == IDTCH4) THEN 
            ! Get the 4 major mixing ratios
            CALL GET_SCALED_CH4(MR90S,MR30S,MR30N,MR90N,
     &                          am_I_Root,Input_Opt)
            ! Figure out which box each is in
            DO J=1,JJPAR   
               ! Grid-box latitude index
               YLAT  = GET_YMID( 1,J,1 )
               IF ( YLAT < -30d0 ) THEN
                  GRID_EMIT(J,N) = MR90S
               ELSE IF ( YLAT >= -30d0 .and. YLAT < 0d0  ) THEN
                  GRID_EMIT(J,N) = MR30S
               ELSE IF ( YLAT >=   0d0 .and. YLAT < 30d0 ) THEN
                  GRID_EMIT(J,N) = MR30N
               ELSE
                  GRID_EMIT(J,N) = MR90N
               ENDIF
            ENDDO
         ELSEIF (TARG_IDT == IDTOCS) THEN 
            ! Simple boundary condition
            GRID_EMIT(:,N) = OCS_MR
         ELSEIF (TARG_IDT == IDTH2) THEN 
            ! Simple boundary condition
            GRID_EMIT(:,N) = H2_MR
         ELSEIF (TARG_IDT == IDTCH3Br) THEN
            ! 4 latitude bands
            CALL GET_SCALED_CH3Br(MR90S,MR55S,MR55N,MR90N,
     &                                am_I_Root,Input_Opt)
            DO J=1,JJPAR   
               ! Grid-box latitude index
               YLAT  = GET_YMID( 1,J,1 )
               IF ( YLAT < -55d0 ) THEN
                  GRID_EMIT(J,N) = MR90S
               ELSE IF ( YLAT >= -55d0 .and. YLAT < 0d0  ) THEN
                  GRID_EMIT(J,N) = MR55S
               ELSE IF ( YLAT >=   0d0 .and. YLAT < 55d0 ) THEN
                  GRID_EMIT(J,N) = MR55N
               ELSE
                  GRID_EMIT(J,N) = MR90N
               ENDIF
            ENDDO
         ELSE

            LUCNAME = TRACER_NAME(TARG_IDT)
            ! Standard 1:1
            GRID_EMIT(:,N) = 
     &           GET_LLS_MR(TRIM(LUCNAME),NCYR,NCMO,Input_Opt,am_I_Root)
         ENDIF
      ENDDO

      END SUBROUTINE READ_SFC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_scaled_ch3br
!
! !DESCRIPTION: Subroutine GET\_SCALED\_CH3BR retrieves zonal mean surface mixing
!  ratios for 4 regions (90S-55S,55S-0,0-55N,55N-90N)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_SCALED_CH3Br(MR90S,MR55S,MR55N,MR90N,am_I_Root,
     &                                       Input_Opt)
!
! !USES:
!
      USE TIME_MOD,              ONLY : GET_YEAR
      USE FUTURE_EMISSIONS_MOD,  ONLY : GET_FUTURE_YEAR
      USE GIGC_Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT)         :: MR90S,MR55S,MR55N,MR90N
!
! !REVISION HISTORY: 
!  18 Dec 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CH3Br values ( from pptv )
      REAL*8, PARAMETER :: gt55N     = 8.35d-12 
      REAL*8, PARAMETER :: gt0_lt55N = 8.27d-12
      REAL*8, PARAMETER :: lt0_gt55S = 6.94d-12
      REAL*8, PARAMETER :: lt55S     = 6.522d-12
!
! !LOCAL VARIABLES:
!
      REAL*8                :: CH3Br_SCALE, LOWLATSUM, UCX_CH3Br
      INTEGER, SAVE         :: LAST_YEAR=-1
      INTEGER               :: CH3Br_YEAR
      REAL*8, SAVE          :: LAST_MR(4)
      REAL*8, SAVE          :: GC_CH3Br

      !=================================================================
      ! GET_SCALED_CH3Br begins here!
      !=================================================================

      ! Calculate RCP adjustment factor
      CALL CHECK_LLSYEAR(am_I_Root,Input_Opt)
      CH3Br_YEAR = LLSYEAR + OFFSETYEAR

      IF (CH3Br_YEAR.ne.LAST_YEAR) THEN
         IF (LAST_YEAR.lt.0) THEN
            ! Calculate mean adjustment
            ! This performs area weighting for 4 segments divided at
            ! +90, +55, 0, -55 and -90
            LOWLATSUM = gt0_lt55N + lt0_gt55S
            GC_CH3Br = (DSIN(55.d0)*(gt55N+lt55S-LOWLATSUM)+LOWLATSUM)/2.d0
         ENDIF
         UCX_CH3Br = GET_LLS_MR('CH3Br',CH3Br_YEAR,1,
     &                                       Input_Opt,am_I_Root)
         CH3Br_SCALE = UCX_CH3Br/GC_CH3Br

         LAST_MR(1) = lt55s*CH3Br_SCALE
         LAST_MR(2) = lt0_gt55S*CH3Br_SCALE
         LAST_MR(3) = gt0_lt55N*CH3Br_SCALE
         LAST_MR(4) = gt55N*CH3Br_SCALE

         ! Debug
         IF (am_I_Root.and.Input_Opt%LPRT) THEN
            WRITE(6,'(a,e10.4,a)') 'Scaled CH3Br N90: ',
     &                                              LAST_MR(1), ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH3Br N55: ',
     &                                              LAST_MR(2), ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH3Br S55: ',
     &                                              LAST_MR(3), ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH3Br S90: ',
     &                                              LAST_MR(4), ' (v/v)'
         ENDIF

         LAST_YEAR = CH3Br_YEAR
      ENDIF

      MR90S = LAST_MR(1)
      MR55S = LAST_MR(2)
      MR55N = LAST_MR(3)
      MR90N = LAST_MR(4)

      END SUBROUTINE GET_SCALED_CH3Br
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_scaled_ch4
!
! !DESCRIPTION: Subroutine GET\_SCALED\_CH4 retrieves zonal mean surface mixing
!  ratios for 4 equal-area regions (90S-30S,30S-0,0-30N,30N-90N)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_SCALED_CH4(MR90S,MR30S,MR30N,MR90N,am_I_Root,
     &                                     Input_Opt)
!
! !USES:
!
      USE TIME_MOD,              ONLY : GET_YEAR, GET_MONTH
      USE FUTURE_EMISSIONS_MOD,  ONLY : GET_FUTURE_YEAR
      USE GIGC_Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT)         :: MR90S,MR30S,MR30N,MR90N
!
! !REVISION HISTORY: 
!  18 Dec 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: BASE_YEAR, CH4_TARG_YEAR
      INTEGER, PARAMETER    :: REF_YEAR=2007
      REAL*8                :: CH4_SCALE
      INTEGER, SAVE         :: LAST_YEAR=-1
      REAL*8, SAVE          :: LAST_MR(4)
      REAL*8, SAVE          :: BASE_MR(4)
      INTEGER               :: CH4_YEAR   ! Target year

      !=================================================================
      ! GET_SCALED_CH4 begins here!
      !=================================================================

      ! Adapted from GEOS-Chem routine
      CALL CHECK_LLSYEAR(am_I_Root,Input_Opt)
      CH4_YEAR = LLSYEAR + OFFSETYEAR

      ! Pick the CH4 concentration [ppbv] for the proper lat bin
      ! CH4 values are read in "chemdr.f"
      IF ( Input_Opt%LFUTURE ) THEN
         CH4_TARG_YEAR = GET_FUTURE_YEAR()
         BASE_YEAR = CH4_TARG_YEAR
      ELSE
         CH4_TARG_YEAR = LLSYEAR + OFFSETYEAR
         IF (Input_Opt%LRCPLLS) THEN
            BASE_YEAR = 2007
         ELSE
            BASE_YEAR = MIN(REF_YEAR,CH4_TARG_YEAR)
         ENDIF
      ENDIF

      IF (LAST_YEAR .ne. CH4_TARG_YEAR) THEN
         IF (LAST_YEAR.lt.0) THEN
            ! Get CH4 [ppbv] in 4 latitude bins for 2007 (last year with
            ! gradient information)
            CALL GET_GLOBAL_CH4( BASE_YEAR,.TRUE., MR90S,
     &                              MR30S,   MR30N, MR90N,
     &                              am_I_Root, Input_Opt )
            BASE_MR(1) = MR90S
            BASE_MR(2) = MR30S
            BASE_MR(3) = MR30N
            BASE_MR(4) = MR90N
         ELSE
            MR90S = BASE_MR(1)
            MR30S = BASE_MR(2)
            MR30N = BASE_MR(3)
            MR90N = BASE_MR(4)
         ENDIF

         ! Each zone has the same area
         CH4_SCALE = GET_LLS_MR('CH4',LLSYEAR+OFFSETYEAR,1,
     &                 Input_Opt,am_I_Root)/
     &                 (0.25d0*(MR90S+MR30S+MR30N+MR90N))

         ! Convert from [ppbv] to [vv] and scale by year
         ! Using just one month's data to be consistent with
         ! base GEOS-Chem approach to CH4 (1-yr resolution)

         MR90S = MR90S * CH4_SCALE
         MR30S = MR30S * CH4_SCALE
         MR30N = MR30N * CH4_SCALE
         MR90N = MR90N * CH4_SCALE

         ! Debug
         IF (am_I_Root.and.Input_Opt%LPRT) THEN
            WRITE(6,'(a,e10.4,a)') 'Scaled CH4 N90: ', MR90N, ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH4 N30: ', MR30N, ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH4 S30: ', MR30S, ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH4 S90: ', MR90S, ' (v/v)'
         ENDIF

         ! Save for later calls
         LAST_YEAR = CH4_TARG_YEAR
         LAST_MR(1) = MR90S
         LAST_MR(2) = MR30S
         LAST_MR(3) = MR30N
         LAST_MR(4) = MR90N
      ELSE
         ! Re-use data
         MR90S = LAST_MR(1)
         MR30S = LAST_MR(2)
         MR30N = LAST_MR(3)
         MR90N = LAST_MR(4)
      ENDIF

      END SUBROUTINE GET_SCALED_CH4
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lls_scalefac
!
! !DESCRIPTION: Functoin LLS\_SCALEFAC returns (RCP or WMO) the
!  ratio of the projected surface MR in a target year to that in a
!  reference year.
!\\
!\\
! !INTERFACE:
!
      FUNCTION LLS_SCALEFAC(MS_TRAC,MS_INITYR,MS_TARGMO,
     &                 Input_Opt,am_I_Root) RESULT(SCALEFAC)
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input options
      CHARACTER(*),INTENT(IN)     :: MS_TRAC
      INTEGER,INTENT(IN)          :: MS_INITYR  ! Baseline year
      INTEGER,INTENT(IN)          :: MS_TARGMO  ! Target month
      LOGICAL,INTENT(IN)          :: am_I_Root
!
! !OUTPUT VARIABLES:
!
      REAL*8                      :: SCALEFAC   ! Scaling factor
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Initial version
!  27 Mar 2014 - S. D. Eastham - Revised for RCP scenarios
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                     :: MS_FUTRYR  ! Future year

      !=================================================================
      ! LLS_SCALEFAC begins here!
      !=================================================================
       
      MS_FUTRYR = LLSYEAR + OFFSETYEAR
      SCALEFAC = 
     & GET_LLS_MR(MS_TRAC,MS_FUTRYR,MS_TARGMO,Input_Opt,am_I_Root)/
     & GET_LLS_MR(MS_TRAC,MS_INITYR,MS_TARGMO,Input_Opt,am_I_Root)

      END FUNCTION LLS_SCALEFAC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_llsyear
!
! !DESCRIPTION: Subroutine CHECK\_LLSYEAR updates LLSYEAR and checks if a
!  new month has occurred, necessitating a new surface read-in.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_LLSYEAR( am_I_Root, Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : GET_MONTH
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      LOGICAL, INTENT(IN)           :: am_I_Root
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  15 Jul 2014 - S. D. Eastham - Moved out of READ_SFC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255) :: DBGMSG
      !=================================================================
      ! CHECK_LLSYEAR begins here!
      !=================================================================

      ! If it's a new month, read in new fields
      IF (SFC_MONTH/=GET_MONTH()) THEN
         IF ((Input_Opt%LRCPLLS).and.(Input_Opt%LFIXLLS)) THEN
            OFFSETYEAR = 0
         ELSEIF (GET_MONTH().lt.SFC_MONTH) THEN
            OFFSETYEAR = OFFSETYEAR + 1
         ENDIF
         SFC_MONTH = GET_MONTH()
         IF (Input_Opt%LPRT.and.am_I_Root) THEN
            WRITE(DBGMSG,'(a,a,1x,I4)') ' ### LLS: ',
     &         'Using reference year', (LLSYEAR+OFFSETYEAR)
            CALL DEBUG_MSG( DBGMSG )
         ENDIF
            
         CALL READ_SFC( am_I_Root, Input_Opt )
      END IF

      END SUBROUTINE CHECK_LLSYEAR
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_basic
!
! !DESCRIPTION: Subroutine EMISS\_BASIC sets surface mixing ratios of N2O,
!  OCS and ozone-depleting substances covered by the Montreal protocol.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_BASIC( am_I_Root, Input_Opt, State_Met,
     &                        State_Chm  )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_MONTH
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT VARIABLES:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Initial version
!  11 Feb 2014 - R. Yantosca   - Rewrote DO loop to be more efficient
!  14 Feb 2014 - R. Yantosca   - Now parallelize main DO loop
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm 
!                                objects via the arg list
!  27 Mar 2014 - S. D. Eastham - Consolidated into new module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, N, N_EMITTED
      REAL*8             :: N2O_MULT,OCS_MULT,CFC_MULT,HCFC_MULT
      REAL*8             :: DT_YR
      CHARACTER(LEN=255) :: OUTMSG
      LOGICAL            ::NEWMONTH

      ! Local variables for quantities from Input_Opt
      LOGICAL            :: LPRT
      INTEGER            :: N_TRACERS
      REAL*8             :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL*8, POINTER    :: STT(:,:,:,:)

      !=================================================================
      ! EMISS_BASIC begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LPRT      = Input_Opt%LPRT
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      ! If not correct hour block and day read in new fields
      ! If it's a new month, read in new fields
      CALL CHECK_LLSYEAR( am_I_Root, Input_Opt )
  
       ! Loop over emitted species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, N_EMITTED )
!$OMP+SCHEDULE( DYNAMIC )
       DO N = 1, N_EMIT

          ! Get the tracer # of the emitted species
          N_EMITTED = IDT_EMIT(N)

          ! Convert units if the tracer is valid
          IF ( N_EMITTED /= 0 ) THEN
        
             ! Loop over surface grmid boxes
             DO J = 1, JJPAR
             DO I = 1, IIPAR
           
                ! Set kg/box using vol/vol (see CONVERT_UNITS in DAO_MOD)
                STT(I,J,1,N_EMITTED) = GRID_EMIT(J,N)
     &                               * State_Met%AD(I,J,1)   
     &                               / TCVV(N_EMITTED)

             ENDDO
             ENDDO
          ENDIF
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      IF ( LPRT ) CALL DEBUG_MSG( '### LLS: Emissions complete' )

      END SUBROUTINE EMISS_BASIC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_lls_mr
!
! !DESCRIPTION: Function GET\_LLS\_MR is a wrapper function which
! retrieves either the WMO or RCP estimate for a surface mixing ratio
! depending on user input switches.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_LLS_MR (TRAC_NAME,TARG_YR_NC,TARG_MO_NC,
     &                  Input_Opt,am_I_Root) RESULT (MEAN_MIXRATIO)
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : TRANUC
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      REAL*8                        :: MEAN_MIXRATIO
      CHARACTER*(*),INTENT(IN)      :: TRAC_NAME
      INTEGER,INTENT(IN)            :: TARG_YR_NC
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      INTEGER,INTENT(IN)            :: TARG_MO_NC
      LOGICAL,INTENT(IN)            :: am_I_Root
!
! !REVISION HISTORY: 
!  27 Mar 2014 - S. D. Eastham - Initial version
!  15 Jul 2014 - S. D. Eastham - Updated for v10-01c
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: TARG_TRAC, DBGMSG

      !=================================================================
      ! GET_LLS_MR begins here!
      !=================================================================

      ! Get a full upper-case version of TRAC_NAME
      TARG_TRAC = TRAC_NAME
      CALL TRANUC( TARG_TRAC )

      ! Are we using RCP or WMO?
      ! In each case, translate cases which are named differently in
      ! GEOS-Chem and in the NetCDF
      IF (Input_Opt%LRCPLLS) THEN
         SELECT CASE ( TRIM( TARG_TRAC ) )
            CASE ( 'CFC11' )
               TARG_TRAC = 'CFC-11'
            CASE ( 'CFC12' )
               TARG_TRAC = 'CFC-12'
            CASE ( 'CFC113' )
               TARG_TRAC = 'CFC-113'
            CASE ( 'CFC114' )
               TARG_TRAC = 'CFC-114'
            CASE ( 'CFC115' )
               TARG_TRAC = 'CFC-115'
            CASE ( 'HCFC22' )
               TARG_TRAC = 'HCFC-22'
            CASE ( 'HCFC141B' )
               TARG_TRAC = 'HCFC-141b'
            CASE ( 'HCFC142B' )
               TARG_TRAC = 'HCFC-142b'
            CASE ( 'CCL4' )
               TARG_TRAC = 'CCl4'
            CASE ( 'H1211' )
               TARG_TRAC = 'H-1211'
            CASE ( 'H1301' )
               TARG_TRAC = 'H-1301'
            CASE ( 'H2402' )
               TARG_TRAC = 'H-2402'
            CASE ( 'C2BR2F4' )
               ! Alias
               TARG_TRAC = 'H-2402'
            CASE ( 'CH3BR' )
               TARG_TRAC = 'CH3Br'
            CASE ( 'CH3CL' )
               TARG_TRAC = 'CH3Cl'
            CASE ( 'CH3CCL3' )
               TARG_TRAC = 'CH3CCl3'
            CASE DEFAULT
               ! Restore case to be the same as before
               TARG_TRAC = TRAC_NAME
         END SELECT
         MEAN_MIXRATIO = GET_RCP_NCDF(TARG_TRAC,TARG_YR_NC,
     &                        TARG_MO_NC,Input_Opt,am_I_Root)
      ELSE
!         WRITE(6,'(a,a)') 'Normal   : ',TRIM(TRAC_NAME)
!         WRITE(6,'(a,a)') 'Uppercase: ',TRIM(TARG_TRAC)
         SELECT CASE ( TRIM( TARG_TRAC ) )
            CASE ( 'CFC-113' )
               TARG_TRAC = 'CFC113'
            CASE ( 'CFC-114' )
               TARG_TRAC = 'CFC114'
            CASE ( 'CFC-115' )
               TARG_TRAC = 'CFC115'
            CASE ( 'HCFC-141B' )
               TARG_TRAC = 'HCFC141b'
            CASE ( 'HCFC-142B' )
               TARG_TRAC = 'HCFC142b'
            CASE ( 'H24O2' )
               TARG_TRAC = 'C2BR2F4'
            CASE ( 'CCL4' )
               TARG_TRAC = 'CCl4'
            CASE ( 'CH3BR' )
               TARG_TRAC = 'CH3Br'
            CASE ( 'CH3CL' )
               TARG_TRAC = 'CH3Cl'
            CASE ( 'CH3CCL3' )
               TARG_TRAC = 'CH3CCl3'
            CASE ( 'H2402' )
               TARG_TRAC = 'C2BR2F4'
            CASE DEFAULT
               ! Restore case to be the same as before
               TARG_TRAC = TRAC_NAME
         END SELECT
         MEAN_MIXRATIO = GET_WMO_NCDF(TARG_TRAC,TARG_YR_NC,
     &             TARG_MO_NC,Input_Opt,am_I_Root)
      ENDIF

      IF (Input_Opt%LPRT.and.am_I_Root) THEN
         WRITE(DBGMSG,'(a,a,a,E14.6,a)') 'Reading ', TRIM(TRAC_NAME), 
     &     ' MR: ', MEAN_MIXRATIO, ' (v/v)'
         CALL DEBUG_MSG(TRIM(DBGMSG))
      ENDIF

      ! Return to calling program
      END FUNCTION GET_LLS_MR
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_rcp_ncdf
!
! !DESCRIPTION: Function GET\_RCP\_NCDF gets mean surface mixing ratio
!  of a given species based on the current RCP scenario.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_RCP_NCDF (TRAC_NAME,TARG_YR_NC,TARG_MO_NC,Input_Opt,
     &                          am_I_Root) RESULT (MEAN_MIXRATIO)
!
! !USES:
!
      USE FILE_MOD,           ONLY : IOERROR
      USE ERROR_MOD,          ONLY : ERROR_STOP, ALLOC_ERR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      REAL*8                         :: MEAN_MIXRATIO
      CHARACTER*(*),INTENT(IN)       :: TRAC_NAME
      INTEGER,INTENT(IN)             :: TARG_YR_NC
      INTEGER,INTENT(IN)             :: TARG_MO_NC 
      TYPE(OptInput), INTENT(IN)     :: Input_Opt   ! Input options
      LOGICAL,INTENT(IN)             :: am_I_Root
!
! !REVISION HISTORY: 
!  27 Mar 2014 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: FILENAME, DBGMSG
      ! For NetCDF
      INTEGER                        :: fId
      INTEGER                        :: st1d(1), ct1d(1)   ! Start + count, for 1D arrays 
      REAL*8                         :: ncTemp(1)

      !=================================================================
      ! GET_RCP_NCDF begins here!
      !=================================================================

      ! Read direct from master file
      FILENAME = TRIM(RCP_FILE_NCDF)

      IF (Input_Opt%LPRT .and. am_I_Root) THEN
         WRITE(DBGMSG,'(a,a,a,a)') ' ### UCX: Reading ', 
     &        TRIM(TRAC_NAME), ' from ', TRIM( FILENAME )
         CALL DEBUG_MSG( TRIM(DBGMSG) )
      ENDIF

      CALL NcOp_Rd (fId,TRIM(FILENAME))

      ! Start and count indices - data is only yearly
      ! May attempt interpolation at a later date
      st1d(1) = TARG_YR_NC-1764
      ct1d(1) = 1

      ! Data is stored as v/v, so no adjustment necessary
      ! Seems to be having difficulties with overloaded function
      CALL NcRd( ncTemp, fId, TRIM(TRAC_NAME), st1d, ct1d )
      MEAN_MIXRATIO = ncTemp(1)
      CALL NcCl( fId )

      ! Return to calling program
      END FUNCTION GET_RCP_NCDF
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_wmo_ncdf
!
! !DESCRIPTION: Function GET\_WMO\_NCDF gets mean surface mixing ratio
!  of a given species based on forecasts from WMO estimates.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_WMO_NCDF (TRAC_NAME,TARG_YR_NC,TARG_MO_NC,Input_Opt,
     &                     am_I_Root)  RESULT (MEAN_MIXRATIO)
!
! !USES:
!
      USE FILE_MOD,           ONLY : IOERROR
      USE ERROR_MOD,          ONLY : ERROR_STOP, ALLOC_ERR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      REAL*8                         :: MEAN_MIXRATIO
      CHARACTER*(*),INTENT(IN)       :: TRAC_NAME
      INTEGER,INTENT(IN)             :: TARG_YR_NC
      INTEGER,INTENT(IN)             :: TARG_MO_NC
      TYPE(OptInput), INTENT(IN)     :: Input_Opt   ! Input options
      LOGICAL,INTENT(IN)             :: am_I_Root
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!  10 Oct 2013 - S. D. Eastham - Split into two routines (GET and SET)
!  04 Mar 2014 - S. D. Eastham - Added NetCDF compatibility
!  27 Mar 2014 - S. D. Eastham - Restructured into separate module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: FILENAME, DBGMSG
      ! For NetCDF
      INTEGER                        :: fId
      INTEGER                        :: st1d(1), ct1d(1)   ! Start + count, for 1D arrays 
      REAL*8                         :: ncTemp(1)

      !=================================================================
      ! GET_WMO_NCDF begins here!
      !=================================================================

      ! Read direct from master file
      FILENAME = TRIM(WMO_FILE_NCDF)

      IF (Input_Opt%LPRT .and. am_I_Root) THEN
         WRITE(DBGMSG,'(a,a,a,a)') ' ### UCX: Reading ', 
     &        TRIM(TRAC_NAME), ' from ', TRIM( FILENAME )
         CALL DEBUG_MSG( TRIM(DBGMSG) )
      ENDIF

      ! By the end of this, need to have data in the form
      ! of an array size JJPARx51 [LAT x LEV]

      ! UCX_MR:      JJPARx51
      ! UCX_MR_IN:   19x51
      CALL NcOp_Rd (fId,TRIM(FILENAME))

      ! Start and count indices
      st1d(1) = (12*(TARG_YR_NC-1959)) + TARG_MO_NC
      ct1d(1) = 1

      ! Data is stored as v/v, so no adjustment necessary
      ! Seems to be having difficulties with overloaded function
      CALL NcRd( ncTemp, fId, TRIM(TRAC_NAME), st1d, ct1d )
      MEAN_MIXRATIO = ncTemp(1)
      CALL NcCl( fId )

      ! Return to calling program
      END FUNCTION GET_WMO_NCDF
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lls
!
! !DESCRIPTION: Subroutine INIT\_LLS initializes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_LLS( am_I_Root, Input_Opt )
!
! !USES:
!
      USE TRACERID_MOD
      USE TIME_MOD,           ONLY : GET_YEAR
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE DIRECTORY_MOD,      ONLY : DATA_DIR_1x1
      USE RCP_MOD,            ONLY : RCPNAME, RCPYEAR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)        :: am_I_Root  ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
!
! !OUTPUT VARIABLES:
!
!
! !REVISION HISTORY: 
!  27 Mar 2014 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N, AS
      CHARACTER(LEN=255) :: DBGMSG
      LOGICAL :: LPRT, LCH4SBC, LCFCEMIS, LFUTURECFC, LCLEMIS
      LOGICAL :: LOCSEMIS, LN2OEMIS, LBREMIS

      !=================================================================
      ! INIT_LLS begins here!
      !=================================================================

      ! Copy logicals from input options structure
      LPRT         = Input_Opt%LPRT
      LFUTURECFC   = Input_Opt%LFUTURECFC
      LCH4SBC      = Input_Opt%LCH4SBC
      LCFCEMIS     = Input_Opt%LCFCEMIS
      LCLEMIS      = Input_Opt%LCLEMIS
      LOCSEMIS     = Input_Opt%LOCSEMIS
      LN2OEMIS     = Input_Opt%LN2OEMIS
      LBREMIS      = Input_Opt%LBREMIS

      WRITE( 6,'(a)') REPEAT( '=', 79 )
      WRITE( 6,'(a)') 'L O N G - L I V E D   S P E C I E S'
      WRITE( 6,'(a)') 'Routines written by SEBASTIAN D. EASTHAM'
      WRITE( 6,'(a)') REPEAT( '=', 79 )

      WRITE(WMO_FILE_NCDF,'(a,a)') TRIM(DATA_DIR_1x1), 
     &   '/UCX_201403/Init1D/Montreal.nc'
      IF (Input_Opt%LRCPLLS) THEN
         WRITE(RCP_FILE_NCDF,'(a,a,a,a)') TRIM(DATA_DIR_1x1), 
     &      '/UCX_201403/Init1D/',TRIM(RCPNAME),'.nc'
         ! Overwrite LLSYEAR with RCPYEAR
         LLSYEAR = RCPYEAR
      ELSE
         IF (.not.LFUTURECFC) THEN
            LLSYEAR = GET_YEAR()
         ENDIF
         RCP_FILE_NCDF=''
      ENDIF
  
      IF ( LPRT ) THEN
         WRITE(DBGMSG,'(a,a)') 
     &     '### UCX: Setting long-lived species BCs from ',
     &     TRIM(WMO_FILE_NCDF)
         CALL DEBUG_MSG( TRIM(DBGMSG) )
      ENDIF

      ! Set "current" month to impossible value
      SFC_MONTH = 0
      OFFSETYEAR = 0

      ! Determine which tracers are to be emitted
      N_EMIT = 0
      IDT_EMIT(:) = 0

      IF (LCH4SBC) THEN
         ! Simple CH4 fixed surface emissions
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH4
      ENDIF 

      IF (LCFCEMIS) THEN
         ! R-11/12/113/114/115
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFC11

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFC12

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFCX

         ! HCFC-22/123/141b/142b
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTHCFC22

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTHCFCX

         ! Halons (H-1301/2402/1211)
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH1301

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH1211

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH2402
      ENDIF

      IF (LCLEMIS) THEN
         ! CCl4
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCCl4

         ! CH3Cl
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3Cl

         ! CH3CCl3
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3CCl3
      ENDIF

      IF (LOCSEMIS) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTOCS
      ENDIF

      IF (LN2OEMIS) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTN2O
      ENDIF 

      IF (LBREMIS) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3Br
      ENDIF

      IF (IDTH2.gt.0) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH2
      ENDIF

      IF (N_EMIT.gt.0) THEN
         ALLOCATE( GRID_EMIT( JJPAR, N_EMIT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GRID_EMIT' )
         GRID_EMIT = 0d0
      ENDIF

      END SUBROUTINE INIT_LLS

!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_lls
!
! !DESCRIPTION: Subroutine CLEANUP\_LLS deallocates module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_LLS
!
! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_LLS begins here!
      !=================================================================

      IF ( ALLOCATED( GRID_EMIT  ) ) DEALLOCATE( GRID_EMIT  )

      END SUBROUTINE CLEANUP_LLS
!
!EOC
      END MODULE LONGLIVED_MOD
