!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: aeic_mod
!
! !DESCRIPTION: Module AEIC\_MOD contains variables and routines for aircraft
!  flight emissions into the chemistry and transport grids. (sde, 12/14/12)
!\\
!\\
! !INTERFACE: 
!
      MODULE AEIC_MOD
!
! !USES:
!
      IMPLICIT NONE
#     include "netcdf.inc" ! Needed for netCDF libraries
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: 
!
      PUBLIC EMIT_AEIC
      PUBLIC CLEANUP_AEIC
!
! 
! !REMARKS:
!  The AEIC fields are stored on 1x1 grids.
!  These fields will be interpolated onto the current horizontal grid.
!                                                                             .
!  References:
!  ============================================================================
!  (1) Stettler, M.E.J., S. Eastham, S.R.H. Barrett, “Air quality and public
!       health impacts of UK airports. Part I: Emissions,” Atmos. Env., 2011.
!       DOI: 10.1016/j.atmosenv.2011.07.012
! 
! !REVISION HISTORY: 
!  14 Dec 2012 - S. Eastham  - Adapted from AIRCRAFT_NOX_MOD. Now emits
!                              aircraft NOx, CO, HC, SO2, SO4, OC and BC
!  01 Aug 2013 - M. Sulprizio- Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Number of levels on the 1x1 grid
      ! Standard resolution (~3 GB of data)
      INTEGER, PARAMETER   :: NLVL_HD     = 5
      REAL*8,  PARAMETER   :: DEPTHLVL_HD = 60.96d0
      INTEGER, PARAMETER   :: NLVL_SD     = 44
      REAL*8,  PARAMETER   :: DEPTHLVL_SD = 304.8d0
      INTEGER, PARAMETER   :: NLVLAEIC    = 49
 
      ! Indices for emitted species
      INTEGER, PARAMETER   :: FB_INDEX    = 1
      INTEGER, PARAMETER   :: HC_INDEX    = 2
      INTEGER, PARAMETER   :: CO_INDEX    = 3
      INTEGER, PARAMETER   :: NOx_INDEX   = 4
      INTEGER, PARAMETER   :: BC_INDEX    = 5
      INTEGER, PARAMETER   :: OC_INDEX    = 6
  
      ! THC to TOG conversion factor for aircraft emissions
      REAL,    PARAMETER   :: THC2TOG     = 1.16d0

      ! HC speciation is included below
      REAL,    PARAMETER   :: ACF_ACET    = 0.003693477d0
      REAL,    PARAMETER   :: ACF_ALD2    = 0.042718224d0
      REAL,    PARAMETER   :: ACF_ALK4    = 0.213791063d0
      REAL,    PARAMETER   :: ACF_C2H6    = 0.005214505d0
      REAL,    PARAMETER   :: ACF_C3H8    = 0.000780871d0
      REAL,    PARAMETER   :: ACF_CH2O    = 0.123081099d0
      REAL,    PARAMETER   :: ACF_PRPE    = 0.178041756d0
      REAL,    PARAMETER   :: ACF_MACR    = 0.005362609d0
      REAL,    PARAMETER   :: ACF_RCHO    = 0.036769436d0
      ! Note not all aircraft hydrocarbon species modeled in GEOS-Chem

      ! Fuel sulfur properties
      REAL                 :: FSC                    ! Fraction by mass
      REAL                 :: SULFCONV               ! Conversion efficiency
 
!
! !PRIVATE TYPES:
!
      ! Emissions on native 1x1 grid
      REAL*8,  ALLOCATABLE :: AEIC_1x1(:,:,:,:)

      ! Regridded emissions
      REAL*8,  ALLOCATABLE :: AEIC_REGRID(:,:,:,:)
 
      ! Aircraft emissions on GEOS-CHEM grid
      REAL*8,  ALLOCATABLE :: AEIC_EMIT(:,:,:,:)

      ! LTO ceiling level (3000 ft = 914.4 m)
      REAL*8               :: LTOCEIL     = 9.144d2

      ! For AEDT emissions
      ! NAIR     - Maximum number for aircraft emissions in SAGE  
#if defined( GRIDREDUCED )
      INTEGER, PARAMETER :: NAIR = 40
#else
      INTEGER, PARAMETER :: NAIR = 70
#endif

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_aedt_file
!
! !DESCRIPTION: Subroutine READ\_AEDT\_FILE reads fuel burn and
!  emissions from AEDT text files (daily average) (SDE, 2013-11-07)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_AEDT_FILE( INMONTH, INDAY, State_Met, Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY   : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE FILE_MOD,       ONLY : IU_FILE, IOERROR
      USE TIME_MOD,       ONLY : GET_YEAR!, GET_HOUR
      USE TIME_MOD,       ONLY : GET_MONTH
      USE TIME_MOD,       ONLY : GET_TS_EMIS!, GET_HOUR
      USE ERROR_MOD,      ONLY : ERROR_STOP, ALLOC_ERR
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1 
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE CMN_SIZE_MOD
      USE CMN_GCTM_MOD, ONLY : Rdg0
 
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
      INTEGER,INTENT(IN) :: INMONTH, INDAY 
      CHARACTER(LEN=2)   :: MONTH_STR, DAY_STR
 
      ! From FAA emissions database
      INTEGER               :: J, I, K, N, IOS
      INTEGER               :: J0, I0
      INTEGER               :: IREF, JREF, L
      REAL*8                :: FB, CO, HC, NOx, PMNV, PMFO, MTEMP
      REAL*8                :: PLOW, PHIGH, XSUM, PAIR1, PAIR2, TMP

      ! SDE 06/24/12: New emissions code
      REAL*8                :: cellfrac,maxalt
      REAL*8                :: emzdelta,emzmin,emzmax,ltofrac
      REAL*8                :: ltoceil,prevalt,pnext,pseudoground
      REAL*8                :: actmax, actmin
      REAL*8,PARAMETER      :: psealev = 1013.25d0
      REAL*8, PARAMETER     :: ltoalt = 914.4d0
      LOGICAL               :: evallto
      INTEGER               :: AS
      REAL*8, dimension(:), allocatable :: emiss_temp,emiss_mini
      
      CHARACTER(LEN=255)    :: FILENAME
      CHARACTER(LEN=255)    :: EMISS_DIR
      REAL*8                :: LTO_FB_MULT, CRUISE_FB_MULT
      REAL*8,DIMENSION(6)   :: LTOSCALE,CRUISESCALE,EMITSCALE
      
      !=================================================================
      ! READ_AEDT_FILE begins here!
      !=================================================================
!
      ! Re-zero output grid
      AEIC_EMIT = 0d0

      ! Construct filename
      WRITE( MONTH_STR, '(i0.2)' )  INMONTH
      WRITE( DAY_STR, '(i0.2)' )  INDAY

      ! SDE: If working on a leap year, repeat Feb 28th for Feb 29th
      IF (( INDAY .EQ. 29) .AND. ( INMONTH .EQ. 2 )) THEN
         WRITE( DAY_STR, '(i0.2)' ) 28
      ENDIF

      ! Allocate temporary emission vector
      ! Fuel burn, CO, HC, NOx, PMNV, PMFO, NO2, HONO
      ALLOCATE( emiss_temp(8) , STAT=AS) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'emiss_temp' )
      ALLOCATE( emiss_mini(8) , STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'emiss_mini' )
      emiss_temp = 0d0
      emiss_mini = 0d0

      EMISS_DIR = Input_Opt%AEDT_DIR
!      WRITE(EMISS_DIR,'(a)') '/net/d03/data/seb_data/AEDT_Regrid_4x5/'
      FILENAME = TRIM(EMISS_DIR) // '/LAE_ACEmit'
      FILENAME = TRIM(FILENAME) // '_2006_' 
     &    // TRIM(MONTH_STR) // '_' // TRIM(DAY_STR)
     &    // '.txt'
      PRINT *, ''
      PRINT *,'************************************************'
      PRINT *,'Reading aircraft emissions: ' // TRIM(FILENAME)

      ! Do we include LTO and/or cruise?
      IF ( .NOT. Input_Opt%LLTO_EMIS) THEN
         LTO_FB_MULT = 0.0d0
      ENDIF
      IF ( .NOT. Input_Opt%LCRUISE_EMIS ) THEN
         CRUISE_FB_MULT = 0.0d0
      END IF ! Else all emissions included

      ! Now open file
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS)
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'RCE:1' )
 
      DO 
         !  Read aircraft emissions
         !LAT_IDX, LON_IDX, ALT_MIN, ALT_MAX,
         !FB, CO, HC, NOx, PMNV, PMFO
         READ( IU_FILE, *, IOSTAT=IOS ) J,
     &   I, emzmin, emzmax,
     &   emiss_temp(1), emiss_temp(2),
     &   emiss_temp(3), emiss_temp(4),
     &   emiss_temp(5), emiss_temp(6)
         ! Convert emissions g/s to kg/s
         emiss_temp(2:6) = emiss_temp(2:6) * 1.0d-3
         
         ! IOS < 0 is end of file
         ! IOS > 0 is an I/O error
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'RCE:2' )

         ! Determine emission cells 
         ! Altitude readings
         maxalt = 0.0d0             ! Output ceiling
         prevalt = 0.0d0            ! Output base
         ! LTO ceiling is local ground plus 3000 ft
         actmin = 0.0d0               ! Effective base
         actmax = 0.0d0               ! Effective ceiling

         ! Fractional multipliers
         ltofrac = 0.0d0  ! Emission fraction in LTO
         cellfrac = 0.0d0 ! Fraction of emission in cell
        
         ! Calculate LTO fraction?
         evallto = .true.

         ! LTO ceiling taken as 3000 ft (914.4 m)
         ! above local ground; need to define in
         ! reference to sea level
         pnext = GET_PEDGE(I,J,1)
         pseudoground = (Rdg0 * State_Met%T(I,J,1) * LOG(psealev/pnext))
         if (psealev .le. pnext) then
            ltoceil = ltoalt
         else
            ltoceil =  pseudoground + ltoalt
         endif

         ! Move the base of the emission to be at least at
         ! ground altitude
         if (emzmin .le. pseudoground) then
            emzmin = pseudoground
            if (emzmax .le. pseudoground) then
               emzmax = emzmin + (0.1*State_Met%BXHEIGHT(I,J,1))
            endif
         endif            

         emzdelta = emzmax - emzmin
         if (emzmin .ge. ltoceil) then
            ! 100% cruise emissioni
            evallto = .false.
            ltofrac = 0.0d0
         elseif (emzmax .le. ltoceil) then
            ! 100% LTO emission
            evallto = .false.
            ltofrac = 1.0d0
         else
            evallto = .true.
         endif
 
         maxalt = pseudoground
         loop_lmax: do K=1,NAIR
            ! Cell limits: prevalt -> maxalt
            prevalt = maxalt
            maxalt = prevalt + State_Met%BXHEIGHT(i,j,k)
            ! Check if emission cell below sea/ground level
            if (maxalt .le. max(0.0d0,pseudoground)) then
               cycle loop_lmax
            endif
            ! Reached emission cell?
            if (emzmin .le. maxalt) then
               ! Clear temporary emissions
               emiss_mini = 0d0
               ! Establish temporary emission limits
               actmin = max(emzmin,prevalt)
               actmax = min(emzmax,maxalt)
               ! Calculate fraction of emission to apply
               cellfrac = (actmax-actmin)/emzdelta
               if (evallto) then
                  if (maxalt .le. ltoceil) then
                     ! Cell entirely within LTO
                     ltofrac = 1.0d0
                  else
                     ! Border cell discovered
                     ! Base is LTO, top is not
                     ! Stop LTO evaluation after this
                     evallto = .false.
                     ! ltofrac is the LTO fraction of the emission
                     ltofrac = ltoceil - actmin
                     ltofrac = ltofrac/(actmax-actmin)
                  endif
               endif
   
               ! Calculate LTO emissions
               emiss_mini(:) = emiss_temp(:) * ltofrac * 
     &            Input_Opt%LTO_FB_MULT
               emiss_mini(2) = emiss_mini(2) * Input_Opt%LTO_CO_MULT
               emiss_mini(3) = emiss_mini(3) * Input_Opt%LTO_HC_MULT
               emiss_mini(4) = emiss_mini(4) * Input_Opt%LTO_NOx_MULT
               ! Keep BC and OC as calculated for now
               emiss_mini(7) = 0d0 ! Ignore NO2 fraction
               emiss_mini(8) = 0d0 ! Ignore H2O
               if (actmax .ge. ltoceil) then
                  ! Add cruise emissions
                  emiss_mini(1:4) = emiss_mini(1:4) + ((1.0d0-ltofrac) 
     &               * emiss_temp(1:4) * Input_Opt%CRUISE_FB_MULT)
                  emiss_mini(5) = emiss_mini(5) + ((1.0d0-ltofrac) 
     &              * emiss_temp(1) * Input_Opt%CRUISE_BC_EI * 1.0d-3
     &              * Input_Opt%CRUISE_FB_MULT)
                  emiss_mini(6) = emiss_mini(6) + ((1.0d0-ltofrac) 
     &              * emiss_temp(1) * Input_Opt%CRUISE_OC_EI * 1.0d-3
     &              * Input_Opt%CRUISE_FB_MULT)
                  ! Ignore H2O/NO2 fractions
               endif
               emiss_mini(:) = emiss_mini(:) * cellfrac
              
               AEIC_EMIT(I,J,K,FB_INDEX) = AEIC_EMIT(I,J,K,FB_INDEX) +
     &                                       emiss_mini(1)
               AEIC_EMIT(I,J,K,CO_INDEX) = AEIC_EMIT(I,J,K,CO_INDEX) +
     &                                       emiss_mini(2)
               AEIC_EMIT(I,J,K,HC_INDEX) = AEIC_EMIT(I,J,K,HC_INDEX) +
     &                                       emiss_mini(3)
               AEIC_EMIT(I,J,K,NOx_INDEX) = AEIC_EMIT(I,J,K,NOx_INDEX) +
     &                                       emiss_mini(4)
               AEIC_EMIT(I,J,K,BC_INDEX) = AEIC_EMIT(I,J,K,BC_INDEX) +
     &                                       emiss_mini(5)
               AEIC_EMIT(I,J,K,OC_INDEX) = AEIC_EMIT(I,J,K,OC_INDEX) +
     &                                       emiss_mini(6)
            endif
            if (emzmax .le. maxalt) then
               ! Last cell
               exit loop_lmax
            endif
         enddo loop_lmax
      ENDDO
      ! Get nested-grid offsets
!     I0 = GET_XOFFSET()
!     J0 = GET_YOFFSET()
      
      CLOSE(IU_FILE)
  
      IF ( ALLOCATED( emiss_temp )) THEN
         DEALLOCATE( emiss_temp )
      ENDIF
      IF ( ALLOCATED( emiss_mini )) THEN
         DEALLOCATE( emiss_mini )
      ENDIF
      
      ! Apply exclusion zones
      IF (Input_Opt%EXC_CODE .eq. 0) THEN
         ! Do nothing
      ELSE IF ( Input_Opt%EXC_CODE .eq. 1 ) THEN
         ! Only use one zone
         AEIC_EMIT(1:(Input_Opt%EXC_LL_I-1),:,:,:) = 0.0d0
         AEIC_EMIT((Input_Opt%EXC_UR_I+1):IIPAR,:,:,:) = 0.0d0
         AEIC_EMIT(:,1:(Input_Opt%EXC_LL_J-1),:,:) = 0.0d0
         AEIC_EMIT(:,(Input_Opt%EXC_UR_J+1):JJPAR,:,:) = 0.0d0
      ELSE IF (Input_Opt%EXC_CODE .eq. -1 ) THEN
         ! Exclude one zone
         AEIC_EMIT(Input_Opt%EXC_LL_I:Input_Opt%EXC_UR_I,
     &         Input_Opt%EXC_LL_J:Input_Opt%EXC_UR_J,:,:) = 0.0d0
      ELSE
         ! Invalid
         CALL ERROR_STOP( 
     &      'Aircraft emissions exlusion flag invalid!',
     &      'READ_AEDT_FILE ("aeic_mod.F")' )
      END IF

      ! Return to calling program
      END SUBROUTINE READ_AEDT_FILE
!
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_aeic_file
!
! !DESCRIPTION: Subroutine READ\_AEIC\_FILE reads fuel burn and
!  emissions from AEIC netCDF files (monthly average). (sde, 12/14/12)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_AEIC_FILE( INMONTH, State_Met, Input_Opt )
!
! !USES:
!
      ! NcdfUtil modules for netCDF I/O
      USE m_netcdf_io_open         ! netCDF open
      USE m_netcdf_io_get_dimlen   ! netCDF dimension queries
      USE m_netcdf_io_read         ! netCDF data reads
      USE m_netcdf_io_close        ! netCDF close

      USE DIRECTORY_MOD,      ONLY : DATA_DIR_1x1
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE REGRID_A2A_MOD,     ONLY : DO_REGRID_A2A
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE PRESSURE_MOD,       ONLY : GET_AP, GET_BP

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_GCTM_MOD,       ONLY : Rdg0
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: INMONTH     ! Current month number (1-12)
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY: 
!  14 Dec 2012 - S. Eastham  - Adapted from READ_NONERUP_VOLC
!  31 Jul 2013 - M. Sulprizio- Now pass met fields using State_Met
!  01 Aug 2013 - M. Sulprizio- Update to read emissions from netCDF file
!                              created by Christoph Keller. He combined original
!                              monthly files into one file and made file
!                              COARDS compliant.
!  26 Aug 2013 - R. Yantosca - Avoid array temporaries in call to NcRd
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8, TARGET       :: ACARRAY_IN(360,180,4)
      REAL*8, TARGET       :: ACARRAY_OUT(IIPAR,JJPAR)
      INTEGER              :: I,J,L, AEIC_FID, EMIT_INDEX
      CHARACTER(LEN=255)   :: AEIC_FILE,MONTH_STR
      CHARACTER(LEN=255)   :: LLFILENAME, EMIT_NAME
      CHARACTER(LEN=255)   :: AEIC_GRIDSTR
      REAL*8, POINTER      :: OUTGRID(:,:) => NULL()
      REAL*8, POINTER      :: INGRID(:,:)  => NULL()
      REAL*8               :: PBASE, PCEIL
      REAL*8               :: CURRBASE, CURRCEIL
      REAL*8,PARAMETER     :: PSEALEV=1013.25d0
      REAL*8,DIMENSION(6)  :: LTOSCALE,CRUISESCALE,EMITSCALE
      REAL*8               :: TEMP ! Local temperature
      INTEGER              :: L_AEIC
      REAL*8               :: FRAC_AEIC,AEIC_CEIL,AEIC_BASE,ACTIVEDEPTH
      REAL*8               :: DEPTHLVLAEIC, LTOFRAC, PLUSALT, MINALT
      LOGICAL              :: AEIC_DONE
      INTEGER              :: st4d(4), ct4d(4)
      LOGICAL              :: INCRUISE

      !=================================================================
      ! READ_AEIC begins here!
      !=================================================================

      ! Retrieve values from options structure and convert from km to m
      PLUSALT = Input_Opt%ALT_OFFSET * 1.0d3
      MINALT = (Input_Opt%ALT_FLR * 1.0d3) + PLUSALT

      ! Filename
      AEIC_FILE = TRIM( DATA_DIR_1x1 ) // '/AEIC_201301/' //
     &           'AEIC.nc'

      ! Echo output
      WRITE( 6, 100 ) TRIM( AEIC_FILE )
 100  FORMAT( '     - READ_AEIC: Reading ', a )

      ! Read 3D NetCDF file containing emissions data for this month
      ! Open file
      CALL NcOp_Rd( AEIC_FID, TRIM(AEIC_FILE))

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1 ) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      DO EMIT_INDEX = 1, 4
  
         SELECT CASE (EMIT_INDEX)
            CASE (FB_INDEX)
               EMIT_NAME = 'fuelBurn'
            CASE (CO_INDEX)
               EMIT_NAME = 'CO'
            CASE (HC_INDEX)
               EMIT_NAME = 'HC'
            CASE (NOx_INDEX)
               EMIT_NAME = 'NOx'
            CASE DEFAULT
               CALL ERROR_STOP( 'Unrecognized aircraft emission!',
     &               'READ_AEIC (aeic_mod.F90)' )
         END SELECT
 
         st4d = (/   1,   1,        1 , INMONTH /)
         ct4d = (/ 360, 180, NLVLAEIC,        1 /)
         CALL NcRd( AEIC_1x1(:,:,:,EMIT_INDEX),
     &              AEIC_FID,TRIM(EMIT_NAME), st4d, ct4d )
 
      ENDDO

      ! Close NetCDF file
      CALL NcCl(AEIC_FID)

      ! Split out LTO and cruise emissions
      LTOSCALE(FB_INDEX) = Input_Opt%LTO_FB_MULT
      LTOSCALE(HC_INDEX) = Input_Opt%LTO_FB_MULT*
     &      Input_Opt%LTO_HC_MULT*1.0d-3
      LTOSCALE(NOx_INDEX)= Input_Opt%LTO_FB_MULT*
     &      Input_Opt%LTO_NOX_MULT*1.0d-3
      LTOSCALE(CO_INDEX) = Input_Opt%LTO_FB_MULT*
     &      Input_Opt%LTO_CO_MULT*1.0d-3
      LTOSCALE(BC_INDEX) = Input_Opt%LTO_FB_MULT*
     &      Input_Opt%LTO_BC_EI*1.0d-3
      LTOSCALE(OC_INDEX) = Input_Opt%LTO_FB_MULT*
     &      Input_Opt%LTO_OC_EI*1.0d-3
      CRUISESCALE(FB_INDEX) = Input_Opt%CRUISE_FB_MULT
      CRUISESCALE(HC_INDEX) = Input_Opt%CRUISE_FB_MULT*
     &      Input_Opt%CRUISE_HC_MULT*1.0d-3
      CRUISESCALE(NOx_INDEX)= Input_Opt%CRUISE_FB_MULT*
     &      Input_Opt%CRUISE_NOX_MULT*1.0d-3
      CRUISESCALE(CO_INDEX) = Input_Opt%CRUISE_FB_MULT*
     &      Input_Opt%CRUISE_CO_MULT*1.0d-3
      CRUISESCALE(BC_INDEX) = Input_Opt%CRUISE_FB_MULT*
     &      Input_Opt%CRUISE_BC_EI*1.0d-3
      CRUISESCALE(OC_INDEX) = Input_Opt%CRUISE_FB_MULT*
     &      Input_Opt%CRUISE_OC_EI*1.0d-3

      ! First step: rescale grid horizontally at each level
      DO L=1,NLVLAEIC
         ACARRAY_IN = AEIC_1x1(:,:,L,:)
         DO EMIT_INDEX=1,4
            ! Point to input & output arrays
            INGRID  => ACARRAY_IN(:,:,EMIT_INDEX) !Units kg/s
            OUTGRID => ACARRAY_OUT

            ! Regrid
            CALL DO_REGRID_A2A( LLFILENAME, 360,     180,
     &                          INGRID,     OUTGRID, IS_MASS=1,
     &                          netCDF=.TRUE.                   )

            AEIC_REGRID(:,:,L,EMIT_INDEX) = OUTGRID

            ! Free pointers
            NULLIFY( INGRID, OUTGRID )
         ENDDO
      ENDDO

      ! Re-zero output grid
      AEIC_EMIT = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, EMIT_INDEX )
!$OMP+PRIVATE( PBASE, PCEIL, TEMP, CURRBASE, CURRCEIL, AEIC_DONE )
!$OMP+PRIVATE( ACARRAY_IN, DEPTHLVLAEIC, AEIC_CEIL, AEIC_BASE )
!$OMP+PRIVATE( FRAC_AEIC, ACTIVEDEPTH, L_AEIC, EMITSCALE )
!$OMP+PRIVATE( LTOFRAC, INCRUISE )
      DO J=1,JGLOB
      DO I=1,IGLOB
         CURRCEIL    = 0d0
         PCEIL       = PSEALEV
         L_AEIC      = 1
         INCRUISE    = .FALSE.
         AEIC_CEIL   = PLUSALT
         ACARRAY_IN  = 0.0d0
         AEIC_DONE   = .FALSE.
         EMITSCALE   = LTOSCALE

         DO L=1,LLPAR
            ! Are we in LTO or cruise (or fraction?)
            CURRBASE = CURRCEIL
            PBASE    = PCEIL

            ! Calculate ceiling of next box on standard grid
            PCEIL    = GET_AP(L+1) + (GET_BP(L+1) * PSEALEV)

            ! Temperature [K]
            TEMP     = State_Met%T(I,J,L)

            CURRCEIL = CURRBASE + (Rdg0*TEMP*LOG(PBASE/PCEIL))
 
            IF ((.not. INCRUISE).and.(LTOCEIL .lt. CURRCEIL)) THEN
               ! Transition box
               INCRUISE = .TRUE.
               LTOFRAC = (LTOCEIL-CURRBASE)/(CURRCEIL-CURRBASE)
               EMITSCALE = LTOSCALE*LTOFRAC
               EMITSCALE = EMITSCALE+(CRUISESCALE*(1.0d0-LTOFRAC))
            ELSEIF (INCRUISE) THEN
               EMITSCALE = CRUISESCALE
            ELSE
               EMITSCALE = LTOSCALE
            ENDIF

            ! Scan through layers until we hit a GC boundary
            AEIC_DONE  = ((AEIC_CEIL.le.CURRCEIL) .and. 
     &                    (L_AEIC.gt.NLVLAEIC))
            ACARRAY_IN = 0d0
            
            DO WHILE (.not. AEIC_DONE)
               ! Calculate fraction of AEIC level within GC level
               IF (L_AEIC .le. NLVL_HD) THEN
                  ! In finely resolved area
                  DEPTHLVLAEIC = DEPTHLVL_HD
                  AEIC_CEIL = (REAL(L_AEIC)*DEPTHLVLAEIC)
               ELSE
                  ! Standard resolution region
                  DEPTHLVLAEIC = DEPTHLVL_SD
                  AEIC_CEIL = (REAL(NLVL_HD)*DEPTHLVL_HD) +
     &                        (REAL(L_AEIC-NLVL_HD)*DEPTHLVL_SD)
               ENDIF

               AEIC_BASE = AEIC_CEIL-DEPTHLVLAEIC
               ! Establish base level for input cell
               IF (AEIC_CEIL .gt. MINALT) THEN
                   IF (AEIC_BASE .gt. MINALT) THEN
                      FRAC_AEIC = MAX(CURRBASE,AEIC_BASE)
                      ACTIVEDEPTH = DEPTHLVLAEIC
                   ELSE
                      FRAC_AEIC = MAX(CURRBASE,AEIC_BASE,MINALT)
                      ACTIVEDEPTH=AEIC_CEIL-MAX(MINALT,AEIC_BASE)
                   ENDIF

                   IF (AEIC_CEIL.lt.CURRCEIL) THEN
                      ! Switch input level
                      FRAC_AEIC = (AEIC_CEIL-FRAC_AEIC)/ACTIVEDEPTH
                   ELSE
                      ! Switch output level
                      FRAC_AEIC = (CURRCEIL-FRAC_AEIC)/ACTIVEDEPTH
                   ENDIF
               ELSE
                   ! Below cutoff
                   FRAC_AEIC = 0.0d0
               ENDIF
                      

               ! Add to output grid
               DO EMIT_INDEX=1,4
                  AEIC_EMIT(I,J,L,EMIT_INDEX) = 
     &                 AEIC_EMIT(I,J,L,EMIT_INDEX) + 
     &                (AEIC_REGRID(I,J,L_AEIC,EMIT_INDEX)*FRAC_AEIC)
               ENDDO

               IF (AEIC_CEIL.lt.CURRCEIL) THEN
                  L_AEIC = L_AEIC + 1
                  IF (L_AEIC.gt.NLVLAEIC) THEN
                     AEIC_DONE = .TRUE.
                  ENDIF
               ELSE
                  ! Change GC level now
                  AEIC_DONE = .TRUE.
               ENDIF
            ENDDO

            ! Copy OC/BC from fuelburn, apply EIs in cruise/LTO
            AEIC_EMIT(I,J,L,OC_INDEX) = AEIC_EMIT(I,J,L,FB_INDEX)
            AEIC_EMIT(I,J,L,BC_INDEX) = AEIC_EMIT(I,J,L,FB_INDEX)

            ! Store emissions
            AEIC_EMIT(I,J,L,:) = EMITSCALE*AEIC_EMIT(I,J,L,:)
            
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
  
      ! Apply exclusion zones
      IF (Input_Opt%EXC_CODE .eq. 0) THEN
         ! Do nothing
      ELSE IF ( Input_Opt%EXC_CODE .eq. 1 ) THEN
         ! Only use one zone
         AEIC_EMIT(1:(Input_Opt%EXC_LL_I-1),:,:,:) = 0.0d0
         AEIC_EMIT((Input_Opt%EXC_UR_I+1):IIPAR,:,:,:) = 0.0d0
         AEIC_EMIT(:,1:(Input_Opt%EXC_LL_J-1),:,:) = 0.0d0
         AEIC_EMIT(:,(Input_Opt%EXC_UR_J+1):JJPAR,:,:) = 0.0d0
      ELSE IF (Input_Opt%EXC_CODE .eq. -1 ) THEN
         ! Exclude one zone
         AEIC_EMIT(Input_Opt%EXC_LL_I:Input_Opt%EXC_UR_I,
     &         Input_Opt%EXC_LL_J:Input_Opt%EXC_UR_J,:,:) = 0.0d0
      ELSE
         ! Invalid
         CALL ERROR_STOP( 
     &      'Aircraft emissions exlusion flag invalid!',
     &      'READ_AEIC_FILE ("aeic_mod.F")' )
      END IF

      END SUBROUTINE READ_AEIC_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emit_aeic
!
! !DESCRIPTION: Subroutine EMIT\_AEIC interpolates AEIC aircraft emissions from
!  the native grid onto the given GEOS-CHEM grid. (sde, 12/14/12)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMIT_AEIC( Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY   : GET_CHEMGRID_LEVEL
      USE DIAG_MOD,           ONLY   : AD32_ac,     AD13_SO2_ac
      USE GIGC_Input_Opt_Mod, ONLY   : OptInput
      USE GIGC_State_Chm_Mod, ONLY   : ChmState
      USE GIGC_State_Met_Mod, ONLY   : MetState
      USE GRID_MOD,           ONLY   : GET_XOFFSET, GET_YOFFSET
      USE GRID_MOD,           ONLY   : GET_AREA_CM2
      USE PRESSURE_MOD,       ONLY   : GET_AP,      GET_BP
      USE TIME_MOD,           ONLY   : GET_MONTH,   GET_TS_EMIS
      USE TIME_MOD,           ONLY   : GET_DAY
      USE TRACERID_MOD

      USE CMN_SIZE_MOD               ! Size parameters
      USE CMN_DIAG_MOD               ! Diagnostic switches
      USE CMN_MOD                    ! PTOP, SIGE, AVP
!
! !INPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY:
!  31 Jul 2013 - M. Sulprizio- Now add aircraft NOx and SO2 emissions to
!                              AD32 and AD13 arrays for diagnostics
!  01 Aug 2013 - M. Sulprizio- Added ProTeX headers
!  13 Dec 2013 - M. Sulprizio- Only add SO2 aircraft emissions to AD13_SO2_ac
!                              if level is <= LD13 in order to avoid array-
!                              out-of-bounds error (J. Fisher)
!  07 Feb 2014 - R. Yantosca - Bug fix: Only save into LLCHEM levels of the
!                              array AD32_nc.  This is set for UCX now.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I,  J,  IREF,    JREF,  L
      INTEGER          :: I0, J0, MAXLEV
      REAL*8           :: DTSRCE, TMPMULT, TEMP
      REAL*8           :: AREA_CM2
      INTEGER, SAVE    :: LAST_MONTH = 0
      INTEGER, SAVE    :: LAST_DAY   = 0
      LOGICAL, SAVE    :: FIRST      = .TRUE.
      REAL*8,  SAVE    :: FBSUM      = 0d0

      ! AEDT modifications
      LOGICAL          :: PRINT_AEIC
      INTEGER          :: CDAY, CMONTH

      ! For fields from Input_Opt
      LOGICAL          :: LVARTROP
      REAL*8           :: XNUMOL(Input_Opt%N_TRACERS)

      ! For fields from State_Chm
      REAL*8, POINTER  :: STT(:,:,:,:)

      !=================================================================
      ! EMIT_AEIC begins here!
      !=================================================================

      ! Copy logical fields from INPUT_OPT to local variables
      XNUMOL    = Input_Opt%XNUMOL

      ! Initialize fields from State_Chm
      STT      => State_Chm%Tracers

      IF ( FIRST ) THEN
         FIRST = .FALSE.

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   )
         WRITE( 6, 110   )
         WRITE( 6, 120   )
         WRITE( 6, 130   )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! FORMAT strings
 100     FORMAT( 'A I R C R A F T   E M I S S I O N S'   )
 110     FORMAT( 'Routines originally by STEVEN BARRETT'    )
 120     FORMAT( 'Edited by Jamin Koo and S. D. Eastham'   )
 130     FORMAT( 'Last Modification Date: 12/19/12'       )

         CALL INIT_AEIC

#if defined( UCX )
         IF ( IDTH2O .ne. 0 ) THEN
            WRITE(6,*) ' - H2O EI          : ',1260.d0
         ELSE
            WRITE(6,*) 'Aircraft H2O emissions ignored.'
         ENDIF
#endif
         ! Get conversion coefficients from inputs
         FSC = Input_Opt%FSC
         SULFCONV = Input_Opt%SULFCONV

      END IF

      ! Read aircraft emissions
      PRINT_AEIC=.FALSE.
      CDAY = GET_DAY()
      CMONTH = GET_MONTH()
      IF ((Input_Opt%LAEDT) .AND. (LAST_DAY/=CDAY)) THEN
          LAST_DAY = CDAY
          LAST_MONTH = CMONTH
          CALL READ_AEDT_FILE(LAST_MONTH,LAST_DAY,State_Met,Input_Opt)
          PRINT_AEIC=.TRUE.
         WRITE(6,'(a,I02,a,I02)') 'Aircraft emissions (AEDT) for 2006-',
     &      LAST_MONTH,'-',LAST_DAY
         IF ((LAST_MONTH.eq.1).and.(LAST_DAY.eq.1)) THEN
            FBSUM = 0d0
         ENDIF
      ELSEIF (LAST_MONTH/=CMONTH) THEN
         LAST_MONTH = CMONTH
         CALL READ_AEIC_FILE( LAST_MONTH, State_Met, Input_Opt )

         WRITE(6,*) 'Aircraft emissions (AEIC) for month ', LAST_MONTH
         PRINT_AEIC=.TRUE.
         IF (LAST_MONTH .eq. 1) THEN
            FBSUM = 0d0
         ENDIF
      ELSE
         PRINT_AEIC=.FALSE.
      ENDIF

      IF (PRINT_AEIC) THEN
         WRITE(6,*) 'Aircraft fuelburn: ',
     &       SUM(AEIC_EMIT(:,:,:,FB_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft NOx     : ',
     &       SUM(AEIC_EMIT(:,:,:,NOx_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft CO      : ',
     &       SUM(AEIC_EMIT(:,:,:,CO_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft BC/PMNV : ',
     &       SUM(AEIC_EMIT(:,:,:,BC_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft OC/PMV  : ',
     &       SUM(AEIC_EMIT(:,:,:,OC_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft HC      : ',
     &       SUM(AEIC_EMIT(:,:,:,HC_INDEX)), ' kg/s'
      ENDIF

      ! Emissions timestep
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      ! Debug
      !WRITE(6,*) 'PRESO4: ', SUM(STT(:,:,:,IDTSO4))
      !WRITE(6,*) 'PRENOx: ', SUM(STT(:,:,:,IDTNO))

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, MAXLEV, TMPMULT, AREA_CM2 )
      ! Loop over grid boxes
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         MAXLEV = GET_CHEMGRID_LEVEL( I, J, State_Met)

         ! First handle fuelburn-multiplied emissinos
         TMPMULT = (FSC*DTSRCE/32.0d0)

         STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4) + 
     &          (AEIC_EMIT(I,J,L,FB_INDEX)*96.0d0*TMPMULT
     &           *SULFCONV)
         STT(I,J,L,IDTSO2) = STT(I,J,L,IDTSO2) +
     &          ((AEIC_EMIT(I,J,L,FB_INDEX)*64.0d0*TMPMULT)
     &           *(1.d0-SULFCONV))
         STT(I,J,L,IDTNO) = STT(I,J,L,IDTNO) +
     &           (AEIC_EMIT(I,J,L,NOx_INDEX)*DTSRCE)
         STT(I,J,L,IDTCO) = STT(I,J,L,IDTCO) +
     &           (AEIC_EMIT(I,J,L,CO_INDEX)*DTSRCE)
         STT(I,J,L,IDTOCPI) = STT(I,J,L,IDTOCPI) +
     &           (AEIC_EMIT(I,J,L,OC_INDEX)*DTSRCE)
         STT(I,J,L,IDTBCPI) = STT(I,J,L,IDTBCPI) +
     &           (AEIC_EMIT(I,J,L,BC_INDEX)*DTSRCE)
#if defined( UCX )
         IF ( IDTH2O .ne. 0 ) THEN
            STT(I,J,L,IDTH2O) = STT(I,J,L,IDTH2O) +
     &              (1.26d0*AEIC_EMIT(I,J,L,FB_INDEX)*DTSRCE)
         ENDIF
#endif

         ! Only emit HC within chemgrid to avoid buildup
         IF ( L .le. MAXLEV ) THEN
            TMPMULT = THC2TOG*DTSRCE*AEIC_EMIT(I,J,L,HC_INDEX)

            STT(I,J,L,IDTMACR) = STT(I,J,L,IDTMACR) +
     &               (TMPMULT*ACF_MACR)
            STT(I,J,L,IDTRCHO) = STT(I,J,L,IDTRCHO) +
     &               (TMPMULT*ACF_RCHO)
            STT(I,J,L,IDTACET) = STT(I,J,L,IDTACET) +
     &               (TMPMULT*ACF_ACET)
            STT(I,J,L,IDTALD2) = STT(I,J,L,IDTALD2) +
     &               (TMPMULT*ACF_ALD2)
            STT(I,J,L,IDTALK4) = STT(I,J,L,IDTALK4) +
     &               (TMPMULT*ACF_ALK4)
            STT(I,J,L,IDTC2H6) = STT(I,J,L,IDTC2H6) +
     &               (TMPMULT*ACF_C2H6)
            STT(I,J,L,IDTC3H8) = STT(I,J,L,IDTC3H8) +
     &               (TMPMULT*ACF_C3H8)
            STT(I,J,L,IDTCH2O) = STT(I,J,L,IDTCH2O) +
     &               (TMPMULT*ACF_CH2O)
            STT(I,J,L,IDTPRPE) = STT(I,J,L,IDTPRPE) +
     &               (TMPMULT*ACF_PRPE)
         ENDIF

         ! Surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( I, J, L )

         ! ND32 -- NOx aircraft emissions in [molec/cm2]
         IF ( ND32 > 0 .and. L <= LLCHEM ) THEN
            AD32_ac(I,J,L) = AD32_ac(I,J,L) +
     &                       (AEIC_EMIT(I,J,L,NOx_INDEX)  *
     &                       XNUMOL(IDTNO2) / AREA_CM2 )
         ENDIF

         ! ND13 -- SO2 aircraft emissions in [kg S/box/timestep]
         IF ( ND13 > 0 .and. L <= LD13 ) THEN 
            AD13_SO2_ac(I,J,L) = AD13_SO2_ac(I,J,L) +
     &                           ((AEIC_EMIT(I,J,L,FB_INDEX)*FSC*DTSRCE*
     &                           32.0d0/64.0d0)*(1.d0-SULFCONV))
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Debug
      !WRITE(6,*) 'POSTSO4:  ', SUM(STT(:,:,:,IDTSO4))
      !WRITE(6,*) 'POSTNOx:  ', SUM(STT(:,:,:,IDTNO))
      !FBSUM = FBSUM + (DTSRCE*SUM(AEIC_EMIT(:,:,:,FB_INDEX)))
      !WRITE(6,*) 'Fuel burn in current year: ', FBSUM*1.d-9, 'Tg'

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE EMIT_AEIC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_aeic
!
! !DESCRIPTION: Subroutine INIT\_AEIC allocates and initializes module
!  variables. (sde, 12/16/12)
!\\
!\\
! !INTERFACE:

      SUBROUTINE INIT_AEIC
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD
! 
! !REVISION HISTORY:
!  01 Aug 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_AEIC begins here!
      !=================================================================

      ! AEIC_1x1 is used temporarily for input grids
      ALLOCATE( AEIC_1x1( 360, 180, NLVLAEIC, 4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEIC_1x1' )
      AEIC_1x1 = 0d0

      ALLOCATE( AEIC_REGRID( IIPAR, JJPAR, NLVLAEIC, 4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEIC_REGRID' )
      AEIC_REGRID = 0d0
  
      ! AEIC_EMIT holds the output grids
      ALLOCATE( AEIC_EMIT( IIPAR, JJPAR, LLPAR, 6 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEIC_EMIT' )
      AEIC_EMIT = 0d0

      END SUBROUTINE INIT_AEIC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_aeic
!
! !DESCRIPTION: Subroutine CLEANUP\_AEIC deallocates module variables.
!  (sde, 12/16/12)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_AEIC
! 
! !REVISION HISTORY:
!  01 Aug 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_AEIC begins here!
      !=================================================================

      IF ( ALLOCATED( AEIC_1x1    ) ) DEALLOCATE( AEIC_1x1    )
      IF ( ALLOCATED( AEIC_EMIT   ) ) DEALLOCATE( AEIC_EMIT   )

      END SUBROUTINE CLEANUP_AEIC
!EOC
      END MODULE AEIC_MOD
