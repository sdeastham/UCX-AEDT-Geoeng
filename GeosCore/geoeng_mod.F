!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: geoeng_mod
!     
! !DESCRIPTION: Module GEOENG\_MOD contains routines and variables which
!  are associated with SO2 or SO4 geoengineering emissions
!\\   
!\\   
! !INTERFACE: 
!
      MODULE GEOENG_MOD
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
      REAL*8,             PUBLIC :: GE_TGS         ! TgS/year to deliver
      REAL*8,             PUBLIC :: GE_PBASE       ! Geoengineering delivery base (hPa)
      REAL*8,             PUBLIC :: GE_PCEIL       ! Geoengineering delivery ceiling (hPa)
      REAL*8,             PUBLIC :: GE_MINLAT      ! Lower latitude edge (deg)
      REAL*8,             PUBLIC :: GE_MAXLAT      ! Upper latitude edge (deg)

      REAL*8,             PUBLIC :: GE_TCONST      ! Surface temperature offset (K)
      REAL*8,             PUBLIC :: GE_TGMFACT     ! Factor of GeoMIP T data
      REAL*8,             PUBLIC :: GE_LAPSEGMFACT ! Factor of GeoMIP lapse data
      REAL*8,             PUBLIC :: GE_QQCONST     ! Precipitation factor
      REAL*8,             PUBLIC :: GE_QQFACTOR    ! Final precip factor
      REAL*8,             PUBLIC :: GE_QQGMFACT    ! Factor of GeoMIP precip data
      CHARACTER(LEN=255), PUBLIC :: GE_ROOT_NCDF   ! GeoMIP NetCDF folder
      LOGICAL,            PUBLIC :: QQADJUST       ! Are we adjusting precip
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: EMIT_GEOENG
      PUBLIC  :: OFFSET_MET_FIELDS
      PUBLIC  :: INIT_GEOENG
      PUBLIC  :: CLEANUP_GEOENG
!
! !REVISION HISTORY: 
!  11 Nov 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! MODULE PARAMETERS
      !
      ! AVOGADRO           : Avogadro's number (#/mole)
      ! R_UNIV             : Universal gas constant (J/mol K)
      ! GE_NLAT            : Number of latitudes in NetCDF data
      ! GE_NLEV            : Number of P levels in NetCDF data
      !
      !=================================================================

      REAL*8,  PARAMETER            :: AVOGADRO=6.022d23
      REAL*8, PARAMETER             :: R_UNIV = 8.31448d0
      INTEGER,PARAMETER             :: GE_NLEV = 36

!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Arrays
      !
      ! GE_INJECT          : Injection mass (kg/box) 
      ! GE_TOFFSET         : Temperature offset (K)
      ! GE_LMIN            : Minimum geoengineering height (level)
      ! GE_LMAX            : Maximum geoengineering height (level)
      ! GE_LATS            : Latitude edges of 1D data (deg)
      ! GE_LEVS            : Level centers of 1D data (hPa)
      ! GE_REGRID          : Regridding matrix
      ! 
      !=================================================================

      ! Arrays
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE   :: GE_INJECT
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE   :: GE_TOFFSET
      INTEGER,DIMENSION(:,:),ALLOCATABLE    :: GE_LMIN, GE_LMAX
      REAL*8,DIMENSION(:),ALLOCATABLE       :: GE_LATS
      REAL*8,DIMENSION(GE_NLEV)             :: GE_LEVS
      REAL*8,DIMENSION(:,:),ALLOCATABLE     :: GE_REGRID
      LOGICAL                               :: LPRT_GE ! Laziness

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
! !IROUTINE: emit_geoeng
!
! !DESCRIPTION: Subroutine EMIT_GEOENG is used to model sulfate geoengineering
!  using H2SO4.           
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMIT_GEOENG( Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY   : ChmState
      USE TRACERID_MOD, ONLY : IDTSO4
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE GRID_MOD,     ONLY : GET_YEDGE, GET_AREA_M2
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput),INTENT(IN) :: Input_Opt
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  28 Jul 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                   :: KG_PER_STEP
      REAL*8                   :: DTEMIS
      REAL*8                   :: CELL_AREA, CELL_VOL
      REAL*8                   :: CELL_MINLAT, CELL_MAXLAT
      REAL*8                   :: CELL_PCEIL, CELL_PBASE
      REAL*8                   :: CELL_ALTFRAC, CELL_LATFRAC
      REAL*8                   :: GE_GLOBALVOL
      LOGICAL, SAVE            :: FIRST=.TRUE.
      INTEGER, SAVE            :: JMIN, JMAX
      INTEGER                  :: I, J, L
      CHARACTER(LEN=255)       :: OUTMSG

      ! For fields from State_Chm
      REAL*8, POINTER  :: STT(:,:,:,:)

      !=================================================================
      ! EMIT_GEOENG begins here!
      !=================================================================

      ! Initialize fields from State_Chm
      STT      => State_Chm%Tracers

      ! Timestep length in days
      DTEMIS = GET_TS_EMIS() / (60d0 * 24d0)
      IF (FIRST) THEN
         IF (LPRT_GE) THEN
            CALL DEBUG_MSG( '### GEOENG: Initializing geoeng module')
         ENDIF
         ! Zero out injection mass and limits
         GE_INJECT = 0d0
         GE_LMIN = LLPAR
         GE_LMAX = 1

         ! Calculate how much we need to emit per timestep
         ! Assume that we want to emit an average of X TgS per 365.25 days
         KG_PER_STEP = 1.d9 * GE_TGS * DTEMIS / 365.25d0
         ! Convert from TgS to TgSO4
         KG_PER_STEP = KG_PER_STEP * 96d0 / 32d0
         ! Calculate the total volume into which we are injecting
         DO J=1,JJPAR
         CELL_MINLAT = GET_YEDGE(1,J,1)
         CELL_MAXLAT = GET_YEDGE(1,J+1,1)
         IF (CELL_MAXLAT<GE_MINLAT) THEN
            JMIN = MAX(J,JMIN)
            CYCLE
         ENDIF
         IF (CELL_MINLAT>GE_MAXLAT) THEN
            EXIT
         ENDIF
         JMAX = J
         ! Calculate fraction of cell within limits
         CELL_LATFRAC = MIN(CELL_MAXLAT,GE_MAXLAT) - 
     &                       MAX(CELL_MINLAT,GE_MINLAT)
         CELL_LATFRAC = CELL_LATFRAC/(CELL_MAXLAT-CELL_MINLAT)
         CELL_AREA    = GET_AREA_M2(1,J,1)
         DO I=1,IIPAR
         DO L=1,LLPAR
         CELL_PBASE = GET_PEDGE(I,J,L)
         CELL_PCEIL = GET_PEDGE(I,J+1,L)
         IF (CELL_PCEIL.gt.GE_PBASE) THEN
            GE_LMIN(I,J) = L
            CYCLE
         ENDIF
         IF (CELL_PBASE.lt.GE_PCEIL) THEN
            GE_LMAX(I,J) = L
            EXIT
         ENDIF
         ! As for latitude, calculate limits
         ! Are we using the entire cell?
         IF ((CELL_PBASE.gt.GE_PBASE).or.(CELL_PCEIL.lt.GE_PCEIL)) THEN
            CELL_ALTFRAC = MIN(CELL_PBASE,GE_PBASE)/
     &                     MAX(CELL_PCEIL,GE_PCEIL)
            CELL_ALTFRAC = LOG(CELL_ALTFRAC)/LOG(CELL_PBASE/CELL_PCEIL)
         ELSE
            CELL_ALTFRAC = 1.d0
         ENDIF
         CELL_VOL = CELL_ALTFRAC*CELL_LATFRAC*CELL_AREA*
     &          State_Met%BXHEIGHT(I,J,L)
         GE_GLOBALVOL = GE_GLOBALVOL+CELL_VOL
         GE_INJECT(I,J,L) = CELL_VOL
         ENDDO
         ENDDO
         ENDDO
         ! Divide cell volumes by global volume to get fraction
         ! For each target cell, calculate the volume
         DO J=JMIN,JMAX
         DO I=1,IIPAR
         DO L=GE_LMIN(I,J),GE_LMAX(I,J)
         GE_INJECT(I,J,L) = KG_PER_STEP*GE_INJECT(I,J,L)/GE_GLOBALVOL
         ENDDO
         ENDDO
         ENDDO
         FIRST = .FALSE.
         IF (LPRT_GE) THEN
            WRITE(OUTMSG,'(a,a,1x,E10.4)') 'GE_INJECT grid calculated.',
     &          ' Annual mass: ', SUM(GE_INJECT)*365.25/DTEMIS
            CALL DEBUG_MSG(TRIM(OUTMSG))
         ENDIF
      ENDIF

      ! Simply apply the precalculated matrix
      ! SDE DEBUG
      GE_GLOBALVOL = SUM(SUM(SUM(STT(:,:,:,IDTSO4),3),2),1)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
      DO J=JMIN,JMAX
      DO I=1,IIPAR
      DO L=GE_LMIN(I,J),GE_LMAX(I,J)
      STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4) + GE_INJECT(I,J,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF (LPRT_GE) THEN
          WRITE(OUTMSG,'(a,1x,E10.4,1x,a)') 'Geoengineering emission: ',
     & ((1.d0/3.d0)*(SUM(SUM(SUM(STT(:,:,:,IDTSO4),3),2),1)
     &   -GE_GLOBALVOL)), ' TgS as SO4'
          CALL DEBUG_MSG(TRIM(OUTMSG))
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE EMIT_GEOENG
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: offset_met_fields
!
! !DESCRIPTION: Subroutine OFFSET\_MET\_FIELDS applies offsets to met fields
!  to crudely simulate the effects of climate change when appropriate met
!  forecasts are not available
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OFFSET_MET_FIELDS(FIELDTYPE,State_Met,Input_Opt)
!
! !USES:
!
!      USE DAO_MOD,     ONLY : T, TMPU1, TMPU2, TSKIN, TS
      USE PRESSURE_MOD,ONLY : GET_PCENTER, GET_PEDGE
!      USE LOGICAL_MOD, ONLY : LCC_SFCTGEOMIP, LCC_SFCTCONST
!      USE LOGICAL_MOD, ONLY : LCC_QQGEOMIP,   LCC_QQCONST
!      USE LOGICAL_MOD, ONLY : LCC_LAPSEGEOMIP,LCC_LAPSECONST
!      USE LOGICAL_MOD, ONLY : LCC_TDELTA,LCC_QQDELTA,LCC_TLAPSE
      USE TIME_MOD,    ONLY : GET_MONTH
      USE ERROR_MOD,   ONLY : ALLOC_ERR, ERROR_STOP
!      USE DAO_MOD,     ONLY : FROCEAN
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: FIELDTYPE  ! Type of field (A3, A6, I6, A1, I3)
      TYPE(OptInput),INTENT(IN) :: Input_Opt
      TYPE(MetState),INTENT(INOUT) :: State_Met
!
! !REVISION HISTORY: 
!  15 Oct 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER,SAVE                       :: LASTMIN=-1
      INTEGER,SAVE                       :: LASTMONTH=0
      INTEGER                            :: I, J, L, AS, CURRMONTH
      REAL*8,PARAMETER                   :: PMAX = 50.d0
      REAL*8,PARAMETER                   :: PMIN = 300.d0
      LOGICAL,SAVE                       :: FIRST=.TRUE.
      REAL*8, SAVE                       :: INTERPFAC, TMAX, TMIN
      REAL*8                             :: TLOCAL, PLOCAL, QQLOCAL
      REAL*8                             :: QQLATLAND, QQLATOCEAN
      REAL*8                             :: QQOCEAN
      REAL*8                             :: PLAST,PNEXT,PINTERP,TINTERP
      REAL*8                             :: PLASTGM,PNEXTGM
      INTEGER                            :: PLEVGM, GE_NC_MAXLEV
      LOGICAL                            :: EXTRAPGM
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:):: DTARRAY
      CHARACTER(LEN=255)                 :: DBGMSG
      LOGICAL :: LCC_SFCTGEOMIP, LCC_SFCTCONST
      LOGICAL :: LCC_QQGEOMIP,   LCC_QQCONST
      LOGICAL :: LCC_LAPSEGEOMIP,LCC_LAPSECONST
      LOGICAL :: LCC_TDELTA,LCC_QQDELTA,LCC_TLAPSE

      !=================================================================
      ! OFFSET_MET_FIELDS begins here!
      !=================================================================

      ! Get variables from input options
      LCC_SFCTGEOMIP = Input_Opt%LCC_SFCTGEOMIP
      LCC_TDELTA = Input_Opt%LCC_TDELTA
      LCC_TLAPSE = Input_Opt%LCC_TLAPSE
      LCC_SFCTCONST = Input_Opt%LCC_SFCTCONST
      LCC_QQDELTA = Input_Opt%LCC_QQDELTA
      LCC_QQGEOMIP = Input_Opt%LCC_QQGEOMIP
      LCC_LAPSEGEOMIP = Input_Opt%LCC_LAPSEGEOMIP
      LCC_QQCONST = Input_Opt%LCC_QQCONST
      LCC_LAPSECONST = Input_Opt%LCC_LAPSECONST

      ! Do we need to apply any offsets?
      IF (.not.(LCC_TDELTA.or.LCC_QQDELTA)) THEN
          FIRST = .FALSE.
          RETURN
      ENDIF

      ! Do we need to recalculate GE_TOFFSET etc.?
      IF (FIRST) THEN
         ! Allocate temperature offset array
         IF (LCC_TDELTA) THEN
            ALLOCATE( GE_TOFFSET( IIPAR, JJPAR, LLPAR ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'GE_TOFFSET' )
            GE_TOFFSET = 0d0
         ENDIF
         IF (LCC_QQDELTA) THEN
            GE_QQFACTOR = 1d0
         ENDIF
         IF (LCC_SFCTCONST) THEN
            GE_TOFFSET(:,:,1) = GE_TCONST
         ENDIF
         IF (LCC_LAPSECONST) THEN
            DO L=1,LLPAR
               GE_TOFFSET(:,:,L) = GE_TOFFSET(:,:,1)
            ENDDO
         ENDIF
         IF (LCC_QQCONST) THEN
            GE_QQFACTOR = 1.d0 + GE_QQCONST
         ENDIF
      ENDIF
      IF (LCC_QQGEOMIP.or.LCC_SFCTGEOMIP) THEN
         CURRMONTH = GET_MONTH()
         IF (CURRMONTH.ne.LASTMONTH) THEN
            ! Load new monthly data from NetCDF
            LASTMONTH = CURRMONTH
            IF (LCC_SFCTGEOMIP) THEN
               ALLOCATE(DTARRAY(IIPAR,JJPAR,GE_NLEV), STAT=AS)
               IF (AS /= 0 ) CALL ALLOC_ERR( 'DTARRAY' )
               DTARRAY = 0d0
               DTARRAY = GET_GE_AIRTEMP(CURRMONTH,GE_NC_MAXLEV)
               ! The NetCDF data may not extend far enough
               GE_NC_MAXLEV = MIN(GE_NC_MAXLEV,GE_NLEV)
               ! Do surface now
               GE_TOFFSET = 0d0
               GE_TOFFSET(:,:,1) = DTARRAY(:,:,1)*GE_TGMFACT
               IF (LCC_LAPSEGEOMIP) THEN
                  ! Just use current pressure
                  DO I=1,IIPAR
                  DO J=1,JJPAR
                  PNEXT = GET_PEDGE(I,J,1)
                  PLASTGM = -1.d0
                  PNEXTGM = GE_LEVS(1)
                  PLEVGM = 1
                  EXTRAPGM = .FALSE.
                  DO L=2,LLPAR
                     IF (.not.EXTRAPGM) THEN
                        PLAST = PNEXT
                        PNEXT = GET_PCENTER(I,J,L)
                        DO WHILE (PNEXT.lt.PNEXTGM)
                           PLASTGM = PNEXTGM
                           PLEVGM = PLEVGM + 1
                           IF (PLEVGM.gt.GE_NC_MAXLEV) THEN
                              PNEXTGM = 0d0
                              EXTRAPGM = .TRUE.
                           ELSE
                              PNEXTGM = GE_LEVS(PLEVGM)
                           ENDIF
                        ENDDO
                        IF (EXTRAPGM) THEN
                           TINTERP = DTARRAY(I,J,GE_NC_MAXLEV)
                        ELSEIF (PLEVGM.eq.1) THEN
                           TINTERP = DTARRAY(I,J,1)
                        ELSE
                           PINTERP = (PNEXT-PNEXTGM)/(PLASTGM-PNEXTGM)
                           TINTERP = DTARRAY(I,J,PLEVGM-1) + 
     &            (PINTERP*(DTARRAY(I,J,PLEVGM)-DTARRAY(I,J,PLEVGM-1)))
                        ENDIF
                     ENDIF
                     GE_TOFFSET(I,J,L) = TINTERP*GE_TGMFACT
                  ENDDO
                  ENDDO
                  ENDDO
               ELSEIF (LCC_LAPSECONST) THEN
                  DO L=2,LLPAR
                     GE_TOFFSET(:,:,L) = GE_TOFFSET(:,:,1)
                  ENDDO
               ENDIF
               DEALLOCATE(DTARRAY)
            ENDIF
            IF (LCC_QQGEOMIP) THEN
                ! The value returned by GET_GE_PRECIP is the final value
                ! which we expect to multiply the precipitation
                ! Need to convert to a "deviation" before we adjust
                ! to allow for the relative injection mass
                ! e.g. start as 0.95x normal precipitation, convert to
                ! a -0.05 adjustment
                GE_QQFACTOR = GET_GE_PRECIP(CURRMONTH)-1.d0
                IF (LPRT_GE) THEN
                    WRITE(DBGMSG,'(a,a,F6.2,a)') ' ### GEOENG: GeoMIP ',
     &              'base factor: ', 1.d2*GE_QQFACTOR, '%'
                    CALL DEBUG_MSG( TRIM(DBGMSG) )
                ENDIF
                ! Multiply by relative size (e.g. 1/5 of original
                ! mass therefore get -0.01) and add to 1 to get
                ! the final factor (0.99 in this case)
                GE_QQFACTOR = 1.d0 + (GE_QQGMFACT*GE_QQFACTOR)
                IF (LPRT_GE) THEN
                    WRITE(DBGMSG,'(a,a,F6.2)') ' ### GEOENG: ',
     &         'TgS scaling factor: ', GE_QQGMFACT
                    CALL DEBUG_MSG( TRIM(DBGMSG) )
                    WRITE(DBGMSG,'(a,a,F6.2,a)') ' ### GEOENG: ',
     &         'Adjusting precip to ', 1.d2*GE_QQFACTOR, '%'
                    CALL DEBUG_MSG( TRIM(DBGMSG) )
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      ! Can't do precipitation with MERRA or GEOS-FP as they rely on
      ! fields which cannot easily be co-perturbed
#if defined( MERRA ) || defined( GEOS_57 )
      IF (LCC_QQDELTA) THEN
         CALL ERROR_STOP('Cannot adjust MERRA or GEOS-FP precip!',
     &                  'OFFSET_MET_FIELDS(geoeng_mod.F)')
      ENDIF
#endif

      ! Run through available field types
      SELECT CASE ( FIELDTYPE )
         CASE ( 1 )
            IF (LCC_TDELTA) THEN
            ! A3
#if defined( GEOS_4 ) || defined( GEOS_5 ) || defined( GCAP )
               State_Met%TS(:,:) = State_Met%TS(:,:) + GE_TOFFSET(:,:,1)
#if !defined( GCAP )
               State_Met%TSKIN(:,:) = State_Met%TSKIN(:,:) +
     &                              GE_TOFFSET(:,:,1)
#elif defined( MERRA )
               State_Met%T(:,:,:) = State_Met%T(:,:,:) +
     &                              GE_TOFFSET(:,:,:)
#endif
#endif
            ENDIF
         CASE ( 2 )
#if defined( GEOS_4 ) || defined( GEOS_5 ) || defined( GCAP )
            ! A6
            IF (LCC_TDELTA) THEN
               State_Met%T(:,:,:) = State_Met%T(:,:,:)+GE_TOFFSET(:,:,:)
            ENDIF
            !=================================================
            ! Now adjust precip within precip calculation
            ! (SDE 2014-06-21)
            !=================================================
            !IF (LCC_QQDELTA) THEN
            !   DO L=1,LLPAR
            !      MOISTQ(L,:,:) = MOISTQ(L,:,:)*GE_QQFACTOR(:,:)
            !   ENDDO
            !ENDIF
            !=================================================
#endif
         CASE ( 3 )
            ! I6
         CASE ( 4 )
            ! A1
            IF (LCC_TDELTA) THEN
#if defined( MERRA ) || defined( GEOS_57 )
               State_Met%TS(:,:) = State_Met%TS(:,:) + GE_TOFFSET(:,:,1)
               State_Met%TSKIN(:,:) = State_Met%TSKIN(:,:) +
     &                                GE_TOFFSET(:,:,1)
#endif
            ENDIF
         CASE ( 5 )
            ! I3
            IF (LCC_TDELTA) THEN
#if defined( GEOS_57 )
               State_Met%TMPU1(:,:,:) = State_Met%TMPU1(:,:,:) + 
     &                                         GE_TOFFSET(:,:,:)
               State_Met%TMPU2(:,:,:) = State_Met%TMPU2(:,:,:) + 
     &                                             GE_TOFFSET(:,:,:)
#endif
            ENDIF
         CASE DEFAULT
            ! Do nothing
      END SELECT

      FIRST = .FALSE.

      END SUBROUTINE OFFSET_MET_FIELDS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ge_airtemp
!
! !DESCRIPTION: Function GET\_GE\_AIRTEMP loads temperature
!  delta fields from a NetCDF file.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GE_AIRTEMP (TARG_MO_NC,OUTMAX) RESULT (OUTARRAY)
!
! !USES:
!
      USE FILE_MOD,       ONLY : IOERROR
      USE ERROR_MOD,      ONLY : ERROR_STOP, ALLOC_ERR
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
!
! !INPUT PARAMETERS:
!
      REAL*8,DIMENSION(IIPAR,JJPAR,GE_NLEV) :: OUTARRAY
      INTEGER,INTENT(OUT)                   :: OUTMAX
      INTEGER,INTENT(IN)                    :: TARG_MO_NC
!
! !REVISION HISTORY: 
!  19 Apr 2014 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)                :: FILENAME, DBGMSG
      ! For NetCDF
      INTEGER                           :: fId
      INTEGER                           :: st4d(4), ct4d(4)   ! Start + count, for 4D arrays 
      REAL*8,DIMENSION(360,181,GE_NLEV) :: INARRAY

      !=================================================================
      ! GET_GE_AIRTEMP begins here!
      !=================================================================

      ! Read direct from master file
      WRITE(FILENAME,'(a,a)') TRIM(GE_ROOT_NCDF),'/GeoMIPTDelta.nc'

      ! In 4D, lon, lat, pressure level, then month
      ! Data archived on generic 1x1 grid
      st4d = (/1,1,1,TARG_MO_NC/)
      ct4d = (/360,181,GE_NLEV,1/)

      IF (LPRT_GE) THEN
         WRITE(DBGMSG,'(a,a,a)') ' ### GEOENG: Reading ', 
     &        'temperature delta from ', TRIM( FILENAME )
         CALL DEBUG_MSG( TRIM(DBGMSG) )
      ENDIF

      CALL NcOp_Rd (fId,TRIM(FILENAME))
      CALL NcRd( INARRAY, fId, 'dT', st4d, ct4d )
      CALL NcRd( OUTMAX, fId, 'maxLev' )
      CALL NcCl( fId )

      ! Regrid to local grid
      CALL DO_REGRID_1x1(GE_NLEV,'K',INARRAY,OUTARRAY)

      ! Return to calling program
      END FUNCTION GET_GE_AIRTEMP
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ge_precip
!
! !DESCRIPTION: Function GET\_GE\_PRECIP loads precipitation factors
!  from a NetCDF file.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GE_PRECIP (TARG_MO_NC) RESULT (OUTVAL)
!
! !USES:
!
      USE FILE_MOD,  ONLY : IOERROR
      USE ERROR_MOD, ONLY : ERROR_STOP, ALLOC_ERR
!
! !INPUT PARAMETERS:
!
      REAL*8                            :: OUTVAL  
      INTEGER,INTENT(IN)                :: TARG_MO_NC
!
! !REVISION HISTORY: 
!  19 Apr 2014 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)                :: FILENAME, DBGMSG
      ! For NetCDF
      INTEGER                           :: fId
      INTEGER                           :: st1d(1), ct1d(1)   ! Start + count, for 1D arrays 
      REAL*8                            :: INARRAY(1)

      !=================================================================
      ! GET_GE_PRECIP begins here!
      !=================================================================

      ! Read direct from master file
      WRITE(FILENAME,'(a,a)') TRIM(GE_ROOT_NCDF),'/GeoMIPQQDelta.nc'

      ! In 1D - just month
      st1d(1) = TARG_MO_NC
      ct1d(1) = 1

      IF (LPRT_GE) THEN
         WRITE(DBGMSG,'(a,a,a)') ' ### GEOENG: Reading ', 
     &        'precip factors from ', TRIM( FILENAME )
         CALL DEBUG_MSG( TRIM(DBGMSG) )
      ENDIF

      ! Assume that we want the global delta for now
      CALL NcOp_Rd (fId,TRIM(FILENAME))
      CALL NcRd( INARRAY, fId, 'dQQ', st1d, ct1d )
      CALL NcCl( fId )

      OUTVAL = INARRAY(1)

      ! Return to calling program
      END FUNCTION GET_GE_PRECIP
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_geoeng
!
! !DESCRIPTION: Subroutine INIT\_GEOENG initializes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GEOENG( am_I_Root, Input_Opt )
!
! !USES:
!
      USE TIME_MOD,     ONLY : GET_YEAR
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root  ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt
!
! !OUTPUT VARIABLES:
!
!
! !REVISION HISTORY: 
!  11 Nov 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_GEOENG begins here!
      !=================================================================

      WRITE( 6,'(a)') REPEAT( '=', 79 )
      WRITE( 6,'(a)') 'G E O E N G I N E E R I N G'
      WRITE( 6,'(a)') 'Routines written by SEBASTIAN D. EASTHAM'
      WRITE( 6,'(a)') 'Last Modification Date: 2014-07-15'
      WRITE( 6,'(a)') REPEAT( '=', 79 )

      ! Get LPRT
      LPRT_GE = Input_Opt%LPRT

      IF (Input_Opt%LGESULF) THEN
         ! Allocate geoengineering injection array
         ALLOCATE( GE_INJECT( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GE_INJECT' )
         GE_INJECT = 0d0

         ! Allocate arrays of altitude limits
         ALLOCATE( GE_LMIN( IIPAR,JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GE_LMIN' )
         GE_LMIN = LLPAR

         ALLOCATE( GE_LMAX( IIPAR,JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GE_LMAX' )
         GE_LMAX = 1
      ENDIF

      ! Level center data
      GE_LEVS(1:36) = (/1013.25d0,  1000.d0,    925.d0,     850.d0,
     &                  700.d0,     600.d0,     500.d0,     400.d0,
     &                  300.d0,     250.d0,     200.d0,     170.d0,
     &                  150.d0,     130.d0,     115.d0,     100.d0,
     &                  90.d0,      80.d0,      70.d0,      50.d0,
     &                  30.d0,      20.d0,      15.d0,      10.d0,
     &                  7.d0,       5.d0,       3.d0,       2.d0,
     &                  1.5d0,      1.d0 ,      5.d-1,      3.d-1,
     &                  2.d-1,      1.d-1,      5.d-2,      3.d-2  /)

      ! Safety factors
      IF (Input_Opt%LCC_SFCTCONST) GE_TGMFACT=0.d0
      IF (Input_Opt%LCC_SFCTGEOMIP) GE_TCONST=0.d0

      IF (Input_Opt%LCC_LAPSECONST) GE_LAPSEGMFACT=0.d0

      IF (Input_Opt%LCC_QQCONST) GE_QQGMFACT=0.d0
      IF (Input_Opt%LCC_QQGEOMIP) GE_QQCONST=0.d0

      QQADJUST=(.not.(Input_Opt%LCC_QQGEOMIP.or.Input_Opt%LCC_QQCONST))
      IF (.not.QQADJUST) THEN
         GE_QQFACTOR = 1.0d0
      ENDIF

      END SUBROUTINE INIT_GEOENG
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_geoeng
!
! !DESCRIPTION: Subroutine CLEANUP\_GEOENG deallocates module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GEOENG
!
! !REVISION HISTORY: 
!  11 Nov 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_GEOENGAT begins here!
      !=================================================================

      IF ( ALLOCATED( GE_INJECT   ) ) DEALLOCATE( GE_INJECT   )
      IF ( ALLOCATED( GE_TOFFSET  ) ) DEALLOCATE( GE_TOFFSET  )
      IF ( ALLOCATED( GE_LMIN     ) ) DEALLOCATE( GE_LMIN     )
      IF ( ALLOCATED( GE_LMAX     ) ) DEALLOCATE( GE_LMAX     )
      IF ( ALLOCATED( GE_REGRID   ) ) DEALLOCATE( GE_REGRID   )
      
      END SUBROUTINE CLEANUP_GEOENG
!
!EOC
      END MODULE GEOENG_MOD
