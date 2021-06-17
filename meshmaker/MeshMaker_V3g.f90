!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>       T_MeshMaker.f95: Code unit including the main program         >
!>                      and all related routines                       >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      MODULE MeshMaker_Data
!
         SAVE
!
! ----------
! ...... Double precision arrays
! ----------
!
         REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: V,  A,  D1, D2, RC
         REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: DX, DY, DZ, X,  Y,  Z,  Xr,  DRm, DRp
         REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: DX_adj, DY_adj
!
         REAL(KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: DLXYZ
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND = 8) :: X_ref, Y_ref, Z_ref, factor, deg
!
! ----------
! ...... Integer arrays
! ----------
!
         INTEGER, ALLOCATABLE, DIMENSION(:) :: ijk_location
!
         INTEGER(KIND = 2), ALLOCATABLE, DIMENSION(:  ) :: jmin, jmax
!
         INTEGER(KIND = 2), ALLOCATABLE, DIMENSION(:,:) :: RefSurfLoc, kmin, kmax
!
! ----------
! ...... Integer parameters
! ----------
!
         INTEGER, PARAMETER :: MESH_Unit = 12, MINC_Unit = 14, InactElem_Unit = 15
         INTEGER, PARAMETER :: Elem_File_Unit     = 20, Conx_File_Unit     = 21, RefSurfaceCoord_Unit  = 22
         INTEGER, PARAMETER :: TempElem_File_Unit = 25, TempConx_File_Unit = 26, Elem_Coordinates_Unit = 30
         INTEGER, PARAMETER :: Corner_Unit = 100, InactElement_Unit = 101
!
! ----------
! ...... Common integer variables
! ----------
!
         INTEGER :: ElemNameLength
!
         INTEGER :: Max_NumElem
!
         INTEGER :: NXmax, NYmax, NZmax, NLayer
         INTEGER :: LLoop_begin, LLoop_end, max_offset, letter_shift = 0
!
! ----------
! ...... Character variables
! ----------
!
         CHARACTER(LEN = 1)  :: level_of_grid_generation = 'F'
         CHARACTER(LEN = 3)  :: FormatType, LengthUnits, coordinate_order
         CHARACTER(LEN = 8)  :: grid_numbering_system = 'standard'
         CHARACTER(LEN = 11) :: coordinates
!
         CHARACTER(LEN = 8), PARAMETER :: TempElem_file_name = 'TempElem', TempConx_file_name = 'TempConx'
!
         CHARACTER(LEN = :), ALLOCATABLE :: mesh_file_name, elem_file_name, conx_file_name, corners_file_name, IE_corners_file_name
!
! ----------
! ...... Character arrays
! ----------
!
         CHARACTER(LEN = :), ALLOCATABLE, DIMENSION(:) :: name
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: areas_for_HeatExch_Solution = .TRUE., media_by_number = .FALSE.
         LOGICAL :: partial_processing, VTK_output, renumbered_inactive_elements, old_style_inactive_elements
         LOGICAL :: ensure_grid_continuity = .FALSE., smooth_continuous_surface = .FALSE., active_boundaries = .FALSE.
!
! ----------
! ...... Logical arrays
! ----------
!
         LOGICAL, ALLOCATABLE, DIMENSION(:)   :: excluded, fractured
!
!
      END MODULE MeshMaker_Data
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      PROGRAM T_MeshMaker
!
!        Current version developed by G. Moridis in February/March 2015
!
!        New capabilities:
!    (a) Namelist-based data entry
!    (b) Definition of heterogeneous regions/subdomains (regularly/irregularly-shaped)
!    (c) Definition of regularly/irregularly-shaped boundaries
!    (d) Ability to define regularly/irregularly-shaped exclusion and inclusion zones
!
! ...... Modules to be used
!
         USE MeshMaker_Data
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     MeshMaker Main Program: Calls several lower-level routines      *
!*                    which execute computations                       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: AxesOrigin_X = 0.0d0, AxesOrigin_R = 0.0d0, AxesOrigin_Y = 0.0d0, AxesOrigin_Z = 0.0d0
      REAL(KIND = 8) :: inclination_angle = 0.0d0
!
! -------
! ... Integer Variables
! -------
!
      INTEGER :: ierG
!
      INTEGER :: MaxNum_X_Subdivisions  = 0, MaxNum_R_Subdivisions = 0, MaxNum_Y_Subdivisions = 0, MaxNum_Th_Subdivisions = 0, MaxNum_Z_Subdivisions = 0
      INTEGER :: ElemName_NumCharacters = 5
      INTEGER :: beginning_of_partial_LLoop, end_of_partial_LLoop, NLMax
!
! -------
! ... Character Variables
! -------
!
      CHARACTER(LEN = 120) :: TITLE
      CHARACTER(LEN = 11)  :: coordinate_system  = '           '
      CHARACTER(LEN = 1)   :: addition_axis
      CHARACTER(LEN = 3)   :: output_file_format = 'old'
      CHARACTER(LEN = 4)   :: length_units = '  '
      CHARACTER(LEN = 6)   :: header
      CHARACTER(LEN = 14)  :: name_of_elem_file, name_of_conx_file, name_of_mesh_file, name_of_corners_file
!
! -------
! ... Logical Variables
! -------
!
      LOGICAL :: addition_at_MinAxis, addition_at_MaxAxis
!
      LOGICAL :: First_Call = .TRUE., exists = .FALSE., Flag_MINC, special_features
!
! -------
! ... Saving Variables
! -------
!
      SAVE :: First_Call
!
! -------
! ... Namelists
! -------
!
      NAMELIST/ Grid_Specifications / coordinate_system, output_file_format, length_units, grid_numbering_system, ElemName_NumCharacters,                        &
     &                                MaxNum_X_Subdivisions, MaxNum_R_Subdivisions, MaxNum_Y_Subdivisions, MaxNum_Th_Subdivisions, MaxNum_Z_Subdivisions,        &
     &                                AxesOrigin_X, AxesOrigin_R, AxesOrigin_Y, AxesOrigin_Z, coordinate_order, inclination_angle, areas_for_HeatExch_Solution,  &
     &                                media_by_number, name_of_mesh_file, level_of_grid_generation, special_features
!
      NAMELIST/ Special_Grid_Features / renumbered_inactive_elements, old_style_inactive_elements, VTK_output, name_of_corners_file, letter_shift,  &
     &                                  addition_at_MinAxis, addition_at_MaxAxis, addition_axis, ensure_grid_continuity, smooth_continuous_surface, &
     &                                  active_boundaries
!
      NAMELIST/ Partial_Grid_Specifications / beginning_of_partial_LLoop, end_of_partial_LLoop, name_of_elem_file, name_of_conx_file
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of TOUGH+
!
      IF (First_Call) THEN
         WRITE(*,6000)
         First_Call = .FALSE.
      END IF
!
! -------
! ... Reading the heading/title of the input data file
! -------
!
      READ (UNIT = *, FMT = 6002) TITLE
      WRITE(UNIT = *, FMT = 6003) TITLE
!
! ... Initialization
!
!
      coordinates = '           '
      FormatType  = 'old'
      LengthUnits = '  '
      ElemNameLength = 5
      addition_axis  = ' '
!
      level_of_grid_generation = 'F'
!
      coordinate_order     = '   '
      name_of_mesh_file    = '            '
      name_of_elem_file    = '            '
      name_of_conx_file    = '            '
      name_of_corners_file = '            '
!
      LLoop_begin  = 0
      LLoop_end    = 0
      letter_shift = 0
!
      beginning_of_partial_LLoop = 0
      end_of_partial_LLoop       = 0
!
      special_features   = .FALSE.
      partial_processing = .FALSE.
      VTK_output         = .FALSE.
      renumbered_inactive_elements = .TRUE.
      old_style_inactive_elements  = .FALSE.
!
      ensure_grid_continuity    = .FALSE.
      smooth_continuous_surface = .FALSE.
!
      addition_at_MinAxis = .FALSE.
      addition_at_MaxAxis = .FALSE.
!
      NXmax = 0
      NYmax = 0
      NZmax = 0
!
      X_ref = 0.0d0
      Y_ref = 0.0d0
      Z_ref = 0.0d0
      deg   = 0.0d0
!
! -------
! ... Read the data describing the grid characteristics
! -------
!
 1000 READ (UNIT = *, FMT = '(A6)', IOSTAT = ierG ) header
!
      IF(ierG /= 0) THEN
         WRITE(UNIT = *, FMT = 6050)
         STOP
      END IF
!
      IF( header(1:3) == '   ' ) THEN
         GO TO 1000
      ELSE IF( header == '>>>GEN' ) THEN
         CONTINUE
      ELSE
         WRITE(UNIT = *, FMT = 6080) header(1:6)
         STOP
      END IF
!
      READ (UNIT = *, NML = Grid_Specifications, IOSTAT = ierG )
!
! ... Stop if there is a problem reading the namelist
!
      IF(ierG /= 0) THEN
         WRITE(UNIT = *, FMT = 6100) 'Grid_Specifications'
         STOP
      END IF
!
! -------
! ... Read the data describing special features of the grid
! -------
!
      IF( special_features ) READ (UNIT = *, NML = Special_Grid_Features, IOSTAT = ierG )
!
! ... Stop if there is a problem reading the namelist
!
      IF(ierG /= 0) THEN
         WRITE(UNIT = *, FMT = 6100) 'Special_Grid_Features'
         STOP
      END IF
!
! -------
! ... Read the data describing specifications for partial grid creation
! -------
!
      IF( (level_of_grid_generation == 'P') .OR. (level_of_grid_generation == 'p') ) THEN
         READ (UNIT = *, NML = Partial_Grid_Specifications, IOSTAT = ierG )
      END IF
!
! ... Stop if there is a problem reading the namelist
!
      IF(ierG /= 0) THEN
         WRITE(UNIT = *, FMT = 6100) 'Partial_Grid_Specifications'
         STOP
      END IF
!
      SELECT CASE (level_of_grid_generation)
      CASE('P','p')
         CONTINUE
      CASE DEFAULT
         level_of_grid_generation = 'F'
      END SELECT
!
! >>> Check the coordinate system
!
      SELECT CASE (coordinate_system(1:3))
      CASE('CAR', 'car', 'Car')
!
         coordinates = 'cartesian'
         NXmax = MaxNum_X_Subdivisions
         NYmax = MaxNum_Y_Subdivisions
         NZmax = MaxNum_Z_Subdivisions
!
         X_ref = AxesOrigin_X
         Y_ref = AxesOrigin_Y
         Z_ref = AxesOrigin_Z
!
         deg = inclination_angle
!
      CASE('TRA', 'Tra', 'tra')
!
         coordinates = 'transformed'
         NXmax = MaxNum_R_Subdivisions
         NYmax = MaxNum_Th_Subdivisions
         NZmax = MaxNum_Z_Subdivisions
!
         X_ref = AxesOrigin_R
         Y_ref = 0.0d0
         Z_ref = AxesOrigin_Z
!
         deg = inclination_angle
!
      CASE('CYL', 'cyl', 'Cyl')
!
         level_of_grid_generation = 'F'
         coordinates = 'cylindrical'
         NXmax = MaxNum_R_Subdivisions
         NYmax = 1
         NZmax = MaxNum_Z_Subdivisions
!
         Z_ref = AxesOrigin_Z
!
      CASE DEFAULT
         WRITE(UNIT = *, FMT = 6101) coordinate_system
         STOP
      END SELECT
!
! >>>
! >>> Open the output files - ELEMENTS and CONNECTIONS next
! >>>
!
      IF( (level_of_grid_generation == 'P') .OR. (level_of_grid_generation == 'p') .OR. &
     &    (ensure_grid_continuity .AND. (level_of_grid_generation == 'F' .OR. level_of_grid_generation == 'f') ) ) &
     &THEN
!
         IF_ElemName: IF( name_of_elem_file == '            ' ) THEN
            ALLOCATE( CHARACTER(LEN=8) :: elem_file_name )
            elem_file_name = 'ELEMENTS'
         ELSE
            ierG = LEN( TRIM(ADJUSTL(name_of_elem_file)) )
            ALLOCATE( CHARACTER(LEN=ierG) :: elem_file_name )
            elem_file_name = TRIM(ADJUSTL(name_of_elem_file))
         END IF IF_ElemName
!
         IF_ConxName: IF( name_of_conx_file == '            ' ) THEN
            ALLOCATE( CHARACTER(LEN=11) :: conx_file_name )
            conx_file_name = 'CONNECTIONS'
         ELSE
            ierG = LEN( TRIM(ADJUSTL(name_of_conx_file)) )
            ALLOCATE( CHARACTER(LEN=ierG) :: conx_file_name )
            conx_file_name = TRIM(ADJUSTL(name_of_conx_file))
         END IF IF_ConxName
!
! ...... Inquire whether the element file exists; open accordingly
!
         INQUIRE( FILE = elem_file_name, EXIST = exists )
         IF_Elem: IF(exists) THEN
                     WRITE(UNIT = *, FMT = 6004) elem_file_name
                     OPEN( UNIT = Elem_File_Unit, FILE = elem_file_name, STATUS = 'OLD')
                  ELSE
                     WRITE(UNIT = *, FMT = 6005) elem_file_name
                     OPEN( UNIT = Elem_File_Unit, FILE = elem_file_name, STATUS = 'NEW')
                  END IF IF_Elem
         REWIND (UNIT = Elem_File_Unit)
!
         INQUIRE( FILE = conx_file_name, EXIST = exists )
         IF_Conx: IF(exists) THEN
                     WRITE(UNIT = *, FMT = 6004) conx_file_name
                     OPEN( UNIT = Conx_File_Unit, FILE = conx_file_name, STATUS = 'OLD')
                  ELSE
                     WRITE(UNIT = *, FMT = 6005) conx_file_name
                     OPEN( UNIT = Conx_File_Unit, FILE = conx_file_name, STATUS = 'NEW')
                  END IF IF_Conx
         REWIND (UNIT = Conx_File_Unit)
!
      END IF

!
      IF( (level_of_grid_generation == 'P' .OR. level_of_grid_generation == 'p') .AND. ensure_grid_continuity ) THEN
!
         OPEN( UNIT = TempElem_File_Unit, FILE = TempElem_file_name )
         OPEN( UNIT = TempConx_File_Unit, FILE = TempConx_file_name )
!
         REWIND (UNIT = TempElem_File_Unit)
         REWIND (UNIT = TempConx_File_Unit)
!
      END IF
!
!
     IF( level_of_grid_generation == 'F' .OR. level_of_grid_generation == 'f' ) THEN
!
         IF_MeshName: IF( name_of_mesh_file == '            ' ) THEN
            ALLOCATE( CHARACTER(LEN=4) :: mesh_file_name )
            mesh_file_name = 'MESH'
         ELSE
            ierG = LEN( TRIM(ADJUSTL(name_of_mesh_file)) )
            ALLOCATE( CHARACTER(LEN=ierG) :: mesh_file_name )
            mesh_file_name = TRIM(ADJUSTL(name_of_mesh_file))
         END IF IF_MeshName
!
! ...... Inquire whether the mesh file exists; open accordingly
!
         INQUIRE( FILE = mesh_file_name, EXIST = exists )
         IF_MESH: IF(exists) THEN
                     WRITE(UNIT = *, FMT = 6004) mesh_file_name
                     OPEN( UNIT = MESH_Unit, FILE = mesh_file_name, STATUS = 'OLD')
                  ELSE
                     WRITE(UNIT = *, FMT = 6005) mesh_file_name
                     OPEN( UNIT = MESH_Unit, FILE = mesh_file_name, STATUS = 'NEW')
                  END IF IF_MESH
         REWIND (UNIT = MESH_Unit)
!
         elem_file_name = '            '
         conx_file_name = '            '
!
     END IF
!
      IF( coordinate_system(1:2) == 'CY' .OR. coordinate_system(1:2) == 'Cy' .OR. coordinate_system(1:2) == 'cy' ) GO TO 1500
!
! ... Old style inactive elements
!
 1500 IF( old_style_inactive_elements )  THEN
         renumbered_inactive_elements = .TRUE.
         OPEN(  UNIT = InactElem_Unit, FILE = 'INACT_ELEM' )
         WRITE( UNIT = InactElem_Unit, FMT  = '(A)' ) 'ina00                                                                                 '
      END IF
!
      IF( active_boundaries )  THEN
         renumbered_inactive_elements = .FALSE.
         old_style_inactive_elements  = .FALSE.
      END IF
!
! >>>
!
      NLMax = MAX( NXmax, NYMax, NZMax )
!
      IF( ( beginning_of_partial_LLoop == 1 .AND. end_of_partial_LLoop <  NLMax ) .OR. &
     &    ( beginning_of_partial_LLoop >  1 .AND. end_of_partial_LLoop >= NLMax ) .OR. &
     &    ( beginning_of_partial_LLoop >  1 .AND. end_of_partial_LLoop <  NLMax ) )    &
     &THEN
!
         IF( coordinate_system(1:3) /= 'CAR' .AND. coordinate_system(1:3) /= 'Car' .AND. coordinate_system(1:3) /= 'car' .AND. &
             coordinate_system(1:3) /= 'TRA' .AND. coordinate_system(1:3) /= 'Tra' .AND. coordinate_system(1:3) /= 'tra' )     &
     &   THEN
            WRITE(*, 6901)
            STOP
         END IF
!
         IF( level_of_grid_generation == 'F' .OR. level_of_grid_generation == 'f' ) THEN
            WRITE(*, 6902)
            STOP
         END IF
!
         LLoop_begin = beginning_of_partial_LLoop
         LLoop_end   = end_of_partial_LLoop
         partial_processing = .TRUE.
!
      ELSE IF( beginning_of_partial_LLoop == 0 .AND. end_of_partial_LLoop == 0 ) THEN
         partial_processing = .FALSE.
      ELSE IF( beginning_of_partial_LLoop == 1 .AND. end_of_partial_LLoop >= NLMax ) THEN
         partial_processing = .FALSE.
      ELSE
         WRITE(*, 6903) beginning_of_partial_LLoop, end_of_partial_LLoop
         STOP
      END IF
!
! >>>
!
      SELECT CASE (length_units)
      CASE('m', 'M')
         factor = 1.0d0
      CASE('ft', 'FT', 'Ft')
         factor = 3.048d-1
      CASE('in', 'IN', 'In')
         factor = 2.54d-2
      CASE('km', 'KM', 'Km')
         factor = 1.0d3
      CASE('cm', 'CM', 'Cm')
         factor = 1.0d-2
      CASE('mi', 'MI', 'Mi')
         factor = 1.609344d3
      CASE DEFAULT
         WRITE(UNIT = *, FMT = 6102) 'length', LengthUnits
         STOP
      END SELECT
!
! >>> Check the coordinate order
!
      SELECT CASE (coordinate_order(1:3))
      CASE( 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX', 'xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx' )
         CONTINUE
      CASE DEFAULT
         coordinate_order = '   '
      END SELECT
!
! >>>
!
      X_ref = X_ref * factor
      Y_ref = Y_ref * factor
      Z_ref = Z_ref * factor
!
      LengthUnits = length_units
!
      SELECT CASE (output_file_format)
      CASE('OLD', 'Old', 'old')
         FormatType  = output_file_format
      CASE('NEW', 'New', 'new')
         FormatType  = output_file_format
      CASE DEFAULT
         WRITE(UNIT = *, FMT = 6103) output_file_format
         STOP
      END SELECT
!
! >>>
!
      SELECT CASE (grid_numbering_system(1:2))
      CASE('ST', 'St', 'st')
         IF( coordinates(1:2) == 'cy' ) THEN
            IF(NXMax > NZMax) THEN
               grid_numbering_system = 'column'
            ELSE
               grid_numbering_system = 'layer'
            END IF
         END IF
      CASE('LA', 'La', 'la')
         IF( coordinates(1:2) == 'ca' .OR. coordinates(1:2) == 'tr' ) grid_numbering_system = 'standard'
      CASE('CO', 'Co', 'co')
         IF( coordinates(1:2) == 'ca' .OR. coordinates(1:2) == 'tr' ) grid_numbering_system = 'standard'
      CASE DEFAULT
         grid_numbering_system = 'standard'
         IF(NXMax > NZMax) THEN
            grid_numbering_system = 'column'
         ELSE
            grid_numbering_system = 'layer'
         END IF
      END SELECT
!
! >>>
!
      IF( ElemName_NumCharacters /= 5 .AND. ElemName_NumCharacters /= 8 ) THEN
         WRITE(UNIT = *, FMT = 6105) ElemName_NumCharacters
         STOP
      ELSE
         ElemNameLength = ElemName_NumCharacters
      END IF
!
      Max_NumElem = NXmax * NYmax * NZmax
!
!     IF( Max_NumElem <= 200000 .OR. coordinate_system(1:3) == 'CYL' .OR. coordinate_system(1:3) == 'Cyl' .OR. coordinate_system(1:3) == 'cyl' ) THEN
         ALLOCATE(CHARACTER(LEN=ElemNameLength) :: name(Max_NumElem), STAT = ierG)
!     END IF
!
! >>>
! >>> Inquire whether files CORNERS and IE_CORNERS exists; open accordingly
! >>>
!
      IF_VTK: IF ( VTK_output ) THEN
!
         IF_VTKName: IF( name_of_corners_file == '            ' ) THEN
            ALLOCATE( CHARACTER(LEN=7) :: corners_file_name )
            corners_file_name = 'CORNERS'
         ELSE
            ierG = LEN( TRIM(ADJUSTL(name_of_corners_file)) )
            ALLOCATE( CHARACTER(LEN=ierG) :: corners_file_name )
            corners_file_name = TRIM(ADJUSTL(name_of_corners_file))
         END IF IF_VTKName
!
         INQUIRE( FILE = corners_file_name, EXIST = exists )
!
         IF_Corner: IF(exists) THEN
                       WRITE(UNIT = *, FMT = 6004) corners_file_name
                       OPEN( UNIT = Corner_Unit, FILE = corners_file_name, STATUS = 'OLD')
                    ELSE
                       WRITE(UNIT = *, FMT = 6005) corners_file_name
                       OPEN( UNIT = Corner_Unit, FILE = corners_file_name, STATUS = 'NEW')
                    END IF IF_Corner
!
         REWIND (UNIT = Corner_Unit)
!
         IF( renumbered_inactive_elements ) THEN
!
            IF_IEName: IF( name_of_corners_file == '            ' ) THEN
               ALLOCATE( CHARACTER(LEN=10) :: IE_corners_file_name )
               IE_corners_file_name = 'IE_CORNERS'
            ELSE
               ierG = 3 + LEN( TRIM(ADJUSTL(name_of_corners_file)) )
               ALLOCATE( CHARACTER(LEN=11) :: IE_corners_file_name )
               IE_corners_file_name = 'IE_'//TRIM(ADJUSTL(name_of_corners_file))
            END IF IF_IEName
!
            INQUIRE( FILE = IE_corners_file_name, EXIST = exists )
!
            IF_IECorner: IF(exists) THEN
                          WRITE(UNIT = *, FMT = 6004) IE_corners_file_name
                          OPEN( UNIT = InactElement_Unit, FILE = IE_corners_file_name, STATUS = 'OLD')
                       ELSE
                          WRITE(UNIT = *, FMT = 6005) 'IE_CORNERS'
                          OPEN( UNIT = InactElement_Unit, FILE = IE_corners_file_name, STATUS = 'NEW')
                       END IF IF_IECorner
!
            REWIND (UNIT = InactElement_Unit)
!
         END IF
!
      END IF IF_VTK
!
! >>>>>>
! >>>>>> For thin areal grids
! >>>>>>
!
      IF(ensure_grid_continuity) THEN
!
! ...... Read inputs and define the Reference Surface
!
         CALL Define_Reference_Surface
!
      END IF
!
! -------
! ... End of inputs
! -------
!
      READ( UNIT = *, FMT = '(A3)', IOSTAT = ierG ) header
!
      IF( header(1:3) /= '<<<' ) THEN
         WRITE( UNIT = *, FMT = 6515 ) header
         STOP
      END IF
!
      WRITE( UNIT = *, FMT = 6006 )
!
!
!***********************************************************************
!*                                                                     *
!*                  CALL THE MESH GENERATION PACKAGE                   *
!*                                                                     *
!***********************************************************************
!
!
      CALL Mesh_Maker(Flag_MINC)
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
!6000 FORMAT('MeshMaker         1.0   10 October   2006',//)
!6000 FORMAT('MeshMaker         1.5   15 August    2014',//)
 6000 FORMAT(  'MeshMaker 2.0 .......................... 24 December  2014',//)
!
 6002 FORMAT(A120)
 6003 FORMAT(/1X,131('='),//,5X,'PROBLEM TITLE:  ',A120,//)


!
 6004 FORMAT(' FILE <',A,'> EXISTS ==> OPEN AS AN OLD FILE')
 6005 FORMAT(' FILE <',A,'> DOES NOT EXIST ==> OPEN AS A NEW FILE')
!
 6006 FORMAT(//,' END OF MeshMaker INPUT JOB',/)
!
 6050 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Main Program <T_MeshMaker>: There is a problem reading the heading ">>>GENERAL" of the 1st data block ',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6080 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Main Program <T_MeshMaker>: The heading of the 1st data block "',A,'" different from the expected ">>>GEN" ',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6100 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Main Program <T_MeshMaker>: There is a problem reading the namelist <',A,'> ',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6101 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system of the grid is ',A,': Unknown/Unavailable option'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T10,'>>>  The units of length <',A,'> = "',A,'" are not an available option',//,  &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',      &
     &       //,20('ERROR-'))
!
 6103 FORMAT(//,20('ERROR-'),//,   &
     &       T10,'>>>  The format of the output file <> = "',A,'" are not an available option',//,  &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',      &
     &       //,20('ERROR-'))
!
 6105 FORMAT(//,20('ERROR-'),//,   &
     &       T10,'>>>  The number of characters in the element name <ElemName_NumCharacters> = ',i1,' is not an available option',//,  &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',      &
     &       //,20('ERROR-'))
!
 6515 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Main Program "T_MeshMaker": The ending descriptor "',A3,'" of dataset <GENERAL> is not the required "<<<"'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6901 FORMAT(//,20('ERROR-'),//,   &
     &       T10,'>>>  It is possible to have partial grid processing only when the coordinate system is Cartesian !!!',//,  &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',      &
     &       //,20('ERROR-'))
!
 6902 FORMAT(//,20('ERROR-'),//,   &
     &       T10,'>>>  It is not possible to have partial grid processing when the variable <level_of_grid_generation> = "F" or "f" !!!',//,  &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',      &
     &       //,20('ERROR-'))
!
 6903 FORMAT(//,20('ERROR-'),//,   &
     &       T10,'>>>  It beginning (=',i5.5,') and the end (=',i5.5,') of the loop for partial grid processing are invalid !!!',//,  &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',      &
     &       //,20('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of T_MeshMaker
!
!
      STOP
!
!
!
      END PROGRAM T_MeshMaker
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************




!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
!
      MODULE Grid_Generation_Parameters
!
         SAVE
!
! ----------
! ...... Double precision parameters
! ----------
!
         REAL(KIND = 8), PARAMETER :: pi  = 3.14159265358979324D0
!
! ----------
! ...... Character parameter arrays
! ----------
!
         CHARACTER(LEN = 1), PARAMETER, DIMENSION(36) :: NA = (/ '0','1','2','3','4','5','6','7','8','9',   &
     &                                                           'A','B','C','D','E','F','G','H','I','J',   &
     &                                                           'K','L','M','N','O','P','Q','R','S','T',   &
     &                                                           'U','V','W','X','Y','Z'  /)
!
         CHARACTER(LEN = 5), PARAMETER, DIMENSION(6) :: TYPOS = (/ 'ONE-D','TWO-D','THRED','STANA','STANB','STANT' /)
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER :: NSKIP, NELEMT
!
! ----------
! ...... CHARACTER variables
! ----------
!
         CHARACTER(LEN = 5) :: dominant_medium = '    1'
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: Flag_HetRegions, Flag_Boundaries, Flag_ExclZones, Flag_InclZones
!
! >>>>>
! >>>>>
! >>>>>
! >>>>>
! >>>>>
!
         CONTAINS
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Logarithmic_Discretization( axis, D_max, D_last, D_Log, Nb, IncrFact, DL_number, Dd )
!
!
!
         IMPLICIT NONE
!
! ----------
! ...... Double precision arrays
! ----------
!
         REAL(KIND = 8), INTENT(OUT), DIMENSION(Nb+1:) :: Dd
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND =  8), INTENT(IN) :: D_max, D_last, D_Log
!
         REAL(KIND = 16) :: am, rp, rp1, xm, xn, test, F0x, F1x, F2x, er, IncrFact
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER, INTENT(IN) :: Nb, DL_number
!
         INTEGER :: i, k, np1
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: First_call = .TRUE.
!
! ----------
! ...... Character variables
! ----------
!
         CHARACTER(LEN = 1), INTENT(IN) :: axis
!
! ----------
! ...... Saving variables
! ----------
!
         SAVE First_call
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Logarithmic_Discretization>
!
         IF(First_call) THEN
            WRITE( UNIT = *, FMT = 6000 )
            First_call   = .FALSE.
         END IF
!
! >>>>>>
! >>>>>> Determining basic parameters
! >>>>>>
!
         am  = (D_max - D_last) / D_Log
         rp  = DBLE(DL_number - 1)
         np1 = DL_number
         rp1 = DBLE(np1)
!
         IF( am == rp ) THEN
            xm = 1.1e0_16
            GO TO 1000
         ELSE
            test = am / rp
            IF( test > 1 ) THEN
               xm = MIN( test, 1.5e0_16 )
            ELSE
               xm = MAX( test, 0.8e0_16 )
            END IF
         END IF
!
! >>>>>>
! >>>>>> Find multiplying factor by the modified Newton method
! >>>>>>
!
         DO_Fact: DO k = 1,10000
!
            F0x = xm**np1 - am * xm + am - 1.0e0_16
            F1x = rp1 * (xm**DL_number) - am
            F2x = rp * rp1 * (xm**(DL_number-1))
!
            xn = xm - 2.0e0_16 * F0x * F1x / (2.0e0_16 * F1x * F1x - F0x * F2x )
            WRITE(*,*) 'xm, xn = ',xm,xn
            er = ABS( 1.0d2 * (xm - xn) / xm )
            xm = xn
!
            IF_check: IF( er <= 1.0d-13 ) THEN
               EXIT
            ELSE
               IF( k == 10000 ) THEN
                  WRITE(UNIT = *, FMT = 6010) er
                  STOP
               END IF
            END IF IF_check
!
         END DO DO_Fact
!
  1000   IncrFact = xn
!
         WRITE( UNIT = *, FMT = 6020) k, 'D'//axis, IncrFact
!
         Dd(Nb+1) = D_Log
         DO i=2,DL_number
            Dd(Nb+i) = Dd(Nb+i-1) * IncrFact
         END DO
!
         RETURN
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Logarithmic_Discretization 1.0 ......... 18 January   2018',6X,'Discretize a given length by logarithmically-distributed subdivisions')
!
 6010 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Logarithmic_Discretization>: The increment factor has not converged after 10,000 iterations. The final error is',ES11.4,'%',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6020 FORMAT(/,T2,' After ',I7.7,' iterations, the converged value of the ',A,'-associated increment factor <IncrFact> =',ES20.14E1)
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Logarithmic_Discretization>
!
         RETURN
!
      END SUBROUTINE Logarithmic_Discretization
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      END MODULE Grid_Generation_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Het_Region_Variables
!
!
         SAVE
!
! ----------
! ...... Derived-type variables
! ----------
!
         TYPE Het_Region
!
! ......... For regularly-shaped regions
!
            CHARACTER(LEN = 5)  :: name                       ! Name of the heterogeneous region
!
            CHARACTER(LEN = 11) :: shape                      ! Shape of the heterogeneous region
!
            CHARACTER(LEN = 2 ) :: units                      ! Length units of the heterogeneous region
!
            REAL(KIND = 8), DIMENSION(3) :: LMin, LMax        ! Min and max extent of a cartesian region in the principal directions
!
            REAL(KIND = 8), DIMENSION(3) :: CylBase1Coord, CylBase2Coord     ! Coordinates (X,Y,Z or R,Z) of the centers of the two bases of a cylindrical region
!
            REAL(KIND = 8) :: CylRmin,       CylRmax           ! Minimum and maximum radius of a cylindrical region
            REAL(KIND = 8) :: CylBase1R,     CylBase2R         ! Radii of the top and bottom bases
!
            REAL(KIND = 8), DIMENSION(3) :: EllBase1Coord, EllBase2Coord
!
            REAL(KIND = 8) :: LAxisAngle, LAxis, SAxis
            REAL(KIND = 8) :: Base1LAxis, Base2LAxis, Base1SAxis, Base2SAxis
!
            REAL(KIND = 8), DIMENSION(3) :: SphereCenterCoord                ! Coordinates (X,Y,Z or R,Z) of the center of a spherical region
            REAL(KIND = 8)               :: SphRmin, SphRmax                 ! Minimum and maximum radius of a spherical region
!
            REAL(KIND = 8) :: emissivity
!
! ......... For irregularly-shaped regions
!
            CHARACTER(LEN = 1)  :: DependentVar                     ! Identifiers of the surface equations
!
            CHARACTER(LEN = 1)  :: VertOrStrat1, VertOrStrat2, RefSurf1, RefSurf2
!
            CHARACTER(LEN = 2)  :: EllipsePlane
!
            CHARACTER(LEN = 3)  :: id                               ! Identifier of a boundary or marked region
!
            CHARACTER(LEN = 10) :: TypeEqu1, TypeEqu2               ! Identifiers of the surface equations
!
            INTEGER :: OrderEquSurf1, OrderEquSurf2                 ! Order of the shape equations
!
            INTEGER :: IntTableNum1, IntTableNum2                   ! Numbers of the tables from which interpolations are to be computed
!
            REAL(KIND = 8) :: expon1, expon2, sign1, sign2
!
            REAL(KIND = 8) :: Thick0, Thick1, Thick2
!
            REAL(KIND = 8), POINTER, DIMENSION(:) :: Equ1CoeffA, Equ1CoeffB, Equ1CoeffC, Equ1CoeffD
!
            REAL(KIND = 8), POINTER, DIMENSION(:) :: Equ2CoeffA, Equ2CoeffB, Equ2CoeffC, Equ2CoeffD
!
            REAL(KIND = 8), DIMENSION(3) :: L1Shift, L2Shift           ! Shifts in the computation of the bounding surface equations
!
         END TYPE Het_Region
!
! ----------
! ...... Derived-type variables: Interpolation table definition
! ----------
!
         TYPE Tabular
!
! ......... Double precision arrays
!
            REAL(KIND = 8), POINTER, DIMENSION(:) :: dv, iv1, iv2
!
! ......... Double precision variable
!
            REAL(KIND = 8) :: SearchRadius
!
! ......... Character parameter arrays
!
            CHARACTER(LEN = 8) :: FileName
!
! ......... Integer variables
!
            INTEGER :: NumDataPoints, UnitNum, NumSearchPoints
!
         END TYPE Tabular
!
! ----------
! ...... Derived-type arrays
! ----------
!
         TYPE(Het_Region), ALLOCATABLE, DIMENSION(:) :: Region, bound, ExclZone, InclZone
!
         TYPE(Tabular), ALLOCATABLE, DIMENSION(:) :: IntTable, temp_IntTable
!
! ----------
! ...... Derived-type variables
! ----------
!
         TYPE(Het_Region) :: RefSurface
!
! ----------
! ...... Real variables
! ----------
!
         REAL(KIND = 8) :: dominant_medium_emissivity = 0.0d0
!
! ----------
! ...... Integer arrays
! ----------
!
         INTEGER, ALLOCATABLE, DIMENSION(:) :: MediaSequenceNumber
         INTEGER, ALLOCATABLE, DIMENSION(:) :: InclZone_SequenceNumber
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER :: Num_HetRegions, TotNum_HetRegions, Num_ExclZones, TotNum_ExclZones, Num_InclZones, TotNum_InclZones, Num_Boundaries
!
         INTEGER :: Num_IntTables = 0, Num_Emitting_Units
!
         INTEGER, ALLOCATABLE, DIMENSION(:) :: IntTable_UnitNum, temp_IntDataSet_num, RegionTable_UnitNum
!
!
!
      END MODULE Het_Region_Variables
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Subdomain_Bounding_Surfaces
!
!
!
         USE MeshMaker_Data, ONLY: coordinates
!
         USE Het_Region_Variables
!
! >>>>>
! >>>>>
! >>>>>
! >>>>>
! >>>>>
!
         CONTAINS
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Surface_Coordinate_Limits( subdomain, table_1, table_2, XYZ_coordinates, RZ_coordinates,  &
     &                                         use_InterpData, completed_interpolation, available_InterpData, &
     &                                         XYZ_limits, RZ_limits, find_RefSurface, success )
!
            USE Utility_Functions
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*             Routine determining the media region where              *
!*                     a particular cell belongs                       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... User-defined variable
! -------------
!
            TYPE(Het_Region) :: subdomain
!
            TYPE(Tabular), OPTIONAL, INTENT(IN) :: table_1, table_2
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), OPTIONAL, INTENT(IN), DIMENSION(3) :: XYZ_coordinates
            REAL(KIND = 8), OPTIONAL, INTENT(IN), DIMENSION(2) ::  RZ_coordinates
!
            REAL(KIND = 8), OPTIONAL, INTENT(OUT), DIMENSION(3,2) :: XYZ_limits
            REAL(KIND = 8), OPTIONAL, INTENT(OUT), DIMENSION(2,2) ::  RZ_limits
!
            REAL(KIND = 8), OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables,3) :: available_InterpData
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: value1A, value2A, value1B, value2B, arg1, arg2, dZdX, dZdY
!
            REAL(KIND = 8) :: function_1p, function_2p, function_v, conv_factor
!
            REAL(KIND = 8), DIMENSION(2) :: val
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: i, m, dv = 0, iv_mod, dv_mod, m1, m2, NumData
!
            INTEGER, DIMENSION(2) :: iv = 0
!
! -------------
! ......... Logical variables
! -------------
!
            LOGICAL, OPTIONAL, INTENT(IN)  :: use_InterpData, find_RefSurface
            LOGICAL, INTENT(OUT) :: success
!
            LOGICAL, OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables) :: completed_interpolation
!
            LOGICAL :: First_call = .TRUE., valid, do_not_do
!
! -------------
! ......... Saving variables
! -------------
!
            SAVE First_call
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of <Surface_Coordinate_Limits>
!
!
            IF(First_call) THEN
               WRITE( UNIT = *, FMT = 6000 )
               First_call = .FALSE.
            END IF
!
! ......... Initialization and assignment
!
            success   = .FALSE.
            do_not_do = .FALSE.
!
            function_1p = 0.0d0
            function_2p = 0.0d0
!
!***********************************************************************
!*                                                                     *
!*          Determine the dependent and independent variables          *
!*                                                                     *
!***********************************************************************
!
            IF_Coord: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                    coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
! >>>>>>>>>>>>
!
               IF_DepVarCart: IF( subdomain%DependentVar == 'Z' .OR. subdomain%DependentVar == 'z' ) THEN
!
                  dv     = 3
                  iv(1)  = 1
                  iv(2)  = 2
                  val(1) = XYZ_coordinates(1)
                  val(2) = XYZ_coordinates(2)
!
               ELSE IF( subdomain%DependentVar == 'Y' .OR. subdomain%DependentVar == 'y' ) THEN
!
                  dv     = 2
                  iv(1)  = 1
                  iv(2)  = 3
                  val(1) = XYZ_coordinates(1)
                  val(2) = XYZ_coordinates(3)
!
               ELSE IF( subdomain%DependentVar == 'X' .OR. subdomain%DependentVar == 'x' ) THEN
!
                  dv     = 1
                  iv(1)  = 2
                  iv(2)  = 3
                  val(1) = XYZ_coordinates(2)
                  val(2) = XYZ_coordinates(3)
!
               END IF IF_DepVarCart
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
            ELSE
!
               IF_DepVarCyl: IF( subdomain%DependentVar == 'Z' .OR. subdomain%DependentVar == 'z' ) THEN
!
                  dv     = 3
                  iv(1)  = 1
                  val(1) = RZ_coordinates(1)
!
               ELSE IF( subdomain%DependentVar == 'R' .OR. subdomain%DependentVar == 'r' ) THEN
!
                  dv     = 1
                  iv(1)  = 3
                  val(1) = RZ_coordinates(2)
!
               END IF IF_DepVarCyl
!
            END IF IF_Coord
!
!
!***********************************************************************
!*                                                                     *
!*          Compute the limits of the independent variables            *
!*                                                                     *
!***********************************************************************
!
!
            IF_Coord2: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                     coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
               IF( XYZ_coordinates(iv(1)) < subdomain%LMin(iv(1)) .OR. XYZ_coordinates(iv(2)) < subdomain%LMin(iv(2)) .OR. &
     &             XYZ_coordinates(iv(1)) > subdomain%LMax(iv(1)) .OR. XYZ_coordinates(iv(2)) > subdomain%LMax(iv(2)) )    &
               THEN
                  RETURN
               END IF
!
               XYZ_limits( iv(1), 1 ) = subdomain%LMin(iv(1))
               XYZ_limits( iv(1), 2 ) = subdomain%LMax(iv(1))
!
               XYZ_limits( iv(2), 1 ) = subdomain%LMin(iv(2))
               XYZ_limits( iv(2), 2 ) = subdomain%LMax(iv(2))
!
               success = .TRUE.
!
            ELSE
!
               IF( iv(1) == 3 ) THEN
                  iv_mod = 2
               ELSE
                  iv_mod = 1
               END IF
!
               IF( RZ_coordinates(iv_mod) < subdomain%LMin(iv(1)) .OR. RZ_coordinates(iv_mod) > subdomain%LMax(iv(1)) ) THEN
                  RETURN
               END IF
!
               RZ_limits( iv_mod, 1 ) = subdomain%Lmin(iv(1))
               RZ_limits( iv_mod, 2 ) = subdomain%Lmax(iv(1))
!
               success = .TRUE.
!
            END IF IF_Coord2
!
!
!***********************************************************************
!*                                                                     *
!*          Compute the limits of the dependent variables              *
!*                                                                     *
!***********************************************************************
!
!
            IF( subdomain%units == 'IN' .OR. subdomain%units == 'In' .OR. subdomain%units == 'in') THEN
               conv_factor = 2.54d-2
            ELSE IF( subdomain%units == 'FT' .OR. subdomain%units == 'ft' .OR. subdomain%units == 'Ft') THEN
               conv_factor = 3.038d-1
            ELSE IF( subdomain%units == 'KM' .OR. subdomain%units == 'Km' .OR. subdomain%units == 'km') THEN
               conv_factor = 1.0d3
            ELSE IF( subdomain%units == 'MM' .OR. subdomain%units == 'Mm' .OR. subdomain%units == 'mm') THEN
               conv_factor = 1.0d-3
            ELSE
               conv_factor = 1.0d0
            END IF
!
            IF_Coord3: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                     coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For the FIRST surface
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( subdomain%TypeEqu1(1:4) == 'POLY' .OR. subdomain%TypeEqu1(1:4) == 'Poly' .OR. subdomain%TypeEqu1(1:4) == 'poly' ) THEN
!
! ************
! ************ For a POLYNOMIAL function
! ************
!
! ............... First, compute the polynomial of the 1st surface
!
                  m       = subdomain%OrderEquSurf1
                  arg1    = val(1) + subdomain%L1Shift(iv(1))
                  value1A = Polynomial(argument = arg1, n_poly = m, A = subdomain%Equ1CoeffA )
!
! ............... Then, compute the polynomial of the 2nd independent variable
!
                  arg2    = val(2) + subdomain%L1Shift(iv(2))
                  value1B = Polynomial(argument = arg2, n_poly = m, A = subdomain%Equ1CoeffB )
!
! ............... Finaly, compute the cross products
!
                  value1B = value1A + value1B
                  IF( m > 1) THEN
                     DO i = 1, m-1
                        value1B = value1B + subdomain%Equ1CoeffC(i) * ( arg1 ** (m-i) ) * ( arg2 ** (i) )
                     END DO
                  END IF
!
                  function_1p = subdomain%sign1 * conv_factor * value1B + subdomain%L1Shift(dv)
!
! ************
! ************ For a POWER function
! ************
!
               ELSE IF( subdomain%TypeEqu1(1:4) == 'POWE' .OR. subdomain%TypeEqu1(1:4) == 'Powe' .OR. subdomain%TypeEqu1(1:4) == 'powe' ) THEN
!
! ............... First, compute the powers of the first independent variable
!
                  m       = subdomain%OrderEquSurf1
                  arg1    = val(1) + subdomain%L1Shift(iv(1))
                  value1A = PowerSeries(argument = arg1, n_power = m, A = subdomain%Equ1CoeffA, B = subdomain%Equ1CoeffB, acceptable = valid )
!
! ............... Then, compute the powers of the 2nd independent variable
!
                  arg2    = val(2) + subdomain%L1Shift(iv(2))
                  value1B = PowerSeries(argument = arg2, n_power = m, A = subdomain%Equ1CoeffC, B = subdomain%Equ1CoeffD, acceptable = valid )
!
                  value1B     = value1A + value1B
                  function_1p = subdomain%sign1 * conv_factor * value1B + subdomain%L1Shift(dv)
!
! ************
! ************ For INTERPOLATION from a table
! ************
!
               ELSE IF( subdomain%TypeEqu1(1:4) == 'INTE' .OR. subdomain%TypeEqu1(1:4) == 'Inte' .OR. subdomain%TypeEqu1(1:4) == 'inte' ) THEN
!
                  m1 = subdomain%IntTableNum1
!
! >>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>> If there is no specific instruction to use the earlier-determined interpolations
! >>>>>>>>>>>>>>>
!
                  IF( PRESENT(find_RefSurface) ) THEN
                     IF( find_RefSurface ) do_not_do = .TRUE.
                  END IF
!
                  IF( .NOT. PRESENT(use_InterpData) .OR. do_not_do ) THEN
!
                     arg1 = val(1) + subdomain%L1Shift(iv(1))
                     arg2 = val(2) + subdomain%L1Shift(iv(2))
!
! .................. Then, determine the elevation of the 1st surface using table interpolation
!
                     NumData = IntTable(m1)%NumDataPoints
!
                     CALL Triangular_Interpolation( X0 = arg1,          &
     &                                              Y0 = arg2,          &
     &                                              Z0 = function_1p,   &
     &                                              R  = table_1%SearchRadius, &
     &                                              NumSearchPoints = table_1%NumSearchPoints, &
     &                                              TableNumber     = m1,              &
     &                                              NumData         = NumData,         &
     &                                              XD = IntTable(m1)%iv1(1:NumData),  &
     &                                              YD = IntTable(m1)%iv2(1:NumData),  &
     &                                              ZD = IntTable(m1)%dv( 1:NumData),  &
     &                                              Z_deriv = .TRUE., &
     &                                              dZdX = dZdX,      &
     &                                              dZdY = dZdY )

!
                     GO TO 1000
!
                  END IF
!
! >>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>> If the interpolated values have been computed earlier
! >>>>>>>>>>>>>>>
!
                  IF( completed_interpolation(m1) .AND. use_InterpData ) THEN
                     function_1p = available_InterpData(m1,1)
                     dZdX        = available_InterpData(m1,2)
                     dZdY        = available_InterpData(m1,3)
                  ELSE
!
! .................. Otherwise, first identify the table to be used for interpolation
!
                     arg1 = val(1) + subdomain%L1Shift(iv(1))
                     arg2 = val(2) + subdomain%L1Shift(iv(2))
!
! .................. Then, determine the elevation of the 1st surface using table interpolation
!
                     NumData = IntTable(m1)%NumDataPoints
!
                     CALL Triangular_Interpolation( X0 = arg1,          &
     &                                              Y0 = arg2,          &
     &                                              Z0 = function_1p,   &
     &                                              R  = table_1%SearchRadius, &
     &                                              NumSearchPoints = table_1%NumSearchPoints, &
     &                                              TableNumber     = m1,              &
     &                                              NumData         = NumData,         &
     &                                              XD = IntTable(m1)%iv1(1:NumData),  &
     &                                              YD = IntTable(m1)%iv2(1:NumData),  &
     &                                              ZD = IntTable(m1)%dv( 1:NumData),  &
     &                                              Z_deriv = .TRUE., &
     &                                              dZdX = dZdX,      &
     &                                              dZdY = dZdY )
!
                     IF( use_InterpData ) THEN
                        available_InterpData(m1,1)  = function_1p
                        available_InterpData(m1,2)  = dZdX
                        available_InterpData(m1,3)  = dZdY
                        completed_interpolation(m1) = .TRUE.
                     END IF
!
                  END IF
!
! ************
! ************ For FIXED width
! ************
!
               ELSE IF( subdomain%TypeEqu1(1:4) == 'FIXE' .OR. subdomain%TypeEqu1(1:4) == 'Fixe' .OR. subdomain%TypeEqu1(1:4) == 'fixe' ) THEN
!
! ............... When both surfaces are at a fixed distance from a reference surface
!
                  IF_Fixed: IF( subdomain%TypeEqu2(1:4) == 'FIXE' .OR. subdomain%TypeEqu2(1:4) == 'Fixe' .OR. subdomain%TypeEqu2(1:4) == 'fixe' ) THEN
!
                     m1 = subdomain%IntTableNum1
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>> If there is no specific instruction to use the earlier-determined interpolations
! >>>>>>>>>>>>>>>>>>
!
                     IF( .NOT. PRESENT(use_InterpData) ) THEN
!
                        arg1 = val(1) + subdomain%L1Shift(iv(1))
                        arg2 = val(2) + subdomain%L1Shift(iv(2))
!
! ..................... Then, determine the elevation of the 1st surface using table interpolation
!
                        NumData = IntTable(m1)%NumDataPoints
!
                        CALL Triangular_Interpolation( X0 = arg1,         &
     &                                                 Y0 = arg2,         &
     &                                                 Z0 = function_v,   &
     &                                                 R  = table_1%SearchRadius, &
     &                                                 NumSearchPoints = table_1%NumSearchPoints, &
     &                                                 TableNumber     = m1,              &
     &                                                 NumData         = NumData,         &
     &                                                 XD = IntTable(m1)%iv1(1:NumData),  &
     &                                                 YD = IntTable(m1)%iv2(1:NumData),  &
     &                                                 ZD = IntTable(m1)%dv( 1:NumData),  &
     &                                                 Z_deriv = .TRUE., &
     &                                                 dZdX = dZdX,      &
     &                                                 dZdY = dZdY )
!
                        GO TO 500
!
                     END IF
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>> If the interpolated values have been computed earlier
! >>>>>>>>>>>>>>>>>>
!
                     IF( completed_interpolation(m1) .AND. use_InterpData ) THEN
                        function_v = available_InterpData(m1,1)
                        dZdX       = available_InterpData(m1,2)
                        dZdY       = available_InterpData(m1,3)
                     ELSE
!
! ..................... Otherwise, first identify the table to be used for interpolation
!
                        arg1 = val(1) + subdomain%L1Shift(iv(1))
                        arg2 = val(2) + subdomain%L1Shift(iv(2))
!
                        NumData = IntTable(m1)%NumDataPoints
!
! ..................... Then, determine the elevation of the reference surface using table interpolation
!
                        CALL Triangular_Interpolation( X0 = arg1,          &
     &                                                 Y0 = arg2,          &
     &                                                 Z0 = function_v,    &
     &                                                 R  = table_1%SearchRadius, &
     &                                                 NumSearchPoints = table_1%NumSearchPoints, &
     &                                                 TableNumber     = m1,              &
     &                                                 NumData         = NumData,         &
     &                                                 XD = IntTable(m1)%iv1(1:NumData),  &
     &                                                 YD = IntTable(m1)%iv2(1:NumData),  &
     &                                                 ZD = IntTable(m1)%dv( 1:NumData) , &
     &                                                 Z_deriv = .TRUE., &
     &                                                 dZdX = dZdX,      &
     &                                                 dZdY = dZdY )
!
                        IF( use_InterpData ) THEN
                           available_InterpData(m1,1)  = function_v
                           available_InterpData(m1,2)  = dZdX
                           available_InterpData(m1,3)  = dZdY
                           completed_interpolation(m1) = .TRUE.
                        END IF
!
                     END IF
!
! .................. Compute the elevations of the two sufraces
!
  500                IF( subdomain%RefSurf1 == '0' .AND. subdomain%RefSurf2 == '0' ) THEN
!
                        IF( subdomain%VertOrStrat1 == 'V' .OR. subdomain%VertOrStrat1 == 'v' ) THEN
                           function_1p = function_v + subdomain%thick1
                        ELSE
                           function_1p = function_v + subdomain%thick1 / COS( ATAN(dZdX) ) / COS( ATAN(dZdY) )
                        END IF
!
                        IF( subdomain%VertOrStrat2 == 'V' .OR. subdomain%VertOrStrat2 == 'v' ) THEN
                           function_2p = function_v + subdomain%thick2
                        ELSE
                           function_2p = function_v + subdomain%thick2 / COS( ATAN(dZdX) ) / COS( ATAN(dZdY) )
                        END IF
!
                     ELSE IF( subdomain%RefSurf1 == '0' .AND. subdomain%RefSurf2 == '1' ) THEN
!
                        IF( subdomain%VertOrStrat1 == 'V' .OR. subdomain%VertOrStrat1 == 'v' ) THEN
                           function_1p = function_v + subdomain%thick1
                        ELSE IF( subdomain%VertOrStrat1 == 'C' .OR. subdomain%VertOrStrat1 == 'c' ) THEN
                           function_1p = function_v + subdomain%thick0 / COS( ATAN(dZdX) ) / COS( ATAN(dZdY) ) + subdomain%thick1
                        ELSE
                           function_1p = function_v + subdomain%thick1 / COS( ATAN(dZdX) ) / COS( ATAN(dZdY) )
                        END IF
!
                        IF( subdomain%VertOrStrat2 == 'V' .OR. subdomain%VertOrStrat2 == 'v' ) THEN
                           function_2p = function_1p + subdomain%thick2
                        ELSE
                           function_2p = function_1p + subdomain%thick2 * COS( ATAN(dZdX) ) * COS( ATAN(dZdY) )
                        END IF
!
                     ELSE IF( subdomain%RefSurf1 == '1' .AND. subdomain%RefSurf2 == '0' ) THEN
!
                        IF( subdomain%VertOrStrat2 == 'V' .OR. subdomain%VertOrStrat2 == 'v' ) THEN
                           function_2p = function_v + subdomain%thick2
                        ELSE IF( subdomain%VertOrStrat2 == 'C' .OR. subdomain%VertOrStrat2 == 'c' ) THEN
                           function_2p = function_v + subdomain%thick0 / COS( ATAN(dZdX) ) / COS( ATAN(dZdY) ) + subdomain%thick2
                        ELSE
                           function_2p = function_v + subdomain%thick2 / COS( ATAN(dZdX) ) / COS( ATAN(dZdY) )
                        END IF
!
                        IF( subdomain%VertOrStrat1 == 'V' .OR. subdomain%VertOrStrat1 == 'v' ) THEN
                           function_1p = function_2p + subdomain%thick1
                        ELSE
                           function_1p = function_2p + subdomain%thick1 / COS( ATAN(dZdX) ) / COS( ATAN(dZdY) )
                        END IF
!
                     END IF
!
                     GO TO 2000
!
                  END IF IF_Fixed
!
               END IF
!
! ............ For combination polynomial/power functions
!
               IF( subdomain%TypeEqu1(5:6) == '/P' .OR. subdomain%TypeEqu1(5:6) == '/p' ) THEN
                  function_1p = subdomain%L1Shift(dv) + subdomain%sign1 * conv_factor * ( value1B )**subdomain%expon1
               END IF
!
! ............ For combination polynomial/exponential functions
!
               IF( subdomain%TypeEqu1(5:6) == '/E' .OR. subdomain%TypeEqu1(5:6) == '/e' ) THEN
                  function_1p = subdomain%L1Shift(dv) + EXP( subdomain%sign1 * conv_factor * value1B )
               END IF
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> REPEAT for the 2nd surface
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
 1000          IF( do_not_do ) THEN      ! ... The second surface is not needed when determining the reference surface for thin areal grids
                  function_2p = 0.0d0
                  GO TO 2000
               END IF
!
               IF( subdomain%TypeEqu2(1:4) == 'POLY' .OR. subdomain%TypeEqu2(1:4) == 'Poly' .OR. subdomain%TypeEqu2(1:4) == 'poly' ) THEN
!
! ************
! ************ For a POLYNOMIAL function
! ************
!
! ............... First, compute the polynomial of the 2nd surface
!
                  m       = subdomain%OrderEquSurf2
                  arg1    = val(1) + subdomain%L2Shift(iv(1))
                  value2A = Polynomial(argument = arg1, n_poly = m, A = subdomain%Equ2CoeffA )
!
! ............... Then, compute the polynomial of the 2nd independent variable
!
                  arg2    = val(2) + subdomain%L2Shift(iv(2))
                  value2B = Polynomial(argument = arg2, n_poly = m, A = subdomain%Equ2CoeffB )
!
! ............... Finaly, compute the cross products
!
                  value2B = value2A + value2B
                  IF( m > 1) THEN
                     DO i = 1, m-1
                        value2B = value2B + subdomain%Equ2CoeffC(i) * ( arg1 ** (m-i) ) * ( arg2 ** (i) )
                     END DO
                  END IF
!
                  function_2p = subdomain%sign2 * conv_factor * value2B + subdomain%L2Shift(dv)
!
                  IF( ( subdomain%TypeEqu1(1:4) == 'FIXE' .OR. subdomain%TypeEqu1(1:4) == 'Fixe' .OR. subdomain%TypeEqu1(1:4) == 'fixe' ) .AND. &
     &                ( subdomain%VertOrStrat1 == 'V' .OR. subdomain%VertOrStrat1 == 'v' ) ) &
     &            THEN
                        function_1p = function_2p + subdomain%thick1
                  END IF
!
! ************
! ************ For a POWER function
! ************
!
               ELSE IF( subdomain%TypeEqu2(1:4) == 'POWE' .OR. subdomain%TypeEqu2(1:4) == 'Powe' .OR. subdomain%TypeEqu2(1:4) == 'powe' ) THEN
!
! ............... First, compute the powers of the first independent variable
!
                  m       = subdomain%OrderEquSurf2
                  arg1    = val(1) + subdomain%L2Shift(iv(1))
                  value2A = PowerSeries(argument = arg1, n_power = m, A = subdomain%Equ2CoeffA, B = subdomain%Equ2CoeffB, acceptable = valid )
!
! ............... Then, compute the powers of the 2nd independent variable
!
                  arg2    = val(2) + subdomain%L2Shift(iv(2))
                  value2B = PowerSeries(argument = arg2, n_power = m, A = subdomain%Equ2CoeffC, B = subdomain%Equ2CoeffD, acceptable = valid )
!
                  value2B     = value2A + value2B
                  function_2p = subdomain%sign2 * conv_factor * value2B + subdomain%L2Shift(dv)
!
                  IF( ( subdomain%TypeEqu1(1:4) == 'FIXE' .OR. subdomain%TypeEqu1(1:4) == 'Fixe' .OR. subdomain%TypeEqu1(1:4) == 'fixe' ) .AND. &
     &                ( subdomain%VertOrStrat1 == 'V' .OR. subdomain%VertOrStrat1 == 'v' ) ) &
     &            THEN
                        function_1p = function_2p + subdomain%thick1
                  END IF
!
! ************
! ************ For INTERPOLATION from a table
! ************
!
               ELSE IF( subdomain%TypeEqu2(1:4) == 'INTE' .OR. subdomain%TypeEqu2(1:4) == 'Inte' .OR. subdomain%TypeEqu2(1:4) == 'inte' ) THEN
!
                  m2   = subdomain%IntTableNum2
!
! >>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>> If there is no specific instruction to use the earlier-determined interpolations
! >>>>>>>>>>>>>>>
!
                  IF( .NOT. PRESENT(use_InterpData) ) THEN
!
                     arg1 = val(1) + subdomain%L2Shift(iv(1))
                     arg2 = val(2) + subdomain%L2Shift(iv(2))
!
! .................. Then, determine the elevation of the 1st surface using table interpolation
!
                     NumData = IntTable(m2)%NumDataPoints
!
                     CALL Triangular_Interpolation( X0 = arg1,          &
     &                                              Y0 = arg2,          &
     &                                              Z0 = function_2p,   &
     &                                              R  = table_2%SearchRadius, &
     &                                              NumSearchPoints = table_2%NumSearchPoints, &
     &                                              TableNumber     = m2,              &
     &                                              NumData         = NumData,         &
     &                                              XD = IntTable(m2)%iv1(1:NumData),  &
     &                                              YD = IntTable(m2)%iv2(1:NumData),  &
     &                                              ZD = IntTable(m2)%dv( 1:NumData),  &
     &                                              Z_deriv = .TRUE., &
     &                                              dZdX = dZdX,      &
     &                                              dZdY = dZdY )

!
                     GO TO 2000
!
                  END IF
!
! >>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>> If the interpolated values have been computed earlier
! >>>>>>>>>>>>>>>
!
                  IF( completed_interpolation(m2) .AND. use_InterpData ) THEN
                     function_2p = available_InterpData(m2,1)
                     dZdX        = available_InterpData(m2,2)
                     dZdY        = available_InterpData(m2,3)
                  ELSE
!
! .................. Otherwise, first identify the table to be used for interpolation
!
                     arg1 = val(1) + subdomain%L2Shift(iv(1))
                     arg2 = val(2) + subdomain%L2Shift(iv(2))
!
! .................. Then, determine the elevation of the 2nd surface using table interpolation
!
                     NumData = IntTable(m2)%NumDataPoints
!
                     CALL Triangular_Interpolation( X0 = arg1,          &
     &                                              Y0 = arg2,          &
     &                                              Z0 = function_2p,   &
     &                                              R  = table_2%SearchRadius, &
     &                                              NumSearchPoints = table_2%NumSearchPoints, &
     &                                              TableNumber     = m2,              &
     &                                              NumData         = NumData,         &
     &                                              XD = IntTable(m2)%iv1(1:NumData),  &
     &                                              YD = IntTable(m2)%iv2(1:NumData),  &
     &                                              ZD = IntTable(m2)%dv( 1:NumData),  &
     &                                              Z_deriv = .TRUE., &
     &                                              dZdX = dZdX,      &
     &                                              dZdY = dZdY )
!
                     IF( use_InterpData ) THEN
                        available_InterpData(m2,1)  = function_2p
                        available_InterpData(m2,2)  = dZdX
                        available_InterpData(m2,3)  = dZdY
                        completed_interpolation(m2) = .TRUE.
                     END IF
!
                  END IF
!
! ............... When the 1st surface is at a fixed distance from the 2nd
!
                  IF( ( subdomain%TypeEqu1(1:4) == 'FIXE' .OR. subdomain%TypeEqu1(1:4) == 'Fixe' .OR. subdomain%TypeEqu1(1:4) == 'fixe' ) ) THEN
!
                     IF( subdomain%VertOrStrat1 == 'S' .OR. subdomain%VertOrStrat1 == 's' ) THEN
                        function_1p = function_2p + subdomain%thick1 / COS( ATAN(dZdX) ) / COS( ATAN(dZdY) )
                     ELSE IF( subdomain%VertOrStrat1 == 'V' .OR. subdomain%VertOrStrat1 == 'v' ) THEN
                        function_1p = function_2p + subdomain%thick1
                     END IF
!
                     GO TO 2000
!
                  END IF
!
! ************
! ************ For FIXED width
! ************
!
               ELSE IF( subdomain%TypeEqu2(1:4) == 'FIXE' .OR. subdomain%TypeEqu2(1:4) == 'Fixe' .OR. subdomain%TypeEqu2(1:4) == 'fixe' ) THEN
!
! ............... When the 2nd surface is at a fixed distance from the 1st
!
                  IF( subdomain%VertOrStrat2 == 'S' .OR. subdomain%VertOrStrat2 == 's' ) THEN
                     function_2p = function_1p + subdomain%thick2 / COS( ATAN(dZdX) ) / COS( ATAN(dZdY) )  ! Adjustment for steepness
                  ELSE IF( subdomain%VertOrStrat2 == 'V' .OR. subdomain%VertOrStrat2 == 'v' ) THEN
                     function_2p = function_1p + subdomain%thick2
                  END IF
!
!
               END IF
!
! ............ For combination polynomial/power functions
!
               IF( subdomain%TypeEqu2(5:6) == '/P' .OR. subdomain%TypeEqu2(5:6) == '/p' ) THEN
                  function_2p = subdomain%L2Shift(dv) + subdomain%sign2 * conv_factor * ( value2B )**subdomain%expon2
               END IF
!
! ............ For combination polynomial/exponential functions
!
               IF( subdomain%TypeEqu2(5:6) == '/E' .OR. subdomain%TypeEqu2(5:6) == '/e' ) THEN
                  function_2p = subdomain%L2Shift(dv) + EXP( subdomain%sign2 * conv_factor * value2B )
               END IF
!
 2000          XYZ_limits( dv, 1 ) = MIN( function_1p, function_2p )
               XYZ_limits( dv, 2 ) = MAX( function_1p, function_2p )
!
! >>>>>>>>>
!
            ELSE
!
! >>>>>>>>>
!
! ............ First surface
!
               IF( subdomain%TypeEqu1(1:4) == 'POLY' .OR. subdomain%TypeEqu1(1:4) == 'Poly' .OR. subdomain%TypeEqu1(1:4) == 'poly' ) THEN
!
                  m       = subdomain%OrderEquSurf1
                  arg1    = val(1) + subdomain%L1Shift(iv(1))
                  value1A = Polynomial(argument = arg1, n_poly = m, A = subdomain%Equ1CoeffA )
                  function_1p = subdomain%L1Shift(dv) + subdomain%sign1 * conv_factor * value1A
!
               ELSE IF( subdomain%TypeEqu1(1:4) == 'POWE' .OR. subdomain%TypeEqu1(1:4) == 'Powe' .OR. subdomain%TypeEqu1(1:4) == 'powe' ) THEN
!
                  m       = subdomain%OrderEquSurf1
                  arg1    = val(1) + subdomain%L1Shift(iv(1))
                  value1A = PowerSeries(argument = arg1, n_power = m, A = subdomain%Equ1CoeffA, B = subdomain%Equ1CoeffB, acceptable = valid )
                  function_1p = subdomain%L1Shift(dv) + subdomain%sign1 * conv_factor * value1A
!
               END IF
!
! ............ For combination polynomial/power functions
!
               IF( subdomain%TypeEqu1(5:6) == '/P' .OR. subdomain%TypeEqu1(5:6) == '/p' ) THEN
                  function_1p = subdomain%L1Shift(dv) + subdomain%sign1 * conv_factor * ( value1A )**subdomain%expon1
               END IF
!
! ............ For combination polynomial/exponential functions
!
               IF( subdomain%TypeEqu1(5:6) == '/E' .OR. subdomain%TypeEqu1(5:6) == '/e' ) THEN
                  function_1p = subdomain%L1Shift(dv) + EXP( subdomain%sign1 * conv_factor * value1A )
               END IF
!
! ............ Second surface
!
               IF( subdomain%TypeEqu2(1:4) == 'POLY' .OR. subdomain%TypeEqu2(1:4) == 'Poly' .OR. subdomain%TypeEqu2(1:4) == 'poly' ) THEN
!
                  m       = subdomain%OrderEquSurf2
                  arg1    = val(1) - subdomain%L2Shift(iv(1))
                  value2A = Polynomial(argument = arg1, n_poly = m, A = subdomain%Equ2CoeffA )
                  function_2p = subdomain%L2Shift(dv) + subdomain%sign2 * conv_factor * value2A
!
               ELSE IF( subdomain%TypeEqu2(1:4) == 'POWE' .OR. subdomain%TypeEqu2(1:4) == 'Powe' .OR. subdomain%TypeEqu2(1:4) == 'powe' ) THEN
!
                  m       = subdomain%OrderEquSurf2
                  arg1    = val(1) - subdomain%L2Shift(iv(1))
                  value2A = PowerSeries(argument = arg1, n_power = m, A = subdomain%Equ2CoeffA, B = subdomain%Equ2CoeffB, acceptable = valid )
                  function_2p = subdomain%L2Shift(dv) + subdomain%sign2 * conv_factor * value2A
!
               END IF
!
! ............ For combination polynomial/power functions
!
               IF( subdomain%TypeEqu2(5:6) == '/P' .OR. subdomain%TypeEqu2(5:6) == '/p' ) THEN
                  function_2p = subdomain%L2Shift(dv) + subdomain%sign2 * conv_factor * ( value2A )**subdomain%expon2
               END IF
!
! ............ For combination polynomial/exponential functions
!
               IF( subdomain%TypeEqu2(5:6) == '/E' .OR. subdomain%TypeEqu2(5:6) == '/e' ) THEN
                  function_2p = subdomain%L2Shift(dv) + EXP( subdomain%sign2 * conv_factor * value2A )
               END IF
!
               IF( dv == 3 ) THEN
                  dv_mod = 2
               ELSE
                  dv_mod = 1
               END IF
!
               RZ_limits( dv_mod, 1 ) = MIN( function_1p, function_2p )
               RZ_limits( dv_mod, 2 ) = MAX( function_1p, function_2p )
!
            END IF IF_Coord3
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Surface_Coordinate_Limits 1.0 .......... 20 April     2015',6X,'Computing the coordinate ranges between bounding region surfaces')
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Surface_Coordinate_Limits>
!
!
            RETURN
!
!
!
         END SUBROUTINE Surface_Coordinate_Limits
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Triangular_Interpolation(X0, Y0, Z0, R, NumSearchPoints, TableNumber, NumData, XD, YD, ZD, Z_deriv, dZdX, dZdY, Y_deriv, dYdX )
!
            USE Array_Ranking_Procedures
!
            IMPLICIT NONE
!
! ......... Integer variables
!
            INTEGER, INTENT(IN) :: TableNumber, NumData, NumSearchPoints
!
            INTEGER :: i, third, num_within_R, num_checked_points
            INTEGER :: T1, T2, T3
!
! ......... Integer arrays
!
            INTEGER, DIMENSION(NumData) :: DistRank
!
            INTEGER, DIMENSION(20) :: num_SearchPoints
!
! ......... Double precision variables
!
            REAL(KIND = 8), INTENT(IN)  :: X0, Y0, R
!
            REAL(KIND = 8), INTENT(OUT) :: Z0
!
            REAL(KIND = 8), INTENT(OUT), OPTIONAL :: dZdX, dZdY, dYdX
!
            REAL(KIND = 8) :: R2, XX, YY, RT, Dmax
            REAL(KIND = 8) :: A, s, t, PP1, PP2, PP3, DET, AA1, AA2, AA3
!
! ......... Double precision arrays
!
            REAL(KIND = 8), DIMENSION(20) :: Xmin, Xmax, Ymin, Ymax
!
            REAL(KIND = 8), DIMENSION(NumData) :: distance
!
            REAL(KIND = 8), INTENT(IN), DIMENSION(NumData) :: XD, YD, ZD
!
! ......... Logical variables
!
            LOGICAL, DIMENSION(20) :: limits_determined
!
            LOGICAL, INTENT(IN), OPTIONAL :: Z_deriv, Y_deriv
!
            LOGICAL :: First_call = .TRUE., extrapolation, last
!
! -------------
! ......... Saving variables
! -------------
!
            SAVE First_call, Xmin, Xmax, Ymin, Ymax
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of <Triangular_Interpolation>
!
!
            IF(First_call) THEN
               WRITE( UNIT = *, FMT = 6000 )
               First_call   = .FALSE.
               Xmin = 0.0d0
               Xmax = 0.0d0
               Ymin = 0.0d0
               Ymax = 0.0d0
               limits_determined = .FALSE.
            END IF
!
            IF( NumData < 3 ) THEN
               WRITE(UNIT = *, FMT = 6050 )
               STOP
            END IF
!
! ......... Initialization
!
            distance = 0.0d0
!
            IF( .NOT. limits_determined(TableNumber) ) THEN
!
               Xmin(TableNumber) = MINVAL(XD)
               Xmax(TableNumber) = MAXVAL(XD)
               Ymin(TableNumber) = MINVAL(YD)
               Ymax(TableNumber) = MAXVAL(YD)
!
               IF( NumSearchPoints <= 30 ) THEN
                  num_SearchPoints(TableNumber) = 30
               ELSE
                  num_SearchPoints(TableNumber) = NumSearchPoints
               END IF
!
               limits_determined(TableNumber) = .TRUE.
!
            END IF
!
            extrapolation = .FALSE.
            IF( X0 < Xmin(TableNumber) .OR. X0 > Xmax(TableNumber) .OR. Y0 < Ymin(TableNumber) .OR. Y0 > Ymax(TableNumber) ) THEN
               extrapolation = .TRUE.
            END IF
!
! ......... The search radius is specified
!
            RT = R
!
            XX   = MAX( ABS(X0 - Xmin(TableNumber)), ABS(X0 - Xmax(TableNumber)) )
            YY   = MAX( ABS(Y0 - Ymin(TableNumber)), ABS(Y0 - Ymax(TableNumber)) )
            Dmax = SQRT( XX*XX + YY*YY )
!
            IF( extrapolation ) THEN
               RT = Dmax
            ELSE
               IF( RT == 0.0d0 ) RT = 5.0d-2 * Dmax
            END IF
!
            R2 = RT * RT
!
! ......... Determine the vector of distances from the table points
!
            FORALL (i = 1:NumData) distance(i) = (XD(i) - X0) * (XD(i) - X0) + (YD(i) - Y0) * (YD(i) - Y0)
!
! ......... Rank the <distance> array
!
            CALL Rank_Array( distance, DistRank )
!
! ......... Find the triangle within which lies the interpolation point
!
            T1 = DistRank(1)
            T2 = DistRank(2)
!
            num_within_R = COUNT(distance <= R2)
            num_checked_points = MIN( num_within_R, num_SearchPoints(TableNumber) )
!
            third = 3
            last  = .FALSE.
!
 1000       T3 = DistRank(third)
!
! ......... The case of extrapolation
!
            IF( extrapolation ) GO TO 2000
!
! ......... First compute the area of the potential triangle surrounding the point
!
            A = XD(T2) * YD(T3) - XD(T3) * YD(T2) + YD(T1) * ( XD(T3) - XD(T2) ) + XD(T1) * ( YD(T2) - YD(T3) )
!
! ......... Check if the 3 points are co-linear
!
            IF ( ABS(A) < 1.0d-4 ) THEN
               IF( third < NumData ) THEN
                  third = third + 1
                  GO TO 1000
               ELSE
                  WRITE( UNIT = *, FMT = 6055 )
                  STOP
               END IF
            END IF
!
            IF( last ) GO TO 2000
!
! ......... Use barycentric coordinates to compute the checking criteria for being an internal point
!
            s = ( YD(T1) * XD(T3) - XD(T1) * YD(T3) + X0 * ( YD(T3) - YD(T1) ) + Y0 * ( XD(T1) - XD(T3) ) ) / A
!
            t = ( YD(T2) * XD(T1) - XD(T2) * YD(T1) + X0 * ( YD(T1) - YD(T2) ) + Y0 * ( XD(T2) - XD(T1) ) ) / A
!
! ......... Check whether the point is internal
!
            IF ( s >= 0.0d0 .AND. t >= 0.0d0 .AND. (1.0d0 - s - t ) >= 0.0d0 ) THEN
               CONTINUE
            ELSE
               IF( third < num_checked_points ) THEN
                  third = third + 1
                  GO TO 1000
               ELSE IF( third == num_within_R .AND. num_checked_points == num_within_R ) THEN
                  IF( num_SearchPoints(TableNumber) > num_within_R ) THEN
                     RT = 1.25d0 * RT
                     R2 = RT * RT
                     num_within_R = MIN( COUNT(distance <= R2), num_SearchPoints(TableNumber) )
                     num_checked_points = MIN( num_within_R, num_SearchPoints(TableNumber) )
                     IF( num_within_R <= num_SearchPoints(TableNumber) .AND. third < num_SearchPoints(TableNumber) ) THEN
                        third = third + 1
                        GO TO 1000
                     END IF
                  ELSE
                     third = 3
                     last  = .TRUE.
                     GO TO 1000
                  END IF
               ELSE IF( third >= num_SearchPoints(TableNumber) ) THEN
                  third = 3
                  last  = .TRUE.
                  GO TO 1000
               END IF
            END IF
!
 2000       PP1 = XD(T2)*YD(T3) - XD(T3)*YD(T2)
            PP2 = XD(T3)*YD(T1) - XD(T1)*YD(T3)
            PP3 = XD(T1)*YD(T2) - XD(T2)*YD(T1)
!
            DET = PP1 + PP2 + PP3
!
            IF_Dtm: IF(ABS(DET) <= 1.0d-7) THEN    ! ... Co-linearity of points!
               IF( third < num_checked_points ) THEN
                  third = third + 1
                  GO TO 1000
               ELSE
                  WRITE(UNIT = *, FMT = 6060) X0, Y0
                  STOP
               END IF
            END IF IF_Dtm
!
            AA1 = PP1 * ZD(T1) + PP2 * ZD(T2) + PP3 * ZD(T3)
            AA2 = ( YD(T2) - YD(T3)) * ZD(T1) + (YD(T3) - YD(T1)) * ZD(T2) + (YD(T1) - YD(T2) ) * ZD(T3)
            AA3 = ( XD(T3) - XD(T2)) * ZD(T1) + (XD(T1) - XD(T3)) * ZD(T2) + (XD(T2) - XD(T1) ) * ZD(T3)
!
            Z0 = ( AA1 + AA2 * X0 + AA3 * Y0 ) / DET
!
! >>>>>>>>>
! >>>>>>>>> Z-derivatives
! >>>>>>>>>
!
            IF( PRESENT(Z_deriv) ) THEN
               IF( Z_deriv ) THEN
                  dZdX = AA2 / DET
                  dZdY = AA3 / DET
               END IF
            END IF
!
! >>>>>>>>>
! >>>>>>>>> Y-derivatives
! >>>>>>>>>
!
            IF( PRESENT(Y_deriv) ) THEN
!
               IF( Y_deriv ) THEN
!
                  PP1 = XD(T2)*YD(T3) - XD(T3)*YD(T2)
                  PP2 = XD(T3)*YD(T1) - XD(T1)*YD(T3)
                  PP3 = XD(T1)*YD(T2) - XD(T2)*YD(T1)
!
                  DET = PP1 + PP2 + PP3
!
                  dYdX = 0.0d0
!
               END IF
!
            END IF
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Triangular_Interpolation 1.0 ............. 6 June     2015',6X,'Interpolation over triangular surfaces defined by irregularly-spaced data')
!
 6050 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Triangular_Interpolation>: The interpolation table has fewer than 3 entries',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6055 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Triangular_Interpolation>: The table includes only co-linear data points',/, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6060 FORMAT(//,25('ERROR-'),//, &
     &       T10,'>>>  Procedure <Triangular_Interpolation>: Unable to determine the triangle that includes the point (dV1=',ES11.4,',dV2=',ES11.4,')',/, &
     &       T10,'>>>                         Likely cause : Improper definition of the exclusion zones',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,25('ERROR-'))
!
            RETURN
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Triangular_Interpolation>
!
!
         END SUBROUTINE Triangular_Interpolation
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      END MODULE Subdomain_Bounding_Surfaces
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      MODULE Het_Region_Definition
!
         USE Het_Region_Variables
!
         SAVE
!
! >>>>>
! >>>>>
! >>>>>
! >>>>>
! >>>>>
!
         CONTAINS
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         CHARACTER(LEN = 5) FUNCTION Media_Region( coordinates, dominant_medium, X, Y, Z, media_by_number, emissivity, &
     &                                             use_InterpData, completed_interpolation, available_InterpData )
!
            USE Subdomain_Bounding_Surfaces, ONLY: Surface_Coordinate_Limits
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*             Routine determining the media region where              *
!*                     a particular cell belongs                       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... User-defined variables
! -------------
!
            TYPE(Tabular) :: table_1, table_2
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: phi, theta, omega, omega1, omega2, height
!
            REAL(KIND = 8), OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables,3) :: available_InterpData
!
            REAL(KIND = 8), DIMENSION(3,2) :: XYZ_limits
            REAL(KIND = 8), DIMENSION(2,2) :: RZ_limits
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: X_1p, Y_1p, Z_1p, X_2p, Z_2p, Rtest, long_axis, short_axis, radius, Rmin, RMax, Delta
!
            REAL(KIND = 8), INTENT(IN)  :: X, Y, Z
            REAL(KIND = 8), INTENT(OUT) :: emissivity
!
! -------------
! ......... Double precision parameters
! -------------
!
            REAL(KIND = 8), PARAMETER :: pi = 3.14159265358979324D0
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: n, N_0, N_1, N_2, N_3, N_4, N_5, m1, m2, ier
!
! -------------
! ......... CHARACTER variables
! -------------
!
            CHARACTER(LEN =  5), INTENT(IN) :: dominant_medium
            CHARACTER(LEN = 11), INTENT(IN) :: coordinates
!
! -------------
! ......... Logical variables
! -------------
!
            LOGICAL, INTENT(IN) :: media_by_number
!
            LOGICAL, OPTIONAL, INTENT(IN) :: use_InterpData
!
            LOGICAL, OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables) :: completed_interpolation
!
            LOGICAL :: First_call = .TRUE., correct_region, interp1, interp2
!
! -------------
! ......... Saving variables
! -------------
!
            SAVE First_call, phi, theta, omega, omega1, omega2, height
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of <Media_Region>
!
!
            IF(First_call) THEN
               WRITE(UNIT = *, FMT = 6000)
            ELSE
               GO TO 1000
            END IF
!
            ALLOCATE( phi(Num_HetRegions),   theta(Num_HetRegions),  height(Num_HetRegions), &
     &                omega(Num_HetRegions), omega1(Num_HetRegions), omega2(Num_HetRegions), STAT=ier )
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6850) 'phi, theta, height, omega, omega1, omega2', 'Media_Region'
               STOP
            END IF
!
! ......... Initialization
!
            phi    = 0.0d0
            theta  = 0.0d0
            height = 0.0d0
!
            omega  = 0.0d0
            omega1 = 0.0d0
            omega2 = 0.0d0
!
            emissivity = 0.0d0
!
!
!***********************************************************************
!*                                                                     *
!*          Baseline computations for cylindrical shapes in            *
!*          cartesian coordinate systems                               *
!*                                                                     *
!***********************************************************************
!
!
            IF_Initial: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                      coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
               DO_Regions: DO n=2,TotNum_HetRegions
!
! ............... Cartesian grid system, cylindrically-shaped media regions
!
                  IF_Shape: IF( (Region(n)%shape(1:2) == 'CY' .OR. Region(n)%shape(1:2) == 'cy' .OR. &
     &                           Region(n)%shape(1:2) == 'Cy') )  &
     &            THEN
!
                     CALL Determine_Orientation_Parameters( base1_center_coord = Region(n)%CylBase1Coord, &
     &                                                      base2_center_coord = Region(n)%CylBase2Coord, &
     &                                                      radius1 = Region(n)%CylBase1R, radius2 = Region(n)%CylBase2R, &
     &                                                      phi = phi(n), theta = theta(n), omega = omega(n), height = height(n) )
!
! ............... Cartesian grid system, elliptically-shaped media regions
!
                  ELSE IF( (Region(n)%shape(1:2) == 'EL' .OR. Region(n)%shape(1:2) == 'el' .OR. &
     &                      Region(n)%shape(1:2) == 'El') )  &
     &            THEN
!
                     CALL Determine_Orientation_Parameters( base1_center_coord = Region(n)%EllBase1Coord, &
     &                                                      base2_center_coord = Region(n)%EllBase2Coord, &
     &                                                      LAxis1 = Region(n)%Base1LAxis, SAxis1 = Region(n)%Base1SAxis, &
     &                                                      LAxis2 = Region(n)%Base2LAxis, SAxis2 = Region(n)%Base2SAxis, &
     &                                                      phi = phi(n), theta = theta(n), omega1 = omega1(n), omega2 = omega2(n), height = height(n) )
!
                  END IF IF_Shape
!
               END DO DO_Regions
!
            END IF IF_Initial
!
            First_call = .FALSE.
!
! -------------
! ......... Initialization
! -------------
!
 1000       Media_Region = dominant_medium  ! Initial assignment to the first domain (considered dominant by default)
            emissivity   = dominant_medium_emissivity

!
            IF(media_by_number) THEN
               IF( dominant_medium(1:1) == 'F' .OR. dominant_medium(1:1) == 'f' ) THEN
                  Media_Region = ' 0001'
               ELSE
                  Media_Region = '-0001'
               END IF
            END IF
!
!***********************************************************************
!*                                                                     *
!*          For cartesian coordinate systems                           *
!*                                                                     *
!***********************************************************************
!
!
            IF_Assign: IF(coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                    coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
               DO_Cartesian: DO n=2,TotNum_HetRegions
!
! ............... Cartesian grid system, Cartesian-coordinate-defined or periodic media regions
!
                  IF_RShape: IF( Region(n)%shape(1:4) == 'RECT' .OR. Region(n)%shape(1:4) == 'Rect' .OR. Region(n)%shape(1:4) == 'rect' .OR. &
     &                           Region(n)%shape(1:2) == 'PE'   .OR. Region(n)%shape(1:2) == 'Pe'   .OR. Region(n)%shape(1:2) == 'pe' )      &
                  THEN
!
                     IF( X >= Region(n)%LMin(1) .AND. X <= Region(n)%LMax(1) .AND.   &
     &                   Y >= Region(n)%LMin(2) .AND. Y <= Region(n)%LMax(2) .AND.   &
     &                   Z >= Region(n)%LMin(3) .AND. Z <= Region(n)%LMax(3) )       &
     &               THEN
!
                        Media_Region = Region(n)%name
                        emissivity   = Region(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = MediaSequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( Media_Region(1:1) == 'F' .OR. Media_Region(1:1) == 'f' ) THEN
                              Media_Region = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              Media_Region = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, cylindrically-shaped media regions
!
                  ELSE IF( (Region(n)%shape(1:4) == 'CYLI' .OR. Region(n)%shape(1:4) == 'cyli' .OR. Region(n)%shape(1:4) == 'Cyli') ) THEN
!
! .................. The coordinates in the 1st rotated system
!
                     X_1p =  (X - Region(n)%CylBase1Coord(1)) * COS(phi(n)) + (Y - Region(n)%CylBase1Coord(2)) * SIN(phi(n))
                     Y_1p = -(X - Region(n)%CylBase1Coord(1)) * SIN(phi(n)) + (Y - Region(n)%CylBase1Coord(2)) * COS(phi(n))
                     Z_1p =   Z - Region(n)%CylBase1Coord(3)
!
! .................. The coordinates in the 2nd rotated system
!
                     X_2p =  X_1p * COS(theta(n)) + Z_1p * SIN(theta(n))
                     Z_2p = -X_1p * SIN(theta(n)) + Z_1p * COS(theta(n))
!                    Y_2p =  Y_1p
!
                     IF( Region(n)%CylBase1R /= 0.0d0 ) radius = MAX(Region(n)%CylBase1R, Region(n)%CylBase2R) - X_2p * omega(n)
!
                     Rtest = SQRT( Z_2p * Z_2p + Y_1p * Y_1p )
!
                     IF( ( Region(n)%CylBase1R == 0.0d0 .AND. X_2p >= 0.0d0 .AND. X_2p <= height(n) .AND. Rtest >= Region(n)%CylRmin .AND. Rtest <= Region(n)%CylRmax ) .OR. &
     &                   ( Region(n)%CylBase1R /= 0.0d0 .AND. X_2p >= 0.0d0 .AND. X_2p <= height(n) .AND. Rtest <= radius ) ) &
     &               THEN
!
                        Media_Region = Region(n)%name
                        emissivity   = Region(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = MediaSequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( Media_Region(1:1) == 'F' .OR. Media_Region(1:1) == 'f' ) THEN
                              Media_Region = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              Media_Region = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, elliptical media regions
!
                  ELSE IF( Region(n)%shape(1:4) == 'ELLI' .OR. Region(n)%shape(1:4) == 'Elli' .OR. Region(n)%shape(1:4) == 'elli' ) THEN
!
! .................. The angle between the original x- or y-axis and the long axis of the ellipse (rotated system)
!
                     phi(n) = Region(n)%LAxisAngle
!
                     SELECT CASE( Region(n)%EllipsePlane )
!
                     CASE( 'XY', 'xy', 'Xy', 'xY' )
!
                        IF( Region(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = Z - Region(n)%EllBase1Coord(3)
                           long_axis  = Region(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = Region(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = Region(n)%LAxis
                           short_axis = Region(n)%SAxis
                        END IF
!
                        X_1p =  (X - Region(n)%EllBase1Coord(1)) * COS(phi(n)) + (Y - Region(n)%EllBase1Coord(2)) * SIN(phi(n))
                        Y_1p = -(X - Region(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Y - Region(n)%EllBase1Coord(2)) * COS(phi(n))
!
                        Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Y_1p * Y_1p / ( short_axis * short_axis )
!
                     CASE( 'XZ', 'xz', 'Xz', 'xZ' )
!
                        IF( Region(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = Y - Region(n)%EllBase1Coord(2)
                           long_axis  = Region(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = Region(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = Region(n)%LAxis
                           short_axis = Region(n)%SAxis
                        END IF
!
                        X_1p =  (X - Region(n)%EllBase1Coord(1)) * COS(phi(n)) + (Z - Region(n)%EllBase1Coord(3)) * SIN(phi(n))
                        Z_1p = -(X - Region(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Z - Region(n)%EllBase1Coord(3)) * COS(phi(n))
!
                        Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     CASE( 'YZ', 'yz', 'Yz', 'yZ' )
!
                        IF( Region(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = X - Region(n)%EllBase1Coord(1)
                           long_axis  = Region(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = Region(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = Region(n)%LAxis
                           short_axis = Region(n)%SAxis
                        END IF
!
                        Y_1p =  (Y - Region(n)%EllBase1Coord(2)) * COS(phi(n)) + (Z - Region(n)%EllBase1Coord(3)) * SIN(phi(n))
                        Z_1p = -(Y - Region(n)%EllBase1Coord(2)) * SIN(phi(n)) + (Z - Region(n)%EllBase1Coord(3)) * COS(phi(n))
!
                        Rtest = Y_1p * Y_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     END SELECT
!
                     IF( Delta >= 0.0d0 .AND. Delta <= height(n) .AND. Rtest <= 2.5d-1 ) THEN
!
                        Media_Region = Region(n)%name
                        emissivity   = Region(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = MediaSequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( Media_Region(1:1) == 'F' .OR. Media_Region(1:1) == 'f' ) THEN
                              Media_Region = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              Media_Region = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, spherically-shaped media regions
!
                  ELSE IF( (Region(n)%shape(1:4) == 'SPHE' .OR. Region(n)%shape(1:4) == 'sphe' .OR. Region(n)%shape(1:4) == 'Sphe') ) THEN
!
                     X_1p = (X - Region(n)%SphereCenterCoord(1))
                     Y_1p = (Y - Region(n)%SphereCenterCoord(2))
                     Z_1p = (Z - Region(n)%SphereCenterCoord(3))
!
                     Rtest = SQRT( X_1p * X_1p + Y_1p * Y_1p + Z_1p * Z_1p )
!
                     IF( Rtest >= Region(n)%SphRmin .AND. Rtest <= Region(n)%SphRmax ) THEN
!
                        Media_Region = Region(n)%name
                        emissivity   = Region(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = MediaSequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( Media_Region(1:1) == 'F' .OR. Media_Region(1:1) == 'f' ) THEN
                              Media_Region = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              Media_Region = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, irregularly-shaped media regions
!
                  ELSE IF( Region(n)%shape(1:4) == 'IRRE' .OR. Region(n)%shape(1:4) == 'irre' .OR. Region(n)%shape(1:4) == 'Irre' ) THEN
!
                     m1      = Region(n)%IntTableNum1
                     interp1 = .FALSE.
                     IF( m1 > 0 ) THEN
                        table_1 = IntTable(m1)
                        IF( completed_interpolation(m1) ) interp1 = .TRUE.
                     END IF
!
                     m2      = Region(n)%IntTableNum2
                     interp2 = .FALSE.
                     IF( m2 > 0 ) THEN
                        table_2 = IntTable(m2)
                        IF( completed_interpolation(m2) ) interp2 = .TRUE.
                     END IF
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
!
                     IF_tables: IF( m1 == 0 .AND. m2 == 0 ) THEN
!
                        CALL Surface_Coordinate_Limits( subdomain        = Region(n),     &
     &                                                  XYZ_coordinates  = (/ X, Y, Z /), &
     &                                                  XYZ_limits       = XYZ_limits,    &
     &                                                  success          = correct_region )
!
                     ELSE
!
                        IF( (PRESENT(use_InterpData)) .AND. (interp1 .OR. interp2) ) THEN
!
                           CALL Surface_Coordinate_Limits( subdomain        = Region(n),      &
     &                                                     table_1          = table_1,        &
     &                                                     table_2          = table_2,        &
     &                                                     XYZ_coordinates  = (/ X, Y, Z /),  &
     &                                                     XYZ_limits       = XYZ_limits,     &
     &                                                     success          = correct_region, &
     &                                                     use_InterpData   = use_InterpData, &
     &                                                     completed_interpolation = completed_interpolation, &
     &                                                     available_InterpData    = available_InterpData )
!
                        ELSE
!
                           CALL Surface_Coordinate_Limits( subdomain        = Region(n),     &
     &                                                     table_1          = table_1,       &
     &                                                     table_2          = table_2,       &
     &                                                     XYZ_coordinates  = (/ X, Y, Z /), &
     &                                                     XYZ_limits       = XYZ_limits,    &
     &                                                     success          = correct_region )
!
                        END IF
!
                     END IF IF_tables
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
!
                     IF( .NOT. correct_region ) CYCLE DO_Cartesian
!
                     IF( X >= XYZ_limits(1,1) .AND. X <= XYZ_limits(1,2) .AND.   &
     &                   Y >= XYZ_limits(2,1) .AND. Y <= XYZ_limits(2,2) .AND.   &
     &                   Z >= XYZ_limits(3,1) .AND. Z <= XYZ_limits(3,2) )       &
     &               THEN
!
                        Media_Region = Region(n)%name
                        emissivity   = Region(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = MediaSequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( Media_Region(1:1) == 'F' .OR. Media_Region(1:1) == 'f' ) THEN
                              Media_Region = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              Media_Region = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
!                       RETURN
!
                     END IF
!
! ............... Cartesian grid system, unknown/unavailable shape of boundary
!
                  ELSE
!
! .................. Error
!
                     WRITE(*,6100) n, Region(n)%shape(1:4)
                     STOP
!
!
                  END IF IF_RShape
!
               END DO DO_Cartesian
!
!***********************************************************************
!*                                                                     *
!*          For cylindrical coordinate systems                         *
!*                                                                     *
!***********************************************************************
!
            ELSE IF(coordinates(1:4) == 'CYLI' .OR. coordinates(1:4) == 'cyli' .OR. coordinates(1:4) == 'Cyli') THEN
!
               DO_Cylindrical: DO n=2,TotNum_HetRegions
!
! ............... Cylindrical grid system, cylindrically-shaped or periodic media regions
!
                  IF_HRegions: IF( Region(n)%shape(1:4) == 'CYLI' .OR. Region(n)%shape(1:4) == 'cyli' .OR. Region(n)%shape(1:4) == 'Cyli' .OR. &
     &                             Region(n)%shape(1:2) == 'PE'   .OR. Region(n)%shape(1:2) == 'Pe'   .OR. Region(n)%shape(1:2) == 'pe' )      &
                  THEN
!
                     IF( Region(n)%CylBase1R /= 0.0d0 ) radius = MAX(Region(n)%CylBase1R, Region(n)%CylBase2R) - Z * omega(n)
!
                     RMax = MAX( Region(n)%CylBase1Coord(3), Region(n)%CylBase2Coord(3) )
                     RMin = MIN( Region(n)%CylBase1Coord(3), Region(n)%CylBase2Coord(3) )
!
                     IF( ( Region(n)%CylBase1R == 0.0d0 .AND. X >= Region(n)%CylRmin .AND. X <= Region(n)%CylRmax .AND. Z >= RMin .AND. Z <= RMax ) .OR.  &
     &                   ( Region(n)%CylBase1R /= 0.0d0 .AND. X <= radius .AND. Z >= RMin .AND. Z <= RMax ) )  &
     &               THEN
!
                        Media_Region = Region(n)%name
                        emissivity   = Region(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = MediaSequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( Media_Region(1:1) == 'F' .OR. Media_Region(1:1) == 'f' ) THEN
                              Media_Region = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              Media_Region = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, elliptical media regions
!
                  ELSE IF( Region(n)%shape(1:4) == 'ELLI' .OR. Region(n)%shape(1:4) == 'Elli' .OR. Region(n)%shape(1:4) == 'elli' ) THEN
!
! .................. The angle between the original angle-axis and the long axis of the ellipse (rotated system)
!
                     phi(n) = Region(n)%LAxisAngle
!
                     IF( Region(n)%Base1LAxis /= 0.0d0 ) THEN
                        Delta      = Z - Region(n)%EllBase1Coord(3)
                        long_axis  = Region(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                        short_axis = Region(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                     ELSE
                        long_axis  = Region(n)%LAxis
                        short_axis = Region(n)%SAxis
                     END IF
!
                     X_1p =  (X - Region(n)%EllBase1Coord(1)) * COS(phi(n)) + (Z - Region(n)%EllBase1Coord(3)) * SIN(phi(n))
                     Z_1p = -(X - Region(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Z - Region(n)%EllBase1Coord(3)) * COS(phi(n))
!
                     Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     IF( Delta >= 0.0d0 .AND. Delta <= height(n) .AND. Rtest <= 2.5d-1 ) THEN
!
                        Media_Region = Region(n)%name
                        emissivity   = Region(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = MediaSequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( Media_Region(1:1) == 'F' .OR. Media_Region(1:1) == 'f' ) THEN
                              Media_Region = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              Media_Region = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cylindrical grid system, spherically-shaped media regions
!
                  ELSE IF( (Region(n)%shape(1:4) == 'SPHE' .OR. Region(n)%shape(1:4) == 'sphe' .OR. Region(n)%shape(1:4) == 'Sphe') ) THEN
!
                     Z_1p = (Z - Region(n)%SphereCenterCoord(3))
!
                     Rtest = SQRT( X * X + Z_1p * Z_1p )
!
!
                     IF( Rtest >= Region(n)%SphRmin .AND. Rtest <= Region(n)%SphRmax ) THEN
!
                        Media_Region = Region(n)%name
                        emissivity   = Region(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = MediaSequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( Media_Region(1:1) == 'F' .OR. Media_Region(1:1) == 'f' ) THEN
                              Media_Region = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              Media_Region = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cylindrical grid system, irregularly-shaped exclusion regions
!
                  ELSE IF( (Region(n)%shape(1:4) == 'IRRE' .OR. Region(n)%shape(1:4) == 'irre' .OR. Region(n)%shape(1:4) == 'Irre') ) THEN
!
                     CALL Surface_Coordinate_Limits( subdomain        =  Region(n),    &
     &                                               RZ_coordinates   = (/ X, Z /),    &
     &                                               RZ_limits        = RZ_limits,     &
     &                                               success          = correct_region )
!
                     IF( .NOT. correct_region ) CYCLE DO_Cylindrical
!
                     IF( X >= RZ_limits(1,1) .AND. X <= RZ_limits(1,2) .AND. Z >= RZ_limits(2,1) .AND. Z <= RZ_limits(2,2) ) THEN
!
                        Media_Region = Region(n)%name
                        emissivity   = Region(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = MediaSequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( Media_Region(1:1) == 'F' .OR. Media_Region(1:1) == 'f' ) THEN
                              Media_Region = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              Media_Region = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
!                       RETURN
!
                     END IF
!
! ............... Cylindrical grid system, unknown/unavailable shape of exclusion region
!
                  ELSE
!
! .................. Error
!
                     WRITE(*,6102) n, Region(n)%shape
                     STOP
!
                  END IF IF_HRegions
!
               END DO DO_Cylindrical
!
!
!
            END IF IF_Assign
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
!6000 FORMAT(/,'Media_Region      1.0   30 August    2005',6X,'Assigning cells to regions of heterogeneous media')
!6000 FORMAT(/,'Media_Region      1.1    3 August    2014',6X,'Assigning cells to regions of heterogeneous media')
 6000 FORMAT(/,'Media_Region 2.0 ....................... 27 March     2015',6X,'Assigning cells to regions of heterogeneous media')
!
 6100 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system of the grid is "Cartesian", ', &
     &          'but the shape of the region #',I3,', is "',A,'": Unknown/Unavailable'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system of the grid is "cylindrical", ', &
     &          'but the shape of the region #',I3,', is "',A,'": Unknown/Unavailable'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6850 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Media_Region
!
!
            RETURN
!
         END FUNCTION Media_Region
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Boundary_Status( coordinates, X, Y, Z, bound_id, medium, media_by_number, emissivity, &
     &                               use_InterpData, completed_interpolation, available_InterpData )
!
            USE Subdomain_Bounding_Surfaces, ONLY: Surface_Coordinate_Limits
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      Routine determining whether a cell belongs to a boundary       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... User-defined variables
! -------------
!
            TYPE(Tabular) :: table_1, table_2
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: phi, theta, omega, omega1, omega2, height
!
            REAL(KIND = 8), OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables,3) :: available_InterpData
!
            REAL(KIND = 8), DIMENSION(3,2) :: XYZ_limits
            REAL(KIND = 8), DIMENSION(2,2) :: RZ_limits
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: X_1p, Y_1p, Z_1p, X_2p, Z_2p, Rtest, long_axis, short_axis, radius, Rmin, Rmax, Delta
!
            REAL(KIND = 8), INTENT(IN)  :: X, Y, Z
            REAL(KIND = 8), INTENT(OUT) :: emissivity
!
! -------------
! ......... Double precision parameters
! -------------
!
            REAL(KIND = 8), PARAMETER :: pi = 3.14159265358979324D0
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: n, N_0, N_1, N_2, N_3, N_4, N_5, m1, m2, ier
!
! -------------
! ......... CHARACTER variables
! -------------
!
            CHARACTER(LEN =  5) :: medium
!
            CHARACTER(LEN =  3), INTENT(OUT) :: bound_id
            CHARACTER(LEN = 11), INTENT(IN)  :: coordinates
!
! -------------
! ......... Logical variables
! -------------
!
            LOGICAL, INTENT(IN) :: media_by_number
!
            LOGICAL, OPTIONAL, INTENT(IN) :: use_InterpData
!
            LOGICAL, OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables) :: completed_interpolation
!
            LOGICAL :: First_call = .TRUE., correct_boundary, interp1, interp2
!
! -------------
! ......... Saving variables
! -------------
!
            SAVE First_call, phi, theta, omega, omega1, omega2, height
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Boundary_Status>
!
!
            IF(First_call) THEN
               WRITE(UNIT = *, FMT = 6000)
            ELSE
               GO TO 1000
            END IF
!
            ALLOCATE( phi(Num_Boundaries),   theta(Num_Boundaries),  height(Num_Boundaries), &
     &                omega(Num_Boundaries), omega1(Num_Boundaries), omega2(Num_Boundaries), STAT=ier )
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6850) 'phi, theta, height, omega, omega1, omega2', 'Boundary_Status'
               STOP
            END IF
!
! ......... Initialization
!
            phi    = 0.0d0
            theta  = 0.0d0
            height = 0.0d0
!
            omega  = 0.0d0
            omega1 = 0.0d0
            omega2 = 0.0d0
!
            emissivity = 0.0d0
!
!
!***********************************************************************
!*                                                                     *
!*          Baseline computations for cylindrical shapes in            *
!*          cartesian coordinate systems                               *
!*                                                                     *
!***********************************************************************
!
!
            IF_Initial: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                      coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
               DO_Boundaries: DO n=1,Num_Boundaries
!
! ............... Cartesian grid system, cylindrically-shaped media regions
!
                  IF_Shape: IF( (bound(n)%shape(1:2) == 'CY' .OR. bound(n)%shape(1:2) == 'cy' .OR. &
     &                           bound(n)%shape(1:2) == 'Cy') )  &
     &            THEN
!
                     CALL Determine_Orientation_Parameters( base1_center_coord = bound(n)%CylBase1Coord, &
     &                                                      base2_center_coord = bound(n)%CylBase2Coord, &
     &                                                      radius1 = bound(n)%CylBase1R, radius2 = bound(n)%CylBase2R, &
     &                                                      phi = phi(n), theta = theta(n), omega = omega(n), height = height(n) )
!
! ............... Cartesian grid system, elliptically-shaped media regions
!
                  ELSE IF( (bound(n)%shape(1:2) == 'EL' .OR. bound(n)%shape(1:2) == 'El' .OR. &
     &                      bound(n)%shape(1:2) == 'el') )  &
     &            THEN
!
                     CALL Determine_Orientation_Parameters( base1_center_coord = bound(n)%EllBase1Coord, &
     &                                                      base2_center_coord = bound(n)%EllBase2Coord, &
     &                                                      LAxis1 = bound(n)%Base1LAxis, SAxis1 = bound(n)%Base1SAxis, &
     &                                                      LAxis2 = bound(n)%Base2LAxis, SAxis2 = bound(n)%Base2SAxis, &
     &                                                      phi = phi(n), theta = theta(n), omega1 = omega1(n), omega2 = omega2(n), height = height(n) )
!
                  END IF IF_Shape
!
               END DO DO_Boundaries
!
            END IF IF_Initial
!
            First_call = .FALSE.
!
! -------------
! ......... Initialization
! -------------
!
 1000       bound_id  = '   '    ! Initial assignment to active status
!
!
!***********************************************************************
!*                                                                     *
!*          For cartesian coordinate systems                           *
!*                                                                     *
!***********************************************************************
!
!
            IF_Assign: IF(coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                    coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
               DO_Cartesian: DO n=1,Num_Boundaries
!
! ............... Cartesian grid system, Cartesian-coordinate-defined boundaries
!
                  IF_BShape: IF( bound(n)%shape(1:4) == 'RECT' .OR. bound(n)%shape(1:4) == 'Rect' .OR. bound(n)%shape(1:4) == 'rect') THEN
!
                     IF( X >= bound(n)%LMin(1) .AND. X <= bound(n)%LMax(1) .AND.   &
     &                   Y >= bound(n)%LMin(2) .AND. Y <= bound(n)%LMax(2) .AND.   &
     &                   Z >= bound(n)%LMin(3) .AND. Z <= bound(n)%LMax(3) )       &
     &               THEN
!
                        IF(bound(n)%name /= '*****') medium = bound(n)%name
                        emissivity = bound(n)%emissivity
                        bound_id   = bound(n)%id
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = Num_HetRegions + n, &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, cylindrically-shaped boundaries
!
                  ELSE IF( bound(n)%shape(1:4) == 'CYLI' .OR. bound(n)%shape(1:4) == 'cyli' .OR. &
     &                     bound(n)%shape(1:4) == 'Cyli' ) &
     &            THEN
!
! .................. The coordinates in the 1st rotated system
!
                     X_1p =  (X - bound(n)%CylBase1Coord(1)) * COS(phi(n)) + (Y - bound(n)%CylBase1Coord(2)) * SIN(phi(n))
                     Y_1p = -(X - bound(n)%CylBase1Coord(1)) * SIN(phi(n)) + (Y - bound(n)%CylBase1Coord(2)) * COS(phi(n))
                     Z_1p =   Z - bound(n)%CylBase1Coord(3)
!
! .................. The coordinates in the 1st rotated system
!
                     X_2p =  X_1p * COS(theta(n)) + Z_1p * SIN(theta(n))
                     Z_2p = -X_1p * SIN(theta(n)) + Z_1p * COS(theta(n))
!                    Y_2p =  Y_1p
!
                     IF( bound(n)%CylBase1R /= 0.0d0 ) radius = MAX(bound(n)%CylBase1R, bound(n)%CylBase2R) - Z_2p * omega(n)
!
                     Rtest = SQRT( Z_2p * Z_2p + Y_1p * Y_1p )
!
                     IF( ( bound(n)%CylBase1R == 0.0d0 .AND. X_2p >= 0.0d0 .AND. X_2p <= height(n) .AND. Rtest >= bound(n)%CylRmin .AND. Rtest <= bound(n)%CylRmax ) .OR. &
     &                   ( bound(n)%CylBase1R /= 0.0d0 .AND. X_2p >= 0.0d0 .AND. X_2p <= height(n) .AND. Rtest <= radius ) ) &
     &               THEN
!
                        IF(bound(n)%name /= '*****') medium = bound(n)%name
                        emissivity = bound(n)%emissivity
                        bound_id   = bound(n)%id
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = Num_HetRegions + n, &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, elliptical media regions
!
                  ELSE IF( bound(n)%shape(1:4) == 'ELLI' .OR. bound(n)%shape(1:4) == 'Elli' .OR. &
     &                     bound(n)%shape(1:4) == 'elli' ) &
     &            THEN
!
! .................. The angle between the original x- or y-axis and the long axis of the ellipse (rotated system)
!
                     phi(n) = bound(n)%LAxisAngle
!
                     SELECT CASE( bound(n)%EllipsePlane )
!
                     CASE( 'XY', 'xy', 'Xy', 'xY' )
!
                        IF( bound(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = Z - bound(n)%EllBase1Coord(3)
                           long_axis  = bound(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = bound(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = bound(n)%LAxis
                           short_axis = bound(n)%SAxis
                        END IF
!
                        X_1p =  (X - bound(n)%EllBase1Coord(1)) * COS(phi(n)) + (Y - bound(n)%EllBase1Coord(2)) * SIN(phi(n))
                        Y_1p = -(X - bound(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Y - bound(n)%EllBase1Coord(2)) * COS(phi(n))
!
                        Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Y_1p * Y_1p / ( short_axis * short_axis )
!
                     CASE( 'XZ', 'xz', 'Xz', 'xZ' )
!
                        IF( bound(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = Y - bound(n)%EllBase1Coord(2)
                           long_axis  = bound(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = bound(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = bound(n)%LAxis
                           short_axis = bound(n)%SAxis
                        END IF
!
                        X_1p =  (X - bound(n)%EllBase1Coord(1)) * COS(phi(n)) + (Z - bound(n)%EllBase1Coord(3)) * SIN(phi(n))
                        Z_1p = -(X - bound(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Z - bound(n)%EllBase1Coord(3)) * COS(phi(n))
!
                        Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     CASE( 'YZ', 'yz', 'Yz', 'yZ' )
!
                        IF( bound(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = X - bound(n)%EllBase1Coord(1)
                           long_axis  = bound(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = bound(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = bound(n)%LAxis
                           short_axis = bound(n)%SAxis
                        END IF
!
                        Y_1p =  (Y - bound(n)%EllBase1Coord(2)) * COS(phi(n)) + (Z - bound(n)%EllBase1Coord(3)) * SIN(phi(n))
                        Z_1p = -(Y - bound(n)%EllBase1Coord(2)) * SIN(phi(n)) + (Z - bound(n)%EllBase1Coord(3)) * COS(phi(n))
!
                        Rtest = Y_1p * Y_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     END SELECT
!
                     IF( Delta >= 0.0d0 .AND. Delta <= height(n) .AND. Rtest <= 2.5d-1 ) THEN
!
                        IF(bound(n)%name /= '*****') medium = bound(n)%name
                        emissivity = bound(n)%emissivity
                        bound_id   = bound(n)%id
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = Num_HetRegions + n, &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, spherically-shaped boundaries
!
                  ELSE IF( bound(n)%shape(1:4) == 'SPHE' .OR. bound(n)%shape(1:4) == 'sphe' .OR. bound(n)%shape(1:4) == 'Sphe') THEN
!
                     X_1p = (X - bound(n)%SphereCenterCoord(1))
                     Y_1p = (Y - bound(n)%SphereCenterCoord(2))
                     Z_1p = (Z - bound(n)%SphereCenterCoord(3))
!
                     Rtest = SQRT( X_1p * X_1p + Y_1p * Y_1p + Z_1p * Z_1p )
!
                     IF( Rtest >= bound(n)%SphRmin .AND. Rtest <= bound(n)%SphRmax ) THEN
!
                        IF(bound(n)%name /= '*****') medium = bound(n)%name
                        emissivity = bound(n)%emissivity
                        bound_id   = bound(n)%id
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = Num_HetRegions + n, &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, irregularly-shaped boundaries
!
                  ELSE IF( (bound(n)%shape(1:4) == 'IRRE' .OR. bound(n)%shape(1:4) == 'irre' .OR. bound(n)%shape(1:4) == 'Irre') ) THEN
!
                     m1      = bound(n)%IntTableNum1
                     interp1 = .FALSE.
                     IF( m1 > 0 ) THEN
                        table_1 = IntTable(m1)
                        IF( completed_interpolation(m1) ) interp1 = .TRUE.
                     END IF
!
                     m2      = bound(n)%IntTableNum2
                     interp2 = .FALSE.
                     IF( m2 > 0 ) THEN
                        table_2 = IntTable(m2)
                        IF( completed_interpolation(m2) ) interp2 = .TRUE.
                     END IF
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
!
                     IF_tables: IF( m1 == 0 .AND. m2 == 0 ) THEN
!
                        CALL Surface_Coordinate_Limits( subdomain        = bound(n),        &
     &                                                  XYZ_coordinates  = (/ X, Y, Z /),   &
     &                                                  XYZ_limits       = XYZ_limits,      &
     &                                                  success          = correct_boundary )
!
                     ELSE
!
                        IF( PRESENT(use_InterpData) .AND. (interp1 .OR. interp2) ) THEN
!
                           CALL Surface_Coordinate_Limits( subdomain        = bound(n),         &
     &                                                     table_1          = table_1,          &
     &                                                     table_2          = table_2,          &
     &                                                     XYZ_coordinates  = (/ X, Y, Z /),    &
     &                                                     XYZ_limits       = XYZ_limits,       &
     &                                                     success          = correct_boundary, &
     &                                                     use_InterpData   = use_InterpData,   &
     &                                                     completed_interpolation = completed_interpolation, &
     &                                                     available_InterpData    = available_InterpData )
!
                        ELSE
!
                           CALL Surface_Coordinate_Limits( subdomain        = bound(n),        &
     &                                                     table_1          = table_1,         &
     &                                                     table_2          = table_2,         &
     &                                                     XYZ_coordinates  = (/ X, Y, Z /),   &
     &                                                     XYZ_limits       = XYZ_limits,      &
     &                                                     success          = correct_boundary )
!
                        END IF
!
                     END IF IF_tables
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
!
                     IF( .NOT. correct_boundary ) CYCLE DO_Cartesian
!
                     IF( X >= XYZ_limits(1,1) .AND. X <= XYZ_limits(1,2) .AND.   &
     &                   Y >= XYZ_limits(2,1) .AND. Y <= XYZ_limits(2,2) .AND.   &
     &                   Z >= XYZ_limits(3,1) .AND. Z <= XYZ_limits(3,2) )       &
     &               THEN
!
                        IF(bound(n)%name /= '*****') medium = bound(n)%name
                        emissivity = bound(n)%emissivity
                        bound_id   = bound(n)%id
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = Num_HetRegions + n, &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, unknown/unavailable shape of boundary
!
                  ELSE
!
! .................. Error
!
                     WRITE(*,6100) n, bound(n)%shape(1:4)
                     STOP
!
                  END IF IF_BShape
!
               END DO DO_Cartesian
!
!***********************************************************************
!*                                                                     *
!*          For cylindrical coordinate systems                         *
!*                                                                     *
!***********************************************************************
!
            ELSE IF(coordinates(1:4) == 'CYLI' .OR. coordinates(1:4) == 'cyli' .OR. coordinates(1:4) == 'Cyli') THEN
!
               DO_Cylindrical: DO n=1,Num_Boundaries
!
! ............... Cylindrical grid system, cylindrical-coordinate-defined media boundaries
!
                  IF_BRegions: IF( bound(n)%shape(1:4) == 'CYLI' .OR. bound(n)%shape(1:4) == 'cyli' .OR. bound(n)%shape(1:4) == 'Cyli' ) THEN
!
                     IF( bound(n)%CylBase1R /= 0.0d0 ) radius = MAX(bound(n)%CylBase1R, bound(n)%CylBase2R) - Z * omega(n)
!
                     RMax = MAX( bound(n)%CylBase1Coord(3), bound(n)%CylBase2Coord(3) )
                     RMin = MIN( bound(n)%CylBase1Coord(3), bound(n)%CylBase2Coord(3) )
!
                     IF( ( bound(n)%CylBase1R == 0.0d0 .AND. X >= bound(n)%CylRmin .AND. X <= bound(n)%CylRmax .AND. Z >= RMin .AND. Z <= RMax ) .OR.  &
     &                   ( bound(n)%CylBase1R /= 0.0d0 .AND. X <= radius .AND. Z >= RMin .AND. Z <= RMax ) )         &
     &               THEN
!
                        IF(bound(n)%name /= '*****') medium = bound(n)%name
                        emissivity = bound(n)%emissivity
                        bound_id   = bound(n)%id
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = Num_HetRegions + n, &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cartesian grid system, elliptical boundaries
!
                  ELSE IF( bound(n)%shape(1:4) == 'ELLI' .OR. bound(n)%shape(1:4) == 'Elli' .OR. &
     &                     bound(n)%shape(1:4) == 'elli' ) &
     &            THEN
!
! .................. The angle between the original angle-axis and the long axis of the ellipse (rotated system)
!
                     phi(n) = bound(n)%LAxisAngle
!
                     IF( bound(n)%Base1LAxis /= 0.0d0 ) THEN
                        Delta      = Z - bound(n)%EllBase1Coord(3)
                        long_axis  = bound(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                        short_axis = bound(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                     ELSE
                        long_axis  = bound(n)%LAxis
                        short_axis = bound(n)%SAxis
                     END IF
!
                     X_1p =  (X - bound(n)%EllBase1Coord(1)) * COS(phi(n)) + (Z - bound(n)%EllBase1Coord(3)) * SIN(phi(n))
                     Z_1p = -(X - bound(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Z - bound(n)%EllBase1Coord(3)) * COS(phi(n))
!
                     Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     IF( Delta >= 0.0d0 .AND. Delta <= height(n) .AND. Rtest <= 2.5d-1 ) THEN
!
                        IF(bound(n)%name /= '*****') medium = bound(n)%name
                        emissivity = bound(n)%emissivity
                        bound_id   = bound(n)%id
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = Num_HetRegions + n, &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cylindrical grid system, spherically-shaped media boundaries
!
                  ELSE IF( bound(n)%shape(1:4) == 'SPHE' .OR. bound(n)%shape(1:4) == 'sphe' .OR. &
     &                     bound(n)%shape(1:4) == 'Sphe' ) &
     &            THEN
!
                     X_1p = (X - bound(n)%SphereCenterCoord(1))
                     Z_1p = (Z - bound(n)%SphereCenterCoord(3))
!
                     Rtest = SQRT( X_1p * X_1p + Z_1p * Z_1p )
!
                     IF( Rtest >= bound(n)%SphRmin .AND. Rtest <= bound(n)%SphRmax ) THEN
!
                        IF(bound(n)%name /= '*****') medium = bound(n)%name
                        emissivity = bound(n)%emissivity
                        bound_id   = bound(n)%id
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = Num_HetRegions + n, &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cylindrical grid system, irregularly-shaped boundaries
!
                  ELSE IF( bound(n)%shape(1:4) == 'IRRE' .OR. bound(n)%shape(1:4) == 'irre' .OR. &
     &                     bound(n)%shape(1:4) == 'Irre' ) &
     &            THEN
!
                     CALL Surface_Coordinate_Limits( subdomain        = bound(n),        &
     &                                               RZ_coordinates   = (/ X, Z /),      &
     &                                               RZ_limits        = RZ_limits,       &
     &                                               success          = correct_boundary )

!
                     IF( .NOT. correct_boundary ) CYCLE DO_Cylindrical
!
                     IF( X >= RZ_limits(1,1) .AND. X <= RZ_limits(1,2) .AND. Z >= RZ_limits(2,1) .AND. Z <= RZ_limits(2,2) ) THEN
!
                        IF(bound(n)%name /= '*****') medium = bound(n)%name
                        emissivity = bound(n)%emissivity
                        bound_id   = bound(n)%id
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = Num_HetRegions + n, &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                     END IF
!
! ............... Cylindrical grid system, unknown/unavailable shape of boundaries
!
                  ELSE
!
! .................. Error
!
                     WRITE(*,6102) n, bound(n)%shape(1:4)
                     STOP
!
                  END IF IF_BRegions
!
               END DO DO_Cylindrical
!
!
!
            END IF IF_Assign
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
!6000 FORMAT(/,'Boundary_Status   1.0   31 August    2005',6X,'Determinimg whether a cell belongs to a boundary')
!6000 FORMAT(/,'Boundary_Status   1.1    3 August    2014',6X,'Determinimg whether a cell belongs to a boundary')
 6000 FORMAT(/,'Boundary_Status 2.0 .................... 29 March     2015',6X,'Determinimg whether a cell belongs to a boundary')
!
 6100 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system of the grid is "Cartesian", ', &
     &          'but the shape of the boundary region #',I3,', is "',A,'": Unknown/Unavailable'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system of the grid is "cylindrical", ', &
     &          'but the shape of the boundary region #',I3,', is "',A,'": Unknown/Unavailable'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6850 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Boundary_Status
!
!
            RETURN
!
         END SUBROUTINE Boundary_Status
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         LOGICAL FUNCTION Excluded_Element( coordinates, X, Y, Z, exclusion_zone_type, use_InterpData, &
     &                                      completed_interpolation, available_InterpData, find_RefSurface )
!
            USE Subdomain_Bounding_Surfaces, ONLY: Surface_Coordinate_Limits
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*             Routine determining the media region where              *
!*                     a particular cell belongs                       *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... User-defined variables
! -------------
!
            TYPE(Tabular) :: table_1, table_2
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: phi, theta, omega, omega1, omega2, height
!
            REAL(KIND = 8), OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables,3) :: available_InterpData
!
            REAL(KIND = 8), DIMENSION(3,2) :: XYZ_limits
            REAL(KIND = 8), DIMENSION(2,2) :: RZ_limits
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: X_1p, Y_1p, Z_1p, X_2p, Z_2p, Rtest, radius, long_axis, short_axis, Rmax, Rmin, Delta, ZLimit_min, ZLimit_max
!
            REAL(KIND = 8), INTENT(IN) :: X, Y, Z
!
! -------------
! ......... Double precision parameters
! -------------
!
            REAL(KIND = 8), PARAMETER :: pi = 3.14159265358979324D0
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: n, m1, m2, ier
!
! -------------
! ......... CHARACTER variables
! -------------
!
            CHARACTER(LEN =  3), INTENT(OUT) :: exclusion_zone_type
            CHARACTER(LEN = 11), INTENT(IN)  :: coordinates
!
! -------------
! ......... Logical variables
! -------------
!
            LOGICAL, OPTIONAL, INTENT(IN) :: use_InterpData, find_RefSurface
!
            LOGICAL, OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables) :: completed_interpolation
!
            LOGICAL :: First_call = .TRUE., correct_exclusion, interp1, interp2
!
! -------------
! ......... Saving variables
! -------------
!
            SAVE First_call, phi, theta, omega, omega1, omega2, height
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of <Excluded_Element>
!
!
            IF(First_call) THEN
               WRITE(UNIT = *, FMT = 6000)
            ELSE
               GO TO 1000
            END IF
!
            ALLOCATE( phi(Num_ExclZones),   theta(Num_ExclZones),  height(Num_ExclZones), &
     &                omega(Num_ExclZones), omega1(Num_ExclZones), omega2(Num_ExclZones), STAT=ier )
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6850) 'phi, theta, height, omega, omega1, omega2', 'Excluded_Element'
               STOP
            END IF
!
! ......... Initialization
!
            phi    = 0.0d0
            theta  = 0.0d0
            height = 0.0d0
!
            omega  = 0.0d0
            omega1 = 0.0d0
            omega2 = 0.0d0
!
!
!***********************************************************************
!*                                                                     *
!*          Baseline computations for cylindrical shapes in            *
!*          cartesian coordinate systems                               *
!*                                                                     *
!***********************************************************************
!
!
            IF_Initial: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                      coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
               DO_Exclusions: DO n=1,Num_ExclZones
!
! ............... Cartesian grid system, cylindrically-shaped media regions
!
                  IF_Shape: IF( (ExclZone(n)%shape(1:2) == 'CY' .OR. ExclZone(n)%shape(1:2) == 'cy' .OR. &
     &                           ExclZone(n)%shape(1:2) == 'Cy') )  &
     &            THEN
!
                     CALL Determine_Orientation_Parameters( base1_center_coord = ExclZone(n)%CylBase1Coord, &
     &                                                      base2_center_coord = ExclZone(n)%CylBase2Coord, &
     &                                                      radius1 = ExclZone(n)%CylBase1R, radius2 = ExclZone(n)%CylBase2R, &
     &                                                      phi = phi(n), theta = theta(n), omega = omega(n), height = height(n) )
!
! ............... Cartesian grid system, elliptically-shaped media regions
!
                  ELSE IF( (ExclZone(n)%shape(1:2) == 'EL' .OR. ExclZone(n)%shape(1:2) == 'El' .OR. &
     &                      ExclZone(n)%shape(1:2) == 'el') )  &
     &            THEN
!
                     CALL Determine_Orientation_Parameters( base1_center_coord = ExclZone(n)%EllBase1Coord, &
     &                                                      base2_center_coord = ExclZone(n)%EllBase2Coord, &
     &                                                      LAxis1 = ExclZone(n)%Base1LAxis, SAxis1 = ExclZone(n)%Base1SAxis, &
     &                                                      LAxis2 = ExclZone(n)%Base2LAxis, SAxis2 = ExclZone(n)%Base2SAxis, &
     &                                                      phi = phi(n), theta = theta(n), omega1 = omega1(n), omega2 = omega2(n), height = height(n) )
!
                  END IF IF_Shape
!
               END DO DO_Exclusions
!
            END IF IF_Initial
!
            First_call = .FALSE.
!
! -------------
! ......... Initialization
! -------------
!
 1000       Excluded_Element    = .FALSE.  ! Initial assignments
            exclusion_zone_type = '   '
!
!
!***********************************************************************
!*                                                                     *
!*          For cartesian coordinate systems                           *
!*                                                                     *
!***********************************************************************
!
!
            IF_Assign: IF(coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                    coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
               DO_Cartesian: DO n=1,Num_ExclZones
!
! ............... For determination of reference surface in thin areal grids, ensure that the dependent variable is NOT 'Z'
!
                  IF_RefSurf: IF( PRESENT(find_RefSurface) ) THEN
                     IF( find_RefSurface ) THEN
                       IF( ExclZone(n)%DependentVar == 'Z' .OR. ExclZone(n)%DependentVar == 'z' ) CYCLE DO_Cartesian
                       IF( ExclZone(n)%shape(1:1) /= 'I' .AND. ExclZone(n)%shape(1:1) /= 'i' ) CYCLE DO_Cartesian
                     END IF
                  END IF IF_RefSurf
!
! ............... Cartesian grid system, Cartesian-coordinate-defined exclusion zones
!
                  IF_ExclShape: IF( ExclZone(n)%shape(1:4) == 'RECT' .OR. ExclZone(n)%shape(1:4) == 'rect' .OR. ExclZone(n)%shape(1:4) == 'Rect') THEN
!
                     IF( X >= ExclZone(n)%LMin(1) .AND. X <= ExclZone(n)%LMax(1) .AND.   &
     &                   Y >= ExclZone(n)%LMin(2) .AND. Y <= ExclZone(n)%LMax(2) .AND.   &
     &                   Z >= ExclZone(n)%LMin(3) .AND. Z <= ExclZone(n)%LMax(3) ) &
     &               THEN
                        Excluded_Element    = .TRUE.
                        exclusion_zone_type = ExclZone(n)%id
                        RETURN
                     END IF
!
! ............... Cartesian grid system, cylindrically-shaped exclusion zones
!
                  ELSE IF( (ExclZone(n)%shape(1:4) == 'CYLI' .OR. ExclZone(n)%shape(1:4) == 'cyli' .OR. ExclZone(n)%shape(1:4) == 'Cyli') ) THEN
!
! .................. The coordinates in the 1st rotated system
!
                     X_1p =  (X - ExclZone(n)%CylBase1Coord(1)) * COS(phi(n)) + (Y - ExclZone(n)%CylBase1Coord(2)) * SIN(phi(n))
                     Y_1p = -(X - ExclZone(n)%CylBase1Coord(1)) * SIN(phi(n)) + (Y - ExclZone(n)%CylBase1Coord(2)) * COS(phi(n))
                     Z_1p =   Z - ExclZone(n)%CylBase1Coord(3)
!
! .................. The coordinates in the 2nd rotated system
!
                     X_2p =  X_1p * COS(theta(n)) + Z_1p * SIN(theta(n))
                     Z_2p = -X_1p * SIN(theta(n)) + Z_1p * COS(theta(n))
!                    Y_2p =  Y_1p
!
                     IF( ExclZone(n)%CylBase1R /= 0.0d0 ) radius = MAX(ExclZone(n)%CylBase1R, ExclZone(n)%CylBase2R) - X_2p * omega(n)
!
                     Rtest = SQRT( Z_2p * Z_2p + Y_1p * Y_1p )
!
                     IF( ( ExclZone(n)%CylBase1R == 0.0d0 .AND. X_2p >= 0.0d0 .AND. X_2p <= height(n) .AND. Rtest >= ExclZone(n)%CylRmin .AND. Rtest <= ExclZone(n)%CylRmax ) .OR. &
     &                   ( ExclZone(n)%CylBase1R /= 0.0d0 .AND. X_2p >= 0.0d0 .AND. X_2p <= height(n) .AND. Rtest <= radius ) ) &
     &               THEN
                        Excluded_Element    = .TRUE.
                        exclusion_zone_type = ExclZone(n)%id
                        RETURN
                     END IF
!
! ............... Cartesian grid system, elliptical exclusion zones
!
                  ELSE IF( ExclZone(n)%shape(1:4) == 'ELLI' .OR. ExclZone(n)%shape(1:4) == 'Elli' .OR. ExclZone(n)%shape(1:4) == 'elli' ) THEN
!
! .................. The angle between the original x- or y-axis and the long axis of the ellipse (rotated system)
!
                     phi(n) = ExclZone(n)%LAxisAngle
!
                     SELECT CASE( ExclZone(n)%EllipsePlane )
!
                     CASE( 'XY', 'xy', 'Xy', 'xY' )
!
                        IF( ExclZone(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = Z - ExclZone(n)%EllBase1Coord(3)
                           long_axis  = ExclZone(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = ExclZone(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = ExclZone(n)%LAxis
                           short_axis = ExclZone(n)%SAxis
                        END IF
!
                        X_1p =  (X - ExclZone(n)%EllBase1Coord(1)) * COS(phi(n)) + (Y - ExclZone(n)%EllBase1Coord(2)) * SIN(phi(n))
                        Y_1p = -(X - ExclZone(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Y - ExclZone(n)%EllBase1Coord(2)) * COS(phi(n))
!
                        Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Y_1p * Y_1p / ( short_axis * short_axis )
!
                     CASE( 'XZ', 'xz', 'Xz', 'xZ' )
!
                        IF( ExclZone(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = Y - ExclZone(n)%EllBase1Coord(2)
                           long_axis  = ExclZone(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = ExclZone(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = ExclZone(n)%LAxis
                           short_axis = ExclZone(n)%SAxis
                        END IF
!
                        X_1p =  (X - ExclZone(n)%EllBase1Coord(1)) * COS(phi(n)) + (Z - ExclZone(n)%EllBase1Coord(3)) * SIN(phi(n))
                        Z_1p = -(X - ExclZone(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Z - ExclZone(n)%EllBase1Coord(3)) * COS(phi(n))
!
                        Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     CASE( 'YZ', 'yz', 'Yz', 'yZ' )
!
                        IF( ExclZone(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = X - ExclZone(n)%EllBase1Coord(1)
                           long_axis  = ExclZone(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = ExclZone(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = ExclZone(n)%LAxis
                           short_axis = ExclZone(n)%SAxis
                        END IF
!
                        Y_1p =  (Y - ExclZone(n)%EllBase1Coord(2)) * COS(phi(n)) + (Z - ExclZone(n)%EllBase1Coord(3)) * SIN(phi(n))
                        Z_1p = -(Y - ExclZone(n)%EllBase1Coord(2)) * SIN(phi(n)) + (Z - ExclZone(n)%EllBase1Coord(3)) * COS(phi(n))
!
                        Rtest = Y_1p * Y_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     END SELECT
!
                     IF( Delta >= 0.0d0 .AND. Delta <= height(n) .AND. Rtest <= 2.5d-1 ) THEN
                        Excluded_Element    = .TRUE.
                        exclusion_zone_type = ExclZone(n)%id
                        RETURN
                     END IF
!
! ............... Cartesian grid system, spherically-shaped exclusion zones
!
                  ELSE IF( (ExclZone(n)%shape(1:4) == 'SPHE' .OR. ExclZone(n)%shape(1:4) == 'sphe' .OR. ExclZone(n)%shape(1:4) == 'Sphe') ) THEN
!
                     X_1p = (X - ExclZone(n)%SphereCenterCoord(1))
                     Y_1p = (Y - ExclZone(n)%SphereCenterCoord(2))
                     Z_1p = (Z - ExclZone(n)%SphereCenterCoord(3))
!
                     Rtest = SQRT( X_1p * X_1p + Y_1p * Y_1p + Z_1p * Z_1p )
!
                     IF( Rtest >= ExclZone(n)%SphRmin .AND. Rtest <= ExclZone(n)%SphRmax ) THEN
                        Excluded_Element    = .TRUE.
                        exclusion_zone_type = ExclZone(n)%id
                        RETURN
                     END IF
!
! ............... Cartesian grid system, irregularly-shaped exclusion zones
!
                  ELSE IF( (ExclZone(n)%shape(1:4) == 'IRRE' .OR. ExclZone(n)%shape(1:4) == 'irre' .OR. ExclZone(n)%shape(1:4) == 'Irre') ) THEN
!
                     m1      = ExclZone(n)%IntTableNum1
                     interp1 = .FALSE.
                     IF( m1 > 0 ) THEN
                        table_1 = IntTable(m1)
                        IF( completed_interpolation(m1) ) interp1 = .TRUE.
                     END IF
!
                     m2      = ExclZone(n)%IntTableNum2
                     interp2 = .FALSE.
                     IF( m2 > 0 ) THEN
                        table_2 = IntTable(m2)
                        IF( completed_interpolation(m2) ) interp2 = .TRUE.
                     END IF
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
!
                     IF_tables: IF( m1 == 0 .AND. m2 == 0 ) THEN
!
                        CALL Surface_Coordinate_Limits( subdomain        = ExclZone(n),      &
     &                                                  XYZ_coordinates  = (/ X, Y, Z /),    &
     &                                                  XYZ_limits       = XYZ_limits,       &
     &                                                  success          = correct_exclusion )
!
                     ELSE
!
                        IF( (PRESENT(use_InterpData)) .AND. (interp1 .OR. interp2) ) THEN
!
                           CALL Surface_Coordinate_Limits( subdomain        = ExclZone(n),       &
     &                                                     table_1          = table_1,           &
     &                                                     table_2          = table_2,           &
     &                                                     XYZ_coordinates  = (/ X, Y, Z /),     &
     &                                                     XYZ_limits       = XYZ_limits,        &
     &                                                     success          = correct_exclusion, &
     &                                                     use_InterpData   = use_InterpData,    &
     &                                                     completed_interpolation = completed_interpolation, &
     &                                                     available_InterpData    = available_InterpData )
!
                        ELSE
!
                           CALL Surface_Coordinate_Limits( subdomain        = ExclZone(n),      &
     &                                                     table_1          = table_1,          &
     &                                                     table_2          = table_2,          &
     &                                                     XYZ_coordinates  = (/ X, Y, Z /),    &
     &                                                     XYZ_limits       = XYZ_limits,       &
     &                                                     success          = correct_exclusion )
!
                        END IF
!
                     END IF IF_tables
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
!
  200                IF( .NOT. correct_exclusion ) CYCLE DO_Cartesian
!
! .................. For determination of reference surface in thin areal grids, ensure that the dependent variable is NOT 'Z'
!
                     ZLimit_min = XYZ_limits(3,1)
                     ZLimit_max = XYZ_limits(3,2)
!
                     IF( PRESENT(find_RefSurface) ) THEN
                        IF( find_RefSurface ) THEN
                           ZLimit_min = -1.0d9
                           ZLimit_max =  1.0d9
                        END IF
                     END IF
!
                     IF( X >= XYZ_limits(1,1) .AND. X <= XYZ_limits(1,2) .AND.   &
     &                   Y >= XYZ_limits(2,1) .AND. Y <= XYZ_limits(2,2) .AND.   &
     &                   Z >= ZLimit_min      .AND. Z <= ZLimit_max )            &
     &               THEN
                        Excluded_Element    = .TRUE.
                        exclusion_zone_type = ExclZone(n)%id
                        RETURN
                     END IF
!
! ............... Cartesian grid system, unknown/unavailable shape of exclusion region
!
                  ELSE
!
! .................. Error
!
                     WRITE(*,6100) n
                     STOP
!
                  END IF IF_ExclShape
!
               END DO DO_Cartesian
!
!***********************************************************************
!*                                                                     *
!*          For cylindrical coordinate systems                         *
!*                                                                     *
!***********************************************************************
!
            ELSE IF(coordinates(1:4) == 'CYLI' .OR. coordinates(1:4) == 'cyli' .OR. coordinates(1:4) == 'Cyli') THEN
!
               DO_Cylindrical: DO n=1,Num_ExclZones
!
! ............... Cylindrical grid system, cylindrically-shaped media exclusion zones
!
                  IF_ExclShape2: IF( ExclZone(n)%shape(1:4) == 'CYLI' .OR. ExclZone(n)%shape(1:4) == 'cyli' .OR. ExclZone(n)%shape(1:4) == 'Cyli') THEN
!
                     IF( ExclZone(n)%CylBase1R /= 0.0d0 ) radius = MAX(ExclZone(n)%CylBase1R, ExclZone(n)%CylBase2R) - Z * omega(n)
!
                     RMax = MAX( ExclZone(n)%CylBase1Coord(3), ExclZone(n)%CylBase2Coord(3) )
                     RMin = MIN( ExclZone(n)%CylBase1Coord(3), ExclZone(n)%CylBase2Coord(3) )
!
                     IF( ( ExclZone(n)%CylBase1R == 0.0d0 .AND. X >= ExclZone(n)%CylRmin .AND. X <= ExclZone(n)%CylRmax .AND. Z >= RMin .AND. Z <= RMax ) .OR.    &
     &                   ( ExclZone(n)%CylBase1R /= 0.0d0 .AND. X <= radius .AND. Z >= RMin .AND. Z <= RMax ) )  &
     &               THEN
                        Excluded_Element    = .TRUE.
                        exclusion_zone_type = ExclZone(n)%id
                        RETURN
                     END IF
!
! ............... Cartesian grid system, elliptical exclusion zones
!
                  ELSE IF( ExclZone(n)%shape(1:4) == 'ELLI' .OR. ExclZone(n)%shape(1:4) == 'Elli' .OR. ExclZone(n)%shape(1:4) == 'elli' ) THEN
!
! .................. The angle between the original angle-axis and the long axis of the ellipse (rotated system)
!
                     phi(n) = ExclZone(n)%LAxisAngle
!
                     IF( ExclZone(n)%Base1LAxis /= 0.0d0 ) THEN
                        Delta      = Z - ExclZone(n)%EllBase1Coord(3)
                        long_axis  = ExclZone(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                        short_axis = ExclZone(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                     ELSE
                        long_axis  = ExclZone(n)%LAxis
                        short_axis = ExclZone(n)%SAxis
                     END IF
!
                     X_1p =  (X - ExclZone(n)%EllBase1Coord(1)) * COS(phi(n)) + (Z - ExclZone(n)%EllBase1Coord(3)) * SIN(phi(n))
                     Z_1p = -(X - ExclZone(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Z - ExclZone(n)%EllBase1Coord(3)) * COS(phi(n))
!
                     Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     IF( Delta >= 0.0d0 .AND. Delta <= height(n) .AND. Rtest <= 2.5d-1 ) THEN
                        Excluded_Element    = .TRUE.
                        exclusion_zone_type = ExclZone(n)%id
                        RETURN
                     END IF
!
! ............... Cylindrical grid system, spherically-shaped exclusion zones
!
                  ELSE IF( (ExclZone(n)%shape(1:4) == 'SPHE' .OR. ExclZone(n)%shape(1:4) == 'sphe' .OR. ExclZone(n)%shape(1:4) == 'Sphe') ) THEN
!
                     X_1p = (X - ExclZone(n)%SphereCenterCoord(1))
                     Z_1p = (Z - ExclZone(n)%SphereCenterCoord(3))
!
                     Rtest = SQRT( X_1p * X_1p + Z_1p * Z_1p )
!
                     IF( Rtest >= ExclZone(n)%SphRmin .AND. Rtest <= ExclZone(n)%SphRmax ) THEN
                        Excluded_Element    = .TRUE.
                        exclusion_zone_type = ExclZone(n)%id
                        RETURN
                     END IF
!
! ............... Cylindrical grid system, irregularly-shaped exclusion zones
!
                  ELSE IF( (ExclZone(n)%shape(1:4) == 'IRRE' .OR. ExclZone(n)%shape(1:4) == 'irre' .OR. ExclZone(n)%shape(1:4) == 'Irre') ) THEN
!
                     CALL Surface_Coordinate_Limits( subdomain        = ExclZone(n),      &
     &                                               RZ_coordinates   = (/ X, Z /),       &
     &                                               RZ_limits        = RZ_limits ,       &
     &                                               success          = correct_exclusion )

!
                     IF( .NOT. correct_exclusion ) CYCLE DO_Cylindrical
!
                     IF( X >= RZ_limits(1,1) .AND. X <= RZ_limits(1,2) .AND. Z >= RZ_limits(2,1) .AND. Z <= RZ_limits(2,2) ) THEN
                        Excluded_Element    = .TRUE.
                        exclusion_zone_type = ExclZone(n)%id
                        RETURN
                     END IF
!
! ............... Cylindrical grid system, unknown/unavailable shape of exclusion zones
!
                  ELSE
!
! .................. Error
!
                     WRITE(*,6102) n
                     STOP
!
                  END IF IF_ExclShape2
!
               END DO DO_Cylindrical
!
!
!
            END IF IF_Assign
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Excluded_Element 1.0 ................... 29 March     2015',6X,'Checking exclusion regions and removing elements from consideration')
!
 6100 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system of the grid is "Cartesian", ', &
     &          'but the shape of the exclusion region #',I3,', is "',A,'": Unknown/Unavailable'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system of the grid is "cylindrical", ', &
     &          'but the shape of the exclusion region #',I3,', is "',A,'": Unknown/Unavailable'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6850 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Excluded_Element>
!
!
            RETURN
!
         END FUNCTION Excluded_Element
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         LOGICAL FUNCTION Included_Element( coordinates, X, Y, Z, inclusion_id, medium, media_by_number, emissivity, &
     &                                      use_InterpData, completed_interpolation, available_InterpData )
!
            USE Subdomain_Bounding_Surfaces, ONLY: Surface_Coordinate_Limits
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*   Routine determining whether a cell belongs to an inclusion zone   *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... User-defined variables
! -------------
!
            TYPE(Tabular) :: table_1, table_2
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: phi, theta, omega, omega1, omega2, height
!
            REAL(KIND = 8), OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables,3) :: available_InterpData
!
            REAL(KIND = 8), DIMENSION(3,2) :: XYZ_limits
            REAL(KIND = 8), DIMENSION(2,2) :: RZ_limits
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: X_1p, Y_1p, Z_1p, X_2p, Z_2p, Rtest, long_axis, short_axis, radius, Rmin, Rmax, Delta
!
            REAL(KIND = 8), INTENT(IN)  :: X, Y, Z
            REAL(KIND = 8), INTENT(OUT) :: emissivity
!
! -------------
! ......... Double precision parameters
! -------------
!
            REAL(KIND = 8), PARAMETER :: pi = 3.14159265358979324D0
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: n, N_0, N_1, N_2, N_3, N_4, N_5, mbegin, m1, m2, ier
!
! -------------
! ......... CHARACTER variables
! -------------
!
            CHARACTER(LEN =  5) :: medium
!
            CHARACTER(LEN = 11), INTENT(IN)  :: coordinates
            CHARACTER(LEN =  3), INTENT(OUT) :: inclusion_id
!
! -------------
! ......... Logical variables
! -------------
!
            LOGICAL, INTENT(IN) :: media_by_number
!
            LOGICAL, OPTIONAL, INTENT(IN) :: use_InterpData
!
            LOGICAL, OPTIONAL, INTENT(INOUT), DIMENSION(Num_IntTables) :: completed_interpolation
!
            LOGICAL :: First_call = .TRUE., correct_inclusion, interp1, interp2
!
! -------------
! ......... Saving variables
! -------------
!
            SAVE First_call, phi, theta, omega, omega1, omega2, height, mbegin
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Inclusion_Status>
!
!
            IF(First_call) THEN
               WRITE(UNIT = *, FMT = 6000)
               mbegin = Num_HetRegions + Num_Boundaries
            ELSE
               GO TO 1000
            END IF
!
            ALLOCATE( phi(Num_InclZones),   theta(Num_InclZones),  height(Num_InclZones), &
     &                omega(Num_InclZones), omega1(Num_InclZones), omega2(Num_InclZones), STAT=ier )
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6850) 'phi, theta, height, omega, omega1, omega2', 'Inclusion_Status'
               STOP
            END IF
!
! ......... Initialization
!
            phi    = 0.0d0
            theta  = 0.0d0
            height = 0.0d0
!
            omega  = 0.0d0
            omega1 = 0.0d0
            omega2 = 0.0d0
!
            emissivity = 0.0d0
!
!
!***********************************************************************
!*                                                                     *
!*          Baseline computations for cylindrical shapes in            *
!*          cartesian coordinate systems                               *
!*                                                                     *
!***********************************************************************
!
!
            IF_Initial: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                      coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
               DO_InclZones: DO n=1,Num_InclZones
!
! ............... Cartesian grid system, cylindrically-shaped inclusion zones
!
                  IF_Shape: IF( (InclZone(n)%shape(1:2) == 'CY' .OR. InclZone(n)%shape(1:2) == 'cy' .OR. &
     &                           InclZone(n)%shape(1:2) == 'Cy') )  &
     &            THEN
!
                     CALL Determine_Orientation_Parameters( base1_center_coord = InclZone(n)%CylBase1Coord, &
     &                                                      base2_center_coord = InclZone(n)%CylBase2Coord, &
     &                                                      radius1 = InclZone(n)%CylBase1R, radius2 = InclZone(n)%CylBase2R, &
     &                                                      phi = phi(n), theta = theta(n), omega = omega(n), height = height(n) )
!
! ............... Cartesian grid system, elliptically-shaped inclusion zones
!
                  ELSE IF( (InclZone(n)%shape(1:2) == 'EL' .OR. InclZone(n)%shape(1:2) == 'El' .OR. &
     &                      InclZone(n)%shape(1:2) == 'el') )  &
     &            THEN
!
                     CALL Determine_Orientation_Parameters( base1_center_coord = InclZone(n)%EllBase1Coord, &
     &                                                      base2_center_coord = InclZone(n)%EllBase2Coord, &
     &                                                      LAxis1 = InclZone(n)%Base1LAxis, SAxis1 = InclZone(n)%Base1SAxis, &
     &                                                      LAxis2 = InclZone(n)%Base2LAxis, SAxis2 = InclZone(n)%Base2SAxis, &
     &                                                      phi = phi(n), theta = theta(n), omega1 = omega1(n), omega2 = omega2(n), height = height(n) )
!
                  END IF IF_Shape
!
               END DO DO_InclZones
!
            END IF IF_Initial
!
            First_call = .FALSE.
!
! -------------
! ......... Initialization
! -------------
!
 1000       Included_Element = .FALSE.
!
!
!***********************************************************************
!*                                                                     *
!*          For cartesian coordinate systems                           *
!*                                                                     *
!***********************************************************************
!
!
            IF_Assign: IF(coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'cart' .OR. coordinates(1:4) == 'Cart' .OR. &
     &                    coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &      THEN

!
               DO_Cartesian: DO n=1,Num_InclZones
!
! ............... Cartesian grid system, Cartesian-coordinate-defined inclusion zones
!
                  IF_BShape: IF( InclZone(n)%shape(1:4) == 'RECT' .OR. InclZone(n)%shape(1:4) == 'Rect' .OR. InclZone(n)%shape(1:4) == 'rect') THEN
!
                     IF( X >= InclZone(n)%LMin(1) .AND. X <= InclZone(n)%LMax(1) .AND.   &
     &                   Y >= InclZone(n)%LMin(2) .AND. Y <= InclZone(n)%LMax(2) .AND.   &
     &                   Z >= InclZone(n)%LMin(3) .AND. Z <= InclZone(n)%LMax(3) )       &
     &               THEN
!
                        medium = InclZone(n)%name
                        emissivity = InclZone(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = mbegin + InclZone_SequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                        Included_Element = .TRUE.
                        inclusion_id     = InclZone(n)%id
                        RETURN
!
                     END IF
!
! ............... Cartesian grid system, cylindrically-shaped inclusion zones
!
                  ELSE IF( InclZone(n)%shape(1:4) == 'CYLI' .OR. InclZone(n)%shape(1:4) == 'cyli' .OR. &
     &                     InclZone(n)%shape(1:4) == 'Cyli' ) &
     &            THEN
!
! .................. The coordinates in the 1st rotated system
!
                     X_1p =  (X - InclZone(n)%CylBase1Coord(1)) * COS(phi(n)) + (Y - InclZone(n)%CylBase1Coord(2)) * SIN(phi(n))
                     Y_1p = -(X - InclZone(n)%CylBase1Coord(1)) * SIN(phi(n)) + (Y - InclZone(n)%CylBase1Coord(2)) * COS(phi(n))
                     Z_1p =   Z - InclZone(n)%CylBase1Coord(3)
!
! .................. The coordinates in the 1st rotated system
!
                     X_2p =  X_1p * COS(theta(n)) + Z_1p * SIN(theta(n))
                     Z_2p = -X_1p * SIN(theta(n)) + Z_1p * COS(theta(n))
!                    Y_2p =  Y_1p
!
                     IF( InclZone(n)%CylBase1R /= 0.0d0 ) radius = MAX(InclZone(n)%CylBase1R, InclZone(n)%CylBase2R) - Z_2p * omega(n)
!
                     Rtest = SQRT( Z_2p * Z_2p + Y_1p * Y_1p )
!
                     IF( ( InclZone(n)%CylBase1R == 0.0d0 .AND. X_2p >= 0.0d0 .AND. X_2p <= height(n) .AND. Rtest >= InclZone(n)%CylRmin .AND. Rtest <= InclZone(n)%CylRmax ) .OR. &
     &                   ( InclZone(n)%CylBase1R /= 0.0d0 .AND. X_2p >= 0.0d0 .AND. X_2p <= height(n) .AND. Rtest <= radius ) ) &
     &               THEN
!
                        medium = InclZone(n)%name
                        emissivity = InclZone(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = mbegin + InclZone_SequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                        Included_Element = .TRUE.
                        inclusion_id     = InclZone(n)%id
                        RETURN
!
                     END IF
!
! ............... Cartesian grid system, elliptical inclusion zones
!
                  ELSE IF( InclZone(n)%shape(1:4) == 'ELLI' .OR. InclZone(n)%shape(1:4) == 'Elli' .OR. &
     &                     InclZone(n)%shape(1:4) == 'elli' ) &
     &            THEN
!
! .................. The angle between the original x- or y-axis and the long axis of the ellipse (rotated system)
!
                     phi(n) = InclZone(n)%LAxisAngle
!
                     SELECT CASE( InclZone(n)%EllipsePlane )
!
                     CASE( 'XY', 'xy', 'Xy', 'xY' )
!
                        IF( InclZone(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = Z - InclZone(n)%EllBase1Coord(3)
                           long_axis  = InclZone(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = InclZone(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = InclZone(n)%LAxis
                           short_axis = InclZone(n)%SAxis
                        END IF
!
                        X_1p =  (X - InclZone(n)%EllBase1Coord(1)) * COS(phi(n)) + (Y - InclZone(n)%EllBase1Coord(2)) * SIN(phi(n))
                        Y_1p = -(X - InclZone(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Y - InclZone(n)%EllBase1Coord(2)) * COS(phi(n))
!
                        Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Y_1p * Y_1p / ( short_axis * short_axis )
!
                     CASE( 'XZ', 'xz', 'Xz', 'xZ' )
!
                        IF( InclZone(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = Y - InclZone(n)%EllBase1Coord(2)
                           long_axis  = InclZone(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = InclZone(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = InclZone(n)%LAxis
                           short_axis = InclZone(n)%SAxis
                        END IF
!
                        X_1p =  (X - InclZone(n)%EllBase1Coord(1)) * COS(phi(n)) + (Z - InclZone(n)%EllBase1Coord(3)) * SIN(phi(n))
                        Z_1p = -(X - InclZone(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Z - InclZone(n)%EllBase1Coord(3)) * COS(phi(n))
!
                        Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     CASE( 'YZ', 'yz', 'Yz', 'yZ' )
!
                        IF( InclZone(n)%Base1LAxis /= 0.0d0 ) THEN
                           Delta      = X - InclZone(n)%EllBase1Coord(1)
                           long_axis  = InclZone(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                           short_axis = InclZone(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                        ELSE
                           long_axis  = InclZone(n)%LAxis
                           short_axis = InclZone(n)%SAxis
                        END IF
!
                        Y_1p =  (Y - InclZone(n)%EllBase1Coord(2)) * COS(phi(n)) + (Z - InclZone(n)%EllBase1Coord(3)) * SIN(phi(n))
                        Z_1p = -(Y - InclZone(n)%EllBase1Coord(2)) * SIN(phi(n)) + (Z - InclZone(n)%EllBase1Coord(3)) * COS(phi(n))
!
                        Rtest = Y_1p * Y_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     END SELECT
!
                     IF( Delta >= 0.0d0 .AND. Delta <= height(n) .AND. Rtest <= 2.5d-1 ) THEN
!
                        medium = InclZone(n)%name
                        emissivity = InclZone(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = mbegin + InclZone_SequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                        Included_Element = .TRUE.
                        inclusion_id     = InclZone(n)%id
                        RETURN
!
                     END IF
!
! ............... Cartesian grid system, spherically-shaped inclusion zones
!
                  ELSE IF( InclZone(n)%shape(1:4) == 'SPHE' .OR. InclZone(n)%shape(1:4) == 'sphe' .OR. InclZone(n)%shape(1:4) == 'Sphe') THEN
!
                     X_1p = (X - InclZone(n)%SphereCenterCoord(1))
                     Y_1p = (Y - InclZone(n)%SphereCenterCoord(2))
                     Z_1p = (Z - InclZone(n)%SphereCenterCoord(3))
!
                     Rtest = SQRT( X_1p * X_1p + Y_1p * Y_1p + Z_1p * Z_1p )
!
                     IF( Rtest >= InclZone(n)%SphRmin .AND. Rtest <= InclZone(n)%SphRmax ) THEN
!
                        medium = InclZone(n)%name
                        emissivity = InclZone(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = mbegin + InclZone_SequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                        Included_Element = .TRUE.
                        inclusion_id     = InclZone(n)%id
                        RETURN
!
                     END IF
!
! ............... Cartesian grid system, irregularly-shaped inclusion zones
!
                  ELSE IF( (InclZone(n)%shape(1:4) == 'IRRE' .OR. InclZone(n)%shape(1:4) == 'irre' .OR. InclZone(n)%shape(1:4) == 'Irre') ) THEN
!
                     m1      = InclZone(n)%IntTableNum1
                     interp1 = .FALSE.
                     IF( m1 > 0 ) THEN
                        table_1 = IntTable(m1)
                        IF( completed_interpolation(m1) ) interp1 = .TRUE.
                     END IF
!
                     m2      = InclZone(n)%IntTableNum2
                     interp2 = .FALSE.
                     IF( m2 > 0 ) THEN
                        table_2 = IntTable(m2)
                        IF( completed_interpolation(m2) ) interp2 = .TRUE.
                     END IF
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
!
                     IF_tables: IF( m1 == 0 .AND. m2 == 0 ) THEN
!
                        CALL Surface_Coordinate_Limits( subdomain        = InclZone(n),      &
     &                                                  XYZ_coordinates  = (/ X, Y, Z /),    &
     &                                                  XYZ_limits       = XYZ_limits,       &
     &                                                  success          = correct_inclusion )
!
                     ELSE
!
                        IF( (PRESENT(use_InterpData)) .AND. (interp1 .OR. interp2) ) THEN
!
                           CALL Surface_Coordinate_Limits( subdomain        = InclZone(n),       &
     &                                                     table_1          = table_1,           &
     &                                                     table_2          = table_2,           &
     &                                                     XYZ_coordinates  = (/ X, Y, Z /),     &
     &                                                     XYZ_limits       = XYZ_limits,        &
     &                                                     success          = correct_inclusion, &
     &                                                     use_InterpData   = use_InterpData,    &
     &                                                     completed_interpolation = completed_interpolation, &
     &                                                     available_InterpData    = available_InterpData )
!
                        ELSE
!
                           CALL Surface_Coordinate_Limits( subdomain        = InclZone(n),      &
     &                                                     table_1          = table_1,          &
     &                                                     table_2          = table_2,          &
     &                                                     XYZ_coordinates  = (/ X, Y, Z /),    &
     &                                                     XYZ_limits       = XYZ_limits,       &
     &                                                     success          = correct_inclusion )
!
                        END IF
!
                     END IF IF_tables
!
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>
!
                     IF( .NOT. correct_inclusion ) CYCLE DO_Cartesian
!
                     IF( X >= XYZ_limits(1,1) .AND. X <= XYZ_limits(1,2) .AND.   &
     &                   Y >= XYZ_limits(2,1) .AND. Y <= XYZ_limits(2,2) .AND.   &
     &                   Z >= XYZ_limits(3,1) .AND. Z <= XYZ_limits(3,2) )       &
     &               THEN
!
                        medium = InclZone(n)%name
                        emissivity = InclZone(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = mbegin + InclZone_SequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                        Included_Element = .TRUE.
                        inclusion_id     = InclZone(n)%id
                        RETURN
!
                     END IF
!
! ............... Cartesian grid system, unknown/unavailable shape of the inclusion zone
!
                  ELSE
!
! .................. Error
!
                     WRITE(*,6100) n, InclZone(n)%shape(1:4)
                     STOP
!
                  END IF IF_BShape
!
               END DO DO_Cartesian
!
!***********************************************************************
!*                                                                     *
!*          For cylindrical coordinate systems                         *
!*                                                                     *
!***********************************************************************
!
            ELSE IF(coordinates(1:4) == 'CYLI' .OR. coordinates(1:4) == 'cyli' .OR. coordinates(1:4) == 'Cyli') THEN
!
               DO_Cylindrical: DO n=1,Num_InclZones
!
! ............... Cylindrical grid system, cylindrical-coordinate-defined media inclusion zones
!
                  IF_BRegions: IF( InclZone(n)%shape(1:4) == 'CYLI' .OR. InclZone(n)%shape(1:4) == 'cyli' .OR. InclZone(n)%shape(1:4) == 'Cyli' ) THEN
!
                     IF( InclZone(n)%CylBase1R /= 0.0d0 ) radius = MAX(InclZone(n)%CylBase1R, InclZone(n)%CylBase2R) - Z * omega(n)
!
                     RMax = MAX( InclZone(n)%CylBase1Coord(3), InclZone(n)%CylBase2Coord(3) )
                     RMin = MIN( InclZone(n)%CylBase1Coord(3), InclZone(n)%CylBase2Coord(3) )
!
                     IF( ( InclZone(n)%CylBase1R == 0.0d0 .AND. X >= InclZone(n)%CylRmin .AND. X <= InclZone(n)%CylRmax .AND. Z >= RMin .AND. Z <= RMax ) .OR.   &
     &                   ( InclZone(n)%CylBase1R /= 0.0d0 .AND. X <= radius .AND. Z >= RMin .AND. Z <= RMax ) ) &
     &               THEN
!
                        medium = InclZone(n)%name
                        emissivity = InclZone(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = mbegin + InclZone_SequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                        Included_Element = .TRUE.
                        inclusion_id     = InclZone(n)%id
                        RETURN
!
                     END IF
!
! ............... Cartesian grid system, elliptical inclusion zones
!
                  ELSE IF( InclZone(n)%shape(1:4) == 'ELLI' .OR. InclZone(n)%shape(1:4) == 'Elli' .OR. &
     &                     InclZone(n)%shape(1:4) == 'elli' ) &
     &            THEN
!
! .................. The angle between the original angle-axis and the long axis of the ellipse (rotated system)
!
                     phi(n) = InclZone(n)%LAxisAngle
!
                     IF( InclZone(n)%Base1LAxis /= 0.0d0 ) THEN
                        Delta      = Z - InclZone(n)%EllBase1Coord(3)
                        long_axis  = InclZone(n)%Base1LAxis - 2.0d0 * Delta * omega1(n)
                        short_axis = InclZone(n)%Base1SAxis - 2.0d0 * Delta * omega2(n)
                     ELSE
                        long_axis  = InclZone(n)%LAxis
                        short_axis = InclZone(n)%SAxis
                     END IF
!
                     X_1p =  (X - InclZone(n)%EllBase1Coord(1)) * COS(phi(n)) + (Z - InclZone(n)%EllBase1Coord(3)) * SIN(phi(n))
                     Z_1p = -(X - InclZone(n)%EllBase1Coord(1)) * SIN(phi(n)) + (Z - InclZone(n)%EllBase1Coord(3)) * COS(phi(n))
!
                     Rtest = X_1p * X_1p / ( long_axis * long_axis ) + Z_1p * Z_1p / ( short_axis * short_axis )
!
                     IF( Delta >= 0.0d0 .AND. Delta <= height(n) .AND. Rtest <= 2.5d-1 ) THEN
!
                        medium = InclZone(n)%name
                        emissivity = InclZone(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = mbegin + InclZone_SequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                        Included_Element = .TRUE.
                        inclusion_id     = InclZone(n)%id
                        RETURN
!
                     END IF
!
! ............... Cylindrical grid system, spherically-shaped media inclusion zones
!
                  ELSE IF( InclZone(n)%shape(1:4) == 'SPHE' .OR. InclZone(n)%shape(1:4) == 'sphe' .OR. &
     &                     InclZone(n)%shape(1:4) == 'Sphe' ) &
     &            THEN
!
                     X_1p = (X - InclZone(n)%SphereCenterCoord(1))
                     Z_1p = (Z - InclZone(n)%SphereCenterCoord(3))
!
                     Rtest = SQRT( X_1p * X_1p + Z_1p * Z_1p )
!
                     IF( Rtest >= InclZone(n)%SphRmin .AND. Rtest <= InclZone(n)%SphRmax ) THEN
!
                        medium = InclZone(n)%name
                        emissivity = InclZone(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = mbegin + InclZone_SequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                        Included_Element = .TRUE.
                        inclusion_id     = InclZone(n)%id
                        RETURN
!
                     END IF
!
! ............... Cylindrical grid system, irregularly-shaped inclusion zones
!
                  ELSE IF( InclZone(n)%shape(1:4) == 'IRRE' .OR. InclZone(n)%shape(1:4) == 'irre' .OR. &
     &                     InclZone(n)%shape(1:4) == 'Irre' ) &
     &            THEN
!
                     CALL Surface_Coordinate_Limits( subdomain        = InclZone(n),      &
     &                                               RZ_coordinates   = (/ X, Z /),       &
     &                                               RZ_limits        = RZ_limits,        &
     &                                               success          = correct_inclusion )
!
                     IF( .NOT. correct_inclusion ) CYCLE DO_Cylindrical

!
                     IF( X >= RZ_limits(1,1) .AND. X <= RZ_limits(1,2) .AND. Z >= RZ_limits(2,1) .AND. Z <= RZ_limits(2,2) ) THEN
!
                        medium = InclZone(n)%name
                        emissivity = InclZone(n)%emissivity
!
! ..................... Listing by material number
!
                        IF(media_by_number) THEN
!
                           CALL Determine_MINC_Media_Number( M = mbegin + InclZone_SequenceNumber(n), &
     &                                                       N_0 = N_0, N_1 = N_1, N_2 = N_2, N_3 = N_3, N_4 = N_4, N_5 = N_5 )
!
                           IF( medium(1:1) == 'F' .OR. medium(1:1) == 'f' ) THEN
                              medium = ' '//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           ELSE
                              medium = '-'//CHAR(N_0+48)//CHAR(N_2+48)//CHAR(N_4+48)//CHAR(N_5+48)
                           END IF
!
                        END IF
!
                        Included_Element = .TRUE.
                        inclusion_id     = InclZone(n)%id
                        RETURN
!
                     END IF
!
! ............... Cylindrical grid system, unknown/unavailable shape of inclusion zone
!
                  ELSE
!
! .................. Error
!
                     WRITE(*,6102) n, InclZone(n)%shape(1:4)
                     STOP
!
                  END IF IF_BRegions
!
               END DO DO_Cylindrical
!
!
!
            END IF IF_Assign
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Included_Element 1.0 .................... 9 April     2015',6X,'Determinimg whether a cell belongs to a given inclusion zone')
!
 6100 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system of the grid is "Cartesian", ', &
     &          'but the shape of the inclusion zone #',I3,', is "',A,'": Unknown/Unavailable'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The coordinate system of the grid is "cylindrical", ', &
     &          'but the shape of the inclusion zone #',I3,', is "',A,'": Unknown/Unavailable'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6850 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Included_Element>
!
!
            RETURN
!
         END FUNCTION Included_Element
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Determine_Orientation_Parameters( base1_center_coord, base2_center_coord,           &
     &                                                radius1, radius2, LAxis1, SAxis1, LAxis2, SAxis2, &
     &                                                phi, theta, omega, omega1, omega2, height )
!
            IMPLICIT NONE
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), INTENT(IN), DIMENSION(3) :: base1_center_coord, base2_center_coord
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8), OPTIONAL :: SAxis1, SAxis2
!
            REAL(KIND = 8), INTENT(IN), OPTIONAL  :: radius1, radius2, LAxis1, LAxis2
!
            REAL(KIND = 8), INTENT(OUT), OPTIONAL :: omega, omega1, omega2
!
            REAL(KIND = 8), INTENT(OUT) :: phi, theta, height
!
            REAL(KIND = 8) :: DeltX, DeltY, DeltZ, length_xy, RadDif, LAxDif, SAxDif
!
! -------------
! ......... Double precision parameters
! -------------
!
            REAL(KIND = 8), PARAMETER :: pi = 3.14159265358979324D0
!
! -------------
! ......... Logical variables
! -------------
!
            LOGICAL :: First_call = .TRUE.
!
! -------------
! ......... Saving variables
! -------------
!
            SAVE First_call
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of <Determine_Orientation_Parameters>
!
!
            IF(First_call) THEN
               WRITE(UNIT = *, FMT = 6000)
               First_call = .FALSE.
            END IF
!
!***********************************************************************
!*                                                                     *
!*          Determine parameters for cylindrical subdomains            *
!*                                                                     *
!***********************************************************************
!
            DeltX = base2_center_coord(1) - base1_center_coord(1)
            DeltY = base2_center_coord(2) - base1_center_coord(2)
            DeltZ = base2_center_coord(3) - base1_center_coord(3)
!
            height    = SQRT( DeltX * DeltX + DeltY * DeltY + DeltZ * DeltZ )
            length_xy = SQRT( DeltX * DeltX + DeltY * DeltY )
!
! ......... The angle of rotation around the z-axis
!
            IF( ABS(DeltX) < 1.0d-7 .AND. ABS(DeltY) < 1.0d-7 ) THEN
               phi = 0.0d0
            ELSE IF( ABS(DeltY) < 1.0d-7 .AND. ABS(DeltX) > 1.0d-7  ) THEN
               phi = 0.0d0
            ELSE IF( ABS(DeltX) < 1.0d-7 .AND. ABS(DeltY) > 1.0d-7 .AND. DeltY > 0.0d0 ) THEN
               phi = 5.0d-1 * PI
            ELSE IF( ABS(DeltX) < 1.0d-7 .AND. ABS(DeltY) > 1.0d-7 .AND. DeltY < 0.0d0 ) THEN
               phi =-5.0d-1 * PI
            ELSE
               phi = ATAN( DeltY / DeltX )
            END IF
!
! ......... The 2nd angle of rotation
!
            IF( length_xy < 1.0d-7 .AND. ABS(DeltZ) < 1.0d-7 ) THEN
               theta = 0.0d0
            ELSE IF( ABS(DeltZ) < 1.0d-7 .AND. length_xy > 1.0d-7 ) THEN
               theta = 0.0d0
            ELSE IF( length_xy < 1.0d-7 .AND. ABS(DeltZ) > 1.0d-7 .AND. DeltZ > 0.0d0) THEN
               theta = 5.0d-1 * PI
            ELSE IF( length_xy < 1.0d-7 .AND. ABS(DeltZ) > 1.0d-7 .AND. DeltZ < 0.0d0) THEN
               theta =-5.0d-1 * PI
            ELSE
               theta = ATAN( DeltZ / length_xy )
            END IF
!
!***********************************************************************
!*                                                                     *
!*          Determine parameters for partial cone                      *
!*                                                                     *
!***********************************************************************
!
            IF_PCone: IF( present(radius1) .AND. present(radius2) ) THEN
!
               IF( radius1 < 1.0d-7 .AND. radius2 < 1.0d-7 ) GO TO 1000
!
               RadDif = ABS(radius1 - radius2)
!
               IF( MAX( radius1, radius2 ) > 0.0d0 .AND. RadDif > 0.0d0 ) THEN
                  omega = RadDif / height
               END IF
!
               RETURN
!
            END IF IF_PCone
!
!***********************************************************************
!*                                                                     *
!*          Determine parameters for partial elliptical cone           *
!*                                                                     *
!***********************************************************************
!
 1000       IF_PElCone: IF( present(LAxis1) .AND. present(SAxis1) .AND. present(LAxis2) .AND. present(SAxis2) ) THEN
!
               IF( LAxis1 < 1.0d-7 .AND. LAxis2 < 1.0d-7 .AND. SAxis1 < 1.0d-7 .AND. SAxis2 < 1.0d-7 ) RETURN
!
               IF( LAxis2 > 0.0d0 .AND. SAxis2 > 0.0d0 .AND. LAxis1 > 0.0d0 .AND. SAxis1 <= 1.0d-7 ) THEN
                  SAxis1 = SAxis2 * ( Laxis1 / Laxis2 )
               END IF
!
               LAxDif = 5.0d-1 * ( LAxis1 - LAxis2 )
!
               IF( LAxis1 > 0.0d0 .AND. SAxis1 > 0.0d0 .AND. LAxis2 > 0.0d0 .AND. SAxis2 <= 1.0d-7 ) THEN
                  SAxis2 = SAxis1 * ( Laxis2 / Laxis1 )
               END IF
!
               SAxDif = 5.0d-1 * ( SAxis1 - SAxis2 )
!
               omega1 = LAxDif / height
               omega2 = SAxDif / height
!
               RETURN
!
            END IF IF_PElCone
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Determine_Orientation_Parameters 1.0 .... 8 April     2015',6X,'Determine orientation parameters and lengths in cylindrical and elliptical domains')
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Determine_Orientation_Parameters>
!
!
            RETURN
!
         END SUBROUTINE Determine_Orientation_Parameters
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Read_Tabular_Data( table_number, file_name, number_of_rows, number_of_columns, read_data_by_row,      &
     &                                 read_data_format, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2 , &
     &                                 interpolation_search_radius, num_interrogated_points  )
!
            IMPLICIT NONE
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: table
            REAL(KIND = 8), INTENT(IN),  OPTIONAL     :: interpolation_search_radius
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER, INTENT(IN) :: table_number
!
            INTEGER, INTENT(IN) :: number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2
!
            INTEGER, INTENT(IN), OPTIONAL :: num_interrogated_points
!
            INTEGER :: i, j, m, ier, U_limit
!
! -------------
! ......... Character variables
! -------------
!
            CHARACTER(LEN =  8), INTENT(IN) :: file_name
            CHARACTER(LEN = 50), INTENT(IN) :: read_data_format
!
! -------------
! ......... Logical variables
! -------------
!
            LOGICAL, INTENT(IN) :: read_data_by_row
!
            LOGICAL :: First_call = .TRUE., file_exists
!
! -------------
! ......... Saving variables
! -------------
!
            SAVE First_call
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of <Read_Tabular_Data>
!
!
            IF(First_call) THEN
               WRITE(UNIT = *, FMT = 6000)
               First_call = .FALSE.
            END IF
!
!***********************************************************************
!*                                                                     *
!*           Determine length of table and allocate memory             *
!*                                                                     *
!***********************************************************************
!
            m = table_number
!
            IntTable(m)%FileName = file_name
            IntTable(m)%UnitNum  = 50 + m
!
            IntTable(m)%NumDataPoints = number_of_rows
!
            IF( PRESENT(interpolation_search_radius) ) IntTable(m)%SearchRadius = interpolation_search_radius
!
            IF( PRESENT(num_interrogated_points) ) IntTable(m)%NumSearchPoints = num_interrogated_points
!
            ALLOCATE( IntTable(m)%dv(number_of_rows), IntTable(m)%iv1(number_of_rows), IntTable(m)%iv2(number_of_rows), STAT=ier )
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6100) 'IntTable(m)%dv; iv1; iv2'
               STOP
            END IF
!
! ......... Reading the tabular data for interpolation Table 1
!
!
            INQUIRE(FILE = ADJUSTL(TRIM(file_name)), EXIST = file_exists)
!
            IF(file_exists) THEN
               OPEN( UNIT = IntTable(m)%UnitNum, FILE = ADJUSTL(TRIM(file_name)) )
            ELSE
               WRITE( UNIT = *, FMT = 6105 ) ADJUSTL(TRIM(file_name))
               STOP
            END IF
!
            REWIND ( UNIT = IntTable(m)%UnitNum )
!
! >>>>>>>>>
!
            U_limit = MIN( number_of_columns, MAX( RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2 ) )
!
            IF_ByRow: IF( read_data_by_row ) THEN
!
               IF_Format1: IF( read_data_format == '(*)' .OR. read_data_format == '*' ) THEN
!
                  ALLOCATE( table(U_limit) )
!
                  DO j = 1, number_of_rows
                     READ( UNIT = IntTable(m)%UnitNum, FMT = * )  ( table(i), i=1,U_limit )
                     IntTable(m)%dv(j)  = table(RowCol_DepVariable)
                     IntTable(m)%iv1(j) = table(RowCol_IndVariable_1)
                     IntTable(m)%iv2(j) = table(RowCol_IndVariable_2)
                  END DO
!
               ELSE
!
                  DO j = 1, number_of_rows
                     READ( UNIT = IntTable(m)%UnitNum, FMT = read_data_format ) ( table(i), i=1,U_limit )
                     IntTable(m)%dv(j)  = table(RowCol_DepVariable)
                     IntTable(m)%iv1(j) = table(RowCol_IndVariable_1)
                     IntTable(m)%iv2(j) = table(RowCol_IndVariable_2)
                  END DO
!
               END IF IF_Format1
!
            ELSE
!
               IF_Format2: IF( read_data_format == '(*)' .OR. read_data_format == '*' ) THEN
!
                  ALLOCATE( table(number_of_rows) )
!
                  DO j = 1,U_limit
                     READ( UNIT = IntTable(m)%UnitNum, FMT = * )  ( table(i), i=1,number_of_rows )
                     IF(j == RowCol_DepVariable)   IntTable(m)%dv(1:number_of_rows)  = table(1:number_of_rows)
                     IF(j == RowCol_IndVariable_1) IntTable(m)%iv1(1:number_of_rows) = table(1:number_of_rows)
                     IF(j == RowCol_IndVariable_2) IntTable(m)%iv2(1:number_of_rows) = table(1:number_of_rows)
                  END DO
!
               ELSE
!
                  DO j = 1,U_limit
                     READ( UNIT = IntTable(m)%UnitNum, FMT = read_data_format )  ( table(i), i=1,number_of_rows )
                     IF(j == RowCol_DepVariable)   IntTable(m)%dv(1:number_of_rows)  = table(1:number_of_rows)
                     IF(j == RowCol_IndVariable_1) IntTable(m)%iv1(1:number_of_rows) = table(1:number_of_rows)
                     IF(j == RowCol_IndVariable_2) IntTable(m)%iv2(1:number_of_rows) = table(1:number_of_rows)
                  END DO
!
               END IF IF_Format2
!
            END IF IF_ByRow
!
            DEALLOCATE (table)
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Read_Tabular_Data 1.0 .................. 28 April     2015',6X,'Read the tabular data that are used for interpolation')
!
 6100 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <Read_Tabular_Data> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6105 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'In subroutine <Read_Tabular_Data>, the file <',A,'> with tabular data for interpolation does NOT exist',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Read_Tabular_Data>
!
!
            RETURN
!
         END SUBROUTINE Read_Tabular_Data
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Determine_MINC_Media_Number( M, N_0, N_1, N_2, N_3, N_4, N_5 )
!
            IMPLICIT NONE
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER, INTENT(IN)  :: M
!
            INTEGER, INTENT(OUT) :: N_0, N_1, N_2, N_3, N_4, N_5
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of <Determine_MINC_Media_Number>
!
!
            N_0 = M / 1000
            N_1 = MOD(M,1000)
            N_2 = N_1 /100
            N_3 = MOD(N_1,100)
            N_4 = N_3 / 10
            N_5 = MOD(N_3,10)
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Determine_MINC_Media_Number>
!
!
            RETURN
!
         END SUBROUTINE Determine_MINC_Media_Number
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      END MODULE Het_Region_Definition
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE Mesh_Maker(Flag_MINC)
!
         USE MeshMaker_Data
         USE Grid_Generation_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*            Executive routine for generating TOUGH+ grid             *
!*        and assigning initial property distribution by region        *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 6) :: KeyWord
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: exists
      LOGICAL, INTENT(OUT) :: Flag_MINC
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Mesh_Maker
!
!
      WRITE(*,6000)
!
! ----------
! ... Initialization
! ----------
!
      Flag_MINC = .FALSE.        !  Indicates non-MINC type of MESH
!
      Flag_HetRegions = .FALSE.  !  Indicates uniform geologic media distribution
      Flag_Boundaries = .FALSE.  !  Indicates no-flow boundaries
      Flag_ExclZones  = .FALSE.  !  Indicates regularly-shaped domain
!
! >>>
! >>>
! >>>
!
 1000 READ(UNIT = *, FMT = 5001) KeyWord
!
!
!
      CASE_KeyWord: SELECT CASE(KeyWord(1:6))
!
!
!***********************************************************************
!*                                                                     *
!*    KeyWord = 'HET': Assign heterogeneous property distributions     *
!*                     by region (geometry-based)                      *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('>>>HET','>>>Het','>>>het','>>>REG','>>>Reg','>>>reg')
!
         WRITE(*,6015)
!
         CALL Define_Heterogeneous_Regions
!
         Flag_HetRegions = .TRUE.  !  Indicates heterogeneous distribution of geologic media
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*    KeyWord = 'EXC': Describe exclusion zones of the original grid   *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('>>>EXC','>>>Exc','>>>exc')
!
         WRITE(*,6016)
!
         Flag_ExclZones = .TRUE.  !  Indicates irregularly-shaped domain
!
         CALL Define_Exclusion_Zones
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*    KeyWord = 'INC': Describe inclusion zones of the original grid   *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('>>>INC','>>>Inc','>>>inc')
!
         WRITE(*,6019)
!
         Flag_InclZones = .TRUE.  !  Indicates irregularly-shaped domain
!
         CALL Define_Inclusion_Zones
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*    KeyWord = 'BOU': Assign boundary subdomains                      *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('>>>BOU','>>>Bou','>>>bou')
!
         WRITE(*,6018)
!
         CALL Define_Boundaries
!
         Flag_Boundaries = .TRUE.  !  Indicates boundaries of time-variant or invariant conditions
!
         GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*    KeyWord = 'MINC': GENERATE GRID FOR FRACTURED MEDIA              *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('>>>MIN','>>>Min','>>>min')
!
         Flag_MINC = .TRUE.
!
         INQUIRE(FILE='MINC', EXIST = exists)
!
         IF(exists) THEN
            WRITE(*,6021)
            OPEN(UNIT = MINC_Unit, FILE = 'MINC', STATUS = 'OLD')
         ELSE
            WRITE(*,6022)
            OPEN(UNIT = MINC_Unit, FILE = 'MINC', STATUS = 'NEW')
         END IF
!
         CALL MINC
!
         RETURN
!
!
!***********************************************************************
!*                                                                     *
!*    KeyWord = 'GRID': READ THE GRID SPECIFICS                        *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('>>>DIS','>>>Dis','>>>dis')
!
!        IF( .NOT. Flag_HetRegions ) THEN
!           WRITE( UNIT = *, FMT = 6200 )
!           STOP
!        END IF
!
         GO TO 2000
!
!
!***********************************************************************
!*                                                                     *
!*    KeyWord = ' ': IGNORE BLANKS BETWEEN DATA BLOCKS                 *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('      ')
!
      GO TO 1000
!
!
!***********************************************************************
!*                                                                     *
!*    Default case: Unknown keyword - The simulation is aborted        *
!*                                                                     *
!***********************************************************************
!
!
      CASE DEFAULT
!
! ...... Printing the unknown keyword, and then stopping
!
         WRITE(*,6100) KeyWord
         STOP
!
      END SELECT CASE_KeyWord
!
! >>>
! >>> Read the grid data
! >>>
!
 2000 CASE_GridNum: SELECT CASE(grid_numbering_system(1:2))
!
!
!***********************************************************************
!*                                                                     *
!*    GENERATE CYLINDRICAL GRID BY COLUMNS                             *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('CO','Co','co')
!
         WRITE(*,6012)
!
         CALL RZ2D
!
         CALL WRZ2D(1)                                      ! Grid by columns
!
!
!***********************************************************************
!*                                                                     *
!*    GENERATE CYLINDRICAL GRID BY LAYERS                              *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('LA','La','la')
!
         WRITE(*,6012)
!
         CALL RZ2D
!
         CALL WRZ2D(2)                                      ! Grid by layers
!
!
!***********************************************************************
!*                                                                     *
!*    GENERATE CARTESIAN CYLINDRICAL GRID BY LAYERS                    *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('ST','St','st')
!
         WRITE(*,6014)
!
         CALL GXYZ
!
!
!***********************************************************************
!*                                                                     *
!*    KeyWord = '  ' or '<<<': END OF THE DATA SET                     *
!*                                                                     *
!***********************************************************************
!
!
      CASE ('<<<')
!
         RETURN
!
      END SELECT CASE_GridNum
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
!6000 FORMAT(/,'Mesh_Maker        1.0   29 August    2004',6X,'EXECUTIVE ROUTINE FOR INTERNAL MESH GENERATION')
!6000 FORMAT(/,'Mesh_Maker        1.5    8 August    2014',6X,'EXECUTIVE ROUTINE FOR INTERNAL MESH GENERATION')
 6000 FORMAT(/,'Mesh_Maker 2.0 ......................... 26 December  2014',6X,'EXECUTIVE ROUTINE FOR INTERNAL MESH GENERATION')
!
 5001 FORMAT(A6)
!
 6012 FORMAT(/,' ',131('*'),   &
     &       /,' *',20X,'"Mesh_Maker" - RZ2D: GENERATE 2-D R-Z MESH',T132,'*',   &
     &       /,' ',131('*'))
 6014 format(/,' ',131('*'),   &
     &       /,' *',20X,'"Mesh_Maker" - XYZ: GENERATE 1, 2, OR 3-D CARTESIAN MESH',T132,'*',   &
     &       /,' ',131('*'))
 6015 format(/,' ',131('*'),   &
     &       /,' *',20X,'"Mesh_Maker" - HET: Define regions of heterogeneous media in the grid',T132,'*',   &
     &       /,' ',131('*'))
 6016 format(/,' ',131('*'),   &
     &       /,' *',20X,'"Mesh_Maker" - EXC: Describe irregularly-shaped grid by defining grid exclusion zones',T132,'*',   &
     &       /,' ',131('*'))
 6018 format(/,' ',131('*'),   &
     &       /,' *',20X,'"Mesh_Maker" - BOU: Define boundary regions/subdomains of the grid',T132,'*',   &
     &       /,' ',131('*'))
 6019 format(/,' ',131('*'),   &
     &       /,' *',20X,'"Mesh_Maker" - INC: Describe irregularly-shaped grid by defining grid inclusion zones',T132,'*',   &
     &       /,' ',131('*'))
!
 6020 format(/,' ',131('*'),   &
     &       /,' *',18X,'"Mesh_Maker" - MINC: GENERATE MULTIPLE INTERACTING CONTINUA MESH FOR FRACTURED MEDIUM',T132,'*',   &
     &       /,' ',131('*')/)
!
 6021 FORMAT(' FILE "MINC" EXISTS ==> OPEN AS AN OLD FILE')
 6022 FORMAT(' FILE "MINC" DOES NOT EXIST ==> OPEN AS A NEW FILE')
!
 6024 FORMAT( /,' ',131('*'),   &
     &       //,' MESH GENERATION COMPLETE ==> EXIT FROM THE "Mesh_Maker" FACILITY',/)
!
 6100 FORMAT(' "Mesh_Maker" HAS READ THE UNKNOWN KEYWORD ',A5,'  ==> >>>>> STOP EXECUTION')
!
 6200 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Mesh_Maker>: There is no data block describing homogeneous or heterogeneous media',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Mesh_Maker
!
!
      RETURN
!
      END SUBROUTINE Mesh_Maker
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Define_Reference_Surface
!
            USE MeshMaker_Data, ONLY: coordinates, NXMax, NYMax, NZMax, max_offset
            USE Het_Region_Definition
            USE Grid_Generation_Parameters
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*   Routine defining RefSurfaces of heterogeneous media in the grid   *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: X_min, Y_min, X_max, Y_max, X_shift, Y_shift
!
            REAL(KIND = 8) :: exponent, sign1, RefSurface_offset, interpolation_search_radius
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), DIMENSION(0:10) :: RefSurface_EquCoeff_A, RefSurface_EquCoeff_B, RefSurface_EquCoeff_C
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: ier, n, i1, m
!
            INTEGER :: equation_order_of_RefSurface, num_interrogated_points
!
            INTEGER :: number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2
!
! -------------
! ......... CHARACTER variables
! -------------
!
            CHARACTER(LEN =  1) :: vertical_or_stratigraphic
            CHARACTER(LEN =  2) :: length_units
            CHARACTER(LEN =  5) :: RefSurface_name
            CHARACTER(LEN =  8) :: RefSurface_Data_FileName
            CHARACTER(LEN = 12) :: RefSurface_shape, type_of_equation
            CHARACTER(LEN = 50) :: read_data_format
!
            CHARACTER(LEN = 1) :: dependent_variable_of_RefSurface
!
! -------------
! ......... LOGICAL variables
! -------------
!
            LOGICAL :: read_data_by_row
!
! -------
! ... Namelists
! -------
!
      NAMELIST/ RefSurface_GeneralInfo / RefSurface_name, RefSurface_shape, length_units, max_offset
!
      NAMELIST/ Rectangular_RefSurface / X_min, Y_min, X_max, Y_max
!
      NAMELIST/ Irregular_RefSurface   / dependent_variable_of_RefSurface, type_of_equation, X_min, Y_min, X_max, Y_max,  &
     &                                   RefSurface_Data_FileName, vertical_or_stratigraphic, RefSurface_offset
!
      NAMELIST/ Irregular_RefSurface_Equation / equation_order_of_RefSurface, X_shift, Y_shift, exponent, sign1, &
     &                                          RefSurface_EquCoeff_A, RefSurface_EquCoeff_B, RefSurface_EquCoeff_C
!
      NAMELIST/ Irregular_RefSurface_Table / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                       read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Define_Reference_Surface>
!
!
            WRITE(*,6000)
!
!***********************************************************************
!*                                                                     *
!*         Read the specifics of the Reference Surface                 *
!*                                                                     *
!***********************************************************************
!
!
! ............ Initializations - Namelist components
!
               RefSurface_shape = '           '
               RefSurface_name  = '     '
               length_units = 'm'
!
               RefSurface%name  = '     '
               RefSurface%shape = '           '
               RefSurface%units = 'm'
!
               max_offset = 0
!
! -----------------
! ............ Read general info on the heterogeneous RefSurface
! -----------------
!
               READ (UNIT = *, NML = RefSurface_GeneralInfo, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6600) '<RefSurface_GeneralInfo>'
                  STOP
               END IF
!
               SELECT CASE( RefSurface_shape(1:1) )
               CASE('R', 'r', 'I', 'i' )
                  CONTINUE
               CASE DEFAULT
                  WRITE( UNIT = *, FMT = 6010 ) RefSurface_shape
                  STOP
               END SELECT
!
               IF( coordinates(1:2) /= 'CA' .AND. coordinates(1:2) /= 'Ca' .AND. coordinates(1:2) /= 'ca' .AND.  &
     &             coordinates(1:2) /= 'TR' .AND. coordinates(1:2) /= 'Tr' .AND. coordinates(1:2) /= 'tr' )      &
     &         THEN
                  WRITE(*,6520)
                  STOP
               END IF
!
! ............ Assignment of the namelist values
!
               RefSurface%name  = RefSurface_name
               RefSurface%shape = RefSurface_shape
               RefSurface%units = length_units
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    Reading the shape-specific data from the remaining namelists     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
               SELECT CASE(RefSurface_shape(1:1))
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Rectangular Reference Surface
! >>>>>>>>>>>>
!
               CASE( 'R', 'r' )
!
! ............... Initializations
!
                  X_min = 0.0d0
                  Y_min = 0.0d0
                  X_max = 0.0d0
                  Y_max = 0.0d0
!
                  RefSurface%LMin = 0.0d0     ! ... Whole array operations
                  RefSurface%LMax = 0.0d0
!
! ............... Reading the rectangle-related data using a namelist
!
                  READ (UNIT = *, NML = Rectangular_RefSurface, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6600) '<Rectangular_RefSurface>'
                     STOP
                  END IF
!
                  RefSurface%LMin(1) = X_min
                  RefSurface%LMin(2) = Y_min
!
                  RefSurface%LMax(1) = X_max
                  RefSurface%LMax(2) = Y_max
!
! ............... Ensuring that the LMin of the range is not >= of the LMax
!
                  IF( RefSurface%LMin(1) >= RefSurface%LMax(1) ) THEN
                     WRITE(*,6050) RefSurface%LMin(1), RefSurface%LMax(1)
                     STOP
                  END IF
!
                  IF( RefSurface%LMin(2) >= RefSurface%LMax(2) ) THEN
                     WRITE(*,6051) RefSurface%LMin(2), RefSurface%LMax(2)
                     STOP
                  END IF
!
                  GO TO 1500
!
               END SELECT
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Irregular RefSurfaces
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               type_of_equation = '           '  ! ... Initialization
!
               X_min = 0.0d0
               Y_min = 0.0d0
               X_max = 0.0d0
               Y_max = 0.0d0
!
               X_shift = 0.0d0
               Y_shift = 0.0d0
!
               dependent_variable_of_RefSurface = ' '
               equation_order_of_RefSurface     = 0
!
               RefSurface_Data_FileName = '        '
!
               sign1    = 1.0d0
               exponent = 0.0d0
               RefSurface_offset = 0.0d0
!
               interpolation_search_radius = 0.0d0
!
               vertical_or_stratigraphic = ' '
!
! ............ Initializations - User-defined variable components
!
               RefSurface%TypeEqu1 = '           '
               RefSurface%TypeEqu2 = '           '
!
               RefSurface%DependentVar = ' '
!
               RefSurface%OrderEquSurf1 = 0
               RefSurface%OrderEquSurf2 = 0
!
               RefSurface%LMin = 0.0d0     ! ... Whole array operations
               RefSurface%LMax = 0.0d0
!
               RefSurface%L1Shift = 0.0d0
               RefSurface%L2Shift = 0.0d0
!
               RefSurface%sign1  = 1.0d0
               RefSurface%sign2  = 1.0d0
               RefSurface%expon1 = 0.0d0
               RefSurface%expon2 = 0.0d0
!
               RefSurface%IntTableNum1 = 0
               RefSurface%IntTableNum2 = 0
!
               RefSurface%thick0 = 0.0d0
               RefSurface%thick1 = 0.0d0
               RefSurface%thick2 = 0.0d0
!
               RefSurface%VertOrStrat1 = ' '
               RefSurface%VertOrStrat2 = ' '
               RefSurface%RefSurf1     = ' '
               RefSurface%RefSurf2     = ' '
!
! ............ Reading the irregular zone data using a namelist
!
               READ (UNIT = *, NML = Irregular_RefSurface, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6600) '<Irregular_RefSurface>'
                  STOP
               END IF
!
! ............ Checking the dependent variable
!
               IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &             coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &         THEN

!
                  IF( dependent_variable_of_RefSurface /= 'Z' .AND. dependent_variable_of_RefSurface /= 'z' ) THEN
                     WRITE( UNIT = *, FMT = 6608 ) 'cartesian',  'Z or z'
                     STOP
                  END IF
!
               END IF
!
! ............ Assignment
!
               RefSurface%LMin(1) = X_min
               RefSurface%LMin(2) = Y_min
!
               RefSurface%LMax(1) = X_max
               RefSurface%LMax(2) = Y_max
!
               RefSurface%DependentVar = dependent_variable_of_RefSurface
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Setting the type of equation of the reference surface
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( type_of_equation(1:4) == 'POLY' .OR. type_of_equation(1:4) == 'Poly' .OR. type_of_equation(1:4) == 'poly' ) THEN
!
! ............ Polynomial equation
!
                  type_of_equation(1:4) = 'Poly'
                  IF( type_of_equation(5:6) == '/E' .OR. type_of_equation(5:6) == '/e' ) THEN
                     type_of_equation(5:6) = '/E'
                  ELSE IF( type_of_equation(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) type_of_equation(1:6)
                     STOP
                  END IF
!
! ............ Power equation
!
               ELSE IF( type_of_equation(1:4) == 'POWE' .OR. type_of_equation(1:4) == 'Powe' .OR. type_of_equation(1:4) == 'powe' ) THEN
!
                  type_of_equation(1:4) = 'Powe'
                  IF( type_of_equation(5:6) == '/E' .OR. type_of_equation(5:6) == '/e' ) THEN
                     type_of_equation(5:6) = '/E'
                  ELSE IF( type_of_equation(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) type_of_equation(1:6)
                     STOP
                  END IF
!
! ............ No equation: interpolation from a grid/tabular data set
!
               ELSE IF( type_of_equation(1:4) == 'INTE' .OR. type_of_equation(1:4) == 'Inte' .OR. type_of_equation(1:4) == 'inte' ) THEN
!
                  ALLOCATE( IntTable(1), temp_IntTable(1), STAT=ier )
!
! ............ Unknown option: ERROR!
!
               ELSE
                  WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation(1:6)
                  STOP
               END IF
!
               RefSurface%TypeEqu1 = type_of_equation
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For Cartesian coodinates
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( ( NXMax >  1 .AND. NYMax == 1 .AND. NZMax == 1 ) .OR. ( NXMax == 1 .AND. NYMax >  1 .AND. NZMax == 1 ) .OR.  &
     &             ( NXMax == 1 .AND. NYMax == 1 .AND. NZMax >  1 ) )    &
     &         THEN
                  WRITE( UNIT = *, FMT = 6606 ) 'cartesian'
                  STOP
               END IF
!
               IF_Inte: IF( type_of_equation(1:4) == 'INTE' .OR. type_of_equation(1:4) == 'Inte' .OR. type_of_equation(1:4) == 'inte' ) THEN
!
                  type_of_equation(1:4) = 'Inte'
!
                  m = 1
                  Num_IntTables           = m
                  RefSurface%IntTableNum1 = Num_IntTables
!
! ............... Reading info describing the structure of interpolation Table 1
!
                  READ (UNIT = *, NML = Irregular_RefSurface_Table, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                    WRITE(UNIT = *, FMT = 6600) '<Irregular_RefSurface_Table>'
                     STOP
                  END IF
!
                  CALL Read_Tabular_Data( table_number = m, file_name = RefSurface_Data_FileName, number_of_rows = number_of_rows, &
     &                                    number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,                  &
     &                                    read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                    RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                    interpolation_search_radius = interpolation_search_radius, &
     &                                    num_interrogated_points     = num_interrogated_points )
!
                  RefSurface%thick0 = RefSurface_offset
!
                  GO TO 1500
!
               END IF IF_Inte
!
! ............ Initializations
!
               RefSurface_EquCoeff_A = 0.0d0
               RefSurface_EquCoeff_B = 0.0d0
               RefSurface_EquCoeff_C = 0.0d0
!
! ............ Reading the data for surface 1
!
               READ (UNIT = *, NML = Irregular_RefSurface_Equation, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6600) '<Irregular_RefSurface_Equation>'
                  STOP
               END IF
!
               RefSurface%OrderEquSurf1 = equation_order_of_RefSurface
!
               i1 = equation_order_of_RefSurface
!
               ALLOCATE ( RefSurface%Equ1CoeffA(0:i1), RefSurface%Equ1CoeffB(0:i1) )
!
               RefSurface%Equ1CoeffA = 0.0d0
               RefSurface%Equ1CoeffB = 0.0d0
!
               IF_Surf1: IF( RefSurface%TypeEqu1(1:4) == 'Powe' ) THEN
                  ALLOCATE ( RefSurface%Equ1CoeffC(0:i1) )
                  RefSurface%Equ1CoeffC = 0.0d0
               ELSE
                  IF(i1 > 1 ) THEN
                     ALLOCATE ( RefSurface%Equ1CoeffC(1:Max(1,i1)) )
                     RefSurface%Equ1CoeffC = 0.0d0
                  END IF
               END IF IF_Surf1
!
! ............ Assignments
!
               RefSurface%sign1  = sign1
               RefSurface%expon1 = exponent
!
               RefSurface%Equ1CoeffA(0:i1) = RefSurface_EquCoeff_A(0:i1)
               RefSurface%Equ1CoeffB(0:i1) = RefSurface_EquCoeff_B(0:i1)
!
               IF( RefSurface%TypeEqu1(1:4) == 'Powe' ) THEN
                  RefSurface%Equ1CoeffC(0:i1) = RefSurface_EquCoeff_C(0:i1)
               ELSE
                  IF(i1 > 1) RefSurface%Equ1CoeffC(1:i1-1) = RefSurface_EquCoeff_C(1:i1-1)
               END IF
!
               RefSurface%L1Shift(1) = X_shift
               RefSurface%L1Shift(2) = Y_shift
!
!***********************************************************************
!*                                                                     *
!*     Convert the Reference Surface quantities into SI units (m)      *
!*                                                                     *
!***********************************************************************
!
! ......... Conversion into METERS from INCHES
!
 1500       IF( RefSurface%units == 'IN' .OR. RefSurface%units == 'in' .OR. RefSurface%units == 'In' ) THEN
!
               RefSurface%LMin(1:3) = RefSurface%LMin(1:3) * 2.54d-2
               RefSurface%LMax(1:3) = RefSurface%LMax(1:3) * 2.54d-2
!
               RefSurface%CylBase1Coord(1:3) = RefSurface%CylBase1Coord(1:3) * 2.54d-2
!
               RefSurface%L1Shift(1:3) = RefSurface%L1Shift(1:3) * 2.54d-2
!
               RefSurface%thick0 = RefSurface%thick0 * 2.54d-2
!
               RETURN
!
            END IF
!
! ......... Conversion into METERS from FEET
!
            IF( RefSurface%units == 'FT' .OR. RefSurface%units == 'ft' .OR. RefSurface%units == 'Ft' ) THEN
!
               RefSurface%LMin(1:3) = RefSurface%LMin(1:3) * 3.038d-1
               RefSurface%LMax(1:3) = RefSurface%LMax(1:3) * 3.038d-1
!
               RefSurface%L1Shift(1:3) = RefSurface%L1Shift(1:3) * 3.038d-1
!
               RefSurface%thick0 = RefSurface%thick0 * 3.038d-1
!
               RETURN
!
            END IF
!
! ......... Conversion into METERS from KM
!
            IF( RefSurface%units == 'KM' .OR. RefSurface%units == 'km' .OR. RefSurface%units == 'Km' ) THEN
!
               RefSurface%LMin(1:3) = RefSurface%LMin(1:3) * 1.0d3
               RefSurface%LMax(1:3) = RefSurface%LMax(1:3) * 1.0d3
!
               RefSurface%L1Shift(1:3) = RefSurface%L1Shift(1:3) * 1.0d3
!
               RefSurface%thick0 = RefSurface%thick0 * 1.0d3
!
               RETURN
!
            END IF
!
! ......... Conversion into METERS from MM
!
            IF( RefSurface%units == 'MM' .OR. RefSurface%units == 'mm' .OR. RefSurface%units == 'Mm' ) THEN
!
               RefSurface%LMin(1:3) = RefSurface%LMin(1:3) * 1.0d-3
               RefSurface%LMax(1:3) = RefSurface%LMax(1:3) * 1.0d-3
!
               RefSurface%L1Shift(1:3) = RefSurface%L1Shift(1:3) * 1.0d-3
!
               RefSurface%thick0 = RefSurface%thick0 * 1.0d-3
!
               RETURN
!
            END IF
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Define_Reference_Surface 2.0 ........... 18 January   2015',6X,'Defining the Reference Surface for thin areal grids')
!
 6010 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The shape of the Reference Surface "',A,'": Unknown/Unavailable option'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6050 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input X_min =',ES12.5,' of the Reference Surface is larger than (or equal to) the X_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6051 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input Y_min =',ES12.5,' of the Reference Surface is larger than (or equal to) the Y_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6101 FORMAT(T5,'Memory allocation to arrays <',A,'> in subroutine <Define_Reference_Surface> was successful')
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <Define_Reference_Surface> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6105 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The subdomain shape (read by "RefSurface_shape" = ',a11,' in subroutine <Define_Reference_Surface>) is unavailable',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6401 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Reference_Surface>: The Reference Surface has an irregular shape that is described by ',/, &
     &       T10,'                                               an unknown/unavailable type of equation (= <',A,'>)',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6402 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Reference_Surface: The RefSurface #',i3.3,' is defined by two irregular surfaces of type "FIXED" but the ' ,/,  &
     &       T10,'                                              name of the needed reference interpolation file is not defined',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6515 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Reference_Surface>: The dataset/namelist <',A,'> must be ended by the "<<<" descriptor'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6520 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Reference_Surface>: It is not possible to have a Reference Surface in anything but Cartesian coordinates ',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6600 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Reference_Surface>: There is a problem reading the namelist <',A,'>',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6604 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Reference_Surface>: The orders of the equation describing the the Reference Surface is < 0',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6606 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Reference_Surface>: This ',A,' system has 1 active dimension - it is not possible for the ', /, &
     &       T10,'                                           the Reference Surface to be irregularly shaped',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6608 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Reference_Surface>: In this ',A,' system, one or more of the dependent variables defining the ', /,&
     &       T10,'                                           Reference Surface is not among the possible options (',A,')',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Define_Reference_Surface>
!
!
            RETURN
!
         END SUBROUTINE Define_Reference_Surface
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Define_Heterogeneous_Regions
!
            USE MeshMaker_Data, ONLY: coordinates, NXMax, NYMax, NZMax
            USE Het_Region_Definition
            USE Grid_Generation_Parameters
            USE Utility_Functions
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    Routine defining regions of heterogeneous media in the grid      *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: D1, D2, D3
!
            REAL(KIND = 8) :: X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max
            REAL(KIND = 8) :: D1_min, D1_max, D2_min, D2_max, D3_min, D3_max
!
            REAL(KIND = 8) :: sphere_Rmax, sphere_Rmin, cylinder_Rmin, cylinder_Rmax, cylinder_radius1, cylinder_radius2
            REAL(KIND = 8) :: Base1_Long_Axis, Base2_Long_Axis, Base1_Short_Axis, Base2_Short_Axis, Long_Axis, Short_Axis, ellipse_long_axis_angle
!
            REAL(KIND = 8) :: X1_shift, Y1_shift, Z1_shift, R1_shift, X2_shift, Y2_shift, Z2_shift, R2_shift
!
            REAL(KIND = 8) :: exponent1, sign1, exponent2, sign2, emissivity
!
            REAL(KIND = 8) :: location_of_1st_periodic_occurrence, width_of_periodic_region, period_of_occurrence
            REAL(KIND = 8) :: period, L_first, L_width, thickness0, thickness1, thickness2, interpolation_search_radius
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), DIMENSION(3) :: sphere_center_coordinates, CylBase1_center_coordinates, CylBase2_center_coordinates
            REAL(KIND = 8), DIMENSION(3) :: Base1_center_coordinates,  Base2_center_coordinates
!
            REAL(KIND = 8), DIMENSION(0:10) :: BoundSurface1_EquCoeff_A, BoundSurface1_EquCoeff_B, BoundSurface1_EquCoeff_C, BoundSurface1_EquCoeff_D
            REAL(KIND = 8), DIMENSION(0:10) :: BoundSurface2_EquCoeff_A, BoundSurface2_EquCoeff_B, BoundSurface2_EquCoeff_C, BoundSurface2_EquCoeff_D
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: ier, n, i1, i2, j, k, m, AxisNum, n_table, mn, mk1, mk2, mk3
!
            INTEGER :: number_of_regions, NumTables_with_regions, TotNum_tabular_subdomains, NumRegions_this_table
            INTEGER :: number_of_periodic_regions, total_number_periodic_subdomains, number_of_periodic_occurrences
!
            INTEGER :: equation_order_of_bounding_surface1, equation_order_of_bounding_surface2, num_interrogated_points
!
            INTEGER :: number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2
!
            INTEGER, ALLOCATABLE, DIMENSION(:) :: frequency_of_occurrences
            INTEGER, ALLOCATABLE, DIMENSION(:) :: number_of_regions_in_table
!
! -------------
! ......... CHARACTER variables
! -------------
!
            CHARACTER(LEN =  1) :: axis_of_periodicity, vertical_or_stratigraphic1, vertical_or_stratigraphic2, reference_surface1, reference_surface2
            CHARACTER(LEN =  2) :: length_units, plane_of_ellipse_bases
            CHARACTER(LEN =  3) :: mn_num
            CHARACTER(LEN =  5) :: region_name
            CHARACTER(LEN =  6) :: first_part
            CHARACTER(LEN =  8) :: interpolation_data_file_name1, interpolation_data_file_name2, ref_interpolation_data_file
            CHARACTER(LEN = 12) :: region_shape, type_of_equation1, type_of_equation2
            CHARACTER(LEN = 20) :: table_name
            CHARACTER(LEN = 50) :: read_data_format
!
            CHARACTER(LEN = 1) :: dependent_variable_of_surfaces, direction
!
! -------------
! ......... LOGICAL variables
! -------------
!
            LOGICAL :: IntDat_match, read_data_by_row, read_regions_from_table, active_table
!
! -------
! ... Namelists
! -------
!
      NAMELIST/ Heterogeneous_Regions / number_of_regions, number_of_periodic_regions, total_number_periodic_subdomains, &
     &                                  NumTables_with_regions, TotNum_tabular_subdomains, dominant_medium, dominant_medium_emissivity
!
      NAMELIST/ HetRegion_GeneralInfo / region_name, region_shape, length_units, emissivity, read_regions_from_table, &
     &                                  table_name, NumRegions_this_table
!
      NAMELIST/ Rectangular_HetRegion / X_min, Y_min, Z_min, X_max, Y_max, Z_max
!
      NAMELIST/ Periodic_HetRegion    / X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max,        &
     &                                  number_of_periodic_occurrences, axis_of_periodicity,           &          ! Needed for periodic shapes
     &                                  location_of_1st_periodic_occurrence, width_of_periodic_region, &
     &                                  period_of_occurrence
!
      NAMELIST/ Cylindrical_HetRegion / CylBase1_center_coordinates, CylBase2_center_coordinates, &
     &                                  cylinder_Rmin, cylinder_Rmax, cylinder_radius1, cylinder_radius2
!
      NAMELIST/ Spherical_HetRegion   / sphere_center_coordinates, sphere_Rmax, sphere_Rmin
!
      NAMELIST/ Elliptical_HetRegion  / Base1_center_coordinates, Base2_center_coordinates, ellipse_long_axis_angle, plane_of_ellipse_bases, &
     &                                  Base1_Long_Axis, Base2_Long_Axis, Base1_Short_Axis, Base2_Short_Axis, Long_Axis, Short_Axis
!
      NAMELIST/ Irregular_HetRegion   /  dependent_variable_of_surfaces, type_of_equation1, type_of_equation2, &
     &                                   X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max,               &
     &                                   interpolation_data_file_name1, interpolation_data_file_name2,         &
     &                                   vertical_or_stratigraphic1, vertical_or_stratigraphic2, reference_surface1, reference_surface2, &
     &                                   thickness0, thickness1, thickness2, ref_interpolation_data_file
!
      NAMELIST/ Irregular_HetRegion_Surf1 / equation_order_of_bounding_surface1,                      &
     &                                      X1_shift, Y1_shift, Z1_shift, R1_shift, exponent1, sign1, &
     &                                      BoundSurface1_EquCoeff_A, BoundSurface1_EquCoeff_B, BoundSurface1_EquCoeff_C, BoundSurface1_EquCoeff_D
!
      NAMELIST/ Irregular_HetRegion_Surf2 / equation_order_of_bounding_surface2,                      &
     &                                      X2_shift, Y2_shift, Z2_shift, R2_shift, exponent2, sign2, &
     &                                      BoundSurface2_EquCoeff_A, BoundSurface2_EquCoeff_B, BoundSurface2_EquCoeff_C, BoundSurface2_EquCoeff_D
!
      NAMELIST/ Irregular_HetRegion_IntTable1 / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                          read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
      NAMELIST/ Irregular_HetRegion_IntTable2 / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                          read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
      NAMELIST/ Irregular_HetRegion_RefTable  / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                          read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Define_Heterogeneous_Regions>
!
!
            WRITE(*,6000)
!
! ......... Initialization
!
            NumTables_with_regions          = 0
            number_of_periodic_regions       = 0
            total_number_periodic_subdomains = 0
            TotNum_tabular_subdomains    = 0
!
            dominant_medium_emissivity = 0.0d0
!
! ------------
! ......... Read the number of heterogeneous domains
! ------------
!
            READ (UNIT = *, NML = Heterogeneous_Regions, IOSTAT = ier )
!
! ......... Stop if there is a problem reading the namelist
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6600) '<Heterogeneous_Regions>'
               STOP
            END IF
!
! ......... If a set of regions is read through a table/file ...
!
            IF( NumTables_with_regions > 0 ) THEN
!
               ALLOCATE( RegionTable_UnitNum(NumTables_with_regions) )
               ALLOCATE( number_of_regions_in_table(NumTables_with_regions) )
!
               DO i1 = 1,NumTables_with_regions
                  RegionTable_UnitNum = 200 + i1 - 1
               END DO
!
            END IF
!
            ALLOCATE( frequency_of_occurrences(number_of_periodic_regions) )
!
            frequency_of_occurrences = 0 ! ... Initialization
!
            Num_HetRegions    = number_of_regions
!
            TotNum_HetRegions =   number_of_regions &
     &                          - number_of_periodic_regions + total_number_periodic_subdomains &
     &                          - NumTables_with_regions     + TotNum_tabular_subdomains
!
! ......... Storing temporarily table names, before re-dimensioning the grid/interpolation data sets
!
            i1 = Num_IntTables
            i2 = Num_IntTables + 2 * TotNum_HetRegions
!
            IF( i1 > 0 ) temp_IntTable(1:i1) = IntTable(1:i1)
!
            IF( ALLOCATED(IntTable) ) DEALLOCATE( IntTable )
!
            ALLOCATE( IntTable(i2), STAT=ier )
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'IntTable'
            ELSE
               WRITE(*,6102) 'IntTable'
               STOP
            END IF
!
            IF( i1 > 0 ) IntTable(1:i1) = temp_IntTable
!
            IF( ALLOCATED(temp_IntTable) ) DEALLOCATE( temp_IntTable )
            ALLOCATE( temp_IntTable(i2), STAT=ier )
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'temp_IntTable'
            ELSE
               WRITE(*,6102) 'temp_IntTable'
               STOP
            END IF
!
! ------------
! ......... Allocate memory to the arrays defining the heterogeneous regions
! ------------
!
            ALLOCATE(Region(1:TotNum_HetRegions), MediaSequenceNumber(1:TotNum_HetRegions), STAT=ier)
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'Region, MediaSequenceNumber'
            ELSE
               WRITE(*,6102) 'Region, MediaSequenceNumber'
               STOP
            END IF
!
!***********************************************************************
!*                                                                     *
!*         Read the specifics of the remaining media regions           *
!*                                                                     *
!***********************************************************************
!
            k = 1
            m = 0   ! ... Counter of number of periodic regions
!
            n_table = 0   ! ... Counter of number of region sets read from tables
!
            DO_Regions: DO n = 2, Num_HetRegions + 1
!
! -----------------
! ............ Read heterogeneous region specifics for Num_HetRegions > 1
! -----------------
!
               IF( n == Num_HetRegions + 1 ) THEN
!
                  READ(UNIT = *, FMT = '(A3)') first_part
!
                  IF( first_part /= '<<<' ) THEN
                     WRITE(UNIT = *, FMT = 6515) first_part
                     STOP
                  ELSE
                     EXIT DO_Regions
                  END IF
!
               END IF
!
! ............ Advance the counter
!
               k = k + 1
!
! ............ Initializations - Namelist components
!
               region_shape = '           '
               region_name  = '     '
               length_units = 'm'
!
               emissivity = 0.0d0
!
               read_regions_from_table = .FALSE.
               active_table = .FALSE.

!
!  ........... Initialization
!
               Region(k)%name  = '     '
               Region(k)%shape = '           '
               Region(k)%units = 'm'
               Region(k)%emissivity = 0.0d0
!
! -----------------
! ............ Read general info on the heterogeneous region
! -----------------
!
               READ (UNIT = *, NML = HetRegion_GeneralInfo, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6601) '<HetRegion_GeneralInfo>', n
                  STOP
               END IF
!
               SELECT CASE( region_shape(1:1) )
               CASE('R', 'r', 'C', 'c', 'S', 's', 'E', 'e', 'P', 'p', 'I', 'i' )
                  CONTINUE
               CASE DEFAULT
                  WRITE( UNIT = *, FMT = 6010 ) n, region_shape
                  STOP
               END SELECT
!
               IF( ( coordinates(1:2) == 'CY' .OR. coordinates(1:4) == 'Cy' .OR. coordinates(1:4) == 'cy' ) .AND. &
     &             ( region_shape(1:1) == 'R' .OR. region_shape(1:1) == 'r' ) )  &
     &         THEN
                  WRITE(*,6520) region_shape
                  STOP
               END IF
!
! ............ Assignment of the namelist values
!
               Region(k)%name  = region_name
               Region(k)%shape = region_shape
               Region(k)%units = length_units
!
               Region(k)%emissivity = emissivity
!
               IF_Table1: IF( read_regions_from_table ) THEN
!
                  n_table = n_table + 1
!
                  OPEN(   UNIT = RegionTable_UnitNum(n_table), FILE = table_name, STATUS = 'OLD' )
                  REWIND( UNIT = RegionTable_UnitNum(n_table) )
!
                  mn = 0                          ! ... Counter of table entries
!
                  active_table = .TRUE.
!
                  number_of_regions_in_table(n_table) = NumRegions_this_table
!
               END IF IF_Table1
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    Reading the shape-specific data from the remaining namelists     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
               SELECT CASE(region_shape(1:1))
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Rectangular regions
! >>>>>>>>>>>>
!
               CASE( 'R', 'r' )
!
! ............... Initializations
!
  500             X_min = 0.0d0
                  Y_min = 0.0d0
                  Z_min = 0.0d0
                  X_max = 0.0d0
                  Y_max = 0.0d0
                  Z_max = 0.0d0
!
                  Region(k)%LMin = 0.0d0     ! ... Whole array operations
                  Region(k)%LMax = 0.0d0
!
! ............... Reading the rectangle-related data using a namelist
!
                  IF_Table2: IF( .NOT. read_regions_from_table ) THEN
!
                     READ (UNIT = *, NML = Rectangular_HetRegion, IOSTAT = ier )
!
                     IF(ier /= 0) THEN
                        WRITE(UNIT = *, FMT = 6601) '<Rectangular_HetRegion>', n
                        STOP
                     END IF
!
                  ELSE IF( read_regions_from_table .AND. active_table) THEN
!
                     mn = mn+1
!
! .................. Read 1st line of table file: determine the variable and invariable coordinates
!
                     IF_TEntry1: IF( mn == 1 ) THEN
!
                        READ (UNIT = RegionTable_UnitNum(n_table), FMT = *, IOSTAT = ier ) Direction, D2_min, D2_max, D3_min, D3_max
!
                        IF(ier /= 0) THEN
                           WRITE(UNIT = *, FMT = 6610) table_name, n
                           STOP
                        END IF
!
                        IF( direction == 'X' ) THEN
                           Y_min = D2_min
                           Y_max = D2_max
                           Z_min = D3_min
                           Z_max = D3_max
                        ELSE IF( direction == 'Y' ) THEN
                           X_min = D2_min
                           X_max = D2_max
                           Z_min = D3_min
                           Z_max = D3_max
                        ELSE IF( direction == 'Z' ) THEN
                           X_min = D2_min
                           X_max = D2_max
                           Y_min = D3_min
                           Y_max = D3_max
                        ELSE
                           WRITE( UNIT = *, FMT = 6050 ) direction, table_name
                           STOP
                        END IF
!
                     END IF IF_TEntry1
!
                     READ (UNIT = RegionTable_UnitNum(n_table), FMT = *, IOSTAT = ier ) D1_min, D1_max
!
                     IF(ier /= 0) THEN
                        WRITE(UNIT = *, FMT = 6612) table_name, n
                        STOP
                     END IF
!
                     IF( direction == 'X' ) THEN
                        X_min = D1_min
                        X_max = D1_max
                     ELSE IF( direction == 'Y' ) THEN
                        Y_min = D1_min
                        Y_max = D1_max
                     ELSE IF( direction == 'Z' ) THEN
                        Z_min = D1_min
                        Z_max = D1_max
                     END IF
!
                     mk1 = mn / 100
                     mn_num(1:1) = N_Character(mK1)
!
                     mk1 = MOD( mn, 100 )
                     IF( mk1 == 0 ) THEN
                        mn_num(2:2) = N_Character(0)
                        mn_num(3:3) = N_Character(0)
                     ELSE
                        mk2 = mk1 / 10
                        mk3 = MOD( mk1, 10 )
                        IF( mk3 == 0 ) THEN
                           mn_num(2:2) = N_Character(mk2)
                           mn_num(3:3) = N_Character(0)
                        ELSE
                           mn_num(2:2) = N_Character(mk2)
                           mn_num(3:3) = N_Character(mk3)
                        END IF
                     END IF
!
                  END IF IF_Table2
!
                  MediaSequenceNumber(k) = n
!
                  Region(k)%LMin(1) = X_min
                  Region(k)%LMin(2) = Y_min
                  Region(k)%LMin(3) = Z_min
!
                  Region(k)%LMax(1) = X_max
                  Region(k)%LMax(2) = Y_max
                  Region(k)%LMax(3) = Z_max
!
! ............... Ensuring that the LMin of the range is not >= of the LMax
!
                  IF( Region(k)%LMin(1) >= Region(k)%LMax(1) ) THEN
                     WRITE(*,6052) 'X',Region(k)%LMin(1), n, 'X',Region(k)%LMax(1)
                     STOP
                  END IF
!
                  IF( Region(k)%LMin(2) >= Region(k)%LMax(2) ) THEN
                     WRITE(*,6052) 'Y',Region(k)%LMin(2), n, 'Y',Region(k)%LMax(2)
                     STOP
                  END IF
!
                  IF( Region(k)%LMin(3) >= Region(k)%LMax(3) ) THEN
                     WRITE(*,6052)  'Z',Region(k)%LMin(3), n, 'Z',Region(k)%LMax(3)
                     STOP
                  END IF
!
! ............... If reading regions from table
!
                  IF( read_regions_from_table ) THEN
                     IF( mn < number_of_regions_in_table(n_table) ) THEN
                        k  = k+1
                        GO TO 500
                     ELSE
                        CLOSE( UNIT = RegionTable_UnitNum(n_table) )
                     END IF
                  END IF
!
                  CYCLE DO_Regions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Cylindrical regions
! >>>>>>>>>>>>
!
               CASE( 'C', 'c' )
!
! ............... Initializations
!
  600             CylBase1_center_coordinates = 0.0d0  ! ... Whole array operation
                  CylBase2_center_coordinates = 0.0d0
!
                  cylinder_Rmin    = 0.0d0
                  cylinder_Rmax    = 0.0d0
                  cylinder_radius1 = 0.0d0
                  cylinder_radius2 = 0.0d0
!
                  Region(k)%CylBase1Coord = 0.0d0    ! ... Whole array operation
                  Region(k)%CylBase2Coord = 0.0d0
!
                  Region(k)%CylRmin   = 0.0d0
                  Region(k)%CylRmax   = 0.0d0
                  Region(k)%CylBase1R = 0.0d0
                  Region(k)%CylBase2R = 0.0d0
!
! ............... Reading the cylinder-related data using a namelist
!
                  IF_Table3: IF( .NOT. read_regions_from_table ) THEN
!
                     READ (UNIT = *, NML = Cylindrical_HetRegion, IOSTAT = ier )
!
                     IF(ier /= 0) THEN
                        WRITE(UNIT = *, FMT = 6601) 'Cylindrical_HetRegion', n
                        STOP
                     END IF
!
                  ELSE IF( read_regions_from_table .AND. active_table) THEN
!
                     mn = mn+1
!
! ................... Read 1st line of table file: determine the variable and invariable coordinates
!
                     IF_TEntry2: IF( mn == 1 ) THEN
!
                        READ (UNIT = RegionTable_UnitNum(n_table), FMT = * ) direction, D2_min, D2_max, D3_min, D3_max
!
                        IF(ier /= 0) THEN
                           WRITE(UNIT = *, FMT = 6610) table_name, n
                           STOP
                        END IF
!
                     END IF IF_TEntry2
!
                     IF( direction == 'R' ) THEN
                        CylBase1_center_coordinates = (/ D3_min, D3_max, D2_min /)
                        CylBase2_center_coordinates = (/ D3_min, D3_max, D2_max /)
                     ELSE IF( direction == 'Z' ) THEN
                        Cylinder_Rmin = D2_min
                        Cylinder_Rmax = D2_max
                     ELSE
                        WRITE( UNIT = *, FMT = 6050 ) direction, table_name
                        STOP
                     END IF
!
                     READ (UNIT = RegionTable_UnitNum(n_table), FMT = * ) D1_min, D1_max
!
                     IF(ier /= 0) THEN
                        WRITE(UNIT = *, FMT = 6612) table_name, n
                        STOP
                     END IF
!
                     IF( direction == 'R' ) THEN
                        Cylinder_Rmin = D1_min
                        Cylinder_Rmax = D1_max
                     ELSE IF( direction == 'Z' ) THEN
                        CylBase1_center_coordinates = (/ D3_min, D3_max, D1_min /)
                        CylBase2_center_coordinates = (/ D3_min, D3_max, D1_max /)
                     END IF
!
                     mk1 = mn / 100
                     mn_num(1:1) = N_Character(mK1)
!
                     mk1 = MOD( mn, 100 )
                     IF( mk1 == 0 ) THEN
                        mn_num(2:2) = N_Character(0)
                        mn_num(3:3) = N_Character(0)
                     ELSE
                        mk2 = mk1 / 10
                        mk3 = MOD( mk1, 10 )
                        IF( mk3 == 0 ) THEN
                           mn_num(2:2) = N_Character(mk2)
                           mn_num(3:3) = N_Character(0)
                        ELSE
                           mn_num(2:2) = N_Character(mk2)
                           mn_num(3:3) = N_Character(mk3)
                        END IF
                     END IF
!
                     Region(k)%name  = region_name(1:2)//mn_num
                     Region(k)%shape = region_shape
                     Region(k)%units = length_units
                     Region(k)%emissivity = emissivity
!
                  END IF IF_Table3
!
                  MediaSequenceNumber(k) = n
!
! ............... For Cartesian coordinates
!
                  IF_Coord: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                          coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN
!
                     Region(k)%CylBase1Coord = CylBase1_center_coordinates
                     Region(k)%CylBase2Coord = CylBase2_center_coordinates
!
                     IF( ( cylinder_Rmin > 1.0d-7 .OR. cylinder_Rmax > 1.0d-7 ) .AND. ( cylinder_radius1 < 1.0d-7 .AND. cylinder_radius2 < 1.0d-7 ) ) THEN
                        Region(k)%CylRmin = cylinder_Rmin
                        Region(k)%CylRmax = cylinder_Rmax
                     ELSE IF( ( cylinder_Rmin < 1.0d-7 .AND. cylinder_Rmax < 1.0d-7 ) .AND. ( cylinder_radius1 > 1.0d-7 .OR. cylinder_radius2 > 1.0d-7 ) ) THEN
                        Region(k)%CylBase1R = cylinder_radius1
                        Region(k)%CylBase2R = cylinder_radius2
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'cylindrical',n
                        STOP
                     END IF
!
                     IF( cylinder_radius1 /= 0.0d0 ) CYCLE DO_Regions
!
! .................. Ensuring that the two bases of the cylindrical region do not coincide
!
                     D1 = Region(k)%CylBase1Coord(1) - Region(k)%CylBase2Coord(1)
                     D2 = Region(k)%CylBase1Coord(2) - Region(k)%CylBase2Coord(2)
                     D3 = Region(k)%CylBase1Coord(3) - Region(k)%CylBase2Coord(3)
!
                     IF( ABS(D1) < 1.0d-7 .AND. ABS(D2) < 1.0d-7 .AND. ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! .................. Ensuring that the Rmin of the cylindrical region is not >= of the Rmax
!
                     IF( Region(k)%CylRmin >= Region(k)%CylRmax ) THEN
                        WRITE(*,6054) Region(k)%CylRmin, n, Region(k)%CylRmax
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     Region(k)%CylBase1Coord(3) = CylBase1_center_coordinates(3)
                     Region(k)%CylBase2Coord(3) = CylBase2_center_coordinates(3)
!
                     IF( ( cylinder_Rmax > 1.0d-7 ) .AND. ( cylinder_radius1 < 1.0d-7 .AND. cylinder_radius2 < 1.0d-7 ) ) THEN
                        Region(k)%CylRmin = cylinder_Rmin
                        Region(k)%CylRmax = cylinder_Rmax
                     ELSE IF( ( cylinder_Rmax < 1.0d-7 ) .AND. ( cylinder_radius1 > 1.0d-7 .AND. cylinder_radius2 > 1.0d-7 ) ) THEN
                        Region(k)%CylBase1R = cylinder_radius1
                        Region(k)%CylBase2R = cylinder_radius2
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'cylindrical',n
                        STOP
                     END IF
!
                     IF( cylinder_radius1 /= 0.0d0 ) CYCLE DO_Regions
!
! .................. Ensuring that the two bases of the cylindrical region do not coincide
!
                     D3 = Region(k)%CylBase1Coord(3) - Region(k)%CylBase2Coord(3)
!
                     IF( ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! .................. Ensuring that the Rmin of the cylindrical region is not >= of the Rmax
!
                     IF( Region(k)%CylRmin >= Region(k)%CylRmax ) THEN
                       WRITE(*,6054) Region(k)%CylRmin, n, Region(k)%CylRmax
                        STOP
                     END IF
!
                  END IF IF_Coord
!
! ............... If reading regions from table
!
                  IF( read_regions_from_table ) THEN
                     IF( mn < number_of_regions_in_table(n_table) ) THEN
                        k  = k+1
                        GO TO 600
                     ELSE
                        CLOSE( UNIT = RegionTable_UnitNum(n_table) )
                     END IF
                  END IF
!
                  CYCLE DO_Regions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Periodic regions
! >>>>>>>>>>>>
!
               CASE( 'P', 'p' )
!
! ............... Initializations
!
                  X_min = 0.0d0
                  R_min = 0.0d0
                  Y_min = 0.0d0
                  Z_min = 0.0d0
                  X_max = 0.0d0
                  R_max = 0.0d0
                  Y_max = 0.0d0
                  Z_max = 0.0d0
!
                  region_shape = 'periodic'
!
                  axis_of_periodicity                 = ' '
                  number_of_periodic_occurrences      = 0
                  location_of_1st_periodic_occurrence = 0.0d0
                  width_of_periodic_region = 0.0d0
                  period_of_occurrence     = 0.0d0
!
! ............... Reading the periodic-related data using a namelist
!
                  IF( .NOT. read_regions_from_table ) THEN
!
                     READ (UNIT = *, NML = Periodic_HetRegion, IOSTAT = ier )
!
                     IF(ier /= 0) THEN
                        WRITE(UNIT = *, FMT = 6601) '<Periodic_HetRegion>', n
                        STOP
                     END IF
!
                  END IF
!
                  m = m + 1
!
                  IF ( axis_of_periodicity == 'X' .AND. axis_of_periodicity == 'x' .AND.  &
     &                 axis_of_periodicity == 'R' .AND. axis_of_periodicity == 'r' .AND.  &
     &                 axis_of_periodicity == 'Y' .AND. axis_of_periodicity == 'y' .AND.  &
     &                 axis_of_periodicity == 'Z' .AND. axis_of_periodicity == 'z' )   &
     &            THEN
                     WRITE( UNIT = *, FMT = 6702 ) n, axis_of_periodicity
                     STOP
                  END IF
!
                  L_first = location_of_1st_periodic_occurrence
                  L_width = width_of_periodic_region
                  period  = period_of_occurrence
!
                  frequency_of_occurrences(m) = number_of_periodic_occurrences
!
                  IF_Period: IF (axis_of_periodicity == 'X' .OR. axis_of_periodicity == 'x' ) THEN
!
                     AxisNum = 1
                     IF( coordinates(1:4) == 'CYLI' .OR. coordinates(1:4) == 'Cyli' .OR. coordinates(1:4) == 'cyli' ) THEN
                        WRITE( UNIT = *, FMT = 6701 ) n, coordinates, axis_of_periodicity
                        STOP
                     END IF
!
                  ELSE IF (axis_of_periodicity == 'R' .OR. axis_of_periodicity == 'r' ) THEN
!
                     AxisNum = 1
                     IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                   coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &               THEN
                        WRITE( UNIT = *, FMT = 6701 ) n, coordinates, axis_of_periodicity
                        STOP
                     END IF
!
                  ELSE IF (axis_of_periodicity == 'Y' .OR. axis_of_periodicity == 'y' ) THEN
!
                     AxisNum = 2
!
                  ELSE IF (axis_of_periodicity == 'Z' .OR. axis_of_periodicity == 'z' ) THEN
!
                     AxisNum = 3
!
                  END IF IF_Period
!
                  j = 0
!
 1000             IF( AxisNum == 1 ) THEN
                     X_min = L_first + j * period
                     X_max = X_min + L_width
                     R_min = X_min
                     R_max = X_max
                  ELSE IF( AxisNum == 2 ) THEN
                     Y_min = L_first + j * period
                     Y_max = Y_min + L_width
                  ELSE
                     Z_min = L_first + j * period
                     Z_max = Z_min + L_width
                  END IF
!
                  MediaSequenceNumber(k) = n
!
                  Region(k)%name  = region_name
                  Region(k)%shape = region_shape
                  Region(k)%units = length_units
!
                  Region(k)%LMin = 0.0d0     ! ... Whole array operations
                  Region(k)%LMax = 0.0d0
!
                  Region(k)%LMin(1) = X_min
                  Region(k)%LMin(2) = Y_min
                  Region(k)%LMin(3) = Z_min
!
                  Region(k)%LMax(1) = X_max
                  Region(k)%LMax(2) = Y_max
                  Region(k)%LMax(3) = Z_max
!
! ............... Ensuring that the LMin of the range is not >= of the LMax
!
                  IF( Region(k)%LMin(1) >= Region(k)%LMax(1) ) THEN
                     WRITE(*,6052) 'X',Region(k)%LMin(1), n, 'X',Region(k)%LMax(1)
                     STOP
                  END IF
!
                  IF( Region(k)%LMin(2) >= Region(k)%LMax(2) ) THEN
                     WRITE(*,6052) 'Y',Region(k)%LMin(2), n, 'Y',Region(k)%LMax(2)
                     STOP
                  END IF
!
                  IF( Region(k)%LMin(3) >= Region(k)%LMax(3) ) THEN
                     WRITE(*,6052) 'Z',Region(k)%LMin(3), n, 'Z',Region(k)%LMax(3)
                     STOP
                  END IF
!
                  IF ( j < frequency_of_occurrences(m) - 1 ) THEN
                     j = j + 1
                     k = k + 1
                     GO TO 1000
                  END IF
!
                  CYCLE DO_Regions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Elliptical regions
! >>>>>>>>>>>>
!
               CASE( 'E', 'e' )
!
! ............... Initializations
!
                  Base1_center_coordinates = 0.0d0  ! ... Whole array operation
                  Base2_center_coordinates = 0.0d0
                  ellipse_long_axis_angle  = 0.0d0
!
                  plane_of_ellipse_bases = '  '
!
                  Base1_Long_Axis  = 0.0d0
                  Base2_Long_Axis  = 0.0d0
                  Base1_Short_Axis = 0.0d0
                  Base2_Short_Axis = 0.0d0
                  Long_Axis  = 0.0d0
                  Short_Axis = 0.0d0
!
                  Region(k)%EllBase1Coord = 0.0d0    ! ... Whole array operation
                  Region(k)%EllBase1Coord = 0.0d0
                  Region(k)%EllipsePlane  = '  '
!
                  Region(k)%LAxisAngle = 0.0d0
                  Region(k)%Base1LAxis = 0.0d0
                  Region(k)%Base2LAxis = 0.0d0
                  Region(k)%Base1SAxis = 0.0d0
                  Region(k)%Base2SAxis = 0.0d0
                  Region(k)%LAxis = 0.0d0
                  Region(k)%SAxis = 0.0d0
!
! ............... Reading the cylinder-related data using a namelist
!
                  IF( .NOT. read_regions_from_table ) THEN
!
                     READ (UNIT = *, NML = Elliptical_HetRegion, IOSTAT = ier )
!
                     IF(ier /= 0) THEN
                        WRITE(UNIT = *, FMT = 6601) 'Elliptical_HetRegion', n
                        STOP
                     END IF
!
                  END IF
!
! ............... Ensuring that there is no conflict between the long and short axes lengths
!
                  IF( Base1_Short_Axis > Base1_Long_Axis ) THEN
                     WRITE(*,6056) Base1_Short_Axis, n, Base1_Long_Axis
                     STOP
                  END IF
!
                  IF( Base2_Short_Axis > Base2_Long_Axis ) THEN
                     WRITE(*,6058) Base2_Short_Axis, n, Base2_Long_Axis
                     STOP
                  END IF
!
                  SELECT CASE(plane_of_ellipse_bases)
                  CASE('XY','Xy','xY','xy','XZ','Xz','xZ','xz','YZ','Yz','yZ','yz')
                     CONTINUE
                  CASE DEFAULT
                     WRITE(UNIT = *, FMT = 6059) plane_of_ellipse_bases, n
                     STOP
                  END SELECT
!
                  MediaSequenceNumber(k) = n
!
! ............... For Cartesian coordinates
!
                  IF_Coord2: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                           coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     Region(k)%EllBase1Coord = Base1_center_coordinates
                     Region(k)%EllBase2Coord = Base2_center_coordinates
                     Region(k)%LAxisAngle    = PI * ellipse_long_axis_angle / 1.8d2
                     Region(k)%EllipsePlane  = plane_of_ellipse_bases
!
                     IF( ( Long_Axis > 1.0d-7 .AND. Short_Axis > 1.0d-7 ) .AND. &
     &                   ( Base1_Long_Axis < 1.0d-7 .AND. Base1_Long_Axis < 1.0d-7 .AND. Base1_Short_Axis < 1.0d-7 .AND. Base2_Short_Axis < 1.0d-7 ) ) &
     &               THEN
                        Region(k)%LAxis = Long_Axis
                        Region(k)%SAxis = Short_Axis
                     ELSE IF( ( Long_Axis < 1.0d-7 .AND. Short_Axis < 1.0d-7 ) .AND. &
     &                        ( Base1_Long_Axis > 1.0d-7 .AND. Base1_Long_Axis > 1.0d-7 .AND. Base1_Short_Axis > 1.0d-7 .AND. Base2_Short_Axis > 1.0d-7 ) ) &
     &               THEN
                        Region(k)%Base1LAxis = Base1_Long_Axis
                        Region(k)%Base2LAxis = Base2_Long_Axis
                        Region(k)%Base1SAxis = Base1_Short_Axis
                        Region(k)%Base2SAxis = Base2_Short_Axis
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'elliptical',n
                        STOP
                     END IF
!
! .................. Ensuring that the two bases of the elliptical region do not coincide
!
                     D1 = Region(k)%EllBase1Coord(1) - Region(k)%EllBase2Coord(1)
                     D2 = Region(k)%EllBase1Coord(2) - Region(k)%EllBase2Coord(2)
                     D3 = Region(k)%EllBase1Coord(3) - Region(k)%EllBase2Coord(3)
!
                     IF( ABS(D1) < 1.0d-7 .AND. ABS(D2) < 1.0d-7 .AND. ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     Region(k)%EllBase1Coord(3) = Base1_center_coordinates(3)
                     Region(k)%EllBase2Coord(3) = Base2_center_coordinates(3)
                     Region(k)%EllipsePlane     = '  '
                     Region(k)%LAxisAngle       = 0.0d0
!
                     IF( ( Long_Axis > 1.0d-7 .AND. Short_Axis > 1.0d-7 ) .AND. &
     &                   ( Base1_Long_Axis < 1.0d-7 .AND. Base1_Long_Axis < 1.0d-7 .AND. Base1_Short_Axis < 1.0d-7 .AND. Base2_Short_Axis < 1.0d-7 ) ) &
     &               THEN
                        Region(k)%LAxis = Long_Axis
                        Region(k)%SAxis = Short_Axis
                     ELSE IF( ( Long_Axis < 1.0d-7 .AND. Short_Axis < 1.0d-7 ) .AND. &
     &                        ( Base1_Long_Axis > 1.0d-7 .AND. Base1_Long_Axis > 1.0d-7 .AND. Base1_Short_Axis > 1.0d-7 .AND. Base2_Short_Axis > 1.0d-7 ) ) &
     &               THEN
                        Region(k)%Base1LAxis = Base1_Long_Axis
                        Region(k)%Base2LAxis = Base2_Long_Axis
                        Region(k)%Base1SAxis = Base1_Short_Axis
                        Region(k)%Base2SAxis = Base2_Short_Axis
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'elliptical',n
                        STOP
                     END IF
!
! .................. Ensuring that the two bases of the cylindrical region do not coincide
!
                     D3 = Region(k)%CylBase1Coord(3) - Region(k)%CylBase2Coord(3)
!
                     IF( ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
                  END IF IF_Coord2
!
                  CYCLE DO_Regions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Spherical regions
! >>>>>>>>>>>>
!
               CASE( 'S', 's' )
!
! ............... Initializations
!
                  sphere_center_coordinates = 0.0d0    ! ... Whole array operation
                  sphere_Rmax = 0.0d0
                  sphere_Rmin = 0.0d0
!
                  Region(k)%SphereCenterCoord = 0.0d0  ! ... Whole array operation
                  Region(k)%SphRmin = 0.0d0
                  Region(k)%SphRmax = 0.0d0
!
! ............... Reading the sphere-related data using a namelist
!
                  IF( .NOT. read_regions_from_table ) THEN
!
                     READ (UNIT = *, NML = Spherical_HetRegion, IOSTAT = ier )
!
                     IF(ier /= 0) THEN
                        WRITE(UNIT = *, FMT = 6601) 'Spherical_HetRegion', n
                        STOP
                     END IF
!
                  END IF
!
                  MediaSequenceNumber(k) = n
!
                  Region(k)%SphRmin = sphere_Rmin
                  Region(k)%SphRmax = sphere_Rmax
!
! ............... For Cartesian coordinates
!
                  IF_Coord3: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                           coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     Region(k)%SphereCenterCoord = sphere_center_coordinates
!
! .................. Ensuring that the Rmin of the spherical region is not >= of the Rmax
!
                     IF( Region(k)%SphRmin >= Region(k)%SphRmax ) THEN
                        WRITE(*,6055) Region(k)%SphRmin, n, Region(k)%CylRmax
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     Region(k)%SphereCenterCoord(3) = sphere_center_coordinates(3)
!
! .................. Ensuring that the Rmin of the spherical region is not >= of the Rmax
!
                     IF( Region(k)%SphRmin >= Region(k)%SphRmax ) THEN
                        WRITE(*,6055) Region(k)%SphRmin, n, Region(k)%CylRmax
                        STOP
                     END IF
!
                  END IF IF_Coord3
!
                  CYCLE DO_Regions
!
               END SELECT
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Irregular regions
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               type_of_equation1 = '           '  ! ... Initialization
               type_of_equation2 = '           '
!
               X_min = 0.0d0
               Y_min = 0.0d0
               Z_min = 0.0d0
               X_max = 0.0d0
               Y_max = 0.0d0
               Z_max = 0.0d0
!
               X1_shift = 0.0d0
               Y1_shift = 0.0d0
               Z1_shift = 0.0d0
!
               X2_shift = 0.0d0
               Y2_shift = 0.0d0
               Z2_shift = 0.0d0
!
               R_min = 0.0d0
               R_max = 0.0d0
!
               R1_shift = 0.0d0
               R2_shift = 0.0d0
!
               dependent_variable_of_surfaces = ' '
               equation_order_of_bounding_surface1 = 0
               equation_order_of_bounding_surface2 = 0
!
               interpolation_data_file_name1 = '        '
               interpolation_data_file_name2 = '        '
               ref_interpolation_data_file   = '        '
!
               sign1 = 1.0d0
               sign2 = 1.0d0
               exponent1 = 0.0d0
               exponent2 = 0.0d0
!
               thickness0 = 0.0d0
               thickness1 = 0.0d0
               thickness2 = 0.0d0
               interpolation_search_radius = 0.0d0
!
               vertical_or_stratigraphic1 = ' '
               vertical_or_stratigraphic2 = ' '
!
               reference_surface1  = ' '
               reference_surface2  = ' '
!
! ............ Initializations - Array components
!
               Region(k)%TypeEqu1 = '           '
               Region(k)%TypeEqu2 = '           '
!
               Region(k)%DependentVar = ' '
!
               Region(k)%OrderEquSurf1 = 0
               Region(k)%OrderEquSurf2 = 0
!
               Region(k)%LMin = 0.0d0     ! ... Whole array operations
               Region(k)%LMax = 0.0d0
!
               Region(k)%L1Shift = 0.0d0
               Region(k)%L2Shift = 0.0d0
!
               Region(k)%sign1  = 1.0d0
               Region(k)%sign2  = 1.0d0
               Region(k)%expon1 = 0.0d0
               Region(k)%expon2 = 0.0d0
!
               Region(k)%IntTableNum1 = 0
               Region(k)%IntTableNum2 = 0
!
               Region(k)%thick0 = 0.0d0
               Region(k)%thick1 = 0.0d0
               Region(k)%thick2 = 0.0d0
!
               Region(k)%VertOrStrat1 = ' '
               Region(k)%VertOrStrat2 = ' '
               Region(k)%RefSurf1     = ' '
               Region(k)%RefSurf2     = ' '
!
! ............ Reading the irregular zone data using a namelist
!
               IF( .NOT. read_regions_from_table ) THEN
!
                  READ (UNIT = *, NML = Irregular_HetRegion, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_HetRegion>', n
                     STOP
                  END IF
!
               END IF
!
! ............ Checking the dependent variable
!
               IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &             coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &         THEN

!
                  IF( ( dependent_variable_of_surfaces /= 'X' .AND. dependent_variable_of_surfaces /= 'x'  .AND.  &
     &                  dependent_variable_of_surfaces /= 'Y' .AND. dependent_variable_of_surfaces /= 'y'  .AND.  &
     &                  dependent_variable_of_surfaces /= 'Z' .AND. dependent_variable_of_surfaces /= 'z' ) )     &
     &            THEN
                     WRITE( UNIT = *, FMT = 6608 ) 'cartesian', n, 'X, Y and Z'
                     STOP
                  END IF
!
               ELSE
!
                  IF( ( dependent_variable_of_surfaces /= 'R' .AND. dependent_variable_of_surfaces /= 'r'   .AND. &
     &                  dependent_variable_of_surfaces /= 'Z' .AND. dependent_variable_of_surfaces /= 'z' ) )     &
     &            THEN
                     WRITE( UNIT = *, FMT = 6608 ) 'cylindrical', n, 'R and Z'
                     STOP
                  END IF
!
               END IF
!
! ............ Assignment
!
               Region(k)%LMin(1) = X_min
               Region(k)%LMin(2) = Y_min
               Region(k)%LMin(3) = Z_min
!
               Region(k)%LMax(1) = X_max
               Region(k)%LMax(2) = Y_max
               Region(k)%LMax(3) = Z_max
!
               Region(k)%DependentVar = dependent_variable_of_surfaces
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Setting the type of equation of surface 1
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( type_of_equation1(1:4) == 'POLY' .OR. type_of_equation1(1:4) == 'Poly' .OR. type_of_equation1(1:4) == 'poly' ) THEN
!
! ............ Polynomial equation
!
                  type_of_equation1(1:4) = 'Poly'
                  IF( type_of_equation1(5:6) == '/E' .OR. type_of_equation1(5:6) == '/e' ) THEN
                     type_of_equation1(5:6) = '/E'
                  ELSE IF( type_of_equation1(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                     STOP
                  END IF
!
! ............ Power equation
!
               ELSE IF( type_of_equation1(1:4) == 'POWE' .OR. type_of_equation1(1:4) == 'Powe' .OR. type_of_equation1(1:4) == 'powe' ) THEN
!
                  type_of_equation1(1:4) = 'Powe'
                  IF( type_of_equation1(5:6) == '/E' .OR. type_of_equation1(5:6) == '/e' ) THEN
                     type_of_equation1(5:6) = '/E'
                  ELSE IF( type_of_equation1(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                     STOP
                  END IF
!
! ............ No equation: interpolation from a grid/tabular data set
!
               ELSE IF( type_of_equation1(1:4) == 'INTE' .OR. type_of_equation1(1:4) == 'Inte' .OR. type_of_equation1(1:4) == 'inte' ) THEN
!
                  CONTINUE
!
! ............ For Fixed distance surfaces
!
               ELSE IF( type_of_equation1(1:4) == 'FIXE' .OR. type_of_equation1(1:4) == 'Fixe' .OR. type_of_equation1(1:4) == 'fixe' ) THEN
!
! ............... Making sure that only one Fixed Width surface is defined
!
                  IF( ( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) .AND. &
     &                ( ref_interpolation_data_file(1:3) == '   ' ) )  &
     &            THEN
                     WRITE( UNIT = *, FMT = 6402 ) n
                     STOP
                  END IF
!
                  CONTINUE
!
! ............ Unknown option: ERROR!
!
               ELSE
                  WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                  STOP
               END IF
!
               Region(k)%TypeEqu1 = type_of_equation1
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Setting the type of equation of surface 2
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( type_of_equation2(1:4) == 'POLY' .OR. type_of_equation2(1:4) == 'Poly' .OR. type_of_equation2(1:4) == 'poly' ) THEN
!
! ............ Polynomial equation
!
                  type_of_equation2(1:4) = 'Poly'
                  IF( type_of_equation2(5:6) == '/E' .OR. type_of_equation2(5:6) == '/e' ) THEN
                     type_of_equation2(5:6) = '/E'
                  ELSE IF( type_of_equation2(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                     STOP
                  END IF
!
! ............ Power equation
!
               ELSE IF( type_of_equation2(1:4) == 'POWE' .OR. type_of_equation2(1:4) == 'Powe' .OR. type_of_equation2(1:4) == 'powe' ) THEN
                  type_of_equation2(1:4) = 'Powe'
                  IF( type_of_equation2(5:6) == '/E' .OR. type_of_equation2(5:6) == '/e' ) THEN
                     type_of_equation2(5:6) = '/E'
                  ELSE IF( type_of_equation2(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                     STOP
                  END IF
!
! ............ No equation: interpolation from a grid/tabular data set
!
               ELSE IF( type_of_equation2(1:4) == 'INTE' .OR. type_of_equation2(1:4) == 'Inte' .OR. type_of_equation2(1:4) == 'inte' ) THEN
!
                  CONTINUE
!
! ............ No equation: Fixed distance from another surface
!
               ELSE IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) THEN
!
                  CONTINUE
!
! ............ Unknown option: ERROR!
!
               ELSE
                  WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                  STOP
               END IF
!
               Region(k)%TypeEqu2 = type_of_equation2
!
               MediaSequenceNumber(k) = n
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For Cartesian coodinates
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF_Coord4: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                        coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &         THEN

!
                  IF( ( NXMax >  1 .AND. NYMax == 1 .AND. NZMax == 1 ) .OR. ( NXMax == 1 .AND. NYMax >  1 .AND. NZMax == 1 ) .OR.  &
     &                ( NXMax == 1 .AND. NYMax == 1 .AND. NZMax >  1 ) )    &
     &            THEN
                     WRITE( UNIT = *, FMT = 6606 ) 'cartesian', n
                     STOP
                  END IF
!
                  IF_Inte: IF( type_of_equation1(1:4) == 'INTE' .OR. type_of_equation1(1:4) == 'Inte' .OR. type_of_equation1(1:4) == 'inte' ) THEN
!
                     type_of_equation1(1:4) = 'Inte'
!
                     IntDat_match = .FALSE.
                     IF( Num_IntTables > 0 ) THEN
                        DO_m1: DO m = 1,Num_IntTables
                           IF( interpolation_data_file_name1 == IntTable(m)%FileName ) THEN
                              Region(k)%IntTableNum1 = m
                              IntDat_match           = .TRUE.
                              EXIT
                           END IF
                        END DO DO_m1
                     END IF
!
                     IF( .NOT. IntDat_match ) THEN
!
                        m = Num_IntTables + 1
                        Num_IntTables          = m
                        Region(k)%IntTableNum1 = Num_IntTables
!
! ..................... Reading info describing the structure of interpolation Table 1
!
                        READ (UNIT = *, NML = Irregular_HetRegion_IntTable1, IOSTAT = ier )
!
                        IF(ier /= 0) THEN
                          WRITE(UNIT = *, FMT = 6601) '<Irregular_HetRegion_IntTable1>', n
                           STOP
                        END IF
!
                        CALL Read_Tabular_Data( table_number = m, file_name = interpolation_data_file_name1, number_of_rows = number_of_rows, &
     &                                          number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,                  &
     &                                          read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                          RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                          interpolation_search_radius = interpolation_search_radius, &
     &                                          num_interrogated_points     = num_interrogated_points )
!
                     END IF
!
                     GO TO 1500
!
                  ELSE IF( type_of_equation1(1:4) == 'FIXE' .OR. type_of_equation1(1:4) == 'Fixe' .OR. type_of_equation1(1:4) == 'fixe' )  THEN
!
                     Region(k)%thick0 = thickness0
                     Region(k)%thick1 = thickness1
!
                     Region(k)%VertOrStrat1 = vertical_or_stratigraphic1
                     Region(k)%RefSurf1     = reference_surface1
!
                     IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) THEN
!
                        type_of_equation1(1:4) = 'Fixe'
                        type_of_equation2(1:4) = 'Fixe'
!
                        IntDat_match = .FALSE.
!
                        IF( Num_IntTables > 0 ) THEN
                           DO_m2: DO m = 1,Num_IntTables
                              IF( ref_interpolation_data_file == IntTable(m)%FileName ) THEN
                                 Region(k)%IntTableNum1 = m
                                 IntDat_match           = .TRUE.
                                 EXIT
                              END IF
                           END DO DO_m2
                        END IF
!
                        IF( .NOT. IntDat_match ) THEN
!
                           m = Num_IntTables + 1
                           Num_IntTables          = m
                           Region(k)%IntTableNum1 = Num_IntTables
!
! ........................ Reading info describing the structure of the reference
!
                           READ (UNIT = *, NML = Irregular_HetRegion_RefTable, IOSTAT = ier )
!
                           IF(ier /= 0) THEN
                              WRITE(UNIT = *, FMT = 6601) '<Irregular_HetRegion_RefTable>', n
                              STOP
                           END IF
!
                           CALL Read_Tabular_Data( table_number = m, file_name = ref_interpolation_data_file, number_of_rows = number_of_rows, &
     &                                             number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                             read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                             RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                             interpolation_search_radius = interpolation_search_radius,  &
     &                                             num_interrogated_points     = num_interrogated_points )
!
                        END IF
!
                     END IF
!
                     GO TO 1500
!
                  END IF IF_Inte
!
! ............... Initializations
!
                  BoundSurface1_EquCoeff_A = 0.0d0
                  BoundSurface1_EquCoeff_B = 0.0d0
                  BoundSurface1_EquCoeff_C = 0.0d0
                  BoundSurface1_EquCoeff_D = 0.0d0
!
! ............... Reading the data for surface 1
!
                  READ (UNIT = *, NML = Irregular_HetRegion_Surf1, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_HetRegion_Surf1>',n
                     STOP
                  END IF
!
                  Region(k)%OrderEquSurf1 = equation_order_of_bounding_surface1
!
                  i1 = equation_order_of_bounding_surface1
!
                  ALLOCATE ( Region(k)%Equ1CoeffA(0:i1), Region(k)%Equ1CoeffB(0:i1) )
!
                  Region(k)%Equ1CoeffA = 0.0d0
                  Region(k)%Equ1CoeffB = 0.0d0
!
                  IF_Surf1: IF( Region(k)%TypeEqu1(1:4) == 'Powe' ) THEN
                     ALLOCATE ( Region(k)%Equ1CoeffC(0:i1), Region(k)%Equ1CoeffD(0:i1) )
                     Region(k)%Equ1CoeffC = 0.0d0
                     Region(k)%Equ1CoeffD = 0.0d0
                  ELSE
                     IF(i1 > 1 ) THEN
                        ALLOCATE ( Region(k)%Equ1CoeffC(1:Max(1,i1)) )
                        Region(k)%Equ1CoeffC = 0.0d0
                     END IF
                  END IF IF_Surf1
!
! ............... Assignments
!
                  Region(k)%sign1  = sign1
                  Region(k)%expon1 = exponent1
!
                  Region(k)%Equ1CoeffA(0:i1) = BoundSurface1_EquCoeff_A(0:i1)
                  Region(k)%Equ1CoeffB(0:i1) = BoundSurface1_EquCoeff_B(0:i1)
!
                  IF( Region(k)%TypeEqu1(1:4) == 'Powe' ) THEN
                     Region(k)%Equ1CoeffC(0:i1) = BoundSurface1_EquCoeff_C(0:i1)
                     Region(k)%Equ1CoeffD(0:i1) = BoundSurface1_EquCoeff_D(0:i1)
                  ELSE
                     IF(i1 > 1) Region(k)%Equ1CoeffC(1:i1-1) = BoundSurface1_EquCoeff_C(1:i1-1)
                  END IF
!
                  Region(k)%L1Shift(1) = X1_shift
                  Region(k)%L1Shift(2) = Y1_shift
                  Region(k)%L1Shift(3) = Z1_shift
!
! >>>>>
!
 1500             IF_Inte2: IF( type_of_equation2(1:4) == 'INTE' .OR. type_of_equation2(1:4) == 'Inte' .OR. type_of_equation2(1:4) == 'inte' ) THEN
!
                     type_of_equation2(1:4) = 'Inte'
!
                     IntDat_match = .FALSE.
                     IF( Num_IntTables > 0 ) THEN
                         DO_m3: DO m = 1,Num_IntTables
                            IF( interpolation_data_file_name2 == IntTable(m)%FileName ) THEN
                               Region(k)%IntTableNum2 = m
                               IntDat_match           = .TRUE.
                               EXIT
                            END IF
                         END DO DO_m3
                     END IF
!
                     IF( .NOT. IntDat_match ) THEN
!
                        m = Num_IntTables + 1
                        Num_IntTables = m
                        Region(k)%IntTableNum2 = Num_IntTables
!
! ..................... Reading info describing the structure of interpolation Table 2
!
                        READ (UNIT = *, NML = Irregular_HetRegion_IntTable2, IOSTAT = ier )
!
                        IF(ier /= 0) THEN
                          WRITE(UNIT = *, FMT = 6601) '<Irregular_HetRegion_IntTable2>', n
                           STOP
                        END IF
!
                        CALL Read_Tabular_Data( table_number = m, file_name = interpolation_data_file_name2, number_of_rows = number_of_rows, &
     &                                          number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                          read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                          RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                          interpolation_search_radius = interpolation_search_radius,  &
     &                                          num_interrogated_points     = num_interrogated_points )
!
                     END IF
!
                     CYCLE DO_Regions
!
                  ELSE IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' )  THEN
!
                     Region(k)%thick2 = thickness2
!
                     Region(k)%VertOrStrat2 = vertical_or_stratigraphic2
                     Region(k)%RefSurf2     = reference_surface2
!
                     CYCLE DO_Regions
!
                  END IF IF_Inte2
!
! ............... Initializations
!
                  BoundSurface2_EquCoeff_A = 0.0d0
                  BoundSurface2_EquCoeff_B = 0.0d0
                  BoundSurface2_EquCoeff_C = 0.0d0
                  BoundSurface2_EquCoeff_D = 0.0d0
!
! ............... Reading the data for surface 2
!
                  READ (UNIT = *, NML = Irregular_HetRegion_Surf2, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_HetRegion_Surf2>',n
                     STOP
                  END IF
!
! ............... Checking consistency of orders of equation
!
                  IF( ( equation_order_of_bounding_surface1 < 0 ) .OR. ( equation_order_of_bounding_surface2 < 0 ) ) THEN
                     WRITE( UNIT = *, FMT = 6604 ) n
                     STOP
                  END IF
!
                  Region(k)%OrderEquSurf2 = equation_order_of_bounding_surface2
!
                  i2 = equation_order_of_bounding_surface2
!
                  ALLOCATE ( Region(k)%Equ2CoeffA(0:i2), Region(k)%Equ2CoeffB(0:i2) )
!
                  Region(k)%Equ2CoeffA = 0.0d0
                  Region(k)%Equ2CoeffB = 0.0d0
!
                  IF_Surf2: IF( Region(k)%TypeEqu2(1:4) == 'Powe' ) THEN
                     ALLOCATE ( Region(k)%Equ2CoeffC(0:i2), Region(k)%Equ2CoeffD(0:i2) )
                     Region(k)%Equ2CoeffC = 0.0d0
                     Region(k)%Equ2CoeffD = 0.0d0
                  ELSE
                     IF(i2 > 1 ) THEN
                        ALLOCATE ( Region(k)%Equ2CoeffC(1:Max(1,i2)) )
                        Region(k)%Equ2CoeffC = 0.0d0
                     END IF
                  END IF IF_Surf2
!
! ............... Assignments
!
                  Region(k)%sign2  = sign2
                  Region(k)%expon2 = exponent2
!
                  Region(k)%Equ2CoeffA(0:i2) = BoundSurface2_EquCoeff_A(0:i2)
                  Region(k)%Equ2CoeffB(0:i2) = BoundSurface2_EquCoeff_B(0:i2)
!
                  IF( Region(k)%TypeEqu2(1:4) == 'Powe' ) THEN
                     Region(k)%Equ2CoeffC(0:i2) = BoundSurface2_EquCoeff_C(0:i2)
                     Region(k)%Equ2CoeffD(0:i2) = BoundSurface2_EquCoeff_D(0:i2)
                  ELSE
                     IF(i2 > 1) Region(k)%Equ2CoeffC(1:i2-1) = BoundSurface2_EquCoeff_C(1:i2-1)
                  END IF
!
                  Region(k)%L2Shift(1) = X2_shift
                  Region(k)%L2Shift(2) = Y2_shift
                  Region(k)%L2Shift(3) = Z2_shift
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For Cylindrical coodinates
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               ELSE
!
! ............... Irregular shapes
!
                  IF( NXMax >  1 .AND. NZMax == 1 ) THEN
                     WRITE( UNIT = *, FMT = 6606 ) 'radial', n
                     STOP
                  END IF
!
! ............... Initializations
!
                  BoundSurface1_EquCoeff_A = 0.0d0
                  BoundSurface1_EquCoeff_B = 0.0d0
                  BoundSurface2_EquCoeff_A = 0.0d0
                  BoundSurface2_EquCoeff_B = 0.0d0
!
! ............... Reading the data for surface 1
!
                  READ (UNIT = *, NML = Irregular_HetRegion_Surf1, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_HetRegion_Surf1>', n
                     STOP
                  END IF
!
! ............... Reading the data for surface 2
!
                  READ (UNIT = *, NML = Irregular_HetRegion_Surf2, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_HetRegion_Surf2>', n
                     STOP
                  END IF
!
! ............... Checking consistency of orders of equation
!
                  IF( ( equation_order_of_bounding_surface1 < 0 ) .OR. ( equation_order_of_bounding_surface2 < 0 ) ) THEN
                     WRITE( UNIT = *, FMT = 6604 ) n
                     STOP
                  END IF
!
                  Region(k)%OrderEquSurf1 = equation_order_of_bounding_surface1
                  Region(k)%OrderEquSurf2 = equation_order_of_bounding_surface2
!
                  i1 = equation_order_of_bounding_surface1
                  ALLOCATE ( Region(k)%Equ1CoeffA(0:i1) )
!
                  Region(k)%Equ1CoeffA = 0.0d0
!
                  IF( Region(k)%TypeEqu1(1:4) == 'Powe' ) THEN
                     ALLOCATE ( Region(k)%Equ1CoeffB(0:i1) )
                     Region(k)%Equ1CoeffB = 0.0d0
                  END IF
!
                  i2 = equation_order_of_bounding_surface2
                  ALLOCATE ( Region(k)%Equ2CoeffA(0:i2) )
!
                  Region(k)%Equ2CoeffA = 0.0d0
!
                  IF( Region(k)%TypeEqu2(1:4) == 'Powe' ) THEN
                     ALLOCATE ( Region(k)%Equ2CoeffB(0:i2) )
                     Region(k)%Equ2CoeffB = 0.0d0
                  END IF
!
! ............... Assignments
!
                  Region(k)%LMin(1) = R_min
                  Region(k)%LMin(3) = Z_min
!
                  Region(k)%LMax(1) = R_max
                  Region(k)%LMax(3) = Z_max
!
                  Region(k)%Equ1CoeffA(0:i1) = BoundSurface1_EquCoeff_A(0:i1)
                  Region(k)%Equ2CoeffA(0:i2) = BoundSurface2_EquCoeff_A(0:i2)
!
                  IF( Region(k)%TypeEqu1(1:4) == 'Powe' ) Region(k)%Equ1CoeffB(0:i1) = BoundSurface1_EquCoeff_B(0:i1)
                  IF( Region(k)%TypeEqu2(1:4) == 'Powe' ) Region(k)%Equ2CoeffB(0:i2) = BoundSurface2_EquCoeff_B(0:i2)
!
                  Region(k)%L1Shift(1) = R1_shift
                  Region(k)%L1Shift(3) = Z1_shift
!
                  Region(k)%L2Shift(1) = R2_shift
                  Region(k)%L2Shift(3) = Z2_shift
!
               END IF IF_Coord4
!
!
!
            END DO DO_Regions
!
! <<<<<<<<<
! <<<<<<<<<
! <<<<<<<<<
!
            IF( k /= TotNum_HetRegions ) THEN
               WRITE(UNIT = *, FMT = 6705) k, TotNum_HetRegions, number_of_regions, number_of_periodic_regions, total_number_periodic_subdomains
               STOP
            END IF
!
!***********************************************************************
!*                                                                     *
!*          Convert the region ranges into SI units (m)                *
!*                                                                     *
!***********************************************************************
!
! ......... Conversion into METERS from INCHES
!
            FORALL (n=2:TotNum_HetRegions, Region(k)%units == 'IN' .OR. Region(k)%units == 'in' .OR. Region(k)%units == 'In')
!
               Region(n)%LMin(1:3) = Region(n)%LMin(1:3) * 2.54d-2
               Region(n)%LMax(1:3) = Region(n)%LMax(1:3) * 2.54d-2
!
               Region(n)%CylBase1Coord(1:3) = Region(n)%CylBase1Coord(1:3) * 2.54d-2
               Region(n)%CylBase2Coord(1:3) = Region(n)%CylBase2Coord(1:3) * 2.54d-2
!
               Region(n)%CylRmin = Region(n)%CylRmin * 2.54d-2
               Region(n)%CylRmax = Region(n)%CylRmax * 2.54d-2
!
               Region(n)%SphereCenterCoord(1:3) = Region(n)%SphereCenterCoord(1:3) * 2.54d-2
!
               Region(n)%SphRmin = Region(n)%SphRmin * 2.54d-2
               Region(n)%SphRmax = Region(n)%SphRmax * 2.54d-2
!
               Region(n)%L1Shift(1:3) = Region(n)%L1Shift(1:3) * 2.54d-2
               Region(n)%L2Shift(1:3) = Region(n)%L2Shift(1:3) * 2.54d-2
!
               Region(n)%EllBase1Coord(1:3) = Region(n)%EllBase1Coord(1:3) * 2.54d-2
               Region(n)%EllBase1Coord(1:3) = Region(n)%EllBase1Coord(1:3) * 2.54d-2
!
               Region(n)%LAxis = Region(n)%LAxis * 2.54d-2
               Region(n)%SAxis = Region(n)%SAxis * 2.54d-2
!
               Region(n)%Base1LAxis = Region(n)%Base1LAxis * 2.54d-2
               Region(n)%Base1SAxis = Region(n)%Base1SAxis * 2.54d-2
               Region(n)%Base2LAxis = Region(n)%Base2LAxis * 2.54d-2
               Region(n)%Base2SAxis = Region(n)%Base2SAxis * 2.54d-2
!
               Region(n)%thick0 = Region(n)%thick0 * 2.54d-2
               Region(n)%thick1 = Region(n)%thick1 * 2.54d-2
               Region(n)%thick2 = Region(n)%thick2 * 2.54d-2
!
            END FORALL
!
! ......... Conversion into METERS from FEET
!
            FORALL (n=2:TotNum_HetRegions, Region(n)%units == 'FT' .OR. Region(n)%units == 'ft' .OR. Region(n)%units == 'Ft')
!
               Region(n)%LMin(1:3) = Region(n)%LMin(1:3) * 3.038d-1
               Region(n)%LMax(1:3) = Region(n)%LMax(1:3) * 3.038d-1
!
               Region(n)%CylBase1Coord(1:3) = Region(n)%CylBase1Coord(1:3) * 3.038d-1
               Region(n)%CylBase2Coord(1:3) = Region(n)%CylBase2Coord(1:3) * 3.038d-1
!
               Region(n)%CylRmin = Region(n)%CylRmin * 3.038d-1
               Region(n)%CylRmax = Region(n)%CylRmax * 3.038d-1
!
               Region(n)%SphereCenterCoord(1:3) = Region(n)%SphereCenterCoord(1:3) * 3.038d-1
!
               Region(n)%SphRmin = Region(n)%SphRmin * 3.038d-1
               Region(n)%SphRmax = Region(n)%SphRmax * 3.038d-1
!
               Region(n)%L1Shift(1:3) = Region(n)%L1Shift(1:3) * 3.038d-1
               Region(n)%L2Shift(1:3) = Region(n)%L2Shift(1:3) * 3.038d-1
!
               Region(n)%EllBase1Coord(1:3) = Region(n)%EllBase1Coord(1:3) * 3.038d-1
               Region(n)%EllBase1Coord(1:3) = Region(n)%EllBase1Coord(1:3) * 3.038d-1
!
               Region(n)%LAxis = Region(n)%LAxis * 3.038d-1
               Region(n)%SAxis = Region(n)%SAxis * 3.038d-1
!
               Region(n)%Base1LAxis = Region(n)%Base1LAxis * 3.038d-1
               Region(n)%Base1SAxis = Region(n)%Base1SAxis * 3.038d-1
               Region(n)%Base2LAxis = Region(n)%Base2LAxis * 3.038d-1
               Region(n)%Base2SAxis = Region(n)%Base2SAxis * 3.038d-1
!
               Region(n)%thick0 = Region(n)%thick0 * 3.038d-1
               Region(n)%thick1 = Region(n)%thick1 * 3.038d-1
               Region(n)%thick2 = Region(n)%thick2 * 3.038d-1
!
            END FORALL
!
! ......... Conversion into METERS from KM
!
            FORALL (n=2:TotNum_HetRegions, Region(n)%units == 'KM' .OR. Region(n)%units == 'km' .OR. Region(n)%units == 'Km')
!
               Region(n)%LMin(1:3) = Region(n)%LMin(1:3) * 1.0d3
               Region(n)%LMax(1:3) = Region(n)%LMax(1:3) * 1.0d3
!
               Region(n)%CylBase1Coord(1:3) = Region(n)%CylBase1Coord(1:3) * 1.0d3
               Region(n)%CylBase2Coord(1:3) = Region(n)%CylBase2Coord(1:3) * 1.0d3
!
               Region(n)%CylRmin = Region(n)%CylRmin * 1.0d3
               Region(n)%CylRmax = Region(n)%CylRmax * 1.0d3
!
               Region(n)%SphereCenterCoord(1:3) = Region(n)%SphereCenterCoord(1:3) * 1.0d3
!
               Region(n)%SphRmin = Region(n)%SphRmin * 1.0d3
               Region(n)%SphRmax = Region(n)%SphRmax * 1.0d3
!
               Region(n)%L1Shift(1:3) = Region(n)%L1Shift(1:3) * 1.0d3
               Region(n)%L2Shift(1:3) = Region(n)%L2Shift(1:3) * 1.0d3
!
               Region(n)%EllBase1Coord(1:3) = Region(n)%EllBase1Coord(1:3) * 1.0d3
               Region(n)%EllBase1Coord(1:3) = Region(n)%EllBase1Coord(1:3) * 1.0d3
!
               Region(n)%LAxis = Region(n)%LAxis * 1.0d3
               Region(n)%SAxis = Region(n)%SAxis * 1.0d3
!
               Region(n)%Base1LAxis = Region(n)%Base1LAxis * 1.0d3
               Region(n)%Base1SAxis = Region(n)%Base1SAxis * 1.0d3
               Region(n)%Base2LAxis = Region(n)%Base2LAxis * 1.0d3
               Region(n)%Base2SAxis = Region(n)%Base2SAxis * 1.0d3
!
               Region(n)%thick0 = Region(n)%thick0 * 1.0d3
               Region(n)%thick1 = Region(n)%thick1 * 1.0d3
               Region(n)%thick2 = Region(n)%thick2 * 1.0d3
!
            END FORALL
!
! ......... Conversion into METERS from MM
!
            FORALL (n=2:TotNum_HetRegions, Region(n)%units == 'MM' .OR. Region(n)%units == 'mm' .OR. Region(n)%units == 'Mm')
!
               Region(n)%LMin(1:3) = Region(n)%LMin(1:3) * 1.0d-3
               Region(n)%LMax(1:3) = Region(n)%LMax(1:3) * 1.0d-3
!
               Region(n)%CylBase1Coord(1:3) = Region(n)%CylBase1Coord(1:3) * 1.0d-3
               Region(n)%CylBase2Coord(1:3) = Region(n)%CylBase2Coord(1:3) * 1.0d-3
!
               Region(n)%CylRmin = Region(n)%CylRmin * 1.0d-3
               Region(n)%CylRmax = Region(n)%CylRmax * 1.0d-3
!
               Region(n)%SphereCenterCoord(1:3) = Region(n)%SphereCenterCoord(1:3) * 1.0d-3
!
               Region(n)%SphRmin = Region(n)%SphRmin * 1.0d-3
               Region(n)%SphRmax = Region(n)%SphRmax * 1.0d-3
!
               Region(n)%L1Shift(1:3) = Region(n)%L1Shift(1:3) * 1.0d-3
               Region(n)%L2Shift(1:3) = Region(n)%L2Shift(1:3) * 1.0d-3
!
               Region(n)%EllBase1Coord(1:3) = Region(n)%EllBase1Coord(1:3) * 1.0d-3
               Region(n)%EllBase1Coord(1:3) = Region(n)%EllBase1Coord(1:3) * 1.0d-3
!
               Region(n)%LAxis = Region(n)%LAxis * 1.0d-3
               Region(n)%SAxis = Region(n)%SAxis * 1.0d-3
!
               Region(n)%Base1LAxis = Region(n)%Base1LAxis * 1.0d-3
               Region(n)%Base1SAxis = Region(n)%Base1SAxis * 1.0d-3
               Region(n)%Base2LAxis = Region(n)%Base2LAxis * 1.0d-3
               Region(n)%Base2SAxis = Region(n)%Base2SAxis * 1.0d-3
!
               Region(n)%thick0 = Region(n)%thick0 * 1.0d-3
               Region(n)%thick1 = Region(n)%thick1 * 1.0d-3
               Region(n)%thick2 = Region(n)%thick2 * 1.0d-3
!
            END FORALL
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Define_Heterogeneous_Regions 1.0 ....... 23 September 2019',6X,'Defining regions of heterogeneous media in the grid')
!
 6010 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The shape of the heterogeneous region # ',I3.3,' is "',A,'": Unknown/Unavailable option'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6050 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The direction <',A,'> of the increasing coordinate in the table of regions in file <',A,'> does not exist',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6052 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input ',A,'_min =',ES12.5,' of the media region # ',I3.3,' is larger than (or equal to) the ',A,'_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6053 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'In the cylindrical media region # ',I3.3,', the centers of the two bases coincide ',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6054 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <cylinder_Rmin> =',ES12.5,' of the cylindrical media region # ',I3.3,  &
     &          ' is larger than (or equal to) the <cylinder_Rmax> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6055 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <sphere_Rmin> =',ES12.5,' of the spherical media region # ',I3.3,  &
     &          ' is larger than (or equal to) the <sphere_Rmax> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6056 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <Base1_Short_Axis> =',ES12.5,' of the elliptical media region # ',I3.3,  &
     &          ' is larger than (or equal to) the <Base1_Long_Axis> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6058 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <Base2_Short_Axis> =',ES12.5,' of the elliptical media region # ',I3.3,  &
     &          ' is larger than (or equal to) the <Base2_Long_Axis> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6059 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <plane_of_ellipse_bases> =',A,' of the elliptical media region # ',I3.3,' is NOT an available option',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6101 FORMAT(T5,'Memory allocation to array <',A,'> in subroutine <Define_Heterogeneous_Regions> was successful')
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to array <',A,'> in subroutine <Define_Heterogeneous_Regions> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6105 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The subdomain shape (read by "region_shape" = ',a11,' in subroutine <Define_Heterogeneous_Regions>) is unavailable',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6401 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: The region #',i3.3,' has an irregular shape that is described by ',/, &
     &       T10,'                                               an unknown/unavailable type of equation (= <',A,'>)',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6402 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions: The region #',i3.3,' is defined by two irregular surfaces of type "FIXED" but the ' ,/,  &
     &       T10,'                                              name of the needed reference interpolation file is not defined',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6515 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: The dataset/namelist <',A,'> must be ended by the "<<<" descriptor'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6520 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: It is not possible to have a rectangular region in a cylindrical coordinate system ',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6600 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: There is a problem reading the namelist <',A,'>',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6601 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: There is a problem reading the namelist <',A,'> in medium #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6602 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: There is a conflict between the various radii options that describe the ', /,&
     &       T10,'                                               ',A,' shape of region #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6604 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: At least one the orders of the equations describing the bounding surfaces 1 and/or 2 ', /,&
     &       T10,'                                               of region #',i3.3,' is < 0',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6605 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: The region #',i3.3,' has an irregular shape but the types of the ',/, &
     &       T10,'     equations describing the bounding surfaces 1 and/or 2 of the region are non-blanks',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6606 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: This ',A,' system has 1 active dimension - it is not possible for the ', /, &
     &       T10,'                                               region #',i3.3,' to be irregularly shaped',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6608 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: In this ',A,' system, one or more of the dependent variables defining the bounding surfaces 1 and 2', /,&
     &       T10,'                                               of region #',i3.3,' is not among the possible options (',A,')',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6610 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: There is a problem reading the 1st data line in the region table <',A,'> in medium #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6612 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: There is a problem reading the coordinate data in the region table <',A,'> in medium #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6701 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: The shape of region #',i3.3,' is "PERIODIC", but the coordinate system = "',A,'"' /,&
     &       T10,'                                               conflicts with the <axis_of_perodicity> = "',A1,'"',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6702 FORMAT(//,20('ERROR-'),//,   &
     &       T10,'>>>  Procedure <Define_Heterogeneous_Regions>: The shape of region #',i3.3,' is "PERIODIC", but the axis of periodicity = "',A1,'"' /,&
     &       T10,'                                               is not among the available options ("X,", "R", "Y", or "Z")',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6705 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'>>>  Procedure <Define_Heterogeneous_Regions>: The total number of computed subdomains = ',i3.3,' does not match the number <TotNum_HetRegions> that is estimated from:',/, &
     &       T5,'                                               TotNum_HetRegions(=',i3.3,') = number_of_regions(=',i3.3,') - number_of_periodic_regions(=',i3.3,') + total_number_periodic_subdomains(=',i3.3,')',/, &
     &       T5,'                                               Check all the related inputs', //, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Define_Heterogeneous_Regions>
!
!
            RETURN
!
         END SUBROUTINE Define_Heterogeneous_Regions
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Define_Boundaries
!
            USE MeshMaker_Data, ONLY: coordinates, NXMax, NYMax, NZMax
            USE Het_Region_Definition
            USE Grid_Generation_Parameters, ONLY: pi
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*            Routine defining the boundaries in the grid              *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: D1, D2, D3
!
            REAL(KIND = 8) :: X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max, emissivity
!
            REAL(KIND = 8) :: sphere_Rmax, sphere_Rmin, cylinder_Rmin, cylinder_Rmax, cylinder_radius1, cylinder_radius2
            REAL(KIND = 8) :: Base1_Long_Axis, Base2_Long_Axis, Base1_Short_Axis, Base2_Short_Axis, Long_Axis, Short_Axis, ellipse_long_axis_angle
!
            REAL(KIND = 8) :: X1_shift, Y1_shift, Z1_shift, R1_shift, X2_shift, Y2_shift, Z2_shift, R2_shift
!
            REAL(KIND = 8) :: exponent1, sign1, exponent2, sign2, thickness0, thickness1, thickness2, interpolation_search_radius
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), DIMENSION(3) :: sphere_center_coordinates, CylBase1_center_coordinates, CylBase2_center_coordinates
            REAL(KIND = 8), DIMENSION(3) :: Base1_center_coordinates,  Base2_center_coordinates
!
            REAL(KIND = 8), DIMENSION(0:10) :: BoundSurface1_EquCoeff_A, BoundSurface1_EquCoeff_B, BoundSurface1_EquCoeff_C, BoundSurface1_EquCoeff_D
            REAL(KIND = 8), DIMENSION(0:10) :: BoundSurface2_EquCoeff_A, BoundSurface2_EquCoeff_B, BoundSurface2_EquCoeff_C, BoundSurface2_EquCoeff_D
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: ier, m, n, i1, i2
!
            INTEGER :: number_of_boundaries
!
            INTEGER :: equation_order_of_bounding_surface1, equation_order_of_bounding_surface2, num_interrogated_points
!
            INTEGER :: number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2
!
! -------------
! ......... CHARACTER variables
! -------------
!
            CHARACTER(LEN =  1) :: vertical_or_stratigraphic1, vertical_or_stratigraphic2, reference_surface1, reference_surface2
            CHARACTER(LEN =  2) :: length_units, plane_of_ellipse_bases
            CHARACTER(LEN =  3) :: boundary_type
            CHARACTER(LEN =  5) :: boundary_name
            CHARACTER(LEN =  6) :: first_part
            CHARACTER(LEN =  8) :: interpolation_data_file_name1, interpolation_data_file_name2, ref_interpolation_data_file
            CHARACTER(LEN = 12) :: boundary_shape, type_of_equation1, type_of_equation2
            CHARACTER(LEN = 50) :: read_data_format
!
            CHARACTER(LEN = 1) :: dependent_variable_of_surfaces
!
! -------------
! ......... LOGICAL variables
! -------------
!
            LOGICAL :: IntDat_match, read_data_by_row
!
! -------
! ... Namelists
! -------
!
      NAMELIST/ Boundary_Regions / number_of_boundaries
!
      NAMELIST/ Boundary_GeneralInfo / boundary_name, boundary_type, boundary_shape, length_units, emissivity
!
      NAMELIST/ Rectangular_Boundary / X_min, Y_min, Z_min, X_max, Y_max, Z_max
!
      NAMELIST/ Cylindrical_Boundary / CylBase1_center_coordinates, CylBase2_center_coordinates, &
     &                                 cylinder_Rmin, cylinder_Rmax, cylinder_radius1, cylinder_radius2
!
      NAMELIST/ Spherical_Boundary   / sphere_center_coordinates, sphere_Rmax, sphere_Rmin
!
      NAMELIST/ Elliptical_Boundary  / Base1_center_coordinates, Base2_center_coordinates, ellipse_long_axis_angle, plane_of_ellipse_bases, &
     &                                 Base1_Long_Axis, Base2_Long_Axis, Base1_Short_Axis, Base2_Short_Axis, Long_Axis, Short_Axis
!
      NAMELIST/ Irregular_Boundary   /  dependent_variable_of_surfaces, type_of_equation1, type_of_equation2, &
     &                                  X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max,               &
     &                                  interpolation_data_file_name1, interpolation_data_file_name2,         &
     &                                  vertical_or_stratigraphic1, vertical_or_stratigraphic2, reference_surface1, reference_surface2, &
     &                                  thickness0, thickness1, thickness2, ref_interpolation_data_file

!
      NAMELIST/ Irregular_Boundary_Surf1 / equation_order_of_bounding_surface1,                      &
     &                                     X1_shift, Y1_shift, Z1_shift, R1_shift, exponent1, sign1, &
     &                                     BoundSurface1_EquCoeff_A, BoundSurface1_EquCoeff_B, BoundSurface1_EquCoeff_C, BoundSurface1_EquCoeff_D
!
      NAMELIST/ Irregular_Boundary_Surf2 / equation_order_of_bounding_surface2,                      &
     &                                     X2_shift, Y2_shift, Z2_shift, R2_shift, exponent2, sign2, &
     &                                     BoundSurface2_EquCoeff_A, BoundSurface2_EquCoeff_B, BoundSurface2_EquCoeff_C, BoundSurface2_EquCoeff_D
!
      NAMELIST/ Irregular_Boundary_IntTable1 / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                         read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
      NAMELIST/ Irregular_Boundary_IntTable2 / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                         read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
      NAMELIST/ Irregular_Boundary_RefTable  / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                         read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Define_Boundaries>
!
!
            WRITE( UNIT = *, FMT = 6000 )
!
! ------------
! ......... Read the number of heterogeneous domains
! ------------
!
            READ (UNIT = *, NML = Boundary_Regions, IOSTAT = ier )
!
! ......... Stop if there is a problem reading the namelist
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6600) '<Boundary_Regions>'
               STOP
            END IF
!
            IF( number_of_boundaries == 0 ) RETURN
!
            Num_Boundaries = number_of_boundaries
!
! ------------
! ......... Allocate memory to the arrays defining the boundaries
! ------------
!
            ALLOCATE(bound(Num_Boundaries), STAT=ier)
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'bound'
            ELSE
               WRITE(*,6102) 'bound'
               STOP
            END IF
!
! ......... Storing temporarily table names, before re-dimensioning the grid/interpolation data sets
!
            i1 = Num_IntTables
            i2 = Num_IntTables + 2 * Num_Boundaries
!
            IF( i1 > 0 ) temp_IntTable(1:i1) = IntTable(1:i1)
!
            DEALLOCATE( IntTable )
            ALLOCATE( IntTable(i2), STAT=ier )
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'IntTable'
            ELSE
               WRITE(*,6102) 'IntTable'
               STOP
            END IF
!
            IF( i1 > 0 ) IntTable(1:i1) = temp_IntTable
!
            DEALLOCATE( temp_IntTable )
            ALLOCATE( temp_IntTable(i2), STAT=ier )
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'temp_IntTable'
            ELSE
               WRITE(*,6102) 'temp_IntTable'
               STOP
            END IF
!
!***********************************************************************
!*                                                                     *
!*        Read the specifics of the remaining media boundaries         *
!*                                                                     *
!***********************************************************************
!
            DO_Boundaries: DO n = 1, Num_Boundaries + 1
!
! -----------------
! ............ Read boundary specifics for Num_Boundaries > 1
! -----------------
!
               IF( n == Num_Boundaries + 1 ) THEN
!
                  READ(UNIT = *, FMT = '(A3)') first_part
!
                  IF( first_part /= '<<<' ) THEN
                     WRITE(UNIT = *, FMT = 6515) first_part
                     STOP
                  ELSE
                     EXIT DO_Boundaries
                  END IF
!
               END IF
!
! ............ Initializations - Namelist components
!
               boundary_shape = '           '
               boundary_name  = '     '
               boundary_type  = '   '
               length_units   = 'm'
               emissivity     = 0.0d0
!
               bound(n)%id    = '   '
               bound(n)%name  = '     '
               bound(n)%shape = '           '
               bound(n)%units = 'm'
!
               bound(n)%emissivity = 0.0d0

!
! -----------------
! ............ Read general info on the boundary
! -----------------
!
               READ (UNIT = *, NML = Boundary_GeneralInfo, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6601) '<Boundary_GeneralInfo>',n
                  STOP
               END IF
!
               SELECT CASE( boundary_shape(1:1) )
               CASE('R', 'r', 'C', 'c', 'S', 's', 'E', 'e', 'I', 'i' )
                  CONTINUE
               CASE DEFAULT
                  WRITE( UNIT = *, FMT = 6010 ) n, boundary_shape
                  STOP
               END SELECT
!
               IF( ( coordinates(1:2) == 'CY' .OR. coordinates(1:4) == 'Cy' .OR. coordinates(1:4) == 'cy' ) .AND. &
     &             ( boundary_shape(1:1) == 'R' .OR. boundary_shape(1:1) == 'r' ) )  &
     &         THEN
                  WRITE(*,6520) boundary_shape
                  STOP
               END IF
!
! ............ Assignment of the namelist values
!
               bound(n)%name  = boundary_name
               bound(n)%shape = boundary_shape
               bound(n)%units = length_units
!
               bound(n)%emissivity = emissivity
!
! ............ Setting the boundary type
!
               bound(n)%id = boundary_type
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    Reading the shape-specific data from the remaining namelists     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
               SELECT CASE(boundary_shape(1:1))
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Rectangular boundary
! >>>>>>>>>>>>
!
               CASE( 'R', 'r' )
!
! ............... Initializations
!
                  X_min = 0.0d0
                  Y_min = 0.0d0
                  Z_min = 0.0d0
                  X_max = 0.0d0
                  Y_max = 0.0d0
                  Z_max = 0.0d0
!
                  bound(n)%LMin = 0.0d0     ! ... Whole array operations
                  bound(n)%LMax = 0.0d0
!
! ............... Reading the rectangle-related data using a namelist
!
                  READ (UNIT = *, NML = Rectangular_Boundary, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Rectangular_Boundary>', n
                     STOP
                  END IF
!
                  bound(n)%LMin(1) = X_min
                  bound(n)%LMin(2) = Y_min
                  bound(n)%LMin(3) = Z_min
!
                  bound(n)%LMax(1) = X_max
                  bound(n)%LMax(2) = Y_max
                  bound(n)%LMax(3) = Z_max
!
! ............... Ensuring that the LMin of the range is not >= of the LMax
!
                  IF( bound(n)%LMin(1) >= bound(n)%LMax(1) ) THEN
                     WRITE(*,6050) bound(n)%LMin(1), n, bound(n)%LMax(1)
                     STOP
                  END IF
!
                  IF( bound(n)%LMin(2) >= bound(n)%LMax(2) ) THEN
                     WRITE(*,6051) bound(n)%LMin(2), n, bound(n)%LMax(2)
                     STOP
                  END IF
!
                  IF( bound(n)%LMin(3) >= bound(n)%LMax(3) ) THEN
                     WRITE(*,6052)  bound(n)%LMin(3), n, bound(n)%LMax(3)
                     STOP
                  END IF
!
                  CYCLE DO_Boundaries
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Cylindrical boundary
! >>>>>>>>>>>>
!
               CASE( 'C', 'c' )
!
! ............... Initializations
!
                  CylBase1_center_coordinates = 0.0d0  ! ... Whole array operation
                  CylBase2_center_coordinates = 0.0d0
                  cylinder_Rmin    = 0.0d0
                  cylinder_Rmax    = 0.0d0
                  cylinder_radius1 = 0.0d0
                  cylinder_radius2 = 0.0d0
!
                  bound(n)%CylBase1Coord = 0.0d0    ! ... Whole array operation
                  bound(n)%CylBase1Coord = 0.0d0
                  bound(n)%CylRmin   = 0.0d0
                  bound(n)%CylRmax   = 0.0d0
                  bound(n)%CylBase1R = 0.0d0
                  bound(n)%CylBase2R = 0.0d0
!
! ............... Reading the cylinder-related data using a namelist
!
                  READ (UNIT = *, NML = Cylindrical_Boundary, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) 'Cylindrical_Boundary', n
                     STOP
                  END IF
!
! ............... For Cartesian coordinates
!
                  IF_Coord: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                          coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     bound(n)%CylBase1Coord = CylBase1_center_coordinates
                     bound(n)%CylBase2Coord = CylBase2_center_coordinates
!
                     IF( ( cylinder_Rmin > 1.0d-7 .OR. cylinder_Rmax > 1.0d-7 ) .AND. ( cylinder_radius1 < 1.0d-7 .AND. cylinder_radius2 < 1.0d-7 ) ) THEN
                        bound(n)%CylRmin = cylinder_Rmin
                        bound(n)%CylRmax = cylinder_Rmax
                     ELSE IF( ( cylinder_Rmin < 1.0d-7 .AND. cylinder_Rmax < 1.0d-7 ) .AND. ( cylinder_radius1 > 1.0d-7 .OR. cylinder_radius2 > 1.0d-7 ) ) THEN
                        bound(n)%CylBase1R = cylinder_radius1
                        bound(n)%CylBase2R = cylinder_radius2
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'cylindrical',n
                        STOP
                     END IF
!
                     IF( cylinder_radius1 /= 0.0d0 ) CYCLE DO_Boundaries
!
! .................. Ensuring that the two bases of the cylindrical boundary do not coincide
!
                     D1 = bound(n)%CylBase1Coord(1) - bound(n)%CylBase2Coord(1)
                     D2 = bound(n)%CylBase1Coord(2) - bound(n)%CylBase2Coord(2)
                     D3 = bound(n)%CylBase1Coord(3) - bound(n)%CylBase2Coord(3)
!
                     IF( ABS(D1) < 1.0d-7 .AND. ABS(D2) < 1.0d-7 .AND. ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! .................. Ensuring that the Rmin of the cylindrical boundary is not >= of the Rmax
!
                     IF( bound(n)%CylRmin >= bound(n)%CylRmax ) THEN
                       WRITE(*,6054) bound(n)%CylRmin, n, bound(n)%CylRmax
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     bound(n)%CylBase1Coord(3) = CylBase1_center_coordinates(3)
                     bound(n)%CylBase2Coord(3) = CylBase2_center_coordinates(3)
!
                     IF( ( cylinder_Rmax > 1.0d-7 ) .AND. ( cylinder_radius1 < 1.0d-7 .AND. cylinder_radius2 < 1.0d-7 ) ) THEN
                        bound(n)%CylRmin = cylinder_Rmin
                        bound(n)%CylRmax = cylinder_Rmax
                     ELSE IF( ( cylinder_Rmax < 1.0d-7 ) .AND. ( cylinder_radius1 > 1.0d-7 .AND. cylinder_radius2 > 1.0d-7 ) ) THEN
                        bound(n)%CylBase1R = cylinder_radius1
                        bound(n)%CylBase2R = cylinder_radius2
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'cylindrical',n
                        STOP
                     END IF
!
                     IF( cylinder_radius1 /= 0.0d0 ) CYCLE DO_Boundaries
!
! .................. Ensuring that the two bases of the cylindrical boundary do not coincide
!
                     D3 = bound(n)%CylBase1Coord(3) - bound(n)%CylBase2Coord(3)
!
                     IF( ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! .................. Ensuring that the Rmin of the cylindrical boundary is not >= of the Rmax
!
                     IF( bound(n)%CylRmin >= bound(n)%CylRmax ) THEN
                       WRITE(*,6054) bound(n)%CylRmin, n, bound(n)%CylRmax
                        STOP
                     END IF
!
                  END IF IF_Coord
!
                  CYCLE DO_Boundaries
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Elliptical boundary
! >>>>>>>>>>>>
!
               CASE( 'E', 'e' )
!
! ............... Initializations
!
                  Base1_center_coordinates = 0.0d0  ! ... Whole array operation
                  Base2_center_coordinates = 0.0d0
                  ellipse_long_axis_angle  = 0.0d0
!
                  plane_of_ellipse_bases = '  '
!
                  Base1_Long_Axis  = 0.0d0
                  Base2_Long_Axis  = 0.0d0
                  Base1_Short_Axis = 0.0d0
                  Base2_Short_Axis = 0.0d0
                  Long_Axis  = 0.0d0
                  Short_Axis = 0.0d0
!
                  bound(n)%EllBase1Coord = 0.0d0    ! ... Whole array operation
                  bound(n)%EllBase1Coord = 0.0d0
                  bound(n)%EllipsePlane  = '  '
!
                  bound(n)%LAxisAngle = 0.0d0
                  bound(n)%Base1LAxis = 0.0d0
                  bound(n)%Base2LAxis = 0.0d0
                  bound(n)%Base1SAxis = 0.0d0
                  bound(n)%Base2SAxis = 0.0d0
                  bound(n)%LAxis = 0.0d0
                  bound(n)%SAxis = 0.0d0
!
! ............... Reading the cylinder-related data using a namelist
!
                  READ (UNIT = *, NML = Elliptical_Boundary, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) 'Elliptical_Boundary', n
                     STOP
                  END IF
!
! ............... Ensuring that there is no conflict between the long and short axes lengths
!
                  IF( Base1_Short_Axis > Base1_Long_Axis ) THEN
                     WRITE(*,6056) Base1_Short_Axis, n, Base1_Long_Axis
                     STOP
                  END IF
!
                  IF( Base2_Short_Axis > Base2_Long_Axis ) THEN
                     WRITE(*,6058) Base2_Short_Axis, n, Base2_Long_Axis
                     STOP
                  END IF
!
                  SELECT CASE(plane_of_ellipse_bases)
                  CASE('XY','Xy','xY','xy','XZ','Xz','xZ','xz','YZ','Yz','yZ','yz')
                     CONTINUE
                  CASE DEFAULT
                     WRITE(UNIT = *, FMT = 6059) plane_of_ellipse_bases, n
                     STOP
                  END SELECT
!
! ............... For Cartesian coordinates
!
                  IF_Coord2: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                           coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     bound(n)%EllBase1Coord = Base1_center_coordinates
                     bound(n)%EllBase2Coord = Base2_center_coordinates
                     bound(n)%LAxisAngle    = pi * ellipse_long_axis_angle / 1.8d2
                     bound(n)%EllipsePlane  = plane_of_ellipse_bases
!
                     IF( ( Long_Axis > 1.0d-7 .AND. Short_Axis > 1.0d-7 ) .AND. &
     &                   ( Base1_Long_Axis < 1.0d-7 .AND. Base1_Long_Axis < 1.0d-7 .AND. Base1_Short_Axis < 1.0d-7 .AND. Base2_Short_Axis < 1.0d-7 ) ) &
     &               THEN
                        bound(n)%LAxis = Long_Axis
                        bound(n)%SAxis = Short_Axis
                     ELSE IF( ( Long_Axis < 1.0d-7 .AND. Short_Axis < 1.0d-7 ) .AND. &
     &                        ( Base1_Long_Axis > 1.0d-7 .AND. Base1_Long_Axis > 1.0d-7 .AND. Base1_Short_Axis > 1.0d-7 .AND. Base2_Short_Axis > 1.0d-7 ) ) &
     &               THEN
                        bound(n)%Base1LAxis = Base1_Long_Axis
                        bound(n)%Base2LAxis = Base2_Long_Axis
                        bound(n)%Base1SAxis = Base1_Short_Axis
                        bound(n)%Base2SAxis = Base2_Short_Axis
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'elliptical',n
                        STOP
                     END IF
!
! .................. Ensuring that the two bases of the elliptical boundary do not coincide
!
                     D1 = bound(n)%EllBase1Coord(1) - bound(n)%EllBase2Coord(1)
                     D2 = bound(n)%EllBase1Coord(2) - bound(n)%EllBase2Coord(2)
                     D3 = bound(n)%EllBase1Coord(3) - bound(n)%EllBase2Coord(3)
!
                     IF( ABS(D1) < 1.0d-7 .AND. ABS(D2) < 1.0d-7 .AND. ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     bound(n)%EllBase1Coord(3) = Base1_center_coordinates(3)
                     bound(n)%EllBase2Coord(3) = Base2_center_coordinates(3)
                     bound(n)%EllipsePlane     = '  '
                     bound(n)%LAxisAngle       = 0.0d0
!
                     IF( ( Long_Axis > 1.0d-7 .AND. Short_Axis > 1.0d-7 ) .AND. &
     &                   ( Base1_Long_Axis < 1.0d-7 .AND. Base1_Long_Axis < 1.0d-7 .AND. Base1_Short_Axis < 1.0d-7 .AND. Base2_Short_Axis < 1.0d-7 ) ) &
     &               THEN
                        bound(n)%LAxis = Long_Axis
                        bound(n)%SAxis = Short_Axis
                     ELSE IF( ( Long_Axis < 1.0d-7 .AND. Short_Axis < 1.0d-7 ) .AND. &
     &                        ( Base1_Long_Axis > 1.0d-7 .AND. Base1_Long_Axis > 1.0d-7 .AND. Base1_Short_Axis > 1.0d-7 .AND. Base2_Short_Axis > 1.0d-7 ) ) &
     &               THEN
                        bound(n)%Base1LAxis = Base1_Long_Axis
                        bound(n)%Base2LAxis = Base2_Long_Axis
                        bound(n)%Base1SAxis = Base1_Short_Axis
                        bound(n)%Base2SAxis = Base2_Short_Axis
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'elliptical',n
                        STOP
                     END IF
!
! .................. Ensuring that the two bases of the cylindrical boundary do not coincide
!
                     D3 = bound(n)%CylBase1Coord(3) - bound(n)%CylBase2Coord(3)
!
                     IF( ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
                  END IF IF_Coord2
!
                  CYCLE DO_Boundaries
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Spherical boundary
! >>>>>>>>>>>>
!
               CASE( 'S', 's' )
!
! ............... Initializations
!
                  sphere_center_coordinates = 0.0d0    ! ... Whole array operation
                  sphere_Rmax = 0.0d0
                  sphere_Rmin = 0.0d0
!
                  bound(n)%SphereCenterCoord = 0.0d0  ! ... Whole array operation
                  bound(n)%SphRmin = 0.0d0
                  bound(n)%SphRmax = 0.0d0
!
! ............... Reading the sphere-related data using a namelist
!
                  READ (UNIT = *, NML = Spherical_Boundary, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) 'Spherical_Boundary', n
                     STOP
                  END IF
!
                  bound(n)%SphRmin = sphere_Rmin
                  bound(n)%SphRmax = sphere_Rmax
!
! ............... For Cartesian coordinates
!
                  IF_Coord3: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                           coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     bound(n)%SphereCenterCoord = sphere_center_coordinates
!
! .................. Ensuring that the Rmin of the spherical boundary is not >= of the Rmax
!
                     IF( bound(n)%SphRmin >= bound(n)%SphRmax ) THEN
                        WRITE(*,6055) bound(n)%SphRmin, n, bound(n)%CylRmax
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     bound(n)%SphereCenterCoord(3) = sphere_center_coordinates(3)
!
! .................. Ensuring that the Rmin of the spherical boundary is not >= of the Rmax
!
                     IF( bound(n)%SphRmin >= bound(n)%SphRmax ) THEN
                        WRITE(*,6055) bound(n)%SphRmin, n, bound(n)%CylRmax
                        STOP
                     END IF
!
                  END IF IF_Coord3
!
                  CYCLE DO_Boundaries
!
               END SELECT
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Irregular boundary
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               type_of_equation1 = '           '  ! ... Initialization
               type_of_equation2 = '           '
!
               X_min = 0.0d0
               Y_min = 0.0d0
               Z_min = 0.0d0
               X_max = 0.0d0
               Y_max = 0.0d0
               Z_max = 0.0d0
!
               X1_shift = 0.0d0
               Y1_shift = 0.0d0
               Z1_shift = 0.0d0
!
               X2_shift = 0.0d0
               Y2_shift = 0.0d0
               Z2_shift = 0.0d0
!
               R_min = 0.0d0
               R_max = 0.0d0
!
               R1_shift = 0.0d0
               R2_shift = 0.0d0
!
               dependent_variable_of_surfaces = ' '
               equation_order_of_bounding_surface1 = 0
               equation_order_of_bounding_surface2 = 0
!
               dependent_variable_of_surfaces = ' '
               equation_order_of_bounding_surface1 = 0
               equation_order_of_bounding_surface2 = 0
!
               interpolation_data_file_name1 = '        '
               interpolation_data_file_name2 = '        '
               ref_interpolation_data_file   = '        '
!
               sign1 = 1.0d0
               sign2 = 1.0d0
               exponent1 = 0.0d0
               exponent2 = 0.0d0
!
               thickness0 = 0.0d0
               thickness1 = 0.0d0
               thickness2 = 0.0d0
!
               vertical_or_stratigraphic1 = ' '
               vertical_or_stratigraphic2 = ' '
!
               reference_surface1  = ' '
               reference_surface2  = ' '
!
! ............ Initializations - Array components
!
               bound(n)%TypeEqu1 = '           '
               bound(n)%TypeEqu2 = '           '
!
               bound(n)%DependentVar = ' '
!
               bound(n)%OrderEquSurf1 = 0
               bound(n)%OrderEquSurf2 = 0
!
               bound(n)%LMin = 0.0d0     ! ... Whole array operations
               bound(n)%LMax = 0.0d0
!
               bound(n)%L1Shift = 0.0d0
               bound(n)%L2Shift = 0.0d0
!
               bound(n)%sign1  = 1.0d0
               bound(n)%sign2  = 1.0d0
               bound(n)%expon1 = 0.0d0
               bound(n)%expon2 = 0.0d0
!
               bound(n)%IntTableNum1 = 0
               bound(n)%IntTableNum2 = 0
!
               bound(n)%thick0 = 0.0d0
               bound(n)%thick1 = 0.0d0
               bound(n)%thick2 = 0.0d0
!
               bound(n)%VertOrStrat1 = ' '
               bound(n)%VertOrStrat2 = ' '
               bound(n)%RefSurf1     = ' '
               bound(n)%RefSurf2     = ' '
!
! ............ Reading the irregular zone data using a namelist
!
               READ (UNIT = *, NML = Irregular_Boundary, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6601) '<Irregular_Boundary>', n
                  STOP
               END IF
!
! ............ Checking the dependent variable
!
               IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &             coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &         THEN

!
                  IF( ( dependent_variable_of_surfaces /= 'X' .AND. dependent_variable_of_surfaces /= 'x'  .AND.  &
     &                  dependent_variable_of_surfaces /= 'Y' .AND. dependent_variable_of_surfaces /= 'y'  .AND.  &
     &                  dependent_variable_of_surfaces /= 'Z' .AND. dependent_variable_of_surfaces /= 'z' ) )     &
     &            THEN
                     WRITE( UNIT = *, FMT = 6608 ) 'cartesian', n, 'X, Y and Z'
                     STOP
                  END IF
!
               ELSE
!
                  IF( ( dependent_variable_of_surfaces /= 'R' .AND. dependent_variable_of_surfaces /= 'r'   .AND. &
     &                  dependent_variable_of_surfaces /= 'Z' .AND. dependent_variable_of_surfaces /= 'z' ) )     &
     &            THEN
                     WRITE( UNIT = *, FMT = 6608 ) 'cylindrical', n, 'R and Z'
                     STOP
                  END IF
!
               END IF
!
! ............ Assignment
!
               bound(n)%LMin(1) = X_min
               bound(n)%LMin(2) = Y_min
               bound(n)%LMin(3) = Z_min
!
               bound(n)%LMax(1) = X_max
               bound(n)%LMax(2) = Y_max
               bound(n)%LMax(3) = Z_max
!
               bound(n)%DependentVar = dependent_variable_of_surfaces
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Setting the type of equation of surface 1
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( type_of_equation1(1:4) == 'POLY' .OR. type_of_equation1(1:4) == 'Poly' .OR. type_of_equation1(1:4) == 'poly' ) THEN
!
! ............ Polynomial equation
!
                  type_of_equation1(1:4) = 'Poly'
                  IF( type_of_equation1(5:6) == '/E' .OR. type_of_equation1(5:6) == '/e' ) THEN
                     type_of_equation1(5:6) = '/E'
                  ELSE IF( type_of_equation1(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                     STOP
                  END IF
!
! ............ Power equation
!
               ELSE IF( type_of_equation1(1:4) == 'POWE' .OR. type_of_equation1(1:4) == 'Powe' .OR. type_of_equation1(1:4) == 'powe' ) THEN
                  type_of_equation1(1:4) = 'Powe'
                  IF( type_of_equation1(5:6) == '/E' .OR. type_of_equation1(5:6) == '/e' ) THEN
                     type_of_equation1(5:6) = '/E'
                  ELSE IF( type_of_equation1(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                     STOP
                  END IF
!
! ............ No equation: interpolation from a grid/tabular data set
!
               ELSE IF( type_of_equation1(1:4) == 'INTE' .OR. type_of_equation1(1:4) == 'Inte' .OR. type_of_equation1(1:4) == 'inte' ) THEN
!
                  CONTINUE
!
! ............ For Fixed Width surfaces
!
               ELSE IF( type_of_equation1(1:4) == 'FIXE' .OR. type_of_equation1(1:4) == 'Fixe' .OR. type_of_equation1(1:4) == 'fixe' ) THEN
!
! ............... Making sure that only one Fixed Width surface is defined
!
                  IF( ( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) .AND. &
     &                ( ref_interpolation_data_file(1:3) == '   ' ) )  &
     &            THEN
                     WRITE( UNIT = *, FMT = 6402 ) n
                     STOP
                  END IF
!
                  CONTINUE
!
! ............ Unknown option: ERROR!
!
               ELSE
                  WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                  STOP
               END IF
!
               bound(n)%TypeEqu1 = type_of_equation1
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Setting the type of equation of surface 2
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( type_of_equation2(1:4) == 'POLY' .OR. type_of_equation2(1:4) == 'Poly' .OR. type_of_equation2(1:4) == 'poly' ) THEN
!
! ............ Polynomial equation
!
                  type_of_equation2(1:4) = 'Poly'
                  IF( type_of_equation2(5:6) == '/E' .OR. type_of_equation2(5:6) == '/e' ) THEN
                     type_of_equation2(5:6) = '/E'
                  ELSE IF( type_of_equation2(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                     STOP
                  END IF
!
! ............ Power equation
!
               ELSE IF( type_of_equation2(1:4) == 'POWE' .OR. type_of_equation2(1:4) == 'Powe' .OR. type_of_equation2(1:4) == 'powe' ) THEN
                  type_of_equation2(1:4) = 'Powe'
                  IF( type_of_equation2(5:6) == '/E' .OR. type_of_equation2(5:6) == '/e' ) THEN
                     type_of_equation2(5:6) = '/E'
                  ELSE IF( type_of_equation2(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                     STOP
                  END IF
!
! ............ No equation: interpolation from a grid/tabular data set
!
               ELSE IF( type_of_equation2(1:4) == 'INTE' .OR. type_of_equation2(1:4) == 'Inte' .OR. type_of_equation2(1:4) == 'inte' ) THEN
!
                  CONTINUE
!
! ............ No equation: Fixed distance from another surface
!
               ELSE IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) THEN
!
                  CONTINUE
!
! ............ Unknown option: ERROR!
!
               ELSE
                  WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                  STOP
               END IF
!
               bound(n)%TypeEqu2 = type_of_equation2
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For Cartesian coodinates
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF_Coord4: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                        coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &         THEN

!
                  IF( ( NXMax >  1 .AND. NYMax == 1 .AND. NZMax == 1 ) .OR. ( NXMax == 1 .AND. NYMax >  1 .AND. NZMax == 1 ) .OR.  &
     &                ( NXMax == 1 .AND. NYMax == 1 .AND. NZMax >  1 ) )    &
     &            THEN
                     WRITE( UNIT = *, FMT = 6606 ) 'cartesian', n
                     STOP
                  END IF
!
                  IF_Inte: IF( type_of_equation1(1:4) == 'INTE' .OR. type_of_equation1(1:4) == 'Inte' .OR. type_of_equation1(1:4) == 'inte' ) THEN
!
                     type_of_equation1(1:4) = 'Inte'
!
                     IntDat_match = .FALSE.
                     IF( Num_IntTables > 0 ) THEN
                        DO_m1: DO m = 1,Num_IntTables
                           IF( interpolation_data_file_name1 == IntTable(m)%FileName ) THEN
                              bound(n)%IntTableNum1 = m
                              IntDat_match          = .TRUE.
                              EXIT
                           END IF
                        END DO DO_m1
                     END IF
!
                     IF( .NOT. IntDat_match ) THEN
!
                        m = Num_IntTables + 1
                        Num_IntTables         = m
                        bound(n)%IntTableNum1 = Num_IntTables
!
! ..................... Reading info describing the structure of interpolation Table 1
!
                        READ (UNIT = *, NML = Irregular_Boundary_IntTable1, IOSTAT = ier )
!
                        IF(ier /= 0) THEN
                          WRITE(UNIT = *, FMT = 6601) '<Irregular_Boundary_IntTable1>', n
                           STOP
                        END IF
!
                        CALL Read_Tabular_Data( table_number = m, file_name = interpolation_data_file_name1, number_of_rows = number_of_rows, &
     &                                          number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                          read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                          RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                          interpolation_search_radius = interpolation_search_radius, &
     &                                          num_interrogated_points     = num_interrogated_points )
!
                     END IF
!
                     GO TO 1500
!
                  ELSE IF( type_of_equation1(1:4) == 'FIXE' .OR. type_of_equation1(1:4) == 'Fixe' .OR. type_of_equation1(1:4) == 'fixe' )  THEN
!
                     bound(n)%thick0 = thickness0
                     bound(n)%thick1 = thickness1
!
                     bound(n)%VertOrStrat1 = vertical_or_stratigraphic1
                     bound(n)%RefSurf1     = reference_surface1
!
                     IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) THEN
!
                        type_of_equation1(1:4) = 'Fixe'
                        type_of_equation2(1:4) = 'Fixe'
!
                        IntDat_match = .FALSE.
!
                        IF( Num_IntTables > 0 ) THEN
                           DO_m2: DO m = 1,Num_IntTables
                              IF( ref_interpolation_data_file == IntTable(m)%FileName ) THEN
                                 bound(n)%IntTableNum1 = m
                                 IntDat_match = .TRUE.
                                 EXIT
                              END IF
                           END DO DO_m2
                        END IF
!
                        IF( .NOT. IntDat_match ) THEN
!
                           m = Num_IntTables + 1
                           Num_IntTables = m
                           bound(n)%IntTableNum1 = Num_IntTables
!
! ........................ Reading info describing the structure of the reference
!
                           READ (UNIT = *, NML = Irregular_Boundary_RefTable, IOSTAT = ier )
!
                           IF(ier /= 0) THEN
                              WRITE(UNIT = *, FMT = 6601) '<Irregular_Boundary_RefTable>', n
                              STOP
                           END IF
!
                           CALL Read_Tabular_Data( table_number = m, file_name = ref_interpolation_data_file, number_of_rows = number_of_rows, &
     &                                             number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                             read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                             RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                             interpolation_search_radius = interpolation_search_radius,  &
     &                                             num_interrogated_points     = num_interrogated_points )
!
                        END IF
!
                     END IF
!
                     GO TO 1500
!
                  END IF IF_Inte

!
! ............... Initializations
!
                  BoundSurface1_EquCoeff_A = 0.0d0
                  BoundSurface1_EquCoeff_B = 0.0d0
                  BoundSurface1_EquCoeff_C = 0.0d0
                  BoundSurface1_EquCoeff_D = 0.0d0
!
! ............... Reading the data for surface 1
!
                  READ (UNIT = *, NML = Irregular_Boundary_Surf1, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_Boundary_Surf1>',n
                     STOP
                  END IF
!
                  bound(n)%OrderEquSurf1 = equation_order_of_bounding_surface1
!
                  i1 = equation_order_of_bounding_surface1
!
                  ALLOCATE ( bound(n)%Equ1CoeffA(0:i1), bound(n)%Equ1CoeffB(0:i1) )
!
                  bound(n)%Equ1CoeffA = 0.0d0
                  bound(n)%Equ1CoeffB = 0.0d0
!
                  IF_Surf1: IF( bound(n)%TypeEqu1(1:4) == 'Powe' ) THEN
                     ALLOCATE ( bound(n)%Equ1CoeffC(0:i1), bound(n)%Equ1CoeffD(0:i1) )
                     bound(n)%Equ1CoeffC = 0.0d0
                     bound(n)%Equ1CoeffD = 0.0d0
                  ELSE
                     IF(i1 > 1 ) THEN
                        ALLOCATE ( bound(n)%Equ1CoeffC(1:Max(1,i1)) )
                        bound(n)%Equ1CoeffC = 0.0d0
                     END IF
!
                  END IF IF_Surf1
!
! ............... Assignments
!
                  bound(n)%sign1  = sign1
                  bound(n)%expon1 = exponent1
!
                  bound(n)%Equ1CoeffA(0:i1) = BoundSurface1_EquCoeff_A(0:i1)
                  bound(n)%Equ1CoeffB(0:i1) = BoundSurface1_EquCoeff_B(0:i1)
!
                  IF( bound(n)%TypeEqu1(1:4) == 'Powe' ) THEN
                     bound(n)%Equ1CoeffC(0:i1) = BoundSurface1_EquCoeff_C(0:i1)
                     bound(n)%Equ1CoeffD(0:i1) = BoundSurface1_EquCoeff_D(0:i1)
                  ELSE
                     IF(i1 > 1) bound(n)%Equ1CoeffC(1:i1-1) = BoundSurface1_EquCoeff_C(1:i1-1)
                  END IF
!
                  bound(n)%L1Shift(1) = X1_shift
                  bound(n)%L1Shift(2) = Y1_shift
                  bound(n)%L1Shift(3) = Z1_shift
!
! >>>>>
!
 1500             IF_Inte2: IF( type_of_equation2(1:4) == 'INTE' .OR. type_of_equation2(1:4) == 'Inte' .OR. type_of_equation2(1:4) == 'inte' ) THEN
!
                     type_of_equation2(1:4) = 'Inte'
!
                     IntDat_match = .FALSE.
                     IF( Num_IntTables > 0 ) THEN
                        DO_m3: DO m = 1,Num_IntTables
                           IF( interpolation_data_file_name2 == IntTable(m)%FileName ) THEN
                              bound(n)%IntTableNum2 = m
                              IntDat_match          = .TRUE.
                              EXIT
                           END IF
                        END DO DO_m3
                     END IF
!
                     IF( .NOT. IntDat_match ) THEN
!
                        m = Num_IntTables + 1
                        Num_IntTables = m
                        bound(n)%IntTableNum2 = Num_IntTables
!
! ..................... Reading info describing the structure of interpolation Table 2
!
                        READ (UNIT = *, NML = Irregular_Boundary_IntTable2, IOSTAT = ier )
!
                        IF(ier /= 0) THEN
                          WRITE(UNIT = *, FMT = 6601) '<Irregular_Boundary_IntTable2>', n
                           STOP
                        END IF
!
                        CALL Read_Tabular_Data( table_number = m, file_name = interpolation_data_file_name2, number_of_rows = number_of_rows, &
     &                                          number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                          read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                          RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                          interpolation_search_radius = interpolation_search_radius,  &
     &                                          num_interrogated_points     = num_interrogated_points )
!
                     END IF
!
                     CYCLE DO_Boundaries
!
                  ELSE IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' )  THEN
!
                     bound(n)%thick2 = thickness2
!
                     bound(n)%VertOrStrat2 = vertical_or_stratigraphic2
                     bound(n)%RefSurf2     = reference_surface2
!
                     CYCLE DO_Boundaries
!
                  END IF IF_Inte2
!
! ............... Initializations
!
                  BoundSurface2_EquCoeff_A = 0.0d0
                  BoundSurface2_EquCoeff_B = 0.0d0
                  BoundSurface2_EquCoeff_C = 0.0d0
                  BoundSurface2_EquCoeff_D = 0.0d0
!
! ............... Reading the data for surface 2
!
                  READ (UNIT = *, NML = Irregular_Boundary_Surf2, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_Boundary_Surf2>',n
                     STOP
                  END IF
!
! ............... Checking consistency of orders of equation
!
                  IF( ( equation_order_of_bounding_surface1 < 0 ) .OR. ( equation_order_of_bounding_surface2 < 0 ) ) THEN
                     WRITE( UNIT = *, FMT = 6604 ) n
                     STOP
                  END IF
!
                  bound(n)%OrderEquSurf2 = equation_order_of_bounding_surface2
!
                  i2 = equation_order_of_bounding_surface2
!
                  ALLOCATE ( bound(n)%Equ2CoeffA(0:i2), bound(n)%Equ2CoeffB(0:i2) )
!
                  bound(n)%Equ2CoeffA = 0.0d0
                  bound(n)%Equ2CoeffB = 0.0d0
!
                  IF_Surf2: IF( bound(n)%TypeEqu2(1:4) == 'Powe' ) THEN
                     ALLOCATE ( bound(n)%Equ2CoeffC(0:i2), bound(n)%Equ2CoeffD(0:i2) )
                     bound(n)%Equ2CoeffC = 0.0d0
                     bound(n)%Equ2CoeffD = 0.0d0
                  ELSE
                     IF(i2 > 1 ) THEN
                        ALLOCATE ( bound(n)%Equ2CoeffC(1:Max(1,i2)) )
                        bound(n)%Equ2CoeffC = 0.0d0
                     END IF
                  END IF IF_Surf2
!
! ............... Assignments
!
                  bound(n)%sign2  = sign2
                  bound(n)%expon2 = exponent2
!
                  bound(n)%Equ2CoeffA(0:i2) = BoundSurface2_EquCoeff_A(0:i2)
                  bound(n)%Equ2CoeffB(0:i2) = BoundSurface2_EquCoeff_B(0:i2)
!
                  IF( bound(n)%TypeEqu2(1:4) == 'Powe' ) THEN
                     bound(n)%Equ2CoeffC(0:i2) = BoundSurface2_EquCoeff_C(0:i2)
                     bound(n)%Equ2CoeffD(0:i2) = BoundSurface2_EquCoeff_D(0:i2)
                  ELSE
                     IF(i2 > 1) bound(n)%Equ2CoeffC(1:i2-1) = BoundSurface2_EquCoeff_C(1:i2-1)
                  END IF
!
                  bound(n)%L2Shift(1) = X2_shift
                  bound(n)%L2Shift(2) = Y2_shift
                  bound(n)%L2Shift(3) = Z2_shift
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For Cylindrical coodinates
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               ELSE
!
! ............... Irregular shapes
!
                  IF( NXMax >  1 .AND. NZMax == 1 ) THEN
                     WRITE( UNIT = *, FMT = 6606 ) 'radial', n
                     STOP
                  END IF
!
! ............... Initializations
!
                  BoundSurface1_EquCoeff_A = 0.0d0
                  BoundSurface1_EquCoeff_B = 0.0d0
                  BoundSurface2_EquCoeff_A = 0.0d0
                  BoundSurface2_EquCoeff_B = 0.0d0
!
! ............... Reading the data for surface 1
!
                  READ (UNIT = *, NML = Irregular_Boundary_Surf1, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_Boundary_Surf1>',n
                     STOP
                  END IF
!
! ............... Reading the data for surface 2
!
                  READ (UNIT = *, NML = Irregular_Boundary_Surf2, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_Boundary_Surf2>',n
                     STOP
                   END IF
!
! ............... Checking consistency of orders of equation
!
                  IF( ( equation_order_of_bounding_surface1 < 0 ) .OR. ( equation_order_of_bounding_surface2 < 0 ) ) THEN
                     WRITE( UNIT = *, FMT = 6604 ) n
                     STOP
                  END IF
!
                  bound(n)%DependentVar = dependent_variable_of_surfaces
!
                  bound(n)%OrderEquSurf1 = equation_order_of_bounding_surface1
                  bound(n)%OrderEquSurf2 = equation_order_of_bounding_surface2
!
                  i1 = equation_order_of_bounding_surface1
                  ALLOCATE ( bound(n)%Equ1CoeffA(0:i1) )
!
                  bound(n)%Equ1CoeffA = 0.0d0
!
                  IF( bound(n)%TypeEqu1(1:4) == 'Powe' ) THEN
                     ALLOCATE ( bound(n)%Equ1CoeffB(0:i1) )
                     bound(n)%Equ1CoeffB = 0.0d0
                  END IF
!
                  i2 = equation_order_of_bounding_surface2
                  ALLOCATE ( bound(n)%Equ2CoeffA(0:i2) )
!
                  bound(n)%Equ2CoeffA = 0.0d0
!
                  IF( bound(n)%TypeEqu2(1:4) == 'Powe' ) THEN
                     ALLOCATE ( bound(n)%Equ2CoeffB(0:i2) )
                     bound(n)%Equ2CoeffB = 0.0d0
                  END IF
!
! ............... Assignments
!
                  bound(n)%LMin(1) = R_min
                  bound(n)%LMin(3) = Z_min
!
                  bound(n)%LMax(1) = R_max
                  bound(n)%LMax(3) = Z_max
!
                  bound(n)%Equ1CoeffA(0:i1) = BoundSurface1_EquCoeff_A(0:i1)
                  bound(n)%Equ2CoeffA(0:i2) = BoundSurface2_EquCoeff_A(0:i2)
!
                  IF( bound(n)%TypeEqu1(1:4) == 'Powe' ) bound(n)%Equ1CoeffB(0:i1) = BoundSurface1_EquCoeff_B(0:i1)
                  IF( bound(n)%TypeEqu2(1:4) == 'Powe' ) bound(n)%Equ2CoeffB(0:i2) = BoundSurface2_EquCoeff_B(0:i2)
!
                  bound(n)%L1Shift(1) = R1_shift
                  bound(n)%L1Shift(3) = Z1_shift
!
                  bound(n)%L2Shift(1) = R2_shift
                  bound(n)%L2Shift(3) = Z2_shift
!
               END IF IF_Coord4
!
! <<<<<<<<<
! <<<<<<<<<
! <<<<<<<<<
!
            END DO DO_Boundaries
!
!***********************************************************************
!*                                                                     *
!*           Convert the boundary ranges into SI units (m)             *
!*                                                                     *
!***********************************************************************
!
! ......... Conversion into METERS from INCHES
!
            FORALL (n=2:Num_Boundaries, bound(n)%units == 'IN' .OR. bound(n)%units == 'in' .OR. bound(n)%units == 'In')
!
               bound(n)%LMin(1:3) = bound(n)%LMin(1:3) * 2.54d-2
               bound(n)%LMax(1:3) = bound(n)%LMax(1:3) * 2.54d-2
!
               bound(n)%CylBase1Coord(1:3) = bound(n)%CylBase1Coord(1:3) * 2.54d-2
               bound(n)%CylBase2Coord(1:3) = bound(n)%CylBase2Coord(1:3) * 2.54d-2
!
               bound(n)%CylRmin = bound(n)%CylRmin * 2.54d-2
               bound(n)%CylRmax = bound(n)%CylRmax * 2.54d-2
!
               bound(n)%SphereCenterCoord(1:3) = bound(n)%SphereCenterCoord(1:3) * 2.54d-2
!
               bound(n)%SphRmin = bound(n)%SphRmin * 2.54d-2
               bound(n)%SphRmax = bound(n)%SphRmax * 2.54d-2
!
               bound(n)%L1Shift(1:3) = bound(n)%L1Shift(1:3) * 2.54d-2
               bound(n)%L2Shift(1:3) = bound(n)%L2Shift(1:3) * 2.54d-2
!
               bound(n)%EllBase1Coord(1:3) = bound(n)%EllBase1Coord(1:3) * 2.54d-2
               bound(n)%EllBase1Coord(1:3) = bound(n)%EllBase1Coord(1:3) * 2.54d-2
!
               bound(n)%LAxis = bound(n)%LAxis * 2.54d-2
               bound(n)%SAxis = bound(n)%SAxis * 2.54d-2
!
               bound(n)%Base1LAxis = bound(n)%Base1LAxis * 2.54d-2
               bound(n)%Base1SAxis = bound(n)%Base1SAxis * 2.54d-2
               bound(n)%Base2LAxis = bound(n)%Base2LAxis * 2.54d-2
               bound(n)%Base2SAxis = bound(n)%Base2SAxis * 2.54d-2
!
               bound(n)%thick0 = bound(n)%thick0 * 2.54d-2
               bound(n)%thick1 = bound(n)%thick1 * 2.54d-2
               bound(n)%thick2 = bound(n)%thick2 * 2.54d-2
!
            END FORALL
!
! ......... Conversion into METERS from FEET
!
            FORALL (n=2:Num_Boundaries, bound(n)%units == 'FT' .OR. bound(n)%units == 'ft' .OR. bound(n)%units == 'Ft')
!
               bound(n)%LMin(1:3) = bound(n)%LMin(1:3) * 3.038d-1
               bound(n)%LMax(1:3) = bound(n)%LMax(1:3) * 3.038d-1
!
               bound(n)%CylBase1Coord(1:3) = bound(n)%CylBase1Coord(1:3) * 3.038d-1
               bound(n)%CylBase2Coord(1:3) = bound(n)%CylBase2Coord(1:3) * 3.038d-1
!
               bound(n)%CylRmin = bound(n)%CylRmin * 3.038d-1
               bound(n)%CylRmax = bound(n)%CylRmax * 3.038d-1
!
               bound(n)%SphereCenterCoord(1:3) = bound(n)%SphereCenterCoord(1:3) * 3.038d-1
!
               bound(n)%SphRmin = bound(n)%SphRmin * 3.038d-1
               bound(n)%SphRmax = bound(n)%SphRmax * 3.038d-1
!
               bound(n)%L1Shift(1:3) = bound(n)%L1Shift(1:3) * 3.038d-1
               bound(n)%L2Shift(1:3) = bound(n)%L2Shift(1:3) * 3.038d-1
!
               bound(n)%EllBase1Coord(1:3) = bound(n)%EllBase1Coord(1:3) * 3.038d-1
               bound(n)%EllBase1Coord(1:3) = bound(n)%EllBase1Coord(1:3) * 3.038d-1
!
               bound(n)%LAxis = bound(n)%LAxis * 3.038d-1
               bound(n)%SAxis = bound(n)%SAxis * 3.038d-1
!
               bound(n)%Base1LAxis = bound(n)%Base1LAxis * 3.038d-1
               bound(n)%Base1SAxis = bound(n)%Base1SAxis * 3.038d-1
               bound(n)%Base2LAxis = bound(n)%Base2LAxis * 3.038d-1
               bound(n)%Base2SAxis = bound(n)%Base2SAxis * 3.038d-1
!
               bound(n)%thick0 = bound(n)%thick0 * 3.038d-1
               bound(n)%thick1 = bound(n)%thick1 * 3.038d-1
               bound(n)%thick2 = bound(n)%thick2 * 3.038d-1
!
            END FORALL
!
! ......... Conversion into METERS from KM
!
            FORALL (n=2:Num_Boundaries, bound(n)%units == 'KM' .OR. bound(n)%units == 'km' .OR. bound(n)%units == 'Km')
!
               bound(n)%LMin(1:3) = bound(n)%LMin(1:3) * 1.0d3
               bound(n)%LMax(1:3) = bound(n)%LMax(1:3) * 1.0d3
!
               bound(n)%CylBase1Coord(1:3) = bound(n)%CylBase1Coord(1:3) * 1.0d3
               bound(n)%CylBase2Coord(1:3) = bound(n)%CylBase2Coord(1:3) * 1.0d3
!
               bound(n)%CylRmin = bound(n)%CylRmin * 1.0d3
               bound(n)%CylRmax = bound(n)%CylRmax * 1.0d3
!
               bound(n)%SphereCenterCoord(1:3) = bound(n)%SphereCenterCoord(1:3) * 1.0d3
!
               bound(n)%SphRmin = bound(n)%SphRmin * 1.0d3
               bound(n)%SphRmax = bound(n)%SphRmax * 1.0d3
!
               bound(n)%L1Shift(1:3) = bound(n)%L1Shift(1:3) * 1.0d3
               bound(n)%L2Shift(1:3) = bound(n)%L2Shift(1:3) * 1.0d3
!
               bound(n)%EllBase1Coord(1:3) = bound(n)%EllBase1Coord(1:3) * 1.0d3
               bound(n)%EllBase1Coord(1:3) = bound(n)%EllBase1Coord(1:3) * 1.0d3
!
               bound(n)%LAxis = bound(n)%LAxis * 1.0d3
               bound(n)%SAxis = bound(n)%SAxis * 1.0d3
!
               bound(n)%Base1LAxis = bound(n)%Base1LAxis * 1.0d3
               bound(n)%Base1SAxis = bound(n)%Base1SAxis * 1.0d3
               bound(n)%Base2LAxis = bound(n)%Base2LAxis * 1.0d3
               bound(n)%Base2SAxis = bound(n)%Base2SAxis * 1.0d3
!
               bound(n)%thick0 = bound(n)%thick0 * 1.0d3
               bound(n)%thick1 = bound(n)%thick1 * 1.0d3
               bound(n)%thick2 = bound(n)%thick2 * 1.0d3
!
            END FORALL
!
! ......... Conversion into METERS from MM
!
            FORALL (n=2:Num_Boundaries, bound(n)%units == 'MM' .OR. bound(n)%units == 'mm' .OR. bound(n)%units == 'Mm')
!
               bound(n)%LMin(1:3) = bound(n)%LMin(1:3) * 1.0d-3
               bound(n)%LMax(1:3) = bound(n)%LMax(1:3) * 1.0d-3
!
               bound(n)%CylBase1Coord(1:3) = bound(n)%CylBase1Coord(1:3) * 1.0d-3
               bound(n)%CylBase2Coord(1:3) = bound(n)%CylBase2Coord(1:3) * 1.0d-3
!
               bound(n)%CylRmin = bound(n)%CylRmin * 1.0d-3
               bound(n)%CylRmax = bound(n)%CylRmax * 1.0d-3
!
               bound(n)%SphereCenterCoord(1:3) = bound(n)%SphereCenterCoord(1:3) * 1.0d-3
!
               bound(n)%SphRmin = bound(n)%SphRmin * 1.0d-3
               bound(n)%SphRmax = bound(n)%SphRmax * 1.0d-3
!
               bound(n)%L1Shift(1:3) = bound(n)%L1Shift(1:3) * 1.0d-3
               bound(n)%L2Shift(1:3) = bound(n)%L2Shift(1:3) * 1.0d-3
!
               bound(n)%EllBase1Coord(1:3) = bound(n)%EllBase1Coord(1:3) * 1.0d-3
               bound(n)%EllBase1Coord(1:3) = bound(n)%EllBase1Coord(1:3) * 1.0d-3
!
               bound(n)%LAxis = bound(n)%LAxis * 1.0d-3
               bound(n)%SAxis = bound(n)%SAxis * 1.0d-3
!
               bound(n)%Base1LAxis = bound(n)%Base1LAxis * 1.0d-3
               bound(n)%Base1SAxis = bound(n)%Base1SAxis * 1.0d-3
               bound(n)%Base2LAxis = bound(n)%Base2LAxis * 1.0d-3
               bound(n)%Base2SAxis = bound(n)%Base2SAxis * 1.0d-3
!
               bound(n)%thick0 = bound(n)%thick0 * 1.0d-3
               bound(n)%thick1 = bound(n)%thick1 * 1.0d-3
               bound(n)%thick2 = bound(n)%thick2 * 1.0d-3
!
            END FORALL
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Define_Boundaries 1.0 .................. 18 January   2015',6X,'Defining the boundaries of the grid')
!
 6010 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The shape of the boundary # ',I3.3,' is "',A,'": Unknown/Unavailable option'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6050 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input X_min =',ES12.5,' of the boundary # ',I3.3,' is larger than (or equal to) the X_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6051 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input Y_min =',ES12.5,' of the boundary # ',I3.3,' is larger than (or equal to) the Y_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6052 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input Z_min =',ES12.5,' of the boundary # ',I3.3,' is larger than (or equal to) the Z_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6053 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'In the cylindrical boundary # ',I3.3,', the centers of the two bases coincide ',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6054 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <cylinder_Rmin> =',ES12.5,' of the cylindrical boundary # ',I3.3,  &
     &          ' is larger than (or equal to) the <cylinder_Rmax> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6055 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <sphere_Rmin> =',ES12.5,' of the spherical boundary # ',I3.3,  &
     &          ' is larger than (or equal to) the <sphere_Rmax> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6056 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <Base1_Short_Axis> =',ES12.5,' of the elliptical boundary # ',I3.3,  &
     &          ' is larger than (or equal to) the <Base1_Long_Axis> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6058 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <Base2_Short_Axis> =',ES12.5,' of the elliptical boundary # ',I3.3,  &
     &          ' is larger than (or equal to) the <Base2_Long_Axis> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6059 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <plane_of_ellipse_bases> =',A,' of the elliptical boundary # ',I3.3,' is NOT an available option',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6101 FORMAT(T5,'Memory allocation to arrays <',A,'> in subroutine <Define_Boundaries> was successful')
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <Define_Boundaries> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6105 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The subdomain shape (read by "boundary_shape" = ',a11,' in subroutine <Define_Boundaries>) is unavailable',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6401 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: The boundary #',i3.3,' has an irregular shape that is described by ',/, &
     &       T10,'                                    an unknown/unavailable type of equation (= "',A,'")',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6402 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: The boundary #',i3.3,' is defined by two irregular surfaces of type "FIXED" but the ' ,/,  &
     &       T10,'                                    name of the needed reference interpolation file is not defined',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6515 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: The dataset/namelist <',A,'> must be ended by the "<<<" descriptor'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6520 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: It is not possible to have a rectangular boundary in a cylindrical coordinate system ',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6600 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: There is a problem reading the namelist <',A,'>',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6601 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: There is a problem reading the namelist <',A,'> in boundary #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6602 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: There is a conflict between the various radii options that describe the ', /,&
     &       T10,'                                    ',A,' shape of boundary #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6604 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: At least one the orders of the equations describing the bounding surfaces 1 and/or 2 ', /,&
     &       T10,'                                    of boundary #',i3.3,' is < 0',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6605 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: The boundary #',i3.3,' has an irregular shape but the types of the ',/, &
     &       T10,'     equations describing the bounding surfaces 1 and/or 2 of the boundary are non-blanks',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6606 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: This ',A,' system has 1 active dimension - it is not possible for the ', /, &
     &       T10,'                                    boundary #',i3.3,' to be irregularly shaped',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6608 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Boundaries>: In this ',A,' system, one or more of the dependent variables defining the bounding surfaces 1 and 2', /,&
     &       T10,'                                    of boundary #',i3.3,' is not among the possible options (',A,')',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Define_Boundaries>
!
!
            RETURN
!
         END SUBROUTINE Define_Boundaries
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Define_Exclusion_Zones
!
            USE MeshMaker_Data, ONLY: coordinates, NXMax, NYMax, NZMax
            USE Het_Region_Definition
            USE Grid_Generation_Parameters, ONLY: pi
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          Routine defining the excluded zones in the grid            *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: D1, D2, D3
!
            REAL(KIND = 8) :: X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max
!
            REAL(KIND = 8) :: sphere_Rmax, sphere_Rmin, cylinder_Rmin, cylinder_Rmax, cylinder_radius1, cylinder_radius2
            REAL(KIND = 8) :: Base1_Long_Axis, Base2_Long_Axis, Base1_Short_Axis, Base2_Short_Axis, Long_Axis, Short_Axis, ellipse_long_axis_angle
!
            REAL(KIND = 8) :: X1_shift, Y1_shift, Z1_shift, R1_shift, X2_shift, Y2_shift, Z2_shift, R2_shift
!
            REAL(KIND = 8) :: exponent1, sign1, exponent2, sign2
!
            REAL(KIND = 8) :: location_of_1st_periodic_occurrence, width_of_periodic_ExclZone, period_of_occurrence
            REAL(KIND = 8) :: period, L_first, L_width, thickness0, thickness1, thickness2, interpolation_search_radius
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), DIMENSION(3) :: sphere_center_coordinates, CylBase1_center_coordinates, CylBase2_center_coordinates
            REAL(KIND = 8), DIMENSION(3) :: Base1_center_coordinates,  Base2_center_coordinates
!
            REAL(KIND = 8), DIMENSION(0:10) :: BoundSurface1_EquCoeff_A, BoundSurface1_EquCoeff_B, BoundSurface1_EquCoeff_C, BoundSurface1_EquCoeff_D
            REAL(KIND = 8), DIMENSION(0:10) :: BoundSurface2_EquCoeff_A, BoundSurface2_EquCoeff_B, BoundSurface2_EquCoeff_C, BoundSurface2_EquCoeff_D
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: ier, n, m, i1, i2, j, k, AxisNum
!
            INTEGER :: number_of_exclusion_zones, number_of_periodic_ExclZones, total_number_periodic_ExclZones, number_of_periodic_occurrences
!
            INTEGER :: equation_order_of_bounding_surface1, equation_order_of_bounding_surface2, num_interrogated_points
!
            INTEGER :: number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2
!
            INTEGER, ALLOCATABLE, DIMENSION(:) :: frequency_of_occurrences
!
! -------------
! ......... CHARACTER variables
! -------------
!
            CHARACTER(LEN =  1) :: axis_of_periodicity, vertical_or_stratigraphic1, vertical_or_stratigraphic2, reference_surface1, reference_surface2
            CHARACTER(LEN =  2) :: length_units, plane_of_ellipse_bases
            CHARACTER(LEN =  3) :: exclusion_zone_type
            CHARACTER(LEN =  6) :: first_part
            CHARACTER(LEN =  8) :: interpolation_data_file_name1, interpolation_data_file_name2, ref_interpolation_data_file
            CHARACTER(LEN = 12) :: exclusion_zone_shape, type_of_equation1, type_of_equation2
            CHARACTER(LEN = 50) :: read_data_format
!
            CHARACTER(LEN = 1) :: dependent_variable_of_surfaces
!
! -------------
! ......... LOGICAL variables
! -------------
!
            LOGICAL :: IntDat_match, read_data_by_row
!
! -------
! ... Namelists
! -------
!
      NAMELIST/ Exclusion_Zones / number_of_exclusion_zones, number_of_periodic_ExclZones, total_number_periodic_ExclZones
!
      NAMELIST/ ExclZone_GeneralInfo / exclusion_zone_shape, exclusion_zone_type, length_units
!
      NAMELIST/ Rectangular_ExclZone / X_min, Y_min, Z_min, X_max, Y_max, Z_max
!
      NAMELIST/ Periodic_ExclZone    / X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max,          &
     &                                 number_of_periodic_occurrences, axis_of_periodicity,             &
     &                                 location_of_1st_periodic_occurrence, width_of_periodic_ExclZone, &
     &                                 period_of_occurrence
!
      NAMELIST/ Cylindrical_ExclZone / CylBase1_center_coordinates, CylBase2_center_coordinates, &
     &                                 cylinder_Rmin, cylinder_Rmax, cylinder_radius1, cylinder_radius2
!
      NAMELIST/ Spherical_ExclZone   / sphere_center_coordinates, sphere_Rmax, sphere_Rmin
!
      NAMELIST/ Elliptical_ExclZone  / Base1_center_coordinates, Base2_center_coordinates, ellipse_long_axis_angle, plane_of_ellipse_bases, &
     &                                 Base1_Long_Axis, Base2_Long_Axis, Base1_Short_Axis, Base2_Short_Axis, Long_Axis, Short_Axis
!
      NAMELIST/ Irregular_ExclZone   /  dependent_variable_of_surfaces, type_of_equation1, type_of_equation2, &
     &                                  X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max,               &
     &                                  interpolation_data_file_name1, interpolation_data_file_name2,         &
     &                                  vertical_or_stratigraphic1, vertical_or_stratigraphic2, reference_surface1, reference_surface2, &
     &                                  thickness0, thickness1, thickness2, ref_interpolation_data_file
!
      NAMELIST/ Irregular_ExclZone_Surf1 / equation_order_of_bounding_surface1,                      &
     &                                     X1_shift, Y1_shift, Z1_shift, R1_shift, exponent1, sign1, &
     &                                     BoundSurface1_EquCoeff_A, BoundSurface1_EquCoeff_B, BoundSurface1_EquCoeff_C, BoundSurface1_EquCoeff_D
!
      NAMELIST/ Irregular_ExclZone_Surf2 / equation_order_of_bounding_surface2,                      &
     &                                     X2_shift, Y2_shift, Z2_shift, R2_shift, exponent2, sign2, &
     &                                     BoundSurface2_EquCoeff_A, BoundSurface2_EquCoeff_B, BoundSurface2_EquCoeff_C, BoundSurface2_EquCoeff_D
!
      NAMELIST/ Irregular_ExclZone_IntTable1 / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                         read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
      NAMELIST/ Irregular_ExclZone_IntTable2 / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                         read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
      NAMELIST/ Irregular_ExclZone_RefTable  / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                         read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Define_Exclusion_Zones>
!
!
            WRITE(*,6000)
!
! ......... Initialization
!
            number_of_periodic_ExclZones    = 0
            total_number_periodic_ExclZones = 0
!
! ------------
! ......... Read the number of heterogeneous domains
! ------------
!
            READ (UNIT = *, NML = Exclusion_Zones, IOSTAT = ier )
!
! ......... Stop if there is a problem reading the namelist
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6600) '<Exclusion_Zones>'
               STOP
            END IF
!
            ALLOCATE( frequency_of_occurrences(number_of_periodic_ExclZones) )
!
            frequency_of_occurrences = 0    ! ... Initialization
!
            Num_ExclZones    = number_of_exclusion_zones
            TotNum_ExclZones = number_of_exclusion_zones - number_of_periodic_ExclZones + total_number_periodic_ExclZones
!
            IF(Num_ExclZones <= 0) RETURN
!
! ------------
! ......... Allocate memory to the arrays defining the exclusion zones
! ------------
!
            ALLOCATE(ExclZone(1:TotNum_ExclZones), STAT=ier)
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'ExclZone'
            ELSE
               WRITE(*,6102) 'ExclZone'
               STOP
            END IF
!
! ......... Storing temporarily table names, before re-dimensioning the grid/interpolation data sets
!
            i1 = Num_IntTables
            i2 = Num_IntTables + 2 * Num_ExclZones
!
            IF( i1 > 0 ) temp_IntTable(1:i1) = IntTable(1:i1)
!
            DEALLOCATE( IntTable )
            ALLOCATE( IntTable(i2), STAT=ier )
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'IntTable'
            ELSE
               WRITE(*,6102) 'IntTable'
               STOP
            END IF
!
            IF( i1 > 0 ) IntTable(1:i1) = temp_IntTable
!
            DEALLOCATE( temp_IntTable )
            ALLOCATE( temp_IntTable(i2), STAT=ier )
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'temp_IntTable'
            ELSE
               WRITE(*,6102) 'temp_IntTable'
               STOP
            END IF
!
!***********************************************************************
!*                                                                     *
!*      Read the specifics of the remaining media exclusion zones      *
!*                                                                     *
!***********************************************************************
!
            k = 0
            m = 0
!
            DO_Exclusions: DO n = 1, Num_ExclZones + 1
!
! -----------------
! ............ Read exclusion zone specifics for Num_ExclZones > 1
! -----------------
!
               IF( n == Num_ExclZones + 1 ) THEN
!
                  READ(UNIT = *, FMT = '(A3)') first_part
!
                  IF( first_part /= '<<<' ) THEN
                     WRITE(UNIT = *, FMT = 6515) first_part
                     STOP
                  ELSE
                     EXIT DO_Exclusions
                  END IF
!
               END IF
!
! ............ Advance the counter
!
               k = k + 1
!
! ............ Initializations - Namelist components
!
               exclusion_zone_shape = '           '
               exclusion_zone_type  = '   '
               length_units = 'm'
!
               ExclZone(k)%shape = '           '
               ExclZone(k)%units = 'm'
               ExclZone(k)%id    = '    '
!
! -----------------
! ............ Read general info on the exclusion zone
! -----------------
!
               READ (UNIT = *, NML = ExclZone_GeneralInfo, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6601) '<ExclZone_GeneralInfo>',n
                  STOP
               END IF
!
               SELECT CASE( exclusion_zone_shape(1:1) )
               CASE('R', 'r', 'C', 'c', 'S', 's', 'E', 'e', 'P', 'p', 'I', 'i' )
                  CONTINUE
               CASE DEFAULT
                  WRITE( UNIT = *, FMT = 6010 ) n, exclusion_zone_shape
                  STOP
               END SELECT
!
               IF( ( coordinates(1:2) == 'CY' .OR. coordinates(1:4) == 'Cy' .OR. coordinates(1:4) == 'cy' ) .AND. &
     &             ( exclusion_zone_shape(1:1) == 'R' .OR. exclusion_zone_shape(1:1) == 'r' ) )  &
     &         THEN
                  WRITE(*,6520) exclusion_zone_shape
                  STOP
               END IF
!
! ............ Assignment of the namelist values
!
               ExclZone(k)%shape = exclusion_zone_shape
               ExclZone(k)%id    = exclusion_zone_type
               ExclZone(k)%units = length_units
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    Reading the shape-specific data from the remaining namelists     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
               SELECT CASE(exclusion_zone_shape(1:1))
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Rectangular exclusion zones
! >>>>>>>>>>>>
!
               CASE( 'R', 'r' )
!
! ............... Initializations
!
                  X_min = 0.0d0
                  Y_min = 0.0d0
                  Z_min = 0.0d0
                  X_max = 0.0d0
                  Y_max = 0.0d0
                  Z_max = 0.0d0
!
                  ExclZone(k)%LMin = 0.0d0     ! ... Whole array operations
                  ExclZone(k)%LMax = 0.0d0
!
! ............... Reading the rectangle-related data using a namelist
!
                  READ (UNIT = *, NML = Rectangular_ExclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Rectangular_ExclZone>', n
                     STOP
                  END IF
!
                  ExclZone(k)%LMin(1) = X_min
                  ExclZone(k)%LMin(2) = Y_min
                  ExclZone(k)%LMin(3) = Z_min
!
                  ExclZone(k)%LMax(1) = X_max
                  ExclZone(k)%LMax(2) = Y_max
                  ExclZone(k)%LMax(3) = Z_max
!
! ............... Ensuring that the LMin of the range is not >= of the LMax
!
                  IF( ExclZone(k)%LMin(1) >= ExclZone(k)%LMax(1) ) THEN
                     WRITE(*,6050) ExclZone(k)%LMin(1), n, ExclZone(k)%LMax(1)
                     STOP
                  END IF
!
                  IF( ExclZone(k)%LMin(2) >= ExclZone(k)%LMax(2) ) THEN
                     WRITE(*,6051) ExclZone(k)%LMin(2), n, ExclZone(k)%LMax(2)
                     STOP
                  END IF
!
                  IF( ExclZone(k)%LMin(3) >= ExclZone(k)%LMax(3) ) THEN
                     WRITE(*,6052)  ExclZone(k)%LMin(3), n, ExclZone(k)%LMax(3)
                     STOP
                  END IF
!
                  CYCLE DO_Exclusions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Cylindrical exclusion zones
! >>>>>>>>>>>>
!
               CASE( 'C', 'c' )
!
! ............... Initializations
!
                  CylBase1_center_coordinates = 0.0d0  ! ... Whole array operation
                  CylBase2_center_coordinates = 0.0d0
                  cylinder_Rmin    = 0.0d0
                  cylinder_Rmax    = 0.0d0
                  cylinder_radius1 = 0.0d0
                  cylinder_radius2 = 0.0d0
!
                  ExclZone(k)%CylBase1Coord = 0.0d0    ! ... Whole array operation
                  ExclZone(k)%CylBase1Coord = 0.0d0
                  ExclZone(k)%CylRmin   = 0.0d0
                  ExclZone(k)%CylRmax   = 0.0d0
                  ExclZone(k)%CylBase1R = 0.0d0
                  ExclZone(k)%CylBase2R = 0.0d0
!
! ............... Reading the cylinder-related data using a namelist
!
                  READ (UNIT = *, NML = Cylindrical_ExclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) 'Cylindrical_ExclZone', n
                     STOP
                  END IF
!
! ............... For Cartesian coordinates
!
                  IF_Coord: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                          coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     ExclZone(k)%CylBase1Coord = CylBase1_center_coordinates
                     ExclZone(k)%CylBase2Coord = CylBase2_center_coordinates
!
                     IF( ( cylinder_Rmin > 1.0d-7 .OR. cylinder_Rmax > 1.0d-7 ) .AND. ( cylinder_radius1 < 1.0d-7 .AND. cylinder_radius2 < 1.0d-7 ) ) THEN
                        ExclZone(k)%CylRmin = cylinder_Rmin
                        ExclZone(k)%CylRmax = cylinder_Rmax
                     ELSE IF( ( cylinder_Rmin < 1.0d-7 .AND. cylinder_Rmax < 1.0d-7 ) .AND. ( cylinder_radius1 > 1.0d-7 .OR. cylinder_radius2 > 1.0d-7 ) ) THEN
                        ExclZone(k)%CylBase1R = cylinder_radius1
                        ExclZone(k)%CylBase2R = cylinder_radius2
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'cylindrical',n
                        STOP
                     END IF
!
                     IF( cylinder_radius1 /= 0.0d0 ) CYCLE DO_Exclusions
!
! .................. Ensuring that the two bases of the cylindrical exclusion zone do not coincide
!
                     D1 = ExclZone(k)%CylBase1Coord(1) - ExclZone(k)%CylBase2Coord(1)
                     D2 = ExclZone(k)%CylBase1Coord(2) - ExclZone(k)%CylBase2Coord(2)
                     D3 = ExclZone(k)%CylBase1Coord(3) - ExclZone(k)%CylBase2Coord(3)
!
                     IF( ABS(D1) < 1.0d-7 .AND. ABS(D2) < 1.0d-7 .AND. ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! .................. Ensuring that the Rmin of the cylindrical exclusion zone is not >= of the Rmax
!
                     IF( ExclZone(k)%CylRmin >= ExclZone(k)%CylRmax ) THEN
                       WRITE(*,6054) ExclZone(k)%CylRmin, n, ExclZone(k)%CylRmax
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     ExclZone(k)%CylBase1Coord(3) = CylBase1_center_coordinates(3)
                     ExclZone(k)%CylBase2Coord(3) = CylBase2_center_coordinates(3)
!
                     IF( ( cylinder_Rmax > 1.0d-7 ) .AND. ( cylinder_radius1 < 1.0d-7 .AND. cylinder_radius2 < 1.0d-7 ) ) THEN
                        ExclZone(k)%CylRmin = cylinder_Rmin
                        ExclZone(k)%CylRmax = cylinder_Rmax
                     ELSE IF( ( cylinder_Rmax < 1.0d-7 ) .AND. ( cylinder_radius1 > 1.0d-7 .AND. cylinder_radius2 > 1.0d-7 ) ) THEN
                        ExclZone(k)%CylBase1R = cylinder_radius1
                        ExclZone(k)%CylBase2R = cylinder_radius2
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'cylindrical',n
                        STOP
                     END IF
!
                     IF( cylinder_radius1 /= 0.0d0 ) CYCLE DO_Exclusions
!
! .................. Ensuring that the two bases of the cylindrical exclusion zone do not coincide
!
                     D3 = ExclZone(k)%CylBase1Coord(3) - ExclZone(k)%CylBase2Coord(3)
!
                     IF( ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! .................. Ensuring that the Rmin of the cylindrical exclusion zone is not >= of the Rmax
!
                     IF( ExclZone(k)%CylRmin >= ExclZone(k)%CylRmax ) THEN
                       WRITE(*,6054) ExclZone(k)%CylRmin, n, ExclZone(k)%CylRmax
                        STOP
                     END IF
!
                  END IF IF_Coord
!
                  CYCLE DO_Exclusions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Periodic exclusion zones
! >>>>>>>>>>>>
!
               CASE( 'P', 'p' )
!
! ............... Initializations
!
                  X_min = 0.0d0
                  R_min = 0.0d0
                  Y_min = 0.0d0
                  Z_min = 0.0d0
                  X_max = 0.0d0
                  R_max = 0.0d0
                  Y_max = 0.0d0
                  Z_max = 0.0d0
!
                  exclusion_zone_shape = 'periodic'
!
                  axis_of_periodicity                 = ' '
                  number_of_periodic_occurrences      = 0
                  location_of_1st_periodic_occurrence = 0.0d0
                  width_of_periodic_ExclZone          = 0.0d0
                  period_of_occurrence                = 0.0d0
!
! ............... Reading the rectangle-related data using a namelist
!
                  READ (UNIT = *, NML = Periodic_ExclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Periodic_ExclZone>', n
                     STOP
                  END IF
!
                  m = m + 1
!
                  IF ( axis_of_periodicity == 'X' .AND. axis_of_periodicity == 'x' .AND.  &
     &                 axis_of_periodicity == 'R' .AND. axis_of_periodicity == 'r' .AND.  &
     &                 axis_of_periodicity == 'Y' .AND. axis_of_periodicity == 'y' .AND.  &
     &                 axis_of_periodicity == 'Z' .AND. axis_of_periodicity == 'z' )   &
     &            THEN
                     WRITE( UNIT = *, FMT = 6702 ) n, axis_of_periodicity
                     STOP
                  END IF
!
                  L_first = location_of_1st_periodic_occurrence
                  L_width = width_of_periodic_ExclZone
                  period  = period_of_occurrence
!
                  frequency_of_occurrences(m) = number_of_periodic_occurrences
!
                  IF_Period: IF (axis_of_periodicity == 'X' .OR. axis_of_periodicity == 'x' ) THEN
!
                     AxisNum = 1
                     IF( coordinates(1:4) == 'CYLI' .OR. coordinates(1:4) == 'Cyli' .OR. coordinates(1:4) == 'cyli' ) THEN
                        WRITE( UNIT = *, FMT = 6701 ) n, coordinates, axis_of_periodicity
                        STOP
                     END IF
!
                  ELSE IF (axis_of_periodicity == 'R' .OR. axis_of_periodicity == 'r' ) THEN
!
                     AxisNum = 1
                     IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                   coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &               THEN

                        WRITE( UNIT = *, FMT = 6701 ) n, coordinates, axis_of_periodicity
                        STOP
                     END IF
!
                  ELSE IF (axis_of_periodicity == 'Y' .OR. axis_of_periodicity == 'y' ) THEN
!
                     AxisNum = 2
!
                  ELSE IF (axis_of_periodicity == 'Z' .OR. axis_of_periodicity == 'z' ) THEN
!
                     AxisNum = 3
!
                  END IF IF_Period
!
                  j = 0
!
 1000             IF( AxisNum == 1 ) THEN
                     X_min = L_first + j * period
                     X_max = X_min + L_width
                     R_min = X_min
                     R_max = X_max
                  ELSE IF( AxisNum == 2 ) THEN
                     Y_min = L_first + j * period
                     Y_max = Y_min + L_width
                  ELSE
                     Z_min = L_first + j * period
                     Z_max = Z_min + L_width
                  END IF
!
                  ExclZone(k)%shape = exclusion_zone_shape
                  ExclZone(k)%units = length_units
!
                  ExclZone(k)%LMin = 0.0d0     ! ... Whole array operations
                  ExclZone(k)%LMax = 0.0d0
!
                  ExclZone(k)%LMin(1) = X_min
                  ExclZone(k)%LMin(2) = Y_min
                  ExclZone(k)%LMin(3) = Z_min
!
                  ExclZone(k)%LMax(1) = X_max
                  ExclZone(k)%LMax(2) = Y_max
                  ExclZone(k)%LMax(3) = Z_max
!
! ............... Ensuring that the LMin of the range is not >= of the LMax
!
                  IF( ExclZone(k)%LMin(1) >= ExclZone(k)%LMax(1) ) THEN
                     WRITE(*,6050) ExclZone(k)%LMin(1), n, ExclZone(k)%LMax(1)
                     STOP
                  END IF
!
                  IF( ExclZone(k)%LMin(2) >= ExclZone(k)%LMax(2) ) THEN
                     WRITE(*,6051) ExclZone(k)%LMin(2), n, ExclZone(k)%LMax(2)
                     STOP
                  END IF
!
                  IF( ExclZone(k)%LMin(3) >= ExclZone(k)%LMax(3) ) THEN
                     WRITE(*,6052)  ExclZone(k)%LMin(3), n, ExclZone(k)%LMax(3)
                     STOP
                  END IF
!
                  IF ( j < frequency_of_occurrences(m) - 1 ) THEN
                     j = j + 1
                     k = k + 1
                     GO TO 1000
                  END IF
!
                  CYCLE DO_Exclusions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Elliptical exclusion zones
! >>>>>>>>>>>>
!
               CASE( 'E', 'e' )
!
! ............... Initializations
!
                  Base1_center_coordinates = 0.0d0  ! ... Whole array operation
                  Base2_center_coordinates = 0.0d0
                  ellipse_long_axis_angle  = 0.0d0
!
                  plane_of_ellipse_bases = '  '
!
                  Base1_Long_Axis  = 0.0d0
                  Base2_Long_Axis  = 0.0d0
                  Base1_Short_Axis = 0.0d0
                  Base2_Short_Axis = 0.0d0
                  Long_Axis  = 0.0d0
                  Short_Axis = 0.0d0
!
                  ExclZone(k)%EllBase1Coord = 0.0d0    ! ... Whole array operation
                  ExclZone(k)%EllBase1Coord = 0.0d0
                  ExclZone(k)%EllipsePlane  = '  '
!
                  ExclZone(k)%LAxisAngle = 0.0d0
                  ExclZone(k)%Base1LAxis = 0.0d0
                  ExclZone(k)%Base2LAxis = 0.0d0
                  ExclZone(k)%Base1SAxis = 0.0d0
                  ExclZone(k)%Base2SAxis = 0.0d0
                  ExclZone(k)%LAxis = 0.0d0
                  ExclZone(k)%SAxis = 0.0d0
!
! ............... Reading the cylinder-related data using a namelist
!
                  READ (UNIT = *, NML = Elliptical_ExclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) 'Elliptical_ExclZone', n
                     STOP
                  END IF
!
! ............... Ensuring that there is no conflict between the long and short axes lengths
!
                  IF( Base1_Short_Axis > Base1_Long_Axis ) THEN
                     WRITE(*,6056) Base1_Short_Axis, n, Base1_Long_Axis
                     STOP
                  END IF
!
                  IF( Base2_Short_Axis > Base2_Long_Axis ) THEN
                     WRITE(*,6058) Base2_Short_Axis, n, Base2_Long_Axis
                     STOP
                  END IF
!
                  SELECT CASE(plane_of_ellipse_bases)
                  CASE('XY','Xy','xY','xy','XZ','Xz','xZ','xz','YZ','Yz','yZ','yz')
                     CONTINUE
                  CASE DEFAULT
                     WRITE(UNIT = *, FMT = 6059) plane_of_ellipse_bases, n
                     STOP
                  END SELECT
!
! ............... For Cartesian coordinates
!
                  IF_Coord2: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                           coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     ExclZone(k)%EllBase1Coord = Base1_center_coordinates
                     ExclZone(k)%EllBase2Coord = Base2_center_coordinates
                     ExclZone(k)%LAxisAngle    = pi * ellipse_long_axis_angle / 1.8d2
                     ExclZone(k)%EllipsePlane  = plane_of_ellipse_bases
!
                     IF( ( Long_Axis > 1.0d-7 .AND. Short_Axis > 1.0d-7 ) .AND. &
     &                   ( Base1_Long_Axis < 1.0d-7 .AND. Base1_Long_Axis < 1.0d-7 .AND. Base1_Short_Axis < 1.0d-7 .AND. Base2_Short_Axis < 1.0d-7 ) ) &
     &               THEN
                        ExclZone(k)%LAxis = Long_Axis
                        ExclZone(k)%SAxis = Short_Axis
                     ELSE IF( ( Long_Axis < 1.0d-7 .AND. Short_Axis < 1.0d-7 ) .AND. &
     &                        ( Base1_Long_Axis > 1.0d-7 .AND. Base1_Long_Axis > 1.0d-7 .AND. Base1_Short_Axis > 1.0d-7 .AND. Base2_Short_Axis > 1.0d-7 ) ) &
     &               THEN
                        ExclZone(k)%Base1LAxis = Base1_Long_Axis
                        ExclZone(k)%Base2LAxis = Base2_Long_Axis
                        ExclZone(k)%Base1SAxis = Base1_Short_Axis
                        ExclZone(k)%Base2SAxis = Base2_Short_Axis
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'elliptical',n
                        STOP
                     END IF
!
! .................. Ensuring that the two bases of the elliptical exclusion zone do not coincide
!
                     D1 = ExclZone(k)%EllBase1Coord(1) - ExclZone(k)%EllBase2Coord(1)
                     D2 = ExclZone(k)%EllBase1Coord(2) - ExclZone(k)%EllBase2Coord(2)
                     D3 = ExclZone(k)%EllBase1Coord(3) - ExclZone(k)%EllBase2Coord(3)
!
                     IF( ABS(D1) < 1.0d-7 .AND. ABS(D2) < 1.0d-7 .AND. ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     ExclZone(k)%EllBase1Coord(3) = Base1_center_coordinates(3)
                     ExclZone(k)%EllBase2Coord(3) = Base2_center_coordinates(3)
                     ExclZone(k)%EllipsePlane     = '  '
                     ExclZone(k)%LAxisAngle       = 0.0d0
!
                     IF( ( Long_Axis > 1.0d-7 .AND. Short_Axis > 1.0d-7 ) .AND. &
     &                   ( Base1_Long_Axis < 1.0d-7 .AND. Base1_Long_Axis < 1.0d-7 .AND. Base1_Short_Axis < 1.0d-7 .AND. Base2_Short_Axis < 1.0d-7 ) ) &
     &               THEN
                        ExclZone(k)%LAxis = Long_Axis
                        ExclZone(k)%SAxis = Short_Axis
                     ELSE IF( ( Long_Axis < 1.0d-7 .AND. Short_Axis < 1.0d-7 ) .AND. &
     &                        ( Base1_Long_Axis > 1.0d-7 .AND. Base1_Long_Axis > 1.0d-7 .AND. Base1_Short_Axis > 1.0d-7 .AND. Base2_Short_Axis > 1.0d-7 ) ) &
     &               THEN
                        ExclZone(k)%Base1LAxis = Base1_Long_Axis
                        ExclZone(k)%Base2LAxis = Base2_Long_Axis
                        ExclZone(k)%Base1SAxis = Base1_Short_Axis
                        ExclZone(k)%Base2SAxis = Base2_Short_Axis
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'elliptical',n
                        STOP
                     END IF
!
! .................. Ensuring that the two bases of the cylindrical exclusion zone do not coincide
!
                     D3 = ExclZone(k)%CylBase1Coord(3) - ExclZone(k)%CylBase2Coord(3)
!
                     IF( ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
                  END IF IF_Coord2
!
                  CYCLE DO_Exclusions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Spherical exclusion zones
! >>>>>>>>>>>>
!
               CASE( 'S', 's' )
!
! ............... Initializations
!
                  sphere_center_coordinates = 0.0d0    ! ... Whole array operation
                  sphere_Rmax = 0.0d0
                  sphere_Rmin = 0.0d0
!
                  ExclZone(k)%SphereCenterCoord = 0.0d0  ! ... Whole array operation
                  ExclZone(k)%SphRmin = 0.0d0
                  ExclZone(k)%SphRmax = 0.0d0
!
! ............... Reading the sphere-related data using a namelist
!
                  READ (UNIT = *, NML = Spherical_ExclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) 'Spherical_ExclZone', n
                     STOP
                  END IF
!
                  ExclZone(k)%SphRmin = sphere_Rmin
                  ExclZone(k)%SphRmax = sphere_Rmax
!
! ............... For Cartesian coordinates
!
                  IF_Coord3: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                           coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     ExclZone(k)%SphereCenterCoord = sphere_center_coordinates
!
! .................. Ensuring that the Rmin of the spherical exclusion zone is not >= of the Rmax
!
                     IF( ExclZone(k)%SphRmin >= ExclZone(k)%SphRmax ) THEN
                        WRITE(*,6055) ExclZone(k)%SphRmin, n, ExclZone(k)%CylRmax
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     ExclZone(k)%SphereCenterCoord(3) = sphere_center_coordinates(3)
!
! .................. Ensuring that the Rmin of the spherical exclusion zone is not >= of the Rmax
!
                     IF( ExclZone(k)%SphRmin >= ExclZone(k)%SphRmax ) THEN
                        WRITE(*,6055) ExclZone(k)%SphRmin, n, ExclZone(k)%CylRmax
                        STOP
                     END IF
!
                  END IF IF_Coord3
!
                  CYCLE DO_Exclusions
!
               END SELECT
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Irregular exclusion zones
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               type_of_equation1 = '           '  ! ... Initialization
               type_of_equation2 = '           '
!
               X_min = 0.0d0
               Y_min = 0.0d0
               Z_min = 0.0d0
               X_max = 0.0d0
               Y_max = 0.0d0
               Z_max = 0.0d0
!
               X1_shift = 0.0d0
               Y1_shift = 0.0d0
               Z1_shift = 0.0d0
!
               X2_shift = 0.0d0
               Y2_shift = 0.0d0
               Z2_shift = 0.0d0
!
               R_min = 0.0d0
               R_max = 0.0d0
!
               R1_shift = 0.0d0
               R2_shift = 0.0d0
!
               dependent_variable_of_surfaces = ' '
               equation_order_of_bounding_surface1 = 0
               equation_order_of_bounding_surface2 = 0
!
               interpolation_data_file_name1 = '        '
               interpolation_data_file_name2 = '        '
               ref_interpolation_data_file   = '        '
!
               sign1 = 1.0d0
               sign2 = 1.0d0
               exponent1 = 0.0d0
               exponent2 = 0.0d0
!
               thickness0 = 0.0d0
               thickness1 = 0.0d0
               thickness2 = 0.0d0
!
               vertical_or_stratigraphic1 = ' '
               vertical_or_stratigraphic2 = ' '
!
               reference_surface1  = ' '
               reference_surface2  = ' '
!
! ............ Initializations - Array components
!
               ExclZone(k)%TypeEqu1 = '           '
               ExclZone(k)%TypeEqu2 = '           '
!
               ExclZone(k)%DependentVar = ' '
!
               ExclZone(k)%OrderEquSurf1 = 0
               ExclZone(k)%OrderEquSurf2 = 0
!
               ExclZone(k)%LMin = 0.0d0     ! ... Whole array operations
               ExclZone(k)%LMax = 0.0d0
!
               ExclZone(k)%LMin = 0.0d0
               ExclZone(k)%LMax = 0.0d0
!
               ExclZone(k)%L1Shift = 0.0d0
               ExclZone(k)%L2Shift = 0.0d0
!
               ExclZone(k)%sign1  = 1.0d0
               ExclZone(k)%sign2  = 1.0d0
               ExclZone(k)%expon1 = 0.0d0
               ExclZone(k)%expon2 = 0.0d0
!
               ExclZone(k)%IntTableNum1 = 0
               ExclZone(k)%IntTableNum2 = 0
!
               ExclZone(k)%thick0 = 0.0d0
               ExclZone(k)%thick1 = 0.0d0
               ExclZone(k)%thick2 = 0.0d0
!
               ExclZone(k)%VertOrStrat1 = ' '
               ExclZone(k)%VertOrStrat2 = ' '
               ExclZone(k)%RefSurf1     = ' '
               ExclZone(k)%RefSurf2     = ' '
!
! ............ Reading the irregular zone data using a namelist
!
               READ (UNIT = *, NML = Irregular_ExclZone, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6601) '<Irregular_ExclZone>', n
                  STOP
               END IF
!
! ............ Checking the dependent variable
!
               IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &             coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &         THEN

!
                  IF( ( dependent_variable_of_surfaces /= 'X' .AND. dependent_variable_of_surfaces /= 'x'  .AND.  &
     &                  dependent_variable_of_surfaces /= 'Y' .AND. dependent_variable_of_surfaces /= 'y'  .AND.  &
     &                  dependent_variable_of_surfaces /= 'Z' .AND. dependent_variable_of_surfaces /= 'z' ) )     &
     &            THEN
                     WRITE( UNIT = *, FMT = 6608 ) 'cartesian', n, 'X, Y and Z'
                     STOP
                  END IF
!
               ELSE
!
                  IF( ( dependent_variable_of_surfaces /= 'R' .AND. dependent_variable_of_surfaces /= 'r'   .AND. &
     &                  dependent_variable_of_surfaces /= 'Z' .AND. dependent_variable_of_surfaces /= 'z' ) )     &
     &            THEN
                     WRITE( UNIT = *, FMT = 6608 ) 'cylindrical', n, 'R and Z'
                     STOP
                  END IF
!
               END IF
!
! ............ Assignment
!
               ExclZone(k)%LMin(1) = X_min
               ExclZone(k)%LMin(2) = Y_min
               ExclZone(k)%LMin(3) = Z_min
!
               ExclZone(k)%LMax(1) = X_max
               ExclZone(k)%LMax(2) = Y_max
               ExclZone(k)%LMax(3) = Z_max
!
               ExclZone(k)%DependentVar = dependent_variable_of_surfaces
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Setting the type of equation of surface 1
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( type_of_equation1(1:4) == 'POLY' .OR. type_of_equation1(1:4) == 'Poly' .OR. type_of_equation1(1:4) == 'poly' ) THEN
!
! ............ Polynomial equation
!
                  type_of_equation1(1:4) = 'Poly'
                  IF( type_of_equation1(5:6) == '/E' .OR. type_of_equation1(5:6) == '/e' ) THEN
                     type_of_equation1(5:6) = '/E'
                  ELSE IF( type_of_equation1(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                     STOP
                  END IF
!
! ............ Power equation
!
               ELSE IF( type_of_equation1(1:4) == 'POWE' .OR. type_of_equation1(1:4) == 'Powe' .OR. type_of_equation1(1:4) == 'powe' ) THEN
                  type_of_equation1(1:4) = 'Powe'
                  IF( type_of_equation1(5:6) == '/E' .OR. type_of_equation1(5:6) == '/e' ) THEN
                     type_of_equation1(5:6) = '/E'
                  ELSE IF( type_of_equation1(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                     STOP
                  END IF
!
! ............ No equation: interpolation from a grid/tabular data set
!
               ELSE IF( type_of_equation1(1:4) == 'INTE' .OR. type_of_equation1(1:4) == 'Inte' .OR. type_of_equation1(1:4) == 'inte' ) THEN
!
                  CONTINUE
!
! ............ For Fixed Width surfaces
!
               ELSE IF( type_of_equation1(1:4) == 'FIXE' .OR. type_of_equation1(1:4) == 'Fixe' .OR. type_of_equation1(1:4) == 'fixe' ) THEN
!
! ............... Making sure that only one Fixed Width surface is defined
!
                  IF( ( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) .AND. &
     &                ( ref_interpolation_data_file(1:3) == '   ' ) )  &
     &            THEN
                     WRITE( UNIT = *, FMT = 6402 ) n
                     STOP
                  END IF
!
                  CONTINUE
!
! ............ Unknown option: ERROR!
!
               ELSE
                  WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                  STOP
               END IF
!
               ExclZone(k)%TypeEqu1 = type_of_equation1
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Setting the type of equation of surface 2
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( type_of_equation2(1:4) == 'POLY' .OR. type_of_equation2(1:4) == 'Poly' .OR. type_of_equation2(1:4) == 'poly' ) THEN
!
! ............ Polynomial equation
!
                  type_of_equation2(1:4) = 'Poly'
                  IF( type_of_equation2(5:6) == '/E' .OR. type_of_equation2(5:6) == '/e' ) THEN
                     type_of_equation2(5:6) = '/E'
                  ELSE IF( type_of_equation2(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                     STOP
                  END IF
!
! ............ Power equation
!
               ELSE IF( type_of_equation2(1:4) == 'POWE' .OR. type_of_equation2(1:4) == 'Powe' .OR. type_of_equation2(1:4) == 'powe' ) THEN
                  type_of_equation2(1:4) = 'Powe'
                  IF( type_of_equation2(5:6) == '/E' .OR. type_of_equation2(5:6) == '/e' ) THEN
                     type_of_equation2(5:6) = '/E'
                  ELSE IF( type_of_equation2(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                     STOP
                  END IF
!
! ............ No equation: interpolation from a grid/tabular data set
!
               ELSE IF( type_of_equation2(1:4) == 'INTE' .OR. type_of_equation2(1:4) == 'Inte' .OR. type_of_equation2(1:4) == 'inte' ) THEN
!
                  CONTINUE
!
! ............ No equation: Fixed distance from another surface
!
               ELSE IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) THEN
!
                  CONTINUE
!
! ............ Unknown option: ERROR!
!
               ELSE
                  WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                  STOP
               END IF
!
               ExclZone(k)%TypeEqu2 = type_of_equation2
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For Cartesian coodinates
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF_Coord4: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                        coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &         THEN

!
                  IF( ( NXMax >  1 .AND. NYMax == 1 .AND. NZMax == 1 ) .OR. ( NXMax == 1 .AND. NYMax >  1 .AND. NZMax == 1 ) .OR.  &
     &                ( NXMax == 1 .AND. NYMax == 1 .AND. NZMax >  1 ) )    &
     &            THEN
                     WRITE( UNIT = *, FMT = 6606 ) 'cartesian', n
                     STOP
                  END IF
!
                  IF_Inte: IF( type_of_equation1(1:4) == 'INTE' .OR. type_of_equation1(1:4) == 'Inte' .OR. type_of_equation1(1:4) == 'inte' ) THEN
!
                     type_of_equation1(1:4) = 'Inte'
!
                     IntDat_match = .FALSE.
                     IF( Num_IntTables > 0 ) THEN
                         DO_m1: DO m = 1,Num_IntTables
                            IF( interpolation_data_file_name1 == IntTable(m)%FileName ) THEN
                               ExclZone(k)%IntTableNum1 = m
                               IntDat_match             = .TRUE.
                               EXIT
                            END IF
                         END DO DO_m1
                     END IF
!
                     IF( .NOT. IntDat_match ) THEN
!
                        m = Num_IntTables + 1
                        Num_IntTables            = m
                        ExclZone(k)%IntTableNum1 = Num_IntTables
!
! ..................... Reading info describing the structure of interpolation Table 1
!
                        READ (UNIT = *, NML = Irregular_ExclZone_IntTable1, IOSTAT = ier )
!
                        IF(ier /= 0) THEN
                          WRITE(UNIT = *, FMT = 6601) '<Irregular_ExclZone_IntTable1>', n
                           STOP
                        END IF
!
                        CALL Read_Tabular_Data( table_number = m, file_name = interpolation_data_file_name1, number_of_rows = number_of_rows, &
     &                                          number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                          read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                          RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                          interpolation_search_radius = interpolation_search_radius,  &
     &                                          num_interrogated_points     = num_interrogated_points )

!
                     END IF
!
                     GO TO 1500
!
                  ELSE IF( type_of_equation1(1:4) == 'FIXE' .OR. type_of_equation1(1:4) == 'Fixe' .OR. type_of_equation1(1:4) == 'fixe' )  THEN
!
                     ExclZone(k)%thick0 = thickness0
                     ExclZone(k)%thick1 = thickness1
!
                     ExclZone(k)%VertOrStrat1 = vertical_or_stratigraphic1
                     ExclZone(k)%RefSurf1     = reference_surface1
!
                     IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) THEN
!
                        type_of_equation1(1:4) = 'Fixe'
                        type_of_equation2(1:4) = 'Fixe'
!
                        IntDat_match = .FALSE.
!
                        IF( Num_IntTables > 0 ) THEN
                           DO_m2: DO m = 1,Num_IntTables
                              IF( ref_interpolation_data_file == IntTable(m)%FileName ) THEN
                                 ExclZone(k)%IntTableNum1 = m
                                 IntDat_match           = .TRUE.
                                 EXIT
                              END IF
                           END DO DO_m2
                        END IF
!
                        IF( .NOT. IntDat_match ) THEN
!
                           m = Num_IntTables + 1
                           Num_IntTables          = m
                           ExclZone(k)%IntTableNum1 = Num_IntTables
!
! ........................ Reading info describing the structure of the reference
!
                           READ (UNIT = *, NML = Irregular_ExclZone_RefTable, IOSTAT = ier )
!
                           IF(ier /= 0) THEN
                              WRITE(UNIT = *, FMT = 6601) '<Irregular_ExclZone_RefTable>', n
                              STOP
                           END IF
!
                           CALL Read_Tabular_Data( table_number = m, file_name = ref_interpolation_data_file, number_of_rows = number_of_rows, &
     &                                             number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                             read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                             RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                             interpolation_search_radius = interpolation_search_radius,  &
     &                                             num_interrogated_points     = num_interrogated_points )

!
                        END IF
!
                     END IF
!
                     GO TO 1500
!
                  END IF IF_Inte
!
! ............... Initializations
!
                  BoundSurface1_EquCoeff_A = 0.0d0
                  BoundSurface1_EquCoeff_B = 0.0d0
                  BoundSurface1_EquCoeff_C = 0.0d0
                  BoundSurface1_EquCoeff_D = 0.0d0
!
! ............... Reading the data for surface 1
!
                  READ (UNIT = *, NML = Irregular_ExclZone_Surf1, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_ExclZone_Surf1>',n
                     STOP
                  END IF
!
                  ExclZone(k)%OrderEquSurf1 = equation_order_of_bounding_surface1
!
                  i1 = equation_order_of_bounding_surface1
!
                  ALLOCATE ( ExclZone(k)%Equ1CoeffA(0:i1), ExclZone(k)%Equ1CoeffB(0:i1) )
!
                  ExclZone(k)%Equ1CoeffA = 0.0d0
                  ExclZone(k)%Equ1CoeffB = 0.0d0
!
                  IF_Surf1: IF( ExclZone(k)%TypeEqu1(1:4) == 'Powe' ) THEN
                     ALLOCATE ( ExclZone(k)%Equ1CoeffC(0:i1), ExclZone(k)%Equ1CoeffD(0:i1) )
                     ExclZone(k)%Equ1CoeffC = 0.0d0
                     ExclZone(k)%Equ1CoeffD = 0.0d0
                  ELSE
                     IF(i1 > 1 ) THEN
                        ALLOCATE ( ExclZone(k)%Equ1CoeffC(1:Max(1,i1)) )
                        ExclZone(k)%Equ1CoeffC = 0.0d0
                     END IF
                  END IF IF_Surf1
!
! ............... Assignments
!
                  ExclZone(k)%sign1  = sign1
                  ExclZone(k)%expon1 = exponent1
!
                  ExclZone(k)%Equ1CoeffA(0:i1) = BoundSurface1_EquCoeff_A(0:i1)
                  ExclZone(k)%Equ1CoeffB(0:i1) = BoundSurface1_EquCoeff_B(0:i1)
!
                  IF( ExclZone(k)%TypeEqu1(1:4) == 'Powe' ) THEN
                     ExclZone(k)%Equ1CoeffC(0:i1) = BoundSurface1_EquCoeff_C(0:i1)
                     ExclZone(k)%Equ1CoeffD(0:i1) = BoundSurface1_EquCoeff_D(0:i1)
                  ELSE
                     IF(i1 > 1) ExclZone(k)%Equ1CoeffC(1:i1-1) = BoundSurface1_EquCoeff_C(1:i1-1)
                  END IF
!
                  ExclZone(k)%L1Shift(1) = X1_shift
                  ExclZone(k)%L1Shift(2) = Y1_shift
                  ExclZone(k)%L1Shift(3) = Z1_shift
!
! >>>>>
!
 1500             IF_Inte2: IF( type_of_equation2(1:4) == 'INTE' .OR. type_of_equation2(1:4) == 'Inte' .OR. type_of_equation2(1:4) == 'inte' ) THEN
!
                     type_of_equation2(1:4) = 'Inte'
!
                     IntDat_match = .FALSE.
                     IF( Num_IntTables > 0 ) THEN
                         DO_m3: DO m = 1,Num_IntTables
                            IF( interpolation_data_file_name2 == IntTable(m)%FileName ) THEN
                               ExclZone(k)%IntTableNum2 = m
                               IntDat_match             = .TRUE.
                               EXIT
                            END IF
                         END DO DO_m3
                     END IF
!
                     IF( .NOT. IntDat_match ) THEN
!
                        m = Num_IntTables + 1
                        Num_IntTables = m
                        ExclZone(k)%IntTableNum2 = Num_IntTables
!
! ..................... Reading info describing the structure of interpolation Table 2
!
                        READ (UNIT = *, NML = Irregular_ExclZone_IntTable2, IOSTAT = ier )
!
                        IF(ier /= 0) THEN
                          WRITE(UNIT = *, FMT = 6601) '<Irregular_ExclZone_IntTable2>', n
                           STOP
                        END IF
!
                        CALL Read_Tabular_Data( table_number = m, file_name = interpolation_data_file_name2, number_of_rows = number_of_rows, &
     &                                          number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                          read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                          RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                          interpolation_search_radius = interpolation_search_radius,  &
     &                                          num_interrogated_points     = num_interrogated_points )

!
                     END IF
!
                     CYCLE DO_Exclusions
!
                  ELSE IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' )  THEN
!
                     ExclZone(k)%thick2 = thickness2
!
                     ExclZone(k)%VertOrStrat2 = vertical_or_stratigraphic2
                     ExclZone(k)%RefSurf2     = reference_surface2
!
                     CYCLE DO_Exclusions
!
                  END IF IF_Inte2
!
! ............... Initializations
!
                  BoundSurface2_EquCoeff_A = 0.0d0
                  BoundSurface2_EquCoeff_B = 0.0d0
                  BoundSurface2_EquCoeff_C = 0.0d0
                  BoundSurface2_EquCoeff_D = 0.0d0
!
! ............... Reading the data for surface 2
!
                  READ (UNIT = *, NML = Irregular_ExclZone_Surf2, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_ExclZone_Surf2>',n
                     STOP
                  END IF
!
! .................. Checking consistency of orders of equation
!
                  IF( ( equation_order_of_bounding_surface1 < 0 ) .OR. ( equation_order_of_bounding_surface2 < 0 ) ) THEN
                     WRITE( UNIT = *, FMT = 6604 ) n
                     STOP
                  END IF
!
                  ExclZone(k)%OrderEquSurf2 = equation_order_of_bounding_surface2
!
                  i2 = equation_order_of_bounding_surface2
!
                  ALLOCATE ( ExclZone(k)%Equ2CoeffA(0:i2), ExclZone(k)%Equ2CoeffB(0:i2) )
!
                  ExclZone(k)%Equ2CoeffA = 0.0d0
                  ExclZone(k)%Equ2CoeffB = 0.0d0
!
                  IF_Surf2: IF( ExclZone(k)%TypeEqu2(1:4) == 'Powe' ) THEN
                     ALLOCATE ( ExclZone(k)%Equ2CoeffC(0:i2), ExclZone(k)%Equ2CoeffD(0:i2) )
                     ExclZone(k)%Equ2CoeffC = 0.0d0
                     ExclZone(k)%Equ2CoeffD = 0.0d0
                  ELSE
                     IF(i2 > 1 ) THEN
                        ALLOCATE ( ExclZone(k)%Equ2CoeffC(1:Max(1,i2)) )
                        ExclZone(k)%Equ2CoeffC = 0.0d0
                     END IF
                  END IF IF_Surf2
!
! ............... Assignments
!
                  ExclZone(k)%sign2  = sign2
                  ExclZone(k)%expon2 = exponent2
!
                  ExclZone(k)%Equ2CoeffA(0:i2) = BoundSurface2_EquCoeff_A(0:i2)
                  ExclZone(k)%Equ2CoeffB(0:i2) = BoundSurface2_EquCoeff_B(0:i2)
!
                  IF( ExclZone(k)%TypeEqu2(1:4) == 'Powe' ) THEN
                     ExclZone(k)%Equ2CoeffC(0:i2) = BoundSurface2_EquCoeff_C(0:i2)
                     ExclZone(k)%Equ2CoeffD(0:i2) = BoundSurface2_EquCoeff_D(0:i2)
                  ELSE
                     IF(i2 > 1) ExclZone(k)%Equ2CoeffC(1:i2-1) = BoundSurface2_EquCoeff_C(1:i2-1)
                  END IF
!
                  ExclZone(k)%L2Shift(1) = X2_shift
                  ExclZone(k)%L2Shift(2) = Y2_shift
                  ExclZone(k)%L2Shift(3) = Z2_shift
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For Cylindrical coodinates
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               ELSE
!
! ............... Irregular shapes
!
                  IF( NXMax >  1 .AND. NZMax == 1 ) THEN
                     WRITE( UNIT = *, FMT = 6606 ) 'radial', n
                     STOP
                  END IF
!
! ............... Initializations
!
                  BoundSurface1_EquCoeff_A = 0.0d0
                  BoundSurface1_EquCoeff_B = 0.0d0
                  BoundSurface2_EquCoeff_A = 0.0d0
                  BoundSurface2_EquCoeff_B = 0.0d0
!
! ............... Reading the data for surface 1
!
                  READ (UNIT = *, NML = Irregular_ExclZone_Surf1, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_ExclZone_Surf1>',n
                     STOP
                  END IF
!
! ............... Reading the data for surface 2
!
                  READ (UNIT = *, NML = Irregular_ExclZone_Surf2, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_ExclZone_Surf2>',n
                     STOP
                  END IF
!
! ............... Checking consistency of orders of equation
!
                  IF( ( equation_order_of_bounding_surface1 < 0 ) .OR. ( equation_order_of_bounding_surface2 < 0 ) ) THEN
                     WRITE( UNIT = *, FMT = 6604 ) n
                     STOP
                  END IF
!
                  ExclZone(k)%DependentVar = dependent_variable_of_surfaces
!
                  ExclZone(k)%OrderEquSurf1 = equation_order_of_bounding_surface1
                  ExclZone(k)%OrderEquSurf2 = equation_order_of_bounding_surface2
!
                  i1 = equation_order_of_bounding_surface1
                  ALLOCATE ( ExclZone(k)%Equ1CoeffA(0:i1) )
!
                  ExclZone(k)%Equ1CoeffA = 0.0d0
!
                  IF( ExclZone(k)%TypeEqu1(1:4) == 'Powe' ) THEN
                     ALLOCATE ( ExclZone(k)%Equ1CoeffB(0:i1) )
                     ExclZone(k)%Equ1CoeffB = 0.0d0
                  END IF
!
                  i2 = equation_order_of_bounding_surface2
                  ALLOCATE ( ExclZone(k)%Equ2CoeffA(0:i2) )
!
                  ExclZone(k)%Equ2CoeffA = 0.0d0
!
                  IF( ExclZone(k)%TypeEqu2(1:4) == 'Powe' ) THEN
                     ALLOCATE ( ExclZone(k)%Equ2CoeffB(0:i2) )
                     ExclZone(k)%Equ2CoeffB = 0.0d0
                  END IF
!
! ............... Assignments
!
                  ExclZone(k)%LMin(1) = R_min
                  ExclZone(k)%LMin(3) = Z_min
!
                  ExclZone(k)%LMax(1) = R_max
                  ExclZone(k)%LMax(3) = Z_max
!
                  ExclZone(k)%Equ1CoeffA(0:i1) = BoundSurface1_EquCoeff_A(0:i1)
                  ExclZone(k)%Equ2CoeffA(0:i2) = BoundSurface2_EquCoeff_A(0:i2)
!
                  IF( ExclZone(k)%TypeEqu1(1:4) == 'Powe' ) ExclZone(k)%Equ1CoeffB(0:i1) = BoundSurface1_EquCoeff_B(0:i1)
                  IF( ExclZone(k)%TypeEqu2(1:4) == 'Powe' ) ExclZone(k)%Equ2CoeffB(0:i2) = BoundSurface2_EquCoeff_B(0:i2)
!
                  ExclZone(k)%L1Shift(1) = R1_shift
                  ExclZone(k)%L1Shift(3) = Z1_shift
!
                  ExclZone(k)%L2Shift(1) = R2_shift
                  ExclZone(k)%L2Shift(3) = Z2_shift
!
               END IF IF_Coord4
!
!
!
            END DO DO_Exclusions
!
! <<<<<<<<<
! <<<<<<<<<
! <<<<<<<<<
!
            IF( k /= TotNum_ExclZones ) THEN
               WRITE(UNIT = *, FMT = 6705) k, TotNum_ExclZones, number_of_exclusion_zones, number_of_periodic_ExclZones, total_number_periodic_ExclZones
               STOP
            END IF
!
!***********************************************************************
!*                                                                     *
!*        Convert the exclusion zone ranges into SI units (m)          *
!*                                                                     *
!***********************************************************************
!
! ......... Conversion into METERS from INCHES
!
            FORALL (n=2:TotNum_ExclZones, ExclZone(k)%units == 'IN' .OR. ExclZone(k)%units == 'in' .OR. ExclZone(k)%units == 'In')
!
               ExclZone(n)%LMin(1:3) = ExclZone(n)%LMin(1:3) * 2.54d-2
               ExclZone(n)%LMax(1:3) = ExclZone(n)%LMax(1:3) * 2.54d-2
!
               ExclZone(n)%CylBase1Coord(1:3) = ExclZone(n)%CylBase1Coord(1:3) * 2.54d-2
               ExclZone(n)%CylBase2Coord(1:3) = ExclZone(n)%CylBase2Coord(1:3) * 2.54d-2
!
               ExclZone(n)%CylRmin = ExclZone(n)%CylRmin * 2.54d-2
               ExclZone(n)%CylRmax = ExclZone(n)%CylRmax * 2.54d-2
!
               ExclZone(n)%SphereCenterCoord(1:3) = ExclZone(n)%SphereCenterCoord(1:3) * 2.54d-2
!
               ExclZone(n)%SphRmin = ExclZone(n)%SphRmin * 2.54d-2
               ExclZone(n)%SphRmax = ExclZone(n)%SphRmax * 2.54d-2
!
               ExclZone(n)%L1Shift(1:3) = ExclZone(n)%L1Shift(1:3) * 2.54d-2
               ExclZone(n)%L2Shift(1:3) = ExclZone(n)%L2Shift(1:3) * 2.54d-2
!
               ExclZone(n)%EllBase1Coord(1:3) = ExclZone(n)%EllBase1Coord(1:3) * 2.54d-2
               ExclZone(n)%EllBase1Coord(1:3) = ExclZone(n)%EllBase1Coord(1:3) * 2.54d-2
!
               ExclZone(n)%LAxis = ExclZone(n)%LAxis * 2.54d-2
               ExclZone(n)%SAxis = ExclZone(n)%SAxis * 2.54d-2
!
               ExclZone(n)%Base1LAxis = ExclZone(n)%Base1LAxis * 2.54d-2
               ExclZone(n)%Base1SAxis = ExclZone(n)%Base1SAxis * 2.54d-2
               ExclZone(n)%Base2LAxis = ExclZone(n)%Base2LAxis * 2.54d-2
               ExclZone(n)%Base2SAxis = ExclZone(n)%Base2SAxis * 2.54d-2
!
               ExclZone(n)%thick0 = ExclZone(n)%thick0 * 2.54d-2
               ExclZone(n)%thick1 = ExclZone(n)%thick1 * 2.54d-2
               ExclZone(n)%thick2 = ExclZone(n)%thick2 * 2.54d-2
!
            END FORALL
!
! ......... Conversion into METERS from FEET
!
            FORALL (n=2:TotNum_ExclZones, ExclZone(n)%units == 'FT' .OR. ExclZone(n)%units == 'ft' .OR. ExclZone(n)%units == 'Ft')
!
               ExclZone(n)%LMin(1:3) = ExclZone(n)%LMin(1:3) * 3.038d-1
               ExclZone(n)%LMax(1:3) = ExclZone(n)%LMax(1:3) * 3.038d-1
!
               ExclZone(n)%CylBase1Coord(1:3) = ExclZone(n)%CylBase1Coord(1:3) * 3.038d-1
               ExclZone(n)%CylBase2Coord(1:3) = ExclZone(n)%CylBase2Coord(1:3) * 3.038d-1
!
               ExclZone(n)%CylRmin = ExclZone(n)%CylRmin * 3.038d-1
               ExclZone(n)%CylRmax = ExclZone(n)%CylRmax * 3.038d-1
!
               ExclZone(n)%SphereCenterCoord(1:3) = ExclZone(n)%SphereCenterCoord(1:3) * 3.038d-1
!
               ExclZone(n)%SphRmin = ExclZone(n)%SphRmin * 3.038d-1
               ExclZone(n)%SphRmax = ExclZone(n)%SphRmax * 3.038d-1
!
               ExclZone(n)%L1Shift(1:3) = ExclZone(n)%L1Shift(1:3) * 3.038d-1
               ExclZone(n)%L2Shift(1:3) = ExclZone(n)%L2Shift(1:3) * 3.038d-1
!
               ExclZone(n)%EllBase1Coord(1:3) = ExclZone(n)%EllBase1Coord(1:3) * 3.038d-1
               ExclZone(n)%EllBase1Coord(1:3) = ExclZone(n)%EllBase1Coord(1:3) * 3.038d-1
!
               ExclZone(n)%LAxis = ExclZone(n)%LAxis * 3.038d-1
               ExclZone(n)%SAxis = ExclZone(n)%SAxis * 3.038d-1
!
               ExclZone(n)%Base1LAxis = ExclZone(n)%Base1LAxis * 3.038d-1
               ExclZone(n)%Base1SAxis = ExclZone(n)%Base1SAxis * 3.038d-1
               ExclZone(n)%Base2LAxis = ExclZone(n)%Base2LAxis * 3.038d-1
               ExclZone(n)%Base2SAxis = ExclZone(n)%Base2SAxis * 3.038d-1
!
               ExclZone(n)%thick0 = ExclZone(n)%thick0 * 3.038d-1
               ExclZone(n)%thick1 = ExclZone(n)%thick1 * 3.038d-1
               ExclZone(n)%thick2 = ExclZone(n)%thick2 * 3.038d-1
!
            END FORALL
!
! ......... Conversion into METERS from KM
!
            FORALL (n=2:TotNum_ExclZones, ExclZone(n)%units == 'KM' .OR. ExclZone(n)%units == 'km' .OR. ExclZone(n)%units == 'Km')
!
               ExclZone(n)%LMin(1:3) = ExclZone(n)%LMin(1:3) * 1.0d3
               ExclZone(n)%LMax(1:3) = ExclZone(n)%LMax(1:3) * 1.0d3
!
               ExclZone(n)%CylBase1Coord(1:3) = ExclZone(n)%CylBase1Coord(1:3) * 1.0d3
               ExclZone(n)%CylBase2Coord(1:3) = ExclZone(n)%CylBase2Coord(1:3) * 1.0d3
!
               ExclZone(n)%CylRmin = ExclZone(n)%CylRmin * 1.0d3
               ExclZone(n)%CylRmax = ExclZone(n)%CylRmax * 1.0d3
!
               ExclZone(n)%SphereCenterCoord(1:3) = ExclZone(n)%SphereCenterCoord(1:3) * 1.0d3
!
               ExclZone(n)%SphRmin = ExclZone(n)%SphRmin * 1.0d3
               ExclZone(n)%SphRmax = ExclZone(n)%SphRmax * 1.0d3
!
               ExclZone(n)%L1Shift(1:3) = ExclZone(n)%L1Shift(1:3) * 1.0d3
               ExclZone(n)%L2Shift(1:3) = ExclZone(n)%L2Shift(1:3) * 1.0d3
!
               ExclZone(n)%EllBase1Coord(1:3) = ExclZone(n)%EllBase1Coord(1:3) * 1.0d3
               ExclZone(n)%EllBase1Coord(1:3) = ExclZone(n)%EllBase1Coord(1:3) * 1.0d3
!
               ExclZone(n)%LAxis = ExclZone(n)%LAxis * 1.0d3
               ExclZone(n)%SAxis = ExclZone(n)%SAxis * 1.0d3
!
               ExclZone(n)%Base1LAxis = ExclZone(n)%Base1LAxis * 1.0d3
               ExclZone(n)%Base1SAxis = ExclZone(n)%Base1SAxis * 1.0d3
               ExclZone(n)%Base2LAxis = ExclZone(n)%Base2LAxis * 1.0d3
               ExclZone(n)%Base2SAxis = ExclZone(n)%Base2SAxis * 1.0d3
!
               ExclZone(n)%thick0 = ExclZone(n)%thick0 * 1.0d3
               ExclZone(n)%thick1 = ExclZone(n)%thick1 * 1.0d3
               ExclZone(n)%thick2 = ExclZone(n)%thick2 * 1.0d3
!
            END FORALL
!
! ......... Conversion into METERS from MM
!
            FORALL (n=2:TotNum_ExclZones, ExclZone(n)%units == 'MM' .OR. ExclZone(n)%units == 'mm' .OR. ExclZone(n)%units == 'Mm')
!
               ExclZone(n)%LMin(1:3) = ExclZone(n)%LMin(1:3) * 1.0d-3
               ExclZone(n)%LMax(1:3) = ExclZone(n)%LMax(1:3) * 1.0d-3
!
               ExclZone(n)%CylBase1Coord(1:3) = ExclZone(n)%CylBase1Coord(1:3) * 1.0d-3
               ExclZone(n)%CylBase2Coord(1:3) = ExclZone(n)%CylBase2Coord(1:3) * 1.0d-3
!
               ExclZone(n)%CylRmin = ExclZone(n)%CylRmin * 1.0d-3
               ExclZone(n)%CylRmax = ExclZone(n)%CylRmax * 1.0d-3
!
               ExclZone(n)%SphereCenterCoord(1:3) = ExclZone(n)%SphereCenterCoord(1:3) * 1.0d-3
!
               ExclZone(n)%SphRmin = ExclZone(n)%SphRmin * 1.0d-3
               ExclZone(n)%SphRmax = ExclZone(n)%SphRmax * 1.0d-3
!
               ExclZone(n)%L1Shift(1:3) = ExclZone(n)%L1Shift(1:3) * 1.0d-3
               ExclZone(n)%L2Shift(1:3) = ExclZone(n)%L2Shift(1:3) * 1.0d-3
!
               ExclZone(n)%EllBase1Coord(1:3) = ExclZone(n)%EllBase1Coord(1:3) * 1.0d-3
               ExclZone(n)%EllBase1Coord(1:3) = ExclZone(n)%EllBase1Coord(1:3) * 1.0d-3
!
               ExclZone(n)%LAxis = ExclZone(n)%LAxis * 1.0d-3
               ExclZone(n)%SAxis = ExclZone(n)%SAxis * 1.0d-3
!
               ExclZone(n)%Base1LAxis = ExclZone(n)%Base1LAxis * 1.0d-3
               ExclZone(n)%Base1SAxis = ExclZone(n)%Base1SAxis * 1.0d-3
               ExclZone(n)%Base2LAxis = ExclZone(n)%Base2LAxis * 1.0d-3
               ExclZone(n)%Base2SAxis = ExclZone(n)%Base2SAxis * 1.0d-3
!
               ExclZone(n)%thick0 = ExclZone(n)%thick0 * 1.0d-3
               ExclZone(n)%thick1 = ExclZone(n)%thick1 * 1.0d-3
               ExclZone(n)%thick2 = ExclZone(n)%thick2 * 1.0d-3
!
            END FORALL
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Define_Exclusion_Zones 1.0 ............. 18 January   2015',6X,'Defining exclusion zones in the grid')
!
 6010 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The shape of the exclusion zone # ',I3.3,' is "',A,'": Unknown/Unavailable option'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6050 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input X_min =',ES12.5,' of the exclusion zone # ',I3.3,' is larger than (or equal to) the X_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6051 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input Y_min =',ES12.5,' of the exclusion zone # ',I3.3,' is larger than (or equal to) the Y_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6052 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input Z_min =',ES12.5,' of the exclusion zone # ',I3.3,' is larger than (or equal to) the Z_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6053 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'In the cylindrical exclusion zone # ',I3.3,', the centers of the two bases coincide ',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6054 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <cylinder_Rmin> =',ES12.5,' of the cylindrical exclusion zone # ',I3.3,  &
     &          ' is larger than (or equal to) the <cylinder_Rmax> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6055 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <sphere_Rmin> =',ES12.5,' of the spherical exclusion zone # ',I3.3,  &
     &          ' is larger than (or equal to) the <sphere_Rmax> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6056 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <Base1_Short_Axis> =',ES12.5,' of the elliptical exclusion zone # ',I3.3,  &
     &          ' is larger than (or equal to) the <Base1_Long_Axis> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6058 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <Base2_Short_Axis> =',ES12.5,' of the elliptical exclusion zone # ',I3.3,  &
     &          ' is larger than (or equal to) the <Base2_Long_Axis> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6059 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <plane_of_ellipse_bases> =',A,' of the elliptical exclusion zone # ',I3.3,' is NOT an available option',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6101 FORMAT(T5,'Memory allocation to array <',A,'> in subroutine <Define_Exclusion_Zones> was successful')
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <Define_Exclusion_Zones> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6105 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The subdomain shape (read by "exclusion_zone_shape" = ',a11,' in subroutine <Define_Exclusion_Zones>) is unavailable',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6401 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: The exclusion zone #',i3.3,' has an irregular shape that is described by ',/, &
     &       T10,'                                         an unknown/unavailable type of equation (= "',A,'")',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6402 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: The exclusion zone #',i3.3,' is defined by two irregular surfaces of type "FIXED" but the ' ,/,  &
     &       T10,'                                         name of the needed reference interpolation file is not defined',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED     !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6515 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: The dataset/namelist "',A,'" must be ended by the "<<<" descriptor'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6520 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: It is not possible to have a rectangular exclusion zone in a cylindrical coordinate system ',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6600 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: There is a problem reading the namelist <',A,'>',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6601 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: There is a problem reading the namelist ',A,' in exclusion zone #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6602 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: There is a conflict between the various radii options that describe the ', /,&
     &       T10,'                                        ',A,' shape of exclusion zone #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6604 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: At least one the orders of the equations describing the bounding surfaces 1 and/or 2 ', /,&
     &       T10,'                                         of exclusion zone #',i3.3,' is < 0',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6605 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: The exclusion zone #',i3.3,' has an irregular shape but the types of the ',/, &
     &       T10,'     equations describing the bounding surfaces 1 and/or 2 of the exclusion zone are non-blanks',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6606 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: This ',A,' system has 1 active dimension - it is not possible for the ', /, &
     &       T10,'                                         exclusion zone #',i3.3,' to be irregularly shaped',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6608 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: In this ',A,' system, one or more of the dependent variables defining the bounding surfaces 1 and 2', /,&
     &       T10,'                                         of exclusion zone #',i3.3,' is not among the possible options (',A,')',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6701 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: The shape of exclusion zone #',i3.3,' is "PERIODIC", but the coordinate system = "',A,'"' /,&
     &       T10,'                                         conflicts with the <axis_of_perodicity> = "',A1,'"',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6702 FORMAT(//,20('ERROR-'),//,   &
     &       T10,'>>>  Procedure <Define_Exclusion_Zones>: The shape of exclusion zone #',i3.3,' is "PERIODIC", but the axis of periodicity = "',A1,'"' /,&
     &       T10,'                                         is not among the available options ("X,", "R", "Y", or "Z")',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6705 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'>>>  Procedure <Define_Exclusion_Zones>: The total number of computed exclusion zones = ',i3.3,' does not match the number <TotNum_ExclZones> that is estimated from:',/, &
     &       T5,'                                         TotNum_ExclZones(=',i3.3,') = number_of_exclusion_zones(=',i3.3,') - number_of_periodic_ExclZones(=',i3.3,') + total_number_periodic_ExclZones(=',i3.3,')',/, &
     &       T5,'                                         Check all the related inputs', //, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Define_Exclusion_Zones>
!
!
            RETURN
!
         END SUBROUTINE Define_Exclusion_Zones
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
         SUBROUTINE Define_Inclusion_Zones
!
            USE MeshMaker_Data, ONLY: coordinates, NXMax, NYMax, NZMax
            USE Het_Region_Definition
            USE Grid_Generation_Parameters, ONLY: pi
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          Routine defining the included zones in the grid            *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
!
! -------------
! ......... Double precision variables
! -------------
!
            REAL(KIND = 8) :: D1, D2, D3
!
            REAL(KIND = 8) :: X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max
!
            REAL(KIND = 8) :: sphere_Rmax, sphere_Rmin, cylinder_Rmin, cylinder_Rmax, cylinder_radius1, cylinder_radius2
            REAL(KIND = 8) :: Base1_Long_Axis, Base2_Long_Axis, Base1_Short_Axis, Base2_Short_Axis, Long_Axis, Short_Axis, ellipse_long_axis_angle
!
            REAL(KIND = 8) :: X1_shift, Y1_shift, Z1_shift, R1_shift, X2_shift, Y2_shift, Z2_shift, R2_shift
!
            REAL(KIND = 8) :: exponent1, sign1, exponent2, sign2, emissivity
!
            REAL(KIND = 8) :: location_of_1st_periodic_occurrence, width_of_periodic_InclZone, period_of_occurrence
            REAL(KIND = 8) :: period, L_first, L_width, thickness0, thickness1, thickness2, interpolation_search_radius
!
! -------------
! ......... Double precision arrays
! -------------
!
            REAL(KIND = 8), DIMENSION(3) :: sphere_center_coordinates, CylBase1_center_coordinates, CylBase2_center_coordinates
            REAL(KIND = 8), DIMENSION(3) :: Base1_center_coordinates,  Base2_center_coordinates
!
            REAL(KIND = 8), DIMENSION(0:10) :: BoundSurface1_EquCoeff_A, BoundSurface1_EquCoeff_B, BoundSurface1_EquCoeff_C, BoundSurface1_EquCoeff_D
            REAL(KIND = 8), DIMENSION(0:10) :: BoundSurface2_EquCoeff_A, BoundSurface2_EquCoeff_B, BoundSurface2_EquCoeff_C, BoundSurface2_EquCoeff_D
!
! -------------
! ......... Integer variables
! -------------
!
            INTEGER :: ier, n, m, i1, i2, j, k, AxisNum
!
            INTEGER :: number_of_inclusion_zones, number_of_periodic_InclZones, total_number_periodic_InclZones, number_of_periodic_occurrences
!
            INTEGER :: equation_order_of_bounding_surface1, equation_order_of_bounding_surface2, num_interrogated_points
!
            INTEGER :: number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2
!
            INTEGER, ALLOCATABLE, DIMENSION(:) :: frequency_of_occurrences, temp
!
! -------------
! ......... CHARACTER variables
! -------------
!
            CHARACTER(LEN =  1) :: axis_of_periodicity, vertical_or_stratigraphic1, vertical_or_stratigraphic2, reference_surface1, reference_surface2
            CHARACTER(LEN =  2) :: length_units, plane_of_ellipse_bases
            CHARACTER(LEN =  3) :: inclusion_zone_type
            CHARACTER(LEN =  5) :: inclusion_zone_name
            CHARACTER(LEN =  6) :: first_part
            CHARACTER(LEN =  8) :: interpolation_data_file_name1, interpolation_data_file_name2, ref_interpolation_data_file
            CHARACTER(LEN = 12) :: inclusion_zone_shape, type_of_equation1, type_of_equation2
            CHARACTER(LEN = 50) :: read_data_format
!
            CHARACTER(LEN = 1) :: dependent_variable_of_surfaces
!
! -------------
! ......... LOGICAL variables
! -------------
!
            LOGICAL :: IntDat_match, read_data_by_row
!
! -------
! ... Namelists
! -------
!
      NAMELIST/ Inclusion_Zones / number_of_inclusion_zones, number_of_periodic_InclZones, total_number_periodic_InclZones
!
      NAMELIST/ InclZone_GeneralInfo / inclusion_zone_name, inclusion_zone_shape, inclusion_zone_type, length_units, emissivity
!
      NAMELIST/ Rectangular_InclZone / X_min, Y_min, Z_min, X_max, Y_max, Z_max
!
      NAMELIST/ Periodic_InclZone    / X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max,           &
     &                                 number_of_periodic_occurrences, axis_of_periodicity,             &
     &                                 location_of_1st_periodic_occurrence, width_of_periodic_InclZone, &
     &                                 period_of_occurrence
!
      NAMELIST/ Cylindrical_InclZone / CylBase1_center_coordinates, CylBase2_center_coordinates, &
     &                                 cylinder_Rmin, cylinder_Rmax, cylinder_radius1, cylinder_radius2
!
      NAMELIST/ Spherical_InclZone   / sphere_center_coordinates, sphere_Rmax, sphere_Rmin
!
      NAMELIST/ Elliptical_InclZone  / Base1_center_coordinates, Base2_center_coordinates, ellipse_long_axis_angle, plane_of_ellipse_bases, &
     &                                 Base1_Long_Axis, Base2_Long_Axis, Base1_Short_Axis, Base2_Short_Axis, Long_Axis, Short_Axis
!
      NAMELIST/ Irregular_InclZone   /  dependent_variable_of_surfaces, type_of_equation1, type_of_equation2, &
     &                                  X_min, R_min, Y_min, Z_min, X_max, R_max, Y_max, Z_max,               &
     &                                  interpolation_data_file_name1, interpolation_data_file_name2,         &
     &                                  vertical_or_stratigraphic1, vertical_or_stratigraphic2, reference_surface1, reference_surface2, &
     &                                  thickness0, thickness1, thickness2, ref_interpolation_data_file

!
      NAMELIST/ Irregular_InclZone_Surf1 / equation_order_of_bounding_surface1,                      &
     &                                     X1_shift, Y1_shift, Z1_shift, R1_shift, exponent1, sign1, &
     &                                     BoundSurface1_EquCoeff_A, BoundSurface1_EquCoeff_B, BoundSurface1_EquCoeff_C, BoundSurface1_EquCoeff_D
!
      NAMELIST/ Irregular_InclZone_Surf2 / equation_order_of_bounding_surface2,                      &
     &                                     X2_shift, Y2_shift, Z2_shift, R2_shift, exponent2, sign2, &
     &                                     BoundSurface2_EquCoeff_A, BoundSurface2_EquCoeff_B, BoundSurface2_EquCoeff_C, BoundSurface2_EquCoeff_D
!
      NAMELIST/ Irregular_InclZone_IntTable1 / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                         read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
      NAMELIST/ Irregular_InclZone_IntTable2 / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                         read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
      NAMELIST/ Irregular_InclZone_RefTable  / number_of_rows, number_of_columns, RowCol_DepVariable, RowCol_IndVariable_1, RowCol_IndVariable_2, &
     &                                         read_data_by_row, read_data_format, interpolation_search_radius, num_interrogated_points
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of <Define_Inclusion_Zones>
!
!
            WRITE(*,6000)
!
! ......... Initialization
!
            number_of_periodic_InclZones    = 0
            total_number_periodic_InclZones = 0
!
! ------------
! ......... Read the number of heterogeneous domains
! ------------
!
            READ (UNIT = *, NML = Inclusion_Zones, IOSTAT = ier )
!
! ......... Stop if there is a problem reading the namelist
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6600) '<Inclusion_Zones>'
               STOP
            END IF
!
            ALLOCATE( frequency_of_occurrences(number_of_periodic_InclZones) )
!
            frequency_of_occurrences = 0    ! ... Initialization
!
            Num_InclZones    = number_of_inclusion_zones
            TotNum_InclZones = number_of_inclusion_zones - number_of_periodic_InclZones + total_number_periodic_InclZones
!
            IF(Num_InclZones <= 0) RETURN
!
! ......... Re-dimension the <MediaSequenceNumber> array; needed if the inclusions increase the number of media
!
            ALLOCATE( temp(TotNum_HetRegions) )
!
            temp = MediaSequenceNumber ! ... Array operation
!
            DEALLOCATE( MediaSequenceNumber )
            ALLOCATE  ( MediaSequenceNumber( TotNum_HetRegions + TotNum_InclZones ) )
!
            MediaSequenceNumber(1:TotNum_HetRegions) = temp(1:TotNum_HetRegions)
!
            DEALLOCATE( temp )
!
! ------------
! ......... Allocate memory to the arrays defining the inclusion zones
! ------------
!
            ALLOCATE(InclZone(1:TotNum_InclZones), InclZone_SequenceNumber(1:TotNum_InclZones), STAT=ier)
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'InclZone'
            ELSE
               WRITE(*,6102) 'InclZone'
               STOP
            END IF
!
! ......... Storing temporarily table names, before re-dimensioning the grid/interpolation data sets
!
            i1 = Num_IntTables
            i2 = Num_IntTables + 2 * Num_InclZones
!
            IF( i1 > 0 ) temp_IntTable(1:i1) = IntTable(1:i1)
!
            DEALLOCATE( IntTable )
            ALLOCATE( IntTable(i2), STAT=ier )
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'IntTable'
            ELSE
               WRITE(*,6102) 'IntTable'
               STOP
            END IF
!
            IF( i1 > 0 ) IntTable(1:i1) = temp_IntTable
!
            DEALLOCATE( temp_IntTable )
            ALLOCATE( temp_IntTable(i2), STAT=ier )
!
! ......... Print explanatory comments
!
            IF(ier == 0) THEN
               WRITE(*,6101) 'temp_IntTable'
            ELSE
               WRITE(*,6102) 'temp_IntTable'
               STOP
            END IF
!
!***********************************************************************
!*                                                                     *
!*      Read the specifics of the remaining media inclusion zones      *
!*                                                                     *
!***********************************************************************
!
            k = 0
            m = 0
!
            DO_Inclusions: DO n = 1, Num_InclZones + 1
!
! -----------------
! ............ Read inclusion zone specifics for Num_InclZones > 1
! -----------------
!
               IF( n == Num_InclZones + 1 ) THEN
!
                  READ(UNIT = *, FMT = '(A3)') first_part
!
                  IF( first_part /= '<<<' ) THEN
                     WRITE(UNIT = *, FMT = 6515) first_part
                     STOP
                  ELSE
                     EXIT DO_Inclusions
                  END IF
!
               END IF
!
! ............ Advance the counter
!
               k = k + 1
!
! ............ Initializations - Namelist components
!
               inclusion_zone_name  = '     '
               inclusion_zone_shape = '           '
               inclusion_zone_type  = '   '
!
               length_units = 'm'
               emissivity   = 0.0d0
!
               InclZone(k)%name  = '     '
               InclZone(k)%shape = '           '
               InclZone(k)%id    = '   '
               InclZone(k)%units = 'm'
!
               InclZone(k)%emissivity = 0.0d0
!
! -----------------
! ............ Read general info on the inclusion zone
! -----------------
!
               READ (UNIT = *, NML = InclZone_GeneralInfo, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6601) '<InclZone_GeneralInfo>',n
                  STOP
               END IF
!
               SELECT CASE( inclusion_zone_shape(1:1) )
               CASE('R', 'r', 'C', 'c', 'S', 's', 'E', 'e', 'P', 'p', 'I', 'i' )
                  CONTINUE
               CASE DEFAULT
                  WRITE( UNIT = *, FMT = 6010 ) n, inclusion_zone_shape
                  STOP
               END SELECT
!
               IF( ( coordinates(1:2) == 'CY' .OR. coordinates(1:4) == 'Cy' .OR. coordinates(1:4) == 'cy' ) .AND. &
     &             ( inclusion_zone_shape(1:1) == 'R' .OR. inclusion_zone_shape(1:1) == 'r' ) )  &
     &         THEN
                  WRITE(*,6520) inclusion_zone_shape
                  STOP
               END IF
!
! ............ Assignment of the namelist values
!
               InclZone(k)%name   = inclusion_zone_name
               InclZone(k)%id     = inclusion_zone_type
               InclZone(k)%shape  = inclusion_zone_shape
               InclZone(k)%units  = length_units
!
               InclZone(k)%emissivity = emissivity
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    Reading the shape-specific data from the remaining namelists     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
               SELECT CASE(inclusion_zone_shape(1:1))
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Rectangular inclusion zones
! >>>>>>>>>>>>
!
               CASE( 'R', 'r' )
!
! ............... Initializations
!
                  X_min = 0.0d0
                  Y_min = 0.0d0
                  Z_min = 0.0d0
                  X_max = 0.0d0
                  Y_max = 0.0d0
                  Z_max = 0.0d0
!
                  InclZone(k)%LMin = 0.0d0     ! ... Whole array operations
                  InclZone(k)%LMax = 0.0d0
!
! ............... Reading the rectangle-related data using a namelist
!
                  READ (UNIT = *, NML = Rectangular_InclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Rectangular_InclZone>', n
                     STOP
                  END IF
!
                  InclZone_SequenceNumber(k) = n
!
                  InclZone(k)%LMin(1) = X_min
                  InclZone(k)%LMin(2) = Y_min
                  InclZone(k)%LMin(3) = Z_min
!
                  InclZone(k)%LMax(1) = X_max
                  InclZone(k)%LMax(2) = Y_max
                  InclZone(k)%LMax(3) = Z_max
!
! ............... Ensuring that the LMin of the range is not >= of the LMax
!
                  IF( InclZone(k)%LMin(1) >= InclZone(k)%LMax(1) ) THEN
                     WRITE(*,6050) InclZone(k)%LMin(1), n, InclZone(k)%LMax(1)
                     STOP
                  END IF
!
                  IF( InclZone(k)%LMin(2) >= InclZone(k)%LMax(2) ) THEN
                     WRITE(*,6051) InclZone(k)%LMin(2), n, InclZone(k)%LMax(2)
                     STOP
                  END IF
!
                  IF( InclZone(k)%LMin(3) >= InclZone(k)%LMax(3) ) THEN
                     WRITE(*,6052)  InclZone(k)%LMin(3), n, InclZone(k)%LMax(3)
                     STOP
                  END IF
!
                  CYCLE DO_Inclusions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Cylindrical inclusion zones
! >>>>>>>>>>>>
!
               CASE( 'C', 'c' )
!
! ............... Initializations
!
                  CylBase1_center_coordinates = 0.0d0  ! ... Whole array operation
                  CylBase2_center_coordinates = 0.0d0
                  cylinder_Rmin    = 0.0d0
                  cylinder_Rmax    = 0.0d0
                  cylinder_radius1 = 0.0d0
                  cylinder_radius2 = 0.0d0
!
                  InclZone(k)%CylBase1Coord = 0.0d0    ! ... Whole array operation
                  InclZone(k)%CylBase1Coord = 0.0d0
                  InclZone(k)%CylRmin   = 0.0d0
                  InclZone(k)%CylRmax   = 0.0d0
                  InclZone(k)%CylBase1R = 0.0d0
                  InclZone(k)%CylBase2R = 0.0d0
!
! ............... Reading the cylinder-related data using a namelist
!
                  READ (UNIT = *, NML = Cylindrical_InclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) 'Cylindrical_InclZone', n
                     STOP
                  END IF
!
                  InclZone_SequenceNumber(k) = n
!
! ............... For Cartesian coordinates
!
                  IF_Coord: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                          coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     InclZone(k)%CylBase1Coord = CylBase1_center_coordinates
                     InclZone(k)%CylBase2Coord = CylBase2_center_coordinates
!
                     IF( ( cylinder_Rmin > 1.0d-7 .OR. cylinder_Rmax > 1.0d-7 ) .AND. ( cylinder_radius1 < 1.0d-7 .AND. cylinder_radius2 < 1.0d-7 ) ) THEN
                        InclZone(k)%CylRmin = cylinder_Rmin
                        InclZone(k)%CylRmax = cylinder_Rmax
                     ELSE IF( ( cylinder_Rmin < 1.0d-7 .AND. cylinder_Rmax < 1.0d-7 ) .AND. ( cylinder_radius1 > 1.0d-7 .OR. cylinder_radius2 > 1.0d-7 ) ) THEN
                        InclZone(k)%CylBase1R = cylinder_radius1
                        InclZone(k)%CylBase2R = cylinder_radius2
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'cylindrical',n
                        STOP
                     END IF
!
                     IF( cylinder_radius1 /= 0.0d0 ) CYCLE DO_Inclusions
!
! .................. Ensuring that the two bases of the cylindrical inclusion zone do not coincide
!
                     D1 = InclZone(k)%CylBase1Coord(1) - InclZone(k)%CylBase2Coord(1)
                     D2 = InclZone(k)%CylBase1Coord(2) - InclZone(k)%CylBase2Coord(2)
                     D3 = InclZone(k)%CylBase1Coord(3) - InclZone(k)%CylBase2Coord(3)
!
                     IF( ABS(D1) < 1.0d-7 .AND. ABS(D2) < 1.0d-7 .AND. ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! .................. Ensuring that the Rmin of the cylindrical inclusion zone is not >= of the Rmax
!
                     IF( InclZone(k)%CylRmin >= InclZone(k)%CylRmax ) THEN
                       WRITE(*,6054) InclZone(k)%CylRmin, n, InclZone(k)%CylRmax
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     InclZone(k)%CylBase1Coord(3) = CylBase1_center_coordinates(3)
                     InclZone(k)%CylBase2Coord(3) = CylBase2_center_coordinates(3)
!
                     IF( ( cylinder_Rmax > 1.0d-7 ) .AND. ( cylinder_radius1 < 1.0d-7 .AND. cylinder_radius2 < 1.0d-7 ) ) THEN
                        InclZone(k)%CylRmin = cylinder_Rmin
                        InclZone(k)%CylRmax = cylinder_Rmax
                     ELSE IF( ( cylinder_Rmax < 1.0d-7 ) .AND. ( cylinder_radius1 > 1.0d-7 .AND. cylinder_radius2 > 1.0d-7 ) ) THEN
                        InclZone(k)%CylBase1R = cylinder_radius1
                        InclZone(k)%CylBase2R = cylinder_radius2
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'cylindrical',n
                        STOP
                     END IF
!
                     IF( cylinder_radius1 /= 0.0d0 ) CYCLE DO_Inclusions
!
! .................. Ensuring that the two bases of the cylindrical inclusion zone do not coincide
!
                     D3 = InclZone(k)%CylBase1Coord(3) - InclZone(k)%CylBase2Coord(3)
!
                     IF( ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! .................. Ensuring that the Rmin of the cylindrical inclusion zone is not >= of the Rmax
!
                     IF( InclZone(k)%CylRmin >= InclZone(k)%CylRmax ) THEN
                       WRITE(*,6054) InclZone(k)%CylRmin, n, InclZone(k)%CylRmax
                        STOP
                     END IF
!
                  END IF IF_Coord
!
                  CYCLE DO_Inclusions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Periodic inclusion zones
! >>>>>>>>>>>>
!
               CASE( 'P', 'p' )
!
! ............... Initializations
!
                  X_min = 0.0d0
                  R_min = 0.0d0
                  Y_min = 0.0d0
                  Z_min = 0.0d0
                  X_max = 0.0d0
                  R_max = 0.0d0
                  Y_max = 0.0d0
                  Z_max = 0.0d0
!
                  inclusion_zone_shape = 'periodic'
!
                  axis_of_periodicity                 = ' '
                  number_of_periodic_occurrences      = 0
                  location_of_1st_periodic_occurrence = 0.0d0
                  width_of_periodic_InclZone          = 0.0d0
                  period_of_occurrence                = 0.0d0
!
! ............... Reading the rectangle-related data using a namelist
!
                  READ (UNIT = *, NML = Periodic_InclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Periodic_InclZone>', n
                     STOP
                  END IF
!
                  m = m + 1
!
                  IF ( axis_of_periodicity == 'X' .AND. axis_of_periodicity == 'x' .AND.  &
     &                 axis_of_periodicity == 'R' .AND. axis_of_periodicity == 'r' .AND.  &
     &                 axis_of_periodicity == 'Y' .AND. axis_of_periodicity == 'y' .AND.  &
     &                 axis_of_periodicity == 'Z' .AND. axis_of_periodicity == 'z' )   &
     &            THEN
                     WRITE( UNIT = *, FMT = 6702 ) n, axis_of_periodicity
                     STOP
                  END IF
!
                  L_first = location_of_1st_periodic_occurrence
                  L_width = width_of_periodic_InclZone
                  period  = period_of_occurrence
!
                  frequency_of_occurrences(m) = number_of_periodic_occurrences
!
                  IF_Period: IF (axis_of_periodicity == 'X' .OR. axis_of_periodicity == 'x' ) THEN
!
                     AxisNum = 1
                     IF( coordinates(1:4) == 'CYLI' .OR. coordinates(1:4) == 'Cyli' .OR. coordinates(1:4) == 'cyli' ) THEN
                        WRITE( UNIT = *, FMT = 6701 ) n, coordinates, axis_of_periodicity
                        STOP
                     END IF
!
                  ELSE IF (axis_of_periodicity == 'R' .OR. axis_of_periodicity == 'r' ) THEN
!
                     AxisNum = 1
                     IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                   coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &               THEN

                        WRITE( UNIT = *, FMT = 6701 ) n, coordinates, axis_of_periodicity
                        STOP
                     END IF
!
                  ELSE IF (axis_of_periodicity == 'Y' .OR. axis_of_periodicity == 'y' ) THEN
!
                     AxisNum = 2
!
                  ELSE IF (axis_of_periodicity == 'Z' .OR. axis_of_periodicity == 'z' ) THEN
!
                     AxisNum = 3
!
                  END IF IF_Period
!
                  j = 0
!
 1000             IF( AxisNum == 1 ) THEN
                     X_min = L_first + j * period
                     X_max = X_min + L_width
                     R_min = X_min
                     R_max = X_max
                  ELSE IF( AxisNum == 2 ) THEN
                     Y_min = L_first + j * period
                     Y_max = Y_min + L_width
                  ELSE
                     Z_min = L_first + j * period
                     Z_max = Z_min + L_width
                  END IF
!
                  InclZone_SequenceNumber(k) = n
!
                  InclZone(k)%name  = inclusion_zone_name
                  InclZone(k)%shape = inclusion_zone_shape
                  InclZone(k)%id    = inclusion_zone_type
                  InclZone(k)%units = length_units
!
                  InclZone(k)%LMin = 0.0d0     ! ... Whole array operations
                  InclZone(k)%LMax = 0.0d0
!
                  InclZone(k)%LMin(1) = X_min
                  InclZone(k)%LMin(2) = Y_min
                  InclZone(k)%LMin(3) = Z_min
!
                  InclZone(k)%LMax(1) = X_max
                  InclZone(k)%LMax(2) = Y_max
                  InclZone(k)%LMax(3) = Z_max
!
! ............... Ensuring that the LMin of the range is not >= of the LMax
!
                  IF( InclZone(k)%LMin(1) >= InclZone(k)%LMax(1) ) THEN
                     WRITE(*,6050) InclZone(k)%LMin(1), n, InclZone(k)%LMax(1)
                     STOP
                  END IF
!
                  IF( InclZone(k)%LMin(2) >= InclZone(k)%LMax(2) ) THEN
                     WRITE(*,6051) InclZone(k)%LMin(2), n, InclZone(k)%LMax(2)
                     STOP
                  END IF
!
                  IF( InclZone(k)%LMin(3) >= InclZone(k)%LMax(3) ) THEN
                     WRITE(*,6052)  InclZone(k)%LMin(3), n, InclZone(k)%LMax(3)
                     STOP
                  END IF
!
                  IF ( j < frequency_of_occurrences(m) - 1 ) THEN
                     j = j + 1
                     k = k + 1
                     GO TO 1000
                  END IF
!
                  CYCLE DO_Inclusions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Elliptical inclusion zones
! >>>>>>>>>>>>
!
               CASE( 'E', 'e' )
!
! ............... Initializations
!
                  Base1_center_coordinates = 0.0d0  ! ... Whole array operation
                  Base2_center_coordinates = 0.0d0
                  ellipse_long_axis_angle  = 0.0d0
!
                  plane_of_ellipse_bases = '  '
!
                  Base1_Long_Axis  = 0.0d0
                  Base2_Long_Axis  = 0.0d0
                  Base1_Short_Axis = 0.0d0
                  Base2_Short_Axis = 0.0d0
                  Long_Axis  = 0.0d0
                  Short_Axis = 0.0d0
!
                  InclZone(k)%EllBase1Coord = 0.0d0    ! ... Whole array operation
                  InclZone(k)%EllBase1Coord = 0.0d0
                  InclZone(k)%EllipsePlane  = '  '
!
                  InclZone(k)%LAxisAngle = 0.0d0
                  InclZone(k)%Base1LAxis = 0.0d0
                  InclZone(k)%Base2LAxis = 0.0d0
                  InclZone(k)%Base1SAxis = 0.0d0
                  InclZone(k)%Base2SAxis = 0.0d0
                  InclZone(k)%LAxis = 0.0d0
                  InclZone(k)%SAxis = 0.0d0
!
! ............... Reading the cylinder-related data using a namelist
!
                  READ (UNIT = *, NML = Elliptical_InclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) 'Elliptical_InclZone', n
                     STOP
                  END IF
!
! ............... Ensuring that there is no conflict between the long and short axes lengths
!
                  IF( Base1_Short_Axis > Base1_Long_Axis ) THEN
                     WRITE(*,6056) Base1_Short_Axis, n, Base1_Long_Axis
                     STOP
                  END IF
!
                  IF( Base2_Short_Axis > Base2_Long_Axis ) THEN
                     WRITE(*,6058) Base2_Short_Axis, n, Base2_Long_Axis
                     STOP
                  END IF
!
                  SELECT CASE(plane_of_ellipse_bases)
                  CASE('XY','Xy','xY','xy','XZ','Xz','xZ','xz','YZ','Yz','yZ','yz')
                     CONTINUE
                  CASE DEFAULT
                     WRITE(UNIT = *, FMT = 6059) plane_of_ellipse_bases, n
                     STOP
                  END SELECT
!
                  InclZone_SequenceNumber(k) = n
!
! ............... For Cartesian coordinates
!
                  IF_Coord2: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                           coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     InclZone(k)%EllBase1Coord = Base1_center_coordinates
                     InclZone(k)%EllBase2Coord = Base2_center_coordinates
                     InclZone(k)%LAxisAngle    = pi * ellipse_long_axis_angle / 1.8d2
                     InclZone(k)%EllipsePlane  = plane_of_ellipse_bases
!
                     IF( ( Long_Axis > 1.0d-7 .AND. Short_Axis > 1.0d-7 ) .AND. &
     &                   ( Base1_Long_Axis < 1.0d-7 .AND. Base1_Long_Axis < 1.0d-7 .AND. Base1_Short_Axis < 1.0d-7 .AND. Base2_Short_Axis < 1.0d-7 ) ) &
     &               THEN
                        InclZone(k)%LAxis = Long_Axis
                        InclZone(k)%SAxis = Short_Axis
                     ELSE IF( ( Long_Axis < 1.0d-7 .AND. Short_Axis < 1.0d-7 ) .AND. &
     &                        ( Base1_Long_Axis > 1.0d-7 .AND. Base1_Long_Axis > 1.0d-7 .AND. Base1_Short_Axis > 1.0d-7 .AND. Base2_Short_Axis > 1.0d-7 ) ) &
     &               THEN
                        InclZone(k)%Base1LAxis = Base1_Long_Axis
                        InclZone(k)%Base2LAxis = Base2_Long_Axis
                        InclZone(k)%Base1SAxis = Base1_Short_Axis
                        InclZone(k)%Base2SAxis = Base2_Short_Axis
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'elliptical',n
                        STOP
                     END IF
!
! .................. Ensuring that the two bases of the elliptical inclusion zone do not coincide
!
                     D1 = InclZone(k)%EllBase1Coord(1) - InclZone(k)%EllBase2Coord(1)
                     D2 = InclZone(k)%EllBase1Coord(2) - InclZone(k)%EllBase2Coord(2)
                     D3 = InclZone(k)%EllBase1Coord(3) - InclZone(k)%EllBase2Coord(3)
!
                     IF( ABS(D1) < 1.0d-7 .AND. ABS(D2) < 1.0d-7 .AND. ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     InclZone(k)%EllBase1Coord(3) = Base1_center_coordinates(3)
                     InclZone(k)%EllBase2Coord(3) = Base2_center_coordinates(3)
                     InclZone(k)%EllipsePlane     = '  '
                     InclZone(k)%LAxisAngle       = 0.0d0
!
                     IF( ( Long_Axis > 1.0d-7 .AND. Short_Axis > 1.0d-7 ) .AND. &
     &                   ( Base1_Long_Axis < 1.0d-7 .AND. Base1_Long_Axis < 1.0d-7 .AND. Base1_Short_Axis < 1.0d-7 .AND. Base2_Short_Axis < 1.0d-7 ) ) &
     &               THEN
                        InclZone(k)%LAxis = Long_Axis
                        InclZone(k)%SAxis = Short_Axis
                     ELSE IF( ( Long_Axis < 1.0d-7 .AND. Short_Axis < 1.0d-7 ) .AND. &
     &                        ( Base1_Long_Axis > 1.0d-7 .AND. Base1_Long_Axis > 1.0d-7 .AND. Base1_Short_Axis > 1.0d-7 .AND. Base2_Short_Axis > 1.0d-7 ) ) &
     &               THEN
                        InclZone(k)%Base1LAxis = Base1_Long_Axis
                        InclZone(k)%Base2LAxis = Base2_Long_Axis
                        InclZone(k)%Base1SAxis = Base1_Short_Axis
                        InclZone(k)%Base2SAxis = Base2_Short_Axis
                     ELSE
                        WRITE( UNIT = *, FMT = 6602) 'elliptical',n
                        STOP
                     END IF
!
! .................. Ensuring that the two bases of the cylindrical inclusion zone do not coincide
!
                     D3 = InclZone(k)%CylBase1Coord(3) - InclZone(k)%CylBase2Coord(3)
!
                     IF( ABS(D3) < 1.0d-7 ) THEN
                        WRITE(*,6053) n
                        STOP
                     END IF
!
                  END IF IF_Coord2
!
                  CYCLE DO_Inclusions
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Spherical inclusion zones
! >>>>>>>>>>>>
!
               CASE( 'S', 's' )
!
! ............... Initializations
!
                  sphere_center_coordinates = 0.0d0    ! ... Whole array operation
                  sphere_Rmax = 0.0d0
                  sphere_Rmin = 0.0d0
!
                  InclZone(k)%SphereCenterCoord = 0.0d0  ! ... Whole array operation
                  InclZone(k)%SphRmin = 0.0d0
                  InclZone(k)%SphRmax = 0.0d0
!
! ............... Reading the sphere-related data using a namelist
!
                  READ (UNIT = *, NML = Spherical_InclZone, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) 'Spherical_InclZone', n
                     STOP
                  END IF
!
                  InclZone_SequenceNumber(k) = n
!
                  InclZone(k)%SphRmin = sphere_Rmin
                  InclZone(k)%SphRmax = sphere_Rmax
!
! ............... For Cartesian coordinates
!
                  IF_Coord3: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                           coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &            THEN

!
                     InclZone(k)%SphereCenterCoord = sphere_center_coordinates
!
! .................. Ensuring that the Rmin of the spherical inclusion zone is not >= of the Rmax
!
                     IF( InclZone(k)%SphRmin >= InclZone(k)%SphRmax ) THEN
                        WRITE(*,6055) InclZone(k)%SphRmin, n, InclZone(k)%CylRmax
                        STOP
                     END IF
!
! ............... For cylindrical coordinates
!
                  ELSE
!
                     InclZone(k)%SphereCenterCoord(3) = sphere_center_coordinates(3)
!
! .................. Ensuring that the Rmin of the spherical inclusion zone is not >= of the Rmax
!
                     IF( InclZone(k)%SphRmin >= InclZone(k)%SphRmax ) THEN
                        WRITE(*,6055) InclZone(k)%SphRmin, n, InclZone(k)%CylRmax
                        STOP
                     END IF
!
                  END IF IF_Coord3
!
                  CYCLE DO_Inclusions
!
               END SELECT
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Irregular inclusion zones
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               type_of_equation1 = '           '  ! ... Initialization
               type_of_equation2 = '           '
!
               X_min = 0.0d0
               Y_min = 0.0d0
               Z_min = 0.0d0
               X_max = 0.0d0
               Y_max = 0.0d0
               Z_max = 0.0d0
!
               X1_shift = 0.0d0
               Y1_shift = 0.0d0
               Z1_shift = 0.0d0
!
               X2_shift = 0.0d0
               Y2_shift = 0.0d0
               Z2_shift = 0.0d0
!
               R_min = 0.0d0
               R_max = 0.0d0
!
               R1_shift = 0.0d0
               R2_shift = 0.0d0
!
               dependent_variable_of_surfaces = ' '
               equation_order_of_bounding_surface1 = 0
               equation_order_of_bounding_surface2 = 0
!
               interpolation_data_file_name1 = '        '
               interpolation_data_file_name2 = '        '
               ref_interpolation_data_file   = '        '
!
               sign1 = 1.0d0
               sign2 = 1.0d0
               exponent1 = 0.0d0
               exponent2 = 0.0d0
!
               thickness0 = 0.0d0
               thickness1 = 0.0d0
               thickness2 = 0.0d0
!
               vertical_or_stratigraphic1 = ' '
               vertical_or_stratigraphic2 = ' '
!
               reference_surface1  = ' '
               reference_surface2  = ' '
!
! ............ Initializations - Array components
!
               InclZone(k)%TypeEqu1 = '           '
               InclZone(k)%TypeEqu2 = '           '
!
               InclZone(k)%DependentVar = ' '
!
               InclZone(k)%OrderEquSurf1 = 0
               InclZone(k)%OrderEquSurf2 = 0
!
               InclZone(k)%LMin = 0.0d0     ! ... Whole array operations
               InclZone(k)%LMax = 0.0d0
!
               InclZone(k)%LMin = 0.0d0
               InclZone(k)%LMax = 0.0d0
!
               InclZone(k)%L1Shift = 0.0d0
               InclZone(k)%L2Shift = 0.0d0
!
               InclZone(k)%sign1  = 1.0d0
               InclZone(k)%sign2  = 1.0d0
               InclZone(k)%expon1 = 0.0d0
               InclZone(k)%expon2 = 0.0d0
!
               InclZone(k)%thick0 = 0.0d0
               InclZone(k)%thick1 = 0.0d0
               InclZone(k)%thick2 = 0.0d0
!
               InclZone(k)%VertOrStrat1 = ' '
               InclZone(k)%VertOrStrat2 = ' '
               InclZone(k)%RefSurf1     = ' '
               InclZone(k)%RefSurf2     = ' '
!
! ............ Reading the irregular zone data using a namelist
!
               READ (UNIT = *, NML = Irregular_InclZone, IOSTAT = ier )
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6601) '<Irregular_InclZone>', n
                  STOP
               END IF
!
! ............ Checking the dependent variable
!
               IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &             coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &         THEN

!
                  IF( ( dependent_variable_of_surfaces /= 'X' .AND. dependent_variable_of_surfaces /= 'x'  .AND.  &
     &                  dependent_variable_of_surfaces /= 'Y' .AND. dependent_variable_of_surfaces /= 'y'  .AND.  &
     &                  dependent_variable_of_surfaces /= 'Z' .AND. dependent_variable_of_surfaces /= 'z' ) )     &
     &            THEN
                     WRITE( UNIT = *, FMT = 6608 ) 'cartesian', n, 'X, Y and Z'
                     STOP
                  END IF
!
               ELSE
!
                  IF( ( dependent_variable_of_surfaces /= 'R' .AND. dependent_variable_of_surfaces /= 'r'   .AND. &
     &                  dependent_variable_of_surfaces /= 'Z' .AND. dependent_variable_of_surfaces /= 'z' ) )     &
     &            THEN
                     WRITE( UNIT = *, FMT = 6608 ) 'cylindrical', n, 'R and Z'
                     STOP
                  END IF
!
               END IF
!
! ............ Assignment
!
               InclZone(k)%LMin(1) = X_min
               InclZone(k)%LMin(2) = Y_min
               InclZone(k)%LMin(3) = Z_min
!
               InclZone(k)%LMax(1) = X_max
               InclZone(k)%LMax(2) = Y_max
               InclZone(k)%LMax(3) = Z_max
!
               InclZone(k)%DependentVar = dependent_variable_of_surfaces
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Setting the type of equation of surface 1
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( type_of_equation1(1:4) == 'POLY' .OR. type_of_equation1(1:4) == 'Poly' .OR. type_of_equation1(1:4) == 'poly' ) THEN
!
! ............ Polynomial equation
!
                  type_of_equation1(1:4) = 'Poly'
                  IF( type_of_equation1(5:6) == '/E' .OR. type_of_equation1(5:6) == '/e' ) THEN
                     type_of_equation1(5:6) = '/E'
                  ELSE IF( type_of_equation1(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                     STOP
                  END IF
!
! ............ Power equation
!
               ELSE IF( type_of_equation1(1:4) == 'POWE' .OR. type_of_equation1(1:4) == 'Powe' .OR. type_of_equation1(1:4) == 'powe' ) THEN
!
                  type_of_equation1(1:4) = 'Powe'
                  IF( type_of_equation1(5:6) == '/E' .OR. type_of_equation1(5:6) == '/e' ) THEN
                     type_of_equation1(5:6) = '/E'
                  ELSE IF( type_of_equation1(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                     STOP
                  END IF
!
! ............ No equation: interpolation from a grid/tabular data set
!
               ELSE IF( type_of_equation1(1:4) == 'INTE' .OR. type_of_equation1(1:4) == 'Inte' .OR. type_of_equation1(1:4) == 'inte' ) THEN
!
                  CONTINUE
!
! ............ For Fixed width surfaces
!
               ELSE IF( type_of_equation1(1:4) == 'FIXE' .OR. type_of_equation1(1:4) == 'Fixe' .OR. type_of_equation1(1:4) == 'fixe' ) THEN
!
! ............... Making sure that only one Fixed Width surface is defined
!
                  IF( ( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) .AND. &
     &                ( ref_interpolation_data_file(1:3) == '   ' ) )  &
     &            THEN
                     WRITE( UNIT = *, FMT = 6402 ) n
                     STOP
                  END IF
!
                  CONTINUE
!
! ............ Unknown option: ERROR!
!
               ELSE
                  WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation1(1:6)
                  STOP
               END IF
!
               InclZone(k)%TypeEqu1 = type_of_equation1
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Setting the type of equation of surface 2
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF( type_of_equation2(1:4) == 'POLY' .OR. type_of_equation2(1:4) == 'Poly' .OR. type_of_equation2(1:4) == 'poly' ) THEN
!
! ............ Polynomial equation
!
                  type_of_equation2(1:4) = 'Poly'
                  IF( type_of_equation2(5:6) == '/E' .OR. type_of_equation2(5:6) == '/e' ) THEN
                     type_of_equation2(5:6) = '/E'
                  ELSE IF( type_of_equation2(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                     STOP
                  END IF
!
! ............ Power equation
!
               ELSE IF( type_of_equation2(1:4) == 'POWE' .OR. type_of_equation2(1:4) == 'Powe' .OR. type_of_equation2(1:4) == 'powe' ) THEN
                  type_of_equation2(1:4) = 'Powe'
                  IF( type_of_equation2(5:6) == '/E' .OR. type_of_equation2(5:6) == '/e' ) THEN
                     type_of_equation2(5:6) = '/E'
                  ELSE IF( type_of_equation2(5:6) == '  ' ) THEN
                     CONTINUE
                  ELSE
                     WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                     STOP
                  END IF
!
! ............ No equation: interpolation from a grid/tabular data set
!
               ELSE IF( type_of_equation2(1:4) == 'INTE' .OR. type_of_equation2(1:4) == 'Inte' .OR. type_of_equation2(1:4) == 'inte' ) THEN
!
                  CONTINUE
!
! ............ For Fixed width surfaces
!
               ELSE IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) THEN
!
                  CONTINUE
!
! ............ Unknown option: ERROR!
!
               ELSE
                  WRITE( UNIT = *, FMT = 6401 ) n, type_of_equation2(1:6)
                  STOP
               END IF
!
               InclZone(k)%TypeEqu2 = type_of_equation2
!
               MediaSequenceNumber( TotNum_HetRegions + k ) = n
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For Cartesian coodinates
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               IF_Coord4: IF( coordinates(1:4) == 'CART' .OR. coordinates(1:4) == 'Cart' .OR. coordinates(1:4) == 'cart' .OR. &
     &                        coordinates(1:4) == 'TRAN' .OR. coordinates(1:4) == 'Tran' .OR. coordinates(1:4) == 'tran' )    &
     &         THEN

!
                  IF( ( NXMax >  1 .AND. NYMax == 1 .AND. NZMax == 1 ) .OR. ( NXMax == 1 .AND. NYMax >  1 .AND. NZMax == 1 ) .OR.  &
     &                ( NXMax == 1 .AND. NYMax == 1 .AND. NZMax >  1 ) )    &
     &            THEN
                     WRITE( UNIT = *, FMT = 6606 ) 'cartesian', n
                      STOP
                  END IF
!
                  IF_Inte: IF( type_of_equation1(1:4) == 'INTE' .OR. type_of_equation1(1:4) == 'Inte' .OR. type_of_equation1(1:4) == 'inte' )  THEN
!
                     type_of_equation1(1:4) = 'Inte'
!
                     IntDat_match = .FALSE.
                     IF( Num_IntTables > 0 ) THEN
                         DO_m1: DO m = 1,Num_IntTables
                            IF( interpolation_data_file_name1 == IntTable(m)%FileName ) THEN
                               InclZone(k)%IntTableNum1 = m
                               IntDat_match             = .TRUE.
                               EXIT
                            END IF
                         END DO DO_m1
                     END IF
!
                     IF( .NOT. IntDat_match ) THEN
!
                        m = Num_IntTables + 1
                        Num_IntTables            = m
                        InclZone(k)%IntTableNum1 = Num_IntTables
!
! ..................... Reading info describing the structure of interpolation Table 1
!
                        READ (UNIT = *, NML = Irregular_InclZone_IntTable1, IOSTAT = ier )
!
                        IF(ier /= 0) THEN
                          WRITE(UNIT = *, FMT = 6601) '<Irregular_InclZone_IntTable1>', n
                           STOP
                        END IF
!
                        CALL Read_Tabular_Data( table_number = m, file_name = interpolation_data_file_name1, number_of_rows = number_of_rows, &
     &                                          number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                          read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                          RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                          interpolation_search_radius = interpolation_search_radius,  &
     &                                          num_interrogated_points     = num_interrogated_points )

!
                     END IF
!
                     GO TO 1500
!
                  ELSE IF( type_of_equation1(1:4) == 'FIXE' .OR. type_of_equation1(1:4) == 'Fixe' .OR. type_of_equation1(1:4) == 'fixe' )  THEN
!
                     InclZone(k)%thick0 = thickness0
                     InclZone(k)%thick1 = thickness1
!
                     InclZone(k)%VertOrStrat1 = vertical_or_stratigraphic1
                     InclZone(k)%RefSurf1    = reference_surface1
!
                     IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' ) THEN
!
                        type_of_equation1(1:4) = 'Fixe'
                        type_of_equation2(1:4) = 'Fixe'
!
                        IntDat_match = .FALSE.
!
                        IF( Num_IntTables > 0 ) THEN
                           DO_m2: DO m = 1,Num_IntTables
                              IF( ref_interpolation_data_file == IntTable(m)%FileName ) THEN
                                 InclZone(k)%IntTableNum1 = m
                                 IntDat_match           = .TRUE.
                                 EXIT
                              END IF
                           END DO DO_m2
                        END IF
!
                        IF( .NOT. IntDat_match ) THEN
!
                           m = Num_IntTables + 1
                           Num_IntTables          = m
                           InclZone(k)%IntTableNum1 = Num_IntTables
!
! ........................ Reading info describing the structure of the reference
!
                           READ (UNIT = *, NML = Irregular_InclZone_RefTable, IOSTAT = ier )
!
                           IF(ier /= 0) THEN
                              WRITE(UNIT = *, FMT = 6601) '<Irregular_InclZone_RefTable>', n
                              STOP
                           END IF
!
                           CALL Read_Tabular_Data( table_number = m, file_name = ref_interpolation_data_file, number_of_rows = number_of_rows, &
     &                                             number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                             read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                             RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                             interpolation_search_radius = interpolation_search_radius,  &
     &                                             num_interrogated_points     = num_interrogated_points )

!
                        END IF
!
                     END IF
!
                     GO TO 1500
!
                  END IF IF_Inte

!
! ............... Initializations
!
                  BoundSurface1_EquCoeff_A = 0.0d0
                  BoundSurface1_EquCoeff_B = 0.0d0
                  BoundSurface1_EquCoeff_C = 0.0d0
                  BoundSurface1_EquCoeff_D = 0.0d0
!
! ............... Reading the data for surface 1
!
                  READ (UNIT = *, NML = Irregular_InclZone_Surf1, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_InclZone_Surf1>',n
                     STOP
                  END IF
!
                  InclZone(k)%OrderEquSurf1 = equation_order_of_bounding_surface1
!
                  i1 = equation_order_of_bounding_surface1
!
                  ALLOCATE ( InclZone(k)%Equ1CoeffA(0:i1), InclZone(k)%Equ1CoeffB(0:i1) )
!
                  InclZone(k)%Equ1CoeffA = 0.0d0
                  InclZone(k)%Equ1CoeffB = 0.0d0
!
                  IF_Surf1: IF( InclZone(k)%TypeEqu1(1:4) == 'Powe' ) THEN
                     ALLOCATE ( InclZone(k)%Equ1CoeffC(0:i1), InclZone(k)%Equ1CoeffD(0:i1) )
                     InclZone(k)%Equ1CoeffC = 0.0d0
                     InclZone(k)%Equ1CoeffD = 0.0d0
                  ELSE
                     IF(i1 > 1 ) THEN
                        ALLOCATE ( InclZone(k)%Equ1CoeffC(1:Max(1,i1)) )
                        InclZone(k)%Equ1CoeffC = 0.0d0
                     END IF
                  END IF IF_Surf1
!
! ............... Assignments
!
                  InclZone(k)%sign1  = sign1
                  InclZone(k)%expon1 = exponent1
!
                  InclZone(k)%Equ1CoeffA(0:i1) = BoundSurface1_EquCoeff_A(0:i1)
                  InclZone(k)%Equ1CoeffB(0:i1) = BoundSurface1_EquCoeff_B(0:i1)
!
                  IF( InclZone(k)%TypeEqu1(1:4) == 'Powe' ) THEN
                     InclZone(k)%Equ1CoeffC(0:i1) = BoundSurface1_EquCoeff_C(0:i1)
                     InclZone(k)%Equ1CoeffD(0:i1) = BoundSurface1_EquCoeff_D(0:i1)
                  ELSE
                     IF(i1 > 1) InclZone(k)%Equ1CoeffC(1:i1-1) = BoundSurface1_EquCoeff_C(1:i1-1)
                  END IF
!
                  InclZone(k)%L1Shift(1) = X1_shift
                  InclZone(k)%L1Shift(2) = Y1_shift
                  InclZone(k)%L1Shift(3) = Z1_shift
!
! >>>>>
!
 1500             IF_Inte2: IF( type_of_equation2(1:4) == 'INTE' .OR. type_of_equation2(1:4) == 'Inte' .OR. type_of_equation2(1:4) == 'inte' ) THEN
!
                     type_of_equation2(1:4) = 'Inte'
!
                     IntDat_match = .FALSE.
                     IF( Num_IntTables > 0 ) THEN
                        DO_m3: DO m = 1,Num_IntTables
                           IF( interpolation_data_file_name2 == IntTable(m)%FileName ) THEN
                              InclZone(k)%IntTableNum2 = m
                              IntDat_match             = .TRUE.
                              EXIT
                           END IF
                        END DO DO_m3
                     END IF
!
                     IF( .NOT. IntDat_match ) THEN
!
                        m = Num_IntTables + 1
                        Num_IntTables = m
                        InclZone(k)%IntTableNum2 = Num_IntTables
!
! ..................... Reading info describing the structure of interpolation Table 2
!
                        READ (UNIT = *, NML = Irregular_InclZone_IntTable2, IOSTAT = ier )
!
                        IF(ier /= 0) THEN
                           WRITE(UNIT = *, FMT = 6601) '<Irregular_InclZone_IntTable2>', n
                           STOP
                        END IF
!
                        CALL Read_Tabular_Data( table_number = m, file_name = interpolation_data_file_name2, number_of_rows = number_of_rows, &
     &                                          number_of_columns = number_of_columns, read_data_by_row  = read_data_by_row,              &
     &                                          read_data_format  = read_data_format,  RowCol_DepVariable   = RowCol_DepVariable,         &
     &                                          RowCol_IndVariable_1 = RowCol_IndVariable_1, RowCol_IndVariable_2 = RowCol_IndVariable_2, &
     &                                          interpolation_search_radius = interpolation_search_radius,  &
     &                                          num_interrogated_points     = num_interrogated_points )

!
                     END IF
!
                     CYCLE DO_Inclusions
!
                  ELSE IF( type_of_equation2(1:4) == 'FIXE' .OR. type_of_equation2(1:4) == 'Fixe' .OR. type_of_equation2(1:4) == 'fixe' )  THEN
!
                     InclZone(k)%thick2 = thickness2
!
                     InclZone(k)%VertOrStrat2 = vertical_or_stratigraphic2
                     InclZone(k)%RefSurf2    = reference_surface2
!
                     CYCLE DO_Inclusions
!
                  END IF IF_Inte2
!
! ............... Initializations
!
                  BoundSurface2_EquCoeff_A = 0.0d0
                  BoundSurface2_EquCoeff_B = 0.0d0
                  BoundSurface2_EquCoeff_C = 0.0d0
                  BoundSurface2_EquCoeff_D = 0.0d0
!
! ............... Reading the data for surface 2
!
                  READ (UNIT = *, NML = Irregular_InclZone_Surf2, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_InclZone_Surf2>', n
                     STOP
                  END IF
!
! .................. Checking consistency of orders of equation
!
                  IF( ( equation_order_of_bounding_surface1 < 0 ) .OR. ( equation_order_of_bounding_surface2 < 0 ) ) THEN
                     WRITE( UNIT = *, FMT = 6604 ) n
                     STOP
                  END IF
!
                  InclZone(k)%OrderEquSurf2 = equation_order_of_bounding_surface2
!
                  i2 = equation_order_of_bounding_surface2
!
                  ALLOCATE ( InclZone(k)%Equ2CoeffA(0:i2), InclZone(k)%Equ2CoeffB(0:i2) )
!
                  InclZone(k)%Equ2CoeffA = 0.0d0
                  InclZone(k)%Equ2CoeffB = 0.0d0
!
                  IF_Surf2: IF( InclZone(k)%TypeEqu2(1:4) == 'Powe' ) THEN
                     ALLOCATE ( InclZone(k)%Equ2CoeffC(0:i2), InclZone(k)%Equ2CoeffD(0:i2) )
                     InclZone(k)%Equ2CoeffC = 0.0d0
                     InclZone(k)%Equ2CoeffD = 0.0d0
                  ELSE
                     IF(i2 > 1 ) THEN
                        ALLOCATE ( InclZone(k)%Equ2CoeffC(1:Max(1,i2)) )
                        InclZone(k)%Equ2CoeffC = 0.0d0
                     END IF
                  END IF IF_Surf2
!
! ............... Assignments
!
                  InclZone(k)%sign2  = sign2
                  InclZone(k)%expon2 = exponent2
!
                  InclZone(k)%Equ2CoeffA(0:i2) = BoundSurface2_EquCoeff_A(0:i2)
                  InclZone(k)%Equ2CoeffB(0:i2) = BoundSurface2_EquCoeff_B(0:i2)
!
                  IF( InclZone(k)%TypeEqu2(1:4) == 'Powe' ) THEN
                     InclZone(k)%Equ2CoeffC(0:i2) = BoundSurface2_EquCoeff_C(0:i2)
                     InclZone(k)%Equ2CoeffD(0:i2) = BoundSurface2_EquCoeff_D(0:i2)
                  ELSE
                     IF(i2 > 1) InclZone(k)%Equ2CoeffC(1:i2-1) = BoundSurface2_EquCoeff_C(1:i2-1)
                  END IF
!
                  InclZone(k)%L2Shift(1) = X2_shift
                  InclZone(k)%L2Shift(2) = Y2_shift
                  InclZone(k)%L2Shift(3) = Z2_shift
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> For Cylindrical coodinates
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               ELSE
!
! ............... Irregular shapes
!
                  IF( NXMax >  1 .AND. NZMax == 1 ) THEN
                     WRITE( UNIT = *, FMT = 6606 ) 'radial', n
                     STOP
                  END IF
!
! ............... Initializations
!
                  BoundSurface1_EquCoeff_A = 0.0d0
                  BoundSurface1_EquCoeff_B = 0.0d0
                  BoundSurface2_EquCoeff_A = 0.0d0
                  BoundSurface2_EquCoeff_B = 0.0d0
!
! ............... Reading the data for surface 1
!
                  READ (UNIT = *, NML = Irregular_InclZone_Surf1, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_InclZone_Surf1>', n
                     STOP
                  END IF
!
! ............... Reading the data for surface 2
!
                  READ (UNIT = *, NML = Irregular_InclZone_Surf2, IOSTAT = ier )
!
                  IF(ier /= 0) THEN
                     WRITE(UNIT = *, FMT = 6601) '<Irregular_InclZone_Surf2>', n
                     STOP
                  END IF
!
! ............... Checking consistency of orders of equation
!
                  IF( ( equation_order_of_bounding_surface1 < 0 ) .OR. ( equation_order_of_bounding_surface2 < 0 ) ) THEN
                     WRITE( UNIT = *, FMT = 6604 ) n
                     STOP
                  END IF
!
                  InclZone(k)%OrderEquSurf1 = equation_order_of_bounding_surface1
                  InclZone(k)%OrderEquSurf2 = equation_order_of_bounding_surface2
!
                  i1 = equation_order_of_bounding_surface1
                  ALLOCATE ( InclZone(k)%Equ1CoeffA(0:i1) )
!
                  InclZone(k)%Equ1CoeffA = 0.0d0
!
                  IF( InclZone(k)%TypeEqu1(1:4) == 'Powe' ) THEN
                     ALLOCATE ( InclZone(k)%Equ1CoeffB(0:i1) )
                     InclZone(k)%Equ1CoeffB = 0.0d0
                  END IF
!
                  i2 = equation_order_of_bounding_surface2
                  ALLOCATE ( InclZone(k)%Equ2CoeffA(0:i2) )
!
                  InclZone(k)%Equ2CoeffA = 0.0d0
!
                  IF( InclZone(k)%TypeEqu2(1:4) == 'Powe' ) THEN
                     ALLOCATE ( InclZone(k)%Equ2CoeffB(0:i2) )
                     InclZone(k)%Equ2CoeffB = 0.0d0
                  END IF
!
! ............... Assignments
!
                  InclZone(k)%LMin(1) = R_min
                  InclZone(k)%LMin(3) = Z_min
!
                  InclZone(k)%LMax(1) = R_max
                  InclZone(k)%LMax(3) = Z_max
!
                  InclZone(k)%Equ1CoeffA(0:i1) = BoundSurface1_EquCoeff_A(0:i1)
                  InclZone(k)%Equ2CoeffA(0:i2) = BoundSurface2_EquCoeff_A(0:i2)
!
                  IF( InclZone(k)%TypeEqu1(1:4) == 'Powe' ) InclZone(k)%Equ1CoeffB(0:i1) = BoundSurface1_EquCoeff_B(0:i1)
                  IF( InclZone(k)%TypeEqu2(1:4) == 'Powe' ) InclZone(k)%Equ2CoeffB(0:i2) = BoundSurface2_EquCoeff_B(0:i2)
!
                  InclZone(k)%L1Shift(1) = R1_shift
                  InclZone(k)%L1Shift(3) = Z1_shift
!
                  InclZone(k)%L2Shift(1) = R2_shift
                  InclZone(k)%L2Shift(3) = Z2_shift
!
               END IF IF_Coord4
!
!
!
            END DO DO_Inclusions
!
! <<<<<<<<<
! <<<<<<<<<
! <<<<<<<<<
!
            IF( k /= TotNum_InclZones ) THEN
               WRITE(UNIT = *, FMT = 6705) k, TotNum_InclZones, number_of_inclusion_zones, number_of_periodic_InclZones, total_number_periodic_InclZones
               STOP
            END IF
!
!***********************************************************************
!*                                                                     *
!*        Convert the inclusion zone ranges into SI units (m)          *
!*                                                                     *
!***********************************************************************
!
! ......... Conversion into METERS from INCHES
!
            FORALL (n=2:TotNum_InclZones, InclZone(k)%units == 'IN' .OR. InclZone(k)%units == 'in' .OR. InclZone(k)%units == 'In')
!
               InclZone(n)%LMin(1:3) = InclZone(n)%LMin(1:3) * 2.54d-2
               InclZone(n)%LMax(1:3) = InclZone(n)%LMax(1:3) * 2.54d-2
!
               InclZone(n)%CylBase1Coord(1:3) = InclZone(n)%CylBase1Coord(1:3) * 2.54d-2
               InclZone(n)%CylBase2Coord(1:3) = InclZone(n)%CylBase2Coord(1:3) * 2.54d-2
!
               InclZone(n)%CylRmin = InclZone(n)%CylRmin * 2.54d-2
               InclZone(n)%CylRmax = InclZone(n)%CylRmax * 2.54d-2
!
               InclZone(n)%SphereCenterCoord(1:3) = InclZone(n)%SphereCenterCoord(1:3) * 2.54d-2
!
               InclZone(n)%SphRmin = InclZone(n)%SphRmin * 2.54d-2
               InclZone(n)%SphRmax = InclZone(n)%SphRmax * 2.54d-2
!
               InclZone(n)%L1Shift(1:3) = InclZone(n)%L1Shift(1:3) * 2.54d-2
               InclZone(n)%L2Shift(1:3) = InclZone(n)%L2Shift(1:3) * 2.54d-2
!
               InclZone(n)%EllBase1Coord(1:3) = InclZone(n)%EllBase1Coord(1:3) * 2.54d-2
               InclZone(n)%EllBase1Coord(1:3) = InclZone(n)%EllBase1Coord(1:3) * 2.54d-2
!
               InclZone(n)%LAxis = InclZone(n)%LAxis * 2.54d-2
               InclZone(n)%SAxis = InclZone(n)%SAxis * 2.54d-2
!
               InclZone(n)%Base1LAxis = InclZone(n)%Base1LAxis * 2.54d-2
               InclZone(n)%Base1SAxis = InclZone(n)%Base1SAxis * 2.54d-2
               InclZone(n)%Base2LAxis = InclZone(n)%Base2LAxis * 2.54d-2
               InclZone(n)%Base2SAxis = InclZone(n)%Base2SAxis * 2.54d-2
!
               InclZone(n)%thick0 = InclZone(n)%thick0 * 2.54d-2
               InclZone(n)%thick1 = InclZone(n)%thick1 * 2.54d-2
               InclZone(n)%thick2 = InclZone(n)%thick2 * 2.54d-2
!
            END FORALL
!
! ......... Conversion into METERS from FEET
!
            FORALL (n=2:TotNum_InclZones, InclZone(n)%units == 'FT' .OR. InclZone(n)%units == 'ft' .OR. InclZone(n)%units == 'Ft')
!
               InclZone(n)%LMin(1:3) = InclZone(n)%LMin(1:3) * 3.038d-1
               InclZone(n)%LMax(1:3) = InclZone(n)%LMax(1:3) * 3.038d-1
!
               InclZone(n)%CylBase1Coord(1:3) = InclZone(n)%CylBase1Coord(1:3) * 3.038d-1
               InclZone(n)%CylBase2Coord(1:3) = InclZone(n)%CylBase2Coord(1:3) * 3.038d-1
!
               InclZone(n)%CylRmin = InclZone(n)%CylRmin * 3.038d-1
               InclZone(n)%CylRmax = InclZone(n)%CylRmax * 3.038d-1
!
               InclZone(n)%SphereCenterCoord(1:3) = InclZone(n)%SphereCenterCoord(1:3) * 3.038d-1
!
               InclZone(n)%SphRmin = InclZone(n)%SphRmin * 3.038d-1
               InclZone(n)%SphRmax = InclZone(n)%SphRmax * 3.038d-1
!
               InclZone(n)%L1Shift(1:3) = InclZone(n)%L1Shift(1:3) * 3.038d-1
               InclZone(n)%L2Shift(1:3) = InclZone(n)%L2Shift(1:3) * 3.038d-1
!
               InclZone(n)%EllBase1Coord(1:3) = InclZone(n)%EllBase1Coord(1:3) * 3.038d-1
               InclZone(n)%EllBase1Coord(1:3) = InclZone(n)%EllBase1Coord(1:3) * 3.038d-1
!
               InclZone(n)%LAxis = InclZone(n)%LAxis * 3.038d-1
               InclZone(n)%SAxis = InclZone(n)%SAxis * 3.038d-1
!
               InclZone(n)%Base1LAxis = InclZone(n)%Base1LAxis * 3.038d-1
               InclZone(n)%Base1SAxis = InclZone(n)%Base1SAxis * 3.038d-1
               InclZone(n)%Base2LAxis = InclZone(n)%Base2LAxis * 3.038d-1
               InclZone(n)%Base2SAxis = InclZone(n)%Base2SAxis * 3.038d-1
!
               InclZone(n)%thick0 = InclZone(n)%thick0 * 3.038d-1
               InclZone(n)%thick1 = InclZone(n)%thick1 * 3.038d-1
               InclZone(n)%thick2 = InclZone(n)%thick2 * 3.038d-1
!
            END FORALL
!
! ......... Conversion into METERS from KM
!
            FORALL (n=2:TotNum_InclZones, InclZone(n)%units == 'KM' .OR. InclZone(n)%units == 'km' .OR. InclZone(n)%units == 'Km')
!
               InclZone(n)%LMin(1:3) = InclZone(n)%LMin(1:3) * 1.0d3
               InclZone(n)%LMax(1:3) = InclZone(n)%LMax(1:3) * 1.0d3
!
               InclZone(n)%CylBase1Coord(1:3) = InclZone(n)%CylBase1Coord(1:3) * 1.0d3
               InclZone(n)%CylBase2Coord(1:3) = InclZone(n)%CylBase2Coord(1:3) * 1.0d3
!
               InclZone(n)%CylRmin = InclZone(n)%CylRmin * 1.0d3
               InclZone(n)%CylRmax = InclZone(n)%CylRmax * 1.0d3
!
               InclZone(n)%SphereCenterCoord(1:3) = InclZone(n)%SphereCenterCoord(1:3) * 1.0d3
!
               InclZone(n)%SphRmin = InclZone(n)%SphRmin * 1.0d3
               InclZone(n)%SphRmax = InclZone(n)%SphRmax * 1.0d3
!
               InclZone(n)%L1Shift(1:3) = InclZone(n)%L1Shift(1:3) * 1.0d3
               InclZone(n)%L2Shift(1:3) = InclZone(n)%L2Shift(1:3) * 1.0d3
!
               InclZone(n)%EllBase1Coord(1:3) = InclZone(n)%EllBase1Coord(1:3) * 1.0d3
               InclZone(n)%EllBase1Coord(1:3) = InclZone(n)%EllBase1Coord(1:3) * 1.0d3
!
               InclZone(n)%LAxis = InclZone(n)%LAxis * 1.0d3
               InclZone(n)%SAxis = InclZone(n)%SAxis * 1.0d3
!
               InclZone(n)%Base1LAxis = InclZone(n)%Base1LAxis * 1.0d3
               InclZone(n)%Base1SAxis = InclZone(n)%Base1SAxis * 1.0d3
               InclZone(n)%Base2LAxis = InclZone(n)%Base2LAxis * 1.0d3
               InclZone(n)%Base2SAxis = InclZone(n)%Base2SAxis * 1.0d3
!
               InclZone(n)%thick0 = InclZone(n)%thick0 * 1.0d3
               InclZone(n)%thick1 = InclZone(n)%thick1 * 1.0d3
               InclZone(n)%thick2 = InclZone(n)%thick2 * 1.0d3
!
            END FORALL
!
! ......... Conversion into METERS from MM
!
            FORALL (n=2:TotNum_InclZones, InclZone(n)%units == 'MM' .OR. InclZone(n)%units == 'mm' .OR. InclZone(n)%units == 'Mm')
!
               InclZone(n)%LMin(1:3) = InclZone(n)%LMin(1:3) * 1.0d-3
               InclZone(n)%LMax(1:3) = InclZone(n)%LMax(1:3) * 1.0d-3
!
               InclZone(n)%CylBase1Coord(1:3) = InclZone(n)%CylBase1Coord(1:3) * 1.0d-3
               InclZone(n)%CylBase2Coord(1:3) = InclZone(n)%CylBase2Coord(1:3) * 1.0d-3
!
               InclZone(n)%CylRmin = InclZone(n)%CylRmin * 1.0d-3
               InclZone(n)%CylRmax = InclZone(n)%CylRmax * 1.0d-3
!
               InclZone(n)%SphereCenterCoord(1:3) = InclZone(n)%SphereCenterCoord(1:3) * 1.0d-3
!
               InclZone(n)%SphRmin = InclZone(n)%SphRmin * 1.0d-3
               InclZone(n)%SphRmax = InclZone(n)%SphRmax * 1.0d-3
!
               InclZone(n)%L1Shift(1:3) = InclZone(n)%L1Shift(1:3) * 1.0d-3
               InclZone(n)%L2Shift(1:3) = InclZone(n)%L2Shift(1:3) * 1.0d-3
!
               InclZone(n)%EllBase1Coord(1:3) = InclZone(n)%EllBase1Coord(1:3) * 1.0d-3
               InclZone(n)%EllBase1Coord(1:3) = InclZone(n)%EllBase1Coord(1:3) * 1.0d-3
!
               InclZone(n)%LAxis = InclZone(n)%LAxis * 1.0d-3
               InclZone(n)%SAxis = InclZone(n)%SAxis * 1.0d-3
!
               InclZone(n)%Base1LAxis = InclZone(n)%Base1LAxis * 1.0d-3
               InclZone(n)%Base1SAxis = InclZone(n)%Base1SAxis * 1.0d-3
               InclZone(n)%Base2LAxis = InclZone(n)%Base2LAxis * 1.0d-3
               InclZone(n)%Base2SAxis = InclZone(n)%Base2SAxis * 1.0d-3
!
               InclZone(n)%thick0 = InclZone(n)%thick0 * 1.0d-3
               InclZone(n)%thick1 = InclZone(n)%thick1 * 1.0d-3
               InclZone(n)%thick2 = InclZone(n)%thick2 * 1.0d-3
!
            END FORALL
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Define_Inclusion_Zones 1.0 .............. 9 April     2015',6X,'Defining inclusion zones in the grid')
!
 6010 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The shape of the inclusion zone # ',I3.3,' is "',A,'": Unknown/Unavailable option'//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6050 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input X_min =',ES12.5,' of the inclusion zone # ',I3.3,' is larger than (or equal to) the X_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6051 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input Y_min =',ES12.5,' of the inclusion zone # ',I3.3,' is larger than (or equal to) the Y_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6052 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input Z_min =',ES12.5,' of the inclusion zone # ',I3.3,' is larger than (or equal to) the Z_max =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6053 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'In the cylindrical inclusion zone # ',I3.3,', the centers of the two bases coincide ',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6054 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <cylinder_Rmin> =',ES12.5,' of the cylindrical inclusion zone # ',I3.3,  &
     &          ' is larger than (or equal to) the <cylinder_Rmax> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6055 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <sphere_Rmin> =',ES12.5,' of the spherical inclusion zone # ',I3.3,  &
     &          ' is larger than (or equal to) the <sphere_Rmax> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6056 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <Base1_Short_Axis> =',ES12.5,' of the elliptical inclusion zone # ',I3.3,  &
     &          ' is larger than (or equal to) the <Base1_Long_Axis> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6058 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <Base2_Short_Axis> =',ES12.5,' of the elliptical inclusion zone # ',I3.3,  &
     &          ' is larger than (or equal to) the <Base2_Long_Axis> =',ES12.5,//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6059 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The input <plane_of_ellipse_bases> =',A,' of the elliptical inclusion zone # ',I3.3,' is NOT an available option',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6101 FORMAT(T5,'Memory allocation to array <',A,'> in subroutine <Define_Inclusion_Zones> was successful')
 6102 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to array <',A,'> in subroutine <Define_Inclusion_Zones> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6105 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The subdomain shape (read by <inclusion_zone_shape> = ',a11,' in subroutine "Define_Inclusion_Zones") is unavailable',//, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6401 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: The inclusion zone #',i3.3,' has an irregular shape that is described by ',/, &
     &       T10,'                                         an unknown/unavailable type of equation (= "',A,'")',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6402 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: The inclusion zone #',i3.3,' is defined by two irregular surfaces of type "FIXED" but the ' ,/,  &
     &       T10,'                                         name of the needed reference interpolation file is not defined',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED     !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6515 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: The dataset/namelist "',A,'" must be ended by the "<<<" descriptor'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6520 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: It is not possible to have a rectangular inclusion zone in a cylindrical coordinate system ',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6600 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: There is a problem reading the namelist <',A,'>',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6601 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: There is a problem reading the namelist <',A,'> in inclusion zone #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6602 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: There is a conflict between the various radii options that describe the ', /,&
     &       T10,'                                        ',A,' shape of inclusion zone #',i3.3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6604 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure "Define_Inclusion_Zones": At least one the orders of the equations describing the bounding surfaces 1 and/or 2 ', /,&
     &       T10,'                                         of inclusion zone #',i3.3,' is < 0',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6605 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: The inclusion zone #',i3.3,' has an irregular shape but the types of the ',/, &
     &       T10,'     equations describing the bounding surfaces 1 and/or 2 of the inclusion zone are non-blanks',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6606 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: This ',A,' system has 1 active dimension - it is not possible for the ', /, &
     &       T10,'                                         inclusion zone #',i3.3,' to be irregularly shaped',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6608 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: In this ',A,' system, one or more of the dependent variables defining the bounding surfaces 1 and 2', /,&
     &       T10,'                                         of inclusion zone #',i3.3,' is not among the possible options (',A,')',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6701 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Define_Inclusion_Zones>: The shape of inclusion zone #',i3.3,' is "PERIODIC", but the coordinate system = "',A,'"' /,&
     &       T10,'                                         conflicts with the <axis_of_perodicity> = "',A1,'"',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6702 FORMAT(//,20('ERROR-'),//,   &
     &       T10,'>>>  Procedure "Define_Inclusion_Zones": The shape of inclusion zone #',i3.3,' is "PERIODIC", but the axis of periodicity = "',A1,'"' /,&
     &       T10,'                                         is not among the available options ("X,", "R", "Y", or "Z")',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6705 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'>>>  Procedure <Define_Inclusion_Zones>: The total number of computed inclusion zones = ',i3.3,' does not match the number <TotNum_InclZones> that is estimated from:',/, &
     &       T5,'                                         TotNum_InclZones(=',i3.3,') = number_of_inclusion_zones(=',i3.3,') - number_of_periodic_InclZones(=',i3.3,') + total_number_periodic_InclZones(=',i3.3,')',/, &
     &       T5,'                                         Check all the related inputs', //, &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of <Define_Inclusion_Zones>
!
!
            RETURN
!
         END SUBROUTINE Define_Inclusion_Zones
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE RZ2D
!
         USE MeshMaker_Data, DR => DX
         USE Grid_Generation_Parameters
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND=8) :: Delta_R, Delta_Z, R_max, Z_max, DR_log, DZ_log, Z_last
!
      REAL(KIND = 16) :: IncrFact
!
! -------
! ... Character variables
! -------
!
      CHARACTER( LEN =  6 ) :: first_part
      CHARACTER( LEN = 15 ) :: option
      CHARACTER( LEN = 22 ) :: discretization_descriptor
!
      CHARACTER( LEN = 120) :: read_format
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: NR = 0, Nrb, NRe, NZ = 0, NZo, NZ1, ICALL = 0, error_code, i, n, n1
!
      INTEGER :: num_DR_subsets, num_DZ_subsets
      INTEGER :: number_of_DRs,  number_of_DZs, number_of_radii
!
! -------
! ... Saving variables
! -------
!
      SAVE ICALL, NR, NZ
!
! -------
! ... Namelists
! -------
!
      NAMELIST/ DR_Subsets / num_DR_subsets
      NAMELIST/ DZ_Subsets / num_DZ_subsets
!
      NAMELIST/ DR_Data / number_of_DRs, Delta_R, read_format, R_max, DR_Log, option, number_of_radii
      NAMELIST/ DZ_Data / number_of_DZs, Delta_Z, read_format, Z_max, DZ_Log, option
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of <RZ2D>
!
!
      ICALL=ICALL+1
      IF(ICALL == 1) WRITE(*,6000)
!6000 FORMAT('RZ2D              1.0   11 FEBRUARY  1992',6X,'CALCULATE 2-D R-Z MESH FROM INPUT DATA')
!6000 FORMAT('RZ2D              1.5   24 September 2013',6X,'CALCULATE 2-D R-Z MESH FROM INPUT DATA')
 6000 FORMAT(/,'RZ2D 2.0 ............................... 15 January   2015',6X,'CALCULATE 2-D R-Z MESH FROM INPUT DATA')
!
      ALLOCATE ( V(NXmax), A(NXmax), D1(NXmax), D2(NXmax), RC(NXmax+1), DR(NXmax), DZ(NZMax), excluded( NXmax * NZmax ), STAT=error_code )
!
      IF(error_code /= 0) THEN
         WRITE(UNIT = *, FMT = 6850) 'V, A, D1, D2, RC, DR, DZ, excluded', 'RZ2D'
         STOP
      END IF
!
      excluded = .FALSE.
!
! ... Initialization
!
      RC = 0.0d0 ! ... Whole array operation
!
      NR = 0
!
!
!
      DO_NumActDir: DO n = 1, 2
!
! >>>>>>
! >>>>>>
! >>>>>>
!
         READ( UNIT = *, FMT = '(a22)', ADVANCE = 'NO', IOSTAT = error_code ) discretization_descriptor
!
         IF(error_code /= 0) THEN
            WRITE(UNIT = *, FMT = 6500)
            STOP
         END IF
!
         SELECT CASE (discretization_descriptor)
!
         CASE(':::>>>R-DISCRETIZATION', ':::>>>R-Discretization', ':::>>>R-discretization', ':::>>>r-discretization')
!
            IF_R: IF( n /= 1 ) THEN
               WRITE(UNIT = *, FMT = 6501) discretization_descriptor, ':::>>>R-DISCRETIZATION'
               STOP
            ELSE
!
               num_DR_subsets = 0
               READ (UNIT = *, NML = DR_Subsets, IOSTAT = error_code )
!
! ............ Stop if there is a problem reading the namelist
!
               IF(error_code /= 0) THEN
                  WRITE(UNIT = *, FMT = 6601) '<DR_Subsets>'
                  STOP
               END IF
!
               DO_RDiscr: DO n1=1,num_DR_subsets + 1
!
                  IF( n1 == num_DR_subsets + 1) THEN
                     READ( UNIT = *, FMT = '(a6)' ) first_part
                     IF( first_part /= ':::<<<' ) THEN
                        WRITE(UNIT = *, FMT = 6515) discretization_descriptor
                        STOP
                     ELSE
                        CYCLE DO_NumActDir
                     END IF
                  END IF
!
! ............... Initializations
!
                  number_of_radii = 0
                  number_of_DRs   = 0
                  DR_Log          = 0.0d0
                  Delta_R         = 0.0d0
                  R_max           = 0.0d0
                  option      = '  '
                  read_format = '  '
!
! ............... Reading the namelists
!
                  READ (UNIT = *, NML = DR_Data, IOSTAT = error_code )
!
                  IF(error_code /= 0) THEN
                     WRITE(UNIT = *, FMT = 6602) '<D'//discretization_descriptor(7:7)//'_Data>',n1
                     STOP
                  END IF
!
                  SELECT_R: SELECT CASE(option(1:1))
!
! ............... Various discretization options
!
                  CASE('L', 'l')
!
                     IF( NR+number_of_DRs > NXMax ) THEN
                        WRITE(UNIT = *, FMT = 6520) 'NR', NR+number_of_DRs, 'MaxNum_R_Subdivisions', NXMax
                        STOP
                     END IF
!
                     R_max = R_max * factor
!
                     IF(DR_Log /= 0.0d0) THEN
                        DR_Log = DR_Log * factor
                     ELSE
                        IF( NR >= 2 ) THEN
                           DR_Log = (RC(NR) - RC(NR-1)) * factor
                        ELSE
                           WRITE(UNIT = *, FMT = 6508)
                           STOP
                        END IF
                     END IF
!
                     IF(NR == 1) RC(1) = 1.0d-5
!
                     NRb = NR+1
                     NRe = NR+number_of_DRs
!
                     WRITE( UNIT = *, FMT = 109)
!
                     CALL Logarithmic_Discretization( axis   = 'R',     D_max    = R_max,     &
    &                                                 D_last = RC(NR),  D_Log    = DR_Log,    &
    &                                                 Nb     = NRb,     IncrFact = IncrFact,  &
    &                                                 DL_number = number_of_DRs,  Dd = DR(NRb:NRe) )
!
                     DO i=1,number_of_DRs
                        RC(NR+i) = RC(NR+i-1) + DR(NR+i)
                     END DO
!
                     NR = NRe
!
                     IF( ABS((RC(NR) - R_max)/R_max) > 1.0d-2 ) THEN
                        WRITE(UNIT = *, FMT = 6045) RC(NR), R_max, IncrFact
                        STOP
                     END IF
!
                  CASE('R', 'r')
!
                     IF( NR+number_of_radii-1 > NXMax ) THEN
                        WRITE(UNIT = *, FMT = 6520) 'NR', NR+number_of_radii-1, 'MaxNum_R_Subdivisions', NXMax
                        STOP
                     END IF
!
                     IF( number_of_radii == 0 ) THEN
                        WRITE(UNIT = *, FMT = 6510)
                        STOP
                     END IF
!
                     IF( TRIM(ADJUSTL(read_format)) == '(*)' .OR. TRIM(ADJUSTL(read_format)) == '*' ) THEN
                        READ( UNIT = *, FMT = * ) (RC(i), i = NR+1, NR+number_of_radii)
                     ELSE IF( read_format(1:2) == '  ' ) THEN
                        WRITE( UNIT = *, FMT = 6521 ) 'R', 'Radii', read_format, n1
                        STOP
                     ELSE
                        READ( UNIT = *, FMT = read_format ) ( RC(i), i = NR+1, NR+number_of_radii )
                     END IF
!
                     NR1 = NR + 1
                     NR  = NR + number_of_radii
!
                     RC(NR1:NR) = RC(NR1:NR) * factor  ! ... Whole array operation
!
                  CASE('U', 'u', 'E', 'e', 'V', 'v')
!
                     IF( NR+number_of_DRs > NXMax ) THEN
                        WRITE(UNIT = *, FMT = 6520) 'NR', number_of_DRs, 'MaxNum_R_Subdivisions', NXMax
                        STOP
                     END IF
!
                     Case_DR: SELECT CASE(option(1:1))
                     CASE('U', 'u', 'E', 'e')
                        DR(1:number_of_DRs) = Delta_R
                     CASE('V', 'v')
                        IF( TRIM(ADJUSTL(read_format)) == '(*)' .OR. TRIM(ADJUSTL(read_format)) == '*' ) THEN
                           READ( UNIT = *, FMT = * ) (DR(i), i = 1,number_of_DRs)
                        ELSE IF( read_format(1:2) == '  ' ) THEN
                           WRITE( UNIT = *, FMT = 6521 ) 'R', 'Variable', read_format, n1
                           STOP
                        ELSE
                           READ( UNIT = *, FMT = read_format ) (DR(i), i = 1,number_of_DRs)
                        END IF
                     END SELECT Case_DR
!
                     IF(NR == 0) THEN
                        NR    = 1
                        RC(1) = 1.0d-5
                     END IF
!
                     DO i = 1,number_of_DRs
                        RC(NR+i) = RC(NR+i-1) + DR(i)
                     END DO
!
                     NR1 = NR + 1
                     NR  = NR + number_of_DRs
!
                     RC(NR1:NR) = RC(NR1:NR) * factor  ! ... Whole array operation
!
                  CASE DEFAULT
!
                     WRITE(UNIT = *, FMT = 6505) option
                     STOP
!
                  END SELECT SELECT_R
!
               END DO DO_RDiscr
!
            END IF IF_R
!
         CASE(':::>>>Z-DISCRETIZATION',':::>>>Z-Discretization', ':::>>>Z-discretization', ':::>>>z-discretization')
!
            IF_Z: IF( n /= 2 ) THEN
               WRITE(UNIT = *, FMT = 6501) discretization_descriptor, ':::>>>Z-DISCRETIZATION'
               STOP
            ELSE
!
               num_DZ_subsets = 0
               READ (UNIT = *, NML = DZ_Subsets, IOSTAT = error_code )
!
! ............ Stop if there is a problem reading the namelist
!
               IF(error_code /= 0) THEN
                  WRITE(UNIT = *, FMT = 6601) '<DZ_Subsets>'
                  STOP
               END IF
!
               DO_ZDiscr: DO n1=1,num_DZ_subsets + 1
!
                  IF( n1 == num_DZ_subsets + 1 ) THEN
                     READ( UNIT = *, FMT = '(a6)' ) first_part
                     IF( (n1 == num_DZ_subsets + 1) .AND. (first_part /= ':::<<<') ) THEN
                        WRITE(UNIT = *, FMT = 6515) discretization_descriptor
                        STOP
                     ELSE
                        EXIT DO_NumActDir
                     END IF
                  END IF
!
! ............... Initializations
!
                  number_of_DZs   = 0
                  DZ_Log          = 0.0d0
                  Delta_Z         = 0.0d0
                  Z_max           = 0.0d0
                  option      = '  '
                  read_format = '  '
!
! ............... Reading the namelists
!
                  READ (UNIT = *, NML = DZ_Data, IOSTAT = error_code )
!
                  IF(error_code /= 0) THEN
                     WRITE(UNIT = *, FMT = 6602) '<D'//discretization_descriptor(7:7)//'_Data>',n1
                     STOP
                  END IF
!
                  NZo = NZ
                  NZ1 = NZ + 1
                  NZ  = NZ + number_of_DZs
!
                  IF( NZ > NZMax ) THEN
                     WRITE(UNIT = *, FMT = 6520) 'NZ', NZ, 'MaxNum_Z_Subdivisions', NZMax
                     STOP
                  END IF
!
                  Case_DZ: SELECT CASE(option(1:1))
                  CASE('L', 'l')
 !
                    IF_DZ0: IF( ABS(DZ_Log) <= 1.0d-7 ) THEN
!
                        IF( NZo > 1 ) THEN
                           DZ_Log = DZ(NZo)
                           Z_last = SUM(DZ(1:NZo))
                        ELSE IF( NZo == 1 ) THEN
                           DZ_Log = DZ(1)
                           Z_last = DZ(1)
                        ELSE IF( NZo == 0 ) THEN
                           WRITE( UNIT = *, FMT = 6050) 'Z'
                           STOP
                        END IF
!
                     ELSE
!
                        IF( NZo > 1 ) THEN
                           Z_last = SUM(DZ(1:NZo))
                        ELSE IF( NZo == 1 ) THEN
                           Z_last = DZ(1)
                        ELSE IF( NZo == 0 ) THEN
                           Z_last = 0.0d0
                        END IF
!
                     END IF IF_DZ0
!
                     CALL Logarithmic_Discretization( axis   = 'Z',     D_max    = Z_max,     &
    &                                                 D_last = Z_last,  D_Log    = DZ_Log,   &
    &                                                 Nb     = NZo,     IncrFact = IncrFact,  &
    &                                                 DL_number = number_of_DZs,  Dd = DZ(NZ1:NZ) )
!
                  CASE('U', 'u', 'E', 'e')
                     DZ(NZ1:NZ) = Delta_Z
                  CASE DEFAULT
                     IF( TRIM(ADJUSTL(read_format)) == '(*)' .OR. TRIM(ADJUSTL(read_format)) == '*' ) THEN
                        READ( UNIT = *, FMT = * ) (DZ(i), i = NZ1,NZ)
                     ELSE IF( read_format(1:2) == '  ' ) THEN
                        WRITE( UNIT = *, FMT = 6521 ) 'Z', 'Variable', read_format, n1
                        STOP
                     ELSE
                        READ( UNIT = *, FMT = read_format ) (DZ(i), i = NZ1,NZ)
                     END IF
                  END SELECT Case_DZ
!
               END DO DO_ZDiscr
!
               DZ = DZ * factor  ! ... Whole array operation
!
            END IF IF_Z
!
         CASE DEFAULT
!
            WRITE(UNIT = *, FMT = 6505) discretization_descriptor
            STOP
!
         END SELECT
!
! >>>>>>
! >>>>>>
! >>>>>>
!
      END DO DO_NumActDir
!
! ... Converting all lengths to meters
!
  109 FORMAT(/' ***** FIND APPROPRIATE FACTOR FOR INCREASING RADIAL DISTANCES *****')
!
   42 FORMAT(' AT ITERATION',I3,', XL=',ES20.14E1,', XR =',ES20.14E1, ', XM =',ES12.6E1,', SMID=',ES12.6E1,', AM=',ES12.6E1)
   44 FORMAT(' CONVERGENCE FAILURE IN SUBROUTINE RZ2D',/,' STOP EXECUTION')
   30 FORMAT(/,T2,' After ',I4.4,' iterations, the converged value of the ',A,'-associated increment factor <IncrFact> =',ES20.14E1)
!
!
!
 6044 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The increment factor has not converged after 50 iterations. The final error is',ES11.4,'%',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6045 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The last radius of the geometric system =',ES12.5,' differs from the specified <R_max> =',ES12.5,//, &
     &       T10,'>>>  There is something wrong in the estimation of the increment factor <IncrFact> =',ES20.14E1,//   &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6050 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>: The D',A,' distribution is logarithmic, but D',A,'_Log = 0 and this is the 1st subdivision',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6500 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: There is a problem reading the discretization axis descriptor'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6501 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The discretization axis descriptor = "',A,'" is in conflict with the expected value "',A,'"'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6502 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The indicator of discretization data = "',A,'" is in conflict with the expected value ":::>>>'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6505 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The discretization axis descriptor = "',A,'" is not an available option'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6508 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The R-discretization option is LOGARITHMIC, but <Delta_R> = 0 and no radii '/, &
     &       T10,'                       have been defined previously from which to extract a non-zero <Delta_R>'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6510 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The R-discretization option is RADII, but the <number_of_radii_R> = 0 and no radii '//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6515 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The dataset/namelist "',A,'" must be ended by the ":::<<<" descriptor'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6520 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The number of subdivisions ',A' = ',i5.5,' > the maximum declared number ',A,' = ',i5.5//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6521 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The ',A,'-discretization option is "',A,'" but the value of <read_format>=',A,' in dataset #',i3,' is invalid',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6601 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: There is a problem reading the namelist ',A,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6602 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: There is a problem reading the namelist ',A,' in dataset #',i3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6850 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!
!
!-----NOW ASSIGN ALL GEOMETRY DATA FOR A LAYER OF UNIT THICKNESS.
!
!=====PRINT ALL GEOMETRIC QUANTITIES=====
!
!
      WRITE(UNIT = *, FMT = 19)
   19 FORMAT(//, ' * * * * * M E S H  G E O M E T R Y * * * * *',10X,   &
     &           'VOLUME AND AREA ARE GIVEN FOR HEIGHT = 1 METER',//,   &
     &           ' ELEM     REL           RCON           D1             D2            V              A'/)
!
!
!
      NELEMT = NR-1
      NLayer = NZ
!
! ----------------
! ........ KP - 2/11/92
! ........ Begin change to avoid "funny" multiply-add-related
! ........ behavior on IBM RISC6000
! ----------------
!
      DO 15 I=1,NELEMT
         A(I)  = 2.0d0*PI*RC(I+1)
         REL   = SQRT(RC(I+1)*RC(I))
!
         D1(I) = REL     - RC(I)
         D2(I) = RC(I+1) - REL
!        D(I)  = (RC(I+1) - RC(I))/2.0D0
         v(i)  = pi*(rc(i+1)+rc(i))*(rc(i+1)-rc(i))
   15 CONTINUE
!
      DO I=1,NELEMT
!        REL = RC(I+1) - D(I)
         REL = SQRT(RC(I+1)*RC(I))
!
         WRITE(*,18) I,REL,RC(I+1),D1(I),D2(I),V(I),A(I)
      END DO
   18 FORMAT(' ',I4,6(2X,ES12.5))
!
!
! ----------------
! ........ KP - 2/11/92
! ........ End change
! ----------------
!
!
      REWIND (UNIT = MESH_Unit)
      WRITE(*,20)
   20 FORMAT(/' WRITE FILE "MESH" FROM INPUT DATA')
!
! - GJM: Begin modification for 8-character elements
!
      IF_Format: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
         WRITE(MESH_Unit, FMT = '(T1, A)') '>>>ELEMENTS: Cylindrical coordinates'
         WRITE(MESH_Unit, FMT = '(T1, A)') '&Units  length_units = "m" /'
      ELSE
!
         IF(ElemNameLength == 8) THEN
            WRITE(MESH_Unit,5003) NELEMT,RC(1),RC(NELEMT+1)
         ELSE
            WRITE(MESH_Unit,3) NELEMT,RC(1),RC(NELEMT+1)
         END IF
!
      END IF IF_Format
!
!
 5003 FORMAT('ELEMEext2 --- ',I5,2(2x,ES12.5))
    3 FORMAT('ELEME --- ',I5,2(2x,ES12.5))
!
! - GJM: End modification for 8-character elements
!
      RETURN
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of RZ2D
!
!
      STOP
      END SUBROUTINE RZ2D
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE WRZ2D(K)
!
         USE Grid_Generation_Parameters
         USE Het_Region_Definition
!
         USE MeshMaker_Data, H => DZ
!
!
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: emiss
!
      REAL(KIND = 8) :: emissivity
!
      CHARACTER(LEN=1) :: NII, NII1, NOV, NOV1
      CHARACTER(LEN=3) :: bound_id, exclusion_zone_type
      CHARACTER(LEN=5) :: DOM, medium
!
      INTEGER :: ICALL = 0, ij, ij8
!
      SAVE ICALL
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of WRZ2D
!
!
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(*,6000)
!6000 FORMAT(6X,'WRZ2D    1.0      26 MARCH     1991',6X,
!6000 FORMAT('WRZ2D             1.1   20 September 2000',6X,'WRITE DATA FOR 2-D R-Z MESH ON FILE *MESH*')
!6000 FORMAT('WRZ2D             1.5   11 August    2014',6X,'WRITE DATA FOR 2-D R-Z MESH ON FILE *MESH*')
 6000 FORMAT(/,'WRZ2D 2.0 .............................. 16 January   2015',6X,'WRITE DATA FOR 2-D R-Z MESH ON FILE *MESH*')
!
      NSKIP1=1
      NSKIP2=2
!
! ... Allocate memory and initialize
!
      ALLOCATE( emiss(NXMax*NZMax) )
!
      emiss = 0.0d0
!
      IF(K == 2) GO TO 1000
!
!-----COME HERE TO GENERATE MESH BY VERTICAL COLUMNS
!
      DO_ILoop1: DO I8 = NSKIP1-1, NELEMT-1
!
         I = I8+1
!
! ----------------
! ........ GJM - 9/19/2000
! ........ Begin change to assign correct location for average gridblock pressure
! ----------------
!
!        REL = RC(I+1)-D(I)
         REL = SQRT(RC(I+1)*RC(I))
!
! ...... GJM - 9/19/2000
!
!        REL = SQRT( 5.0d-1*(RC(I+1)*RC(I+1) - RC(I)*RC(I))/ DLOG((RC(I+1)/RC(I))))
!
! ----------------
! ........ GJM - 9/19/2000
! ........ End change
! ----------------
!
         ZJ  = - H(1)/2.0d0
!
! - GJM: Begin modification for 8-character elements
!
         IF(ElemNameLength == 8) THEN

            MI  = MOD(I8,1000)
            IF(MOD(I,1000) == 0) THEN
               II = I/1000
            ELSE
               II = I/1000+1
            END IF
            NII = NA(II+10)
!
            DO_JLoop1: DO J=1,NLayer
!
               IF(MOD(J,1000) == 0) THEN
                  JJ = J/1000
               ELSE
                  JJ = J/1000 + 1
               END IF
!
               MJ  = MOD(J-1,1000)
!              NOV = NA((J-1)/1000+11)
               NOV = NA(24-(J-1)/1000)
!
               IF(J > 1) ZJ = ZJ - H(J)/2.0d0 - H(J-1)/2.0d0
!
! ............ Determine if the element is to be excluded
!
               IF( Excluded_Element( coordinates, REL, 0.0d0, Z_ref+ZJ, exclusion_zone_type ) ) THEN
                  excluded( i + (j-1) * NXMax ) = .TRUE.
                  CYCLE DO_JLoop1
               END IF
!
               IF( areas_for_HeatExch_Solution ) THEN
                  AIJ = 0.0d0
                  IF(J == 1 .OR. J == NLayer) AIJ = AIJ + V(I)
               END IF
!
               VIJ = V(I)*H(J)
!
! ----------------
! ............ When the domain is homogeneous, assign the standard domain number/name
! ----------------
!
               IF_Media1: IF(Flag_HetRegions .EQV. .FALSE.) THEN
!
                  DOM = '    1'
                  emissivity = dominant_medium_emissivity
!
! ----------------
! ............ When the domain is heterogeneous ...
! ----------------
!
               ELSE
!
! ............... Assign the dominant domain number/name
!
                  IF_HetRegions1: IF(TotNum_HetRegions == 1) THEN
!
                     DOM = dominant_medium
                     emissivity = dominant_medium_emissivity
!
! ............... Assign the appropriate heterogeneous domain number/name
!
                  ELSE
                     DOM = Media_Region(coordinates,dominant_medium, REL, Z_ref+ZJ, Z_ref+ZJ, media_by_number, emissivity )
                  END IF IF_HetRegions1
!
               END IF IF_Media1
!
! ----------------
! ............ Determine whether this is a boundary cell
! ----------------
!
               IF(Flag_Boundaries) THEN
!
                  medium = DOM
!
                  CALL Boundary_Status(coordinates, REL, Z_ref+ZJ, Z_ref+ZJ, bound_id, medium, media_by_number, emissivity )
!
                  IF(bound_id(1:1) == 'I' .OR. bound_id(1:1) == 'V') THEN
                     DOM = medium
                  END IF
!
               ELSE
                  bound_id = '   '
               END IF
!
 5401 FORMAT(2(A1,I3.3),'  &Elem  V=',ES17.10,', MedNum=',A,', r=',ES17.10,', z=',ES17.10,', act="',A,'" /')
 5400 FORMAT(2(A1,I3.3),'  &Elem  V=',ES17.10,', MedNum=',A,', r=',ES17.10,', z=',ES17.10,' /')
!
 5501 FORMAT(2(A1,I3.3),'  &Elem  V=',ES17.10,', MedNam="',A,'", r=',ES17.10,', z=',ES17.10,', act="',A,'" /')
 5500 FORMAT(2(A1,I3.3),'  &Elem  V=',ES17.10,', MedNam="',A,'", r=',ES17.10,', z=',ES17.10,' /')
!
               IF( active_boundaries .AND. (bound_id(1:1) == 'I' .OR. bound_id(1:1) == 'V') ) THEN
                  VIJ = 1.0d50
                  bound_id(1:1) = '*'
               END IF
!
               IF_Format: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
                  IF(bound_id(1:1) == ' ') THEN
                     IF(DOM == '    1') THEN
                        WRITE(MESH_Unit, FMT = 5400) NOV,MJ,NII,MI,VIJ,DOM,REL,Z_ref+ZJ
                     ELSE
                        WRITE(MESH_Unit, FMT = 5500) NOV,MJ,NII,MI,VIJ,DOM,REL,Z_ref+ZJ
                     END IF
                  ELSE
                     IF(DOM == '    1') THEN
                        WRITE(MESH_Unit, FMT = 5401) NOV,MJ,NII,MI,VIJ, DOM, REL, Z_ref+ZJ, bound_id
                     ELSE
                        WRITE(MESH_Unit, FMT = 5501) NOV,MJ,NII,MI,VIJ, DOM, REL, Z_ref+ZJ, bound_id
                     END IF
                  END IF
               ELSE
                  IF( areas_for_HeatExch_Solution ) THEN
                     WRITE(MESH_Unit,5005) NII, MI, NOV, MJ, DOM, VIJ, AIJ, REL, Z_ref+ZJ, bound_id
                  ELSE
                     WRITE(MESH_Unit,5815) NII, MI, NOV, MJ, DOM, VIJ, REL, Z_ref+ZJ, bound_id
                  END IF
               END IF IF_Format
!
               ij = i + (j-1) * NXMax
               emiss(ij) = emissivity
!
            END DO DO_JLoop1
!
!
!
         ELSE
!
!
            MI = MOD(I8,100)
!
            IF(MOD(I,100) == 0) THEN
               II = I/100
            ELSE
               II = I/100+1
            END IF
!
            NII = NA(II)
!
!
            DO_JLoop2: DO J=1,NLayer
!
               MJ  = MOD(J-1,36)+1
               NOV = NA((J-1)/36+11)
!
               IF(J > 1) ZJ = ZJ - H(J)/2.0d0 - H(J-1)/2.0d0
!
! ............ Determine if the element is to be excluded
!
               IF( Excluded_Element( coordinates, REL, 0.0d0, Z_ref+ZJ, exclusion_zone_type ) ) THEN
                  excluded( i + (j-1) * NXMax ) = .TRUE.
                  CYCLE DO_JLoop2
               END IF
!
               IF( areas_for_HeatExch_Solution ) THEN
                  AIJ = 0.0d0
                  IF(J == 1 .OR. J == NLayer) AIJ = AIJ+V(I)
               END IF
!
               VIJ = V(I)*H(J)
!
! ----------------
! ............ When the domain is homogeneous, assign the standard domain number/name
! ----------------
!
               IF_Media2: IF(Flag_HetRegions .EQV. .FALSE.) THEN
!
                  DOM = '    1'
                  emissivity = dominant_medium_emissivity
!
! ----------------
! ............ When the domain is heterogeneous ...
! ----------------
!
               ELSE
!
! ............... Assign the dominant domain number/name
!
                  IF_HetRegions2: IF(TotNum_HetRegions == 1) THEN
!
                     DOM = dominant_medium
                     emissivity = dominant_medium_emissivity
!
! ............... Assign the appropriate heterogeneous domain number/name
!
                  ELSE
                     DOM = Media_Region( coordinates, dominant_medium, REL, Z_ref+ZJ, Z_ref+ZJ, media_by_number, emissivity )
                  END IF IF_HetRegions2
!
               END IF IF_Media2
!
! ----------------
! ............ Determine whether this is a boundary cell
! ----------------
!
               IF(Flag_Boundaries) THEN
!
                  medium = DOM
!
                  CALL Boundary_Status(coordinates, REL, Z_ref+ZJ, Z_ref+ZJ, bound_id, medium, media_by_number, emissivity )
!
                  IF(bound_id(1:1) == 'I' .OR. bound_id(1:1) == 'V') THEN
                     DOM = medium
                  END IF
!
               ELSE
                  bound_id = '   '
               END IF
!
 5403 FORMAT(3A1,I2.2,'  &Elem  V=',ES17.10,', MedNum=',A,', r=',ES17.10,', z=',ES17.10,', act="',A,'" /')
 5402 FORMAT(3A1,I2.2,'  &Elem  V=',ES17.10,', MedNum=',A,', r=',ES17.10,', z=',ES17.10,' /')
!
 5503 FORMAT(3A1,I2.2,'  &Elem  V=',ES17.10,', MedNam="',A,'", r=',ES17.10,', z=',ES17.10,', act="',A,'" /')
 5502 FORMAT(3A1,I2.2,'  &Elem  V=',ES17.10,', MedNam="',A,'", r=',ES17.10,', z=',ES17.10,' /')
!
               IF( active_boundaries .AND. (bound_id(1:1) == 'I' .OR. bound_id(1:1) == 'V') ) THEN
                  VIJ = 1.0d50
                  bound_id(1:1) = '*'
               END IF
!
!
               IF_Format2: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
                  IF(bound_id(1:1) == ' ') THEN
                     IF(DOM == '    1') THEN
                        WRITE(MESH_Unit, FMT = 5402) NOV,NA(MJ),NII,MI,VIJ,DOM,REL,Z_ref+ZJ
                     ELSE
                        WRITE(MESH_Unit, FMT = 5502) NOV,NA(MJ),NII,MI,VIJ,DOM,REL,Z_ref+ZJ
                     END IF
                  ELSE
                     IF(DOM == '    1') THEN
                        WRITE(MESH_Unit, FMT = 5403) NOV,NA(MJ),NII,MI,VIJ,DOM,REL,Z_ref+ZJ,bound_id
                     ELSE
                        WRITE(MESH_Unit, FMT = 5503) NOV,NA(MJ),NII,MI,VIJ,DOM,REL,Z_ref+ZJ,bound_id
                     END IF
                  END IF
               ELSE
                  IF( areas_for_HeatExch_Solution ) THEN
                     WRITE(MESH_Unit,5) NOV,NA(MJ),NII,MI,DOM,VIJ,AIJ,REL,Z_ref+ZJ,bound_id
                  ELSE
                     WRITE(MESH_Unit,5805) NOV,NA(MJ),NII,MI,DOM,VIJ,REL,Z_ref+ZJ,bound_id
                  END IF
               END IF IF_Format2
!
               ij = i + (j-1) * NXMax
               emiss(ij) = emissivity
!
            END DO DO_JLoop2
!
         END IF
!
      END DO DO_ILoop1
!
    5 FORMAT(   3A1,I2.2,10X,a5,2(ES10.4),10X,ES10.4,10X,ES10.4E1,1X,A3 )
 5805 FORMAT(   3A1,I2.2,10X,a5,1(ES10.4),20X,ES10.4,10X,ES10.4E1,1X,A3 )
 5005 FORMAT( 2(A1,I3.3), 7X,a5,2(ES10.4),10X,ES10.4,10X,ES10.4E1,1X,A3 )
 5815 FORMAT( 2(A1,I3.3), 7X,a5,1(ES10.4),20X,ES10.4,10X,ES10.4E1,1X,A3 )
!
!! - GJM: End modification for 8-character elements
!
!
!     MJ RUNS FROM 1 TO 36 FOR CONSECUTIVE LAYERS.
!     NOV RUNS THROUGH A, B, C, ... FOR THE 1ST, 2ND, 3RD, ... SET
!     OF 36 LAYERS.
!
!-----CONVENTION FOR ELEMENT NAMES
!     NOV = A, B, C, .... LABELS 1ST, 2ND, 3RD, ... SET OF 37 LAYERS.
!     NA = 0, 1, 2, 3, ... 9, A, B, C, ... LABELS LAYERS.
!     NII = 0, 1, 2, 3, ... AFTER 100, 200, 300, ... RADIAL ELEMENTS.
!     MI = 0, 1, 2, 3, ... 99; NUMBERS RADIAL ELEMENTS IN EXCESS OF
!                              FULL HUNDREDS.
!
      GO TO 2000
!
! >>>>>>>
! >>>>>>>
! >>>>>>>
!
 1000 CONTINUE
!-----COME HERE TO GENERATE MESH BY HORIZONTAL LAYERS
!
      ZJ=-H(1)/2.0d0
!
!
      DO_JLoop3: DO  J = 1, NLayer
!
         IF(ElemNameLength == 8) THEN
!
            IF(MOD(J,1000) == 0) THEN
               JJ = J/1000
            ELSE
               JJ = J/1000 + 1
            END IF
!
            MJ  = MOD(J-1,1000)
!           NOV = NA((J-1)/1000+11)
            NOV = NA(24-(J-1)/1000)
!
         ELSE
!
            MJ  = MOD(J-1,36)+1
            NOV = NA((J-1)/36+11)
!
         END IF
!
         IF(J > 1) ZJ = ZJ - H(J)/2.0d0 - H(J-1)/2.0d0
!
!        DO 3 I = NSKIP1,NELEMT
!
         DO_ILoop3: DO I8 = NSKIP1-1,NELEMT-1
!           MI  = MOD(I,100)
            I   = I8+1
!
            IF(ElemNameLength == 8) THEN
               MI  = MOD(I8,1000)
               IF(MOD(I,1000) == 0) THEN
                  II = I/1000
               ELSE
                  II = I/1000+1
               END IF
               NII = NA(II+10)
            ELSE
               MI  = MOD(I8,100)
               IF(MOD(I,100) == 0) THEN
                  II = I/100
               ELSE
                  II = I/100+1
               END IF
               NII = NA(II)
            END IF
!           IF(II.GE.1) NII = NA(II)
!
! ----------------
! ........ GJM - 9/19/2000
! ........ Begin change to assign correct location for average gridblock pressure
! ----------------
!
!           REL = RC(I+1)-D(I)
            REL = SQRT(RC(I+1)*RC(I))
!
! ......... Determine if the element is to be excluded
!
            IF( Excluded_Element( coordinates, REL, 0.0d0, Z_ref+ZJ, exclusion_zone_type ) )  THEN
               excluded( i + (j-1) * NXMax ) = .TRUE.
               CYCLE DO_ILoop3
            END IF
!
! ----------------
! ........ GJM - 9/19/2000
! ........ End change
! ----------------
!
            IF( areas_for_HeatExch_Solution ) THEN
               AIJ = 0.0d0
               IF(J == 1 .OR. J == NLayer) AIJ = AIJ+V(I)
            END IF
!
            VIJ = V(I)*H(J)
!
! - GJM: Begin modification for 8-character elements
!
!
! -------------
! ......... When the domain is homogeneous, assign the standard domain number/name
! -------------
!
            IF_Media3: IF(Flag_HetRegions .EQV. .FALSE.) THEN
!
               DOM = '    1'
               emissivity = dominant_medium_emissivity
!
! ......... When the domain is heterogeneous ...
!
            ELSE
!
! ............ Assign the dominant domain number/name
!
               IF_HetRegions3: IF(TotNum_HetRegions == 1) THEN
!
                  DOM = dominant_medium
                  emissivity = dominant_medium_emissivity
!
! ............ Assign the appropriate heterogeneous domain number/name
!
               ELSE
                  DOM = Media_Region(coordinates,dominant_medium,REL,Z_ref+ZJ,Z_ref+ZJ, media_by_number, emissivity )
               END IF IF_HetRegions3
!
            END IF IF_Media3
!
! -------------
! ............ Determine whether this is a boundary cell
! -------------
!
            IF(Flag_Boundaries) THEN
!
               medium = DOM
!
               CALL Boundary_Status(coordinates,REL,Z_ref+ZJ,Z_ref+ZJ,bound_id,medium, media_by_number, emissivity )
!
               IF(bound_id(1:1) == 'I' .OR. bound_id(1:1) == 'V') THEN
                  DOM = medium
               END IF
!
            ELSE
               bound_id = '   '
            END IF
!
!
               IF( active_boundaries .AND. (bound_id(1:1) == 'I' .OR. bound_id(1:1) == 'V') ) THEN
                  VIJ = 1.0d50
                  bound_id(1:1) = '*'
               END IF
!
!
            IF_Format3: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
!
               IF(bound_id(1:1) == ' ') THEN
                  IF(ElemNameLength == 8) THEN
                     IF(DOM == '    1') THEN
                        WRITE(MESH_Unit, FMT = 5400) NII,MI,NOV,MJ,VIJ,DOM,REL,Z_ref+ZJ
                     ELSE
                        WRITE(MESH_Unit, FMT = 5500) NII,MI,NOV,MJ,VIJ,DOM,REL,Z_ref+ZJ
                     END IF
                  ELSE
                     IF(DOM == '    1') THEN
                        WRITE(MESH_Unit, FMT = 5402) NII,MI,NOV,NA(MJ),VIJ,DOM,REL,Z_ref+ZJ
                     ELSE
                        WRITE(MESH_Unit, FMT = 5502) NII,MI,NOV,NA(MJ),VIJ,DOM,REL,Z_ref+ZJ
                     END IF
                  END IF
               ELSE
                  IF(ElemNameLength == 8) THEN
                     IF(DOM == '    1') THEN
                        WRITE(MESH_Unit, FMT = 5401) NOV,MJ,NII,MI,VIJ,DOM,REL,Z_ref+ZJ,bound_id
                     ELSE
                        WRITE(MESH_Unit, FMT = 5501) NOV,MJ,NII,MI,VIJ,DOM,REL,Z_ref+ZJ,bound_id
                     END IF
                  ELSE
                     IF(DOM == '    1') THEN
                        WRITE(MESH_Unit, FMT = 5403) NOV,NA(MJ),NII,MI,VIJ,DOM,REL,Z_ref+ZJ,bound_id
                     ELSE
                        WRITE(MESH_Unit, FMT = 5503) NOV,NA(MJ),NII,MI,VIJ,DOM,REL,Z_ref+ZJ,bound_id
                     END IF
                  END IF
               END IF
!
            ELSE
!
               IF(ElemNameLength == 8) THEN
                  IF( areas_for_HeatExch_Solution ) THEN
                     WRITE(MESH_Unit,5005) NII,MI,NOV,MJ,DOM,VIJ,AIJ,REL,Z_ref+ZJ,bound_id
                  ELSE
                     WRITE(MESH_Unit,5815) NII,MI,NOV,MJ,DOM,VIJ,REL,Z_ref+ZJ,bound_id
                  END IF
               ELSE

                  IF( areas_for_HeatExch_Solution ) THEN
                     WRITE(MESH_Unit,5) NOV,NA(MJ),NII,MI,DOM,VIJ,AIJ,REL,Z_ref+ZJ,bound_id
                  ELSE
                     WRITE(MESH_Unit,5805) NOV,NA(MJ),NII,MI,DOM,VIJ,REL,Z_ref+ZJ,bound_id
                  END IF
               END IF
!
            END IF IF_Format3
!
            ij = i + (j-1) * NXMax
            emiss(ij) = emissivity
!
! - GJM: End modification for 8-character elements
!
         END DO DO_ILoop3
!
      END DO DO_JLoop3
!
!     MJ RUNS FROM 1 TO 35 FOR CONSECUTIVE LAYERS.
!     NOV RUNS THROUGH A, B, C, ... FOR THE 1ST, 2ND, 3RD, ... SET
!     OF 35 LAYERS.
!
 2000 CONTINUE
!
!
      IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
         WRITE(MESH_Unit, FMT = 5551)
       ELSE
         WRITE(MESH_Unit, FMT = 6)
      END IF
!
    6 FORMAT('     ')
 5551 FORMAT(T1,'<<<End of the ELEMENTS data block',/,'     ')
 5552 FORMAT(T1,'<<<End of the CONNECTIONS data block',/,'     ')
!
!
!
      IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
         WRITE(MESH_Unit, FMT = '(T1, A)') '>>>CONNECTIONS: Cartesian coordinates'
         WRITE(MESH_Unit, FMT = '(T1, A)') '&Units  length_units = "m" /'
      ELSE
         WRITE(MESH_Unit, FMT = 7)
      END IF !
!
    7 FORMAT('CONNE')
!
!-----ASSIGN HORIZONTAL CONNECTIONS.
!
!     DO 8 I=NSKIP2,NELEMT
      DO_ILoop4: DO I8=NSKIP2-1,NELEMT-1
!
         I = I8+1
!
         IF(ElemNameLength == 8) THEN
            MI  = MOD(I8,1000)
            MI1 = MOD(I8-1,1000)
            IF(MOD(I,1000) == 0) THEN
               II  = I/1000
               II1 = I8/1000
            ELSE
               II  = I/1000+1
               II1 = I8/1000+1
            END IF
            NII  = NA(II+10)
            NII1 = NA(II1+10)
         ELSE
            MI  = MOD(I8,  100)
            MI1 = MOD(I8-1,100)
!
            IF(MOD(I,100) == 0) THEN
               II = I/100
            ELSE
               II = I/100+1
            END IF
!
            IF(MOD(I8,100) == 0) THEN
               II1 = I8/100
            ELSE
               II1 = I8/100+1
            END IF
!
            NII  = NA(II)
            NII1 = NA(II1)
         END IF
!
! - GJM: Begin modification for 8-character elements
!
         IF(ElemNameLength == 8) THEN
!
            DO_JLoop4: DO J=1,NLayer
!
               ij  = i + (j-1) * NXMax
               ij8 = ij - 1
!
! ............ Determine if any of the connection elements are excluded
!
               IF( excluded(i + (j-1) * NXMax) .OR. excluded(i8 + (j-1) * NXMax) ) THEN
                  CYCLE DO_JLoop4
               END IF
!
               IF(MOD(J,1000) == 0) THEN
                  JJ = J/1000
               ELSE
                  JJ = J/1000 + 1
               END IF
!
               MJ  = MOD(J-1,1000)
               NOV = NA(24-(J-1)/1000)
!
               AIJ = A(I-1)*H(J)
!
               IF( emiss(ij) /= emiss(ij8) ) emiss(ij) = 0.0d0
!
               IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
                  WRITE(MESH_Unit,5520) NII1,MI1,NOV,MJ,NII,MI,NOV,MJ, AIJ, D2(I-1),D1(I),1,emiss(ij)
               ELSE
                  WRITE(UNIT = MESH_Unit, FMT = 5009, ADVANCE = 'NO' ) NII1,MI1,NOV,MJ,NII,MI,NOV,MJ, D2(I-1), D1(I), AIJ, 0.0e0
                  IF( ABS(emiss(ij)) > 0.0e0 ) THEN
                     WRITE(UNIT = MESH_Unit, FMT = '(ES10.4E1)' ) emiss(ij)
                  ELSE
                     WRITE(UNIT = MESH_Unit, FMT = '(A)' ) '          '
                  END IF
               END IF
!
            END DO DO_JLoop4
!
         ELSE
!
            DO_JLoop5: DO J=1,NLayer
!
               ij  = i + (j-1) * NXMax
               ij8 = ij - 1
!
! ............ Determine if any of the connection elements are excluded
!
               IF( excluded(ij) .OR. excluded(ij8) ) THEN
                  CYCLE DO_JLoop5
               END IF
!
               MJ  = MOD(J-1,36)+1
               NOV = NA((J-1)/36+11)
               AIJ = A(I-1)*H(J)
!
               IF( emiss(ij) /= emiss(ij8) ) emiss(ij) = 0.0d0
!
               IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
                  WRITE(MESH_Unit,5530) NOV,NA(MJ),NII1,MI1, NOV,NA(MJ),NII,MI, AIJ, D2(I-1),D1(I),1,emiss(ij)
               ELSE
                  WRITE(UNIT = MESH_Unit, FMT = 9, ADVANCE = 'NO' ) NOV,NA(MJ),NII1,MI1, NOV,NA(MJ),NII,MI,D2(I-1),D1(I),AIJ, 0.0e0
                  IF( ABS(emiss(ij)) > 0.0e0 ) THEN
                     WRITE(UNIT = MESH_Unit, FMT = '(ES10.4E1)' ) emiss(ij)
                  ELSE
                     WRITE(UNIT = MESH_Unit, FMT = '(A)' ) '          '
                  END IF
               END IF
!
            END DO DO_JLoop5
!
         END IF
!
      END DO DO_ILoop4
!
 5530 FORMAT(2( 3(a1),i2.2 ),'  &Conx  A=',ES17.10,', d1=',ES15.8,', d2=',ES15.8,', dir=',i1,' /')
 5520 FORMAT(4( a1,i3.3 ),'  &Conx  A=',ES17.10,', d1=',ES15.8,', d2=',ES15.8,', dir=',i1,', emsv=',ES15.8,' /')
!
 5009 FORMAT(4( a1,i3.3 ),13X,'1',4(ES10.4E1))
    9 FORMAT(2( 3(a1),i2.2 ),19X,'1',4(ES10.4E1))
!
! - GJM: End modification for 8-character elements
!
!     MJ RUNS FROM 1 TO 36 FOR CONSECUTIVE LAYERS.
!     NOV RUNS THROUGH A, B, C, ... FOR THE 1ST, 2ND, 3RD, ... SET
!     OF 35 LAYERS.
!
!
      IF(NLayer <= 1) GO TO 2500
!
!-----NOW FOR THE VERTICAL CONNECTIONS.
!
!     DO 40 I=NSKIP1,NELEMT
      DO_ILoop5: DO I8=NSKIP1-1,NELEMT-1

         I = I8+1
!
         IF(ElemNameLength == 8) THEN
            MI  = MOD(I8,1000)
            IF(MOD(I,1000) == 0) THEN
               II  = I/1000
            ELSE
               II  = I/1000+1
            END IF
            NII  = NA(II+10)
         ELSE
            MI  = MOD(I8,  100)
            IF(MOD(I,100) == 0) THEN
               II  = I/100
            ELSE
               II  = I/100+1
            END IF
            NII  = NA(II)
         END IF
!
         AIJ = V(I)
!
! - GJM: Begin modification for 8-character elements
!
         IF(ElemNameLength == 8) THEN
            DO_JLoop6: DO J=2,NLayer
!
               J1 = J-1
!
               ij  = i + (j-1) * NXMax
               ij8 = ij - NXMax
!
! ............ Determine if any of the connection elements are excluded
!
               IF( excluded(ij) .OR. excluded(ij8)) THEN
                  CYCLE DO_JLoop6
               END IF

               IF(MOD(J1,1000) == 0) THEN
                  JJ = J1/1000
               ELSE
                  JJ = J1/1000 + 1
               END IF
!
               MJ1  = MOD(J1-1,1000)
!              NOV1 = NA((J1-1)/1000+11)
               NOV1 = NA(24-(J1-1)/1000)
!
               IF(MOD(J,1000) == 0) THEN
                  JJ = J/1000
               ELSE
                  JJ = J/1000 + 1
               END IF
!
               MJ  = MOD(J1,1000)
!              NOV = NA((J1)/1000+11)
               NOV = NA(24-(J1)/1000)
!
               D1U  = H(J1)/2.0d0
               D2U  = H(J)/2.0d0
!
               IF( emiss(ij) /= emiss(ij8) ) emiss(ij) = 0.0d0
!
               IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
                  WRITE(MESH_Unit,5521) NOV1, MJ1,NII,MI, NOV, MJ, NII, MI, AIJ, D1U, D2U, 3, 1.0d0, emiss(ij)
               ELSE
                  WRITE(UNIT = MESH_Unit, FMT = 5041, ADVANCE = 'NO' ) NII, MI, NOV1,MJ1,NII,MI,NOV, MJ, D1U, D2U, AIJ, 1.0e0
                  IF( ABS(emiss(ij)) > 0.0e0 ) THEN
                     WRITE(UNIT = MESH_Unit, FMT = '(ES10.4E1)' ) emiss(ij)
                  ELSE
                     WRITE(UNIT = MESH_Unit, FMT = '(A)' ) '          '
                  END IF
               END IF
!
            END DO DO_JLoop6
!
         ELSE
!
            DO_JLoop7: DO J=2,NLayer
!
               J1 = J-1
!
               ij  = i + (j-1) * NXMax
               ij8 = ij - NXMax
!
! ............ Determine if any of the connection elements are excluded
!
               IF( excluded(ij) .OR. excluded(ij8)) THEN
                  CYCLE DO_JLoop7
               END IF
!
               MJ   = MOD(J1,36)+1
               NOV  = NA((J1)/36+11)
!
               MJ1  = MOD(J1-1,36)+1
               NOV1 = NA((J1-1)/36+11)
!
               D1U  = H(J1)/2.0d0
               D2U  = H(J)/2.0d0
!
               IF( emiss(ij) /= emiss(ij8) ) emiss(ij) = 0.0d0
!
               IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
                  WRITE(MESH_Unit,5531) NOV1,NA(MJ1),NII,MI, NOV,NA(MJ),NII,MI,AIJ,D1U,D2U,3,1.0e0, emiss(ij)
               ELSE
                  WRITE(UNIT = MESH_Unit, FMT = 41, ADVANCE = 'NO' ) NOV1,NA(MJ1),NII,MI,NOV,NA(MJ),NII,MI,D1U,D2U,AIJ,1.0e0
                  IF( ABS(emiss(ij)) > 0.0e0 ) THEN
                     WRITE(UNIT = MESH_Unit, FMT = '(ES10.4E1)' ) emiss(ij)
                  ELSE
                     WRITE(UNIT = MESH_Unit, FMT = '(A)' ) '          '
                  END IF
               END IF
!
            END DO DO_JLoop7
!
         END IF
!
      END DO DO_ILoop5
!
!
 5531 FORMAT(2( 3(a1),i2.2 ),'  &Conx  A=',ES17.10,', d1=',1ES15.8,', d2=',1ES15.8,', dir=',i1,', beta=',ES17.10,' /')
 5521 FORMAT(4( a1,i3.3 ),'  &Conx  A=',ES17.10,', d1=',1ES15.8,', d2=',1ES15.8,', dir=',i1,', beta=',ES17.10,', emsv=',ES15.8,' /')
!
 5041 FORMAT(4( a1,i3.3 ),13X,'3',4(ES10.4E1))
   41 FORMAT(2( 3(a1),i2.2 ),19X,'3',4(ES10.4E1))
!
! - GJM: End modification for 8-character elements
!
 2500 CONTINUE
!
!
!
      IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
         WRITE(MESH_Unit, FMT = 5555)
       ELSE
         WRITE(MESH_Unit, FMT = 6)
      END IF
!
 5555 FORMAT(T1,'<<<End of the CONNECTIONS data block',/,'     ')
!
!
!
      REWIND (UNIT = MESH_Unit)
      READ(MESH_Unit,43) DOM
!
      L=0
!
! - GJM: Begin modification for 8-character elements
!
      DO 44 I = NSKIP1,NELEMT
         DO 44 J = 1,NLayer
            L = L+1
            IF(ElemNameLength == 8) THEN
               READ(MESH_Unit,5043) name(L)
            ELSE
               READ(MESH_Unit,43) name(L)(1:5)
            END IF
   44 CONTINUE
!
 5043 FORMAT(A8)
   43 FORMAT(A5)
!
! - GJM: End modification for 8-character elements
!
      CALL PRZ2D(K,NSKIP1,NELEMT,NLayer)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of WRZ2D
!
!
      RETURN
!
      END SUBROUTINE WRZ2D
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE PRZ2D(KK,NSKIP1,NELEMT,NLAY)
!
         USE MeshMaker_Data
         USE Grid_Generation_Parameters, ONLY: NA
!
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      SAVE ICALL
      DATA ICALL/0/
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of PRZ2D
!
!
!
!     FIND LOCATION OF ELEMENT NAME WITH R-INDEX I, Z-INDEX K.
!     MESH ORGANIZED "BY COLUMNS" (KK = 1)
!
      NIK1(I,K)=(I-1)*NLAY+K
!
!     MESH ORGANIZED "BY LAYERS" (KK = 2)
!
      NIK2(I,K)=(K-1)*(NELEMT-NSKIP1+1)+I
!
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(UNIT = *, FMT = 6000)
! 899 FORMAT('PRZ2D             1.0   27 MARCH     1991',6X,'MAKE STRUCTURED PRINTOUT OF 2-D R-Z MESH')
 6000 FORMAT(/,'PRZ2D 1.5 .............................. 29 July      2014',6X,'MAKE STRUCTURED PRINTOUT OF 2-D R-Z MESH')
!
      NRZ = NELEMT*NLAY
      WRITE(*,31) NELEMT,NLAY,NRZ
   31 FORMAT(/,' ',131('*')/' *',20X,'2-D R-Z MESH WITH NR*NLAY = ',I4,' *',I4,'  =',I7,' GRID BLOCKS',52X,'*'/' ',131('*'))
!
      WRITE(*,18) NLAY,NELEMT
   18 FORMAT(' *',129X,'*'/   &
     &   ' *',20X,'THE MESH WILL BE PRINTED AS VERTICAL SLICES',66X,   &
     &   '*'/' *',129X,'*'/   &
     &   ' *',20X,'LAYERS GO FROM K = 1 TO K = NLAY =',I4,71X,'*'/   &
     &   ' *',129X,'*'/   &
     &   ' *',20X,'RADIAL GRID BLOCKS GO IN COLUMNS FROM I = 1 TO',   &
     &   ' I = NR =',I4,50X,'*'/   &
     &   ' *',129X,'*'/' ',131('*')/)
!
!
! - GJM: Begin modification for 8-character elements
!
      IF(ElemNameLength == 8) THEN
         WRITE(*,5004)
         DO K=1,NLAY
            IF(KK == 1) WRITE(*,5006) K,(name(NIK1(I,K)), I=NSKIP1,NELEMT)
            IF(KK == 2) WRITE(*,5006) K,(name(NIK2(I,K)), I=NSKIP1,NELEMT)
         END DO
      ELSE
         WRITE(*,4)
         DO K=1,NLAY
            IF(KK == 1) WRITE(*,6) K,(name(NIK1(I,K))(1:5), I=NSKIP1,NELEMT)
            IF(KK == 2) WRITE(*,6) K,(name(NIK2(I,K))(1:5), I=NSKIP1,NELEMT)
         END DO
      END IF
!
 5004 FORMAT('   COLUMN I =    1        2        3        4        ',   &
     &       '5        6        7        8        9       10',   &
     &       '       11       12'/' LAYERS')
    4 FORMAT('   COLUMN I =  1     2     3     4     5     6     ',   &
     &       '7     8     9    10    11    12    13    14    15    ',   &
     &       '16    17    18    19    20'/' LAYERS')
 5006 FORMAT('  K = ',I4,2X,12(1X,A8),/,(12X,12(1X,A8)),/,(12X,12(1X,A8)),/,(12X,12(1X,A8)))
    6 FORMAT('  K = ',I4,2X,20(1X,A5)/(12X,20(1X,A5)))
!
! - Describe the discretization in the Z-direction
!
      WRITE(*,44)
!
      WRITE(*,43) NLAY,(NA(I),I=2,10),(DZ(I),I=1,NLAY)
!
      WRITE(*,443) NLAY,(NA(I),I=2,10),(Z_ref - SUM(DZ(1:I)),I=1,NLAY)
!
      WRITE(*,444)
      WRITE(*,445) Z_ref-DZ(1), Z_ref, 1
      DO i = 2,NLAY
         WRITE(*,445) Z_ref-SUM(DZ(1:i)), Z_ref-SUM(DZ(1:i-1)), i
      END DO
!
      WRITE(*,44)
!
   43 FORMAT(//'   DELTA-Z  (',I6,' INCREMENTS)'//13X,8(A1,11X),A1,10X,'10'//(6X,10(1X,ES11.4)))
   44 FORMAT(/' ',131('*'))
!
  443 FORMAT(//'   Z  (',I6,' COORDINATES)'//13X,8(A1,13X),A1,12X,'10'//(6X,10(1X,ES13.6)))
  444 FORMAT(//'   Z  LAYERS - Beginning from the top, moving downward:')
  445 FORMAT(T5,2(1X,ES15.8),3x,'! Layer ',i5.5)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of PRZ2D
!
!
      RETURN
      END SUBROUTINE PRZ2D
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE GXYZ
!
         USE Grid_Generation_Parameters
         USE Het_Region_Definition
!
         USE MeshMaker_Data
!
!        USE IFPORT
!
      IMPLICIT NONE
!
! -------
! ... Double precision arrays
! -------
!
      REAL(KIND = 8), DIMENSION(3) :: DDXYZ, RRR, Lref
      REAL(KIND = 8), DIMENSION(2) :: D1d, D2d, D3d
      REAL(KIND = 8), DIMENSION(2) :: coeff
!
      REAL(KIND = 8), DIMENSION(Num_IntTables,3) :: available_InterpData
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: volume, emissivity, emissivity_p
      REAL(KIND = 8) :: RAD, BET, BETA, phi_x, phi_y, beta_x, beta_y, beta_z
      REAL(KIND = 8) :: CBi, CMi, CSm
      REAL(KIND = 8) :: Z_00, Z_xp, Z_yp, X_00, X_xp, Y_00, Y_yp, angle
!
! -------
! ... Integer arrays
! -------
!
      INTEGER, DIMENSION(3) :: NLXYZ
      INTEGER, DIMENSION(2) :: output_file_unit
!
! -------
! ... Integer variables
! -------
!
      INTEGER :: grid_type
!
      INTEGER :: NX = 0, NY = 0, NZ = 0, NN, NXNY, ier, error_code, n, n1, m
      INTEGER :: I, J, K, L, ILOOP, coord_order, jbegin, jend, kbegin, kend, imin, imax
      INTEGER :: itrue, jtrue, ktrue, ijk, ipjk, ijpk, ijk_min, ijk_max, ijk_min_IE, ijk_max_IE
      INTEGER :: kx_offset, ky_offset, k_xp, k_yp, ijk_array_size, Max_ijk_location
!
      INTEGER :: NumElem, NumInactElem, NumTotalElem, NumElem_Bound1, NumElem_Bound2, NumConx
!
! -------
! ... Character variables
! -------
!
      CHARACTER( LEN = ElemNameLength ) :: elem_name
!
      CHARACTER(LEN = 1) :: index
!
      CHARACTER(LEN = 2) :: connectivity
      CHARACTER(LEN = 3) :: exclusion_zone_type
      CHARACTER(LEN = 3) :: subdomain_id
!
      CHARACTER(LEN =  5) :: DOM, medium
      CHARACTER(LEN =  6) :: first_part
      CHARACTER(LEN = 23) :: discretization_descriptor
!
      CHARACTER(LEN =  90) :: line
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: valid_element, outside_I_Range_limits
!
      LOGICAL :: accepted_ijk, accepted_ipjk, accepted_ijpk, accepted_ijkp, accepted_ip_true, accepted_jp_true, accepted_kp_true, excluded_ijk
!
      LOGICAL :: is_boundary, is_inclusion, discontinue_smoothing, exists
!
      LOGICAL, DIMENSION(Num_IntTables) :: completed_interpolation
!
      LOGICAL :: First_call = .TRUE.
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call, NXNY, grid_type
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of GXYZ
!
!
      IF(First_call) THEN
!
         WRITE( UNIT = *, FMT = 6000 )
!
         First_call = .FALSE.
!
         SELECT CASE(coordinates(1:2))
         CASE('CA', 'Ca', 'ca')
            grid_type = 1
         CASE('TR', 'Tr', 'tr')
            grid_type = 2
         END SELECT
!
      END IF
!
!
!
      ALLOCATE( DX(NXmax), DY(NYmax), DZ(NZmax), DLXYZ(MAX(NXmax,NYmax,NZmax),3), STAT = error_code )
!
      IF(error_code /= 0) THEN
         WRITE(UNIT = *, FMT = 6850) 'DX, DY, DZ, DLXYZ', 'GXYZ'
         STOP
      END IF
!
!
      OPEN( UNIT = Elem_Coordinates_Unit, FILE = 'ElemCoordinates', IOSTAT = error_code )
!
! ... Open file to store temporarily the element numbers and types - to be used for connection creation
!
      IF(error_code /= 0) THEN
         WRITE(UNIT = *, FMT = 6602) '<ElemCoordinates>'
         STOP
      END IF
!
!
!
      DO_NumActDir: DO n = 1, 3
!
! >>>>>>
! >>>>>> Discretize the axes of the Cartesian system or full cylindrical system
! >>>>>>
!
         READ( UNIT = *, FMT = '(a22)', ADVANCE = 'NO', IOSTAT = error_code ) discretization_descriptor
!
         IF(error_code /= 0) THEN
            WRITE(UNIT = *, FMT = 6500)
            STOP
         END IF
!
         SELECT CASE (discretization_descriptor(1:11))
!
         CASE(':::>>>X-DIS', ':::>>>X-Dis', ':::>>>X-dis', ':::>>>x-dis')
!
            IF( grid_type == 2 ) THEN
               WRITE(*,6501) 'X','Transformed'
               STOP
            END IF
!
            CALL Discretize_Axis( n = n, direction = 'X', DL = DX, NL = NX )
!
         CASE(':::>>>R-DIS', ':::>>>R-Dis', ':::>>>R-dis', ':::>>>r-dis')
!
            IF( grid_type == 1 ) THEN
               WRITE(*,6501) 'R','Cartesian'
               STOP
            END IF
!
            CALL Discretize_Axis( n = n, direction = 'R', DL = DX, NL = NX )
!
         CASE(':::>>>Y-DIS',':::>>>Y-Dis', ':::>>>Y-dis', ':::>>>y-dis')
!
            IF( grid_type == 2 ) THEN
               WRITE(*,6501) 'Y','Transformed'
               STOP
            END IF
!
            CALL Discretize_Axis( n = n, direction = 'Y', DL = DY, NL = NY )
!
         CASE(':::>>>Th-DI',':::>>>Th-Di', ':::>>>Th-di', ':::>>>th-di')
!
            IF( grid_type == 1 ) THEN
               WRITE(*,6501) 'Th','Cartesian'
               STOP
            END IF
!
            CALL Discretize_Axis( n = n, direction = 'T', DL = DY, NL = NY )
!
         CASE(':::>>>Z-DIS',':::>>>Z-Dis', ':::>>>Z-dis', ':::>>>z-dis')
!
            CALL Discretize_Axis( n = n, direction = 'Z', DL = DZ, NL = NZ )
!
         CASE DEFAULT
!
            WRITE(UNIT = *, FMT = 6505) TRIM(ADJUSTL(discretization_descriptor))
            STOP
!
         END SELECT
!
! >>>>>>
! >>>>>>
! >>>>>>
!
      END DO DO_NumActDir
!
! >>>>>>>>
! >>>>>>>>
! >>>
! >>> For grid continuity (if specified), determine location of a reference surface (top or bottom)
! >>>
! >>>>>>>>
! >>>>>>>>
!
      IF( ensure_grid_continuity ) THEN
!
         CALL Determine_Reference_Surface( grid_type )
!
      END IF
!
! <<<
! <<<
! <<<
!
      RAD  = PI*DEG/180.
      BET  = SIN(RAD)
      BETA = COS(RAD)
!
! ...
! ... MAKE PRINTOUT OF MESH SPECIFICATIONS.
! ...
!
      NN = NX * NY * NZ
!
      IF( grid_type == 2 ) THEN
!
         IF( .NOT. ensure_grid_continuity) CALL Cylindrical_to_Cartesian
!
         WRITE(*,440) NX,NY,NZ,NN
!
         WRITE(*,46)  NX,(NA(I),I=2,10),(DX(I),I=1,NX)
!
         WRITE(*,446) NX,(NA(I),I=2,10),(X_ref + SUM(DX(1:I)),I=1,NX)
!
         WRITE(*,48)  NY,(NA(I),I=2,10),(DY(I),I=1,NY)
!
         WRITE(*,448) NY,(NA(I),I=2,10),( SUM(DY(1:I)),I=1,NY)
!
      ELSE
!
         WRITE(*,40) NX,NY,NZ,NN,DEG
!
         WRITE(*,41) NX,(NA(I),I=2,10),(DX(I),I=1,NX)
!
         WRITE(*,441) NX,(NA(I),I=2,10),(X_ref + SUM(DX(1:I)),I=1,NX)
!
         WRITE(*,42) NY,(NA(I),I=2,10),(DY(I),I=1,NY)
!
         WRITE(*,442) NY,(NA(I),I=2,10),(Y_ref + SUM(DY(1:I)),I=1,NY)
!
      END IF
!
      WRITE(*,43) NZ,(NA(I),I=2,10),(DZ(I),I=1,NZ)
!
      WRITE(*,443) NZ,(NA(I),I=2,10),(Z_ref - SUM(DZ(1:I)),I=1,NZ)
!
      WRITE(*,444)
      WRITE(*,445) Z_ref-DZ(1), Z_ref, 1
      DO i = 2,NZ
         WRITE(*,445) Z_ref-SUM(DZ(1:i)), Z_ref-SUM(DZ(1:i-1)), i
      END DO
!
      WRITE(*,44)
!
! >>>>>>>>
! >>>>>>>>
! >>>
! >>> Determine optimal numbering to minimize matrix bandwidth
! >>>
! >>>>>>>>
! >>>>>>>>
!
      CALL Determine_Optimal_Ordering
!
! <<<
! <<<
! <<<
!
! ... Counter initialization
!
      NumElem      = 0  ! ... Initializing the number of active elements
      NumInactElem = 0  ! ... Initializing the number of inactive elements
      NumConx      = 0  ! ... Initializing the number of connections
      NumTotalElem = 0  ! ... Initializing the total number of elements
!
      NumElem_Bound1 = 0  ! ... Initializing the number of boundary elements at L = Lmin
      NumElem_Bound2 = 0  ! ... Initializing the number of boundary elements at L = Lmax
!
      NXNY = NX * NY
!
! ----------------
! ........ GJM - 9/19/2000
! ........ End addition for optimal numbering to minimize bandwidth
! ----------------
!
      SELECT CASE (level_of_grid_generation)
      CASE('F', 'f')
!
         output_file_unit(1) = MESH_unit
         output_file_unit(2) = MESH_unit
!
         IF( ensure_grid_continuity ) THEN
            output_file_unit(1) = elem_file_unit
            output_file_unit(2) = conx_file_unit
         END IF
!
      CASE('P', 'p')
!
         IF( ensure_grid_continuity ) THEN
            output_file_unit(1) = TempElem_file_unit
            output_file_unit(2) = TempConx_file_unit
         ELSE
            output_file_unit(1) = elem_file_unit
            output_file_unit(2) = conx_file_unit
         END IF
!
      END SELECT
!
      REWIND ( UNIT = output_file_unit(1) )
      REWIND ( UNIT = output_file_unit(2) )
!
      ILOOP = 0
!
   30 CONTINUE
!
      ILOOP = ILOOP+1
!
! >>>
! >>> The Element loop
! >>>
!
      IF_Loop1: IF(ILOOP == 1) THEN
!
         IF_Format0: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
            WRITE(output_file_unit(1), FMT = '(T1, A)') '>>>ELEMENTS: Cartesian coordinates'
            WRITE(output_file_unit(1), FMT = '(T1, A)') '&Units  length_units = "m" /'
         ELSE
!
            IF(ElemNameLength == 8) THEN
               IF(ILOOP == 1) THEN
                  WRITE(UNIT = output_file_unit(1), FMT = 5020)
                  IF( (elem_file_unit /= output_file_unit(1)) .AND. (level_of_grid_generation == 'P' .OR. level_of_grid_generation == 'p') ) WRITE(UNIT = elem_file_unit, FMT = 5020)
               END IF
            ELSE
               IF(ILOOP == 1) THEN
                  WRITE(UNIT = output_file_unit(1), FMT = 20)
                  IF( (elem_file_unit /= output_file_unit(1)) .AND. (level_of_grid_generation == 'P' .OR. level_of_grid_generation == 'p') ) WRITE(UNIT = elem_file_unit, FMT = 20)
               END IF
            END IF
!
         END IF IF_Format0
!
         IF(partial_processing) THEN
            IF( LLoop_begin == 1 ) THEN
               imin = 1
            ELSE
               imin = LLoop_begin - 1
            END IF
            imax = LLoop_end
         ELSE
            imin = 1
            imax = NLXYZ(3)
         END IF
!
         ALLOCATE( jmin(imin:imax), jmax(imin:imax), kmin(imin:imax,NLXYZ(2)), kmax(imin:imax,NLXYZ(2)), STAT=ier )
!
         IF(ier /= 0) THEN
            WRITE(UNIT = *, FMT = 6850) 'jmin, jmax, kmin, kmax', 'GXYZ'
            STOP
         END IF
!
         jmin = NLXYZ(2)
         jmax = 1
!
         kmin = NLXYZ(1)
         kmax = 1
!
      END IF IF_Loop1
!
! >>>
! >>> The Connection loop
! >>>
!
      IF_Loop2: IF(ILOOP == 2) THEN
!
         n1 = 0
!
         IF_Format2: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
            WRITE(UNIT = output_file_unit(2), FMT = '(T1, A)') '>>>CONNECTIONS: Cartesian coordinates'
            WRITE(UNIT = output_file_unit(2), FMT = '(T1, A)') '&Units  length_units = "m" /'
         ELSE
            WRITE(UNIT = output_file_unit(2), FMT = 21)
         END IF IF_Format2
!
! ...... Determine the size of the temporary arrays
!
         IF( ensure_grid_continuity .AND. LLoop_end /= NLXYZ(3) ) THEN
            ijk_array_size = NumTotalElem + NLXYZ(2) * NLXYZ(1)
         ELSE
            ijk_array_size = NumTotalElem
         END IF
!
         IF_DXY_adj: IF( ensure_grid_continuity ) THEN
!
            ALLOCATE( DX_adj(ijk_array_size), DY_adj(ijk_array_size), STAT=ier )
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6850) 'DX_adj, DY_adj', 'GXYZ'
               STOP
            END IF
!
         END IF IF_DXY_adj
!
         ALLOCATE( ijk_location(ijk_array_size), STAT=ier )
!
         IF(ier /= 0) THEN
            WRITE(UNIT = *, FMT = 6100) 'ijk_location'
            STOP
         END IF
!
! ...... Initializing
!
         ijk_location = 0
!
! ...... Transferring data from the <ElemCoordinates> file into arrays
!
         REWIND( UNIT = Elem_Coordinates_Unit)
!
         READ(   UNIT = Elem_Coordinates_Unit, FMT = 6901 ) ( ijk_location(n), n = 1, NumTotalElem )
!
         CLOSE ( UNIT = Elem_Coordinates_Unit, STATUS = 'DELETE' )
!
         Max_ijk_location = MAXVAL(ijk_location)
!
      END IF IF_Loop2
!
!
! ----------------
! ........ GJM - 9/19/2000
! ........ Begin addition
! ........   (a) for optimal numbering to minimize bandwidth
! ........   (b) to introduce different numbering/naming system to allow the use of the MINC facility for large grids
! ----------------
!
! >>>
! >>> The longest coordinate
! >>>
!
      CBi = Lref(3)
!
      DO_IMax: DO I=1,NLXYZ(3)
!
         IF( partial_processing .AND. ( I > LLoop_end ) ) EXIT DO_IMax
!
         CALL Coordinate_Processing( n = I, num = 3, type_of_grid = grid_type, DNd = D1d, CNi = CBi )
!
         IF( partial_processing .AND. ( I < LLoop_begin - 1 ) ) CYCLE DO_IMax
!
         IF( ILOOP == 1) THEN
            jbegin = 1
            jend   = NLXYZ(2)
         ELSE
            IF( jmin(i) >= jmax(i) ) THEN
               jbegin = 1
               jend   = NLXYZ(2)
            ELSE
               jbegin = jmin(i)
               jend   = jmax(i)
            END IF
         END IF
!
! >>>>>>
! >>>>>> The mid-sized coordinate
! >>>>>>
!
         CMi = Lref(2)
!
         DO_JMed: DO J = 1, NLXYZ(2)
!
            IF( partial_processing .AND. ( J > jend ) ) EXIT DO_JMed
!
            CALL Coordinate_Processing( n = J, num = 2, type_of_grid = grid_type, DNd = D2d, CNi = CMi )
!
            IF( partial_processing .AND. ( J < jbegin ) ) CYCLE DO_JMed
!
            IF( ILOOP == 1) THEN
               kbegin = 1
               kend   = NLXYZ(1)
            ELSE
               IF( kmin(i,j) >= kmax(i,j) ) THEN
                  kbegin = 1
                  kend   = NLXYZ(1)
               ELSE
                  kbegin = kmin(i,j)
                  kend   = kmax(i,j)
               END IF
            END IF
!
!
! >>>>>>
! >>>>>> The shortest coordinate
! >>>>>>
!
            CSm = Lref(1)
!
            DO_KMin: DO K= 1, NLXYZ(1)
!
               IF( partial_processing .AND. ( K > kend ) ) EXIT DO_KMin
!
               CALL Coordinate_Processing( n = K, num = 1, type_of_grid = grid_type, DNd = D3d, CNi = CSm )
!
               IF( partial_processing .AND. ( K < kbegin ) ) CYCLE DO_Kmin
!
! >>>>>>>>>>>>
! >>>>>>>>>>>> Determine the relationship between the optimized and the standard ordering schemes
! >>>>>>>>>>>>
!
               CALL Ordering_Relationship( coord_order = coord_order, type_of_grid = grid_type )
!
! <<<<<<<<<<<<
! <<<<<<<<<<<<
! <<<<<<<<<<<<
!
               subdomain_id = '   '
               completed_interpolation = .FALSE.
               available_InterpData    =  0.0d0   ! Array operation
!
               accepted_ijk  = .FALSE.
!
! ............ Determine the status of the element
!
               CALL Determine_Element_Status( Xloc = DDXYZ(1), Yloc = DDXYZ(2), Zloc = DDXYZ(3), &
     &                                        is_accepted  = accepted_ijk, &
     &                                        is_boundary  = is_boundary,  &
     &                                        is_inclusion = is_inclusion, &
     &                                        discontinue_smoothing = discontinue_smoothing, &
     &                                        medium_type  = medium, emissivity = emissivity )
!
               IF( accepted_ijk ) THEN
!
                  index = 'A'
                  IF( is_boundary  ) index = 'B'
                  IF( is_inclusion ) index = 'i'
                  IF( discontinue_smoothing ) index = 'D'
                  connectivity = medium(4:5)
!
               END IF
!
!>>>>>>>>>>>>>
!>>>>>>>>>>>>>
!>>>>>>>>>>>>> The element computations
!>>>>>>>>>>>>>
!>>>>>>>>>>>>>
!
               IF_LoopNum: IF( ILOOP == 1 ) THEN
!
                  IF_acc: IF( accepted_ijk ) THEN
!
                     IF( j < jmin(i) ) jmin(i) = j
                     IF( j > jmax(i) ) jmax(i) = j
!
                     IF( k < kmin(i,j) ) kmin(i,j) = k
                     IF( k > kmax(i,j) ) kmax(i,j) = k
!
                     IF( partial_processing ) THEN
                        IF( LLoop_begin == 1 .AND. I >= 1 .AND. I <= LLoop_end ) THEN
                           NumTotalElem  = NumTotalElem + 1
                           valid_element = .TRUE.
                        ELSE IF( I == LLoop_begin - 1 .AND. LLoop_begin > 1 ) THEN
                           NumElem_Bound1 = NumElem_Bound1 + 1
                           NumTotalElem   = NumTotalElem + 1
                           valid_element    = .FALSE.
                        ELSE IF( I >= LLoop_begin .AND. LLoop_begin > 1 .AND. I <= LLoop_end ) THEN
                           NumTotalElem  = NumTotalElem + 1
                           valid_element = .TRUE.
                        END IF
                     ELSE
                        NumTotalElem  = NumTotalElem + 1
                        valid_element = .TRUE.
                     END IF
!
                     IF_valid: IF( valid_element ) THEN
!
! ..................... First add to the element list
!
                        CALL Create_Element_List( i = I, j = J, k = K, medium_type = medium, volume = volume, elem_name = elem_name )
!
! ..................... Write element data in the VTK format for visualization
!
                        IF( VTK_output ) CALL Write_VTK_Data_File( itrue = itrue, jtrue = jtrue, ktrue = ktrue, elem_name = elem_name )
!
                     END IF IF_valid
!
! .................. Store element locations in a temporary file
!
                     ijk = itrue + (jtrue - 1) * NX + ( ktrue - 1) * NXNY
!
                     m = MOD(NumTotalElem,40)
                     IF( m > 0 .AND. m <= 39 ) THEN
                        WRITE( UNIT = Elem_Coordinates_Unit, FMT = 6900, ADVANCE = 'NO' ) ijk
                     ELSE
                        WRITE( UNIT = Elem_Coordinates_Unit, FMT = 6900 ) ijk
                     END IF
!
                  END IF IF_acc
!
                  CYCLE DO_KMin
!
!>>>>>>>>>>>>>
!>>>>>>>>>>>>>
!>>>>>>>>>>>>> The connection computations
!>>>>>>>>>>>>>
!>>>>>>>>>>>>>
!
               ELSE IF(ILOOP == 2) THEN
!
                  accepted_ip_true = .FALSE.
                  accepted_jp_true = .FALSE.
                  accepted_kp_true = .FALSE.
!
                  accepted_ipjk = .FALSE.
                  accepted_ijpk = .FALSE.
                  accepted_ijkp = .FALSE.
!
                  IF_accepted_ijk: IF( accepted_ijk ) THEN
!
                     coeff = 1.0d0  ! ... Array operation
!
                     IF_cont0: IF( ensure_grid_continuity ) THEN
!
                        Z_00 = Z(ktrue) + 5.0d-1 * DZ(ktrue)
                        X_00 = X(itrue)
                        Y_00 = Y(jtrue)
!
                        IF( itrue < NX ) THEN
                           kx_offset = RefSurfLoc( itrue+1, jtrue ) - RefSurfLoc( itrue, jtrue )
                        ELSE
                           kx_offset = 0
                        END IF
!
                        IF( jtrue < NY ) THEN
                           ky_offset = RefSurfLoc( itrue, jtrue+1 ) - RefSurfLoc( itrue, jtrue )
                        ELSE
                           ky_offset = 0
                        END IF
!
                        IF( itrue == NX ) THEN
                           k_xp = ktrue
                           X_xp = X(itrue)
                        ELSE
                           k_xp = ktrue + kx_offset
                           IF( k_xp <= 0 .OR. k_xp > NZ ) THEN
                              k_xp = ktrue
                              kx_offset = 0
                           END IF
                           X_xp = X(itrue+1)
                        END IF
!
                        IF( jtrue == NY ) THEN
                           k_yp = ktrue
                           Y_yp = Y(jtrue)
                        ELSE
                           k_yp = ktrue + ky_offset
                           IF( k_yp <= 0 .OR. k_yp > NZ ) THEN
                              k_yp = ktrue
                              ky_offset = 0
                           END IF
                           Y_yp = Y(jtrue+1)
                        END IF
!
                        IF_smooth: IF( .NOT. smooth_continuous_surface ) THEN
                           IF( ABS(kx_offset) <= max_offset ) THEN
                              kx_offset = 0
                              k_xp      = ktrue
                           END IF
                           IF( ABS(ky_offset) <= max_offset ) THEN
                              ky_offset = 0
                              k_yp      = ktrue
                           END IF
                        END IF IF_smooth
!
                        Z_xp = Z(k_xp) + 5.0d-1 * DZ(k_xp)
                        Z_yp = Z(k_yp) + 5.0d-1 * DZ(k_yp)
!
                        IF_Xp: IF( ABS(Z_xp - Z_00) > 1.0d-8 ) THEN
                           phi_x = ATAN ( (X_xp - X_00) / (Z_xp - Z_00) )
                           IF( phi_x > 0.0d0 ) phi_x = PI - phi_x
                        ELSE
                           phi_x = 5.0d-1 * PI
                        END IF IF_Xp
!
                        IF_Yp: IF( ABS(Z_yp - Z_00) > 1.0d-8 ) THEN
                           IF(grid_type == 2) THEN
                              phi_y = ATAN ( 2.0d0 * X_00 * SIN( 5.0d-1*(Y_yp -Y_00) ) / (Z_yp - Z_00) )
                           ELSE
                              phi_y = ATAN ( (Y_yp - Y_00) / (Z_yp - Z_00) )
                           END IF
                           IF( phi_y > 0.0d0 ) phi_y = PI - phi_y
                        ELSE
                           phi_y = 5.0d-1 * PI
                        END IF IF_Yp
!
                     ELSE
!
                        kx_offset = 0
                        ky_offset = 0
!
                        IF( itrue < NX ) THEN
                           IF( grid_type == 2 ) THEN
                              X_xp = X(itrue+1) * COS( Y(jtrue) )
                           ELSE
                              X_xp = DDXYZ(1) + 5.0d-1 * ( DX(itrue) + DX(itrue+1) ) ! X(itrue+1)
                           END IF
                        END IF
!
                        IF( jtrue < NY ) THEN
                           IF( grid_type == 2 ) THEN
                              Y_yp = X(itrue) * SIN( Y(jtrue + 1) )
                           ELSE
                              Y_yp = DDXYZ(2) + 5.0d-1 * ( DY(jtrue) + DY(jtrue+1) ) ! Y(jtrue+1)
                           END IF
                        END IF
!
                        beta_x = BET
                        beta_y = 0.0d0
                        beta_z = BETA
!
                     END IF IF_cont0
!
!
! .................. The TRUE (i+1,j,k)-element
!
!
                     IF( NX == 1 .OR. itrue == NX ) GO TO 1000
!
                     IF_cont_ip: IF( ensure_grid_continuity ) THEN
!
! ..................... For transformed coordinates
!
                        IF (grid_type == 2 ) THEN
                           angle = Y_00
                           Y_00  = X_xp * SIN(angle)
                           X_xp  = X_xp * COS(angle)
                        END IF
!
                        CALL Determine_Element_Status( Xloc = X_xp, Yloc = Y_00, Zloc = Z_xp, &
     &                                                 is_accepted  = accepted_ip_true,       &
     &                                                 is_boundary  = is_boundary,            &
     &                                                 is_inclusion = is_inclusion,           &
     &                                                 discontinue_smoothing = discontinue_smoothing, &
     &                                                 medium_type  = medium, emissivity = emissivity_p )
!
                        IF_ip_true: IF( accepted_ip_true ) THEN
!
                           IF_index_1a: IF( index /= 'D' .AND. (.NOT. discontinue_smoothing) ) THEN
!
                              IF( index == 'B' .AND. is_boundary ) THEN
                                 accepted_ip_true = .FALSE.
                              ELSE
                                 coeff(1) = 1.0d0 / ABS( SIN(phi_x) )
                                 ipjk = itrue + 1 + (jtrue - 1) * NX + ( k_xp - 1 ) * NXNY
                                 IF( ipjk > Max_ijk_location ) THEN
                                    NumElem_Bound2 = NumElem_Bound2 + 1
                                    NumTotalElem   = NumTotalElem + 1
                                    ijk_location(NumTotalElem) = ipjk
                                 END IF
                              END IF
!
                           ELSE
!
                              CALL Determine_Element_Status( Xloc = X_xp, Yloc = Y_00, Zloc = Z_00, &
     &                                                       is_accepted  = accepted_ip_true,       &
     &                                                       is_boundary  = is_boundary,            &
     &                                                       is_inclusion = is_inclusion,           &
     &                                                       discontinue_smoothing = discontinue_smoothing, &
     &                                                       medium_type  = medium, emissivity = emissivity_p )
!
                              IF( accepted_ip_true ) THEN
                                 phi_x     = 5.0d-1 * pi
                                 coeff(1)  = 1.0d0
                                 kx_offset = 0
                              END IF
!
                           END IF IF_index_1a
!
                        ELSE
!
                           IF_index_1b: IF( index == 'D' .OR. (discontinue_smoothing) ) THEN
!
                              CALL Determine_Element_Status( Xloc = X_xp, Yloc = Y_00, Zloc = Z_00, &
     &                                                       is_accepted  = accepted_ip_true,       &
     &                                                       is_boundary  = is_boundary,            &
     &                                                       is_inclusion = is_inclusion,           &
     &                                                       discontinue_smoothing = discontinue_smoothing, &
     &                                                       medium_type  = medium, emissivity = emissivity_p )
!
                              IF( accepted_ip_true ) THEN
                                 phi_x     = 5.0d-1 * pi
                                 coeff(1)  = 1.0d0
                                 kx_offset = 0
                              END IF
!
                           END IF IF_index_1b
!
                        END IF IF_ip_true
!
                     ELSE
!
                        CALL Determine_Element_Status( Xloc = X_xp, Yloc = DDXYZ(2), Zloc = DDXYZ(3), &
     &                                                 is_accepted  = accepted_ip_true,   &
     &                                                 is_boundary  = is_boundary,        &
     &                                                 is_inclusion = is_inclusion,       &
     &                                                 discontinue_smoothing = discontinue_smoothing, &
     &                                                 medium_type  = medium, emissivity = emissivity_p )
!
                        IF( index == 'B' .AND. is_boundary ) accepted_ip_true = .FALSE.
!
                     END IF IF_cont_ip
!
                     IF( (connectivity == medium(4:5)) .AND. (connectivity == '-x' .OR. connectivity == '-X') ) accepted_ip_true = .FALSE.
!
 1000                SELECT CASE(coord_order)
                     CASE(1,2)
                        accepted_ipjk = accepted_ip_true
                     CASE(3,5)
                        accepted_ijpk = accepted_ip_true
                     Case(4,6)
                        accepted_ijkp = accepted_ip_true
                     END SELECT
!
!
! .................. The TRUE (i,j+1,k)-element
!
!
                     IF( NY == 1 .OR. jtrue == NY ) GO TO 1500
!
                     IF_cont_jp: IF( ensure_grid_continuity ) THEN
!
! ..................... For transformed coordinates
!
                        IF (grid_type == 2 ) THEN
                           angle = Y_yp
                           Y_yp  = X(itrue) * SIN(angle)
                           X_00  = X(itrue) * COS(angle)
                        END IF
!
                        CALL Determine_Element_Status( Xloc = X_00, Yloc = Y_yp, Zloc = Z_yp, &
     &                                                 is_accepted  = accepted_jp_true,       &
     &                                                 is_boundary  = is_boundary,            &
     &                                                 is_inclusion = is_inclusion,           &
     &                                                 discontinue_smoothing = discontinue_smoothing, &
     &                                                 medium_type  = medium, emissivity = emissivity_p )
!
                        IF_jp_true: IF( accepted_jp_true ) THEN
!
                           IF_index_2a: IF( index /= 'D' .AND. (.NOT. discontinue_smoothing) ) THEN
!
                              IF( index == 'B' .AND. is_boundary ) THEN
                                 accepted_jp_true = .FALSE.
                              ELSE
                                 coeff(2) = 1.0d0 / ABS( SIN(phi_y) )
                                 ijpk = itrue + jtrue * NX + ( k_yp - 1 ) * NXNY
                                 IF( ijpk > Max_ijk_location ) THEN
                                    NumElem_Bound2 = NumElem_Bound2 + 1
                                    NumTotalElem   = NumTotalElem + 1
                                    ijk_location(NumTotalElem) = ijpk
                                 END IF
                              END IF
!
                           ELSE
!
                              CALL Determine_Element_Status( Xloc = X_00, Yloc = Y_yp, Zloc = Z_00, &
     &                                                       is_accepted  = accepted_jp_true,       &
     &                                                       is_boundary  = is_boundary,            &
     &                                                       is_inclusion = is_inclusion,           &
     &                                                       discontinue_smoothing = discontinue_smoothing, &
     &                                                       medium_type  = medium, emissivity = emissivity_p )
!
                              IF( accepted_jp_true ) THEN
                                 phi_y     = 5.0d-1 * pi
                                 coeff(2)  = 1.0d0
                                 ky_offset = 0
                              END IF
!
                           END IF IF_index_2a
!
                        ELSE
!
                           IF_index_2b: IF( index == 'D' .OR. (discontinue_smoothing) ) THEN
!
                              CALL Determine_Element_Status( Xloc = X_00, Yloc = Y_yp, Zloc = Z_00, &
     &                                                       is_accepted  = accepted_jp_true,       &
     &                                                       is_boundary  = is_boundary,            &
     &                                                       is_inclusion = is_inclusion,           &
     &                                                       discontinue_smoothing = discontinue_smoothing, &
     &                                                       medium_type  = medium, emissivity = emissivity_p )
!
                              IF( accepted_jp_true ) THEN
                                 phi_y     = 5.0d-1 * pi
                                 coeff(2)  = 1.0d0
                                 ky_offset = 0
                              END IF
!
                           END IF IF_index_2b
!
                        END IF IF_jp_true
!
                     ELSE
!
                        CALL Determine_Element_Status( Xloc = DDXYZ(1), Yloc = Y_yp, Zloc = DDXYZ(3), &
     &                                                 is_accepted  = accepted_jp_true,   &
     &                                                 is_boundary  = is_boundary,        &
     &                                                 is_inclusion = is_inclusion,       &
     &                                                 discontinue_smoothing = discontinue_smoothing, &
     &                                                 medium_type  = medium, emissivity = emissivity_p )
!
                        IF( index == 'B' .AND. is_boundary ) accepted_jp_true = .FALSE.
!
                     END IF IF_cont_jp
!
                     IF( (connectivity == medium(4:5)) .AND. (connectivity == '-y' .OR. connectivity == '-Y') ) accepted_jp_true = .FALSE.
!
 1500                SELECT CASE(coord_order)
                     CASE(3,4)
                        accepted_ipjk = accepted_jp_true
                     CASE(1,6)
                        accepted_ijpk = accepted_jp_true
                     Case(2,5)
                        accepted_ijkp = accepted_jp_true
                     END SELECT
!
!
! .................. The TRUE (i,j,k+1)-element
!
!
                     IF( NZ == 1 .OR. ktrue == NZ ) GO TO 2000
!
                     CALL Determine_Element_Status( Xloc = DDXYZ(1), Yloc = DDXYZ(2), Zloc = DDXYZ(3) - 5.0d-1 * ( DZ(ktrue) + DZ(ktrue+1) ), &
     &                                              is_accepted  = accepted_kp_true, &
     &                                              is_boundary  = is_boundary,      &
     &                                              is_inclusion = is_inclusion,     &
     &                                              discontinue_smoothing = discontinue_smoothing, &
     &                                              medium_type  = medium, emissivity = emissivity_p )
!
                     IF( index == 'B' .AND. is_boundary ) accepted_kp_true = .FALSE.
!
                     IF( (connectivity == medium(4:5)) .AND. (connectivity == '-z' .OR. connectivity == '-Z') ) accepted_kp_true = .FALSE.
!
                     SELECT CASE(coord_order)
                     CASE(5,6)
                        accepted_ipjk = accepted_kp_true
                     CASE(2,4)
                        accepted_ijpk = accepted_kp_true
                     Case(1,3)
                        accepted_ijkp = accepted_kp_true
                     END SELECT
!
                  ELSE
!
                     CYCLE DO_KMin
!
                  END IF IF_accepted_ijk
!
               END IF IF_LoopNum
!
 2000          IF( accepted_ijk .AND. (.NOT. accepted_ip_true) .AND. (.NOT. accepted_jp_true) .AND. (.NOT. accepted_kp_true) ) THEN
                  CYCLE DO_KMin
               END IF
!
! >>>>>>>>>>>>
! >>>>>>>>>>>>
! >>>>>>>>>>>> Come here for connection processing and printing
! >>>>>>>>>>>>
! >>>>>>>>>>>>
!
               CALL Create_Connection_List( i = I, j = J, k = K, accepted_ijk = accepted_ijk, &
     &                                      accepted_ijk_ext = (/ accepted_ipjk, accepted_ijpk, accepted_ijkp /), &
     &                                      outside_I_Range_limits = outside_I_Range_limits )
!
               IF( outside_I_Range_limits .AND. (.NOT. ensure_grid_continuity) ) CYCLE DO_IMax
!
            END DO DO_KMin
!
         END DO DO_JMed
!
      END DO DO_IMax
!
!
! ----------------
! ........ GJM - 9/19/2000
! ........ End addition
! ........   (a) for optimal numbering to minimize bandwidth
! ........   (b) to introduce diffenet numbering/naming system to allow
! ........       the use of the MINC facility for large grids
! ----------------
!
!
!
      IF_Format6: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
!
         IF(ILOOP == 1) THEN
            WRITE(output_file_unit(1), FMT = 5551)
!
            SELECT CASE (level_of_grid_generation)
            CASE('P', 'p')
               WRITE(output_file_unit(2), FMT = 5552)
            END SELECT
         END IF
!
         IF(ILOOP == 2) WRITE(output_file_unit(1), FMT = 5552)
!
      ELSE
!
         IF( VTK_output .AND. ILOOP == 1 ) THEN
            WRITE( UNIT = Corner_Unit, FMT = 6905 ) NumElem, ijk_min, ijk_max
            IF( renumbered_inactive_elements ) WRITE( UNIT = InactElement_Unit, FMT = 6905 ) NumInactElem, ijk_min_IE, ijk_max_IE
         END IF
!
! >>>>>>
! >>>>>> Modifying the ELEMENT data block for old-style inactive element
! >>>>>>
!
         IF_Style: IF( old_style_inactive_elements ) THEN
!
            IF_contin: IF( ensure_grid_continuity ) THEN
!
               IF_ElemAdj: IF( ILOOP == 2 ) THEN
!
! ............... Resize the element volumes (for forced areal grid continuity )
!
                  IF( ensure_grid_continuity ) CALL Adjust_Element_Volumes
                  GO TO 2500
!
               ELSE IF( (ILOOP == 1) .AND. (level_of_grid_generation == 'F' .OR. level_of_grid_generation == 'f') ) THEN
!
                  CONTINUE
!
               ELSE
!
                  GO TO 3000
!
               END IF IF_ElemAdj
!
! ............ Go to the top of the <INACT_ELEM> file
!
!              ier = FSEEK( InactElem_Unit,  0,  0 )
               CALL FSEEK( InactElem_Unit,  0,  0 )
!
! ............ Stop if there is a problem
!
               IF(ier /= 0) THEN
                  WRITE(UNIT = *, FMT = 6800)
                  STOP
               END IF
!
               DO i = 1, NumInactElem + 1
                  READ ( UNIT = InactElem_Unit,      FMT = '(A)' ) line
                  WRITE( UNIT = output_file_unit(1), FMT = '(A)' ) line
               END DO
!
! ............ Close and delete the <INACT_ELEM> file
!
               CLOSE ( UNIT = InactElem_Unit, STATUS = 'DELETE' )
!
            END IF IF_contin
!
         ELSE
!
            IF( (ILOOP == 2) .AND. ensure_grid_continuity ) CALL Adjust_Element_Volumes
!
         END IF IF_Style
!
         IF( ILOOP == 1 ) THEN
            INQUIRE( FILE = TempElem_file_name, EXIST = exists )
            IF( exists .AND. ( TempElem_file_unit == output_file_unit(1) ) ) THEN
               WRITE( output_file_unit(1), FMT = 22 )
            ELSE
               INQUIRE( FILE = elem_file_name, EXIST = exists )
               IF( exists .AND. ( TempElem_file_unit /= output_file_unit(1) ) ) WRITE( output_file_unit(1), FMT = 22 )
            END IF
         END IF
!
 2500    IF( ILOOP == 2 ) THEN
            INQUIRE( FILE = TempConx_file_name, EXIST = exists )
            IF( exists .AND. ( TempConx_file_unit == output_file_unit(2) ) ) THEN
               WRITE( output_file_unit(2), FMT = 22 )
            ELSE
               INQUIRE( FILE = conx_file_name, EXIST = exists )
               IF( exists .AND. ( TempConx_file_unit /= output_file_unit(2) ) ) WRITE( output_file_unit(2), FMT = 22 )
            END IF
         END IF
!
         INQUIRE( FILE = 'MESH', EXIST = exists )
         IF( exists .AND. ( level_of_grid_generation == 'F' .OR. level_of_grid_generation == 'f' ) .AND. (.NOT. ensure_grid_continuity) ) THEN
            WRITE( MESH_unit, FMT = 22 )
            WRITE( MESH_unit, FMT = 22 )
         END IF
!
      END IF IF_Format6
!
! >>>
! >>>
! >>>
!
 3000 IF( ILOOP == 1 ) GO TO 30
!
! >>>
! >>>
! >>>
!
      IF( Max_NumElem >= 200000 .OR. Flag_ExclZones .OR. Flag_InclZones) RETURN
!
      REWIND (UNIT = output_file_unit(1))
      READ(output_file_unit(1),50) DOM
!
      L=0
!
! ... Begin modification for 8-character elements
!
      IF(ElemNameLength == 8) THEN
         DO I=1,NX
            DO J=1,NY
               DO K=1,NZ
                  L = L+1
                  READ(output_file_unit(1),5050) name(L)
               END DO
            END DO
         END DO
      ELSE
         DO I=1,NX
            DO J=1,NY
               DO K=1,NZ
                  L = L+1
                  READ(output_file_unit(1),50) name(L)(1:5)
               END DO
            END DO
         END DO
      END IF
!
 5050 FORMAT(A8)
   50 FORMAT(A5)
!
! ... End modification for 8-character elements
!
      CALL PCAR(NLXYZ(3),NLXYZ(2),NLXYZ(1),NX,NY,NZ,coord_order)
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
!6000 FORMAT(/6X,'GXYZ     1.0      18 MARCH     1991',6X,'GENERATE 1, 2, OR 3-D CARTESIAN MESH')
!6000 FORMAT(/,  'GXYZ     2.0      18 September 2000',6X,'GENERATE 1, 2, OR 3-D CARTESIAN MESH')
 6000 FORMAT(/,'GXYZ 2.0 ............................... 08 July      2016',6X,'GENERATE 1, 2, OR 3-D CARTESIAN OR HALF-CYLINDRICAL MESH')
!
 5010 FORMAT(A8, 7X,A5,2ES10.4E1,10X,3ES10.4E1,1X,A3)
 5510 FORMAT(A8, 7X,A5,1ES10.4E1,20X,3ES10.4E1,1X,A3)
   10 FORMAT(A5,10X,A5,2ES10.4E1,10X,3ES10.4E1,1X,A3)
 5610 FORMAT(A5,10X,A5,1ES10.4E1,20X,3ES10.4E1,1X,A3)
 5011 FORMAT(2(A8),13X,A1,4ES10.4E1)
   11 FORMAT(2(A5),19X,A1,4ES10.4E1)
!
!
!
   20 FORMAT('ELEME')
   21 FORMAT('CONNE')
   22 FORMAT(120(' '))
 5020 FORMAT('ELEMEext2')
!
   40 FORMAT(/,' GENERATE CARTESIAN MESH WITH Nx*Ny*Nz = ',   &
     &I5,' *',I5,' *',I5,'  =  ',I9,'  GRID BLOCKS'/1X,131('*')//   &
     &'   THE X-AXIS IS HORIZONTAL; THE Y-AXIS IS ROTATED BY ',ES11.4,' DEGREES AGAINST THE HORIZONTAL'//   &
     &'   THE GRID INCREMENTS ARE AS FOLLOWS')
!
  440 FORMAT(/,' GENERATE 3D CYLINDRICAL MESH WITH Nr*NTh*Nz = ',   &
     &I5,' *',I5,' *',I5,'  =  ',I9,'  GRID BLOCKS'/1X,131('*')//   &
     &'   THE R- AND Th-AXES ARE HORIZONTAL; USED FOR INCLINED SINGLE_WELL PROBLEMS, DESCRIBES AN ANGLE Theta = Pi'//   &
     &'   THE GRID INCREMENTS ARE AS FOLLOWS')
!
   41 FORMAT(//'   DELTA-X  (',I6,' INCREMENTS)'//13X,8(A1,11X),A1,10X,'10'//(6X,10(1X,ES11.4)))
   42 FORMAT(//'   DELTA-Y  (',I6,' INCREMENTS)'//13X,8(A1,11X),A1,10X,'10'//(6X,10(1X,ES11.4)))
   43 FORMAT(//'   DELTA-Z  (',I6,' INCREMENTS)'//13X,8(A1,11X),A1,10X,'10'//(6X,10(1X,ES11.4)))
   44 FORMAT(/' ',131('*'))
!
   46 FORMAT(//'   DELTA-R  (',I6,' INCREMENTS)'//13X,8(A1,11X),A1,10X,'10'//(6X,10(1X,ES11.4)))
   48 FORMAT(//'   DELTA-Th (',I6,' INCREMENTS)'//13X,8(A1,11X),A1,10X,'10'//(6X,10(1X,ES11.4)))
!
  441 FORMAT(//'   X  (',I6,' COORDINATES)'//13X,8(A1,13X),A1,12X,'10'//(6X,10(1X,ES13.6)))
  442 FORMAT(//'   Y  (',I6,' COORDINATES)'//13X,8(A1,13X),A1,12X,'10'//(6X,10(1X,ES13.6)))
  443 FORMAT(//'   Z  (',I6,' COORDINATES)'//13X,8(A1,13X),A1,12X,'10'//(6X,10(1X,ES13.6)))
  444 FORMAT(//'   Z  LAYERS - Beginning from the top, moving downward:')
  445 FORMAT(T5,2(1X,ES15.8),3x,'! Layer ',i5.5)
!
  446 FORMAT(//'   R  (',I6,' COORDINATES)',//,13X,8(A1,13X),A1,12X,'10'//(6X,10(1X,ES13.6)))
  448 FORMAT(//'   Th (',I6,' COORDINATES) - in Radians',//,13X,8(A1,13X),A1,12X,'10'//(6X,10(1X,ES13.6)))
!
 5551 FORMAT(T1,'<<<END of ELEMENTS data block',/,'     ')
 5552 FORMAT(T1,'<<<END of CONNECTIONS data block',/,'     ')
!
 6100 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <Read_Tabular_Data> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6500 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <GXYZ>: There is a problem reading the discretization axis descriptor'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6501 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <GXYZ>: ',A,'-discretization is specified in a ',A,' system',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6505 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <GXYZ>: The discretization axis descriptor = "',A,'" is not an available option',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6602 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <GXYZ>: There is a problem opening the file ',A,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6800 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <GXYZ>: There is a problem going to the top of the <INACT_ELEM> file',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6850 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6900 FORMAT(i10,1x)
 6901 FORMAT( 40(i10,1x) )
!
 6902 FORMAT(ES14.7,1X)
 6903 FORMAT( 30(ES14.7,1X) )
!
 6905 FORMAT( '>>>',3(2x,i10) )
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of GXYZ
!
!
      RETURN
!
!
!***********************************************************************
!*                                                                     *
!*                        INTERNAL PROCEDURES                          *
!*                                                                     *
!***********************************************************************
!
!
      CONTAINS
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Discretize_Axis( n, direction, DL, NL )
!
!
      IMPLICIT NONE
!
! -------
! ... Double precision arrays
! -------
!
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: radii
!
! -------
! ... Double precision arrays
! -------
!
      REAL(KIND = 8), INTENT(INOUT), DIMENSION(:) :: DL
!
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: Lxyz
!
! -------
! ... Double precision variables
! -------
!
      REAL(KIND = 8) :: Delta_L, L_max, DL_Log, L_last
      REAL(KIND = 8) :: Delta_X, X_max, DX_Log
      REAL(KIND = 8) :: Delta_Y, Y_max, DY_Log
      REAL(KIND = 8) :: Delta_Z, Z_max, DZ_Log
!
      REAL(KIND = 8) :: Delta_R, R_max, DR_Log
      REAL(KIND = 8) :: Delta_Th, Th_max, Rdev
!
      REAL(KIND = 16) :: IncrFact
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN)  :: n
      INTEGER, INTENT(OUT) :: NL
!
      INTEGER :: num_DL_subsets, num_DX_subsets, num_DY_subsets, num_DZ_subsets, num_DR_subsets, num_DTh_subsets
      INTEGER :: number_of_DLs,  number_of_DXs,  number_of_DYs,  number_of_DZs,  number_of_DRs,  number_of_DThs, number_of_radii, error
!
      INTEGER :: i, n1, NLMax, NL1, NLo, NLb, mx, imax
!
! -------
! ... Character variables
! -------
!
      CHARACTER(LEN = 1), INTENT(IN) :: direction
!
      CHARACTER(LEN = 2)   :: direction2
      CHARACTER(LEN = 15)  :: option
      CHARACTER(LEN = 120) :: read_format
!
! -------
! ... Logical variables
! -------
!
      LOGICAL :: First_call = .TRUE.
!
! -------
! ... Saving variables
! -------
!
      SAVE First_call
!
! -------
! ... Namelists
! -------
!
      NAMELIST/ DX_Subsets / num_DX_subsets
      NAMELIST/ DY_Subsets / num_DY_subsets
      NAMELIST/ DZ_Subsets / num_DZ_subsets
!
      NAMELIST/ DR_Subsets  / num_DR_subsets
      NAMELIST/ DTh_Subsets / num_DTh_subsets
!
      NAMELIST/ DX_Data / number_of_DXs, Delta_X, read_format, X_max, DX_Log, option
      NAMELIST/ DY_Data / number_of_DYs, Delta_Y, read_format, Y_max, DY_Log, option
      NAMELIST/ DZ_Data / number_of_DZs, Delta_Z, read_format, Z_max, DZ_Log, option
!
      NAMELIST/ DR_Data  / number_of_DRs,  Delta_R,  read_format, R_max,  option, DR_Log, number_of_radii
      NAMELIST/ DTh_Data / number_of_DThs, Delta_Th, read_format, Th_max, option
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Discretize_Axis>
!
!
      IF(First_call) THEN
!
         WRITE( UNIT = *, FMT = 6000 )
         First_call = .FALSE.
!
         IF ( coordinates(1:1) == 'T' .OR. coordinates(1:1) == 't' ) THEN
            ALLOCATE( radii(0:NXMax) )
            radii(0) = X_ref
         END IF
!
      END IF
!
! ... Initialization
!
      mx    = 1
      NL    = 0
      NLMax = SIZE(DL)
!
      direction2(1:1) = direction
!
! ... Posible directions
!
      SELECT CASE(direction)
!
      CASE('X','x')
!
         IF( n /= 1 ) THEN
            WRITE(UNIT = *, FMT = 6501) discretization_descriptor, ':::>>>X-DISCRETIZATION'
            STOP
         END IF
!
         READ (UNIT = *, NML = DX_Subsets, IOSTAT = error )
         num_DL_subsets = num_DX_subsets
!
      CASE('R','r')
!
         IF( n /= 1 ) THEN
            WRITE(UNIT = *, FMT = 6501) discretization_descriptor, ':::>>>R-DISCRETIZATION'
            STOP
         END IF
!
         READ (UNIT = *, NML = DR_Subsets, IOSTAT = error )
         num_DL_subsets = num_DR_subsets
!
      CASE('Y','y')
!
         IF( n /= 2 ) THEN
            WRITE(UNIT = *, FMT = 6501) discretization_descriptor, ':::>>>Y-DISCRETIZATION'
            STOP
         END IF
!
         READ (UNIT = *, NML = DY_Subsets, IOSTAT = error )
         num_DL_subsets = num_DY_subsets
!
      CASE('T','t')
!
         IF( n /= 2 ) THEN
            WRITE(UNIT = *, FMT = 6501) discretization_descriptor, ':::>>>Th-DISCRETIZATION'
            STOP
         END IF
!
         mx = 2
         direction2(2:2) = 'h'
         READ (UNIT = *, NML = DTh_Subsets, IOSTAT = error )
         num_DL_subsets = num_DTh_subsets
!
      CASE('Z','z')
!
         IF( n /= 3 ) THEN
            WRITE(UNIT = *, FMT = 6501) discretization_descriptor, ':::>>>Z-DISCRETIZATION'
            STOP
         END IF
!
         READ (UNIT = *, NML = DZ_Subsets, IOSTAT = error )
         num_DL_subsets = num_DZ_subsets
!
      END SELECT
!
! ... Stop if there is a problem reading the namelists
!
      IF(error /= 0) THEN
         WRITE(UNIT = *, FMT = 6601) '<D'//direction2(1:mx)//'_Subsets>'
         STOP
      END IF
!
! ... Allocate distance array
!
      ALLOCATE( Lxyz(NLMax), STAT = error )
!
      IF(error /= 0) THEN
         WRITE(UNIT = *, FMT = 6850) direction//'-associated Lxyz', 'Discretize_Axis'
         STOP
      END IF
!
! ... Array initialization
!
      DL = 0.0d0
!
! ... Read the namelists describing axis discretization
!
      Rdev = 0.0d0
!
      DO_LDiscr: DO n1=1,num_DL_subsets + 1
!
! ...... Inquiring about the end of the dataset
!
         IF( n1 == num_DL_subsets + 1 ) THEN
            READ( UNIT = *, FMT = '(a6)' ) first_part
            IF( first_part /= ':::<<<' ) THEN
               WRITE(UNIT = *, FMT = 6515) discretization_descriptor
               STOP
            ELSE
               DEALLOCATE ( Lxyz )
               RETURN
            END IF
         END IF
!
! ...... Initialization of data for the logarithmic DX distribution
!
         mx = 1
         read_format = '  '
         option      = '  '
         direction2(1:1) = direction
!
         SELECT CASE(direction)
         CASE('X','x')
!
            number_of_DXs  = 0
            X_max          = 0.0d0
            DX_Log         = 0.0d0
            Delta_X        = 0.0d0
!
            READ (UNIT = *, NML = DX_Data, IOSTAT = error )
!
            IF( number_of_DXs < 0 ) THEN
               WRITE(UNIT = *, FMT = 6603) '<D'//direction2(1:mx)//'_Data>', n1, 'number_of_DXs', number_of_DXs
            END IF
!
         CASE('R','r')
!
            number_of_DRs  = 0
            R_max          = 0.0d0
            DR_Log         = 0.0d0
            Delta_R        = 0.0d0
 !
            READ (UNIT = *, NML = DR_Data, IOSTAT = error )
!
            IF( number_of_DRs < 0 ) THEN
               WRITE(UNIT = *, FMT = 6603) '<D'//direction2(1:mx)//'_Data>', n1, 'number_of_DRs', number_of_DRs
            END IF
!
         CASE('Y','y')
!
            number_of_DYs  = 0
            Y_max          = 0.0d0
            DY_Log         = 0.0d0
            Delta_Y        = 0.0d0
!
            READ (UNIT = *, NML = DY_Data, IOSTAT = error )
!
            IF( number_of_DYs < 0 ) THEN
               WRITE(UNIT = *, FMT = 6603) '<D'//direction2(1:mx)//'_Data>', n1, 'number_of_DYs', number_of_DYs
            END IF
!
         CASE('T','t')
!
            mx = 2
            direction2(2:2) = 'h'
            number_of_DThs  = 0
            Th_max          = 0.0d0
            Delta_Th        = 0.0d0
!
            READ (UNIT = *, NML = DTh_Data, IOSTAT = error )
!
            IF( number_of_DThs < 0 ) THEN
               WRITE(UNIT = *, FMT = 6603) '<D'//direction2(1:mx)//'_Data>', n1, 'number_of_DThs', number_of_DThs
            END IF
!
! ......... Convert degrees to radiants
!
            Delta_Th = Delta_Th * PI / 1.8d2
!
         CASE('Z','z')
!
            number_of_DZs  = 0
            Z_max          = 0.0d0
            DZ_Log         = 0.0d0
            Delta_Z        = 0.0d0
!
            READ (UNIT = *, NML = DZ_Data, IOSTAT = error )
!
            IF( number_of_DZs < 0 ) THEN
               WRITE(UNIT = *, FMT = 6603) '<D'//direction2(1:mx)//'_Data>', n1, 'number_of_DZs', number_of_DZs
            END IF
!
         END SELECT
!
         IF(error /= 0) THEN
            WRITE(UNIT = *, FMT = 6602) '<D'//direction2(1:mx)//'_Data>',n1
            STOP
         END IF
!
         SELECT CASE(option(1:1))
         CASE('E','e','U','u','V','v','L','l','R','r')
            CONTINUE
         CASE DEFAULT
            WRITE(UNIT = *, FMT = 6605) '<D'//direction2(1:mx)//'_Data>', n1, option
            STOP
         END SELECT
!
! ...... Assignment to working variables/arrays
!
         SELECT CASE(direction)
         CASE('X','x')
!
            number_of_DLs  = number_of_DXs
            L_max          = X_max
            DL_Log         = DX_Log
            Delta_L        = Delta_X
!
            DL = DX   ! ... Array operation
!
         CASE('R','r')
!
            number_of_DLs  = number_of_DRs
            L_max          = R_max
            DL_Log         = DR_Log
            Delta_L        = Delta_R
!
            DL = DX   ! ... Array operation
!
         CASE('Y','y')
!
            number_of_DLs  = number_of_DYs
            L_max          = Y_max
            DL_Log         = DY_Log
            Delta_L        = Delta_Y
!
            DL = DY   ! ... Array operation
!
         CASE('T','t')
!
            number_of_DLs  = number_of_DThs
            L_max          = Th_max
            Delta_L        = Delta_Th
!
            DL = DY   ! ... Array operation
!
         CASE('Z','z')
!
            number_of_DLs  = number_of_DZs
            L_max          = Z_max
            DL_Log         = DZ_Log
            Delta_L        = Delta_Z
!
            DL = DZ   ! ... Array operation
!
         END SELECT
!
         NLo = NL
         NL1 = NL + 1
         NL  = NL + number_of_DLs
!
         IF( NL > NLMax ) THEN
            WRITE(UNIT = *, FMT = 6520) 'N'//direction2(1:mx), NL+number_of_DLs, 'MaxNum_'//direction//'_Subdivisions', NLMax
            STOP
         END IF
!
! >>>>>>
! >>>>>> The case of logarithmic DL distribution
! >>>>>>
!
         IF_Log: IF( option(1:1) == 'L' .OR. option(1:1) == 'l' ) THEN
!
            IF_DL0: IF( ABS(DL_Log) <= 1.0e-7_16 ) THEN
!
                WRITE(*,*)  'Route 1', DL_Log
               IF( NLo > 1 ) THEN
                  DL_Log = DL(NLo)
                  L_last = Lxyz(NLo)
               ELSE IF( NLo == 1 ) THEN
                  DL_Log = DL(1)
                  L_last = DL(1)
               ELSE IF( NLo == 0 ) THEN
                  WRITE( UNIT = *, FMT = 6050) direction
                  STOP
               END IF
!
            ELSE
!
                WRITE(*,*)  'Route 2', DL_Log
               IF( NLo > 1 ) THEN
                  L_last = Lxyz(NLo)
               ELSE IF( NLo == 1 ) THEN
                  L_last = DL(1)
               ELSE IF( NLo == 0 ) THEN
                  L_last = 0.0e0_16
               END IF
!
            END IF IF_DL0
!
            WRITE(*,*)  'Parameters : ', L_max, L_Last, DL_Log, NLo, number_of_DLs, NL1, NL
            CALL Logarithmic_Discretization( axis   = direction, D_max    = L_max,     &
    &                                        D_last = L_Last,    D_Log    = DL_Log,    &
    &                                        Nb     = NLo,       IncrFact = IncrFact,  &
    &                                        DL_number = number_of_DLs,  Dd = DL(NL1: NL) )
!
            IF( direction == 'R' ) THEN
               DO i=1,number_of_DLs
                  radii(NLo+i) = radii(NLo+i-1) + DL(i)
               END DO
            END IF
!
! >>>>>>>>>
! >>>>>>>>> For non-logarithmic distributions
! >>>>>>>>>
!
         ELSE
!
            Case_NonLog: SELECT CASE(option(1:1))
            CASE('E','e','U','u')
!
               DL(NL1:NL) = Delta_L
!
               IF( direction == 'R' ) THEN
                  DO i=1,number_of_DLs
                     radii(NLo+i) = radii(NLo+i-1) + DL(i)
                  END DO
               END IF
!
            CASE('V','v')
!
               IF( TRIM(ADJUSTL(read_format)) == '(*)' .OR. TRIM(ADJUSTL(read_format)) == '*' ) THEN
                   READ( UNIT = *, FMT = * ) (DL(i),i=NL1,NL)
               ELSE IF( read_format(1:2) == '  ' ) THEN
                   WRITE( UNIT = *, FMT = 6521 ) direction2(1:mx), 'Variable', read_format, n1
                   STOP
               ELSE
                   READ( UNIT = *, FMT = read_format ) (DL(i),i=NL1,NL)
               END IF
!
               IF( direction == 'R' ) THEN
                  DO i=1,number_of_DLs
                     radii(NLo+i) = radii(NLo+i-1) + DL(i)
                  END DO
               END IF
!
            CASE('R','r')
!
               IF( direction == 'R' ) THEN
!
                  IF( TRIM(ADJUSTL(read_format)) == '(*)' .OR. TRIM(ADJUSTL(read_format)) == '*' ) THEN
                     IF( NLo == 0 ) THEN
                        READ( UNIT = *, FMT = * ) ( radii(i-1), i = NLo+1,NLo+number_of_radii )
                     ELSE
                        READ( UNIT = *, FMT = * ) ( radii(i), i = NLo+1,NLo+number_of_radii )
                     END IF
                  ELSE IF( read_format(1:2) == '  ' ) THEN
                     WRITE( UNIT = *, FMT = 6521 ) direction, 'Radii', read_format, n1
                     STOP
                  ELSE
                      READ( UNIT = *, FMT = read_format ) ( radii(i), i= NLo+1,NLo+number_of_radii )
                  END IF
!
                  NL = NLo + number_of_radii
!
                  IF( NLo == 0 ) THEN
                     X_ref = radii(0)
                     Rdev  = radii(0)
                     imax  = number_of_radii - 1
                  ELSE
                     Rdev = 0.0d0
                     imax = number_of_radii
                  END IF
!
                  NL = NLo + imax
!
                  DO i=1,imax
!
                     DL(NLo+i) = radii(NLo+i) - radii(NLo+i-1)
!
                     IF( DL(NLo+i) < 0.0d0 ) THEN
                        WRITE( UNIT = *, FMT = 6052) NLo+i, DL(NLo+i)
                        STOP
                     END IF
!
                  END DO
!
               END IF
!
            END SELECT Case_NonLog
!
         END IF IF_Log
!
! ...... Scaling
!
         DL(NL1:NL) = DL(NL1:NL)
!
!
        IF( NL1 == 1 ) THEN
            Lxyz(1) = DL(1) + Rdev
            NLb     = NL1 + 1
         ELSE
            NLb = NL1
         END IF
!
         DO i = NLb, NL
            Lxyz(i) = Lxyz(i-1) + DL(i)
         END DO
!
      END DO DO_LDiscr
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Discretize_Axis 1.0 .................... 08 July      2016',6X,'Discretize the axes of a Cartesian system for grid creation')
!
 6050 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>: The D',A,' distribution is logarithmic, but D',A,'_Log = 0 and this is the 1st subdivision',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6052 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>: The radial increment DR(',i3,') =',ES12.5,'<0: wrong values of the radii defining the increments',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6501 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>: The discretization axis descriptor = "',A,'" is in conflict with the expected value "',A,'"',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6515 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>: The dataset/namelist "',A,'" must be ended by the ":::<<<" descriptor',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6520 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>: The number of subdivisions ',A' = ',i5.5,' > the maximum declared number ',A,' = ',i5.5,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6521 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <RZ2D>: The ',A,'-discretization option is "',A,'" but the value of <read_format>=',A,' in dataset #',i3,' is invalid',//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6601 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>: There is a problem reading the namelist ',A,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6602 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>: There is a problem reading the namelist ',A,' in dataset #',i3,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6603 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>, namelist <',A,'>, dataset #',i3,': ',/,&
     &                'The value of the variable <',A,'> =',i4,'<0: NOT POSSIBLE', //, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6605 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure <Discretize_Axis>, namelist <',A,'>, dataset #',i3,': ',/,&
     &                'The value of the variable <option> = "',A,'" is not an acceptable option', //, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6850 FORMAT(//,22('ERROR-'),//,   &
     &       T10,'>>>  Procedure <Discretize_Axis>: Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
 6855 FORMAT(//,22('ERROR-'),//,   &
     &       T10,'>>>  Procedure <Discretize_Axis>: The angle units in the cylindrical system <',A,'> are of an unknown type',//,   &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Discretize_Axis>
!
!
         RETURN
!
      END SUBROUTINE Discretize_Axis
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Determine_Optimal_Ordering
!
!
      IMPLICIT NONE
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Determine_Optimal_Ordering>
!
      WRITE( UNIT = *, FMT = 6000 )
!
!
!
      IF_order: IF( coordinate_order == '   ' ) THEN
!
         IF_standard: IF( NX >= NY .AND. NX >= NZ ) THEN
            IF(NY >= NZ) THEN
               coord_order = 1
            ELSE
               coord_order = 2
            END IF
         ELSE IF( NY >= NX .AND. NY >= NZ ) THEN
            IF(NX >= NZ) THEN
               coord_order = 3
            ELSE
               coord_order = 4
            END IF
         ELSE IF( NZ >= NX .AND. NZ >= NY ) THEN
            IF(NX >= NY) THEN
               coord_order = 5
            ELSE
               coord_order = 6
            END IF
         END IF IF_standard
!
      ELSE
!
         CASE_specified: SELECT CASE (coordinate_order(1:3))
         CASE( 'XYZ', 'xyz' )
            coord_order = 1
         CASE( 'XZY', 'xzy' )
            coord_order = 2
         CASE( 'YXZ', 'yxz' )
            coord_order = 3
         CASE( 'YZX', 'yzx' )
            coord_order = 4
         CASE( 'ZXY', 'zxy' )
            coord_order = 5
         CASE( 'ZYX', 'zyx' )
            coord_order = 6
         END SELECT CASE_specified
!
      END IF IF_order
!
!
!
      IF( coord_order == 1 .OR. coord_order == 2 ) THEN
!
         NLXYZ(3)      = NX
         Lref(3)       = X_ref
         DLXYZ(1:NX,3) = DX(1:NX)
!
         IF( coord_order == 1 ) THEN
!
            NLXYZ(2)      = NY
            Lref(2)       = Y_ref
            DLXYZ(1:NY,2) = DY(1:NY)
!
            NLXYZ(1)      = NZ
            Lref(1)       = Z_ref
            DLXYZ(1:NZ,1) = DZ(1:NZ)
!
            GO TO 1000
!
         ELSE
!
            NLXYZ(2)      = NZ
            Lref(2)       = Z_ref
            DLXYZ(1:NZ,2) = DZ(1:NZ)
!
            NLXYZ(1)      = NY
            Lref(1)       = Y_ref
            DLXYZ(1:NY,1) = DY(1:NY)
!
            GO TO 1000
!
         END IF
!
      END IF
!
!
!
      IF( coord_order == 3 .OR. coord_order == 4 ) THEN
!
         NLXYZ(3)      = NY
         Lref(3)       = Y_ref
         DLXYZ(1:NY,3) = DY(1:NY)
!
         IF( coord_order == 3 ) THEN
!
            NLXYZ(2)      = NX
            Lref(2)       = X_ref
            DLXYZ(1:NX,2) = DX(1:NX)
!
            NLXYZ(1)      = NZ
            Lref(1)       = Z_ref
            DLXYZ(1:NZ,1) = DZ(1:NZ)
!
            GO TO 1000
!
         ELSE
!
            NLXYZ(2)      = NZ
            Lref(2)       = Z_ref
            DLXYZ(1:NZ,2) = DZ(1:NZ)
!
            NLXYZ(1)      = NX
            Lref(1)       = X_ref
            DLXYZ(1:NX,1) = DX(1:NX)
!
            GO TO 1000
!
         END IF
!
      END IF
!
!
!
      IF( coord_order == 5 .OR. coord_order == 6 ) THEN
!
         NLXYZ(3)      = NZ
         Lref(3)       = Z_ref
         DLXYZ(1:NZ,3) = DZ(1:NZ)
!
         IF( coord_order == 5 ) THEN
!
            NLXYZ(2)      = NX
            Lref(2)       = X_ref
            DLXYZ(1:NX,2) = DX(1:NX)
!
            NLXYZ(1)      = NY
            Lref(1)       = Y_ref
            DLXYZ(1:NY,1) = DY(1:NY)
!
            GO TO 1000
!
         ELSE
!
            NLXYZ(2)      = NY
            Lref(2)       = Y_ref
            DLXYZ(1:NY,2) = DY(1:NY)
!
            NLXYZ(1)      = NX
            Lref(1)       = X_ref
            DLXYZ(1:NX,1) = DX(1:NX)
!
         END IF
!
      END IF
!
!
!
 1000 IF((coord_order == 1) .OR. (coord_order == 3)) THEN
         RRR(1) =-1.0d0
      ELSE
         RRR(1) = 1.0d0
      END IF
!
      IF((coord_order == 2) .OR. (coord_order == 4)) THEN
         RRR(2) =-1.0d0
      ELSE
         RRR(2) = 1.0d0
      END IF
!
      IF((coord_order == 5) .OR. (coord_order == 6)) THEN
         RRR(3) =-1.0d0
      ELSE
         RRR(3) = 1.0d0
      END IF
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Determine_Optimal_Ordering 1.0 ......... 17 September 2015',6X,'Determine the optimal grid numbering to minimize the matrix bandwidth')
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Determine_Optimal_Ordering>
!
!
         RETURN
!
      END SUBROUTINE Determine_Optimal_Ordering
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Ordering_Relationship( coord_order, type_of_grid )
!
!
         IMPLICIT NONE
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER, INTENT(IN) :: coord_order, type_of_grid
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: First_call = .TRUE.
!
! ----------
! ...... Saving variables
! ----------
!
         SAVE First_call
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Ordering_Relationship>
!
         IF(First_call) THEN
            WRITE( UNIT = *, FMT = 6000 )
            First_call = .FALSE.
         END IF
!
         SELECT CASE(coord_order)
         CASE(1)
!
            DDXYZ(1) = CBi
            DDXYZ(2) = CMi
            DDXYZ(3) = CSm
!
            itrue = i
            jtrue = j
            ktrue = k
!
         CASE(2)
!
            DDXYZ(1) = CBi
            DDXYZ(2) = CSm
            DDXYZ(3) = CMi
!
            itrue = i
            jtrue = k
            ktrue = j
!
         CASE(3)
!
            DDXYZ(1) = CMi
            DDXYZ(2) = CBi
            DDXYZ(3) = CSm
!
            itrue = j
            jtrue = i
            ktrue = k
!
         CASE(4)
!
            DDXYZ(1) = CSm
            DDXYZ(2) = CBi
            DDXYZ(3) = CMi
!
            itrue = k
            jtrue = i
            ktrue = j
!
         CASE(5)
!
            DDXYZ(1) = CMi
            DDXYZ(2) = CSm
            DDXYZ(3) = CBi
!
            itrue = j
            jtrue = k
            ktrue = i
!
         CASE(6)
!
            DDXYZ(1) = CSm
            DDXYZ(2) = CMi
            DDXYZ(3) = CBi
!
            itrue = k
            jtrue = j
            ktrue = i
!
         END SELECT
!
         IF( type_of_grid == 1 ) RETURN
!
 1000    CONTINUE
!
         SELECT CASE(coord_order)
         CASE(1)
!
            DDXYZ(1) = CBi * COS(CMi)
            DDXYZ(2) = CBi * SIN(CMi)
!
         CASE(2)
!
            DDXYZ(1) = CBi * COS(CSm)
            DDXYZ(2) = CBi * SIN(CSm)
!
         CASE(3)
!
            DDXYZ(1) = CMi * COS(CBi)
            DDXYZ(2) = CMi * SIN(CBi)
!
         CASE(4)
!
            DDXYZ(1) = CSm * COS(CBi)
            DDXYZ(2) = CSm * SIN(CBi)
!
         CASE(5)
!
            DDXYZ(1) = CMi * COS(CSm)
            DDXYZ(2) = CMi * SIN(CSm)
!
         CASE(6)
!
            DDXYZ(1) = CSm * COS(CMi)
            DDXYZ(2) = CSm * SIN(CMi)
!
         END SELECT
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Ordering_Relationship 1.0 .............. 17 September 2015',6X,'Determine the relationship between the optimized and the standard ordering')
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Ordering_Relationship>
!
!
         RETURN
!
      END SUBROUTINE Ordering_Relationship
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Cylindrical_to_Cartesian
!
! ...... Integer variables
!
         INTEGER :: n
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Cylindrical_to_Cartesian>
!
         WRITE(UNIT = *, FMT = 6000)
!
         ALLOCATE( Xr(0:NX), X(1:NX), Y(1:NY), DRm(1:NX), DRp(1:NX), STAT=ier )
!
         IF(ier /= 0) THEN
            WRITE(UNIT = *, FMT = 6850) 'Xr, X, Y', 'Cylindrical_to_Cartesian'
            STOP
         END IF
!
! ...... The R-coordinates
!
         Xr(0) = X_ref
!
         DO n = 1, NX
            Xr(n) = Xr(n-1) + DX(n)
         END DO
!
         DO n = 1, NX
            X(n)   = SQRT( 5.0d-1 * ( Xr(n) * Xr(n) + Xr(n-1) * Xr(n-1) ) )
            DRm(n) = X(n)  - Xr(n-1)
            DRp(n) = Xr(n) - X(n)
         END DO
!
! ...... The Th-coordinates
!
         Y(1) = 5.0d-1 * DY(1)
!
         DO n = 2, NY
            Y(n) = Y(n-1) + 5.0d-1 * ( DY(n-1) + DY(n) )
         END DO
!
         IF( ABS(SUM(DY(1:NY)) - PI) > 1.0d-6) THEN
            WRITE(UNIT = *, FMT = 6200) SUM(DY(1:NY)), PI
            STOP
         END IF
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Cylindrical_to_Cartesian 1.0 ........... 25 November  2017',6X,'Transform cylindrical coordinates to Cartesian ones')
!
 6200 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'The sum of the angles in the cylindrical coordinate system is',ES15.8,' instead of PI =',ES15.8,//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 6850 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Cylindrical_to_Cartesian>
!
         RETURN
!
      END SUBROUTINE Cylindrical_to_Cartesian
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Coordinate_Processing( n, num, type_of_grid, DNd, CNi )
!
! ...... Double precision arrays
!
         REAL(KIND = 8), DIMENSION(2) :: DNd
!
! ...... Double precision variables
!
         REAL(KIND = 8), INTENT(OUT) :: CNi
!
! ...... Integer variables
!
         INTEGER, INTENT(IN) :: n, num, type_of_grid
!
         INTEGER :: np, nm
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: First_call = .TRUE.
!
! ----------
! ...... Saving variables
! ----------
!
         SAVE First_call
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Coordinate_Processing>
!
         IF(First_call) THEN
            WRITE( UNIT = *, FMT = 6000 )
            First_call   = .FALSE.
         END IF
!
         np = n+1
         nm = n-1
!
         IF( type_of_grid == 2 ) GO TO 1000
!
         DNd(1) = DLXYZ(n,num)/2.0d0
!
         IF( NLXYZ(num) > n ) THEN
            DNd(2) = DLXYZ(np,num)/2.0d0
         ELSE
            DNd(2) = 0.0d0
         END IF
!
         IF(n == 1) THEN
            CNi = CNi + RRR(num)*DLXYZ(1,num)/2.0d0
         ELSE
            CNi = CNi + RRR(num) * ( DLXYZ(n,num) + DLXYZ(nm,num) ) / 2.0d0
         END IF
!
         RETURN
!
 1000    IF( (num == 3 .AND. (coord_order == 1 .OR. coord_order == 2) ) .OR. &
     &       (num == 2 .AND. (coord_order == 3 .OR. coord_order == 5) ) .OR. &
     &       (num == 1 .AND. (coord_order == 4 .OR. coord_order == 6) ) )    &
     &   THEN
!
            DNd(1) = DRp(n)
!
            IF( NLXYZ(num) > n ) THEN
               DNd(2) = DRm(np)
            ELSE
               DNd(2) = 0.0d0
            END IF
!
            CNi = X(n)
!
         ELSE
!
            DNd(1) = DLXYZ(n,num)/2.0d0
!
            IF( NLXYZ(num) > n ) THEN
               DNd(2) = DLXYZ(np,num)/2.0d0
            ELSE
               DNd(2) = 0.0d0
            END IF
!
            IF(n == 1) THEN
               CNi = CNi + RRR(num)*DLXYZ(1,num)/2.0d0
            ELSE
               CNi = CNi + RRR(num) * ( DLXYZ(n,num) + DLXYZ(nm,num) ) / 2.0d0
            END IF
!
         END IF
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Coordinate_Processing 1.0 .............. 25 November  2017',6X,'Determining element dimensions and edge locations along coordinates')
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Coordinate_Processing>
!
         RETURN
!
      END SUBROUTINE Coordinate_Processing
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Determine_Reference_Surface(grid_type)
!
! -------------
! ......... Modules to be used
! -------------
!
         USE Het_Region_Variables, ONLY: RefSurface, Tabular
         USE Subdomain_Bounding_Surfaces, ONLY: Surface_Coordinate_Limits
!
         IMPLICIT NONE
!
! ...... User-defined variables
!
         TYPE(Tabular) :: table_1, table_2
!
! ...... Double precision arrays
!
         REAL(KIND = 8), DIMENSION(3,2) :: XYZ_limits
!
! ...... Double precision variables
!
         REAL(KIND = 8) :: ZZ_RefSurf, XX, YY
!
! ...... Integer variables
!
         INTEGER, INTENT(IN) :: grid_type
!
! ...... Integer variables
!
         INTEGER :: i, j, n, k
!
! ...... Logical variables
!
         LOGICAL :: find_RefSurface, correct_RefSurface
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Determine_Reference_Surface>
!
         WRITE(UNIT = *, FMT = 6000)
!
! >>>>>>
! >>>>>> For standard Cartesian systems
! >>>>>>
!
         IF( grid_type == 2 ) THEN
!
            CALL Cylindrical_to_Cartesian
!
            ALLOCATE( Z(0:NZ), RefSurfLoc(NX,NY), STAT=ier )
!
            IF(ier /= 0) THEN
               WRITE(UNIT = *, FMT = 6850) 'Z, RefSurfLoc', 'Determine_Reference_Surface'
               STOP
            END IF
!
            GO TO 1000
         END IF
!
         ALLOCATE( X(1:NX), Y(1:NY), Z(0:NZ), RefSurfLoc(NX,NY), STAT=ier )
!
         IF(ier /= 0) THEN
            WRITE(UNIT = *, FMT = 6850) 'X, Y, Z, RefSurfLoc', 'Determine_Reference_Surface'
            STOP
         END IF
!
! ...... The X-coordinates
!
         X(1) = X_ref + 5.0d-1 * DX(1)
!
         DO n = 2, NX
            X(n) = X(n-1) + 5.0d-1 * ( DX(n-1) + DX(n) )
         END DO
!
! ...... The Y-coordinates
!
         Y(1) = Y_ref + 5.0d-1 * DY(1)
!
         DO n = 2, NY
            Y(n) = Y(n-1) + 5.0d-1 * ( DY(n-1) + DY(n) )
         END DO
!
! ...... The Z-coordinates
!
 1000    Z(0) = Z_ref
         Z(1) = Z_ref - DZ(1)
!
         DO n = 2, NZ
            Z(n) = Z(n-1) - DZ(n)
         END DO
!
! ...... Determine the elevation of the reference surface
!
         find_RefSurface = .TRUE.
!
         n = 0
!
         DO_XGrid: DO i = 1, NX
!
            DO_YGrid: DO j = 1, NY
!
! ............ No vertical exclusions
!
               IF( grid_type == 2 ) THEN
                  XX = X(i) * COS(Y(j))
                  YY = X(i) * SIN(Y(j))
               ELSE
                  XX = X(i)
                  YY = Y(j)
               END IF
!
               excluded_ijk = Excluded_Element( coordinates = coordinates,                  &
     &                                          X = XX, Y = YY, Z = 0.0d0,                  &
     &                                          exclusion_zone_type = exclusion_zone_type,  &
     &                                          find_RefSurface     = find_RefSurface  )
!
               IF( excluded_ijk ) CYCLE DO_YGrid
!
               ZZ_RefSurf = 0.0d0
!
               IF( RefSurface%IntTableNum1 > 0 ) table_1 = IntTable(1)
!
               CALL Surface_Coordinate_Limits( subdomain        = RefSurface,              &
     &                                         table_1          = table_1,                 &
     &                                         table_2          = table_2,                 &
     &                                         XYZ_coordinates  = (/ XX, YY, 0.0d0 /),     &
     &                                         XYZ_limits       = XYZ_limits,              &
     &                                         find_RefSurface  = find_RefSurface,         &
     &                                         success          = correct_RefSurface )
!
               ZZ_RefSurf = XYZ_limits(3,1) + RefSurface%thick0  ! ... The reference surface offset by a fixed distance
!
               DO k = 1, NZ
!
                  IF( ZZ_RefSurf <= Z(k-1) .AND. ZZ_RefSurf > Z(k) ) THEN
                     RefSurfLoc(i,j) = k
                     n = n + 1
                     EXIT
                  END IF
!
               END DO
!
            END DO DO_YGrid
!
         END DO DO_XGrid
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Determine_Reference_Surface 1.0 ........ 30 August    2015',6X,'Determine the I,J,K coordinates of the Reference Surface')
!
 6850 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
 7000 FORMAT( 3(i6,2x) )

!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Determine_Reference_Surface>
!
         RETURN
!
      END SUBROUTINE Determine_Reference_Surface
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      SUBROUTINE Determine_Element_Status( XLoc, Yloc, Zloc, is_accepted, is_boundary, is_inclusion, discontinue_smoothing, medium_type, emissivity )
!
!
!
         IMPLICIT NONE
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND=8), INTENT(IN)  :: XLoc, Yloc, Zloc
         REAL(KIND=8), INTENT(OUT) :: emissivity
!
! ----------
! ...... Character variables
! ----------
!
         CHARACTER(LEN=5), INTENT(OUT) :: medium_type
!
         CHARACTER(LEN=1) :: boundary_id
         CHARACTER(LEN=5) :: rock
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL, INTENT(OUT) :: is_accepted, is_boundary, is_inclusion, discontinue_smoothing
!
         LOGICAL :: excluded_ijk, included_ijk
!
         LOGICAL :: First_call = .TRUE.
!
! ----------
! ...... Saving variables
! ----------
!
         SAVE First_call
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Determine_Element_Status>
!
         IF(First_call) THEN
            WRITE( UNIT = *, FMT = 6000 )
            First_call   = .FALSE.
         END IF
!
! >>>>>>
! >>>>>> Initialization
! >>>>>>
!
         is_accepted  = .FALSE.
         is_boundary  = .FALSE.
         is_inclusion = .FALSE.
         included_ijk = .FALSE.
!
         discontinue_smoothing = .FALSE.
!
         emissivity = dominant_medium_emissivity
!
! >>>>>>
! >>>>>> For standard cases
! >>>>>>
!
         IF( Num_ExclZones == 0 ) THEN
!
            excluded_ijk = .FALSE.
            included_ijk = .FALSE.
            GO TO 500
!
         END IF
!
! >>>>>>
! >>>>>> Treatment for exclusion and inclusion zones
! >>>>>>
!
         excluded_ijk = Excluded_Element( coordinates = coordinates, X = Xloc, Y = Yloc, Z = Zloc,            &
     &                                    exclusion_zone_type = exclusion_zone_type, use_InterpData = .TRUE., &
     &                                    completed_interpolation = completed_interpolation, &
     &                                    available_InterpData    = available_InterpData )
!
         IF( excluded_ijk ) THEN
!
            IF( Num_InclZones == 0 .OR. exclusion_zone_type == 'ABS' .OR. exclusion_zone_type == 'Abs' .OR. exclusion_zone_type == 'abs' ) THEN
               RETURN
            ELSE
!
               included_ijk = Included_Element( coordinates = coordinates, X = Xloc, Y = Yloc, Z = Zloc,       &
     &                                          inclusion_id = subdomain_id, medium = medium_type, emissivity = emissivity, &
     &                                          media_by_number = media_by_number,  use_InterpData = .TRUE.,   &
     &                                          completed_interpolation = completed_interpolation, &
     &                                          available_InterpData    = available_InterpData )
!
               IF( .NOT. included_ijk ) THEN
                  RETURN
               ELSE
                  is_inclusion = .TRUE.
                  GO TO 1000
               END IF
!
            END IF
!
         END IF
!
! >>>>>>
! >>>>>> Assigning properties
! >>>>>>
!
! ...... When the domain is homogeneous, assign the standard domain number/name
!
  500    IF_Media: IF(Flag_HetRegions .EQV. .FALSE.) THEN
!
            medium_type = '    1'
            emissivity  = dominant_medium_emissivity
!
! ----------
! ...... When the domain is heterogeneous ...
! ----------
!
         ELSE
!
! ......... Assign the dominant domain number/name
!
            IF_HetRegions: IF(TotNum_HetRegions == 1) THEN
!
               medium_type = dominant_medium
               emissivity  = dominant_medium_emissivity
!
! ......... Assign the appropriate heterogeneous domain number/name
!
            ELSE
               medium_type = Media_Region( coordinates = coordinates,  dominant_medium = dominant_medium, &
     &                                     X = Xloc, Y = Yloc, Z = Zloc,    emissivity = emissivity,      &
     &                                     media_by_number = media_by_number, use_InterpData = .TRUE.,    &
     &                                     completed_interpolation = completed_interpolation,             &
     &                                     available_InterpData    = available_InterpData )
            END IF IF_HetRegions
!
         END IF IF_Media
!
! ----------
! ...... Determine whether this is a boundary cell
! ----------
!
         IF_Bound: IF(Flag_Boundaries) THEN
!
            rock = medium_type
!
            CALL Boundary_Status( coordinates = coordinates, X = Xloc, Y = Yloc, Z = Zloc,           &
     &                            bound_id = subdomain_id,   medium = rock, emissivity = emissivity, &
     &                            media_by_number = media_by_number, use_InterpData = .TRUE., &
     &                            completed_interpolation = completed_interpolation,          &
     &                            available_InterpData    = available_InterpData )
!
            IF_id: IF( subdomain_id(1:1) == 'V' .OR. subdomain_id(1:1) == 'I' ) THEN
!
               is_boundary = .TRUE.
               medium_type = rock
               boundary_id = subdomain_id
!
               IF_Inclu: IF( Num_InclZones > 0 .AND. subdomain_id /= 'III' ) THEN
!
                  included_ijk = Included_Element( coordinates = coordinates, X = Xloc, Y = Yloc, Z = Zloc,           &
     &                                             inclusion_id  = subdomain_id, medium = medium_type, emissivity = emissivity, &
     &                                             media_by_number = media_by_number,  use_InterpData = .TRUE.,       &
     &                                             completed_interpolation = completed_interpolation,                 &
     &                                             available_InterpData    = available_InterpData )
!
                  IF( included_ijk ) THEN
                     medium_type  = medium_type
                     is_inclusion = .TRUE.
                     subdomain_id = '=>i'
                  END IF
!
               END IF IF_Inclu
!
            END IF IF_id
!
         END IF IF_Bound
!
! >>>>>>
! >>>>>> Accepting this element (main body, boundary or inclusion)
! >>>>>>
!
 1000    IF( medium_type(5:5) == '*' ) discontinue_smoothing = .TRUE.
!
         is_accepted = .TRUE.
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Determine_Element_Status 1.0 ........... 31 August    2015',6X,'Determine element occurrence and type at any location')
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Determine_Element_Status>
!
         RETURN
!
      END SUBROUTINE Determine_Element_Status
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Create_Element_List( i, j, k, medium_type, volume, elem_name )
!
!
!
         IMPLICIT NONE
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND = 8), INTENT(OUT) :: volume
!
         REAL(KIND = 8) :: AIJ, DD1, DD3
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER, INTENT(IN) :: i, j, k
!
! ----------
! ...... Character variables
! ----------
!
         CHARACTER( LEN = 5 ), INTENT(IN)  :: medium_type
!
         CHARACTER( LEN = ElemNameLength ), INTENT(OUT) :: elem_name
!
         CHARACTER( LEN = 5 ) :: ELNAME,  Name_of_5Character_Element
         CHARACTER( LEN = 8 ) :: ELNAME8, Name_of_8Character_Element
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: First_call = .TRUE.
!
! ----------
! ...... Saving variables
! ----------
!
         SAVE First_call
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Create_Element_List>
!
         IF(First_call) THEN
            WRITE( UNIT = *, FMT = 6000 )
            First_call   = .FALSE.
         END IF
!
         AIJ = 0.0d0  ! Heat exchange area initialization
!
! >>>>>>
! >>>>>> Computing volumes and surface areas - Cartesian coordinates (grid_type = 1)
! >>>>>>
!
         IF_reg: IF( grid_type /= 2 ) THEN
!
            volume = DLXYZ(I,3)*DLXYZ(J,2)*DLXYZ(K,1)
!
            IF( areas_for_HeatExch_Solution ) THEN
!
               IF((coord_order == 1) .OR. (coord_order == 3)) THEN
                  IF( k == 1 .OR. k == NZ )  AIJ = DLXYZ(j,2)*DLXYZ(i,3)
               END IF
!
               IF((coord_order == 2) .OR. (coord_order == 4)) THEN
                  IF( j == 1 .OR. j == NZ )  AIJ = DLXYZ(k,1)*DLXYZ(i,3)
               END IF
!
               IF((coord_order == 5) .OR. (coord_order == 6)) THEN
                  IF( i == 1 .OR. i == NZ ) AIJ = DLXYZ(k,1)*DLXYZ(j,2)
               END IF
!
            END IF
!
            GO TO 1000
!
         END IF IF_reg
!
! >>>>>>
! >>>>>> Computing volumes and surface areas - Transformed cylindrical coordinates (grid_type = 2)
! >>>>>>
!
         SELECT CASE(coord_order)
         CASE(1)
!
            volume = 5.0d-1 * ( Xr(i) * Xr(i) - Xr(i-1) * Xr(i-1) ) * DLXYZ(j,2) * DLXYZ(k,1)
!
            IF( areas_for_HeatExch_Solution ) THEN
               IF( k == 1 .OR. k == NZ )  AIJ = 5.0d-1 * ( Xr(i) * Xr(i) - Xr(i-1) * Xr(i-1) ) * DLXYZ(j,2)
            END IF
!
         CASE(2)
!
            volume = 5.0d-1 * ( Xr(i) * Xr(i) - Xr(i-1) * Xr(i-1) ) * DLXYZ(j,2) * DLXYZ(k,1)
!
            IF( areas_for_HeatExch_Solution ) THEN
               IF( j == 1 .OR. j == NZ )  AIJ = 5.0d-1 * ( Xr(i) * Xr(i) - Xr(i-1) * Xr(i-1) ) * DLXYZ(k,1)
            END IF
!
         CASE(3)
!
            volume = 5.0d-1 * ( Xr(j) * Xr(j) - Xr(j-1) * Xr(j-1) ) * DLXYZ(i,3) * DLXYZ(k,1)
!
            IF( areas_for_HeatExch_Solution ) THEN
               IF( k == 1 .OR. k == NZ )  AIJ = 5.0d-1 * ( Xr(j) * Xr(j) - Xr(j-1) * Xr(j-1) ) * DLXYZ(i,3)
            END IF
!
         CASE(4)
!
            volume = 5.0d-1 * ( Xr(k) * Xr(k) - Xr(k-1) * Xr(k-1) ) * DLXYZ(i,3) * DLXYZ(j,2)
!
            IF( areas_for_HeatExch_Solution ) THEN
               IF( j == 1 .OR. j == NZ )  AIJ = 5.0d-1 * ( Xr(k) * Xr(k) - Xr(k-1) * Xr(k-1) ) * DLXYZ(i,3)
            END IF
!
         CASE(5)
!
            volume = 5.0d-1 * ( Xr(j) * Xr(j) - Xr(j-1) * Xr(j-1) ) * DLXYZ(i,3) * DLXYZ(k,1)
!
            IF( areas_for_HeatExch_Solution ) THEN
               IF( i == 1 .OR. i == NZ )  AIJ = 5.0d-1 * ( Xr(j) * Xr(j) - Xr(j-1) * Xr(j-1) ) * DLXYZ(k,1)
            END IF
!
         CASE(6)
!
            volume = 5.0d-1 * ( Xr(k) * Xr(k) - Xr(k-1) * Xr(k-1) ) * DLXYZ(i,3) * DLXYZ(j,2)
!
            IF( areas_for_HeatExch_Solution ) THEN
               IF( i == 1 .OR. i == NZ )  AIJ = 5.0d-1 * ( Xr(k) * Xr(k) - Xr(k-1) * Xr(k-1) ) * DLXYZ(j,2)
            END IF
!
         END SELECT
!
         IF( (subdomain_id(1:1) == 'I' .OR. subdomain_id(1:1) == 'V') .AND. active_boundaries ) THEN
            volume = 1.0d50
            subdomain_id(1:1) = '*'
         END IF
!
! ...... Adjust global coordinates for inclined systems
!
 1000    IF( .NOT. ensure_grid_continuity ) THEN
            DD1 = DDXYZ(1) * BETA + DDXYZ(3) * BET
            DD3 = DDXYZ(3) * BETA - DDXYZ(1) * BET
         ELSE
            DD1 = DDXYZ(1)
            DD3 = DDXYZ(3)
         END IF
!
! >>>>>>
! >>>>>> Naming elements and printing element list
! >>>>>>
!
         IF(ElemNameLength == 8) THEN
            ELNAME8   = Name_of_8Character_Element(NLXYZ(1),NLXYZ(2),NLXYZ(3),i,j,k)
            elem_name = ELNAME8
         ELSE
            ELNAME    = Name_of_5Character_Element(NLXYZ(1),NLXYZ(2),NLXYZ(3),i,j,k)
            elem_name = ELNAME
         END IF
!
! ...... Print the element list
!
         IF_Format: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
!
            IF_Act: IF(subdomain_id(1:1) == ' ') THEN
!
               IF_8Char_A: IF(ElemNameLength == 8) THEN
!
                  IF(DOM == '    1') THEN
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5400) ELNAME8,Volume, medium_type, DD1, DDXYZ(2), DD3, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5400) ELNAME8,Volume, medium_type, DD1, DDXYZ(2), DD3
                     END IF
                  ELSE
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5500) ELNAME8,Volume, medium_type, DD1, DDXYZ(2), DD3, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5500) ELNAME8,Volume, medium_type, DD1, DDXYZ(2), DD3
                     END IF
                  END IF
!
               ELSE
!
                  IF(medium_type == '    1') THEN
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5402) ELNAME,Volume, medium_type, DD1, DDXYZ(2), DD3, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5402) ELNAME,Volume, medium_type, DD1, DDXYZ(2), DD3
                     END IF
                  ELSE
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5502) ELNAME,Volume, medium_type, DD1, DDXYZ(2), DD3, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5502) ELNAME,Volume, medium_type, DD1, DDXYZ(2), DD3
                     END IF
                  END IF
!
               END IF IF_8Char_A
!
            ELSE
!
               IF_8Char_B: IF(ElemNameLength == 8) THEN
!
                  IF(medium_type == '    1') THEN
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5401) ELNAME8,Volume, medium_type, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5401) ELNAME8,Volume, medium_type, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  ELSE
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5501) ELNAME8,Volume, medium_type, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5501) ELNAME8,Volume, medium_type, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  END IF
!
               ELSE
!
                  IF(medium_type == '    1') THEN
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5403) ELNAME,Volume, medium_type, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5403) ELNAME,Volume, medium_type, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  ELSE
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5503) ELNAME,Volume, medium_type, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5503) ELNAME,Volume, medium_type, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  END IF
!
               END IF IF_8Char_B
!
            END IF IF_Act
!
         ELSE
!
            IF_Name8: IF(ElemNameLength == 8) THEN
!
               IF_heat1: IF( areas_for_HeatExch_Solution ) THEN
!
                  IF_act1: IF( old_style_inactive_elements .AND. (subdomain_id(1:1) == 'I' .OR. subdomain_id(1:1) == 'V') ) THEN
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5010) ELNAME8, medium_type, Volume, AIJ, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(InactElem_Unit, FMT = 5010) ELNAME8, medium_type, Volume, AIJ, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  ELSE
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5010) ELNAME8, medium_type, Volume, AIJ, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5010) ELNAME8, medium_type, Volume, AIJ, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  END IF IF_act1
!
               ELSE
!
                  IF_act2: IF( old_style_inactive_elements .AND. (subdomain_id(1:1) == 'I' .OR. subdomain_id(1:1) == 'V') ) THEN
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5510) ELNAME8, medium_type, Volume, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(InactElem_Unit, FMT = 5510) ELNAME8, medium_type, Volume, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  ELSE
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5510) ELNAME8, medium_type, Volume, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5510) ELNAME8, medium_type, Volume, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  END IF IF_act2
!
               END IF IF_heat1
!
            ELSE
!
               IF_heat2: IF( areas_for_HeatExch_Solution ) THEN
!
                  IF_act3: IF( old_style_inactive_elements .AND. (subdomain_id(1:1) == 'I' .OR. subdomain_id(1:1) == 'V') ) THEN
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 10) ELNAME, medium_type, Volume, AIJ, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(InactElem_Unit, FMT = 10) ELNAME, medium_type, Volume, AIJ, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  ELSE
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 10) ELNAME, medium_type, Volume, AIJ, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 10) ELNAME, medium_type, Volume, AIJ, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  END IF IF_act3
!
               ELSE
!
                  IF_act4: IF( old_style_inactive_elements .AND. (subdomain_id(1:1) == 'I' .OR. subdomain_id(1:1) == 'V') ) THEN
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5610) ELNAME, medium_type, Volume, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(InactElem_Unit, FMT = 5610) ELNAME, medium_type, Volume, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  ELSE
                     IF( ensure_grid_continuity ) THEN
                        WRITE(output_file_unit(1), FMT = 5610) ELNAME, medium_type, Volume, DD1, DDXYZ(2), DD3, subdomain_id, itrue, jtrue, ktrue
                     ELSE
                        WRITE(output_file_unit(1), FMT = 5610) ELNAME, medium_type, Volume, DD1, DDXYZ(2), DD3, subdomain_id
                     END IF
                  END IF IF_act4
!
               END IF IF_heat2
!
            END IF IF_Name8
!
         END IF IF_Format
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Create_Element_List 1.0 ................ 31 August    2015',6X,'Create the ELEMENT data block for TOUGH2/TOUGH+')
!
 5010 FORMAT(A8, 7X,A5,1ES10.4E2,1ES10.4E1,10X,3ES10.4E1,1X,A3,3(1x,i5))
 5510 FORMAT(A8, 7X,A5,1ES10.4E2,20X,3ES10.4E1,1X,A3,3(1x,i5))
   10 FORMAT(A5,10X,A5,ES10.4E2,ES10.4E1,10X,3ES10.4E1,1X,A3,3(1x,i5))
 5610 FORMAT(A5,10X,A5,1ES10.4E2,20X,3ES10.4E1,1X,A3,3(1x,i5))
!
 5400 FORMAT(A,'  &Elem  V=',ES17.10,', MedNum=',A5,', x=',ES17.10,', y=',ES17.10,', z=',ES17.10,', i=',i5,', j=',i5,'k=',i5,'/')
 5401 FORMAT(A,'  &Elem  V=',ES17.10,', MedNum=',A5,', x=',ES17.10,', y=',ES17.10,', z=',ES17.10,', act="',A,'", i=',i5,', j=',i5,'k=',i5,'/')
!
 5500 FORMAT(A,'  &Elem  V=',ES17.10,', MedNam="',A5,'", x=',ES17.10,', y=',ES17.10,', z=',ES17.10,', i=',i5,', j=',i5,'k=',i5,'/')
 5501 FORMAT(A,'  &Elem  V=',ES17.10,', MedNam="',A5,'", x=',ES17.10,', y=',ES17.10,', z=',ES17.10,', act="',A,'", i=',i5,', j=',i5,'k=',i5,'/')
!
 5402 FORMAT(A,'  &Elem  V=',ES17.10,', MedNum=',A5,', x=',ES17.10,', y=',ES17.10,', z=',ES17.10,', i=',i5,', j=',i5,'k=',i5,'/')
 5403 FORMAT(A,'  &Elem  V=',ES17.10,', MedNum=',A5,', x=',ES17.10,', y=',ES17.10,', z=',ES17.10,', act="',A,'", i=',i5,', j=',i5,'k=',i5,'/')
!
 5502 FORMAT(A,'  &Elem  V=',ES17.10,', MedNam="',A5,'", x=',ES17.10,', y=',ES17.10,', z=',ES17.10,', i=',i5,', j=',i5,'k=',i5,'/')
 5503 FORMAT(A,'  &Elem  V=',ES17.10,', MedNam="',A5,'", x=',ES17.10,', y=',ES17.10,', z=',ES17.10,', act="',A,'", i=',i5,', j=',i5,'k=',i5,'/')
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Create_Element_List>
!
         RETURN
!
      END SUBROUTINE Create_Element_List
!
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Write_VTK_Data_File( itrue, jtrue, ktrue, elem_name )
!
!
!
         IMPLICIT NONE
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND = 8) :: XXX_1, XXX_2, YYY_1, YYY_2, ZZZ_1, ZZZ_2
         REAL(KIND = 8) :: XX1Z1, XX1Z2, XX2Z1, XX2Z2, ZZ1X1, ZZ1X2, ZZ2X1, ZZ2X2
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER, INTENT(IN) :: itrue, jtrue, ktrue
!
         INTEGER :: ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8
         INTEGER :: NXplus1, NXYext
!
! ...... Character variables
!
         CHARACTER( LEN = ElemNameLength ), INTENT(IN) :: elem_name
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: First_call = .TRUE.
!
! ----------
! ...... Saving variables
! ----------
!
         SAVE First_call, NXplus1, NXYext
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Write_VTK_Data_File>
!
         IF(First_call) THEN
!
            WRITE( UNIT = *, FMT = 6000 )
            First_call   = .FALSE.
!
            NXplus1 = NX + 1
            NXYext  = NXplus1 * ( NY + 1 )
!
         END IF
!
! >>>>>>
! >>>>>> Writing data in VTK files and the ELEME block in MESH
! >>>>>>
!
         IF( itrue == 1 ) THEN
            XXX_1 = X_ref
            XXX_2 = X_ref + DX(1)
         ELSE
            XXX_1 = X_ref + SUM( DX(1:itrue - 1) )
            XXX_2 = X_ref + SUM( DX(1:itrue) )
         END IF
!
         IF( jtrue == 1 ) THEN
            YYY_1 = Y_ref
            YYY_2 = Y_ref + DY(1)
         ELSE
            YYY_1 = Y_ref + SUM( DY(1:jtrue - 1) )
            YYY_2 = Y_ref + SUM( DY(1:jtrue) )
         END IF
!
         IF( ktrue == 1 ) THEN
            ZZZ_2 = Z_ref
            ZZZ_1 = Z_ref - DZ(1)
         ELSE
            ZZZ_2 = Z_ref - SUM( DZ(1:ktrue - 1) )
            ZZZ_1 = Z_ref - SUM( DZ(1:ktrue) )
         END IF
!
! ...... Adjust global coordinates for inclined systems
!
         IF( .NOT. ensure_grid_continuity ) THEN
!
            XX1Z1 = XXX_1 * BETA + ZZZ_1 * BET
            XX1Z2 = XXX_1 * BETA + ZZZ_2 * BET
            XX2Z1 = XXX_2 * BETA + ZZZ_1 * BET
            XX2Z2 = XXX_2 * BETA + ZZZ_2 * BET
!
            ZZ1X1 = ZZZ_1 * BETA - XXX_1 * BET
            ZZ1X2 = ZZZ_1 * BETA - XXX_2 * BET
            ZZ2X1 = ZZZ_2 * BETA - XXX_1 * BET
            ZZ2X2 = ZZZ_2 * BETA - XXX_2 * BET
!
         ELSE
!
            XX1Z1 = XXX_1
            XX1Z2 = XXX_1
            XX2Z1 = XXX_2
            XX2Z2 = XXX_2
!
            ZZ1X1 = ZZZ_1
            ZZ1X2 = ZZZ_1
            ZZ2X1 = ZZZ_2
            ZZ2X2 = ZZZ_2
!
         END IF
!
! :::::::::::::
!
         ijk_1 = itrue + (jtrue - 1) * NXplus1 + ktrue * NXYext
         ijk_2 = itrue + 1 + (jtrue - 1) * NXplus1 + ktrue * NXYext
!
         ijk_3 = itrue + jtrue * NXplus1 + ktrue * NXYext
         ijk_4 = itrue + 1 + jtrue * NXplus1 + ktrue * NXYext
!
         ijk_5 = itrue + (jtrue - 1) * NXplus1 + ( ktrue - 1) * NXYext
         ijk_6 = itrue + 1 + (jtrue - 1) * NXplus1 + ( ktrue - 1) * NXYext
!
         ijk_7 = itrue + jtrue * NXplus1 + ( ktrue - 1) * NXYext
         ijk_8 = itrue + 1 + jtrue * NXplus1 + ( ktrue - 1) * NXYext
!
! :::::::::::::
!
         IF( .NOT. renumbered_inactive_elements ) THEN
!
            NumElem = NumElem + 1
!
            IF( NumElem == 1 ) THEN
               ijk_min = MIN( ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
               ijk_max = MAX( ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
            ELSE
               ijk_min = MIN( ijk_min, ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
               ijk_max = MAX( ijk_max, ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
            END IF
!
            WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX1Z1, YYY_1, ZZ1X1, ijk_1, elem_name, NumElem
            WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX2Z1, YYY_1, ZZ1X2, ijk_2
            WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX2Z1, YYY_2, ZZ1X2, ijk_4
            WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX1Z1, YYY_2, ZZ1X1, ijk_3
            WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX1Z2, YYY_1, ZZ2X1, ijk_5
            WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX2Z2, YYY_1, ZZ2X2, ijk_6
            WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX2Z2, YYY_2, ZZ2X2, ijk_8
            WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX1Z2, YYY_2, ZZ2X1, ijk_7
!
         ELSE
!
            IF( subdomain_id(1:1) == 'I' .OR. subdomain_id(1:1) == 'V' ) THEN
!
               NumInactElem = NumInactElem + 1
!
               IF( NumInactElem == 1 ) THEN
                  ijk_min_IE = MIN( ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
                  ijk_max_IE = MAX( ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
               ELSE
                  ijk_min_IE = MIN( ijk_min, ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
                  ijk_max_IE = MAX( ijk_max, ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
               END IF
!
               WRITE( UNIT = InactElement_Unit, FMT = 6100 ) XX1Z1, YYY_1, ZZ1X1, ijk_1, elem_name, NumInactElem
               WRITE( UNIT = InactElement_Unit, FMT = 6100 ) XX2Z1, YYY_1, ZZ1X2, ijk_2
               WRITE( UNIT = InactElement_Unit, FMT = 6100 ) XX2Z1, YYY_2, ZZ1X2, ijk_4
               WRITE( UNIT = InactElement_Unit, FMT = 6100 ) XX1Z1, YYY_2, ZZ1X1, ijk_3
               WRITE( UNIT = InactElement_Unit, FMT = 6100 ) XX1Z2, YYY_1, ZZ2X1, ijk_5
               WRITE( UNIT = InactElement_Unit, FMT = 6100 ) XX2Z2, YYY_1, ZZ2X2, ijk_6
               WRITE( UNIT = InactElement_Unit, FMT = 6100 ) XX2Z2, YYY_2, ZZ2X2, ijk_8
               WRITE( UNIT = InactElement_Unit, FMT = 6100 ) XX1Z2, YYY_2, ZZ2X1, ijk_7
!
            ELSE
!
               NumElem = NumElem + 1
!
               IF( NumElem == 1 ) THEN
                  ijk_min = MIN( ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
                  ijk_max = MAX( ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
               ELSE
                  ijk_min = MIN( ijk_min, ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
                  ijk_max = MAX( ijk_max, ijk_1, ijk_2, ijk_3, ijk_4, ijk_5, ijk_6, ijk_7, ijk_8 )
               END IF
!
               WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX1Z1, YYY_1, ZZ1X1, ijk_1, elem_name, NumElem
               WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX2Z1, YYY_1, ZZ1X2, ijk_2
               WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX2Z1, YYY_2, ZZ1X2, ijk_4
               WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX1Z1, YYY_2, ZZ1X1, ijk_3
               WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX1Z2, YYY_1, ZZ2X1, ijk_5
               WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX2Z2, YYY_1, ZZ2X2, ijk_6
               WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX2Z2, YYY_2, ZZ2X2, ijk_8
               WRITE( UNIT = Corner_Unit, FMT = 6100 ) XX1Z2, YYY_2, ZZ2X1, ijk_7
!
            END IF
!
         END IF
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Write_VTK_Data_File 1.0 ................ 10 September 2015',6X,'Write data using the appropriate format for VTK-based visualization')
!
 6100 FORMAT( 3(2x,ES12.5E2,2x),i10,2x,A,2x,i8 )
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Write_VTK_Data_File>
!
         RETURN
!
      END SUBROUTINE Write_VTK_Data_File
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      SUBROUTINE Create_Connection_List( i, j, k, accepted_ijk, accepted_ijk_ext, outside_I_Range_limits )
!
!
!
         IMPLICIT NONE
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND = 8) :: A12, A13, A23, Fc, Fr, Fcr
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER, INTENT(IN) :: i, j, k
!
         INTEGER :: first, second, third, i1, j1, k1, i2, j2, n1, ijk, next1, next2, next3
!
! ----------
! ...... Character variables
! ----------
!
         CHARACTER( LEN = 5 ) :: ELNAME,  ELNAM1,  Name_of_5Character_Element
         CHARACTER( LEN = 8 ) :: ELNAME8, ELNAM18, Name_of_8Character_Element
!
! ----------
! ...... Logical arrays
! ----------
!
         LOGICAL, INTENT(IN), DIMENSION(3) :: accepted_ijk_ext
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL, INTENT(IN) :: accepted_ijk
!
         LOGICAL, INTENT(OUT) :: outside_I_Range_limits
!
         LOGICAL :: First_call = .TRUE., printout
!
! ----------
! ...... Saving variables
! ----------
!
         SAVE First_call, next1, next2, next3
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Create_Connection_List>
!
         IF_1st: IF(First_call) THEN
!
            WRITE( UNIT = *, FMT = 6000 )
            First_call   = .FALSE.
!
            next1 = 1
            next2 = 1
            next3 = 1
!
! ......... Initializing the DX_adj and DY_adj to be used for the storage of the adjusted DX, DY
!
            IF_cont1: IF( ensure_grid_continuity ) THEN
!
               DX_adj = 0.0d0
               DY_adj = 0.0d0
!
! ............ Assigning initial DX values
!
               DO_n1: DO n1 = 1, NumTotalElem
!
                  ijk = ijk_location(n1)
!
                  k1 = ijk/NXNY
                  j2 = MOD(ijk,NXNY)
!
                  IF(j2 == 0 ) THEN
                     j1 = NY
                     i1 = NX
                  ELSE
                     k1 = k1 + 1
                     j1 = j2/NX
                     i2 = MOD(j2,NX)
                     IF( i2 == 0 ) THEN
                        i1 = NX
                     ELSE
                        j1 = j1 + 1
                        i1 = i2
                     END IF
                  END IF
!
                  IF( grid_type == 1 ) THEN
                     DX_adj(n1) = DX(i1)
                     DY_adj(n1) = DY(j1)
                  ELSE
                     DX_adj(n1) = DX(i1)
                     DY_adj(n1) = 5.0d-1 * ( Xr(i1-1) + Xr(i1) ) * DY(j1)
                  END IF
!
               END DO DO_n1
!
            END IF IF_cont1
!
         END IF IF_1st
!
! >>>>>>
! >>>>>> Computing connection areas
! >>>>>>
!
         outside_I_Range_limits = .FALSE.
!
! ...... For Cartesian systems
!
         IF_Type: IF( grid_type == 1 ) THEN
!
            A12 = DLXYZ(j,2)*DLXYZ(k,1)
            A13 = DLXYZ(i,3)*DLXYZ(k,1)
            A23 = DLXYZ(i,3)*DLXYZ(j,2)
!
! ...... For transformed cylindrical systems
!
         ELSE
!
            SELECT CASE(coord_order)
            CASE(1)
!
               A12 = Xr(i) * DLXYZ(j,2) * DLXYZ(k,1)
               A13 = DX(i) * DLXYZ(k,1)
               A23 = 5.0d-1 * ( Xr(i) * Xr(i) - Xr(i-1) * Xr(i-1) ) * DLXYZ(j,2)
!
            CASE(2)
!
               A12 = Xr(i) * DLXYZ(j,2) * DLXYZ(k,1)
               A13 = 5.0d-1 * ( Xr(i) * Xr(i) - Xr(i-1) * Xr(i-1) ) * DLXYZ(k,1)
               A23 = DX(i) * DLXYZ(j,2)
!
            CASE(3)
!
               A12 = DX(j) * DLXYZ(k,1)
               A13 = Xr(j) * DLXYZ(i,3) * DLXYZ(k,1)
               A23 = 5.0d-1 * ( Xr(j) * Xr(j) - Xr(j-1) * Xr(j-1) ) * DLXYZ(i,3)
!
            CASE(4)
!
               A12 = DX(k) * DLXYZ(j,2)
               A13 = 5.0d-1 * ( Xr(k) * Xr(k) - Xr(k-1) * Xr(k-1) ) * DLXYZ(i,3)
               A23 = Xr(k) * DLXYZ(i,3) * DLXYZ(j,2)
!
            CASE(5)
!
               A12 = 5.0d-1 * ( Xr(j) * Xr(j) - Xr(j-1) * Xr(j-1) ) * DLXYZ(k,1)
               A13 = Xr(j) * DLXYZ(i,3) * DLXYZ(k,1)
               A23 = DX(j) * DLXYZ(i,3)
!
            CASE(6)
!
               A12 = 5.0d-1 * ( Xr(k) * Xr(k) - Xr(k-1) * Xr(k-1) ) * DLXYZ(j,2)
               A13 = DX(k) * DLXYZ(j,2)
               A23 = Xr(k) * DLXYZ(i,3) * DLXYZ(j,2)
!
            END SELECT
!
         END IF IF_Type
!
! >>>>>>
! >>>>>> Checking the (i,j,k)-element
! >>>>>>
!
         IF( .NOT. accepted_ijk ) THEN
            RETURN
         ELSE
            IF(ElemNameLength == 8) THEN
               ELNAME8 = Name_of_8Character_Element(NLXYZ(1),NLXYZ(2),NLXYZ(3),i,j,k)
            ELSE
               ELNAME = Name_of_5Character_Element(NLXYZ(1),NLXYZ(2),NLXYZ(3),i,j,k)
            END IF
         END IF
!
! >>>>>>
! >>>>>> Checking the (i+1,j,k)-element
! >>>>>>
!
         IF_Icheck: IF(I < NLXYZ(3)) THEN
!
            IF( .NOT. accepted_ijk_ext(1) ) THEN
               GO TO 1000
            ELSE
!
               Fc = 1.0d0
               Fr = 1.0d0
               beta_x = 0.0d0
               beta_y = 0.0d0
               beta_z = BETA
               IF( ABS(beta_z) < 1.0d-9 ) beta_z = 0.0d0
!
               IF( coord_order == 1 ) THEN
!
                  first  = i+1
                  second = j
                  third  = k + kx_offset
!
                  IF_cont2: IF( ensure_grid_continuity .AND. kx_offset /= 0 ) THEN
!
                     Fc     = coeff(1)
                     beta_x = COS(phi_x)
                     IF( ABS(beta_x) < 1.0d-9 ) beta_x = 0.0d0
!
                     ijk = i + (j - 1) * NX + (k - 1) * NXNY
                     DO n1 = next1, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(i) + Fc * D1d(1)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRp(i) + Fc * D1d(1)
                           END IF
                           next1 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + 1 + kx_offset * NXNY
                     DO n1 = next1, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(i+1) + Fc * D1d(2)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRm(i+1) + Fc * D1d(2)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  ELSE
!
                     beta_x = BET
!
                  END IF IF_cont2
!
               ELSE IF( coord_order == 2 ) THEN
!
                  first  = i+1
                  second = j + kx_offset
                  third  = k
!
                  IF_cont3: IF( ensure_grid_continuity .AND. kx_offset /= 0 ) THEN
!
                     Fc     = coeff(1)
                     beta_x = COS(phi_x)
                     IF( ABS(beta_x) < 1.0d-9 ) beta_x = 0.0d0
!
                     ijk = i + (k - 1) * NX + (j - 1) * NXNY
                     DO n1 = next1, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(i) + Fc * D1d(1)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRp(i) + Fc * D1d(1)
                           END IF
                           next1 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + 1 + kx_offset * NXNY
                     DO n1 = next1, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(i+1) + Fc * D1d(2)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRm(i+1) + Fc * D1d(2)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  ELSE
!
                     beta_x = BET
!
                  END IF IF_cont3
!
               ELSE IF( coord_order == 3 ) THEN
!
                  first  = i+1
                  second = j
                  third  = k + ky_offset
!
                  IF_cont4: IF( ensure_grid_continuity ) THEN
!
                     IF( ky_offset == 0 ) THEN
                        IF( grid_type == 2 ) Fr = X(j)
                        GO TO 500
                     END IF
!
                     Fc = coeff(2)
                     beta_y = COS(phi_y)
                     IF( ABS(beta_y) < 1.0d-9 ) beta_y = 0.0d0
!
                     ijk = j + (i - 1) * NX + (k - 1) * NXNY
                     DO n1 = next1, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(i) + Fc * D1d(1)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(j-1) + Xr(j) ) * ( 5.0d-1 * DY(i) - Fc * D1d(1) )
                              Fr = X(j)
                           END IF
                           next1 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + NX + ky_offset * NXNY
                     DO n1 = next1, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(i+1) + Fc * D1d(2)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(j-1) + Xr(j) ) * ( 5.0d-1 * DY(i+1) - Fc * D1d(2) )
                              Fr = X(j)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  END IF IF_cont4
!
               ELSE IF( coord_order == 4 ) THEN
!
                  first  = i+1
                  second = j + ky_offset
                  third  = k
!
                  IF_cont5: IF( ensure_grid_continuity ) THEN
!
                     IF( ky_offset == 0 ) THEN
                        IF( grid_type == 2 ) Fr = X(k)
                        GO TO 500
                     END IF
!
                     Fc = coeff(2)
                     beta_y = COS(phi_y)
                     IF( ABS(beta_y) < 1.0d-9 ) beta_y = 0.0d0
!
                     ijk = k + (i - 1) * NX + (j - 1) * NXNY
                     DO n1 = next1, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(i) + Fc * D1d(1)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(k-1) + Xr(k) ) * ( 5.0d-1 * DY(i) - Fc * D1d(1) )
                              Fr = X(k)
                           END IF
                           next1 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + NX + ky_offset * NXNY
                     DO n1 = next1, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(i+1) + Fc * D1d(2)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(k-1) + Xr(k) ) * ( 5.0d-1 * DY(i+1) - Fc * D1d(2) )
                              Fr = X(k)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  END IF IF_cont5
!
               ELSE
!
                  Fc = 1.0d0
                  Fr = 1.0d0
                  beta_x = BET
                  beta_z = BETA
                  IF( ABS(beta_z) < 1.0d-9 ) beta_z = 0.0d0
                  first  = i+1
                  second = j
                  third  = k
!
               END IF
!
! ............ For partial processing of connections: need to compute connection to I-1 to avoid discontinuities
!
  500          IF( partial_processing .AND. (I == LLoop_begin - 1) .AND. (LLoop_begin /= 1) ) THEN
                  outside_I_Range_limits = .TRUE.
                  IF( ensure_grid_continuity ) THEN
                     GO TO 1000
                  ELSE
                     RETURN
                  END IF
               END IF
!
               NumConx = NumConx + 1
!
               IF(ElemNameLength == 8) THEN
                  ELNAM18 = Name_of_8Character_Element( NLXYZ(1), NLXYZ(2), NLXYZ(3), first, second, third )
               ELSE
                  ELNAM1 = Name_of_5Character_Element( NLXYZ(1), NLXYZ(2), NLXYZ(3), first, second, third )
               END IF
!
               IF( emissivity /= emissivity_p ) emissivity = 0.0d0
!
            END IF
!
! ......... Writing the connection information
!
            IF( ABS(beta_x) < 1.0d-10 ) beta_x = 0.0d0
            IF( ABS(beta_y) < 1.0d-10 ) beta_y = 0.0d0
            IF( ABS(beta_z) < 1.0d-10 ) beta_z = 0.0d0
!
            Fcr = Fc*Fr
!
            IF_Format3: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
               IF(ElemNameLength == 8) THEN
                  IF((coord_order == 1) .OR. (coord_order == 2)) WRITE(output_file_unit(2),5521) ELNAME8,ELNAM18,A12,Fcr*D1d(1),Fcr*D1d(2),NA(2),beta_x,emissivity
                  IF((coord_order == 3) .OR. (coord_order == 4)) WRITE(output_file_unit(2),5521) ELNAME8,ELNAM18,A12,Fcr*D1d(1),Fcr*D1d(2),NA(3),beta_y,emissivity
                  IF((coord_order == 5) .OR. (coord_order == 6)) WRITE(output_file_unit(2),5521) ELNAME8,ELNAM18,A12,Fcr*D1d(1),Fcr*D1d(2),NA(4),beta_z,emissivity
               ELSE
                  IF((coord_order == 1) .OR. (coord_order == 2)) WRITE(output_file_unit(2),5531) ELNAME,ELNAM1,A12,Fcr*D1d(1),Fcr*D1d(2),NA(2),beta_x,emissivity
                  IF((coord_order == 3) .OR. (coord_order == 4)) WRITE(output_file_unit(2),5531) ELNAME,ELNAM1,A12,Fcr*D1d(1),Fcr*D1d(2),NA(3),beta_y,emissivity
                  IF((coord_order == 5) .OR. (coord_order == 6)) WRITE(output_file_unit(2),5531) ELNAME,ELNAM1,A12,Fcr*D1d(1),Fcr*D1d(2),NA(4),beta_z,emissivity
               END IF
            ELSE
               IF(ElemNameLength == 8) THEN
                  IF((coord_order == 1) .OR. (coord_order == 2)) WRITE(UNIT = output_file_unit(2), FMT = 5011, ADVANCE = 'NO') ELNAME8,ELNAM18,NA(2),Fcr*D1d(1),Fcr*D1d(2),A12,beta_x
                  IF((coord_order == 3) .OR. (coord_order == 4)) WRITE(UNIT = output_file_unit(2), FMT = 5011, ADVANCE = 'NO') ELNAME8,ELNAM18,NA(3),Fcr*D1d(1),Fcr*D1d(2),A12,beta_y
                  IF((coord_order == 5) .OR. (coord_order == 6)) WRITE(UNIT = output_file_unit(2), FMT = 5011, ADVANCE = 'NO') ELNAME8,ELNAM18,NA(4),Fcr*D1d(1),Fcr*D1d(2),A12,beta_z
                  IF( emissivity > 0.0e0 ) THEN
                     WRITE(UNIT = output_file_unit(2), FMT = '(ES10.4E1)' ) emissivity
                  ELSE
                     WRITE(UNIT = output_file_unit(2), FMT = '(A)' ) '          '
                  END IF
               ELSE
                  IF((coord_order == 1) .OR. (coord_order == 2)) WRITE(UNIT = output_file_unit(2), FMT = 11, ADVANCE = 'NO') ELNAME,ELNAM1,NA(2),Fcr*D1d(1),Fcr*D1d(2),A12,beta_x
                  IF((coord_order == 3) .OR. (coord_order == 4)) WRITE(UNIT = output_file_unit(2), FMT = 11, ADVANCE = 'NO') ELNAME,ELNAM1,NA(3),Fcr*D1d(1),Fcr*D1d(2),A12,beta_y
                  IF((coord_order == 5) .OR. (coord_order == 6)) WRITE(UNIT = output_file_unit(2), FMT = 11, ADVANCE = 'NO') ELNAME,ELNAM1,NA(4),Fcr*D1d(1),Fcr*D1d(2),A12,beta_z
                  IF( emissivity > 0.0e0 ) THEN
                     WRITE(UNIT = output_file_unit(2), FMT = '(ES10.4E1)' ) emissivity
                  ELSE
                     WRITE(UNIT = output_file_unit(2), FMT = '(A)' ) '          '
                  END IF
               END IF
!
            END IF IF_Format3
!
         END IF IF_Icheck
!
! >>>>>>
! >>>>>> Checking the (i,j+1,k)-element
! >>>>>>
!
 1000    IF_Jcheck: IF(J < NLXYZ(2)) THEN
!
            printout = .TRUE.
!
            IF( .NOT. accepted_ijk_ext(2) ) THEN
               GO TO 2000
            ELSE
!
               Fc = 1.0d0
               Fr = 1.0d0
               beta_x = 0.0d0
               beta_y = 0.0d0
               beta_z = BETA
               IF( ABS(beta_z) < 1.0d-9 ) beta_z = 0.0d0
!
               IF( coord_order == 3 ) THEN
!
                  first  = i
                  second = j+1
                  third  = k + kx_offset
!
                  IF_cont6: IF( ensure_grid_continuity ) THEN
!
                     IF( partial_processing ) THEN
                        IF( i < LLoop_begin .OR. i > LLoop_end ) RETURN
                     END IF
!
                     IF( X(j) < X(jmin(i)) ) printout = .FALSE.
!
                     IF( kx_offset == 0 ) GO TO 1500
!
                     Fc = coeff(1)
                     beta_x = COS(phi_x)
                     IF( ABS(beta_x) < 1.0d-9 ) beta_x = 0.0d0
!
                     ijk = j + (i - 1) * NX + (k - 1) * NXNY
                     DO n1 = next2, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(j) + Fc * D2d(1)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRp(j) + Fc * D2d(1)
                           END IF
                           next2 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + 1 + kx_offset * NXNY
                     DO n1 = next2, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(j+1) + Fc * D2d(2)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRm(j+1) + Fc * D2d(2)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  ELSE
!
                     beta_x = BET
!
                  END IF IF_cont6
!
               ELSE IF( coord_order == 5 ) THEN
!
                  first  = i + kx_offset
                  second = j+1
                  third  = k
!
                  IF_cont7: IF( ensure_grid_continuity ) THEN
!
                     IF( partial_processing ) THEN
                        IF( i < LLoop_begin .OR. i > LLoop_end ) RETURN
                     END IF
!
                     IF( X(j) < X(jmin(i)) ) printout = .FALSE.
!
                     IF( kx_offset == 0 ) GO TO 1500
!
                     Fc = coeff(1)
                     beta_x = COS(phi_x)
                     IF( ABS(beta_x) < 1.0d-9 ) beta_x = 0.0d0
!
                     ijk = j + (k - 1) * NX + (i - 1) * NXNY
                     DO n1 = next2, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(j) + Fc * D2d(1)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRp(j) + Fc * D2d(1)
                           END IF
                           next2 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + 1 + kx_offset * NXNY
                     DO n1 = next2, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(j+1) + Fc * D2d(2)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRm(j+1) + Fc * D2d(2)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  ELSE
!
                     beta_x = BET
!
                  END IF IF_cont7
!
               ELSE IF( coord_order == 1 ) THEN
!
                  first  = i
                  second = j+1
                  third  = k + ky_offset
!
                  IF_cont8: IF( ensure_grid_continuity ) THEN
!
                     IF( partial_processing ) THEN
                        IF( i < LLoop_begin .OR. i > LLoop_end ) RETURN
                     END IF
!
                     IF( Y(j) < Y(jmin(i)) ) printout = .FALSE.
!
                     IF( ky_offset == 0 ) THEN
                        IF( grid_type == 2 ) Fr = X(i)
                        GO TO 1500
                     END IF
!
                     Fc = coeff(2)
                     beta_y = COS(phi_y)
                     IF( ABS(beta_y) < 1.0d-9 ) beta_y = 0.0d0
!
                     ijk = i + (j - 1) * NX + (k - 1) * NXNY
                     DO n1 = next2, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(j) + Fc * D2d(1)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(i-1) + Xr(i) ) * ( 5.0d-1 * DY(j) - Fc * D2d(1) )
                              Fr = X(i)
                           END IF
                           next2 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + NX + ky_offset * NXNY
                     DO n1 = next2, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(j+1) + Fc * D2d(2)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(i-1) + Xr(i) ) * ( 5.0d-1 * DY(j+1) - Fc * D2d(2)  )
                              Fr = X(i)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  END IF IF_cont8
!
               ELSE IF( coord_order == 6 ) THEN
!
                  first  = i + ky_offset
                  second = j+1
                  third  = k
!
                  IF_cont9: IF( ensure_grid_continuity ) THEN
!
                     IF( partial_processing ) THEN
                        IF( i < LLoop_begin .OR. i > LLoop_end ) RETURN
                     END IF
!
                     IF( Y(j) < Y(jmin(i)) ) printout = .FALSE.
!
                     IF( ky_offset == 0 ) THEN
                        IF( grid_type == 2 ) Fr = X(k)
                        GO TO 1500
                     END IF
!
                     Fc = coeff(2)
                     beta_y = COS(phi_y)
                     IF( ABS(beta_y) < 1.0d-9 ) beta_y = 0.0d0
!
                     ijk = k + (j - 1) * NX + (i - 1) * NXNY
                     DO n1 = next2, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(j) + Fc * D2d(1)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(k-1) + Xr(k) ) * ( 5.0d-1 * DY(j) - Fc * D2d(1) )
                              Fr = X(k)
                           END IF
                           next2 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = k + NX + ky_offset * NXNY
                     DO n1 = next2, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(j+1) + Fc * D2d(2)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(k-1) + Xr(k) ) * ( 5.0d-1 * DY(j+1) - Fc * D2d(2) )
                              Fr = X(k)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  END IF IF_cont9
!
               ELSE
!
                  Fr = 1.0d0
                  Fc = 1.0d0
                  beta_x = BET
                  beta_z = BETA
                  IF( ABS(beta_z) < 1.0d-9 ) beta_z = 0.0d0
                  first  = i
                  second = j+1
                  third  = k
!
               END IF
!
 1500          IF( outside_I_Range_limits .OR. (.NOT. printout) ) GO TO 2000
!
               NumConx = NumConx + 1
!
               IF(ElemNameLength == 8) THEN
                  ELNAM18 = Name_of_8Character_Element( NLXYZ(1), NLXYZ(2), NLXYZ(3), first, second, third )
               ELSE
                  ELNAM1 = Name_of_5Character_Element( NLXYZ(1), NLXYZ(2), NLXYZ(3), first, second, third )
               END IF
!
               IF( emissivity /= emissivity_p ) emissivity = 0.0d0
!
            END IF
!
            IF( ABS(beta_x) < 1.0d-10 ) beta_x = 0.0d0
            IF( ABS(beta_y) < 1.0d-10 ) beta_y = 0.0d0
            IF( ABS(beta_z) < 1.0d-10 ) beta_z = 0.0d0
!
            Fcr = Fc*Fr
!
            IF_Format4: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
!
               IF(ElemNameLength == 8) THEN
                  IF((coord_order == 3) .OR. (coord_order == 5)) WRITE(output_file_unit(2),5520) ELNAME8,ELNAM18,A13,Fcr*D2d(1),Fcr*D2d(2),NA(2),beta_x,emissivity
                  IF((coord_order == 1) .OR. (coord_order == 6)) WRITE(output_file_unit(2),5521) ELNAME8,ELNAM18,A13,Fcr*D2d(1),Fcr*D2d(2),NA(3),beta_y,emissivity
                  IF((coord_order == 2) .OR. (coord_order == 4)) WRITE(output_file_unit(2),5521) ELNAME8,ELNAM18,A13,Fcr*D2d(1),Fcr*D2d(2),NA(4),beta_z,emissivity
               ELSE
                  IF((coord_order == 3) .OR. (coord_order == 5)) WRITE(output_file_unit(2),5530) ELNAME,ELNAM1,A13,Fcr*D2d(1),Fcr*D2d(2),NA(2),beta_x,emissivity
                  IF((coord_order == 1) .OR. (coord_order == 6)) WRITE(output_file_unit(2),5531) ELNAME,ELNAM1,A13,Fcr*D2d(1),Fcr*D2d(2),NA(3),beta_y,emissivity
                  IF((coord_order == 2) .OR. (coord_order == 4)) WRITE(output_file_unit(2),5531) ELNAME,ELNAM1,A13,Fcr*D2d(1),Fcr*D2d(2),NA(4),beta_z,emissivity
               END IF
 !
            ELSE
!
               IF(ElemNameLength == 8) THEN
                  IF((coord_order == 3) .OR. (coord_order == 5)) WRITE(UNIT = output_file_unit(2), FMT = 5011, ADVANCE = 'NO') ELNAME8,ELNAM18,NA(2),Fcr*D2d(1),Fcr*D2d(2),A13,beta_x
                  IF((coord_order == 1) .OR. (coord_order == 6)) WRITE(UNIT = output_file_unit(2), FMT = 5011, ADVANCE = 'NO') ELNAME8,ELNAM18,NA(3),Fcr*D2d(1),Fcr*D2d(2),A13,beta_y
                  IF((coord_order == 2) .OR. (coord_order == 4)) WRITE(UNIT = output_file_unit(2), FMT = 5011, ADVANCE = 'NO') ELNAME8,ELNAM18,NA(4),Fcr*D2d(1),Fcr*D2d(2),A13,beta_z
                  IF( emissivity > 0.0e0 ) THEN
                     WRITE(UNIT = output_file_unit(2), FMT = '(ES10.4E1)' ) emissivity
                  ELSE
                     WRITE(UNIT = output_file_unit(2), FMT = '(A)' ) '          '
                  END IF
               ELSE
                  IF((coord_order == 3) .OR. (coord_order == 5)) WRITE(UNIT = output_file_unit(2), FMT = 11, ADVANCE = 'NO') ELNAME,ELNAM1,NA(2),Fcr*D2d(1),Fcr*D2d(2),A13,beta_x
                  IF((coord_order == 1) .OR. (coord_order == 6)) WRITE(UNIT = output_file_unit(2), FMT = 11, ADVANCE = 'NO') ELNAME,ELNAM1,NA(3),Fcr*D2d(1),Fcr*D2d(2),A13,beta_y
                  IF((coord_order == 2) .OR. (coord_order == 4)) WRITE(UNIT = output_file_unit(2), FMT = 11, ADVANCE = 'NO') ELNAME,ELNAM1,NA(4),Fcr*D2d(1),Fcr*D2d(2),A13,beta_z
                  IF( emissivity > 0.0e0 ) THEN
                     WRITE(UNIT = output_file_unit(2), FMT = '(ES10.4E1)' ) emissivity
                  ELSE
                     WRITE(UNIT = output_file_unit(2), FMT = '(A)' ) '          '
                  END IF
               END IF
!
            END IF IF_Format4
!
         END IF IF_Jcheck
!
! >>>>>>
! >>>>>> Checking the (i,j,k+1)-element
! >>>>>>
!
 2000    IF_Kcheck: IF(K < NLXYZ(1)) THEN
!
            printout = .TRUE.
!
            IF( .NOT. accepted_ijk_ext(3) ) THEN
               RETURN
            ELSE
!
               Fc = 1.0d0
               Fr = 1.0d0
               beta_x = 0.0d0
               beta_y = 0.0d0
               beta_z = BETA
               IF( ABS(beta_z) < 1.0d-9 ) beta_z = 0.0d0
!
               IF( coord_order == 4 ) THEN
!
                  first  = i
                  second = j + kx_offset
                  third  = k + 1
!
                  IF_contA: IF( ensure_grid_continuity ) THEN
!
                     IF( partial_processing ) THEN
                        IF( i < LLoop_begin .OR. i > LLoop_end ) RETURN
                     END IF
!
                     IF( X(k) < X(kmin(i,j)) ) printout = .FALSE.
!
                     IF( kx_offset == 0 ) GO TO 2500
!
                     Fc = coeff(1)
                     beta_x = COS(phi_x)
                     IF( ABS(beta_x) < 1.0d-9 ) beta_x = 0.0d0
!
                     ijk = k + (i - 1) * NX + (j - 1) * NXNY
                     DO n1 = next3, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(k) + Fc * D3d(1)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRp(k) + Fc * D3d(1)
                           END IF
                           next3 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + 1 + kx_offset * NXNY
                     DO n1 = next3, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(k+1) + Fc * D3d(2)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRm(k+1) + Fc * D3d(2)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  ELSE
!
                     beta_x = BET
!
                  END IF IF_contA
!
               ELSE IF( coord_order == 6 ) THEN
!
                  first  = i + kx_offset
                  second = j
                  third  = k + 1
!
                  IF_contB: IF( ensure_grid_continuity ) THEN
!
                     IF( partial_processing ) THEN
                        IF( i < LLoop_begin .OR. i > LLoop_end ) RETURN
                     END IF
!
                     IF( X(k) < X(kmin(i,j)) ) printout = .FALSE.
!
                     IF( kx_offset == 0 ) GO TO 2500
!
                     Fc = coeff(1)
                     beta_x = COS(phi_x)
                     IF( ABS(beta_x) < 1.0d-9 ) beta_x = 0.0d0
!
!
                     ijk = k + (j - 1) * NX + (i - 1) * NXNY
                     DO n1 = next3, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(k) + Fc * D3d(1)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRp(k) + Fc * D3d(1)
                           END IF
                           next3 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + 1 + kx_offset  * NXNY
                     DO n1 = next3, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DX_adj(n1) = DX_adj(n1) - 5.0d-1 * DX(k+1) + Fc * D3d(2)
                           ELSE
                              DX_adj(n1) = DX_adj(n1) - DRm(k+1) + Fc * D3d(2)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  ELSE
!
                     beta_x = BET
!
                  END IF IF_contB
!
               ELSE IF( coord_order == 2 ) THEN
!
                  first  = i
                  second = j + ky_offset
                  third  = k + 1
!
                  IF_contC: IF( ensure_grid_continuity ) THEN
!
                     IF( partial_processing ) THEN
                        IF( i < LLoop_begin .OR. i > LLoop_end ) RETURN
                     END IF
!
                     IF( Y(k) < Y(kmin(i,j)) ) printout = .FALSE.
!
                     IF( ky_offset == 0 ) THEN
                        IF( grid_type == 2 ) Fr = X(i)
                        GO TO 2500
                     END IF
!
                     Fc = coeff(2)
                     beta_y = COS(phi_y)
                     IF( ABS(beta_y) < 1.0d-9 ) beta_y = 0.0d0
!
                     ijk = i + (k - 1) * NX + (j - 1) * NXNY
                     DO n1 = next3, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(k) + Fc * D3d(1)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(i-1) + Xr(i) ) * ( 5.0d-1 * DY(k) - Fc * D3d(1) )
                              Fr = X(i)
                           END IF
                           next3 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + NX + ky_offset * NXNY
                     DO n1 = next3, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(k+1) + Fc * D3d(2)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(i-1) + Xr(i) ) * ( 5.0d-1 * DY(k+1) - Fc * D3d(2) )
                              Fr = X(i)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  END IF IF_contC
!
               ELSE IF( coord_order == 5 ) THEN
!
                  first  = i + ky_offset
                  second = j
                  third  = k + 1
!
                  IF_contD: IF( ensure_grid_continuity ) THEN
!
                     IF( partial_processing ) THEN
                        IF( i < LLoop_begin .OR. i > LLoop_end ) RETURN
                     END IF
!
                     IF( Y(k) < Y(kmin(i,j)) ) printout = .FALSE.
!
                     IF( ky_offset == 0 ) THEN
                        IF( grid_type == 2 ) Fr = X(j)
                        GO TO 2500
                     END IF
!
                     Fc = coeff(2)
                     beta_y = COS(phi_y)
                     IF( ABS(beta_y) < 1.0d-9 ) beta_y = 0.0d0
!
                     ijk = j + (k - 1) * NX + (i - 1) * NXNY
                     DO n1 = next3, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(k) + Fc * D3d(1)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(j-1) + Xr(j) ) * ( 5.0d-1 * DY(k) - Fc * D3d(1) )
                              Fr = X(j)
                           END IF
                           next3 = n1 + 1
                           EXIT
                        END IF
                     END DO
!
                     ijk = ijk + NX + ky_offset  * NXNY
                     DO n1 = next3, NumTotalElem
                        IF( ijk_location(n1) == ijk ) THEN
                           IF( grid_type == 1 ) THEN
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * DY(k+1) + Fc * D3d(2)
                           ELSE
                              DY_adj(n1) = DY_adj(n1) - 5.0d-1 * ( Xr(j-1) + Xr(j) ) * ( 5.0d-1 * DY(k+1) - Fc * D3d(2) )
                              Fr = X(j)
                           END IF
                           EXIT
                        END IF
                     END DO
!
                  END IF IF_contD
!
               ELSE
!
                  Fc = 1.0d0
                  Fr = 1.0d0
                  beta_x = BET
                  beta_z = BETA
                  IF( ABS(beta_z) < 1.0d-9 ) beta_z = 0.0d0
                  first  = i
                  second = j
                  third  = k+1
!
               END IF
!
 2500          IF( outside_I_Range_limits .OR. (.NOT. printout) ) RETURN
!
               NumConx = NumConx + 1
!
               IF(ElemNameLength == 8) THEN
                  ELNAM18 = Name_of_8Character_Element( NLXYZ(1), NLXYZ(2), NLXYZ(3), first, second, third )
               ELSE
                  ELNAM1 = Name_of_5Character_Element(NLXYZ(1), NLXYZ(2), NLXYZ(3), first, second, third )
               END IF
!
               IF( emissivity /= emissivity_p ) emissivity = 0.0d0
!
            END IF
!
            IF( ABS(beta_x) < 1.0d-10 ) beta_x = 0.0d0
            IF( ABS(beta_y) < 1.0d-10 ) beta_y = 0.0d0
            IF( ABS(beta_z) < 1.0d-10 ) beta_z = 0.0d0
!
            Fcr = Fc*Fr
!
            IF_Format5: IF(FormatType == 'New' .OR. FormatType == 'new' .OR. FormatType == 'NEW') THEN
!
               IF(ElemNameLength == 8) THEN
                  IF((coord_order == 4) .OR. (coord_order == 6)) WRITE(output_file_unit(2),5520) ELNAME8,ELNAM18,A23,Fcr*D3d(1),Fcr*D3d(2),NA(2),beta_x,emissivity
                  IF((coord_order == 2) .OR. (coord_order == 5)) WRITE(output_file_unit(2),5521) ELNAME8,ELNAM18,A23,Fcr*D3d(1),Fcr*D3d(2),NA(3),beta_y,emissivity
                  IF((coord_order == 1) .OR. (coord_order == 3)) WRITE(output_file_unit(2),5521) ELNAME8,ELNAM18,A23,Fcr*D3d(1),Fcr*D3d(2),NA(4),beta_z,emissivity
               ELSE
                  IF((coord_order == 4) .OR. (coord_order == 6)) WRITE(output_file_unit(2),5530) ELNAME,ELNAM1,A23,Fcr*D3d(1),Fcr*D3d(2),NA(2),beta_x,emissivity
                  IF((coord_order == 2) .OR. (coord_order == 5)) WRITE(output_file_unit(2),5531) ELNAME,ELNAM1,A23,Fcr*D3d(1),Fcr*D3d(2),NA(3),beta_y,emissivity
                  IF((coord_order == 1) .OR. (coord_order == 3)) WRITE(output_file_unit(2),5531) ELNAME,ELNAM1,A23,Fcr*D3d(1),Fcr*D3d(2),NA(4),beta_z,emissivity
               END IF
!
            ELSE
!
               IF(ElemNameLength == 8) THEN
                  IF((coord_order == 4) .OR. (coord_order == 6)) WRITE(UNIT = output_file_unit(2), FMT = 5011, ADVANCE = 'NO') ELNAME8,ELNAM18,NA(2),Fcr*D3d(1),Fcr*D3d(2),A23,beta_x
                  IF((coord_order == 2) .OR. (coord_order == 5)) WRITE(UNIT = output_file_unit(2), FMT = 5011, ADVANCE = 'NO') ELNAME8,ELNAM18,NA(3),Fcr*D3d(1),Fcr*D3d(2),A23,beta_y
                  IF((coord_order == 1) .OR. (coord_order == 3)) WRITE(UNIT = output_file_unit(2), FMT = 5011, ADVANCE = 'NO') ELNAME8,ELNAM18,NA(4),Fcr*D3d(1),Fcr*D3d(2),A23,beta_z
                  IF( emissivity > 0.0e0 ) THEN
                     WRITE(UNIT = output_file_unit(2), FMT = '(ES10.4E1)' ) emissivity
                  ELSE
                     WRITE(UNIT = output_file_unit(2), FMT = '(A)' ) '          '
                  END IF
               ELSE
                  IF((coord_order == 4) .OR. (coord_order == 6)) WRITE(UNIT = output_file_unit(2), FMT = 11, ADVANCE = 'NO') ELNAME,ELNAM1,NA(2),Fcr*D3d(1),Fcr*D3d(2),A23,beta_x
                  IF((coord_order == 2) .OR. (coord_order == 5)) WRITE(UNIT = output_file_unit(2), FMT = 11, ADVANCE = 'NO') ELNAME,ELNAM1,NA(3),Fcr*D3d(1),Fcr*D3d(2),A23,beta_y
                  IF((coord_order == 1) .OR. (coord_order == 3)) WRITE(UNIT = output_file_unit(2), FMT = 11, ADVANCE = 'NO') ELNAME,ELNAM1,NA(4),Fcr*D3d(1),Fcr*D3d(2),A23,beta_z
                  IF( emissivity > 0.0e0 ) THEN
                     WRITE(UNIT = output_file_unit(2), FMT = '(ES10.4E1)' ) emissivity
                  ELSE
                     WRITE(UNIT = output_file_unit(2), FMT = '(A)' ) '          '
                  END IF
               END IF
!
            END IF IF_Format5
!
         END IF IF_Kcheck
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Create_Connection_List 1.0 ............. 31 August    2015',6X,'Create the CONNECTION data block for TOUGH2/TOUGH+')
!
 5011 FORMAT(2(A8),13X,A1,4ES10.4E1)
   11 FORMAT(2(A5),19X,A1,4ES10.4E1)
!
 5530 FORMAT(A,A,'  &ConX  A=',ES17.10,', d1=',ES15.8,', d2=',ES15.8,', dir=',A1,', emsv=',ES15.8,' /')
 5531 FORMAT(A,A,'  &Conx  A=',ES17.10,', d1=',ES15.8,', d2=',ES15.8,', dir=',A1,', beta=',ES15.8,', emsv=',ES15.8,' /')
!
 5520 FORMAT(A,A,'  &ConX  A=',ES17.10,', d1=',ES15.8,', d2=',ES15.8,', dir=',A1,', emsv=',ES15.8,' /')
 5521 FORMAT(A,A,'  &Conx  A=',ES17.10,', d1=',ES15.8,', d2=',ES15.8,', dir=',A1,', beta=',ES15.8,', emsv=',ES15.8,' /')
!
 6850 FORMAT(//,20('ERROR-'),//,   &
     &       T5,'Memory allocation to arrays <',A,'> in subroutine <',A,'> was unsuccessful',//,   &
     &       T5,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,20('ERROR-'))
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of <Create_Connection_List>
!
         RETURN
!
      END SUBROUTINE Create_Connection_List
!
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      SUBROUTINE Adjust_Element_Volumes
!
!
         IMPLICIT NONE
!
! ----------
! ...... Double precision variables
! ----------
!
         REAL(KIND = 8) :: volume
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER :: unit_elem, unit_conx
!
! ----------
! ...... Character variables
! ----------
!
         CHARACTER(LEN =  5) :: header
         CHARACTER(LEN = 20) :: part1
         CHARACTER(LEN = 54) :: part2
         CHARACTER(LEN = 80) :: connection_line
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Beginning of <Adjust_Element_Volumes>
!
         WRITE( UNIT = *, FMT = 6000 )
!
!
! >>>>>>
! >>>>>> Compute the modified volume
! >>>>>>
!
!
         REWIND ( UNIT = output_file_unit(1) )
!
         READ   ( UNIT = output_file_unit(1), FMT = '(A5)' ) header
!
         IF( level_of_grid_generation == 'F' .OR. level_of_grid_generation == 'f' ) THEN
            unit_elem = MESH_Unit
            unit_conx = MESH_Unit
         ELSE
            unit_elem = elem_file_unit
            unit_conx = conx_file_unit
         END IF
!
         REWIND ( UNIT = unit_elem )
         WRITE  ( UNIT = unit_elem, FMT = '(A8)' ) 'ELEMENTS'
!
         DO n = NumElem_Bound1 + 1, NumTotalElem - NumElem_Bound2
!
            READ ( UNIT = output_file_unit(1), FMT = 5000 ) part1, volume, part2, itrue, jtrue, ktrue
!
            volume = DX_adj(n) * DY_adj(n) * DZ(ktrue)
!
! ......... Modify the ELEMENT data block in the appropriate file
!
            WRITE ( UNIT = unit_elem, FMT = 5000 ) part1, volume, part2
!
         END DO
!
         WRITE ( UNIT = unit_elem, FMT = '(A10)' ) '          '
!
         CLOSE ( UNIT = output_file_unit(1), STATUS = 'DELETE' )
!
         REWIND( UNIT = output_file_unit(2) )
!
         READ  ( UNIT = output_file_unit(2), FMT = '(A5)' ) header
!
! ...... Adding the CONNECTIONS block to the the appropriate file (MESH or CONNECTIONS)
!
         WRITE ( UNIT = unit_conx, FMT = '(A11)' ) 'CONNECTIONS'
!
! ...... Adding the CONNECTIONS block to the MESH file
!
         DO n = 1, NumConx
            READ ( UNIT = output_file_unit(2), FMT = '(A80)' ) connection_line
            WRITE( UNIT = unit_conx,           FMT = '(A80)' ) connection_line
         END DO
!
         WRITE ( UNIT = unit_conx, FMT = '(A20)' ) '                    '
!
         CLOSE ( UNIT = output_file_unit(2), STATUS = 'DELETE' )
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
!
 6000 FORMAT(/,'Adjust_Element_Volumes 1.0 ............. 29 September 2015',6X,'Modify the volumes for a system with no areal grid discontinuity')
!
 5000 FORMAT(A20,ES10.4E2,A54,3(1x,i5))
!
!
!
         RETURN
!
      END SUBROUTINE Adjust_Element_Volumes
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
      END SUBROUTINE GXYZ
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE PCAR(NX,NY,NZ,NXA,NYA,NZA,coord_order)
!
         USE MeshMaker_Data
!
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! -------
! ... Integer variables
! -------
!
      INTEGER, INTENT(IN) :: coord_order
!
      INTEGER :: ICALL = 0
!
! -------
! ... CHARACTER arrays
! -------
!
      CHARACTER(LEN = 1), DIMENSION(6), PARAMETER :: L1 = (/ 'K','J','K','I','J','I' /), &
     &                                               L2 = (/ 'J','K','I','K','I','J' /), &
     &                                               L3 = (/ 'I','I','J','J','K','K' /)
!
      CHARACTER(LEN = 2), DIMENSION(6), PARAMETER :: M1 = (/ 'NZ','NY','NZ','NX','NY','NX' /), &
     &                                               M2 = (/ 'NY','NZ','NX','NZ','NX','NY' /), &
     &                                               M3 = (/ 'NX','NX','NY','NY','NZ','NZ' /)
!
! -------
! ... Saving variables
! -------
!
      SAVE ICALL
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of PCAR
!
!
!.... Local function to find location of element name with X-Index I, Y-Index J, And Z-Index K
!
      NLOC(I,J,K)=(I-1)*NYA*NZA+(J-1)*NZA+K
!
      ICALL=ICALL+1
      IF(ICALL == 1) WRITE(UNIT = *, FMT = 6000)
!
! 899 FORMAT(6X,'PCAR     1.0      25 MARCH     1991',6X,
! 899 FORMAT('PCAR        1.1   19 September    2000',6X,'MAKE STRUCTURED PRINTOUT OF CARTESIAN MESH')
 6000 FORMAT(/,'PCAR 1.5 ................................ 2 August    2014',6X,'MAKE STRUCTURED PRINTOUT OF CARTESIAN MESH')
!
!
!
      WRITE( UNIT = *, FMT = 31) NXA,NYA,NZA
!
   31 FORMAT(/,' ',131('*'),   &
     &       /,' *',20X,'CARTESIAN MESH WITH NX*NY*NZ = ',I6,' *',I6,' *',I6,'  GRID BLOCKS',43X,'*',   &
     &       /,' ',131('*'))
!
      N1 = NZ
      N2 = NY
      N3 = NX
!
      WRITE( UNIT = *, FMT = 18) L1(coord_order), L1(coord_order), M1(coord_order), N1
   18 FORMAT(' *',129X,'*',/,' *',20X,'THE MESH WILL BE PRINTED AS SLICES FOR ',A1,' = 1 TO ',A1,' = ',A2,' =',I6,47X,'*')
!
      WRITE( UNIT = *, FMT = 19) L2(coord_order), L2(coord_order), M2(coord_order), N2
   19 FORMAT(' *',129X,'*',/,' *',20X,'IN EACH MESH SLICE, ROWS WILL GO FROM  ',A1,' = 1 TO ',A1,' = ',A2,' =',I6,47X,'*')
!
      WRITE( UNIT = *, FMT = 20) L3(coord_order), L3(coord_order),M3(coord_order), N3
   20 FORMAT(' *',129X,'*',/,' *',20X,'IN EACH ROW, COLUMNS WILL GO FROM      ',A1,' = 1 TO ',A1,' = ',A2,' =',I6,47X,'*',/,' *',129X,'*')
!
! ... GJM: Begin modification for 8-character elements
!
      IF(ElemNameLength == 8) THEN
!
         DO_K1: DO K = 1, N1
!
            WRITE( UNIT = *, FMT = 3 )   L1(coord_order), K
            WRITE( UNIT = *, FMT = 5004) L3(coord_order)
!
            SELECT CASE(coord_order)
            CASE(1)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 5006 ) L2(coord_order),J,(name(NLOC(I,J,K)),I=1,N3)
               END DO
            CASE(2)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 5006 ) L2(coord_order),J,(name(NLOC(I,K,J)),I=1,N3)
               END DO
            CASE(3)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 5006 ) L2(coord_order),J,(name(NLOC(J,I,K)),I=1,N3)
               END DO
            CASE(4)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 5006 ) L2(coord_order),J,(name(NLOC(K,I,J)),I=1,N3)
               END DO
            CASE(5)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 5006 ) L2(coord_order),J,(name(NLOC(J,K,I)),I=1,N3)
               END DO
            CASE(6)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 5006 ) L2(coord_order),J,(name(NLOC(K,J,I)),I=1,N3)
               END DO
            END SELECT
!
         END DO DO_K1
!
      ELSE
!
!
         DO_K2: DO K = 1, N1
!
            WRITE( UNIT = *, FMT = 3 ) L1(coord_order), K
            WRITE( UNIT = *, FMT = 4 ) L3(coord_order)
!
            SELECT CASE(coord_order)
            CASE(1)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 6 ) L2(coord_order),J,(name(NLOC(I,J,K))(1:5),I=1,N3)
               END DO
            CASE(2)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 6 ) L2(coord_order),J,(name(NLOC(I,K,J))(1:5),I=1,N3)
               END DO
            CASE(3)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 6 ) L2(coord_order),J,(name(NLOC(J,I,K))(1:5),I=1,N3)
               END DO
            CASE(4)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 6 ) L2(coord_order),J,(name(NLOC(K,I,J))(1:5),I=1,N3)
               END DO
            CASE(5)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 6 ) L2(coord_order),J,(name(NLOC(J,K,I))(1:5),I=1,N3)
               END DO
            CASE(6)
               DO j = 1, N2
                  WRITE( UNIT = *, FMT = 6 ) L2(coord_order),J,(name(NLOC(K,J,I))(1:5),I=1,N3)
               END DO
            END SELECT
!
         END DO DO_K2
!
      END IF
!
!    MBK: Modifications end here (Oct 4, 2004)
!
    3 FORMAT(' ',131('*'),//,' SLICE WITH ',A1,' =',I6/)
!
 5004 FORMAT('   COLUMN ',A1,' =      1        2        3        4        5        6        7        8        9       10       11       12',/,' LAYERS')
    4 FORMAT('   COLUMN ',A1,' =    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20',/,' ROWS')
!
 5006 FORMAT('  ',A1,' = ',I6,2X,12(1X,A8),/,(14X,12(1X,A8)),/,(14X,12(1X,A8)),/,(14X,12(1X,A8)))
    6 FORMAT('  ',A1,' = ',I6,2X,20(1X,A5)/(14X,20(1X,A5)))
!
! - GJM: End modification for 8-character elements
!
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of PCAR
!
!
      RETURN
      END SUBROUTINE PCAR
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE MINC
!
         USE MeshMaker_Data
!
!-----THIS IS A MODIFIED AND ENHANCED VERSION OF PROGRAM MeshMaker (extension of GMIC)
!
!     Initially:
!     GMINC IMPLEMENTS THE METHOD OF
!                  MULTIPLE INTERACTING CONTINUA (MINC)
!     AS DEVELOPED BY PRUESS AND NARASIMHAN.
!
!     REFERENCES:
!
!     (1) K. PRUESS AND T.N. NARASIMHAN, A PRACTICAL METHOD FOR
!         MODELING FLUID AND HEAT FLOW IN FRACTURED POROUS MEDIA,
!         PAPER SPE-10509, PRESENTED AT THE SIXTH SPE-SYMPOSIUM ON
!         RESERVOIR SIMULATION, NEW ORLEANS, LA; (FEBRUARY 1982);
!         ALSO: SOCIETY OF PETROLEUM ENGINEERS JOURNAL, VOL. 25, NO.1,
!               PP. 14-26, 1985.
!
!     (2) K. PRUESS AND T.N. NARASIMHAN, ON FLUID RESERVES AND THE
!         PRODUCTION OF SUPERHEATED STEAM FROM FRACTURED, VAPOR-
!         DOMINATED GEOTHERMAL RESERVOIRS, J. GEOPHYS. RES. 87 (B11),
!         9329-9339, 1982.
!
!     (3) K. PRUESS AND K. KARASAKI, PROXIMITY FUNCTIONS FOR MODELING
!         FLUID AND HEAT FLOW IN RESERVOIRS WITH STOCHASTIC FRACTURE
!         DISTRIBUTIONS, PAPER PRESENTED AT EIGTH STANFORD WORKSHOP
!         ON GEOTHERMAL RESERVOIR ENGINEERING, STANFORD, CA.
!         (DECEMBER 1982).
!
!     (4) K. PRUESS, GMINC - A MESH GENERATOR FOR FLOW SIMULATIONS IN
!         FRACTURED RESERVOIRS, LAWRENCE BERKELEY LABORATORY REPORT
!         LBL-15227, 1983.
!
!-----   GMINC GENERATES ONE-, TWO-, OR THREE-DIMENSIONAL MESHES
!        FOR FLOW SIMULATIONS IN FRACTURED POROUS MEDIA.
!
!        Current version developed by G. Moridis in February/March 2015
!
!        New capabilities:
!    (a) Namelist-based data entry
!    (b) Definition of heterogeneous regions/subdomains
!    (c) Ability to define interacting fractured and unfractured media
!
!
!***********************************************************************
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      SAVE ICALL
      DATA ICALL/0/
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of MINC
!
!
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(UNIT = *, FMT = 6000)
 6000 FORMAT(/,'MINC 1.0 ............................... 22 JANUARY   1990',6X,'EXECUTIVE ROUTINE FOR MAKING A "SECONDARY" FRACTURED-POROUS MEDIUM MESH')
!
!-----READ DATA ON MINC-PARTITIONING.
!
      CALL PART
!
! ... CALL TO GEOM INSERTED 11-19-84, TO OBTAIN ALL VOL(M) BEFORE PROCESSING CONNECTIONS.
!
      CALL GEOM
!
!-----READ ELEMENT DATA FROM FILE *MESH* AND PROCESS SEQUENTIALLY.
!
      CALL MINCME
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of MINC
!
!
      RETURN
!
      END SUBROUTINE MINC
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE MINC_Parameters
!
!
         SAVE
!
! ----------
! ...... Double precision arrays
! ----------
!
         REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: VOL, Am, Dm
!
         REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: PAR
!
! ----------
! ...... Integer variables
! ----------
!
         INTEGER :: J, L, NVOL
         INTEGER :: num_media, num_FracMedia
!
! ----------
! ...... Character variables and arrays
! ----------
!
         CHARACTER(LEN = 14) :: matrix_to_matrix_flow
!
! ----------
! ...... Logical variables
! ----------
!
         LOGICAL :: uniform_matrix_properties, fracture_is_1st_continuum
!
         LOGICAL, ALLOCATABLE, DIMENSION(:) :: fractured_medium
!
!
      END MODULE MINC_Parameters
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE PART
!
!
         USE MeshMaker_Data
         USE Grid_Generation_Parameters
         USE MINC_Parameters
!
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      REAL(KIND = 8), DIMENSION(7) :: fracture_DistParameters
!
      INTEGER :: number_of_MINCs, number_of_specified_VolFractions
      INTEGER :: fractured_medium_number, nfm
!
      CHARACTER(LEN=  3) :: indicator
      CHARACTER(LEN=  5) :: proximity_function_type
      CHARACTER(LEN=100) :: VolFraction_data_format
!
      CHARACTER(LEN=1), PARAMETER, DIMENSION(6) :: MM_flow = (/ 'V', 'v', 'A', 'a', 'F', 'f' /)
!
      LOGICAL :: First_call = .TRUE.
!
      SAVE First_call
!
! -------
! ... Namelists
! -------
!
      NAMELIST/ General_MINC_Data / number_of_media,                  number_of_fractured_media,  number_of_MINCs,       &
                                    number_of_specified_VolFractions, fracture_is_1st_continuum,  matrix_to_matrix_flow
!
      NAMELIST/ Medium_MINC_Data / fractured_medium_number,   proximity_function_type, fracture_DistParameters,  &
     &                             uniform_matrix_properties, VolFraction_data_format
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of PART
!
!
      IF (First_call) THEN
         First_call = .FALSE.
         WRITE( UNIT = *, FMT = 6000 )
      END IF
!
! 899 FORMAT('PART              1.0   22 JANUARY   1990',6X,'READ SPECIFICATIONS OF MINC-PARTITIONING FROM FILE "INPUT"')
 6000 FORMAT(/,'PART 2.0 ..............................  25 FEBRUARY  2015',6X,'READ SPECIFICATIONS OF MINC-PARTITIONING FROM FILE "INPUT"')
!
! ... Initialization
!
      fracture_is_1st_continuum = .TRUE.


!
      READ( UNIT = *, NML = General_MINC_Data, IOSTAT = ier )
! !
      IF(ier /= 0) THEN
         WRITE(UNIT = *, FMT = 6001) '<General_MINC_Data>'
         STOP
      END IF
!

      num_media     = number_of_media
      num_FracMedia = number_of_fractured_media
!
      J    = number_of_MINCs
      NVOL = number_of_specified_VolFractions
!
      ALLOCATE( VOL(num_media, J), Am(num_media, J), Dm(num_media, J), PAR(num_media, 7) )
      ALLOCATE( fractured_medium(num_media), fractured(0:Max_NumElem) )
!
      fractured_medium = .FALSE.   ! ... Whole array operation
!
      DO_BigLoop: DO ixx = 1, num_FracMedia
!
         fractured_medium_number   = 0
         proximity_function_type   = ' '
         uniform_matrix_properties = .TRUE.
         fracture_DistParameters   = 0.0d0
!
         READ( UNIT = *, NML = Medium_MINC_Data, IOSTAT = ier )
!
         IF(ier /= 0) THEN
            WRITE(UNIT = *, FMT = 6001) '<Medium_MINC_Data>'
            STOP
         END IF
!
 6001 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure "PART": There is a problem reading the namelist ',A,//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
         nfm = fractured_medium_number
         IF( nfm == 0 .OR. nfm > num_media ) THEN
            WRITE(*, FMT = 6002) nfm
            STOP
         END IF
!
         fractured_medium(nfm) = .TRUE.
!
6002 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure "PART", Namelist <Medium_MINC_Data>:',&
     &       T10,'!!!!!!! The variable < fractured_medium_number> =',i3,' is outside the possible range',/, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))

!
         IF( ALL(TYPOS(1:6) /= proximity_function_type) ) THEN
           WRITE( UNIT = *, FMT = 6004 ) proximity_function_type
            STOP
         ELSE
            DO L = 1,6
               IF( proximity_function_type == TYPOS(L) ) EXIT
            END DO
         END IF
!
6004 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure "PART", Namelist <Medium_MINC_Data>:',&
     &       T10,'!!!!!!! HAVE READ UNKNOWN PROXIMITY FUNCTION IDENTIFIER  "',A5,'"',/, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))

!
         PAR(nfm,1:7) = fracture_DistParameters(1:7)
!
         IF( TRIM(ADJUSTL(VolFraction_data_format)) == '(*)' .OR. TRIM(ADJUSTL(VolFraction_data_format)) == '*' ) THEN
             IF( fracture_is_1st_continuum ) THEN
                READ( UNIT = *, FMT = * ) ( Vol(nfm, i), i = 1,NVOL )
             ELSE
                READ( UNIT = *, FMT = * ) ( Vol(nfm, J+1-i), i = 1,NVOL )
             END IF
         ELSE
             IF( fracture_is_1st_continuum ) THEN
                READ( UNIT = *, FMT = VolFraction_data_format ) ( Vol(nfm, i), i = 1,NVOL )
             ELSE
                READ( UNIT = *, FMT = VolFraction_data_format ) ( Vol(nfm, J+1-i), i = 1,NVOL )
             END IF
         END IF
!
!----- END OF MINC-DATA ------------------------------------------------
!
         IF( (L == 2 .OR. L == 3) .AND. PAR(nfm,2) == 0.0e0 ) PAR(nfm,2) = PAR(nfm,1)
         IF( (L == 3) .AND. PAR(nfm,3) == 0.0e0 ) PAR(nfm,3) = PAR(nfm,2)
!
!-----IDENTIFY CHOICE MADE FOR GLOBAL MATRIX-MATRIX FLOW.
!
         IF( ALL( MM_flow /= matrix_to_matrix_flow(1:1) ) ) matrix_to_matrix_flow(1:1) = ' '
!
         WRITE( UNIT = *, FMT = 5010 ) matrix_to_matrix_flow
 5010 FORMAT(/' CHOICE OF MATRIX-MATRIX FLOW HANDLING: "',A5,'"')
!
         WRITE( UNIT = *, FMT = 5011 )
 5011 FORMAT(/'                       THE OPTIONS ARE: "     "',   &
     &        ' (DEFAULT):       NO GLOBAL MATRIX-MATRIX FLOW; GLOBAL FLOW ONLY THROUGH FRACTURES'/   &
     &    40X,'"Vertical":       GLOBAL MATRIX-MATRIX FLOW IN VERTICAL DIRECTION ONLY'/   &
     &    40X,'"All_Directions": GLOBAL MATRIX-MATRIX FLOW IN ALL DIRECTIONS'/)
!
         IF(i == 2 .AND. J /= 2) WRITE( UNIT = *, FMT = 5012 )
 5012 FORMAT(' !!!!! WARNING !!!!! THE "ALL" OPTION SHOULD ONLY',   &
     &       ' BE USED WITH TWO CONTINUA, WHERE IT AMOUNTS TO A DUAL-PERMEABILITY TREATMENT')
!
      END DO DO_BigLoop
!
      READ( UNIT = *, FMT = '(A3)' ) indicator
!
      IF( indicator /= '<<<' ) THEN
         WRITE( UNIT = *, FMT = 6003 ) indicator
         STOP
      END IF
!
 6003 FORMAT(//,22('ERROR-'),//, &
     &       T10,'>>>  Procedure "PART", Datablock "MINC": The ending descriptor of datablock = "',A,'" is not the required "<<<'//, &
     &       T10,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',   &
     &       //,22('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of PART
!
!
      RETURN
      END SUBROUTINE PART
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE GEOM
!
!
         USE MeshMaker_Data
         USE Grid_Generation_Parameters
         USE MINC_Parameters
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Xm
!
      REAL(KIND=8), PARAMETER :: DELTA = 1.0d-8
!
      REAL(KIND=8) :: XRm
!
      INTEGER :: ICALL = 0
!
      SAVE ICALL
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of GEOM
!
!
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(*,899)
  899 FORMAT(/,'GEOM 1.0 ................................ 1 MAY       1991',6X,'CALCULATE GEOMETRY PARAMETERS OF SECONDARY (FRACTURED-POROUS) MESH')
!
      ALLOCATE( Xm(num_media, J) )
!
!-----COME HERE TO ASSIGN EQUAL VOLUMINA TO SUBDIVISIONS WHICH HAVE NOT
!     BEEN EXPLICITLY SPECIFIED-----
!
      CONTINUE
!
! ----------------
! ........ GJM - 9/19/2000
! ........ Begin correction
! ........ Avoid problem of stalling when location does not have the desired value
! ----------------
!
      DO_k: DO k=1,num_media
!
         IF( .NOT. fractured_medium(k) ) CYCLE
!
         IF(NVOL >= J) GO TO 3
!
         VEX = 0.0d0
         VF  = 0.0d0
!
         DO M=1,NVOL
            IF( fracture_is_1st_continuum ) THEN
               VEX = VEX + VOL(k,m)
            ELSE
               VEX = VEX + VOL(k,J+1-m)
            END IF
         END DO
!
! ----------------
! ........ GJM - 9/19/2000
! ........ End correction
! ----------------
!
!     VEX IS THE TOTAL EXPLICITLY ASSIGNED VOLUME FRACTION.
!
      IF(VEX  > 1.0d0) GO TO 10
      IF(VEX == 1.0d0)   GO TO 3
!
!-----VF IS THE VOLUME FRACTION FOR PARTITIONS WHICH ARE NOT EXPLICITLY ASSIGNED.
!
      VF    = (1.0d0-VEX)/FLOAT(J-NVOL)
      NVOL1 = NVOL+1
!
      DO 5 M=NVOL1,J
         IF( fracture_is_1st_continuum ) THEN
            VOL(k,m) = VF
         ELSE
            VOL(k,J+1-m) = VF
         END IF
    5 CONTINUE
      GO TO 3
!
   10 CONTINUE
!-----COME HERE IF EXPLICITLY ASSIGNED VOLUMINA EXCEED 100%-----
      WRITE(*,11) VEX
   11 FORMAT(' PROGRAM STOPS BECAUSE TOTAL VOLUME VEX = ',ES12.5,' > 100%  ---  NEED TO CORRECT INPUT DATA')
      STOP
!
    3 CONTINUE
!
!-----NOW FIND DISTANCES FROM FRACTURES WHICH CORRESPOND TO
!     DESIRED VOLUME FRACTIONS.
!     INDEXING STARTS AT THE OUTSIDE; I.E. *1* IS THE OUTERMOST
!     VOLUME ELEMENT, AND *J* IS THE INNERMOST ONE.
!
!-----INITIALIZE TOTAL VOLUME FRACTION.
      TVOL = 0.0d0
!
!-----FIRST INTERFACE WILL BE AT FRACTURE FACE.
      Xm(k,1) = 0.0d0
      Dm(k,1) = 0.0d0
      Am(k,1) = (1.0d0 - VOL(k,1)) * PROX(k,1.0d-10)/1.0d-10
!
!
!-----INITIALIZE SEARCH INTERVAL.
!
      XL  = 0.0d0
      XRm = VOL(k,2)/Am(k,1)
!
!
      DO 30 M=2,J
!
!-----COMPUTE TOTAL FRACTION OF MATRIX VOLUME.
         TVOL = TVOL+VOL(k,M)/(1.0d0-VOL(k,1))
!
         IF(M == J) TVOL = 1.0d0-1.0d-9
!
         CALL INVER( k, TVOL, XMID, XL, XRm )
!
         Xm(k,M) = XMID
!
         XMD     = XMID*DELTA
         Am(k,M) = (1.0d0-VOL(k,1)) * (PROX(k,XMID+XMD)-PROX(k,XMID-XMD)) / (2.0d0*XMD)
!
         Dm(k,M) = (Xm(k,M) - Xm(k,M-1)) / 2.0d0
!
!-----PUT LEFT END OF NEXT ITERATION INTERVAL AT PRESENT X.
         XL = XMID
!
   30 CONTINUE
!
!
!-----COME HERE TO COMPUTE A QUASI-STEADY VALUE FOR INNERMOST
!     NODAL DISTANCE.
!
      GOTO (41,42,43,44,45,46,47,48,49,50),L
!
   41 CONTINUE
!
!----- ONE-D CASE.
!
      Dm(k,J) = (PAR(k,1)-2.0d0*Xm(k,J-1))/6.0d0
      GO TO 40
!
   42 CONTINUE
!
!----- TWO-D CASE.
!
      U  = PAR(k,1) - 2.0d0*Xm(k,J-1)
      Vm = PAR(k,2) - 2.0d0*Xm(k,J-1)
!
      Dm(k,J) = U*Vm/(4.0d0*(U+Vm))
      GO TO 40
!
   43 CONTINUE
!
!----- THRED CASE.
!
      U  = PAR(k,1)-2.0d0*Xm(k,J-1)
      Vm = PAR(k,2)-2.0d0*Xm(k,J-1)
      W  = PAR(k,3)-2.0d0*Xm(k,J-1)
!
      Dm(k,J) = 3.0d0*U*Vm*W/(1.0d1*(U*Vm+Vm*W+U*W))
      GO TO 40
!
   44 CONTINUE
   45 CONTINUE
   46 CONTINUE
   47 CONTINUE
   48 CONTINUE
   49 CONTINUE
   50 CONTINUE
!
      Dm(k,J) = (Xm(k,J)-Xm(k,J-1))/5.0d0
!
   40 CONTINUE
!
!-----PRINT OUT GEOMETRY DATA.
!
      WRITE(*, 27)
      WRITE(*, 23) k
      WRITE(*, 24)
      WRITE(*, 25) VOL(k,1), Dm(k,1)
      WRITE(*, 26) Am(k,1),  Xm(k,1)
!
   23 FORMAT('  FRACTURED MEDIUM #',i3.3,/, &
     &       '  CONTINUUM     IDENTIFIER       VOLUME      NODAL DISTANCE      INTERFACE AREA   INTERFACE DISTANCE')
   24 FORMAT(84X,'FROM FRACTURES',/)
   25 FORMAT('  1-FRACTURES      * *    ',2(4X,E12.5))
   26 FORMAT(66X,E12.5,7X,E12.5)
   27 FORMAT(//' ==================== GEOMETRY DATA, NORMALIZED TO A DOMAIN OF UNIT VOLUME =========================',//)
!
      DO 100 M = 2, J
         WRITE(*, 101) M, NA(M), VOL(k,M), Dm(k,M)
         IF(M /= J) WRITE(*, 102) Am(k,M), Xm(k,M)
  100 CONTINUE
!
!
!
      END DO DO_k
!
!
!
  101 FORMAT(' ',I2,'-MATRIX',9X,'*',A1,'*',8X,E12.5,4X,E12.5)
  102 FORMAT(66X,E12.5,7X,E12.5)
!
      WRITE(*,103)
  103 FORMAT(/,1X,99('='))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of GEOM
!
!
      RETURN
      END SUBROUTINE GEOM
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      REAL(KIND = 8) FUNCTION PROX(k,X)
!
         USE MINC_Parameters, ONLY: L, PAR
!
!
!-----THE PROXIMITY FUNCTION PROX(k,X) REPRESENTS THE FRACTION OF
!     MATRIX VOLUME [VM=(1.-VOL(1))*V0 WITHIN A DOMAIN V0] WHICH
!     IS WITHIN A DISTANCE X FROM THE FRACTURES.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      SAVE ICALL,A,B,C,D
!     NOW ASSIGN DATA FOR STANFORD LARGE RESERVOIR MODEL.
      DATA A,B,C,D/.263398,.190754,.2032,.191262/
      DATA ICALL/0/
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of PROX
!
!
      ICALL=ICALL+1
      IF(ICALL == 1) WRITE(UNIT = *, FMT = 6000)
 6000 FORMAT(/,'PROX 1.0 ............................... 22 JANUARY   1990',6X,'CALCULATE PROXIMITY FUNCTIONS FOR DIFFERENT MATRIX BLOCK SHAPES')
!
      GOTO(1,2,3,4,4,4,1,1,1,1),L
!
    1 CONTINUE
!
!----- ONE-D CASE.
!
      PROX = 2.0d0*X/PAR(k,1)
      IF(X >= PAR(k,1)/2.0d0) PROX = 1.0d0
      RETURN
!
    2 CONTINUE
!
!----- TWO-D CASE.
!     THE MATRIX BLOCKS HAVE THICKNESS OF PAR(k,1) AND PAR(k,2),
!     RESPECTIVELY, MEASURED PERPENDICULAR TO THE FRACTURES.
!     THE PROXIMITY FUNCTION IS VALID FOR ARBITRARY ANGLE
!     BETWEEN THE FRACTURE SETS.
!
      PROX =  2.0d0 * (PAR(k,1)+PAR(k,2)) * X / (PAR(k,1)*PAR(k,2)) - 4.0d0*X*X/(PAR(k,1)*PAR(k,2))
!
      IF(X >= PAR(k,1)/2.00d0 .OR. X >= PAR(k,2)/2.0d0) PROX = 1.0d0
      RETURN
!
    3 CONTINUE
!
!----- THREE DIMENSIONAL CASE.
!
      U    = 2.0d0*X/PAR(k,1)
      V    = 2.0d0*X/PAR(k,2)
      W    = 2.0d0*X/PAR(k,3)
      PROX = U*V*W-(U*V+U*W+V*W)+U+V+W
      IF(U >= 1.0d0 .OR. V >= 1.0d0 .OR. W >= 1.0d0) PROX = 1.0d0
      RETURN
!
    4 CONTINUE
!
!
!***** MATRIX OF STANFORD LARGE RESERVOIR MODEL *****
!
!
!     RECTANGULAR BLOCKS IN LAYERS B1,B2,M1,M2,T1.
      VR = 8.0d0*X**3-(8.0d0*B+4.0d0*A)*X**2+(4.0d0*A*B+2.0d0*B**2)*X
!
      IF(X >= B/2.0d0) VR = A*B*B
!
!     TRIANGULAR BLOCKS IN LAYERS B1,B2,M1,M2,T1.
      VT =  (6.0d0+4.0d0*SQRT(2.0d0))*X**3 - (A*(6.0d0+4.0d0*SQRT(2.0d0))/2.0d0   &
     &     + 2.0d0*B*(2.0d0+SQRT(2.0d0)))*X**2 + (A*B*(2.0d0+SQRT(2.0d0))+B*B)*X
!
      IF(X >= B/(2.0d0+SQRT(2.0d0))) VT = A*B*B/2.0d0
!
!     RECTANGULAR BLOCKS IN LAYER T2.
      VRT2 = 8.0d0*X**3-(8.0d0*D+4.0d0*C)*X**2+(4.0d0*C*D+2.0d0*D**2)*X
!
      IF(X >= D/2.0d0) VRT2 = C*D*D
!
!     TRIANGULAR BLOCKS IN LAYER T2.
      VTT2 = (6.0d0+4.0d0*SQRT(2.0d0))*X**3 - (C*(6.0d0+4.0d0*SQRT(2.0d0))/2.0d0   &
     &      + 2.0d0*D*(2.0d0+SQRT(2.0d0)))*X**2 + (C*D*(2.0d0+SQRT(2.0d0))+D*D)*X
!
      IF(X >= D/(2.0d0+SQRT(2.0d0))) VTT2 = C*D*D/2.0d0
!
      IF(L.EQ.4) GOTO 14
      IF(L.EQ.5) GOTO 15
      IF(L.EQ.6) GOTO 16
!
!***** NOW COMPUTE TOTAL MATRIX VOLUME WITHIN DISTANCE X.
   14 V = 5.0d0*(5.0d0*VR+4.0d0*VT) + 5.0d0*VRT2+4.0d0*VTT2
!
!-----AVERAGE PROXIMITY FUNCTION FOR ENTIRE ROCK LOADING.
!
      VTOT =35.0d0*A*B**2+7.0d0*C*D**2
!     VOLUME FRACTION.
      PROX = V/VTOT
      RETURN
!
   15 PROX = (5.0d0*VR+4.*VT)/(7.0d0*A*B*B)
!-----PROXIMITY FUNCTION FOR FIVE BOTTOM LAYERS.
      RETURN
!
   16 PROX = (5.0d0*VRT2+4.*VTT2)/(7.0d0*C*D*D)
!-----PROXIMITY FUNCTION FOR TOP LAYER.
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of PROX
!
!
      RETURN
!
!
      END FUNCTION PROX
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE INVER(k, F,X,XL,XR)
!
!
!===== THIS ROUTINE INVERTS THE PROXIMITY FUNCTION, TO GIVE A
!     DISTANCE *X* FROM FRACTURE FACES FOR A DESIRED FRACTION *F* OF
!     MATRIX VOLUME.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      REAL(KIND=8), PARAMETER :: TOL = 1.0d-10
!
      INTEGER :: ICALL = 0
!
      SAVE ICALL
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of INVER
!
!
      ICALL = ICALL+1
      IF(ICALL == 1) WRITE(UNIT = *, FMT = 6000)
 6000 FORMAT(/,'INVER 1.0 .............................. 22 JANUARY   1990',6X,'INVERT A MONOTONIC FUNCTION THROUGH NESTED BISECTIONS')
!
!-----CHECK AND ADJUST UPPER LIMIT OF SEARCH INTERVAL.
!
   22 FR = PROX(k, XR)
      IF(FR > F) GO TO 20
      XR = 2.0d0*XR
      GO TO 22
!
!-----PERFORM ITERATIVE BISECTING, TO OBTAIN A SEQUENCE OF NESTED
!     INTERVALS CONTAINING THE DESIRED POINT, X.
!
   20 XMID = (XR+XL)/2.0d0
!
      IF(XR-XL <= TOL*XR) GO TO 21
      FMID = PROX(k, XMID)
      IF(FMID <= F) XL=XMID
      IF(FMID >= F) XR=XMID
      GOTO 20
!
   21 CONTINUE
!
!-----COME HERE FOR CONVERGENCE.
!
      X = XMID
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of INVER
!
!
      RETURN
      END SUBROUTINE INVER
!
!
!
!***************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***************************************************************************
!
!
!
      SUBROUTINE MINCME
!
!
         USE MeshMaker_Data, Xa => X, Ya => Y, Za => Z
         USE Grid_Generation_Parameters
!
         USE MINC_Parameters
!
!
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!=====THIS ROUTINE WORKS SEQUENTIALLY THROUGH THE ELEMENTS OF THE
!     PRIMARY MESH, ASSIGNING ALL SECONDARY ELEMENTS AND INTRA-BLOCK
!     CONNECTIONS.
!
!
      CHARACTER(LEN=1) :: EL1, EM1
      CHARACTER(LEN=4) :: EL2, EM2
      CHARACTER(LEN=5) :: ELE, EC1, EC2, ELREF
!
      CHARACTER(LEN=5) :: header(2)
      CHARACTER(LEN=3) :: activity
      CHARACTER(LEN=1) :: first
!
      INTEGER :: NELP,NELAP,NCONP
      INTEGER :: ICALL = 0
!
      LOGICAL :: volume_flag
!
! - GJM: Begin modification for 8-character elements
!
      CHARACTER(LEN=8) :: ELE_8, EREF_8, EC1_8, EC2_8
      CHARACTER EL1_3*3,EL2_5*5,EM1_3*3,EM2_5*5
!
      CHARACTER MESHX*4
!
! - GJM: End modification for 8-character elements
!
      SAVE ICALL
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of MINCME
!
!
      ICALL=ICALL+1
!     IF(ICALL.EQ.1) WRITE(*,899)
! 899 FORMAT('MINCME            1.0   14 FEBRUARY  1990',6X,'PROCESS PRIMARY MESH FROM FILE *MESH*, AND WRITE SECONDARY MESH ON FILE *MINC*')
!
      IF(ICALL == 1) WRITE(*,6000)
 6000 FORMAT(/,'MINCME 2.0 .............................. 4 MARCH     2015',6X,'PROCESS PRIMARY MESH FROM FILE *MESH*, AND WRITE SECONDARY MESH ON FILE *MINC*')
!
      REWIND (UNIT = MESH_Unit)
      REWIND (UNIT = MINC_Unit)
!
      WRITE(*,2)
    2 FORMAT(/' READ PRIMARY MESH FROM FILE *MESH*')
!
!*****READ ELEMENT DATA FROM FILE *MESH*.*******************************
!
      READ(MESH_Unit,1) (header(i),I=1,2)
    1 FORMAT(2A5)
!
      IF(header(1) /= 'ELEME') THEN
         WRITE(*,4) header(1)
         STOP
      END IF
!
    4 FORMAT(' HAVE READ UNKNOWN BLOCK LABEL "',A5,'" ON FILE *MESH* WHILE EXPECTING THE ELEMENT BLOCK', &
     &       ' --- STOP EXECUTION ---')
!
! - GJM: Begin modification for 8-character elements
!
      IF(header(2)(1:4) == 'ext2') THEN
         IF( ElemNameLength /= 8 ) THEN
            WRITE(*,6010) MESHX, ElemNameLength
            STOP
         END IF
         MESHX = 'ext2'
         WRITE(MINC_Unit, 5006)
      ELSE
         IF(ElemNameLength /= 5) THEN
            WRITE(*,6000) MESHX, ElemNameLength
            STOP
         END IF
         WRITE(MINC_Unit,6)
      END IF
!
 5006 FORMAT('ELEMEext2')
    6 FORMAT('ELEME')
!
 6010 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',/,   &
     &         T2,'The mesh type inicator MESHX = ',a4,' is in conflict with the # of characters in the element name <ElemNameLength> = ',i1,   &
     &       /,T32,'               CORRECT AND TRY AGAIN',    &
     &       //,20('ERROR-'))
!
! - GJM: End modification for 8-character elements
!
!
      NEL  = 0
      NELA = 0
!
      Nel_2nd  = 0
      NelA_2nd = 0
!
      volume_flag = .FALSE.
      num_dom     = 0
!
      fractured = .FALSE.  ! ... Array operation
!
    9 CONTINUE
!
! - GJM: Begin modification for 8-character elements
!
      IF(MESHX.EQ.'ext2') THEN
         READ(MESH_Unit,5010) EL1_3, EL2_5, MA, VOLX, AHTX, X, Y, Z, activity
         IF(EL1_3.EQ.'   '.AND.EL2_5.EQ.'     ') GO TO 40
      ELSE
         READ(MESH_Unit,10) EL1, EL2, MA, VOLX, AHTX, X, Y, Z, activity
         IF(EL1.EQ.' '.AND.EL2.EQ.'    ') GO TO 40
      END IF
!
 5010 FORMAT(A3,A5, 7X,I5,ES10.4E1,ES10.4E1,10X,3F10.3,1x,a3)
   10 FORMAT(A1,A4,10X,I5,ES10.4E1,ES10.4E1,10X,3ES10.4E1,1x,a3)
!
! - GJM: End modification for 8-character elements
!
         NEL = NEL+1
!
         IF( activity(1:1) == 'I' .OR. activity(1:1) == 'V' ) THEN
            MA      = ABS(MA)
            MA1     = (MA-1)*J+1
            Nel_2nd = Nel_2nd + 1
            GO TO 8
         END IF
!
         IF( (activity(1:1) /= 'I' .OR. activity(1:1) /= 'V') .AND. VOLX <= 0.0e0 ) THEN
            MA      = ABS(MA)
            MA1     = (MA-1)*J+1
            Nel_2nd = Nel_2nd + 1
            GO TO 8
         END IF
!
         NELA =  NELA + 1
!
! ...... For domains that will remain unfractured
!
         IF( MA < 0 ) THEN
            Nel_2nd  = Nel_2nd + 1
            Nela_2nd = Nela_2nd + 1
            MA       = ABS(MA)
            MA1      = (MA-1)*J+1
            GO TO 8
         END IF
!
!-----COME HERE FOR ACTIVE ELEMENTS
!
         VVV = VOL(MA,1)*VOLX
         MA1 = (MA-1)*J+1
         AH  = AHTX
!
         SELECT CASE (matrix_to_matrix_flow(1:1))
         CASE('V','v')
            AH = VOL(MA,1)*AHTX
         END SELECT
!
         Nel_2nd  = Nel_2nd + 1
         Nela_2nd = Nela_2nd + 1
!
! - GJM: Begin modification for 8-character elements
!
         IF(MESHX == 'ext2') THEN
!
            i_CHAR = ICHAR(EL1_3(1:1))
            IF( i_CHAR >= 65 .OR. i_CHAR <= 90 ) THEN
               first = CHAR(i_CHAR + 32)
            ELSE IF( i_CHAR >= 97 .OR. i_CHAR <= 122 ) THEN
               first = CHAR(i_CHAR - 32)
            ELSE IF( i_CHAR == ICHAR(' ') ) THEN
               first = 'f'
            END IF
!
            ELE_8 = first//EL1_3(2:3)//EL2_5
            WRITE(MINC_Unit,5007) ELE_8, MA1, VVV, AH, X, Y, Z, activity
!
         ELSE
!
            i_CHAR = ICHAR(EL1)
            IF( i_CHAR >= 65 .OR. i_CHAR <= 90 ) THEN
               first = CHAR(i_CHAR + 32)
            ELSE IF( i_CHAR >= 97 .OR. i_CHAR <= 122 ) THEN
               first = CHAR(i_CHAR - 32)
            ELSE IF( i_CHAR == ICHAR(' ') ) THEN
               first = 'f'
            END IF
!
            ELE = first//EL2
            WRITE(MINC_Unit,7) ELE, MA1, VVV, AH, X, Y, Z, activity
!
         END IF
!
 5007 FORMAT(A8, 7X,I5,2ES10.4E1,10X,3F10.3,1x,A3)
    7 FORMAT(A5,10X,I5,2ES10.4E1,10X,3ES10.4E1,1x,A3)
!
! - GJM: End modification for 8-character elements
!
!
!     STORE FRACTURE ELEMENT NAME FOR LATER CROSS REFERENCING
!
! ----------------
! ........ GJM - 2/14/2015
! ----------------
!
         IF(MESHX == 'ext2') THEN
            name(NEL) = ELE_8
         ELSE
            name(NEL) = ELE
         END IF
!
         fractured(NEL) = .TRUE.
!
! ----------------
! ........ GJM - 2/14/2015
! ........ End modification
! ----------------
!
!-----FOR EACH PRIMARY ELEMENT, ASSIGN *J* SECONDARY ELEMENTS.
!
!
! - GJM: Begin modification for 8-character elements
!
         Nel_2nd  = Nel_2nd + 1
         Nela_2nd = Nela_2nd + 1
!
         IF(MESHX.EQ.'ext2') THEN
!
            DO M=2,J
!
               VVV = VOL(MA,M)*VOLX
!
               IF( .NOT. uniform_matrix_properties ) THEN
                  MA2 = (MA-1)*J + m
               ELSE
                  MA2 = (MA-1)*J + 2
               END IF
               AH  = 0.0e0
!
               IF( matrix_to_matrix_flow(1:1) == 'V' .OR. matrix_to_matrix_flow(1:1) == 'v' ) AH = VOL(MA,M)*AHTX

               IF( J > 2 ) EL1_3(1:1) = NA(m)
               WRITE(MINC_Unit,5010) EL1_3, EL2_5, MA2, VVV, AH, X, Y, Z, activity
!
            END DO
!
         ELSE
!
            DO M=2,J
!
               ELE = EL1//EL2
               VVV = VOL(MA,M)*VOLX
!
               IF( .NOT. uniform_matrix_properties ) THEN
                  MA2 = (MA-1)*J + m
               ELSE
                  MA2 = (MA-1)*J + 2
               END IF
!
               AH  = 0.0e0
!
               IF( matrix_to_matrix_flow(1:1) == 'V' .OR. matrix_to_matrix_flow(1:1) == 'v' ) AH = VOL(MA,M)*AHTX
!
               IF( J > 2 ) ELE = NA(m)//EL2
               WRITE(MINC_Unit,7) ELE, MA2, VVV, AH, X, Y, Z, activity
!
            END DO
!
         END IF
!
! - GJM: End modification for 8-character elements
!
!
         GO TO 9
!
    8    CONTINUE
!-----COME HERE FOR INACTIVE ELEMENTS. THEY WILL NOT BE SUBJECTED TO
!     THE MINC PROCEDURE, BUT WILL BE HANDLED AS A SINGLE CONTINUUM.
!
! - GJM: Begin modification for 8-character elements
!
         IF(MESHX == 'ext2') THEN
            ELE_8 = EL1_3//EL2_5
            WRITE(MINC_Unit,5007) ELE_8, MA1, VOLX, AHTX, X, Y, Z, activity
         ELSE
            ELE = EL1//EL2
            WRITE(MINC_Unit,7) ELE, MA1, VOLX, AHTX, X, Y, Z, activity
         END IF
!
! - GJM: End modification for 8-character elements
!
!
         GO TO 9
!
   40 CONTINUE
 8100 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',/,   &
     &         T29,'The array dimension MNEL is insufficient for the MINC needs',   &
     &       /,T33,'               CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
!
!
!
      WRITE(MINC_Unit,103)
      WRITE(MINC_Unit,104)
  103 FORMAT('     ')
  104 FORMAT('CONNE')
!
!-----NOW REDEFINE ELEMENT COUNTERS.
      NBC=0
      IF(NELA.NE.0) NBC  = NEL-NELA
      IF(NELA.EQ.0) NELA = NEL
!     PARAMETERS FOR PRIMARY MESH
      NELP  = NEL
      NELAP = NELA
!     PARAMETERS FOR SECONDARY MESH
      NELA = J*NELA
      NEL  = NELA + NBC
!
!
!****READ CONNECTION DATA.*********************************************
!
      N     = 0
      NCONP = 0
!
   45 READ(MESH_Unit,1) header(1)
!
      IF( header(1)(1:1) == ' ') THEN
         GO TO 45
      ELSE IF(header(1) /= 'CONNE') THEN
         WRITE(*,1201) header(1)
         STOP
      END IF
!
 1201 FORMAT(' HAVE READ UNKNOWN BLOCK LABEL "',A5,'" ON FILE *MESH* WHILE EXPECTING THE CONNECTION BLOCK', &
     &       ' --- STOP EXECUTION ---')
!
      SELECT CASE (matrix_to_matrix_flow(1:1))
      CASE('A', 'a', 'V', 'v')
         J_Max = J
      CASE DEFAULT
         J_Max = 2
      END SELECT
!
 1200 CONTINUE
!
! - GJM: Begin modification for 8-character elements
!
      IF(MESHX == 'ext2') THEN
         READ(MESH_Unit,20) EL1_3,EL2_5,EM1_3,EM2_5, ISOT, Dist1, Dist2, AREAX, BETAX
         IF(EL1_3 == '   ' .AND. EL2_5 == '     ') GO TO 1400
         IF(EL1_3(1:1).EQ.'+') GO TO 1400
      ELSE
         READ(MESH_Unit,20) EL1, EL2, EM1, EM2, ISOT, Dist1, Dist2, AREAX, BETAX
         IF(EL1 == '<') GOTO 1400
         IF(EL1 == ' ' .AND. EL2 == '    ') GOTO 1400
         IF(EL1 == '+') GO TO 1400
      END IF
!
 5020 FORMAT(2(A3,A5), 9X,I5,4ES10.4E1)
   20 FORMAT(2(A1,A4),15X,I5,4ES10.4E1)
!
! - GJM: End modification for 8-character elements
!
      N     = N+1
      NCONP = NCONP+1
!
!-----NOW DETERMINE ACTIVITY STATUS OF ELEMENTS AT THIS CONNECTION.
!     ASSIGN THE "WOULD-BE" FRACTURE ELEMENTS, AND SEE WHETHER THEY
!     APPEAR IN THE LIST OF ACTIVE ELEMENTS
!
      IF(MESHX == 'ext2') THEN
!
         i_CHAR = ICHAR(EL1_3(1:1))
         IF( i_CHAR >= 65 .OR. i_CHAR <= 90 ) THEN
            first = CHAR(i_CHAR + 32)
         ELSE IF( i_CHAR >= 97 .OR. i_CHAR <= 122 ) THEN
            first = CHAR(i_CHAR - 32)
         ELSE IF( i_CHAR == ICHAR(' ') ) THEN
            first = 'f'
         END IF
!
         EC1_8 = first//EL1_3(2:3)//EL2_5
!
         i_CHAR = ICHAR(EM1_3(1:1))
         IF( i_CHAR >= 65 .OR. i_CHAR <= 90 ) THEN
            first = CHAR(i_CHAR + 32)
         ELSE IF( i_CHAR >= 97 .OR. i_CHAR <= 122 ) THEN
            first = CHAR(i_CHAR - 32)
         ELSE IF( i_CHAR == ICHAR(' ') ) THEN
            first = 'f'
         END IF
!
         EC2_8 = first//EM1_3(2:3)//EM2_5
!
      ELSE
!
         i_CHAR = ICHAR(EL1)
         IF( i_CHAR >= 65 .OR. i_CHAR <= 90 ) THEN
            first = CHAR(i_CHAR + 32)
         ELSE IF( i_CHAR >= 97 .OR. i_CHAR <= 122 ) THEN
            first = CHAR(i_CHAR - 32)
         ELSE IF( i_CHAR == ICHAR(' ') ) THEN
            first = 'f'
         END IF
!
         EC1 = first//EL2
!
         i_CHAR = ICHAR(EM1)
         IF( i_CHAR >= 65 .OR. i_CHAR <= 90 ) THEN
            first = CHAR(i_CHAR + 32)
         ELSE IF( i_CHAR >= 97 .OR. i_CHAR <= 122 ) THEN
            first = CHAR(i_CHAR - 32)
         ELSE IF( i_CHAR == ICHAR(' ') ) THEN
            first = 'f'
         END IF
!
         EC2 = first//EM2
!
      END IF
!
!     INITIALIZE ACTIVITY INDEX
!
      IAC = 0
      N1  = 0
      N2  = 0
!
! ----------------
! ........ GJM - 2/14/2015
! ----------------
!
      IF(MESHX == 'ext2') THEN
!
         DO_Loop1: DO i=1,NELP
!
            EREF_8 = name(i)
!
            IF(EREF_8 /= EC1_8) GO TO 516
            IAC = IAC+1
            N1  = I
  516       IF(EREF_8 /= EC2_8) GO TO 517
            IAC = IAC+1
            N2  = I
  517       CONTINUE
            IF(IAC == 2) GO TO 18
!
         END DO DO_Loop1
!
         IF( .NOT. fractured(N1) ) EC1_8 = EL1_3//EL2_5
         IF( .NOT. fractured(N2) ) EC2_8 = EM1_3//EM2_5
!
      ELSE
!
         DO_Loop2: DO i=1,NELP
!
            ELREF = name(i)(1:5)
!
            IF(ELREF /= EC1) GO TO 16
            IAC = IAC+1
            N1  = I
   16       IF(ELREF /= EC2) GO TO 17
            IAC = IAC+1
            N2  = I
   17       CONTINUE
!
            IF(IAC == 2) GO TO 18
!
         END DO DO_Loop2
!
         IF( .NOT. fractured(N1) ) EC1 = EL1//EL2
         IF( .NOT. fractured(N2) ) EC2 = EM1//EM2
!
      END IF
!
! ----------------
! ........ GJM - 2/14/2015
! ........ End modification
! ----------------
!
   18 CONTINUE
!
!-----GENERATE GLOBAL FRACTURE CONNECTION DATA.
!
      IF(MESHX == 'ext2') THEN
         WRITE(MINC_Unit,5105) EC1_8,EC2_8,ISOT,Dist1,Dist2,AREAX,BETAX
      ELSE
         WRITE(MINC_Unit,105) EC1,EC2,ISOT,Dist1,Dist2,AREAX,BETAX
      END IF
!
 5105 FORMAT(2A8, 9X,I5,4ES10.4E1)
  105 FORMAT(2A5,15X,I5,4ES10.4E1)
!
!
!
      SELECT CASE (matrix_to_matrix_flow(1:1))
      CASE(' ')
         GO TO 1200
      CASE('A', 'a', 'F', 'f')
         CONTINUE
      CASE('V', 'v')
         IF( ABS(BETAX) == 1.0e0 ) THEN
            CONTINUE
         ELSE
            GO TO 1200
         END IF
      CASE DEFAULT
         GO TO 1200
      END SELECT
!
! ... Check if both are unfractured media
!
      IF( ( .NOT. fractured(N1) ) .AND. ( .NOT. fractured(N2) ) ) GO TO 1200
!
! ... ASSIGN GLOBAL CONNECTIONS BETWEEN MATRIX CONTINUA
!
! - GJM: Begin modification for 8-character elements
!
      IF(MESHX == 'ext2') THEN
!
         DO m = 2,J_Max
!
            N = N+1
!
            SELECT CASE (matrix_to_matrix_flow(1:1))
            CASE('A', 'a', 'F', 'f')
               AX = AREAX
            CASE('V', 'v')
               AX = AREAX * VOL(1,m)
            END SELECT
!
            IF( J > 2 ) THEN
               IF( fractured(N1) ) EC1_8 = NA(m)//EL1_3(2:3)//EL2_5
               IF( fractured(N2) ) EC2_8 = NA(m)//EM1_3(2:3)//EM2_5
            ELSE
               IF( fractured(N1) ) EC1_8 = EL1_3//EL2_5
               IF( fractured(N2) ) EC2_8 = EM1_3//EM2_5
            END IF
!
            WRITE(MINC_Unit,5105) EC1_8, EC2_8, ISOT, Dist1, Dist2, AX, BETAX
!
         END DO
!
      ELSE
!
         DO m = 2,J_Max
!
            N = N+1
!
            SELECT CASE (matrix_to_matrix_flow(1:1))
            CASE('A', 'a', 'F', 'f')
               AX = AREAX
            CASE('V', 'v')
               AX = AREAX * VOL(1,m)
            END SELECT
!
            IF( J > 2 ) THEN
               IF( fractured(N1) ) EC1 = NA(m)//EL2
               IF( fractured(N2) ) EC2 = NA(m)//EM2
            ELSE
               IF( fractured(N1) ) EC1 = EL1//EL2
               IF( fractured(N2) ) EC2 = EM1//EM2
            END IF
!
            WRITE(MINC_Unit,105) EC1, EC2, ISOT, Dist1, Dist2, AX, BETAX
!
         END DO
!
      END IF
!
! - GJM: End modification for 8-character elements
!
!
      GO TO 1200
!
!-----END OF CONNECTION DATA.-------------------------------------------
!
 1400 CONTINUE
      WRITE(*,5) NELP,NELAP,NCONP
    5 FORMAT('            THE PRIMARY MESH HAS',I7,' ELEMENTS (',I7,' ACTIVE) AND ',I8,' CONNECTIONS (INTERFACES) BETWEEN THEM')
!
!-----NOW ASSIGN INTRA-BLOCK CONNECTION DATA.  LOOP AGAIN OVER ELEMENTS
!
      REWIND (UNIT = MESH_Unit)
      READ(MESH_Unit,1) (header(i),i=1,2)
      I=0
!
! - GJM: Begin modification for 8-character elements
!
  110 IF(MESHX == 'ext2') THEN
         READ(MESH_Unit,5010) EL1_3, EL2_5, MA, VOLX, AHTX
         IF(EL1_3 == '   ' .AND. EL2_5 == '     ') GO TO 111
      ELSE
         READ(MESH_Unit,10) EL1, EL2, MA, VOLX, AHTX
         IF(EL1 == ' ' .AND. EL2 == '    ') GO TO 111
      END IF
!
! - GJM: End modification for 8-character elements
!
      I = I+1
!
!-----FOR INACTIVE ELEMENTS, DO NOT ASSIGN INTRABLOCK CONNECTIONS
!
      IF( .NOT. fractured(i) ) GO TO 110
!
!-----COME HERE TO ASSIGN INTRABLOCK CONNECTIONS FOR ACTIVE ELEMENTS ONLY
!
      IF(MESHX == 'ext2') THEN
         i_CHAR = ICHAR(EL1_3(1:1))
      ELSE
         i_CHAR = ICHAR(EL1)
      END IF
!
      IF( i_CHAR >= 65 .OR. i_CHAR <= 90 ) THEN
         first = CHAR(i_CHAR + 32)
      ELSE IF( i_CHAR >= 97 .OR. i_CHAR <= 122 ) THEN
         first = CHAR(i_CHAR - 32)
      ELSE IF( i_CHAR == ICHAR(' ') ) THEN
         first = 'f'
      END IF
!
      BETAX = 0.0d0
!
! - GJM: Begin modification for 8-character elements
!


      IF(MESHX == 'ext2') THEN
!
         DO m = 2, J
            N     = N+1
            AREAX = VOLX * Am(MA,m-1)
            IF( J > 2 ) EL1_3(1:1) = NA(m)
            WRITE(MINC_Unit,5104) first,      EL1_3(2:2),EL1_3(3:3),EL2_5,   &
     &                            EL1_3(1:1), EL1_3(2:2),EL1_3(3:3),EL2_5, Dm(MA,m-1), Dm(MA,m), AREAX
            IF(J > 2) first = NA(m)
         END DO
!
      ELSE
!
         DO m = 2, J
            N     = N+1
            AREAX = VOLX * Am(MA,m-1)
            IF( J > 2 ) EL1 = NA(m)
            WRITE(MINC_Unit,102) first, EL2, EL1, EL2, Dm(MA,m-1), Dm(MA,m), AREAX
            IF(J > 2) first = NA(m)
         END DO
!
      END IF
!
 5104 FORMAT(2(A1,A1,A1,A5),13X,'1',3ES10.4E1)
 5102 FORMAT(2(A3,A5),13X,'1',3ES10.4E1)
  102 FORMAT(2(A1,A4),19X,'1',3ES10.4E1)
!
! - GJM: Begin modification for 8-character elements
!
      GO TO 110
!
  111 CONTINUE
      NCON=N
      WRITE(*,113) NEL,NELA,NCON
  113 FORMAT(/' WRITE SECONDARY MESH ON FILE *MINC*',/,   &
     &        '          THE SECONDARY MESH HAS',I7,' ELEMENTS (',I7,' ACTIVE) AND ',I8,' CONNECTIONS (INTERFACES) BETWEEN THEM')
!
      WRITE(MINC_Unit,103)
      ENDFILE (UNIT = MINC_Unit)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of MINCME
!
!
      RETURN
!
      END SUBROUTINE MINCME
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      CHARACTER*5 FUNCTION Name_of_5Character_Element(NLXYZ1,NLXYZ2,NLXYZ3,i,j,k)
!
         USE Grid_Generation_Parameters
         USE MeshMaker_Data, ONLY: letter_shift
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NLXYZ1, NLXYZ2, NLXYZ3, i, j, k
!
      INTEGER :: nnn111, nnn222, nnn333, nnn444, nnn555
      INTEGER :: iii222, iii333, mmm444, n1_adj
!
      CHARACTER(LEN=1) :: L1, L2, L3, L4, L5
!
      LOGICAL :: First_Call = .TRUE.
!
      SAVE First_Call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of <Name_of_5Character_Element>
!
!
      IF(First_Call) THEN
         First_Call = .FALSE.
         WRITE(*,6000)
      END IF
!
!*********************************************************************
!*                                                                   *
!*                       DETERMINATION PROCESS                       *
!*                                                                   *
!*********************************************************************
!
      IF( NLXYZ1 > 1 ) THEN
!
         IF(NLXYZ1 > 100) THEN
!
            nnn333 = 1+((k-1)/100)
            IF(nnn333 > 36) THEN
               WRITE(*,8100)
               STOP
            ELSE
               iii222 = (j-1)+(i-1)*NLXYZ2
               nnn222 = 1+MOD(iii222,36)
               nnn111 = 1+(iii222/36)
            END IF
!
         ELSE
!
            iii333 = (j-1)+(i-1)*NLXYZ2
            nnn333 = 1+MOD(iii333,36)
            iii222 = 1+(iii333/36)
!
            IF(iii222 == 36) THEN
               nnn222 = iii222
               nnn111 = nnn222/36
            ELSE IF(iii222 < 36) THEN
               nnn222 = iii222
               nnn111 = 1+nnn222/36
            ELSE
               IF(MOD(iii222,36) /= 0) THEN
                  nnn222 = MOD(iii222,36)
                  nnn111 = 1+(iii222/36)
               ELSE
                  nnn222 = 36
                  nnn111 = iii222/36
               END IF
            END IF
!
         END IF
!
         mmm444 = mod(k-1,100)
         nnn444 = mmm444/10
         nnn555 = mod(mmm444,10)
!
      ELSE
!
         IF(NLXYZ2 /= 1) THEN
!
            IF(NLXYZ2 > 100) THEN
!
               nnn333 = 1+((j-1)/100)
!
               IF(nnn333 > 36) THEN
                  nnn333 = 1+MOD(nnn333,36)
                  iii333 = 1+nnn333/36
                  IF(iii333 > 36) THEN
                     nnn222 = 1+MOD(iii333,36)
                     nnn111 = 1+nnn222/36
                  ELSE
                     nnn222 = iii333
                     nnn111 = 1
                  END IF
               ELSE
                  nnn222 = 1+MOD((i-1),36)
                  nnn111 = 1+((i-1)/36)
               END IF
!
            ELSE
!
               nnn333 = 1+MOD((i-1),36)
               nnn222 = 1+((i-1)/36)
               iii222 = MOD(nnn222,36)
               IF(nnn222 > 36) THEN
                  IF( iii222 /= 0 ) THEN
                     nnn111 = 1 + nnn222/36
                     nnn222 = iii222
                  ELSE
                     nnn111 = 1 + nnn222/36
                     nnn222 = 1 + iii222
                  END IF
               ELSE
                  nnn111 = 1 + nnn222/36
                  IF(iii222 == 0 ) THEN
                     nnn111 = nnn222/36
                  END IF
               END IF
!
            END IF
!
            mmm444 = mod(j-1,100)
            nnn444 = mmm444/10
            nnn555 = mod(mmm444,10)
!
         ELSE
!
            IF(NLXYZ3 > 100) THEN
!
               nnn333 = 1+((i-1)/100)
               IF(nnn333 > 36) THEN
                  nnn333 = 1+MOD(nnn333,36)
                  iii333 = 1+nnn333/36
                  IF(iii333 > 36) THEN
                     nnn222 = 1+MOD(iii333,36)
                     nnn111 = 1+nnn222/36
                  ELSE
                     nnn222 = iii333
                     nnn111 = 1
                  END IF
               ELSE
                  nnn222 = 1
                  nnn111 = 1
               END IF
!
            ELSE
!
               nnn333 = 1
               nnn222 = 1
               nnn111 = 1
!
            END IF
!
            mmm444 = mod(i-1,100)
            nnn444 = mmm444/10
            nnn555 = mod(mmm444,10)
!
         END IF
!
      END IF
!
!***********************************************************************
!*                                                                     *
!*                     RETURNING THE onoma VALUE                       *
!*                                                                     *
!***********************************************************************
!
      n1_adj = nnn111 + 10 + letter_shift
!
      IF(n1_adj <= 36) THEN
         L1 = NA(n1_adj)
      ELSE
         L1 = NA(n1_adj-36)
      END IF
!
      L2 = NA(nnn222)
      L3 = NA(nnn333)
      L5 = NA(nnn555+1)
!
      IF(nnn444 >= 1) THEN
         L4 = NA(nnn444+1)
      ELSE
         L4 = '0'
      END IF
!
 1000 Name_of_5Character_Element = L1//L2//L3//L4//L5
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
 6000 FORMAT(/,'Name_of_5Character_Element 1.0 ......... 18 September 2000',6X,'Function returning 5-character names of the grid elements')
!
 8100 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',/,   &
     &         T19,'The number of elements along the shortest ',   &
     &             'dimension of a 3-D grid exceeds 3600',   &
     &       /,T33,'               CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of <Name_of_5Character_Element>
!
!
      RETURN
!
!
      END FUNCTION Name_of_5Character_Element
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      CHARACTER*8 FUNCTION Name_of_8Character_Element(NLXYZ1,NLXYZ2,NLXYZ3,i,j,k)
!
         USE Grid_Generation_Parameters
         USE MeshMaker_Data, ONLY: letter_shift
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NLXYZ1, NLXYZ2, NLXYZ3, i, j, k
!
      INTEGER :: i2, i3, i4, i5, i6, i7
      INTEGER :: n1, n2, n3, n4, n5, n6, n7, n8
      INTEGER :: m1, m2, m3, n8m1, n8m2, n8m3, ns, nsm1, nsm2, i8m2, ism1, n1_adj
!
      CHARACTER(LEN=1), DIMENSION(8) :: L
!
      LOGICAL :: First_Call = .TRUE.
!
      SAVE First_Call
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of <Name_of_8Character_Element>
!
!
      IF(First_Call) THEN
         First_Call = .FALSE.
         WRITE(*,6000)
      END IF
!
!*********************************************************************
!*                                                                   *
!*                            1-D SYSTEMS                            *
!*                                                                   *
!*********************************************************************
!
      IF( NLXYZ1 == 1 .AND. NLXYZ2 == 1 .AND. NLXYZ3 > 1 ) THEN
!
         L(1) = 'A'
!
         n2   = i / 1000000
         L(2) = NA( n2 + 1 )
!
         i3   = MOD(i,1000000)
         n3   = i3 / 100000
         L(3) = NA( n3 + 1 )
!
         i4   = MOD(i3,100000)
         n4   = i4 / 10000
         L(4) = NA( n4 + 1 )
!
         i5   = MOD(i4,10000)
         n5   = i5 / 1000
         L(5) = NA( n5 + 1 )
!
         i6   = MOD(i5,1000)
         n6   = i6 / 100
         L(6) = NA( n6 + 1 )
!
         i7   = MOD(i6,100)
         n7   = i7 / 10
         L(7) = NA( n7 + 1 )
!
         n8   = MOD(i7,10)
         L(8) = NA( n8 + 1 )
!
         n1_adj = ICHAR(L(1)) + letter_shift
         L(1)   = CHAR(n1_adj)
!
         Name_of_8Character_Element = L(1)//L(2)//L(3)//L(4)//L(5)//L(6)//L(7)//L(8)
!
         RETURN
!
      END IF
!
!*********************************************************************
!*                                                                   *
!*                            2-D SYSTEMS                            *
!*                                                                   *
!*********************************************************************
!
      IF( NLXYZ1 == 1 .AND. NLXYZ2 > 1 .AND. NLXYZ3 > 1 ) THEN
!
         IF( NLXYZ2 <= 100 ) THEN
!
            L(5) = 'A'
            L(6) = '0'
!
            n7   = j / 10
            IF( MOD(j,10) == 0 ) n7 = n7 - 1
            L(7) = NA( n7 + 1)
!
            n8 = MOD(j,10)
            IF(n8 == 0) n8 = 10
            L(8) = NA( n8 )
!
         ELSE
!
            n5 = j / 3600
            IF(MOD(j,3600) == 0) n5 = n5 - 1
            L(5) = NA( n5 + 1 )
!
            i6 = MOD(j,3600)
            IF(MOD(i6,100) /= 0) THEN
               n6 = i6 / 100
            ELSE
               n6 = i6 / 100 - 1
            END IF
            IF( i6 == 0 ) n6 = 35
            L(6) = NA( n6 + 1 )
!
            i7 = MOD(i6,100)
            IF(MOD(i7,10) /= 0) THEN
               n7 = i7 / 10
            ELSE
               n7 = i7 / 10 - 1
            END IF
            IF( i7 == 0 ) n7 = 9
            L(7) = NA( n7 + 1 )
!
            n8 = MOD(i7,10)
            IF(n8 == 0) n8 = 10
            L(8) = NA( n8 )
!
         END IF
!
         IF( NLXYZ3 <= 100 ) THEN
!
            L(1) = 'A'
            L(2) = '0'
!
            n3 = i / 100
            IF(MOD(i,10) == 0) n3 = n3 - 1
            L(3) = NA( n3 + 1 )
!
            n4 = MOD(i,10)
            IF(n4 == 0) n4 = 10
            L(4) = NA( n4 )
!
         ELSE
!
            n1 = i / 3600
            IF(MOD(i,3600) == 0) n1 = n1 - 1
            IF(n1 <= 26) THEN
               L(1) = NA( n1 + 11 )
            ELSE
               L(1) = NA( n1 - 25 )
            END IF
!
            i2 = MOD(i,3600)
            IF(MOD(i2,100) /= 0) THEN
               n2 = i2 / 100
            ELSE
               n2 = i2 / 100 - 1
            END IF
            IF( i2 == 0 ) n2 = 35
            L(2) = NA( n2 + 1 )
!
            i3 = MOD(i2,100)
            IF(MOD(i3,10) /= 0) THEN
               n3 = i3 / 10
            ELSE
               n3 = i3 / 10 - 1
            END IF
            IF( i3 == 0 ) n3 = 9
            L(3) = NA( n3 + 1 )
!
            n4 = MOD(i3,10)
            IF(n4 == 0) n4 = 10
            L(4) = NA( n4 )
!
         END IF
!
         n1_adj = ICHAR(L(1)) + letter_shift
         L(1)   = CHAR(n1_adj)
!
         Name_of_8Character_Element = L(1)//L(2)//L(3)//L(4)//L(5)//L(6)//L(7)//L(8)
!
         RETURN
!
      END IF
!
!*********************************************************************
!*                                                                   *
!*                            3-D SYSTEMS                            *
!*                                                                   *
!*********************************************************************
!
      IF( NLXYZ1 <= 100 ) THEN
!
         m1 = 2
         n7 = k / 10
         IF( MOD(k,10) == 0 ) n7 = n7 - 1
         L(7) = NA( n7 + 1)
!
         n8 = MOD(k,10)
         IF(n8 == 0) n8 = 10
         L(8) = NA( n8 )
!
      ELSE IF( NLXYZ1 <= 3600 .AND. NLXYZ1 > 100 ) THEN
!
         m1 = 3
         n6 = k / 100
         IF( MOD(k,100) == 0 ) n6 = n6 - 1
         L(6) = NA( n6 + 1 )
!
         i7 = MOD(k,100)
         IF(MOD(i7,10) /= 0) THEN
            n7 = i7 / 10
         ELSE
            n7 = i7 / 10 - 1
         END IF
         IF( i7 == 0 ) n7 = 9
         L(7) = NA( n7 + 1 )
!
         n8 = MOD(i7,10)
         IF(n8 == 0) n8 = 10
         L(8) = NA( n8 )
!
      END IF
!
! >>>
! >>>
! >>>
!
      IF( NLXYZ2 <= 36 ) THEN
!
         m2 = 1
!
         n8m1 = MOD(j,36)
         IF(n8m1 == 0) n8m1 = 36
         L(8-m1) = NA( n8m1 )

!
      ELSE IF( NLXYZ2 <= 1296 .AND. NLXYZ2 > 36 ) THEN
!
         m2 = 2
!
         n8m2 = j / 36
         IF( MOD(j,36) == 0 ) n8m2 = n8m2 - 1
         L(7-m1) = NA( n8m2 + 1)
!
         n8m1 = MOD(j,36)
         IF(n8m1 == 0) n8m1 = 36
         L(8-m1) = NA( n8m1 )
!
      ELSE IF( NLXYZ2 <= 46656 .AND. NLXYZ2 > 1296 ) THEN
!
         m2 = 3
!
         n8m3 = j / 1296
         IF( MOD(j,1296) == 0 ) n8m3 = n8m3 - 1
         L(6-m1) = NA( n8m3 + 1 )
!
         i8m2 = MOD(j,1296)
         IF(MOD(i8m2,36) /= 0) THEN
            n8m2 = i8m2 / 36
         ELSE
            n8m2 = i8m2 / 36 - 1
         END IF
         IF( i8m2 == 0 ) n8m2 = 35
         L(6-m1) = NA( n8m2 + 1 )
!
         n8m1 = MOD(i8m2,36)
         IF( n8m1 == 0 ) n8m1 = 36
         L(5-m1) = NA( n8m1 )
!
      END IF
!
! >>>
! >>>
! >>>
!
      m3 = 8 - m1 - m2
!
      IF( NLXYZ3 <= 36 ) THEN
!
         ns = MOD(i,36)
         IF(ns == 0) ns = 36
         L(m3) = NA( ns )
!
         L(1) = 'A'
         L(2:m3-1) = '0'
!
      ELSE IF( NLXYZ3 <= 1296 .AND. NLXYZ3 > 36 ) THEN
!
         nsm1 = i / 36
         IF( MOD(i,36) == 0 ) nsm1 = nsm1 - 1
         L(m3-1) = NA( nsm1 + 1)
!
         ns = MOD(i,36)
         IF(ns == 0) ns = 36
         L(m3) = NA( ns )
!
         IF( m3 >  2) L(1) = 'A'
         IF( m3 >= 4) L(2:m3-2) = '0'
!
      ELSE IF( NLXYZ3 <= 46656 .AND. NLXYZ3 > 1296 ) THEN
!
         nsm2 = i / 1296
!
         IF( MOD(i,1296) == 0 ) nsm2 = nsm2 - 1
         IF( m3 == 3 ) THEN
            IF(nsm2 <= 26) THEN
               L(m3-2) = NA( nsm2 + 11 )
            ELSE
               L(m3-2) = NA( nsm2 - 25 )
            END IF
         ELSE
            L(m3-2) = NA( nsm2 + 1)
         END IF
!
         ism1 = MOD(i,1296)
         IF(MOD(ism1,36) /= 0) THEN
            nsm1 = ism1 / 36
         ELSE
            nsm1 = ism1 / 36 - 1
         END IF
         IF( ism1 == 0 ) nsm1 = 35
         L(m3-1) = NA( nsm1 + 1 )
!
         ns = MOD(ism1,36)
         IF( ns == 0 ) ns = 36
         L(m3) = NA( ns )
!
         IF( m3 >  3) L(1) = 'A'
         IF( m3 >= 5) L(2:m3-3) = '0'
!
      END IF
!
      n1_adj = ICHAR(L(1)) + letter_shift
      L(1)   = CHAR(n1_adj)
!
      Name_of_8Character_Element = L(1)//L(2)//L(3)//L(4)//L(5)//L(6)//L(7)//L(8)
!
      RETURN
!
!
!***********************************************************************
!*                                                                     *
!*                        F  O  R  M  A  T  S                          *
!*                                                                     *
!***********************************************************************
!
 6000 FORMAT(/,'Name_of_8Character_Element 1.0 ......... 14 December  2014',6X,'Function returning 8-character names of the grid elements')
!
 8100 FORMAT(//,20('ERROR-'),//,T33,   &
     &             '       S I M U L A T I O N   A B O R T E D',/,   &
     &         T19,'The number of elements along the shortest ',   &
     &             'dimension of a 3-D grid exceeds 3600',   &
     &       /,T33,'               CORRECT AND TRY AGAIN',   &
     &       //,20('ERROR-'))
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of <Name_of_8Character_Element>
!
!
      RETURN
!
      END FUNCTION Name_of_8Character_Element
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!     END MODULE Grid_Generator
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
