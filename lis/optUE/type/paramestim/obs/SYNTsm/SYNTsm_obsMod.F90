!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: SYNTsm_obsMod
! 
! !DESCRIPTION:  This routine contains interfaces and subroutines to
!   handle SYNTHETIC observations of soil moisture
!
!   
! !REVISION HISTORY: 
! 25 May 2023 Sara Modanesi;   Initial Specification. Synthetic hourly observation 
! 
module SYNTsm_obsMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SYNTsm_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SYNTsm_obs_struc

  type, public ::  SYNTsm_obs_data_dec

     integer             :: smcField
     integer             :: nc,nr
     real,    allocatable    :: smobs(:,:)
     real,    allocatable    :: smtime(:,:)
     logical                 :: startMode 
     type(ESMF_Time)             :: startTime
     type(ESMF_TimeInterval)     :: timestep

  end type SYNTsm_obs_data_dec

  type(SYNTsm_obs_data_dec), allocatable :: SYNTsm_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: SYNTsm_obs_setup
! \label{SYNTsm_obs_setup}
! 
! !INTERFACE: 
  subroutine SYNTsm_obs_setup(Obs_State)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
    use map_utils, only : latlon_to_ij 
    use LIS_timeMgrMod
    use netcdf

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obs_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   The arguments are: 
!   \begin{description}
!    \item[Obs\_State]   observation state object 
!   \end{description}
!EOP
    integer                   ::  n 
    type(ESMF_ArraySpec)      ::  realarrspec !,intarrspec,pertArrSpec
    type(ESMF_Field)          ::  obsField(LIS_rc%nnest)
!    character*100             ::  obsAttribFile(LIS_rc%nnest)
    character*100          ::  obsdir
    character*80           ::  filename

    integer                 :: ftn
    real                   :: dx, dy
    integer                :: NX, NY
    integer                :: ncid
    integer                :: ios
!    character*100          :: infile
    character*100          :: xname, yname
    integer                 :: status


    allocate(SYNTsm_obs_struc(LIS_rc%nnest))

    write(LIS_logunit,*) '[INFO] Setting up SYNT sm data reader....'

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SYNTHETIC soil moisture data directory:",&
         rc=status)

    call ESMF_ConfigGetattribute(LIS_config,obsdir,&
            rc=status)
    call LIS_verify(status,'SYNTHETIC soil moisture data directory: not defined')

    do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            obsdir, rc=status)
       call LIS_verify(status)
       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Update Time",&
            -99.0, rc=status)
       call LIS_verify(status)
    enddo   
 
    write(LIS_logunit,*)'[INFO] read SYNTHETIC soil moisture data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations. 
!   The array size is LIS_rc%ngrid(n). 
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest

       obsField = ESMF_FieldCreate(arrayspec=realarrspec, &
            grid=LIS_vecGrid(n), &
            name="SYNT_sm", rc=status)
       call LIS_verify(status, 'Error in ESMF_FieldCreate: SYNT_sm ')


       call ESMF_StateAdd(Obs_State(n),(/obsField/),rc=status)
       call LIS_verify(status, 'Error in ESMF_StateAdd: obsField')

    enddo
    write(LIS_logunit,*) '[INFO] created the States to hold the SYNTHETIC soil moisture data'
    
!-------------------------------------------------------------
! set up the SYNTHETIC SM domain %and interpolation weights. 
!-------------------------------------------------------------

    ! get nx ny dimensions
    ! 2015 01 01 as example file. If call file gives error stops run
    call SYNTsmobs_filename(filename,obsdir,&
                     2015,1,1,1)
    ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)
    call LIS_verify(ios,'Error reading in SYNTHETIC sm data dimensions: Error opening file'// filename)
    ios = nf90_inquire_dimension(ncid,1,yname,NY)
    ios = nf90_inquire_dimension(ncid,2,xname,NX)
    ios = nf90_close(ncid)

    do n=1,LIS_rc%nnest
       SYNTsm_obs_struc(n)%nc = NX
       SYNTsm_obs_struc(n)%nr = NY

       allocate(SYNTsm_obs_struc(n)%smobs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(SYNTsm_obs_struc(n)%smtime(&
            LIS_rc%lnc(n), LIS_rc%lnr(n)))
       SYNTsm_obs_struc(n)%smobs = LIS_rc%udef
       SYNTsm_obs_struc(n)%smtime = -1 !check

       call LIS_registerAlarm("SYNT sm read alarm",&
            3600.0, 3600.0) !from 86400

       SYNTsm_obs_struc(n)%startMode = .true.
    enddo
  end subroutine SYNTsm_obs_setup
  
end module SYNTsm_obsMod
