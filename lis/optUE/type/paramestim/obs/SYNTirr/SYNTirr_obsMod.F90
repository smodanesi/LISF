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
! !MODULE: SYNTirr_obsMod
! 
! !DESCRIPTION: This routine contains interfaces and subroutines to
!   handle SYNTHETIC IRRIGATION observatio
!
!   
! !REVISION HISTORY: 
!25 May 2023 Sara Modanesi;   Initial Specification. Synthetic observation of
!irrigation every hour
! 
module SYNTirr_obsMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SYNTirr_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SYNTirr_obs_struc

  type, public ::  SYNTirr_obs_data_dec

     integer             :: irrigRateField
     integer             :: nc,nr
     real,    allocatable    :: irr(:,:)
     real,    allocatable    :: irrtime(:,:)
     logical                 :: startMode 

  end type SYNTirr_obs_data_dec

  type(SYNTirr_obs_data_dec), allocatable :: SYNTirr_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: SYNTirr_obs_setup
! \label{SYNTirr_obs_setup}
! 
! !INTERFACE: 
  subroutine SYNTirr_obs_setup(Obs_State)
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
    character*100          :: xname, yname
    integer                 :: status


    allocate(SYNTirr_obs_struc(LIS_rc%nnest))

    write(LIS_logunit,*) '[INFO] Setting up SYNT irr data reader....'

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"SYNTHETIC irrigation data directory:",&
         rc=status)

    call ESMF_ConfigGetattribute(LIS_config,obsdir,&
            rc=status)
    call LIS_verify(status,'SYNTHETIC irrigation data directory: not defined')

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
 
    write(LIS_logunit,*)'[INFO] read SYNTHETIC irrigation data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations. 
!   The array size is LIS_rc%ngrid(n). 
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest

       obsField = ESMF_FieldCreate(arrayspec=realarrspec, &
            grid=LIS_vecGrid(n), &
            name="SYNT_irr", rc=status)
       call LIS_verify(status, 'Error in ESMF_FieldCreate: SYNT_irr ')


       call ESMF_StateAdd(Obs_State(n),(/obsField/),rc=status)
       call LIS_verify(status, 'Error in ESMF_StateAdd: obsField')

    enddo
    write(LIS_logunit,*) '[INFO] created the States to hold the SYNTHETIC irrigation data'
    
!-------------------------------------------------------------
! set up the SYNTHETIC  IRR domain 
!-------------------------------------------------------------

    ! get nx ny dimensions
    ! 2015 01 01 as example file
    call SYNTirrobs_filename(filename,obsdir,&
                     2016,1,1,1)
    ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)
    call LIS_verify(ios,'Error reading in SYNTHETIC IRR data dimensions: Error opening file'// filename)
    ios = nf90_inquire_dimension(ncid,1,yname,NY)
    ios = nf90_inquire_dimension(ncid,2,xname,NX)
    ios = nf90_close(ncid)

    do n=1,LIS_rc%nnest
       SYNTirr_obs_struc(n)%nc = NX
       SYNTirr_obs_struc(n)%nr = NY

       allocate(SYNTirr_obs_struc(n)%irr(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(SYNTirr_obs_struc(n)%irrtime(&
            LIS_rc%lnc(n), LIS_rc%lnr(n)))
       SYNTirr_obs_struc(n)%irr = LIS_rc%udef
       SYNTirr_obs_struc(n)%irrtime = -1

       call LIS_registerAlarm("SYNT irr read alarm",&
            3600.0, 3600.0)

       SYNTirr_obs_struc(n)%startMode = .true.
    enddo
  end subroutine SYNTirr_obs_setup
  
end module SYNTirr_obsMod
