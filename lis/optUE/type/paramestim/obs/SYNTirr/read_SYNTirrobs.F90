!-----------------------BEGIN NOTICE -- DO NOT
!EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT
!EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_SYNTirrobs
! \label{read_SYNTirrobs}
!
! !REVISION HISTORY:
!  25 May 2023  Sara Modanesi;   Initial Specification
!  21 Aug 2023 Sara Modanesi; added configuration to run the code with hourly
!  data

! !INTERFACE: 
subroutine read_SYNTirrobs(Obj_Space) 
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use map_utils
  use LIS_fileIOMod
  use SYNTirr_obsMod

  implicit none
! !ARGUMENTS: 
  integer             :: n
  type(ESMF_State)    :: Obj_Space

  type(ESMF_Field)              :: irrigRateField
  logical                       :: data_update
  logical                       :: alarmCheck
  logical                       :: file_exists !,data_upd
  integer                       :: c,r !, p, t
  character*100                 :: obsdir
  character*80                  :: IRR_filename
  integer                       :: ftn
  real                          :: lon, lhour
  real                          :: gmt
  real                          :: dt
  integer                       :: zone
  integer                       :: grid_index
  real,             pointer     :: irrigrate(:)
  integer                       :: status !, iret, ierr
  integer                       :: days(12)
  
  n=1
  
  call ESMF_AttributeGet(Obj_Space,"Data Directory",& 
       obsdir, rc=status)
  call LIS_verify(status,'Error in ESMF_AttributeGet: Data Directory')
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status, 'Error in ESMF_AttributeGet: Data Update Status')

!-------------------------------------------------------------------------
!   Read the data 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SYNT irr read alarm")
  
  if(alarmCheck.or.SYNTirr_obs_struc(n)%startMode) then
     SYNTirr_obs_struc(n)%startMode = .false.
     SYNTirr_obs_struc(n)%irr= LIS_rc%udef
     SYNTirr_obs_struc(n)%irrtime = -1

     call SYNTirrobs_filename(IRR_filename,obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da, LIS_rc%hr)

     inquire(file=IRR_filename,exist=file_exists)
     if(file_exists) then

        write(LIS_logunit,*)  '[INFO] Reading SYNT IRR data', IRR_filename
        call read_SYNTirrobs_data(n, IRR_filename,SYNTirr_obs_struc(n)%irr) !removed k index

!-------------------------------------------------------------------------
! Store the IRR time
!-------------------------------------------------------------------------

        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then
                 grid_index = c+(r-1)*LIS_rc%lnc(n)
                 lon = LIS_domain(n)%lon(grid_index)

                 !lhour = 9.0 !from 6.0 am to 10.0 am
                 !call LIS_localtime2gmt(gmt,lon,lhour,zone)
                 !SYNTirr_obs_struc(n)%irrtime(c,r) = gmt
                 SYNTirr_obs_struc(n)%irrtime(c,r) = LIS_rc%hr
              endif
           enddo
        enddo
     endif
  endif
  call ESMF_StateGet(Obj_Space,"SYNT_irr",irrigRateField,&
          rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: SYNT_irr')

  call ESMF_FieldGet(irrigRateField,localDE=0,farrayPtr=irrigrate,rc=status)
  call LIS_verify(status, 'Error in ESMF_FieldGet: irrigRateField')

  irrigrate = LIS_rc%udef
!-------------------------------------------------------------------------
!  Update the OBJ_Space (OBJECTIVE SPACE) by subsetting to the local grid time  
!-------------------------------------------------------------------------     

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then
           grid_index = c+(r-1)*LIS_rc%lnc(n)

           dt = (LIS_rc%gmt - SYNTirr_obs_struc(n)%irrtime(c,r))*3600.0
           lon = LIS_domain(n)%lon(grid_index)

           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then
              irrigrate(LIS_domain(n)%gindex(c,r)) = &
                   SYNTirr_obs_struc(n)%irr(c,r)
           endif
        endif
     enddo
  enddo
  
  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status, 'Error in ESMF_AttributeSet: Data Update Status')

end subroutine read_SYNTirrobs

!BOP
!
! !ROUTINE: read_SYNTirrobs_data
! \label{read_SYNTirrobs_data}
!
! !INTERFACE:
subroutine read_SYNTirrobs_data(n, fname, irr_ip)   !(n, k, fname, irr_ip)
!
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use map_utils,    only : latlon_to_ij
  use SYNTirr_obsMod, only : SYNTirr_obs_struc

  implicit none
!
! !INPUT PARAMETERS:
!
  integer                       :: n
!  integer                       :: k
  character (len=*)             :: fname
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine reads the SYNTHETIC IRRIG. NETCDF files
!
!  The arguments are:
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SYNT IRRIG. file
!  \item[irr\_ip]   irrigation data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
  real                        :: irrigation(SYNTirr_obs_struc(n)%nr,SYNTirr_obs_struc(n)%nc)
  real                        :: lat_nc(SYNTirr_obs_struc(n)%nr)
  real                        :: lat_nc_fold(SYNTirr_obs_struc(n)%nr)
  real                        :: lon_nc(SYNTirr_obs_struc(n)%nc)
  real                        :: irr_ip(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                     :: nirr_ip(LIS_rc%lnc(n),LIS_rc%lnr(n))
  logical                     :: file_exists
  integer                     :: c,r,i,j
  integer                     :: stn_col,stn_row
  real                        :: col,row
  integer                     :: nid
  integer                     :: irrId,latId,lonId
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=fname, exist=file_exists)
  if(file_exists) then

     write(LIS_logunit,*) 'Reading ',trim(fname)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(fname))

     ! variables
     ios = nf90_inq_varid(nid, 'irr',irrid)
     call LIS_verify(ios, 'Error nf90_inq_varid: SYNT IRR data')

     ios = nf90_inq_varid(nid, 'lat',latid)
     call LIS_verify(ios, 'Error nf90_inq_varid: latitude data')

     ios = nf90_inq_varid(nid, 'lon',lonid)
     call LIS_verify(ios, 'Error nf90_inq_varid: longitude data')

     !values
     ios = nf90_get_var(nid, irrid, irrigation)
     call LIS_verify(ios, 'Error nf90_get_var: irrigation')

     ios = nf90_get_var(nid, latid, lat_nc_fold)
     call LIS_verify(ios, 'Error nf90_get_var: lat')

     !flipped lat variable. SM: LIS lat is from lower to larger lat
     do i=1,size(lat_nc_fold)
       lat_nc(i)=lat_nc_fold(size(lat_nc_fold)+1-i)
     end do
     ios = nf90_get_var(nid, lonid, lon_nc)
     call LIS_verify(ios, 'Error nf90_get_var: lon')

     ! close file
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))

!    ! Initialize 
     irr_ip = 0
     nirr_ip = 0

     ! Interpolate the data by averaging 
     do i=1,SYNTirr_obs_struc(n)%nr
        do j=1,SYNTirr_obs_struc(n)%nc

           call latlon_to_ij(LIS_domain(n)%lisproj,&
                !lat_nc(i),lon_nc(j),col,row)
                lat_nc(SYNTirr_obs_struc(n)%nr-(i-1)),lon_nc(j),col,row)
           stn_col = nint(col)
           stn_row = nint(row)

           if(irrigation(i,j).ge.-999.and.&
                stn_col.gt.0.and.stn_col.le.LIS_rc%lnc(n).and.&
                stn_row.gt.0.and.stn_row.le.LIS_rc%lnr(n)) then
              ! SM: part of the code from Sentinel-1 DA
              irr_ip(stn_col,stn_row) = irr_ip(stn_col,stn_row) + irrigation(i,j)
              nirr_ip(stn_col,stn_row) = nirr_ip(stn_col,stn_row) + 1
           endif
        enddo
     enddo
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(nirr_ip(c,r).ne.0) then
              ! average
              irr_ip(c,r) = irr_ip(c,r)/nirr_ip(c,r)
           else
              irr_ip(c,r) = LIS_rc%udef
           endif
        enddo
     enddo
  endif
!     ! Fix stripes in obs caused by regridding (missing obs for LIS grid cell)
!     ----REMOVED-----LOOK AT THE S1 READER IN CASE YOU NEED IT

#endif

end subroutine read_SYNTirrobs_data



! ! INTERFACE:
subroutine SYNTirrobs_filename(filename, ndir, yr, mo, da,hr)

  implicit none
! !ARGUMENTS: 
  character*80      :: filename
  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped SYNT IRR filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the SYNT IRR filename
!  \item[ndir] name of the SYNT IRR root directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  filename = trim(ndir)//'/SYNT_IRR_'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//'.nc'

end subroutine SYNTirrobs_filename
  
