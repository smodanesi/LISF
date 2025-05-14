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
! !ROUTINE: read_SATlaiobs
! \label{read_SATlaiobs}
!
! !REVISION HISTORY:
!  21 Jun 2023  Sara Modanesi;   Initial Specification
!  14 May 2025  Sara Modanesi;   changed specification from SYNTlai to SATlai to avoid confusion and run opt. with satellite LAI 
! !INTERFACE: 
subroutine read_SATlaiobs(Obj_Space) 
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use map_utils
  use LIS_fileIOMod
  use SATlai_obsMod

  implicit none
! !ARGUMENTS: 
  integer             :: n
  type(ESMF_State)    :: Obj_Space
  type(ESMF_Field)              :: laiField
  real,             pointer     :: lai(:)
  logical                       :: data_update
  logical                       :: alarmCheck
  logical                       :: file_exists !,data_upd
  integer                       :: c,r !, p, t
  character*100                 :: obsdir
  character*80                  :: fname
  integer                       :: ftn
  real                          :: lon, lhour
  real                          :: gmt
  real                          :: dt
  integer                       :: zone
  integer                       :: grid_index
  integer                       :: status !, iret, ierr
  !integer                       :: days(12)
  
  n=1
  
  call ESMF_AttributeGet(Obj_Space,"Data Directory",& 
       obsdir, rc=status)
  call LIS_verify(status,'Error in ESMF_AttributeGet: Data Directory')
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status, 'Error in ESMF_AttributeGet: Data Update Status')

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SAT lai read alarm")
  
  if(alarmCheck.or.SATlai_obs_struc(n)%startMode) then
     SATlai_obs_struc(n)%startMode = .false.
     SATlai_obs_struc(n)%laiobs= LIS_rc%udef
     SATlai_obs_struc(n)%laitime = -1

     call SATlaiobs_filename(fname,obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da, LIS_rc%hr)

     inquire(file=fname,exist=file_exists)
     if(file_exists) then

        write(LIS_logunit,*)  '[INFO] Reading SAT lai data', fname
        call read_SATlaiobs_data(n, fname,SATlai_obs_struc(n)%laiobs)

!-------------------------------------------------------------------------
! Store the LAI time (modified from Sentinel1 backscatter DA)
!-------------------------------------------------------------------------

        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then
                 grid_index = c+(r-1)*LIS_rc%lnc(n)
                 lon = LIS_domain(n)%lon(grid_index)

                 !lhour = 11.0 
                 !call LIS_localtime2gmt(gmt,lon,lhour,zone)
                 !SATlai_obs_struc(n)%laitime(c,r) = gmt
                 !or another option
                 SATlai_obs_struc(n)%laitime(c,r) = LIS_rc%hr
              endif
           enddo
        enddo
     endif
  endif
  call ESMF_StateGet(Obj_Space,"SAT_lai",laiField,&
          rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: SAT_lai')

  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status, 'Error in ESMF_FieldGet: laiField')

  lai = LIS_rc%udef
!-------------------------------------------------------------------------
!  Update the OBJ_Space (OBJECTIVE SPACE) by subsetting to the local grid time  
!-------------------------------------------------------------------------     

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then
           grid_index = c+(r-1)*LIS_rc%lnc(n)

           dt = (LIS_rc%gmt - SATlai_obs_struc(n)%laitime(c,r))*3600.0
           lon = LIS_domain(n)%lon(grid_index)

           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then
              lai(LIS_domain(n)%gindex(c,r)) = &
                   SATlai_obs_struc(n)%laiobs(c,r)
           endif
        endif
     enddo
  enddo
  
  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status, 'Error in ESMF_AttributeSet: Data Update Status')

end subroutine read_SATlaiobs

!BOP
!
! !ROUTINE: read_SATlaiobs_data
! \label{read_SATlaiobs_data}
!
! !INTERFACE:
subroutine read_SATlaiobs_data(n, fname, lai_ip)
!
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use map_utils,    only : latlon_to_ij
  use SATlai_obsMod, only : SATlai_obs_struc

  implicit none
!
! !INPUT PARAMETERS:
!
  integer                       :: n
  character (len=*)             :: fname
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine reads the Synthetic LAI  NETCDF files
!
!  The arguments are:
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SAT lai file
!  \item[lai\_ip]   lai data simulated in LIS (same LIS domain)
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
  real                        :: leafareaindex(SATlai_obs_struc(n)%nr,SATlai_obs_struc(n)%nc)
  real                        :: lat_nc(SATlai_obs_struc(n)%nr)
  real                        :: lat_nc_fold(SATlai_obs_struc(n)%nr)
  real                        :: lon_nc(SATlai_obs_struc(n)%nc)
  real                        :: lai_ip(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                     :: nlai_ip(LIS_rc%lnc(n),LIS_rc%lnr(n))
  logical                     :: file_exists
  integer                     :: c,r,i,j
  integer                     :: stn_col,stn_row
  real                        :: col,row
  integer                     :: nid
  integer                     :: laiId,latId,lonId
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=fname, exist=file_exists)
  if(file_exists) then

     write(LIS_logunit,*) 'Reading ',trim(fname)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(fname))

     ! variables
     ios = nf90_inq_varid(nid, 'lai',laiid)
     call LIS_verify(ios, 'Error nf90_inq_varid: SAT lai data')

     ios = nf90_inq_varid(nid, 'lat',latid)
     call LIS_verify(ios, 'Error nf90_inq_varid: latitude data')

     ios = nf90_inq_varid(nid, 'lon',lonid)
     call LIS_verify(ios, 'Error nf90_inq_varid: longitude data')

     !values
     ios = nf90_get_var(nid, laiid, leafareaindex )
     call LIS_verify(ios, 'Error nf90_get_var: leafareaindex')

     ios = nf90_get_var(nid, latid, lat_nc_fold)
     call LIS_verify(ios, 'Error nf90_get_var: lat')

     !flipped lat variable. LAI: LIS lat is from lower to larger lat
     do i=1,size(lat_nc_fold)
       lat_nc(i)=lat_nc_fold(size(lat_nc_fold)+1-i)
     end do
     ios = nf90_get_var(nid, lonid, lon_nc)
     call LIS_verify(ios, 'Error nf90_get_var: lon')

     ! close file
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))

!    ! Initialize 
     lai_ip = 0
     nlai_ip = 0

     ! Interpolate the data by averaging 
     do i=1,SATlai_obs_struc(n)%nr
        do j=1,SATlai_obs_struc(n)%nc

           call latlon_to_ij(LIS_domain(n)%lisproj,&
                !lat_nc(i),lon_nc(j),col,row)
                lat_nc(SATlai_obs_struc(n)%nr-(i-1)),lon_nc(j),col,row)
           stn_col = nint(col)
           stn_row = nint(row)

           if(leafareaindex(i,j).ge.-999.and.&
                stn_col.gt.0.and.stn_col.le.LIS_rc%lnc(n).and.&
                stn_row.gt.0.and.stn_row.le.LIS_rc%lnr(n)) then
              ! SM: part of the code from Sentinel-1 DA
              lai_ip(stn_col,stn_row) = lai_ip(stn_col,stn_row) + leafareaindex(i,j)
              nlai_ip(stn_col,stn_row) = nlai_ip(stn_col,stn_row) + 1
           endif
        enddo
     enddo
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(nlai_ip(c,r).ne.0) then
              ! average
              lai_ip(c,r) = lai_ip(c,r)/nlai_ip(c,r)
           else
              lai_ip(c,r) = LIS_rc%udef
           endif
        enddo
     enddo
  endif
!     ! Fix stripes in obs caused by regridding (missing obs for LIS grid cell)
!     ----REMOVED-----LOOK AT THE S1 READER IN CASE YOU NEED IT

#endif

end subroutine read_SATlaiobs_data



! ! INTERFACE:
subroutine SATlaiobs_filename(filename, ndir, yr, mo, da, hr)

  implicit none
! !ARGUMENTS: 
  character*80      :: filename
  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped SATELLITE leafareaindex filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the SAT lai filename
!  \item[ndir] name of the SAT lai root directory
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
  filename =trim(ndir)//'/SAT_LAI_'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//'.nc'

end subroutine SATlaiobs_filename
     
