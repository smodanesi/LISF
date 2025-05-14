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
! !ROUTINE: reset_SATlaiobs
! \label{reset_SATlaiobs}
!
! !REVISION HISTORY:
!  21 June 2023: Sara Modanesi; Initial Specification
!  14 May 2025 : Sara Modanesi; changed specifications from SYNTlai to SATlai (avoid confusion and calibrate with satellite data)
! !INTERFACE: 
subroutine reset_SATlaiobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,        only : LIS_rc
  use SATlai_obsMod,       only : SATlai_obs_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  resets the satellite leaf area index data structure for parameter
!  optimization
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP

  integer            :: n

end subroutine reset_SATlaiobs

