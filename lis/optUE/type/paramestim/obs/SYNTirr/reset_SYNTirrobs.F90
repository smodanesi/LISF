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
! !ROUTINE: reset_SYNTirrobs
! \label{reset_SYNTirrobs}
!
! !REVISION HISTORY:
!  7 Nov 2022: Sara Modanesi; Initial Specification
!
! !INTERFACE: 
subroutine reset_SYNTirrobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,        only : LIS_rc
  use SYNTirr_obsMod,       only : SYNTirr_obs_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  resets the SYNT irr data structures for parameter
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

end subroutine reset_SYNTirrobs

