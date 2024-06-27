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
! !ROUTINE: NoahMP36_setupobspred_SYNTirrobs
!  \label{NoahMP36_setupobspred_SYNTirrobs}
!
! !REVISION HISTORY:
! 07 Nov 2020: Sara Modanesi; Initial Specification
!
! !INTERFACE:
subroutine NoahMP36_setupobspred_SYNTirrobs(OBSPred)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_vecPatch
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: OBSPred
!
! !DESCRIPTION:
!  
!  This routine creates an entry in the Obs pred object used for 
!  parameter estimation
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: irrigRateField
  integer                :: status

  n = 1
  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

  irrigRateField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecPatch(n,LIS_rc%lsm_index), &
       name="SYNT_irr", rc=status)
  call LIS_verify(status)
  
  call ESMF_StateAdd(OBSPred,(/irrigRateField/),rc=status)
  call LIS_verify(status)
end subroutine NoahMP36_setupobspred_SYNTirrobs

