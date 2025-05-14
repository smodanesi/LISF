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
! !ROUTINE: NoahMP36_setupobspred_SATlaiobs
!  \label{NoahMP36_setupobspred_SATlaiobs}
!
! !REVISION HISTORY:
! 21 JUN 2023: Sara Modanesi; Initial Specification modified based on the SMAPsm dir
! 14 May 2025: Sara Modanesi; changed specification from SYNTlai to SATlai to avoid confusion and calibrate with sat obs
! !INTERFACE:
subroutine NoahMP36_setupobspred_SATlaiobs(OBSPred)
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
  type(ESMF_Field)       :: laiField
  integer                :: status

  n = 1
  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

  laiField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecPatch(n,LIS_rc%lsm_index), &
       name="SAT_lai", rc=status)
  call LIS_verify(status)
  
  call ESMF_StateAdd(OBSPred,(/laiField/),rc=status)
  call LIS_verify(status)

end subroutine NoahMP36_setupobspred_SATlaiobs

