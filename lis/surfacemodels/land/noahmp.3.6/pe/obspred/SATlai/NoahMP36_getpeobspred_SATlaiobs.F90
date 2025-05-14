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
! !ROUTINE: NoahMP36_getpeobspred_SATlaiobs
!  \label{NoahMP36_getpeobspred_SATlaiobs}
!
! !REVISION HISTORY:
! 21 Jun 2023: Sara Modanesi; Initial Specification modified based on the SMAPsm dir
! 14 May 2025: Sara Modanesi; changed from SYNTlai to SATlai to allow confusion and calibrate with sat dat 
! !INTERFACE:
subroutine NoahMP36_getpeobspred_SATlaiobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use NoahMP36_lsmMod !, only : NoahMP36_struc
  use noahmp36_lsmMod
  use LIS_logMod,       only : LIS_verify, LIS_logunit

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: Obj_Func
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: laiField
  real, pointer          :: lai(:)
  integer                :: t
  integer                :: i
  integer                :: status

  n = 1

  call ESMF_StateGet(Obj_Func,"SAT_lai",laiField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lai(t) = NoahMP36_struc(n)%noahmp36(t)%lai
  enddo

end subroutine NoahMP36_getpeobspred_SATlaiobs



