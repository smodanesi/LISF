!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_getpeobspred_SATlaiobs
!  \label{NoahMP401_getpeobspred_SATlaiobs}
!
! !REVISION HISTORY:
! 15 May 2025: Sara Modanesi; Initial Specification modified based on Noah-MP3.6 
! !INTERFACE:
subroutine NoahMP401_getpeobspred_SATlaiobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use NoahMP401_lsmMod !, only : NoahMP36_struc
  use noahmp401_lsmMod
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
     lai(t) = NoahMP401_struc(n)%noahmp401(t)%lai
  enddo

end subroutine NoahMP401_getpeobspred_SATlaiobs



