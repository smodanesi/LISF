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
! !ROUTINE: NoahMP36_getpeobspred_SYNTsmobs
!  \label{NoahMP36_getpeobspred_SYNTsmobs}
!
! !REVISION HISTORY:
! 25 May 2023: Sara Modanesi; Initial Specification modified from SMAPsm dir
!
! !INTERFACE:
subroutine NoahMP36_getpeobspred_SYNTsmobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_soilsMod,  only : LIS_soils
  use NoahMP36_lsmMod, only : NoahMP36_struc
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
  type(ESMF_Field)       :: smcField
  real, pointer          :: smc(:)
  integer                :: t
  integer                :: i
  integer                :: status

!  write(LIS_logunit,*) '[INFO] Here 1 in NoahMP36_getpeobspred_SMAPsmobs '

  n = 1

  call ESMF_StateGet(Obj_Func,"SYNT_sm",smcField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=smc,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     smc(t) = NoahMP36_struc(n)%noahmp36(t)%smc(1)
  enddo

end subroutine NoahMP36_getpeobspred_SYNTsmobs



