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
! !ROUTINE: NoahMP36_getpeobspred_SYNTirrobs
!  \label{NoahMP36_getpeobspred_SYNTirrobs}
!
! !REVISION HISTORY:
! 07 Nov 2022: Sara Modanesi; Initial Specifications
! !INTERFACE:
subroutine NoahMP36_getpeobspred_SYNTirrobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use NoahMP36_lsmMod
  use LIS_irrigationMod
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
  type(ESMF_Field)       :: irrigRateField
  !real, pointer          :: irrigRate(:) !changed from irrigrate to irrigRate
  real, pointer          :: irrigcum(:)
  integer                :: t
  integer                :: i
  integer                :: status
  ! SM 21/09/2023 add variable needed to create cumulated sum over specific tme
  ! interval
  real                   :: IRR_movsum
  integer                :: itemp

!  write(LIS_logunit,*) '[INFO] Here 1 in NoahMP36_getpeobspred_SYNTirr '

  n = 1

  call ESMF_StateGet(Obj_Func,"SYNT_irr",irrigRateField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(irrigRateField,localDE=0,farrayPtr=irrigcum,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     do itemp=NOAHMP36_struc(n)%Tcal_windowsize,2,-1
        NOAHMP36_struc(n)%noahmp36(t)%IRR_ac_antecedent(itemp)=NOAHMP36_struc(n)%noahmp36(t)%IRR_ac_antecedent(itemp-1) 
     enddo
     NOAHMP36_struc(n)%noahmp36(t)%IRR_ac_antecedent(1)=NOAHMP36_struc(n)%noahmp36(t)%irr
     IRR_movsum = 0.0
     do itemp=1,NOAHMP36_struc(n)%Tcal_windowsize
        IRR_movsum=IRR_movsum+NOAHMP36_struc(n)%noahmp36(t)%IRR_ac_antecedent(itemp)
     enddo

     irrigcum(t) = IRR_movsum
  enddo
!  write(LIS_logunit,*) '[INFO] Finished NoahMP36_getpeobspred_SYNTirr, irr  = ', irr

end subroutine NoahMP36_getpeobspred_SYNTirrobs



