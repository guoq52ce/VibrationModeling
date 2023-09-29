!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USR0

      USE USR
!**************************************qg***************************************************
      USE GET_STL_DATA_MOD
      use stl_preproc_des
!**************************************qg***************************************************      


      IMPLICIT NONE

!**************************************qg***************************************************
!  Define local variables here

      CALL GET_STL_DATA_IS
      CALL COPY_INIT_STL_DATA

! Get group ID of the stl files and asign a BC_ID to them
! The BC_ID must be diffferent from existing BCs
      is_1_group_id = GET_GROUP_ID('is_0001.stl')

! Not needed       
      CALL ASSIGN_BC_ID_TO_STL_GROUP(is_1_group_id, 11)

! Bin the STL to the DES grid (brute force, this could be optimizedto only rebin
! the moving facets).
      CALL BIN_FACETS_TO_DG
!**************************************qg***************************************************

      RETURN
      END SUBROUTINE USR0
