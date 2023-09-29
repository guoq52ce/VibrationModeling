      MODULE usr
!
!       Declare the user-defined namelist variables (usrnlst.inc) in this module.
!       Also Include user-defined variables in this module.  To access the
!       variables from a subroutine add the statement "Use usr".  If allocatable
!       arrays are defined in this module allocate them in usr0.  To turn on the
!       user defined subroutines (usr0, usr1, and usr2) set keyword CALL_USR to true.

!**************************************qg***************************************************
            USE GET_STL_DATA_MOD
            use stl_preproc_des
            IMPLICIT NONE
            integer :: is_1_group_id
!**************************************qg***************************************************            

      END MODULE usr
