subroutine setup_num_levels(icv_fine, jcv_fine, num_level)

! argument variables
integer, intent(in):: icv_fine, jcv_fine
integer, intent(out):: num_level

! local variables

integer:: i, ncv

ncv = jcv_fine
i = 1

do

    if( mod(ncv,2) /= 0) exit
    if( ncv/2 < 16 ) exit

    ncv = ncv/2
    i = i+1

end do

num_level = i

end subroutine setup_num_levels
