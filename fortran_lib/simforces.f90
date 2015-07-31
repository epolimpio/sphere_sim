subroutine calc_force_elastic(n, pos, dir_vec, f_tot, v0, k, sigma)
!	Calculate the pair forces in the case
!	of only elastic (repulsive) forces
!
!	
	implicit none
	integer(kind=4), intent(in) :: n
	real(kind=8), intent(in) :: v0
	real(kind=8), intent(in) :: k
	real(kind=8), intent(in) :: sigma
	real(kind=8), intent(in), dimension(3,n) :: pos
	real(kind=8), intent(in), dimension(3,n) :: dir_vec
	real(kind=8), intent(out), dimension(3,n) :: f_tot
	real(kind=8), dimension(3) :: rij
	real(kind=8) :: mod_rij
	real(kind=8) :: f
	real(kind=8), dimension(3) :: Fij
	integer(kind=4) :: i
	integer(kind=4) :: j

	f_tot = v0*dir_vec
	do i = 1,n
		do j = i+1, n
			rij(:) = pos(:,i) - pos(:,j)
			mod_rij = sqrt(sum(rij**2))
			f = 2*sigma - mod_rij
			if (f .gt. 0) then
				Fij = k*rij/mod_rij*f
				f_tot(:,i) = f_tot(:,i) + Fij(:)
				f_tot(:,j) = f_tot(:,j) - Fij(:)
			end if
		end do
	end do

end