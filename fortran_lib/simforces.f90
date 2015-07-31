subroutine calc_force_elastic(n, pos, dir_vec, v0, f_tot)
!	Calculate the pair forces in the case of only elastic (repulsive) forces
!	
!	INPUT: n-> number of particles
!	INPUT: pos-> (3,n)-array with the coordinates of all the particles
!	INPUT: dir_vec-> (3,n)-array with the direction movement of the particles
!	INPUT: v0-> scalar with the active force parameter
!
!	OUTPUT: f_tot-> (3,n)-array with the coordinates of the total force

	implicit none
	integer(kind=4), intent(in) :: n
	real(kind=8), intent(in) :: v0
	real(kind=8), intent(in), dimension(3,n) :: pos
	real(kind=8), intent(in), dimension(3,n) :: dir_vec
	real(kind=8), intent(out), dimension(3,n) :: f_tot
	real(kind=8), dimension(3) :: rij
	real(kind=8) :: mod_rij
	real(kind=8) :: f
	real(kind=8), dimension(3) :: Fij
	integer(kind=4) :: i
	integer(kind=4) :: j

!	Calculate the active force
	f_tot = v0*dir_vec
!	Sum the elastic force for each pair
	do i = 1,n
		do j = i+1, n
			rij(:) = pos(:,i) - pos(:,j)
			mod_rij = sqrt(sum(rij**2))
			f = 2 - mod_rij
			if (f .gt. 0) then
				Fij = rij/mod_rij*f
				f_tot(:,i) = f_tot(:,i) + Fij(:)
				f_tot(:,j) = f_tot(:,j) - Fij(:)
			end if
		end do
	end do

end subroutine calc_force_elastic

subroutine get_all_pairs(n_tri, list, pairs)

	implicit none
	integer(kind=4), intent(in) :: n_tri
	integer(kind=4), intent(in), dimension(n_tri, 3) :: list
	integer(kind=4), intent(out), dimension(3*n_tri, 2) :: pairs
	integer(kind=4) :: i, index

!	For each triangle write each of the 3 pairs
	do i=1,n_tri
		index = 3*(i-1) + 1
		pairs(index,1) = list(i,1)
		pairs(index,2) = list(i,2)
		pairs(index+1,1) = list(i,1)
		pairs(index+1,2) = list(i,3)
		pairs(index+2,1) = list(i,2)
		pairs(index+2,2) = list(i,3)
	end do

end subroutine get_all_pairs

function cross(a, b)
! Calculate the cross product of two vectors a and b in 3 dimensions
  
  integer, dimension(3) :: cross
  integer, dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)

end function cross

