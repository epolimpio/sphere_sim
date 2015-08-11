subroutine calc_force_elastic(n, pos, dir_vec, v0, f_tot, stress)
!   Calculate the pair forces in the case of only elastic (repulsive) forces
!   
!   INPUT: n-> number of particles
!   INPUT: pos-> (3,n)-array with the coordinates of all the particles
!   INPUT: dir_vec-> (3,n)-array with the direction movement of the particles
!   INPUT: v0-> scalar with the active force parameter
!
!   OUTPUT: f_tot-> (3,n)-array with the coordinates of the total force

    implicit none
    integer(kind=4), intent(in) :: n
    
    real(kind=8), intent(in) :: v0
    real(kind=8), intent(in), dimension(3,n) :: pos
    real(kind=8), intent(in), dimension(3,n) :: dir_vec
    real(kind=8), intent(out), dimension(3,n) :: f_tot
    real(kind=8), intent(out), dimension(9,n) :: stress
    real(kind=8), dimension(3) :: rij
    real(kind=8) :: mod_rij, f
    real(kind=8), dimension(3) :: Fij
    integer(kind=4) :: i, j

!   Start stress matrix
    stress = 0

!   Calculate the active force
    f_tot = v0*dir_vec

!   Sum the elastic force for each pair
    do i = 1,n
        do j = i+1, n
            rij(:) = pos(:,i) - pos(:,j)
            mod_rij = sqrt(dot_product(rij,rij))
            f = 2 - mod_rij
            if (f .gt. 0) then
                Fij = rij/mod_rij*f
                f_tot(:,i) = f_tot(:,i) + Fij(:)
                call add_stress(pos(:,i), rij, Fij, stress(:,i))
                f_tot(:,j) = f_tot(:,j) - Fij(:)
                call add_stress(pos(:,j), -rij, -Fij, stress(:,j))
            end if
        end do
    end do

end subroutine calc_force_elastic

subroutine calc_force_hooke_break(n, n_tri, pos, dir_vec, v0, anisotropy, max_dist, list, f_tot, stress)
!   Calculate the pair forces in the case of only elastic (repulsive) forces
!   
!   INPUT: n-> number of particles
!   INPUT: pos-> (3,n)-array with the coordinates of all the particles
!   INPUT: dir_vec-> (3,n)-array with the direction movement of the particles
!   INPUT: v0-> scalar with the active force parameter
!   INPUT: anisotropy-> ratio of the force in pulling and pushing of the spring (k_pull/k_push)
!   INPUT: max_dist-> distance on which force is zero. If -1 it is neglected
!
!   OUTPUT: f_tot-> (3,n)-array with the coordinates of the total force
!   OUTPUT: stress -> (9,n)-array with all the stress components in apherical coordinates

    implicit none
    integer(kind=4), intent(in) :: n_tri
    integer(kind=4), intent(in) :: n
    
    real(kind=8), intent(in) :: v0
    real(kind=8), intent(in) :: anisotropy
    real(kind=8), intent(in) :: max_dist
    real(kind=8), intent(in), dimension(3,n) :: pos
    real(kind=8), intent(in), dimension(3,n) :: dir_vec
    integer(kind=4), intent(in), dimension(n_tri, 3) :: list
    real(kind=8), intent(out), dimension(3,n) :: f_tot
    real(kind=8), intent(out), dimension(9,n) :: stress
    integer(kind=4), dimension(3*n_tri,2) :: pairs
    real(kind=8), dimension(3) :: rij
    real(kind=8), dimension(3) :: Fij
    real(kind=8) :: mod_rij, f, k
    integer(kind=4) :: i, i1, i2

!   Start stress matrix
    stress = 0

!   Calculate the active force
    f_tot = v0*dir_vec

!   Get the list of all pairs
    call get_all_pairs(n_tri, list, pairs)

!   Sum the Hooke force for each pair
    do i=1,3*n_tri
        i1 = pairs(i,1)
        i2 = pairs(i,2)
        rij(:) = pos(:,i1) - pos(:,i2)
        mod_rij = sqrt(dot_product(rij,rij))
        f = 2-mod_rij
        if (f .lt. 0) then
!       pushing
            k = 1
        else
!       pulling
            k = anisotropy
        end if          
!       check if the distance is above max_dist
        if ((max_dist .gt. 2) .and. (mod_rij .gt. max_dist))  then
            k = 0
        end if
!       correct f and check if it is above 0
        f = k*f
        if (f .gt. 1e-9) then                           
            Fij = rij/mod_rij*(2-mod_rij)
            f_tot(:,i1) = f_tot(:,i1) + Fij(:)
            call add_stress(pos(:,i1), rij, Fij, stress(:,i1))
            f_tot(:,i2) = f_tot(:,i2) - Fij(:)
            call add_stress(pos(:,i2), -rij, -Fij, stress(:,i2))
        end if
    end do
    
end subroutine calc_force_hooke_break

subroutine calc_force_hooke(n, n_tri, pos, dir_vec, v0, list, f_tot, stress)
!   Calculate the pair forces in the case of only elastic (repulsive) forces
!   
!   INPUT: n-> number of particles
!   INPUT: n_tri-> number of triangles
!   INPUT: pos-> (3,n)-array with the coordinates of all the particles
!   INPUT: list-> (n_tri,3)-array with the list of the nodes of the triangles
!   INPUT: dir_vec-> (3,n)-array with the direction movement of the particles
!   INPUT: v0-> scalar with the active force parameter
!
!   OUTPUT: f_tot-> (3,n)-array with the coordinates of the total force

    implicit none
    integer(kind=4), intent(in) :: n_tri
    integer(kind=4), intent(in) :: n
    
    real(kind=8), intent(in) :: v0
    real(kind=8), intent(in), dimension(3,n) :: pos
    real(kind=8), intent(in), dimension(3,n) :: dir_vec
    integer(kind=4), intent(in), dimension(n_tri, 3) :: list
    real(kind=8), intent(out), dimension(3,n) :: f_tot
    real(kind=8), intent(out), dimension(9,n) :: stress
    integer(kind=4), dimension(3*n_tri,2) :: pairs
    real(kind=8), dimension(3) :: rij
    real(kind=8), dimension(3) :: Fij
    real(kind=8) :: mod_rij
    integer(kind=4) :: i, i1, i2

!   Start stress matrix
    stress = 0

!   Calculate the active force
    f_tot = v0*dir_vec

!   Get the list of all pairs
    call get_all_pairs(n_tri, list, pairs)

!   Sum the Hooke force for each pair
    do i=1,3*n_tri
        i1 = pairs(i,1)
        i2 = pairs(i,2)
        rij(:) = pos(:,i1) - pos(:,i2) 
        mod_rij = sqrt(dot_product(rij,rij))
        Fij = rij/mod_rij*(2-mod_rij)
        f_tot(:,i1) = f_tot(:,i1) + Fij(:)
        call add_stress(pos(:,i1), rij, Fij, stress(:,i1))
        f_tot(:,i2) = f_tot(:,i2) - Fij(:)
        call add_stress(pos(:,i2), -rij, -Fij, stress(:,i2))
    end do

end subroutine calc_force_hooke

subroutine add_stress(r, rij, Fij, stress)

    real(kind=8), intent(in), dimension(3) :: r
    real(kind=8), intent(in), dimension(3) :: rij
    real(kind=8), intent(in), dimension(3) :: Fij
    real(kind=8), intent(inout), dimension(9) :: stress
    
    real(kind=8), dimension(3) :: e_r
    real(kind=8), dimension(3) :: e_phi
    real(kind=8), dimension(3) :: e_theta
    real(kind=8) :: Fr, rr
    real(kind=8) :: Fth, rth
    real(kind=8) :: Fphi, rphi

    call get_spherical_unit_vec(r, e_r, e_phi, e_theta)

    ! Calculate the components in spherical coordinates
    Fr = dot_product(e_r, Fij)
    Fth = dot_product(e_theta, Fij)
    Fphi = dot_product(e_phi, Fij)
    rr = dot_product(e_r, rij)
    rth = dot_product(e_theta, rij)
    rphi = dot_product(e_phi, rij)

    ! Add the nine components of the stress
    stress(1) = stress(1) + rr*Fr
    stress(2) = stress(2) + rth*Fr
    stress(3) = stress(3) + rphi*Fr
    stress(4) = stress(4) + rr*Fth
    stress(5) = stress(5) + rth*Fth
    stress(6) = stress(6) + rphi*Fth
    stress(7) = stress(7) + rr*Fphi
    stress(8) = stress(8) + rth*Fphi
    stress(9) = stress(9) + rphi*Fphi

end subroutine add_stress

subroutine get_all_pairs(n_tri, list, pairs)

    implicit none
    integer(kind=4), intent(in) :: n_tri
    
    integer(kind=4), intent(in), dimension(n_tri, 3) :: list
    integer(kind=4), intent(out), dimension(3*n_tri, 2) :: pairs
    
    integer(kind=4) :: i, index

!   For each triangle write each of the 3 pairs
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

subroutine calc_pairs_dist(n_tri, n, pos, list, pairs_dist, pairs)

    implicit none
    integer(kind=4), intent(in) :: n_tri
    integer(kind=4), intent(in) :: n
    
    real(kind=8), intent(in), dimension(3,n) :: pos
    integer(kind=4), intent(in), dimension(n_tri,3) :: list
    real(kind=8), intent(out), dimension(3*n_tri) :: pairs_dist
    integer(kind=4), intent(out), dimension(3*n_tri,2) :: pairs
    
    real(kind=8), dimension(3) :: rij
    integer(kind=4) :: i

    call get_all_pairs(n_tri, list, pairs)

    do i=1,3*n_tri
        rij(:) = pos(:,pairs(i,1)) - pos(:,pairs(i,2)) 
        pairs_dist(i) = sqrt(dot_product(rij,rij))
    end do

end subroutine calc_pairs_dist

subroutine get_spherical_unit_vec(r, e_r, e_phi, e_theta)

    real(kind=8), intent(in), dimension(3) :: r
    real(kind=8), intent(out), dimension(3) :: e_r
    real(kind=8), intent(out), dimension(3) :: e_phi
    real(kind=8), intent(out), dimension(3) :: e_theta
    
    real(kind=8) :: mod_r, theta, phi

    mod_r = sqrt(dot_product(r,r))
    e_r = r/mod_r
    theta = acos(r(3)/mod_r)
    phi = atan2(r(2), r(1))

    ! Calculate e_phi
    e_phi(1) = -sin(phi)
    e_phi(2) = cos(phi)
    e_phi(3) = 0

    ! Calculate e_theta
    e_theta(1) = cos(theta)*cos(phi)
    e_theta(2) = cos(theta)*sin(phi)
    e_theta(3) = -sin(theta)

end subroutine get_spherical_unit_vec

function cross(a, b)
! Calculate the cross product of two vectors a and b in 3 dimensions
  
  integer, dimension(3) :: cross
  integer, dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)

end function cross


