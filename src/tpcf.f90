module tpcf
use, intrinsic :: ISO_FORTRAN_ENV
!use opm_lib
use kdtree_module
use kmeans_module
use utils

implicit none

private
public shear_2pcf, kappa_2pcf, simple_bootstrap_foreground,  pair_weighted_block_bootstrap, deallocate_lensing_survey, allocate_shear_2pcf_results, deallocate_shear_2pcf_results

!public simple_bootstrap_background, fixed_block_bootstrap, randomize_shear

public shear_2pcf_results, marks_field_collection, shape_catalog, kappa_map, lens_catalog, lensing_survey_field, lensing_survey, equal_weight_blocks, rectangular_blocks, voronoi_blocks, survey_marked_point_block_bootstrap, survey_marked_point_bootstrap

type shear_2pcf_results
    integer :: n_bin
    real, allocatable, dimension(:) :: theta, mean_theta, g_t, g_x, weights, marks_weight
    real, allocatable, dimension(:,:) :: xi_E_marks, xi_B_marks, K_marks, w_marks
    integer(kind=8), allocatable, dimension(:) :: n_pairs
end type shear_2pcf_results

type shape_catalog
    integer :: n
    real, allocatable, dimension(:) :: x, y, z, e1, e2, m, w
end type shape_catalog

type kappa_map
    integer :: n
    real, allocatable, dimension(:) :: x, y, kappa_E, kappa_B, w
end type kappa_map

type lens_catalog
    integer :: n, n_bin
    real, allocatable, dimension(:) :: x, y, z, val, w, marks_weight
end type lens_catalog

type lensing_survey_field
    type(shape_catalog) :: shapes
    type(lens_catalog) :: lenses
    type(kappa_map) :: convergence
end type lensing_survey_field

type lensing_survey
    integer :: n_field
    type(lensing_survey_field), allocatable, dimension(:) :: fields
end type lensing_survey

type marks_field
    integer :: n, n_bin
    real, allocatable, dimension(:) :: x, y, w
    real, allocatable, dimension(:,:) :: xi_E_marks, xi_B_marks, K_marks, w_marks
    integer, allocatable, dimension(:) :: block_indicies
end type marks_field

type marks_field_collection
    integer :: n_field, n_marks
    type(marks_field), allocatable, dimension(:) :: fields
end type marks_field_collection

integer, parameter :: equal_weight_blocks = 1, rectangular_blocks = 2, voronoi_blocks = 3

contains

    subroutine deallocate_lensing_survey(survey)
        !Arguments
        type(lensing_survey), intent(inout) :: survey
        !Variables
        integer :: i

        if(allocated(survey%fields)) then
            do i=1, size(survey%fields)
                if(allocated(survey%fields(i)%shapes%x)) deallocate(survey%fields(i)%shapes%x)
                if(allocated(survey%fields(i)%shapes%y)) deallocate(survey%fields(i)%shapes%y)
                if(allocated(survey%fields(i)%shapes%z)) deallocate(survey%fields(i)%shapes%z)
                if(allocated(survey%fields(i)%shapes%e1)) deallocate(survey%fields(i)%shapes%e1)
                if(allocated(survey%fields(i)%shapes%e2)) deallocate(survey%fields(i)%shapes%e2)
                if(allocated(survey%fields(i)%shapes%m)) deallocate(survey%fields(i)%shapes%m)
                if(allocated(survey%fields(i)%shapes%w)) deallocate(survey%fields(i)%shapes%w)

                if(allocated(survey%fields(i)%convergence%x)) deallocate(survey%fields(i)%convergence%x)
                if(allocated(survey%fields(i)%convergence%y)) deallocate(survey%fields(i)%convergence%y)
                if(allocated(survey%fields(i)%convergence%kappa_E)) deallocate(survey%fields(i)%convergence%kappa_E)
                if(allocated(survey%fields(i)%convergence%kappa_B)) deallocate(survey%fields(i)%convergence%kappa_B)
                if(allocated(survey%fields(i)%convergence%w)) deallocate(survey%fields(i)%convergence%w)

                if(allocated(survey%fields(i)%lenses%x)) deallocate(survey%fields(i)%lenses%x)
                if(allocated(survey%fields(i)%lenses%y)) deallocate(survey%fields(i)%lenses%y)
                if(allocated(survey%fields(i)%lenses%z)) deallocate(survey%fields(i)%lenses%z)
                if(allocated(survey%fields(i)%lenses%val)) deallocate(survey%fields(i)%lenses%val)
                if(allocated(survey%fields(i)%lenses%marks_weight)) deallocate(survey%fields(i)%lenses%marks_weight)
                if(allocated(survey%fields(i)%lenses%w)) deallocate(survey%fields(i)%lenses%w)
            end do
            deallocate(survey%fields)
        end if
    end subroutine deallocate_lensing_survey

    subroutine allocate_shear_2pcf_results(results, n_bin, n_A)
        !Arguments
        type(shear_2pcf_results), intent(out) :: results
        integer, intent(in) :: n_bin, n_A

        if(.not. allocated(results%theta)) allocate(results%theta(n_bin), source=0.0)
        if(.not. allocated(results%mean_theta)) allocate(results%mean_theta(n_bin), source=0.0)
        if(.not. allocated(results%g_t)) allocate(results%g_t(n_bin), source=0.0)
        if(.not. allocated(results%g_x)) allocate(results%g_x(n_bin), source=0.0)
        if(.not. allocated(results%weights)) allocate(results%weights(n_bin), source=0.0)
        if(.not. allocated(results%n_pairs)) allocate(results%n_pairs(n_bin), source=int(0, kind=8))
        if(.not. allocated(results%marks_weight)) allocate(results%marks_weight(n_A), source=0.0)
        if(.not. allocated(results%xi_E_marks) .and. n_A > 0) allocate(results%xi_E_marks(n_bin, n_A), source=0.0)
        if(.not. allocated(results%xi_B_marks) .and. n_A > 0) allocate(results%xi_B_marks(n_bin, n_A), source=0.0)
        if(.not. allocated(results%K_marks) .and. n_A > 0) allocate(results%K_marks(n_bin, n_A), source=0.0)
        if(.not. allocated(results%w_marks) .and. n_A > 0) allocate(results%w_marks(n_bin, n_A), source=0.0)
    end subroutine allocate_shear_2pcf_results

    subroutine deallocate_shear_2pcf_results(results)
        !Arguments
        type(shear_2pcf_results), intent(out) :: results

        if(allocated(results%theta)) deallocate(results%theta)
        if(allocated(results%mean_theta)) deallocate(results%mean_theta)
        if(allocated(results%g_t)) deallocate(results%g_t)
        if(allocated(results%g_x)) deallocate(results%g_x)
        if(allocated(results%weights)) deallocate(results%weights)
        if(allocated(results%n_pairs)) deallocate(results%n_pairs)
        if(allocated(results%marks_weight)) deallocate(results%marks_weight)
        if(allocated(results%xi_E_marks)) deallocate(results%xi_E_marks)
        if(allocated(results%xi_B_marks)) deallocate(results%xi_B_marks)
        if(allocated(results%K_marks)) deallocate(results%K_marks)
        if(allocated(results%w_marks)) deallocate(results%w_marks)
    end subroutine deallocate_shear_2pcf_results

    subroutine shear_2pcf(field, &
                          bins, &
                          results, &
                          spherical_coords, &
                          num_threads, use_kd_tree, verbose)
        !Arguments
        type(lensing_survey_field), intent(in) :: field
        type(binning), intent(in) :: bins
        type(shear_2pcf_results) :: results
        logical, intent(in) :: spherical_coords, use_kd_tree, verbose
        integer, intent(in) :: num_threads
        !Variables
        logical :: calculate_marks

        call allocate_shear_2pcf_results(results, bins%n_bin, field%lenses%n)

        results%n_bin = bins%n_bin
        results%theta = bins%bin_centers

        call tangential_shear_tpcf(field%lenses%x, field%lenses%y, field%lenses%z, field%lenses%val, field%lenses%w,&
                                   field%shapes%x, field%shapes%y, field%shapes%z, field%shapes%e1, field%shapes%e2, field%shapes%w, field%shapes%m, &
                                   bins, &
                                   results%mean_theta, results%g_t, results%g_x, results%weights, results%n_pairs, &
                                   spherical_coords, num_threads, use_kd_tree, verbose, &
                                    calculate_marks, results%xi_E_marks, results%xi_B_marks, results%K_marks, results%w_marks, results%marks_weight)

    end subroutine shear_2pcf

    subroutine kappa_2pcf(field, &
                          bins, &
                          results, &
                          spherical_coords, &
                          num_threads, use_kd_tree, verbose)
        !Arguments
        type(lensing_survey_field), intent(in) :: field
        type(binning), intent(in) :: bins
        type(shear_2pcf_results) :: results
        logical, intent(in) :: spherical_coords, use_kd_tree, verbose
        integer, intent(in) :: num_threads
        !Variables
        logical :: calculate_marks

        call allocate_shear_2pcf_results(results, bins%n_bin, field%lenses%n)

        results%n_bin = bins%n_bin
        results%theta = bins%bin_centers

        call scalar_kappa_tpcf(field%lenses%x, field%lenses%y, field%lenses%val, field%lenses%w, &
                               field%convergence%x, field%convergence%y, field%convergence%kappa_E, field%convergence%kappa_B, field%convergence%w, &
                               bins, &
                               results%mean_theta, results%g_t, results%g_x, results%weights, results%n_pairs, &
                               spherical_coords, num_threads, use_kd_tree, verbose, &
                               calculate_marks, results%xi_E_marks, results%xi_B_marks, results%K_marks, results%w_marks, results%marks_weight)

    end subroutine kappa_2pcf

    subroutine tangential_shear_tpcf(x_A, y_A, z_A, val_A, weight_A,&
                                     x_B, y_B, z_B, e1_B, e2_B, weight_B, m_B, &
                                     bins, &
                                     avg_theta, g_t, g_x, weights, n_pairs, &
                                     spherical_coords, &
                                     num_threads, use_kd_tree, verbose, &
                                     calculate_marks, xi_E_marks, xi_B_marks, K_marks, w_marks, marks_weight)
        !Arguments
        real, dimension(:), intent(in) :: x_A, y_A, z_A, val_A, weight_A, x_B, y_B, z_B, e1_B, e2_B, weight_B, m_B
        type(binning), intent(in) :: bins
        real, dimension(bins%n_bin), intent(out) :: avg_theta, g_t, g_x, weights
        integer(kind=8), dimension(bins%n_bin), intent(out) :: n_pairs
        logical, intent(in) :: spherical_coords, use_kd_tree, verbose, calculate_marks
        integer, intent(in) :: num_threads
        real, dimension(bins%n_bin, size(x_A)), intent(out) :: xi_E_marks, xi_B_marks, K_marks, w_marks
        real, dimension(size(x_A)), intent(out) :: marks_weight
        !Variables
        real, dimension(bins%n_bin) :: K
        real :: distance, delta_x, delta_y, cos2phi, sin2phi, x_foreground, y_foreground, z_foreground, value_foreground, weight_foreground, x_background, y_background, z_background, e1, e2, w, m, x_range, y_range, point_density
        real, dimension(2) :: coord_A, cos_sin_theta
        integer :: n_A, n_B, i, j, index, bin_index, progress_counter, last_progress_point, approx_points_in_theta_max_ball, n_nearest_found, size_kd_tree_results

        real, allocatable, dimension(:,:) :: tree_data
        type(kdtree) :: kd_tree
        integer, allocatable, dimension(:) :: kd_tree_results

        n_A = size(x_A)
        n_B = size(x_B)

        x_range = maxval(x_B) - minval(x_B)
        y_range = maxval(y_B) - minval(y_B)
        point_density = n_B/(x_range*y_range)
        approx_points_in_theta_max_ball = pi*bins%x_max**2 * point_density
        size_kd_tree_results = approx_points_in_theta_max_ball*10.0

        if(verbose) then
            if(use_kd_tree) then
                print *,"Using kd tree."
            else
                print *,"Not using kd tree."
            end if

            print *, "Number of foreground objects:", n_A
            print *, "Number of background objects:", n_B
            print *, "Expected average number of background objects within theta_max: ", approx_points_in_theta_max_ball
            print "(A, F7.2)", "Range of background x coordinates:", x_range
            print "(A, F7.2)", "Range of background y coordinates:", y_range
        end if

        progress_counter = 0
        last_progress_point = 0

        avg_theta = 0.0
        g_t = 0.0
        g_x = 0.0
        weights = 0.0
        n_pairs = 0
        K = 0.0
        xi_E_marks = 0.0
        xi_B_marks = 0.0
        K_marks = 0.0
        w_marks = 0.0
        marks_weight = 0.0

        if(use_kd_tree) then
            tree_data = reshape([x_B, y_B], shape=[n_B, 2])
            if(verbose) print *, "Building tree"
            call grow_tree(kd_tree, tree_data, leafsize=10)
            if(verbose) print *, "Finished building tree"
        end if

!$OMP parallel default(none) num_threads(num_threads) shared(n_A, n_B, use_kd_tree, verbose, progress_counter, last_progress_point, tree_data, kd_tree, size_kd_tree_results, bins, x_A, y_A, z_A, val_A, weight_A, x_B, y_B, z_B, e1_B, e2_B, weight_B, m_B, spherical_coords, calculate_marks, xi_E_marks, xi_B_marks, K_marks, w_marks, marks_weight) private(i, j, x_foreground, y_foreground, z_foreground, value_foreground, weight_foreground, coord_A, x_background, y_background, z_background, n_nearest_found, kd_tree_results, index, delta_x, delta_y, distance, w, m, bin_index, cos2phi, sin2phi, cos_sin_theta, e1, e2) reduction(+: avg_theta, g_t, g_x, weights, n_pairs, K)

        if(use_kd_tree) allocate(kd_tree_results(size_kd_tree_results))

        !$OMP do schedule(auto) 
        do i=1,n_A
            !$OMP atomic
            progress_counter = progress_counter + 1
            !$OMP atomic
            last_progress_point = last_progress_point + 1
            !$OMP critical
            if(verbose .and. 1.0*last_progress_point/n_A >= 0.1) then
                print "(A, F5.1, A)", "Progress: ", 100.0*progress_counter/n_A, "%."
                last_progress_point = 0
            end if
            !$OMP end critical

            x_foreground = x_A(i)
            y_foreground = y_A(i)
            z_foreground = z_A(i)
            value_foreground = val_A(i)
            weight_foreground = weight_A(i)
            coord_A = [x_foreground, y_foreground]

            if(weight_foreground == 0.0) cycle

            if(use_kd_tree) then
                if(.not. spherical_coords) then
                    call search_within_r(kd_tree, coord_A, bins%x_max, kd_tree_results, n_nearest_found, .false.)
                else if(spherical_coords) then
                    call search_within_r_spherical(kd_tree, coord_A, bins%x_max, kd_tree_results, n_nearest_found, .false.)
                end if
            else
                n_nearest_found = n_B
            end if
            
            do j=1,n_nearest_found
                if(use_kd_tree) then
                    index = kd_tree_results(j)
                else
                    index = j
                end if

                z_background = z_B(index)

                if((z_foreground < 0.0) .or. (z_background < 0.0)) then
                    !Masked object
                    cycle
                end if

                x_background = x_B(index)
                y_background = y_B(index)

                if(spherical_coords) then
                    distance = distance_sphere(x_foreground, y_foreground, x_background, y_background)
                else
                    delta_x = x_foreground - x_background
                    delta_y = y_foreground - y_background
                    distance = sqrt(delta_x**2 + delta_y**2)
                end if

                if(z_foreground > z_background) then
                    cycle
                end if

                bin_index = find_bin_index(distance, bins)
                if(bin_index < 1) then
                    !distance doesn't fall into any of the bins
                    cycle
                end if

                w = weight_B(index)
                m = m_B(index)
                e1 = e1_B(index)
                e2 = e2_B(index)

                if(spherical_coords) then
                    cos_sin_theta = angles_sphere(x_foreground, y_foreground, &
                                                  x_background, y_background, distance)
                    cos2phi = 2.0*cos_sin_theta(1)**2 - 1.0
                    sin2phi = 2.0*cos_sin_theta(1)*cos_sin_theta(2)
                else
                    cos2phi = 2.0*(delta_x/distance)**2 - 1.0
                    sin2phi = 2.0*delta_x*delta_y/distance**2
                end if

                avg_theta(bin_index) = avg_theta(bin_index) + distance*w*weight_foreground
                g_t(bin_index) = g_t(bin_index) - (e1*cos2phi + e2*sin2phi)*w*value_foreground*weight_foreground
                g_x(bin_index) = g_x(bin_index) - (-e1*sin2phi + e2*cos2phi )*w*value_foreground*weight_foreground
                weights(bin_index) = weights(bin_index) + w*weight_foreground
                n_pairs(bin_index) = n_pairs(bin_index) + 1
                K(bin_index) = K(bin_index) + w*(1.0+m)*weight_foreground

                xi_E_marks(bin_index, i) = xi_E_marks(bin_index, i) - (e1*cos2phi + e2*sin2phi)*w*value_foreground*weight_foreground
                xi_B_marks(bin_index, i) = xi_B_marks(bin_index, i) - (-e1*sin2phi + e2*cos2phi )*w*value_foreground*weight_foreground
                w_marks(bin_index, i) = w_marks(bin_index, i) + w*weight_foreground
                K_marks(bin_index, i) = K_marks(bin_index, i) + w*(1.0+m)*weight_foreground
                marks_weight(i) = marks_weight(i) + w*weight_foreground
            end do
        end do
        !$OMP end do nowait
        if(use_kd_tree) deallocate(kd_tree_results)
!$OMP end parallel
        if(use_kd_tree) then
            call burn_tree(kd_tree)
            deallocate(tree_data)
        end if

        avg_theta = avg_theta/weights
        where(K > 0.0)
            g_t = g_t/K
            g_x = g_x/K
        else where
            g_t = 0.0
            g_x = 0.0
        end where

    end subroutine tangential_shear_tpcf

    subroutine scalar_kappa_tpcf(x_A, y_A, val_A, w_A,&
                                 x_B, y_B, kappa_E, kappa_B, w_B, &
                                 bins, &
                                 avg_theta, g_t, g_x, weights, n_pairs, &
                                 spherical_coords, &
                                 num_threads, use_kd_tree, verbose, &
                                 calculate_marks, xi_E_marks, xi_B_marks, K_marks, w_marks, marks_weight)
        !Arguments
        real, dimension(:), intent(in) :: x_A, y_A, val_A, w_A, x_B, y_B, kappa_E, kappa_B, w_B
        type(binning), intent(in) :: bins
        real, dimension(bins%n_bin), intent(out) :: avg_theta, g_t, g_x, weights
        integer(kind=8), dimension(bins%n_bin), intent(out) :: n_pairs
        logical, intent(in) :: spherical_coords, use_kd_tree, verbose, calculate_marks
        integer, intent(in) :: num_threads
        real, dimension(bins%n_bin,size(x_A)), intent(out) :: xi_E_marks, xi_B_marks, K_marks, w_marks
        real, dimension(size(x_A)), intent(out) :: marks_weight
        !Variables
        real, dimension(bins%n_bin) :: K
        real :: distance, delta_x, delta_y, x_foreground, y_foreground, value_foreground, weight_foreground, x_background, y_background, k_E, k_B, w, x_range, y_range, point_density
        real, dimension(2) :: coord_A
        integer :: n_A, n_B, i, j, index, bin_index, progress_counter, last_progress_point, approx_points_in_theta_max_ball, n_nearest_found, size_kd_tree_results

        real, allocatable, dimension(:,:) :: tree_data
        type(kdtree) :: kd_tree
        integer, allocatable, dimension(:) :: kd_tree_results

        n_A = size(x_A)
        n_B = size(x_B)

        x_range = maxval(x_B) - minval(x_B)
        y_range = maxval(y_B) - minval(y_B)
        point_density = n_B/(x_range*y_range)
        approx_points_in_theta_max_ball = pi*bins%x_max**2 * point_density
        size_kd_tree_results = approx_points_in_theta_max_ball*10.0

        if(verbose) then
            if(use_kd_tree) then
                print *,"Using kd tree."
            else
                print *,"Not using kd tree."
            end if

            print *, "Number of foreground objects:", n_A
            print *, "Number of background objects:", n_B
            print *, "Expected average number of background objects within theta_max: ", approx_points_in_theta_max_ball
            print "(A, F7.2)", "Range of background x coordinates:", x_range
            print "(A, F7.2)", "Range of background y coordinates:", y_range
        end if

        progress_counter = 0
        last_progress_point = 0

        avg_theta = 0.0
        g_t = 0.0
        g_x = 0.0
        weights = 0.0
        n_pairs = 0
        K = 0.0
        xi_E_marks = 0.0
        xi_B_marks = 0.0
        K_marks = 0.0
        w_marks = 0.0
        marks_weight = 0.0

        if(use_kd_tree) then
            tree_data = reshape([x_B, y_B], shape=[n_B, 2])
            if(verbose) print *, "Building tree"
            call grow_tree(kd_tree, tree_data, leafsize=10)
            if(verbose) print *, "Finished building tree"
        end if

!$OMP parallel default(none) num_threads(num_threads) shared(n_A, n_B, use_kd_tree, verbose, progress_counter, last_progress_point, tree_data, kd_tree, size_kd_tree_results, bins, x_A, y_A, val_A, w_A, x_B, y_B, kappa_E, kappa_B, w_B, spherical_coords, calculate_marks, xi_E_marks, xi_B_marks, K_marks, w_marks, marks_weight) private(i, j, x_foreground, y_foreground, value_foreground, weight_foreground, coord_A, x_background, y_background, n_nearest_found, kd_tree_results, index, delta_x, delta_y, distance, w, bin_index, k_E, k_B) reduction(+: avg_theta, g_t, g_x, weights, n_pairs, K)

        if(use_kd_tree) allocate(kd_tree_results(size_kd_tree_results))

        !$OMP do schedule(auto) 
        do i=1,n_A
            !$OMP atomic
            progress_counter = progress_counter + 1
            !$OMP atomic
            last_progress_point = last_progress_point + 1
            !$OMP critical
            if(verbose .and. 1.0*last_progress_point/n_A >= 0.1) then
                print "(A, F5.1, A)", "Progress: ", 100.0*progress_counter/n_A, "%."
                last_progress_point = 0
            end if
            !$OMP end critical

            x_foreground = x_A(i)
            y_foreground = y_A(i)
            value_foreground = val_A(i)
            weight_foreground = w_A(i)
            coord_A = [x_foreground, y_foreground]

            if(weight_foreground == 0.0) cycle

            if(use_kd_tree) then
                if(.not. spherical_coords) then
                    call search_within_r(kd_tree, coord_A, bins%x_max, kd_tree_results, n_nearest_found, .false.)
                else if(spherical_coords) then
                    call search_within_r_spherical(kd_tree, coord_A, bins%x_max, kd_tree_results, n_nearest_found, .false.)
                end if
            else
                n_nearest_found = n_B
            end if
            
            do j=1,n_nearest_found
                if(use_kd_tree) then
                    index = kd_tree_results(j)
                else
                    index = j
                end if

                x_background = x_B(index)
                y_background = y_B(index)

                if(spherical_coords) then
                    distance = distance_sphere(x_foreground, y_foreground, x_background, y_background)
                else
                    delta_x = x_foreground - x_background
                    delta_y = y_foreground - y_background
                    distance = sqrt(delta_x**2 + delta_y**2)
                end if

                bin_index = find_bin_index(distance, bins)
                if(bin_index < 1) then
                    !distance doesn't fall into any of the bins
                    cycle
                end if

                w = w_B(index)
                k_E = kappa_E(index)
                k_B = kappa_B(index)

                avg_theta(bin_index) = avg_theta(bin_index) + distance*w*weight_foreground
                g_t(bin_index) = g_t(bin_index) + k_E*w*value_foreground*weight_foreground
                g_x(bin_index) = g_x(bin_index) + k_B*w*value_foreground*weight_foreground
                weights(bin_index) = weights(bin_index) + w*weight_foreground
                n_pairs(bin_index) = n_pairs(bin_index) + 1
                K(bin_index) = K(bin_index) + w*weight_foreground

                xi_E_marks(bin_index, i) = xi_E_marks(bin_index, i) + k_E*w*value_foreground*weight_foreground
                xi_B_marks(bin_index, i) = xi_B_marks(bin_index, i) + k_B*w*value_foreground*weight_foreground
                w_marks(bin_index, i) = w_marks(bin_index, i) + w*weight_foreground
                K_marks(bin_index, i) = K_marks(bin_index, i) + w*weight_foreground
                marks_weight(i) = marks_weight(i) + w*weight_foreground
            end do
        end do
        !$OMP end do nowait
        if(use_kd_tree) deallocate(kd_tree_results)
!$OMP end parallel
        if(use_kd_tree) then
            call burn_tree(kd_tree)
            deallocate(tree_data)
        end if

        avg_theta = avg_theta/weights
        where(K > 0.0)
            g_t = g_t/K
            g_x = g_x/K
        else where
            g_t = 0.0
            g_x = 0.0
        end where

    end subroutine scalar_kappa_tpcf

    subroutine simple_bootstrap_foreground(field, &
                                           bins, &
                                           n_bootstrap_samples, &
                                           results, &
                                           spherical_coords, &
                                           n_threads, use_kd_tree, verbose, &
                                           use_marks, xi_E_marks, xi_B_marks, K_marks, marks_weight)
        !Arguments
        type(lensing_survey_field), intent(in) :: field
        type(binning), intent(in) :: bins
        integer, intent(in) :: n_bootstrap_samples
        type(shear_2pcf_results), dimension(n_bootstrap_samples), intent(out) :: results
        logical, intent(in) :: spherical_coords
        integer, intent(in) :: n_threads
        logical, intent(in) :: use_kd_tree, verbose, use_marks
        real, optional, dimension(:,:), intent(in) :: xi_E_marks, xi_B_marks, K_marks
        real, optional, dimension(:), intent(in) :: marks_weight
        !Variables
        integer :: i, j, n_A, n_unique, index, object_weight
        integer, allocatable, dimension(:) :: indicies, index_counter
        real, allocatable, dimension(:) :: K
        type(lensing_survey_field) :: field_bs

        if(use_marks .and. .not. (present(xi_E_marks) .and. present(xi_B_marks) .and. present(K_marks) .and. present(marks_weight))) then
            print *, "Did not supply data needed for marked point bootstrap. Exiting"
            return
        end if

        n_A = field%lenses%n

        allocate(field_bs%shapes%x, mold=field%shapes%x)
        allocate(field_bs%shapes%y, mold=field%shapes%y)
        allocate(field_bs%shapes%z, mold=field%shapes%z)
        allocate(field_bs%shapes%e1, mold=field%shapes%e1)
        allocate(field_bs%shapes%e2, mold=field%shapes%e2)
        allocate(field_bs%shapes%w, mold=field%shapes%w)
        allocate(field_bs%shapes%m, mold=field%shapes%m)

        allocate(indicies(n_A), index_counter(n_A), K(bins%n_bin))
        do i=1, n_bootstrap_samples
            call random_number_int(array=indicies, max=n_A, min=1)
            index_counter = 0.0
            do j=1, n_A
                index = indicies(j)
                index_counter(index) = index_counter(index) + 1
            end do

            if(use_marks) then
                call allocate_shear_2pcf_results(results(i), bins%n_bin, 0)
                results(i)%g_t = 0.0
                results(i)%g_x = 0.0
                results(i)%n_pairs = 0
                results(i)%weights = 0.0
                K = 0.0
                do j=1, n_A
                    object_weight = index_counter(j)
                    if(object_weight > 0) then
                        results(i)%g_t = results(i)%g_t + xi_E_marks(:,j)*object_weight
                        results(i)%g_x = results(i)%g_x + xi_B_marks(:,j)*object_weight
                        results(i)%weights = results(i)%weights + marks_weight(j)*object_weight
                        K(:) = K(:) + K_marks(:,j)*object_weight
                    end if
                end do
                where(K > 0.0)
                    results(i)%g_t = results(i)%g_t/K
                    results(i)%g_x = results(i)%g_x/K
                else where
                    results(i)%g_t = 0.0
                    results(i)%g_x = 0.0
                end where
                cycle
            end if

            field_bs%lenses%n = count(index_counter > 0)
            field_bs%lenses%x = pack(field%lenses%x, mask=index_counter > 0)
            field_bs%lenses%y = pack(field%lenses%y, mask=index_counter > 0)
            field_bs%lenses%z = pack(field%lenses%z, mask=index_counter > 0)
            field_bs%lenses%val = pack(field%lenses%val, mask=index_counter > 0)
            field_bs%lenses%w = pack(field%lenses%w*index_counter, mask=index_counter > 0)

            if(verbose) print *, "Running tpcf with randomized foreground catalog nr", i
            call shear_2pcf(field_bs, &
                            bins, &
                            results(i), &
                            spherical_coords, n_threads, use_kd_tree, .false.)
        end do
        deallocate(indicies, index_counter, K)
    end subroutine simple_bootstrap_foreground

!    subroutine simple_bootstrap_background(x_A, y_A, z_A, val_A, weight_A, &
!                                           x_B, y_B, z_B, e1_B, e2_B, weight_B, m_B, &
!                                           bins, &
!                                           n_bootstrap_samples, &
!                                           results, &
!                                           spherical_coords, &
!                                           n_threads, use_kd_tree, verbose)
!        !Arguments
!        real, intent(in) :: x_A(:), y_A(:), z_A(:), val_A(:), weight_A(:), x_B(:), y_B(:), z_B(:), e1_B(:), e2_B(:), weight_B(:), m_B(:)
!        type(binning), intent(in) :: bins
!        integer, intent(in) :: n_bootstrap_samples
!        type(shear_2pcf_results), dimension(n_bootstrap_samples), intent(out) :: results
!        logical, intent(in) :: spherical_coords
!        integer, intent(in) :: n_threads
!        logical, intent(in) :: use_kd_tree, verbose
!        !Variables
!        integer :: i, j, k, n_B, n_unique, index
!        integer, allocatable :: indicies(:), index_counter(:)
!        real, allocatable :: x_bs(:), y_bs(:), z_bs(:), e1_bs(:), e2_bs(:), weight_bs(:), m_bs(:)
!
!        n_B = size(x_B)
!
!        allocate(indicies(n_B), index_counter(n_B))
!        do i=1, n_bootstrap_samples
!            call random_number_int(array=indicies, max=n_B, min=1)
!            index_counter = 0
!            do j=1, n_B
!                index = indicies(j)
!                index_counter(index) = index_counter(index) + 1
!            end do
!            n_unique = count(index_counter > 0)
!            allocate(x_bs(n_unique), y_bs(n_unique), z_bs(n_unique), e1_bs(n_unique), e2_bs(n_unique), weight_bs(n_unique), m_bs(n_unique))
!            k = 1
!            do j=1, n_B
!                if(index_counter(j) > 0) then
!                    x_bs(k) = x_B(j)
!                    y_bs(k) = y_B(j)
!                    z_bs(k) = z_B(j)
!                    e1_bs(k) = e1_B(j)
!                    e2_bs(k) = e2_B(j)
!                    m_bs(k) = m_B(j)
!                    weight_bs(k) = weight_B(j)*index_counter(j)
!                    k = k + 1
!                    if(k > n_unique) then
!                        exit
!                    end if
!                else
!                    cycle
!                end if
!            end do
!            if(verbose) print *, "Running tpcf with randomized background catalog nr", i
!            call shear_2pcf(x_A, y_A, z_A, val_A, weight_A, &
!                            x_bs, y_bs, z_bs, e1_bs, e2_bs, weight_bs, m_bs, &
!                            bins, &
!                            results(i), &
!                            spherical_coords, n_threads, use_kd_tree, .false.)
!            deallocate(x_bs, y_bs, z_bs, e1_bs, e2_bs, weight_bs, m_bs)
!        end do
!        deallocate(indicies, index_counter)
!    end subroutine simple_bootstrap_background
!
!    subroutine fixed_block_bootstrap(x_A, y_A, z_A, val_A, weight_A, &
!                                     x_B, y_B, z_B, e1_B, e2_B, weight_B, m_B, &
!                                     bins, &
!                                     n_x_blocks, n_y_blocks, n_bootstrap_samples, supersampling_factor, &
!                                     results, &
!                                     spherical_coords, &
!                                     n_threads, use_kd_tree, verbose, &
!                                     use_marks, xi_E_marks, xi_B_marks, K_marks, marks_weight)
!        !Arguments
!        real, intent(in) :: x_A(:), y_A(:), z_A(:), val_A(:), weight_A(:), x_B(:), y_B(:), z_B(:), e1_B(:), e2_B(:), weight_B(:), m_B(:)
!        type(binning), intent(in) :: bins
!        integer, intent(in) :: n_bootstrap_samples, n_x_blocks, n_y_blocks
!        real, intent(in) :: supersampling_factor
!        type(shear_2pcf_results), dimension(n_bootstrap_samples), intent(out) :: results
!        logical, intent(in) :: spherical_coords
!        integer, intent(in) :: n_threads
!        logical, intent(in) :: use_kd_tree, verbose, use_marks
!        real, optional, dimension(bins%n_bin,size(x_A)), intent(in) :: xi_E_marks, xi_B_marks, K_marks
!        real, optional, dimension(size(x_A)), intent(in) :: marks_weight
!        !Variables
!        integer :: i, j, n_A, n_B, n_foreground_bs, n_background_bs, x_block_idx, y_block_idx, n_blocks, n_block_samples
!        type(binning) :: x_block_bins, y_block_bins
!        real :: x_min, x_max, y_min, y_max
!        integer, allocatable :: foreground_block_indices(:), background_block_indices(:), block_indicies(:), index_counter(:)
!        integer :: foreground_block_n(n_x_blocks*n_y_blocks), background_block_n(n_x_blocks*n_y_blocks)
!        real :: block_weights(n_x_blocks*n_y_blocks), block_weight
!        real, allocatable :: x_A_bs(:), y_A_bs(:), z_A_bs(:), val_bs(:), w_A_bs(:), K(:)
!        real, allocatable :: x_B_bs(:), y_B_bs(:), z_B_bs(:), e1_bs(:), e2_bs(:), w_B_bs(:), m_bs(:)
!
!        if(use_marks .and. .not. (present(xi_E_marks) .and. present(xi_B_marks) .and. present(K_marks) .and. present(marks_weight))) then
!            print *, "Did not supply data needed for marked point bootstrap. Exiting"
!            return
!        end if
!
!        n_blocks = n_x_blocks*n_y_blocks
!        if(n_blocks <= 1) then
!            print *, "At least 2 blocks are required. Exiting."
!            return
!        end if
!
!        n_block_samples = int(n_blocks*supersampling_factor)
!        n_A = size(x_A)
!        n_B = size(x_B)
!        allocate(foreground_block_indices(n_A), background_block_indices(n_B), block_indicies(n_block_samples), K(bins%n_bin), index_counter(n_A))
!
!        x_min = minval(x_A)
!        x_max = maxval(x_A)
!        y_min = minval(y_A)
!        y_max = maxval(y_A)
!        call create_bins(x_min=x_min, x_max=x_max, n_bin=n_x_blocks, spacing=lin_binning, bins=x_block_bins)
!        call create_bins(x_min=y_min, x_max=y_max, n_bin=n_y_blocks, spacing=lin_binning, bins=y_block_bins)
!
!        !Assign index of block that contains object to object
!        do i=1,n_A
!            x_block_idx = find_bin_index(x_A(i), x_block_bins)
!            y_block_idx = find_bin_index(y_A(i), y_block_bins)
!            foreground_block_indices(i) = (y_block_idx-1)*n_y_blocks + x_block_idx
!        end do
!        do i=1,n_B
!            if(x_B(i) < x_min) then
!                x_block_idx = 1
!            else if(x_B(i) >= x_max) then
!                x_block_idx = n_x_blocks
!            else
!                x_block_idx = find_bin_index(x_B(i), x_block_bins)
!            endif
!            if(y_B(i) < y_min) then
!                y_block_idx = 1
!            else if(y_B(i) >= y_max) then
!                y_block_idx = n_y_blocks
!            else
!                y_block_idx = find_bin_index(y_B(i), y_block_bins)
!            endif
!            background_block_indices(i) = (y_block_idx-1)*n_y_blocks + x_block_idx
!        end do
!
!        !Count how many objects are within each block
!        do i=1,n_blocks
!            foreground_block_n(i) = count(foreground_block_indices == i)
!            background_block_n(i) = count(background_block_indices == i)
!        end do
!
!        do i=1,n_bootstrap_samples
!            call random_number_int(array=block_indicies, max=n_blocks, min=1)
!
!            !Count number of unique objects within the drawn blocks
!            n_foreground_bs = 0
!            n_background_bs = 0
!            do j=1,n_blocks
!                block_weights(j) = count(block_indicies == j)
!                if(block_weights(j) > 0) then
!                    n_foreground_bs = n_foreground_bs + foreground_block_n(j)
!                    n_background_bs = n_background_bs + background_block_n(j)
!                end if
!            end do
!
!            if(use_marks) then
!                call allocate_shear_2pcf_results(results(i), bins%n_bin, 0)
!                results(i)%g_t = 0.0
!                results(i)%g_x = 0.0
!                results(i)%weights = 0.0
!                results(i)%n_pairs = 0
!                K = 0.0
!                do j=1, n_A
!                    block_weight = block_weights(foreground_block_indices(j))
!                    if(block_weight > 0) then
!                        results(i)%g_t = results(i)%g_t + xi_E_marks(:,j)*block_weight
!                        results(i)%g_x = results(i)%g_x + xi_B_marks(:,j)*block_weight
!                        results(i)%weights = results(i)%weights + marks_weight(j)*block_weight
!                        K(:) = K(:) + K_marks(:,j)*block_weight
!                    end if
!                end do
!                results(i)%g_t = results(i)%g_t/K
!                results(i)%g_x = results(i)%g_x/K
!                cycle
!            end if
!
!            forall(j=1:n_A) index_counter(j) = block_weights(foreground_block_indices(j))
!            x_A_bs = pack(x_A, mask=index_counter > 0)
!            y_A_bs = pack(y_A, mask=index_counter > 0)
!            z_A_bs = pack(z_A, mask=index_counter > 0)
!            val_bs = pack(val_A, mask=index_counter > 0)
!            w_A_bs = pack(weight_A*index_counter, mask=index_counter > 0)
!
!            forall(j=1:n_B) index_counter(j) = block_weights(background_block_indices(j))
!            x_B_bs = pack(x_B, mask=index_counter > 0)
!            y_B_bs = pack(y_B, mask=index_counter > 0)
!            z_B_bs = pack(z_B, mask=index_counter > 0)
!            e1_bs = pack(e1_B, mask=index_counter > 0)
!            e2_bs = pack(e2_B, mask=index_counter > 0)
!            w_B_bs = pack(weight_B*index_counter, mask=index_counter > 0)
!            m_bs = pack(m_B, mask=index_counter > 0)
!
!            if(verbose) print *, "Running tpcf with randomized blocks, run", i
!            call shear_2pcf(x_A_bs, y_A_bs, z_A_bs, val_bs, w_A_bs, &
!                            x_B_bs, y_B_bs, z_B_bs, e1_bs, e2_bs, w_B_bs, m_bs, &
!                            bins, &
!                            results(i), &
!                            spherical_coords, n_threads, use_kd_tree, .false.)
!
!            deallocate(x_A_bs, y_A_bs, z_A_bs, val_bs, w_A_bs)
!            deallocate(x_B_bs, y_B_bs, z_B_bs, e1_bs, e2_bs, w_B_bs, m_bs)
!        end do
!
!        deallocate(foreground_block_indices, background_block_indices, block_indicies, K, index_counter)
!    end subroutine fixed_block_bootstrap

    subroutine pair_weighted_block_bootstrap(x_A, y_A, z_A, val_A, weight_A, &
                                             x_B, y_B, z_B, e1_B, e2_B, weight_B, m_B, &
                                             bins, &
                                             n_blocks, n_bootstrap_samples, supersampling_factor, recalculate_blocks, &
                                             results, &
                                             spherical_coords, &
                                             n_threads, use_kd_tree, verbose, &
                                             xi_E_marks, xi_B_marks, K_marks, marks_weight)
        !Arguments
        real, dimension(:), intent(in) :: x_A, y_A, z_A, val_A, weight_A, x_B, y_B, z_B, e1_B, e2_B, weight_B, m_B
        type(binning), intent(in) :: bins
        integer, intent(in) :: n_bootstrap_samples, n_blocks
        real, intent(in) :: supersampling_factor
        logical, intent(in) :: recalculate_blocks
        type(shear_2pcf_results), dimension(n_bootstrap_samples), intent(inout) :: results
        logical, intent(in) :: spherical_coords
        integer, intent(in) :: n_threads
        logical, intent(in) :: use_kd_tree, verbose
        real, dimension(bins%n_bin,size(x_A)), intent(in) :: xi_E_marks, xi_B_marks, K_marks
        real, dimension(size(x_A)), intent(in) :: marks_weight
        !Variables
        integer :: i, j, n_A, n_B, n_block_samples, block_count
        integer, allocatable, dimension(:) :: block_indicies, drawn_blocks, drawn_block_counts
        real, allocatable, dimension(:) :: pair_weights, K(:)
        real, allocatable, dimension(:,:) :: block_centroids
        procedure(distance_function_interface), pointer :: distance_func_ptr

        if(spherical_coords) then
            distance_func_ptr => distance_sphere
        else
            distance_func_ptr => distance_euclidean
        end if

        if(n_blocks <= 1) then
            print *, "At least 2 blocks are required. Exiting."
            return
        end if

        n_A = size(x_A)
        n_B = size(x_B)
        if(n_A /= count(z_A >= 0.0)) then
            print *, "Found masked regions in data. Exiting."
            return
        end if

        allocate(pair_weights(n_A), block_indicies(n_A), block_centroids(n_A, 2), K(bins%n_bin))
        pair_weights = 0.01 + marks_weight/maxval(marks_weight)

        if(.not. recalculate_blocks) then
            call partition_data(data_points=reshape([x_A, y_A], shape=[n_A, 2]), weights=pair_weights, &
                                                    k=n_blocks, cluster_indicies=block_indicies, cluster_centroids=block_centroids, &
                                                    distance_func=distance_func_ptr, &
                                                    equal_size=.true., compute_initial_centroids=kmeans_compute_random_initial_centroids, &
                                                    weight_equalizer_scaling=3.0, n_threads=n_threads, iterations=1000, reassignment_tol=0.01, weight_tol=0.05)
        end if

        n_block_samples = int(n_blocks*supersampling_factor)
        allocate(drawn_blocks(n_block_samples), drawn_block_counts(n_blocks))

        do i=1,n_bootstrap_samples
            if(recalculate_blocks) then
                print *, "Computing blocks, n = ", i
                call partition_data(data_points=reshape([x_A, y_A], shape=[n_A, 2]), weights=pair_weights, &
                                                        k=n_blocks, cluster_indicies=block_indicies, cluster_centroids=block_centroids, &
                                                        distance_func=distance_func_ptr, &
                                                        equal_size=.true., compute_initial_centroids=kmeans_compute_random_initial_centroids, &
                                                        weight_equalizer_scaling=3.0, n_threads=n_threads, iterations=1000, reassignment_tol=0.01, weight_tol=0.05)
            end if
            call random_number_int(array=drawn_blocks, max=n_blocks, min=1)

            !Count number of unique objects within the drawn blocks
            do j=1,n_blocks
                drawn_block_counts(j) = count(drawn_blocks == j)
            end do

            call allocate_shear_2pcf_results(results(i), bins%n_bin, 0)
            results(i)%g_t = 0.0
            results(i)%g_x = 0.0
            results(i)%n_pairs = 0
            results(i)%weights = 0.0
            K = 0.0
            do j=1,n_A
                block_count = drawn_block_counts(block_indicies(j))
                if(block_count > 0) then
                    results(i)%g_t = results(i)%g_t + xi_E_marks(:,j)*block_count
                    results(i)%g_x = results(i)%g_x + xi_B_marks(:,j)*block_count
                    results(i)%weights = results(i)%weights + marks_weight(j)*block_count
                    K(:) = K(:) + K_marks(:,j)*block_count
                end if
            end do
            results(i)%g_t = results(i)%g_t/K
            results(i)%g_x = results(i)%g_x/K
        end do

        deallocate(pair_weights, block_indicies, block_centroids, drawn_blocks, drawn_block_counts, K)
    end subroutine pair_weighted_block_bootstrap

!    subroutine randomize_shear(x_A, y_A, z_A, val_A, weight_A, &
!                               x_B, y_B, z_B, e1_B, e2_B, weight_B, m_B, &
!                               bins, &
!                               n_bootstrap_samples, &
!                               results, &
!                               spherical_coords, &
!                               n_threads, use_kd_tree, verbose)
!        !Arguments
!        real, intent(in) :: x_A(:), y_A(:), z_A(:), val_A(:), weight_A(:), x_B(:), y_B(:), z_B(:), e1_B(:), e2_B(:), weight_B(:), m_B(:)
!        type(binning), intent(in) :: bins
!        integer, intent(in) :: n_bootstrap_samples
!        type(shear_2pcf_results), dimension(n_bootstrap_samples), intent(out) :: results
!        logical, intent(in) :: spherical_coords
!        integer, intent(in) :: n_threads
!        logical, intent(in) :: use_kd_tree, verbose
!        !Variables
!        integer :: i, n_B
!        real, allocatable :: alpha(:), e1_bs(:), e2_bs(:)
!
!        n_B = size(x_B)
!
!        allocate(alpha(n_B),e1_bs(n_B), e2_bs(n_B))
!        do i=1, n_bootstrap_samples
!            call random_number(alpha)
!            alpha = 2.0*pi*alpha
!            e1_bs = cos(2.0*alpha)*e1_B - sin(2.0*alpha)*e2_B
!            e2_bs = sin(2.0*alpha)*e1_B + cos(2.0*alpha)*e2_B
!            if(verbose) print *, "Running tpcf with randomized shear catalog nr", i
!            call shear_2pcf(x_A, y_A, z_A, val_A, weight_A, &
!                            x_B, y_B, z_B, e1_bs, e2_bs, weight_B, m_B, &
!                            bins, &
!                            results(i), &
!                            spherical_coords, n_threads, use_kd_tree, .false.)
!        end do
!        deallocate(alpha, e1_bs, e2_bs)
!    end subroutine randomize_shear

    subroutine survey_marked_point_block_bootstrap( survey, bins, &
                                                    n_bootstrap_samples, n_blocks, n_x_blocks, n_y_blocks, &
                                                    clustering_scheme, supersampling_factor, &
                                                    spherical_coords, n_threads, &
                                                    results, block_output_prefix)
        !Arguments
        type(marks_field_collection), intent(inout) :: survey
        type(binning), intent(in) :: bins
        integer, intent(in) :: n_bootstrap_samples, n_blocks, n_x_blocks, n_y_blocks, clustering_scheme, n_threads
        real, intent(in) :: supersampling_factor
        logical, intent(in) :: spherical_coords
        type(shear_2pcf_results), dimension(n_bootstrap_samples), intent(inout) :: results
        character(len=256), optional, intent(in) :: block_output_prefix
        !Variables
        integer :: i, j, n_objects, x_block_idx, y_block_idx, block_index_offset, n_blocks_total, block_index, n_block_samples, n_blocks_field
        type(binning) :: x_block_bins, y_block_bins
        real :: x_min, x_max, y_min, y_max, weight_per_block
        real, allocatable, dimension(:) :: block_weight, normalized_weights, field_weight, K, g_t, g_x
        real, allocatable, dimension(:,:) :: block_g_t, block_g_x, block_K, block_centroids
        integer, allocatable, dimension(:) :: drawn_block_indicies, drawn_block_multiplicity
        procedure(distance_function_interface), pointer :: distance_func_ptr

        logical :: output_blocks
        character(len=256) :: block_output_filename
        integer :: file_unit, iostat

        output_blocks = .false.
        if(present(block_output_prefix) .and. block_output_prefix /= "") output_blocks = .true.

        if(.not. (n_blocks > 1 .or. n_x_blocks*n_y_blocks > 1)) then
            print *, "At least 2 blocks are required. Exiting."
            return
        end if

        !Average weight per block (use weight from marks)
        allocate(field_weight(survey%n_field), K(bins%n_bin))
        forall(i=1:survey%n_field) field_weight(i) = sum(survey%fields(i)%w)

        !Assign blocks
        if(clustering_scheme == rectangular_blocks) then
            print *, "Forming rectangular blocks."
            n_blocks_total = survey%n_field * n_x_blocks*n_y_blocks
            weight_per_block = sum(field_weight)/n_blocks_total

            allocate(block_weight(n_blocks_total), block_g_t(bins%n_bin, n_blocks_total), block_g_x(bins%n_bin, n_blocks_total), block_K(bins%n_bin, n_blocks_total))
            block_index_offset = 0
            block_g_t = 0.0
            block_g_x = 0.0
            block_weight = 0.0
            block_K = 0.0

            do i=1, survey%n_field
                n_objects = survey%fields(i)%n

                if(.not. allocated(survey%fields(i)%block_indicies)) allocate(survey%fields(i)%block_indicies(n_objects))

                x_min = minval(survey%fields(i)%x)
                x_max = maxval(survey%fields(i)%x)
                y_min = minval(survey%fields(i)%y)
                y_max = maxval(survey%fields(i)%y)
                call create_bins(x_min=x_min, x_max=x_max, n_bin=n_x_blocks, spacing=lin_binning, bins=x_block_bins)
                call create_bins(x_min=y_min, x_max=y_max, n_bin=n_y_blocks, spacing=lin_binning, bins=y_block_bins)

                !Assign index of block that contains object
                do j=1,n_objects
                    x_block_idx = find_bin_index(survey%fields(i)%x(j), x_block_bins)
                    y_block_idx = find_bin_index(survey%fields(i)%y(j), y_block_bins)
                    block_index = (y_block_idx-1)*n_y_blocks + x_block_idx + block_index_offset
                    survey%fields(i)%block_indicies(j) = block_index

                    !Calcualte weight, g_t, g_x contribution to block
                    block_weight(block_index) = block_weight(block_index) + survey%fields(i)%w(j)
                    block_g_t(:,block_index) = block_g_t(:,block_index) + survey%fields(i)%xi_E_marks(:,j)
                    block_g_x(:,block_index) = block_g_x(:,block_index) + survey%fields(i)%xi_B_marks(:,j)
                    block_K(:,block_index) = block_K(:,block_index) + survey%fields(i)%K_marks(:,j)
                end do
                block_index_offset = block_index_offset + n_x_blocks*n_y_blocks
            end do
            if(count(block_weight == 0.0) > 0) then
                !Remove empty blocks
                print *, "Empty blocks: ", count(block_weight == 0.0)
                n_blocks_total = count(block_weight > 0)
                block_g_t = reshape(pack(block_g_t, mask=spread(block_weight > 0, dim=1, ncopies=bins%n_bin)), shape=[bins%n_bin, n_blocks_total])
                block_g_x = reshape(pack(block_g_x, mask=spread(block_weight > 0, dim=1, ncopies=bins%n_bin)), shape=[bins%n_bin, n_blocks_total])
                block_K = reshape(pack(block_K, mask=spread(block_weight > 0, dim=1, ncopies=bins%n_bin)), shape=[bins%n_bin, n_blocks_total])
                block_weight = pack(block_weight, mask=block_weight > 0)
            end if
            if(count(block_K == 0.0) > 0) then
                print *, "Empty bins: ", count(block_K == 0.0)
                where(block_K == 0)
                    block_g_t = 0.0
                    block_g_x = 0.0
                end where
            end if
        else if(clustering_scheme == equal_weight_blocks .or. clustering_scheme == voronoi_blocks) then
            !=====Equal weight blocks======
            print *, "Forming equal weight blocks."

            if(spherical_coords) then
                distance_func_ptr => distance_sphere
            else
                distance_func_ptr => distance_euclidean
            end if

            weight_per_block = sum(field_weight)/n_blocks
            n_blocks_total = 0
            do i=1,survey%n_field
                n_blocks_total = n_blocks_total + sum(survey%fields(i)%w)/weight_per_block + 1
            end do
            print "(A, I5, A)", "Forming ", n_blocks_total, " blocks."

            allocate(block_weight(n_blocks_total), block_g_t(bins%n_bin, n_blocks_total), block_g_x(bins%n_bin, n_blocks_total), block_K(bins%n_bin, n_blocks_total))
            block_index_offset = 0
            block_g_t = 0.0
            block_g_x = 0.0
            block_weight = 0.0
            block_K = 0.0

            do i=1, survey%n_field
                n_blocks_field = field_weight(i)/weight_per_block + 1
                print "(A, I4, A, I4, A)", "Partitioning field nr. ", i, " into ", n_blocks_field, " blocks."

                n_objects = survey%fields(i)%n
                if(.not. allocated(survey%fields(i)%block_indicies)) allocate(survey%fields(i)%block_indicies(n_objects))

                allocate(normalized_weights(n_objects), block_centroids(n_blocks_field, 2))
                normalized_weights = 0.01 + survey%fields(i)%w/maxval(survey%fields(i)%w)

                if(clustering_scheme == equal_weight_blocks) then
                    call partition_data(data_points=reshape([survey%fields(i)%x, survey%fields(i)%y], shape=[n_objects, 2]), &
                                        weights=normalized_weights, &
                                        k=n_blocks_field, &
                                        cluster_indicies=survey%fields(i)%block_indicies, cluster_centroids=block_centroids, &
                                        distance_func=distance_func_ptr, &
                                        equal_size=.true., compute_initial_centroids=kmeans_compute_fixed_initial_centroids, &
                                        weight_equalizer_scaling=5.0, n_threads=n_threads, iterations=1000, &
                                        reassignment_tol=0.01, weight_tol=0.05)
                else
                    !Voronoi blocks
                    call partition_data(data_points=reshape([survey%fields(i)%x, survey%fields(i)%y], shape=[n_objects, 2]), &
                                        weights=normalized_weights, &
                                        k=n_blocks_field, &
                                        cluster_indicies=survey%fields(i)%block_indicies, cluster_centroids=block_centroids, &
                                        distance_func=distance_func_ptr, &
                                        equal_size=.false., compute_initial_centroids=kmeans_compute_fixed_initial_centroids, &
                                        weight_equalizer_scaling=0.0, n_threads=n_threads, iterations=1000, &
                                        reassignment_tol=0.0, weight_tol=0.0)
                end if

                survey%fields(i)%block_indicies = survey%fields(i)%block_indicies + block_index_offset
                block_index_offset = block_index_offset + n_blocks_field

                do j=1,n_objects
                    block_index = survey%fields(i)%block_indicies(j)

                    block_g_t(:,block_index) = block_g_t(:,block_index) + survey%fields(i)%xi_E_marks(:,j)
                    block_g_x(:,block_index) = block_g_x(:,block_index) + survey%fields(i)%xi_B_marks(:,j)
                    block_weight(block_index) = block_weight(block_index) + survey%fields(i)%w(j)
                    block_K(:,block_index) = block_K(:,block_index) + survey%fields(i)%K_marks(:,j)
                end do

                deallocate(normalized_weights, block_centroids)
            end do
        end if

        print "(A, F5.2)", "Maximum relative block weight devation across fields: ", (maxval(block_weight)-minval(block_weight))/weight_per_block
        print "(A, ES9.2)", "Relative block weight std: ", sqrt(variance(block_weight))/weight_per_block

        if(output_blocks) then
            allocate(g_t(bins%n_bin), g_x(bins%n_bin))
            print *, "Writing block g_t and g_x to files in ", trim(block_output_prefix)
            do i=1,n_blocks_total
                where(block_K(:,i) > 0.0)
                    g_t = block_g_t(:,i)/block_K(:,i)
                    g_x = block_g_x(:,i)/block_K(:,i)
                else where
                    g_t = 0.0
                    g_x = 0.0
                end where
                write(block_output_filename, fmt="(A,I0,A)") trim(block_output_prefix) // "block_", i, ".dat"
                open(newunit=file_unit, file=block_output_filename, iostat=iostat, status="replace")
                call assert_iostat(iostat, block_output_filename)
                write(file_unit, fmt="(es, es, es, I)") ((g_t(j), g_x(j), block_weight(i), 0), j=1,bins%n_bin)
                close(file_unit)
            end do
        end if

        print *, "Performing bootstrap resamplings."
        n_block_samples = int(n_blocks_total*supersampling_factor)
        allocate(drawn_block_indicies(n_block_samples), drawn_block_multiplicity(n_blocks_total))
        do i=1,n_bootstrap_samples
            call random_number_int(array=drawn_block_indicies, max=n_blocks_total, min=1)

            call allocate_shear_2pcf_results(results(i), bins%n_bin, 0)
            results(i)%g_t = 0.0
            results(i)%g_x = 0.0
            results(i)%n_pairs = 0
            results(i)%weights = 0.0
            K(:) = 0.0
            do j=1,n_blocks_total
                drawn_block_multiplicity(j) = count(drawn_block_indicies == j)
                results(i)%g_t(:) = results(i)%g_t(:) + block_g_t(:,j)*drawn_block_multiplicity(j)
                results(i)%g_x(:) = results(i)%g_x(:) + block_g_x(:,j)*drawn_block_multiplicity(j)
                results(i)%weights(:) = results(i)%weights(:) + block_weight(j)*drawn_block_multiplicity(j)
                K(:) = K(:) + block_K(:,j)*drawn_block_multiplicity(j)
            end do
            where(K > 0.0)
                results(i)%g_t = results(i)%g_t/K
                results(i)%g_x = results(i)%g_x/K
            else where
                results(i)%g_t = 0.0
                results(i)%g_x = 0.0
            end where
        end do
    end subroutine survey_marked_point_block_bootstrap

    subroutine survey_marked_point_bootstrap(survey, bins, &
                                             n_bootstrap_samples, supersampling_factor, &
                                             spherical_coords, n_threads, &
                                             results)
        !Arguments
        type(marks_field_collection), intent(inout) :: survey
        type(binning), intent(in) :: bins
        integer, intent(in) :: n_bootstrap_samples, n_threads
        real, intent(in) :: supersampling_factor
        logical, intent(in) :: spherical_coords
        type(shear_2pcf_results), dimension(n_bootstrap_samples), intent(inout) :: results
        !Variables
        integer :: i, j, k_idx, n_objects, n_samples, field_index, object_index
        real, allocatable, dimension(:) :: K
        integer, allocatable, dimension(:) :: drawn_indicies, drawn_multiplicity
        integer, allocatable, dimension(:,:) :: field_index_lookup
        procedure(distance_function_interface), pointer :: distance_func_ptr

        allocate(K(bins%n_bin))

        n_objects = survey%n_marks
        allocate(field_index_lookup(2, n_objects))

        j = 1
        k_idx = 0
        do i=1,n_objects
            if(i > k_idx + survey%fields(j)%n) then
                k_idx = k_idx + survey%fields(j)%n
                j = j + 1
            end if

            field_index_lookup(1,i) = j
            field_index_lookup(2,i) = i-k_idx
        end do

        print *, "Performing bootstrap resamplings."
        n_samples = int(n_objects*supersampling_factor)
        allocate(drawn_indicies(n_samples), drawn_multiplicity(n_objects))
!$OMP PARALLEL DO default(none) num_threads(n_threads) shared(n_bootstrap_samples, n_objects, n_samples, results, bins, field_index_lookup, survey) private(i, j, drawn_indicies, drawn_multiplicity, K, field_index, object_index)
        do i=1,n_bootstrap_samples
            call random_number_int(array=drawn_indicies, max=n_objects, min=1)

            call allocate_shear_2pcf_results(results(i), bins%n_bin, 0)
            results(i)%g_t = 0.0
            results(i)%g_x = 0.0
            results(i)%n_pairs = 0
            results(i)%weights = 0.0
            K(:) = 0.0
            drawn_multiplicity = 0
            do j=1,n_samples
                object_index = drawn_indicies(j)
                drawn_multiplicity(object_index) = drawn_multiplicity(object_index) + 1
            end do

            do j=1,n_objects
                field_index = field_index_lookup(1,j)
                object_index = field_index_lookup(2,j)
                results(i)%g_t(:) = results(i)%g_t(:) + survey%fields(field_index)%xi_E_marks(:,object_index)*drawn_multiplicity(j)
                results(i)%g_x(:) = results(i)%g_x(:) + survey%fields(field_index)%xi_B_marks(:,object_index)*drawn_multiplicity(j)
                results(i)%weights(:) = results(i)%weights(:) + survey%fields(field_index)%w(object_index)*drawn_multiplicity(j)
                K(:) = K(:) + survey%fields(field_index)%K_marks(:,object_index)*drawn_multiplicity(j)
            end do
            where(K > 0.0)
                results(i)%g_t = results(i)%g_t/K
                results(i)%g_x = results(i)%g_x/K
            else where
                results(i)%g_t = 0.0
                results(i)%g_x = 0.0
            end where
        end do
!$OMP END PARALLEL DO
    end subroutine survey_marked_point_bootstrap

end module tpcf