module kmeans_module
    use utils
    implicit none

    public partition_data, kmeans_compute_random_initial_centroids, kmeans_compute_fixed_initial_centroids

    private
    integer, parameter :: kmeans_compute_random_initial_centroids = 1, kmeans_compute_fixed_initial_centroids = 2

    contains

    subroutine partition_data(data_points, weights, k, cluster_indicies, cluster_centroids, distance_func, equal_size, compute_initial_centroids, n_threads, weight_equalizer_scaling, iterations, reassignment_tol, weight_tol)
        !Arguments
        real, dimension(:,:), intent(in) :: data_points
        real, dimension(:), intent(in) :: weights
        integer, intent(in) :: k
        integer, dimension(:), intent(out) :: cluster_indicies
        real, dimension(k,2), intent(inout) :: cluster_centroids
        procedure(distance_function_interface) :: distance_func
        logical, intent(in) :: equal_size
        integer, intent(in) :: compute_initial_centroids, n_threads
        integer, optional, intent(in) :: iterations
        real, optional, intent(in) :: weight_equalizer_scaling, reassignment_tol, weight_tol
        !Variables
        integer :: n_data, n_dim, i, j, index, n_iterations, n_reassignments, max_iterations, cluster_index, old_cluster_index, destination_cluster_index, data_index, n_side, row
        real :: distance, min_distance, rand, mean_cluster_weight, distance_weight, distance_weight_scaling, max_weight_deviation, weight_tolerance, reassignment_tolerance, x_min, x_max, y_min, y_max, delta_x, delta_y, data_diameter
        integer, dimension(k) :: inital_point_indicies, cluster_size_order
        integer, allocatable, dimension(:) :: candidate_indicies, candidate_order, cluster_member_indicies, transfer_indicies
        real, dimension(k) :: cluster_weights
        real, allocatable, dimension(:) :: candidate_memberships

        max_iterations = 1000
        if(present(iterations)) max_iterations = iterations
        reassignment_tolerance=0.0
        if(present(reassignment_tol)) reassignment_tolerance = reassignment_tol

        n_data = size(data_points, dim=1)
        n_dim = size(data_points, dim=2)

        if(k <= 1) then
            print *, "At least 2 blocks are required. Exiting"
            return
        end if
        if(n_dim /= 2) then
            print *, "Only 2 dimensional data are supported for now. Exiting."
            return
        end if

        if(equal_size) then
            allocate(candidate_indicies(n_data), candidate_memberships(n_data), candidate_order(n_data))
            candidate_memberships = 0.0
        end if
        weight_tolerance=0.0
        if(present(weight_tol)) weight_tolerance = weight_tol
        if(present(weight_equalizer_scaling)) then
            distance_weight_scaling = weight_equalizer_scaling
        else
            distance_weight_scaling = 0.0
        end if
        mean_cluster_weight = sum(weights)/k

        x_min = minval(data_points(:,1))
        x_max = maxval(data_points(:,1))
        y_min = minval(data_points(:,2))
        y_max = maxval(data_points(:,2))
        data_diameter = sqrt((x_max-x_min)**2 + (y_max-y_min)**2)
        if(compute_initial_centroids == kmeans_compute_random_initial_centroids) then
            !Randomly set initial centroids by choosing k points from the data set
            inital_point_indicies = 0
            do j=1,k
                call random_number(rand)
                index = int(rand*n_data) + 1
                do while(count(inital_point_indicies == index) > 0)
                    call random_number(rand)
                    index = int(rand*n_data) + 1
                end do
                inital_point_indicies(j) = index
                cluster_centroids(j,:) = data_points(index,:)
            end do
            !Assign data points to random clusters
            call random_number_int(array=cluster_indicies, min=1, max=k)
            do j=1,k
                cluster_weights(j) = sum(weights, mask=cluster_indicies == j)
            end do
        else if(compute_initial_centroids == kmeans_compute_fixed_initial_centroids) then
            !Algorithm: if field is wider than high, divide field in n_side rows
            !x positions are then derived from dividing n_side*width_of_field by n_blocks_field
            n_side = sqrt(real(k))
            if(x_max - x_min > y_max - y_min) then
                delta_y = (y_max - y_min)/n_side
                delta_x = (x_max - x_min)*n_side/k
                do j=1, k
                    row = int((delta_x/2.0 + delta_x*(j-1))/(x_max - x_min))
                    cluster_centroids(j,1) = x_min + mod(delta_x/2.0 + delta_x*(j-1), x_max - x_min)
                    cluster_centroids(j,2) = y_min + delta_y/2.0 + delta_y*row
                end do
            else
                delta_x = (x_max - x_min)/n_side
                delta_y = (y_max - y_min)*n_side/k
                do j=1, k
                    row = int((delta_y/2.0 + delta_y*(j-1))/(y_max - y_min))
                    cluster_centroids(j,2) = y_min + mod(delta_y/2.0 + delta_y*(j-1), y_max - y_min)
                    cluster_centroids(j,1) = x_min + delta_x/2.0 + delta_x*row
                end do
            end if
            cluster_indicies = [(i, i=1,n_data)]
            cluster_indicies = real(cluster_indicies-1)/n_data*k + 1
            !cluster_indicies = mod(cluster_indicies, k) + 1
            do j=1,k
                cluster_weights(j) = sum(weights, mask=cluster_indicies == j)
            end do
        end if

        !TODO: Do separate initialization step

        n_iterations = 0
        !k-means algorithm, stop after iteration limit is reached or no new assignments have taken place
        do while(n_iterations <= max_iterations)
            n_iterations = n_iterations + 1
            n_reassignments = 0
            !Assign data point to clusters
            call assign_clusters(reassignment_tolerance)
            !print *, "Reassignments: ", n_reassignments
            !Find new centroids
            call calculate_centroids()

            max_weight_deviation = (maxval(cluster_weights) - minval(cluster_weights))/mean_cluster_weight

            if(n_reassignments == 0) then
                !No new assignments necessary, thus the solution has been found
                print *, "Number of iterations necessary to find solution:", n_iterations
                exit
            else if(real(n_reassignments)/n_data < reassignment_tolerance) then
                print *, "Reached reassignment tolerance. Iterations needed:", n_iterations
                exit
            end if
        end do

        if(n_iterations > max_iterations .and. n_reassignments > 0) then
            print *, "Maximum number of iterations reached. Number of reassignments in last pass:", n_reassignments
        end if

        if(max_weight_deviation > weight_tolerance) then
            print *, "Did not reach weight deviation tolerance. Maximum size deviation:", max_weight_deviation
        end if

        !Try to equalize cluster weights
        if(equal_size .and. max_weight_deviation > weight_tolerance) then
            call postprocess_weights(weight_tolerance)
            max_weight_deviation = (maxval(cluster_weights) - minval(cluster_weights))/mean_cluster_weight
            if(max_weight_deviation > weight_tolerance) then
                print *, "Did not reach weight deviation tolerance. Maximum size deviation:", max_weight_deviation
            end if
        end if

        if(equal_size) then
            deallocate(candidate_indicies, candidate_memberships)
            deallocate(candidate_order)
        end if

        contains
        subroutine assign_clusters(reassignment_tolerance)
            !Arguments
            real :: reassignment_tolerance

!$OMP PARALLEL DO default(none) num_threads(n_threads) shared(k, n_data, n_iterations, data_points, weights, cluster_centroids, distance_weight_scaling, mean_cluster_weight, equal_size, cluster_indicies, candidate_memberships, candidate_indicies, cluster_weights, n_reassignments, data_diameter) private(i, j, min_distance, cluster_index, old_cluster_index, distance, distance_weight)
            do i=1,n_data
                min_distance = data_diameter*10.0
                cluster_index = 0
                do j=1,k
                    if(n_iterations > 1 .and. distance_weight_scaling > 0.0) then
                        !Increase distances to overpopulated clusters
                        distance_weight = 1.0 + (cluster_weights(j)/mean_cluster_weight - 1.0)*distance_weight_scaling
                    else
                        !Find approximate clusters first
                        distance_weight = 1.0
                    end if
                    distance = distance_weight * distance_func(data_points(i,1), data_points(i,2), cluster_centroids(j,1), cluster_centroids(j,2))
                    if(equal_size .and. n_iterations > 1 .and. distance > 0.0) then
                        if(distance < min_distance) then
                            candidate_memberships(i) = min_distance/distance
                            candidate_indicies(i) = cluster_index
                            min_distance = distance
                            cluster_index = j
                        else if(distance < candidate_memberships(i)*min_distance) then
                            candidate_memberships(i) = distance/min_distance
                            candidate_indicies(i) = j
                        end if
                    else
                        if(distance < min_distance) then
                            min_distance = distance
                            cluster_index = j
                        end if
                    end if
                end do
                if(cluster_index /= cluster_indicies(i)) then
                    old_cluster_index = cluster_indicies(i)
                    cluster_indicies(i) = cluster_index
                    !$OMP ATOMIC UPDATE
                    cluster_weights(old_cluster_index) = cluster_weights(old_cluster_index) - weights(i)
                    !$OMP ATOMIC UPDATE
                    cluster_weights(cluster_index) = cluster_weights(cluster_index) + weights(i)
                    !$OMP ATOMIC UPDATE
                    n_reassignments = n_reassignments + 1
                end if
            end do
!$OMP END PARALLEL DO
        end subroutine assign_clusters

        subroutine calculate_centroids()
            cluster_centroids = 0.0
            do i=1,n_data
                cluster_index = cluster_indicies(i)
                cluster_centroids(cluster_index,:) = cluster_centroids(cluster_index,:) + data_points(i,1:n_dim)*weights(i)
            end do
            do j=1,k
                cluster_centroids(j,:) = cluster_centroids(j,:)/cluster_weights(j)
            end do
        end subroutine calculate_centroids

        subroutine postprocess_weights(weight_tolerance)
            !Arguments
            real :: weight_tolerance
            !Variables
            integer :: n_postprocess_reassignments

            n_postprocess_reassignments = 0
            print *, "Equalizing cluster weights."

            call qsort_arg(cluster_weights, cluster_size_order, reverse=.true.)
            call qsort_arg(candidate_memberships, candidate_order)

!            print *, mean_cluster_size
!            print *, cluster_size
!            print *, cluster_size_order

            do j=1,k
                cluster_index = cluster_size_order(j)
                cluster_member_indicies = pack(candidate_order, cluster_indicies(candidate_order) == cluster_index)
                transfer_indicies = pack(candidate_indicies(candidate_order), cluster_indicies(candidate_order) == cluster_index)
                !print *, "Cluster size: ", cluster_size(cluster_index)
                do i=1,size(cluster_member_indicies)
                    data_index = cluster_member_indicies(i)
                    destination_cluster_index = transfer_indicies(i)
                    !Check if there is a candidate cluster
                    if(destination_cluster_index == 0) cycle
                    !Only assign to lower weighted cluster in the list
                    if(count(cluster_size_order(1:j) == destination_cluster_index) > 0) then
                        cycle
                    end if
                    cluster_indicies(data_index) = destination_cluster_index
                    cluster_weights(cluster_index) = cluster_weights(cluster_index) - weights(data_index)
                    cluster_weights(destination_cluster_index) = cluster_weights(destination_cluster_index) + weights(data_index)
                    n_postprocess_reassignments = n_postprocess_reassignments + 1
                    if(cluster_weights(cluster_index) <= mean_cluster_weight*(1.0+weight_tolerance)) then
                        exit
                    end if
                end do
                !print *, "Adjusted cluster size: ", cluster_size(cluster_index)
            end do
            print *, "Number of reassignments in postprocessing: ", n_postprocess_reassignments
        end subroutine postprocess_weights
    end subroutine partition_data

end module kmeans_module