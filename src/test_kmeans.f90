program test_kmeans
use kmeans_module
use utils
implicit none

integer , parameter :: N = 100000, k = 15, n_masks = 4

real, dimension(N,2) :: points
real, dimension(k,2) :: centroids
logical, dimension(N) :: point_mask
real, allocatable, dimension(:,:) :: masked_points, pixels
real, allocatable, dimension(:) :: weights, x, y

integer, allocatable, dimension(:) :: cluster_indicies
character(len=256) :: debug_output_filename, map_filename, marks_filename
integer :: fileunit, i, j, n_points, n_seed, n_lines, iostat
integer :: dummy_int

real, dimension(n_masks, 3) :: masks
logical :: equal_size

integer, allocatable, dimension(:) :: seed
equal_size = .true.
masks = transpose(reshape([0.54, 0.37, 0.21, 0.23, 0.91, 0.15, 0.8, 0.12, 0.2, 0.24, 0.05, 0.13], shape(transpose(masks))))

!call random_seed()
!call random_seed(size=n_seed)
!allocate(seed(n_seed))
!call random_seed(get=seed)

!map_filename = "/Users/yooken/Research/tpcf/data/CDE0133_gamma_map.dat"
!marks_filename = "/Users/yooken/Research/tpcf/data/CDE0133_marks.dat"
!n_lines = count_number_of_lines(map_filename)
!n_points = count_number_of_lines(marks_filename)
!
!allocate(pixels(n_lines, 4), weights(n_points), cluster_indicies(n_points))
!
!open(newunit=fileunit, file=map_filename, status='old', iostat=iostat)
!read(unit=fileunit, fmt=*, iostat=iostat) ((pixels(i,1), pixels(i,2), pixels(i,3), pixels(i,4)), i=1, n_lines)
!close(fileunit)
!
!open(newunit=fileunit, file=marks_filename, status='old', iostat=iostat)
!read(unit=fileunit, fmt=*, iostat=iostat) ((weights(i), dummy_int), i=1, n_points)
!close(fileunit)
!
!if( n_points /= count(pixels(:,3) >= 0.0)) then
!    print *, "size mismatch:", n_points, count(pixels(:,3) >= 0.0)
!    stop
!end if
!
!print *, count(weights > 0.0)
!
!!masked_points = pack(pixels(:,1:2), mask=spread(pixels(:,3) >= 0.0, dim=2, ncopies=2))
!x = pack(pixels(:,1), mask=pixels(:,3) >= 0.0)
!y = pack(pixels(:,2), mask=pixels(:,3) >= 0.0)

!======== real marks ======
!marks_filename = "/Users/yooken/Research/tpcf/data/lenses_groups.dat"
!n_points = count_number_of_lines(marks_filename)
!
!allocate(x(n_points), y(n_points), weights(n_points), cluster_indicies(n_points))
!
!open(newunit=fileunit, file=marks_filename, status='old', iostat=iostat)
!read(unit=fileunit, fmt=*, iostat=iostat) ((x(i), y(i), weights(i)), i=1, n_points)
!close(fileunit)
!
!x = pack(x, mask=weights > 0.0)
!y = pack(y, mask=weights > 0.0)
!n_points = size(x)
!weights = pack(weights, mask=weights > 0.0)
!
!masked_points = reshape([x,y], shape=[n_points, 2])
!weights = 0.01 + weights/maxval(weights)
!======= end real marks ======

!====== random points ========
call random_number(centroids)
centroids(:,1) = 30.0 + 10.0*centroids(:,1)
centroids(:,2) = -12.0 + 10.0*centroids(:,2)

call random_number(points)

point_mask = .true.
do i=1,n_masks
    do j=1,N
        if(distance_euclidean(points(j,1), points(j,2), masks(i,1), masks(i,2)) < masks(i,3)) then
            point_mask(j) = .false.
        end if
    end do
end do
!point_mask = .true.
n_points = count(point_mask)
print *, "n_points: ", n_points

masked_points = reshape(pack(points, mask=spread(source=point_mask, dim=2, ncopies=2)), shape=[n_points,2])

allocate(weights(n_points), cluster_indicies(n_points))
call random_number(weights)
do i=1,n_points
    weights(i) = weights(i)*exp(-distance_euclidean(masked_points(i,1), masked_points(i,2), 0.5, 0.5)**2/0.1)
end do
weights = 1.0

!====== end random points ========

!call random_seed(put=seed)
!call partition_data(masked_points, weights, k, cluster_indicies, centroids, distance_euclidean, equal_size=.true., &
!                    compute_initial_centroids=kmeans_compute_fixed_initial_centroids, weight_equalizer_scaling=5.0, n_threads=4, iterations=1000, reassignment_tol=0.001, weight_tol=0.05)
!
call partition_data(masked_points, weights, k, cluster_indicies, centroids, distance_euclidean, equal_size=.false., &
                    compute_initial_centroids=kmeans_compute_fixed_initial_centroids, weight_equalizer_scaling=0.0, n_threads=4, iterations=1000, reassignment_tol=0.0, weight_tol=0.0)


open(newunit=fileunit, file="/Users/yooken/Research/tpcf/kmeans_debug_data.txt", status="replace")
write(fileunit, fmt="(F, F, F, I)") ((masked_points(i,:), weights(i), cluster_indicies(i)), i=1,n_points)
close(fileunit)
open(newunit=fileunit, file="/Users/yooken/Research/tpcf/kmeans_debug_centroids.txt", status="replace")
write(fileunit, fmt="(F, F)") ((centroids(i,:)), i=1,k)
close(fileunit)

!do j=0,30
!    call random_seed(put=seed)
!
!    call partition_data(masked_points, k, cluster_indicies, centroids, distance_euclidean, .false., equal_size, j, size_tol=0.05)
!
!    write(debug_output_filename, fmt="(A,I0,A)") "/Users/yooken/Research/tpcf/kmeans_debug_data_step_", j, ".txt"
!    open(newunit=fileunit, file=debug_output_filename, status="replace")
!    write(fileunit, fmt="(F, F, I)") ((masked_points(i,:), cluster_indicies(i)), i=1,n_points)
!    close(fileunit)
!    write(debug_output_filename, fmt="(A,I0,A)") "/Users/yooken/Research/tpcf/kmeans_debug_centroids_step_", j, ".txt"
!    open(newunit=fileunit, file=debug_output_filename, status="replace")
!    write(fileunit, fmt="(F, F)") ((centroids(i,:)), i=1,k)
!    close(fileunit)
!end do

end program