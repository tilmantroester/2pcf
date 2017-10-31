program test_kdtree
use kdtree_module
use utils
implicit none

type(kdtree) :: tree

integer , parameter :: N = 100000

real(fpp), allocatable, dimension(:,:) :: points
real(fpp), dimension(2) :: query
integer, allocatable, dimension(:) :: results, results_new, results_bruteforce, sort_array, sort_array_bf
type(search_record) :: dr
integer :: n_found, n_found_bruteforce, i
real(fpp) :: r, dummy
integer :: fileunit, n_trails
logical :: spherical_coords = .true.

integer(kind=8) :: clock_count_start, clock_count_end, elapsed_time
real :: clock_rate, average_time

allocate(points(N,2), results(N), results_new(N), results_bruteforce(N), sort_array(N), sort_array_bf(N))

query = [2.5, 0.25]
r = 1.1

call random_seed()

call random_number(points)

points(:,1) = points(:,1)*2.0*pi
points(:,2) = -pi/2.0 + points(:,2)*pi

!open(newunit=fileunit, file="../../../data/background.dat", status='old')
!read(unit=fileunit, fmt=*) ((points(i,1), points(i,2), dummy, dummy, dummy, dummy), i=1, N)
!close(fileunit)

call grow_tree(tree, points, leafsize=10)
allocate(dr%included_nodes(tree%n_nodes))

open(newunit=fileunit, file="/Users/yooken/Research/tpcf/tree_debug_nodes.txt", status="replace")
call write_node_data_to_file(tree, tree%root_ptr, fileunit)
close(fileunit)
open(newunit=fileunit, file="/Users/yooken/Research/tpcf/tree_debug_idx.txt", status="replace")
write(fileunit, fmt="(I)") (tree%indicies(i), i=1,N)
close(fileunit)
open(newunit=fileunit, file="/Users/yooken/Research/tpcf/tree_debug_pos.txt", status="replace")
write(fileunit, fmt="(2F15.3)") (points(i,:), i=1,N)
close(fileunit)

print *, "Total number of children:", tree%n_children
print *, "Number of nodes:", tree%n_nodes
print *, "Number of midpoint slides:", tree%n_mid_point_slides

if(spherical_coords) then
!    n_trails = 1000
!    elapsed_time = 0
!    do i=1,n_trails
!        call system_clock(count=clock_count_start)
!        call search_within_r_spherical(tree, query, r, results, n_found, .false.)
!        call system_clock(count=clock_count_end, count_rate=clock_rate)
!        elapsed_time = elapsed_time + clock_count_end-clock_count_start
!    end do
!    print *, "Run time old kdtree ", elapsed_time/clock_rate/n_trails
!    elapsed_time = 0
!    do i=1,n_trails
!        call system_clock(count=clock_count_start)
!        call kdtree_estimator(tree, query, r, results, n_found, .false.)
!        call system_clock(count=clock_count_end, count_rate=clock_rate)
!        elapsed_time = elapsed_time + clock_count_end-clock_count_start
!    end do
!    print *, "Run time new kdtree ", elapsed_time/clock_rate/n_trails

    call search_within_r_spherical(tree, query, r, results, n_found, .false., dr)
!    call kdtree_estimator(tree, query, r, results_new, n_found, .false., dr)

    call search_within_r_spherical(tree, query, r, results_bruteforce, n_found_bruteforce, .true.)
else
    call search_within_r(tree, query, r, results, n_found, .false., dr)
    !call search_within_r(tree, query, r, results_bruteforce, n_found_bruteforce, .true.)
end if

print *, "Number of points found within", r, " of", query," :", n_found
!print *, "Number of points found within", r, " of", query," :(brute force)", n_found_bruteforce

open(newunit=fileunit, file="/Users/yooken/Research/tpcf/tree_debug_search.txt", status="replace")
call write_search_debug_data_to_file(tree, dr, fileunit)
close(fileunit)
open(newunit=fileunit, file="/Users/yooken/Research/tpcf/tree_debug_search_idx.txt", status="replace")
write(fileunit, fmt="(I)") (results(i), i=1,n_found)
close(fileunit)

!if(all(results == results_new)) then
!    print *, "Results between old and new are identical."
!end if
call qsort_arg_int(results, sort_array)
call qsort_arg_int(results_bruteforce, sort_array_bf)
results = results(sort_array)
results_bruteforce = results_bruteforce(sort_array_bf)

if(all(results == results_bruteforce)) then
    print *, "Results are identical."
else
    print *, "tree search is broken."
    do i=1, N
        if(results(i) /= 0) then
            write(*, fmt="(I, 2F6.3, A, F6.3)") results(i), points(results(i), 1), points(results(i), 2), " distance: ", sqrt(sum((points(results(i),:)-query(:))**2))
        end if
    end do
    print *, "vs"
    do i=1,N
        if(results_bruteforce(i) /= 0) then
            write(*, fmt="(I, 2F6.3, A, F6.3)")results_bruteforce(i), points(results_bruteforce(i), 1), points(results_bruteforce(i), 2), " distance: ", sqrt(sum((points(results_bruteforce(i),:)-query(:))**2))
        end if
    end do
end if



call burn_tree(tree)

end program