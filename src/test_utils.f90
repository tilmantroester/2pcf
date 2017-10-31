program test_utils
use utils
implicit none

!real :: mu
!integer :: i
!
!mu = 7.0
!
!do i=1,100
!    write(*, fmt="(A, I0)") "test:", poisson(mu)
!end do

integer :: i, j, n_bin, n_trails
integer :: int_array(6), order_array(6), count_array(20)
real :: real_array(6), array_copy(10000), m
real :: phi1, phi2, lambda1, lambda2, d, theta_min, theta_max, phi_rot, lambda_rot
real, dimension(10000) :: phi1_array, lambda1_array, phi2_array, lambda2_array, distance_array, lon_reach1_array, lon_reach2_array
real, dimension(10000, 2) :: cos_sin_array1, cos_sin_array2

real, dimension(:), allocatable :: big_real_array

type(binning) :: bins

integer(kind=8) :: clock_count_start, clock_count_end, elapsed_time
real :: clock_rate, average_time

d = 0.2
phi1 = 0.5
lambda1 = 0.5

print *, max_longitude_difference(lambda1, phi1, d)

n_bin = 10
theta_min = 1.0/60.0
theta_max = 20.0/60.0

theta_min = deg2rad(theta_min)
theta_max = deg2rad(theta_max)

call create_bins(x_min=theta_min, x_max=theta_max, n_bin=n_bin, spacing=lin_binning, bins=bins)

!print *, bins%bin_edges
!print *, bins%bin_centers

call random_seed()

call random_number_int(array=int_array, max=10, min=0)
print *, int_array
call qsort_arg_int(int_array, order_array)
print *, "sorted array: ", int_array(order_array)

count_array = 0
do i=1, 10000
    call random_unique_number_int(int_array, max=20, min=1)
    count_array(int_array) = count_array(int_array) + 1
end do
print *, "counts ", count_array


call random_number(real_array)
real_array = real_array*10

n_trails = 10
elapsed_time = 0
do i=1,n_trails
    call random_number(real_array)

    call system_clock(count=clock_count_start)
    call qsort_arg(real_array, order_array)
    !write(*,fmt="(10F3.2)") real_array(order_array)
    call system_clock(count=clock_count_end, count_rate=clock_rate)
    elapsed_time = elapsed_time + clock_count_end-clock_count_start
end do
average_time = real(elapsed_time)/n_trails
print *, "Time elapsed for copy qsort: ", average_time/clock_rate


print *, "Variance: ", variance(real_array)

phi1 = 0.3
lambda1 = pi*1.5
phi_rot = pi/2
lambda_rot = pi*0.6
call rotation_sphere(lambda1, phi1, lambda_rot, phi_rot, lambda2, phi2)
print "(A, F5.2, A, F5.2, A, F5.2, A)", "lambda rot: ", lambda_rot/pi, " : ", lambda1/pi, " -> ", lambda2/pi
print "(A, F5.2, A, F5.2, A, F5.2, A)", "phi rot: ", phi_rot/pi, " : ", phi1/pi, " -> ", phi2/pi

print *, "Checking pointing angles"
call random_number(phi1_array)
call random_number(phi2_array)
call random_number(lambda1_array)
call random_number(lambda2_array)
phi1_array = phi1_array * 2.0*pi
phi2_array = phi2_array * 2.0*pi
lambda1_array = lambda1_array *pi - pi/2.0
lambda2_array = lambda2_array *pi - pi/2.0

distance_array = distance_sphere(lambda1_array, phi1_array, lambda2_array, phi2_array)
n_trails = 100
elapsed_time = 0
do j=1,n_trails
    do i=1,size(distance_array)
        call system_clock(count=clock_count_start)
        cos_sin_array1(i,:) = angles_sphere(lambda1_array(i), phi1_array(i), lambda2_array(i), phi2_array(i), distance_array(i))
        call system_clock(count=clock_count_end, count_rate=clock_rate)
        elapsed_time = elapsed_time + clock_count_end-clock_count_start
    end do
end do
average_time = real(elapsed_time)/n_trails
print *, "Time elapsed for old pointing angle: ", average_time/clock_rate
elapsed_time = 0
do j=1,n_trails
    do i=1,size(distance_array)
        call system_clock(count=clock_count_start)
        cos_sin_array2(i,:) = pointing_angle_cos_sin_sphere(lambda1_array(i), phi1_array(i), lambda2_array(i), phi2_array(i), distance_array(i))
        call system_clock(count=clock_count_end, count_rate=clock_rate)
        elapsed_time = elapsed_time + clock_count_end-clock_count_start
    end do
end do
average_time = real(elapsed_time)/n_trails
print *, "Time elapsed for new pointing angle: ", average_time/clock_rate

lon_reach1_array = max_longitude_difference(lambda1_array, phi1_array*0.1, distance_array(1)*0.1)
lon_reach2_array = longitude_reach_sphere(phi1_array*0.1, distance_array(1)*0.1)

print *, sum(abs(cos_sin_array1-cos_sin_array2))/size(distance_array)
print *, sum(abs(lon_reach1_array-lon_reach2_array))/size(distance_array)
print *, sum(distance_array)/size(distance_array)

print *, distance_to_meridian_sphere(pi, -0.3, 0.0)
print *, distance_to_meridian_sphere(3.0, 0.0, 2.0*pi)

print "(3(A, F4.2))", "Mean : ", mean(phi1_array), " Median : ", median(phi1_array), " Variance : ", variance(phi1_array)

real_array = [(i, i=1,size(real_array))]
real_array = [1.0, 1.0, 1.0, 2.0, 3.0, 2.0]

print "(3(A, F4.2))", "Mean : ", mean(real_array), " Median : ", median(real_array), " Variance : ", variance(real_array)

allocate(big_real_array(300000))
call random_number(big_real_array)
print *, big_real_array(:10)
print "(A, F4.2)", "Median : ", median(big_real_array)

n_trails = 1
elapsed_time = 0
do j=1,n_trails
    call system_clock(count=clock_count_start)
    m = median(big_real_array)
    call system_clock(count=clock_count_end, count_rate=clock_rate)
    elapsed_time = elapsed_time + clock_count_end-clock_count_start
end do
average_time = real(elapsed_time)/n_trails
print *, "Time elapsed for median: ", average_time/clock_rate

end program test_utils