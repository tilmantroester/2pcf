module utils
    use, intrinsic :: ISO_FORTRAN_ENV
    implicit none
    
    real, parameter :: pi = 3.141592653589793
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: lin_binning = 0, log_binning = 1, ihs_binning = 2, sqrt_binning = 3

    private
    public pi, assert_iostat, count_number_of_lines
    public lin_binning, log_binning, ihs_binning, sqrt_binning
    public binning, find_bin_index, create_bins, destroy_bins
    public qsort_arg, qsort_arg_int, qselect
    public distance_function_interface, distance_euclidean, deg2rad, rad2deg, distance_sphere, max_longitude_difference, angles_sphere, rotation_sphere, bearing_sphere, pointing_angle_cos_sin_sphere, distance_to_meridian_sphere, longitude_reach_sphere
    public random_number_int, random_unique_number_int, poisson
    public mean, median, variance, covariance_matrix

    type binning
        real :: x_min, x_max, bin_size = 0.0
        integer :: n_bin = 0
        integer :: spacing
        real :: ihs_theta
        real, allocatable :: bin_edges(:), bin_centers(:), bin_centers_eff(:)
    end type binning

    abstract interface
        function distance_function_interface(x1, y1, x2, y2) result(d)
            implicit none
            real, intent(in) :: x1, y1, x2, y2
            real :: d
        end function
    end interface

    contains

    subroutine assert_iostat(iostat, filename)
        !Arguments
        integer, intent(in) :: iostat
        character(len=*), intent(in) :: filename

        if(iostat /= 0) then
            print "(A, A, A, I5)", "Failed to open ", trim(filename), ". Iostat: ", iostat
            error stop
        end if
    end subroutine assert_iostat

    function count_number_of_lines(filename) result(n)
        implicit none
        character(len=*), intent(in) :: filename
        integer :: n, file_unit, iostat
        character :: c

        open(newunit=file_unit, file=filename, status='old', iostat=iostat)
        if(iostat > 0) then
            print *, "Failed to open file ", filename
            n = -1
            return
        end if

        n = 0
        do
            read(unit=file_unit, fmt=*, iostat=iostat) c
            if(iostat == 0) then
                if(c == '#') then
                    print *, "Deteced comment lines in file ", filename, ". Check your data file."
                    cycle
                end if

                n = n + 1
            else if(iostat < 0) then
                !print *, "Reached end of file."
                exit
            else
                print *, "Error reading file ", filename
                n = -1
                return
            end if
        end do

        close(file_unit)
    end function

    pure elemental function isinh_transform(x, p) result(y)
        !Argument
        real, intent(in) :: x, p
        !Return value
        real :: y

        y = log(p*x + sqrt((p*x)**2 + 1))/p
    end function isinh_transform

    pure elemental function sinh_transform(y, p) result(x)
        !Argument
        real, intent(in) :: y, p
        !Return value
        real :: x

        x = sinh(p*y)/p
    end function sinh_transform

    function find_bin_index(x, bins) result(i)
        real, intent(in) :: x
        type(binning), intent(in) :: bins
        !Variables
        integer :: i
        real :: bin_size, x_trans, x_min_trans

        if(x < bins%x_min .or. x > bins%x_max) then
            i = -1
            return
        else if(x == bins%x_max) then
            i = bins%n_bin
            return
        end if

        if(bins%spacing == lin_binning) then
            i = int((x - bins%x_min)/bins%bin_size) + 1
        else if(bins%spacing == log_binning) then
            x_trans = log10(x)
            x_min_trans = log10(bins%x_min)
            i = int((x_trans - x_min_trans)/bins%bin_size) + 1
        else if(bins%spacing == ihs_binning) then
            x_trans = isinh_transform(x, bins%ihs_theta)
            x_min_trans = isinh_transform(bins%x_min, bins%ihs_theta)
            i = int((x_trans - x_min_trans)/bins%bin_size) + 1
        else if(bins%spacing == sqrt_binning) then
            x_trans = sqrt(x)
            x_min_trans = sqrt(bins%x_min)
            i = int((x_trans - x_min_trans)/bins%bin_size) + 1
        end if

        if(i > bins%n_bin) then
            i = -1
            print *, "Something went wrong during binning."
            print "(A, F, A, F)", "x : ", x, "x_max : ", bins%x_max
        end if
    end function

    subroutine create_bins(x_min, x_max, bin_size, n_bin, p, spacing, bins)
        real, optional, intent(in) :: x_min, x_max, bin_size
        integer, optional, intent(in) :: n_bin
        real, optional, intent(in) :: p
        integer, intent(in) :: spacing
        type(binning), intent(out) :: bins
        !Variables
        integer :: i
        real :: theta

        bins%spacing = spacing
        if(.not. (spacing == lin_binning .or. &
                  spacing == log_binning .or. &
                  spacing == ihs_binning .or. &
                  spacing == sqrt_binning)) then
            print *, "Bin spacing unknown:", spacing
            error stop
        end if
        if(spacing == ihs_binning) then
            theta = 1.0
            if(present(p)) theta = p
            bins%ihs_theta = theta
        end if

        if(present(n_bin)) then
            bins%n_bin = n_bin
        else if(present(x_min) .and. present(x_max) .and. present(bin_size) .and. .not. present(n_bin)) then
            if(spacing == lin_binning) then
                bins%n_bin = (x_max - x_min)/bin_size
            else if(spacing == log_binning) then
                bins%n_bin = (log10(x_max) - log10(x_min))/bin_size
            else if(spacing == ihs_binning) then
                bins%n_bin = (isinh_transform(x_max, theta) - isinh_transform(x_min, theta))/bin_size
            else if(spacing == sqrt_binning) then
                bins%n_bin = (sqrt(x_max) - sqrt(x_min))/bin_size
            end if
        else
            print *, "Binning is ill-defined."
            error stop
        end if

        allocate(bins%bin_edges(bins%n_bin+1), bins%bin_centers(bins%n_bin))

        if(present(x_min) .and. present(x_max) .and. present(n_bin) .and. .not. present(bin_size)) then
            bins%x_min = x_min
            bins%x_max = x_max
            if(spacing == lin_binning) then
                bins%bin_size = (x_max - x_min)/n_bin
            else if(spacing == log_binning) then
                bins%bin_size = (log10(x_max) - log10(x_min))/n_bin
            else if(spacing == ihs_binning) then
                bins%bin_size = (isinh_transform(x_max, theta) - isinh_transform(x_min, theta))/n_bin
            else if(spacing == sqrt_binning) then
                bins%bin_size = (sqrt(x_max) - sqrt(x_min))/n_bin
            end if
        else if(present(x_min) .and. present(bin_size) .and. present(n_bin) .and. .not. present(x_max)) then
            bins%x_min = x_min
            bins%bin_size = bin_size
            if(spacing == lin_binning) then
                bins%x_max = x_min +  bin_size*n_bin
            else if(spacing == log_binning) then
                bins%x_max = 10**(log10(x_min) +  bin_size*n_bin)
            else if(spacing == ihs_binning) then
                bins%x_max = sinh_transform(isinh_transform(x_min, theta) +  bin_size*n_bin, theta)
            else if(spacing == sqrt_binning) then
                bins%x_max = (sqrt(x_min) +  bin_size*n_bin)**2
            end if
        else
            print *, "Binning is ill-defined."
            error stop
        end if

        bins%bin_edges(bins%n_bin+1) = bins%x_max
        if(spacing == lin_binning) then
            do i=1, bins%n_bin
                bins%bin_edges(i) = bins%x_min + bins%bin_size * (i-1)
            end do
        else if(spacing == log_binning) then
            do i=1, bins%n_bin
                bins%bin_edges(i) = 10.0**(log10(bins%x_min) + bins%bin_size * (i-1))
            end do
        else if(spacing == ihs_binning) then
            do i=1, bins%n_bin
                bins%bin_edges(i) = sinh_transform(isinh_transform(bins%x_min, theta) + bins%bin_size * (i-1), theta)
            end do
        else if(spacing == sqrt_binning) then
            do i=1, bins%n_bin
                bins%bin_edges(i) = (sqrt(bins%x_min) + bins%bin_size * (i-1))**2
            end do
        end if
        bins%bin_centers = bins%bin_edges(1:bins%n_bin) + (bins%bin_edges(2:bins%n_bin+1) - bins%bin_edges(1:bins%n_bin))/2.0
    end subroutine create_bins

    subroutine destroy_bins(bins)
        type(binning) :: bins

        if(allocated(bins%bin_edges)) deallocate(bins%bin_edges)
        if(allocated(bins%bin_centers)) deallocate(bins%bin_centers)
    end subroutine destroy_bins

    subroutine qsort_arg(array, order, reverse)
        !Arguments
        real, dimension(:), intent(in) :: array
        integer, dimension(:), intent(out) :: order
        logical, optional, intent(in) :: reverse
        !Variables
        real, allocatable, dimension(:) :: work_array
        integer :: i

        do i=1,size(array)
            order(i) = i
        end do

        allocate(work_array(size(array)), source=array)

        call qsort_func(1, size(array))

        if(present(reverse)) then
            if(reverse) order = order(size(array):1:-1)
        end if

        deallocate(work_array)
        
        contains
        recursive subroutine qsort_func(l, r)
            integer, intent(in) :: l, r
            integer :: i, j, tmp_idx
            real :: tmp, pivot

            pivot = work_array((l+r)/2)
            i = l-1
            j = r+1

            do
                do
                    j = j-1
                    if (work_array(j) <= pivot) exit
                end do
                do
                    i = i+1
                    if (work_array(i) >= pivot) exit
                end do
                if (i < j) then
                    tmp = work_array(i)
                    work_array(i) = work_array(j)
                    work_array(j) = tmp
                    tmp_idx = order(i)
                    order(i) = order(j)
                    order(j) = tmp_idx
                else if(i == j) then
                    i = i + 1
                    exit
                else
                    exit
                endif
            end do

            if(l < j) call qsort_func(l, j)
            if(i < r) call qsort_func(i, r)
        end subroutine qsort_func
    end subroutine qsort_arg

    subroutine qsort_arg_int(array, order, reverse)
        !Arguments
        integer, dimension(:), intent(in) :: array
        integer, dimension(:), intent(out) :: order
        logical, optional, intent(in) :: reverse
        !Variables
        integer, allocatable, dimension(:) :: work_array
        integer :: i

        do i=1,size(array)
            order(i) = i
        end do

        allocate(work_array(size(array)), source=array)

        call qsort_func(1, size(array))

        if(present(reverse)) then
            if(reverse) order = order(size(array):1:-1)
        end if

        deallocate(work_array)
        
        contains
        recursive subroutine qsort_func(l, r)
            integer, intent(in) :: l, r
            integer :: i, j, tmp_idx
            integer :: tmp, pivot

            pivot = work_array((l+r)/2)
            i = l-1
            j = r+1

            do
                do
                    j = j-1
                    if (work_array(j) <= pivot) exit
                end do
                do
                    i = i+1
                    if (work_array(i) >= pivot) exit
                end do
                if (i < j) then
                    tmp = work_array(i)
                    work_array(i) = work_array(j)
                    work_array(j) = tmp
                    tmp_idx = order(i)
                    order(i) = order(j)
                    order(j) = tmp_idx
                else if(i == j) then
                    i = i + 1
                    exit
                else
                    exit
                endif
            end do

            if(l < j) call qsort_func(l, j)
            if(i < r) call qsort_func(i, r)
        end subroutine qsort_func
    end subroutine qsort_arg_int

    pure function qselect(array, k) result(k_value)
        !Arguments
        real, dimension(:), intent(in) :: array
        integer, intent(in) :: k
        !Return value
        real :: k_value
        !Variables
        real, dimension(:), allocatable :: work_array

        allocate(work_array(size(array)), source=array)
        call quickselect(work_array, 1, size(array), k, k_value)

        contains
            pure recursive subroutine quickselect(array, l, r, k, k_value)
                !Arguments
                real, dimension(:), intent(inout) :: array
                integer, intent(in) :: l, r, k
                real, intent(out) :: k_value
                !Variables
                integer :: pivot_index

                do while(.true.)
                    if(l == r) then
                        k_value = array(l)
                        exit
                    end if
                    call partition(array, l, r, int((l+r)/2), pivot_index)
                    if(k == pivot_index) then
                        k_value = array(k)
                        exit
                    else if( k < pivot_index) then
                        call quickselect(array, l, r-1, k, k_value)
                        exit
                    else
                        call quickselect(array, l+1, r, k, k_value)
                        exit
                    end if
                end do
            end subroutine

            pure subroutine swap(array, i, j)
                !Arguments
                real, dimension(:), intent(inout) :: array
                integer, intent(in) :: i, j
                !Variables
                real :: tmp

                tmp = array(j)
                array(j) = array(i)
                array(i) = tmp
            end subroutine swap

            pure subroutine partition(array, l, r, pivot_index, k)
                !Arguments
                real, dimension(:), intent(inout) :: array
                integer, intent(in) :: l, r, pivot_index
                integer, intent(out) :: k
                !Variables
                real :: pivot_value
                integer :: i

                pivot_value = array(pivot_index)
                call swap(array, pivot_index, r)
                k = l
                do i=l, r-1
                    if(array(i) < pivot_value) then
                        call swap(array, i, k)
                        k = k+1
                    end if
                end do
                call swap(array, r, k)
            end subroutine partition
    end function qselect

    function poisson(mu) result(x)
        !Arguments
        real, intent(in) :: mu
        !Result
        integer :: x
        !Variables
        real :: sum, prod, u

        x = 0
        prod = exp(-mu)
        sum = prod
        call random_number(u)
        do while (u > sum)
            x = x + 1
            prod = prod * mu/x
            sum = sum + prod
        end do
    end function poisson

    pure elemental function distance_euclidean(x1, y1, x2, y2) result(d)
        !Arguments
        real, intent(in) :: x1, y1, x2, y2
        !Result
        real :: d

        d = sqrt((x1-x2)**2 + (y1-y2)**2)
    end function distance_euclidean

    pure elemental function deg2rad(deg) result(rad)
        !Argument
        real, intent(in) :: deg
        !Result
        real :: rad

        rad = deg/180.0*pi
    end function deg2rad

    pure elemental function rad2deg(rad) result(deg)
        !Argument
        real, intent(in) :: rad
        !Result
        real :: deg

        deg = rad*180.0/pi
    end function rad2deg

    pure elemental subroutine rotation_sphere(lambda, phi, lambda_rot, phi_rot, lambda_out, phi_out)
        !Arguments
        real, intent(in) :: lambda, phi, lambda_rot, phi_rot
        real, intent(out) :: lambda_out, phi_out
        !Varaibles
        real :: sigma, theta, theta_rot, x, y, z, tan_lambda_out, cos_phi_out

        !Azimuthal angle: phi, ploar angle: theta
        sigma = lambda + lambda_rot
        theta = pi/2.0 - phi
        theta_rot = -phi_rot
        x = cos(sigma)*cos(theta_rot)*sin(theta) + cos(theta)*sin(theta_rot)
        y = sin(sigma)*sin(theta)
        z = -cos(sigma)*sin(theta_rot)*sin(theta) + cos(theta)*cos(theta_rot)
        lambda_out = atan(y/x)
        if(x < 0.0) then
            lambda_out = lambda_out + pi
        end if

        if(z >= 1.0) then
            phi_out = pi/2.0
        else if(z <= -1.0) then
            phi_out = -pi/2.0
        else
            phi_out = pi/2.0 - acos(z)
        end if
    end subroutine rotation_sphere


    elemental function distance_sphere(lambda1, phi1, lambda2, phi2) result(d)
        !Arguments
        real, intent(in) :: phi1, lambda1, phi2, lambda2
        !Result
        real :: d
        !Variables
        real :: cos_d

        cos_d = sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda2-lambda1)
        if(cos_d >= 1.0) then
            d = 0.0
        else if(cos_d <= -1.0) then
            d = pi
        else
            d = acos(cos_d)
        end if
    end function distance_sphere

    elemental function max_longitude_difference(lambda1, phi1, d) result(Delta_lambda)
        !Arguments
        real, intent(in) :: phi1, lambda1, d
        !Result
        real :: Delta_lambda
        !Variables
        real :: cos_d, sin_phi1, sin_phi2, cos_phi2, cos_Delta_lambda

        cos_d = cos(d)
        sin_phi1 = sin(phi1)
        sin_phi2 = sin_phi1/cos_d
        cos_phi2 = sqrt(1.0 - sin_phi2**2)
        cos_Delta_lambda = (cos_d - sin_phi1*sin_phi2)/(cos(phi1)*cos_phi2)
        if(cos_Delta_lambda >= 1.0) then
            Delta_lambda = 0.0
        else if(cos_d <= -1.0) then
            Delta_lambda = pi
        else
            Delta_lambda = acos(cos_Delta_lambda)
        end if
    end function max_longitude_difference

    elemental function bearing_sphere(lambda1, phi1, lambda2, phi2) result(alpha)
        !Arguments
        real, intent(in) :: phi1, lambda1, phi2, lambda2
        !Result
        real :: alpha
        !Variables
        real :: x, y, cos_phi2

        cos_phi2 = cos(phi2)
        y = cos_phi2 * sin(lambda2-lambda1)
        x = cos(phi1)*sin(phi2) - sin(phi1)*cos_phi2*cos(lambda2-lambda1)
        alpha = atan2(y,x)
    end function bearing_sphere

    pure function pointing_angle_cos_sin_sphere(lambda1, phi1, lambda2, phi2, distance) result(cos_sin)
        !Arguments
        real, intent(in) :: phi1, lambda1, phi2, lambda2, distance
        !Result
        real :: cos_sin(2)
        !Variables
        real :: sin_d

        sin_d = sin(distance)
        cos_sin(1) = sin(lambda2-lambda1)*cos(phi2)/sin_d
        cos_sin(2) = (sin(phi2)-sin(phi1)*cos(distance))/(cos(phi1)*sin_d)
    end function pointing_angle_cos_sin_sphere

    elemental function distance_to_meridian_sphere(lambda1, phi1, lambda2) result(distance)
        !Arguments
        real, intent(in) :: phi1, lambda1, lambda2
        !Result
        real :: distance
        !Variables
        real :: sin_distance

        sin_distance = abs(sin(lambda2-lambda1)*cos(phi1))
        if(abs(lambda2-lambda1) <= pi/2.0) then
            !Meridian is closer than poles
            if(sin_distance >= 1.0) then
                distance = pi/2.0
            else
                distance = asin(sin_distance)
            end if
        else
            !Poles are closer
            distance = pi/2.0 - abs(phi1)
        end if
    end function distance_to_meridian_sphere

    elemental function longitude_reach_sphere(phi1, distance) result(Delta_lambda)
        !Arguments
        real, intent(in) :: phi1, distance
        !Result
        real :: Delta_lambda
        !Variables
        real :: sin_Delta_lambda

        sin_Delta_lambda = abs(sin(distance)/cos(phi1))
        if(sin_Delta_lambda >= 1.0) then
            Delta_lambda = pi/2.0
        else
            Delta_lambda = asin(sin_Delta_lambda)
        end if
    end function longitude_reach_sphere


    pure function angles_sphere(lambda1, phi1, lambda2, phi2, distance) result(cos_sin)
        !Arguments
        real, intent(in) :: phi1, lambda1, phi2, lambda2
        real, optional, intent(in) :: distance
        !Result
        real :: cos_sin(2)
        !Variables
        real :: d, sin_d, cos_phi2

        if(present(distance)) then
            d = distance
        else
            d = distance_sphere (lambda1, phi1, lambda2, phi2)
        end if

        sin_d = sin(d)
        cos_phi2 = cos(phi2)
        !cos theta
        cos_sin(1) = cos_phi2/sin_d * sin(lambda2-lambda1)
        !sin theta
        cos_sin(2) = (cos(phi1)*sin(phi2) - sin(phi1)*cos_phi2*cos(lambda2-lambda1))/sin_d
    end function angles_sphere

    subroutine random_number_int(array, max, min)
        !Arguments
        integer, intent(out) :: array(:)
        integer, intent(in) :: max
        integer, optional, intent(in) :: min
        !Variables
        integer :: i, min_number
        real :: rand

        if(present(min)) then
            min_number = min
        else
            min_number = 0
        endif

        do i=1,size(array)
            call random_number(rand)
            array(i) = int((max-min_number+1)*rand) + min_number
        end do
    end subroutine random_number_int

    subroutine random_unique_number_int(array, max, min)
        !Arguments
        integer, intent(out) :: array(:)
        integer, intent(in) :: max
        integer, optional, intent(in) :: min
        !Variables
        integer :: min_number, N, n_sample, m, t
        real :: u

        if(present(min)) then
            min_number = min
        else
            min_number = 0
        endif

        !Algorithm 3.4.2S from Knuth's book Seminumeric Algorithms
        N = max + 1 - min_number
        n_sample = size(array)
        if(n_sample > N - 1) then
            print "(A, I, A, I, A, I, A)", "Cannot draw ", n_sample, " unique numbers from the interval [", min_number, ",", max,"]."
            return
        end if
        m = 0
        t = 0

        do while(m < n_sample)
            call random_number(u)
            if( (N-t)*u >= n_sample-m) then
                t = t+1
            else
                array(m+1) = t
                t = t+1
                m = m+1
            end if
        end do
        array = array + min_number
    end subroutine random_unique_number_int

    pure function mean(x) result(m)
        !Arguments
        real, dimension(:), intent(in) :: x
        !Returns
        real :: m
        !Variables
        real :: n

        n = size(x)
        m = sum(x)/n
    end function mean

    pure function median(x) result(m)
        !Arguments
        real, dimension(:), intent(in) :: x
        !Returns
        real :: m
        !Variables
        integer :: n

        n = size(x)
        if(modulo(n, 2) == 0) then
            m = (qselect(x, n/2+1) + qselect(x, n/2))/2
        else
            m = qselect(x, n/2+1)
        end if
    end function median

    pure function variance(x) result(var)
        !Arguments
        real, dimension(:), intent(in) :: x
        !Returns
        real :: var
        !Variables
        real :: n, mean

        n = size(x)
        if(n <= 1) then
            var = 0.0
            return
        end if
        mean = sum(x)/n
        var = sum((x-mean)**2)/(n-1)
    end function variance

    subroutine covariance_matrix(x, cov)
        !Arguments
        real, intent(in) :: x(:,:)
        real, intent(out) :: cov(:,:)
        !Variables
        integer :: i, j, n_samples, n_data
        real :: mean_i, mean_j

        n_samples = size(x, dim=1)
        n_data = size(x, dim=2)

        do i=1,n_data
            mean_i = sum(x(:,i))/n_samples
            do j=1,i
                mean_j = sum(x(:,j))/n_samples
                cov(i,j) = sum((x(:,i)-mean_i)*(x(:,j)-mean_j))/(n_samples-1)
                cov(j,i) = cov(i,j)
            end do
        end do
    end subroutine covariance_matrix

end module utils