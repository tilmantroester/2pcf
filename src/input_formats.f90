module file_formats
    use utils
    use tpcf
    implicit none

    public load_standard_shape_catalog, load_standard_kappa_map, load_standard_lens_catalog, load_KiDS_shear_catalog, load_KiDS_lens_catalog, load_survey_filelist, load_survey_fits_filelist, save_marks, load_survey_marks, load_marks_fits_file, create_marks_fits_file, standard_scalar_foreground_data, standard_shape_data, standard_kappa_data

    private

    integer, parameter :: standard_scalar_foreground_data = 1, standard_shape_data = 2, standard_kappa_data = 3
    character(len=*), dimension(3), parameter :: data_type_name = ["lenses", "shapes", "convergence"]
    contains

    subroutine assert_fits_status(iounit, status)
        !Argument
        integer, intent(in) :: iounit, status
        !Variable
        character(len=80) :: comment
        integer :: tmp_status
        tmp_status = 0
        if (status /= 0) then
            print "(A, I4)", "FITS error status: ", status
            comment = "x"
            do while(trim(comment) /= "")
                call ftgmsg(comment)
                print "(A)", comment
            end do
            !Close fits file
            call ftclos(iounit, tmp_status)
            !Deallocate io unit
            call ftfiou(iounit, tmp_status)
            error stop
        end if
    end subroutine assert_fits_status

    function trim_fits_strings(str) result(trimmed)
        !Arguments
        character(len=*), intent(in) :: str
        !Return value
        character(len=len_trim(str(2:len_trim(str)-1))) :: trimmed

        trimmed = trim(str(2:len_trim(str)-1))
        return
    end function trim_fits_strings


    subroutine load_lensing_fits_catalog(filename, field, apply_mask, verbose)
        !Arguments
        character(len=256), intent(in) :: filename
        type(lensing_survey_field), intent(inout) :: field
        logical, intent(in) :: apply_mask
        logical, optional, intent(in) :: verbose
        !Variables
        logical :: verbose_output
        integer :: status, fits_iounit, rwmode, hdutype, blocksize, n, repeat, width, data_type
        character(len=80) :: survey_name, field_name, comment
        real :: nullval, mag_min, mag_max, scale
        logical :: anyf
        integer(kind=4) :: int32
        real(kind=4) :: float32

        verbose_output = .false.
        if(present(verbose)) verbose_output = verbose

        status = 0
        nullval = 0.0

        !Allocate io unit
        call ftgiou(fits_iounit, status)
        !Open fits file
        rwmode = 0
        call ftdkopn(fits_iounit, filename, rwmode, blocksize, status)
        call assert_fits_status(fits_iounit, status)
        !Read primary header information
        call ftgkey(fits_iounit, "SURVEY", survey_name, comment, status)
        call ftgkey(fits_iounit, "FIELD", field_name, comment, status)
        call ftgkyj(fits_iounit, "DATATYPE", int32, comment, status)
        data_type = int32
        call assert_fits_status(fits_iounit, status)
        call ftgkye(fits_iounit, "MAGMIN", float32, comment, status)
        if(status == 202) then
            mag_min = -100.0
            status = 0
        else
            mag_min = float32
        end if
        call ftgkye(fits_iounit, "MAGMAX", float32, comment, status)
        if(status == 202) then
            mag_max = -100.0
            status = 0
        else
            mag_max = float32
        end if
        call ftgkye(fits_iounit, "SCALE", float32, comment, status)
        if(status == 202) then
            scale = -1.0
            status = 0
        else
            scale = float32
        end if

        !Select data HDU
        call ftmahd(fits_iounit, 2, hdutype, status)
        call assert_fits_status(fits_iounit, status)
        !Read data HDU header information
        call FTGNRW(fits_iounit, n, status)
        call assert_fits_status(fits_iounit, status)

        if(verbose_output) print "(7A, I7)", "Loading ", trim(data_type_name(data_type)), ". Field ", trim_fits_strings(field_name), " of ", trim_fits_strings(survey_name), ". Number of objects: ", n

        if(data_type == standard_kappa_data) then
            !kappa
            field%convergence%n = n
            allocate(field%convergence%x(n), field%convergence%y(n), field%convergence%kappa_E(n), field%convergence%kappa_B(n), field%convergence%w(n))

            call read_fits_column_real(fits_iounit, "x", field%convergence%x)
            call read_fits_column_real(fits_iounit, "y", field%convergence%y)
            call read_fits_column_real(fits_iounit, "kappa_E", field%convergence%kappa_E)
            call read_fits_column_real(fits_iounit, "kappa_B", field%convergence%kappa_B)
            call read_fits_column_real(fits_iounit, "w", field%convergence%w)
            call assert_fits_status(fits_iounit, status)
            if(apply_mask) then
                field%convergence%x = pack(field%convergence%x, mask=field%convergence%w>0.0)
                field%convergence%y = pack(field%convergence%y, mask=field%convergence%w>0.0)
                field%convergence%kappa_E = pack(field%convergence%kappa_E, mask=field%convergence%w>0.0)
                field%convergence%kappa_B = pack(field%convergence%kappa_B, mask=field%convergence%w>0.0)
                field%convergence%w = pack(field%convergence%w, mask=field%convergence%w>0.0)
                field%convergence%n = size(field%convergence%x)
            end if
        else if(data_type == standard_shape_data) then
            !Shapes
            field%shapes%n = n
            allocate(field%shapes%x(n), field%shapes%y(n), field%shapes%z(n), field%shapes%e1(n), field%shapes%e2(n), field%shapes%m(n), field%shapes%w(n))

            call read_fits_column_real(fits_iounit, "x", field%shapes%x)
            call read_fits_column_real(fits_iounit, "y", field%shapes%y)
            call read_fits_column_real(fits_iounit, "z", field%shapes%z)
            call read_fits_column_real(fits_iounit, "e1", field%shapes%e1)
            call read_fits_column_real(fits_iounit, "e2", field%shapes%e2)
            call read_fits_column_real(fits_iounit, "m", field%shapes%m)
            call read_fits_column_real(fits_iounit, "w", field%shapes%w)
            call assert_fits_status(fits_iounit, status)
            if(apply_mask) then
                field%shapes%x = pack(field%shapes%x, mask=field%shapes%w>0.0)
                field%shapes%y = pack(field%shapes%y, mask=field%shapes%w>0.0)
                field%shapes%z = pack(field%shapes%z, mask=field%shapes%w>0.0)
                field%shapes%e1 = pack(field%shapes%e1, mask=field%shapes%w>0.0)
                field%shapes%e2 = pack(field%shapes%e2, mask=field%shapes%w>0.0)
                field%shapes%m = pack(field%shapes%m, mask=field%shapes%w>0.0)
                field%shapes%w = pack(field%shapes%w, mask=field%shapes%w>0.0)
                field%shapes%n = size(field%shapes%x)
            end if
        else if(data_type == standard_scalar_foreground_data) then
            !lenses
            field%lenses%n = n
            allocate(field%lenses%x(n), field%lenses%y(n), field%lenses%z(n), field%lenses%val(n), field%lenses%w(n))

            call read_fits_column_real(fits_iounit, "x", field%lenses%x)
            call read_fits_column_real(fits_iounit, "y", field%lenses%y)
            call read_fits_column_real(fits_iounit, "z", field%lenses%z)
            call read_fits_column_real(fits_iounit, "val", field%lenses%val)
            call read_fits_column_real(fits_iounit, "w", field%lenses%w)
            call assert_fits_status(fits_iounit, status)
            if(apply_mask) then
                field%lenses%x = pack(field%lenses%x, mask=field%lenses%w>0.0)
                field%lenses%y = pack(field%lenses%y, mask=field%lenses%w>0.0)
                field%lenses%z = pack(field%lenses%z, mask=field%lenses%w>0.0)
                field%lenses%val = pack(field%lenses%val, mask=field%lenses%w>0.0)
                field%lenses%w = pack(field%lenses%w, mask=field%lenses%w>0.0)

                field%lenses%x = pack(field%lenses%x, mask=field%lenses%z>=0.0)
                field%lenses%y = pack(field%lenses%y, mask=field%lenses%z>=0.0)
                field%lenses%z = pack(field%lenses%z, mask=field%lenses%z>=0.0)
                field%lenses%val = pack(field%lenses%val, mask=field%lenses%z>=0.0)
                field%lenses%w = pack(field%lenses%w, mask=field%lenses%z>=0.0)
                field%lenses%n = size(field%lenses%x)
            end if
        end if
        !Close fits file
        call ftclos(fits_iounit, status)
        !Deallocate io unit
        call ftfiou(fits_iounit, status)
    end subroutine load_lensing_fits_catalog

    subroutine read_fits_column_real(fits_iounit, colname, array)
        !Arguments
        integer, intent(in) :: fits_iounit
        character(len=*), intent(in) :: colname
        real, dimension(:), intent(out) :: array
        !Variables
        integer :: n, status, colnum, datacode, repeat, width
        real :: nullval
        logical :: anyf
        real(kind=4), dimension(:), allocatable :: single_precision_tmp
        real(kind=8), dimension(:), allocatable :: double_precision_tmp

        status = 0
        n = size(array)
        !Get column number
        call ftgcno(fits_iounit, .false., colname, colnum, status)
        call assert_fits_status(fits_iounit, status)
        !Get column format
        call FTGTCL(fits_iounit, colnum, datacode, repeat, width, status)
        call assert_fits_status(fits_iounit, status)

        !Read column data
        if(datacode == 42 .and. kind(array) == 4) then
            call FTGCVE(fits_iounit, colnum, 1, 1, n, nullval, array, anyf, status)
        else if(datacode == 82 .and. kind(array) == 8) then
            call FTGCVD(fits_iounit, colnum, 1, 1, n, nullval, array, anyf, status)
        else if(datacode == 42 .and. kind(array) == 8) then
            allocate(single_precision_tmp(n))
            call FTGCVE(fits_iounit, colnum, 1, 1, n, nullval, single_precision_tmp, anyf, status)
            array = single_precision_tmp
            deallocate(single_precision_tmp)
        else if(datacode == 82 .and. kind(array) == 4) then
            allocate(double_precision_tmp(n))
            call FTGCVD(fits_iounit, colnum, 1, 1, n, nullval, double_precision_tmp, anyf, status)
            array = double_precision_tmp
            deallocate(double_precision_tmp)
        else
            print *, "Wrong data format: ", datacode
            error stop
        end if
    end subroutine read_fits_column_real

    subroutine load_standard_shape_catalog(filename, field, verbose)
        !Arguments
        character(len=256), intent(in) :: filename
        type(lensing_survey_field), intent(inout) :: field
        logical, optional, intent(in) :: verbose
        !Variables
        integer :: i, file_unit, iostat, n
        logical :: verbose_output

        verbose_output = .false.
        if(present(verbose)) verbose_output = verbose

        n = count_number_of_lines(filename)
        if(verbose_output) print "(A, I7)", "Number of objects in file " // trim(filename) // ": ", n

        allocate(field%shapes%x(n), field%shapes%y(n), field%shapes%z(n), field%shapes%e1(n), field%shapes%e2(n), field%shapes%m(n), field%shapes%w(n))

        if(verbose_output) print "(A)", "Loading file " // trim(filename)
        open(newunit=file_unit, file=filename, iostat=iostat, status="old")
        call assert_iostat(iostat, filename)
        read(unit=file_unit, fmt=*) ((field%shapes%x(i), field%shapes%y(i), field%shapes%z(i), field%shapes%e1(i), field%shapes%e2(i), field%shapes%w(i), field%shapes%m(i)), i=1,n)
        close(unit=file_unit)

        if(verbose_output) print "(A, I7)", "Number of objects with zero weight: ", count(field%shapes%w <= 0.0)
        field%shapes%x = pack(field%shapes%x, mask=field%shapes%w>0.0)
        field%shapes%y = pack(field%shapes%y, mask=field%shapes%w>0.0)
        field%shapes%z = pack(field%shapes%z, mask=field%shapes%w>0.0)
        field%shapes%e1 = pack(field%shapes%e1, mask=field%shapes%w>0.0)
        field%shapes%e2 = pack(field%shapes%e2, mask=field%shapes%w>0.0)
        field%shapes%m = pack(field%shapes%m, mask=field%shapes%w>0.0)
        field%shapes%w = pack(field%shapes%w, mask=field%shapes%w>0.0)
        field%shapes%n = size(field%shapes%x)
    end subroutine load_standard_shape_catalog

    subroutine load_standard_kappa_map(filename, field, verbose)
        !Arguments
        character(len=256), intent(in) :: filename
        type(lensing_survey_field), intent(inout) :: field
        logical, optional, intent(in) :: verbose
        !Variables
        integer :: i, file_unit, iostat, n
        logical :: verbose_output

        verbose_output = .false.
        if(present(verbose)) verbose_output = verbose

        n = count_number_of_lines(filename)
        if(verbose_output) print "(A, I7)", "Number of pixels in file " // trim(filename) // ": ", n

        allocate(field%convergence%x(n), field%convergence%y(n), field%convergence%kappa_E(n), field%convergence%kappa_B(n), field%convergence%w(n))

        if(verbose_output) print "(A)", "Loading file " // trim(filename)
        open(newunit=file_unit, file=filename, iostat=iostat, status="old")
        call assert_iostat(iostat, filename)
        read(unit=file_unit, fmt=*) ((field%convergence%x(i), field%convergence%y(i), field%convergence%kappa_E(i), field%convergence%kappa_B(i), field%convergence%w(i)), i=1,n)
        close(unit=file_unit)

        if(verbose_output) print "(A, I7)", "Number of pixels with zero weight: ", count(field%convergence%w <= 0.0)
        field%convergence%x = pack(field%convergence%x, mask=field%convergence%w>0.0)
        field%convergence%y = pack(field%convergence%y, mask=field%convergence%w>0.0)
        field%convergence%kappa_E = pack(field%convergence%kappa_E, mask=field%convergence%w>0.0)
        field%convergence%kappa_B = pack(field%convergence%kappa_B, mask=field%convergence%w>0.0)
        field%convergence%w = pack(field%convergence%w, mask=field%convergence%w>0.0)
        field%convergence%n = size(field%convergence%x)
    end subroutine load_standard_kappa_map

    subroutine load_standard_lens_catalog(filename, field, apply_masks, verbose)
        !Arguments
        character(len=256), intent(in) :: filename
        type(lensing_survey_field), intent(inout) :: field
        logical, intent(in) :: apply_masks
        logical, optional, intent(in) :: verbose
        !Variables
        integer :: i, file_unit, iostat, n
        logical :: verbose_output

        verbose_output = .false.
        if(present(verbose)) verbose_output = verbose

        n = count_number_of_lines(filename)
        if(verbose_output) print "(A, I7)", "Number of objects in file " // trim(filename) // ": ", n

        allocate(field%lenses%x(n), field%lenses%y(n), field%lenses%z(n), field%lenses%val(n), field%lenses%w(n))

        if(verbose_output) print "(A)", "Loading file " // trim(filename)
        open(newunit=file_unit, file=filename, iostat=iostat, status="old")
        call assert_iostat(iostat, filename)
        read(unit=file_unit, fmt=*) ((field%lenses%x(i), field%lenses%y(i), field%lenses%z(i), field%lenses%val(i), field%lenses%w(i)), i=1,n)
        close(unit=file_unit)

        if(apply_masks) then
            if(verbose_output) print "(A, I7)", "Number of objects with zero weight: ", count(field%lenses%w <= 0.0)
            field%lenses%x = pack(field%lenses%x, mask=field%lenses%w>0.0)
            field%lenses%y = pack(field%lenses%y, mask=field%lenses%w>0.0)
            field%lenses%z = pack(field%lenses%z, mask=field%lenses%w>0.0)
            field%lenses%val = pack(field%lenses%val, mask=field%lenses%w>0.0)
            field%lenses%w = pack(field%lenses%w, mask=field%lenses%w>0.0)

            if(verbose_output) print "(A, I7)", "Number of objects with negative redshift: ", count(field%lenses%z < 0.0)
            field%lenses%x = pack(field%lenses%x, mask=field%lenses%z>=0.0)
            field%lenses%y = pack(field%lenses%y, mask=field%lenses%z>=0.0)
            field%lenses%z = pack(field%lenses%z, mask=field%lenses%z>=0.0)
            field%lenses%val = pack(field%lenses%val, mask=field%lenses%z>=0.0)
            field%lenses%w = pack(field%lenses%w, mask=field%lenses%z>=0.0)
        end if
        field%lenses%n = size(field%lenses%x)
    end subroutine load_standard_lens_catalog

    subroutine load_KiDS_shear_catalog(filename1, filename2, z, n_sample, field, verbose)
        !Arguments
        character(len=256), intent(in) :: filename1, filename2
        real, intent(in) :: z
        integer :: n_sample
        type(lensing_survey_field), intent(inout) :: field
        logical, optional, intent(in) :: verbose
        !Variables
        integer :: i, file_unit1, file_unit2, iostat, n
        real(kind=4), allocatable, dimension(:) :: x, y, e1, e2
        integer, allocatable, dimension(:) :: sample_pos, sort_array
        logical :: verbose_output
        integer, parameter :: n_side = 6000, n_pixel = n_side**2

        verbose_output = .false.
        if(present(verbose)) verbose_output = verbose

        n = n_sample
        if(n_sample <= 0 .or. n_sample > n_pixel) then
            print "(A)", "Invalid sample size. Loading full file."
            n = n_pixel
        end if

        allocate(sample_pos(n), sort_array(n), x(n), y(n), e1(n), e2(n))

        if(n < n_pixel) then
            call random_unique_number_int(array=sample_pos, max=n_pixel, min=1)
            call qsort_arg_int(sample_pos, sort_array)
            sample_pos = sample_pos(sort_array)
        else
            sample_pos = [(i, i=1,n_pixel)]
        end if

        if(verbose_output) print "(A)", "Loading file " // trim(filename1)
        open(newunit=file_unit1, file=filename1, iostat=iostat, status="old", form="unformatted", access="stream")
        call assert_iostat(iostat, filename1)
        do i=1,n
            read(unit=file_unit1, pos=(sample_pos(i)-1)*4 + 1) e1(i)
        end do
        close(unit=file_unit1)
        if(verbose_output) print "(A)", "Loading file " // trim(filename2)
        open(newunit=file_unit2, file=filename2, iostat=iostat, status="old", form="unformatted", access="stream")
        call assert_iostat(iostat, filename2)
        do i=1,n
            read(unit=file_unit2, pos=(sample_pos(i)-1)*4 + 1) e2(i)
        end do
        close(unit=file_unit2)

        where(mod(sample_pos, n_side) /= 0)
            x = mod(sample_pos, n_side)
            y = sample_pos/n_side + 1
        elsewhere
            x = n_side
            y = sample_pos/n_side
        end where

        allocate(field%shapes%x(n), field%shapes%y(n), field%shapes%z(n), field%shapes%e1(n), field%shapes%e2(n), field%shapes%m(n), field%shapes%w(n))

        field%shapes%x = x
        field%shapes%y = y
        field%shapes%z = z
        field%shapes%e1 = e1
        field%shapes%e2 = e2
        field%shapes%m = 0.0
        field%shapes%w = 1.0
        field%shapes%n = n

        deallocate(sample_pos, sort_array, x, y, e1, e2)
    end subroutine load_KiDS_shear_catalog

    subroutine load_KiDS_lens_catalog(filename, z_slice, n_sample, field, verbose)
        !Arguments
        character(len=256), intent(in) :: filename
        real, intent(in) :: z_slice
        integer :: n_sample
        type(lensing_survey_field), intent(inout) :: field
        logical, optional, intent(in) :: verbose
        !Variables
        integer :: i, file_unit, iostat, n, file_size, n_objects
        real(kind=4), allocatable, dimension(:) :: x, y, z, tmp_array
        integer, allocatable, dimension(:) :: sample_pos, sort_array
        logical :: verbose_output
        integer, parameter :: n_variables = 28

        verbose_output = .false.
        if(present(verbose)) verbose_output = verbose

        open(newunit=file_unit, file=filename, iostat=iostat, status="old", form="unformatted", access="stream")
        call assert_iostat(iostat, filename)
        inquire(unit=file_unit, size=file_size)
        n_objects = file_size/(4*n_variables)
        if(mod(file_size, 4*n_variables) /= 0) then
            print "(A, I3, A)", "Cannot load sets of ", n_variables, " variables. Exiting."
            close(unit=file_unit)
            return
        end if

        if(verbose_output) print "(A, I7)", "Number of objects in file " // trim(filename) // ": ", n_objects

        n = n_sample
        if(n_sample <= 0 .or. n_sample > n_objects) then
            print "(A)", "Invalid sample size. Loading full file."
            n = n_objects
        end if

        allocate(sample_pos(n), sort_array(n), x(n), y(n), z(n), tmp_array(n_variables))

        if(n < n_objects) then
            call random_unique_number_int(array=sample_pos, max=n_objects, min=1)
            call qsort_arg_int(sample_pos, sort_array)
            sample_pos = sample_pos(sort_array)
        else
            sample_pos = [(i, i=1,n_objects)]
        end if
        
        if(verbose_output) print "(A)", "Loading file " // trim(filename)
        do i=1,n
            !read(unit=file_unit, pos=(sample_pos(i)-1)*n_variables*4 + 1) tmp_array
            read(unit=file_unit) tmp_array

            x(i) = tmp_array(1)
            y(i) = tmp_array(2)
            z(i) = tmp_array(3)
        end do
        close(unit=file_unit)

        allocate(field%lenses%x(n), field%lenses%y(n), field%lenses%z(n), field%lenses%val(n), field%lenses%w(n))

        field%lenses%x = x
        field%lenses%y = y
        field%lenses%z = z_slice
        field%lenses%val = 1.0
        field%lenses%w = 1.0
        field%lenses%n = n

        deallocate(sample_pos, sort_array, x, y, z)
    end subroutine load_KiDS_lens_catalog

    subroutine load_survey_filelist(survey_filelist, survey, data_type, apply_masks, verbose)
        !Arguments
        character(len=256), dimension(:), intent(in) :: survey_filelist
        type(lensing_survey), intent(inout) :: survey
        integer, intent(in) :: data_type
        logical, optional, intent(in) :: apply_masks, verbose
        !Variables
        integer :: i
        logical :: apply_masks_arg

        apply_masks_arg = .true.
        if(present(apply_masks)) apply_masks_arg=apply_masks

        if(allocated(survey%fields) .and. size(survey%fields) /= survey%n_field) then
            print *, "Number of fields does not match. Aborting."
            stop
        else if(.not. allocated(survey%fields)) then
            allocate(survey%fields(survey%n_field))
        end if

        do i=1, survey%n_field
            if(data_type == standard_scalar_foreground_data) then
                call load_standard_lens_catalog(survey_filelist(i), survey%fields(i), apply_masks=apply_masks_arg, verbose=verbose)
            else if(data_type == standard_shape_data) then
                call load_standard_shape_catalog(survey_filelist(i), survey%fields(i), verbose)
            else if(data_type == standard_kappa_data) then
                call load_standard_kappa_map(survey_filelist(i), survey%fields(i), verbose)
            end if
        end do
    end subroutine load_survey_filelist

    subroutine load_survey_fits_filelist(survey_filelist, survey, apply_masks, verbose)
        !Arguments
        character(len=256), dimension(:), intent(in) :: survey_filelist
        type(lensing_survey), intent(inout) :: survey
        logical, optional, intent(in) :: apply_masks, verbose
        !Variables
        integer :: i
        logical :: apply_masks_arg

        apply_masks_arg = .true.
        if(present(apply_masks)) apply_masks_arg=apply_masks

        if(allocated(survey%fields) .and. size(survey%fields) /= survey%n_field) then
            print *, "Number of fields does not match. Aborting."
            stop
        else if(.not. allocated(survey%fields)) then
            allocate(survey%fields(survey%n_field))
        end if

        do i=1, survey%n_field
            call load_lensing_fits_catalog(survey_filelist(i), survey%fields(i), apply_mask=apply_masks_arg, verbose=verbose)
        end do
    end subroutine load_survey_fits_filelist
    
    subroutine save_marks(filename, survey, results)
        !Arguments
        character(len=256), intent(in) :: filename
        type(lensing_survey), intent(in) :: survey
        type(shear_2pcf_results), dimension(:), intent(in) :: results
        !Variables
        integer(kind=8), parameter :: header_size = 1024
        integer(kind=8) :: n_objects, n_bin
        integer :: i, file_unit, iostat, write_pos

        open(newunit=file_unit, file=filename, iostat=iostat, form="unformatted", access="stream", status="replace")
        call assert_iostat(iostat, filename)
        write(unit=file_unit) header_size
        write(unit=file_unit, pos=header_size-8+1) header_size
        inquire(unit=file_unit, size=write_pos)
        if(write_pos /= header_size) then
            print *, "Failed to create header. Aborting."
            stop
        end if

        write(unit=file_unit, pos=9) survey%n_field

        do i=1,survey%n_field
            n_objects = survey%fields(i)%lenses%n
            n_bin = size(results(i)%xi_E_marks, dim=1)

            inquire(unit=file_unit, size=write_pos)
            write_pos = write_pos + 1

            write(unit=file_unit, pos=write_pos) n_objects, n_bin
            write(unit=file_unit) survey%fields(i)%lenses%x, survey%fields(i)%lenses%y, results(i)%marks_weight
            write(unit=file_unit) results(i)%xi_E_marks, results(i)%xi_B_marks, results(i)%K_marks
        end do

        close(unit=file_unit)
    end subroutine save_marks

    subroutine load_survey_marks(filename, survey)
        !Arguments
        character(len=256), intent(in) :: filename
        type(marks_field_collection), intent(out) :: survey
        !Variables
        integer(kind=8) :: n_field, n_objects, n_bin, header_size, read_pos, dummy
        integer :: i, file_unit, iostat

        open(file=filename, newunit=file_unit, iostat=iostat, form="unformatted", access="stream", status="old")
        call assert_iostat(iostat, filename)
        read(unit=file_unit) header_size, n_field
        read(unit=file_unit, pos=header_size-8+1) dummy
        if(dummy /= header_size) then
            print *, "Header is corrupted. Aborting."
            stop
        end if
        survey%n_field = n_field
        print "(A, I4)", "Number of fields: ", n_field
        if(n_field > 20) then
            print *, "Too many fields. Aborting."
            stop
        end if
        allocate(survey%fields(n_field))
        do i=1,n_field
            read(unit=file_unit) n_objects, n_bin
            print "(A, I4, A, I8)", "Number of objects in field ", i, ": ", n_objects
            if(n_objects > 50000000) then
                print *, "Too many objects. Aborting."
                stop
            end if
            survey%fields(i)%n = n_objects
            survey%fields(i)%n_bin = n_bin
            allocate(survey%fields(i)%x(n_objects), survey%fields(i)%y(n_objects), survey%fields(i)%w(n_objects))
            allocate(survey%fields(i)%xi_E_marks(n_bin, n_objects), survey%fields(i)%xi_B_marks(n_bin, n_objects), survey%fields(i)%K_marks(n_bin, n_objects))
            read(unit=file_unit) survey%fields(i)%x, survey%fields(i)%y, survey%fields(i)%w
            read(unit=file_unit) survey%fields(i)%xi_E_marks, survey%fields(i)%xi_B_marks, survey%fields(i)%K_marks
        end do
        close(file_unit)
    end subroutine load_survey_marks

    subroutine create_marks_fits_file(filename, survey, results, bins, foreground_filelist, double_precision)
        !Arguments
        character(len=*), intent(in) :: filename
        type(lensing_survey), intent(in) :: survey
        type(shear_2pcf_results), dimension(:), intent(in) :: results
        type(binning), intent(in) :: bins
        character(len=256), dimension(:), intent(in) :: foreground_filelist
        logical, intent(in) :: double_precision
        integer :: i, status, fits_iounit, exists, blocksize, nelements, n_objects
        character(len=80) :: field_name, marks_type_format, data_type_format
        real :: nullval

        status = 0
        nullval = 0.0

        !Allocate io unit
        call ftgiou(fits_iounit, status)

        call FTEXIST(filename, exists, status)
        if(exists == 1) then
            print *, "Deleting old marks file: ", trim(filename)
            call ftdkopn(fits_iounit, filename, 0, blocksize, status)
            status = 0
            call FTDELT(fits_iounit, status)
            !Allocate new io unit
            call ftgiou(fits_iounit, status)
            call assert_fits_status(fits_iounit, status)
        end if

        !Open fits file
        call ftdkinit(fits_iounit, filename, blocksize, status)
        call assert_fits_status(fits_iounit, status)

        !Set primary header
        call FTPHPS(fits_iounit, 8, 0, 0, status)
        call assert_fits_status(fits_iounit, status)
        !Write survey information to primary header
        call FTPKYS(fits_iounit, "SURVEY", "", "", status)
        call FTPKYJ(fits_iounit, "NFIELD", int(survey%n_field, kind=4), "Number of fields", status)
        call FTPKYJ(fits_iounit, "NBIN", int(bins%n_bin, kind=4), "Number of angular bins", status)
        call FTPKYj(fits_iounit, "BINTYPE", bins%spacing, "Bin spacing. lin = 0, log = 1, ihs = 2", status)
        call FTPKYE(fits_iounit, "THETAMIN", real(bins%x_min, kind=4), 7, "Innermost bin edge", status)
        call FTPKYE(fits_iounit, "THETAMAX", real(bins%x_max, kind=4), 7, "Outermost bin edge", status)
        if(bins%spacing == ihs_binning) call FTPKYE(fits_iounit, "IHSTHETA", real(bins%ihs_theta, kind=4), 7, "IHS binning theta parameter", status)
        call FTPKNE(fits_iounit, "BINCE", 1, bins%n_bin, real(bins%bin_centers, kind=4), 7, "Bin center&&" , status)
        call FTPKNE(fits_iounit, "BINED", 1, bins%n_bin+1, real(bins%bin_edges, kind=4), 7, "Bin edge&&" , status)
        call FTPKNE(fits_iounit, "BINTH", 1, bins%n_bin, real(bins%bin_centers_eff, kind=4), 7, "Mean theta&&" , status)
        call assert_fits_status(fits_iounit, status)

        if(double_precision) then
            write(marks_type_format, fmt="(I, A)") bins%n_bin, "D"
            data_type_format = "D"
        else
            write(marks_type_format, fmt="(I, A)") bins%n_bin, "E"
            data_type_format = "E"
        end if

        do i=1,survey%n_field
            n_objects = survey%fields(i)%lenses%n

            write(field_name, fmt="(A,I3)") "FIELD", i

            !Create data HDU
            call FTIBIN(fits_iounit, n_objects, 1, "x", data_type_format, "", field_name, 0, status)
            call assert_fits_status(fits_iounit, status)
            !Insert column
            call FTICOL(fits_iounit, 2, "y", data_type_format, status)
            call FTICOL(fits_iounit, 3, "w", data_type_format, status)
            call FTICOL(fits_iounit, 4, "xi_E_MARKS", marks_type_format, status)
            call FTICOL(fits_iounit, 5, "xi_B_MARKS", marks_type_format, status)
            call FTICOL(fits_iounit, 6, "K_MARKS", marks_type_format, status)
            call assert_fits_status(fits_iounit, status)
            !Write information to header
            call FTPKLS(fits_iounit, "LENSFILE", foreground_filelist(i), "", status)
            !Include LONGSTR information
            call FTPLSW(fits_iounit, status)
            
            !Write data
            nelements = size(results(i)%xi_E_marks)
            call FTPCLD(fits_iounit, 1, 1, 1, n_objects, real(survey%fields(i)%lenses%x, kind=8), status)
            call FTPCLD(fits_iounit, 2, 1, 1, n_objects, real(survey%fields(i)%lenses%y, kind=8), status)
            if(double_precision) then
                call FTPCLD(fits_iounit, 3, 1, 1, n_objects, real(results(i)%marks_weight, kind=8), status)
                call FTPCLD(fits_iounit, 4, 1, 1, nelements, real(results(i)%xi_E_marks, kind=8), status)
                call FTPCLD(fits_iounit, 5, 1, 1, nelements, real(results(i)%xi_B_marks, kind=8), status)
                call FTPCLD(fits_iounit, 6, 1, 1, nelements, real(results(i)%K_marks, kind=8), status)
            else
                call FTPCLE(fits_iounit, 3, 1, 1, n_objects, real(results(i)%marks_weight, kind=4), status)
                call FTPCLE(fits_iounit, 4, 1, 1, nelements, real(results(i)%xi_E_marks, kind=4), status)
                call FTPCLE(fits_iounit, 5, 1, 1, nelements, real(results(i)%xi_B_marks, kind=4), status)
                call FTPCLE(fits_iounit, 6, 1, 1, nelements, real(results(i)%K_marks, kind=4), status)
            end if
            call assert_fits_status(fits_iounit, status)
        end do

        !Close fits file
        call ftclos(fits_iounit, status)
        !Deallocate io unit
        call ftfiou(fits_iounit, status)
    end subroutine create_marks_fits_file

    subroutine load_marks_fits_file(filename, marks, bins, verbose)
        !Arguments
        character(len=*), intent(in) :: filename
        type(marks_field_collection), intent(out) :: marks
        type(binning), intent(inout) :: bins
        logical, optional, intent(in) :: verbose
        integer :: i, status, fits_iounit, rwmode, hdutype, hdunum, blocksize, nelements, n_objects, n_bin, bin_spacing
        character(len=256) :: survey_name, field_name, comment
        real :: nullval, theta_min, theta_max, ihs_theta
        logical :: verbose_output
        integer(kind=4) :: int32
        real(kind=4) :: float32
        real(kind=8) :: float64
        real(kind=4), allocatable, dimension(:) :: float32_array

        verbose_output = .false.
        if(present(verbose)) verbose_output = verbose
        status = 0
        nullval = 0.0
        ihs_theta = 1.0

        !Allocate io unit
        call ftgiou(fits_iounit, status)

        !Open fits file
        rwmode = 0
        call ftdkopn(fits_iounit, filename, rwmode, blocksize, status)
        call assert_fits_status(fits_iounit, status)

        !Read primary header information
        call ftgkey(fits_iounit, "SURVEY", survey_name, comment, status)

        call ftgkyj(fits_iounit, "NFIELD", int32, comment, status)
        marks%n_field = int32
        call FTTHDU(fits_iounit, hdunum, status)
        call assert_fits_status(fits_iounit, status)
        if(hdunum-1 /= marks%n_field) then
            print *, "Number of HDUs does not match number of fields. Aborting."
            error stop
        end if

        call FTGKYJ(fits_iounit, "NBIN", int32, comment, status)
        n_bin = int32
        call FTGKYj(fits_iounit, "BINTYPE", int32, comment, status)
        bin_spacing = int32
        if(bin_spacing == ihs_binning) then
            call FTGKYE(fits_iounit, "IHSTHETA", float32, comment, status)
            ihs_theta = float32
        end if
        call FTGKYE(fits_iounit, "THETAMIN", float32, comment, status)
        theta_min = float32
        call FTGKYE(fits_iounit, "THETAMAX", float32, comment, status)
        theta_max = float32
        call assert_fits_status(fits_iounit, status)
        if(verbose_output) print "(A, I0.1, A, I1, A, ES8.1, A, ES8.1, A, ES8.1)", "n_bin = ", n_bin, ", bin spacing = ", bin_spacing, ", theta_min = ", theta_min, ", theta_max = ", theta_max, ", ihs_theta = ", ihs_theta
        call create_bins(x_min=theta_min, x_max=theta_max, n_bin=n_bin, p=ihs_theta, spacing=bin_spacing, bins=bins)

        allocate(float32_array(bins%n_bin))
        call FTGKNE(fits_iounit, "BINTH", 1, bins%n_bin, float32_array, nelements, status)
        if(status == 0) then
            allocate(bins%bin_centers_eff(bins%n_bin))
            bins%bin_centers_eff = real(float32_array, kind=kind(bins%bin_centers_eff))
        else
            deallocate(float32_array)
            status = 0
            comment = "x"
            do while(trim(comment) /= "")
                call ftgmsg(comment)
                print "(A)", comment
            end do
        end if

        allocate(marks%fields(marks%n_field))

        do i=1,marks%n_field
            !Select data HDU
            call ftmahd(fits_iounit, i+1, hdutype, status)
            call assert_fits_status(fits_iounit, status)
            !Read data HDU header information
            call ftgkey(fits_iounit, "LENSFILE", field_name, comment, status)
            call FTGNRW(fits_iounit, n_objects, status)
            if(verbose_output) print "(3A, I7)", "Loading lenses from ", trim_fits_strings(field_name), ". Number of objects: ", n_objects
            marks%fields(i)%n = n_objects
            marks%fields(i)%n_bin = bins%n_bin

            allocate(marks%fields(i)%x(n_objects), marks%fields(i)%y(n_objects), marks%fields(i)%w(n_objects))
            allocate(marks%fields(i)%xi_E_marks(bins%n_bin, n_objects), marks%fields(i)%xi_B_marks(bins%n_bin, n_objects), marks%fields(i)%K_marks(bins%n_bin, n_objects))

            call read_fits_column_real(fits_iounit, "x", marks%fields(i)%x)
            call read_fits_column_real(fits_iounit, "y", marks%fields(i)%y)
            call read_fits_column_real(fits_iounit, "w", marks%fields(i)%w)
            call read_fits_marks_column_real(fits_iounit, "xi_E_marks", bins%n_bin, marks%fields(i)%xi_E_marks)
            call read_fits_marks_column_real(fits_iounit, "xi_B_marks", bins%n_bin, marks%fields(i)%xi_B_marks)
            call read_fits_marks_column_real(fits_iounit, "K_marks", bins%n_bin, marks%fields(i)%K_marks)
        end do

        !Close fits file
        call ftclos(fits_iounit, status)
        !Deallocate io unit
        call ftfiou(fits_iounit, status)
    end subroutine load_marks_fits_file

    subroutine read_fits_marks_column_real(fits_iounit, colname, n_bin, array)
        !Arguments
        integer, intent(in) :: fits_iounit
        character(len=*), intent(in) :: colname
        integer, intent(in) :: n_bin
        real, dimension(:,:), intent(out) :: array
        !Variables
        integer :: nelements, status, colnum, datacode, repeat, width
        real :: nullval
        logical :: anyf
        real(kind=4), dimension(:), allocatable :: single_precision_tmp
        real(kind=8), dimension(:), allocatable :: double_precision_tmp

        status = 0
        nelements = size(array)
        !Get column number
        call ftgcno(fits_iounit, .false., colname, colnum, status)
        call assert_fits_status(fits_iounit, status)
        !Get column data format
        call FTGTCL(fits_iounit, colnum, datacode, repeat, width, status)
        call assert_fits_status(fits_iounit, status)
        if(repeat /= n_bin) then
            print *, "Number of bins does not match. Aborting."
            error stop
        end if

        !Read column data
        if(datacode == 42 .and. kind(array) == 4) then
            call FTGCVE(fits_iounit, colnum, 1, 1, nelements, nullval, array, anyf, status)
        else if(datacode == 82 .and. kind(array) == 8) then
            call FTGCVD(fits_iounit, colnum, 1, 1, nelements, nullval, array, anyf, status)
        else if(datacode == 42 .and. kind(array) == 8) then
            allocate(single_precision_tmp(nelements))
            call FTGCVE(fits_iounit, colnum, 1, 1, nelements, nullval, single_precision_tmp, anyf, status)
            array = reshape(single_precision_tmp, shape(array))
            deallocate(single_precision_tmp)
        else if(datacode == 82 .and. kind(array) == 4) then
            allocate(double_precision_tmp(nelements))
            call FTGCVD(fits_iounit, colnum, 1, 1, nelements, nullval, double_precision_tmp, anyf, status)
            array = reshape(double_precision_tmp, shape(array))
            deallocate(double_precision_tmp)
        else
            print *, "Wrong data format: ", datacode
            error stop
        end if
    end subroutine read_fits_marks_column_real

end module file_formats