program test_fits

implicit none

type catalog
    integer :: n
    real, dimension(:), allocatable :: x, y, z, e1, e2, w, m
end type catalog

type marks_field
    integer :: n, n_bin
    real(kind=4), dimension(:,:), allocatable :: g_t_marks
end type marks_field

real(kind=4), allocatable, dimension(:,:) :: array
real(kind=4), allocatable, dimension(:) :: flat_array

type(catalog) :: field
type(marks_field) :: marks_out, marks_in

call readimage("/Users/yooken/Research/tpcf/data/CDE0047_finalkappaE_scale10.fits", array)

flat_array = reshape(array, shape=[size(array)])
print *, shape(array)

print *, array(231, 341)
print *, flat_array(12346)

call read_catalog("/Users/yooken/Research/tpcf/data/W4.fits", field)

marks_out%n = 2
marks_out%n_bin = 5

allocate(marks_out%g_t_marks(marks_out%n_bin, marks_out%n))
marks_out%g_t_marks(:, 1) = [-1, -2, -3, -4, -5]
marks_out%g_t_marks(:, 2) = [6, 7, 8, 9, 10]

print "(10F6.1)", reshape(marks_out%g_t_marks, shape=[size(marks_out%g_t_marks)])

call create_marks_file("/Users/yooken/Research/tpcf/data/test_marks.fits", marks_out, .true.)

call read_marks_file("/Users/yooken/Research/tpcf/data/test_marks.fits", marks_in)

print "(10F6.1)", reshape(marks_in%g_t_marks, shape=[size(marks_in%g_t_marks)])

contains
    subroutine assert_fits_status(status)
        !Argument
        integer, intent(in) :: status
        !Variable
        character(len=80) :: comment
        if (status /= 0) then
            print "(A, I4)", "FITS error status: ", status
            comment = "x"
            do while(trim(comment) /= "")
                call ftgmsg(comment)
                print "(A)", comment
            end do
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
    subroutine read_catalog(filename, field)
        !Arguments
        character(len=*), intent(in) :: filename
        type(catalog), intent(out) :: field
        !Read a FITS image and determine the minimum and maximum pixel value
        integer :: status, fits_iounit, rwmode, hdutype, blocksize, colnum, nrows, datacode, repeat, width
        character(len=80) :: survey_name, field_name, comment
        real :: nullval
        logical :: anyf

        status = 0
        nullval = 0.0

        !Allocate io unit
        call ftgiou(fits_iounit, status)

        !Open fits file
        rwmode = 0
        call ftdkopn(fits_iounit, filename, rwmode, blocksize, status)
        call assert_fits_status(status)


        !Read primary header information
        call ftgkey(fits_iounit, "SURVEY", survey_name, comment, status)
        print *, "Survey name: ", trim(survey_name), ". Comment: ", trim(comment)
        call ftgkey(fits_iounit, "FIELD", field_name, comment, status)
        print *, "Field name: ", trim(field_name), ". Comment: ", trim(comment)

        !Select data HDU
        call ftmahd(fits_iounit, 2, hdutype, status)
        print *, "HDU 2 type:", hdutype

        !Read data HDU header information
        call ftgkyk(fits_iounit, "NAXIS2", field%n, comment, status)
        print *, "Number of objects: ", field%n, ". Comment: ", trim(comment)

        call FTGNRW(fits_iounit, nrows, status)
        if(nrows /= field%n) then
            print *, "Number of rows doesn't match NAXIS2."
            error stop
        end if

        allocate(field%x(field%n))

        !Read data HDU columns
        call ftgcno(fits_iounit, .false., "RA", colnum, status)

        !Get column data format
        call FTGTCL(fits_iounit, colnum, datacode, repeat, width, status)

        !Read column data
        if(datacode == 42) then
            call FTGCVE(fits_iounit, colnum, 1, 1, field%n, nullval, field%x, anyf, status)
        else if(datacode == 82) then
            call FTGCVD(fits_iounit, colnum, 1, 1, field%n, nullval, field%x, anyf, status)
        else
            print *, "Wrong data format."
        end if

        !Close fits file
        call ftclos(fits_iounit, status)
        !Deallocate io unit
        call ftfiou(fits_iounit, status)
    end subroutine read_catalog

subroutine create_marks_file(filename, marks, double_precision)
        !Arguments
        character(len=*), intent(in) :: filename
        type(marks_field), intent(in) :: marks
        logical, intent(in) :: double_precision
        !Read a FITS image and determine the minimum and maximum pixel value
        integer :: status, fits_iounit, rwmode, hdutype, blocksize, colnum, nrows, datacode, nelements, exists
        character(len=80) :: survey_name, field_name, comment, marks_type_format
        real :: nullval
        logical :: anyf
        real(kind=4) :: sp_value
        real(kind=8) :: dp_value
        integer(kind=4) :: sint
        integer(kind=8) :: lint

        status = 0
        nullval = 0.0

        !Allocate io unit
        call ftgiou(fits_iounit, status)

        call FTEXIST(filename, exists, status)
        if(exists == 1) then
            print *, "Deleting ", trim(filename)
            call ftdkopn(fits_iounit, filename, 0, blocksize, status)
            status = 0
            call FTDELT(fits_iounit, status)
            !Allocate new io unit
            call ftgiou(fits_iounit, status)
        end if

        !Open fits file
        rwmode = 1
        call ftdkinit(fits_iounit, trim(filename), blocksize, status)
        call assert_fits_status(status)

        !Set primary header
        call FTPHPS(fits_iounit, 8, 0, 0, status)
        call assert_fits_status(status)
        !Write survey information to primary header
        call FTPKYS(fits_iounit, "SURVEY", "CFHTLenS", "", status)
        call assert_fits_status(status)
        sp_value = 42.0
        dp_value = 16.5
        call FTPKYD(fits_iounit, "DOUBLE", sp_value, 8, "", status)
        call FTPKYE(fits_iounit, "SINGLE", dp_value, 8, "", status)
        sint = 32
        lint = 42
        call FTPKYJ(fits_iounit, "SINT", lint, "", status)
        call FTPKYK(fits_iounit, "LINT", sint, "", status)

        if(double_precision) then
            write(marks_type_format, fmt="(I, A)") size(marks%g_t_marks, dim=1), "D"
        else
            write(marks_type_format, fmt="(I, A)") size(marks%g_t_marks, dim=1), "E"
        end if

        !Create data HDU
        call FTIBIN(fits_iounit, marks%n, 1, "g_t_marks", marks_type_format, "", "", 0, status)

        !Set data HDU header
        !call FTPHBN(unit,nrows,tfields,ttype,tform,tunit,extname,varidat, > status)

        !Insert column
        call FTICOL(fits_iounit, 2, "g_x_marks", marks_type_format, status)

        !Write data
        nelements = size(marks%g_t_marks)
        if(double_precision) then
            call FTPCLD(fits_iounit, 1, 1, 1, nelements, real(marks%g_t_marks, kind=8), status)
        else
            call FTPCLE(fits_iounit, 1, 1, 1, nelements, real(marks%g_t_marks, kind=4), status)
        end if

        !Close fits file
        call ftclos(fits_iounit, status)
        !Deallocate io unit
        call ftfiou(fits_iounit, status)
    end subroutine create_marks_file

    subroutine read_marks_file(filename, marks)
        !Arguments
        character(len=*), intent(in) :: filename
        type(marks_field), intent(out) :: marks
        !Read a FITS image and determine the minimum and maximum pixel value
        integer :: status, fits_iounit, rwmode, hdutype, blocksize, colnum, nrows, datacode, repeat, width, nelements
        character(len=80) :: survey_name, field_name, comment
        real :: nullval
        logical :: anyf
        real(kind=4), dimension(:), allocatable :: single_precision_array
        real(kind=8), dimension(:), allocatable :: double_precision_array
        real(kind=4) :: sp_value
        real(kind=8) :: dp_value
        integer(kind=4) :: sint
        integer(kind=8) :: lint

        status = 0
        nullval = 0.0

        !Allocate io unit
        call ftgiou(fits_iounit, status)

        !Open fits file
        rwmode = 0
        call ftdkopn(fits_iounit, filename, rwmode, blocksize, status)
        if (status > 0) then
            comment = "x"
            do while(trim(comment) /= "")
                call ftgmsg(comment)
                print *, comment
            end do
            error stop
        end if


        !Read primary header information
        call ftgkey(fits_iounit, "SURVEY", survey_name, comment, status)
        print *, "Survey name: ", trim(survey_name), ". Comment: ", trim(comment)
        call ftgkye(fits_iounit, "SINGLE", sp_value, comment, status)
        call ftgkyd(fits_iounit, "DOUBLE", dp_value, comment, status)
        print *, "E to sp: ", sp_value, ". D to dp: ", dp_value
        call ftgkyd(fits_iounit, "DOUBLE", sp_value, comment, status)
        call ftgkye(fits_iounit, "SINGLE", dp_value, comment, status)
        print *, "E to dp: ", dp_value, ". D to sp: ", sp_value

        call ftgkyj(fits_iounit, "SINT", sint, comment, status)
        call ftgkyk(fits_iounit, "LINT", lint, comment, status)
        print *, "J to sint: ", sint, ". K to lint: ", lint
        call ftgkyj(fits_iounit, "SINT", lint, comment, status)
        call ftgkyk(fits_iounit, "LINT", sint, comment, status)
        print *, "J to lint: ", lint, ". K to sint: ", sint

        !Select data HDU
        call ftmahd(fits_iounit, 2, hdutype, status)
        print *, "HDU 2 type (shoud be BINTABLE->2):", hdutype

        !Read data HDU header information
        call ftgkyk(fits_iounit, "NAXIS2", marks%n, comment, status)
        print *, "Number of objects: ", marks%n, ". Comment: ", trim(comment)

        call FTGNRW(fits_iounit, nrows, status)
        if(nrows /= marks%n) then
            print *, "Number of rows doesn't match NAXIS2: ", nrows
            error stop
        end if

        !Read data HDU columns
        call ftgcno(fits_iounit, .false., "g_t_marks", colnum, status)

        !Get column data format
        call FTGTCL(fits_iounit, colnum, datacode, marks%n_bin, width, status)

        allocate(marks%g_t_marks(marks%n_bin, marks%n))
        nelements = size(marks%g_t_marks)

        !Read column data
        if(datacode == 42) then
            allocate(single_precision_array(nelements))
            call FTGCVE(fits_iounit, colnum, 1, 1, nelements, nullval, single_precision_array, anyf, status)
            marks%g_t_marks = reshape(single_precision_array, shape(marks%g_t_marks))
        else if(datacode == 82) then
            allocate(double_precision_array(nelements))
            call FTGCVD(fits_iounit, colnum, 1, 1, nelements, nullval, double_precision_array, anyf, status)
            marks%g_t_marks = reshape(double_precision_array, shape(marks%g_t_marks))
        else
            print *, "Wrong data format: ", datacode
        end if

        !Close fits file
        call ftclos(fits_iounit, status)
        !Deallocate io unit
        call ftfiou(fits_iounit, status)
    end subroutine read_marks_file

    subroutine readimage(filename, data_array)
        !Argument
        character(len=*) :: filename
        real(kind=4), allocatable, dimension(:,:) :: data_array
        !Read a FITS image and determine the minimum and maximum pixel value
        integer :: status, unit, readwrite, blocksize, nfound, group, firstpix, npixels, n_axis
        integer, dimension(2) :: naxes
        character(len=80) :: comment
        real :: nullval
        logical :: anynull

        status = 0

        !Get an unused Logical Unit Number to use to open the FITS file
        call ftgiou(unit, status)

        !open the FITS file previously created by WRITEIMAGE
        readwrite = 0
        call ftopen(unit, filename, readwrite, blocksize, status)
        if (status > 0) then
            comment = "x"
            do while(trim(comment) /= "")
                call ftgmsg(comment)
                print *, comment
            end do
        end if

        call ftgkyj(unit, "NAXIS", n_axis, comment, status)
        if(n_axis /= 2)then
            print *, "Not a 2d image."
            return
        end if

        !determine the size of the image
        call ftgknj(unit, 'NAXIS', 1, 2, naxes, nfound, status)

        !check that it found both NAXIS1 and NAXIS2 keywords
        if(nfound /= 2)then
            print *, 'READIMAGE failed to read the NAXISn keywords.'
            print *, "nfound: ", nfound
            return
        end if

        allocate(data_array(naxes(1), naxes(2)))
        !initialize variables
        npixels = naxes(1)*naxes(2)
        group = 1
        firstpix = 1
        nullval = -999

        call ftgpve(unit, group, firstpix, npixels, nullval, data_array, anynull, status)

        !close the file and free the unit number
        call ftclos(unit, status)
        call ftfiou(unit, status)

        !check for any error, and if so print out error messages
        if (status > 0) print *, "Some error:", status
    end subroutine readimage
end program test_fits