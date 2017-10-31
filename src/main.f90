program tangential_shear
use utils
use tpcf
use file_formats
implicit none

integer, parameter :: simple_foreground_bootstrap = 1, simple_background_bootstrap = 2,block_bootstrap = 3, marked_simple_bootstrap = 4, marked_block_bootstrap = 5, marked_equal_weight_block_bootstrap = 6, randomized_shear_estimator = 7, survey_block_bootstrap = 8, survey_point_bootstrap = 9
integer, parameter :: standard_field_input_format = 1, KiDS_mocks_slice_input_format = 2, precomputed_marks_input_format = 3, standard_survey_input_format = 4
integer, parameter :: tangential_shear_mode = 1, kappa_mode = 2

integer :: iostat, file_unit, i, j, k, n_args, input_format, tpcf_mode, n_field, n_map

character(len=256) :: cmd_option, cmd_argument

character(len=256) :: foreground_map_filename_prefix, map_output_filename_prefix, map_output_filename, foreground_filename, background_filename, output_filename, lens_catalog_filename, gamma1_filename, gamma2_filename

type(lensing_survey) :: survey

type(binning) :: bins
integer :: bin_spacing, n_bin, n_threads, n_lens_samples, n_shape_samples
real :: theta_min, theta_max, ihs_theta, z_lens_slice, z_shape_slice
logical :: use_kd_tree, verbose, left_handed_coord_system, spherical_coords, output_marks

type(shear_2pcf_results), allocatable, dimension(:) :: results
type(shear_2pcf_results) :: stacked_result

!Bootstrap
integer :: n_A_old, n_bootstrap_samples, masked_pixels, n_x_blocks, n_y_blocks, n_blocks, bootstrap_type
real :: supersampling_factor
type(shear_2pcf_results), allocatable, dimension(:) :: bootstrap_results
real, allocatable, dimension(:,:) :: g_t_covariance, g_x_covariance, g_t_bootstrap, g_x_bootstrap
logical(kind=1), allocatable, dimension(:) :: weight_mask
logical(kind=1), allocatable, dimension(:,:) :: marks_weight_mask

character(len=30) :: bootstrap_output_fmt
character(len=256) :: bootstrap_output_prefix, bootstrap_output_filename

character(len=256), allocatable, dimension(:) :: map_field_filelist, foreground_filelist, background_filelist, output_filelist, covariance_filelist, bootstrap_output_prefix_list

!Precomputed marks
character(len=256) :: precomputed_marks_filename, covariance_filename, block_output_prefix, output_marks_filename
!real, allocatable, dimension(:) :: K_sum
type(marks_field_collection) :: raw_marks_survey, marks_survey
integer :: clustering_scheme

input_format = -1
tpcf_mode = tangential_shear_mode

n_field = 0
n_map = 0

foreground_map_filename_prefix = ""
map_output_filename_prefix = ""
map_output_filename = ""

survey%n_field = -1

foreground_filename = ""
background_filename = ""
output_filename = ""
covariance_filename = ""
lens_catalog_filename = ""
gamma1_filename = ""
gamma1_filename = ""


bootstrap_output_prefix = ""
bootstrap_output_filename = ""

bin_spacing = lin_binning
n_bin = -1
theta_min = -1
theta_max = -1
ihs_theta = 1.0
spherical_coords = .false.
left_handed_coord_system = .false.
n_threads = 1
use_kd_tree = .true.
verbose = .false.

n_lens_samples = -1
n_shape_samples = -1
z_lens_slice = -1.0
z_shape_slice = -1.0
output_marks = .false.

bootstrap_type = 0

n_bootstrap_samples = -1

supersampling_factor = 1.0


n_blocks = -1
n_x_blocks = -1
n_y_blocks = -1

block_output_prefix = ""

clustering_scheme = rectangular_blocks


n_args = iargc()
i = 1
do while(i <= n_args)
    call getarg(i, cmd_option)
    select case(cmd_option)
    case("--help", "-h")
        print *, "There is no help, you are on your own."
        stop
    case("--mode")
        i = i + 1
        call getarg(i, cmd_argument)
        select case(cmd_argument)
            case("tangential-shear")
                tpcf_mode = tangential_shear_mode
                print *, "Calculating tangential shear correlation function."
            case("kappa")
                tpcf_mode = kappa_mode
                print *, "Calculating convergence correlation function."
            case default
                print *, "Invalid mode:", cmd_argument
                stop
        end select
    case("--input-format")
        i = i + 1
        call getarg(i, cmd_argument)
        select case(cmd_argument)
            case("standard-field")
                input_format = standard_field_input_format
                print *, "Using standard format for single field."
            case("standard-survey")
                input_format = standard_survey_input_format
                print *, "Using standard format for a survey."
            case("KiDS-mocks-slice")
                input_format = KiDS_mocks_slice_input_format
                print *, "Using KiDS mocks format for single slice."
            case("precomputed-marks")
                input_format = precomputed_marks_input_format
                print *, "Using precomputed marks."
            case default
                print *, "Invalid input format:", cmd_argument
                stop
        end select
    case("--n-map")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_map
        print *, "Number of (foreground) maps:", n_map
    case("--foreground-map-prefix")
        i = i + 1
        call getarg(i, cmd_argument)
        foreground_map_filename_prefix = cmd_argument
        print *, "(Foreground) map file prefix: ", trim(foreground_map_filename_prefix)
    case("--map-output-prefix")
        i = i + 1
        call getarg(i, cmd_argument)
        map_output_filename_prefix = cmd_argument
        print *, "Map output file prefix: ", trim(map_output_filename_prefix)
    case("--n-field")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_field
        print *, "Number of fields:", n_field
    case("--foreground-file")
        i = i + 1
        call getarg(i, cmd_argument)
        foreground_filename = cmd_argument
        print *, "Foreground file: ", trim(foreground_filename)
    case("--foreground-filelist")
        if(n_field <= 0) then
            print *, "Number of fields not specified. --n-field needs to be passed before any file list."
            stop
        end if
        allocate(foreground_filelist(n_field))
        do j=1, n_field
            i = i + 1
            call getarg(i, cmd_argument)
            foreground_filelist(j) = cmd_argument
        end do
    case("--background-file")
        i = i + 1
        call getarg(i, cmd_argument)
        background_filename = cmd_argument
        print *, "Background file: ", trim(background_filename)
    case("--background-filelist")
        if(n_field <= 0) then
            print *, "Number of fields not specified. --n-field needs to be passed before any file list."
            stop
        end if
        allocate(background_filelist(n_field))
        do j=1, n_field
            i = i + 1
            call getarg(i, cmd_argument)
            background_filelist(j) = cmd_argument
        end do
    case("--lens-catalog-file")
        i = i + 1
        call getarg(i, cmd_argument)
        lens_catalog_filename = cmd_argument
        print *, "Lens catalog file: ", trim(lens_catalog_filename)
    case("--shear1-file")
        i = i + 1
        call getarg(i, cmd_argument)
        gamma1_filename = cmd_argument
        print *, "Shear1 field file: ", trim(gamma1_filename)
    case("--shear2-file")
        i = i + 1
        call getarg(i, cmd_argument)
        gamma2_filename = cmd_argument
        print *, "Shear1 field file: ", trim(gamma2_filename)
    case("--output-file")
        i = i + 1
        call getarg(i, cmd_argument)
        output_filename = cmd_argument
        print *, "Output file: ", trim(output_filename)
    case("--output-filelist")
        if(n_field <= 0) then
            print *, "Number of fields not specified. --n-field needs to be passed before any file list."
            stop
        end if
        allocate(output_filelist(n_field))
        do j=1, n_field
            i = i + 1
            call getarg(i, cmd_argument)
            output_filelist(j) = cmd_argument
        end do
    case("--output-marks-file")
        i = i + 1
        call getarg(i, cmd_argument)
        output_marks_filename = cmd_argument
        print *, "Output marks file: ", trim(output_marks_filename)
        output_marks = .true.
    case("--precomputed-marks-file")
        i = i + 1
        call getarg(i, cmd_argument)
        precomputed_marks_filename = cmd_argument
        print *, "Precomputed marks file: ", trim(precomputed_marks_filename)
    case("--covariance-file")
        i = i + 1
        call getarg(i, cmd_argument)
        covariance_filename = cmd_argument
        print *, "Covariance file: ", trim(covariance_filename)
    case("--covariance-filelist")
        if(n_field <= 0) then
            print *, "Number of fields not specified. --n-field needs to be passed before any file list."
            stop
        end if
        allocate(covariance_filelist(n_field))
        do j=1, n_field
            i = i + 1
            call getarg(i, cmd_argument)
            covariance_filelist(j) = cmd_argument
        end do
    case("--block-output-prefix")
        i = i + 1
        call getarg(i, cmd_argument)
        block_output_prefix = cmd_argument
        print *, "Block output filenames: ", trim(block_output_prefix) // "block_(N).dat"
    case("--n-lens-samples")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_lens_samples
        print *, "Number of lens samples: ", n_lens_samples
    case("--n-shape-samples")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_shape_samples
        print *, "Number of shape samples: ", n_shape_samples
    case("--z-lens-slice")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(F)") z_lens_slice
        print *, "Redshift of lens slice: ", z_lens_slice
    case("--z-shape-slice")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(F)") z_shape_slice
        print *, "Redshift of shape slice: ", z_shape_slice
    case("--bootstrap-type")
        i = i + 1
        call getarg(i, cmd_argument)
        select case(cmd_argument)
            case("simple-foreground")
                bootstrap_type = simple_foreground_bootstrap
                print *, "Using simple foreground bootstrap."
            case("simple-background")
                bootstrap_type = simple_background_bootstrap
                print *, "Using simple background bootstrap."
            case("block")
                bootstrap_type = block_bootstrap
                print *, "Using fixed block bootstrap."
            case("block-marks")
                bootstrap_type = marked_block_bootstrap
                print *, "Using marked point block bootstrap."
            case("simple-marks")
                bootstrap_type = marked_simple_bootstrap
                print *, "Using simple marked point bootstrap."
            case("random-shear")
                bootstrap_type = randomized_shear_estimator
                print *, "Using randomized shear error estimate."
            case("equal-weight-blocks")
                bootstrap_type = marked_equal_weight_block_bootstrap
                print *, "Using pair weighted block bootstrap."
            case("field-rectangular-blocks")
                bootstrap_type = survey_block_bootstrap
                clustering_scheme = rectangular_blocks
                print *, "Using multi field rectangular block bootstrap."
            case("field-voronoi-blocks")
                bootstrap_type = survey_block_bootstrap
                clustering_scheme = voronoi_blocks
                print *, "Using multi field voronoi block bootstrap."
            case("field-equal-weight-blocks")
                bootstrap_type = survey_block_bootstrap
                clustering_scheme = equal_weight_blocks
                print *, "Using multi field equal weight block bootstrap."
            case("field-points")
                bootstrap_type = survey_point_bootstrap
                print *, "Using multi field marked point bootstrap."
            case default
                print *, "Invalid bootstrap type:", cmd_argument
                error stop
        end select
    case("--bootstrap-n-blocks")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_blocks
        print *, "Number of total blocks: ", n_blocks
    case("--bootstrap-n-x-blocks")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_x_blocks
        print *, "Number of blocks in x direction: ", n_x_blocks
    case("--bootstrap-n-y-blocks")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_y_blocks
        print *, "Number of blocks in y direction: ", n_y_blocks
    case("--bootstrap-supersampling-factor")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(F)") supersampling_factor
        print "(A, F4.1)", "Supersampling factor: ", supersampling_factor
    case("--bootstrap-n-samples")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_bootstrap_samples
        print *, "Number bootstrap resamplings: ", n_bootstrap_samples
    case("--bootstrap-output-prefix")
        i = i + 1
        call getarg(i, cmd_argument)
        bootstrap_output_prefix = cmd_argument
        print *, "Bootstrap output filenames: ", trim(bootstrap_output_prefix) // "bootstrap_(N).dat"
    case("--n-bin")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_bin
    case("--theta-min")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(F)") theta_min
    case("--theta-max")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(F)") theta_max
    case("--log-spaced")
        i = i + 1
        call getarg(i, cmd_argument)
        select case(cmd_argument)
            case("true", "True", "1")
                bin_spacing = log_binning
            case("false", "False", "0")
                bin_spacing = lin_binning
            case default
                print *, "Invalid log-spaced value:", cmd_argument
                stop
        end select
    case("--bin-spacing")
        i = i + 1
        call getarg(i, cmd_argument)
        select case(cmd_argument)
            case("lin")
                bin_spacing = lin_binning
                print *, "Using linear binning."
            case("log")
                bin_spacing = log_binning
                print *, "Using logarithmic binning."
            case("ihs")
                bin_spacing = ihs_binning
                print *, "Using inverse hyperbolic sine binning."
            case("sqrt")
                bin_spacing = sqrt_binning
                print *, "Using quadratic binning."
            case default
                print *, "Invalid bin spacing:", cmd_argument
                stop
        end select
    case("--ihs-theta")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(F)") ihs_theta
    case("--n-threads")
        i = i + 1
        call getarg(i, cmd_argument)
        read(cmd_argument, fmt="(I)") n_threads
        print *, "Using ", n_threads, " threads."
    case("--left-handed-coordinates")
        left_handed_coord_system = .true.
        print *, "Using a left-handed coordinate system."
    case("--spherical-coordinates")
        spherical_coords = .true.
        print *, "Using spherical coordinates."
    case("--no-kd-tree")
        use_kd_tree = .false.
    case("--verbose")
        verbose = .true.
    case default
        print *, "Unrecognized option ", cmd_option
        print *, "Exiting."
        stop
    end select
    i = i + 1
end do

call validate_options()

!Convert to radians if using spherical coordinates
if(spherical_coords) then
    theta_min = deg2rad(theta_min)
    theta_max = deg2rad(theta_max)
end if

call create_bins(x_min=theta_min, x_max=theta_max, n_bin=n_bin, p=ihs_theta, spacing=bin_spacing, bins=bins)

survey%n_field = n_field

!====================================================
!========Precomputed marks===========================
!====================================================
if(input_format == precomputed_marks_input_format) then
    call random_seed()
    if(verbose) print *, "Loading precomputed marks."
    !call load_survey_marks(precomputed_marks_filename, marks_survey)
    call load_marks_fits_file(precomputed_marks_filename, raw_marks_survey, bins, verbose=.true.)

    if(n_map > 0) then
        if(verbose) print "(A, I4, A)", "Running on ", n_map, " maps."
        allocate(map_field_filelist(n_field))
        if(foreground_map_filename_prefix(len_trim(foreground_map_filename_prefix):len_trim(foreground_map_filename_prefix)) /= "/") then
            foreground_map_filename_prefix(len_trim(foreground_map_filename_prefix)+1:len_trim(foreground_map_filename_prefix)+1) = "/"
        end if
        do i=1, n_map
            do j=1, n_field
                write(map_field_filelist(j), fmt="(A, I0.1, A)") trim(foreground_map_filename_prefix), i-1, "/" // trim(foreground_filelist(j))
            end do

            survey%n_field = n_field
            call load_survey_fits_filelist(map_field_filelist, survey, apply_masks=.false., verbose=.true.)
            call adjust_survey_to_coord_system(survey, bins, spherical_coords, left_handed_coord_system, verbose)
            call calculate_tpcf_from_marks(raw_marks_survey, survey, bins, spherical_coords, left_handed_coord_system, verbose, marks_survey, results, stacked_result)

            write(map_output_filename, fmt="(A, I0.1, A)") trim(map_output_filename_prefix), i-1, "/" // trim(output_filename)

            if(verbose) print *, "Writing correlation function to output file: ", trim(map_output_filename)
            open(newunit=file_unit, file=map_output_filename, iostat=iostat, status="replace")
            call assert_iostat(iostat, map_output_filename)
            write(file_unit, fmt="(es, es, es, es, I, es)", iostat=iostat) ((stacked_result%theta(k), stacked_result%g_t(k), stacked_result%g_x(k), stacked_result%weights(k), stacked_result%n_pairs(k), stacked_result%mean_theta(k)), k=1,bins%n_bin)
            close(file_unit)

            if(allocated(output_filelist)) then
                do j=1,n_field
                    write(map_output_filename, fmt="(A, I0.1, A)") trim(map_output_filename_prefix), i-1, "/" // trim(output_filelist(j))
                    if(verbose) print *, "Writing correlation function to output file: ", trim(map_output_filename)
                    open(newunit=file_unit, file=map_output_filename, iostat=iostat, status="replace")
                    call assert_iostat(iostat, map_output_filename)
                    write(unit=file_unit, fmt="(es, es, es, es, I, es)", iostat=iostat) ((results(j)%theta(k), results(j)%g_t(k), results(j)%g_x(k), results(j)%weights(k), results(j)%n_pairs(k), results(j)%mean_theta(k)), k=1,bins%n_bin)
                    close(unit=file_unit)
                end do
            end if

            call deallocate_lensing_survey(survey)
        end do
    else
        !Single map only
        if(allocated(foreground_filelist)) then
            !load foreground data
            !call load_survey_filelist(foreground_filelist, survey, data_type=standard_scalar_foreground_data, apply_masks=.false., verbose=.true.)
            call load_survey_fits_filelist(foreground_filelist, survey, apply_masks=.false., verbose=.true.)
            call adjust_survey_to_coord_system(survey, bins, spherical_coords, left_handed_coord_system, verbose)
            call calculate_tpcf_from_marks(raw_marks_survey, survey, bins, spherical_coords, left_handed_coord_system, verbose, marks_survey, results, stacked_result)
        end if

        if(verbose) print *, "Writing correlation function to output file: ", trim(output_filename)
        open(newunit=file_unit, file=output_filename, iostat=iostat, status="replace")
        call assert_iostat(iostat, output_filename)
        write(file_unit, fmt="(es, es, es, es, I, es)", iostat=iostat) ((stacked_result%theta(i), stacked_result%g_t(i), stacked_result%g_x(i), stacked_result%weights(i), stacked_result%n_pairs(i), stacked_result%mean_theta(i)), i=1,bins%n_bin)
        close(file_unit)
        if(allocated(output_filelist)) then
            do i=1,marks_survey%n_field
                if(verbose) print *, "Writing correlation function to output file: ", trim(output_filelist(i))
                open(newunit=file_unit, file=output_filelist(i), iostat=iostat, status="replace")
                call assert_iostat(iostat, output_filelist(i))
                write(unit=file_unit, fmt="(es, es, es, es, I, es)", iostat=iostat) ((results(i)%theta(j), results(i)%g_t(j), results(i)%g_x(j), results(i)%weights(j), results(i)%n_pairs(j), results(i)%mean_theta(j)), j=1,bins%n_bin)
                close(unit=file_unit)
            end do
        end if

        if(n_bootstrap_samples > 0) then
            allocate(bootstrap_results(n_bootstrap_samples))

            if(bootstrap_type == survey_block_bootstrap) then
                call survey_marked_point_block_bootstrap(marks_survey, bins, &
                                                         n_bootstrap_samples, n_blocks, n_x_blocks, n_y_blocks, &
                                                         clustering_scheme, supersampling_factor, &
                                                         spherical_coords, n_threads, &
                                                         bootstrap_results, block_output_prefix)
            else if(bootstrap_type == survey_point_bootstrap) then
                call survey_marked_point_bootstrap(marks_survey, bins, &
                                                   n_bootstrap_samples, supersampling_factor, &
                                                   spherical_coords, n_threads, &
                                                   bootstrap_results)
            end if

            allocate(g_t_covariance(bins%n_bin, bins%n_bin), g_x_covariance(bins%n_bin, bins%n_bin), g_t_bootstrap(n_bootstrap_samples, bins%n_bin), g_x_bootstrap(n_bootstrap_samples, bins%n_bin))

            do i=1,n_bootstrap_samples
                bootstrap_results(i)%theta = bins%bin_centers
                bootstrap_results(i)%mean_theta = 0.0
                if(allocated(bins%bin_centers_eff)) then
                    bootstrap_results(i)%mean_theta = bins%bin_centers_eff
                end if
                g_t_bootstrap(i,:) = bootstrap_results(i)%g_t(:)
                g_x_bootstrap(i,:) = bootstrap_results(i)%g_x(:)
            end do

            call covariance_matrix(g_t_bootstrap, g_t_covariance)
            call covariance_matrix(g_x_bootstrap, g_x_covariance)

            if(bootstrap_output_prefix /= "") then
                print *, "Writing bootstrap samples to file."
                do i=1,n_bootstrap_samples
                    write(bootstrap_output_filename, fmt="(A,I0,A)") trim(bootstrap_output_prefix) // "bootstrap_", i-1, ".dat"
                    open(newunit=file_unit, file=bootstrap_output_filename, iostat=iostat, status="replace")
                    call assert_iostat(iostat, bootstrap_output_filename)
                    write(file_unit, fmt="(es, es, es, es, I, es)", iostat=iostat) ((bootstrap_results(i)%theta(j), bootstrap_results(i)%g_t(j), bootstrap_results(i)%g_x(j), bootstrap_results(i)%weights(j), bootstrap_results(i)%n_pairs(j), bootstrap_results(i)%mean_theta(j)), j=1,bins%n_bin)
                    close(file_unit)
                end do
            end if

            if(covariance_filename /= "") then
                open(newunit=file_unit, file=covariance_filename , iostat=iostat, status="replace")
                call assert_iostat(iostat, covariance_filename)
                write(bootstrap_output_fmt, fmt="(A, I, A)") "(", bins%n_bin, "es)"
                write(unit=file_unit, fmt=bootstrap_output_fmt, iostat=iostat) ((g_t_covariance(i,:)), i=1,bins%n_bin)
                write(unit=file_unit, fmt=bootstrap_output_fmt, iostat=iostat) ((g_x_covariance(i,:)), i=1,bins%n_bin)
                close(unit=file_unit)
            end if

            deallocate(g_t_covariance, g_x_covariance, g_t_bootstrap, g_x_bootstrap)
            deallocate(bootstrap_results)
        end if
    end if
    stop
end if

!====================================================
!========Tangential shear============================
!====================================================
if(tpcf_mode == tangential_shear_mode) then
    !Load data
    if(input_format == standard_field_input_format) then
        allocate(survey%fields(1), results(1))
        survey%n_field = 1
        call load_standard_lens_catalog(foreground_filename, survey%fields(1), apply_masks=.true., verbose=.true.)
        call load_standard_shape_catalog(background_filename, survey%fields(1), verbose=.true.)
    else if(input_format == KiDS_mocks_slice_input_format) then
        allocate(survey%fields(1))
        survey%n_field = 1
        call load_KiDS_lens_catalog(lens_catalog_filename, z_lens_slice, n_lens_samples, survey%fields(1), verbose=.true.)
        call load_KiDS_shear_catalog(gamma1_filename, gamma2_filename, z_shape_slice, n_shape_samples, survey%fields(1), verbose=.true.)

        print "(A, F7.1, A, F7.1)", "Lens x range:", minval(survey%fields(1)%lenses%x), " - ", maxval(survey%fields(1)%lenses%x)
        print "(A, F7.1, A, F7.1)", "Lens y range:", minval(survey%fields(1)%lenses%y), " - ", maxval(survey%fields(1)%lenses%y)
        print "(A, F7.1, A, F7.1)", "Lens z range:", minval(survey%fields(1)%lenses%z), " - ", maxval(survey%fields(1)%lenses%z)
        print "(A, F7.1, A, F7.1)", "Shape x range:", minval(survey%fields(1)%shapes%x), " - ", maxval(survey%fields(1)%shapes%x)
        print "(A, F7.1, A, F7.1)", "Shape y range:", minval(survey%fields(1)%shapes%y), " - ", maxval(survey%fields(1)%shapes%y)

    else if(input_format == standard_survey_input_format) then
        !call load_survey_filelist(foreground_filelist, survey, data_type=standard_scalar_foreground_data, verbose=.true.)
        !call load_survey_filelist(background_filelist, survey, data_type=standard_shape_data, verbose=.true.)
        call load_survey_fits_filelist(foreground_filelist, survey, apply_masks=.true., verbose=.true.)
        call load_survey_fits_filelist(background_filelist, survey, apply_masks=.true., verbose=.true.)
        allocate(results(survey%n_field))
    end if

    do i=1, survey%n_field
        !Convert everything to radians if using spherical coordinates
        if(spherical_coords) then
            survey%fields(i)%lenses%x = deg2rad(survey%fields(i)%lenses%x)
            survey%fields(i)%lenses%y = deg2rad(survey%fields(i)%lenses%y)
            survey%fields(i)%shapes%x = deg2rad(survey%fields(i)%shapes%x)
            survey%fields(i)%shapes%y = deg2rad(survey%fields(i)%shapes%y)

            call rotate_data_from_poles_and_meridian(survey%fields(i)%lenses%x, survey%fields(i)%lenses%y, bins%x_max, verbose, survey%fields(i)%shapes%x, survey%fields(i)%shapes%y)
        end if
        if(left_handed_coord_system) then
            survey%fields(i)%lenses%x = -survey%fields(i)%lenses%x
            survey%fields(i)%shapes%x = -survey%fields(i)%shapes%x
        end if

        print *, "Run tangential shear tpcf on field", i
        call shear_2pcf(survey%fields(i), &
                        bins, &
                        results(i), &
                        spherical_coords, &
                        n_threads, use_kd_tree, verbose)
    end do


!====================================================
!========Convergence=================================
!====================================================
else if(tpcf_mode == kappa_mode) then
!Load data
    if(input_format == standard_field_input_format) then
        allocate(survey%fields(1), results(1))
        survey%n_field = 1
        call load_standard_lens_catalog(foreground_filename, survey%fields(1), apply_masks=.true., verbose=.true.)
        call load_standard_kappa_map(background_filename, survey%fields(1), verbose=.true.)

    else if(input_format == standard_survey_input_format) then
        !call load_survey_filelist(foreground_filelist, survey, data_type=standard_scalar_foreground_data, verbose=.true.)
        !call load_survey_filelist(background_filelist, survey, data_type=standard_kappa_data, verbose=.true.)
        call load_survey_fits_filelist(foreground_filelist, survey, apply_masks=.true., verbose=.true.)
        call load_survey_fits_filelist(background_filelist, survey, apply_masks=.true., verbose=.true.)
        allocate(results(survey%n_field))
    end if

    do i=1, survey%n_field
        !Convert everything to radians if using spherical coordinates
        if(spherical_coords) then
            survey%fields(i)%lenses%x = deg2rad(survey%fields(i)%lenses%x)
            survey%fields(i)%lenses%y = deg2rad(survey%fields(i)%lenses%y)
            survey%fields(i)%convergence%x = deg2rad(survey%fields(i)%convergence%x)
            survey%fields(i)%convergence%y = deg2rad(survey%fields(i)%convergence%y)

            call rotate_data_from_poles_and_meridian(survey%fields(i)%lenses%x, survey%fields(i)%lenses%y, bins%x_max, verbose, survey%fields(i)%convergence%x, survey%fields(i)%convergence%y)
        end if
        if(left_handed_coord_system) then
            survey%fields(i)%lenses%x = -survey%fields(i)%lenses%x
            survey%fields(i)%convergence%x = -survey%fields(i)%convergence%x
        end if

        print *, "Run convergence tpcf on field", i
        call kappa_2pcf(survey%fields(i), &
                        bins, &
                        results(i), &
                        spherical_coords, &
                        n_threads, use_kd_tree, verbose)
    end do
end if

!If output filelist is specified, write individual field results to disk
if(allocated(output_filelist)) then
    do i=1,survey%n_field
        open(newunit=file_unit, file=output_filelist(i), iostat=iostat, status="replace")
        call assert_iostat(iostat, output_filelist(i))
        if(spherical_coords) then
            !Convert back to deggrees
            write(unit=file_unit, fmt="(es, es, es, es, I, es)", iostat=iostat) ((rad2deg(results(i)%theta(j)), results(i)%g_t(j), results(i)%g_x(j), results(i)%weights(j), results(i)%n_pairs(j), rad2deg(results(i)%mean_theta(j))), j=1,bins%n_bin)
        else
            write(unit=file_unit, fmt="(es, es, es, es, I, es)", iostat=iostat) ((results(i)%theta(j), results(i)%g_t(j), results(i)%g_x(j), results(i)%weights(j), results(i)%n_pairs(j), results(i)%mean_theta(j)), j=1,bins%n_bin)
        end if
        close(unit=file_unit)
    end do
end if

call allocate_shear_2pcf_results(stacked_result, bins%n_bin, 0)
stacked_result%theta = bins%bin_centers
do i=1,survey%n_field
    stacked_result%g_t = stacked_result%g_t + results(i)%weights * results(i)%g_t
    stacked_result%g_x = stacked_result%g_x + results(i)%weights * results(i)%g_x
    stacked_result%weights = stacked_result%weights + results(i)%weights
    stacked_result%n_pairs = stacked_result%n_pairs + results(i)%n_pairs
    stacked_result%mean_theta = stacked_result%mean_theta + results(i)%weights * results(i)%mean_theta
end do
stacked_result%g_t = stacked_result%g_t/stacked_result%weights
stacked_result%g_x = stacked_result%g_x/stacked_result%weights
stacked_result%mean_theta = stacked_result%mean_theta/stacked_result%weights
if(spherical_coords) then
    !Convert back to deggrees
    stacked_result%theta = rad2deg(stacked_result%theta)
    stacked_result%mean_theta = rad2deg(stacked_result%mean_theta)
end if
allocate(bins%bin_centers_eff, source=stacked_result%mean_theta)

print *, "Finished running tpcf. Writing to output file. ", trim(output_filename)

open(newunit=file_unit, file=output_filename, iostat=iostat, status="replace")
call assert_iostat(iostat, output_filename)
write(file_unit, fmt="(es, es, es, es, I, es)", iostat=iostat) ((stacked_result%theta(i), stacked_result%g_t(i), stacked_result%g_x(i), stacked_result%weights(i), stacked_result%n_pairs(i), stacked_result%mean_theta(i)), i=1,bins%n_bin)
close(file_unit)

if(output_marks) then
    !call save_marks(output_marks_filename, survey, results)
    call create_marks_fits_file(output_marks_filename, survey, results, bins, foreground_filelist, double_precision=.false.)
end if

do i=1,survey%n_field
    allocate(weight_mask(survey%fields(i)%lenses%n))
    weight_mask = results(i)%marks_weight > 0.0

    !Remove foreground objects that do not contribute to the correlation function (marked weights are zero)
    n_A_old = survey%fields(i)%lenses%n
    survey%fields(i)%lenses%n = count(weight_mask)
    print "(A, I6, A, F6.2, A)", "Masking foreground objects with zero marks weight. n masked: ", n_A_old - survey%fields(i)%lenses%n, ". Fraction remaining: ", real(survey%fields(i)%lenses%n)/n_A_old*100.0, "%."
    survey%fields(i)%lenses%x = pack(survey%fields(i)%lenses%x, mask=weight_mask)
    survey%fields(i)%lenses%y = pack(survey%fields(i)%lenses%y, mask=weight_mask)
    survey%fields(i)%lenses%z = pack(survey%fields(i)%lenses%z, mask=weight_mask)
    survey%fields(i)%lenses%val = pack(survey%fields(i)%lenses%val, mask=weight_mask)
    survey%fields(i)%lenses%w = pack(survey%fields(i)%lenses%w, mask=weight_mask)

    deallocate(weight_mask)
end do

!Bootstrapping
do i=1, survey%n_field
    if(tpcf_mode == tangential_shear_mode) then
        if(n_bootstrap_samples > 0) then
            call random_seed()
            allocate(bootstrap_results(n_bootstrap_samples))
            do j=1,n_bootstrap_samples
                call allocate_shear_2pcf_results(bootstrap_results(j), bins%n_bin, survey%fields(i)%lenses%n)
            end do

            if(bootstrap_type == simple_foreground_bootstrap) then
                call simple_bootstrap_foreground(survey%fields(i), &
                                                 bins, &
                                                 n_bootstrap_samples, &
                                                 bootstrap_results, &
                                                 spherical_coords, n_threads, use_kd_tree, verbose, .false.)
    !        else if(bootstrap_type == block_bootstrap) then
    !            call fixed_block_bootstrap(field%lenses%x, field%lenses%y, field%lenses%z, field%lenses%val, field%lenses%w, &
    !                                       field%shapes%x, field%shapes%y, field%shapes%z, field%shapes%e1, field%shapes%e2, field%shapes%w, field%shapes%m, &
    !                                       bins, &
    !                                       n_x_blocks, n_y_blocks, n_bootstrap_samples, supersampling_factor, &
    !                                       bootstrap_results, &
    !                                       spherical_coords, n_threads, use_kd_tree, verbose, .false.)
    !        else if(bootstrap_type == marked_block_bootstrap) then
    !            call fixed_block_bootstrap(field%lenses%x, field%lenses%y, field%lenses%z, field%lenses%val, field%lenses%w, &
    !                                       field%shapes%x, field%shapes%y, field%shapes%z, field%shapes%e1, field%shapes%e2, field%shapes%w, field%shapes%m, &
    !                                       bins, &
    !                                       n_x_blocks, n_y_blocks, n_bootstrap_samples, supersampling_factor, &
    !                                       bootstrap_results, &
    !                                       spherical_coords, n_threads, use_kd_tree, verbose, &
    !                                       .true., results%xi_E_marks, results%xi_B_marks, results%K_marks, results%marks_weight)
            else if(bootstrap_type == marked_simple_bootstrap) then
                call simple_bootstrap_foreground(survey%fields(i), &
                                                 bins, &
                                                 n_bootstrap_samples, &
                                                 bootstrap_results, &
                                                 spherical_coords, n_threads, use_kd_tree, verbose, &
                                                 .true., results(i)%xi_E_marks, results(i)%xi_B_marks, results(i)%K_marks, results(i)%marks_weight)
    !        else if(bootstrap_type == marked_equal_weight_block_bootstrap) then
    !            call pair_weighted_block_bootstrap(field%lenses%x, field%lenses%y, field%lenses%z, field%lenses%val, field%lenses%w, &
    !                                               field%shapes%x, field%shapes%y, field%shapes%z, field%shapes%e1, field%shapes%e2, field%shapes%w, field%shapes%m, &
    !                                               bins, &
    !                                               n_blocks, n_bootstrap_samples, supersampling_factor, .false., &
    !                                               bootstrap_results, &
    !                                               spherical_coords, n_threads, use_kd_tree, verbose, &
    !                                               results%xi_E_marks, results%xi_B_marks, results%K_marks, results%marks_weight)
    !        else if(bootstrap_type == simple_background_bootstrap) then
    !            call simple_bootstrap_background(field%lenses%x, field%lenses%y, field%lenses%z, field%lenses%val, field%lenses%w, &
    !                                             field%shapes%x, field%shapes%y, field%shapes%z, field%shapes%e1, field%shapes%e2, field%shapes%w, field%shapes%m, &
    !                                             bins, &
    !                                             n_bootstrap_samples, &
    !                                             bootstrap_results, &
    !                                             spherical_coords, n_threads, use_kd_tree, verbose)
    !        else if(bootstrap_type == randomized_shear_estimator) then
    !            call randomize_shear(field%lenses%x, field%lenses%y, field%lenses%z, field%lenses%val, field%lenses%w, &
    !                                 field%shapes%x, field%shapes%y, field%shapes%z, field%shapes%e1, field%shapes%e2, field%shapes%w, field%shapes%m, &
    !                                 bins, &
    !                                 n_bootstrap_samples, &
    !                                 bootstrap_results, &
    !                                 spherical_coords, n_threads, use_kd_tree, verbose)
            else
                print *, "No bootstrap type specified. Exiting."
                stop
            end if
            if(bootstrap_type == marked_simple_bootstrap .or. bootstrap_type == marked_block_bootstrap .or. bootstrap_type == marked_equal_weight_block_bootstrap) then
                do j=1,n_bootstrap_samples
                    bootstrap_results(j)%theta = results(i)%theta
                    bootstrap_results(j)%mean_theta = 0.0
                end do
            end if
            if(allocated(bootstrap_output_prefix_list)) then
                print *, "Writing bootstrap samples to file."
                do j=1, n_bootstrap_samples
                    write(bootstrap_output_filename, fmt="(A,I0,A)") trim(bootstrap_output_prefix_list(i)) // "bootstrap_", j, ".dat"
                    open(newunit=file_unit, file=bootstrap_output_filename, iostat=iostat, status="replace")
                    call assert_iostat(iostat, bootstrap_output_filename)
                    write(file_unit, fmt="(es, es, es, es, I, es)", iostat=iostat) ((bootstrap_results(j)%theta(k), bootstrap_results(j)%g_t(k), bootstrap_results(j)%g_x(k), bootstrap_results(j)%weights(k), bootstrap_results(j)%n_pairs(k), bootstrap_results(j)%mean_theta(k)), k=1,bins%n_bin)
                    close(file_unit)
                end do
            end if
            if(allocated(covariance_filelist)) then
                print *, "Calculating covariance."
                allocate(g_t_covariance(bins%n_bin, bins%n_bin), g_x_covariance(bins%n_bin, bins%n_bin), g_t_bootstrap(n_bootstrap_samples, bins%n_bin), g_x_bootstrap(n_bootstrap_samples, bins%n_bin))

                do j=1,n_bootstrap_samples
                    g_t_bootstrap(j,:) = bootstrap_results(j)%g_t(:)
                    g_x_bootstrap(j,:) = bootstrap_results(j)%g_x(:)
                end do
                call covariance_matrix(g_t_bootstrap, g_t_covariance)
                call covariance_matrix(g_x_bootstrap, g_x_covariance)

                open(newunit=file_unit, file=covariance_filelist(i) , iostat=iostat, status="replace")
                call assert_iostat(iostat, covariance_filelist(i))
                write(bootstrap_output_fmt, fmt="(A, I, A)") "(", bins%n_bin, "es)"
                write(unit=file_unit, fmt=bootstrap_output_fmt, iostat=iostat) ((g_t_covariance(j,:)), j=1,bins%n_bin)
                write(unit=file_unit, fmt=bootstrap_output_fmt, iostat=iostat) ((g_x_covariance(j,:)), j=1,bins%n_bin)
                close(unit=file_unit)
            end if

            deallocate(bootstrap_results)
        end if
    end if
end do

call destroy_bins(bins)

contains
    subroutine validate_options()
        if(input_format < 0) then
            print *, "Input format not specified."
            stop
        else if(tpcf_mode < 0) then
            print *, "Mode not specified."
            stop
        end if
        if(input_format == standard_field_input_format) then
            if(foreground_filename == "") then
                print *, "Foreground filename not provided."
                stop
            else if(background_filename == "") then
                print *, "Background filename not provided."
                stop
            end if
        else if(input_format == KiDS_mocks_slice_input_format) then
            if(lens_catalog_filename == "") then
                print *, "Lens catalog filename not provided."
                stop
            else if(gamma1_filename == "") then
                print *, "Gamm1 filename not provided."
                stop
            else if(gamma2_filename == "") then
                print *, "Gamm2 filename not provided."
                stop
            else if(z_lens_slice < 0.0 .or. z_shape_slice <= 0.0) then
                print *, "Slice redshift not specifed."
                stop
            end if
        else if(input_format == precomputed_marks_input_format) then
            if(precomputed_marks_filename == "") then
                print *, "Precomputed marks filename not provided."
                stop
            end if
        else if(input_format == standard_survey_input_format) then
            if(.not. allocated(foreground_filelist)) then
                print *, "No foreground data provided."
                stop
            else if(.not. allocated(background_filelist)) then
                print *, "No background data provided."
                stop
            end if
        end if

        if(theta_min < 0.0 .or. theta_max < 0.0) then
            print *, "Angular range not specified."
            stop
        end if

        if(n_bin < 0) then
            print *, "Number of bins not specified."
            stop
        end if

        if(output_filename == "") then
            print *, "Output filename not specified."
            stop
        end if

        if(bootstrap_type == simple_foreground_bootstrap) then
            if(n_bootstrap_samples < 0) then
                print *, "Number of bootstrap resamplings not specified."
                stop
            else if(bootstrap_output_prefix == "" .and. covariance_filename == "") then
                print *, "Bootstrap output path not specified."
                stop
            end if
        else if(bootstrap_type == simple_foreground_bootstrap) then
            if(n_bootstrap_samples < 0) then
                print *, "Number of bootstrap resamplings not specified."
                stop
            else if(bootstrap_output_prefix == "" .and. covariance_filename == "") then
                print *, "Bootstrap output path not specified."
                stop
            end if
        else if(bootstrap_type == block_bootstrap) then
            if(n_bootstrap_samples < 0) then
                print *, "Number of bootstrap resamplings not specified."
                stop
            else if(n_x_blocks < 0 .or. n_y_blocks < 0) then
                print *, "Number of blocks not specified."
                stop
            else if(bootstrap_output_prefix == "" .and. covariance_filename == "") then
                print *, "Bootstrap output path not specified."
                stop
            end if
        else if(bootstrap_type == marked_block_bootstrap) then
            if(n_bootstrap_samples < 0) then
                print *, "Number of bootstrap resamplings not specified."
                stop
            else if(n_x_blocks < 0 .or. n_y_blocks < 0) then
                print *, "Number of blocks not specified."
                stop
            else if(bootstrap_output_prefix == "" .and. covariance_filename == "") then
                print *, "Bootstrap output path not specified."
                stop
            end if
        else if(bootstrap_type == marked_simple_bootstrap) then
            if(n_bootstrap_samples < 0) then
                print *, "Number of bootstrap resamplings not specified."
                stop
            else if(bootstrap_output_prefix == "" .and. covariance_filename == "") then
                print *, "Bootstrap output path not specified."
                stop
            end if
        else if(bootstrap_type == marked_equal_weight_block_bootstrap) then
            if(n_bootstrap_samples < 0) then
                print *, "Number of bootstrap resamplings not specified."
                stop
            else if(n_blocks < 0) then
                print *, "Number of blocks not specified."
                stop
            else if(bootstrap_output_prefix == "" .and. covariance_filename == "") then
                print *, "Bootstrap output path not specified."
                stop
            end if
        end if
    end subroutine validate_options

    subroutine rotate_data_from_poles_and_meridian(x_A, y_A, r_max, verbose, x_B, y_B)
        !Arguments
        real, dimension(:), intent(inout) :: x_A, y_A
        real, intent(in) :: r_max
        logical, intent(in) :: verbose
        real, dimension(:), optional, intent(inout) :: x_B, y_B
        !Variables
        real :: x_min, x_max, y_min, y_max, field_center_x, field_center_y, x_median, y_median
        integer :: n_A
        logical :: rotate_field

        rotate_field = .false.
        n_A = size(x_A)

        x_min = minval(x_A)
        x_max = maxval(x_A)
        y_min = minval(y_A)
        y_max = maxval(y_A)
        !Check for poles in search region
        if(y_max + r_max >= pi/2.0) then
            rotate_field = .true.
            if(verbose) print *, "Search region includes north pole."
        else if(y_min - r_max <= -pi/2.0) then
            rotate_field = .true.
            if(verbose) print *, "Search region includes south pole."
        end if
        !Check for prime meridian in x_A
        if(x_max - x_min > pi) then
            if(verbose) print "(A, F7.2, A, F7.2)", "RA max: ", rad2deg(x_max), " min: ", rad2deg(x_min)
            !Field lies on both sides of prime meridian or field is larger than 180 degrees (which I assume it's not).
            !Make x_A values continuous over prime meridian
            where(x_A > pi) x_A = x_A - 2*pi
            x_min = minval(x_A)
            x_max = maxval(x_A)
            rotate_field = .true.
            if(verbose) print *, "Search region includes prime meridian."
            if(verbose) print *, "Making x coordinates continuous over prime meridian."
        else if( distance_to_meridian_sphere(x_max, y_max, 0.0) <= r_max .or. &
            distance_to_meridian_sphere(x_min, y_max, 0.0) <= r_max .or. &
            distance_to_meridian_sphere(x_max, y_min, 0.0) <= r_max .or. &
            distance_to_meridian_sphere(x_min, y_min, 0.0) <= r_max ) then
            rotate_field = .true.
            if(verbose) print *, "Search region includes prime meridian."
        end if
        !Check for prime meridian in x_B
        if(present(x_B)) then
            if(maxval(x_B) - minval(x_B) > pi) then
                !Make x_B values continuous over prime meridian
                where(x_B > pi) x_B = x_B - 2*pi
            end if
        end if
        if(rotate_field) then
            x_median = median(x_A)
            y_median = median(y_A)
            field_center_x = x_median
            field_center_y = y_median
            if(verbose) print "(A, F7.2, A, F7.2)", "Field center RA: ", rad2deg(field_center_x), " Dec: ", rad2deg(field_center_y)
            if(verbose) print "(A, F7.2, A, F7.2)", "Field size RA: ", rad2deg(x_max-x_min), " Dec: ", rad2deg(y_max-y_min)
            if(verbose) print "(A, F7.2, A, F7.2)", "Rotating field by RA: ", rad2deg(pi - field_center_x), " Dec: ", -rad2deg(field_center_y)
            call rotation_sphere(x_A, y_A, pi-field_center_x, -field_center_y, x_A, y_A)
            if(present(x_B) .and. present(y_B)) call rotation_sphere(x_B, y_B, pi-field_center_x, -field_center_y, x_B, y_B)
        end if
    end subroutine rotate_data_from_poles_and_meridian

    subroutine calculate_tpcf_from_marks(marks, foreground_survey, bins, spherical_coords, left_handed_coord_system, verbose, marks_survey, results, stacked_result)
        !Arguments
        type(marks_field_collection), intent(in) :: marks
        type(lensing_survey), intent(in) :: foreground_survey
        type(binning), intent(in) :: bins
        logical, intent(in) :: spherical_coords, left_handed_coord_system, verbose
        type(marks_field_collection), intent(out) :: marks_survey
        type(shear_2pcf_results), allocatable, dimension(:), intent(out) :: results
        type(shear_2pcf_results), intent(out) :: stacked_result
        !Variables
        logical(kind=1), allocatable, dimension(:) :: weight_mask
        logical(kind=1), allocatable, dimension(:,:) :: marks_weight_mask
        real, allocatable, dimension(:) :: K_sum
        integer :: i

        if(marks%n_field /= foreground_survey%n_field) then
            print *, "Number of fields does not match. Aborting."
            stop
        end if
        do i=1, marks%n_field
            if(marks%fields(i)%n /= foreground_survey%fields(i)%lenses%n) then
                print "(A, I3, A)", "Number of lenses do not match in field ", i, ". Aborting."
                stop
            else if(any(abs(marks%fields(i)%x - foreground_survey%fields(i)%lenses%x) > 1.0e-5)) then
                print "(A, I3, A)", "X-coordinates of lenses do not match in field ", i, ". Aborting."
                print *, "Marks x:", marks%fields(i)%x(1:5)
                print *, "loaded x:", foreground_survey%fields(i)%lenses%x(1:5)
                print *, "Total difference: ", sum(abs(marks%fields(i)%x - foreground_survey%fields(i)%lenses%x))
                stop
            else if(any(abs(marks%fields(i)%y - foreground_survey%fields(i)%lenses%y) > 1.0e-5)) then
                print "(A, I3, A)", "Y-coordinates of lenses do not match in field ", i, ". Aborting."
                stop
            end if
        end do

        print *, "Copy marks."
        marks_survey%n_field = marks%n_field
        marks_survey%n_marks = marks%n_marks
        allocate(marks_survey%fields(marks_survey%n_field))
        do i=1,marks_survey%n_field
            marks_survey%fields(i)%n = marks%fields(i)%n
            marks_survey%fields(i)%n_bin = marks%fields(i)%n_bin

            allocate(marks_survey%fields(i)%x, source=marks%fields(i)%x)
            allocate(marks_survey%fields(i)%y, source=marks%fields(i)%y)
            allocate(marks_survey%fields(i)%w, source=marks%fields(i)%w)

            allocate(marks_survey%fields(i)%xi_E_marks, source=marks%fields(i)%xi_E_marks)
            allocate(marks_survey%fields(i)%xi_B_marks, source=marks%fields(i)%xi_B_marks)
            allocate(marks_survey%fields(i)%K_marks, source=marks%fields(i)%K_marks)
        end do

        do i=1, marks_survey%n_field
            do j=1,marks_survey%fields(i)%n
                marks_survey%fields(i)%xi_E_marks(:,j) = marks_survey%fields(i)%xi_E_marks(:,j) * foreground_survey%fields(i)%lenses%val(j)*survey%fields(i)%lenses%w(j)
                marks_survey%fields(i)%xi_B_marks(:,j) = marks_survey%fields(i)%xi_B_marks(:,j) * foreground_survey%fields(i)%lenses%val(j)*survey%fields(i)%lenses%w(j)
                marks_survey%fields(i)%K_marks(:,j) = marks_survey%fields(i)%K_marks(:,j) * foreground_survey%fields(i)%lenses%w(j)
            end do
            marks_survey%fields(i)%w = marks_survey%fields(i)%w * foreground_survey%fields(i)%lenses%w
        end do

        marks_survey%n_marks = 0
        do i=1, marks_survey%n_field
            !Remove masked points
            allocate(weight_mask(marks_survey%fields(i)%n), marks_weight_mask(marks_survey%fields(i)%n, bins%n_bin))
            weight_mask = marks_survey%fields(i)%w > 0.0
            marks_survey%fields(i)%n = count(weight_mask)

            marks_weight_mask = spread(weight_mask, dim=1, ncopies=bins%n_bin)
            marks_survey%fields(i)%xi_E_marks = reshape(pack(marks_survey%fields(i)%xi_E_marks, mask=marks_weight_mask), shape=[bins%n_bin, marks_survey%fields(i)%n])
            marks_survey%fields(i)%xi_B_marks = reshape(pack(marks_survey%fields(i)%xi_B_marks, mask=marks_weight_mask), shape=[bins%n_bin, marks_survey%fields(i)%n])
            marks_survey%fields(i)%K_marks = reshape(pack(marks_survey%fields(i)%K_marks, mask=marks_weight_mask), shape=[bins%n_bin, marks_survey%fields(i)%n])
            !marks_survey%fields(i)%w_marks = reshape(pack(marks_survey%fields(i)%w_marks, mask=marks_weight_mask), shape=[bins%n_bin, marks_survey%fields(i)%n])
            marks_survey%fields(i)%x = pack(marks_survey%fields(i)%x, mask=weight_mask)
            marks_survey%fields(i)%y = pack(marks_survey%fields(i)%y, mask=weight_mask)
            marks_survey%fields(i)%w = pack(marks_survey%fields(i)%w, mask=weight_mask)

            deallocate(weight_mask, marks_weight_mask)
            print "(A, I3, A, I7)", "Field ", i, " remaining points: ", marks_survey%fields(i)%n
            marks_survey%n_marks = marks_survey%n_marks + marks_survey%fields(i)%n
        end do
        print "(A, I7)", "Total remaining points: ", marks_survey%n_marks

        allocate(results(marks_survey%n_field), K_sum(bins%n_bin))
        call allocate_shear_2pcf_results(stacked_result, bins%n_bin, 0)
        stacked_result%theta = bins%bin_centers
        stacked_result%mean_theta = 0.0
        if(allocated(bins%bin_centers_eff)) then
            stacked_result%mean_theta = bins%bin_centers_eff
        end if
        stacked_result%n_pairs = 0.0
        stacked_result%g_t = 0.0
        stacked_result%g_t = 0.0
        K_sum = 0.0
        do i=1,marks_survey%n_field
            call allocate_shear_2pcf_results(results(i), bins%n_bin, 0)
            results(i)%theta = bins%bin_centers
            results(i)%g_t = sum(marks_survey%fields(i)%xi_E_marks, dim=2)/sum(marks_survey%fields(i)%K_marks, dim=2)
            results(i)%g_x = sum(marks_survey%fields(i)%xi_B_marks, dim=2)/sum(marks_survey%fields(i)%K_marks, dim=2)
            results(i)%weights = sum(marks_survey%fields(i)%w)
            results(i)%n_pairs = 0
            results(i)%mean_theta = 0.0
            if(allocated(bins%bin_centers_eff)) then
                results(i)%mean_theta = bins%bin_centers_eff
            end if
            if(spherical_coords) then
                !Convert back to deggrees
                results(i)%theta = rad2deg(results(i)%theta)
                results(i)%mean_theta = rad2deg(results(i)%mean_theta)
            end if

            stacked_result%g_t = stacked_result%g_t + sum(marks_survey%fields(i)%xi_E_marks, dim=2)
            stacked_result%g_x = stacked_result%g_x + sum(marks_survey%fields(i)%xi_B_marks, dim=2)
            K_sum = K_sum + sum(marks_survey%fields(i)%K_marks, dim=2)
            stacked_result%weights = stacked_result%weights + results(i)%weights
        end do

        stacked_result%g_t = stacked_result%g_t/K_sum
        stacked_result%g_x = stacked_result%g_x/K_sum
        if(spherical_coords) stacked_result%theta = rad2deg(stacked_result%theta)
        if(spherical_coords) stacked_result%mean_theta = rad2deg(stacked_result%mean_theta)

    end subroutine calculate_tpcf_from_marks

    subroutine adjust_survey_to_coord_system(survey, bins, spherical_coords, left_handed_coord_system, verbose)
        !Arguments
        type(lensing_survey), intent(inout) :: survey
        type(binning), intent(in) :: bins
        logical, intent(in) :: spherical_coords, left_handed_coord_system, verbose
        !Variables
        integer :: i

        do i=1,survey%n_field
            !Convert everything to radians if using spherical coordinates
            if(spherical_coords) then
                survey%fields(i)%lenses%x = deg2rad(survey%fields(i)%lenses%x)
                survey%fields(i)%lenses%y = deg2rad(survey%fields(i)%lenses%y)
                call rotate_data_from_poles_and_meridian(survey%fields(i)%lenses%x, survey%fields(i)%lenses%y, bins%x_max, verbose)
            end if
            if(left_handed_coord_system) then
                survey%fields(i)%lenses%x = -survey%fields(i)%lenses%x
            end if
        end do
    end subroutine adjust_survey_to_coord_system

end program tangential_shear