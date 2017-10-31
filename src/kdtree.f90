module kdtree_module
    use constants
    use utils
    
    implicit none
    !Floating point precision
    !integer, parameter :: fpp = kind(1.0d0)
    integer, parameter :: fpp = kind(1.0)

    type tree_node
        integer :: split_dim
        integer :: l_idx, u_idx
        real(fpp) :: split_value
        real(fpp), allocatable :: bbox(:,:)
        type(tree_node), pointer :: left_ptr, right_ptr
    end type tree_node

    type kdtree
        integer :: n = 0, n_dim = 0, leafsize = 10, max_nodes = 0
        real(fpp), pointer :: data_ptr(:,:)
        integer, allocatable :: indicies(:)
        type(tree_node), pointer :: root_ptr
        !Debug information
        integer :: n_nodes = 0, n_mid_point_slides = 0, n_children = 0
    end type kdtree

    type node_pointer
        type(tree_node), pointer :: ptr
    end type node_pointer

    type search_record
        integer :: n_included_nodes = 0
        type(node_pointer), allocatable :: included_nodes(:)
    end type search_record

    public fpp
    public tree_node, kdtree, search_record
    public grow_tree, burn_tree, search_within_r, search_within_r_spherical, write_node_data_to_file, write_search_debug_data_to_file
    private

    contains

    subroutine burn_tree(tree)
        implicit none
        !Arguments
        type(kdtree) :: tree

        if(associated(tree%root_ptr)) then
            call destroy_node(tree%root_ptr)
        end if
        deallocate(tree%indicies)
        nullify(tree%data_ptr)

        contains
        recursive subroutine destroy_node(node)
            implicit none
            !Arguments
            type(tree_node), pointer :: node

            if(associated(node%left_ptr)) then
                call destroy_node(node%left_ptr)
            end if
            if(associated(node%right_ptr)) then
                call destroy_node(node%right_ptr)
            end if

            deallocate(node%bbox)
            deallocate(node)
            nullify(node)
        end subroutine destroy_node
    end subroutine burn_tree

    subroutine grow_tree(tree, input_data, leafsize)
        implicit none
        !Arguments
        type(kdtree) :: tree
        real(fpp), target, intent(in) :: input_data(:,:)
        integer, optional, intent(in) :: leafsize
        !Variables
        integer :: i
        real(fpp) :: bbox(size(input_data, 2), 2)

        if(present(leafsize)) then
            tree%leafsize = leafsize
        end if

        tree%n = size(input_data, 1)
        tree%n_dim = size(input_data, 2)
        tree%data_ptr => input_data

        tree%max_nodes = tree%n * log(real(tree%n))

        allocate(tree%indicies(size(input_data, dim=1)))
        tree%indicies = [(i, i=1,tree%n)]

        bbox(:,1) = minval(tree%data_ptr, dim=1)
        bbox(:,2) = maxval(tree%data_ptr, dim=1)

        tree%root_ptr => build_tree_from_index_range(tree, 1, tree%n, bbox)
        return
        
        contains
        recursive function build_tree_from_index_range(tree, l, u, bbox) result(node)
            implicit none
            !Arguments
            type(kdtree), intent(inout) :: tree
            integer, intent(in) :: l, u
            real(fpp), intent(in) :: bbox(tree%n_dim,2)
            !Result
            type(tree_node), pointer :: node
            !Variables
            real(fpp) :: left_bbox(tree%n_dim,2), right_bbox(tree%n_dim,2)
            real(fpp) :: mid_point, min_value
            integer :: i, l_tmp, u_tmp, min_index, tmp_index
            integer :: tmp_split_dim(1)

            tree%n_nodes = tree%n_nodes + 1

            if(tree%n_nodes > tree%max_nodes) then
                print *, "Reached maximum number of nodes. Check for overlapping points. Aborting."
                stop
            end if

            allocate(node)
            allocate(node%bbox(tree%n_dim,2))
            node%bbox = bbox
            node%l_idx = l
            node%u_idx = u

            if(u-l+1 <= tree%leafsize) then
                !Terminal node
                nullify(node%left_ptr, node%right_ptr)
                if(u-l+1 > 0) then
                    !Calculate bounding box for leaf node
                    node%bbox(:,1) = minval(tree%data_ptr(tree%indicies(l:u), :), dim=1)
                    node%bbox(:,2) = maxval(tree%data_ptr(tree%indicies(l:u), :), dim=1)
                    !Debug information
                    tree%n_children = tree%n_children + u-l+1
                else
                    !Empty node
                    node%l_idx = 0
                    node%u_idx = 0
                end if
                !print *, "Hit terminal node. Number of childern: ", u-l+1
                return
            end if

            tmp_split_dim = maxloc(bbox(:,2) - bbox(:,1))
            node%split_dim = tmp_split_dim(1)

            mid_point = (bbox(node%split_dim,1) + bbox(node%split_dim,2))/2.0

            !Shuffle indicies. Adapted from Matt Kennel's kdtree2
            l_tmp = l
            u_tmp = u
            min_value = tree%data_ptr(tree%indicies(l_tmp), node%split_dim)
            min_index = l_tmp
            do while(l_tmp < u_tmp)
                if(tree%data_ptr(tree%indicies(l_tmp), node%split_dim) < min_value) then
                    min_value = tree%data_ptr(tree%indicies(l_tmp), node%split_dim)
                    min_index = l_tmp
                end if
                if(tree%data_ptr(tree%indicies(l_tmp), node%split_dim) <= mid_point) then
                    l_tmp = l_tmp + 1
                else
                    tmp_index = tree%indicies(l_tmp)
                    tree%indicies(l_tmp) = tree%indicies(u_tmp)
                    tree%indicies(u_tmp) = tmp_index
                    u_tmp = u_tmp - 1
                end if
            end do
            !Check final value, as l_tmp = u_tmp has not been checked yet
            if(tree%data_ptr(tree%indicies(l_tmp), node%split_dim) > mid_point) then
                l_tmp = l_tmp - 1
            end if

!            !Slide midpoint if either side is empty
!            if(l_tmp < l .or. l_tmp == u) then
!                mid_point = min_value
!                l_tmp = min_index
!                !Debugging
!                tree%n_mid_point_slides = tree%n_mid_point_slides + 1
!            end if

            node%split_value = mid_point
            left_bbox = bbox
            right_bbox = bbox
            left_bbox(node%split_dim, 2) = node%split_value
            right_bbox(node%split_dim, 1) = node%split_value

            !Check for empty nodes. If empty, choose same bbox as sibling node. Somewhat of a hack. Should implement proper midpoint sliding
            if(l > l_tmp) then
                left_bbox = right_bbox
            else if(l_tmp+1 > u) then
                right_bbox = left_bbox
            end if

            node%left_ptr => build_tree_from_index_range(tree, l, l_tmp, left_bbox)
            node%right_ptr => build_tree_from_index_range(tree, l_tmp+1, u, right_bbox)

            !Update bounding box
            node%bbox(:,1) = min(node%left_ptr%bbox(:,1), node%right_ptr%bbox(:,1))
            node%bbox(:,2) = max(node%left_ptr%bbox(:,2), node%right_ptr%bbox(:,2))
        end function build_tree_from_index_range
    end subroutine grow_tree

    subroutine search_within_r(tree, x, r, results, n_found, brute_force, debug_record)
        !Arguments
        type(kdtree), intent(in) :: tree
        real(fpp), intent(in) :: x(:)
        real(fpp), intent(in) :: r
        integer, intent(out) :: results(:)
        integer, intent(out) :: N_found
        logical, intent(in), optional :: brute_force
        type(search_record), optional :: debug_record
        !Variables
        logical :: brute_force_search

        if(.not. present(brute_force)) then
            brute_force_search = .true.
        else
            brute_force_search = brute_force
        end if

        n_found = 0

        if(brute_force_search) then
            call search_index_range_within_r(tree, x, r, tree%root_ptr%l_idx, tree%root_ptr%u_idx, results, n_found)
        else
            call search_node_within_r(tree, tree%root_ptr, x, r, results, n_found, debug_record)
        end if

        contains
        recursive subroutine search_node_within_r(tree, node, x, r, results, n_found, debug_record)
            !Arguments
            type(kdtree), intent(in) :: tree
            type(tree_node), pointer, intent(in) :: node
            real(fpp), intent(in) :: x(:)
            real(fpp), intent(in) :: r
            integer, intent(out) :: results(:)
            integer, intent(inout) :: n_found
            type(search_record), optional :: debug_record
            !Variables
            real(fpp) :: distance_to_split
            integer :: i

            if(present(debug_record)) then
                debug_record%n_included_nodes = debug_record%n_included_nodes + 1
                debug_record%included_nodes(debug_record%n_included_nodes)%ptr => node
            end if

            if(node%l_idx == 0) then
                !Empty node
                return
            end if

            !Check if bounding box is completely outside of r
            if(min_distance_to_box(node%bbox, x) > r) then
                !print *, "distance to box: ", min_distance_to_box(node%bbox, x)
                return
            end if

            !Check if bounding box is completely within r
            if(max_distance_to_box(node%bbox, x) <= r) then
                if(n_found + node%u_idx -node%l_idx+1 > size(results)) then
                    print *, "Result array to small. Found more than", n_found + node%u_idx -node%l_idx+1
                    stop
                end if
                do i=node%l_idx,node%u_idx
                    n_found = n_found + 1
                    results(n_found) = tree%indicies(i)
                end do
                return
            end if

            if(.not. associated(node%left_ptr) .and. .not. associated(node%right_ptr)) then
                !Terminal node
                call search_index_range_within_r(tree, x, r, node%l_idx, node%u_idx, results, n_found)
                return
            end if

            distance_to_split = x(node%split_dim) - node%split_value

            if(distance_to_split**2 > r**2) then
                !only one side
                if(distance_to_split < 0) then
                    !left
                    call search_node_within_r(tree, node%left_ptr, x, r, results, n_found, debug_record)
                else
                    !right
                    call search_node_within_r(tree, node%right_ptr, x, r, results, n_found, debug_record)
                end if
            else
                !Search both child nodes
                call search_node_within_r(tree, node%left_ptr, x, r, results, n_found, debug_record)
                call search_node_within_r(tree, node%right_ptr, x, r, results, n_found, debug_record)
            end if
        end subroutine search_node_within_r

        subroutine search_index_range_within_r(tree, x, r, l, u, results, n_found)
            !Arguments
            type(kdtree), intent(in) :: tree
            real(fpp), intent(in) :: x(:)
            real(fpp), intent(in) :: r
            integer, intent(in) :: l, u
            integer, intent(out) :: results(:)
            integer, intent(inout) :: n_found
            !Variables
            integer :: i
            !print *, tree%indicies(l), tree%indicies(u)
            do i=l,u
                if(sum((tree%data_ptr(tree%indicies(i),:) - x(:))**2) <= r**2) then
                    n_found = n_found + 1
                    if(n_found > size(results)) then
                        print *, "Result array to small. Found at least", n_found
                        return
                    end if
                    results(n_found) = tree%indicies(i)
                end if
            end do
        end subroutine search_index_range_within_r

        function max_distance_to_box(bbox, x) result(distance)
            !Arguments
            real(fpp), intent(in) :: bbox(:,:)
            real(fpp), intent(in) :: x(:)
            !Result
            real(fpp) :: distance
            !Variables
            integer :: i, n_dim
            real(fpp) :: dl, dr

            distance = 0.0
            n_dim = size(bbox, 1)
            do i=1,n_dim
                dl = (bbox(i,1) - x(i))**2
                dr = (bbox(i,2) - x(i))**2
                if(dl > dr) then
                    distance = distance + dl
                else
                    distance = distance + dr
                end if
            end do
            distance = sqrt(distance)
        end function max_distance_to_box

        function min_distance_to_box(bbox, x) result(distance)
            !Arguments
            real(fpp), intent(in) :: bbox(:,:)
            real(fpp), intent(in) :: x(:)
            !Result
            real(fpp) :: distance
            !Variables
            integer :: i, n_dim
            real(fpp) :: dl, dr, dl2, dr2

            distance = 0.0
            n_dim = size(bbox, 1)
            do i=1,n_dim
                dl = (bbox(i,1) - x(i))
                dr = (bbox(i,2) - x(i))
                dl2 = dl**2
                dr2 = dr**2
                if(dl<0.0 .and. dr>0.0) then
                    cycle
                end if
                if(dl2 > dr2) then
                    distance = distance + dr2
                else
                    distance = distance + dl2
                end if
            end do
            distance = sqrt(distance)
        end function min_distance_to_box
    end subroutine search_within_r

    subroutine search_within_r_spherical(tree, x, r, results, n_found, brute_force, debug_record)
        !Arguments
        type(kdtree), intent(in) :: tree
        real(fpp), intent(in) :: x(:)
        real(fpp), intent(in) :: r
        integer, intent(out) :: results(:)
        integer, intent(out) :: N_found
        logical, intent(in), optional :: brute_force
        type(search_record), optional :: debug_record
        !Variables
        logical :: brute_force_search
        real :: longitude_reach

        if(.not. present(brute_force)) then
            brute_force_search = .true.
        else
            brute_force_search = brute_force
        end if

        n_found = 0

        !Check for poles in search region
        if(x(2) + r >= pi/2.0) then
            !North pole included
            print *, "Search region includes north pole. Aborting."
            stop
        else if(x(2) - r <= -pi/2.0) then
            !South pole included
            print *, "Search region includes south pole. Aborting."
            stop
        end if
        !Check for prime meridian
        if(distance_to_meridian_sphere(x(1), x(2), 0.0) <= r) then
            !Prime meridian included
            print *, "Search region includes prime meridian. Aborting."
            stop
        end if

        if(brute_force_search) then
            call search_index_range_within_r(tree, x, r, tree%root_ptr%l_idx, tree%root_ptr%u_idx, results, n_found)
        else
            call search_node_within_r(tree, tree%root_ptr, x, r, results, n_found, debug_record)
        end if

        contains
        recursive subroutine search_node_within_r(tree, node, x, r, results, n_found, debug_record)
            !Arguments
            type(kdtree), intent(in) :: tree
            type(tree_node), pointer, intent(in) :: node
            real(fpp), intent(in) :: x(:)
            real(fpp), intent(in) :: r
            integer, intent(out) :: results(:)
            integer, intent(inout) :: n_found
            type(search_record), optional :: debug_record
            !Variables
            real(fpp) :: distance_to_split
            integer :: i
            logical :: search_left_node, search_right_node

            if(present(debug_record)) then
                debug_record%n_included_nodes = debug_record%n_included_nodes + 1
                debug_record%included_nodes(debug_record%n_included_nodes)%ptr => node
            end if

            if(node%l_idx == 0) then
                !Empty node
                return
            end if

            !Check if bounding box is completely outside of r
            if(min_distance_to_box(node%bbox, x, r) > r) then
                !print *, "distance to box: ", min_distance_to_box(node%bbox, x)
                return
            end if
            !Check if bounding box is completely within r
            if(max_distance_to_box(node%bbox, x) <= r) then
                if(n_found + node%u_idx -node%l_idx+1 > size(results)) then
                    print *, "Result array to small. Found more than", n_found + node%u_idx -node%l_idx+1
                    stop
                end if
                do i=node%l_idx,node%u_idx
                    n_found = n_found + 1
                    results(n_found) = tree%indicies(i)
                end do
                return
            end if

            if(.not. associated(node%left_ptr) .and. .not. associated(node%right_ptr)) then
                !Terminal node
                call search_index_range_within_r(tree, x, r, node%l_idx, node%u_idx, results, n_found)
                return
            end if

            search_left_node = .true.
            search_right_node = .true.
            if(node%split_dim == 1) then
                !longitude
                if(distance_to_meridian_sphere(x(1), x(2), node%split_value) > r) then
                    if(x(1) < node%split_value) then
                        search_right_node = .false.
                    else
                        search_left_node = .false.
                    end if
                end if
            else if(node%split_dim == 2) then
                !latitude
                if(abs(x(2) - node%split_value) > r) then
                    if(x(2) < node%split_value) then
                        search_right_node = .false.
                    else
                        search_left_node = .false.
                    end if
                end if
            else
                !no support for more dimensions
                print *, "Only 2 dimensions are supported in sqherical coordinates."
                return
            end if

            if(search_left_node) call search_node_within_r(tree, node%left_ptr, x, r, results, n_found, debug_record)
            if(search_right_node) call search_node_within_r(tree, node%right_ptr, x, r, results, n_found, debug_record)

        end subroutine search_node_within_r

        function max_distance_to_box(bbox, x) result(distance)
            !Arguments
            real(fpp), intent(in) :: bbox(:,:)
            real(fpp), intent(in) :: x(:)
            !Result
            real(fpp) :: distance
            !Variables
            integer :: i, j
            real(fpp) :: d

            distance = 0.0
            do i=1,2
                do j=1,2
                    d = distance_sphere(x(1), x(2), bbox(1,i), bbox(2,j))
                    if(d > distance) distance = d
                end do
            end do
        end function max_distance_to_box

        function min_distance_to_box(bbox, x, r) result(distance)
            !Arguments
            real(fpp), intent(in) :: bbox(:,:)
            real(fpp), intent(in) :: x(:)
            real(fpp), intent(in) :: r
            !Result
            real(fpp) :: distance
            !Variables
            integer :: i, j
            real(fpp) :: d, distance_edge, distance_corner
            logical :: within_x, within_y

            within_x = (bbox(1,1) < x(1)) .and. (x(1) < bbox(1,2))
            within_y = (bbox(2,1) < x(2)) .and. (x(2) < bbox(2,2))
            distance = 2.0*pi
            if(within_x .and. within_y) then
                !Within the box
                distance = 0.0
            else if(within_x) then
                !Within x but outside along y.
                distance = minval(abs(x(2) - bbox(2,:)))
            else if(within_y) then
                !Within y but outside along x.
                if(abs(x(1) - bbox(1,1)) < abs(x(1) - bbox(1,2))) then
                    !left
                    distance_edge = distance_to_meridian_sphere(x(1), x(2), bbox(1,1))
                    distance_corner = minval(distance_sphere(x(1), x(2), bbox(1,1), bbox(2,:)))
                else
                    !right
                    distance_edge = distance_to_meridian_sphere(x(1), x(2), bbox(1,2))
                    distance_corner = minval(distance_sphere(x(1), x(2), bbox(1,2), bbox(2,:)))
                end if
                distance = minval([distance_edge, distance_corner])
            else
                do i=1,2
                    do j=1,2
                        d = distance_sphere(x(1), x(2), bbox(1,i), bbox(2,j))
                        if(d < distance) distance = d
                    end do
                end do
            end if
        end function min_distance_to_box

        subroutine search_index_range_within_r(tree, x, r, l, u, results, n_found)
            !Arguments
            type(kdtree), intent(in) :: tree
            real(fpp), intent(in) :: x(:)
            real(fpp), intent(in) :: r
            integer, intent(in) :: l, u
            integer, intent(out) :: results(:)
            integer, intent(inout) :: n_found
            !Variables
            integer :: i
            real :: d
            !print *, tree%indicies(l), tree%indicies(u)
            do i=l,u
                d = distance_sphere(tree%data_ptr(tree%indicies(i),1), &
                                    tree%data_ptr(tree%indicies(i),2), &
                                    x(1), &
                                    x(2))
                if(d <= r) then
                    n_found = n_found + 1
                    if(n_found > size(results)) then
                        print *, "Result array to small. Found at least", n_found
                        return
                    end if
                    results(n_found) = tree%indicies(i)
                end if
            end do
        end subroutine search_index_range_within_r
    end subroutine search_within_r_spherical

    recursive subroutine write_node_data_to_file(tree, node, fileunit)
        !Arguments
        type(kdtree), intent(in) :: tree
        type(tree_node), pointer, intent(in) :: node
        integer, intent(in) :: fileunit
        !Variables
        logical :: terminal

        terminal = .false.
        if((associated(node%left_ptr) .and. associated(node%right_ptr))) then
            terminal = .false.
        else
            terminal = .true.
        end if

        write(unit=fileunit, fmt="(2I, 4F15.3, I)") node%l_idx, node%u_idx, node%bbox(1,1), node%bbox(1,2), node%bbox(2,1), node%bbox(2,2), terminal
        if(.not. terminal) then
            call write_node_data_to_file(tree, node%left_ptr, fileunit)
            call write_node_data_to_file(tree, node%right_ptr, fileunit)
        end if
    end subroutine write_node_data_to_file

    subroutine write_search_debug_data_to_file(tree, debug_record, fileunit)
        !Arguments
        type(kdtree), intent(in) :: tree
        type(search_record), intent(in) :: debug_record
        integer, intent(in) :: fileunit
        !Variables
        integer :: i
        type(tree_node), pointer :: node
        logical :: terminal

        terminal = .false.
        do i=1,debug_record%n_included_nodes
            node => debug_record%included_nodes(i)%ptr
            if((associated(node%left_ptr) .and. associated(node%right_ptr))) then
                terminal = .false.
            else
                terminal = .true.
            end if

            write(unit=fileunit, fmt="(2I, 4F15.3, I)") node%l_idx, node%u_idx, node%bbox(1,1), node%bbox(1,2), node%bbox(2,1), node%bbox(2,2), terminal
        end do
    end subroutine write_search_debug_data_to_file

end module kdtree_module