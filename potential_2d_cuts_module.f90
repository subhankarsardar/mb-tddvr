  ! ============================================================================
! 2D POTENTIAL CUTS MODULE
! Creates 2D potential surfaces along specified mode pairs for all electronic
! surfaces
! Compatible with SO-TDDVR.f90
! Based on: Journal of Molecular Structure 1110 (2016) 32-43
! ============================================================================

 module potential_2d_cuts_module
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  
  private
  
  ! Public types and procedures
  public :: potential_2d_analyzer, dp, ip
  public :: create_2d_potential_analyzer
  
  ! Kind parameters
  integer, parameter :: dp = real64
  integer, parameter :: ip = int32
  
  ! Default mode pairs to analyze
  integer, parameter :: DEFAULT_NPAIRS = 4
  integer, parameter :: DEFAULT_PAIRS(2,4) = reshape([ &
    1, 2, &   ! Mode pair 1: 1-2
    5, 8, &   ! Mode pair 2: 5-8  
    8, 15, &  ! Mode pair 3: 8-15
    5, 18 &   ! Mode pair 4: 5-18
  ], [2, DEFAULT_NPAIRS])
  
  ! Plotting parameters
  integer, parameter :: N_POINTS_2D = 50
  real(dp), parameter :: CUT_RANGE_2D = 4.0_dp
  real(dp), parameter :: SMALL = 1.0e-12_dp
  
  ! Type for 2D potential analyzer
  type :: potential_2d_analyzer
    private
    logical :: enabled = .false.
    logical :: verbose = .true.
    logical :: parallel = .true.
    
    ! Mode pairs to analyze
    integer(ip) :: npairs_to_analyze = DEFAULT_NPAIRS
    integer(ip), allocatable :: mode_pairs(:,:)
    
    ! Output control
    integer :: output_frequency = 100
    character(len=100) :: output_dir = "./potential_2d_cuts/"
    
    ! Internal storage for each pair
    real(dp), allocatable :: potential_2d(:,:,:,:)
    real(dp), allocatable :: coord_x(:,:)
    real(dp), allocatable :: coord_y(:,:)
    real(dp), allocatable :: min_potentials(:,:)
    real(dp), allocatable :: max_potentials(:,:)
    real(dp), allocatable :: saddle_points(:,:,:)
    
    ! File units
    integer :: unit_summary = -1
    integer :: unit_gnuplot = -1
    
    ! Counters and state info
    integer :: cut_counter = 0
    integer(ip) :: nstate = 0
    integer(ip) :: nmode = 0
    
  contains
    procedure :: initialize => init_2d_potential_analyzer
    procedure :: set_mode_pairs
    procedure :: compute_2d_potential_surfaces
    procedure :: find_critical_points
    procedure :: save_2d_potential_data
    procedure :: save_surface_data
    procedure :: generate_contour_plots
    procedure :: generate_surface_plots
    procedure :: save_summary_data
    procedure :: finalize => cleanup_2d_potential_analyzer
    procedure :: is_enabled
    procedure :: enable
    procedure :: disable
    procedure :: set_output_frequency
    procedure :: get_critical_points
  end type potential_2d_analyzer
  
  ! Constructor interface
  interface potential_2d_analyzer
    module procedure create_2d_potential_analyzer
  end interface

  contains

  ! ============================================================================
  ! CONSTRUCTOR
  ! ============================================================================
  function create_2d_potential_analyzer(enabled, verbose, parallel, &
                                       pairs_list, output_freq) result(analyzer)
    logical, intent(in), optional :: enabled, verbose, parallel
    integer(ip), intent(in), optional :: pairs_list(:,:)
    integer, intent(in), optional :: output_freq
    type(potential_2d_analyzer) :: analyzer
    
    integer :: i
    
    if (present(enabled)) analyzer%enabled = enabled
    if (present(verbose)) analyzer%verbose = verbose
    if (present(parallel)) analyzer%parallel = parallel
    if (present(output_freq)) analyzer%output_frequency = output_freq
    
    if (present(pairs_list)) then
      analyzer%npairs_to_analyze = size(pairs_list, 2)
      allocate(analyzer%mode_pairs(2, analyzer%npairs_to_analyze))
      analyzer%mode_pairs = pairs_list
    else
      allocate(analyzer%mode_pairs(2, DEFAULT_NPAIRS))
      analyzer%mode_pairs = DEFAULT_PAIRS
    end if
    
    if (analyzer%verbose .and. analyzer%enabled) then
      print '(a)', "2D Potential cuts analyzer created (enabled)"
      print '(a,i0)', "  Mode pairs to analyze: ", analyzer%npairs_to_analyze
      do i = 1, analyzer%npairs_to_analyze
        print '(a,i0,a,i0,a,i0)', "    Pair ", i, ": Modes ", &
              analyzer%mode_pairs(1,i), "-", analyzer%mode_pairs(2,i)
      end do
      print '(a,i0)', "  Grid points: ", N_POINTS_2D, "x", N_POINTS_2D
      print '(a,f6.2)', "  Cut range: ±", CUT_RANGE_2D
      print '(a,i0)', "  Output frequency: every ", analyzer%output_frequency, " steps"
      if (analyzer%parallel) print '(a)', "  Parallel execution: enabled"
    end if
  end function create_2d_potential_analyzer

  ! ============================================================================
  ! INITIALIZE 2D POTENTIAL ANALYZER
  ! ============================================================================
  subroutine init_2d_potential_analyzer(this, nstate, nmode, freqn, ffn)
    class(potential_2d_analyzer), intent(inout) :: this
    integer(ip), intent(in) :: nstate, nmode
    real(dp), intent(in) :: freqn(:), ffn(:)
    
    integer :: ios, i, p, m1, m2
    character(len=200) :: dir_command
    real(dp) :: dx, dy
    
    if (.not. this%enabled) return
    
    this%nstate = nstate
    this%nmode = nmode
    
    if (this%verbose) then
      print '(a)', "Initializing 2D potential cuts analyzer..."
      print '(a,i0)', "  Number of states: ", nstate
      print '(a,i0)', "  Total modes: ", nmode
    end if
    
    ! Check mode indices are valid
    do i = 1, this%npairs_to_analyze
      m1 = this%mode_pairs(1,i)
      m2 = this%mode_pairs(2,i)
      
      if (m1 < 1 .or. m1 > nmode .or. m2 < 1 .or. m2 > nmode) then
        print '(a,i0,a,i0,a)', "ERROR: Mode pair ", m1, "-", m2, " is out of bounds!"
        this%enabled = .false.
        return
      end if
      
      if (m1 == m2) then
        print '(a,i0,a)', "ERROR: Mode pair has identical modes: ", m1
        this%enabled = .false.
        return
      end if
    end do
    
    ! Allocate storage
    allocate(this%potential_2d(nstate, N_POINTS_2D, N_POINTS_2D, this%npairs_to_analyze))
    allocate(this%coord_x(N_POINTS_2D, this%npairs_to_analyze))
    allocate(this%coord_y(N_POINTS_2D, this%npairs_to_analyze))
    allocate(this%min_potentials(nstate, this%npairs_to_analyze))
    allocate(this%max_potentials(nstate, this%npairs_to_analyze))
    allocate(this%saddle_points(3, nstate, this%npairs_to_analyze))
    
    ! Initialize coordinate grids
    dx = 2.0_dp * CUT_RANGE_2D / real(N_POINTS_2D - 1, dp)
    dy = dx
    
    do i = 1, this%npairs_to_analyze
      do p = 1, N_POINTS_2D
        this%coord_x(p, i) = -CUT_RANGE_2D + (p-1) * dx
        this%coord_y(p, i) = -CUT_RANGE_2D + (p-1) * dy
        
        m1 = this%mode_pairs(1,i)
        m2 = this%mode_pairs(2,i)
        
        if (freqn(m1) > SMALL) then
          this%coord_x(p, i) = this%coord_x(p, i) / sqrt(2.0_dp * freqn(m1))
        end if
        if (freqn(m2) > SMALL) then
          this%coord_y(p, i) = this%coord_y(p, i) / sqrt(2.0_dp * freqn(m2))
        end if
      end do
    end do
    
    this%potential_2d = 0.0_dp
    this%min_potentials = huge(1.0_dp)
    this%max_potentials = -huge(1.0_dp)
    this%saddle_points = 0.0_dp
    this%cut_counter = 0
    
    ! Create output directory
    dir_command = "mkdir -p " // trim(this%output_dir)
    call system(dir_command)
    
    ! Create subdirectories for each mode pair
    do i = 1, this%npairs_to_analyze
      m1 = this%mode_pairs(1,i)
      m2 = this%mode_pairs(2,i)
      write(dir_command, '(a,a,a,i0,a,i0)') &
            "mkdir -p ", trim(this%output_dir), "mode", m1, "_", m2
      call system(dir_command)
    end do
    
    ! Open summary file
    open(newunit=this%unit_summary, &
         file=trim(this%output_dir)//"potential_2d_summary.dat", &
         status='replace', iostat=ios)
    if (ios /= 0) then
      print '(a)', "ERROR: Could not open summary file!"
      this%enabled = .false.
      return
    end if
    
    write(this%unit_summary, '(a)') '# 2D Potential Analysis Summary'
    write(this%unit_summary, '(a)') '# Time(fs), Mode1, Mode2, State, Min_Pot, Max_Pot, Saddle_X, Saddle_Y, Saddle_E'
    
    if (this%verbose) then
      print '(a)', "2D potential cuts analyzer initialized successfully."
      print '(a,a)', "  Output directory: ", trim(this%output_dir)
    end if
  end subroutine init_2d_potential_analyzer

  ! ============================================================================
  ! SET MODE PAIRS
  ! ============================================================================
  subroutine set_mode_pairs(this, pairs_list)
    class(potential_2d_analyzer), intent(inout) :: this
    integer(ip), intent(in) :: pairs_list(:,:)
    
    integer :: i
    
    if (allocated(this%mode_pairs)) deallocate(this%mode_pairs)
    
    this%npairs_to_analyze = size(pairs_list, 2)
    allocate(this%mode_pairs(2, this%npairs_to_analyze))
    this%mode_pairs = pairs_list
    
    if (this%verbose .and. this%enabled) then
      print '(a,i0)', "Updated mode pairs to analyze: ", this%npairs_to_analyze
      do i = 1, this%npairs_to_analyze
        print '(a,i0,a,i0,a,i0)', "  Pair ", i, ": Modes ", &
              this%mode_pairs(1,i), "-", this%mode_pairs(2,i)
      end do
    end if
  end subroutine set_mode_pairs

  ! ============================================================================
  ! COMPUTE 2D POTENTIAL SURFACES
  ! ============================================================================
  subroutine compute_2d_potential_surfaces(this, time, step, potential_func, &
                                          freqn, ffn, center_coords)
    class(potential_2d_analyzer), intent(inout) :: this
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    interface
      function potential_func(coords, nstate) result(v)
        import dp, ip
        real(dp), intent(in) :: coords(:)
        integer(ip), intent(in) :: nstate
        real(dp) :: v(nstate, nstate)
      end function potential_func
    end interface
    real(dp), intent(in) :: freqn(:), ffn(:)
    real(dp), intent(in) :: center_coords(:)
    
    integer(ip) :: i, j, k, s, m1, m2, nstate_local, nmode_local
    real(dp), allocatable :: coords(:)
    real(dp), allocatable :: v_diag(:,:)
    
    if (.not. this%enabled) return
    if (mod(step, this%output_frequency) /= 0) return
    
    nstate_local = this%nstate
    nmode_local = this%nmode
    
    this%cut_counter = this%cut_counter + 1
    
    if (this%verbose .and. mod(step, 200) == 0) then
      print '(a,f8.3,a,i0)', "Computing 2D potential surfaces at t=", time, &
            " fs (step ", step, ")"
    end if
    
    allocate(coords(nmode_local))
    allocate(v_diag(nstate_local, nstate_local))
    
    do i = 1, this%npairs_to_analyze
      m1 = this%mode_pairs(1,i)
      m2 = this%mode_pairs(2,i)
      
      do j = 1, N_POINTS_2D
        do k = 1, N_POINTS_2D
          coords = center_coords
          coords(m1) = this%coord_x(j, i)
          coords(m2) = this%coord_y(k, i)
          
          v_diag = potential_func(coords, nstate_local)
          
          do s = 1, nstate_local
            this%potential_2d(s, j, k, i) = v_diag(s, s)
          end do
        end do
      end do
      
      do s = 1, nstate_local
        this%min_potentials(s, i) = min(this%min_potentials(s, i), &
                                       minval(this%potential_2d(s, :, :, i)))
        this%max_potentials(s, i) = max(this%max_potentials(s, i), &
                                       maxval(this%potential_2d(s, :, :, i)))
      end do
    end do
    
    deallocate(coords, v_diag)
    
    call this%find_critical_points()
    call this%save_2d_potential_data(time, step)
    call this%save_summary_data(time, step)
    
    if (mod(step, this%output_frequency * 5) == 0) then
      call this%generate_contour_plots()
      call this%generate_surface_plots()
    end if
  end subroutine compute_2d_potential_surfaces

  ! ============================================================================
  ! FIND CRITICAL POINTS
  ! ============================================================================
  subroutine find_critical_points(this)
    class(potential_2d_analyzer), intent(inout) :: this
    
    integer(ip) :: i, s, j, k, nstate_local, m1, m2
    real(dp) :: v, v_up, v_down, v_left, v_right
    real(dp) :: v_up_left, v_up_right, v_down_left, v_down_right
    real(dp) :: v_min, v_max, saddle_energy
    integer :: min_j, min_k, max_j, max_k, saddle_j, saddle_k
    integer :: min_loc(2), max_loc(2)
    logical :: found_saddle
    
    if (.not. this%enabled) return
    
    nstate_local = this%nstate
    
    do i = 1, this%npairs_to_analyze
      do s = 1, nstate_local
        v_min = minval(this%potential_2d(s, :, :, i))
        v_max = maxval(this%potential_2d(s, :, :, i))
        
        min_loc = minloc(this%potential_2d(s, :, :, i))
        max_loc = maxloc(this%potential_2d(s, :, :, i))
        
        min_j = min_loc(1)
        min_k = min_loc(2)
        max_j = max_loc(1)
        max_k = max_loc(2)
        
        saddle_energy = 0.0_dp
        saddle_j = 0
        saddle_k = 0
        found_saddle = .false.
        
        do j = 2, N_POINTS_2D-1
          do k = 2, N_POINTS_2D-1
            v = this%potential_2d(s, j, k, i)
            v_up = this%potential_2d(s, j, k+1, i)
            v_down = this%potential_2d(s, j, k-1, i)
            v_left = this%potential_2d(s, j-1, k, i)
            v_right = this%potential_2d(s, j+1, k, i)
            
            if ((v > v_up .and. v > v_down .and. v < v_left .and. v < v_right) .or. &
                (v < v_up .and. v < v_down .and. v > v_left .and. v > v_right)) then
              v_up_left = this%potential_2d(s, j-1, k+1, i)
              v_up_right = this%potential_2d(s, j+1, k+1, i)
              v_down_left = this%potential_2d(s, j-1, k-1, i)
              v_down_right = this%potential_2d(s, j+1, k-1, i)
              
              if (abs(v - v_up_left) > 1.0e-6_dp .or. &
                  abs(v - v_up_right) > 1.0e-6_dp .or. &
                  abs(v - v_down_left) > 1.0e-6_dp .or. &
                  abs(v - v_down_right) > 1.0e-6_dp) then
                saddle_energy = v
                saddle_j = j
                saddle_k = k
                found_saddle = .true.
                exit
              end if
            end if
          end do
          if (found_saddle) exit
        end do
        
        if (found_saddle) then
          this%saddle_points(1, s, i) = this%coord_x(saddle_j, i)
          this%saddle_points(2, s, i) = this%coord_y(saddle_k, i)
          this%saddle_points(3, s, i) = saddle_energy
        else
          this%saddle_points(1, s, i) = 0.5_dp * (this%coord_x(min_j, i) + this%coord_x(max_j, i))
          this%saddle_points(2, s, i) = 0.5_dp * (this%coord_y(min_k, i) + this%coord_y(max_k, i))
          this%saddle_points(3, s, i) = 0.5_dp * (v_min + v_max)
        end if
        
        if (this%verbose .and. mod(s, 2) == 0) then
          m1 = this%mode_pairs(1,i)
          m2 = this%mode_pairs(2,i)
          print '(a,i0,a,i0,a,i0,a,f12.6,a,f8.4,a,f8.4,a)', &
                "Mode ", m1, "-", m2, ", State ", s, &
                ": Min = ", v_min, " at (", this%coord_x(min_j,i), ",", this%coord_y(min_k,i), ")"
          if (found_saddle) then
            print '(a,f12.6,a,f8.4,a,f8.4,a)', & 
                  "Saddle = ", saddle_energy, & 
                  " at (", this%coord_x(saddle_j,i), ",", this%coord_y(saddle_k,i), ")"
          end if
        end if
      end do
    end do
  end subroutine find_critical_points

  ! ============================================================================
  ! SAVE 2D POTENTIAL DATA
  ! ============================================================================
  subroutine save_2d_potential_data(this, time, step)
    class(potential_2d_analyzer), intent(in) :: this
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    
    character(len=200) :: filename, dirname
    integer :: i, j, k, s, unit_file, m1, m2, ios
    integer(ip) :: nstate_local
    character(len=30) :: fmt_str
    
    if (.not. this%enabled) return
    
    nstate_local = this%nstate
    write(fmt_str, '(a,i0,a)') '(2f12.6,', nstate_local, 'f16.8)'
    
    do i = 1, this%npairs_to_analyze
      m1 = this%mode_pairs(1,i)
      m2 = this%mode_pairs(2,i)
      
      write(dirname, '(a,a,i0,a,i0,a)') trim(this%output_dir), "mode", m1, "_", m2, "/"
      
      do s = 1, nstate_local
        write(filename, '(a,a,i0,a,f8.3,a)') &
              trim(dirname), "potential_state", s, "_t", time, ".dat"
        
        open(newunit=unit_file, file=trim(filename), status='replace', iostat=ios)
        if (ios /= 0) cycle
        
        write(unit_file, '(a)') '# 2D Potential Surface'
        write(unit_file, '(a,f12.6,a,i0)') '# Time: ', time, ' fs, Step: ', step
        write(unit_file, '(a,i0,a,i0)') '# Modes: ', m1, ' - ', m2
        write(unit_file, '(a,i0)') '# State: ', s
        write(unit_file, '(a)') '# X_Coord, Y_Coord, Potential'
        
        do j = 1, N_POINTS_2D
          do k = 1, N_POINTS_2D
            write(unit_file, '(3f16.8)') &
                  this%coord_x(j, i), &
                  this%coord_y(k, i), &
                  this%potential_2d(s, j, k, i)
          end do
          write(unit_file, *)
        end do
        
        close(unit_file)
      end do
      
      write(filename, '(a,a,f8.3,a)') &
            trim(dirname), "potential_all_states_t", time, ".dat"
      
      open(newunit=unit_file, file=trim(filename), status='replace', iostat=ios)
      if (ios /= 0) cycle
      
      write(unit_file, '(a)') '# Combined 2D Potential Surfaces for All States'
      write(unit_file, '(a,f12.6)') '# Time: ', time
      write(unit_file, '(a,i0,a,i0)') '# Modes: ', m1, ' - ', m2
      write(unit_file, '(a,i0)') '# Number of states: ', nstate_local
      write(unit_file, '(a)') '# X_Coord, Y_Coord, State1, State2, ...'
      
      do j = 1, N_POINTS_2D
        do k = 1, N_POINTS_2D
          write(unit_file, fmt_str) &
                this%coord_x(j, i), this%coord_y(k, i), &
                (this%potential_2d(s, j, k, i), s=1, min(nstate_local, 6))
        end do
        write(unit_file, *)
      end do
      
      close(unit_file)
    end do
  end subroutine save_2d_potential_data

  ! ============================================================================
  ! SAVE SUMMARY DATA
  ! ============================================================================
  subroutine save_summary_data(this, time, step)
    class(potential_2d_analyzer), intent(in) :: this
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    
    integer(ip) :: i, s, m1, m2, nstate_local
    
    if (.not. this%enabled) return
    if (this%unit_summary < 0) return
    
    nstate_local = this%nstate
    
    do i = 1, this%npairs_to_analyze
      m1 = this%mode_pairs(1,i)
      m2 = this%mode_pairs(2,i)
      
      do s = 1, nstate_local
        write(this%unit_summary, '(f12.6,4i6,5f16.8)') &
              time, m1, m2, s, &
              this%min_potentials(s, i), &
              this%max_potentials(s, i), &
              this%saddle_points(1, s, i), &
              this%saddle_points(2, s, i), &
              this%saddle_points(3, s, i)
      end do
    end do
    
    flush(this%unit_summary)
  end subroutine save_summary_data

  ! ============================================================================
  ! SAVE SURFACE DATA (VTK)
  ! ============================================================================
  subroutine save_surface_data(this)
    class(potential_2d_analyzer), intent(in) :: this
    
    character(len=200) :: filename
    integer :: i, j, k, s, unit_file, m1, m2, ios
    integer(ip) :: nstate_local
    
    if (.not. this%enabled) return
    
    nstate_local = this%nstate
    
    do i = 1, this%npairs_to_analyze
      m1 = this%mode_pairs(1,i)
      m2 = this%mode_pairs(2,i)
      
      write(filename, '(a,a,i0,a,i0,a)') &
            trim(this%output_dir), "surface_mode", m1, "_", m2, ".vtk"
      
      open(newunit=unit_file, file=trim(filename), status='replace', iostat=ios)
      if (ios /= 0) cycle
      
      write(unit_file, '(a)') '# vtk DataFile Version 3.0'
      write(unit_file, '(a)') '2D Potential Surfaces'
      write(unit_file, '(a)') 'ASCII'
      write(unit_file, '(a)') 'DATASET STRUCTURED_GRID'
      write(unit_file, '(a,3i8)') 'DIMENSIONS ', N_POINTS_2D, N_POINTS_2D, 1
      write(unit_file, '(a,i8,a)') 'POINTS ', N_POINTS_2D*N_POINTS_2D, ' float'
      
      do k = 1, N_POINTS_2D
        do j = 1, N_POINTS_2D
          write(unit_file, '(3f12.6)') &
                this%coord_x(j, i), this%coord_y(k, i), 0.0_dp
        end do
      end do
      
      write(unit_file, '(a,i8)') 'POINT_DATA ', N_POINTS_2D*N_POINTS_2D
      
      do s = 1, min(nstate_local, 6)
        write(unit_file, '(a,i0,a)') 'SCALARS State', s, '_Potential float 1'
        write(unit_file, '(a)') 'LOOKUP_TABLE default'
        
        do k = 1, N_POINTS_2D
          do j = 1, N_POINTS_2D
            write(unit_file, '(f16.8)') this%potential_2d(s, j, k, i)
          end do
        end do
      end do
      
      close(unit_file)
    end do
  end subroutine save_surface_data

  ! ============================================================================
  ! GENERATE CONTOUR PLOTS
  ! ============================================================================
  subroutine generate_contour_plots(this)
    class(potential_2d_analyzer), intent(in) :: this
    
    character(len=200) :: scriptname, filename, dirname
    integer :: unit, i, s, m1, m2, ios
    integer(ip) :: nstate_local
    
    if (.not. this%enabled) return
    
    nstate_local = this%nstate
    
    do i = 1, this%npairs_to_analyze
      m1 = this%mode_pairs(1,i)
      m2 = this%mode_pairs(2,i)
      
      write(dirname, '(a,a,i0,a,i0,a)') trim(this%output_dir), "mode", m1, "_", m2, "/"
      write(scriptname, '(a,a)') trim(dirname), "plot_contours.gnu"
      
      open(newunit=unit, file=trim(scriptname), status='replace', iostat=ios)
      if (ios /= 0) cycle
      
      write(unit, '(a)') '# Gnuplot script for 2D potential contour plots'
      write(unit, '(a)') 'set terminal pngcairo enhanced size 1600,1200'
      write(unit, '(a)') 'set output "' // trim(dirname) // 'potential_contours.png"'
      write(unit, '(a)') 'set multiplot layout ' // trim(int_to_str((nstate_local+1)/2)) // ',2'
      write(unit, '(a)') ''
      write(unit, '(a,i0,a,i0,a)') 'set title "Mode ', m1, '-', m2, ' - 2D Potential Contours"'
      write(unit, '(a)') ''
      
      do s = 1, min(nstate_local, 6)
        write(unit, '(a,i0,a)') 'set title "State ', s, '"'
        write(unit, '(a)') 'set xlabel "Mode ' // trim(int_to_str(m1)) // '"'
        write(unit, '(a)') 'set ylabel "Mode ' // trim(int_to_str(m2)) // '"'
        write(unit, '(a)') 'set pm3d map'
        write(unit, '(a)') 'set palette defined (0 "blue", 0.5 "green", 1 "red")'
        write(unit, '(a)') 'set contour base'
        write(unit, '(a)') 'set cntrparam levels 15'
        write(unit, '(a)') 'unset surface'
        write(unit, '(a)') ''
        
        write(filename, '(a,i0,a)') 'potential_state', s, '.dat'
        write(unit, '(3a)') 'splot "', trim(filename), '" u 1:2:3 w l notitle'
        write(unit, '(a)') ''
      end do
      
      write(unit, '(a)') 'unset multiplot'
      close(unit)
    end do
    
  contains
    function int_to_str(i) result(str)
      integer, intent(in) :: i
      character(len=20) :: str
      write(str, '(i0)') i
      str = adjustl(str)
    end function int_to_str
  end subroutine generate_contour_plots

  ! ============================================================================
  ! GENERATE SURFACE PLOTS
  ! ============================================================================
  subroutine generate_surface_plots(this)
    class(potential_2d_analyzer), intent(in) :: this
    
    character(len=200) :: scriptname, dirname, outfile
    integer :: unit, i, m1, m2, ios
    integer(ip) :: nstate_local
    
    if (.not. this%enabled) return
    
    nstate_local = this%nstate
    
    do i = 1, this%npairs_to_analyze
      m1 = this%mode_pairs(1,i)
      m2 = this%mode_pairs(2,i)
      
      write(dirname, '(a,a,i0,a,i0,a)') trim(this%output_dir), "mode", m1, "_", m2, "/"
      write(scriptname, '(a,a)') trim(dirname), "plot_surface.gnu"
      write(outfile, '(a,a)') trim(dirname), "potential_surface.png"
      
      open(newunit=unit, file=trim(scriptname), status='replace', iostat=ios)
      if (ios /= 0) cycle
      
      write(unit, '(a)') '# Gnuplot script for 3D potential surface plots'
      write(unit, '(a)') 'set terminal pngcairo enhanced size 1200,800'
      write(unit, '(a)') 'set output "' // trim(outfile) // '"'
      write(unit, '(a)') 'set title "3D Potential Surface - Mode ' // &
                         trim(int_to_str(m1)) // '-' // trim(int_to_str(m2)) // '"'
      write(unit, '(a)') 'set xlabel "Mode ' // trim(int_to_str(m1)) // '"'
      write(unit, '(a)') 'set ylabel "Mode ' // trim(int_to_str(m2)) // '"'
      write(unit, '(a)') 'set zlabel "Potential (eV)"'
      write(unit, '(a)') 'set pm3d'
      write(unit, '(a)') 'set hidden3d'
      write(unit, '(a)') 'set palette rgb 21,22,23'
      
      if (nstate_local >= 1) then
        write(unit, '(a)') 'splot "potential_state1.dat" u 1:2:3 w pm3d title "State 1"'
      end if
      
      close(unit)
    end do
    
  contains
    function int_to_str(i) result(str)
      integer, intent(in) :: i
      character(len=20) :: str
      write(str, '(i0)') i
      str = adjustl(str)
    end function int_to_str
  end subroutine generate_surface_plots

  ! ============================================================================
  ! CLEANUP
  ! ============================================================================
  subroutine cleanup_2d_potential_analyzer(this)
    class(potential_2d_analyzer), intent(inout) :: this
    
    if (allocated(this%potential_2d)) deallocate(this%potential_2d)
    if (allocated(this%coord_x)) deallocate(this%coord_x)
    if (allocated(this%coord_y)) deallocate(this%coord_y)
    if (allocated(this%min_potentials)) deallocate(this%min_potentials)
    if (allocated(this%max_potentials)) deallocate(this%max_potentials)
    if (allocated(this%saddle_points)) deallocate(this%saddle_points)
    if (allocated(this%mode_pairs)) deallocate(this%mode_pairs)
    
    if (this%unit_summary > 0) then
      close(this%unit_summary)
      this%unit_summary = -1
    end if
    
    if (this%verbose .and. this%enabled) then
      print '(a,i0,a)', "2D potential cuts analyzer cleaned up. Generated ", &
            this%cut_counter, " potential surfaces."
    end if
  end subroutine cleanup_2d_potential_analyzer

  ! ============================================================================
  ! GETTERS AND SETTERS
  ! ============================================================================
  function is_enabled(this) result(enabled)
    class(potential_2d_analyzer), intent(in) :: this
    logical :: enabled
    enabled = this%enabled
  end function is_enabled
  
  subroutine enable(this)
    class(potential_2d_analyzer), intent(inout) :: this
    this%enabled = .true.
    if (this%verbose) print '(a)', "2D potential cuts analyzer enabled."
  end subroutine enable
  
  subroutine disable(this)
    class(potential_2d_analyzer), intent(inout) :: this
    this%enabled = .false.
    if (this%verbose) print '(a)', "2D potential cuts analyzer disabled."
  end subroutine disable
  
  subroutine set_output_frequency(this, freq)
    class(potential_2d_analyzer), intent(inout) :: this
    integer, intent(in) :: freq
    this%output_frequency = freq
    if (this%verbose .and. this%enabled) then
      print '(a,i0)', "Output frequency set to every ", freq, " steps"
    end if
  end subroutine set_output_frequency
  
  function get_critical_points(this) result(points)
    class(potential_2d_analyzer), intent(in) :: this
    real(dp), allocatable :: points(:,:,:)
    
    if (allocated(this%saddle_points)) then
      allocate(points(size(this%saddle_points,1), &
                     size(this%saddle_points,2), &
                     size(this%saddle_points,3)))
      points = this%saddle_points
    else
      allocate(points(0,0,0))
    end if
  end function get_critical_points

 end module potential_2d_cuts_module
