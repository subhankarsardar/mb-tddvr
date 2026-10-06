 ! ============================================================================
! 1D POTENTIAL CUTS MODULE
! Creates 1D potential cuts along specified modes for all electronic surfaces
! Compatible with SO-TDDVR.f90
! Based on: Journal of Molecular Structure 1110 (2016) 32-43
! ============================================================================

 module potential_cuts_module
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  
  private
  
  ! Public types and procedures
  public :: potential_cuts_analyzer, dp, ip
  public :: create_potential_analyzer
  
  ! Kind parameters
  integer, parameter :: dp = real64
  integer, parameter :: ip = int32
  
  ! Default modes to analyze
  integer, parameter :: DEFAULT_NMODES = 6
  integer, parameter :: DEFAULT_MODES(6) = [1, 2, 5, 7, 15, 18]
  
  ! Plotting parameters
  integer, parameter :: N_POINTS = 100
  real(dp), parameter :: CUT_RANGE = 5.0_dp
  real(dp), parameter :: SMALL = 1.0e-12_dp
  
  ! Type for potential cuts analyzer
  type :: potential_cuts_analyzer
    private
    logical :: enabled = .false.
    logical :: verbose = .true.
    logical :: parallel = .true.
    
    ! Modes to analyze
    integer(ip) :: nmodes_to_analyze = DEFAULT_NMODES
    integer(ip), allocatable :: modes_to_analyze(:)
    
    ! Output control
    integer :: output_frequency = 50
    character(len=100) :: output_dir = "./potential_cuts/"
    
    ! Internal storage
    real(dp), allocatable :: potential_cuts(:,:,:)
    real(dp), allocatable :: coordinates(:,:)
    real(dp), allocatable :: times(:)
    real(dp), allocatable :: min_potentials(:,:)
    real(dp), allocatable :: max_potentials(:,:)
    
    ! File units
    integer :: unit_summary = -1
    
    ! Counters
    integer :: cut_counter = 0
    integer :: max_time_steps = 10000
    integer(ip) :: nstate = 0
    integer(ip) :: nmode = 0
    
  contains
    procedure :: initialize => init_potential_analyzer
    procedure :: set_modes_to_analyze
    procedure :: compute_1d_potential_cuts
    procedure :: find_minima_maxima
    procedure :: save_potential_data
    procedure :: save_summary_data
    procedure :: generate_gnuplot_scripts
    procedure :: finalize => cleanup_potential_analyzer
    procedure :: is_enabled
    procedure :: enable
    procedure :: disable
    procedure :: set_output_frequency
    procedure :: get_min_potentials
    procedure :: get_max_potentials
  end type potential_cuts_analyzer
  
  ! Constructor interface
  interface potential_cuts_analyzer
    module procedure create_potential_analyzer
  end interface

 contains

  ! ============================================================================
  ! CONSTRUCTOR
  ! ============================================================================
  function create_potential_analyzer(enabled, verbose, parallel, &
                                    modes_list, output_freq) result(analyzer)
    logical, intent(in), optional :: enabled, verbose, parallel
    integer(ip), intent(in), optional :: modes_list(:)
    integer, intent(in), optional :: output_freq
    type(potential_cuts_analyzer) :: analyzer
    
    if (present(enabled)) analyzer%enabled = enabled
    if (present(verbose)) analyzer%verbose = verbose
    if (present(parallel)) analyzer%parallel = parallel
    if (present(output_freq)) analyzer%output_frequency = output_freq
    
    if (present(modes_list)) then
      analyzer%nmodes_to_analyze = size(modes_list)
      allocate(analyzer%modes_to_analyze(analyzer%nmodes_to_analyze))
      analyzer%modes_to_analyze = modes_list
    else
      allocate(analyzer%modes_to_analyze(DEFAULT_NMODES))
      analyzer%modes_to_analyze = DEFAULT_MODES
    end if
    
    if (analyzer%verbose .and. analyzer%enabled) then
      print '(a)', "Potential cuts analyzer created (enabled)"
      print '(a,i0)', "  Modes to analyze: ", analyzer%nmodes_to_analyze
      print '(a,*(i4))', "  Mode list: ", analyzer%modes_to_analyze
      print '(a,i0)', "  Output frequency: every ", analyzer%output_frequency, " steps"
      if (analyzer%parallel) print '(a)', "  Parallel execution: enabled"
    end if
  end function create_potential_analyzer

  ! ============================================================================
  ! INITIALIZE POTENTIAL ANALYZER
  ! ============================================================================
  subroutine init_potential_analyzer(this, nstate, nmode, freqn, ffn)
    class(potential_cuts_analyzer), intent(inout) :: this
    integer(ip), intent(in) :: nstate, nmode
    real(dp), intent(in) :: freqn(:), ffn(:)
    
    integer :: ios, i, m, j
    character(len=200) :: dir_command
    real(dp) :: dx
    
    if (.not. this%enabled) return
    
    this%nstate = nstate
    this%nmode = nmode
    
    if (this%verbose) then
      print '(a)', "Initializing potential cuts analyzer..."
      print '(a,i0)', "  Number of states: ", nstate
      print '(a,i0)', "  Total modes: ", nmode
    end if
    
    ! Check mode indices are valid
    do i = 1, this%nmodes_to_analyze
      m = this%modes_to_analyze(i)
      if (m < 1 .or. m > nmode) then
        print '(a,i0,a)', "ERROR: Mode ", m, " is out of bounds!"
        this%enabled = .false.
        return
      end if
    end do
    
    ! Allocate storage
    allocate(this%potential_cuts(nstate, N_POINTS, this%nmodes_to_analyze))
    allocate(this%coordinates(N_POINTS, this%nmodes_to_analyze))
    allocate(this%min_potentials(nstate, this%nmodes_to_analyze))
    allocate(this%max_potentials(nstate, this%nmodes_to_analyze))
    allocate(this%times(this%max_time_steps))
    
    ! Initialize coordinate grid
    dx = 2.0_dp * CUT_RANGE / real(N_POINTS - 1, dp)
    do i = 1, this%nmodes_to_analyze
      m = this%modes_to_analyze(i)
      do j = 1, N_POINTS
        this%coordinates(j, i) = -CUT_RANGE + (j-1) * dx
        if (freqn(m) > SMALL) then
          this%coordinates(j, i) = this%coordinates(j, i) / sqrt(2.0_dp * freqn(m))
        end if
      end do
    end do
    
    this%potential_cuts = 0.0_dp
    this%min_potentials = huge(1.0_dp)
    this%max_potentials = -huge(1.0_dp)
    this%times = 0.0_dp
    this%cut_counter = 0
    
    dir_command = "mkdir -p " // trim(this%output_dir)
    call system(dir_command)
    
    open(newunit=this%unit_summary, &
         file=trim(this%output_dir)//"potential_summary.dat", &
         status='replace', iostat=ios)
    if (ios == 0) then
      write(this%unit_summary, '(a)') '# Potential Minima and Maxima Summary'
      write(this%unit_summary, '(a)') '# Time(fs), Mode, State, Min_Potential(eV), Max_Potential(eV), Min_Coord, Max_Coord'
    end if
    
    if (this%verbose) then
      print '(a)', "Potential cuts analyzer initialized successfully."
      print '(a,a)', "  Output directory: ", trim(this%output_dir)
    end if
  end subroutine init_potential_analyzer

  ! ============================================================================
  ! SET MODES TO ANALYZE
  ! ============================================================================
  subroutine set_modes_to_analyze(this, modes_list)
    class(potential_cuts_analyzer), intent(inout) :: this
    integer(ip), intent(in) :: modes_list(:)
    
    if (allocated(this%modes_to_analyze)) deallocate(this%modes_to_analyze)
    
    this%nmodes_to_analyze = size(modes_list)
    allocate(this%modes_to_analyze(this%nmodes_to_analyze))
    this%modes_to_analyze = modes_list
    
    if (this%verbose .and. this%enabled) then
      print '(a,i0)', "Updated modes to analyze: ", this%nmodes_to_analyze
      print '(a,*(i4))', "  New mode list: ", this%modes_to_analyze
    end if
  end subroutine set_modes_to_analyze

  ! ============================================================================
  ! COMPUTE 1D POTENTIAL CUTS
  ! ============================================================================
  subroutine compute_1d_potential_cuts(this, time, step, potential_func, &
                                      freqn, ffn, center_coords)
    class(potential_cuts_analyzer), intent(inout) :: this
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
    
    integer(ip) :: i, j, k, m, nstate_local, nmode_local
    real(dp), allocatable :: coords(:)
    real(dp), allocatable :: v_diag(:,:)
    real(dp) :: v_min, v_max
    integer :: time_idx
    
    if (.not. this%enabled) return
    if (mod(step, this%output_frequency) /= 0) return
    
    nstate_local = this%nstate
    nmode_local = this%nmode
    
    allocate(coords(nmode_local))
    allocate(v_diag(nstate_local, nstate_local))
    
    this%cut_counter = this%cut_counter + 1
    time_idx = this%cut_counter
    
    if (time_idx <= size(this%times)) then
      this%times(time_idx) = time
    end if
    
    if (this%verbose .and. mod(step, 100) == 0) then
      print '(a,f8.3,a,i0)', "Computing 1D potential cuts at t=", time, &
            " fs (step ", step, ")"
    end if
    
    do i = 1, this%nmodes_to_analyze
      m = this%modes_to_analyze(i)
      
      do j = 1, N_POINTS
        coords = center_coords
        coords(m) = this%coordinates(j, i)
        v_diag = potential_func(coords, nstate_local)
        
        do k = 1, nstate_local
          this%potential_cuts(k, j, i) = v_diag(k, k)
        end do
      end do
      
      do k = 1, nstate_local
        v_min = minval(this%potential_cuts(k, :, i))
        v_max = maxval(this%potential_cuts(k, :, i))
        
        if (v_min < this%min_potentials(k, i)) then
          this%min_potentials(k, i) = v_min
        end if
        if (v_max > this%max_potentials(k, i)) then
          this%max_potentials(k, i) = v_max
        end if
      end do
    end do
    
    deallocate(coords, v_diag)
    
    call this%save_potential_data(time, step)
    call this%save_summary_data(time, step)
  end subroutine compute_1d_potential_cuts

  ! ============================================================================
  ! FIND POTENTIAL MINIMA AND MAXIMA
  ! ============================================================================
  subroutine find_minima_maxima(this)
    class(potential_cuts_analyzer), intent(inout) :: this
    
    integer(ip) :: i, k, nstate_local, min_idx, max_idx
    real(dp) :: v_min, v_max
    
    if (.not. this%enabled) return
    
    nstate_local = this%nstate
    
    do i = 1, this%nmodes_to_analyze
      do k = 1, nstate_local
        v_min = minval(this%potential_cuts(k, :, i))
        v_max = maxval(this%potential_cuts(k, :, i))
        
        min_idx = minloc(this%potential_cuts(k, :, i), dim=1)
        max_idx = maxloc(this%potential_cuts(k, :, i), dim=1)
        
        if (v_min < this%min_potentials(k, i)) then
          this%min_potentials(k, i) = v_min
        end if
        if (v_max > this%max_potentials(k, i)) then
          this%max_potentials(k, i) = v_max
        end if
        
        if (this%verbose) then
          print '(a,i0,a,i0,a,f12.6,a,f8.4)', &
                "Mode ", this%modes_to_analyze(i), &
                ", State ", k, &
                ": Min = ", v_min, &
                " at x = ", this%coordinates(min_idx, i)
          print '(a,f12.6,a,f8.4)', &
                "                  Max = ", v_max, &
                " at x = ", this%coordinates(max_idx, i)
        end if
      end do
    end do
  end subroutine find_minima_maxima

  ! ============================================================================
  ! SAVE POTENTIAL DATA
  ! ============================================================================
  subroutine save_potential_data(this, time, step)
    class(potential_cuts_analyzer), intent(in) :: this
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    
    character(len=200) :: filename
    integer :: i, j, k, unit_file, m
    integer(ip) :: nstate_local
    character(len=30) :: fmt_str
    
    if (.not. this%enabled) return
    
    nstate_local = this%nstate
    write(fmt_str, '(a,i0,a)') '(f12.6,', nstate_local, 'f16.8)'
    
    do i = 1, this%nmodes_to_analyze
      m = this%modes_to_analyze(i)
      
      write(filename, '(a,i3.3,a,f8.3,a)') &
            trim(this%output_dir) // "potential_mode", m, &
            "_t", time, ".dat"
      
      open(newunit=unit_file, file=trim(filename), status='replace')
      
      write(unit_file, '(a)') '# 1D Potential Cuts'
      write(unit_file, '(a,f12.6,a,i0)') '# Time: ', time, ' fs, Step: ', step
      write(unit_file, '(a,i0)') '# Mode: ', m
      write(unit_file, '(a,*(i12))') '# Coordinate', (k, k=1, nstate_local)
      
      do j = 1, N_POINTS
        write(unit_file, fmt_str) &
              this%coordinates(j, i), &
              (this%potential_cuts(k, j, i), k=1, nstate_local)
      end do
      
      close(unit_file)
    end do
    
    ! Combined file
    write(filename, '(a,f8.3,a)') &
          trim(this%output_dir) // "potential_all_modes_t", time, ".dat"
    
    open(newunit=unit_file, file=trim(filename), status='replace')
    write(unit_file, '(a)') '# Combined 1D Potential Cuts for All Modes'
    write(unit_file, '(a,f12.6)') '# Time: ', time
    write(unit_file, '(a,i0)') '# Number of states: ', nstate_local
    write(unit_file, '(a)') '# Format: Mode, Coordinate, State1, State2, ...'
    
    do i = 1, this%nmodes_to_analyze
      m = this%modes_to_analyze(i)
      write(unit_file, '(a,i0)') '# Mode ', m
      do j = 1, N_POINTS
        write(unit_file, '(i4,f12.6,6f16.8)') &
              m, this%coordinates(j, i), &
              (this%potential_cuts(k, j, i), k=1, min(nstate_local, 6))
      end do
      write(unit_file, *)
    end do
    
    close(unit_file)
  end subroutine save_potential_data

  ! ============================================================================
  ! SAVE SUMMARY DATA
  ! ============================================================================
  subroutine save_summary_data(this, time, step)
    class(potential_cuts_analyzer), intent(in) :: this
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    
    integer(ip) :: i, k, m, nstate_local, min_idx, max_idx
    real(dp) :: min_coord, max_coord
    
    if (.not. this%enabled) return
    if (this%unit_summary < 0) return
    
    nstate_local = this%nstate
    
    do i = 1, this%nmodes_to_analyze
      m = this%modes_to_analyze(i)
      do k = 1, nstate_local
        min_idx = minloc(this%potential_cuts(k, :, i), dim=1)
        max_idx = maxloc(this%potential_cuts(k, :, i), dim=1)
        min_coord = this%coordinates(min_idx, i)
        max_coord = this%coordinates(max_idx, i)
        
        write(this%unit_summary, '(f12.6,2i6,4f16.8)') &
              time, m, k, &
              this%min_potentials(k, i), &
              this%max_potentials(k, i), &
              min_coord, max_coord
      end do
    end do
    
    flush(this%unit_summary)
  end subroutine save_summary_data

  ! ============================================================================
  ! GENERATE GNUPLOT SCRIPTS
  ! ============================================================================
  subroutine generate_gnuplot_scripts(this)
    class(potential_cuts_analyzer), intent(in) :: this
    
    character(len=200) :: scriptname
    integer :: unit, i, m, idx, nstate_local
    character(len=20) :: state_str
    
    if (.not. this%enabled) return
    
    nstate_local = this%nstate
    
    do i = 1, this%nmodes_to_analyze
      m = this%modes_to_analyze(i)
      
      scriptname = trim(this%output_dir) // "plot_potential_mode" // &
                   trim(int_to_str(m)) // ".gnu"
      
      open(newunit=unit, file=trim(scriptname), status='replace')
      
      write(unit, '(a)') '# Gnuplot script for 1D potential cuts'
      write(unit, '(a)') 'set terminal pngcairo enhanced size 1200,800'
      write(unit, '(a,i0,a)') 'set output "' // trim(this%output_dir) // 'potential_mode', m, '.png"'
      write(unit, '(a)') ''
      write(unit, '(a)') '# Plot settings'
      write(unit, '(a,i0,a)') 'set title "1D Potential Cuts - Mode ', m, '"'
      write(unit, '(a)') 'set xlabel "Coordinate"'
      write(unit, '(a)') 'set ylabel "Potential Energy (eV)"'
      write(unit, '(a)') 'set grid'
      write(unit, '(a)') 'set key outside right top'
      write(unit, '(a)') ''
      write(unit, '(a)') '# Define line styles for different states'
      write(unit, '(a)') 'set style line 1 lc rgb "#FF0000" lw 2'
      write(unit, '(a)') 'set style line 2 lc rgb "#00FF00" lw 2'
      write(unit, '(a)') 'set style line 3 lc rgb "#0000FF" lw 2'
      write(unit, '(a)') 'set style line 4 lc rgb "#FF00FF" lw 2'
      write(unit, '(a)') 'set style line 5 lc rgb "#00FFFF" lw 2'
      write(unit, '(a)') 'set style line 6 lc rgb "#FFA500" lw 2'
      write(unit, '(a)') ''
      
      write(unit, '(a)') '# Plot command'
      write(unit, '(a)') 'plot \'
      
      do idx = 1, nstate_local-1
        write(state_str, '(i0)') idx
        write(unit, '(a,i0,a,i0,a)') '  "' // trim(this%output_dir) // 'potential_mode', m, '.dat" u 1:', &
              idx+1, ' w l ls ', idx, ' t "State ' // trim(state_str) // '", \'
      end do
      
      write(state_str, '(i0)') nstate_local
      write(unit, '(a,i0,a,i0,a)') '  "' // trim(this%output_dir) // 'potential_mode', m, '.dat" u 1:', &
            nstate_local+1, ' w l ls ', nstate_local, ' t "State ' // trim(state_str) // '"'
      
      close(unit)
    end do
    
    ! Script for comparing all modes
    scriptname = trim(this%output_dir) // "plot_all_modes_comparison.gnu"
    open(newunit=unit, file=trim(scriptname), status='replace')
    
    write(unit, '(a)') '# Gnuplot script for comparing all modes'
    write(unit, '(a)') 'set terminal pngcairo enhanced size 1600,1200'
    write(unit, '(a)') 'set output "' // trim(this%output_dir) // 'potential_all_modes_comparison.png"'
    
    if (this%nmodes_to_analyze <= 2) then
      write(unit, '(a)') 'set multiplot layout 1,' // trim(int_to_str(this%nmodes_to_analyze))
    else if (this%nmodes_to_analyze <= 4) then
      write(unit, '(a)') 'set multiplot layout 2,' // trim(int_to_str((this%nmodes_to_analyze+1)/2))
    else
      write(unit, '(a)') 'set multiplot layout 3,3'
    end if
    write(unit, '(a)') ''
    
    do i = 1, this%nmodes_to_analyze
      m = this%modes_to_analyze(i)
      write(unit, '(a,i0,a)') 'set title "Mode ', m, '"'
      write(unit, '(a)') 'set xlabel "Coordinate"'
      write(unit, '(a)') 'set ylabel "Potential (eV)"'
      write(unit, '(a)') 'set grid'
      
      write(unit, '(a,i0,a)') 'plot "' // trim(this%output_dir) // 'potential_mode', m, '.dat" u 1:2 w l t "State 1"'
      
      do idx = 2, min(nstate_local, 6)
        write(state_str, '(i0)') idx
        write(unit, '(a,i0,a,i0,a)') '     "' // trim(this%output_dir) // 'potential_mode', m, '.dat" u 1:', &
              idx+1, ' w l t "State ' // trim(state_str) // '"'
      end do
      write(unit, '(a)') ''
    end do
    
    write(unit, '(a)') 'unset multiplot'
    close(unit)
    
    if (this%verbose) then
      print '(a,a)', "Gnuplot scripts generated in ", trim(this%output_dir)
    end if
    
  contains
    function int_to_str(i) result(str)
      integer, intent(in) :: i
      character(len=20) :: str
      write(str, '(i0)') i
      str = adjustl(str)
    end function int_to_str
  end subroutine generate_gnuplot_scripts

  ! ============================================================================
  ! CLEANUP
  ! ============================================================================
  subroutine cleanup_potential_analyzer(this)
    class(potential_cuts_analyzer), intent(inout) :: this
    
    if (allocated(this%potential_cuts)) deallocate(this%potential_cuts)
    if (allocated(this%coordinates)) deallocate(this%coordinates)
    if (allocated(this%times)) deallocate(this%times)
    if (allocated(this%min_potentials)) deallocate(this%min_potentials)
    if (allocated(this%max_potentials)) deallocate(this%max_potentials)
    if (allocated(this%modes_to_analyze)) deallocate(this%modes_to_analyze)
    
    if (this%unit_summary > 0) then
      close(this%unit_summary)
      this%unit_summary = -1
    end if
    
    if (this%verbose .and. this%enabled) then
      print '(a,i0,a)', "Potential cuts analyzer cleaned up. Generated ", &
            this%cut_counter, " potential cuts."
    end if
  end subroutine cleanup_potential_analyzer

  ! ============================================================================
  ! GETTERS AND SETTERS
  ! ============================================================================
  function is_enabled(this) result(enabled)
    class(potential_cuts_analyzer), intent(in) :: this
    logical :: enabled
    enabled = this%enabled
  end function is_enabled
  
  subroutine enable(this)
    class(potential_cuts_analyzer), intent(inout) :: this
    this%enabled = .true.
    if (this%verbose) print '(a)', "Potential cuts analyzer enabled."
  end subroutine enable
  
  subroutine disable(this)
    class(potential_cuts_analyzer), intent(inout) :: this
    this%enabled = .false.
    if (this%verbose) print '(a)', "Potential cuts analyzer disabled."
  end subroutine disable
  
  subroutine set_output_frequency(this, freq)
    class(potential_cuts_analyzer), intent(inout) :: this
    integer, intent(in) :: freq
    this%output_frequency = freq
    if (this%verbose .and. this%enabled) then
      print '(a,i0)', "Output frequency set to every ", freq, " steps"
    end if
  end subroutine set_output_frequency
  
  function get_min_potentials(this) result(min_pots)
    class(potential_cuts_analyzer), intent(in) :: this
    real(dp), allocatable :: min_pots(:,:)
    
    if (allocated(this%min_potentials)) then
      allocate(min_pots(size(this%min_potentials,1), size(this%min_potentials,2)))
      min_pots = this%min_potentials
    else
      allocate(min_pots(0,0))
    end if
  end function get_min_potentials
  
  function get_max_potentials(this) result(max_pots)
    class(potential_cuts_analyzer), intent(in) :: this
    real(dp), allocatable :: max_pots(:,:)
    
    if (allocated(this%max_potentials)) then
      allocate(max_pots(size(this%max_potentials,1), size(this%max_potentials,2)))
      max_pots = this%max_potentials
    else
      allocate(max_pots(0,0))
    end if
  end function get_max_potentials

 end module potential_cuts_module
