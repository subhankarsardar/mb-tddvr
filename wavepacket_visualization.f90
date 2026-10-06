 ! ============================================================================
! WAVEPACKET VISUALIZATION MODULE FOR TDDVR
! Calculates wavepacket probabilities for all electronic surfaces
! Compatible with SO-TDDVR.f90
! Based on: Journal of Molecular Structure 1110 (2016) 32-43
! ============================================================================

 module wavepacket_visualization
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  
  private
  
  ! Public types and procedures
  public :: wavepacket_visualizer, dp, ip
  public :: create_wavepacket_visualizer
  
  ! Kind parameters
  integer, parameter :: dp = real64
  integer, parameter :: ip = int32
  
  ! Type for wavepacket visualizer
  type :: wavepacket_visualizer
    private
    logical :: enabled = .false.
    logical :: verbose = .true.
    logical :: parallel = .true.
    
    ! System dimensions
    integer(ip) :: nstate = 0
    integer(ip) :: nmode = 0
    integer(ip) :: nn = 0
    integer(ip) :: nnprod = 0
    
    ! Grid dimensions array
    integer(ip), allocatable :: ndim(:)
    
    ! LL array for mode indices
    integer(ip), allocatable :: LL(:)
    
    ! Calculation parameters
    integer :: calculation_frequency = 20
    real(dp) :: time_interval = 1.0_dp
    character(len=100) :: output_dir = "./wavepacket_analysis/"
    
    ! Probability storage
    real(dp), allocatable :: state_probabilities(:,:)
    real(dp), allocatable :: total_probability(:)
    real(dp), allocatable :: probability_times(:)
    integer(ip) :: max_steps = 5000
    integer(ip) :: current_step = 0
    
    ! Mode-specific probability tracking
    logical :: track_mode_probabilities = .false.
    real(dp), allocatable :: mode_probabilities(:,:,:)
    
    ! For 2D slices
    integer(ip) :: mode1 = 1, mode2 = 2
    integer(ip) :: mode3 = 5, mode4 = 15
    real(dp), allocatable :: densities(:,:,:)
    real(dp), allocatable :: grid1(:), grid2(:)
    
    ! File units
    integer :: unit_prob = -1
    integer :: unit_summary = -1
    
    ! State for contour data
    integer :: contour_counter = 0
    
  contains
    procedure :: initialize => init_visualizer
    procedure :: set_LL_array
    procedure :: calculate_probabilities
    procedure :: calculate_mode_probabilities
    procedure :: calculate_total_probability
    procedure :: save_probability_data
    procedure :: save_mode_probability_data
    procedure :: save_summary
    procedure :: finalize => cleanup_visualizer
    procedure :: is_enabled
    procedure :: enable
    procedure :: disable
    procedure :: set_calculation_frequency
    procedure :: set_track_mode_probabilities
    procedure :: get_state_probabilities
    procedure :: get_total_probability
    procedure :: get_probability_times
    procedure :: get_LL_array
    procedure :: visualize_current_step => calculate_current_step_probabilities
    procedure :: extract_2d_slice
    procedure :: save_contour_data
    procedure :: save_contour_data_with_time
  end type wavepacket_visualizer
  
  ! Constructor interface
  interface wavepacket_visualizer
    module procedure create_wavepacket_visualizer
  end interface

 contains

  ! ============================================================================
  ! CONSTRUCTOR
  ! ============================================================================
  function create_wavepacket_visualizer(enabled, verbose, parallel, &
                                        output_freq, mode_pair1, mode_pair2, &
                                        track_modes) result(viz)
    logical, intent(in), optional :: enabled, verbose, parallel, track_modes
    integer, intent(in), optional :: output_freq
    integer(ip), intent(in), optional :: mode_pair1(2), mode_pair2(2)
    type(wavepacket_visualizer) :: viz
    
    if (present(enabled)) viz%enabled = enabled
    if (present(verbose)) viz%verbose = verbose
    if (present(parallel)) viz%parallel = parallel
    if (present(output_freq)) viz%calculation_frequency = output_freq
    if (present(track_modes)) viz%track_mode_probabilities = track_modes
    
    if (present(mode_pair1)) then
      viz%mode1 = mode_pair1(1)
      viz%mode2 = mode_pair1(2)
    end if
    
    if (present(mode_pair2)) then
      viz%mode3 = mode_pair2(1)
      viz%mode4 = mode_pair2(2)
    end if
    
    call viz%set_LL_array()
    
    if (viz%verbose .and. viz%enabled) then
      print '(a)', "Wavepacket visualizer created (enabled)"
      print '(a,i0)', "  Calculation frequency: every ", viz%calculation_frequency, " steps"
      if (viz%parallel) print '(a)', "  Parallel execution: enabled"
    end if
  end function create_wavepacket_visualizer

  ! ============================================================================
  ! INITIALIZE VISUALIZER
  ! ============================================================================
  subroutine init_visualizer(this, nstate, ndim, nn, nmode)
    class(wavepacket_visualizer), intent(inout) :: this
    integer(ip), intent(in) :: nstate, nn, nmode
    integer(ip), intent(in) :: ndim(:)
    
    integer :: ios
    integer(ip) :: i
    character(len=200) :: dir_command
    
    if (.not. this%enabled) return
    
    this%nstate = nstate
    this%nmode = nmode
    this%nn = nn
    allocate(this%ndim(nmode))
    this%ndim = ndim
    
    this%nnprod = 1
    do i = 1, nmode
      this%nnprod = this%nnprod * ndim(i)
    end do
    
    if (.not. allocated(this%LL)) then
      call this%set_LL_array()
    end if
    
    if (this%verbose) then
      print '(a)', "Initializing wavepacket visualizer..."
      print '(a,i0)', "  Number of electronic states: ", nstate
      print '(a,i0)', "  Number of vibrational modes: ", nmode
      print '(a,i0)', "  Total grid points: ", this%nnprod
    end if
    
    this%max_steps = 5000
    
    allocate(this%state_probabilities(nstate, this%max_steps))
    allocate(this%total_probability(this%max_steps))
    allocate(this%probability_times(this%max_steps))
    
    this%state_probabilities = 0.0_dp
    this%total_probability = 0.0_dp
    this%probability_times = 0.0_dp
    this%current_step = 0
    this%contour_counter = 0
    
    ! Allocate 2D slice arrays if needed
    if (this%mode1 > 0 .and. this%mode2 > 0 .and. &
        this%mode1 <= nmode .and. this%mode2 <= nmode) then
      allocate(this%grid1(ndim(this%mode1)))
      allocate(this%grid2(ndim(this%mode2)))
      allocate(this%densities(nstate, ndim(this%mode1), ndim(this%mode2)))
      this%grid1 = 0.0_dp
      this%grid2 = 0.0_dp
      this%densities = 0.0_dp
    end if
    
    if (this%track_mode_probabilities) then
      allocate(this%mode_probabilities(nmode, nstate, this%max_steps))
      this%mode_probabilities = 0.0_dp
    end if
    
    dir_command = "mkdir -p " // trim(this%output_dir)
    call system(dir_command)
    
    if (this%verbose) then
      print '(a)', "Wavepacket visualizer initialized successfully."
      print '(a,a)', "  Output directory: ", trim(this%output_dir)
    end if
  end subroutine init_visualizer

  ! ============================================================================
  ! SET LL ARRAY
  ! ============================================================================
  subroutine set_LL_array(this, custom_LL)
    class(wavepacket_visualizer), intent(inout) :: this
    integer(ip), intent(in), optional :: custom_LL(:)
    
    if (present(custom_LL)) then
      if (allocated(this%LL)) deallocate(this%LL)
      allocate(this%LL(size(custom_LL)))
      this%LL = custom_LL
    else
      if (allocated(this%LL)) deallocate(this%LL)
      allocate(this%LL(43))
      this%LL = [1, 25, 46, 57, 92, 132, 178, 192, 210, 250, &
                305, 350, 400, 450, 500, 550, 600, 650, 701, 755, &
                800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, &
                1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, &
                1800, 1850, 0]
    end if
  end subroutine set_LL_array

  ! ============================================================================
  ! CALCULATE CURRENT STEP PROBABILITIES
  ! ============================================================================
  subroutine calculate_current_step_probabilities(this, time, step, wavefunction, tot, &
                                                 newgrid, ndim, nstate)
    class(wavepacket_visualizer), intent(inout) :: this
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    complex(dp), intent(in) :: wavefunction(:,:)
    integer(ip), intent(in) :: tot(:,:)
    real(dp), intent(in) :: newgrid(:,:)
    integer(ip), intent(in) :: ndim(:)
    integer(ip), intent(in) :: nstate
    
    integer(ip) :: istate
    
    if (.not. this%enabled) return
    
    if (mod(step, this%calculation_frequency) /= 0) return
    
    this%current_step = this%current_step + 1
    
    if (this%current_step > this%max_steps) then
      if (this%verbose) then
        print '(a,i0,a)', "WARNING: Maximum steps (", this%max_steps, ") reached!"
        print '(a)', "         Disabling further probability calculations."
      end if
      this%enabled = .false.
      return
    end if
    
    this%probability_times(this%current_step) = time
    
    ! Calculate probabilities
    call this%calculate_probabilities(time, step, wavefunction)
    
    ! Calculate mode probabilities if enabled
    if (this%track_mode_probabilities) then
      call this%calculate_mode_probabilities(wavefunction, tot, this%current_step)
    end if
    
    ! Calculate 2D slices if requested
    if (allocated(this%densities)) then
      call this%extract_2d_slice(wavefunction, tot, newgrid, ndim, &
                                this%mode1, this%mode2, &
                                this%densities, this%grid1, this%grid2)
      
      do istate = 1, nstate
        call this%save_contour_data_with_time(time, step, &
                                   this%densities(istate,:,:), &
                                   this%grid1, this%grid2, istate)
      end do
    end if
    
    call this%calculate_total_probability(this%current_step)
    
    if (mod(this%current_step, 50) == 0) then
      call this%save_probability_data()
    end if
    
    if (this%verbose .and. mod(step, 100) == 0) then
      print '(a,f8.3,a,i0,a)', "Wavepacket probabilities at t=", time, &
            " fs (step ", step, ")"
      print '(a,6es12.4)', "  State probabilities: ", &
            this%state_probabilities(1:min(6,nstate), this%current_step)
      print '(a,es12.4)', "  Total probability: ", &
            this%total_probability(this%current_step)
    end if
    
  end subroutine calculate_current_step_probabilities

  ! ============================================================================
  ! CALCULATE STATE PROBABILITIES
  ! ============================================================================
  subroutine calculate_probabilities(this, time, step, wavefunction)
    class(wavepacket_visualizer), intent(inout) :: this
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    complex(dp), intent(in) :: wavefunction(:,:)
    
    integer(ip) :: i, k
    real(dp) :: norm_sq
    
    do i = 1, this%nstate
      norm_sq = 0.0_dp
      do k = 1, this%nnprod
        norm_sq = norm_sq + abs(wavefunction(i, k))**2
      end do
      this%state_probabilities(i, this%current_step) = norm_sq
    end do
    
  end subroutine calculate_probabilities

  ! ============================================================================
  ! EXTRACT 2D SLICE
  ! ============================================================================
  subroutine extract_2d_slice(this, wavefunction, tot, newgrid, ndim, &
                             mode_i, mode_j, densities, grid_i, grid_j)
    class(wavepacket_visualizer), intent(in) :: this
    complex(dp), intent(in) :: wavefunction(:,:)
    integer(ip), intent(in) :: tot(:,:)
    real(dp), intent(in) :: newgrid(:,:)
    integer(ip), intent(in) :: ndim(:)
    integer(ip), intent(in) :: mode_i, mode_j
    real(dp), intent(out) :: densities(:,:,:)
    real(dp), intent(out) :: grid_i(:), grid_j(:)
    
    integer(ip) :: nstate, nnprod, nmode
    integer(ip) :: i, j, k, l, idx_i, idx_j
    
    nstate = size(wavefunction, 1)
    nnprod = size(wavefunction, 2)
    nmode = size(tot, 2)
    
    densities = 0.0_dp
    grid_i = 0.0_dp
    grid_j = 0.0_dp
    
    do i = 1, ndim(mode_i)
      grid_i(i) = newgrid(i, mode_i)
    end do
    
    do j = 1, ndim(mode_j)
      grid_j(j) = newgrid(j, mode_j)
    end do
    
    do k = 1, nnprod
      idx_i = tot(k, mode_i)
      idx_j = tot(k, mode_j)
      
      do l = 1, nstate
        densities(l, idx_i, idx_j) = densities(l, idx_i, idx_j) + &
                                    abs(wavefunction(l, k))**2
      end do
    end do
    
  end subroutine extract_2d_slice

  ! ============================================================================
  ! SAVE CONTOUR DATA - ORIGINAL VERSION
  ! ============================================================================
  subroutine save_contour_data(this, densities, grid_i, grid_j, state)
    class(wavepacket_visualizer), intent(inout) :: this
    real(dp), intent(in) :: densities(:,:)
    real(dp), intent(in) :: grid_i(:), grid_j(:)
    integer(ip), intent(in) :: state
    
    character(len=200) :: filename
    integer :: i, j, n_i, n_j, unit
    
    n_i = size(grid_i)
    n_j = size(grid_j)
    
    this%contour_counter = this%contour_counter + 1
    
    write(filename, '(a,i3.3,a,i3.3,a,i1,a,i6.6,a)') &
          trim(this%output_dir) // "wavepacket_mode", this%mode1, &
          "_", this%mode2, "_state", state, "_", this%contour_counter, ".dat"
    
    open(newunit=unit, file=trim(filename), status='replace')
    write(unit, '(a)') '# 2D Wavepacket Density'
    write(unit, '(a,i0,a,i0,a,i0)') '# Modes: ', this%mode1, ' - ', this%mode2, &
          ', State: ', state
    write(unit, '(a)') '# Q1, Q2, Density'
    
    do i = 1, n_i
      do j = 1, n_j
        write(unit, '(3f16.8)') grid_i(i), grid_j(j), densities(i,j)
      end do
      write(unit, *)
    end do
    
    close(unit)
    
  end subroutine save_contour_data

  ! ============================================================================
  ! SAVE CONTOUR DATA WITH TIME
  ! ============================================================================
  subroutine save_contour_data_with_time(this, time, step, densities, &
                                        grid_i, grid_j, state)
    class(wavepacket_visualizer), intent(inout) :: this
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    real(dp), intent(in) :: densities(:,:)
    real(dp), intent(in) :: grid_i(:), grid_j(:)
    integer(ip), intent(in) :: state
    
    character(len=200) :: filename
    integer :: i, j, n_i, n_j, unit
    
    n_i = size(grid_i)
    n_j = size(grid_j)
    
    write(filename, '(a,i3.3,a,i3.3,a,i1,a,f8.3,a)') &
          trim(this%output_dir) // "wavepacket_mode", this%mode1, &
          "_", this%mode2, "_state", state, "_t", time, ".dat"
    
    open(newunit=unit, file=trim(filename), status='replace')
    write(unit, '(a)') '# 2D Wavepacket Density'
    write(unit, '(a,f12.6,a,i0)') '# Time: ', time, ' fs, Step: ', step
    write(unit, '(a,i0,a,i0,a,i0)') '# Modes: ', this%mode1, ' - ', this%mode2, &
          ', State: ', state
    write(unit, '(a)') '# Q1, Q2, Density'
    
    do i = 1, n_i
      do j = 1, n_j
        write(unit, '(3f16.8)') grid_i(i), grid_j(j), densities(i,j)
      end do
      write(unit, *)
    end do
    
    close(unit)
    
  end subroutine save_contour_data_with_time

  ! ============================================================================
  ! CALCULATE MODE-SPECIFIC PROBABILITIES
  ! ============================================================================
  subroutine calculate_mode_probabilities(this, wavefunction, tot, step_idx)
    class(wavepacket_visualizer), intent(inout) :: this
    complex(dp), intent(in) :: wavefunction(:,:)
    integer(ip), intent(in) :: tot(:,:)
    integer(ip), intent(in) :: step_idx
    
    integer(ip) :: i, k, m
    real(dp) :: prob
    
    if (.not. allocated(this%mode_probabilities)) return
    
    this%mode_probabilities(:,:,step_idx) = 0.0_dp
    
    do i = 1, this%nstate
      do k = 1, this%nnprod
        prob = abs(wavefunction(i, k))**2
        do m = 1, this%nmode
          this%mode_probabilities(m, i, step_idx) = &
            this%mode_probabilities(m, i, step_idx) + prob
        end do
      end do
    end do
    
  end subroutine calculate_mode_probabilities

  ! ============================================================================
  ! CALCULATE TOTAL PROBABILITY
  ! ============================================================================
  subroutine calculate_total_probability(this, step_idx)
    class(wavepacket_visualizer), intent(inout) :: this
    integer(ip), intent(in) :: step_idx
    
    integer(ip) :: i
    
    this%total_probability(step_idx) = 0.0_dp
    do i = 1, this%nstate
      this%total_probability(step_idx) = this%total_probability(step_idx) + &
                                        this%state_probabilities(i, step_idx)
    end do
  end subroutine calculate_total_probability

  ! ============================================================================
  ! SAVE PROBABILITY DATA
  ! ============================================================================
  subroutine save_probability_data(this)
    class(wavepacket_visualizer), intent(in) :: this
    
    character(len=200) :: filename
    integer :: unit, i, j
    
    if (.not. this%enabled) return
    if (this%current_step == 0) return
    
    filename = trim(this%output_dir) // "state_probabilities.dat"
    
    open(newunit=unit, file=trim(filename), status='replace')
    
    write(unit, '(a)') '# State Probabilities vs Time'
    write(unit, '(a,i0)') '# Number of states: ', this%nstate
    write(unit, '(a,i0)') '# Number of vibrational modes: ', this%nmode
    write(unit, '(a,i0)') '# Number of time steps calculated: ', this%current_step
    write(unit, '(a)') '# Time (fs), State probabilities, Total probability'
    
    do i = 1, this%current_step
      write(unit, '(f12.6)', advance='no') this%probability_times(i)
      do j = 1, this%nstate
        write(unit, '(es16.8)', advance='no') this%state_probabilities(j, i)
      end do
      write(unit, '(es16.8)') this%total_probability(i)
    end do
    
    close(unit)
    
    ! Save individual state files
    do j = 1, this%nstate
      write(filename, '(a,i3.3,a)') trim(this%output_dir) // &
                                   "state_", j, "_probability.dat"
      open(newunit=unit, file=trim(filename), status='replace')
      
      write(unit, '(a,i0,a)') '# Probability for state ', j
      write(unit, '(a)') '# Time (fs), Probability'
      
      do i = 1, this%current_step
        write(unit, '(2es16.8)') this%probability_times(i), &
                                 this%state_probabilities(j, i)
      end do
      
      close(unit)
    end do
    
    if (this%track_mode_probabilities .and. allocated(this%mode_probabilities)) then
      call this%save_mode_probability_data()
    end if
    
    if (this%verbose) then
      print '(a,a)', "Wavepacket probability data saved to: ", trim(this%output_dir)
    end if
  end subroutine save_probability_data

  ! ============================================================================
  ! SAVE MODE PROBABILITY DATA
  ! ============================================================================
  subroutine save_mode_probability_data(this)
    class(wavepacket_visualizer), intent(in) :: this
    
    character(len=200) :: filename
    integer :: unit, i, j, m
    
    if (.not. allocated(this%mode_probabilities)) return
    
    do j = 1, this%nstate
      write(filename, '(a,i3.3,a)') trim(this%output_dir) // &
                                   "state_", j, "_mode_probabilities.dat"
      
      open(newunit=unit, file=trim(filename), status='replace')
      
      write(unit, '(a,i0,a)') '# Mode probabilities for state ', j
      write(unit, '(a,i0)') '# Number of modes: ', this%nmode
      write(unit, '(a,i0)') '# Number of time steps: ', this%current_step
      write(unit, '(a)') '# Time (fs), Mode probabilities'
      
      do i = 1, this%current_step
        write(unit, '(f12.6)', advance='no') this%probability_times(i)
        do m = 1, this%nmode
          write(unit, '(es16.8)', advance='no') this%mode_probabilities(m, j, i)
        end do
        write(unit, '()')
      end do
      
      close(unit)
    end do
    
  end subroutine save_mode_probability_data

  ! ============================================================================
  ! SAVE SUMMARY
  ! ============================================================================
  subroutine save_summary(this)
    class(wavepacket_visualizer), intent(in) :: this
    
    character(len=200) :: filename
    integer :: unit, i
    real(dp) :: avg_prob, max_prob, min_prob
    
    if (.not. this%enabled) return
    if (this%current_step == 0) return
    
    filename = trim(this%output_dir) // "wavepacket_summary.txt"
    
    open(newunit=unit, file=trim(filename), status='replace')
    
    write(unit, '(a)') '=============================================='
    write(unit, '(a)') 'WAVEPACKET PROBABILITY SUMMARY'
    write(unit, '(a)') '=============================================='
    write(unit, '(a,i0)') 'Number of electronic states: ', this%nstate
    write(unit, '(a,i0)') 'Number of vibrational modes: ', this%nmode
    write(unit, '(a,i0)') 'Time steps calculated: ', this%current_step
    write(unit, '(a)') ''
    write(unit, '(a)') 'STATE-WISE PROBABILITY STATISTICS:'
    write(unit, '(a)') '-----------------------------------'
    
    do i = 1, this%nstate
      avg_prob = sum(this%state_probabilities(i, 1:this%current_step)) / &
                real(this%current_step, dp)
      max_prob = maxval(this%state_probabilities(i, 1:this%current_step))
      min_prob = minval(this%state_probabilities(i, 1:this%current_step))
      
      write(unit, '(a,i3,a)') 'State ', i, ':'
      write(unit, '(a,es16.8)') '  Average: ', avg_prob
      write(unit, '(a,es16.8)') '  Maximum: ', max_prob
      write(unit, '(a,es16.8)') '  Minimum: ', min_prob
      write(unit, '(a)') ''
    end do
    
    ! Total probability statistics
    avg_prob = sum(this%total_probability(1:this%current_step)) / &
              real(this%current_step, dp)
    max_prob = maxval(this%total_probability(1:this%current_step))
    min_prob = minval(this%total_probability(1:this%current_step))
    
    write(unit, '(a)') 'TOTAL PROBABILITY STATISTICS:'
    write(unit, '(a)') '------------------------------'
    write(unit, '(a,es16.8)') 'Average: ', avg_prob
    write(unit, '(a,es16.8)') 'Maximum: ', max_prob
    write(unit, '(a,es16.8)') 'Minimum: ', min_prob
    write(unit, '(a,es16.8)') 'Deviation from 1.0: ', abs(avg_prob - 1.0_dp)
    
    close(unit)
    
    if (this%verbose) then
      print '(a,a)', "Wavepacket summary saved to: ", trim(filename)
    end if
  end subroutine save_summary

  ! ============================================================================
  ! CLEANUP
  ! ============================================================================
  subroutine cleanup_visualizer(this)
    class(wavepacket_visualizer), intent(inout) :: this
    
    if (this%enabled .and. this%current_step > 0) then
      call this%save_probability_data()
      call this%save_summary()
    end if
    
    if (allocated(this%state_probabilities)) deallocate(this%state_probabilities)
    if (allocated(this%total_probability)) deallocate(this%total_probability)
    if (allocated(this%probability_times)) deallocate(this%probability_times)
    if (allocated(this%LL)) deallocate(this%LL)
    if (allocated(this%mode_probabilities)) deallocate(this%mode_probabilities)
    if (allocated(this%densities)) deallocate(this%densities)
    if (allocated(this%grid1)) deallocate(this%grid1)
    if (allocated(this%grid2)) deallocate(this%grid2)
    if (allocated(this%ndim)) deallocate(this%ndim)
    
    if (this%verbose .and. this%enabled) then
      print '(a,i0,a)', "Wavepacket visualizer cleaned up. Calculated ", &
            this%current_step, " time steps."
    end if
  end subroutine cleanup_visualizer

  ! ============================================================================
  ! GETTERS AND SETTERS
  ! ============================================================================
  function is_enabled(this) result(enabled)
    class(wavepacket_visualizer), intent(in) :: this
    logical :: enabled
    enabled = this%enabled
  end function is_enabled
  
  subroutine enable(this)
    class(wavepacket_visualizer), intent(inout) :: this
    this%enabled = .true.
    if (this%verbose) print '(a)', "Wavepacket visualizer enabled."
  end subroutine enable
  
  subroutine disable(this)
    class(wavepacket_visualizer), intent(inout) :: this
    this%enabled = .false.
    if (this%verbose) print '(a)', "Wavepacket visualizer disabled."
  end subroutine disable
  
  subroutine set_calculation_frequency(this, freq)
    class(wavepacket_visualizer), intent(inout) :: this
    integer, intent(in) :: freq
    this%calculation_frequency = freq
    if (this%verbose .and. this%enabled) then
      print '(a,i0)', "Calculation frequency set to every ", freq, " steps"
    end if
  end subroutine set_calculation_frequency
  
  subroutine set_track_mode_probabilities(this, track)
    class(wavepacket_visualizer), intent(inout) :: this
    logical, intent(in) :: track
    this%track_mode_probabilities = track
    if (this%verbose .and. this%enabled) then
      if (track) then
        print '(a)', "Mode probability tracking enabled."
      else
        print '(a)', "Mode probability tracking disabled."
      end if
    end if
  end subroutine set_track_mode_probabilities
  
  function get_state_probabilities(this) result(prob)
    class(wavepacket_visualizer), intent(in) :: this
    real(dp), allocatable :: prob(:,:)
    if (allocated(this%state_probabilities)) then
      allocate(prob(this%nstate, this%current_step))
      prob = this%state_probabilities(:, 1:this%current_step)
    else
      allocate(prob(0, 0))
    end if
  end function get_state_probabilities
  
  function get_total_probability(this) result(prob)
    class(wavepacket_visualizer), intent(in) :: this
    real(dp), allocatable :: prob(:)
    if (allocated(this%total_probability)) then
      allocate(prob(this%current_step))
      prob = this%total_probability(1:this%current_step)
    else
      allocate(prob(0))
    end if
  end function get_total_probability
  
  function get_probability_times(this) result(times)
    class(wavepacket_visualizer), intent(in) :: this
    real(dp), allocatable :: times(:)
    if (allocated(this%probability_times)) then
      allocate(times(this%current_step))
      times = this%probability_times(1:this%current_step)
    else
      allocate(times(0))
    end if
  end function get_probability_times
  
  function get_LL_array(this) result(LL)
    class(wavepacket_visualizer), intent(in) :: this
    integer(ip), allocatable :: LL(:)
    if (allocated(this%LL)) then
      allocate(LL(size(this%LL)))
      LL = this%LL
    else
      allocate(LL(0))
    end if
  end function get_LL_array

  end module wavepacket_visualization
