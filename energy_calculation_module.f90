! ============================================================================
! ENERGY CALCULATION MODULE - BARRELENE VERSION
! Computes total energy, kinetic energy, potential energy, and energy components
! Compatible with TDDVR_barrelene_with_6_switches.f90
! ============================================================================

  module energy_calculation_module
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  
  private
  
  ! Public types and procedures
  public :: energy_calculator, dp, ip
  public :: create_energy_calculator
  
  ! Kind parameters
  integer, parameter :: dp = real64
  integer, parameter :: ip = int32
  
  ! Constants
  real(dp), parameter :: HBAR = 0.06350781278_dp
  real(dp), parameter :: EV_EPS = 0.9648455078_dp
  real(dp), parameter :: SMALL = 1.0e-12_dp
  
  ! Energy components type
  type :: energy_components
    real(dp) :: total_energy = 0.0_dp
    real(dp) :: kinetic_energy = 0.0_dp
    real(dp) :: potential_energy = 0.0_dp
    real(dp) :: classical_energy = 0.0_dp
    real(dp) :: quantum_energy = 0.0_dp
    real(dp) :: electronic_energy = 0.0_dp
    real(dp) :: vibrational_energy = 0.0_dp
    real(dp) :: coupling_energy = 0.0_dp
    real(dp) :: state_energies(6) = 0.0_dp
  end type energy_components
  
  ! Type for energy calculator
  type :: energy_calculator
    private
    logical :: enabled = .false.
    logical :: verbose = .true.
    logical :: parallel = .true.
    logical :: compute_components = .true.
    
    integer :: output_frequency = 10
    character(len=100) :: output_dir = "./energy_analysis/"
    
    type(energy_components), allocatable :: energy_history(:)
    real(dp), allocatable :: times(:)
    real(dp) :: energy_conservation_error = 0.0_dp
    real(dp) :: initial_total_energy = 0.0_dp
    
    integer :: unit_total = 500
    integer :: unit_components = 501
    integer :: unit_conservation = 502
    integer :: unit_state_energies = 503
    
    real(dp), allocatable :: grid_weights(:)
    
    integer(ip) :: nmode = 42
    integer(ip) :: nstate = 6
    integer(ip) :: nprod = 0
    
  contains
    procedure :: initialize => init_energy_calculator
    procedure :: compute_energies
    procedure :: compute_quantum_energy
    procedure :: compute_potential_energy
    procedure :: compute_state_energies
    procedure :: check_energy_conservation
    procedure :: save_energy_data
    procedure :: save_energy_components
    procedure :: save_state_energies
    procedure :: finalize => cleanup_energy_calculator
    procedure :: is_enabled
    procedure :: enable
    procedure :: disable
    procedure :: set_output_frequency
    procedure :: get_energy_conservation_error
    procedure :: print_energy_summary
    procedure :: setup_grid_weights
    procedure :: compute_quantum_energy_barrelene
  end type energy_calculator
  
  ! Constructor interface
  interface energy_calculator
    module procedure create_calculator_instance
  end interface

contains

  ! ============================================================================
  ! CONSTRUCTOR
  ! ============================================================================
  function create_calculator_instance(enabled, verbose, parallel, &
       compute_components, output_freq) result(calc)
    logical, intent(in), optional :: enabled, verbose, parallel
    logical, intent(in), optional :: compute_components
    integer, intent(in), optional :: output_freq
    type(energy_calculator) :: calc
    
    if (present(enabled)) calc%enabled = enabled
    if (present(verbose)) calc%verbose = verbose
    if (present(parallel)) calc%parallel = parallel
    if (present(compute_components)) calc%compute_components = compute_components
    if (present(output_freq)) calc%output_frequency = output_freq
    
    calc%nmode = 42
    calc%nstate = 6
    
    if (calc%verbose .and. calc%enabled) then
      print '(a)', 'Energy calculator created (enabled) for BARRELENE'
      print '(a,i0)', '  Vibrational modes: ', calc%nmode
      print '(a,i0)', '  Electronic states: ', calc%nstate
      print '(a,l)', '  Compute detailed components: ', calc%compute_components
      print '(a,i0)', '  Output frequency: every ', calc%output_frequency, ' steps'
      if (calc%parallel) print '(a)', '  Parallel execution: enabled'
    end if
  end function create_calculator_instance

  ! Alias for the main program
  function create_energy_calculator(enabled, verbose, parallel, &
       compute_components, output_freq) result(calc)
    logical, intent(in), optional :: enabled, verbose, parallel
    logical, intent(in), optional :: compute_components
    integer, intent(in), optional :: output_freq
    type(energy_calculator) :: calc
    
    calc = create_calculator_instance(enabled, verbose, parallel, &
         compute_components, output_freq)
  end function create_energy_calculator

  ! ============================================================================
  ! INITIALIZE ENERGY CALCULATOR
  ! ============================================================================
  subroutine init_energy_calculator(this, tstep, tint)
    class(energy_calculator), intent(inout) :: this
    integer(ip), intent(in) :: tstep
    real(dp), intent(in) :: tint
    
    integer :: ios
    character(len=200) :: dir_command
    integer :: nsteps_to_store
    
    if (.not. this%enabled) return
    
    if (this%verbose) then
      print '(a)', 'Initializing energy calculator for BARRELENE...'
      print '(a,i0)', '  Total time steps: ', tstep
      print '(a,f8.4,a)', '  Time step: ', tint, ' fs'
      print '(a,i0)', '  Vibrational modes: ', this%nmode
      print '(a,i0)', '  Electronic states: ', this%nstate
    end if
    
    nsteps_to_store = tstep / this%output_frequency + 1
    
    allocate(this%energy_history(nsteps_to_store))
    allocate(this%times(nsteps_to_store))
    
    this%times = 0.0_dp
    this%energy_conservation_error = 0.0_dp
    this%initial_total_energy = 0.0_dp
    
    dir_command = 'mkdir -p ' // trim(this%output_dir)
    call system(dir_command)
    
    this%unit_total = 500
    open(unit=this%unit_total, &
         file=trim(this%output_dir)//'total_energy.dat', &
         status='replace', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening total_energy.dat'
      this%enabled = .false.
      return
    end if
    
    write(this%unit_total, '(a)') '# Total Energy vs Time for BARRELENE'
    write(this%unit_total, '(a)') '# Time(fs), Total_Energy(eV), Kinetic_Energy(eV), Potential_Energy(eV)'
    
    if (this%compute_components) then
      this%unit_components = 501
      open(unit=this%unit_components, &
           file=trim(this%output_dir)//'energy_components.dat', &
           status='replace', iostat=ios)
      if (ios == 0) then
        write(this%unit_components, '(a)') '# Energy Components vs Time for BARRELENE'
!        write(this%unit_components, '(a)') '# Time(fs), Classical_E(eV), Quantum_E(eV), & Electronic_E(eV), Vibrational_E(eV), Coupling_E(eV)'
      end if
      
      this%unit_state_energies = 502
      open(unit=this%unit_state_energies, &
           file=trim(this%output_dir)//'state_energies.dat', &
           status='replace', iostat=ios)
      if (ios == 0) then
        write(this%unit_state_energies, '(a)') '# State Energies vs Time for BARRELENE'
        write(this%unit_state_energies, '(a)') '# Time(fs), State1(eV), State2(eV), State3(eV), State4(eV), State5(eV), State6(eV)'
      end if
    end if
    
    this%unit_conservation = 503
    open(unit=this%unit_conservation, &
         file=trim(this%output_dir)//'energy_conservation.dat', &
         status='replace', iostat=ios)
    if (ios == 0) then
      write(this%unit_conservation, '(a)') '# Energy Conservation Check for BARRELENE'
      write(this%unit_conservation, '(a)') '# Time(fs), Total_Energy(eV), Delta_E(eV), Relative_Error'
    end if
    
    if (this%verbose) then
      print '(a)', 'Energy calculator initialized successfully for BARRELENE.'
      print '(a,a)', '  Output directory: ', trim(this%output_dir)
    end if
    
  end subroutine init_energy_calculator

  ! ============================================================================
  ! SETUP GRID WEIGHTS
  ! ============================================================================
  subroutine setup_grid_weights(this, nprod, nndim, basis_type)
    class(energy_calculator), intent(inout) :: this
    integer(ip), intent(in) :: nprod
    integer(ip), intent(in), optional :: nndim(:)
    integer(ip), intent(in), optional :: basis_type(:)
    
    integer :: i
    real(dp) :: total_weight
    
    if (.not. this%enabled) return
    
    this%nprod = nprod
    allocate(this%grid_weights(nprod))
    
    if (present(nndim) .and. present(basis_type)) then
      do i = 1, nprod
        this%grid_weights(i) = 1.0_dp
      end do
      total_weight = real(nprod, dp)
      this%grid_weights = this%grid_weights / total_weight
    else
      this%grid_weights = 1.0_dp / real(nprod, dp)
    end if
    
    if (this%verbose) then
      print '(a,i0,a)', '  Grid weights initialized for BARRELENE: ', nprod, ' points'
    end if
    
  end subroutine setup_grid_weights

  ! ============================================================================
  ! COMPUTE ALL ENERGY COMPONENTS
  ! ============================================================================
  subroutine compute_energies(this, time, step, wavefunction, potential_matrices, &
       totmat, matbar, tot, coord, momen, freqn, ffn, newgrid)
    class(energy_calculator), intent(inout) :: this
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    complex(dp), intent(in) :: wavefunction(:,:)
    real(dp), intent(in) :: potential_matrices(:,:,:)
    real(dp), intent(in) :: totmat(:,:,:)
    real(dp), intent(in) :: matbar(:,:,:,:)
    integer(ip), intent(in) :: tot(:,:)
    real(dp), intent(in) :: coord(:), momen(:)
    real(dp), intent(in) :: freqn(:), ffn(:)
    real(dp), intent(in) :: newgrid(:,:)
    
    type(energy_components) :: energies
    integer(ip) :: nstate, nprod, nmode, nn
    integer(ip) :: storage_index
    real(dp) :: total_prob
    integer :: i, j
    
    if (.not. this%enabled) return
    
    if (mod(step, this%output_frequency) /= 0) return
    
    nstate = size(wavefunction, 1)
    nprod = size(wavefunction, 2)
    nmode = size(coord)
    nn = size(newgrid, 1)
    
    if (nstate /= this%nstate .or. nmode /= this%nmode) then
      print *, 'ERROR: Array dimensions do not match BARRELENE parameters!'
      return
    end if
    
    if (.not. allocated(this%grid_weights)) then
      call this%setup_grid_weights(nprod)
    end if
    
    if (this%verbose .and. mod(step, 100) == 0) then
      print '(a,f8.3,a,i0)', 'Computing BARRELENE energies at t=', time, &
            ' fs (step ', step, ')'
    end if
    
    ! Check normalization
    total_prob = 0.0_dp
    do i = 1, nprod
      do j = 1, nstate
        total_prob = total_prob + abs(wavefunction(j, i))**2
      end do
    end do
    
    if (abs(total_prob - 1.0_dp) > 1.0e-4_dp) then
      print '(a,es16.8)', 'WARNING: Wavefunction probability = ', total_prob
    end if
    
    ! Compute energies
    energies%classical_energy = 0.5_dp * (sum(momen**2) + sum(freqn**2 * coord**2))
    energies%kinetic_energy = 0.5_dp * sum(momen**2)
    
    energies%quantum_energy = this%compute_quantum_energy_barrelene( &
         wavefunction, potential_matrices, totmat, matbar, tot, coord, momen, &
         freqn, ffn, newgrid)
    
    if (this%compute_components) then
      call this%compute_potential_energy(wavefunction, potential_matrices, &
           energies%electronic_energy, energies%vibrational_energy, &
           energies%coupling_energy)
      call this%compute_state_energies(wavefunction, potential_matrices, tot, &
           energies%state_energies)
    end if
    
    energies%potential_energy = energies%classical_energy + energies%quantum_energy
    energies%total_energy = energies%classical_energy + energies%quantum_energy
    
    storage_index = step / this%output_frequency + 1
    if (storage_index <= size(this%energy_history)) then
      this%energy_history(storage_index) = energies
      this%times(storage_index) = time
      
      if (storage_index == 1) then
        this%initial_total_energy = energies%total_energy
        if (this%verbose) then
          print '(a,es16.8,a)', '  Initial total energy: ', &
                this%initial_total_energy, ' eV'
        end if
      end if
    end if
    
    call this%save_energy_data(time, energies)
    
    if (this%compute_components) then
      call this%save_energy_components(time, energies)
      call this%save_state_energies(time, energies)
    end if
    
    call this%check_energy_conservation(time, energies%total_energy)
    
  end subroutine compute_energies

  ! ============================================================================
  ! COMPUTE QUANTUM ENERGY
  ! ============================================================================
  function compute_quantum_energy_barrelene(this, wavefunction, potential_matrices, &
       totmat, matbar, tot, coord, momen, freqn, ffn, newgrid) result(quantum_energy)
    class(energy_calculator), intent(in) :: this
    complex(dp), intent(in) :: wavefunction(:,:)
    real(dp), intent(in) :: potential_matrices(:,:,:)
    real(dp), intent(in) :: totmat(:,:,:)
    real(dp), intent(in) :: matbar(:,:,:,:)
    integer(ip), intent(in) :: tot(:,:)
    real(dp), intent(in) :: coord(:), momen(:)
    real(dp), intent(in) :: freqn(:), ffn(:)
    real(dp), intent(in) :: newgrid(:,:)
    
    real(dp) :: quantum_energy, kinetic_part, potential_part
    integer(ip) :: nstate, nprod, nmode
    integer(ip) :: i, j, k, m, idx
    real(dp), allocatable :: facp(:), facq(:)
    complex(dp) :: psi_conj
    real(dp) :: local_kinetic
    
    nstate = size(wavefunction, 1)
    nprod = size(wavefunction, 2)
    nmode = size(coord)
    
    quantum_energy = 0.0_dp
    kinetic_part = 0.0_dp
    potential_part = 0.0_dp
    
    allocate(facp(nmode), facq(nmode))
    
    facp = 0.5_dp * momen**2
    facq = 0.5_dp * sqrt(HBAR / max(freqn/2.0_dp, SMALL)) * coord * freqn**2
    
    do i = 1, nprod
      do j = 1, nstate
        psi_conj = conjg(wavefunction(j, i))
        
        local_kinetic = 0.0_dp
        do m = 1, nmode
          idx = tot(i, m)
          if (idx > 0 .and. idx <= size(matbar, 2)) then
            local_kinetic = local_kinetic + &
                 (facp(m) + facq(m) * matbar(m, idx, idx, 3) + &
                 0.5_dp * totmat(m, idx, idx))
          end if
        end do
        
        kinetic_part = kinetic_part + &
             real(psi_conj * wavefunction(j, i)) * local_kinetic * &
             this%grid_weights(i)
        
        potential_part = potential_part + &
             real(psi_conj * wavefunction(j, i)) * &
             potential_matrices(j, j, i) * this%grid_weights(i)
        
        do k = 1, nstate
          if (k /= j) then
            potential_part = potential_part + &
                 real(psi_conj * wavefunction(k, i)) * &
                 potential_matrices(j, k, i) * this%grid_weights(i)
          end if
        end do
      end do
    end do
    
    quantum_energy = kinetic_part + potential_part
    
    deallocate(facp, facq)
    
  end function compute_quantum_energy_barrelene

  ! ============================================================================
  ! COMPUTE POTENTIAL ENERGY COMPONENTS
  ! ============================================================================
  subroutine compute_potential_energy(this, wavefunction, potential_matrices, &
       electronic_energy, vibrational_energy, coupling_energy)
    class(energy_calculator), intent(in) :: this
    complex(dp), intent(in) :: wavefunction(:,:)
    real(dp), intent(in) :: potential_matrices(:,:,:)
    real(dp), intent(out) :: electronic_energy, vibrational_energy, coupling_energy
    
    integer(ip) :: nstate, nprod, i, j, k
    complex(dp) :: psi_conj
    
    nstate = size(wavefunction, 1)
    nprod = size(wavefunction, 2)
    
    electronic_energy = 0.0_dp
    vibrational_energy = 0.0_dp
    coupling_energy = 0.0_dp
    
    do i = 1, nprod
      do j = 1, nstate
        psi_conj = conjg(wavefunction(j, i))
        
        if (j == 1) then
          vibrational_energy = vibrational_energy + &
               real(psi_conj * wavefunction(j, i)) * &
               potential_matrices(j, j, i) * 0.7_dp * this%grid_weights(i)
          electronic_energy = electronic_energy + &
               real(psi_conj * wavefunction(j, i)) * &
               potential_matrices(j, j, i) * 0.3_dp * this%grid_weights(i)
        else
          vibrational_energy = vibrational_energy + &
               real(psi_conj * wavefunction(j, i)) * &
               potential_matrices(j, j, i) * 0.3_dp * this%grid_weights(i)
          electronic_energy = electronic_energy + &
               real(psi_conj * wavefunction(j, i)) * &
               potential_matrices(j, j, i) * 0.7_dp * this%grid_weights(i)
        end if
        
        do k = 1, nstate
          if (k /= j) then
            coupling_energy = coupling_energy + &
                 real(psi_conj * wavefunction(k, i)) * &
                 potential_matrices(j, k, i) * this%grid_weights(i)
          end if
        end do
      end do
    end do
    
  end subroutine compute_potential_energy

  ! ============================================================================
  ! COMPUTE STATE-SPECIFIC ENERGIES
  ! ============================================================================
  subroutine compute_state_energies(this, wavefunction, potential_matrices, tot, &
       state_energies)
    class(energy_calculator), intent(in) :: this
    complex(dp), intent(in) :: wavefunction(:,:)
    real(dp), intent(in) :: potential_matrices(:,:,:)
    integer(ip), intent(in) :: tot(:,:)
    real(dp), intent(out) :: state_energies(:)
    
    integer(ip) :: nstate, nprod, i, j, k
    real(dp), allocatable :: state_probability(:)
    complex(dp) :: psi_conj
    
    nstate = size(wavefunction, 1)
    nprod = size(wavefunction, 2)
    
    if (size(state_energies) /= this%nstate) then
      print *, 'ERROR: state_energies array size mismatch for BARRELENE!'
      return
    end if
    
    allocate(state_probability(this%nstate))
    
    state_energies = 0.0_dp
    state_probability = 0.0_dp
    
    do i = 1, nprod
      do j = 1, nstate
        psi_conj = conjg(wavefunction(j, i))
        
        state_probability(j) = state_probability(j) + &
             abs(wavefunction(j, i))**2 * this%grid_weights(i)
        
        state_energies(j) = state_energies(j) + &
             real(psi_conj * wavefunction(j, i)) * &
             potential_matrices(j, j, i) * this%grid_weights(i)
        
        do k = 1, nstate
          if (k /= j) then
            state_energies(j) = state_energies(j) + &
                 real(psi_conj * wavefunction(k, i)) * &
                 potential_matrices(j, k, i) * this%grid_weights(i)
          end if
        end do
      end do
    end do
    
    do j = 1, nstate
      if (state_probability(j) > SMALL) then
        state_energies(j) = state_energies(j) / state_probability(j)
      else
        state_energies(j) = 0.0_dp
      end if
    end do
    
    deallocate(state_probability)
    
  end subroutine compute_state_energies

  ! ============================================================================
  ! CHECK ENERGY CONSERVATION
  ! ============================================================================
  subroutine check_energy_conservation(this, time, current_total_energy)
    class(energy_calculator), intent(inout) :: this
    real(dp), intent(in) :: time, current_total_energy
    
    real(dp) :: delta_energy, relative_error
    
    if (.not. this%enabled) return
    
    delta_energy = current_total_energy - this%initial_total_energy
    
    if (abs(this%initial_total_energy) > 1.0e-6_dp) then
      relative_error = abs(delta_energy) / abs(this%initial_total_energy)
    else
      relative_error = abs(delta_energy)
    end if
    
    this%energy_conservation_error = max(this%energy_conservation_error, relative_error)
    
    write(this%unit_conservation, '(4es16.8)') &
         time, current_total_energy, delta_energy, relative_error
    
    if (relative_error > 1.0e-4_dp .and. this%verbose) then
      print '(a,f8.3,a,es10.2,a,es10.2,a)', &
           'BARRELENE energy conservation: t=', time, &
           ' fs, Delta_E=', delta_energy, ' eV, rel=', relative_error
    end if
    
  end subroutine check_energy_conservation

  ! ============================================================================
  ! SAVE ENERGY DATA
  ! ============================================================================
  subroutine save_energy_data(this, time, energies)
    class(energy_calculator), intent(in) :: this
    real(dp), intent(in) :: time
    type(energy_components), intent(in) :: energies
    
    if (.not. this%enabled) return
    
    write(this%unit_total, '(4es20.10)') &
         time, energies%total_energy, energies%kinetic_energy, energies%potential_energy
    flush(this%unit_total)
    
  end subroutine save_energy_data

  ! ============================================================================
  ! SAVE ENERGY COMPONENTS
  ! ============================================================================
  subroutine save_energy_components(this, time, energies)
    class(energy_calculator), intent(in) :: this
    real(dp), intent(in) :: time
    type(energy_components), intent(in) :: energies
    
    if (.not. this%enabled .or. .not. this%compute_components) return
    
    write(this%unit_components, '(6es20.10)') &
         time, energies%classical_energy, energies%quantum_energy, &
         energies%electronic_energy, energies%vibrational_energy, energies%coupling_energy
    flush(this%unit_components)
    
  end subroutine save_energy_components

  ! ============================================================================
  ! SAVE STATE ENERGIES
  ! ============================================================================
  subroutine save_state_energies(this, time, energies)
    class(energy_calculator), intent(in) :: this
    real(dp), intent(in) :: time
    type(energy_components), intent(in) :: energies
    
    if (.not. this%enabled .or. .not. this%compute_components) return
    
    write(this%unit_state_energies, '(7es20.10)') &
         time, energies%state_energies
    flush(this%unit_state_energies)
    
  end subroutine save_state_energies

  ! ============================================================================
  ! PRINT ENERGY SUMMARY
  ! ============================================================================
  subroutine print_energy_summary(this)
    class(energy_calculator), intent(in) :: this
    
    real(dp) :: avg_total, avg_kinetic, avg_potential
    real(dp) :: std_total, std_kinetic, std_potential
    integer :: npoints, i
    
    if (.not. this%enabled) return
    
    npoints = count(this%times > 0.0_dp)
    
    if (npoints > 0) then
      avg_total = sum(this%energy_history(1:npoints)%total_energy) / npoints
      avg_kinetic = sum(this%energy_history(1:npoints)%kinetic_energy) / npoints
      avg_potential = sum(this%energy_history(1:npoints)%potential_energy) / npoints
      
      std_total = sqrt(sum((this%energy_history(1:npoints)%total_energy - avg_total)**2) / npoints)
      std_kinetic = sqrt(sum((this%energy_history(1:npoints)%kinetic_energy - avg_kinetic)**2) / npoints)
      std_potential = sqrt(sum((this%energy_history(1:npoints)%potential_energy - avg_potential)**2) / npoints)
      
      print '(a)', ''
      print '(a)', 'BARRELENE ENERGY ANALYSIS SUMMARY:'
      print '(a)', '-----------------------------------'
      print '(a,i0)', 'Number of energy points: ', npoints
      print '(a,es16.8,a)', 'Initial total energy: ', this%initial_total_energy, ' eV'
      print '(a)', ''
      print '(a)', 'Average energies (eV):'
      print '(a,es16.8,a,es16.8)', '  Total:     ', avg_total, ' +/- ', std_total
      print '(a,es16.8,a,es16.8)', '  Kinetic:   ', avg_kinetic, ' +/- ', std_kinetic
      print '(a,es16.8,a,es16.8)', '  Potential: ', avg_potential, ' +/- ', std_potential
      print '(a)', ''
      print '(a,es12.4)', 'Maximum energy conservation error: ', this%energy_conservation_error
      
      if (this%compute_components) then
        print '(a)', ''
        print '(a)', 'Final state energies (eV):'
        do i = 1, this%nstate
          print '(a,i0,a,es16.8)', '  State ', i, ': ', &
               this%energy_history(npoints)%state_energies(i)
        end do
      end if
    end if
    
  end subroutine print_energy_summary

  ! ============================================================================
  ! CLEANUP
  ! ============================================================================
  subroutine cleanup_energy_calculator(this)
    class(energy_calculator), intent(inout) :: this
    
    if (allocated(this%energy_history)) deallocate(this%energy_history)
    if (allocated(this%times)) deallocate(this%times)
    if (allocated(this%grid_weights)) deallocate(this%grid_weights)
    
    if (this%unit_total > 0) close(this%unit_total)
    if (this%unit_components > 0) close(this%unit_components)
    if (this%unit_conservation > 0) close(this%unit_conservation)
    if (this%unit_state_energies > 0) close(this%unit_state_energies)
    
    this%unit_total = 500
    this%unit_components = 501
    this%unit_conservation = 502
    this%unit_state_energies = 503
    
    if (this%verbose .and. this%enabled) then
      print '(a)', 'BARRELENE energy calculator cleaned up.'
      if (allocated(this%energy_history)) then
        call this%print_energy_summary()
      end if
    end if
  end subroutine cleanup_energy_calculator

  ! ============================================================================
  ! GETTERS AND SETTERS
  ! ============================================================================
  function is_enabled(this) result(enabled)
    class(energy_calculator), intent(in) :: this
    logical :: enabled
    enabled = this%enabled
  end function is_enabled
  
  subroutine enable(this)
    class(energy_calculator), intent(inout) :: this
    this%enabled = .true.
    if (this%verbose) print '(a)', 'BARRELENE energy calculator enabled.'
  end subroutine enable
  
  subroutine disable(this)
    class(energy_calculator), intent(inout) :: this
    this%enabled = .false.
    if (this%verbose) print '(a)', 'BARRELENE energy calculator disabled.'
  end subroutine disable
  
  subroutine set_output_frequency(this, freq)
    class(energy_calculator), intent(inout) :: this
    integer, intent(in) :: freq
    this%output_frequency = freq
    if (this%verbose .and. this%enabled) then
      print '(a,i0)', 'Output frequency set to every ', freq, ' steps'
    end if
  end subroutine set_output_frequency
  
  function get_energy_conservation_error(this) result(error)
    class(energy_calculator), intent(in) :: this
    real(dp) :: error
    error = this%energy_conservation_error
  end function get_energy_conservation_error

  ! ============================================================================
  ! COMPATIBILITY FUNCTION
  ! ============================================================================
  function compute_quantum_energy(this, wavefunction, potential_matrices, &
       totmat, matbar, tot, coord, momen, freqn) result(quantum_energy)
    class(energy_calculator), intent(in) :: this
    complex(dp), intent(in) :: wavefunction(:,:)
    real(dp), intent(in) :: potential_matrices(:,:,:)
    real(dp), intent(in) :: totmat(:,:,:)
    real(dp), intent(in) :: matbar(:,:,:,:)
    integer(ip), intent(in) :: tot(:,:)
    real(dp), intent(in) :: coord(:), momen(:)
    real(dp), intent(in) :: freqn(:)
    
    real(dp) :: quantum_energy
    real(dp), allocatable :: ffn(:), newgrid(:,:)
    
    allocate(ffn(size(freqn)))
    ffn = sqrt(freqn / HBAR)
    
    allocate(newgrid(1, size(coord)))
    newgrid = 0.0_dp
    
    quantum_energy = this%compute_quantum_energy_barrelene( &
         wavefunction, potential_matrices, totmat, matbar, tot, coord, momen, &
         freqn, ffn, newgrid)
    
    deallocate(ffn, newgrid)
    
  end function compute_quantum_energy

  end module energy_calculation_module
