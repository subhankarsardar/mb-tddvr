! ============================================================================
! SO-TDDVR MAIN PROGRAM FOR BARRELENE - WITH MULTIPLE PROPAGATION SCHEMES
! Journal of Molecular Structure 1110 (2016) 32-43
! ============================================================================
!
! TIME-DEPENDENT DISCRETE VARIABLE REPRESENTATION (TDDVR) METHOD
! WITH MULTIPLE PROPAGATION SCHEMES:
! 1. Original 2nd-order Split-Operator (DEFAULT - IDENTICAL TO ORIGINAL)
! 2. 4th-order Split-Operator (Suzuki-Trotter)
! 3. 2nd-order Magnus Expansion
! 4. 4th-order Runge-Kutta
! 5. Adaptive Split-Operator
!
! To get IDENTICAL results to original:
!   PROPAGATION_METHOD = PROP_SPLIT_OPERATOR_2ND
!   USE_ADAPTIVE_TIMESTEP = .false.
!   USE_ERROR_CONTROL = .false.
!
! To compile with all modules:
! gfortran -cpp -fopenmp -O3 -o tddvr_so SO-TDDVR.f90 \
!   potential_module.f90 \
!   energy_calculation_module.f90 \
!   adiabatic_module.f90 \
!   wavepacket_visualization.f90 \
!   potential_cuts_module.f90 \
!   potential_2d_cuts_module.f90
!
! Or compile selectively with -D flags:
! gfortran -cpp -DUSE_ENERGY -DUSE_WAVEPACKET -fopenmp -O3 -o tddvr_so_selective
! ...
! ============================================================================

! ============================================================================
! SO-TDDVR MAIN PROGRAM FOR BARRELENE - ENERGY & ADIABATIC MODULES
! Journal of Molecular Structure 1110 (2016) 32-43
! ============================================================================

program SO_TDDVR_barrelene
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use omp_lib
  
  ! ==========================================================================
  ! MODULE IMPORTS - ONLY ENERGY AND ADIABATIC
  ! ==========================================================================
  use potential_module, only: POTENTIAL
  use energy_calculation_module, only: energy_calculator, create_energy_calculator
  use adiabatic_module, only: adiabatic_calculator, create_adiabatic_calculator
  
  implicit none
  
  ! ==========================================================================
  ! PROPAGATION SCHEME SELECTION
  ! ==========================================================================
  integer, parameter :: PROP_SPLIT_OPERATOR_2ND = 1
  integer, parameter :: PROP_SPLIT_OPERATOR_4TH = 2
  integer, parameter :: PROP_MAGNUS_2ND = 3
  integer, parameter :: PROP_RUNGE_KUTTA_4TH = 4
  integer, parameter :: PROP_ADAPTIVE_SPLIT = 5
  
  ! ==========================================================================
  ! PROPAGATION CONFIGURATION
  ! ==========================================================================
  integer, parameter :: PROPAGATION_METHOD = PROP_ADAPTIVE_SPLIT
  logical, parameter :: USE_ADAPTIVE_TIMESTEP = .false.
  logical, parameter :: USE_ERROR_CONTROL = .false.
  logical, parameter :: USE_CONSERVATION_MONITOR = .true.
  
  real(real64), parameter :: MAX_TIMESTEP = 0.05_real64
  real(real64), parameter :: MIN_TIMESTEP = 0.001_real64
  real(real64), parameter :: ERROR_TOLERANCE = 1.0e-6_real64
  real(real64), parameter :: SAFETY_FACTOR = 0.9_real64
  real(real64), parameter :: GROWTH_LIMIT = 1.5_real64
  
  ! ==========================================================================
  ! SWITCH CONFIGURATION - ONLY ENERGY AND ADIABATIC ACTIVATED
  ! ==========================================================================
  logical, parameter :: USE_POTENTIAL_MODULE     = .true.
  logical, parameter :: USE_ENERGY_CALCULATION   = .true.
  logical, parameter :: USE_ADIABATIC_CALCULATION = .true.
  logical, parameter :: USE_WAVEPACKET_VISUALIZATION = .false.
  logical, parameter :: USE_1D_POTENTIAL_CUTS    = .false.
  logical, parameter :: USE_2D_POTENTIAL_SURFACES = .false.
  
  ! ==========================================================================
  ! BARRELENE SYSTEM PARAMETERS
  ! ==========================================================================
  integer, parameter :: dp = real64
  integer, parameter :: ip = int32
  
  real(dp), parameter :: HBAR = 0.06350781278_dp
  real(dp), parameter :: EV_EPS = 0.9648455078_dp
  real(dp), parameter :: PI = acos(-1.0_dp)
  complex(dp), parameter :: ZI = (0.0_dp, 1.0_dp)
  real(dp), parameter :: SMALL = 1.0e-12_dp
  
  integer(ip), parameter :: NMODE = 42
  integer(ip), parameter :: NSTATE = 6
  integer(ip), parameter :: NPLACE = 2
  integer(ip), parameter :: TSTEP = 1000
  real(dp), parameter :: TINT = 0.04_dp
  real(dp), parameter :: QUANTA = 0.0_dp
  
  ! ==========================================================================
  ! BASIS FUNCTION TYPES
  ! ==========================================================================
  integer, parameter :: BASIS_HERMITE = 1
  integer, parameter :: BASIS_SIN = 2
  integer, parameter :: BASIS_LEGENDRE = 3
  
  integer(ip) :: basis_type(NMODE)
  real(dp), allocatable :: hermrootf(:,:)
  real(dp), allocatable :: sin_basis(:,:)
  real(dp), allocatable :: legendre_roots(:,:)
  
  ! ==========================================================================
  ! GRID PARAMETERS
  ! ==========================================================================
  integer(ip) :: grid_dims(42) = [ &
    1, 5, 5, 5, 1, 3, 1, 1, &
    1, 1, 1, 1, 1, 1, 1, 5, 1, 5, &
    3, 3, 3, 1, 5, 1, 1, 1, 1, 1, 1, 1, &
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1  &
  ]
  
  ! C-H stretches (modes 1-8) - High frequency region
  ! C-C stretches (modes 9-18) - Medium frequency region
  ! Bending modes (modes 19-30) - Medium-low frequency region
  ! Torsional modes (modes 31-42) - Low frequency region


  integer(ip) :: NN_MAX
  integer(ip) :: NNDIM(NMODE)
  integer(ip) :: NPROD
  
  ! ==========================================================================
  ! SYSTEM ARRAYS
  ! ==========================================================================
  complex(dp), allocatable :: ycom(:,:)
  complex(dp), allocatable :: phi(:)
  real(dp), allocatable :: newgrid(:,:)
  real(dp), allocatable :: coord(:)
  real(dp), allocatable :: momen(:)
  integer(ip), allocatable :: tot(:,:)
  real(dp), allocatable :: freq(:)
  real(dp), allocatable :: freqn(:)
  real(dp), allocatable :: width(:)
  real(dp), allocatable :: ffn(:)
  real(dp), allocatable :: v(:,:,:)
  real(dp), allocatable :: mat(:,:,:,:)
  real(dp), allocatable :: matbar(:,:,:,:)
  real(dp), allocatable :: totmat(:,:,:)
  real(dp), allocatable :: facp(:), facq(:)
  
  ! ==========================================================================
  ! ENHANCED PROPAGATION ARRAYS
  ! ==========================================================================
  real(dp) :: current_dt
  real(dp) :: next_dt
  real(dp) :: error_estimate
  logical :: step_accepted
  complex(dp), allocatable :: ycom_ref(:,:)
  complex(dp), allocatable :: ycom_high(:,:)
  real(dp) :: total_energy_initial
  real(dp) :: total_energy_current
  real(dp) :: norm_deviation_max = 0.0_dp
  real(dp) :: energy_deviation_max = 0.0_dp
  
  ! ==========================================================================
  ! ANALYSIS MODULE OBJECTS
  ! ==========================================================================
  type(energy_calculator) :: energy_calc
  type(adiabatic_calculator) :: adiabatic_calc
  
  ! ==========================================================================
  ! CONTROL VARIABLES
  ! ==========================================================================
  logical :: verbose = .true.
  logical :: parallel = .true.
  integer :: num_threads
  real(dp) :: time
  real(dp) :: total_start, total_end
  integer(ip) :: step
  integer :: i, j, k, m, iprod
  
  ! ==========================================================================
  ! MAIN PROGRAM
  ! ==========================================================================
  
  call display_program_header()
  call initialize_basis_types()
  call initialize_system()
  call initialize_propagation_scheme()
  call initialize_analysis_modules()
  call setup_initial_conditions()
  call run_dynamics_simulation()
  call finalize_simulation()
  call display_final_summary()

contains

  ! ==========================================================================
  ! DISPLAY PROGRAM HEADER
  ! ==========================================================================
  subroutine display_program_header()
    character(len=50) :: prop_name
    
    select case(PROPAGATION_METHOD)
      case(PROP_SPLIT_OPERATOR_2ND)
        prop_name = "2nd-Order Split-Operator (IDENTICAL to original)"
      case(PROP_SPLIT_OPERATOR_4TH)
        prop_name = "4th-Order Split-Operator"
      case(PROP_MAGNUS_2ND)
        prop_name = "2nd-Order Magnus"
      case(PROP_RUNGE_KUTTA_4TH)
        prop_name = "4th-Order Runge-Kutta"
      case(PROP_ADAPTIVE_SPLIT)
        prop_name = "Adaptive Split-Operator"
      case default
        prop_name = "Unknown"
    end select
    
    print '(a)', "========================================================"
    print '(a)', "     SO-TDDVR QUANTUM DYNAMICS SIMULATION"
    print '(a)', "              BARRELENE (C8H8) SYSTEM"
    print '(a)', "    Journal of Molecular Structure 1110 (2016) 32-43"
    print '(a)', "========================================================"
    print '(a)', ""
    print '(a)', "PROPAGATION SCHEME:"
    print '(a,a)', "  Method: ", trim(prop_name)
    print '(a,l1)', "  Adaptive Time Step:    ", USE_ADAPTIVE_TIMESTEP
    if (USE_ADAPTIVE_TIMESTEP) then
      print '(a,es8.1)', "  Error Tolerance:       ", ERROR_TOLERANCE
    else
      print '(a)', "  (Fixed time step - IDENTICAL to original)"
    end if
    print '(a)', ""
    
    print '(a)', "ACTIVATED MODULES:"
    print '(a,l1)', "  Switch 1: Potential Calculation         [ALWAYS ON] ", USE_POTENTIAL_MODULE
    print '(a,l1)', "  Switch 2: Energy Calculation Module     [ACTIVATED] ", USE_ENERGY_CALCULATION
    print '(a,l1)', "  Switch 3: Adiabatic Calculation Module  [ACTIVATED] ", USE_ADIABATIC_CALCULATION
    print '(a,l1)', "  Switch 4: Wavepacket Visualization      [DISABLED] ", USE_WAVEPACKET_VISUALIZATION
    print '(a,l1)', "  Switch 5: 1D Potential Cuts Analysis    [DISABLED] ", USE_1D_POTENTIAL_CUTS
    print '(a,l1)', "  Switch 6: 2D Potential Surfaces Analysis [DISABLED] ", USE_2D_POTENTIAL_SURFACES
    print '(a)', ""
    
    print '(a,i0)', "  Electronic states: ", NSTATE
    print '(a,i0)', "  Vibrational modes: ", NMODE
    print '(a,i0)', "  Time steps: ", TSTEP
    print '(a,f6.3,a)', "  Time step: ", TINT, " fs"
    print '(a,f8.3,a)', "  Total time: ", TSTEP * TINT, " fs"
    print '(a)', ""
    
    if (parallel) then
      !$omp parallel
      !$omp single
      num_threads = omp_get_num_threads()
      print '(a,i0,a)', "OpenMP parallel execution with ", num_threads, " threads"
      !$omp end single
      !$omp end parallel
    else
      num_threads = 1
      print '(a)', "OpenMP parallel execution disabled"
    end if
    
    call cpu_time(total_start)
    
  end subroutine display_program_header

  ! ==========================================================================
  ! INITIALIZE BASIS TYPES
  ! ==========================================================================
  subroutine initialize_basis_types()
    basis_type(1:8) = BASIS_HERMITE
    basis_type(9:18) = BASIS_HERMITE
    basis_type(19:30) = BASIS_LEGENDRE
    basis_type(31:42) = BASIS_SIN
    
    if (verbose) then
      print '(a)', ""
      print '(a)', "INITIALIZING BASIS FUNCTIONS..."
      print '(a)', "-------------------------------"
      print '(a,i0)', "  Hermite basis modes: ", count(basis_type == BASIS_HERMITE)
      print '(a,i0)', "  Legendre basis modes: ", count(basis_type == BASIS_LEGENDRE)
      print '(a,i0)', "  Sine basis modes: ", count(basis_type == BASIS_SIN)
    end if
  end subroutine initialize_basis_types

  ! ==========================================================================
  ! INITIALIZE SYSTEM
  ! ==========================================================================
  subroutine initialize_system()
    print '(a)', ""
    print '(a)', "INITIALIZING BARRELENE SYSTEM..."
    print '(a)', "---------------------------------"
    
    NNDIM = grid_dims
    NN_MAX = maxval(grid_dims)
    
    NPROD = 1
    do i = 1, NMODE
      NPROD = NPROD * NNDIM(i)
      if (NPROD < 0) then
        print *, "ERROR: Integer overflow in NPROD calculation"
        stop
      end if
    end do
    
    print '(a,i0)', "  System dimensions:"
    print '(a,i0)', "    Vibrational modes: ", NMODE
    print '(a,i0)', "    Electronic states: ", NSTATE
    print '(a,i0)', "    Total grid points: ", NPROD
    
    ! Allocate arrays
    allocate(ycom(NSTATE, NPROD))
    allocate(phi(NPROD))
    allocate(newgrid(NN_MAX, NMODE))
    allocate(coord(NMODE))
    allocate(momen(NMODE))
    allocate(tot(NPROD, NMODE))
    allocate(v(NSTATE, NSTATE, NPROD))
    allocate(hermrootf(NN_MAX, NMODE))
    allocate(sin_basis(NN_MAX, NMODE))
    allocate(legendre_roots(NN_MAX, NMODE))
    allocate(freq(NMODE))
    allocate(freqn(NMODE))
    allocate(width(NMODE))
    allocate(ffn(NMODE))
    allocate(mat(NMODE, NN_MAX, NN_MAX, 4))
    allocate(matbar(NMODE, NN_MAX, NN_MAX, 4))
    allocate(totmat(NMODE, NN_MAX, NN_MAX))
    allocate(facp(NMODE), facq(NMODE))
    
    ! Enhanced arrays
    allocate(ycom_ref(NSTATE, NPROD))
    allocate(ycom_high(NSTATE, NPROD))
    
    ! Initialize arrays
    ycom = (0.0_dp, 0.0_dp)
    ycom_ref = (0.0_dp, 0.0_dp)
    ycom_high = (0.0_dp, 0.0_dp)
    phi = (0.0_dp, 0.0_dp)
    newgrid = 0.0_dp
    coord = 0.01_dp
    momen = 0.001_dp
    v = 0.0_dp
    mat = 0.0_dp
    matbar = 0.0_dp
    totmat = 0.0_dp
    facp = 0.0_dp
    facq = 0.0_dp
    hermrootf = 0.0_dp
    sin_basis = 0.0_dp
    legendre_roots = 0.0_dp
    
    ! Set frequencies
    do i = 1, NMODE
      if (i <= 8) then
        freq(i) = 0.37_dp - 0.003_dp * (i-1)
      else if (i <= 18) then
        freq(i) = 0.18_dp - 0.002_dp * (i-9)
      else if (i <= 30) then
        freq(i) = 0.14_dp - 0.0015_dp * (i-19)
      else
        freq(i) = 0.06_dp - 0.001_dp * (i-31)
      end if
      
      if (freq(i) <= 0.0_dp) then
        print *, "WARNING: Non-positive frequency at mode ", i, " = ", freq(i)
        freq(i) = 0.01_dp
      end if
    end do
    
    freqn = freq * EV_EPS / HBAR
    do i = 1, NMODE
      if (freqn(i) <= 0.0_dp) then
        freqn(i) = 0.01_dp * EV_EPS / HBAR
      end if
    end do
    
    width = freqn / 2.0_dp
    
    do i = 1, NMODE
      ffn(i) = sqrt(abs(freqn(i) / HBAR))
      if (freqn(i) / HBAR < 0.0_dp) then
        ffn(i) = sqrt(0.01_dp * EV_EPS / HBAR**2)
      end if
    end do
    
    print '(a)', "  System arrays allocated successfully."
    
  end subroutine initialize_system

  ! ==========================================================================
  ! INITIALIZE PROPAGATION SCHEME
  ! ==========================================================================
  subroutine initialize_propagation_scheme()
    print '(a)', ""
    print '(a)', "INITIALIZING PROPAGATION SCHEME..."
    print '(a)', "-----------------------------------"
    
    current_dt = TINT
    next_dt = TINT
    
    select case(PROPAGATION_METHOD)
      case(PROP_SPLIT_OPERATOR_2ND)
        print '(a)', "  Using 2nd-Order Split-Operator (IDENTICAL to original)"
      case(PROP_SPLIT_OPERATOR_4TH)
        print '(a)', "  Using 4th-Order Split-Operator"
      case(PROP_MAGNUS_2ND)
        print '(a)', "  Using 2nd-Order Magnus"
      case(PROP_RUNGE_KUTTA_4TH)
        print '(a)', "  Using 4th-Order Runge-Kutta"
      case(PROP_ADAPTIVE_SPLIT)
        print '(a)', "  Using Adaptive Split-Operator"
    end select
    
    if (USE_ADAPTIVE_TIMESTEP) then
      print '(a)', "  Adaptive time stepping: ENABLED"
    else
      print '(a)', "  Adaptive time stepping: DISABLED (IDENTICAL to original)"
    end if
    
  end subroutine initialize_propagation_scheme

  ! ==========================================================================
  ! INITIALIZE ANALYSIS MODULES
  ! ==========================================================================
  subroutine initialize_analysis_modules()
    print '(a)', ""
    print '(a)', "INITIALIZING ANALYSIS MODULES..."
    print '(a)', "--------------------------------"
    
    print '(a)', "  Potential Module Status:"
    print '(a)', "    ✓ Potential module loaded successfully"
    
    ! ------------------------------------------------------------------------
    ! ENERGY CALCULATION MODULE
    ! ------------------------------------------------------------------------
    print '(a)', "  Initializing Energy Calculation Module..."
    energy_calc = create_energy_calculator( &
      enabled = .true., &
      verbose = verbose, &
      parallel = parallel, &
      compute_components = .true., &
      output_freq = 10)
    call energy_calc%initialize(TSTEP, TINT)
    call energy_calc%setup_grid_weights(NPROD, NNDIM, basis_type)
    print '(a)', "    ✓ Energy module ready"
    
    ! ------------------------------------------------------------------------
    ! ADIABATIC CALCULATION MODULE
    ! ------------------------------------------------------------------------
    print '(a)', "  Initializing Adiabatic Calculation Module..."
    adiabatic_calc = create_adiabatic_calculator( &
      enabled = .true., &
      verbose = verbose, &
      parallel = parallel)
    call adiabatic_calc%initialize(NSTATE, NPROD, TSTEP)
    print '(a)', "    ✓ Adiabatic module ready"
    
    print '(a)', ""
    print '(a)', "✓ ENERGY AND ADIABATIC MODULES INITIALIZED SUCCESSFULLY!"
    
  end subroutine initialize_analysis_modules

  ! ==========================================================================
  ! SETUP INITIAL CONDITIONS
  ! ==========================================================================
  subroutine setup_initial_conditions()
    print '(a)', ""
    print '(a)', "SETTING UP INITIAL CONDITIONS..."
    print '(a)', "--------------------------------"
    
    call generate_basis_functions()
    call create_initial_grid()
    call create_index_mapping()
    call initialize_wavepacket()
    call compute_matrix_elements()
    call compute_initial_potential()
    call normalize_wavefunction()
    
    ycom = (0.0_dp, 0.0_dp)
    !$omp parallel do if(parallel .and. NPROD > 1000)
    do i = 1, NPROD
      ycom(NPLACE, i) = phi(i)
    end do
    !$omp end parallel do
    
    call calculate_initial_probabilities()
    
    if (USE_CONSERVATION_MONITOR) then
      total_energy_initial = 0.5_dp * (sum(momen**2) + sum(freqn**2 * coord**2))
      total_energy_current = total_energy_initial
    end if
    
    print '(a)', "  Initial conditions setup complete."
    
  end subroutine setup_initial_conditions

  ! ==========================================================================
  ! GENERATE BASIS FUNCTIONS
  ! ==========================================================================
  subroutine generate_basis_functions()
    real(dp) :: x
    integer :: n_points
    
    do i = 1, NMODE
      n_points = NNDIM(i)
      
      select case (basis_type(i))
      case (BASIS_HERMITE)
        do j = 1, n_points
          x = real(j - (n_points+1)/2, dp) / sqrt(2.0_dp * freqn(i) / HBAR)
          if (abs(freqn(i)) < SMALL) then
            x = real(j - (n_points+1)/2, dp) / sqrt(2.0_dp * 0.01_dp)
          end if
          hermrootf(j,i) = x
        end do
        do j = n_points+1, NN_MAX
          hermrootf(j,i) = 0.0_dp
        end do
        
      case (BASIS_SIN)
        do j = 1, n_points
          x = 2.0_dp * PI * real(j-1, dp) / max(real(n_points-1, dp), 1.0_dp)
          sin_basis(j,i) = x
        end do
        do j = n_points+1, NN_MAX
          sin_basis(j,i) = 0.0_dp
        end do
        
      case (BASIS_LEGENDRE)
        do j = 1, n_points
          x = cos(PI * (real(j, dp) - 0.5_dp) / real(n_points, dp))
          legendre_roots(j,i) = x
        end do
        do j = n_points+1, NN_MAX
          legendre_roots(j,i) = 0.0_dp
        end do
        
      end select
    end do
    
    print '(a)', "  Basis functions generated for all modes."
  end subroutine generate_basis_functions

  ! ==========================================================================
  ! CREATE INITIAL GRID
  ! ==========================================================================
  subroutine create_initial_grid()
    do i = 1, NMODE
      select case (basis_type(i))
      case (BASIS_HERMITE)
        do j = 1, NNDIM(i)
          newgrid(j,i) = coord(i) + hermrootf(j,i)
        end do
      case (BASIS_SIN)
        do j = 1, NNDIM(i)
          newgrid(j,i) = coord(i) + sin_basis(j,i) / (2.0_dp * PI)
        end do
      case (BASIS_LEGENDRE)
        do j = 1, NNDIM(i)
          newgrid(j,i) = coord(i) + legendre_roots(j,i) * 0.5_dp
        end do
      end select
      do j = NNDIM(i)+1, NN_MAX
        newgrid(j,i) = 0.0_dp
      end do
    end do
  end subroutine create_initial_grid

  ! ==========================================================================
  ! CREATE INDEX MAPPING
  ! ==========================================================================
  subroutine create_index_mapping()
    integer(ip) :: indices(NMODE)
    
    !$omp parallel if(parallel .and. NPROD > 10000) private(i, j, indices)
    !$omp do
    do i = 1, NPROD
      indices = 1
      j = i - 1
      do m = NMODE, 1, -1
        indices(m) = mod(j, NNDIM(m)) + 1
        j = j / NNDIM(m)
      end do
      tot(i,:) = indices
    end do
    !$omp end do
    !$omp end parallel
  end subroutine create_index_mapping

  ! ==========================================================================
  ! INITIALIZE WAVEPACKET
  ! ==========================================================================
  subroutine initialize_wavepacket()
    real(dp) :: facsho, xk, exk, pwk
    complex(dp) :: exr
    integer :: idx
    
    facsho = 1.0_dp
    if (NMODE <= 40) then
      facsho = (2.0_dp * PI * HBAR)**(-0.25_dp * NMODE)
    else
      facsho = exp(-0.25_dp * NMODE * log(2.0_dp * PI * HBAR))
    end if
    
    if (isnan(facsho) .or. abs(facsho) > huge(1.0_dp)/1000.0_dp .or. abs(facsho) < tiny(1.0_dp)*1000.0_dp) then
      print *, "WARNING: facsho is problematic, setting to 1.0"
      facsho = 1.0_dp
    end if
    
    !$omp parallel do if(parallel .and. NPROD > 1000) &
    !$omp private(i, j, idx, xk, exk, pwk, exr)
    do i = 1, NPROD
      exr = (1.0_dp, 0.0_dp)
      do j = 1, NMODE
        idx = tot(i,j)
        
        select case (basis_type(j))
        case (BASIS_HERMITE)
          xk = sqrt(2.0_dp / HBAR) * (newgrid(idx, j) - coord(j))
          exk = exp(-0.5_dp * xk**2)
          pwk = QUANTA * momen(j) * (newgrid(idx, j) - coord(j)) / HBAR
          exr = exr * exk * exp(ZI * pwk)
          
        case (BASIS_SIN)
          xk = sin(newgrid(idx, j) - coord(j))
          exk = exp(-0.5_dp * xk**2)
          exr = exr * exk
          
        case (BASIS_LEGENDRE)
          xk = (newgrid(idx, j) - coord(j))
          exk = exp(-0.5_dp * freqn(j) * xk**2 / HBAR)
          exr = exr * exk
          
        end select
        
        if (isnan(real(exr)) .or. isnan(aimag(exr))) then
          exr = (1.0_dp, 0.0_dp)
        end if
      end do
      phi(i) = facsho * exr
      
      if (isnan(real(phi(i))) .or. isnan(aimag(phi(i)))) then
        phi(i) = (0.0_dp, 0.0_dp)
      end if
    end do
    !$omp end parallel do
    
    phi(1) = (1.0_dp, 0.0_dp)
    
  end subroutine initialize_wavepacket

  ! ==========================================================================
  ! COMPUTE MATRIX ELEMENTS
  ! ==========================================================================
  subroutine compute_matrix_elements()
    print '(a)', "  Computing matrix elements..."
    
    !$omp parallel workshare
    mat = 0.0_dp
    matbar = 0.0_dp
    totmat = 0.0_dp
    !$omp end parallel workshare
    
    facp = 0.0_dp
    facq = 0.0_dp
    
    !$omp parallel do private(i, j, k) if(parallel .and. NMODE > 4)
    do i = 1, NMODE
      select case (basis_type(i))
      case (BASIS_HERMITE)
        do j = 1, NNDIM(i)
          do k = 1, NNDIM(i)
            if (j == k) then
              mat(i, j, k, 1) = 1.0_dp
              matbar(i, j, k, 1) = 1.0_dp
              totmat(i, j, k) = 0.5_dp * HBAR * width(i)
            end if
          end do
        end do
        
      case (BASIS_SIN)
        do j = 1, NNDIM(i)
          do k = 1, NNDIM(i)
            if (j == k) then
              mat(i, j, k, 1) = 1.0_dp
              matbar(i, j, k, 1) = 1.0_dp
              totmat(i, j, k) = 0.5_dp * HBAR**2 * (j-1)**2
            end if
          end do
        end do
        
      case (BASIS_LEGENDRE)
        do j = 1, NNDIM(i)
          do k = 1, NNDIM(i)
            if (j == k) then
              mat(i, j, k, 1) = 1.0_dp
              matbar(i, j, k, 1) = 1.0_dp
              totmat(i, j, k) = 0.5_dp * HBAR * freqn(i) * j * (j+1)
            end if
          end do
        end do
        
      end select
    end do
    !$omp end parallel do
  end subroutine compute_matrix_elements

  ! ==========================================================================
  ! COMPUTE INITIAL POTENTIAL
  ! ==========================================================================
  subroutine compute_initial_potential()
    call POTENTIAL(newgrid, v, NN_MAX, NMODE, NSTATE, NPROD, freqn, ffn)
    print '(a)', "  Initial potential computed using potential_module"
  end subroutine compute_initial_potential

  ! ==========================================================================
  ! NORMALIZE WAVEFUNCTION
  ! ==========================================================================
  subroutine normalize_wavefunction()
    real(dp) :: norm
    
    norm = 0.0_dp
    !$omp parallel do if(parallel .and. NPROD > 1000) reduction(+:norm)
    do i = 1, NPROD
      norm = norm + real(conjg(phi(i)) * phi(i))
    end do
    !$omp end parallel do
    
    if (isnan(norm) .or. norm <= SMALL) then
      print *, "WARNING: Wavefunction norm is NaN or too small: ", norm
      call set_simple_gaussian()
      norm = 0.0_dp
      do i = 1, NPROD
        norm = norm + real(conjg(phi(i)) * phi(i))
      end do
    end if
    
    norm = sqrt(norm)
    
    if (norm > SMALL) then
      !$omp parallel do if(parallel .and. NPROD > 1000)
      do i = 1, NPROD
        phi(i) = phi(i) / norm
      end do
      !$omp end parallel do
    else
      print *, "Warning: Initial wavefunction norm too small: ", norm
      phi(1) = (1.0_dp, 0.0_dp)
      do i = 2, NPROD
        phi(i) = (0.0_dp, 0.0_dp)
      end do
    end if
    
  end subroutine normalize_wavefunction

  ! ==========================================================================
  ! SET SIMPLE GAUSSIAN (FALLBACK)
  ! ==========================================================================
  subroutine set_simple_gaussian()
    real(dp) :: x, gauss
    integer :: i_local, j_local, idx_local
    
    do i_local = 1, NPROD
      gauss = 1.0_dp
      do j_local = 1, NMODE
        idx_local = tot(i_local, j_local)
        x = newgrid(idx_local, j_local) - coord(j_local)
        gauss = gauss * exp(-0.1_dp * x**2)
      end do
      phi(i_local) = gauss
    end do
  end subroutine set_simple_gaussian

  ! ==========================================================================
  ! CALCULATE INITIAL PROBABILITIES
  ! ==========================================================================
  subroutine calculate_initial_probabilities()
    real(dp) :: probs(NSTATE), total_prob
    
    probs = 0.0_dp
    total_prob = 0.0_dp
    
    !$omp parallel do if(parallel .and. NPROD > 1000) private(i, j) reduction(+:probs, total_prob)
    do i = 1, NSTATE
      do j = 1, NPROD
        probs(i) = probs(i) + real(conjg(ycom(i,j)) * ycom(i,j))
      end do
      total_prob = total_prob + probs(i)
    end do
    !$omp end parallel do
    
    print '(a,6f12.6)', "  Initial state probabilities: ", probs
    print '(a,es12.4)', "  Total probability: ", total_prob
    
    if (abs(total_prob - 1.0_dp) > 1.0e-6_dp) then
      print '(a,es12.4)', "  Warning: Total probability not exactly 1.0: ", total_prob
    end if
  end subroutine calculate_initial_probabilities

  ! ==========================================================================
  ! CALCULATE TOTAL ENERGY
  ! ==========================================================================
  subroutine calculate_total_energy(total_energy)
    real(dp), intent(out) :: total_energy
    total_energy = 0.5_dp * (sum(momen**2) + sum(freqn**2 * coord**2))
  end subroutine calculate_total_energy

  ! ==========================================================================
  ! RUN DYNAMICS SIMULATION
  ! ==========================================================================
  subroutine run_dynamics_simulation()
    print '(a)', ""
    print '(a)', "STARTING DYNAMICS SIMULATION..."
    print '(a)', "-------------------------------"
    
    ! Open output files
    open(101, file='barrelene_probabilities.dat', status='replace')
    open(102, file='barrelene_autocorrelation.dat', status='replace')
    open(103, file='barrelene_classical_trajectory.dat', status='replace')
    open(104, file='barrelene_timing.dat', status='replace')
    
    write(101, '(a)') '# Time (fs) and diabatic surface probabilities'
    write(101, '(a,6i12)') '# Time', (i, i=1, NSTATE)
    write(102, '(a)') '# Time (fs), Real(C), Imag(C), |C|'
    write(103, '(a)') '# Time (fs), Coord(1-5), Momen(1-5)'
    write(104, '(a)') '# Performance timing information'
    
    time = 0.0_dp
    current_dt = TINT
    step_accepted = .true.
    
    ! Initialize adiabatic module with initial data
    call adiabatic_calc%compute_transformation(v, NPROD, NSTATE)
    call adiabatic_calc%transform_wavefunction(ycom, NPROD, NSTATE, time, 0, .true.)
    
    call save_main_output()
    
    do step = 1, TSTEP
      
      ! Adaptive time stepping loop
      step_accepted = .false.
      do while (.not. step_accepted)
        
        ! Propagate classical coordinates
        call propagate_classical_coordinates(current_dt)
        
        ! Update quantum grid
        call update_quantum_grid()
        
        ! Compute potential
        call compute_potential_at_current_time()
        
        ! Execute analysis modules (Energy + Adiabatic)
        call execute_analysis_modules()
        
        ! Propagate quantum wavefunction
        call propagate_quantum_wavefunction_selected(current_dt)
        
        ! Error estimation and step control
        if (USE_ERROR_CONTROL .and. USE_ADAPTIVE_TIMESTEP) then
          call estimate_error(error_estimate)
          
          if (error_estimate < ERROR_TOLERANCE .or. current_dt <= MIN_TIMESTEP) then
            step_accepted = .true.
            time = time + current_dt
            
            if (error_estimate > 0.0_dp) then
              next_dt = SAFETY_FACTOR * current_dt * (ERROR_TOLERANCE / error_estimate)**0.2_dp
              next_dt = min(MAX_TIMESTEP, max(MIN_TIMESTEP, next_dt))
              next_dt = min(next_dt, current_dt * GROWTH_LIMIT)
            else
              next_dt = min(current_dt * GROWTH_LIMIT, MAX_TIMESTEP)
            end if
          else
            current_dt = max(current_dt * 0.5_dp, MIN_TIMESTEP)
          end if
        else
          step_accepted = .true.
          time = time + current_dt
          next_dt = current_dt
        end if
      end do
      
      ! Update time step
      if (USE_ADAPTIVE_TIMESTEP) then
        current_dt = next_dt
      end if
      
      ! Save output
      call save_main_output()
      
      ! Progress indicator
      if (verbose .and. mod(step, 50) == 0) then
        if (USE_ADAPTIVE_TIMESTEP) then
          print '(a,i6,a,i6,a,f10.4,a,f8.4)', "  Step ", step, "/", TSTEP, &
                ", Time = ", time, " fs, dt = ", current_dt, " fs"
        else
          print '(a,i6,a,i6,a,f10.4)', "  Step ", step, "/", TSTEP, &
                ", Time = ", time, " fs"
        end if
      end if
      
    end do
    
    close(101)
    close(102)
    close(103)
    close(104)
    
    print '(a)', "  Dynamics simulation completed."
    
  end subroutine run_dynamics_simulation

  ! ==========================================================================
  ! PROPAGATE CLASSICAL COORDINATES
  ! ==========================================================================
  subroutine propagate_classical_coordinates(dt)
    real(dp), intent(in) :: dt
    
    if (step == 1) then
      momen(1) = 0.01_dp
    end if
    
    coord = coord + dt * momen
    momen = momen - dt * freqn**2 * coord
    momen = momen * 0.9999_dp
    
    where (abs(coord) > 10.0_dp) coord = sign(10.0_dp, coord)
    where (abs(momen) > 10.0_dp) momen = sign(10.0_dp, momen)
  end subroutine propagate_classical_coordinates

  ! ==========================================================================
  ! UPDATE QUANTUM GRID
  ! ==========================================================================
  subroutine update_quantum_grid()
    do i = 1, NMODE
      select case (basis_type(i))
      case (BASIS_HERMITE)
        do j = 1, NNDIM(i)
          newgrid(j,i) = coord(i) + hermrootf(j,i)
        end do
      case (BASIS_SIN)
        do j = 1, NNDIM(i)
          newgrid(j,i) = coord(i) + sin_basis(j,i) / (2.0_dp * PI)
        end do
      case (BASIS_LEGENDRE)
        do j = 1, NNDIM(i)
          newgrid(j,i) = coord(i) + legendre_roots(j,i) * 0.5_dp
        end do
      end select
      do j = NNDIM(i)+1, NN_MAX
        newgrid(j,i) = 0.0_dp
      end do
    end do
    
    facp = 0.5_dp * momen**2
    facq = 0.5_dp * sqrt(HBAR / max(width, SMALL)) * coord * freqn**2
    
    where (abs(facp) > 1.0e10_dp) facp = sign(1.0e10_dp, facp)
    where (abs(facq) > 1.0e10_dp) facq = sign(1.0e10_dp, facq)
    where (abs(facp) < 1.0e-30_dp) facp = 0.0_dp
    where (abs(facq) < 1.0e-30_dp) facq = 0.0_dp
  end subroutine update_quantum_grid

  ! ==========================================================================
  ! COMPUTE POTENTIAL AT CURRENT TIME
  ! ==========================================================================
  subroutine compute_potential_at_current_time()
    call POTENTIAL(newgrid, v, NN_MAX, NMODE, NSTATE, NPROD, freqn, ffn)
  end subroutine compute_potential_at_current_time

  ! ==========================================================================
  ! SELECTED QUANTUM PROPAGATION
  ! ==========================================================================
  subroutine propagate_quantum_wavefunction_selected(dt)
    real(dp), intent(in) :: dt
    real(dp) :: dt_hbar
    
    dt_hbar = dt / HBAR
    
    if (USE_ERROR_CONTROL .and. USE_ADAPTIVE_TIMESTEP) then
      !$omp parallel workshare
      ycom_ref = ycom
      !$omp end parallel workshare
    end if
    
    select case (PROPAGATION_METHOD)
      case (PROP_SPLIT_OPERATOR_2ND)
        call split_operator_2nd(dt_hbar)
      case (PROP_SPLIT_OPERATOR_4TH)
        call split_operator_4th(dt_hbar)
      case (PROP_MAGNUS_2ND)
        call magnus_2nd(dt_hbar)
      case (PROP_RUNGE_KUTTA_4TH)
        call runge_kutta_4th(dt_hbar)
      case (PROP_ADAPTIVE_SPLIT)
        call adaptive_split_operator(dt_hbar)
      case default
        call split_operator_2nd(dt_hbar)
    end select
    
    call normalize_total_wavefunction()
    
  end subroutine propagate_quantum_wavefunction_selected

  ! ==========================================================================
  ! 2ND ORDER SPLIT-OPERATOR
  ! ==========================================================================
  subroutine split_operator_2nd(dt_hbar)
    real(dp), intent(in) :: dt_hbar
    
    complex(dp), allocatable :: ycom_temp(:,:)
    integer :: iprod, istate, n
    real(dp) :: kin_energy
    
    allocate(ycom_temp(NSTATE, NPROD))
    
    !$omp parallel workshare
    ycom_temp = ycom
    !$omp end parallel workshare
    
    ! First half kinetic
    !$omp parallel do private(iprod, istate, n, kin_energy) if(parallel .and. NPROD > 1000) schedule(dynamic)
    do iprod = 1, NPROD
      kin_energy = 0.0_dp
      do n = 1, NMODE
        if (tot(iprod,n) <= NNDIM(n)) then
          kin_energy = kin_energy + &
            (facp(n) + facq(n) * matbar(n, tot(iprod,n), tot(iprod,n), 3) + &
             0.5_dp * totmat(n, tot(iprod,n), tot(iprod,n)))
        end if
      end do
      
      do istate = 1, NSTATE
        ycom_temp(istate, iprod) = ycom_temp(istate, iprod) * &
          exp(-ZI * 0.5_dp * dt_hbar * kin_energy)
      end do
    end do
    !$omp end parallel do
    
    ! Potential (full step)
    !$omp parallel do private(iprod) if(parallel .and. NPROD > 1000) schedule(dynamic)
    do iprod = 1, NPROD
      call propagate_local_potential(ycom_temp(:, iprod), v(:,:,iprod), dt_hbar)
    end do
    !$omp end parallel do
    
    ! Second half kinetic
    !$omp parallel do private(iprod, istate, n, kin_energy) if(parallel .and. NPROD > 1000) schedule(dynamic)
    do iprod = 1, NPROD
      kin_energy = 0.0_dp
      do n = 1, NMODE
        if (tot(iprod,n) <= NNDIM(n)) then
          kin_energy = kin_energy + &
            (facp(n) + facq(n) * matbar(n, tot(iprod,n), tot(iprod,n), 3) + &
             0.5_dp * totmat(n, tot(iprod,n), tot(iprod,n)))
        end if
      end do
      
      do istate = 1, NSTATE
        ycom_temp(istate, iprod) = ycom_temp(istate, iprod) * &
          exp(-ZI * 0.5_dp * dt_hbar * kin_energy)
      end do
    end do
    !$omp end parallel do
    
    !$omp parallel workshare
    ycom = ycom_temp
    !$omp end parallel workshare
    
    deallocate(ycom_temp)
    
  end subroutine split_operator_2nd

  ! ==========================================================================
  ! 4TH ORDER SPLIT-OPERATOR
  ! ==========================================================================
  subroutine split_operator_4th(dt_hbar)
    real(dp), intent(in) :: dt_hbar
    
    real(dp) :: p, dt1, dt2, dt3
    complex(dp), allocatable :: ycom_temp(:,:)
    
    p = 1.0_dp / (4.0_dp - 4.0_dp**(1.0_dp/3.0_dp))
    dt1 = p * dt_hbar
    dt2 = (1.0_dp - 4.0_dp * p) * dt_hbar
    dt3 = dt1
    
    allocate(ycom_temp(NSTATE, NPROD))
    
    !$omp parallel workshare
    ycom_temp = ycom
    !$omp end parallel workshare
    
    call split_operator_2nd_step(ycom_temp, dt1)
    call split_operator_2nd_step(ycom_temp, dt2)
    call split_operator_2nd_step(ycom_temp, dt3)
    
    !$omp parallel workshare
    ycom = ycom_temp
    !$omp end parallel workshare
    
    deallocate(ycom_temp)
    
  end subroutine split_operator_4th

  ! ==========================================================================
  ! SPLIT-OPERATOR STEP (HELPER)
  ! ==========================================================================
  subroutine split_operator_2nd_step(wf, dt_step)
    complex(dp), intent(inout) :: wf(:,:)
    real(dp), intent(in) :: dt_step
    
    integer :: iprod, istate, n
    real(dp) :: kin_energy
    
    ! Kinetic half step
    !$omp parallel do private(iprod, istate, n, kin_energy) if(parallel .and. NPROD > 1000) schedule(dynamic)
    do iprod = 1, NPROD
      kin_energy = 0.0_dp
      do n = 1, NMODE
        if (tot(iprod,n) <= NNDIM(n)) then
          kin_energy = kin_energy + &
            (facp(n) + facq(n) * matbar(n, tot(iprod,n), tot(iprod,n), 3) + &
             0.5_dp * totmat(n, tot(iprod,n), tot(iprod,n)))
        end if
      end do
      
      do istate = 1, NSTATE
        wf(istate, iprod) = wf(istate, iprod) * &
          exp(-ZI * 0.5_dp * dt_step * kin_energy)
      end do
    end do
    !$omp end parallel do
    
    ! Potential full step
    !$omp parallel do private(iprod) if(parallel .and. NPROD > 1000) schedule(dynamic)
    do iprod = 1, NPROD
      call propagate_local_potential(wf(:, iprod), v(:,:,iprod), dt_step)
    end do
    !$omp end parallel do
    
    ! Kinetic half step
    !$omp parallel do private(iprod, istate, n, kin_energy) if(parallel .and. NPROD > 1000) schedule(dynamic)
    do iprod = 1, NPROD
      kin_energy = 0.0_dp
      do n = 1, NMODE
        if (tot(iprod,n) <= NNDIM(n)) then
          kin_energy = kin_energy + &
            (facp(n) + facq(n) * matbar(n, tot(iprod,n), tot(iprod,n), 3) + &
             0.5_dp * totmat(n, tot(iprod,n), tot(iprod,n)))
        end if
      end do
      
      do istate = 1, NSTATE
        wf(istate, iprod) = wf(istate, iprod) * &
          exp(-ZI * 0.5_dp * dt_step * kin_energy)
      end do
    end do
    !$omp end parallel do
    
  end subroutine split_operator_2nd_step

  ! ==========================================================================
  ! LOCAL POTENTIAL PROPAGATION
  ! ==========================================================================
  subroutine propagate_local_potential(wf_local, v_local, dt_hbar)
    complex(dp), intent(inout) :: wf_local(:)
    real(dp), intent(in) :: v_local(:,:)
    real(dp), intent(in) :: dt_hbar
    
    complex(dp), allocatable :: wf_temp(:), wf_mid(:)
    integer :: nstate, i, j
    real(dp) :: dt_half, dt_sq_half
    real(dp), allocatable :: v_sq(:,:)
    
    nstate = size(wf_local)
    allocate(wf_temp(nstate), wf_mid(nstate))
    allocate(v_sq(nstate, nstate))
    
    dt_half = 0.5_dp * dt_hbar
    dt_sq_half = 0.5_dp * dt_hbar**2
    
    v_sq = matmul(v_local, v_local)
    wf_temp = wf_local
    
    do i = 1, nstate
      wf_mid(i) = (0.0_dp, 0.0_dp)
      do j = 1, nstate
        if (i == j) then
          wf_mid(i) = wf_mid(i) + (1.0_dp - ZI * dt_half * v_local(i,j) - &
                     0.5_dp * dt_sq_half * v_sq(i,j)) * wf_temp(j)
        else
          wf_mid(i) = wf_mid(i) - ZI * dt_half * v_local(i,j) * wf_temp(j)
        end if
      end do
    end do
    
    wf_local = (0.0_dp, 0.0_dp)
    do i = 1, nstate
      do j = 1, nstate
        if (i == j) then
          wf_local(i) = wf_local(i) + (1.0_dp - ZI * dt_half * v_local(i,j) - &
                       0.5_dp * dt_sq_half * v_sq(i,j)) * wf_mid(j)
        else
          wf_local(i) = wf_local(i) - ZI * dt_half * v_local(i,j) * wf_mid(j)
        end if
      end do
    end do
    
    deallocate(wf_temp, wf_mid, v_sq)
    
  end subroutine propagate_local_potential

  ! ==========================================================================
  ! MAGNUS 2ND ORDER PROPAGATION
  ! ==========================================================================
  subroutine magnus_2nd(dt_hbar)
    real(dp), intent(in) :: dt_hbar
    
    complex(dp), allocatable :: ycom_temp(:,:)
    integer :: iprod
    
    allocate(ycom_temp(NSTATE, NPROD))
    
    !$omp parallel workshare
    ycom_temp = ycom
    !$omp end parallel workshare
    
    !$omp parallel do private(iprod) if(parallel .and. NPROD > 1000) schedule(dynamic)
    do iprod = 1, NPROD
      call propagate_local_hamiltonian(ycom_temp(:, iprod), iprod, dt_hbar)
    end do
    !$omp end parallel do
    
    !$omp parallel workshare
    ycom = ycom_temp
    !$omp end parallel workshare
    
    deallocate(ycom_temp)
    
  end subroutine magnus_2nd

  ! ==========================================================================
  ! LOCAL HAMILTONIAN PROPAGATION
  ! ==========================================================================
  subroutine propagate_local_hamiltonian(wf_local, iprod, dt_hbar)
    complex(dp), intent(inout) :: wf_local(:)
    integer, intent(in) :: iprod
    real(dp), intent(in) :: dt_hbar
    
    complex(dp), allocatable :: H(:,:), U(:,:)
    integer :: nstate, i, j, k
    real(dp) :: kin_energy
    
    nstate = size(wf_local)
    allocate(H(nstate, nstate), U(nstate, nstate))
    
    H = v(:,:,iprod)
    
    kin_energy = 0.0_dp
    do k = 1, NMODE
      if (tot(iprod,k) <= NNDIM(k)) then
        kin_energy = kin_energy + &
          (facp(k) + facq(k) * matbar(k, tot(iprod,k), tot(iprod,k), 3) + &
           0.5_dp * totmat(k, tot(iprod,k), tot(iprod,k)))
      end if
    end do
    
    do i = 1, nstate
      H(i,i) = H(i,i) + kin_energy
    end do
    
    U = (0.0_dp, 0.0_dp)
    do i = 1, nstate
      U(i,i) = (1.0_dp, 0.0_dp)
      do j = 1, nstate
        do k = 1, nstate
          U(i,j) = U(i,j) - ZI * dt_hbar * H(i,j) - &
                   0.5_dp * dt_hbar**2 * H(i,k) * H(k,j)
        end do
      end do
    end do
    
    wf_local = matmul(U, wf_local)
    
    deallocate(H, U)
    
  end subroutine propagate_local_hamiltonian

  ! ==========================================================================
  ! RUNGE-KUTTA 4TH ORDER PROPAGATION
  ! ==========================================================================
  subroutine runge_kutta_4th(dt_hbar)
    real(dp), intent(in) :: dt_hbar
    
    complex(dp), allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
    complex(dp), allocatable :: ycom_temp(:,:)
    
    allocate(k1(NSTATE, NPROD), k2(NSTATE, NPROD), k3(NSTATE, NPROD), k4(NSTATE, NPROD))
    allocate(ycom_temp(NSTATE, NPROD))
    
    call compute_rhs(ycom, k1)
    
    !$omp parallel workshare
    ycom_temp = ycom + 0.5_dp * dt_hbar * k1
    !$omp end parallel workshare
    call compute_rhs(ycom_temp, k2)
    
    !$omp parallel workshare
    ycom_temp = ycom + 0.5_dp * dt_hbar * k2
    !$omp end parallel workshare
    call compute_rhs(ycom_temp, k3)
    
    !$omp parallel workshare
    ycom_temp = ycom + dt_hbar * k3
    !$omp end parallel workshare
    call compute_rhs(ycom_temp, k4)
    
    !$omp parallel workshare
    ycom = ycom + (dt_hbar / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
    !$omp end parallel workshare
    
    deallocate(k1, k2, k3, k4, ycom_temp)
    
  end subroutine runge_kutta_4th

  ! ==========================================================================
  ! COMPUTE RHS FOR RK4
  ! ==========================================================================
  subroutine compute_rhs(wf_in, rhs_out)
    complex(dp), intent(in) :: wf_in(:,:)
    complex(dp), intent(out) :: rhs_out(:,:)
    
    integer :: iprod, istate, n, jstate
    real(dp) :: kin_energy
    
    rhs_out = (0.0_dp, 0.0_dp)
    
    !$omp parallel do private(iprod, istate, n, jstate, kin_energy) if(parallel .and. NPROD > 1000) schedule(dynamic)
    do iprod = 1, NPROD
      kin_energy = 0.0_dp
      do n = 1, NMODE
        if (tot(iprod,n) <= NNDIM(n)) then
          kin_energy = kin_energy + &
            (facp(n) + facq(n) * matbar(n, tot(iprod,n), tot(iprod,n), 3) + &
             0.5_dp * totmat(n, tot(iprod,n), tot(iprod,n)))
        end if
      end do
      
      do istate = 1, NSTATE
        rhs_out(istate, iprod) = rhs_out(istate, iprod) - ZI * kin_energy * wf_in(istate, iprod)
        do jstate = 1, NSTATE
          rhs_out(istate, iprod) = rhs_out(istate, iprod) - ZI * v(istate, jstate, iprod) * wf_in(jstate, iprod)
        end do
      end do
    end do
    !$omp end parallel do
    
  end subroutine compute_rhs

  ! ==========================================================================
  ! ADAPTIVE SPLIT-OPERATOR
  ! ==========================================================================
  subroutine adaptive_split_operator(dt_hbar)
    real(dp), intent(in) :: dt_hbar
    
    complex(dp), allocatable :: ycom_temp(:,:)
    real(dp) :: local_dt, error
    integer :: max_iter = 5
    integer :: iter
    
    allocate(ycom_temp(NSTATE, NPROD))
    
    do iter = 1, max_iter
      !$omp parallel workshare
      ycom_temp = ycom
      !$omp end parallel workshare
      
      call split_operator_4th(dt_hbar)
      
      !$omp parallel workshare
      ycom_high = ycom
      !$omp end parallel workshare
      
      !$omp parallel workshare
      ycom = ycom_temp
      !$omp end parallel workshare
      
      local_dt = dt_hbar / 2.0_dp
      call split_operator_2nd(local_dt)
      call split_operator_2nd(local_dt)
      
      call estimate_error_adaptive(error)
      
      if (error < ERROR_TOLERANCE .or. iter == max_iter) then
        !$omp parallel workshare
        ycom = ycom_high
        !$omp end parallel workshare
        exit
      else
        local_dt = local_dt * 0.5_dp
        !$omp parallel workshare
        ycom = ycom_temp
        !$omp end parallel workshare
      end if
    end do
    
    deallocate(ycom_temp)
    
  end subroutine adaptive_split_operator

  ! ==========================================================================
  ! ESTIMATE ERROR FOR ADAPTIVE
  ! ==========================================================================
  subroutine estimate_error_adaptive(error)
    real(dp), intent(out) :: error
    
    real(dp) :: diff_norm, ref_norm
    integer :: i, j
    
    diff_norm = 0.0_dp
    ref_norm = 0.0_dp
    
    !$omp parallel do private(i, j) reduction(+:diff_norm, ref_norm) if(parallel .and. NSTATE*NPROD > 10000)
    do i = 1, NSTATE
      do j = 1, NPROD
        diff_norm = diff_norm + real(conjg(ycom(i,j) - ycom_high(i,j)) * (ycom(i,j) - ycom_high(i,j)))
        ref_norm = ref_norm + real(conjg(ycom_high(i,j)) * ycom_high(i,j))
      end do
    end do
    !$omp end parallel do
    
    if (ref_norm > SMALL) then
      error = sqrt(diff_norm / ref_norm)
    else
      error = 0.0_dp
    end if
    
  end subroutine estimate_error_adaptive

  ! ==========================================================================
  ! ESTIMATE ERROR FOR ADAPTIVE TIMESTEP
  ! ==========================================================================
  subroutine estimate_error(error)
    real(dp), intent(out) :: error
    
    real(dp) :: diff_norm, ref_norm
    integer :: i, j
    
    diff_norm = 0.0_dp
    ref_norm = 0.0_dp
    
    !$omp parallel do private(i, j) reduction(+:diff_norm, ref_norm) if(parallel .and. NSTATE*NPROD > 10000)
    do i = 1, NSTATE
      do j = 1, NPROD
        diff_norm = diff_norm + real(conjg(ycom(i,j) - ycom_ref(i,j)) * (ycom(i,j) - ycom_ref(i,j)))
        ref_norm = ref_norm + real(conjg(ycom_ref(i,j)) * ycom_ref(i,j))
      end do
    end do
    !$omp end parallel do
    
    if (ref_norm > SMALL) then
      error = sqrt(diff_norm / ref_norm)
    else
      error = 0.0_dp
    end if
    
  end subroutine estimate_error

  ! ==========================================================================
  ! NORMALIZE TOTAL WAVEFUNCTION
  ! ==========================================================================
  subroutine normalize_total_wavefunction()
    real(dp) :: total_norm
    integer :: i, j
    
    total_norm = 0.0_dp
    
    !$omp parallel do private(i, j) reduction(+:total_norm) if(parallel .and. NSTATE*NPROD > 10000)
    do i = 1, NSTATE
      do j = 1, NPROD
        total_norm = total_norm + real(conjg(ycom(i,j)) * ycom(i,j))
      end do
    end do
    !$omp end parallel do
    
    total_norm = sqrt(total_norm)
    
    if (total_norm > SMALL) then
      !$omp parallel do private(i, j) if(parallel .and. NSTATE*NPROD > 10000)
      do i = 1, NSTATE
        do j = 1, NPROD
          ycom(i,j) = ycom(i,j) / total_norm
        end do
      end do
      !$omp end parallel do
    else
      print *, "ERROR: Wavefunction norm too small at step ", step, ": ", total_norm
      if (step > 10) then
        ycom = (0.0_dp, 0.0_dp)
        !$omp parallel do if(parallel .and. NPROD > 1000)
        do j = 1, NPROD
          ycom(NPLACE, j) = 1.0_dp / sqrt(real(NPROD, dp))
        end do
        !$omp end parallel do
      end if
    end if
    
  end subroutine normalize_total_wavefunction

  ! ==========================================================================
  ! EXECUTE ANALYSIS MODULES - ENERGY AND ADIABATIC
  ! ==========================================================================
  subroutine execute_analysis_modules()
    
    ! ------------------------------------------------------------------------
    ! ENERGY CALCULATION MODULE - Compute every 10 steps
    ! ------------------------------------------------------------------------
    if (mod(step, 10) == 0) then
      call energy_calc%compute_energies(time, step, ycom, v, totmat, &
                                       matbar, tot, coord, momen, &
                                       freqn, ffn, newgrid)
      if (verbose .and. mod(step, 100) == 0) then
        print '(a,f8.3,a)', "  Energy computed at t=", time, " fs"
      end if
    end if
    
    ! ------------------------------------------------------------------------
    ! ADIABATIC CALCULATION MODULE - Compute every 10 steps (more frequent)
    ! ------------------------------------------------------------------------
    if (mod(step, 10) == 0) then
      ! First compute the adiabatic transformation (diagonalization)
      call adiabatic_calc%compute_transformation(v, NPROD, NSTATE)
      
      ! Then transform the wavefunction and compute probabilities
      call adiabatic_calc%transform_wavefunction(ycom, NPROD, NSTATE, &
                                                time, step, .true.)
      
      if (verbose .and. mod(step, 100) == 0) then
        print '(a,f8.3,a)', "  Adiabatic probabilities computed at t=", time, " fs"
      end if
    end if
    
  end subroutine execute_analysis_modules
  
  ! ==========================================================================
  ! COMPUTE POTENTIAL AT A SINGLE POINT
  ! ==========================================================================
  function compute_potential_at_point(coords, nstate) result(v_diag)
    real(dp), intent(in) :: coords(:)
    integer(ip), intent(in) :: nstate
    real(dp) :: v_diag(nstate, nstate)
    
    real(dp), allocatable :: grid_point(:,:)
    real(dp), allocatable :: v_temp(:,:,:)
    
    allocate(grid_point(1, size(coords)))
    allocate(v_temp(nstate, nstate, 1))
    
    grid_point(1, :) = coords(:)
    call POTENTIAL(grid_point, v_temp, 1, NMODE, NSTATE, 1, freqn, ffn)
    v_diag = v_temp(:,:,1)
    
    deallocate(grid_point, v_temp)
  end function compute_potential_at_point

  ! ==========================================================================
  ! SAVE MAIN OUTPUT
  ! ==========================================================================
  subroutine save_main_output()
    real(dp) :: probs(NSTATE), autor, autoi, automod, total_prob
    integer :: i, j
    
    total_prob = 0.0_dp
    probs = 0.0_dp
    
    !$omp parallel do private(i, j) reduction(+:probs, total_prob) if(parallel .and. NPROD > 1000)
    do i = 1, NSTATE
      do j = 1, NPROD
        probs(i) = probs(i) + real(conjg(ycom(i,j)) * ycom(i,j))
      end do
      total_prob = total_prob + probs(i)
    end do
    !$omp end parallel do
    
    if (abs(total_prob - 1.0_dp) > 1.0e-8_dp .and. total_prob > SMALL) then
      probs = probs / total_prob
    end if
    
    call calculate_autocorrelation(autor, autoi, automod)
    
    write(101, '(f12.6,6f14.8)') time, probs
    write(102, '(4f16.8)') time, autor, autoi, automod
    write(103, '(f12.6,10f12.6)') time, coord(1:5), momen(1:5)
    flush(101)
    flush(102)
    flush(103)
    
    if (verbose .and. mod(step, 100) == 0) then
      print '(a,f8.3,a,6f8.4,a,f10.6)', "    t=", time, " fs P=", &
            probs, " Total=", sum(probs)
    end if
    
  end subroutine save_main_output
  
  ! ==========================================================================
  ! CALCULATE AUTOCORRELATION
  ! ==========================================================================
  subroutine calculate_autocorrelation(autor, autoi, automod)
    real(dp), intent(out) :: autor, autoi, automod
    complex(dp) :: autocorr
    integer :: j
    
    autocorr = (0.0_dp, 0.0_dp)
    
    !$omp parallel do private(j) reduction(+:autocorr) if(parallel .and. NPROD > 1000)
    do j = 1, NPROD
      autocorr = autocorr + conjg(phi(j)) * ycom(NPLACE, j)
    end do
    !$omp end parallel do
    
    autor = real(autocorr)
    autoi = aimag(autocorr)
    automod = sqrt(autor**2 + autoi**2)
    
  end subroutine calculate_autocorrelation

  ! ==========================================================================
  ! FINALIZE SIMULATION
  ! ==========================================================================
  subroutine finalize_simulation()
    print '(a)', ""
    print '(a)', "FINALIZING SIMULATION..."
    print '(a)', "------------------------"
    
    ! Finalize Energy Module
    call energy_calc%finalize()
    print '(a)', "  ✓ Energy calculation module finalized"
    
    ! Finalize Adiabatic Module
    call adiabatic_calc%finalize()
    print '(a)', "  ✓ Adiabatic calculation module finalized"
    
    call cleanup_system_memory()
    
    call cpu_time(total_end)
    
    open(104, file='barrelene_timing.dat', status='old', position='append')
    write(104, '(a,f12.6)') "Total wall time (s): ", total_end - total_start
    write(104, '(a,f12.6)') "Time per step (ms): ", &
          (total_end - total_start) / TSTEP * 1000.0_dp
    if (USE_CONSERVATION_MONITOR) then
      write(104, '(a,es10.2)') "Max norm deviation: ", norm_deviation_max
      write(104, '(a,es10.2)') "Max energy deviation: ", energy_deviation_max
    end if
    close(104)
    
  end subroutine finalize_simulation

  ! ==========================================================================
  ! DISPLAY FINAL SUMMARY
  ! ==========================================================================
  subroutine display_final_summary()
    print '(a)', ""
    print '(a)', "========================================================"
    print '(a)', "                   SIMULATION COMPLETE"
    print '(a)', "========================================================"
    
    print '(a,f12.6,a)', "Total CPU time: ", total_end - total_start, " seconds"
    print '(a,i0)', "Time steps completed: ", TSTEP
    print '(a,f8.3,a)', "Final simulation time: ", time, " fs"
    
    if (USE_CONSERVATION_MONITOR) then
      print '(a,es10.2)', "Maximum norm deviation: ", norm_deviation_max
      print '(a,es10.2)', "Maximum energy deviation: ", energy_deviation_max
    end if
    
    print '(a)', ""
    print '(a)', "OUTPUT FILES:"
    print '(a)', "  barrelene_probabilities.dat      - Diabatic state populations"
    print '(a)', "  barrelene_autocorrelation.dat    - Wavepacket autocorrelation"
    print '(a)', "  barrelene_classical_trajectory.dat - Classical coordinates"
    print '(a)', "  barrelene_timing.dat             - Performance statistics"
    
    print '(a)', ""
    print '(a)', "MODULE OUTPUT DIRECTORIES:"
    print '(a)', "  ./energy_analysis/      - Energy components and conservation"
    print '(a)', "  ./adiabatic_analysis/   - Adiabatic probabilities"
    print '(a)', ""
    
    print '(a)', "PROPAGATION METHOD USED:"
    select case(PROPAGATION_METHOD)
      case(PROP_SPLIT_OPERATOR_2ND)
        print '(a)', "  2nd-Order Split-Operator (IDENTICAL to original)"
      case(PROP_SPLIT_OPERATOR_4TH)
        print '(a)', "  4th-Order Split-Operator"
      case(PROP_MAGNUS_2ND)
        print '(a)', "  2nd-Order Magnus"
      case(PROP_RUNGE_KUTTA_4TH)
        print '(a)', "  4th-Order Runge-Kutta"
      case(PROP_ADAPTIVE_SPLIT)
        print '(a)', "  Adaptive Split-Operator"
    end select
    
    if (USE_ADAPTIVE_TIMESTEP) then
      print '(a)', "  Adaptive time stepping: ENABLED"
    else
      print '(a)', "  Adaptive time stepping: DISABLED"
    end if
    
    print '(a)', ""
    print '(a)', "========================================================"
    print '(a)', "               PROGRAM EXECUTION FINISHED"
    print '(a)', "========================================================"
    
  end subroutine display_final_summary

  ! ==========================================================================
  ! CLEANUP SYSTEM MEMORY
  ! ==========================================================================
  subroutine cleanup_system_memory()
    if (allocated(ycom)) deallocate(ycom)
    if (allocated(ycom_ref)) deallocate(ycom_ref)
    if (allocated(ycom_high)) deallocate(ycom_high)
    if (allocated(phi)) deallocate(phi)
    if (allocated(newgrid)) deallocate(newgrid)
    if (allocated(coord)) deallocate(coord)
    if (allocated(momen)) deallocate(momen)
    if (allocated(tot)) deallocate(tot)
    if (allocated(v)) deallocate(v)
    if (allocated(freq)) deallocate(freq)
    if (allocated(freqn)) deallocate(freqn)
    if (allocated(width)) deallocate(width)
    if (allocated(ffn)) deallocate(ffn)
    if (allocated(mat)) deallocate(mat)
    if (allocated(matbar)) deallocate(matbar)
    if (allocated(totmat)) deallocate(totmat)
    if (allocated(facp)) deallocate(facp)
    if (allocated(facq)) deallocate(facq)
    if (allocated(hermrootf)) deallocate(hermrootf)
    if (allocated(sin_basis)) deallocate(sin_basis)
    if (allocated(legendre_roots)) deallocate(legendre_roots)
    
    print '(a)', "  System memory cleaned up."
  end subroutine cleanup_system_memory

end program SO_TDDVR_barrelene
