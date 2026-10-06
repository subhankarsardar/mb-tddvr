! ============================================================================
! ADIABATIC PROBABILITIES CALCULATION MODULE
! Compatible with SO-TDDVR.f90
! ============================================================================

module adiabatic_module
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  
  private
  
  ! Public types and procedures
  public :: adiabatic_calculator, dp, ip
  public :: create_adiabatic_calculator
  
  ! Kind parameters
  integer, parameter :: dp = real64
  integer, parameter :: ip = int32
  
  ! Constants
  real(dp), parameter :: SMALL = 1.0e-12_dp
  real(dp), parameter :: PI = 3.14159265358979323846_dp
  
  ! Type for adiabatic calculator
  type :: adiabatic_calculator
    private
    logical :: enabled = .false.
    logical :: verbose = .true.
    logical :: parallel = .true.
    logical :: store_probabilities = .true.
    
    ! Storage for adiabatic data
    real(dp), allocatable :: energies(:,:)        ! (nstate, nnprod)
    real(dp), allocatable :: transformation(:,:,:) ! (nstate, nstate, nnprod)
    real(dp), allocatable :: probabilities(:,:)   ! (nstate, tstep)
    real(dp), allocatable :: times(:)             ! (tstep)
    
    ! Output file units - use units that don't conflict with main program
    ! (101-104)
    integer :: file_unit_prob = 300
    integer :: file_unit_energies = 301
    integer :: file_unit_trans = 302
    
    ! Parameters
    integer(ip) :: nstate = 0
    integer(ip) :: nnprod = 0
    integer(ip) :: tstep = 0
    integer(ip) :: current_step = 0
    
    ! Output directory
    character(len=100) :: output_dir = "./adiabatic_analysis/"
    
    ! Timing
    real(dp) :: compute_time = 0.0_dp
    logical :: transformation_computed = .false.
    
  contains
    procedure :: initialize => init_adiabatic
    procedure :: compute_transformation
    procedure :: transform_wavefunction
    procedure :: compute_adiabatic_populations
    procedure :: write_probabilities
    procedure :: write_adiabatic_energies
    procedure :: finalize => cleanup_adiabatic
    procedure :: is_enabled
    procedure :: enable
    procedure :: disable
    procedure :: set_output_dir
    procedure :: get_adiabatic_energies
    procedure :: get_transformation
  end type adiabatic_calculator
  
  ! Interface for constructor
  interface adiabatic_calculator
    module procedure create_calculator
  end interface

contains

  ! ============================================================================
  ! CONSTRUCTOR
  ! ============================================================================
  function create_calculator(enabled, verbose, parallel) result(calc)
    logical, intent(in), optional :: enabled, verbose, parallel
    type(adiabatic_calculator) :: calc
    
    if (present(enabled)) calc%enabled = enabled
    if (present(verbose)) calc%verbose = verbose
    if (present(parallel)) calc%parallel = parallel
    
    if (calc%verbose .and. calc%enabled) then
      print '(a)', "Adiabatic calculator created (enabled)"
      if (calc%parallel) then
        print '(a)', "  Parallel execution: enabled"
      end if
    end if
  end function create_calculator

  ! Alias for the main program
  function create_adiabatic_calculator(enabled, verbose, parallel) result(calc)
    logical, intent(in), optional :: enabled, verbose, parallel
    type(adiabatic_calculator) :: calc
    
    calc = create_calculator(enabled, verbose, parallel)
  end function create_adiabatic_calculator

  ! ============================================================================
  ! SET OUTPUT DIRECTORY
  ! ============================================================================
  subroutine set_output_dir(this, dirname)
    class(adiabatic_calculator), intent(inout) :: this
    character(len=*), intent(in) :: dirname
    
    this%output_dir = trim(dirname)
    if (this%verbose) then
      print '(a,a)', "Adiabatic output directory set to: ", trim(dirname)
    end if
  end subroutine set_output_dir

  ! ============================================================================
  ! INITIALIZE ADIABATIC CALCULATOR
  ! ============================================================================
  subroutine init_adiabatic(this, nstate, nnprod, tstep, filename)
    class(adiabatic_calculator), intent(inout) :: this
    integer(ip), intent(in) :: nstate, nnprod, tstep
    character(len=*), intent(in), optional :: filename
    
    integer(ip) :: i, j
    character(len=200) :: dir_command
    character(len=200) :: filepath
    real(dp) :: mem_est
    
    if (.not. this%enabled) return
    
    this%nstate = nstate
    this%nnprod = nnprod
    this%tstep = tstep
    this%current_step = 0
    this%transformation_computed = .false.
    
    if (this%verbose) then
      print '(a)', "Initializing adiabatic calculator..."
      print '(a,i0)', "  Number of states: ", nstate
      print '(a,i0)', "  Grid points: ", nnprod
      print '(a,i0)', "  Time steps: ", tstep
    end if
    
    ! Check memory requirements
    mem_est = real(nstate * nnprod * 8, dp) + &
              real(nstate * nstate * nnprod * 8, dp) + &
              real(nstate * tstep * 8, dp)
    mem_est = mem_est / (1024.0_dp * 1024.0_dp)  ! MB
    
    if (this%verbose) then
      print '(a,f10.2,a)', "  Estimated memory: ", mem_est, " MB"
    end if
    
    ! Create output directory
    dir_command = "mkdir -p " // trim(this%output_dir)
    call system(dir_command)
    
    ! Allocate storage
    allocate(this%energies(nstate, nnprod))
    allocate(this%transformation(nstate, nstate, nnprod))
    if (this%store_probabilities) then
      allocate(this%probabilities(nstate, tstep))
      allocate(this%times(tstep))
    end if
    
    ! Initialize
    this%energies = 0.0_dp
    this%transformation = 0.0_dp
    this%probabilities = 0.0_dp
    this%times = 0.0_dp
    
    ! Setup identity transformation
    do i = 1, nnprod
      do j = 1, nstate
        this%transformation(j, j, i) = 1.0_dp
      end do
    end do
    
    ! Open output files with explicit units (avoid conflict with main program
    ! 101-104)
    if (present(filename)) then
      filepath = trim(this%output_dir)//trim(filename)
    else
      filepath = trim(this%output_dir)//"adiabatic_probabilities.dat"
    end if
    
    open(unit=this%file_unit_prob, file=filepath, status='replace')
    write(this%file_unit_prob, '(a)') '# Time and adiabatic surface probabilities'
    write(this%file_unit_prob, '(a,*(i12))') '# Time', (i, i=1, nstate)
    flush(this%file_unit_prob)
    
    ! Open energies file
    filepath = trim(this%output_dir)//"adiabatic_energies.dat"
    open(unit=this%file_unit_energies, file=filepath, status='replace')
    write(this%file_unit_energies, '(a)') '# Adiabatic energies at each grid point'
    write(this%file_unit_energies, '(a,*(i12))') '# GridPoint', (i, i=1, nstate)
    flush(this%file_unit_energies)
    
    ! Open transformation file
    filepath = trim(this%output_dir)//"adiabatic_transformation.dat"
    open(unit=this%file_unit_trans, file=filepath, status='replace')
    write(this%file_unit_trans, '(a)') '# Adiabatic transformation matrix elements'
    write(this%file_unit_trans, '(a)') '# GridPoint, State_i, State_j, Transform_ij'
    flush(this%file_unit_trans)
    
    if (this%verbose) then
      print '(a)', "Adiabatic calculator initialized successfully."
      print '(a,a)', "  Output directory: ", trim(this%output_dir)
      print '(a,i0)', "  Probability file unit: ", this%file_unit_prob
    end if
  end subroutine init_adiabatic

  ! ============================================================================
  ! COMPUTE ADIABATIC TRANSFORMATION (Jacobi diagonalization)
  ! ============================================================================
  subroutine compute_transformation(this, potential_matrices, nnprod, nstate)
    class(adiabatic_calculator), intent(inout) :: this
    real(dp), intent(in) :: potential_matrices(:,:,:)  ! (nstate, nstate, nnprod)
    integer(ip), intent(in) :: nnprod, nstate
    
    integer(ip) :: i, j, k, iter, max_iter
    real(dp) :: tol, off_norm, c, s, t, tau
    real(dp) :: a_ii, a_jj, a_ij
    real(dp), allocatable :: A(:,:), V(:,:), R(:,:)
    integer :: p, q
    real(dp) :: start_time, end_time
    
    if (.not. this%enabled) return
    
    ! Check dimensions
    if (nstate /= this%nstate .or. nnprod /= this%nnprod) then
      if (this%verbose) then
        print '(a)', "ERROR: Dimension mismatch in compute_transformation"
        print '(a,i0,a,i0)', "  Expected nstate=", this%nstate, " got ", nstate
        print '(a,i0,a,i0)', "  Expected nnprod=", this%nnprod, " got ", nnprod
      end if
      return
    end if
    
    call cpu_time(start_time)
    
    if (this%verbose) then
      print '(a,i0,a)', "Computing adiabatic transformation for ", nnprod, " points"
      print '(a,i0)', "  Number of states: ", nstate
    end if
    
    max_iter = 1000
    tol = 1.0e-10_dp
    
    ! For nstate <= 3, use analytic diagonalization
    if (nstate <= 3) then
      call compute_adiabatic_small(this, potential_matrices, nnprod, nstate)
    else
      ! For larger matrices, use Jacobi
      call compute_adiabatic_jacobi(this, potential_matrices, nnprod, nstate)
    end if
    
    this%transformation_computed = .true.
    
    call cpu_time(end_time)
    this%compute_time = end_time - start_time
    
    ! Write energies to file
    call write_adiabatic_energies(this)
    
    if (this%verbose) then
      print '(a,f10.3,a)', "  Computation time: ", this%compute_time, " seconds"
      print '(a)', "Adiabatic transformation completed."
    end if
  end subroutine compute_transformation

  ! ============================================================================
  ! COMPUTE ADIABATIC USING JACOBI
  ! ============================================================================
  subroutine compute_adiabatic_jacobi(this, potential_matrices, nnprod, nstate)
    class(adiabatic_calculator), intent(inout) :: this
    real(dp), intent(in) :: potential_matrices(:,:,:)
    integer(ip), intent(in) :: nnprod, nstate
    
    integer(ip) :: i, j, k, iter, max_iter
    real(dp) :: tol, off_norm, c, s, t, tau
    real(dp) :: a_ii, a_jj, a_ij
    real(dp), allocatable :: A(:,:), V(:,:), R(:,:)
    integer :: p, q
    
    max_iter = 1000
    tol = 1.0e-10_dp
    
    do i = 1, nnprod
      allocate(A(nstate, nstate), V(nstate, nstate), R(nstate, nstate))
      
      ! Copy potential matrix
      A = potential_matrices(:,:,i)
      
      ! Initialize eigenvectors as identity
      V = 0.0_dp
      do j = 1, nstate
        V(j,j) = 1.0_dp
      end do
      
      ! Jacobi diagonalization
      do iter = 1, max_iter
        ! Find largest off-diagonal element
        p = 1
        q = 2
        off_norm = 0.0_dp
        do j = 1, nstate
          do k = j+1, nstate
            if (abs(A(j,k)) > off_norm) then
              off_norm = abs(A(j,k))
              p = j
              q = k
            end if
          end do
        end do
        
        ! Check convergence
        if (off_norm < tol) exit
        
        ! Compute rotation parameters
        a_ii = A(p,p)
        a_jj = A(q,q)
        a_ij = A(p,q)
        
        if (abs(a_ii - a_jj) < tol) then
          c = sqrt(0.5_dp)
          s = sign(sqrt(0.5_dp), a_ij)
        else
          tau = (a_jj - a_ii) / (2.0_dp * a_ij)
          t = sign(1.0_dp / (abs(tau) + sqrt(1.0_dp + tau**2)), tau)
          c = 1.0_dp / sqrt(1.0_dp + t**2)
          s = t * c
        end if
        
        ! Apply rotation to A
        R = A
        do k = 1, nstate
          if (k /= p .and. k /= q) then
            A(p,k) = c * R(p,k) - s * R(q,k)
            A(k,p) = A(p,k)
            A(q,k) = s * R(p,k) + c * R(q,k)
            A(k,q) = A(q,k)
          end if
        end do
        
        A(p,p) = c**2 * R(p,p) - 2.0_dp * c * s * R(p,q) + s**2 * R(q,q)
        A(q,q) = s**2 * R(p,p) + 2.0_dp * c * s * R(p,q) + c**2 * R(q,q)
        A(p,q) = 0.0_dp
        A(q,p) = 0.0_dp
        
        ! Update eigenvectors
        do k = 1, nstate
          t = V(k,p)
          V(k,p) = c * t - s * V(k,q)
          V(k,q) = s * t + c * V(k,q)
        end do
      end do
      
      ! Store results
      do j = 1, nstate
        this%energies(j,i) = A(j,j)
        this%transformation(:,j,i) = V(:,j)
      end do
      
      deallocate(A, V, R)
    end do
    
  end subroutine compute_adiabatic_jacobi

  ! ============================================================================
  ! COMPUTE ADIABATIC FOR SMALL MATRICES (2x2, 3x3)
  ! ============================================================================
  subroutine compute_adiabatic_small(this, potential_matrices, nnprod, nstate)
    class(adiabatic_calculator), intent(inout) :: this
    real(dp), intent(in) :: potential_matrices(:,:,:)
    integer(ip), intent(in) :: nnprod, nstate
    
    integer(ip) :: i, j, k
    real(dp) :: a, b, c, d, disc, trace
    real(dp) :: norm
    real(dp) :: A3(3,3)
    
    do i = 1, nnprod
      
      select case (nstate)
      case (2)
        ! 2x2 analytic diagonalization
        a = potential_matrices(1,1,i)
        b = potential_matrices(1,2,i)
        c = potential_matrices(2,1,i)
        d = potential_matrices(2,2,i)
        
        trace = a + d
        disc = sqrt(max((a - d)**2 + 4.0_dp * b * c, 0.0_dp))
        
        this%energies(1,i) = 0.5_dp * (trace - disc)
        this%energies(2,i) = 0.5_dp * (trace + disc)
        
        ! Eigenvectors
        if (abs(b) > SMALL .or. abs(c) > SMALL) then
          if (abs(b) >= abs(c)) then
            this%transformation(1,1,i) = (this%energies(1,i) - d) / b
            this%transformation(2,1,i) = 1.0_dp
            this%transformation(1,2,i) = (this%energies(2,i) - d) / b
            this%transformation(2,2,i) = 1.0_dp
          else
            this%transformation(1,1,i) = 1.0_dp
            this%transformation(2,1,i) = (this%energies(1,i) - a) / c
            this%transformation(1,2,i) = 1.0_dp
            this%transformation(2,2,i) = (this%energies(2,i) - a) / c
          end if
        else
          this%transformation(1,1,i) = 1.0_dp
          this%transformation(2,1,i) = 0.0_dp
          this%transformation(1,2,i) = 0.0_dp
          this%transformation(2,2,i) = 1.0_dp
        end if
        
        ! Normalize eigenvectors
        do j = 1, 2
          norm = sqrt(this%transformation(1,j,i)**2 + this%transformation(2,j,i)**2)
          if (norm > SMALL) then
            this%transformation(:,j,i) = this%transformation(:,j,i) / norm
          end if
        end do
        
      case (3)
        ! 3x3: Use Jacobi
        A3 = potential_matrices(:,:,i)
        call jacobi_3x3(A3, this%energies(:,i), this%transformation(:,:,i))
        
      case default
        continue
      end select
      
    end do
    
  end subroutine compute_adiabatic_small

  ! ============================================================================
  ! JACOBI FOR 3x3 MATRICES
  ! ============================================================================
  subroutine jacobi_3x3(A, eigenvalues, eigenvectors)
    real(dp), intent(in) :: A(3,3)
    real(dp), intent(out) :: eigenvalues(3)
    real(dp), intent(out) :: eigenvectors(3,3)
    
    real(dp) :: B(3,3), V(3,3)
    real(dp) :: c, s, t, tau, a_ii, a_jj, a_ij
    real(dp) :: off_norm
    integer :: iter, max_iter, p, q, k
    real(dp) :: tol = 1.0e-10_dp
    integer :: max_iter_local = 100
    integer :: j
    
    B = A
    V = 0.0_dp
    do k = 1, 3
      V(k,k) = 1.0_dp
    end do
    
    do iter = 1, max_iter_local
      p = 1
      q = 2
      off_norm = 0.0_dp
      do k = 1, 3
        do j = k+1, 3
          if (abs(B(k,j)) > off_norm) then
            off_norm = abs(B(k,j))
            p = k
            q = j
          end if
        end do
      end do
      
      if (off_norm < tol) exit
      
      a_ii = B(p,p)
      a_jj = B(q,q)
      a_ij = B(p,q)
      
      if (abs(a_ii - a_jj) < tol) then
        c = sqrt(0.5_dp)
        s = sign(sqrt(0.5_dp), a_ij)
      else
        tau = (a_jj - a_ii) / (2.0_dp * a_ij)
        t = sign(1.0_dp / (abs(tau) + sqrt(1.0_dp + tau**2)), tau)
        c = 1.0_dp / sqrt(1.0_dp + t**2)
        s = t * c
      end if
      
      do k = 1, 3
        if (k /= p .and. k /= q) then
          a_ii = B(p,k)
          a_jj = B(q,k)
          B(p,k) = c * a_ii - s * a_jj
          B(k,p) = B(p,k)
          B(q,k) = s * a_ii + c * a_jj
          B(k,q) = B(q,k)
        end if
      end do
      
      a_ii = B(p,p)
      a_jj = B(q,q)
      a_ij = B(p,q)
      B(p,p) = c**2 * a_ii - 2.0_dp * c * s * a_ij + s**2 * a_jj
      B(q,q) = s**2 * a_ii + 2.0_dp * c * s * a_ij + c**2 * a_jj
      B(p,q) = 0.0_dp
      B(q,p) = 0.0_dp
      
      do k = 1, 3
        a_ii = V(k,p)
        a_jj = V(k,q)
        V(k,p) = c * a_ii - s * a_jj
        V(k,q) = s * a_ii + c * a_jj
      end do
    end do
    
    eigenvalues = [(B(k,k), k=1, 3)]
    eigenvectors = V
    
  end subroutine jacobi_3x3

  ! ============================================================================
  ! TRANSFORM WAVEFUNCTION AND COMPUTE PROBABILITIES
  ! ============================================================================
  subroutine transform_wavefunction(this, diabatic_wf, nnprod, nstate, &
                                   time, step, store_probability)
    class(adiabatic_calculator), intent(inout) :: this
    complex(dp), intent(in) :: diabatic_wf(:,:)  ! (nstate, nnprod)
    integer(ip), intent(in) :: nnprod, nstate
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    logical, intent(in), optional :: store_probability
    
    complex(dp), allocatable :: adiabatic_wf(:,:)
    real(dp), allocatable :: probs(:)
    logical :: store
    integer(ip) :: i, j, k
    real(dp) :: total_prob
    
    if (.not. this%enabled) return
    
    ! Check dimensions
    if (nstate /= this%nstate .or. nnprod /= this%nnprod) then
      if (this%verbose) then
        print '(a)', "WARNING: Dimension mismatch in transform_wavefunction"
      end if
      return
    end if
    
    store = .true.
    if (present(store_probability)) store = store_probability
    
    allocate(adiabatic_wf(nstate, nnprod))
    allocate(probs(nstate))
    probs = 0.0_dp
    
    ! Transform wavefunction from diabatic to adiabatic basis
    do i = 1, nnprod
      do j = 1, nstate
        adiabatic_wf(j, i) = (0.0_dp, 0.0_dp)
        do k = 1, nstate
          ! Use transpose: U^T_{j,k} = U_{k,j}
          adiabatic_wf(j, i) = adiabatic_wf(j, i) + &
               this%transformation(k, j, i) * diabatic_wf(k, i)
        end do
      end do
      
      ! Accumulate probabilities
      do j = 1, nstate
        probs(j) = probs(j) + abs(adiabatic_wf(j, i))**2
      end do
    end do
    
    ! Store probabilities if requested
    if (store .and. step <= size(this%probabilities, 2)) then
      this%probabilities(:,step) = probs
      this%times(step) = time
      this%current_step = step
    end if
    
    ! Write to file - ALWAYS write to ensure data is saved
    ! Use explicit unit 300 (no conflict with main program)
    write(this%file_unit_prob, '(f16.10,*(f16.10))') time, probs
    flush(this%file_unit_prob)
    
    ! Debug output
    if (this%verbose .and. mod(step, 100) == 0) then
      total_prob = sum(probs)
      print '(a,f8.3,a)', "  Adiabatic data written at t=", time, " fs"
      print '(a,6f12.6)', "    Probs: ", probs
      print '(a,f12.6)', "    Total: ", total_prob
    end if
    
    deallocate(adiabatic_wf, probs)
    
  end subroutine transform_wavefunction

  ! ============================================================================
  ! COMPUTE ADIABATIC POPULATIONS
  ! ============================================================================
  subroutine compute_adiabatic_populations(this, diabatic_wf, nnprod, nstate, &
                                         time, step)
    class(adiabatic_calculator), intent(inout) :: this
    complex(dp), intent(in) :: diabatic_wf(:,:)
    integer(ip), intent(in) :: nnprod, nstate
    real(dp), intent(in) :: time
    integer(ip), intent(in) :: step
    
    call this%transform_wavefunction(diabatic_wf, nnprod, nstate, time, step, .true.)
    
  end subroutine compute_adiabatic_populations

  ! ============================================================================
  ! WRITE PROBABILITIES TO FILE
  ! ============================================================================
  subroutine write_probabilities(this, times, filename)
    class(adiabatic_calculator), intent(in) :: this
    real(dp), intent(in) :: times(:)
    character(len=*), intent(in), optional :: filename
    
    integer :: i, nsteps, nstates, unit_out
    character(len=200) :: filepath
    
    if (.not. this%enabled) return
    
    nsteps = min(size(times), size(this%probabilities, 2))
    nstates = size(this%probabilities, 1)
    
    if (nsteps <= 0) then
      if (this%verbose) print '(a)', "No probabilities to write"
      return
    end if
    
    if (present(filename)) then
      filepath = trim(this%output_dir)//trim(filename)
    else
      filepath = trim(this%output_dir)//"final_adiabatic_probs.dat"
    end if
    
    open(newunit=unit_out, file=filepath, status='replace')
    
    write(unit_out, '(a)') '# Final adiabatic probabilities vs time'
    write(unit_out, '(a,*(i12))') '# Time', (i, i=1, nstates)
    
    do i = 1, nsteps
      write(unit_out, '(f16.10,*(f16.10))') times(i), this%probabilities(:,i)
    end do
    
    close(unit_out)
    
    if (this%verbose) then
      print '(a)', "Adiabatic probabilities written to file."
    end if
  end subroutine write_probabilities

  ! ============================================================================
  ! WRITE ADIABATIC ENERGIES
  ! ============================================================================
  subroutine write_adiabatic_energies(this)
    class(adiabatic_calculator), intent(in) :: this
    
    integer(ip) :: i, j, k
    integer :: sample_size
    
    if (.not. this%enabled) return
    
    ! Write all energies
    do i = 1, this%nnprod
      write(this%file_unit_energies, '(i12,*(f16.8))') i, this%energies(:,i)
    end do
    flush(this%file_unit_energies)
    
    ! Write transformation matrix (sample points only)
    sample_size = min(100, this%nnprod)
    
    do i = 1, this%nnprod, max(1, this%nnprod / sample_size)
      do j = 1, this%nstate
        do k = 1, this%nstate
          if (abs(this%transformation(j,k,i)) > 1.0e-6_dp) then
            write(this%file_unit_trans, '(i12,2i12,f16.8)') i, j, k, this%transformation(j,k,i)
          end if
        end do
      end do
    end do
    flush(this%file_unit_trans)
    
  end subroutine write_adiabatic_energies

  ! ============================================================================
  ! GET ADIABATIC ENERGIES
  ! ============================================================================
  function get_adiabatic_energies(this) result(energies)
    class(adiabatic_calculator), intent(in) :: this
    real(dp), allocatable :: energies(:,:)
    
    if (allocated(this%energies)) then
      allocate(energies(size(this%energies,1), size(this%energies,2)))
      energies = this%energies
    end if
  end function get_adiabatic_energies

  ! ============================================================================
  ! GET TRANSFORMATION MATRIX
  ! ============================================================================
  function get_transformation(this) result(trans)
    class(adiabatic_calculator), intent(in) :: this
    real(dp), allocatable :: trans(:,:,:)
    
    if (allocated(this%transformation)) then
      allocate(trans(size(this%transformation,1), &
                     size(this%transformation,2), &
                     size(this%transformation,3)))
      trans = this%transformation
    end if
  end function get_transformation

  ! ============================================================================
  ! CLEANUP
  ! ============================================================================
  subroutine cleanup_adiabatic(this)
    class(adiabatic_calculator), intent(inout) :: this
    
    if (allocated(this%energies)) deallocate(this%energies)
    if (allocated(this%transformation)) deallocate(this%transformation)
    if (allocated(this%probabilities)) deallocate(this%probabilities)
    if (allocated(this%times)) deallocate(this%times)
    
    if (this%file_unit_prob > 0) close(this%file_unit_prob)
    if (this%file_unit_energies > 0) close(this%file_unit_energies)
    if (this%file_unit_trans > 0) close(this%file_unit_trans)
    
    this%file_unit_prob = 300
    this%file_unit_energies = 301
    this%file_unit_trans = 302
    
    if (this%verbose .and. this%enabled) then
      print '(a)', "Adiabatic calculator cleaned up."
      if (this%compute_time > 0.0_dp) then
        print '(a,f10.3,a)', "  Total computation time: ", this%compute_time, " seconds"
      end if
    end if
  end subroutine cleanup_adiabatic

  ! ============================================================================
  ! GETTERS AND SETTERS
  ! ============================================================================
  function is_enabled(this) result(enabled)
    class(adiabatic_calculator), intent(in) :: this
    logical :: enabled
    enabled = this%enabled
  end function is_enabled
  
  subroutine enable(this)
    class(adiabatic_calculator), intent(inout) :: this
    this%enabled = .true.
    if (this%verbose) print '(a)', "Adiabatic calculator enabled."
  end subroutine enable
  
  subroutine disable(this)
    class(adiabatic_calculator), intent(inout) :: this
    this%enabled = .false.
    if (this%verbose) print '(a)', "Adiabatic calculator disabled."
  end subroutine disable

 end module adiabatic_module
