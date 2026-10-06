! ============================================================
! POTENTIAL.F90 - BARRELENE WITH 6 STATES AND 42 MODES
! Based on: Journal of Molecular Structure 1110 (2016) 32e43
! Barrelene radical cation (BI⁺) vibronic Hamiltonian
! nstate = 6, nmode = 42 (full set including degeneracy)
! ============================================================
module potential_module
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  private
  
  ! Public subroutines
  public :: POTENTIAL, potential_
  
  ! Constants - make them module-level so all functions can access
  real(real64), parameter :: EV_TO_AU = 0.0367493d0   ! 1 eV = 0.0367493 Hartree
  real(real64), parameter :: PI = 3.14159265358979323846d0
  
  ! Barrelene has 42 vibrational modes in D₃h symmetry
  ! From Table 2 in paper:
  ! a₁': ν₁-ν₆ (6 modes)
  ! a₁'': ν₇-ν₈ (2 modes)
  ! a₂': ν₉ (1 mode)
  ! a₂'': ν₁₀-ν₁₄ (5 modes)
  ! e': ν₁₅-ν₂₁ (7 modes) - doubly degenerate (x,y)
  ! e'': ν₂₂-ν₂₈ (7 modes) - doubly degenerate (x,y)
  ! Total: 6 + 2 + 1 + 5 + 7×2 + 7×2 = 42 normal coordinates
  
  ! Electronic states (6 total: 4 distinct + degeneracy)
  ! 1: ~X²A₂' (ground)
  ! 2: ~A²E' component 1
  ! 3: ~A²E' component 2
  ! 4: ~B²E'' component 1
  ! 5: ~B²E'' component 2
  ! 6: ~C²A₁'

contains

! ============================================================
! SUBROUTINE: POTENTIAL
! ============================================================
subroutine POTENTIAL(grid, v_potential, nn, nmode, nstate, nnprod, freqn, ffn)

  implicit none
  
  ! Input parameters
  integer(int32), intent(in) :: nn, nmode, nstate, nnprod
  real(real64), intent(in) :: grid(nn, nmode)
  real(real64), intent(in) :: freqn(nmode), ffn(nmode)
  
  ! Output potential matrix
  real(real64), intent(out) :: v_potential(nstate, nstate, nnprod)
  
  ! Local variables
  integer :: i, j, k, m, point
  integer, dimension(:,:), allocatable :: tot
  real(real64), dimension(nstate, nstate) :: vf
  real(real64) :: Q, Qx, Qy
  
  ! Vertical Ionization Energies (VIE) from Table 3 (in eV)
  real(real64) :: VIE_X, VIE_A, VIE_B, VIE_C
  
  ! Initialize with paper values (eV)
  VIE_X = 7.93558d0    ! ~X²A₂'
  VIE_A = 9.46779d0    ! ~A²E'
  VIE_B = 11.40068d0   ! ~B²E''
  VIE_C = 11.66064d0   ! ~C²A₁'
  
  ! Convert to atomic units
  VIE_X = VIE_X * EV_TO_AU
  VIE_A = VIE_A * EV_TO_AU
  VIE_B = VIE_B * EV_TO_AU
  VIE_C = VIE_C * EV_TO_AU

  ! Allocate and compute index mapping
  allocate(tot(nnprod, nmode))
  call packindex_local(tot, nn, nmode, nnprod)

  ! Loop over all grid points
  do point = 1, nnprod

    ! Initialize potential matrix
    vf = 0.0d0
    
    ! ============================================================
    ! DIAGONAL ELEMENTS (Eq. 6-7 in paper)
    ! ============================================================
    
    ! Harmonic potential part (for all states)
    do m = 1, nmode
      Q = grid(tot(point,m), m)
      do i = 1, nstate
        vf(i,i) = vf(i,i) + 0.5d0 * freqn(m)**2 * Q**2
      end do
    end do
    
    ! State 1: ~X²A₂' (Eq. 6)
    vf(1,1) = vf(1,1) + VIE_X + linear_coupling_X(point, tot, grid, nmode) &
                         + quadratic_coupling_X(point, tot, grid, nmode)

    ! State 2: ~A²E' component 1 (Eq. 7, + combination)
    vf(2,2) = vf(2,2) + VIE_A + linear_coupling_A(point, tot, grid, nmode, '+') &
                         + quadratic_coupling_A(point, tot, grid, nmode)

    ! State 3: ~A²E' component 2 (Eq. 7, - combination)
    vf(3,3) = vf(3,3) + VIE_A + linear_coupling_A(point, tot, grid, nmode, '-') &
                         + quadratic_coupling_A(point, tot, grid, nmode)

    ! State 4: ~B²E'' component 1 (Eq. 7, + combination)
    vf(4,4) = vf(4,4) + VIE_B + linear_coupling_B(point, tot, grid, nmode, '+') &
                         + quadratic_coupling_B(point, tot, grid, nmode)

    ! State 5: ~B²E'' component 2 (Eq. 7, - combination)
    vf(5,5) = vf(5,5) + VIE_B + linear_coupling_B(point, tot, grid, nmode, '-') &
                         + quadratic_coupling_B(point, tot, grid, nmode)

    ! State 6: ~C²A₁' (Eq. 6)
    vf(6,6) = vf(6,6) + VIE_C + linear_coupling_C(point, tot, grid, nmode) &
                         + quadratic_coupling_C(point, tot, grid, nmode)

    ! ============================================================
    ! OFF-DIAGONAL ELEMENTS
    ! ============================================================
    
    ! JT coupling within degenerate states (Eq. 8)
    ! ~A state (2-3)
    vf(2,3) = JT_coupling_A(point, tot, grid, nmode)
    vf(3,2) = vf(2,3)
    
    ! ~B state (4-5)
    vf(4,5) = JT_coupling_B(point, tot, grid, nmode)
    vf(5,4) = vf(4,5)
    
    ! PJT couplings (Eq. 9-15 and Table 6)
    
    ! ~X - ~A coupling (Eq. 9-10)
    vf(1,2) = PJT_coupling_XA_x(point, tot, grid, nmode)
    vf(1,3) = PJT_coupling_XA_y(point, tot, grid, nmode)
    vf(2,1) = vf(1,2)
    vf(3,1) = vf(1,3)
    
    ! ~X - ~B coupling (Eq. 11-12)
    vf(1,4) = PJT_coupling_XB_x(point, tot, grid, nmode)
    vf(1,5) = PJT_coupling_XB_y(point, tot, grid, nmode)
    vf(4,1) = vf(1,4)
    vf(5,1) = vf(1,5)
    
    ! ~A - ~B coupling (Eq. 13-15)
    vf(2,4) = PJT_coupling_AB_x(point, tot, grid, nmode) &
              + PJT_coupling_AB_a1(point, tot, grid, nmode)
    vf(2,5) = PJT_coupling_AB_y(point, tot, grid, nmode)
    vf(3,4) = PJT_coupling_AB_y(point, tot, grid, nmode)
    vf(3,5) = -PJT_coupling_AB_x(point, tot, grid, nmode) &
              + PJT_coupling_AB_a1(point, tot, grid, nmode)
    
    ! Symmetric completion
    vf(4,2) = vf(2,4)
    vf(5,2) = vf(2,5)
    vf(4,3) = vf(3,4)
    vf(5,3) = vf(3,5)
    
    ! ~A - ~C coupling (from Table 6)
    vf(2,6) = PJT_coupling_AC(point, tot, grid, nmode)
    vf(3,6) = PJT_coupling_AC(point, tot, grid, nmode)
    vf(6,2) = vf(2,6)
    vf(6,3) = vf(3,6)
    
    ! ~B - ~C coupling (strong PJT from Table 6)
    vf(4,6) = PJT_coupling_BC(point, tot, grid, nmode)
    vf(5,6) = PJT_coupling_BC(point, tot, grid, nmode)
    vf(6,4) = vf(4,6)
    vf(6,5) = vf(5,6)

    ! Store in output array
    v_potential(:,:,point) = vf

  enddo
  
  deallocate(tot)

end subroutine POTENTIAL

! ============================================================
! COMPATIBILITY SUBROUTINE: potential_
! ============================================================
subroutine potential_(grid, v_potential, nn, nmode, nstate, nnprod, freqn, ffn)
  implicit none
  integer(int32), intent(in) :: nn, nmode, nstate, nnprod
  real(real64), intent(in) :: grid(nn, nmode)
  real(real64), intent(in) :: freqn(nmode), ffn(nmode)
  real(real64), intent(out) :: v_potential(nstate, nstate, nnprod)
  
  call POTENTIAL(grid, v_potential, nn, nmode, nstate, nnprod, freqn, ffn)
end subroutine potential_

! ============================================================
! HELPER SUBROUTINE: packindex_local
! ============================================================
subroutine packindex_local(tot, nn, nmode, nnprod)
  implicit none
  integer(int32), intent(in) :: nn, nmode, nnprod
  integer, dimension(nnprod, nmode), intent(out) :: tot
  
  integer :: i, j, prod
  integer, allocatable :: indices(:)
  
  allocate(indices(nmode))
  
  do i = 1, nnprod
    prod = i - 1
    do j = 1, nmode
      indices(j) = mod(prod, nn) + 1
      prod = prod / nn
    end do
    
    do j = 1, nmode
      tot(i,j) = indices(j)
    end do
  end do
  
  deallocate(indices)
end subroutine packindex_local

! ============================================================
! LINEAR COUPLING FOR ~X STATE (Table 3)
! ============================================================
function linear_coupling_X(point, tot, grid, nmode) result(v_lin)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_lin
  
  ! k_i parameters from Table 3 for ~X state (eV)
  ! Mode:    ν₁     ν₂     ν₃     ν₄     ν₅     ν₆
  real(real64) :: kappa_X(6)
  real(real64) :: Q
  integer :: m
  
  kappa_X = (/ 0.0096d0, -0.0167d0, 0.1461d0, 0.0534d0, 0.0173d0, -0.1044d0 /)
  
  v_lin = 0.0d0
  do m = 1, min(6, nmode)
    Q = grid(tot(point,m), m)
    v_lin = v_lin + kappa_X(m) * EV_TO_AU * Q
  end do
end function linear_coupling_X

! ============================================================
! LINEAR COUPLING FOR ~A STATE (Table 3, Eq. 7)
! ============================================================
function linear_coupling_A(point, tot, grid, nmode, sign_char) result(v_lin)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  character(len=1), intent(in) :: sign_char
  real(real64) :: v_lin
  
  ! k_i parameters from Table 3 for ~A state (eV)
  real(real64) :: kappa_A(6)
  ! l_i parameters for JT coupling from Table 4 (eV)
  real(real64) :: lambda_A(7)  ! for e' modes ν₁₅-ν₂₁
  real(real64) :: Q, sign_factor
  integer :: m
  
  kappa_A = (/ 0.0132d0, -0.0178d0, 0.1115d0, 0.0816d0, -0.0401d0, -0.0357d0 /)
  lambda_A = (/ 0.0019d0, 0.1175d0, 0.0712d0, 0.0127d0, 0.0850d0, 0.0192d0, 0.1389d0 /)
  
  if (sign_char == '+') then
    sign_factor = 1.0d0
  else
    sign_factor = -1.0d0
  end if
  
  v_lin = 0.0d0
  
  ! a₁' modes (ν₁-ν₆)
  do m = 1, min(6, nmode)
    Q = grid(tot(point,m), m)
    v_lin = v_lin + kappa_A(m) * EV_TO_AU * Q
  end do
  
  ! JT coupling for e' modes (ν₁₅-ν₂₁)
  do m = 1, min(7, nmode-14)
    if (14+m <= nmode) then
      Q = grid(tot(point,14+m), 14+m)
      v_lin = v_lin + sign_factor * lambda_A(m) * EV_TO_AU * Q
    end if
  end do
end function linear_coupling_A

! ============================================================
! LINEAR COUPLING FOR ~B STATE (Table 3, Eq. 7)
! ============================================================
function linear_coupling_B(point, tot, grid, nmode, sign_char) result(v_lin)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  character(len=1), intent(in) :: sign_char
  real(real64) :: v_lin
  
  ! k_i parameters from Table 3 for ~B state (eV)
  real(real64) :: kappa_B(6)
  ! l_i parameters for JT coupling from Table 4 (eV)
  real(real64) :: lambda_B(7)  ! for e' modes ν₁₅-ν₂₁
  real(real64) :: Q, sign_factor
  integer :: m
  
  kappa_B = (/-0.0846d0, -0.0396d0, -0.2295d0, 0.1266d0, -0.0881d0, 0.0287d0 /)
  lambda_B = (/ 0.0754d0, 0.0993d0, 0.0602d0, 0.0779d0, 0.0901d0, 0.0943d0, 0.0330d0 /)
  
  if (sign_char == '+') then
    sign_factor = 1.0d0
  else
    sign_factor = -1.0d0
  end if
  
  v_lin = 0.0d0
  
  ! a₁' modes (ν₁-ν₆)
  do m = 1, min(6, nmode)
    Q = grid(tot(point,m), m)
    v_lin = v_lin + kappa_B(m) * EV_TO_AU * Q
  end do
  
  ! JT coupling for e' modes (ν₁₅-ν₂₁)
  do m = 1, min(7, nmode-14)
    if (14+m <= nmode) then
      Q = grid(tot(point,14+m), 14+m)
      v_lin = v_lin + sign_factor * lambda_B(m) * EV_TO_AU * Q
    end if
  end do
end function linear_coupling_B

! ============================================================
! LINEAR COUPLING FOR ~C STATE (Table 3)
! ============================================================
function linear_coupling_C(point, tot, grid, nmode) result(v_lin)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_lin
  
  ! k_i parameters from Table 3 for ~C state (eV)
  real(real64) :: kappa_C(6)
  real(real64) :: Q
  integer :: m
  
  kappa_C = (/-0.0334d0, 0.1907d0, 0.2425d0, -0.1407d0, 0.0417d0, 0.0662d0 /)
  
  v_lin = 0.0d0
  do m = 1, min(6, nmode)
    Q = grid(tot(point,m), m)
    v_lin = v_lin + kappa_C(m) * EV_TO_AU * Q
  end do
end function linear_coupling_C

! ============================================================
! QUADRATIC COUPLING FUNCTIONS (From Table 3)
! ============================================================
function quadratic_coupling_X(point, tot, grid, nmode) result(v_quad)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_quad
  
  ! g_i parameters from Table 3 for ~X state (eV)
  real(real64) :: gamma_X(6)
  real(real64) :: Q
  integer :: m
  
  gamma_X = (/ 0.0141d0, -0.0061d0, 0.0362d0, 0.0140d0, 0.0030d0, -0.0099d0 /)
  
  v_quad = 0.0d0
  do m = 1, min(6, nmode)
    Q = grid(tot(point,m), m)
    v_quad = v_quad + 0.5d0 * gamma_X(m) * EV_TO_AU * Q**2
  end do
end function quadratic_coupling_X

function quadratic_coupling_A(point, tot, grid, nmode) result(v_quad)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_quad
  
  ! g_i parameters from Table 3 for ~A state (eV)
  real(real64) :: gamma_A(6)
  real(real64) :: Q
  integer :: m
  
  gamma_A = (/ 0.0115d0, -0.0061d0, -0.0138d0, 0.0193d0, -0.0312d0, -0.0529d0 /)
  
  v_quad = 0.0d0
  do m = 1, min(6, nmode)
    Q = grid(tot(point,m), m)
    v_quad = v_quad + 0.5d0 * gamma_A(m) * EV_TO_AU * Q**2
  end do
end function quadratic_coupling_A

function quadratic_coupling_B(point, tot, grid, nmode) result(v_quad)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_quad
  
  ! g_i parameters from Table 3 for ~B state (eV)
  real(real64) :: gamma_B(6)
  real(real64) :: Q
  integer :: m
  
  gamma_B = (/ -0.0282d0, -0.0161d0, 0.0125d0, 0.0047d0, 0.0058d0, -0.0113d0 /)
  
  v_quad = 0.0d0
  do m = 1, min(6, nmode)
    Q = grid(tot(point,m), m)
    v_quad = v_quad + 0.5d0 * gamma_B(m) * EV_TO_AU * Q**2
  end do
end function quadratic_coupling_B

function quadratic_coupling_C(point, tot, grid, nmode) result(v_quad)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_quad
  
  ! g_i parameters from Table 3 for ~C state (eV)
  real(real64) :: gamma_C(6)
  real(real64) :: Q
  integer :: m
  
  gamma_C = (/ 0.0272d0, 0.0148d0, -0.0630d0, 0.0475d0, 0.0055d0, -0.0328d0 /)
  
  v_quad = 0.0d0
  do m = 1, min(6, nmode)
    Q = grid(tot(point,m), m)
    v_quad = v_quad + 0.5d0 * gamma_C(m) * EV_TO_AU * Q**2
  end do
end function quadratic_coupling_C

! ============================================================
! JT COUPLING FUNCTIONS (Eq. 8)
! ============================================================
function JT_coupling_A(point, tot, grid, nmode) result(v_JT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_JT
  
  ! l_i parameters for ~A state from Table 4 (eV)
  real(real64) :: lambda_A(7)
  real(real64) :: Qx, Qy  ! x and y components for degenerate modes
  integer :: m
  
  lambda_A = (/ 0.0019d0, 0.1175d0, 0.0712d0, 0.0127d0, 0.0850d0, 0.0192d0, 0.1389d0 /)
  
  v_JT = 0.0d0
  ! Eq. 8: Σ l_i^A Q_iy (for e' modes, using Q_iy)
  do m = 1, min(7, nmode-14)
    if (14+m <= nmode) then
      Qy = grid(tot(point,14+m), 14+m)
      v_JT = v_JT + lambda_A(m) * EV_TO_AU * Qy
    end if
  end do
end function JT_coupling_A

function JT_coupling_B(point, tot, grid, nmode) result(v_JT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_JT
  
  ! l_i parameters for ~B state from Table 4 (eV)
  real(real64) :: lambda_B(7)
  real(real64) :: Qx, Qy
  integer :: m
  
  lambda_B = (/ 0.0754d0, 0.0993d0, 0.0602d0, 0.0779d0, 0.0901d0, 0.0943d0, 0.0330d0 /)
  
  v_JT = 0.0d0
  do m = 1, min(7, nmode-14)
    if (14+m <= nmode) then
      Qy = grid(tot(point,14+m), 14+m)
      v_JT = v_JT + lambda_B(m) * EV_TO_AU * Qy
    end if
  end do
end function JT_coupling_B

! ============================================================
! PJT COUPLING FUNCTIONS (Table 6)
! ============================================================
function PJT_coupling_XA_x(point, tot, grid, nmode) result(v_PJT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_PJT
  
  ! ~X - ~A coupling via ν₁₆ and ν₂₁ (Table 6)
  real(real64) :: Q16, Q21
  
  v_PJT = 0.0d0
  
  ! ν₁₆ coupling: 0.1529 eV
  if (nmode >= 16) then
    Q16 = grid(tot(point,16), 16)
    v_PJT = v_PJT + 0.1529d0 * EV_TO_AU * Q16
  end if
  
  ! ν₂₁ coupling: 0.1098 eV
  if (nmode >= 21) then
    Q21 = grid(tot(point,21), 21)
    v_PJT = v_PJT + 0.1098d0 * EV_TO_AU * Q21
  end if
end function PJT_coupling_XA_x

function PJT_coupling_XA_y(point, tot, grid, nmode) result(v_PJT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_PJT
  
  ! Eq. 10: -λ^{X-A} Qy
  ! Using same modes as x-component but with negative sign
  real(real64) :: Q16, Q21
  
  v_PJT = 0.0d0
  
  ! ν₁₆ coupling: -0.1529 eV
  if (nmode >= 16) then
    Q16 = grid(tot(point,16), 16)
    v_PJT = v_PJT - 0.1529d0 * EV_TO_AU * Q16
  end if
  
  ! ν₂₁ coupling: -0.1098 eV
  if (nmode >= 21) then
    Q21 = grid(tot(point,21), 21)
    v_PJT = v_PJT - 0.1098d0 * EV_TO_AU * Q21
  end if
end function PJT_coupling_XA_y

function PJT_coupling_XB_x(point, tot, grid, nmode) result(v_PJT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_PJT
  
  ! ~X - ~B coupling via ν₂₂ and ν₂₈ (Table 6)
  real(real64) :: Q22, Q28
  
  v_PJT = 0.0d0
  
  ! ν₂₂ coupling: -0.0867 eV
  if (nmode >= 22) then
    Q22 = grid(tot(point,22), 22)
    v_PJT = v_PJT - 0.0867d0 * EV_TO_AU * Q22
  end if
  
  ! ν₂₈ coupling: -0.1496 eV
  if (nmode >= 28) then
    Q28 = grid(tot(point,28), 28)
    v_PJT = v_PJT - 0.1496d0 * EV_TO_AU * Q28
  end if
end function PJT_coupling_XB_x

function PJT_coupling_XB_y(point, tot, grid, nmode) result(v_PJT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_PJT
  
  ! Eq. 12: -λ^{X-B} Qy
  real(real64) :: Q22, Q28
  
  v_PJT = 0.0d0
  
  ! ν₂₂ coupling: +0.0867 eV (opposite sign)
  if (nmode >= 22) then
    Q22 = grid(tot(point,22), 22)
    v_PJT = v_PJT + 0.0867d0 * EV_TO_AU * Q22
  end if
  
  ! ν₂₈ coupling: +0.1496 eV (opposite sign)
  if (nmode >= 28) then
    Q28 = grid(tot(point,28), 28)
    v_PJT = v_PJT + 0.1496d0 * EV_TO_AU * Q28
  end if
end function PJT_coupling_XB_y

function PJT_coupling_AB_x(point, tot, grid, nmode) result(v_PJT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_PJT
  
  ! ~A - ~B coupling via e'' modes (Table 6)
  real(real64) :: Q25, Q26, Q27
  
  v_PJT = 0.0d0
  
  ! Using ν₂₅: 0.0606 eV
  if (nmode >= 25) then
    Q25 = grid(tot(point,25), 25)
    v_PJT = v_PJT + 0.0606d0 * EV_TO_AU * Q25
  end if
  
  ! Additional coupling via ν₂₆: -0.1341 eV
  if (nmode >= 26) then
    Q26 = grid(tot(point,26), 26)
    v_PJT = v_PJT - 0.1341d0 * EV_TO_AU * Q26
  end if
  
  ! Additional coupling via ν₂₇: -0.0886 eV
  if (nmode >= 27) then
    Q27 = grid(tot(point,27), 27)
    v_PJT = v_PJT - 0.0886d0 * EV_TO_AU * Q27
  end if
end function PJT_coupling_AB_x

function PJT_coupling_AB_y(point, tot, grid, nmode) result(v_PJT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_PJT
  
  ! ~A - ~B coupling y-component (Table 6)
  real(real64) :: Q25, Q26, Q27
  
  v_PJT = 0.0d0
  
  ! Using ν₂₅: 0.0606 eV
  if (nmode >= 25) then
    Q25 = grid(tot(point,25), 25)
    v_PJT = v_PJT + 0.0606d0 * EV_TO_AU * Q25
  end if
  
  ! ν₂₆: 0.1341 eV (opposite sign from x)
  if (nmode >= 26) then
    Q26 = grid(tot(point,26), 26)
    v_PJT = v_PJT + 0.1341d0 * EV_TO_AU * Q26
  end if
  
  ! ν₂₇: 0.0886 eV (opposite sign from x)
  if (nmode >= 27) then
    Q27 = grid(tot(point,27), 27)
    v_PJT = v_PJT + 0.0886d0 * EV_TO_AU * Q27
  end if
end function PJT_coupling_AB_y

function PJT_coupling_AB_a1(point, tot, grid, nmode) result(v_PJT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_PJT
  
  ! ~A - ~B coupling via a₁'' modes (ν₇ and ν₈)
  real(real64) :: Q7, Q8
  
  v_PJT = 0.0d0
  
  ! ν₇: 0.1106 eV
  if (nmode >= 7) then
    Q7 = grid(tot(point,7), 7)
    v_PJT = v_PJT + 0.1106d0 * EV_TO_AU * Q7
  end if
  
  ! ν₈: 0.0475 eV
  if (nmode >= 8) then
    Q8 = grid(tot(point,8), 8)
    v_PJT = v_PJT + 0.0475d0 * EV_TO_AU * Q8
  end if
end function PJT_coupling_AB_a1

function PJT_coupling_AC(point, tot, grid, nmode) result(v_PJT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_PJT
  
  ! ~A - ~C coupling (from discussion in paper)
  ! Using ν₂₀ as mentioned: 0.0198 eV
  real(real64) :: Q20, Q24
  
  v_PJT = 0.0d0
  
  ! ν₂₀: 0.0198 eV
  if (nmode >= 20) then
    Q20 = grid(tot(point,20), 20)
    v_PJT = v_PJT + 0.0198d0 * EV_TO_AU * Q20
  end if
  
  ! ν₂₄: -0.1421 eV
  if (nmode >= 24) then
    Q24 = grid(tot(point,24), 24)
    v_PJT = v_PJT - 0.1421d0 * EV_TO_AU * Q24
  end if
end function PJT_coupling_AC

function PJT_coupling_BC(point, tot, grid, nmode) result(v_PJT)
  implicit none
  integer, intent(in) :: point, nmode
  integer, dimension(:,:), intent(in) :: tot
  real(real64), dimension(:,:), intent(in) :: grid
  real(real64) :: v_PJT
  
  ! ~B - ~C coupling (strong PJT from Table 6)
  real(real64) :: Q23, Q24
  
  v_PJT = 0.0d0
  
  ! ν₂₃: 0.2293 eV
  if (nmode >= 23) then
    Q23 = grid(tot(point,23), 23)
    v_PJT = v_PJT + 0.2293d0 * EV_TO_AU * Q23
  end if
  
  ! ν₂₄: 0.1421 eV
  if (nmode >= 24) then
    Q24 = grid(tot(point,24), 24)
    v_PJT = v_PJT + 0.1421d0 * EV_TO_AU * Q24
  end if
end function PJT_coupling_BC

end module potential_module