module mod_variables
  use mod_basic_types

  implicit none
  public

  !> Number of flux variables
  integer           :: nwflux = 0

  !> Number of auxiliary variables
  integer           :: nwaux = 0

  !> Number of extra variables
  integer           :: nwextra = 0

  !> Total number of variables
  integer           :: nw = 0

  !> Number of vector variables (used for writing output)
  integer           :: nvector = 0

  !> Indices of vector variables
  integer, dimension(:), allocatable :: iw_vector

  !> Maximum number of variables
  integer, parameter :: max_nw = 50

  !> Primitive variable names
  character(len=name_len) :: prim_wnames(max_nw)

  !> Conservative variable names
  character(len=name_len) :: cons_wnames(max_nw)

  ! Global indices of variables that are often used

  !> Index of the (gas) density
  integer, protected :: iw_rho = -1

  !> Indices of the momentum density
  integer, allocatable, protected :: iw_mom(:)

  !> Index of the energy density
  integer, protected :: iw_e = -1

  !> Indices of the magnetic field components
  integer, allocatable, protected :: iw_mag(:)

contains

  !> Set generic flux variable
  function var_set_fluxvar(name_cons, name_prim, ix) result(iw)
    character(len=*), intent(in)  :: name_cons, name_prim
    integer, intent(in), optional :: ix
    integer                       :: iw

    nwflux = nwflux + 1
    nw     = nw + 1
    iw     = nwflux

    if (.not. present(ix)) then
      prim_wnames(nwflux) = name_cons
      cons_wnames(nwflux) = name_prim
    else
      write(cons_wnames(nwflux),"(A,I0)") name_cons, ix
      write(prim_wnames(nwflux),"(A,I0)") name_prim, ix
    end if
  end function var_set_fluxvar

  !> Set density variable
  function var_set_rho() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nw                  = nw + 1
    iw_rho              = nwflux
    iw                  = nwflux
    prim_wnames(nwflux) = 'rho'
    cons_wnames(nwflux) = 'rho'
  end function var_set_rho

  !> Set momentum variables
  function var_set_momentum(ndir) result(iw)
    integer, intent(in) :: ndir
    integer             :: iw(ndir), idir

    if (allocated(iw_mom)) &
         call mpistop("Error: set_mom was already called")
    allocate(iw_mom(ndir))

    do idir = 1, ndir
      nwflux       = nwflux + 1
      nw           = nw + 1
      iw_mom(idir) = nwflux
      iw(idir)     = nwflux
      write(cons_wnames(nwflux),"(A1,I1)") "m", idir
      write(prim_wnames(nwflux),"(A1,I1)") "v", idir
    end do
  end function var_set_momentum

  !> Set energy variable
  function var_set_energy() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nw                  = nw + 1
    iw_e                = nwflux
    iw                  = nwflux
    cons_wnames(nwflux) = 'e'
    prim_wnames(nwflux) = 'p'
  end function var_set_energy

  !> Set magnetic field variables
  function var_set_bfield(ndir) result(iw)
    integer, intent(in) :: ndir
    integer             :: iw(ndir), idir

    if (allocated(iw_mag)) &
         call mpistop("Error: set_mag was already called")
    allocate(iw_mag(ndir))

    do idir = 1, ndir
      nwflux       = nwflux + 1
      nw           = nw + 1
      iw_mag(idir) = nwflux
      iw(idir)     = nwflux
      write(cons_wnames(nwflux),"(A1,I1)") "b", idir
      write(prim_wnames(nwflux),"(A1,I1)") "b", idir
    end do
  end function var_set_bfield

end module mod_variables