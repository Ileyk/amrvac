module mod_rho_phys

  implicit none
  private

  integer, protected, public          :: rho_       = 1
  double precision, protected, public :: rho_v(^ND) = 1.0d0

  ! Public methods
  public :: rho_phys_init
  public :: rho_get_v

contains

  subroutine rho_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rho_list/ rho_v

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, rho_list, end=111)
111    close(unitpar)
    end do

  end subroutine rho_params_read

  !> Write this module's parameters to a snapsoht
  subroutine rho_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = ^ND
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er
    integer                             :: idim

    call MPI_FILE_WRITE(fh, n_par, ^ND, MPI_INTEGER, st, er)

    do idim=1,ndim
      write(names(idim),'(a,i1)') "v",idim
      values(idim) = rho_v(idim)
    end do
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine rho_write_info

  subroutine rho_phys_init()
    use mod_global_parameters
    use mod_physics

    call rho_params_read(par_files)

    physics_type = "rho"
    phys_energy  = .false.

    rho_ = var_set_rho()

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    phys_get_cmax        => rho_get_cmax
    phys_get_flux        => rho_get_flux
    phys_add_source_geom => rho_add_source_geom
    phys_to_conserved    => rho_to_conserved
    phys_to_primitive    => rho_to_primitive
    phys_get_dt          => rho_get_dt
    phys_write_info      => rho_write_info

  end subroutine rho_phys_init

  subroutine rho_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for rho module)
  end subroutine rho_to_conserved

  subroutine rho_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for rho module)
  end subroutine rho_to_primitive

  subroutine rho_get_v(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(out) :: v(ixI^S)

    v(ixO^S) = rho_v(idim)
  end subroutine rho_get_v

  subroutine rho_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    call rho_get_v(w, x, ixI^L, ixO^L, idim, cmax)

    if (present(cmin)) then
       cmin(ixO^S) = min(cmax(ixO^S), zero)
       cmax(ixO^S) = max(cmax(ixO^S), zero)
    else
       cmax(ixO^S) = abs(cmax(ixO^S))
    end if
  end subroutine rho_get_cmax

  subroutine rho_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble
  end subroutine rho_get_dt

  ! There is nothing to add to the transport flux in the transport equation
  subroutine rho_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    double precision                :: v(ixI^S)

    call rho_get_v(w, x, ixI^L, ixO^L, idim, v)

    f(ixO^S, rho_) = w(ixO^S, rho_) * v(ixO^S)
  end subroutine rho_get_flux

  subroutine rho_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)

    ! Add geometrical source terms to w
    ! There are no geometrical source terms in the transport equation

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

  end subroutine rho_add_source_geom

end module mod_rho_phys