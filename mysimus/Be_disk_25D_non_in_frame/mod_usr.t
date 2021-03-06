module mod_usr
  use mod_hd
  implicit none
  double precision :: q, f
  double precision :: a, Egg, Om0

contains

  subroutine usr_init()
    call set_coordinate_system('spherical_2.5D')

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_source        => pt_grav_source
    usr_get_dt        => get_dt_pt_grav
    usr_special_bc    => specialbound_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    call params_read(par_files)

    call hd_activate()

  end subroutine usr_init

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ q, f

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read

  subroutine initglobaldata_usr

    Egg = ( 0.49*q**(2./3.) ) / ( 0.6*q**(2./3.) + dlog(1.+q**(1./3.)) )
    a  = 1.d0/(f*Egg) ! in units of Rstar

    ! ang. speed of frame in units of sqrt(GM/R3)
    Om0 = dsqrt((f*Egg)**3.d0*(1.d0+q)/q)
    print*, 'Orb. sep. :', a
    print*, 'Om0 :', Om0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ix^S,rho_)   =  (x(ix^S,1)*dsin(x(ix^S,2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ix^S,2))-1.d0)/x(ix^S,1))
    w(ix^S,mom(1)) =  0.d0
    w(ix^S,mom(2)) =  0.d0
    where (w(ix^S,rho_)>small_density)
      w(ix^S,mom(3)) =  w(ix^S,rho_)*dsqrt(1.d0/(x(ix^S,1)*dsin(x(ix^S,2)))-3.5d0*hd_adiab)-w(ix^S,rho_)*x(ix^S,1)*dsin(x(ix^S,2))*Om0
    elsewhere
      w(ix^S,rho_)   =  small_density
      w(ix^S,mom(3)) =  -w(ix^S,rho_)*x(ix^S,1)*dsin(x(ix^S,2))*Om0
    end where

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: cent_str(ixI^S,1:ndir), cor(ixI^S,1:ndir)

    call get_cntfgl_str(ixI^L,ixO^L,x,cent_str)
    call get_coriolis(ixI^L,ixO^L,x,wCT,cor)
    where (wCT(ixO^S,rho_)>10.d0*small_density)
      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) - qdt * wCT(ixO^S,rho_)    / x(ixO^S,1)**2.d0
      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * wCT(ixO^S,rho_  ) * (cent_str(ixO^S,1)+cor(ixO^S,1))
      w(ixO^S,mom(2)) = w(ixO^S,mom(2)) + qdt * wCT(ixO^S,rho_  ) * (cent_str(ixO^S,2)+cor(ixO^S,2))
      w(ixO^S,mom(3)) = w(ixO^S,mom(3)) + qdt * wCT(ixO^S,rho_  ) * (cent_str(ixO^S,3)+cor(ixO^S,3))
    endwhere

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
    double precision, intent(in) :: w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew
    double precision :: force(ixG^S,1:ndir), cent_str(ixG^S,1:ndir), cor(ixG^S,1:ndir)

    dtnew=bigdouble
    force=0.d0
    force(ix^S,1)=1.d0/x(ix^S,1)**2.d0
    call get_cntfgl_str(ixG^L,ix^L,x,cent_str)
    call get_coriolis(ixG^L,ix^L,x,w,cor)
    force(ix^S,1:ndim)=force(ix^S,1:ndim)+cor(ix^S,1:ndim)+cent_str(ix^S,1:ndim)
    ! dtnew=min(dtnew,minval(dsqrt(block%dx(ix^S,1)/(1.d0/x(ix^S,1)**2.d0))))
    dtnew=min(dtnew,minval(dsqrt(block%dx(ix^S,1)/&
      dabs(force(ix^S,1)))),&
                    minval(dsqrt(block%dx(ix^S,1)*block%dx(ix^S,2)/&
      dabs(force(ix^S,2)))))

  end subroutine get_dt_pt_grav

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: rho(ixG^S)

    select case(iB)
    case(1)

    w(ixB^S,rho_)   =  (x(ixB^S,1)*dsin(x(ixB^S,2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ixB^S,2))-1.d0)/x(ixB^S,1))
    w(ixB^S,mom(1)) =  0.d0
    w(ixB^S,mom(2)) =  0.d0
    where (w(ixB^S,rho_)>small_density)
      w(ixB^S,mom(3)) =  w(ixB^S,rho_)*dsqrt(1.d0/(x(ixB^S,1)*dsin(x(ixB^S,2)))-3.5d0*hd_adiab)-w(ixB^S,rho_)*x(ixB^S,1)*dsin(x(ixB^S,2))*Om0
    elsewhere
      w(ixB^S,rho_)   = small_density
      w(ixB^S,mom(3)) = -w(ixB^S,rho_)*x(ixB^S,1)*dsin(x(ixB^S,2))*Om0
    end where

    case(2)

    w(ixB^S,rho_)   =  (x(ixB^S,1)*dsin(x(ixB^S,2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ixB^S,2))-1.d0)/x(ixB^S,1))
    w(ixB^S,mom(1)) =  0.d0
    w(ixB^S,mom(2)) =  0.d0
    where (w(ixB^S,rho_)>small_density)
      w(ixB^S,mom(3)) =  w(ixB^S,rho_)*dsqrt(1.d0/(x(ixB^S,1)*dsin(x(ixB^S,2)))-3.5d0*hd_adiab)-w(ixB^S,rho_)*x(ixB^S,1)*dsin(x(ixB^S,2))*Om0
    elsewhere
      w(ixB^S,rho_)   = small_density
      w(ixB^S,mom(3)) = -w(ixB^S,rho_)*x(ixB^S,1)*dsin(x(ixB^S,2))*Om0
    end where

    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  ! -------------------------------------------------------------------------------
  ! Center @ center of masses
  ! x axis along line joining 2 bodies, from donor to accretor
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_cntfgl_str(ixI^L,ixO^L,x,F_Roche)
  use mod_global_parameters
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in)  :: x(ixI^S,1:ndim)
  double precision, intent(out) :: F_Roche(ixI^S,1:ndir)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! centrifugal centered on star rotating @ Omega orbital so as the star does not spin
  ! in the inertial frame (but it does wobble around the center of mass, hence the need
  ! for the centrifugal term in the next subroutine)
  F_Roche(ixO^S,1)= Om0**2.d0 * x(ixO^S,1)*(dsin(x(ixO^S,2)))**2.d0
  F_Roche(ixO^S,2)= Om0**2.d0 * x(ixO^S,1)*dsin(x(ixO^S,2))*dcos(x(ixO^S,2))
  F_Roche(ixO^S,3)= 0.d0
  end subroutine get_cntfgl_str
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Rotation axis along +z, w/ y axis pointing upward in slices of the orbital plane
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_coriolis(ixI^L,ixO^L,x,wCT,F_Coriolis)
  use mod_global_parameters
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
  double precision, intent(out) :: F_Coriolis(ixI^S,1:ndir)
  integer :: i
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  F_Coriolis(ixO^S,1) = - 2. * Om0 * (-wCT(ixO^S,mom(3))*dsin(x(ixO^S,2)))
  F_Coriolis(ixO^S,2) = - 2. * Om0 * (-wCT(ixO^S,mom(3))*dcos(x(ixO^S,2)))
  F_Coriolis(ixO^S,3) = - 2. * Om0 * ( wCT(ixO^S,mom(2))*dcos(x(ixO^S,2))+wCT(ixO^S,mom(1))*dsin(x(ixO^S,2)) )
  do i=1,ndir
    F_Coriolis(ixO^S,i) = F_Coriolis(ixO^S,i) / wCT(ixO^S,rho_  )
  enddo
  end subroutine get_coriolis
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  integer, intent(in)                :: ixI^L,ixO^L
  double precision, intent(in)       :: x(ixI^S,1:ndim)
  double precision                   :: w(ixI^S,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision :: cent_str(ixI^S,1:ndir), cor(ixI^S,1:ndir), pthermal(ixI^S), gradP(ixI^S)

   w(ixO^S,nw+1) = w(ixO^S,mom(1))/w(ixO^S,rho_)
   w(ixO^S,nw+2) = w(ixO^S,mom(2))/w(ixO^S,rho_)
   w(ixO^S,nw+3) = w(ixO^S,mom(3))/w(ixO^S,rho_)
   w(ixO^S,nw+4) = w(ixO^S,mom(3))/w(ixO^S,rho_)+x(ixO^S,1)*dsin(x(ixO^S,2))*Om0
   w(ixO^S,nw+5) = dsqrt(w(ixO^S,nw+1)**2.d0+w(ixO^S,nw+2)**2.d0+w(ixO^S,nw+3)**2.d0)
   w(ixO^S,nw+6) = - 1.d0 / x(ixO^S,1)**2.d0
   call get_cntfgl_str(ixI^L,ixO^L,x,cent_str)
   w(ixO^S,nw+7) = cent_str(ixO^S,1)
   w(ixO^S,nw+8) = cent_str(ixO^S,2)
   w(ixO^S,nw+9) = cent_str(ixO^S,3)
   call get_coriolis(ixI^L,ixO^L,x,w,cor)
   w(ixO^S,nw+10) = cor(ixO^S,1)
   w(ixO^S,nw+11) = cor(ixO^S,2)
   w(ixO^S,nw+12) = cor(ixO^S,3)
   w(ixO^S,nw+13) = cent_str(ixO^S,1)+cor(ixO^S,1)
   w(ixO^S,nw+14) = cent_str(ixO^S,2)+cor(ixO^S,2)
   w(ixO^S,nw+15) = cent_str(ixO^S,3)+cor(ixO^S,3)
   call hd_get_pthermal(w,x,ixI^L,ixO^L,pthermal)
   call gradient(pthermal,ixI^L,ixO^L,1,gradP)
   w(ixO^S,nw+16) = -gradP(ixO^S)/w(ixO^S, rho_)
   call gradient(pthermal,ixI^L,ixO^L,2,gradP)
   w(ixO^S,nw+17) = -gradP(ixO^S)/w(ixO^S, rho_)
   w(ixO^S,nw+18) = ((w(ixO^S, mom(2))**2.d0+w(ixO^S, mom(3))**2.d0) / w(ixO^S, rho_)**2.d0) / x(ixO^S,1)
   w(ixO^S,nw+19) = w(ixO^S,nw+18) - 1.d0 / x(ixO^S,1)**2.d0 + w(ixO^S,nw+16)
   w(ixO^S,nw+20) = w(ixO^S,nw+19) + w(ixO^S,nw+13)
   w(ixO^S,nw+21) = w(ixO^S,nw+18)+cent_str(ixO^S,1)+cor(ixO^S,1)
   end subroutine specialvar_output
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='vr vt vp vp_inrtl v'
    varnames=trim(varnames)//' g ctrr ctrt ctrp crr'
    varnames=trim(varnames)//' crt crp nninr nnint nninp'
    varnames=trim(varnames)//' _gdP_r_ovr_rh _gdP_t_ovr_rh'
    varnames=trim(varnames)//' cnt_gm_r all_bt_nn_in_r all v2nnin'
  end subroutine specialvarnames_output
  ! -------------------------------------------------------------------------------

end module mod_usr
