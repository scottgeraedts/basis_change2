! z_function_m.f90
! (c) F. Duncan M. Haldane, March 2016
! haldane@princeton.edu
!   v0.1   2016-03-08

!  A collection of code  for the
!  "z-function" (Weiestrass times Gaussian) 
!
! contents
!
!   module z_function_m
!   subroutine set_L(norb,l1,l2)
!   subroutine get_L(norb,l1,l2)
!   subroutine get_g(norb,g1,g2)
!   subroutine get_background_factor(factor)
!   subroutine get_log_lattice_z(norb,m,n,logz,normalization)
!   function lattice_z(norb,m,n,l1,l2,use_table) result(z)
!   subroutine z_function_with_modular_transform(x,y,l1,l2,rationalize,denom,zz,sl2z)
!   subroutine optimize_tau(tau,new_tau,sl2z)
!   subroutine get_z(m,n,zmn)
!   subroutine z_function(x,y,l1,l2,rationalize,denom,zz)
!   subroutine setup_z_function_table
!==================================
module z_function_m
integer, parameter :: dp = kind(1.0d0)


!area = (conjg(l1)*l2 - l1*conjg(l2))/2
!g1 = l2/area
!g2 = -l1/area

! conjg(g1)*l2  + g1*conjg(l2) = 0
! conjg(g2)*l1  + g2*conjg(l1) = 0
! (conjg(g1)*l1 +   g1*conjg(l1))/2 = (-conjg(l2)*l1 + l2*conjg(l1))/(2*area) = 1
! (conjg(g2)*l2 +   g2*conjg(l2))/2  = (conjg(l1)*l2 - l1*conjg(l2))/(2*area) = 1


complex (kind=dp), allocatable :: z_function_table(:,:)
real(kind=dp) :: background_factor
complex (kind=dp), allocatable :: log_z_function_table(:,:)
complex(kind=dp) :: l(2) = (/ (0,0) , (0,0) /)
complex(kind=dp) :: g(2) = (/ (0,0) , (0,0) /)
integer:: norb_z = 0

real(kind=dp) :: pi = 1.0_dp

end module z_function_m

!==================================
subroutine set_L(norb,l1,l2)
  use z_function_m
  implicit none
  integer, intent(in) :: norb
  complex(kind=dp), intent(in) :: l1,l2
  complex (kind=dp) :: area
  integer :: i,j,id, k, norb2

  if(pi == 1_dp) pi = 2*asin(pi)


  l(1) = l1
  l(2) = l2
  norb_z = norb
  area = (conjg(l1)*l2 - l1*conjg(l2))/2
  g(1) = l(2)/area
  g(2) = -l(1)/area
  g = 2*pi*g

  return
end subroutine set_L
!==================================

subroutine get_L(norb,l1,l2)
  use z_function_m
  implicit none
  complex(kind=dp), intent(out) :: l1,l2
  integer, intent(out) :: norb
  l1 = l(1) 
  l2 = l(2) 
  norb = norb_z
  return
end subroutine get_L

!==================================
subroutine get_g(norb,g1,g2)
  use z_function_m
  implicit none
  complex(kind=dp), intent(out) :: g1,g2
  integer, intent(out) :: norb



  g1 = g(1) 
  g2 = g(2) 
  norb = norb_z
  return
end subroutine get_g





!==================================
subroutine get_background_factor(factor)
  use z_function_m
  implicit none
  real(kind=dp), intent(out) :: factor
  factor = background_factor
  return
end subroutine get_background_factor


!==================================
subroutine get_log_lattice_z(norb,m,n,logz,normalization)
  use z_function_m
  implicit none
  integer, intent(in) :: norb, m, n
  complex(kind=dp), intent(out) :: logz
  real(kind=dp), intent(out):: normalization
  !-----------------------------------------
  ! uses the table set up by SETUP_Z_FUNCTION_TABLE
  ! returns the logarithm of "lattice Z-function" LOGZ(M,N)
  ! for the NORB x NORB lattice.
  !
  !  Z(M,N) = NORMALIZATION * exp (LOGZ(M,N))
  !
  !  NORMALIZATION is sum'(j,k /= 0,0) (log(Z(j,k))/(norb**2) 
  !
  ! (NORMALISATION == 0, LOGZ = 1 is returned
  !  if M and N are both multiples of NORB)
  !---------------------------------
  complex(kind=dp) :: tau, new_tau
  integer :: sl2z(2,2), j, k, m1, n1, jk

  
  if (norb /= norb_z) then
     write(6,'("ERROR: subroutine GET_LOG_LATTICE_Z: &
          &no compatible table created by SETUP_Z_FUNCTION_TABLE&
          &was found")')
     write(6,'(" norb, norb(table) = ",i12,12x,i12)') norb,norb_z
     stop
  endif
  
  if(mod(m,norb) == 0 .and. mod(n,norb) == 0) then
     normalization = 0_dp
     logz = cmplx(0,kind=dp)
     return
  endif

  normalization = background_factor
  
  if(abs(m) <= norb .and. abs(n) <= norb) then
     logz = log_z_function_table(m,n)
     return
  else
     m1 = modulo(m,norb)
     if(m1 > norb) m1 = m1 - norb
     n1 = modulo(n,norb)
     if(n1 > norb) n1 = n1 - norb
     j = (m-m1)/norb
     k = (n-n1)/norb
     logz = log_z_function_table(m1,n1)
     jk = modulo(norb*(j + k + j*k) + j*n1 - k*m1,2*norb)
     if (jk > norb) jk = jk - 2*norb
     logz = logz + cmplx(0,jk*pi/norb,kind=dp)
     return
  endif
  return
end subroutine get_log_lattice_z

        

!==================================


function lattice_z(norb,m,n,l1,l2,use_table) result(z)
  use z_function_m
  implicit none
  complex (kind=dp) :: z
  integer, intent(in) :: norb, m, n
  complex (kind=dp), intent(in) :: l1,l2
  logical, intent(in) :: use_table
!-----------------------------------------
! returns the "lattice Z-function" z(M,N)
! for the lattice defined by NORB, L1, L2.
!
! this is sqrt(conjg(L1)/L1) * sigma(z) exp (-conjg(z)*z/(2*norb))
! where z - (M*L1 + N*L2)/NORB
! conjg(L1)*L2 - conjg(L2)*L1 = cmplx(0,2*pi) * NORB
! and sigma(z) = sigma(z | {i*L1 + j*L2})
! is the Weierstrass sigma function for the complex
! lattice with oriented basis (L1,L2)
!---------------------------------
  complex(kind=dp) :: tau, new_tau
  integer :: sl2z(2,2)
  logical, parameter :: rationalize = .true.
  real (kind=dp) :: x,y


  if (use_table) then
     if (norb /= norb_z .or. l1 /= l(1) .or. l2 /=l(2)) then
        write(6,'("ERROR: function LATTICE_Z &
             &called with USE_TABLE = .true. but no compatible &
             &table created by SETUP_Z_FUNCTION_TABLE was found")')
        write(6,'(" norb, norb(table) = ",i12,12x,i12)') norb,norb_z
        write(6,'(" L1 , L1(table)    = ",4e12.3)') l1 ,l(1)
        write(6,'(" L2, L2(table)     = ",4e12.3)') l2 ,l(2)        
        stop
     endif
     if(abs(m) > norb .or. abs(n) > norb) then
        write(6,'("ERROR: function LATTICE_Z called with &
             &USE_TABLE = .true. , but m =",i5," and/or n =",i5,&
             & " outside  range |m|, |n| <= norb =",i12)') norb
        stop
     endif
     z = z_function_table(m,n)
     return
  endif
  
  tau = l2/l1
  call optimize_tau(tau,new_tau,sl2z)
  x = real(m,kind=dp)/norb
  y = real(n,kind=dp)/norb
  call z_function_with_modular_transform(x,y,l1,l2,rationalize,norb,z,sl2z)
  return
end function lattice_z


!==================================
subroutine z_function_with_modular_transform(x,y,l1,l2,rationalize,denom,zz,sl2z)
implicit none
integer, parameter :: dp = kind(1.0d0)
integer, intent(in) :: sl2z(2,2)
real(kind=dp), intent (in) :: x, y
complex (kind=dp), intent(in) :: l1,l2
complex (kind=dp), intent(out) :: zz
logical, intent(in) :: rationalize
integer, intent(in) :: denom
!-----------------------------------------
!
!  use a modular transformation to aid in evaluating
!  the Z_function for a given tau.
!
!  the optimal choice can be  obatined by a call to 
!  subroutine optimize_tau(tau,new_tau,sl2z)
!
! 
! This calls subroutine  z_function(m,n,n_lattice,tau,zz)
!----------------------------------------------
complex (kind=dp) :: new_l1, new_l2, phase_factor
real (kind=dp) :: x1,y1
integer :: m,n

if (sl2z(1,1)*sl2z(2,2) - sl2z(1,2)*sl2z(2,1) /= 1) then
   write(6,'("Z-function-modular-transform: invalid SL2Z:")')
   write(6,'("determinant =/ 1:")')
   write(6,'(2i6)') sl2z(1,:), sl2z(2,:)
endif

if (rationalize) then
   m = nint(x*denom)
   n = nint(y*denom)
   x1 = real(sl2z(2,2)*m - sl2z(2,1)*n,kind=dp)/denom
   y1 = real(-sl2z(1,2)*m + sl2z(1,1)*n,kind=dp)/denom
else
   x1 =  sl2z(2,2)*x  - sl2z(2,1)*y
   y1 = -sl2z(1,2)*x  + sl2z(1,1)*y
endif
   new_l1 = sl2z(1,1)*l1 + sl2z(1,2)*l2 
   new_l2 = sl2z(2,1)*l1 + sl2z(2,2)*l2 

call z_function(x1,y1,new_l1,new_l2,rationalize,denom,zz)

phase_factor = sqrt(new_l1/conjg(new_l1))
zz = phase_factor*zz

return
end subroutine z_function_with_modular_transform


!==================================
subroutine optimize_tau(tau,new_tau,sl2z)
implicit none
integer, parameter :: dp = kind(1.0d0)
complex (kind=dp), intent(in) :: tau
complex (kind=dp), intent(out) :: new_tau
integer, intent(out) :: sl2z(2,2)

!------------------------------------------------------
!  cary out a series of modular transformations to 
!  find the most isotropic value of the
!  modular parameter  new_tau = L2/L1
!  where L1 is the shortest and L2 is the next shortest
!  lattice translation, with  aimag(L2/L1) > 0
!
!  the optimized tau has  
!
!   aimag(tau)**2  >= 1 - real(tau)**2 >=  0.75
!    real(tau)**2  <= 0.25
!  
!  new_tau =  (sl2z(2,1) + sl2z(2,2)*tau)
!              --------------------------
!             (sl2z(1,1) + sl2z(1,2)*tau)
!
!  with det | sl2z | = 1
! 
!----------------------------------------------------------
real (kind=dp) :: test
integer ::  k, mod_t(2,2), mod_s(2,2)
real (kind=dp), parameter :: close_to_one = 0.9999_dp

mod_t(1,1)  =  1
mod_t(1,2)  =  0
mod_t(2,1)  =  0
mod_t(2,2)  =  1

mod_s(1,1)  =  0
mod_s(1,2)  =  1
mod_s(2,1)  = -1
mod_s(2,2)  =  0


if(abs(tau) > 1_dp) then
   new_tau = tau
   sl2z = mod_t
else
   new_tau = -1/tau
   sl2z = mod_s
endif

do
   k = nint(real(new_tau))
   if (k /= 0) then
      new_tau = new_tau - k
      mod_t(2,1) = -k
      sl2z = matmul(mod_t,sl2z)
   endif
   test = abs(new_tau)
   if(test > 1_dp) exit
   new_tau = -1/new_tau
   sl2z = matmul(mod_s,sl2z)
! to avoid cycling when the shortest period is
! degenerate.
   if(test > close_to_one) exit
enddo

return
end subroutine optimize_tau

!==================================
subroutine setup_z_function_table
use z_function_m
implicit none

real (kind=dp) :: area, u, a, pin, background_sum,x,y

integer :: sl2z(2,2), norb,i,j, norb2
complex (kind=dp) :: tau, new_tau
logical, parameter :: rationalize = .true.


norb = norb_z
area = aimag(conjg(l(1))*l(2))
if (nint(area) == 0) then
   write(6,'(" MAKE_Z_FUNCTION_TABLE: geometry is not set")')
   stop
endif

tau = l(2)/l(1)
call  optimize_tau(tau,new_tau,sl2z)  

if (allocated(z_function_table)) deallocate(z_function_table)
allocate (z_function_table(-norb:norb,-norb:norb))
if (allocated(log_z_function_table)) deallocate(log_z_function_table)
allocate (log_z_function_table(-norb:norb,-norb:norb))

do i = -1,1
   do j = -1,1
      z_function_table(i*norb,j*norb) = cmplx(0,kind=dp)
   enddo
enddo
do i  = -norb , norb
   do j = 0, norb
      if(mod(i,norb) == 0 .and. mod(j,norb) == 0) cycle
      if (j == 0  .and. i <=  0) cycle
      x = real(i,kind=dp)/norb
      y = real(j,kind=dp)/norb      
      call z_function_with_modular_transform(x,y,l(1),l(2), &
           & rationalize,norb,z_function_table(i,j),sl2z)
      z_function_table(-i,-j) = -z_function_table(i,j)
   enddo
enddo



! coulomb-compensated log(z_factor)
! express z_function(i,j) =  background_factor * exp (log_z_function(i,j))
! where sum'(x /= 0) log_z_function(x) = 0

do i = -norb,norb
   do j = -norb,norb
      if(mod(i, norb) == 0 .and. mod(j, norb) == 0) then 
         log_z_function_table(i,j) = cmplx(0,kind=dp)
      else
         log_z_function_table(i,j) = log(z_function_table(i,j))
      endif
   enddo
enddo
background_sum = sum (real(log_z_function_table(1:norb,1:norb)))/(norb**2)
log_z_function_table = log_z_function_table - background_sum
background_factor = exp(background_sum)
return
end subroutine setup_z_function_table




!==================================

subroutine get_z(m,n,zmn)
  use z_function_m
  implicit none
  integer, intent(in) :: m,n
  complex(kind=dp), intent(out) :: zmn

! l1 = zmn for (m,n) = (norb,0)                                                                                                                                                  
! l2 = zmn for (m,n) = (0,norb)                                                                                                                                                  
!                                                                                                                                                                                
  zmn = (m*l(1) + n*l(2))/norb_z
  return
end subroutine get_z

!==================================
subroutine z_function(x,y,l1,l2,rationalize,denom,zz)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real (kind=dp), intent(in) :: x,y 
  complex (kind=dp), intent(in) :: l1,l2
  complex (kind=dp), intent(out) :: zz
  logical, intent(in) :: rationalize
  integer, intent(in) :: denom
  !-------------------------------------------------------------------------------
  ! computes the "Z-function"
  ! 
  !  ZZ = sigma( z | L)  * exp ( - (abs(z)**2)/(2*A) )
  !
  !  
  !   z = x*L1 + y*l2
  !
  !  |conjg(L1)*L2 - conjg(L2)*L1)| = 2*pi*A > 0
  !
  ! sigma(z) is the modifed Weierstrass sigma function of the complex lattice
  !  {m*L1 + n*L2}
  !
  ! for A = norb, x = m/norb, y= n/norb, tau = l2/l1 will return  the
  ! same result as lattice_z_function(m,n,norb,tau,zz)
  !-----------------------------------------------------------------------
  !
  ! if RATIONALIZE == .true.  it is assumed that x and y can be replaced by the rational
  ! numbers   (nint(X*DENOM)/DENOM,   (nint(Y*DENOM)/DENOM)
  !
  !--------------------------------------------------------------------------
  complex (kind=dp) :: q, qk, factor, gamma, z_prev, tau
  real (kind=dp) :: theta, scale, x1,y1,a 
  real(kind=dp), save :: pi = 1_dp
  complex (kind=dp), parameter :: one = (1,0)
  integer :: two_n_lattice, n_lattice2,  n1, m1, j, k, jk, m, n

  if(pi == 1_dp) pi = 2*asin(pi)
  

  scale = aimag(conjg(l1)*l2)  ! |scale| = pi*area

  if(scale == 0_dp) then
     write (6,'("Z_FUNCTION: area must not be zero")') 
     stop
  endif
  

  if (rationalize) then
     m = nint(x*denom)
     n = nint(y*denom)
     m1 = modulo(m,denom)
     n1 = modulo(n,denom)
     if(m1 == 0 .and. n1 == 0) then
        zz = cmplx(0,kind=dp)
        return
     endif
     if(m1 > denom/2) m1 = m1 - denom
     if(n1 > denom/2) n1 = n1 - denom
     x1 = real(m1,kind=dp)/denom
     y1 = real(n1,kind=dp)/denom
     j = (m-m1)/denom
     k = (n-n1)/denom
     jk = modulo(j*n1-k*m1,2*denom)
     if (jk > denom) jk = jk - 2*denom
     theta = jk*pi/denom
  else
     a = 1_dp
     x1 = modulo(x,a)
     y1 = modulo(y,a)
     if(x1 == 0_dp .and. y1 == 0_dp) then
        zz = cmplx(0,kind=dp)
        return
     endif
     if (x1 > 0.5_dp) x1 = x1 - 1_dp
     if (y1 > 0.5_dp)  y1 = y1 - 1_dp
     j = nint(x-x1)
     k = nint(y-y1)
     theta = pi*(j*y1 - k*x1)
  endif

  tau = l2/l1  
  scale = abs(scale/aimag(tau))
  zz = sqrt(scale)/cmplx(0,2*pi,kind=dp)
  zz = zz*cmplx(cos(theta),sin(theta),kind=dp)
  if(mod(j,2)/= 0 .or. mod(k,2) /= 0) zz = -zz
  


  gamma =  exp(cmplx(0,pi*y1,kind=dp)*tau)
  theta = pi*x1
  gamma = gamma*cmplx(cos(theta),sin(theta),kind=dp)  
  gamma = gamma - 1/gamma
  
  zz = zz*gamma
  
  theta = pi*y1
  zz = zz *  exp(cmplx(0,y1*theta,kind=dp)*tau)
  theta = x1*theta
  zz = zz * cmplx(cos(theta),sin(theta),kind=dp) 

  q = exp(cmplx(0,pi,kind=dp)*tau)
  qk = one
  do
     qk = qk*q
     factor = (gamma/(qk - (1/qk)))**2
     z_prev = zz
     zz = zz*(one-factor)
     if(real(zz) == real(z_prev) .and. aimag(zz) == aimag(z_prev)) exit
  enddo
  
  return
end subroutine z_function



