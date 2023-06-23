module fun
   use, intrinsic :: iso_c_binding
   use ifcore
   use ifport

   integer(c_int) ne, nt, nz, neqf, neqp, lworkf, lworkp, liworkf, liworkp, nrdf, nrdp, iparf, &
      iparp, ioutf, ioutp, ididf, ididp, itolf, itolp, inharm, aiparp(1), ndtr
   real(c_double) zex_w, zex, dz, tend, dtr(2), q(3), icu(2), th(2), a(2), dcir(2), r(2), f0(3), dt, &
      pitch, f10, f20, f30, p10, p20, p30, rtolf, rtolp, atolf, atolp, rparf, rparp, ftol, ptol, nharm, &
      gamma, ukv, betta, betta2, betta_z, betta_z2, betta_perp, betta_perp2, gmp, w_op_w, w_op, c, e, m, &
      artolp(1), aatolp(1), arparp(1), b(2), w_n(2), ia(2), norm, nfac, &
      dtrb(2), dtrh_w, dtr0
   complex(c_double_complex) fp(2)
   logical(c_bool) wc, fok, lensm, btod, iatoi, inher
   character(len=6) path

   integer(c_int) breaknum(3)
   real(c_double) phitmp0(3), phitmp1(3)
   complex(c_double_complex) fc, fcomp(3)

   integer(c_int), allocatable, target :: idxre(:, :), idxim(:, :), iworkf(:), iworkp(:)
   complex(c_double_complex), allocatable, target :: mean(:), u(:)
   real(c_double), allocatable, target :: tax(:), zax(:), eta(:, :), etag(:, :), w(:, :), f(:, :), w1(:, :), p(:, :), &
                    phi(:, :), phios(:, :), wos(:, :), workf(:), workp(:), cl1(:), lhs1(:), rhs1(:), cl2(:), lhs2(:), rhs2(:)!, u(:)

   include 'z.inc'
   include 're.inc'
   include 'im.inc'

   complex(c_double_complex), parameter :: ic = (0.0d0, 1.0d0)
   real(c_double), parameter :: pi = 2.0D0*dacos(0.0D0)

   private freq_out, zex_w, tend, q, i, th, a, dcir, r, f0, pitch

contains

   subroutine init()
      implicit none

      integer(c_int) ii

      call read_param()

      nharm = dble(inharm)
      gamma = 1.0 + ukv/511.0
      betta = dsqrt(1.0d0 - 1.0d0/(gamma*gamma))
      betta2 = 1.0d0 - 1.0d0/(gamma*gamma)
      betta_z = betta/dsqrt(pitch*pitch + 1.0d0)
      betta_z2 = betta2/(pitch*pitch + 1.0d0)
      betta_perp2 = betta2 - betta_z2
      gmp = 0.048715056967419
      w_op_w = 2*pi*w_op*1e9
      c = 29979245800.0d0
      e = 4.803e-10
      m = 9.1093837015e-28

      if (lensm .eq. .true.) zex_w = betta_perp2/2.0d0/betta_z*w_op_w*zex/nharm/c
      print *, 'Zex = ', zex_w

      if (btod .eq. .true.) then
         w_n(1) = e*B(1)/(m*c)*10000.0d0
         dtr(1) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(1))/(gamma*w_op_w))
         w_n(2) = e*B(2)/(m*c)*10000.0d0
         dtr(2) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(2))/(gamma*w_op_w))
      end if            
      
      call norma(norm)

      if (iatoi .eq. .true.) then         
         nfac = factorial(inharm)         
         icu(1) = 2.35/10000*IA(1)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(1)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
         icu(2) = 2.35/10000*IA(2)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(2)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
      end if      
      
      !print *, 'dtr1 = ', dtr(1), 'dtr2 = ', dtr(2)
      print *, 'icu1 = ', icu(1), 'icu2 = ', icu(2)

      nt = tend/dt + 1
      nz = zex_w/dz + 1

      neqp = 4*ne
      nrdp = 4*ne
      lworkp = 8*neqp + 5*nrdp + 21
      liworkp = nrdp + 21

      neqf = 6
      nrdf = 6
      lworkf = 11*neqf + 8*nrdf + 21
      liworkf = nrdf + 21

      call allocate_arrays()

      f(1, 1) = f10
      f(2, 1) = p10
      f(3, 1) = f20
      f(4, 1) = p20
      f(5, 1) = f30
      f(6, 1) = p30

      do ii = 1, nt
         tax(ii) = (ii - 1)*dt
      end do

      do ii = 1, nz
         zax(ii) = (ii - 1)*dz
      end do

      call calc_u(u, zex_w, nz, zax)

      do ii = 1, 2
         idxre(ii, :) = (/2*(ii - 1)*ne + 1:(2*ii - 1)*ne/)
         idxim(ii, :) = (/(2*ii - 1)*ne + 1:2*ii*ne/)
      end do
      
      if (dtrh_w .ne. 0) then
         ndtr = (dtrb(2) - dtrb(1))/dtrh_w + 1
      else
         ndtr = 1
      end if
      
      dtr0 = dtrb(1)

      !open (1, file='test.dat')
      !do ii = 1, ne
      !    write (1,'(4i)') idxre(1,ii), idxim(1,ii), idxre(2,ii), idxim(2,ii)
      !end do
      !close (1)
      !stop
      
      w(:,:) = 0

   end subroutine init
   
   subroutine init_local(path)
      implicit none

      integer(c_int) ii
      character(*) path

      call read_param_local(path)

      nharm = dble(inharm)
      gamma = 1.0 + ukv/511.0
      betta = dsqrt(1.0d0 - 1.0d0/(gamma*gamma))
      betta2 = 1.0d0 - 1.0d0/(gamma*gamma)
      betta_z = betta/dsqrt(pitch*pitch + 1.0d0)
      betta_z2 = betta2/(pitch*pitch + 1.0d0)
      betta_perp2 = betta2 - betta_z2
      gmp = 0.048715056967419
      w_op_w = 2*pi*w_op*1e9
      c = 29979245800.0d0
      e = 4.803e-10
      m = 9.1093837015e-28

      if (lensm .eq. .true.) zex_w = betta_perp2/2.0d0/betta_z*w_op_w*zex/nharm/c
      print *, 'Zex = ', zex_w

      if (btod .eq. .true.) then
         w_n(1) = e*B(1)/(m*c)*10000.0d0
         dtr(1) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(1))/(gamma*w_op_w))
         w_n(2) = e*B(2)/(m*c)*10000.0d0
         dtr(2) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(2))/(gamma*w_op_w))
      end if                        

      if (iatoi .eq. .true.) then         
         nfac = factorial(inharm)         
         icu(1) = 2.35/10000*IA(1)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(1)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
         icu(2) = 2.35/10000*IA(2)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(2)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
      end if      
      
      print *, 'dtr1 = ', dtr(1), 'dtr2 = ', dtr(2)      

      nt = tend/dt + 1
      nz = zex_w/dz + 1

      neqp = 4*ne
      nrdp = 4*ne
      lworkp = 8*neqp + 5*nrdp + 21
      liworkp = nrdp + 21

      neqf = 6
      nrdf = 6
      lworkf = 11*neqf + 8*nrdf + 21
      liworkf = nrdf + 21      

      f(1, 1) = f10
      f(2, 1) = p10
      f(3, 1) = f20
      f(4, 1) = p20
      f(5, 1) = f30
      f(6, 1) = p30

      do ii = 1, nt
         tax(ii) = (ii - 1)*dt
      end do

      do ii = 1, nz
         zax(ii) = (ii - 1)*dz
      end do

      call calc_u(u, zex_w, nz, zax)

      do ii = 1, 2
         idxre(ii, :) = (/2*(ii - 1)*ne + 1:(2*ii - 1)*ne/)
         idxim(ii, :) = (/(2*ii - 1)*ne + 1:2*ii*ne/)
      end do
      
   end subroutine init_local
   
   subroutine norma(norm)
      use, intrinsic :: iso_c_binding, only: c_double, c_double_complex, c_int
      import, only:rea, ima, betta_perp2, betta_z, inharm, c, w_op_w!, omega
      implicit none

      real(c_double) :: dzz = 0.0280211, norm
      complex(c_double_complex) :: u(663)

      dzz = betta_perp2/2.0d0/betta_z*w_op_w*dzz/inharm/c
      u = dcmplx(rea, ima)

      norm = sum(cdabs(u(:))*cdabs(u(:)))*dzz

      print *, 'N = ', norm

   end subroutine norma
   
   recursive function factorial(p) result(l)
      integer, intent(in) :: p
      integer l
      if (p == 1) then
         l = 1
      else
         l = p*factorial(p - 1)
      end if
   end function

   function squval(zz)

      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) squval
      real(c_double) z, re, im, dz, z1
      integer(c_int) l

      dz = 0.280211
      z = zz/zex_w*185.5d0
      l = z/dz

      if (l .eq. 0) then
         l = 2
      elseif (l .ge. 662) then
         l = 662
      else
         if ((z - za(l)) .gt. 0.5*dz) l = l + 1
      end if

      z1 = za(l)
      z = z - 8.5d0

      re = rea(l - 1) + ((-rea(l - 1) + rea(l))*(dz + z - z1))/dz + ((rea(l - 1)/2.0d0 - rea(l) + &
                                                                      rea(l + 1)/2.0d0)*(z - z1)*(dz + z - z1))/dz/dz
      im = ima(l - 1) + ((-ima(l - 1) + ima(l))*(dz + z - z1))/dz + ((ima(l - 1)/2.0d0 - ima(l) + &
                                                                      ima(l + 1)/2.0d0)*(z - z1)*(dz + z - z1))/dz/dz
      !!!NE RABOTAET
      !re = ((rea(l - 1) - 2*rea(l) + rea(l + 1))*z**2)/(2.*dz**2) &
      !     + (z*(dz*(-rea(l - 1) + rea(l + 1)) - 2*(rea(l - 1) - 2*rea(l) + rea(l + 1))*z1))/(2.*dz**2) + &
      !     -(2*dz**2*rea(l) + dz*(rea(l - 1) - rea(l + 1))*z1 + (rea(l - 1) - 2*rea(l) + rea(l + 1))*z1**2)/(2.*dz**2)
      !im = ((ima(l - 1) - 2*ima(l) + ima(l + 1))*z**2)/(2.*dz**2) + &
      !     (z*(dz*(-ima(l - 1) + ima(l + 1)) - 2*(ima(l - 1) - 2*ima(l) + ima(l + 1))*z1))/(2.*dz**2) + &
      !     -(2*dz**2*ima(l) + dz*(ima(l - 1) - ima(l + 1))*z1 + (ima(l - 1) - 2*ima(l) + ima(l + 1))*z1**2)/(2.*dz**2)

      squval = dcmplx(re, im)

   end function squval

   function uval(zz)

      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) uval
      real(c_double) z, re, im, d
      integer(c_int) l

      z = zz/zex_w*185.5 - 8.5
      l = (z + 8.5)/0.28021 + 1
      d = z - za(l)

      !print *, z, l, d

      if (d .gt. 0.0 .and. l /= 663) then
         re = (rea(l)*za(l + 1) - rea(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (rea(l + 1) - rea(l))/(za(l + 1) - za(l))*z
         im = (ima(l)*za(l + 1) - ima(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (ima(l + 1) - ima(l))/(za(l + 1) - za(l))*z
      else if (d .lt. 0.0 .and. l /= 1) then
         re = (rea(l - 1)*za(l) - rea(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (rea(l) - rea(l - 1))/(za(l) - za(l - 1))*z
         im = (ima(l - 1)*za(l) - ima(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (ima(l) - ima(l - 1))/(za(l) - za(l - 1))*z
      else
         re = rea(l)
         im = ima(l)
      end if

      uval = dcmplx(re, im)

   end function uval

   subroutine allocate_arrays()
      use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int) err_alloc

      allocate (f(6, nt), p(4*ne, nz), u(nz), tax(nt), zax(nz), mean(nz), eta(2, nt), etag(2, nt), w(3, nt - 1), w1(3, nt - 1), &
                idxre(2, ne), idxim(2, ne), workf(lworkf), iworkf(liworkf), workp(lworkp), iworkp(liworkp), &
                wos(3, nt - 1), phi(3, nt), phios(3, nt), cl1(nt), lhs1(nt), rhs1(nt), cl2(nt), lhs2(nt), rhs2(nt), &
                stat=err_alloc)

      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if
   end subroutine allocate_arrays

   subroutine deallocate_arrays()
      use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int) err_dealloc

      deallocate (f, p, u, tax, zax, mean, eta, etag, w, w1, stat=err_dealloc)

      if (err_dealloc /= 0) then
         print *, "deallocation error"
         pause
         stop
      end if
   end subroutine deallocate_arrays   
     
   subroutine read_param_local(path) bind(c, name='read_param_local')
      use, intrinsic :: iso_c_binding
      import
      implicit none

      namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
         dcir1, dcir2, r1, r2, f10, f20, f30, p10, p20, p30, dt, dz, pitch, ftol, ptol, wc, fok, inharm, ukv, &
         w_op, lensm, btod, b1, b2, iatoi, ia1, ia2, &
         dtr1, dtr2

      real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dcir1, dcir2, r1, r2, b1, b2, ia1, ia2, dtr1, dtr2
      character(*) path

      open (unit=1, file=path//'input_fortran.in', status='old', err=101)
      read (unit=1, nml=param, err=102)
      close (unit=1)

      q(1) = q1
      q(2) = q2
      q(3) = q3
      icu(1) = i1
      icu(2) = i2
      th(1) = th1
      th(2) = th2
      a(1) = a1
      a(2) = a2
      dtr(1) = dtr1
      dtr(2) = dtr2
      dcir(1) = dcir1
      dcir(2) = dcir2
      r(1) = r1
      r(2) = r2
      b(1) = b1
      b(2) = b2
      ia(1) = ia1
      ia(2) = ia2

      write (*, nml=param)

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of reading file "input_fortran.in"'; pause; stop
   end subroutine read_param_local   

   subroutine read_param() bind(c, name='read_param')
      use, intrinsic :: iso_c_binding
      import
      implicit none

      namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
         dcir1, dcir2, r1, r2, f10, f20, f30, p10, p20, p30, dt, dz, pitch, ftol, ptol, wc, fok, inharm, ukv, &
         w_op, lensm, btod, b1, b2, iatoi, ia1, ia2, &
         dtrb, dtrh, inher

      real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, b1, b2, ia1, ia2, dtrh

      open (unit=1, file='input_fortran.in', status='old', err=101)
      read (unit=1, nml=param, err=102)
      close (unit=1)

      q(1) = q1
      q(2) = q2
      q(3) = q3
      icu(1) = i1
      icu(2) = i2
      th(1) = th1
      th(2) = th2
      a(1) = a1
      a(2) = a2
      !dtr(1) = dtr1
      !dtr(2) = dtr2
      dcir(1) = dcir1
      dcir(2) = dcir2
      r(1) = r1
      r(2) = r2
      b(1) = b1
      b(2) = b2
      ia(1) = ia1
      ia(2) = ia2
      dtrh_w = dtrh

      !write (*, nml=param)

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of reading file "input_fortran.in"'; pause; stop
   end subroutine read_param

   subroutine write_param(path) bind(c, name='write_param')
      use, intrinsic :: iso_c_binding
      import
      implicit none

      namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
         dcir1, dcir2, r1, r2, f10, f20, f30, p10, p20, p30, dt, dz, pitch, ftol, ptol, wc, fok, inharm, ukv, &
         w_op, lensm, btod, b1, b2, iatoi, ia1, ia2, &
         dtr1, dtr2

      character(*) path
      real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dcir1, dcir2, r1, r2, b1, b2, ia1, ia2, dtr1, dtr2    
      
      dtr1 = dtr(1)
      dtr2 = dtr(2)
      i1 = icu(1)
      i2 = icu(2)
      q1 = q(1)
      q2 = q(2)
      q3 = q(3)
      th1 = th(1)
      th2 = th(2)
      a1 = a(1)
      a2 = a(2)
      dcir1 = dcir(1)
      dcir2 = dcir(2)
      r1 = r(1)
      r2 = r(2)
      b1 = b(1)
      b2 = b(2)
      ia1 = ia(1)
      ia2 = ia(2)      

      open (unit=1, file=path//'input_fortran.in', err=101)     
      write (1, nml=param)
      close (unit=1)

      !write (*, nml=param)

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of reading file "input_fortran.in"'; pause; stop
   end subroutine write_param

   subroutine write_results()

      implicit none

      integer i, j
      
      !if (wc .eq. .true.) then
      !   w(:, 1) = 0
      !   do i = 2, nt
      !      do j = 1, 3
      !         !w(j, i - 1) = dimag(log(f(2*j - 1, i)*cdexp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*cdexp(ic*f(2*j, i - 1)))))/dt
      !         w(j, i) = (f(2*j, i) - f(2*j, i - 1))/dt
      !      end do
      !   end do
      !   print *, 'Frequency calculated from phase. ( WC = ', wc, ')'
      !elseif (wc .eq. .false.) then
         call freq()
         print *, 'Frequency calculated from RHS. ( WC = ', wc, ')'
      !end if

      phi(:, 1) = 0; 
      do i = 2, nt
         do j = 1, 3
            phi(j, i) = phi(j, i - 1) + dimag(log(f(2*j - 1, i)*cdexp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*cdexp(ic*f(2*j, i - 1)))))
         end do
      end do

      breaknum(:) = 0
      fcomp(1) = f(2*1 - 1, 1)*cdexp(ic*f(2*1, 1))
      fcomp(2) = f(2*2 - 1, 1)*cdexp(ic*f(2*2, 1))
      fcomp(3) = f(2*3 - 1, 1)*cdexp(ic*f(2*3, 1))
      phitmp0(:) = datan2(dimag(fcomp(:)), dreal(fcomp(:)))
      !phitmp0(:) = datan2(dimag(f(:, 1)), dreal(f(:, 1)))
      phios(:, 1) = phitmp0(:)
      do i = 2, nt
         do j = 1, 3
            fc = f(2*j - 1, i)*cdexp(ic*f(2*j, i))
            phitmp1(j) = datan2(dimag(fc), dreal(fc))
            if ((phitmp1(j) - phitmp0(j)) .gt. pi) breaknum(j) = breaknum(j) - 1
            if ((phitmp1(j) - phitmp0(j)) .lt. -pi) breaknum(j) = breaknum(j) + 1
            phios(j, i) = phitmp1(j) + 2.*pi*breaknum(j)
            !phios(j, i) = phitmp1(j)
            phitmp0(j) = phitmp1(j)
         end do
      end do

      do i = 1, nt - 1
         do j = 1, 3
            wos(j, i) = (phios(j, i + 1) - phios(j, i))/dt
         end do
      end do

      write (*, '(/)')

      !pause

      open (3, file=path//'cl1.dat')
      do i = 1, nt
         write (3, '(5f12.6,a)') tax(i), cl1(i), lhs1(i), rhs1(i), abs(cl1(i)/lhs1(i))*100,' %'
      end do
      close (3)

      open (3, file=path//'cl2.dat')
      do i = 1, nt
         write (3, '(5f12.6,a)') tax(i), cl2(i), lhs2(i), rhs2(i), abs(cl2(i)/lhs2(i))*100, ' %'
      end do
      close (3)

      open (1, file=path//'F.dat')
      do i = 1, nt
         !write (1, '(4e17.8)') tax(i), dabs(f(1, i)), dabs(f(3, i)), dabs(f(5, i))
         write (1, '(4f12.6)') tax(i), f(1, i), f(3, i), f(5, i)
      end do
      close (1)

      open (13, file=path//'FCMPLX.dat')
      do i = 1, nt
         fcomp(1) = f(2*1 - 1, i)*cdexp(ic*f(2*1, i))
         fcomp(2) = f(2*2 - 1, i)*cdexp(ic*f(2*2, i))
         fcomp(3) = f(2*3 - 1, i)*cdexp(ic*f(2*3, i))
         write (13, '(7f12.6)') tax(i), dreal(fcomp(1)), dimag(fcomp(1)), dreal(fcomp(2)), dimag(fcomp(2)), &
            dreal(fcomp(3)), dimag(fcomp(3))
      end do
      close (13)

      open (2, file=path//'E.dat')
      do i = 1, nt
         write (2, '(5f12.6)') tax(i), eta(1, i), etag(1, i), eta(2, i), etag(2, i)
      end do
      close (2)

      open (3, file=path//'W.dat')
      do i = 1, nt
         write (3, '(4f12.6)') tax(i), w(1, i), w(2, i), w(3, i)
      end do
      close (3)

      open (1, file=path//'P.dat')
      do i = 1, nt
         !write (1, '(4e17.8)') tax(i), phi(1, i), phi(2, i), phi(3, i)
         write (1, '(4f12.6)') tax(i), f(2, i), f(4, i), f(6, i)
      end do
      close (1)

      open (1, file=path//'POS.dat')
      do i = 1, nt
         write (1, '(4f12.6)') tax(i), phios(1, i), phios(2, i), phios(3, i)
      end do
      close (1)

      open (3, file=path//'WOS.dat')
      do i = 1, nt - 1
         write (3, '(4f12.6)') tax(i + 1), wos(1, i), wos(2, i), wos(3, i)
      end do
      close (3)

      call write_param(path)

      return

101   print *, 'error of file open.'
      pause
      stop
102   print *, 'error of file reading.'
      pause
      stop
103   print *, 'error of file writing.'
      pause
      stop
   end subroutine write_results

   subroutine ode4f()
      import
      implicit none

      integer(c_int) i, j, aiparf(1), itp, itf!, aiparp(1)
      real(c_double) :: t, z, artolf(1), aatolf(1), arparf(1), xoutp, xoutf!, artolp(1), aatolp(1), arparp(1)
      real(c_double) yf(6), yp(neqp), pex(neqp)
      logical(4) pressed
      character(1) key
      integer(c_int), parameter :: esc = 27
      common/internp/xoutp, itp
      common/internf/xoutf, itf

      !call d02pvf(neqp, zstart, p(:, 1), zex_w, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrkp, ifail)

      !solve eq. at t=0
      fp(1) = f(1, 1)*cdexp(ic*f(2, 1))
      fp(2) = f(3, 1)*cdexp(ic*f(4, 1))

      call init_dopr_p()
      ioutp = 2
      z = zax(1)
      xoutp = z
      itp = 0
      yp = p(:, 1)            
      iworkp(5) = neqp  
      
      !rparp = 0.0
      !iparp = 0
      !itolp = 0
      !!rtolp = 1.0d-7
      !rtolp = ptol
      !atolp = rtolp
      !ioutp = neqp
      !z = zax(1)
      !xoutp = z
      !itp = 0
      !yp = p(:, 1)
      !iworkp(:) = 0
      !workp(:) = 0.0d0
      !iworkp(5) = neqp
      !if (nz > 0) then
      !   workp(6) = dz
      !end if

      !artolp(1) = rtolp
      !aatolp(1) = atolp
      !arparp(1) = rparp
      !aiparp(1) = iparp

      call dopri5_p(neqp, dpdz, z, yp, zex_w, artolp, aatolp, itolp, soloutp, ioutp, &
                    workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

      p(:, nz) = yp(:)

      eta(:, 1) = eff(p(:, nz))
      etag(:, 1) = pitch**2/(pitch**2 + 1)*eta(:, 1)

      rparf = 0.0
      iparf = 0
      t = tax(1)
      xoutf = t
      itf = 0
      yf = f(:, 1)
      itolf = 0
      !rtolf = 1.0d-7
      rtolf = ftol
      atolf = rtolf
      ioutf = 6

      artolf(1) = rtolf
      aatolf(1) = atolf
      arparf(1) = rparf
      aiparf(1) = iparf

      iworkf(:) = 0
      workf(:) = 0.0d0

      iworkf(5) = 6

      call dopri5_f(6, dfdt, t, yf, tend, artolf, aatolf, itolf, soloutf, ioutf, &
                    workf, lworkf, iworkf, liworkf, arparf, aiparf, ididf)

      do j = 1, neqf
         f(j, nt) = yf(j)
      end do
      call calcpex(f(:, nt), pex, cl1(nt), lhs1(nt), rhs1(nt), cl2(nt), lhs2(nt), rhs2(nt))
      eta(:, nt) = eff(pex)
      !eta(:, nt) = eff(p(:, nz))
      etag(:, nt) = pitch**2/(pitch**2 + 1)*eta(:, nt)

   end subroutine ode4f

   function eff(pex) result(eta)
      use, intrinsic :: iso_c_binding, only: c_double, c_int
      import, only:ne, idxre, idxim

      implicit none

      integer(c_int) i
      real(c_double) eta(2)
      real(c_double), intent(in) :: pex(:)

      do i = 1, 2
         eta(i) = 1 - sum(abs(dcmplx(pex(idxre(i, :)), pex(idxim(i, :))))**2)/ne
      end do
   end function eff

   subroutine dpdz(neqp, z, p, prhs, rparp, iparp)
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
      import, only:ne, zex_w, f, ic, dtr, fok, squval, idxre, idxim, fp, inharm
      implicit none

      real(c_double) z, p(*), prhs(*)

      integer(c_int) i, reidx(ne), imidx(ne), neqp, iparp
      real(c_double) rparp
      complex(c_double_complex) s(ne), ptmp(ne), u

      if (fok .eq. .false.) then
         u = dexp(-3*((z - zex_w/2)/(zex_w/2))**2)
      else
         u = squval(z)
      end if

      do i = 1, 2
         ptmp = dcmplx(p(idxre(i, :)), p(idxim(i, :)))

         !s = ic*(fp(i)*u*dconjg(ptmp) - (dtr(i) + cdabs(ptmp)**2 - 1)*ptmp)
         s = ic*(fp(i)*u*dconjg(ptmp)**(inharm - 1) - (dtr(i) + cdabs(ptmp)**2 - 1)*ptmp)

         prhs(idxre(i, :)) = dreal(s)
         prhs(idxim(i, :)) = dimag(s)
      end do
   end subroutine dpdz

   complex(c_double_complex) function xi(p, num)
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
      import, only:ne, nz, mean, u, dz, idxre, idxim, inharm

      implicit none

      integer(c_int) i, num
      real(c_double), intent(in) :: p(:, :)

      do i = 1, nz
         !mean(i) = sum(dcmplx(p(idxre(num, :), i), p(idxim(num, :), i))**2, 1)/ne
         mean(i) = sum(dcmplx(p(idxre(num, :), i), p(idxim(num, :), i))**inharm, 1)/ne
      end do

      mean = dconjg(u)*mean

      xi = (0.5d0*(mean(1) + mean(nz)) + sum(mean(2:nz - 1)))*dz
   end function

   subroutine init_dopr_p()
      import, only:rparp, iparp, itolp, rtolp, atolp, iworkp, workp, ptol, nz, dz, &
         artolp, aatolp, arparp, aiparp
      implicit none

      rparp = 0.0
      iparp = 0
      itolp = 0
      !rtolp = 1.0d-7
      rtolp = ptol
      atolp = rtolp
      iworkp(:) = 0
      workp(:) = 0.0d0      
      
      artolp(1) = rtolp
      aatolp(1) = atolp
      arparp(1) = rparp
      aiparp(1) = iparp
   end subroutine init_dopr_p

   subroutine dfdt(neqf, t, f, s, rparf, iparf)
      !use odep, only : dopri5
      implicit none

      integer(c_int) :: ii, jj, neqf, iparf, itp!, aiparp(1)
      real(c_double) t, z, f(neqf), yp(neqp), s(neqf), xoutp, &
         x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
         x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
         f1, f2, f3, phi1, phi2, phi3, a1, a2, rparf
      complex(c_double_complex) x1, x2
      common/internp/xoutp, itp

      !call d02pvf(neqp, zstart, p(:, 1), zex_w, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifail)

      fp(1) = f(1)*exp(ic*f(2))
      fp(2) = f(3)*exp(ic*f(4))

      call init_dopr_p()
      ioutp = 2
      z = zax(1)
      xoutp = z
      itp = 0
      yp = p(:, 1)
      iworkp(5) = neqp

      !artolp(1) = rtolp
      !aatolp(1) = atolp
      !arparp(1) = rparp
      !aiparp(1) = iparp

      call dopri5_p(neqp, dpdz, z, yp, zex_w, artolp, aatolp, itolp, soloutp, ioutp, &
                    workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

      p(:, nz) = yp(:)

      !open (1, file='test.dat')
      !do ii = 1, nz
      !    write (1,'(5f17.8)') zax(ii), p(1,ii), p(ne+1,ii), p(10,ii), p(ne+10,ii)
      !end do
      !close (1)
      !stop

      !x1 = xi(p(1:2*ne, :), 1)
      !x2 = xi(p(2*ne + 1:4*ne, :), 1)
      x1 = xi(p, 1)
      x2 = xi(p, 2)

      x1r = dreal(x1)
      x1i = dimag(x1)
      x2r = dreal(x2)
      x2i = dimag(x2)

      f1 = f(1)
      phi1 = f(2)
      f2 = f(3)
      phi2 = f(4)
      f3 = f(5)
      phi3 = f(6)

      q31 = q(3)/q(1)
      i1 = icu(1)
      r1 = r(1)
      th1 = th(1)
      dcir1 = dcir(1)
      cos1 = dcos(f(2))
      sin1 = dsin(f(2))

      q32 = q(3)/q(2)
      i2 = icu(2)
      r2 = r(2)
      th2 = th(2)
      dcir2 = dcir(2)
      cos2 = dcos(f(4))
      sin2 = dsin(f(4))

      q3 = q(3)
      a1 = a(1)
      a2 = a(2)

      s(1) = (-nharm*f1 + i1*(-x1i*cos1 + x1r*sin1) + 2*r1*nharm*f3*dcos(phi3 - phi1 - th1))*q31
      s(2) = -2*dcir1*q3 + (i1/f1*(x1r*cos1 + x1i*sin1) + 2*r1*nharm*(f3/f1)*dsin(phi3 - phi1 - th1))*q31

      s(3) = (-nharm*f2 + i2*(-x2i*cos2 + x2r*sin2) + 2*r2*nharm*f3*dcos(phi3 - phi2 - th2))*q32
      s(4) = -2*dcir2*q3 + (i2/f2*(x2r*cos2 + x2i*sin2) + 2*r2*nharm*(f3/f2)*dsin(phi3 - phi2 - th2))*q32

      s(5) = -f3 + a1*f1*dcos(phi1 - phi3) + a2*f2*dcos(phi2 - phi3)
      s(6) = a1*f1/f3*dsin(phi1 - phi3) + a2*f2/f3*dsin(phi2 - phi3)

   end subroutine dfdt

   subroutine freq()

      implicit none

      integer(c_int) :: ii, iparf, itp!, aiparp(1)
      real(c_double) t, z, yp(neqp), xoutp, &
         x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
         x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
         f1, f2, f3, phi1, phi2, phi3, a1, a2, rparf!, artolp(1), aatolp(1), arparp(1)
      complex(c_double_complex) x1, x2
      common/internp/xoutp, itp

      do ii = 2, nt
         fp(1) = f(1, ii)*cdexp(ic*f(2, ii))
         fp(2) = f(3, ii)*cdexp(ic*f(4, ii))

         call init_dopr_p()
         ioutp = 2
         z = zax(1)
         xoutp = z
         itp = 0
         yp = p(:, 1)
         iworkp(5) = neqp

         !rparp = 0.0
         !iparp = 0
         !itolp = 0
         !!rtolp = 1.0d-7
         !rtolp = ptol
         !atolp = rtolp
         !ioutp = neqp
         !z = zax(1)
         !xoutp = z
         !itp = 0
         !yp = p(:, 1)
         !iworkp(:) = 0
         !workp(:) = 0.0d0
         !iworkp(5) = neqp
         !if (nz > 0) then
         !   workp(6) = dz
         !end if

         !artolp(1) = rtolp
         !aatolp(1) = atolp
         !arparp(1) = rparp
         !aiparp(1) = iparp

         call dopri5_p(neqp, dpdz, z, yp, zex_w, artolp, aatolp, itolp, soloutp, ioutp, &
                       workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

         p(:, nz) = yp(:)

         !open (1, file='test.dat')
         !do ii = 1, nz
         !    write (1,'(5f17.8)') zax(ii), p(1,ii), p(ne+1,ii), p(10,ii), p(ne+10,ii)
         !end do
         !close (1)
         !stop

         !x1 = xi(p(1:2*ne, :), 1)
         !x2 = xi(p(2*ne + 1:4*ne, :), 1)
         x1 = xi(p, 1)
         x2 = xi(p, 2)

         x1r = dreal(x1)
         x1i = dimag(x1)
         x2r = dreal(x2)
         x2i = dimag(x2)

         f1 = f(1, ii)
         phi1 = f(2, ii)
         f2 = f(3, ii)
         phi2 = f(4, ii)
         f3 = f(5, ii)
         phi3 = f(6, ii)

         q31 = q(3)/q(1)
         i1 = icu(1)
         r1 = r(1)
         th1 = th(1)
         dcir1 = dcir(1)
         cos1 = dcos(phi1)
         sin1 = dsin(phi1)

         q32 = q(3)/q(2)
         i2 = icu(2)
         r2 = r(2)
         th2 = th(2)
         dcir2 = dcir(2)
         cos2 = dcos(phi2)
         sin2 = dsin(phi2)

         q3 = q(3)
         a1 = a(1)
         a2 = a(2)

         !s(1) = (-nharm*f1 + i1*(-x1i*cos1 + x1r*sin1) + 2*r1*nharm*f3*dcos(phi3 - phi1 - th1))*q31
         w(1, ii) = -2*dcir1*q3 + (i1/f1*(x1r*cos1 + x1i*sin1) + 2*r1*nharm*(f3/f1)*dsin(phi3 - phi1 - th1))*q31

         !s(3) = (-nharm*f2 + i2*(-x2i*cos2 + x2r*sin2) + 2*r2*nharm*f3*dcos(phi3 - phi2 - th2))*q32
         w(2, ii) = -2*dcir2*q3 + (i2/f2*(x2r*cos2 + x2i*sin2) + 2*r2*nharm*(f3/f2)*dsin(phi3 - phi2 - th2))*q32

         !s(5) = -f3 + a1*f1*dcos(phi1 - phi3) + a2*f2*dcos(phi2 - phi3)
         w(3, ii) = a1*f1/f3*dsin(phi1 - phi3) + a2*f2/f3*dsin(phi2 - phi3)
      end do
   end subroutine freq

   subroutine calc_u(u, zex_w, nz, zax)
      import
      implicit none

      integer(c_int), intent(in) :: nz
      real(c_double), intent(in) :: zex_w, zax(nz)
      complex(c_double_complex), intent(out) :: u(:)

      integer(c_int) i

      if (fok .eq. .false.) then
         do i = 1, nz
            u(i) = dexp(-3*((zax(i) - zex_w/2)/(zex_w/2))**2)
         end do
      else
         do i = 1, nz
            u(i) = squval(zax(i))
         end do
      end if

      !open(1, file = 'test.dat')
      !do i = 1,nz
      !   write(1, '(3f12.6)') zax(i), dreal(u(i)), dimag(u(i))
      !end do
      !close(1)
      !stop

   end subroutine

   subroutine soloutf(nr, xold, x, y, n, con, icomp, nd, rparf, iparf, irtrn)
      implicit none

      interface
         function contd5_f(ii, x, con, icomp, nd)
            implicit double precision(a - h, o - z)
            dimension con(5*nd), icomp(nd)
         end
      end interface

      integer(c_int) nr, n, nd, icomp(nd), iparf, irtrn, j, itf
      real(c_double) xold, x, con(5*nd), rparf, y(neqf), xoutf, pex(neqp), yy(neqf)
      logical(4) pressed
      character(1) key
      integer(c_int), parameter :: esc = 27
      common/internf/xoutf, itf

      if (nr .eq. 1) then
         itf = 1
         do j = 1, neqf
            f(j, itf) = y(j)
         end do
         call calcpex(y, pex, cl1(itf), lhs1(itf), rhs1(itf), cl2(itf), lhs2(itf), rhs2(itf))
         eta(:, itf) = eff(pex)
         !eta(:, itf) = eff(p(:, nz))
         etag(:, itf) = pitch**2/(pitch**2 + 1)*eta(:, itf)
         !write (*, '(a,f10.5,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f6.3,a,f6.3,\,a,a)') 'Time = ', xoutf, &
            write (*, '(a,f8.3,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f9.5,a,f9.5,a,f9.5,a,f5.3,a,f5.3,a,\,a)') 't =', xoutf, &
            '  |F1| = ', abs(f(1, itf)), '  |F2| = ', abs(f(3, itf)), &
            '  |F3| = ', abs(f(5, itf)), '  Eff1 = ', eta(1, itf), '  Eff2 = ', eta(2, itf), &
            '  ph1 = ', f(2,itf), '  ph2 = ', f(4,itf), '  ph3 = ', f(6,itf), &
            '  c1 = ', abs(cl1(itf)/rhs1(itf))*100, ' %  c2 = ', abs(cl2(itf)/rhs2(itf))*100, ' %', char(13)
         xoutf = x + dt
      else
10       continue
         if (x .ge. xoutf) then
            itf = itf + 1
            do j = 1, neqf
               f(j, itf) = contd5_f(j, xoutf, con, icomp, nd)
            end do
            call calcpex(y, pex, cl1(itf), lhs1(itf), rhs1(itf), cl2(itf), lhs2(itf), rhs2(itf))
            eta(:, itf) = eff(pex)
            !eta(:, itf) = eff(p(:, nz))
            etag(:, itf) = pitch**2/(pitch**2 + 1)*eta(:, itf)
            !write (*, '(a,f10.5,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f6.3,a,f6.3,\,a,a)') 'Time = ', xoutf, &
            write (*, '(a,f8.3,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f9.5,a,f9.5,a,f9.5,a,f5.3,a,f5.3,a,\,a)') 't =', xoutf, &
               '  |F1| = ', abs(f(1, itf)), '  |F2| = ', abs(f(3, itf)), &
               '  |F3| = ', abs(f(5, itf)), '  Eff1 = ', eta(1, itf), '  Eff2 = ', eta(2, itf), &
               '  ph1 = ', f(2,itf), '  ph2 = ', f(4,itf), '  ph3 = ', f(6,itf), &
               '  c1 = ', dabs(cl1(itf)/rhs1(itf))*100, ' %  c2 = ', dabs(cl2(itf)/rhs2(itf))*100, ' %', char(13)
            xoutf = xoutf + dt
            goto 10
         end if
      end if

      pressed = peekcharqq()
      if (pressed) then
         key = getcharqq()
         if (ichar(key) .eq. esc) then
            write (*, '(/,a)') 'Quit?'
            key = getcharqq()
            if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
               nt = itf
               irtrn = -1
               !return
            end if
         end if
      end if
      return
   end subroutine soloutf

   subroutine soloutp(nr, xold, x, y, n, con, icomp, nd, rparp, iparp, irtrn)
      implicit none

      interface
         function contd5_p(ii, x, con, icomp, nd)
            implicit double precision(a - h, o - z)
            dimension con(5*nd), icomp(nd)
         end
      end interface

      integer(c_int) nr, n, nd, icomp(nd), iparp, irtrn, j, itp
      real(c_double) xold, x, con(5*nd), rparp, y(neqp), xoutp
      logical(4) pressed
      character(1) key
      integer(c_int), parameter :: esc = 27
      common/internp/xoutp, itp

      if (nr .eq. 1) then
         itp = 1
         do j = 1, neqp
            p(j, itp) = y(j)
         end do
         xoutp = x + dz
      else
10       continue
         if (x .ge. xoutp) then
            itp = itp + 1
            do j = 1, neqp
               p(j, itp) = contd5_p(j, xoutp, con, icomp, nd)
            end do
            xoutp = xoutp + dz
            goto 10
         end if
      end if
      return
   end subroutine soloutp

   subroutine calcpex(f, yp, c1, lhs1, rhs1, c2, lhs2, rhs2)
      implicit none

      real(8), intent(inout) :: c1, lhs1, rhs1, c2, lhs2, rhs2

      real(c_double), intent(in) :: f(neqf)
      real(c_double) z, yp(neqp), xoutp, p2ex_mean(2), p20_mean(2)!, artolp(1), aatolp(1), arparp(1)
      integer(c_int) i, itp!, aiparp(1)
      complex(8) ptmp(ne)
      common/internp/xoutp, itp

      fp(1) = f(1)*exp(ic*f(2))
      fp(2) = f(3)*exp(ic*f(4))
      
      call init_dopr_p()
      ioutp = 0
      z = zax(1)
      xoutp = z
      itp = 0
      yp = p(:, 1)            
      iworkp(5) = 0

      !rparp = 0.0
      !iparp = 0
      !itolp = 0
      !!rtolp = 1.0d-7
      !rtolp = ptol
      !atolp = rtolp
      !ioutp = 0
      !z = zax(1)
      !xoutp = z
      !itp = 0
      !yp = p(:, 1)
      !iworkp(:) = 0
      !workp(:) = 0.0d0
      !!iworkp(5) = neqp
      !if (nz > 0) then
      !   workp(6) = dz
      !end if

      !artolp(1) = rtolp
      !aatolp(1) = atolp
      !arparp(1) = rparp
      !aiparp(1) = iparp

      call dopri5_p(neqp, dpdz, z, yp, zex_w, artolp, aatolp, itolp, solout_fiction, 0, &
                    workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

      do i = 1, 2

         ptmp(:) = dcmplx(p(idxre(i, :), 1), p(idxim(i, :), 1))
         p20_mean(i) = sum(cdabs(ptmp(:)*cdabs(ptmp(:))))/ne

         ptmp(:) = dcmplx(yp(idxre(i, :)), yp(idxim(i, :)))
         p2ex_mean(i) = sum(cdabs(ptmp(:)*cdabs(ptmp(:))))/ne
      end do

      lhs1 = 4*f(1)**2 - 8*r(1)*f(5)*f(1)*dcos(th(1) - f(6) + f(2))
      rhs1 = -icu(1)*(p2ex_mean(1) - p20_mean(1))
      c1 = lhs1 - rhs1

      lhs2 = 4*f(3)**2 - 8*r(2)*f(5)*f(3)*dcos(th(2) - f(6) + f(4))
      rhs2 = -icu(2)*(p2ex_mean(2) - p20_mean(2))
      c2 = lhs2 - rhs2

   end subroutine calcpex

   subroutine solout_fiction
   end subroutine solout_fiction

end module fun
