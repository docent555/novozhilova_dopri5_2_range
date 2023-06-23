program sys15f
   use, intrinsic :: iso_c_binding
   use fun
   use ifport

   implicit none

   integer(c_int) i, j, hours, minutes, seconds, res
   real(c_double) start_time, stop_time, calc_time
   complex(c_double_complex) pc

   call init()

   do i = 1, ne
      p(i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2.0d0*pi))
      p(ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2.0d0*pi))
      p(2*ne + i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2.0d0*pi))
      p(3*ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2.0d0*pi))
   end do

   !open (1, file='p_init.dat')
   !do i = 1, ne
   !   write (1, '(i, 4f17.8)') i, p(i, 1), p(ne + i, 1), p(2*ne + i, 1), p(3*ne + i, 1)
   !end do
   !close (1)

   write (*, '(/)')

   start_time = dclock()
   do i = 1, ndtr
      dtr(1) = dtr0 + (i - 1)*dtrh_w
      dtr(2) = dtr0 + (i - 1)*dtrh_w

      write (path, '(f5.3,a)') dtr(1), '/'
      !print *, 'DTR = ', dtr(1)
      
      res = makedirqq(path)
      
      call write_param(path)
      call init_local(path)

      call ode4f()
      !stop_time = dclock()

      !calc_time = stop_time - start_time
      !
      !hours = calc_time/3600
      !minutes = (calc_time - hours*3600)/60
      !seconds = calc_time - hours*3600 - minutes*60

      !write (*, '(/)')
      !print *, 'Calcualting took:', hours, 'h :', minutes, 'm :', seconds, 's'

      !start_time = dclock()
      write (*, '(/)')
      print *, 'Writing...'
      !write (*, '(/)')

      if (inher .eq. .true.) then
         f10 = f(1, nt)
         f20 = f(3, nt)
         f30 = f(5, nt)
         p10 = mod(f(2, nt), 2*pi)
         p20 = mod(f(4, nt), 2*pi)
         p30 = mod(f(6, nt), 2*pi)
      end if      

      call write_results()
      
      f(1, 1) = f10
      f(2, 1) = p10
      f(3, 1) = f20
      f(4, 1) = p20
      f(5, 1) = f30
      f(6, 1) = p30
   end do
   stop_time = dclock()

   calc_time = stop_time - start_time

   hours = calc_time/3600
   minutes = (calc_time - hours*3600)/60
   seconds = calc_time - hours*3600 - minutes*60

   print *, 'Execution time:', hours, 'h :', minutes, 'm :', seconds, 's'
   write (*, '(/)')

   pause
   stop

!   !do i = 2, nt
!   !   do j = 1, 3
!   !      w(j, i - 1) = (f(2*j, i) - f(2*j, i - 1))/dt
!   !   end do
!   !end do
!
!   call freq()
!
!   phi(:, 1) = 0;
!   do i = 2, nt
!      do j = 1, 3
!         phi(j, i) = phi(j, i - 1) + dimag(log(f(2*j - 1, i)*cdexp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*cdexp(ic*f(2*j, i - 1)))))
!      end do
!   end do
!
!   breaknum(:) = 0
!   fcomp(1) = f(2*1 - 1, 1)*cdexp(ic*f(2*1, 1))
!   fcomp(2) = f(2*2 - 1, 1)*cdexp(ic*f(2*2, 1))
!   fcomp(3) = f(2*3 - 1, 1)*cdexp(ic*f(2*3, 1))
!   phitmp0(:) = datan2(dimag(fcomp(:)), dreal(fcomp(:)))
!   !phitmp0(:) = datan2(dimag(f(:, 1)), dreal(f(:, 1)))
!   phios(:, 1) = phitmp0(:)
!   do i = 2, nt
!      do j = 1, 3
!         fc = f(2*j - 1, i)*cdexp(ic*f(2*j, i))
!         phitmp1(j) = datan2(dimag(fc), dreal(fc))
!         if ((phitmp1(j) - phitmp0(j)) .gt. pi) breaknum(j) = breaknum(j) - 1
!         if ((phitmp1(j) - phitmp0(j)) .lt. -pi) breaknum(j) = breaknum(j) + 1
!         phios(j, i) = phitmp1(j) + 2.*pi*breaknum(j)
!         !phios(j, i) = phitmp1(j)
!         phitmp0(j) = phitmp1(j)
!      end do
!   end do
!
!   do i = 1, nt - 1
!      do j = 1, 3
!         wos(j, i) = (phios(j, i + 1) - phios(j, i))/dt
!      end do
!   end do
!
!   write (*, '(/)')
!
!   pause
!
!   open (3, file='cl1.dat')
!   do i = 1, nt
!      write (3, '(4e17.8)') tax(i), cl1(i), lhs1(i), rhs1(i)
!   end do
!   close (3)
!
!   open (3, file='cl2.dat')
!   do i = 1, nt
!      write (3, '(4e17.8)') tax(i), cl2(i), lhs2(i), rhs2(i)
!   end do
!   close (3)
!
!   open (1, file='F.dat')
!   do i = 1, nt
!      write (1, '(4e17.8)') tax(i), f(1, i), f(3, i), f(5, i)
!   end do
!   close (1)
!
!   open (13, file='FCMPLX.dat')
!   do i = 1, nt
!      fcomp(1) = f(2*1 - 1, i)*cdexp(ic*f(2*1, i))
!      fcomp(2) = f(2*2 - 1, i)*cdexp(ic*f(2*2, i))
!      fcomp(3) = f(2*3 - 1, i)*cdexp(ic*f(2*3, i))
!      write (13, '(7e17.8)') tax(i), dreal(fcomp(1)), dimag(fcomp(1)), dreal(fcomp(2)), dimag(fcomp(2)), &
!         dreal(fcomp(3)), dimag(fcomp(3))
!   end do
!   close (13)
!
!   open (2, file='E.dat')
!   do i = 1, nt
!      write (2, '(5e17.8)') tax(i), eta(1, i), etag(1, i), eta(2, i), etag(2, i)
!   end do
!   close (2)
!
!   open (3, file='W.dat')
!   do i = 1, nt - 1
!      write (3, '(4e17.8)') tax(i + 1), w(1, i), w(2, i), w(3, i)
!   end do
!   close (3)
!
!   open (1, file='P.dat')
!   do i = 1, nt
!      write (1, '(4e17.8)') tax(i), f(2, i), f(4, i), f(6, i)
!   end do
!   close (1)
!
!   open (1, file='POS.dat')
!   do i = 1, nt
!      write (1, '(4e17.8)') tax(i), phios(1, i), phios(2, i), phios(3, i)
!   end do
!   close (1)
!
!   open (3, file='WOS.dat')
!   do i = 1, nt - 1
!      write (3, '(4e17.8)') tax(i + 1), wos(1, i), wos(2, i), wos(3, i)
!   end do
!   close (3)
!
!   call write_param()
!
!   stop
!101 print *, 'error of file open.'
!   pause
!   stop
!102 print *, 'error of file reading.'
!   pause
!   stop
!103 print *, 'error of file writing.'
!   pause
!   stop
end program
