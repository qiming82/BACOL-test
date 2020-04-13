c Copyright (c) 2013, Paul Muir, Jack Pew, Zhi Li
c Paul Muir, Mathematics and Computing Science, Saint Mary's University.
c All rights reserved.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions
c are met:
c * Redistributions of source code must retain the above copyright
c   notice, this list of conditions and the following disclaimer.
c * Redistributions in binary form must reproduce the above copyright
c   notice, this list of conditions and the following disclaimer in the
c   documentation and/or other materials provided with the distribution.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
c HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
c LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
c THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
c (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
c OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c This file contains the BACOLI source code (including the new code
c that implements the interpolation based spatial error estimates.)

c-----------------------------------------------------------------------
c This driver file writes output for each PDE on an almost uniform grid
c of dimension nxout x ntout. tstart is 0, tstop is prompted for and the
c spatial domain is [0,1]. This format can be convenient for plotting.
c A number of parameters are hard coded and should be edited as needed:
c npde, xa, xb, tstart, nintmx, kcol, nintinit, nderiv, nxout, ntout,
c and est.
c
c When prompted to input the tolerance, it is recommended to
c use engineering notation (e.g., '1e-6').
c
c After a successful computation, this program will have written several
c files named X, T, U1, ..., UNPDE, Ux1, ..., UxNPDE, and
c Uxx1, ..., UxxNPDE which contain the x and t coordinates of output,
c the solution approximations, and the first and second spatial
c derivative approximations, respectively.
c A line of a U* file contains values corresponding to a single t value.
c Values of NPDE beyond 9 break this simple driver.
c
c This driver should be linked with with bacoli.f and bacoli-aux.f
c (version 3) and a problem definition file, e.g., burg2.f.
c
c Example Python plotting code:
c
c   import matplotlib as mpl
c   mpl.use('AGG')  # for systems not running a GUI
c   from mpl_toolkits.mplot3d import Axes3D
c   from matplotlib import cm
c   import matplotlib.pyplot as plt
c   import numpy as np
c
c   styling = {
c     'cmap': cm.coolwarm,
c     'linewidth': 0,
c     'antialiased': True
c   }
c
c   ufiles  = ['U1',      'Ux1',       'Uxx1']
c   utitles = ['$u(t,x)$','$u_x(t,x)$','$u_{xx}(t,x)$']
c   udata = [np.loadtxt(f) for f in ufiles]
c   x, t = np.meshgrid(np.loadtxt('X'), np.loadtxt('T'))
c   fig = plt.figure(figsize=(24,6))  # inches
c
c   for i, u, title in zip(range(len(ufiles)), udata, utitles):
c     ax = fig.add_subplot(1, 3, i+1, projection='3d')
c     ax.plot_surface(x, t, u, rstride=1, cstride=1, **styling)
c     ax.set_xlabel('$x$')
c     ax.set_ylabel('$t$')
c     ax.set_zlabel(title)
c     ax.set_title(title)
c   plt.savefig('gridmesh.png')
c
c-----------------------------------------------------------------------
        implicit none
c constants:
        integer                 npde
        parameter              (npde = 3)
c                               number of PDEs in the system.
        double precision        xa
        parameter              (xa = 0d0)
c                               the left spatial boundary point.
        double precision        xb
        parameter              (xb = 10.d0)
c                               the right spatial boundary point.
        double precision        tstart
        parameter              (tstart = 0d0)
c                               the start of the time domain, used to
c                               initialise t0.
        integer                 nintinit
        parameter              (nintinit = 200)
c                               the number of subintervals in the
c                               initial uniform mesh.
        integer                 nintmx
        parameter              (nintmx = 1000)
c                               maximum number of subintervals that
c                               bacoli may use
        integer                 nderiv
        parameter              (nderiv = 2)
c                               number of spatial partial derivatives
c                               to output approximations for, along
c                               with the solution value approximations.
c                               note that the solution approximation is
c                               a piecewise polynomial of degree kcol-1,
c                               so nderiv is limited.
        integer                 nxout
        parameter              (nxout = 201)
c                               number of almost uniformly distributed
c                               points along the spatial domain to
c                               output approximated values. the first
c                               point is not tstart, but close to it.
        integer                 ntout
        parameter              (ntout = 41)
c                               number of uniformly distributed points
c                               along the temporal domain to output
c                               approximated values.
        integer                 est
        parameter              (est = 1)
c                               set est = 1 to use the LOI/LE scheme,
c                               or set est = 2 to use the SCI/ST scheme.
        integer                 kcol
        parameter              (kcol = 4)
c                               kcol is the number of collocation
c                               points to be used in each subinterval.
c                                  3 <= kcol <= 10
c                               A value of 4 is most efficient in
c                               typical cases, and higher when
c                               a sharp tolerance is enforced.
c                               The sci scheme can have issues finding
c                               an initial mesh at kcol = 3,
c                               and the loi scheme can be extremely
c                               inefficient with sharp tolerances at
c                               kcol = 3.
        integer                 maxvec
        parameter              (maxvec = npde*(nintmx*kcol+2))
c                               the maximum size of the vector of
c                               B-spline coefficients, y.
        integer                 lrp
        parameter              (lrp = 113+59*npde+27*nintmx+13*npde**2
     +                              +9*kcol+24*kcol*nintmx
     +                              +6*nintmx*kcol**2
     +                              +27*npde*nintmx*kcol
     +                              +7*nintmx*npde
     +                              +2*npde**2*nintmx*kcol**2
     +                              +4*npde**2*kcol*nintmx
     -                              -15*nintmx+3*kcol
     -                              -8*kcol*nintmx+kcol**2
     -                              -nintmx*kcol**2-3*nintmx*npde+2
     +                         )
c                               length of bacoli's floating point
c                               work array, rpar.
c                               the lines with a `-' continuation
c                               character may be uncommented if only
c                               the loi estimate is used.
c
        integer                 lip
        parameter              (lip = 101+npde*(kcol*nintmx+2))
c                               length of bacoli's integer
c                               work array, ipar.
        integer                 lenwrk
        parameter              (lenwrk = (kcol+2)*(nderiv+1)
     +                                 + kcol*(nintmx+1)+4)
c                               length of the work array used by VALUES
c                               after a successful return from bacoli.
        integer                 lmf
        parameter              (lmf = 12)
c                               lmf is the length of the mflag vector.
c-----------------------------------------------------------------------
c bacoli input:
        double precision        t0
c                               t0 < tout is the initial time for a
c                               given call to bacoli.
        double precision        tout(ntout)
c                               tout(:) are the desired output times.
        double precision        tstop
c                               tout is the final desired output time.
        double precision        atol(npde)
c                               atol is the absolute error tolerance
c                               request and is a scalar quantity if
c                               mflag(2) = 0.
        double precision        rtol(npde)
c                               rtol is the relative error tolerance
c                               request and is a scalar quantity if
c                               mflag(2) = 0.
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
        double precision        x(nintmx+1)
c                               x is the spatial mesh which divides the
c                               interval [xa,xb] as:
c                               xa = x(1) < x(2) < .. < x(nint+1) = xb
        integer                 mflag(lmf)
c                               this vector is used for more specific
c                               settings. see documentation in bacoli.f
        external                f
        external                bndxa
        external                bndxb
        external                uinit
        external                dum
c                               PDE system definition subroutines.
c                               see burg1.f for an example, or the
c                               bacoli.f preamble for descriptions.
c-----------------------------------------------------------------------
c bacoli work storage:
        double precision        rpar(lrp)
c                               rpar is a floating point work array
c                               of size lrp.
        integer                 ipar(lip)
c                               ipar is an integer work array
c                               of size lip.
c-----------------------------------------------------------------------
c bacoli output:
        double precision        y(maxvec), umin
c                               on successful return from bacoli, y is
c                               the vector of time dependent B-spline
c                               coefficients at the current time t0.
        integer                 idid
c                               idid is the bacoli exit status flag
c-----------------------------------------------------------------------
c miscellany:
        double precision        uout(npde,nxout,nderiv+1)
c                               contains solution and derivative
c                               approximations, which are obtained from
c                               the time dependent B-spline coefficients
c                               returned by bacoli (in y) through
c                               calling the subroutine values.
        double precision        valwrk(lenwrk)
c                               valwrk is the work array for values.
        double precision        xout(nxout)
c                               xout is a set of spatial points for
c                               which output is desired.
        integer                 i, j, k, m
c                               loop indices
        character*3             eststr
c                               eststr is printed in the input summary.
        character               ks
        character*32            fname
c                               fname is where output filenames are made
c-----------------------------------------------------------------------
c problem specific variables:
        double precision        eps
        common /burger/         eps
c                               eps is used in the burg1.f, burg2.f and
c                               Cahn_Allen.f problem definition files.
c
        double precision        a1,a2,d1,d2,r,c,n,pe1,pe2,RL,del,epsln
        common /rcdstuff/       a1,a2,d1,d2,r,c,n,pe1,pe2,RL,del,epsln
c                               used in RCDsys.f.
c-----------------------------------------------------------------------
c subroutines called:
c                               bacoli
c                               values
c-----------------------------------------------------------------------
c  Common block initializations

      eps = 1d-3
      del = 0.1d0
      epsln = 1d-7
      RL = 10.d0
      a1 = 1.d0  ! ram
      a2 = 30.d0
      d1 = 1.5d0
      d2 = 1.2d0
      r = 0.00001d0
      c = 1.d0/r
      n = 1.d0
      pe1 = 10000.d0
      pe2 = 10000.d0

c-----------------------------------------------------------------------
c  Get some user input

c      write(*,*) 'Enter tstop, the end of the temporal domain.'
c      read(*,*,err=999) tstop

c      write(*,*) 'Please input the error tolerance.'
c      read(*,*,err=999) atol(1)
c      rtol(1) = atol(1)

      tstop = 40.0
      atol(1) = 0.001d0
      rtol(1) = atol(1)
c-----------------------------------------------------------------------
c  BACOLI parameter initializations

      nint = nintinit
      t0 = tstart

c     Initialize rpar and ipar with zeros.
      do i = 1, lrp
         rpar(i) = 0d0
      end do
      do i = 1, lip
         ipar(i) = 0
      end do

c  define the mesh based on a uniform step size.
      x(1) = xa
      do i = 2, nint
         x(i) = xa + (i-1) * (xb-xa) / nint
      end do
      x(nint+1) = xb

c  define output spacial point sequence, uniformly spaced
c  ----- ###############   ---  CAN I USE THE ADAPTIVE MESH from Calculations??  -----   ############

      xout(1) = xa
      do i = 2, nxout-1
         xout(i) = xa + (i-1) * (xb-xa) / (nxout-1)
      end do
      xout(nxout) = xb

c  define output temporal point sequence, not quite uniformly
c  distributed. tout(1) is only 'close to' tstart so that bacoli
c  may be used to approximate spatial derivatives at all output.
c  this has the downside of skipping interesting behaviour at tstart.
      do i = 2, ntout-1
         tout(i) = tstart + (i-1) * (tstop-tstart) / (ntout-1)
      end do
      tout(1) = (tout(2)-tstart) / 1000
      tout(ntout) = tstop

c  initialize the mflag vector.
      do i = 1, lmf
         mflag(i) = 0
      end do
c  set mflag(5) = 1 for Dirichlet boundary conditions, 0 otherwise
c  recall that Dirichlet boundary conditions involve only function
c  values, i.e., U_x is not involved in bndxa or bndxb.
      mflag(5) = 1

      if (est .eq. 1) then
         mflag(8) = 0
         eststr = 'LOI'
      elseif (est .eq. 2) then
         mflag(8) = 1
         eststr = 'SCI'
      else
         stop
      endif

c-----------------------------------------------------------------------
c  Computation & output section

      write(*,*)     'THE INPUT IS'
      write(*,9000) 'kcol = ', kcol, ', nint0 = ', nint,
     &              ', npde = ', npde, ', tout = ', tstop
      write(*,9001) 'atol = ', atol(1), ', rtol = ', rtol(1),
     &              ',                 ', eststr

c     open files to receive output.
      do k = 1, npde
c        String manipulation sure is awkward in F77.
c        F90/95 intrinsics and/or format specifiers are recommended.
         write(ks,'(i1)') k
         fname = 'U'
         do m = 1, nderiv+1
            if (m .gt. 1) fname(m:m) = 'x'
            fname(m+1:m+1) = ks
            open(unit=m*10+k,file=fname)
         end do
      end do

c     Write the tout and xout sequences for constructing a meshgrid
c     when plotting.
      open(unit=10,file='X')
      do i = 1, nxout
         write(10,*) xout(i)
      end do
      close(10)

      open(unit=10,file='T')
      do i = 1, nxout
         write(10,*) tout(i)
      end do
      close(10)

      Open (1,file='Uval.dat')
      Open (2,file='hmin.dat')

c     Call BACOLI at each tout and save output to files.
      do j = 1, ntout

c        idid has to be set to 1 for continuation calls.
c        On the initial call, it does not matter if it is 1
         idid = 1

         call bacoli(t0, tout(j), atol, rtol, npde, kcol, nintmx,
     &               nint, x, mflag, rpar, lrp, ipar, lip, y, idid,
     &               f, dum, bndxa, dum, bndxb, dum, uinit)
c        Check for errors running BACOLI.
         if (idid .le. 0) goto 700

c        Write approximated solution and derivative values to file.
         call values(kcol, xout, nint, x, npde, nxout, nderiv, uout, y,
     &               valwrk)
c         do k = 1, npde
c            do m = 1, nderiv+1
c               write(m*10+k,*) (uout(k, i, m), i = 1, nxout)
c            end do
c         end do
	 Do i=1,nxout
	   write(1,*) xout(i), uout(1,i,1)
       Enddo
	write(2,*) tout(j), umin
	umin = minval(Uout(1,1:nxout,1))
	if(umin.le.0.0001d0) then
	  stop
	endif

c        mflag(1) must be 0 on the initial call, and 1 on continuation.
         mflag(1) = 1
      end do

c-----------------------------------------------------------------------
c  Sucessful return section
      write(*,9002) 'IDID       = ', idid
      write(*,9002) 'Final nint = ', nint
      write(*,9002) 'nsteps     = ', ipar(10)

  600 continue
c     close files
      do k = 1, npde
         do m = 1, nderiv+1
            close(m*10+k)
         end do
      end do

      stop
c-----------------------------------------------------------------------
c  Unsuccessful return section
  700 continue
      write(*,*) 'CANNOT PROCEED DUE TO ERROR FROM BACOLI.'
      write(*,9003) 'Failure at T = ', t0
      write(*,9002) 'IDID = ', idid
      write(*,9002) 'NINT = ', nint

      goto 600
c-----------------------------------------------------------------------
c  Invalid input error section
  999 continue
      write(*,*) 'Error: Improperly formatted input'

      stop
c-----------------------------------------------------------------------
c  Formats!
 9000 format( a, i2, a, i4, a, i3, a, e7.1 )
 9001 format( 2(a, e7.1), 2a )
 9002 format( a, i10 )
 9003 format( a, f15.8 )
c-----------------------------------------------------------------------
      end

      subroutine dum
c     This is a dummy subroutine to provide bacoli in place of
c     derivf, difbxa and difbxb when using approximated derivatives.
      end
