C The thin film equation (fourth-order PDE).
C 
C-----------------------------------------------------------------------
      SUBROUTINE F(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE 
C       NPDE DIMENSIONAL PARABOLIC PARTIAL DIFFERENTIAL EQUATION 
C                        UT = F(T, X, U, UX, UXX). 
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
        DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE 
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        FVAL(NPDE)
C                               FVAL(1:NPDE) IS THE RIGHT HAND SIDE
C                               VECTOR F(T, X, U, UX, UXX) OF THE PDE.
C-----------------------------------------------------------------------
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL

C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      eps=1.D-7
      FVAL(1) = -3.0*U(1)*U(1)*UX(1)*(UX(1)+UX(2))
     *          - U(1)**3*(UXX(2)+U(2))-Eb*U(2)
      FVAL(2) = (UXX(1) - U(2))/eps
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DERIVF(T, X, U, UX, UXX, DFDU, DFDUX, DFDUXX, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE IS USED TO DEFINE THE INFORMATION ABOUT THE 
C       PDE REQUIRED TO FORM THE ANALYTIC JACOBIAN MATRIX FOR THE DAE
C       OR ODE SYSTEM. ASSUMING THE PDE IS OF THE FORM
C                        UT = F(T, X, U, UX, UXX)
C       THIS ROUTINE RETURNS THE JACOBIANS D(F)/D(U), D(F)/D(UX), AND
C       D(F)/D(UXX).
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
        DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE 
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        DFDU(NPDE,NPDE)
C                               DFDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DFDUX(NPDE,NPDE)
C                               DFDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DFDUXX(NPDE,NPDE)
C                               DFDUXX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SECOND SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C-----------------------------------------------------------------------
C       LOOP INDICES
        INTEGER I, J
C-----------------------------------------------------------------------
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL

C-----------------------------------------------------------------------
C
C     ASSIGN DFDU(1:NPDE,1:NPDE), DFDUX(1:NPDE,1:NPDE), AND
C     DFDUXX(1:NPDE,1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      eps = 1D-7
      DFDU(1,1) = -6.0*U(1)*UX(1)*(UX(1)+UX(2)) 
     +  -3.0*U(1)*U(1)*(UXX(2)+U(2))
      DFDU(1,2) = -U(1)**3-Eb
      DFDU(2,1) = 0.0
      DFDU(2,2) = -1.0/eps
      DFDUX(1,1) = -6.0*U(1)*U(1)*UX(1)
      DFDUX(1,2) = -3.0*U(1)*U(1)*UX(1)
      DFDUX(2,1) = 0.0
      DFDUX(2,2) = 0.0
      DFDUXX(1,1) = 0.0
      DFDUXX(1,2) = -U(1)**3
      DFDUXX(2,1) = 1.0/eps
      DFDUXX(2,2) = 0.0
      
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DIFBXA(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY 
C       CONDITIONS AT THE LEFT SPATIAL END POINT X = XA. FOR THE 
C       BOUNDARY CONDITION EQUATION
C                              B(T, U, UX) = 0
C       THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED 
C       BY THIS ROUTINE.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        DBDU(NPDE,NPDE)
C                               DBDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C                               DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDT(NPDE)
C                               DBDT(I) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO TIME T.
C 
C-----------------------------------------------------------------------
C       LOOP INDICES
        INTEGER I, J
C-----------------------------------------------------------------------
c       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL

C-----------------------------------------------------------------------
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF 
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
        

      DO I=1,2
        DO J=1,2
          DBDU(I,J) = 0.d0
          DBDUX(I,J) = 0.d0
        END DO
        DBDT(I) = 0.d0
        DBDUX(I,I) = 1.d0
      END DO

C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DIFBXB(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY 
C       CONDITIONS AT THE RIGHT SPATIAL END POINT 1 = XB. FOR THE 
C       BOUNDARY CONDITION EQUATION
C                              B(T, U, UX) = 0
C       THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED 
C       BY THIS ROUTINE.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        DBDU(NPDE,NPDE)
C                               DBDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C                               DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDT(NPDE)
C                               DBDT(I) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO TIME T.
C
C-----------------------------------------------------------------------
C       LOOP INDICES
        INTEGER I, J
C-----------------------------------------------------------------------
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL

C-----------------------------------------------------------------------
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF 
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).

      DO I=1,2
        DO J=1,2
          DBDU(I,J) = 0.d0
          DBDUX(I,J) = 0.d0
        END DO
        DBDT(I) = 0.d0
        DBDUX(I,I) = 1.d0
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE BNDXA(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       LEFT SPATIAL END POINT X = XA.
C                           B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,XA).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,XA).
C
C OUTPUT:
        DOUBLE PRECISION        BVAL(NPDE)
C                               BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C                               AT THE LEFT BOUNDARY POINT.
C-----------------------------------------------------------------------
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL

C-----------------------------------------------------------------------
      
      BVAL(1) = UX(1)
      BVAL(2) = UX(2)
c      BVAL(3) = UX(3)
c      BVAL(4) = UX(4)
      
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE BNDXB(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       RIGHT SPATIAL END POINT X = XB.
C                           B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,XB).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,XB).
C
C OUTPUT:
        DOUBLE PRECISION        BVAL(NPDE)
C                               BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C                               AT THE RIGHT BOUNDARY POINT.
C-----------------------------------------------------------------------
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2,Rl

C-----------------------------------------------------------------------

      BVAL(1)= UX(1)
      BVAL(2)= UX(2)
c      BVAL(3)= UX(3)
c      BVAL(4)= UX(4)

C
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE TRUU(T, X, U, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C     THIS FUNCTION PROVIDES THE EXACT SOLUTION OF THE PDE.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
C OUTPUT:
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE EXACT SOLUTION AT THE 
C                               POINT (T,X).
C-----------------------------------------------------------------------
C Ture solution is unknown
      U(1) = 0.d0
      U(2) = 0.d0
c      U(3) = 0.d0
c      U(4) = 0.d0
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE UINIT(X, U, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE IS USED TO RETURN THE NPDE-VECTOR OF INITIAL 
C       CONDITIONS OF THE UNKNOWN FUNCTION AT THE INITIAL TIME T = T0 
C       AT THE SPATIAL COORDINATE X.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        DOUBLE PRECISION        X
C                               THE SPATIAL COORDINATE.
C
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
C OUTPUT:
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS VECTOR OF INITIAL VALUES OF
C                               THE UNKNOWN FUNCTION AT T = T0 AND THE
C                               GIVEN VALUE OF X.
C-----------------------------------------------------------------------
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2,RL

C-----------------------------------------------------------------------
     
      pi = 4.0*ATAN(1.0)
c RL = 10.0
      w = 2.0*pi/RL
      a = 1.0d0
      a1 = 0.5*a
      U(1) = a+a1*dcos(w*X)
      U(2) = -a1*w*w*dcos(w*X) 

C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine header(nout)
c-----------------------------------------------------------------------
c purpose:
c       this subroutine writes a header describing the npde dimensional
c       parabolic partial differential equation
c                        ut = f(t, x, u, ux, uxx).
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 nout
c                               nout is the output unit number.
c-----------------------------------------------------------------------
c constants:
        double precision        t0
        parameter              (t0 = 0.0d0)
c
        double precision        xa
        parameter              (xa = 0.0d0)
c
        double precision        xb
        parameter              (xb = 10.0d0)
c-----------------------------------------------------------------------
c
      write(nout,95) 'The Thin Film system'
      write(nout,95) 'domain:'
      write(nout,96) '   t0 =', t0, ' < t,'
      write(nout,96) '   xa =', xa, ' <= x <= xb =', xb, ','
c
      return
   95 format(a)
   96 format(a,e13.5,a,e13.5,a,e13.5,a,e13.5,a)
      end


