C Copyright (C) 2002, 2010 Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C
C =============================================================================
C
C     This is an example for the usage of IPOPT in double precision.
C     It implements problem 71 from the Hock-Schittkowski test suite:
C
C     min   x1*x4*(x1 + x2 + x3)  +  x3
C     s.t.  x1*x2*x3*x4                   >=  25
C           x1**2 + x2**2 + x3**2 + x4**2  =  40
C           1 <=  x1,x2,x3,x4  <= 5
C
C     Starting point:
C        x = (1, 5, 5, 1)
C
C     Optimal solution:
C        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
C
C =============================================================================
C
C
C =============================================================================
C
C                            Main driver program
C
C =============================================================================
C
      program example
C
      implicit none
C
C     include the Ipopt return codes
C
      include 'IpReturnCodes.inc'
C
C     Size of the problem (number of variables and equality constraints)
C
      integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
      parameter  (N = 4, M = 2, NELE_JAC = 8, NELE_HESS = 10)
      parameter  (IDX_STY = 1 )
C
C     Space for multipliers and constraints
C
      double precision LAM(M)
      double precision G(M)
C
C     Vector of variables
C
      double precision X(N)
C
C     Vector of lower and upper bounds, left and right hand sides, and scaling
C
      double precision X_L(N), X_U(N), X_SCALING(N), Z_L(N), Z_U(N)
      double precision G_L(M), G_U(M), G_SCALING(M)

C     To indicate to IPSETPROBLEMSCALING below that scaling for constraints
C     should be omitted, a NULL pointer can be passed for G_SCALING. With
C     Fortran > 77, this is done by removing the declaration of G_SCALING
C     from the line above and the data G_SCALING line below and enabling
C     the following line instead:
C
C     double precision, allocatable, dimension(:) :: G_SCALING
C
C
C     Private data for evaluation routines
C     This could be used to pass double precision and integer arrays untouched
C     to the evaluation subroutines EVAL_*
C
      double precision DAT(2)
      integer IDAT(1)
C
C     Place for storing the Ipopt Problem Handle
C
CC     for 32 bit platforms
C      integer IPROBLEM
C      integer IPCREATE
C     for 64 bit platforms:
      integer*8 IPROBLEM
      integer*8 IPCREATE
C
      integer IERR
      integer IPSOLVE, IPSETPROBLEMSCALING
      integer IPADDSTROPTION, IPADDNUMOPTION, IPADDINTOPTION
      integer IPOPENOUTPUTFILE
C
      double precision F
      integer i
C
C     The following are the Fortran routines for computing the model
C     functions and their derivatives - their code can be found further
C     down in this file.
C
      external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS
CC
CC     The next is an optional callback method.  It is called once per
CC     iteration.
CC
      external ITER_CB
C
C     Set initial point, bounds, and scaling:
C
      data X   / 1d0, 5d0, 5d0, 1d0 /
      data X_L / 1d0, 1d0, 1d0, 1d0 /
      data X_U / 5d0, 5d0, 5d0, 5d0 /
      data X_SCALING / 1d0, 1d0, 1d0, 1d0 /
C
C     Set bounds and scaling for the constraints
C
      data G_L / 25d0, 40d0 /
      data G_U / 1d40, 40d0 /
      data G_SCALING / 1d0, 1d0 /
C
C     First create a handle for the Ipopt problem (and read the options
C     file)
C
      IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,
     1     IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
      if (IPROBLEM.eq.0) then
         write(*,*) 'Error creating an Ipopt Problem handle.'
         stop
      endif
C
C     Open an output file
C
      IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
      if (IERR.ne.0 ) then
         write(*,*) 'Error opening the Ipopt output file.'
         goto 9000
      endif
C
C     Note: The following options are only examples, they might not be
C           suitable for your optimization problem.
C
C     Set a string option
C
      IERR = IPADDSTROPTION(IPROBLEM, 'mu_strategy', 'adaptive')
      if (IERR.ne.0 ) goto 9990
C
C     Set an integer option
C
      IERR = IPADDINTOPTION(IPROBLEM, 'max_iter', 3000)
      if (IERR.ne.0 ) goto 9990
C
C     Set a double precision option
C
      IERR = IPADDNUMOPTION(IPROBLEM, 'tol', 1.d-7)
      if (IERR.ne.0 ) goto 9990

C
C     Set scaling
C
      IERR = IPSETPROBLEMSCALING(IPROBLEM, 1d0, X_SCALING, G_SCALING)
      if (IERR.ne.0 ) goto 9990

CC
CC     Set a callback function to give you control once per iteration.
CC     You can use it if you want to generate some output, or to stop
CC     the optimization early.
CC     Uncomment this line to print iterate in every iteration and
CC     demonstrate how to stop Ipopt early.
CC
C      call IPSETCALLBACK(IPROBLEM, ITER_CB)

C
C     As a simple example, we pass the constants in the constraints to
C     the EVAL_C routine via the "private" DAT array.
C
      DAT(1) = 0.d0
      DAT(2) = 0.d0
C
C     Call optimization routine
C     Passing IPROBLEM instead of IDAT is a hack to make IPROBLEM
C     accessible within ITER_CB. A proficient Fortran programmer would
C     find a way to store IPROBLEM inside IDAT.
C
      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IPROBLEM, DAT)
C
C     Output:
C
      if( IERR.eq.IP_SOLVE_SUCCEEDED ) then
         write(*,*)
         write(*,*) 'The solution was found.'
         write(*,*)
         write(*,*) 'The final value of the objective function is ',F
         write(*,*)
         write(*,*) 'The optimal values of X are:'
         write(*,*)
         do i = 1, N
            write(*,*) 'X  (',i,') = ',X(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the lower bounds are:'
         write(*,*)
         do i = 1, N
            write(*,*) 'Z_L(',i,') = ',Z_L(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the upper bounds are:'
         write(*,*)
         do i = 1, N
            write(*,*) 'Z_U(',i,') = ',Z_U(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the equality constraints are:'
         write(*,*)
         do i = 1, M
            write(*,*) 'LAM(',i,') = ',LAM(i)
         enddo
         write(*,*)
      else
         write(*,*)
         write(*,*) 'An error occoured.'
         write(*,*) 'The error code is ',IERR
         write(*,*)
      endif
C
 9000 continue
C
C     Clean up
C
      call IPFREE(IPROBLEM)
      stop
C
 9990 continue
      write(*,*) 'Error setting an option'
      goto 9000
      end
C
C =============================================================================
C
C                    Computation of objective function
C
C =============================================================================
C
      subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X
      double precision F, X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      F = X(1)*X(4)*(X(1)+X(2)+X(3)) + X(3)
      IERR = 0
      return
      end
C
C =============================================================================
C
C                Computation of gradient of objective function
C
C =============================================================================
C
      subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X
      double precision GRAD(N), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      GRAD(1) = X(4)*(2d0*X(1)+X(2)+X(3))
      GRAD(2) = X(1)*X(4)
      GRAD(3) = X(1)*X(4) + 1d0
      GRAD(4) = X(1)*(X(1)+X(2)+X(3))
      IERR = 0
      return
      end
C
C =============================================================================
C
C                     Computation of equality constraints
C
C =============================================================================
C
      subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X, M
      double precision G(M), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      G(1) = X(1)*X(2)*X(3)*X(4) - DAT(1)
      G(2) = X(1)**2 + X(2)**2 + X(3)**2 + X(4)**2 - DAT(2)
      IERR = 0
      return
      end
C
C =============================================================================
C
C                Computation of Jacobian of equality constraints
C
C =============================================================================
C
      subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,
     1     IDAT, DAT, IERR)
      integer TASK, N, NEW_X, M, NZ
      double precision X(N), A(NZ)
      integer ACON(NZ), AVAR(NZ), I
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C     structure of Jacobian:
C
      integer AVAR1(8), ACON1(8)
      data  AVAR1 /1, 2, 3, 4, 1, 2, 3, 4/
      data  ACON1 /1, 1, 1, 1, 2, 2, 2, 2/
      save  AVAR1, ACON1
C
      if( TASK.eq.0 ) then
        do I = 1, 8
          AVAR(I) = AVAR1(I)
          ACON(I) = ACON1(I)
        enddo
      else
        A(1) = X(2)*X(3)*X(4)
        A(2) = X(1)*X(3)*X(4)
        A(3) = X(1)*X(2)*X(4)
        A(4) = X(1)*X(2)*X(3)
        A(5) = 2d0*X(1)
        A(6) = 2d0*X(2)
        A(7) = 2d0*X(3)
        A(8) = 2d0*X(4)
      endif
      IERR = 0
      return
      end
C
C =============================================================================
C
C                Computation of Hessian of Lagrangian
C
C =============================================================================
C
      subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,
     1     NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
      implicit none
      integer TASK, N, NEW_X, M, NEW_LAM, NNZH, i
      double precision X(N), OBJFACT, LAM(M), HESS(NNZH)
      integer IRNH(NNZH), ICNH(NNZH)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C     structure of Hessian:
C
      integer IRNH1(10), ICNH1(10)
      data  IRNH1 /1, 2, 2, 3, 3, 3, 4, 4, 4, 4/
      data  ICNH1 /1, 1, 2, 1, 2, 3, 1, 2, 3, 4/
      save  IRNH1, ICNH1

      if( TASK.eq.0 ) then
         do i = 1, 10
            IRNH(i) = IRNH1(i)
            ICNH(i) = ICNH1(i)
         enddo
      else
         do i = 1, 10
            HESS(i) = 0d0
         enddo
C
C     objective function
C
         HESS(1) = OBJFACT * 2d0*X(4)
         HESS(2) = OBJFACT * X(4)
         HESS(4) = OBJFACT * X(4)
         HESS(7) = OBJFACT * (2d0*X(1) + X(2) + X(3))
         HESS(8) = OBJFACT * X(1)
         HESS(9) = OBJFACT * X(1)
C
C     first constraint
C
         HESS(2) = HESS(2) + LAM(1) * X(3)*X(4)
         HESS(4) = HESS(4) + LAM(1) * X(2)*X(4)
         HESS(5) = HESS(5) + LAM(1) * X(1)*X(4)
         HESS(7) = HESS(7) + LAM(1) * X(2)*X(3)
         HESS(8) = HESS(8) + LAM(1) * X(1)*X(3)
         HESS(9) = HESS(9) + LAM(1) * X(1)*X(2)
C
C     second constraint
C
         HESS(1) = HESS(1) + LAM(2) * 2d0
         HESS(3) = HESS(3) + LAM(2) * 2d0
         HESS(6) = HESS(6) + LAM(2) * 2d0
         HESS(10)= HESS(10)+ LAM(2) * 2d0
      endif
      IERR = 0
      return
      end
C
C =============================================================================
C
C                   Callback method called once per iteration
C
C =============================================================================
C
      subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU,
     1     MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,
     2     DAT, ISTOP)
      implicit none
      integer ALG_MODE, ITER_COUNT, LS_TRIAL
      double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
      double precision ALPHA_DU, ALPHA_PR
      double precision DAT(*)
      integer IDAT(*)
      integer ISTOP
      integer IPGETCURRITERATE
      integer IPGETCURRVIOLATIONS
      integer i, IERR
      double precision X(4)
      double precision Z_L(4)
      double precision Z_U(4)
      double precision G(2)
      double precision LAMBDA(2)
      double precision COMPL_X_L(4)
      double precision COMPL_X_U(4)
      double precision GRAD_LAG_X(4)
      double precision CONS_VIOL(2)
      double precision COMPL_G(2)
C
C     Get iterate and print
C
      IERR = IPGETCURRITERATE(IDAT, 0, 1, 1, 1, 1, 4, X, Z_L, Z_U, 2, G,
     1   LAMBDA)
      if ( IERR.eq.0 ) then
         do i = 1, 4
            write(*,*) 'X, Z_L, Z_U (',i,') = ', X(i), Z_L(i), Z_U(i)
         enddo
         do i = 1, 2
            write(*,*) 'G, LAMBDA   (',i,') = ', G(i), LAMBDA(i)
         enddo
      endif

C
C     Get violations and print
C
      IERR = IPGETCURRVIOLATIONS(IDAT, 0, 0, 1, 1, 1, 4,
     1   0, 0, COMPL_X_L, COMPL_X_U, GRAD_LAG_X, 2, CONS_VIOL, COMPL_G)
      if ( IERR.eq.0 ) then
         do i = 1, 4
            write(*,*) 'COMPL_X_L, COMPL_X_U, GRAD_LAG_X (',i,') = ',
     1         COMPL_X_L(i), COMPL_X_U(i), GRAD_LAG_X(i)
         enddo
         do i = 1, 2
            write(*,*) 'CONS_VIOL, COMPL_G (',i,') = ',
     1         CONS_VIOL(i), COMPL_G(i)
         enddo
      endif

C
C     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
C     simple example.
C
      if (INF_PR.le.1D-04) then
         write(*,*) 'Requesting stop due to small inf_pr in iteration ',
     1      ITER_COUNT
         ISTOP = 1
      endif

      return
      end
