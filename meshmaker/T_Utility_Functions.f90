!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>      T_UtilityFunctions.f95: Code unit including the routines       >
!>                describing general utility functions                 >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      MODULE Utility_Functions
! 
         PRIVATE
! 
         PUBLIC :: Polynomial,PowerSeries,ExponSeries,SineSeries,CosSeries,LogSeries,  &
      &            Integral_poly,Integral_PowerSrs,Integral_ExponSrs,                  &
      &            Integral_SineSrs,Integral_CosSrs,Integral_LogSrs,                   &
      &            Error_Functions,Incomplete_Beta,LogBeta,Gamma,LogGamma,             &
      &            SORT, Table_Interpolation, Table_Interval, N_Character
!
!
! 
         CONTAINS
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION Polynomial(argument,n_poly,A)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                 Routine for computing a polynomial                  *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_poly
! 
! -------------
! ......... Integer variables
! -------------
! 
            INTEGER :: i
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                      :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_poly) :: A
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of Q_Polynomial
!
!
            Polynomial = A(n_poly)
!
            DO i = n_poly-1,0,-1                                    
!
               Polynomial = Polynomial*argument + A(i)
!
            END DO 
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of Polynomial
!
!
            RETURN
!
         END FUNCTION Polynomial
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION Integral_poly(argument,n_poly,A)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         Routine for computing the integral of a polynomial          *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_poly
! 
! -------------
! ......... Integer variables
! -------------
! 
            INTEGER :: i
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                      :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_poly) :: A
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Integral_poly
!
!
            Integral_poly = A(n_poly)/(DBLE(n_poly) + 1.0d0)
!
            DO i = n_poly-1,0,-1                                    
!
               Integral_poly = Integral_poly*argument + A(i)/(DBLE(i) + 1.0d0)
!
            END DO
!
            Integral_poly = Integral_poly*argument
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Integral_poly
!
!
            RETURN
!
         END FUNCTION Integral_poly
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION PowerSeries(argument,n_power,A,B,acceptable)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                Routine for computing a power series                 *
!*                                                                     *
!*                     Version 1.0 - June 6, 2005                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_power
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                       :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_power) :: A,B
! 
! -------------
! ......... Logical variable
! -------------
! 
            LOGICAL, INTENT(OUT) :: acceptable
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of PowerSeries
!
!
            IF(argument == 0.0d0 .AND. (ANY(B < 0.0d0))) THEN
               acceptable = .FALSE.               
            ELSE
               acceptable = .TRUE.              
               PowerSeries = SUM(A*(argument**B))
            END IF
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of PowerSeries
!
!
            RETURN
!
         END FUNCTION PowerSeries
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION Integral_PowerSrs(argument,n_power,A,B,acceptable)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        Routine for computing the integral of a power series         *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_power
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                       :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_power) :: A,B
! 
! -------------
! ......... Logical variable
! -------------
! 
            LOGICAL, INTENT(OUT) :: acceptable
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Integral_PowerSrs
!
!
            IF(argument == 0.0d0 .AND. (ANY(B < 0.0d0))) THEN
               acceptable = .FALSE.               
            ELSE
               acceptable = .TRUE. 
!
               Integral_PowerSrs = SUM((A/(B + 1.0d0))*(argument**(B + 1.0d0)), MASK = B .NE. -1.0d0)
               IF(ANY(B == -1.0d0)) Integral_PowerSrs = Integral_PowerSrs + SUM((A*LOG(argument)), MASK = B .EQ. -1.0d0) 
            END IF
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Integral_PowerSrs
!
!
            RETURN
!
         END FUNCTION Integral_PowerSrs
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION ExponSeries(argument,n_expon,A,B,C,acceptable)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*             Routine for computing an exponential series             *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_expon
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                       :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_expon) :: A,B,C
! 
! -------------
! ......... Logical variable
! -------------
! 
            LOGICAL, INTENT(OUT) :: acceptable
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of ExponSeries
!
!
            IF(ANY((B + C*argument) > HUGE(1.0d0))) THEN
               acceptable = .FALSE.   
            ELSE
               acceptable = .TRUE.   
               ExponSeries = SUM(A*(EXP(B + C*argument)))
            END IF
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of ExponSeries
!
!
            RETURN
!
         END FUNCTION ExponSeries
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION Integral_ExponSrs(argument,n_expon,A,B,C,acceptable)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     Routine for computing the integral of an exponential series     *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_expon
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                       :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_expon) :: A,B,C
! 
! -------------
! ......... Logical variable
! -------------
! 
            LOGICAL, INTENT(OUT) :: acceptable
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Integral_ExponSrs
!
!
            IF(ANY((B + C*argument) > HUGE(1.0d0))) THEN
               acceptable = .FALSE.   
            ELSE
               acceptable = .TRUE.   
               Integral_ExponSrs = SUM((A/C)*(EXP(B + C*argument)))
            END IF
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Integral_ExponSrs
!
!
            RETURN
!
         END FUNCTION Integral_ExponSrs
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION SineSeries(argument,n_sine,A,B,C)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                Routine for computing a sine series                  *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_sine
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                      :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_sine) :: A,B,C
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of SineSeries
!
!
            SineSeries = SUM(A*(SIN(B + C*argument)))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of SineSeries
!
!
            RETURN
!
         END FUNCTION SineSeries
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION Integral_SineSrs(argument,n_sine,A,B,C)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         Routine for computing the integral of a sine series         *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_sine
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                      :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_sine) :: A,B,C
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of Integral_SineSrs
!
!
            Integral_SineSrs = SUM((-A/C)*(COS(B + C*argument)))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of Integral_SineSrs
!
!
            RETURN
!
         END FUNCTION Integral_SineSrs
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION CosSeries(argument,n_cos,A,B,C)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*               Routine for computing a cosine series                 *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_cos
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                     :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_cos) :: A,B,C
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of CosSeries
!
!
            CosSeries = SUM(A*(COS(B + C*argument)))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of CosSeries
!
!
            RETURN
!
         END FUNCTION CosSeries
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION Integral_CosSrs(argument,n_cos,A,B,C)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        Routine for computing the integral of a cosine series        *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_cos
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                     :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_cos) :: A,B,C
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Integral_CosSrs
!
!
            Integral_CosSrs = SUM((A/B)*(SIN(B + C*argument)))
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Integral_CosSrs
!
!
            RETURN
!
         END FUNCTION Integral_CosSrs
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION LogSeries(base,argument,n_log,A,B,C,acceptable)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                 Routine for computing a log series                  *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_log,base
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                     :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_log) :: A,B,C
! 
! -------------
! ......... Logical variable
! -------------
! 
            LOGICAL, INTENT(OUT) :: acceptable
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of LogSeries
!
!
            IF_OK: IF(ANY((B + C*argument) == 0.0d0)) THEN
               acceptable = .FALSE.   
            ELSE
!
               acceptable = .TRUE.   
!
               IF(base == 10) THEN
                  LogSeries = SUM(A*(LOG10(B + C*argument)))
               ELSE
                  LogSeries = SUM(A*(LOG(B + C*argument)))
               END IF
!
            END IF IF_OK
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of LogSeries
!
!
            RETURN
!
         END FUNCTION LogSeries
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         REAL(KIND = 8) FUNCTION Integral_LogSrs(base,argument,n_log,A,B,C,acceptable)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                 Routine for computing a log series                  *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(IN) :: n_log,base
! 
! -------------
! ......... Real input variables
! -------------
! 
            REAL(KIND = 8), INTENT(IN)                     :: argument
            REAL(KIND = 8), INTENT(IN), DIMENSION(0:n_log) :: A,B,C
! 
! -------------
! ......... Logical variable
! -------------
! 
            LOGICAL, INTENT(OUT) :: acceptable
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Integral_LogSrs
!
!
            IF_OK: IF(ANY((B + C*argument) == 0.0d0)) THEN
               acceptable = .FALSE.   
            ELSE
!
               acceptable = .TRUE.   
!
               IF(base == 10) THEN
                  Integral_LogSrs = SUM( (A/C)*(B + C*argument)*(LOG10(B + C*argument)) - A*argument ) 
               ELSE
                  Integral_LogSrs = SUM( (A/C)*(B + C*argument)*(LOG(B + C*argument)) - A*argument )
               END IF
!
            END IF IF_OK
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of Integral_LogSrs
!
!
            RETURN
!
         END FUNCTION Integral_LogSrs
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
      REAL(KIND = 8) FUNCTION Error_Functions(ARG,JINT)
!
      IMPLICIT NONE
!     
!
!***********************************************************************
!***********************************************************************
!
! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!   for a real argument  x.  It contains the FUNCTION type
!   subprogram ERFUNC.  The calling
!   statements for the primary entries are:
!
!                   Y=ERFUNC(X,0),
!
!                   Y=ERFUNC(X,1),
!   and
!                   Y=ERFUNC(X,2).
!
!   The routine  ERFUNC  is intended for internal packet use only,
!   all computations within the packet being concentrated in this
!   routine.  The parameter usage is as follows
!
!      Function                     Parameters for ERFUNC
!       call              ARG                  ERFUNC          JINT
!
!     ERFUNC(X,0)   ANY REAL ARGUMENT         ERF(ARG)          0
!     ERFUNC(X,1)   ABS(ARG) .LT. XBIG        ERFC(ARG)         1
!     ERFUNC(X,2)   XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
!
!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmeti!
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XMIN   = the smallest positive floating-point number.
!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!            the negative of the solution to the equation
!            2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!            A conservative value is the largest machine number X
!            such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!            W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!            machine precision.  A conservative value is
!            1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!            of XINF and 1/[sqrt(pi)*XMIN].
!
!   Approximate values for some important machines are:
!
!                          XMIN       XINF        XNEG     XSMALL
!
!  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
!
!
!                          XBIG       XHUGE       XMAX
!
!  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ERFC = 0      for  ARG .GE. XBIG;
!
!                       ERFCX = XINF  for  ARG .LT. XNEG;
!      and
!                       ERFCX = 0     for  ARG .GE. XMAX.
!
!
! Intrinsic functions required are:
!
!     DABS, AINT, DEXP
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: May 16, 2004 (by G. J. Moridis, LBNL)
!
!
!***********************************************************************
!***********************************************************************
! 
! 
! -------
! ... Double precision parameter arrays
! -------
! 
! ... Coefficients for approximation to ERF in first interval
! 
      REAL(KIND = 8), PARAMETER, DIMENSION(5) ::    &
     &      A = (/ 3.16112374387056560D00, 1.13864154151050156D02,   &  
     &             3.77485237685302021D02, 3.20937758913846947D03,   &
     &             1.85777706184603153D-1 /)
! 
      REAL(KIND = 8), PARAMETER, DIMENSION(4) ::    &
     &      B = (/ 2.36012909523441209D01, 2.44024637934444173D02,   &
     &             1.28261652607737228D03, 2.84423683343917062D03 /)
! 
! ... Coefficients for approximation to ERFC in second interval
! 
      REAL(KIND = 8), PARAMETER, DIMENSION(9) ::    &
     &      C = (/ 5.64188496988670089D-1, 8.88314979438837594D0,    &
     &             6.61191906371416295D01, 2.98635138197400131D02,   &
     &             8.81952221241769090D02, 1.71204761263407058D03,   &
     &             2.05107837782607147D03, 1.23033935479799725D03,   &
     &             2.15311535474403846D-8 /)
! 
      REAL(KIND = 8), PARAMETER, DIMENSION(8) ::    &
     &      D = (/ 1.57449261107098347D01, 1.17693950891312499D02,   &
     &             5.37181101862009858D02, 1.62138957456669019D03,   &
     &             3.29079923573345963D03, 4.36261909014324716D03,   &
     &             3.43936767414372164D03, 1.23033935480374942D03 /)
! 
! ... Coefficients for approximation to ERFC in third interval
! 
      REAL(KIND = 8), PARAMETER, DIMENSION(6) ::    &
     &      P = (/ 3.05326634961232344D-1, 3.60344899949804439D-1,   &
     &             1.25781726111229246D-1, 1.60837851487422766D-2,   &
     &             6.58749161529837803D-4, 1.63153871373020978D-2 /)
! 
      REAL(KIND = 8), PARAMETER, DIMENSION(5) ::    &
     &      Q = (/ 2.56852019228982242D00, 1.87295284992346047D00,   &
     &             5.27905102951428412D-1, 6.05183413124413191D-2,   &
     &             2.33520497626869185D-3 /)
! 
! -------
! ... Double precision parameters
! -------
! 
      REAL(KIND = 8), PARAMETER :: ZERO = 0.0d0, HALF   = 5.0d-1,   &
     &                             ONE  = 1.0d0, TWO    = 2.0d0,    &
     &                             FOUR = 4.0d0, SIXTEN = 1.6d1
! 
      REAL(KIND = 8), PARAMETER :: THRESH = 4.6875D-1,   &
     &                             SQRPI  = 5.6418958354775628695D-1
! 
      REAL(KIND = 8), PARAMETER :: XINF   =  1.79D308,   &
     &                             XNEG   = -2.6628D1,   &
     &                             XSMALL =  1.11D-16,   &
     &                             XBIG   =  2.6543D1,   &
     &                             XHUGE  =  6.71D7,     & 
     &                             XMAX   =  2.53D307
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8), INTENT(IN) :: ARG
! 
      REAL(KIND = 8) :: DEL,X,Y,XNUM,XDEN,YSQ,ERFUNC
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER, INTENT(IN) :: JINT
! 
      INTEGER :: i
! 
! 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Error_Functions
! 
! 
      X = ARG
      Y = DABS(X)
! 
! **********************************************************************
! *                                                                    *
! *                Evaluate  ERF  for  |X| <= 0.46875                  *
! *                                                                    *
! **********************************************************************
!
      IF_ArgValue: IF (Y <= THRESH) THEN
!
            YSQ = ZERO
            IF (Y > XSMALL) YSQ = Y * Y
!
            XNUM = A(5)*YSQ
            XDEN = YSQ
!
            DO_Int1: DO i = 1, 3
               XNUM = (XNUM + A(i)) * YSQ
               XDEN = (XDEN + B(i)) * YSQ
            END DO DO_Int1
!
            ERFUNC = X * (XNUM + A(4)) / (XDEN + B(4))
!
            IF (JINT == 1) THEN
               ERFUNC = ONE - ERFUNC
            ELSE IF (JINT == 2) THEN
               ERFUNC = DEXP(YSQ) * ERFUNC
            END IF
! 
            Error_Functions = ERFUNC   ! Assign value to the function
!
            RETURN                     ! Done - Exit the function            
! 
! **********************************************************************
! *                                                                    *
! *             Evaluate  ERFC  for 0.46875 <= |X| <= 4.0              *
! *                                                                    *
! **********************************************************************
!
         ELSE IF (Y <= FOUR .AND. Y > THRESH) THEN
!
            XNUM = C(9)*Y
            XDEN = Y
!
            DO_Int2: DO i = 1, 7
               XNUM = (XNUM + C(i)) * Y
               XDEN = (XDEN + D(i)) * Y
            END DO DO_Int2
!
            ERFUNC = (XNUM + C(8)) / (XDEN + D(8))
!
            IF (JINT /= 2) THEN
               YSQ    = AINT(Y*SIXTEN)/SIXTEN
               DEL    = (Y-YSQ)*(Y+YSQ)
               ERFUNC = DEXP(-YSQ*YSQ) * DEXP(-DEL) * ERFUNC
            END IF
! 
! **********************************************************************
! *                                                                    *
! *                   Evaluate  ERFC  for |X| > 4.0                    *
! *                                                                    *
! **********************************************************************
!
         ELSE
! 
            ERFUNC = ZERO
! 
            IF (Y >= XBIG) THEN
! 
               IF ((JINT /= 2) .OR. (Y >= XMAX)) GO TO 1000
! 
               IF (Y >= XHUGE) THEN
                  ERFUNC = SQRPI / Y
                  GO TO 1000
               END IF
! 
            END IF
! 
            YSQ = ONE / (Y * Y)
            XNUM = P(6)*YSQ
            XDEN = YSQ
! 
            DO_Int3: DO i = 1, 4
               XNUM = (XNUM + P(I)) * YSQ
               XDEN = (XDEN + Q(I)) * YSQ
            END DO DO_Int3
! 
            ERFUNC = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
            ERFUNC = (SQRPI -  ERFUNC) / Y
! 
            IF (JINT /= 2) THEN
               YSQ    = AINT(Y*SIXTEN)/SIXTEN
               DEL    = (Y-YSQ)*(Y+YSQ)
               ERFUNC = DEXP(-YSQ*YSQ) * DEXP(-DEL) * ERFUNC
            END IF
! 
      END IF IF_ArgValue
! 
! **********************************************************************
! *                                                                    *
! *                      Finalize adjustments                          *
! *                                                                    *
! **********************************************************************
!
 1000 IF_JintValue: IF (JINT == 0) THEN         ! Compute erf(x) 
! 
            ERFUNC = (HALF - ERFUNC) + HALF
            IF (X < ZERO) ERFUNC = -ERFUNC
! 
         ELSE IF (JINT == 1) THEN               ! Compute erfc(x)
! 
            IF (X < ZERO) ERFUNC = TWO - ERFUNC
! 
         ELSE                                   ! Compute exp(x*x)*erfc(x)
! 
            IF (X < ZERO) THEN
! 
               IF (X < XNEG) THEN
                     ERFUNC = XINF
                  ELSE
                     YSQ    = AINT(X*SIXTEN)/SIXTEN
                     DEL    = (X-YSQ)*(X+YSQ)
                     Y      = DEXP(YSQ*YSQ) * DEXP(DEL)
                     ERFUNC = (Y+Y) - ERFUNC
               END IF
! 
            END IF
! 
      END IF IF_JintValue
! 
      Error_Functions = ERFUNC
! 
! 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Error_Functions
! 
! 
      RETURN
!
      END FUNCTION Error_Functions
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
      REAL(KIND = 8) FUNCTION Gamma(X)
!
!
!
      IMPLICIT NONE
! 
! -------------
! ... Double precision parameters
! -------------
! 
      REAL(KIND=8), PARAMETER :: PI     = 3.1415926535897932384626434d0,  &
     &                           SQRTPI = 0.9189385332046727417803297d0
! 
! -------------
! ... Numerator and denominator coefficients for
! ... rational minimax approximation in the interval (1,2) 
! -------------
! 
      REAL(KIND=8), DIMENSION(8), PARAMETER :: P =                          &
     &     (/ -1.71618513886549492533811D+0,  2.47656508055759199108314D+1, &
     &        -3.79804256470945635097577D+2,  6.29331155312818442661052D+2, &
     &         8.66966202790413211295064D+2, -3.14512729688483675254357D+4, &
     &        -3.61444134186911729807069D+4,  6.64561438202405440627855D+4 /)
      REAL(KIND=8), DIMENSION(8), PARAMETER :: Q =                          &
     &     (/ -3.08402300119738975254353D+1,  3.15350626979604161529144D+2, &
     &        -1.01515636749021914166146D+3, -3.10777167157231109440444D+3, &
     &         2.25381184209801510330112D+4,  4.75584627752788110767815D+3, &
     &        -1.34659959864969306392456D+5, -1.15132259675553483497211D+5 /)
! 
! -------------
! ... Coefficients for minimax approximation in the interval (12, INF)
! -------------
! 
      REAL(KIND=8), DIMENSION(7), PARAMETER :: C =                          &
     &     (/ -1.910444077728D-03,          8.4171387781295D-04,            &
     &        -5.952379913043012D-04,       7.93650793500350248D-04,        &
     &        -2.777777777777681622553D-03, 8.333333333333333331554247D-02, &
     &         5.7083835261D-03 /)
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: CONV,FACT,SUM,X,XDEN,XNUM,Y,Y1,YSQ,Z,F_Gamma
! 
! -------
! ... Integer ariables
! -------
! 
      INTEGER :: I,N
! 
! -------
! ... Logical variables
! -------
! 
      LOGICAL :: PARITY
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Gamma
!
!
         CONV(I) = DBLE(I)     ! .....  Initializations
         PARITY  = .FALSE.
         FACT    = 1.0d0
! 
         N = 0
         Y = X
! 
!***********************************************************************      
!*                                                                     *
!*                       For negative arguments                        *
!*                                                                     *
!***********************************************************************      
!
         IF_Negative: IF(Y <= 0.0d0) THEN
!
            Y       = -X
            Y1      =  INT(Y)
            F_gamma =  Y - Y1
!
            IF(F_gamma /= 0.0d0) THEN
               IF(Y1 /= INT(Y1*5.0d-1)*2.0d0) PARITY = .TRUE.
               FACT = -PI/SIN(PI*F_gamma)
               Y    =  Y + 1.0d0
            ELSE
               Gamma = HUGE(1.0d0)
               RETURN
            END IF
!
         END IF IF_Negative
! 
!***********************************************************************      
!*                                                                     *
!*                       For positive arguments                        *
!*                                                                     *
!***********************************************************************      
!
         IF_Positive: IF(Y < epsilon(1.0d0)) THEN
! 
!-----------
! ......... For a very small argument
!-----------
! 
            IF(Y >= TINY(1.0d0)) THEN
               F_Gamma = 1.0d0 / Y
            ELSE
               Gamma = HUGE(1.0d0)
               RETURN
            END IF
! 
!-----------
! ......... For an argument Y < 12
!-----------
! 
         ELSE IF (Y < 1.2d1) THEN
! 
            Y1 = Y
! 
! ......... For an argument Y <= 1
! 
            IF(Y < 1.0d0) THEN
               Z = Y
               Y = Y + 1.0d0
            ELSE
! 
! ........... For an argument Y > 1 and Y < 12
! 
              N = INT(Y) - 1
              Y = Y - CONV(N)
              Z = Y - 1.0d0
! 
            END IF
! 
!-----------
! ......... Evaluate approximation for 1.0 < argument < 2.0
!-----------
! 
            XNUM = 0.0d0
            XDEN = 1.0d0
            DO I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
            END DO
            F_Gamma = XNUM / XDEN + 1.0d0
! 
!-----------
! ......... Adjust result for case  0.0 < argument < 1.0
!-----------
! 
            IF(Y1 < Y) THEN
! 
               F_Gamma = F_Gamma / Y1
! 
!-----------
! ......... Adjust result for case  2.0 < argument < 12.0
!-----------
! 
            ELSE IF (Y1 > Y) THEN
               DO I = 1, N
                  F_Gamma = F_Gamma * Y
                  Y       = Y + 1.0d0
               END DO
            END IF
! 
!-----------
! ...... Evaluate for argument >= 12.0 
!-----------
! 
         ELSE
! 
            IF(Y <= 171.624D0) THEN
! 
               YSQ = Y * Y
               SUM = C(7)
               DO I = 1, 6
                  SUM = SUM / YSQ + C(I)
               END DO
! 
               SUM     = SUM/Y - Y + SQRTPI
               SUM     = SUM + (Y-5.0d-1)*DLOG(Y)
               F_Gamma = DEXP(SUM)
! 
            ELSE
               Gamma = HUGE(1.0d0)
               RETURN
            END IF
! 
! 
! 
         END IF IF_Positive
! 
!-----------
! ...... Final adjustments and return
!-----------
! 
         IF (PARITY)        F_Gamma = -F_Gamma
         IF (FACT /= 1.0d0) F_Gamma = FACT / F_Gamma
! 
! 
! 
 1000    Gamma = F_Gamma
!
!
!
         RETURN
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Gamma
!
!
      END FUNCTION Gamma
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!
!
!
      REAL(KIND = 8) FUNCTION D9LGMC (X)
!
!
!
      IMPLICIT NONE
! 
! -------
! ... Double precision parameters
! -------
! 
      REAL(KIND = 8), PARAMETER, DIMENSION(15) :: ALGMCS = (/   &
                         +.1666389480451863247205729650822D+0,  & 
                         -.1384948176067563840732986059135D-4,  & 
                         +.9810825646924729426157171547487D-8,  & 
                         -.1809129475572494194263306266719D-10, & 
                         +.6221098041892605227126015543416D-13, & 
                         -.3399615005417721944303330599666D-15, & 
                         +.2683181998482698748957538846666D-17, & 
                         -.2868042435334643284144622399999D-19, & 
                         +.3962837061046434803679306666666D-21, & 
                         -.6831888753985766870111999999999D-23, & 
                         +.1429227355942498147573333333333D-24, & 
                         -.3547598158101070547199999999999D-26, & 
                         +.1025680058010470912000000000000D-27, & 
                         -.3401102254316748799999999999999D-29, & 
                         +.1276642195630062933333333333333D-30  /)
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: X, XBIG, XMAX
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER :: NALGM
! 
! -------
! ... Logical variables
! -------
! 
      LOGICAL :: FIRST = .TRUE.
! 
! -------
! ... Saving variables
! -------
! 
      SAVE NALGM, XBIG, XMAX, FIRST
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of D9LGMC
!
!
         IF (FIRST) THEN
            NALGM = INITDS (ALGMCS, 15, (EPSILON(1.0d0)/2.0d0) )
            XBIG  = 1.0D0/SQRT(EPSILON(1.0d0)/2.0d0)
            XMAX  = EXP (MIN(LOG(HUGE(1.0d0)/12.D0), -LOG(12.D0*TINY(1.0d0))))
            FIRST = .FALSE.
         ENDIF
!
         IF (X >= XMAX) THEN 
            D9LGMC = 0.0d0
         ELSE 
            D9LGMC = 1.D0/(12.D0*X)
            IF (X < XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS, NALGM) / X
         END IF
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of D9LGMC
!
!
      RETURN
!
      END FUNCTION D9LGMC
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!
!
!
      REAL(KIND = 8) FUNCTION Incomplete_Beta(PIN,QIN,x)
!
!
!
      IMPLICIT NONE
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: X,PIN,QIN,ALNEPS,ALNSML,C,EPS,FINSUM,P,PS,Q,SML,TERM,XB,XI,Y,P1
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER :: n,i,IB
! 
! -------
! ... Logical variables
! -------
! 
      LOGICAL :: FIRST = .TRUE.
! 
! -------
! ... Saving variables
! -------
! 
      SAVE EPS, ALNEPS, SML, ALNSML, FIRST
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of Incomplete_Beta
!
!
      IF (FIRST) THEN
         EPS    = EPSILON(1.0d0)/2.0d0
         ALNEPS = LOG (EPS)
         SML    = TINY(1.0d0)
         ALNSML = LOG (SML)
      ENDIF
      FIRST = .FALSE.
!
! ... Initialization
!
      Y = X
      P = PIN
      Q = QIN
!
      IF( (Q <= P .AND. X < 0.8d0) .OR. (X < 0.2d0) ) THEN
         CONTINUE
      ELSE
         Y = 1.0D0 - Y
         P = QIN
         Q = PIN
      END IF
!
      Incomplete_Beta = 0.0d0
!
! ---------
! ... This is a simple case !
! ---------
!
      IF ((P+Q)*Y/(P+1.D0) < EPS) THEN
!
         XB = P*LOG(MAX(Y,SML)) - LOG(P) - LogBeta(P,Q)
!
         IF (XB > ALNSML .AND. Y /= 0.0D0) Incomplete_Beta = EXP(XB)
         IF (Y /= X .OR. P /= PIN)         Incomplete_Beta = 1.0D0 - Incomplete_Beta
!
         RETURN    ! ... We are done !
!
      END IF
!
! ---------
! ... Evaluate the infinite sum first.  The term will equal: Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I)
! ---------
!
      PS = Q - AINT(Q)
      IF (PS == 0.D0) PS = 1.0D0
      XB = P*LOG(Y) - LogBeta(PS,P) - LOG(P)
!
!
!
      IF_XBRelSize: IF (XB >= ALNSML) THEN
!
         Incomplete_Beta = EXP (XB)
         TERM            = Incomplete_Beta*P
!
         IF_PSNotOne: IF(PS /= 1.0d0) THEN
            N = MAX (ALNEPS/LOG(Y), 4.0D0)
            DO i=1,N
              XI              = DBLE(i)
              TERM            = TERM * (XI - PS)*Y/XI
              Incomplete_Beta = Incomplete_Beta + TERM/(P + XI)
            END DO
         END IF IF_PSNotOne
!
      END IF IF_XBRelSize
!
! ---------
! ... Now attempt to evaluate the finite sum, if possible.
! ---------
!
      IF_QPos: IF (Q > 1.0D0) THEN
!
!
         XB   = P*LOG(Y) + Q*LOG(1.0d0 - Y) - LogBeta(P,Q) - LOG(Q)
         IB   = MAX (XB/ALNSML, 0.0d0)
         TERM = EXP(XB - IB*ALNSML)
!
         C    = 1.0D0/(1.0d0 - Y)
         P1   = Q*C/(P + Q - 1.0d0)
!
         FINSUM = 0.0d0
         N      = INT(Q)
         IF (Q == DBLE(N)) N = N - 1
!
         DO_Loop2: DO I=1,N
!
           IF (P1 <= 1.0d0 .AND. ((TERM/EPS) <= FINSUM)) EXIT DO_Loop2
!
           XI = DBLE(I)
           TERM = (Q - XI + 1.0d0)*C*TERM/(P + Q - XI)
!
           IF (TERM > 1.0D0) THEN
              IB = IB - 1
              TERM = TERM*SML
           END IF
!
           IF (IB == 0) FINSUM = FINSUM + TERM
!
         END DO DO_Loop2
!
         Incomplete_Beta = Incomplete_Beta + FINSUM
!
!
      END IF IF_QPos
!
!
!
      IF ( (Y /= X) .OR. (P /= PIN) ) Incomplete_Beta = 1.0d0 - Incomplete_Beta
!
      Incomplete_Beta = MAX (MIN (Incomplete_Beta, 1.0d0), 0.0d0)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of Incomplete_Beta
!
!
      RETURN
!
      END FUNCTION Incomplete_Beta
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!
!
!
      REAL(KIND = 8) FUNCTION DCSEVL(X, CS, N)
!
!
!
      IMPLICIT NONE
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER :: N,ni
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8)               :: B0,B1,B2,TWOX,X
      REAL(KIND = 8), DIMENSION(n) :: CS(n)
! 
! -------
! ... Logical variables
! -------
! 
      LOGICAL :: FIRST = .TRUE.
! 
! -------
! ... Saving variables
! -------
! 
      SAVE FIRST
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DCSEVL
!
!
         IF(FIRST) FIRST = .FALSE.
!
         B1   = 0.0d0
         B0   = 0.0d0
         TWOX = 2.0d0*X
!
         DO ni = N,1,-1
            B2 = B1
            B1 = B0
            B0 = TWOX*B1 - B2 + CS(ni)
         END DO
!
         DCSEVL = 0.5d0*(B0 - B2)
!
!
!
         RETURN
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DCSEVL
!
!
      END FUNCTION DCSEVL 
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!
!
!
      REAL(KIND = 8) FUNCTION LogBeta (A, B)
!
!
!
      IMPLICIT NONE
! 
! -------
! ... Double precision parameters
! -------
! 
      REAL(KIND=8), PARAMETER :: SQ2PIL = 0.91893853320467274178032973640562D0 
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: A,B,P,Q
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of LogBeta
!
!
      P = MIN(A, B)
      Q = MAX(A, B)
!
! ... P and Q are big 
!
      IF (P >= 10.0d0) THEN
!
         LogBeta = -0.5d0*LOG(Q) + SQ2PIL + D9LGMC(P) + D9LGMC(Q) - D9LGMC(P + Q)  &
     &                           + (P - 0.5d0)*LOG(P/(P + Q)) + Q*DLNREL(-P/(P + Q))
!
! ... P is small, but Q is big 
!
      ELSE IF (Q >= 10.0d0) THEN
!
         LogBeta = LogGamma(P) + D9LGMC(Q) - D9LGMC(P + Q) + P*(1.0d0 - LOG(P+Q)) + (Q - 0.5d0)*DLNREL(-P/(P + Q))
!
! ... P and Q are small 
!
      ELSE  
!
         LogBeta = LOG( GAMMA(P)*(GAMMA(Q)/GAMMA(P + Q)) )
!
      END IF
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of LogBeta
!
!
      RETURN
!
      END FUNCTION LogBeta
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!
!
!
      REAL(KIND = 8) FUNCTION LogGamma(X)
!
!
!
      IMPLICIT NONE
! 
! -------
! ... Double precision parameters
! -------
! 
      REAL(KIND = 8), PARAMETER ::  SQ2PIL = 9.11893853320467274178032973640562D-1
      REAL(KIND = 8), PARAMETER ::  SQPI2L = 2.25791352644727432363097614947441D-1   
      REAL(KIND = 8), PARAMETER ::  PI     = 3.14159265358979323846264338327950D0    
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: X, Y
! 
! -------
! ... Logical variables
! -------
! 
      LOGICAL :: FIRST = .TRUE.
! 
! -------
! ... Saving variables
! -------
! 
      SAVE FIRST
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of LogGamma
!
!
      IF (FIRST) FIRST = .FALSE.
!
      Y = ABS (X)
!
! -------
! ... LOG (ABS (GAMMA(X)) ) for ABS(X) <= 10.0
! -------
!
      IF ( Y <= 10.0D0) THEN
!
         LogGamma = LOG( ABS( GAMMA(X) ) )
!
! -------
! ... LOG ( ABS (GAMMA(X)) ) for ABS(X) > 10.0
! -------
!
      ELSE
!
         IF (X > 0.0d0) THEN
            LogGamma = SQ2PIL + (X - 0.5d0)*LOG(X) - X + D9LGMC(Y)
         ELSE
            LogGamma = SQPI2L + (X - 0.5d0)*LOG(Y) - X - LOG(ABS (SIN(PI*Y))) - D9LGMC(Y)
         END IF
!
      END IF
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of LogGamma
!
!
      RETURN
!
      END FUNCTION LogGamma
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!
!
!
      REAL(KIND = 8) FUNCTION DLNREL(X)
!
!
!
      IMPLICIT NONE
! 
! -------
! ... Double precision parameters
! -------
! 
      REAL(KIND = 8), PARAMETER, DIMENSION(43) :: ALNRCS = (/  &
                       +.10378693562743769800686267719098D+1,  &                  
                       -.13364301504908918098766041553133D+0,  &                 
                       +.19408249135520563357926199374750D-1,  &                 
                       -.30107551127535777690376537776592D-2,  &                 
                       +.48694614797154850090456366509137D-3,  &                 
                       -.81054881893175356066809943008622D-4,  &                 
                       +.13778847799559524782938251496059D-4,  &                 
                       -.23802210894358970251369992914935D-5,  &                 
                       +.41640416213865183476391859901989D-6,  &                 
                       -.73595828378075994984266837031998D-7,  &                 
                       +.13117611876241674949152294345011D-7,  &                 
                       -.23546709317742425136696092330175D-8,  &                 
                       +.42522773276034997775638052962567D-9,  &                 
                       -.77190894134840796826108107493300D-10, &                 
                       +.14075746481359069909215356472191D-10, &                 
                       -.25769072058024680627537078627584D-11, &                 
                       +.47342406666294421849154395005938D-12, &                 
                       -.87249012674742641745301263292675D-13, &                 
                       +.16124614902740551465739833119115D-13, &                 
                       -.29875652015665773006710792416815D-14, &                 
                       +.55480701209082887983041321697279D-15, &                 
                       -.10324619158271569595141333961932D-15, &                 
                       +.19250239203049851177878503244868D-16, &                 
                       -.35955073465265150011189707844266D-17, &                 
                       +.67264542537876857892194574226773D-18, &                 
                       -.12602624168735219252082425637546D-18, &                 
                       +.23644884408606210044916158955519D-19, &                 
                       -.44419377050807936898878389179733D-20, &                 
                       +.83546594464034259016241293994666D-21, &                 
                       -.15731559416479562574899253521066D-21, &                 
                       +.29653128740247422686154369706666D-22, &                 
                       -.55949583481815947292156013226666D-23, &                 
                       +.10566354268835681048187284138666D-23, &                 
                       -.19972483680670204548314999466666D-24, &                 
                       +.37782977818839361421049855999999D-25, &                 
                       -.71531586889081740345038165333333D-26, &                 
                       +.13552488463674213646502024533333D-26, &                 
                       -.25694673048487567430079829333333D-27, &                 
                       +.48747756066216949076459519999999D-28, &                 
                       -.92542112530849715321132373333333D-29, &                 
                       +.17578597841760239233269760000000D-29, &                 
                       -.33410026677731010351377066666666D-30, &                 
                       +.63533936180236187354180266666666D-31 /)                 
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: X
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER :: NLNREL
! 
! -------
! ... Logical variables
! -------
! 
      LOGICAL :: FIRST = .TRUE.
! 
! -------
! ... Saving variables
! -------
! 
      SAVE NLNREL, FIRST
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DLNREL
!
!
      IF (FIRST) THEN
         NLNREL = INITDS (ALNRCS, 43, 1.0d-1*DBLE(EPSILON(1.0d0)/2.0d0))
         FIRST  = .FALSE.
      END IF
!
      IF (ABS(X) <= 0.375d0) THEN
         DLNREL = X*( 1.0d0 - X*DCSEVL(X/0.375d0, ALNRCS, NLNREL) )
      ELSE
         DLNREL = LOG(1.0d0 + X)
      END IF
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of DLNREL
!
!
      RETURN
!
      END FUNCTION DLNREL
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!
!
!
      INTEGER FUNCTION INITDS(OS, nos, ETA)
!
!
!
      IMPLICIT NONE
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER :: nos,i
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: OS(nos)
      REAL(KIND = 8) :: ETA,ERR
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of INITDS
!
!
      ERR = 0.0d0
!
      DO_Loop: DO i = NOS,1,-1
        ERR = ERR + ABS(OS(i))
        IF (ERR > ETA) EXIT DO_Loop
      END DO DO_Loop
!
      INITDS = i
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of INITDS
!
!
      RETURN
!
      END FUNCTION INITDS
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!
!
!
      SUBROUTINE SORT(iar,n)
!
!***********************************************************************
!*                                                                     *
!*             SORTING THE VECTOR iar IN ASCENDING ORDER               *     
!*                  Version 1.00 - January 14, 1998                    *     
!*                                                                     *
!***********************************************************************
!
      IMPLICIT NONE
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER :: k,j,m,n,i5,maxd,iert,itemp
! 
! -------
! ... Integer arrays
! -------
! 
      INTEGER, DIMENSION(n) :: IAR
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of SORT
!
!
      m = n
! 
      DO_Outer: DO i5 = 1,n
! 
         m = m/2
         IF(m == 0) EXIT DO_Outer
! 
         maxd = n-m
! 
         DO_jLoop: DO j = 1,maxd
! 
            DO_kLoop: DO k = j,1,-m
! 
               iert = iar(k+m)
               IF(iert >= iar(k)) EXIT DO_kLoop
! 
               itemp    = iar(k+m)
               iar(k+m) = iar(k)
               iar(k)   = itemp
! 
            END DO DO_kLoop
! 
         END DO DO_jLoop
! 
      END DO DO_Outer
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of SORT
!
!
      RETURN 
!
      END SUBROUTINE SORT
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         SUBROUTINE Table_Interpolation(Xvalue, Yvalue, XX, YY, interval, extrapolate)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                 Routine for computing a polynomial                  *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(INOUT) :: interval
! 
! -------------
! ......... Integer variables
! -------------
! 
            INTEGER :: Lo, Hi, mid, m, N, Nplus1, NumVectors, increment
! 
! -------------
! ......... Double precision arrays
! -------------
! 
            REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: XX
            REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: YY
! 
! -------------
! ......... Double precision variables
! -------------
! 
            REAL(KIND=8), INTENT(IN) :: Xvalue
! 
            REAL(KIND=8), DIMENSION(SIZE(YY, DIM=2)), INTENT(OUT) :: Yvalue
! 
! -------------
! ......... Logical variables
! -------------
! 
            LOGICAL, INTENT(IN) :: extrapolate
! 
            LOGICAL :: ascending = .TRUE.
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of Table_Interpolation
!
! ......... Code for a massively parallel system
!
!           Lo = IMINLOCK(ABS(Xvalue - XX))
!
!           IF( (XX(i) > Xvalue) .EQV. ascending) THEN
!              interval = Lo
!           ELSE
!              interval = Lo - 1
!           END IF
! ......... 
! ......... Initialization
! ......... 
!
            N      = SIZE(XX)
            Nplus1 = N + 1
!
            Lo = interval
            NumVectors = SIZE(YY, DIM=2)
!
! -------------
! ......... 
! ......... Initially, hunt the interval
! ......... 
! -------------
!
            IF_InRange: IF(Lo <= 0 .OR. Lo > N) THEN
!
! ............ Initial guess is outside the table range
!
               Lo = 0
               Hi = Nplus1
!
            ELSE
!
               increment = 1
!
               IF_Hunt: IF( Xvalue >= XX(Lo) ) THEN
!
                  DO_Loop1: DO 
                     Hi = Lo + increment
                     IF( Hi > N ) THEN
                        Hi = Nplus1
                        Exit Do_Loop1
                     ELSE
                        IF( Xvalue < XX(hi) ) EXIT Do_Loop1
                        Lo = Hi
                        increment = increment + increment
                     END IF
                  END DO DO_Loop1
!
               ELSE
!
                  Hi = Lo
!
                  DO_Loop2: DO 
                     Lo = Hi - increment
                     IF( Lo < 1 ) THEN
                        Lo = 0
                        Exit Do_Loop2
                     ELSE
                        IF( Xvalue >= XX(Lo) ) EXIT Do_Loop2
                        Hi = Lo
                        increment = increment + increment
                     END IF
                  END DO DO_Loop2
!
               END IF IF_Hunt
!
            END IF IF_InRange
!
! -------------
! ......... 
! ......... Finally, bracketing and interval determination through bisection
! ......... 
! -------------
!
            DO_Bisection: DO
!
               IF_Interval: IF( (Hi - Lo) <= 1 ) THEN
!
                  IF( Xvalue == XX(N) ) Lo = N - 1
!
                  IF( Xvalue == XX(1) ) Lo = 1
!
                  Exit DO_Bisection
!
               ELSE
!
                  mid = ( Hi + Lo )/2
!
                  IF( Xvalue >= XX(mid) ) THEN
                     Lo = mid
                  ELSE
                     Hi = mid
                  END IF
!
               END IF IF_Interval
!
            END DO DO_Bisection
!
! -------------
! ......... 
! ......... Computing the interpolated value
! ......... 
! -------------
!
            interval = Lo
!
!
!
            DO m = 1,NumVectors
!
               IF( Lo == 0 ) THEN
                  IF (extrapolate) THEN
                     Yvalue(m) = YY(1,m) + (Xvalue - XX(1))*( YY(2,m) - YY(1,m) )/( XX(2) - XX(1) )
                  ELSE
                     Yvalue(m) = YY(1,m)
                  END IF
               ELSE IF ( Lo == N ) THEN
                  IF (extrapolate) THEN
                     Yvalue(m) = YY(N,m) + (Xvalue - XX(N-1))*( YY(N,m) - YY(N-1,m) )/( XX(N) - XX(N-1) )
                  ELSE
                     Yvalue(m) = YY(N,m)
                  END IF
               ELSE
                  Yvalue(m) = YY(Lo,m) + (Xvalue - XX(Lo))*( YY(Hi,m) - YY(Lo,m) )/( XX(Hi) - XX(Lo) )
               END IF
!
            END DO
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of Table_Interpolation
!
!
            RETURN
!
         END SUBROUTINE Table_Interpolation
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
         INTEGER FUNCTION Table_Interval(Xvalue, XX, interval)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*                 Routine for computing a polynomial                  *
!*                                                                     *
!*                     Version 1.0 - June 6, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
            IMPLICIT NONE
! 
! -------------
! ......... Double precision arrays
! -------------
! 
            REAL(KIND=8), DIMENSION(:), INTENT(IN) :: XX 
! 
! -------------
! ......... Double precision variables
! -------------
! 
            REAL(KIND=8), INTENT(IN)  :: Xvalue
! 
! -------------
! ......... Integer input variables
! -------------
! 
            INTEGER, INTENT(INOUT) :: interval
! 
! -------------
! ......... Integer variables
! -------------
! 
            INTEGER :: Lo, Hi, mid, N, increment
! 
! -------------
! ......... Logical variables
! -------------
! 
            LOGICAL :: ascending = .TRUE.
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of Table_Interval
!
! ......... Code for a massively parallel system
!
!           Lo = IMINLOCK(ABS(Xvalue - XX))
!
!           IF( (XX(i) > Xvalue) .EQV. ascending) THEN
!              interval = Lo
!           ELSE
!              interval = Lo - 1
!           END IF
! ......... 
! ......... Initialization
! ......... 
!
            N = SIZE(XX)
!
            Lo = interval
!
! -------------
! ......... 
! ......... Initially, hunt the interval
! ......... 
! -------------
!
            IF_InRange: IF(Lo <= 0 .OR. Lo > N) THEN
!
! ............ Initial guess is outside the table range
!
               Lo = 0
               Hi = N + 1
!
            ELSE
!
               increment = 1
!
               IF_Hunt: IF( Xvalue >= XX(Lo) ) THEN
!
                  DO_Loop1: DO 
                     Hi = Lo + increment
                     IF( Hi > N ) THEN
                        Hi = N + 1
                        Exit Do_Loop1
                     ELSE
                        IF( Xvalue < XX(hi) ) EXIT Do_Loop1
                        Lo = Hi
                        increment = increment + increment
                     END IF
                  END DO DO_Loop1
!
               ELSE
!
                  Hi = Lo
!
                  DO_Loop2: DO 
                     Lo = Hi - increment
                     IF( Lo < 1 ) THEN
                        Lo = 0
                        Exit Do_Loop2
                     ELSE
                        IF( Xvalue >= XX(Lo) ) EXIT Do_Loop2
                        Hi = Lo
                        increment = increment + increment
                     END IF
                  END DO DO_Loop2
!
               END IF IF_Hunt
!
            END IF IF_InRange
!
! -------------
! ......... 
! ......... Finally, bracketing and interval determination through bisection
! ......... 
! -------------
!
            DO_Bisection: DO
!
               IF_Interval: IF( (Hi - Lo) <= 1 ) THEN
!
                  IF( Xvalue == XX(N) ) Lo = N - 1
!
                  IF( Xvalue == XX(1) ) Lo = 1
!
                  Exit DO_Bisection
!
               ELSE
!
                  mid = ( Hi + Lo )/2
!
                  IF( Xvalue >= XX(mid) ) THEN
                     Lo = mid
                  ELSE
                     Hi = mid
                  END IF
!
               END IF IF_Interval
!
            END DO DO_Bisection
!
! -------------
! ......... 
! ......... Computing the interpolated value
! ......... 
! -------------
!
            Table_Interval = Lo
!
!
! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> End of Table_Interval
!
!
            RETURN
!
         END FUNCTION Table_Interval
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
      CHARACTER(LEN=1) FUNCTION N_Character(N)
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         ROUTINE FOR ASSIGNING CCHARACTER VALUES TO NUMBERS          *
!*                          IN ELEMENT NAMES                           *
!*                                                                     *
!*                     Version 1.00, July 14, 2007                     *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE
! 
! ----------
! ... Integer variables
! ----------
! 
      INTEGER, INTENT(IN) :: N
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  BEGIN N_Character
!
!
      SELECT CASE (N)
      CASE (0)
         N_Character = '0'
      CASE (1)
         N_Character = '1'
      CASE (2)
         N_Character = '2'
      CASE (3)
         N_Character = '3'
      CASE (4)
         N_Character = '4'
      CASE (5)
         N_Character = '5'
      CASE (6)
         N_Character = '6'
      CASE (7)
         N_Character = '7'
      CASE (8)
         N_Character = '8'
      CASE (9)
         N_Character = '9'
      END SELECT
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  END N_Character
!
!
      RETURN
!
      END FUNCTION N_Character
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
      END MODULE Utility_Functions
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
      MODULE Array_Ranking_Procedures
!
      PUBLIC :: Rank_Array
      PRIVATE :: Rank_RealKind4_Array, Rank_Integer_Array, Rank_RealKind8_Array
!
      INTERFACE Rank_Array
         MODULE PROCEDURE Rank_RealKind8_Array, Rank_RealKind4_Array, Rank_Integer_Array
      END INTERFACE Rank_Array
!
      CONTAINS

      SUBROUTINE Rank_RealKind8_Array (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Real (kind=8), Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Real (kind=8) :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine Rank_RealKind8_Array
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
Subroutine Rank_RealKind4_Array (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Real :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine Rank_RealKind4_Array
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
Subroutine Rank_Integer_Array (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Integer :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine Rank_Integer_Array
!
!
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
! 
!
!
end module Array_Ranking_Procedures
!
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
