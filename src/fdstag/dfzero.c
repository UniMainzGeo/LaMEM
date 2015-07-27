/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This sofware was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   dfzero.c
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

//---------------------------------------------------------------------------
//..................... DFZERO ROOT FINDING ALGORITHM .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "dfzero.h"
//---------------------------------------------------------------------------
void DFZERO(
	PetscScalar (*F)(PetscScalar, void *), // nonlinear function with parameter and context
	void         *FCTX,                    // pointer to a function evaluation context
	PetscScalar  *_B,                      // left bound of root interval (output)
	PetscScalar  *_C,                      // right bound of root interval (output)
	PetscScalar    R,                      // initial guess
	PetscScalar    RE,                     // relative tolerance
	PetscScalar    AE,                     // absolute tolerance
	PetscInt     *_IFLAG)                  // error code (output)
{
    /*=======================================================================
     *   BEGIN PROLOGUE  DFZERO
     *   DATE WRITTEN          700901   (YYMMDD)
     *   REVISION DATE         861211   (YYMMDD)
     *   CATEGORY NO.  F1B
     *   KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(FZERO-S DFZERO-D),
     *             BISECTION,NONLINEAR,NONLINEAR EQUATIONS,ROOTS,ZEROES,
     *             ZEROS
     *   AUTHOR  SHAMPINE,L.F.,SNLA
     *           WATTS,H.A.,SNLA
     *   PURPOSE  Search for a zero of a function F(X) in a given
     *            interval (B,C).  It is designed primarily for problems
     *            where F(B) and F(C) have opposite signs.
     *   DESCRIPTION
     *
     *       **** Double Precision version of FZERO ****
     *
     *     Based on a method by T. J. Dekker
     *     written by L. F. Shampine and H. A. Watts
     *
     *            DFZERO searches for a zero of a function F(X) between
     *            the given values B and C until the width of the interval
     *            (B,C) has collapsed to within a tolerance specified by
     *            the stopping criterion, DABS(B-C) .LE. 2.*(RW*DABS(B)+AE).
     *            The method used is an efficient combination of bisection
     *            and the secant rule.
     *
     *     Description Of Arguments
     *
     *     F,B,C,R,RE and AE are DOUBLE PRECISION input parameters.
     *     B and C are DOUBLE PRECISION output parameters and IFLAG (flagged
     *        by an * below).
     *
     *        F     - Name of the DOUBLE PRECISION valued external function.
     *                This name must be in an EXTERNAL statement in the
     *                calling program.  F must be a function of one double
     *                precision argument.
     *
     *       *B     - One end of the interval (B,C).  The value returned for
     *                B usually is the better approximation to a zero of F.
     *
     *       *C     - The other end of the interval (B,C)
     *
     *        R     - A (better) guess of a zero of F which could help in
     *                speeding up convergence.  If F(B) and F(R) have
     *                opposite signs, a root will be found in the interval
     *                (B,R); if not, but F(R) and F(C) have opposite
     *                signs, a root will be found in the interval (R,C);
     *                otherwise, the interval (B,C) will be searched for a
     *                possible root.  When no better guess is known, it is
     *                recommended that r be set to B or C; because if R is
     *                not interior to the interval (B,C), it will be ignored.
     *
     *        RE    - Relative error used for RW in the stopping criterion.
     *                If the requested RE is less than machine precision,
     *                then RW is set to approximately machine precision.
     *
     *        AE    - Absolute error used in the stopping criterion.  If the
     *                given interval (B,C) contains the origin, then a
     *                nonzero value should be chosen for AE.
     *
     *       *IFLAG - A status code.  User must check IFLAG after each call.
     *                Control returns to the user from FZERO in all cases.
     *                XERROR does not process diagnostics in these cases.
     *
     *                1  B is within the requested tolerance of a zero.
     *                   The interval (B,C) collapsed to the requested
     *                   tolerance, the function changes sign in (B,C), and
     *                   F(X) decreased in magnitude as (B,C) collapsed.
     *
     *                2  F(B) = 0.  However, the interval (B,C) may not have
     *                   collapsed to the requested tolerance.
     *
     *                3  B may be near a singular point of F(X).
     *                   The interval (B,C) collapsed to the requested
     *                   tolerance and the function changes sign in (B,C), but
     *                   F(X) increased in magnitude as (B,C) collapsed,i.e.
     *                     abs(F(B out)) .GT. max(abs(F(B in)),abs(F(C in)))
     *
     *                4  No change in sign of F(X) was found although the
     *                   interval (B,C) collapsed to the requested tolerance.
     *                   The user must examine this case and decide whether
     *                   B is near a local minimum of F(X), or B is near a
     *                   zero of even multiplicity, or neither of these.
     *
     *                5  Too many (.GT. 500) function evaluations used.
     *
     *   REFERENCES  L. F. SHAMPINE AND H. A. WATTS, *FZERO, A ROOT-SOLVING
     *                 CODE*, SC-TM-70-631, SEPTEMBER 1970.
     *               T. J. DEKKER, *FINDING A ZERO BY MEANS OF SUCCESSIVE
     *                 LINEAR INTERPOLATION*, 'CONSTRUCTIVE ASPECTS OF THE
     *                 FUNDAMENTAL THEOREM OF ALGEBRA', EDITED BY B. DEJON
     *                 P. HENRICI, 1969.
     *
     *   END PROLOGUE DFZERO
     *=====================================================================*/
    PetscInt    IC, IFLAG, KOUNT;
    PetscScalar A, ACBS, ACMB, AW, B, C, CMB, ER, FA, FB, FC, FX, FZ, P, Q, RW, T, TOL, Z;
      // BEGIN BLOCK PERMITTING ...EXITS TO 200
      // ER IS TWO TIMES THE COMPUTER UNIT ROUNDOFF VALUE
      // FOR DOUBLE PRECISION ARITHMETIC SHOULD BE 2^{-52} = 2.22e-16.
      // FIRST EXECUTABLE STATEMENT  DFZERO
         ER = 2.0*DBL_EPSILON;
      // INITIALIZE
         B = *_B;
         C = *_C;
         Z = R;
         if (R  <=  FDMIN(B, C)  ||  R  >=  FDMAX(B, C)) Z = C;
         RW = FDMAX(RE, ER);
         AW = FDMAX(AE, 0.0);
         IC = 0;
         T = Z;
         FZ = F(T, FCTX);
         FC = FZ;
         T = B;
         FB = F(T, FCTX);
         KOUNT = 2;
         if (FDSIGN(1.0, FZ)  ==  FDSIGN(1.0, FB)) goto p10;
            C = Z;
         goto p30;
p10:  // CONTINUE
      //    BEGIN BLOCK PERMITTING ...EXITS TO 20
      //    ...EXIT
               if (Z  ==  C) goto p20;
               T = C;
               FC = F(T, FCTX);
               KOUNT = 3;
      //    ...EXIT
               if (FDSIGN(1.0, FZ)  ==  FDSIGN(1.0, FC)) goto p20;
               B = Z;
               FB = FZ;
p20:  //    CONTINUE
p30:  // CONTINUE
         A = C;
         FA = FC;
         ACBS = PetscAbsScalar(B - C);
         FX = FDMAX(PetscAbsScalar(FB), PetscAbsScalar(FC));
p40:  // CONTINUE
      //    BEGIN BLOCK PERMITTING ...EXITS TO 180
               if (PetscAbsScalar(FC)  >=  PetscAbsScalar(FB)) goto p50;
      //          PERFORM INTERCHANGE
                  A = B;
                  FA = FB;
                  B = C;
                  FB = FC;
                  C = A;
                  FC = FA;
p50:  //       CONTINUE
               CMB = 0.5*(C - B);
               ACMB = PetscAbsScalar(CMB);
               TOL = RW*PetscAbsScalar(B) + AW;
      //       TEST STOPPING CRITERION AND FUNCTION COUNT
               if (ACMB   >  TOL) goto p90;
      //          FINISHED. PROCESS RESULTS FOR PROPER SETTING OF IFLAG
                  if (FDSIGN(1.0, FB)  !=  FDSIGN(1.0, FC)) goto p60;
                     IFLAG = 4;
                  goto p80;
p60:  //          CONTINUE
                  if (PetscAbsScalar(FB)  <=  FX) goto p70;
                     IFLAG = 3;
                  goto p80;
p70:  //          CONTINUE
                     IFLAG = 1;
p80:  //          CONTINUE
      //..........EXIT
                  goto p200;
p90:  //       CONTINUE
               if (FB  !=  0.0) goto p100;
                  IFLAG = 2;
      //..........EXIT
                  goto p200;
p100: //       CONTINUE
               if (KOUNT  <   500) goto p110;
                  IFLAG = 5;
      //..........EXIT
                  goto p200;
p110: //       CONTINUE
      //       CALCULATE NEW ITERATE IMPLICITLY AS B+P/Q
      //       WHERE WE ARRANGE P .GE. 0.
      //       THE IMPLICIT FORM IS USED TO PREVENT OVERFLOW.
               P = (B - A)*FB;
               Q = FA - FB;
               if (P  >=  0.0) goto p120;
                  P = -P;
                  Q = -Q;
p120: //       CONTINUE
      //       UPDATE A AND CHECK FOR SATISFACTORY REDUCTION
      //       IN THE SIZE OF THE BRACKETING INTERVAL.
      //       IF NOT, PERFORM BISECTION.
               A = B;
               FA = FB;
               IC = IC + 1;
               if (IC  <   4) goto p140;
                  if (8.0*ACMB  <   ACBS) goto p130;
      //             USE BISECTION
                     B = 0.5*(C + B);
      //    .........EXIT
                     goto p180;
p130: //          CONTINUE
                  IC = 0;
                  ACBS = ACMB;
p140: //       CONTINUE
      //       TEST FOR TOO SMALL A CHANGE
               if (P  >   PetscAbsScalar(Q)*TOL) goto p150;
      //          INCREMENT BY TOLERANCE
                  B = B + FDSIGN(TOL, CMB);
               goto p170;
p150: //       CONTINUE
      //       ROOT OUGHT TO BE BETWEEN B AND (C+B)/2.
               if (P   <  CMB*Q) goto p160;
      //          USE BISECTION
                  B = 0.5*(C + B);
               goto p170;
p160: //       CONTINUE
      //          USE SECANT RULE
                  B = B + P/Q;
p170: //       CONTINUE
p180: //    CONTINUE
      //    HAVE COMPLETED COMPUTATION FOR NEW ITERATE B
            T = B;
            FB = F(T, FCTX);
            KOUNT = KOUNT + 1;
      //    DECIDE WHETHER NEXT STEP IS INTERPOLATION OR EXTRAPOLATION
            if (FDSIGN(1.0, FB)  !=  FDSIGN(1.0, FC)) goto p190;
               C = A;
               FC = FA;
p190: //    CONTINUE
         goto p40;
p200: // CONTINUE
      // FINALIZE
         *_B = B;
         *_C = C;
         *_IFLAG = IFLAG;
}
/*---------------------------------------------------------------------------
// DFZERO ORIGINAL FORTRAN CODE
//---------------------------------------------------------------------------
SUBROUTINE DFZERO(F,B,C,R,RE,AE,IFLAG)
C***BEGIN PROLOGUE  DFZERO
C***DATE WRITTEN   700901   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  F1B
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(FZERO-S DFZERO-D),
C             BISECTION,NONLINEAR,NONLINEAR EQUATIONS,ROOTS,ZEROES,
C             ZEROS
C***AUTHOR  SHAMPINE,L.F.,SNLA
C           WATTS,H.A.,SNLA
C***PURPOSE  Search for a zero of a function F(X) in a given
C            interval (B,C).  It is designed primarily for problems
C            where F(B) and F(C) have opposite signs.
C***DESCRIPTION
C
C       **** Double Precision version of FZERO ****
C
C     Based on a method by T J Dekker
C     written by L F Shampine and H A Watts
C
C            DFZERO searches for a zero of a function F(X) between
C            the given values B and C until the width of the interval
C            (B,C) has collapsed to within a tolerance specified by
C            the stopping criterion, DABS(B-C) .LE. 2.*(RW*DABS(B)+AE).
C            The method used is an efficient combination of bisection
C            and the secant rule.
C
C     Description Of Arguments
C
C     F,B,C,R,RE and AE are DOUBLE PRECISION input parameters.
C     B and C are DOUBLE PRECISION output parameters and IFLAG (flagged
C        by an * below).
C
C        F     - Name of the DOUBLE PRECISION valued external function.
C                This name must be in an EXTERNAL statement in the
C                calling program.  F must be a function of one double
C                precision argument.
C
C       *B     - One end of the interval (B,C).  The value returned for
C                B usually is the better approximation to a zero of F.
C
C       *C     - The other end of the interval (B,C)
C
C        R     - A (better) guess of a zero of F which could help in
C                speeding up convergence.  If F(B) and F(R) have
C                opposite signs, a root will be found in the interval
C                (B,R); if not, but F(R) and F(C) have opposite
C                signs, a root will be found in the interval (R,C);
C                otherwise, the interval (B,C) will be searched for a
C                possible root.  When no better guess is known, it is
C                recommended that r be set to B or C; because if R is
C                not interior to the interval (B,C), it will be ignored.
C
C        RE    - Relative error used for RW in the stopping criterion.
C                If the requested RE is less than machine precision,
C                then RW is set to approximately machine precision.
C
C        AE    - Absolute error used in the stopping criterion.  If the
C                given interval (B,C) contains the origin, then a
C                nonzero value should be chosen for AE.
C
C       *IFLAG - A status code.  User must check IFLAG after each call.
C                Control returns to the user from FZERO in all cases.
C                XERROR does not process diagnostics in these cases.
C
C                1  B is within the requested tolerance of a zero.
C                   The interval (B,C) collapsed to the requested
C                   tolerance, the function changes sign in (B,C), and
C                   F(X) decreased in magnitude as (B,C) collapsed.
C
C                2  F(B) = 0.  However, the interval (B,C) may not have
C                   collapsed to the requested tolerance.
C
C                3  B may be near a singular point of F(X).
C                   The interval (B,C) collapsed to the requested tol-
C                   erance and the function changes sign in (B,C), but
C                   F(X) increased in magnitude as (B,C) collapsed,i.e.
C                     abs(F(B out)) .GT. max(abs(F(B in)),abs(F(C in)))
C
C                4  No change in sign of F(X) was found although the
C                   interval (B,C) collapsed to the requested tolerance.
C                   The user must examine this case and decide whether
C                   B is near a local minimum of F(X), or B is near a
C                   zero of even multiplicity, or neither of these.
C
C                5  Too many (.GT. 500) function evaluations used.
C***REFERENCES  L. F. SHAMPINE AND H. A. WATTS, *FZERO, A ROOT-SOLVING
C                 CODE*, SC-TM-70-631, SEPTEMBER 1970.
C               T. J. DEKKER, *FINDING A ZERO BY MEANS OF SUCCESSIVE
C                 LINEAR INTERPOLATION*, 'CONSTRUCTIVE ASPECTS OF THE
C                 FUNDAMENTAL THEOREM OF ALGEBRA', EDITED BY B. DEJON
C                 P. HENRICI, 1969.
C***ROUTINES CALLED  D1MACH
C***END PROLOGUE  DFZERO
      INTEGER I, IC, ICNT, IERR, IFLAG, IPASS, IPSS, ITEST(36),
     *     ITMP(15), J, KLUS, KOUNT, KPRINT, LUN, NDEG
      DOUBLE PRECISION A, ACBS, ACMB, AE, AW, B, C,
     *     CMB, D1MACH, DABS,
     *     DMAX1, DMIN1, DSIGN1, DSQRT, ER, F, FA, FB, FC,
     *     FX, FZ, P, Q, R, RE, REL, RW, T, TOL, WI, WORK, WR, Z
C     BEGIN BLOCK PERMITTING ...EXITS TO 200
C          ER IS TWO TIMES THE COMPUTER UNIT ROUNDOFF VALUE WHICH IS
C          DEFINED HERE BY THE FUNCTION D1MACH.
C***FIRST EXECUTABLE STATEMENT  DFZERO
ccc         ER = 2.0D0*D1MACH(4)

ccc   For IEEE 754 double precision arithmetic
ccc   d1mach(4) should be 2^{-53} x 2 = 2.22e-16.
ccc   Added by Steve Verrill on 5/25/02.

         er = 2.0d0*2.22d-16

C
C        INITIALIZE
C
         Z = R
         IF (R .LE. DMIN1(B,C) .OR. R .GE. DMAX1(B,C)) Z = C
         RW = DMAX1(RE,ER)
         AW = DMAX1(AE,0.0D0)
         IC = 0
         T = Z
         FZ = F(T)
         FC = FZ
         T = B
         FB = F(T)
         KOUNT = 2
         IF (DSIGN(1.0D0,FZ) .EQ. DSIGN(1.0D0,FB)) GO TO 10
            C = Z
         GO TO 30
   10    CONTINUE
C           BEGIN BLOCK PERMITTING ...EXITS TO 20
C           ...EXIT
               IF (Z .EQ. C) GO TO 20
               T = C
               FC = F(T)
               KOUNT = 3
C           ...EXIT
               IF (DSIGN(1.0D0,FZ) .EQ. DSIGN(1.0D0,FC)) GO TO 20
               B = Z
               FB = FZ
   20       CONTINUE
   30    CONTINUE
         A = C
         FA = FC
         ACBS = DABS(B-C)
         FX = DMAX1(DABS(FB),DABS(FC))
C
   40    CONTINUE
C           BEGIN BLOCK PERMITTING ...EXITS TO 180
               IF (DABS(FC) .GE. DABS(FB)) GO TO 50
C                 PERFORM INTERCHANGE
                  A = B
                  FA = FB
                  B = C
                  FB = FC
                  C = A
                  FC = FA
   50          CONTINUE
C
               CMB = 0.5D0*(C - B)
               ACMB = DABS(CMB)
               TOL = RW*DABS(B) + AW
C
C              TEST STOPPING CRITERION AND FUNCTION COUNT
C
               IF (ACMB .GT. TOL) GO TO 90
C
C
C                 FINISHED. PROCESS RESULTS FOR PROPER SETTING OF IFLAG
C
                  IF (DSIGN(1.0D0,FB) .NE. DSIGN(1.0D0,FC)) GO TO 60
                     IFLAG = 4
                  GO TO 80
   60             CONTINUE
                  IF (DABS(FB) .LE. FX) GO TO 70
                     IFLAG = 3
                  GO TO 80
   70             CONTINUE
                     IFLAG = 1
   80             CONTINUE
C     ............EXIT
                  GO TO 200
   90          CONTINUE
               IF (FB .NE. 0.0D0) GO TO 100
                  IFLAG = 2
C     ............EXIT
                  GO TO 200
  100          CONTINUE
               IF (KOUNT .LT. 500) GO TO 110
                  IFLAG = 5
C     ............EXIT
                  GO TO 200
  110          CONTINUE
C
C              CALCULATE NEW ITERATE IMPLICITLY AS B+P/Q
C              WHERE WE ARRANGE P .GE. 0.
C              THE IMPLICIT FORM IS USED TO PREVENT OVERFLOW.
C
               P = (B - A)*FB
               Q = FA - FB
               IF (P .GE. 0.0D0) GO TO 120
                  P = -P
                  Q = -Q
  120          CONTINUE
C
C              UPDATE A AND CHECK FOR SATISFACTORY REDUCTION
C              IN THE SIZE OF THE BRACKETING INTERVAL.
C              IF NOT, PERFORM BISECTION.
C
               A = B
               FA = FB
               IC = IC + 1
               IF (IC .LT. 4) GO TO 140
                  IF (8.0D0*ACMB .LT. ACBS) GO TO 130
C
C                    USE BISECTION
C
                     B = 0.5D0*(C + B)
C           .........EXIT
                     GO TO 180
  130             CONTINUE
                  IC = 0
                  ACBS = ACMB
  140          CONTINUE
C
C              TEST FOR TOO SMALL A CHANGE
C
               IF (P .GT. DABS(Q)*TOL) GO TO 150
C
C                 INCREMENT BY TOLERANCE
C
                  B = B + DSIGN(TOL,CMB)
               GO TO 170
  150          CONTINUE
C
C              ROOT OUGHT TO BE BETWEEN B AND (C+B)/2.
C
               IF (P .LT. CMB*Q) GO TO 160
C
C                 USE BISECTION
C
                  B = 0.5D0*(C + B)
               GO TO 170
  160          CONTINUE
C
C                 USE SECANT RULE
C
                  B = B + P/Q
  170          CONTINUE
  180       CONTINUE
C
C           HAVE COMPLETED COMPUTATION FOR NEW ITERATE B
C
            T = B
            FB = F(T)
            KOUNT = KOUNT + 1
C
C           DECIDE WHETHER NEXT STEP IS INTERPOLATION OR EXTRAPOLATION
C
            IF (DSIGN(1.0D0,FB) .NE. DSIGN(1.0D0,FC)) GO TO 190
               C = A
               FC = FA
  190       CONTINUE
         GO TO 40
  200 CONTINUE
      RETURN
      END
//-------------------------------------------------------------------------*/
