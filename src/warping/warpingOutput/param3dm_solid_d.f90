!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.3 (r3163) - 09/25/2009 09:03
!
!  Differentiation of param3dm_solid in forward (tangent) mode:
!   variations  of output variables: s
!   with respect to input variables: xyz
SUBROUTINE PARAM3DM_SOLID_D(il, jl, kl, xyz, xyzd, s, sd)
  USE PRECISION
  IMPLICIT NONE
!     ******************************************************************
!     *   PARAM3DM parameterizes the volume of one block of a multi-   *
!     *   block grid structure by setting up the normalized arc-length *
!     *   increments in all three index directions.                    *
!     *                                                                *
!     *   11/29/95  D.Saunders  Adaptation of PARAMXYZ for specialized *
!     *                         WARP-BLK used by FLO107-MB.            *
!     *   06/19/96      "       Allow for degenerate edges.            *
!     *   12/11/08  C.A.Mader   Converted to *.f90                     *
!     *                                                                *
!     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
!     ******************************************************************
!IMPLICIT REAL*8 (A-H,O-Z) ! Take out when all compilers have a switch
!     Arguments.
! I Grid array dimensions.
  INTEGER(kind=inttype) :: il, jl, kl
! I Grid coordinates
  REAL(kind=realtype) :: xyz(3, 0:il+1, 0:jl+1, 0:kl+1)
  REAL(kind=realtype) :: xyzd(3, 0:il+1, 0:jl+1, 0:kl+1)
! O Normalized arc-lengths:
  REAL(kind=realtype) :: s(3, 0:il+1, 0:jl+1, 0:kl+1)
  REAL(kind=realtype) :: sd(3, 0:il+1, 0:jl+1, 0:kl+1)
!   S(1,1,J,K) = 0.,
!   S(2,I,1,K) = 0.,
!   S(3,I,J,1) = 0.,
!   S(1,IL,J,K) = 1.,etc.
!     Local constants.
  REAL :: one, zero
  PARAMETER (one=1.e+0, zero=0.e+0)
!     Local variables.
  INTEGER(kind=inttype) :: i, j, k
!     Local functions.
  REAL :: deli, delj, delk
  REAL(kind=realtype) :: arg1
  REAL(kind=realtype) :: arg1d
  REAL(kind=realtype) :: result1
  REAL(kind=realtype) :: result1d
  INTRINSIC SQRT
!     Execution.
!     ----------
!     Zero the three low-end faces (or edges if one plane is specified).
  DO k=1,kl
    DO j=1,jl
      sd(1, 1, j, k) = 0.0
      s(1, 1, j, k) = zero
    END DO
    DO i=1,il
      sd(2, i, 1, k) = 0.0
      s(2, i, 1, k) = zero
    END DO
  END DO
  DO j=1,jl
    DO i=1,il
      sd(3, i, j, 1) = 0.0
      s(3, i, j, 1) = zero
    END DO
  END DO
  sd = 0.0
!     Set up the low-end edge lines because they are missed by the
!     following loops over most of the low-end faces:
  DO i=2,il
    arg1d = 2*(xyz(1, i, 1, 1)-xyz(1, i-1, 1, 1))*(xyzd(1, i, 1, 1)-xyzd&
&      (1, i-1, 1, 1)) + 2*(xyz(2, i, 1, 1)-xyz(2, i-1, 1, 1))*(xyzd(2, i&
&      , 1, 1)-xyzd(2, i-1, 1, 1)) + 2*(xyz(3, i, 1, 1)-xyz(3, i-1, 1, 1)&
&      )*(xyzd(3, i, 1, 1)-xyzd(3, i-1, 1, 1))
    arg1 = (xyz(1, i, 1, 1)-xyz(1, i-1, 1, 1))**2 + (xyz(2, i, 1, 1)-xyz&
&      (2, i-1, 1, 1))**2 + (xyz(3, i, 1, 1)-xyz(3, i-1, 1, 1))**2
    IF (arg1 .EQ. 0.0) THEN
      result1d = 0.0
    ELSE
      result1d = arg1d/(2.0*SQRT(arg1))
    END IF
    result1 = SQRT(arg1)
    sd(1, i, 1, 1) = sd(1, i-1, 1, 1) + result1d
    s(1, i, 1, 1) = s(1, i-1, 1, 1) + result1
  END DO
  DO j=2,jl
    arg1d = 2*(xyz(1, 1, j, 1)-xyz(1, 1, j-1, 1))*(xyzd(1, 1, j, 1)-xyzd&
&      (1, 1, j-1, 1)) + 2*(xyz(2, 1, j, 1)-xyz(2, 1, j-1, 1))*(xyzd(2, 1&
&      , j, 1)-xyzd(2, 1, j-1, 1)) + 2*(xyz(3, 1, j, 1)-xyz(3, 1, j-1, 1)&
&      )*(xyzd(3, 1, j, 1)-xyzd(3, 1, j-1, 1))
    arg1 = (xyz(1, 1, j, 1)-xyz(1, 1, j-1, 1))**2 + (xyz(2, 1, j, 1)-xyz&
&      (2, 1, j-1, 1))**2 + (xyz(3, 1, j, 1)-xyz(3, 1, j-1, 1))**2
    IF (arg1 .EQ. 0.0) THEN
      result1d = 0.0
    ELSE
      result1d = arg1d/(2.0*SQRT(arg1))
    END IF
    result1 = SQRT(arg1)
    sd(2, 1, j, 1) = sd(2, 1, j-1, 1) + result1d
    s(2, 1, j, 1) = s(2, 1, j-1, 1) + result1
  END DO
  DO k=2,kl
    arg1d = 2*(xyz(1, 1, 1, k)-xyz(1, 1, 1, k-1))*(xyzd(1, 1, 1, k)-xyzd&
&      (1, 1, 1, k-1)) + 2*(xyz(2, 1, 1, k)-xyz(2, 1, 1, k-1))*(xyzd(2, 1&
&      , 1, k)-xyzd(2, 1, 1, k-1)) + 2*(xyz(3, 1, 1, k)-xyz(3, 1, 1, k-1)&
&      )*(xyzd(3, 1, 1, k)-xyzd(3, 1, 1, k-1))
    arg1 = (xyz(1, 1, 1, k)-xyz(1, 1, 1, k-1))**2 + (xyz(2, 1, 1, k)-xyz&
&      (2, 1, 1, k-1))**2 + (xyz(3, 1, 1, k)-xyz(3, 1, 1, k-1))**2
    IF (arg1 .EQ. 0.0) THEN
      result1d = 0.0
    ELSE
      result1d = arg1d/(2.0*SQRT(arg1))
    END IF
    result1 = SQRT(arg1)
    sd(3, 1, 1, k) = sd(3, 1, 1, k-1) + result1d
    s(3, 1, 1, k) = s(3, 1, 1, k-1) + result1
  END DO
!     Set up the rest of the low-end face lines because they are
!     missed by the the main loop over most of the volume.
  DO k=2,kl
    DO j=2,jl
      arg1d = 2*(xyz(1, 1, j, k)-xyz(1, 1, j-1, k))*(xyzd(1, 1, j, k)-&
&        xyzd(1, 1, j-1, k)) + 2*(xyz(2, 1, j, k)-xyz(2, 1, j-1, k))*(&
&        xyzd(2, 1, j, k)-xyzd(2, 1, j-1, k)) + 2*(xyz(3, 1, j, k)-xyz(3&
&        , 1, j-1, k))*(xyzd(3, 1, j, k)-xyzd(3, 1, j-1, k))
      arg1 = (xyz(1, 1, j, k)-xyz(1, 1, j-1, k))**2 + (xyz(2, 1, j, k)-&
&        xyz(2, 1, j-1, k))**2 + (xyz(3, 1, j, k)-xyz(3, 1, j-1, k))**2
      IF (arg1 .EQ. 0.0) THEN
        result1d = 0.0
      ELSE
        result1d = arg1d/(2.0*SQRT(arg1))
      END IF
      result1 = SQRT(arg1)
      sd(2, 1, j, k) = sd(2, 1, j-1, k) + result1d
      s(2, 1, j, k) = s(2, 1, j-1, k) + result1
      arg1d = 2*(xyz(1, 1, j, k)-xyz(1, 1, j, k-1))*(xyzd(1, 1, j, k)-&
&        xyzd(1, 1, j, k-1)) + 2*(xyz(2, 1, j, k)-xyz(2, 1, j, k-1))*(&
&        xyzd(2, 1, j, k)-xyzd(2, 1, j, k-1)) + 2*(xyz(3, 1, j, k)-xyz(3&
&        , 1, j, k-1))*(xyzd(3, 1, j, k)-xyzd(3, 1, j, k-1))
      arg1 = (xyz(1, 1, j, k)-xyz(1, 1, j, k-1))**2 + (xyz(2, 1, j, k)-&
&        xyz(2, 1, j, k-1))**2 + (xyz(3, 1, j, k)-xyz(3, 1, j, k-1))**2
      IF (arg1 .EQ. 0.0) THEN
        result1d = 0.0
      ELSE
        result1d = arg1d/(2.0*SQRT(arg1))
      END IF
      result1 = SQRT(arg1)
      sd(3, 1, j, k) = sd(3, 1, j, k-1) + result1d
      s(3, 1, j, k) = s(3, 1, j, k-1) + result1
    END DO
    DO i=2,il
      arg1d = 2*(xyz(1, i, 1, k)-xyz(1, i-1, 1, k))*(xyzd(1, i, 1, k)-&
&        xyzd(1, i-1, 1, k)) + 2*(xyz(2, i, 1, k)-xyz(2, i-1, 1, k))*(&
&        xyzd(2, i, 1, k)-xyzd(2, i-1, 1, k)) + 2*(xyz(3, i, 1, k)-xyz(3&
&        , i-1, 1, k))*(xyzd(3, i, 1, k)-xyzd(3, i-1, 1, k))
      arg1 = (xyz(1, i, 1, k)-xyz(1, i-1, 1, k))**2 + (xyz(2, i, 1, k)-&
&        xyz(2, i-1, 1, k))**2 + (xyz(3, i, 1, k)-xyz(3, i-1, 1, k))**2
      IF (arg1 .EQ. 0.0) THEN
        result1d = 0.0
      ELSE
        result1d = arg1d/(2.0*SQRT(arg1))
      END IF
      result1 = SQRT(arg1)
      sd(1, i, 1, k) = sd(1, i-1, 1, k) + result1d
      s(1, i, 1, k) = s(1, i-1, 1, k) + result1
      arg1d = 2*(xyz(1, i, 1, k)-xyz(1, i, 1, k-1))*(xyzd(1, i, 1, k)-&
&        xyzd(1, i, 1, k-1)) + 2*(xyz(2, i, 1, k)-xyz(2, i, 1, k-1))*(&
&        xyzd(2, i, 1, k)-xyzd(2, i, 1, k-1)) + 2*(xyz(3, i, 1, k)-xyz(3&
&        , i, 1, k-1))*(xyzd(3, i, 1, k)-xyzd(3, i, 1, k-1))
      arg1 = (xyz(1, i, 1, k)-xyz(1, i, 1, k-1))**2 + (xyz(2, i, 1, k)-&
&        xyz(2, i, 1, k-1))**2 + (xyz(3, i, 1, k)-xyz(3, i, 1, k-1))**2
      IF (arg1 .EQ. 0.0) THEN
        result1d = 0.0
      ELSE
        result1d = arg1d/(2.0*SQRT(arg1))
      END IF
      result1 = SQRT(arg1)
      sd(3, i, 1, k) = sd(3, i, 1, k-1) + result1d
      s(3, i, 1, k) = s(3, i, 1, k-1) + result1
    END DO
  END DO
  DO j=2,jl
    DO i=2,il
      arg1d = 2*(xyz(1, i, j, 1)-xyz(1, i-1, j, 1))*(xyzd(1, i, j, 1)-&
&        xyzd(1, i-1, j, 1)) + 2*(xyz(2, i, j, 1)-xyz(2, i-1, j, 1))*(&
&        xyzd(2, i, j, 1)-xyzd(2, i-1, j, 1)) + 2*(xyz(3, i, j, 1)-xyz(3&
&        , i-1, j, 1))*(xyzd(3, i, j, 1)-xyzd(3, i-1, j, 1))
      arg1 = (xyz(1, i, j, 1)-xyz(1, i-1, j, 1))**2 + (xyz(2, i, j, 1)-&
&        xyz(2, i-1, j, 1))**2 + (xyz(3, i, j, 1)-xyz(3, i-1, j, 1))**2
      IF (arg1 .EQ. 0.0) THEN
        result1d = 0.0
      ELSE
        result1d = arg1d/(2.0*SQRT(arg1))
      END IF
      result1 = SQRT(arg1)
      sd(1, i, j, 1) = sd(1, i-1, j, 1) + result1d
      s(1, i, j, 1) = s(1, i-1, j, 1) + result1
      arg1d = 2*(xyz(1, i, j, 1)-xyz(1, i, j-1, 1))*(xyzd(1, i, j, 1)-&
&        xyzd(1, i, j-1, 1)) + 2*(xyz(2, i, j, 1)-xyz(2, i, j-1, 1))*(&
&        xyzd(2, i, j, 1)-xyzd(2, i, j-1, 1)) + 2*(xyz(3, i, j, 1)-xyz(3&
&        , i, j-1, 1))*(xyzd(3, i, j, 1)-xyzd(3, i, j-1, 1))
      arg1 = (xyz(1, i, j, 1)-xyz(1, i, j-1, 1))**2 + (xyz(2, i, j, 1)-&
&        xyz(2, i, j-1, 1))**2 + (xyz(3, i, j, 1)-xyz(3, i, j-1, 1))**2
      IF (arg1 .EQ. 0.0) THEN
        result1d = 0.0
      ELSE
        result1d = arg1d/(2.0*SQRT(arg1))
      END IF
      result1 = SQRT(arg1)
      sd(2, i, j, 1) = sd(2, i, j-1, 1) + result1d
      s(2, i, j, 1) = s(2, i, j-1, 1) + result1
    END DO
  END DO
!     Traverse the block just once for all lines except those within
!     the low-end faces.
  DO k=2,kl
    DO j=2,jl
      DO i=2,il
        arg1d = 2*(xyz(1, i, j, k)-xyz(1, i-1, j, k))*(xyzd(1, i, j, k)-&
&          xyzd(1, i-1, j, k)) + 2*(xyz(2, i, j, k)-xyz(2, i-1, j, k))*(&
&          xyzd(2, i, j, k)-xyzd(2, i-1, j, k)) + 2*(xyz(3, i, j, k)-xyz(&
&          3, i-1, j, k))*(xyzd(3, i, j, k)-xyzd(3, i-1, j, k))
        arg1 = (xyz(1, i, j, k)-xyz(1, i-1, j, k))**2 + (xyz(2, i, j, k)&
&          -xyz(2, i-1, j, k))**2 + (xyz(3, i, j, k)-xyz(3, i-1, j, k))**&
&          2
        IF (arg1 .EQ. 0.0) THEN
          result1d = 0.0
        ELSE
          result1d = arg1d/(2.0*SQRT(arg1))
        END IF
        result1 = SQRT(arg1)
        sd(1, i, j, k) = sd(1, i-1, j, k) + result1d
        s(1, i, j, k) = s(1, i-1, j, k) + result1
        arg1d = 2*(xyz(1, i, j, k)-xyz(1, i, j-1, k))*(xyzd(1, i, j, k)-&
&          xyzd(1, i, j-1, k)) + 2*(xyz(2, i, j, k)-xyz(2, i, j-1, k))*(&
&          xyzd(2, i, j, k)-xyzd(2, i, j-1, k)) + 2*(xyz(3, i, j, k)-xyz(&
&          3, i, j-1, k))*(xyzd(3, i, j, k)-xyzd(3, i, j-1, k))
        arg1 = (xyz(1, i, j, k)-xyz(1, i, j-1, k))**2 + (xyz(2, i, j, k)&
&          -xyz(2, i, j-1, k))**2 + (xyz(3, i, j, k)-xyz(3, i, j-1, k))**&
&          2
        IF (arg1 .EQ. 0.0) THEN
          result1d = 0.0
        ELSE
          result1d = arg1d/(2.0*SQRT(arg1))
        END IF
        result1 = SQRT(arg1)
        sd(2, i, j, k) = sd(2, i, j-1, k) + result1d
        s(2, i, j, k) = s(2, i, j-1, k) + result1
        arg1d = 2*(xyz(1, i, j, k)-xyz(1, i, j, k-1))*(xyzd(1, i, j, k)-&
&          xyzd(1, i, j, k-1)) + 2*(xyz(2, i, j, k)-xyz(2, i, j, k-1))*(&
&          xyzd(2, i, j, k)-xyzd(2, i, j, k-1)) + 2*(xyz(3, i, j, k)-xyz(&
&          3, i, j, k-1))*(xyzd(3, i, j, k)-xyzd(3, i, j, k-1))
        arg1 = (xyz(1, i, j, k)-xyz(1, i, j, k-1))**2 + (xyz(2, i, j, k)&
&          -xyz(2, i, j, k-1))**2 + (xyz(3, i, j, k)-xyz(3, i, j, k-1))**&
&          2
        IF (arg1 .EQ. 0.0) THEN
          result1d = 0.0
        ELSE
          result1d = arg1d/(2.0*SQRT(arg1))
        END IF
        result1 = SQRT(arg1)
        sd(3, i, j, k) = sd(3, i, j, k-1) + result1d
        s(3, i, j, k) = s(3, i, j, k-1) + result1
      END DO
    END DO
  END DO
!     Normalizing requires another pass through the volume.
!     Handle lines of zero length first by inserting uniform
!     distributions.  Then the standard normalization can be
!     applied safely everywhere.
  DO k=1,kl
!        Zero-length lines in the I direction?
    DO j=1,jl
      IF (s(1, il, j, k) .EQ. zero) THEN
        DO i=2,il
          sd(1, i, j, k) = 0.0
          s(1, i, j, k) = i - 1
        END DO
      END IF
    END DO
!        Zero-length lines in the J direction?
    DO i=1,il
      IF (s(2, i, jl, k) .EQ. zero) THEN
        DO j=2,jl
          sd(2, i, j, k) = 0.0
          s(2, i, j, k) = j - 1
        END DO
      END IF
    END DO
  END DO
!     Zero-length lines in the K direction?
  DO j=1,jl
    DO i=1,il
      IF (s(3, i, j, kl) .EQ. zero) THEN
        DO k=2,kl
          sd(3, i, j, k) = 0.0
          s(3, i, j, k) = k - 1
        END DO
      END IF
    END DO
  END DO
!     Normalize:
  DO k=1,kl
    DO j=1,jl
      DO i=1,il
        sd(1, i, j, k) = (sd(1, i, j, k)*s(1, il, j, k)-s(1, i, j, k)*sd&
&          (1, il, j, k))/s(1, il, j, k)**2
        s(1, i, j, k) = s(1, i, j, k)/s(1, il, j, k)
        sd(2, i, j, k) = (sd(2, i, j, k)*s(2, i, jl, k)-s(2, i, j, k)*sd&
&          (2, i, jl, k))/s(2, i, jl, k)**2
        s(2, i, j, k) = s(2, i, j, k)/s(2, i, jl, k)
        sd(3, i, j, k) = (sd(3, i, j, k)*s(3, i, j, kl)-s(3, i, j, k)*sd&
&          (3, i, j, kl))/s(3, i, j, kl)**2
        s(3, i, j, k) = s(3, i, j, k)/s(3, i, j, kl)
      END DO
    END DO
  END DO
!     Finally, precise 1s for the three high-end faces:
  DO k=1,kl
    DO j=1,jl
      sd(1, il, j, k) = 0.0
      s(1, il, j, k) = one
    END DO
    DO i=1,il
      sd(2, i, jl, k) = 0.0
      s(2, i, jl, k) = one
    END DO
  END DO
  DO j=1,jl
    DO i=1,il
      sd(3, i, j, kl) = 0.0
      s(3, i, j, kl) = one
    END DO
  END DO
END SUBROUTINE PARAM3DM_SOLID_D