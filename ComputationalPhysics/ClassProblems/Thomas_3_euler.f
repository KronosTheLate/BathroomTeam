      PROGRAM EULER
      IMPLICIT NONE
    
      INTEGER, PARAMETER :: N = 1001
      INTEGER, PARAMETER :: ALPHA = -1
      INTEGER :: I



      REAL :: T(N)
      REAL :: Y(N)

      Y(1) = 1

C CREATE A LINSPACE 0,10,N

      DO I = 1,N
        T(I) = 10 * REAL(I - 1) / REAL(N - 1)
        PRINT *,(T(I))
      END DO

C ADD SOME PADDING BETWEEN THE TWO PRINTS

      DO I = 1,4
        PRINT *,0
      END DO

      PRINT *,Y(1)


C USE EULER'S METHOD

      DO I = 2,N
        Y(I) = Y(I-1) * (1 + (T(I) - T(I-1)) * ALPHA)
        PRINT *,(Y(I))
      END DO


      STOP
      END PROGRAM EULER


