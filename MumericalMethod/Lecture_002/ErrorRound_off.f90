PROGRAM Main
  IMPLICIT NONE
  INTEGER, PARAMETER :: N=100000000
  REAL (KIND=4) :: X
  INTEGER      :: i 
  X=0.0
  DO i=1, N
    X=X+1.0
    IF (I > X) THEN
       PRINT*,X,i
       exit
    END IF
  END DO
  PRINT*,X,X/N
END PROGRAM Main

