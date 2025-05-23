MODULE histogram 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE                     :: hist_angle(:)
  REAL*8,ALLOCATABLE                       :: hist_dist(:)
  REAL*8,PARAMETER                         :: max_hblen=3.0
END MODULE histogram

PROGRAM ZDF
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    if (MOD(frame,stride)==0) CALL ANALYSIS
  END DO

  CALL PRINT_RESULTS

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist
  USE histogram, ONLY : hist_angle, hist_dist
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride, nhist
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(hist_angle(nhist)); hist_angle=0
  ALLOCATE(hist_dist(nhist)); hist_dist=0.0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, nhist 
  USE histogram, ONLY : hist_angle, hist_dist, max_hblen
  IMPLICIT NONE
  INTEGER                    :: i, j, ind
  INTEGER, DIMENSION(2)      :: ind_H 
  REAL*8                     :: d_OOw, Dist, Angle, n, d
  REAL*8, DIMENSION(2)       :: ang_OHO

  n = 0

  DO i = 1,natoms

    !Selecting only Ow
    IF ( atype(i) == 2 )  THEN

      ind_H = 0
      CALL get_H2(i, ind_H)

      DO j = 1,natoms

        !Selecting only Ow
        IF ((atype(j) == 2) .and. (i.ne.j)) THEN

          d_OOw = Dist(i, j)

          IF ( (ind_H(1) /= 0) .and. (ind_H(2) /= 0 )) THEN

            ang_OHO(1) = Angle(pos(:,i), pos(:,j), pos(:,ind_H(1)))
            ang_OHO(2) = Angle(pos(:,i), pos(:,j), pos(:,ind_H(2)))

            !Definition of H-bond - Luzar and Chandler
            IF (d_OOw < 3.5) THEN 
              IF (ang_OHO(1) > 0.8660) THEN
                ind=int(((ang_OHO(1))-0.8660)/0.134*nhist) + 1
                hist_angle(ind) = hist_angle(ind) + 1
                d = Dist(j,ind_H(1))
                ind=int(d/max_hblen*nhist) + 1
                hist_dist(ind) = hist_dist(ind)+1
              ELSE IF (ang_OHO(2) > 0.8660) THEN
                ind=int(((ang_OHO(2))-0.8660)/0.134*nhist) + 1
                hist_angle(ind) = hist_angle(ind) + 1
                d = Dist(j,ind_H(2))
                ind=int(d/max_hblen*nhist) + 1
                hist_dist(ind) = hist_dist(ind)+1
              END IF

            END IF !d_OOw
 
          END IF !ind_H

        END IF !ind_atom

      END DO !j

    END IF !ind_atom

  END DO !i 

END SUBROUTINE ANALYSIS

SUBROUTINE get_H2(ind_O, ind_H)
  USE parameters, ONLY : atype, pos, cutoff_OH, natoms
  IMPLICIT NONE
  INTEGER, DIMENSION(2)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1

  DO WHILE (( i <= natoms) .and. ( k <= 2))

    !Should be the index for H
    IF ( atype(i) == 3 ) THEN

      IF ( Dist(ind_O,i) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

END SUBROUTINE get_H2

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist
  USE histogram, ONLY : hist_angle, hist_dist, max_hblen
  IMPLICIT NONE
  INTEGER :: ihist
  REAL*8 :: dl

  OPEN(unit = 2, file = "HB_angle.dat")
  OPEN(unit = 3, file = "HB_length.dat")

  dl=max_hblen/float(nhist)
  hist_dist = hist_dist/sum(hist_dist)

  DO ihist=1,nhist
    !WRITE(2,fmt='(F12.8,3X,I10)') float(ihist-1)*0.5235/float(nhist)*57.3, hist_angle(ihist)
    WRITE(2,fmt='(F12.8,3X,I10)') (float(ihist-1))*0.134/float(nhist) + 0.8660, hist_angle(ihist)
    WRITE(3,fmt='(F12.8,3X,F12.8)') float(ihist-1)*dl, hist_dist(ihist)
  END DO

  CLOSE(2)
  CLOSE(3)

END SUBROUTINE PRINT_RESULTS
