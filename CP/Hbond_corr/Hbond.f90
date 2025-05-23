!Computes the average and instantaneous Hbond

MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: neq = 12
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: pos
  REAL*8, DIMENSION(9)                     :: box
  REAL*8, PARAMETER                        :: cutoff_OH = 1.3D0 
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE   :: hbond_corr
  INTEGER                                  :: natoms, nframes
  INTEGER*8                                :: n_H, n_H2
  INTEGER                                  :: nequil

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  !Remove equilibration part
  DO frame = 1,nequil*(natoms+2)
    READ(1, *)
  END DO
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL MAKE_HISTOGRAM (frame-nequil)
  END DO

  CALL COARSE_GRAIN
  CALL CORR_FUNCT

  CLOSE(1);CLOSE(2)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  READ(1, *), box
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms))
  ALLOCATE(hbond_corr(neq,nframes-nequil,2))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO
 
  CLOSE(1)

  OPEN(unit = 1,file = file_name)
  OPEN(unit = 2,file = "H_Corr.dat")
 

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(5)               :: junk2

  READ(1,*)
  READ(1,*)

  DO i = 1,natoms
    READ(1, *), junk2, pos(i,:)
  END DO

END SUBROUTINE

SUBROUTINE MAKE_HISTOGRAM (frame)
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j, frame
  INTEGER                    :: index_O2c
  INTEGER, DIMENSION(2)      :: ind_H 
  REAL*8                     :: d_OOw, Dist, Angle, n
  REAL*8, DIMENSION(2)       :: ang_OHO

  n = 0
  index_O2c = 0

  DO i = 1,natoms

    !Selecting only O2C from TiO2
    IF ( ind_atom(i) == 21 )  THEN

      index_O2c = index_O2c + 1

      DO j = 1,natoms

        !Selecting only O from H2O
        IF (ind_atom(j) == 3) THEN

          ind_H = 0
          d_OOw = Dist(pos(i,:), pos(j,:))
          
          !Definition of H-bond - Luzar and Chandler
          IF ( d_OOw < 3.5 ) THEN

            CALL get_H(j, ind_H)

            IF ( (ind_H(1) /= 0) .and. (ind_H(2) /= 0 )) THEN

              ang_OHO(1) = Angle(pos(j,:), pos(i,:), pos(ind_H(1),:))
              ang_OHO(2) = Angle(pos(j,:), pos(i,:), pos(ind_H(2),:))

              !Definition of H-bond - Luzar and Chandler
              !Need to know which H is H-bound to the O2c
              IF ( ang_OHO(1) > 0.8660 ) THEN

                !Assigning the H index to the h-bonded O2C
                IF ( hbond_corr(index_O2C, frame, 1) == 0 ) THEN
                  hbond_corr(index_O2C, frame, 1) = ind_H(1)
                ELSE
                  hbond_corr(index_O2C, frame, 2) = ind_H(1)
                END IF
 
              ELSEIF ( ang_OHO(2) > 0.8660 ) THEN

                !Very ugly, but works
                IF ( hbond_corr(index_O2C, frame, 1) == 0 ) THEN
                  hbond_corr(index_O2C, frame, 1) = ind_H(2)
                ELSE
                  hbond_corr(index_O2C, frame, 2) = ind_H(2)
                END IF

              END IF !d_OOw
 
            END IF !ind_H

          END IF !d_OOw

        END IF !ind_atom

      END DO !j

    END IF !ind_atom

  END DO !i 
 
  !WRITE(2,*), hbond_corr(11,frame,:)


END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE COARSE_GRAIN
  !Make a coase grained function to determine if a O2c is H bonded
  !to a same site (1) or not (0)
  !Recycling hbond_corr array!!!
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i_O2c, j_frame, k_H
  INTEGER                    :: new(2), old(2)
  
  DO i_O2c = 1,neq

    DO j_frame = 2, nframes-nequil
      
      old = hbond_corr(i_O2c, j_frame-1, :)
      new = hbond_corr(i_O2c, j_frame, :)

      !Checking if old H-bond is still unbroken
      !IF ( ( any(new == old) ) .or. any(new == (/ old(2),old(1) /)) ) THEN
      IF ( (( new(1) == old(1))  .and. ( new(1)*old(1) .ne. 0 )) .or. &
           (( new(1) == old(2))  .and. ( new(1)*old(2) .ne. 0 )) .or. &
           (( new(2) == old(1))  .and. ( new(2)*old(1) .ne. 0 )) .or. &
           (( new(2) == old(2))  .and. ( new(2)*old(2) .ne. 0 )) ) THEN

        hbond_corr(i_O2c, j_frame, 1) = 1
        hbond_corr(i_O2c, j_frame-1, 1) = 1

      ELSE
 
        hbond_corr(i_O2c, j_frame, 1) = 0 

      END IF !new, old

    END DO

    hbond_corr(i_O2c, 1, 1) = 1

  END DO

  !print *, hbond_corr(:, :, 1)

END SUBROUTINE COARSE_GRAIN

SUBROUTINE CORR_FUNCT
  !Make the h-bond correlation function from the coarse grained function
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i,j,k, tframe
  INTEGER*8                  :: cf(nframes-nequil)
  
  cf = 0
  tframe = nframes-nequil

  DO i=1,tframe
    DO j=1,tframe-i+1
    
      DO k=1,neq
        cf(i) = cf(i) + hbond_corr(k,j,1)*hbond_corr(k,j+i,1)
      END DO

    END DO !j

    WRITE(2,*), i, DBLE(cf(i)) / ( DBLE(j * cf(1) ) / DBLE(tframe) )  

  END DO !i

END SUBROUTINE CORR_FUNCT

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters
  IMPLICIT NONE
  INTEGER, DIMENSION(2)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1

  DO WHILE (( i <= natoms) .and. ( k <= 2))

    !Should be the index for H
    IF ( ind_atom(i) == 4 ) THEN

      IF ( Dist(pos(ind_O,:),pos(i,:)) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

END SUBROUTINE get_H

REAL*8 FUNCTION Dist(a, b)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(4*i-3) ) * box(4*i-3)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

REAL*8 FUNCTION Dist2(a, b)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(4*i-3) ) * box(4*i-3)

  END DO

  Dist2 = SUM(xyz*xyz)

END FUNCTION Dist2

REAL*8 FUNCTION Angle(C, a1, a2)
  !Angle between C-a1 and C-a2
  !Return value in cos(Angle)
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: C, a1, a2
  REAL*8                     :: Dist2, Dist

  Angle = dot_product((C - a1),(C - a2)) / ( norm2(C - a1) * norm2(C - a2) )
  !Angle = Dist2(C,a1) + Dist2(C,a2) - Dist2(a1,a2)
  !Angle = Angle / ( 2 * Dist(C,a1) * Dist(C,a2) )

END FUNCTION Angle
