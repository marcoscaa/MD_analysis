!Computes the average and instantaneous Hbond

MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER                                 :: neq = 12, iprint=10
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE      :: pos
  DOUBLE PRECISION, DIMENSION(9)                     :: box
  DOUBLE PRECISION, PARAMETER                        :: cutoff_OH = 1.3D0 
  INTEGER, DIMENSION(:), ALLOCATABLE                 :: ind_atom
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE             :: hbond_corr
  INTEGER                                            :: natoms, nframes
  INTEGER*8                                          :: n_H, n_H2
  INTEGER                                            :: nequil, index_equil

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame, index_frame

  CALL INITIALIZE

  !Remove equilibration part
  DO frame = 1,nequil
    CALL READ_FRAME (index_frame) 
  END DO

  index_equil = index_frame

  DO frame = nequil+1,nframes
    CALL READ_FRAME (index_frame)
    CALL MAKE_HISTOGRAM (index_frame)
  END DO

  CALL COARSE_GRAIN (index_frame)
  CALL CORR_FUNCT (index_frame)

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
 
  !For exclusion of duplicated lines
  index_equil=0
  CLOSE(1)

  OPEN(unit = 1,file = file_name)
  OPEN(unit = 2,file = "H_Corr.dat")
 

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME (index_frame)
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, index_frame
  REAL(8)                    :: junk

  READ(1,*), index_frame, junk

  ! This will take care of the duplication
  index_frame = index_frame/iprint - index_equil
  print *, index_frame

  DO i = 1,natoms
    READ(1, fmt='(3(3X,E22.14))'), pos(i,1), pos(i,2), pos(i,3)
  END DO

END SUBROUTINE

SUBROUTINE MAKE_HISTOGRAM (frame)
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: i, j, frame
  INTEGER                              :: index_O2c
  INTEGER, DIMENSION(2)                :: ind_H 
  DOUBLE PRECISION                     :: d_OOw, Dist, Angle, n
  DOUBLE PRECISION, DIMENSION(2)       :: ang_OHO

  n = 0
  index_O2c = 0

  DO i = 1,natoms

    !Selecting only O2C from TiO2
    IF ( ind_atom(i) == 31 )  THEN

      index_O2c = index_O2c + 1

      DO j = 1,natoms

        !Selecting only O from H2O
        IF (ind_atom(j) == 4) THEN

          ind_H = 0
          d_OOw = Dist(i,j)
          
          !Definition of H-bond - Luzar and Chandler
          IF ( d_OOw < 3.5 ) THEN

            CALL get_H(j, ind_H)

            IF ( All(ind_H /= 0) ) THEN

              ang_OHO(1) = Angle(pos,j,i,ind_H(1))
              ang_OHO(2) = Angle(pos,j,i,ind_H(2))

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

SUBROUTINE COARSE_GRAIN (tframe)
  !Make a coase grained function to determine if a O2c is H bonded
  !to a same site (1) or not (0). Could be more sofisticated, to check 
  !if a particular h-bond (not any) is still alive
  !Recycling hbond_corr array!!!
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i_O2c, j_frame, k_H, tframe
  INTEGER                    :: new(2)
  
  DO i_O2c = 1,neq

    DO j_frame = 1, tframe
      
      new = hbond_corr(i_O2c, j_frame, :)

      !Checking if old H-bond is still unbroken
      IF ( any(new .ne. 0 ) ) THEN

        hbond_corr(i_O2c, j_frame, 1) = 1

      ELSE
 
        hbond_corr(i_O2c, j_frame, 1) = 0 

      END IF !new, old

    END DO

    hbond_corr(i_O2c, 1, 1) = 1

  END DO

  !print *, hbond_corr(:, :, 1)

END SUBROUTINE COARSE_GRAIN

SUBROUTINE CORR_FUNCT (tframe)
  !Make the h-bond correlation function from the coarse grained function
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i,j,k, tframe
  INTEGER*8                  :: cf(neq), count_hcorr
  DOUBLE PRECISION           :: average_cf, cf0
  
  !Average sqare of the number of h-bonds 
  cf0 = 0.
  DO k=1,neq 
    cf0 = cf0 + sum(hbond_corr(k,:,1)*hbond_corr(k,:,1))
  END DO
  cf0 = cf0 / dble(tframe*neq)

  !Compute the h-bond itself
  DO i=0,tframe

    cf = 0
    count_hcorr = 0

    DO j=1,tframe-i+1
    
      DO k=1,neq

        !Consider only intact h-bonds
        IF ( ALL( hbond_corr(k,j:j+i,1) == 1 ) ) THEN
          cf(k) = cf(k) + hbond_corr(k,j,1)*hbond_corr(k,j+i,1)
        END IF

      END DO

    END DO !j

    average_cf = sum(cf) / dble(j*neq)

    WRITE(2,*), i, average_cf / cf0  

  END DO !i

END SUBROUTINE CORR_FUNCT

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : natoms, cutoff_OH, ind_atom
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                :: ind_H
  INTEGER                              :: i, k, ind_O
  DOUBLE PRECISION                     :: Dist, dist_OH, mindist(2)

  mindist = 100.d0

  DO i = 1,natoms

    !Should be the index for H
    IF ( ind_atom(i) == 2 ) THEN

      dist_OH = Dist(ind_O,i)

      !Select only H's within a threshold
      IF ( dist_OH < cutoff_OH ) THEN
 
        !Very basic sorting of all OH distances within the threshold
        IF ( dist_OH < mindist(1) ) THEN 

          ind_H(2) = ind_H(1)
          ind_H(1) = i
          mindist(2) = mindist(1)
          mindist(1) = dist_OH

        ELSEIF ( dist_OH < mindist(2) ) THEN
 
          ind_H(2) = i
          mindist(2) = dist_OH

        END IF
 
      END IF
 
    END IF

  END DO

END SUBROUTINE get_H

DOUBLE PRECISION FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  DOUBLE PRECISION                     :: xyz(3)
  INTEGER                              :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos(ind1,i) - pos(ind2,i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(4*i-3) ) * box(4*i-3)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

DOUBLE PRECISION FUNCTION Angle(pos,ind1, ind2, ind3)
  !Angle between particles with index ind1-ind2 and ind1-ind3
  !Return value in cos(Angle)
  USE parameters, ONLY : natoms
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)         :: pos(natoms,3)
  INTEGER, INTENT(IN)        :: ind1, ind2, ind3
  DOUBLE PRECISION                     :: v1(3), v2(3)

  v1(:) = pos(ind2,:) - pos(ind1,:)
  v2(:) = pos(ind3,:) - pos(ind1,:)

  Angle = dot_product(v1,v2) / ( norm2(v1) * norm2(v2) )

END FUNCTION Angle
