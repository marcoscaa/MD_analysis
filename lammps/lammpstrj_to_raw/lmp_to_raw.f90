!Converts Lammps atom file to DeepMD raw format

PROGRAM LMP_to_RAW
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    !CALL ORDERED_POS
    IF (MOD(frame,nskip)==0) THEN
      !CALL SET_CENTER_OF_MASS_TO_ZERO
      CALL PRINT_RAW
    END IF
  END DO

  CALL PRINT_TYPE_RAW

  CLOSE(1);CLOSE(2);CLOSE(3)

END PROGRAM LMP_TO_RAW

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, box,&
                         atype, pos, nskip
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nskip

  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))

  CLOSE(1)
  OPEN(unit  =  1,file  =  pos_file)
  OPEN(unit = 2,file = "coord.raw")
  OPEN(unit = 3,file = "box.raw")
 
END SUBROUTINE INITIALIZE

SUBROUTINE PRINT_RAW
  USE parameters, ONLY : natoms, box, pos
  IMPLICIT NONE
  INTEGER                    :: i, j

  WRITE(2, fmt="(*(F12.7,3X))") ( (pos(j,i), j=1,3), i=1,natoms )
  WRITE(3, fmt="(9(F12.7,3X))") box(1), 0., 0., 0., box(2), 0., 0., 0., box(3) 

END SUBROUTINE

SUBROUTINE ORDERED_POS
  USE parameters, ONLY : natoms, pos, box, atype
  IMPLICIT NONE
  DOUBLE PRECISION           :: postmp(3,natoms), d(4)
  INTEGER                    :: i, iat, ind_H(2)
  INTEGER                    :: iatnew

  iatnew=1

  DO iat =1,natoms

    IF ( atype(iat).eq.1 ) THEN !water oxygen

      CALL get_H( iat, ind_H )

      CALL APPLY_PBC(pos(:,iat), postmp(:,iatnew))
      CALL DISTANCE_VECTOR( pos(:,iat), pos(:,ind_H(1)), d )
      postmp(:,iatnew+1) = postmp(:,iatnew) + d(1:3)
      CALL DISTANCE_VECTOR( pos(:,iat), pos(:,ind_H(2)), d )
      postmp(:,iatnew+2) = postmp(:,iatnew) + d(1:3)

      if (d(4)>1.3) PRINT *, 'Too large OH bond', iat, d(4)

      iatnew = iatnew + 3

    END IF

  END DO

  pos = postmp

END SUBROUTINE ORDERED_POS

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : natoms, cutoff_OH, atype
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                :: ind_H
  INTEGER                              :: i, k, ind_O
  DOUBLE PRECISION                     :: Dist, dist_OH, mindist(2)

  mindist = 100.d0
  ind_H = 0

  DO i = 1,natoms

    !Should be the index for H
    IF ( atype(i) == 2 ) THEN

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

  !For consitency: higher index comes first
  IF ( ind_H(1) < ind_H(2) ) THEN
    i = ind_H(2)
    ind_H(2) = ind_H(1)
    ind_H(1) = i
  END IF

END SUBROUTINE get_H

SUBROUTINE PRINT_TYPE_RAW
  USE parameters, only : atype, natoms
  IMPLICIT NONE
  INTEGER :: iat

  OPEN(unit=4, file='type.raw')
  WRITE(4,fmt='(*(I3,1X))') (atype(iat), iat=1,natoms)
  CLOSE(4)

END SUBROUTINE PRINT_TYPE_RAW
