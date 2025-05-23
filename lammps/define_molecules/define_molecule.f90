
PROGRAM ZDF
  USE parameters, ONLY : nequil, nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL DEFINE_MOLECULE
    CALL PRINT_MOL_INDEX
  END DO

  CLOSE(1);CLOSE(10)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, nattype, rcut, molid, &
                         index_bonded 
  IMPLICIT NONE
  INTEGER                    :: i,j
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, nattype 
  ALLOCATE(rcut(nattype,nattype))
  READ(1, *) ( ( rcut(i,j), j=i,nattype ), i=1,nattype )
  
  do i=1,nattype
    do j=i+1,nattype
      rcut(j,i)=rcut(i,j)
    end do
  end do

  ALLOCATE(pos(3,natoms))
  ALLOCATE(molid(natoms))
  ALLOCATE(index_bonded(natoms,10))
  ALLOCATE(atype(natoms))

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
  OPEN(unit = 10,file = 'mol_index.dat')
 
END SUBROUTINE INITIALIZE

SUBROUTINE PRINT_MOL_INDEX
  USE parameters, ONLY : natoms, atype, molid
  IMPLICIT NONE
  INTEGER :: iat

  DO iat=1,natoms
    WRITE(10,fmt='(3(I10,3X))') iat, atype(iat), molid(iat)
  END DO
  WRITE(10,*) '---'

END SUBROUTINE PRINT_MOL_INDEX
