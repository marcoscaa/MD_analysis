!Set of subroutines used to analyse the lipid trajectory

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: pos(:,:), vel(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: ekin(:,:)
  DOUBLE PRECISION                         :: box(9)
  INTEGER*8, ALLOCATABLE                   :: ind_atom(:), counter(:,:)
  INTEGER,PARAMETER                        :: nbins = 250, iprint=10, sep=10
  INTEGER                                  :: natoms, nframes 
  INTEGER                                  :: nequil, index_equil 
  INTEGER                                  :: nsample ! nframes without duplication 

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  !Remove equilibration part
  CALL REMOVE_EQUIL 

  DO frame = 1,nframes-nequil
    CALL READ_FRAME 
    CALL EQUIPARTITION 
  END DO

  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, file_vel, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, file_vel)
  CALL getarg(3, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  READ(1, *), box

  ALLOCATE( ind_atom(natoms) )

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO

  !For exclusion of duplicated lines
  index_equil=0
  !Upper bond estimate - will be overwritten
  nsample = (nframes-nequil) / sep
  CLOSE(1)

  ALLOCATE(pos(3,nat)) ! coordinates
  ALLOCATE(vel(3,nat)) ! velocities 
  ALLOCATE(ekin(nbins,4)) ! 4 is the number of atomtypes
  ALLOCATE(counter(nbins,4)) ! 4 is the number of atomtypes

  OPEN(unit = 1,file = file_pos)
  OPEN(unit = 2,file = file_vel)
  
END SUBROUTINE INITIALIZE

SUBROUTINE REMOVE_EQUIL
  !Subroutine reads gro files
  USE parameters, ONLY : nequil, natoms
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1, nequil*(natoms+1)

    READ(1,*)
    READ(2,*)

  END DO

END SUBROUTINE

SUBROUTINE READ_FRAME 
  !Subroutine reads gro files
  USE parameters
  IMPLICIT NONE
  REAL*8                     :: junk
  INTEGER                    :: i

  !One could test here if the files match
  READ(1,*)
  READ(2,*)

  DO i = 1,natoms

    READ(1, fmt='(3(3X,E22.14))'), pos(1,i), pos(2,i), pos(3,i)
    READ(2, fmt='(3(3X,E22.14))'), vel(1,i), vel(2,i), vel(3,i)

  END DO

  !Bohr to angstrom
  pos = pos * 0.529177

END SUBROUTINE

SUBROUTINE EQUIPARTITION
  !Compute the instantaneous temperature based on classical equipartition
  !theorem
  USE parameters
  IMPLICIT NONE
  DOUBLE PRECISION           :: kin, z
  INTEGER                    :: iat, ipol

  DO iat = 1, nat

    kin = 0.d0

    DO ipol = 1,3

      kin = kin + vel(ipol,iat)*vel(ipol,iat)

    END DO

    z = pos(3,iat) - nint( pos(3,iat) / box(9) ) * box(9)
    z = z + box(9) / 2.d0
    ind_layer = int( dble(nbins) * z / box(9) ) + 1
    ekin(ind_layer,ind_atom(iat)) = ekin(ind_layer,ind_atom(iat)) + kin 
    counter(ind_layer,ind_atom(iat)) = counter(ind_layer,ind_atom(iat)) + 1

  END DO

END SUBROUTINE EQUIPARTIION

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE parameters
  IMPLICIT NONE
  DOUBLE PRECISION           :: zbin
  INTEGER                    :: i, j

  zbin = box(9) / dble(nbins)

  DO i = 1,nbins

    write(*,fmt = '(E11.4,3X,4(E12.5,3X))'), (dble(i)-0.5d0)*zbin, ekin(i,:) / dble(counter(i,:))

  ENDDO

END SUBROUTINE PRINT_RESULTS

