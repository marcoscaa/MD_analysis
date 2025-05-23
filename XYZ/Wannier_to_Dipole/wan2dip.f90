!This codes reads a EXTXYZ file containing Wannier functions centers,
!and writes an EXTXYZ file containing atomic positions plus 3 additional columns.
!These additional columns are the dipole for each dipole center.
!The current version of the code only works for water, assigning 4 
!Wannier function centers for each O atom. This can be generalized to 
!any atoms by improving the current code.

PROGRAM CONVERT
 USE, intrinsic :: iso_fortran_env, Only : iostat_end
 INTEGER :: iostat

 CALL INITIALIZE

 DO 
   CALL READ_EXTXYZ_IO (iostat)
   IF (iostat == iostat_end ) THEN
     EXIT
   END IF
   CALL ASSIGN_WFC_TO_ATOMS 
   CALL COMPUTE_MOLECULAR_DIPOLE
   CALL WRITE_EXTXYZ_DIPOLE
 END DO

 CALL FINALIZE

END PROGRAM

MODULE wannier
  IMPLICIT NONE 
  REAL*8, ALLOCATABLE           :: dipole(:,:)
END MODULE wannier

SUBROUTINE INITIALIZE
  IMPLICIT NONE
  CHARACTER(100)             :: pos_file
  INTEGER :: i

  CALL getarg(1, pos_file)

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = 'pos_dipole.xyz')

END SUBROUTINE INITIALIZE

SUBROUTINE ASSIGN_WFC_TO_ATOMS
  USE parameters, ONLY : natoms, pos, box, atype
  USE wannier, ONLY : dipole 
  IMPLICIT NONE
  REAL*8                     :: wfc_atom(3,4), d(4)
  INTEGER                    :: i, iat

  if (allocated(dipole)) deallocate(dipole)
  allocate(dipole(3,natoms))
  dipole=0.0

  DO iat = 1,natoms

    IF ( trim(atype(iat))=="O" ) THEN 

      CALL NEAREST_WFCs( iat, wfc_atom, 4 )

      DO i=1,4
        dipole(:,iat) = dipole(:,iat) + wfc_atom(:,i)
      END DO

      !Compute the electric dipole
      dipole(:,iat) = -2.0*dipole(:,iat)

    END IF

  END DO

END SUBROUTINE ASSIGN_WFC_TO_ATOMS

SUBROUTINE COMPUTE_MOLECULAR_DIPOLE
  USE parameters, ONLY : natoms, pos, box, atype, OH_index
  USE wannier, ONLY : dipole 
  IMPLICIT NONE
  REAL*8                     :: d(4) 
  INTEGER                    :: i, iat

  if (allocated(OH_index)) deallocate(OH_index)
  allocate(OH_index(natoms))
  !Do a Voronoi tesselation and assign each H to its 
  !closest O atom
  CALL IDENTIFY_OH_GROUPS

  DO iat = 1,natoms

    IF ( trim(atype(iat))=="H" ) THEN 

      !Compute the ionic dipole 
      CALL DISTANCE_VECTOR_IND(iat,OH_index(iat),d)
      dipole(:,OH_index(iat)) = dipole(:,OH_index(iat)) + d(1:3)

    END IF

  END DO

  !Convert to Debye
  dipole = dipole / 0.2081943

END SUBROUTINE COMPUTE_MOLECULAR_DIPOLE

SUBROUTINE NEAREST_WFCs( iat, wfcs, n )
  !Get the coordinates of the n WFCs closest to atom with index "iat"
  !Assuming the WFCs in the XYZ file have atom type "X"
  USE parameters,ONLY : pos, natoms, box, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat, n
  INTEGER                    :: ind_wf(n), ind_0(n)
  INTEGER                    :: iwf, jwf
  REAL*8                     :: d(4), dist(n), wfcs(3,n)
  REAL*8                     :: wfcs_(3,n)

  !Create index array for Bubble sort
  DO iwf = 1,n
    ind_0(iwf)=iwf
  END DO

  dist=100.d0 
  DO iwf = 1,natoms
    IF (trim(atype(iwf))=="X") THEN
      CALL DISTANCE_VECTOR( pos(:,iwf), pos(:,iat), d )

      IF (d(4)<dist(n)) THEN
        dist(n)=d(4)
        wfcs(:,n)=d(1:3)
        ind_wf = ind_0
        CALL Bubble_Sort(dist,ind_wf,n,n)
        wfcs_=wfcs
        DO jwf = 1,n
          wfcs(:,jwf) = wfcs_(:,ind_wf(jwf))
        END DO
      END IF
    END IF
  END DO

END SUBROUTINE NEAREST_WFCs

SUBROUTINE WRITE_EXTXYZ_DIPOLE
  USE parameters, only : pos, natoms, box, atype
  USE wannier, only: dipole 
  IMPLICIT NONE
  INTEGER                    :: iat, nwf, natoms_

  !First determine the number of atoms (no WFCs)
  nwf=0
  DO iat = 1,natoms
    IF (trim(atype(iat))=="X") nwf = nwf + 1
  END DO

  natoms_ = natoms - nwf

  WRITE(2,*) natoms_
  WRITE(2,fmt="(A10,8(F8.5,1X),F8.5,A1,1X,A41)") 'Lattice="', &
    box(1,1), box(1,2), box(1,3), &
    box(2,1), box(2,2), box(2,3), &
    box(3,1), box(3,2), box(3,3), '"', &
    'Properties=species:S:1:pos:R:3:dipole:R:3'

  DO iat = 1,natoms
    IF (trim(atype(iat)).ne."X") THEN
      WRITE(2, fmt="(A3,6(2X,F12.8))") & 
        atype(iat), pos(1,iat), pos(2,iat), pos(3,iat), &
        dipole(1,iat), dipole(2,iat), dipole(3,iat)
    END IF
  END DO

END SUBROUTINE WRITE_EXTXYZ_DIPOLE

SUBROUTINE FINALIZE
  IMPLICIT NONE
  CLOSE(1);CLOSE(2)
END SUBROUTINE FINALIZE
