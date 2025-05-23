!Set of subroutines used to analyse the lipid trajectory

MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: wfc(:,:) 
  INTEGER                                  :: nwfc
END MODULE histogram

PROGRAM SFG
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL_WFC

  DO frame = 1,nframes
    CALL READ_POS
    CALL READ_WFC
    CALL PRINT_POS_WFC 
  END DO

  CLOSE(1);CLOSE(2);CLOSE(3)

END PROGRAM SFG

SUBROUTINE REMOVE_EQUIL_WFC
  USE parameters, ONLY : natoms, nequil, nframes
  USE histogram,  ONLY : nwfc
  IMPLICIT NONE
  INTEGER                    :: iwfc, iat

  DO iwfc = 1, nequil*(nwfc+1)
    READ(2,*)
  END DO

  DO iat = 1, nequil*(natoms+1)
    READ(1,*)
  END DO

  nframes = nframes-nequil

END SUBROUTINE REMOVE_EQUIL_WFC

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : wfc, nwfc 
  USE parameters, ONLY : natoms, nframes, nwater, &
                         pos, coarse, nlayers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file
  CHARACTER(100)             :: wfc_file 

  CALL getarg(1, pos_file)
  CALL getarg(2, wfc_file)
  CALL getarg(3, index_file)

  CALL READ_INDEX_WFC ( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  IF ( nlayers > 1 ) THEN 
    ALLOCATE(coarse(nframes,nwater)); coarse=0
  END IF
  ALLOCATE(wfc(nwfc,3)); wfc=0
  
  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = wfc_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_INDEX_WFC ( index_file ) 
  USE histogram,  ONLY : nwfc
  USE parameters, ONLY : nlayers, nframes, nequil,&
                         natoms, layers, box, ind_atom
  IMPLICIT NONE
  CHARACTER(100), INTENT(IN) :: index_file
  INTEGER                    :: i,j
  
  OPEN(2, file=index_file)

  READ (2, *) natoms, nframes, nequil, nwfc, nlayers
  READ (2, *) ( ( box(j,i), i=1,3 ), j=1,3 )

  IF ( nlayers > 1 ) THEN
    ALLOCATE( layers(nlayers+1) )
    READ(2, *) layers
  END IF

  ALLOCATE( ind_atom(natoms) ) 

  DO i = 1, natoms
    READ(2,*) ind_atom(i)
  END DO

  CLOSE(2)

END SUBROUTINE READ_INDEX_WFC

SUBROUTINE READ_WFC
  USE histogram, ONLY : wfc, nwfc 
  USE parameters,ONLY : natoms, pos
  IMPLICIT NONE
  INTEGER                    :: iwf, ipol

  READ(2,*)

  DO iwf = 1, nwfc
    READ(2, *), ( wfc(iwf,ipol), ipol=1,3 )
  END DO

  !Bohr to angstrom
  wfc = wfc * 0.529177

END SUBROUTINE

SUBROUTINE PRINT_POS_WFC
  USE parameters, ONLY : pos, natoms, box, ind_atom
  USE histogram,  ONLY : wfc, nwfc
  IMPLICIT NONE
  INTEGER                    :: iat, iwfc, ipol
  DOUBLE PRECISION           :: coord(3)
  CHARACTER(5)               :: atom_type

  write(3,*) natoms+nwfc
  write(3,*) 

  DO iat = 1, natoms
    DO ipol = 1,3
      coord(ipol) = pos(iat,ipol) - nint(pos(iat,ipol)/box(ipol,ipol))*box(ipol,ipol)
      coord(ipol) = coord(ipol) + box(ipol,ipol)/2.
    END DO
    write(3,fmt='(A3, 3(3X,E12.5))') atom_type(ind_atom(iat)), coord
  END DO

  DO iwfc = 1, nwfc
    DO ipol = 1,3
      coord(ipol) = wfc(iwfc,ipol) - nint(wfc(iwfc,ipol)/box(ipol,ipol))*box(ipol,ipol)
      coord(ipol) = coord(ipol) + box(ipol,ipol)/2.
    END DO
    write(3,fmt='(A3, 3(3X,E12.5))') 'X  ', coord
  END DO

END SUBROUTINE PRINT_POS_WFC

CHARACTER(5) FUNCTION atom_type(itype)
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: itype

  SELECT CASE (itype)
    CASE (1)
      atom_type='Ti'
    CASE (2)
      atom_type='H'
    CASE (3)
      atom_type='O'
    CASE (4)
      atom_type='O'
  END SELECT
    
END FUNCTION atom_type
