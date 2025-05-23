!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE wannier
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: wfc(:,:) 
  INTEGER, ALLOCATABLE                     :: nwf_atom(:)
  INTEGER                                  :: nwfc
END MODULE wannier

PROGRAM CONVERT
  USE parameters, ONLY : nframes, nequil
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL READ_CELL
    CALL READ_WFC
    CALL ORDERED_WFC
    print *,frame
  END DO

  CLOSE(1);CLOSE(2);CLOSE(3)
  CLOSE(4)

END PROGRAM CONVERT

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, pos
  USE wannier, ONLY : wfc, nwfc
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos
  CHARACTER(100)             :: index_file, file_cel
  CHARACTER(100)             :: file_wfc

  CALL getarg(1, file_pos)
  CALL getarg(2, file_cel)
  CALL getarg(3, file_wfc)
  CALL getarg(4, index_file)

  !First 2 lines of index file are number of atoms and box
  CALL READ_INDEX_WFC( index_file ) 
  
  ALLOCATE(pos(natoms,3)); pos=0
  ALLOCATE(wfc(nwfc,3)); wfc=0

  OPEN(unit  =  1,file  =  file_pos)
  OPEN(unit  =  2,file  =  file_cel)
  OPEN(unit  =  3,file  =  file_wfc)
  OPEN(unit  =  4,file  =  "wannier.raw")
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters, ONLY : pos, natoms
  IMPLICIT NONE
  INTEGER                    :: i,j

  READ(1,*)

  DO i = 1,natoms
    READ(1, *), pos(i,:)
  END DO

  !Bohr to angstrom
  pos = pos * 0.529177

END SUBROUTINE

SUBROUTINE READ_CELL
  USE parameters, ONLY : box
  IMPLICIT NONE
  INTEGER                    :: i

  READ(2,*)

  DO i = 1,3
    READ(2, *) box(i,1), box(i,2), box(i,3)
  END DO

  !Bohr to angstrom
  box = box * 0.529177

END SUBROUTINE

SUBROUTINE READ_WFC 
  USE wannier, ONLY : nwfc, wfc
  IMPLICIT NONE
  INTEGER          :: iwf

  READ(3,*)

  DO iwf=1,nwfc
    READ(3,*) wfc(iwf,1), wfc(iwf,2), wfc(iwf,3)
  END DO

  !Bohr to angstrom
  wfc = wfc * 0.529177

END SUBROUTINE READ_WFC

SUBROUTINE APPLY_PBC(x_in,x_out)
  USE parameters, ONLY : box 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)        :: x_in(3)
  DOUBLE PRECISION                    :: x_out(3)
  INTEGER                             :: ipol

  DO ipol = 1,3

    x_out(ipol) = x_in(ipol) - nint(x_in(ipol)/box(ipol,ipol))*box(ipol,ipol) 
    x_out(ipol) = x_out(ipol) + box(ipol,ipol)/2.d0

  END DO

END SUBROUTINE APPLY_PBC 

SUBROUTINE ORDERED_WFC
  USE parameters, ONLY : natoms, pos, box
  USE wannier, ONLY : nwfc, wfc, nwf_atom
  IMPLICIT NONE
  DOUBLE PRECISION           :: wfctmp(3,2*nwfc), wfc_atom(3,nwfc), d(4)
  INTEGER                    :: i, iat
  INTEGER                    :: iwf

  iwf=1

  DO iat = 1,natoms

    IF ( nwf_atom(iat)>0 ) THEN 

      CALL NEAREST_WFCs( iat, wfc_atom, nwf_atom(iat) )

      DO i=1,nwf_atom(iat)
        wfctmp(:,iwf)=wfc_atom(:,i) 
        iwf=iwf+1
      END DO

    END IF

  END DO

  CALL Print_Wannier_Centers(wfctmp,iwf-1)

END SUBROUTINE ORDERED_WFC

SUBROUTINE PRINT_Wannier_Centers(WFC,n)
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: n 
  REAL*8, INTENT(IN)         :: wfc(3,n)
  INTEGER                    :: i
  WRITE(4, fmt='(*(E16.7E4,3X))'), ( wfc(1,i), wfc(2,i), wfc(3,i), i=1,n )
END SUBROUTINE

SUBROUTINE NEAREST_WFCs( iat, wfcs, n )
  !Sum the coordinates of the WFCs closest to atom with index "iat"
  USE parameters,ONLY : pos, natoms
  USE wannier
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat, n
  INTEGER                    :: ind_wf(nwfc)
  INTEGER                    :: iwf
  DOUBLE PRECISION           :: d(4,nwfc), dist(nwfc), wfcs(3,n)

  dist=0.d0 
  DO iwf = 1,nwfc
    CALL DISTANCE_VECTOR( pos(iat,:), wfc(iwf,:), d(:,iwf) )
    ind_wf(iwf)=iwf
    dist(iwf)=d(4,iwf)
  END DO

  CALL Bubble_Sort(dist,ind_wf,n,nwfc)

  DO iwf=1,n
    wfcs(:,iwf)=d(1:3,ind_wf(iwf))
  END DO

END SUBROUTINE NEAREST_WFCs

SUBROUTINE READ_INDEX_WFC ( index_file ) 
  USE wannier,  ONLY : nwfc, nwf_atom
  USE parameters, ONLY : nlayers, nframes, nequil,&
                         natoms, layers, box, ind_atom
  IMPLICIT NONE
  CHARACTER(100), INTENT(IN) :: index_file
  INTEGER                    :: i,j
  
  OPEN(20, file=index_file)

  READ (20, *) natoms, nframes, nequil, nwfc, nlayers

  IF ( nlayers > 1 ) THEN
    ALLOCATE( layers(nlayers+1) )
    READ(20, *) layers
  END IF

  ALLOCATE( ind_atom(natoms) ) 
  ALLOCATE( nwf_atom(natoms) )

  DO i = 1, natoms
    READ(20,*) ind_atom(i), nwf_atom(i)
  END DO

  CLOSE(20)

END SUBROUTINE READ_INDEX_WFC

SUBROUTINE Bubble_Sort(a,ind,n,ntot)
  !Sort array a in ascending order. Ind is an index array
  !a will contain the smallest n entries sorted from 1 to n
  INTEGER, INTENT(in) :: n, ntot
  REAL*8, INTENT(inout), DIMENSION(ntot) :: a
  INTEGER, INTENT(inout), DIMENSION(ntot) :: ind
  REAL*8 :: temp
  INTEGER :: i, j, itemp
  LOGICAL :: swapped
 
  !DO j = SIZE(a)-1, 1, -1
  DO j = 1, n 
    swapped = .FALSE.
    !DO i = 1, j
    DO i = ntot-1, j, -1
      IF (a(i) > a(i+1)) THEN
        temp = a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        itemp = ind(i)
        ind(i) = ind(i+1)
        ind(i+1) = itemp
        swapped = .TRUE.
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END SUBROUTINE Bubble_Sort
