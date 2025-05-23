PROGRAM CONVERT
 USE parameters, only: nframes, nequil
 INTEGER :: frame

 CALL INITIALIZE

 DO frame=1,nframes
   CALL READ_RAW_POS_BOX
   CALL READ_RAW_WFC
   IF (frame>nequil) THEN
     CALL ORDERED_WFC
   END IF
 END DO

 CALL FINALIZE

END PROGRAM

MODULE wannier
  IMPLICIT NONE 
  REAL, ALLOCATABLE            :: wfc(:,:) 
  INTEGER, ALLOCATABLE                     :: nwf_atom(:)
  INTEGER                                  :: nwfc
END MODULE wannier

SUBROUTINE INITIALIZE
  USE parameters, only : natoms, nframes, nequil, pos, &
                         atype, ntype
  USE wannier
  IMPLICIT NONE
  CHARACTER(100)             :: index_file
  INTEGER :: i

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nwfc, ntype 
  ALLOCATE(nwf_atom(ntype))
  READ(1,*) nwf_atom
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(wfc(3,nwfc)) ! Wannier centers
  ALLOCATE(atype(natoms))

  !Read atype from type.raw
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  atype=atype+1
  CLOSE(1)

  OPEN(unit = 1,file = 'coord.raw')
  OPEN(unit = 2,file = 'box.raw')
  OPEN(unit = 4,file = 'wannier.raw')
  OPEN(unit = 5,file = 'wannier_ordered.raw')
  OPEN(unit = 6,file = 'dipole.raw')

END SUBROUTINE INITIALIZE

SUBROUTINE READ_RAW_WFC 
  !Read DPMD raw file
  USE wannier, only : wfc, nwfc
  IMPLICIT NONE
  INTEGER                    :: i, j
 
  !Read Wannier centers
  READ(4,*) ( ( wfc(i, j), i=1,3 ), j=1,nwfc )

END SUBROUTINE READ_RAW_WFC

SUBROUTINE ORDERED_WFC
  USE parameters, ONLY : natoms, pos, box, atype
  USE wannier, ONLY : nwfc, wfc, nwf_atom
  IMPLICIT NONE
  REAL           :: wfctmp(3,2*nwfc), wfc_atom(3,nwfc)
  REAL           :: dipole(3,natoms) 
  INTEGER                    :: i, iat, idip
  INTEGER                    :: iwf

  iwf=1
  idip=1
  dipole=0.0
  wfctmp=0.0

  DO iat = 1,natoms

    IF ( nwf_atom(atype(iat))>0 ) THEN 

      CALL NEAREST_WFCs( iat, wfc_atom, nwf_atom(atype(iat)) )

      DO i=1,nwf_atom(atype(iat))
        wfctmp(:,iwf)=wfc_atom(:,i) 
        dipole(:,idip) = dipole(:,idip) + wfc_atom(:,i)
        iwf=iwf+1
      END DO

      !Compute the Wannier centroid
      dipole(:,idip) = dipole(:,idip) / dble(nwf_atom(atype(iat)))
      idip = idip + 1

    END IF

  END DO

  CALL Print_Wannier_Centers(wfctmp,iwf-1)
  CALL Print_Dipole(dipole,natoms,idip-1)

END SUBROUTINE ORDERED_WFC

SUBROUTINE PRINT_Wannier_Centers(WFC,n)
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: n 
  REAL, INTENT(IN)         :: wfc(3,n)
  INTEGER                    :: i
  WRITE(5, fmt='(*(E16.7E4,3X))'), ( wfc(1,i), wfc(2,i), wfc(3,i), i=1,n )
END SUBROUTINE

SUBROUTINE Print_Dipole(dipole,n,ndip)
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: n, ndip
  REAL, INTENT(IN)         :: dipole(3,n)
  INTEGER                    :: i
  WRITE(6, fmt='(*(E16.7E4,3X))'), ( dipole(1,i), dipole(2,i), dipole(3,i), i=1,ndip )
END SUBROUTINE

SUBROUTINE NEAREST_WFCs( iat, wfcs, n )
  !Sum the coordinates of the WFCs closest to atom with index "iat"
  USE parameters,ONLY : pos, natoms, box
  USE wannier
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat, n
  INTEGER                    :: ind_wf(nwfc)
  INTEGER                    :: iwf
  REAL           :: d(4,nwfc), dist(nwfc), wfcs(3,n)

  dist=0.d0 
  DO iwf = 1,nwfc
    CALL DISTANCE_VECTOR( pos(:,iat), wfc(:,iwf), d(:,iwf) )
    ind_wf(iwf)=iwf
    dist(iwf)=d(4,iwf)
  END DO

  CALL Bubble_Sort(dist,ind_wf,n,nwfc)

  DO iwf=1,n
    wfcs(:,iwf)=d(1:3,ind_wf(iwf))
  END DO

END SUBROUTINE NEAREST_WFCs

SUBROUTINE Bubble_Sort(a,ind,n,ntot)
  !Sort array a in ascending order. Ind is an index array
  !a will contain the smallest n entries sorted from 1 to n
  INTEGER, INTENT(in) :: n, ntot
  REAL, INTENT(inout), DIMENSION(ntot) :: a
  INTEGER, INTENT(inout), DIMENSION(ntot) :: ind
  REAL :: temp
  INTEGER :: i, j, itemp
  LOGICAL :: swapped
 
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

SUBROUTINE FINALIZE
  IMPLICIT NONE
  CLOSE(1);CLOSE(2);CLOSE(3)
  CLOSE(5);CLOSE(6)
END SUBROUTINE FINALIZE
