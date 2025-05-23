PROGRAM CONVERT
 USE parameters, only: nframes, nequil
 INTEGER :: frame

 CALL INITIALIZE

 DO frame=1,nframes
   CALL READ_RAW_POS_BOX
   CALL READ_RAW_WFC
   !CALL ORDERED_WFC
   CALL HBOND_WFC
 END DO

 CALL FINALIZE

END PROGRAM

MODULE wannier
  IMPLICIT NONE 
  REAL*8, ALLOCATABLE                        :: wfc(:,:) 
  INTEGER, ALLOCATABLE                     :: nwf_atom(:)
  INTEGER                                  :: nwfc
  INTEGER                                  :: Otype, Htype
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
  READ(1,*) Otype, Htype
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(wfc(3,nwfc)) ! Wannier centers
  ALLOCATE(atype(natoms))

  !Read atype from type.raw
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  CLOSE(1)

  OPEN(unit = 1,file = 'coord.raw')
  OPEN(unit = 2,file = 'box.raw')
  OPEN(unit = 4,file = 'wannier.raw')
  OPEN(unit = 5,file = 'wannier_lone_distances.dat')
  OPEN(unit = 6,file = 'wannier_bond_distances.dat')

  WRITE(5, *) "# Dist_O_Hw Dist_O_WFC Dist_Hw_WFC "
  WRITE(6, *) "# Dist_O_Hw Dist_O_WFC Dist_Occ_Hw "

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
  REAL*8                     :: wfctmp(3,nwfc), wfc_atom(3,nwfc)
  INTEGER                    :: i, iat
  INTEGER                    :: iwf

  iwf=1
  wfctmp=0.0

  DO iat = 1,natoms

    IF ( nwf_atom(atype(iat))>0 ) THEN 

      CALL NEAREST_WFCs( iat, wfc_atom, nwf_atom(atype(iat)) )

      DO i=1,nwf_atom(atype(iat))
        wfctmp(:,iwf)=wfc_atom(:,i) 
        iwf=iwf+1
      END DO

    END IF

  END DO

  wfc = wfctmp

END SUBROUTINE ORDERED_WFC

SUBROUTINE NEAREST_WFCs( iat, wfcs, n )
  !Sum the coordinates of the WFCs closest to atom with index "iat"
  USE parameters,ONLY : pos, natoms, box, boxinv
  USE wannier
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat, n
  INTEGER                    :: ind_wf(nwfc)
  INTEGER                    :: iwf
  REAL*8                     :: d(4,nwfc), dist(nwfc), wfcs(3,n)

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

SUBROUTINE HBOND_WFC
  USE parameters, only: natoms, pos, atype, box, boxinv 
  USE wannier, only: wfc, nwf_atom, nwfc, Otype
  IMPLICIT NONE
  INTEGER                    :: iat, iwf
  INTEGER                    :: ind_hw, ind_O, ind(4)
  INTEGER                    :: get_H_Hbonded_to_wfc
  INTEGER                    :: get_O_Hbonded_to_H, ind_Oacc
  REAL*8                     :: O_wfc(3,nwfc)
  REAL*8                     :: coord_wfc(3) 
  REAL*8                     :: Distxyz
  REAL*8                     :: d_hw_wfc, d_o_wfc, d_o_hw, d_Occ_Hw 

  ind_O = 0

  DO iat = 1, natoms

      !This analysis will focus on water only. Oxygen is the
      !dipole center of water.
      IF (atype(iat)==Otype) THEN
          O_wfc=0.0

          !Find the 4 WFCs closest to O. This will return an 
          !array of size (4,3), sorted by O-W distance
          CALL NEAREST_WFCs( iat, O_wfc, nwf_atom(atype(iat)) )

          !Loop over the lone pairs of WFCs
          DO iwf = 1,2
              !O_wfc has the origin on the O atom. This needs to be shifted
              coord_wfc = pos(:,iat) + O_wfc(:,iwf)
              !Find the index of H H-bonded to O (and the WFC)
              ind_hw = get_H_Hbonded_to_wfc(coord_wfc,pos(:,iat))
              if (ind_hw .ne. 0) then
                  !3 distances used for future analyses
                  d_hw_wfc = Distxyz(pos(:,ind_hw),coord_wfc,box,boxinv)
                  d_o_wfc = norm2(O_wfc(:,iwf))
                  d_o_hw = Distxyz(pos(:,ind_hw),pos(:,iat),box,boxinv)
                  WRITE(5,fmt = "(3(F12.8,3X))") d_o_hw, d_o_wfc, d_hw_wfc
              end if
          END DO

          !Loop over the bond pairs of WFCs
          DO iwf = 3,4
              !O_wfc has the origin on the O atom. This needs to be shifted
              coord_wfc = pos(:,iat) + O_wfc(:,iwf)
              !Find the index of H bonded to O (and the WFC)
              ind_hw = get_H_Hbonded_to_wfc(coord_wfc,pos(:,iat))
              !Find the index of O accepting H-bond from ind_hw
              ind_Oacc = get_O_Hbonded_to_H(pos(:,ind_hw),pos(:,iat))
              if ( (ind_hw .ne. 0) .and. (ind_Oacc .ne. 0) ) then
                  !3 distances used for future analyses
                  d_o_hw = Distxyz(pos(:,iat),pos(:,ind_hw),box,boxinv)
                  d_o_wfc = norm2(O_wfc(:,iwf))
                  d_Occ_Hw = Distxyz(pos(:,ind_hw),pos(:,ind_Oacc),box,boxinv)
                  WRITE(6,fmt = "(3(F12.8,3X))") d_o_hw, d_o_wfc, d_Occ_Hw
              end if
          END DO

      ind_O = ind_O + 1

      END IF
  END DO


END SUBROUTINE HBOND_WFC 

SUBROUTINE Bubble_Sort(a,ind,n,ntot)
  !Sort array a in ascending order. Ind is an index array
  !a will contain the smallest n entries sorted from 1 to n
  INTEGER, INTENT(in) :: n, ntot
  REAL*8, INTENT(inout), DIMENSION(ntot) :: a
  INTEGER, INTENT(inout), DIMENSION(ntot) :: ind
  REAL*8 :: temp
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

INTEGER FUNCTION get_H_Hbonded_to_wfc(coord_wfc,coord_O)
    USE parameters, ONLY: pos, natoms, box, boxinv, atype
    USE wannier, ONLY : Htype
    IMPLICIT NONE
    REAL*8, INTENT(IN)       :: coord_wfc(3), coord_O(3)
    REAL*8                   :: angle_OWH, d_min, d_HW
    REAL*8                   :: Distxyz, Angle, ang_OWH
    INTEGER                  :: ind_min, iH
    LOGICAL                  :: is_hydrogen

    ind_min = 0
    d_min = 100.

    DO iH = 1, natoms
    
        IF (atype(iH)==Htype) THEN

            d_HW = Distxyz(coord_wfc, pos(:,iH), box, boxinv)
           
            IF ( (d_HW < 3.0) .and. (d_HW < d_min) ) THEN
           
                ang_OWH = Angle(coord_O, coord_wfc, pos(:,iH))
           
                !Definition of H-bond based on Wannier function 
                IF (ang_OWH > 0.6427) THEN
                   d_min = d_HW
                   ind_min = iH
                END IF
           
            END IF 

        END IF

    END DO

    get_H_Hbonded_to_wfc = ind_min

END FUNCTION get_H_Hbonded_to_wfc

INTEGER FUNCTION get_O_Hbonded_to_H(coord_H,coord_O)
    USE parameters, ONLY: pos, natoms, box, boxinv, atype
    USE wannier, ONLY : Otype
    IMPLICIT NONE
    REAL*8, INTENT(IN)       :: coord_H(3), coord_O(3)
    REAL*8                   :: angle_OHO, d_min, d_OO
    REAL*8                   :: Distxyz, Angle, ang_OHO
    INTEGER                  :: ind_min, iO

    ind_min = 0
    d_min = 100.

    DO iO = 1, natoms
    
        IF (atype(iO)==Otype) THEN

            d_OO = Distxyz(coord_O, pos(:,iO), box, boxinv)
           
            IF (d_OO < 3.5) THEN
           
                ang_OHO = Angle(coord_O, coord_H, pos(:,iO))
           
                !Definition of H-bond using Luzar-Chandler  
                IF (ang_OHO > 0.8660) THEN
                   d_min = Distxyz(pos(:,iO), coord_H, box, boxinv)
                   ind_min = iO
                END IF
           
            END IF 

        END IF

    END DO

    get_O_Hbonded_to_H = ind_min

END FUNCTION get_O_Hbonded_to_H

SUBROUTINE FINALIZE
  IMPLICIT NONE
  CLOSE(1);CLOSE(2);CLOSE(3)
  CLOSE(5)
END SUBROUTINE FINALIZE
