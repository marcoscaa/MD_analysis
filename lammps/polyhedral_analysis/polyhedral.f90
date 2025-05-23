!Computes the average coordination number of 
!atoms along the axis normal to the surface (z)

MODULE histogram
  IMPLICIT NONE 
  INTEGER                            :: atype1, atype2, ntype1
  INTEGER, PARAMETER                 :: maxcoord=12 !Max coordination number
  INTEGER, PARAMETER                 :: maxneigh=26 !Max polyhedral neighbors
  INTEGER, ALLOCATABLE               :: list_coord(:,:), list_neigh(:,:) 
  INTEGER                            :: hist_poly(3)
  REAL*8                             :: dcut_coord, dcut_neigh  
END MODULE histogram

PROGRAM Polyhedral
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    IF (MOD(frame,stride)==0) CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Polyhedral

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, pos, atype, stride, &
                         nframes, nequil, box
  USE histogram
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, ntype1, nframes, nequil, stride
  READ(1,*) atype1, atype2, dcut_coord, dcut_neigh 
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(list_coord(maxcoord,ntype1));list_coord=0 
  ALLOCATE(list_neigh(maxneigh,ntype1));list_neigh=0
  hist_poly=0

  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE parameters, ONLY : natoms, atype
  USE histogram, ONLY : ntype1, maxcoord, hist_poly, atype1, atype2, &
                        dcut_coord, dcut_neigh, maxneigh
  IMPLICIT NONE
  INTEGER                     :: iat, icenter, ineigh, nshared
  INTEGER                     :: coord_ind(maxcoord,ntype1)
  INTEGER                     :: adj_ind(maxneigh,ntype1)
  INTEGER                     :: center_ind(natoms)

  coord_ind=0
  adj_ind=0
  icenter=1
  center_ind = 0

  DO iat = 1, natoms
    IF (atype(iat)==atype1) THEN
      CALL LIST_COORDINATE(iat,atype2,maxcoord,dcut_coord,coord_ind(:,icenter))
      CALL LIST_COORDINATE(iat,atype1,maxneigh,dcut_neigh,adj_ind(:,icenter))
      center_ind(iat) = icenter
      icenter=icenter+1
    END IF
  END DO

  DO icenter = 1,ntype1
    DO ineigh = 1, maxneigh
      IF ( adj_ind(ineigh,icenter) .ne. 0 ) THEN
        adj_ind(ineigh,icenter) = center_ind(adj_ind(ineigh,icenter))
      END IF
    END DO
  END DO

  DO icenter = 1, ntype1
    DO ineigh = 1, maxneigh
      IF ( adj_ind(ineigh,icenter) .ne. 0 ) THEN
        CALL COUNT_SHARED_VERTICES(icenter,ineigh,coord_ind,adj_ind,nshared)
        hist_poly(nshared) = hist_poly(nshared) + 1
      END IF
    END DO
  END DO

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE LIST_COORDINATE( iat, itype, maxn, d_cut, list_coord )
  !Compute the coordination number of atom with index iat
  USE parameters, ONLY : natoms, pos, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat, itype, maxn
  REAL*8,  INTENT(IN)        :: d_cut
  INTEGER                    :: list_coord(maxn)
  INTEGER                    :: ic, n_coord, iat2
  REAL*8                     :: d, Dist

  n_coord = 0
  ic = 1
  DO iat2 = 1, natoms
    IF ( (iat2.ne.iat).and.(atype(iat2).eq.itype) ) THEN
      d = Dist( iat, iat2 )
      IF ( d < d_cut )  THEN
        list_coord(ic) = iat2
        ic = ic + 1
      END IF
    END IF 
  END DO
 
END SUBROUTINE LIST_COORDINATE 

SUBROUTINE COUNT_SHARED_VERTICES(icenter,ineigh,coord_ind,adj_ind,nshared)
  !Compute the number of common ligands between icenter and its neighbohr ineigh
  USE histogram, ONLY : maxcoord, maxneigh, ntype1
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: icenter, ineigh
  INTEGER, INTENT(IN)        :: coord_ind(maxcoord,ntype1)
  INTEGER, INTENT(IN)        :: adj_ind(maxneigh,ntype1) 
  INTEGER                    :: nshared
  INTEGER                    :: ic1, ic2

  nshared = 0

  DO ic1 = 1, maxcoord
    IF ( coord_ind(ic1,icenter) .ne. 0 ) THEN
      DO ic2 = 1, maxcoord
        IF (coord_ind(ic1,icenter).eq.coord_ind(ic2,adj_ind(ineigh,icenter))) THEN
          nshared = nshared + 1
        END IF
      END DO
    END IF
  END DO

END SUBROUTINE COUNT_SHARED_VERTICES

SUBROUTINE PRINT_RESULTS
  USE histogram, ONLY : hist_poly
  IMPLICIT NONE
  INTEGER                    :: ic
  
  OPEN(unit = 2,file = "hist_polyhedral.dat")

  write(2,*) '# ', 'Polyhedral connection', 'Probability'

  DO ic = 1, 3
    WRITE(2,*) ic, float(hist_poly(ic))/float(sum(hist_poly))
  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
