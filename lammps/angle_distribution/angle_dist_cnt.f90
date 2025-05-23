!(atype2-atype1-atype2) Angle distribution

MODULE histogram
  IMPLICIT NONE 
  INTEGER                            :: atype1, atype2
  INTEGER,PARAMETER                  :: maxcoord=16 !buffer
  INTEGER, ALLOCATABLE               :: hist_angle(:,:) 
  REAL*8                             :: d_cut 
END MODULE histogram

PROGRAM angle_dist 
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM angle_dist

SUBROUTINE INITIALIZE
  USE parameters, only : natoms,nframes,nequil,nhist, pos, atype 
  USE histogram
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nhist
  READ(1,*) atype1, atype2,  d_cut !User defined distance cutoff for bonds
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(hist_angle(nhist,2));hist_angle=0

  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

INTEGER FUNCTION delta(x,x0)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x,x0
  INTEGER :: dout
 
  IF (x<=x0) THEN
    dout=1
  ELSE
    dout=2
  END IF

  delta = dout
END FUNCTION delta

SUBROUTINE MAKE_HISTOGRAM
  USE parameters, ONLY : natoms, atype, nhist, pos
  USE histogram, ONLY : atype1, hist_angle, maxcoord
  IMPLICIT NONE
  INTEGER                     :: iat, n1, n2
  INTEGER                     :: nbonds,ind, ind2, delta
  REAL*8                      :: bond_dir(3,maxcoord), angle
  REAL*8                      :: center(3), Dist_Center, d

  center=(/0.0,0.0,0.0/)

  DO iat = 1, natoms
    IF (atype(iat)==atype1) THEN
      d=Dist_Center(center,pos(:,iat))
      ind2=delta(d,10.d0)
      CALL CONNECT_NEAREST_NEIGHBOHRS(iat,bond_dir,nbonds)
      DO n1 = 1,nbonds
        DO n2 = n1+1,nbonds
          angle = sum(bond_dir(:,n1)*bond_dir(:,n2))
          ind = int(nhist*(angle+1)/2.) + 1
          hist_angle(ind,ind2) = hist_angle(ind,ind2) + 1
        END DO
      END DO
    END IF
  END DO

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist 
  USE histogram, ONLY : hist_angle
  IMPLICIT NONE
  INTEGER                    :: ic
  REAL*8                     :: bin
  
  OPEN(unit = 2,file = "hist_angle_cnt.dat")
  write(2,*) '# ', 'cos(Angle)', 'Count(inner)', 'Count(outer_)'

  bin = 2./real(nhist)

  DO ic = 1, nhist
    WRITE(2,*) bin*(real(ic)-0.5) - 1.0, hist_angle(ic,1), hist_angle(ic,2) 
  END DO

  CLOSE(2)

END SUBROUTINE

SUBROUTINE CONNECT_NEAREST_NEIGHBOHRS( iat,bonds,nbonds )
  !Compute the coordination number of atom with index iat
  USE parameters, ONLY : natoms, atype
  USE histogram, ONLY : atype2, d_cut, maxcoord
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat 
  INTEGER                    :: iat2, nbonds
  REAL*8                     :: d(3),normd,bonds(3,maxcoord)

  nbonds = 0
  bonds=0

  DO iat2 = 1, natoms
    IF ( (iat2.ne.iat).and.(atype(iat2).eq.atype2) ) THEN
      CALL DIST_UNIT_VECTOR(iat,iat2,d,normd) 
      IF ( normd < d_cut ) THEN 
        nbonds = nbonds + 1
        bonds(:,nbonds) = d  
      END IF
    END IF 
  END DO
 
END SUBROUTINE CONNECT_NEAREST_NEIGHBOHRS

