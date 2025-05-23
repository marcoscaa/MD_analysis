!Computes the average coordination number of 
!atoms along the axis normal to the surface (z)

MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, stride
  INTEGER                            :: nframes, nequil, nhist
  INTEGER                            :: ind_rdf(2)
  INTEGER, ALLOCATABLE               :: atype(:), lsi_histogram(:) 
  REAL*8                             :: box(3)
  REAL*8, ALLOCATABLE                :: pos(:,:)
  REAL*8, PARAMETER                  :: rcut=3.7d0,max_lsi=4.d0
END MODULE parameters

PROGRAM RDF
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM
    if (MOD(frame,stride)==0) CALL COMPUTE_LSI
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, lsi_histogram, &
                         nhist, nframes, nequil, &
                         ind_rdf, atype, stride
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil,stride, nhist
  READ(*,*) ind_rdf !Atomic species to be computed by RDF
                    !0 means all species

  !Allocate module arrays
  ALLOCATE(pos(natoms,3))
  ALLOCATE(atype(natoms))
  ALLOCATE(lsi_histogram(nhist)); lsi_histogram=0

  CLOSE(1)
  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE REMOVE_EQUIL
  USE parameters, ONLY : natoms, nframes, nequil
  IMPLICIT NONE
  INTEGER                    :: iframe

  DO iframe = 1,nequil*(natoms+9)
    READ(1,*)
  END DO

  nframes = nframes - nequil

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE READ_ATOM
  !Read LAMMPS atom file
  USE parameters, ONLY : pos, natoms, atype, box
  IMPLICIT NONE
  INTEGER                    :: iat, ind, junk
  REAL*8                     :: box_tmp(2)
 
  DO iat=1,5
    READ(1,*)
  END DO

  !Read Box
  DO iat=1,3
    READ(1,*) box_tmp(1), box_tmp(2)
    box(iat) = box_tmp(2)-box_tmp(1)
  END DO

  READ(1,*)

  DO iat = 1, natoms
    READ(1,*) ind, atype(ind), pos(ind,1), pos(ind,2), pos(ind,3)
    pos(ind,:) = pos(ind,:) * box !Non-reduced coordinates
  END DO

END SUBROUTINE READ_ATOM

SUBROUTINE COMPUTE_LSI
  USE parameters, ONLY : natoms, pos, ind_rdf, box, &
                         nhist, lsi_histogram, atype, &
                         rcut
  IMPLICIT NONE
  INTEGER                     :: iat
  INTEGER                     :: n_neighbors
  REAL*8                      :: lsi_op, distances(natoms)

  lsi_op=0.0
  DO iat = 1, natoms

    IF (atype(iat)==ind_rdf(1)) THEN

      CALL DISTANCES_WITHIN_RCUT(iat,ind_rdf(2),distances,rcut,n_neighbors) 
      if (n_neighbors.ne.0) then
        CALL QUICKSORT(distances,1,n_neighbors+1)
        CALL LSI(distances,n_neighbors+1,natoms,lsi_op)
        CALL MAKE_HISTOGRAM(lsi_histogram,lsi_op)
      end if

    END IF

  END DO

END SUBROUTINE COMPUTE_LSI

SUBROUTINE DISTANCES_WITHIN_RCUT(iat,atype_neigh,distances,rcut,n)
  USE parameters, only : pos,atype,natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iat,atype_neigh
  REAL*8, INTENT(IN) :: rcut
  INTEGER :: n, jat, c
  REAL*8 :: distances(natoms), Dist, d
  REAL*8 :: buffer

  buffer=0.5d0
  distances=1000.
101  n=0;c=0
  DO jat=1,natoms

    if ((atype(jat)==atype_neigh).and.(jat.ne.iat)) then
      d=Dist(iat,jat)
      IF (d<rcut+buffer) THEN
        c=c+1
        distances(c) = Dist(iat,jat)
        IF (d<rcut) n=n+1
      END IF
    end if

  END DO

  IF (n==c) THEN
    buffer=buffer+0.5
    if (buffer>10.0) then
       PRINT *, 'Too large buffer'
       STOP
    end if
    GO TO 101
  END IF

  !IF ((n==0).or.(c==n)) THEN
  IF (c==n) THEN
    PRINT *, 'Error in DISTANCES_WITHIN_RCUT'
    PRINT *, n, c
    PRINT *, 'STOP!'
    STOP
  END IF

END SUBROUTINE DISTANCES_WITHIN_RCUT

SUBROUTINE LSI(distances,n,ntot,lsi_op)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, ntot
  REAL*8, INTENT(IN) :: distances(ntot)
  REAL*8 :: lsi_op, delta(n-1), delta_mean
  INTEGER :: i,j

  delta_mean=0.d0
  DO i=1,n-1
    delta(i) = distances(i+1)-distances(i)
    delta_mean=delta_mean+delta(i)
  END DO

  delta_mean = delta_mean / float(n-1)

  DO i=1,n-1
    lsi_op = lsi_op + (delta(i) - delta_mean)**2
  END DO
  
  lsi_op = lsi_op / float(n-1)

END SUBROUTINE LSI

SUBROUTINE MAKE_HISTOGRAM(lsi_histogram,lsi_op)
  USE parameters, only : nhist, max_lsi
  IMPLICIT NONE
  INTEGER :: lsi_histogram(nhist), ind
  REAL*8,INTENT(IN) :: lsi_op
    
  ind = int(float(nhist)*lsi_op/max_lsi) + 1
  lsi_histogram(ind) = lsi_histogram(ind) + 1

END SUBROUTINE MAKE_HISTOGRAM

! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
RECURSIVE SUBROUTINE QUICKSORT(a, first, last)
  implicit none
  real*8  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
END SUBROUTINE QUICKSORT

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, lsi_histogram, max_lsi
  IMPLICIT NONE
  INTEGER                    :: ihist
  
  OPEN(unit = 2,file = "lsi_distribution.dat")

  write(2,*) '#', 'LSI', 'Count'

  DO ihist = 1, nhist

    WRITE(2,fmt='(F12.8,3X,I10)') max_lsi*float(2*ihist-1)/float(nhist)/2.d0, lsi_histogram(ihist) 

  END DO

  CLOSE(2)

END SUBROUTINE

REAL*8 FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  REAL*8                     :: xyz(3)
  INTEGER                    :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos(ind1,i) - pos(ind2,i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

INTEGER FUNCTION closest_atom(icenter,atyp)
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: icenter, atyp
  INTEGER :: iat, iat_min
  REAL*8 :: d, min_d, Dist

  min_d = 100.
  iat_min = 0

  DO iat=1,natoms

    IF (atype(iat)==atyp) THEN
      d = Dist(iat,icenter)
      IF (d < min_d) THEN
        min_d=d
        iat_min=iat
      END IF
    END IF

  END DO

  IF (iat_min==0) THEN
    PRINT *, "Could not find closest atom to atom", icenter
    PRINT *, "STOP!!!!"
    STOP
  END IF

  closest_atom = iat_min

END FUNCTION closest_atom
