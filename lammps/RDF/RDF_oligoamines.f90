!Computes the average coordination number of 
!atoms along the axis normal to the surface (z)

MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, stride
  INTEGER                            :: nframes, nequil, nhist
  INTEGER                            :: ind_rdf(2)
  INTEGER, ALLOCATABLE               :: NC_coordnumb(:), H_id(:)
  INTEGER, ALLOCATABLE               :: atype(:), gofr(:) 
  REAL*8                             :: box(3), maxr, number_density
  REAL*8, ALLOCATABLE                :: pos(:,:)
END MODULE parameters

PROGRAM RDF
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM
    if (MOD(frame,stride)==0) THEN 
      CALL CLASSIFY_N_AND_H_TYPES
      CALL COMPUTE_RDF_NH
    end if
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, number_density, &
                         nhist, gofr, nframes, nequil, &
                         box, ind_rdf, maxr, atype, stride, &
                         H_id, NC_coordnumb
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil,stride, nhist
  READ(1,*) box !3 dimensions
  READ(*,*) ind_rdf !Atomic species to be computed by RDF
                    !0 means all species

  maxr=minval(box)/2. !Max distance considered in g(r)
  number_density = 0.d0

  !Allocate module arrays
  ALLOCATE(pos(natoms,3))
  ALLOCATE(atype(natoms))
  ALLOCATE(gofr(nhist)); gofr=0
  ALLOCATE(H_id(natoms)); H_id=0
  ALLOCATE(NC_coordnumb(natoms)); NC_coordnumb=0

  !DO i=1,natoms
  !  READ(1,*) atype(i)
  !END DO

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

SUBROUTINE CLASSIFY_N_AND_H_TYPES
  USE parameters, only : H_id, NC_coordnumb, natoms, atype
  IMPLICIT NONE
  INTEGER :: iat, jat
  REAL*8 :: Dist, d
  REAL*8, parameter :: rcut_H = 1.5d0, rcut_CN=1.8d0
  
  H_id=0
  NC_coordnumb=0
  do iat=1,natoms
    IF (atype(iat)==2) then
      do jat=1,natoms
        if (atype(jat)==1) then 
          d=Dist(iat,jat)
          if (d<rcut_CN) NC_coordnumb(iat) = NC_coordnumb(iat) + 1
        end if
      end do
    ELSE IF (atype(iat)==3) then
      do jat=1,natoms
        if ((atype(jat)==2).or.(atype(jat)==4)) then 
          d=Dist(iat,jat)
          !if (d<rcut_H) H_id(iat) = atype(jat)
          if (d<rcut_H) H_id(iat) = jat
        end if
      end do
    END IF
  end do

END SUBROUTINE CLASSIFY_N_AND_H_TYPES

SUBROUTINE COMPUTE_RDF_NH
  USE parameters, ONLY : natoms, pos, ind_rdf, box, &
                         nhist, gofr, atype, maxr, number_density, &
                         H_id, NC_coordnumb
  IMPLICIT NONE
  INTEGER                     :: iat, iat2
  INTEGER                     :: n_coord 
  INTEGER                     :: ind, gofr_tmp(nhist), n_tmp
  INTEGER                     :: n
  REAL*8                      :: Dist, d

  n = 0 !total number of pairs
  DO iat = 1, natoms

    IF ( ( atype(iat) == 2 ) .and. ( (NC_coordnumb(iat) == ind_rdf(1)) .or. &
         ( ind_rdf(1) == -1 ) ) ) THEN

      gofr_tmp=0
! !$omp parallel do private(d,ind) reduction(+:gofr_tmp)
      DO iat2 = 1, natoms

      !IF ( ( atype(iat2) == 3 ) .and. ( H_id(iat2) .eq. ind_rdf(2) ) ) THEN
        IF ( ( ( atype(iat2) == 3 ) .and. ( H_id(iat2) .ne. 0 ) .and. & 
             ( ( NC_coordnumb(H_id(iat2)) .eq. ind_rdf(2) )  .or. &
               ( ind_rdf(2) == -1 ) ) ) ) THEN 

          d = Dist(iat,iat2)
          ind = int(nhist*d/maxr) + 1
          IF (ind<=nhist) gofr_tmp(ind) = gofr_tmp(ind) + 1
          n=n+1

        END IF

      END DO
! !$omp end parallel do
      gofr = gofr + gofr_tmp

    END IF

  END DO

  !n = 0 !total number of pairs
  !do iat=1,natoms
  !  do iat2=1,natoms
  !    if ((atype(iat)==2).and.(atype(iat2)==3)) n=n+1
  !  end do
  !end do

  number_density = number_density + float(n)/product(box)

END SUBROUTINE COMPUTE_RDF_NH

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, gofr, maxr, number_density, &
                         atype, natoms, ind_rdf, nframes
  IMPLICIT NONE
  REAL*8                     :: bin, aver, vol
  REAL*8, PARAMETER          :: FourPioverThree = 4.18879
  INTEGER                    :: ihist, N 
  INTEGER                    :: get_natoms_type
  
  OPEN(unit = 2,file = "rdf.dat")

  write(2,*) '#', 'Z', 'RDF'

  bin =  maxr/float(nhist)
  number_density = number_density !/ float(nframes)

  N = get_natoms_type(atype,natoms,ind_rdf(1)) 

  DO ihist = 1, nhist

    vol = FourPioverThree * ( (float(ihist)*bin)**3 - ((float(ihist)-1)*bin)**3 )
    aver = float(gofr(ihist)) / number_density / vol
    !aver = float(gofr(ihist)) /  vol
    WRITE(2,*) bin*(float(ihist)-0.5), aver 

  END DO

  CLOSE(2)

END SUBROUTINE

INTEGER FUNCTION get_natoms_type(atype,natoms,ind_rdf)
  !Return natoms with type ind_rdf
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: natoms, atype(natoms), ind_rdf
  INTEGER                    :: iat, irdf, N

  N = 0

  DO iat = 1, natoms

    IF ( ( atype(iat).eq.ind_rdf ) .or. &
       ( ind_rdf.eq.0 ) ) N = N + 1 

  END DO

  get_natoms_type = N

END FUNCTION get_natoms_type

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
