PROGRAM RDF
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_EXTXYZ
    if (MOD(frame,stride)==0) THEN
      CALL FIND_CENTER_CNT
      CALL COMPUTE_RDF
    end if
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, number_density, &
                         nhist, gofr, nframes, nequil, &
                         box, ind_rdf, maxr, atype, stride, fixed_coord
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil,stride, nhist, maxr
  READ(*,*) ind_rdf !Atomic species to be computed by RDF
                    !0 means all species
  if(trim(ind_rdf(1))=="fixed") then
    READ(1,*) fixed_coord          
  end if

  number_density = 0.d0

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(gofr(nhist)); gofr=0

  CLOSE(1)
  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE FIND_CENTER_CNT
  USE parameters, only : pos, box, natoms, fixed_coord, atype
  IMPLICIT NONE
  REAL*8 :: sum_coord(3),d(4)
  INTEGER :: i, nc

  nc=0
  sum_coord=0.d0
  DO i = 1,natoms
    IF (trim(atype(i))=="C") THEN
      CALL DISTANCE_VECTOR(pos(:,i),fixed_coord,d)
      sum_coord = sum_coord + d(1:3)
      nc = nc+1
    END IF
  END DO

  fixed_coord = fixed_coord + sum_coord/float(nc)

END SUBROUTINE FIND_CENTER_CNT

SUBROUTINE COMPUTE_RDF
  USE parameters, ONLY : natoms, pos, ind_rdf, box, fixed_coord, &
                         nhist, gofr, atype, maxr, number_density
  IMPLICIT NONE
  INTEGER                     :: iat, iat2
  INTEGER                     :: ind, gofr_tmp(nhist)
  REAL*8                      :: DistXY, d, central_atom(3)

    gofr_tmp=0
    central_atom=fixed_coord !pos(:,iat)

!$omp parallel do private(d,ind) reduction(+:gofr_tmp)
    DO iat2 = 1, natoms

      IF ( ( ( trim(atype(iat2)) == trim(ind_rdf(2)) ) .or. &
          ( trim(ind_rdf(2)) == "All" ) ) .and. ( iat.ne.iat2 ) ) THEN

        d = DistXY(central_atom,pos(:,iat2))
        ind = int(nhist*d/maxr) + 1
        IF (ind<=nhist) gofr_tmp(ind) = gofr_tmp(ind) + 1

      END IF

    END DO
!$omp end parallel do
    gofr = gofr + gofr_tmp
    number_density = number_density + 1.d0

END SUBROUTINE COMPUTE_RDF

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, gofr, maxr, number_density, &
                         atype, natoms, ind_rdf, nframes, box
  IMPLICIT NONE
  REAL*8                     :: bin, aver, vol
  REAL*8, PARAMETER          :: Pi = 3.1415
  INTEGER                    :: ihist 
  INTEGER                    :: get_natoms_type
  
  OPEN(unit = 2,file = "rdf.dat")

  write(2,*) '#', 'Z', 'RDF'

  bin =  maxr/float(nhist)

  DO ihist = 1, nhist

    vol = Pi * ( (float(ihist)*bin)**2 - ((float(ihist)-1)*bin)**2 ) * box(3,3)
    aver = float(gofr(ihist)) / vol / number_density 
    !aver = float(gofr(ihist)) /  vol
    WRITE(2,*) bin*(float(ihist)-0.5), aver 

  END DO

  CLOSE(2)

END SUBROUTINE

