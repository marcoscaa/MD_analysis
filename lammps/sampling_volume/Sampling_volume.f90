!Computes the average coordination number of 
!atoms along the axis normal to the surface (z)

MODULE local
  IMPLICIT NONE 
  REAL*16, ALLOCATABLE               :: NofR(:), NofR2(:) 
  REAL*8                             :: maxr
  INTEGER                            :: nconf
  INTEGER, PARAMETER                 :: ncenter=2000
END MODULE local

PROGRAM RDF
  USE parameters, ONLY : nframes, stride
  USE local, only : nconf
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    if (MOD(frame,stride)==0) then 
      CALL COMPUTE_NofR
      nconf = nconf + 1
    end if
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, &
                         nhist, nframes, nequil, &
                         box, atype, stride
  USE local
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

  maxr=minval(box)/2. !Max distance considered 
  nconf = 0

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(NofR(nhist)); NofR=0
  ALLOCATE(NofR2(nhist)); NofR2=0

  CLOSE(1)
  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_NofR
  USE parameters, ONLY : natoms, pos, box, nhist, atype
  USE local
  IMPLICIT NONE
  INTEGER                     :: iat, ind, ic
  INTEGER*8                   :: N(nhist), Ntmp(nhist), N2(nhist)
  REAL*8                      :: Distance, d, origin(3)

  N = 0;N2=0
  DO ic = 1, ncenter
    CALL SELECT_RANDOM_ORIGIN(origin,box)
    Ntmp=0
    !$omp parallel do private(d,ind) reduction(+:Ntmp)
    DO iat = 1, natoms
      d = Distance(pos(:,iat),origin) 
      ind = int(float(nhist)*d/maxr) + 1
      IF (ind<=nhist) CALL AddONE(Ntmp, ind, nhist)
    END DO
    !$end omp parallel do
    N = N + Ntmp
    N2 = N2 + Ntmp*Ntmp
  END DO

  NofR  = ( NofR*float(nconf) + float(N)/float(ncenter) ) / float(nconf+1)
  NofR2  = ( NofR2*float(nconf) + float(N2)/float(ncenter) ) / float(nconf+1)

  print *, NofR(nhist), NofR2(nhist)

END SUBROUTINE COMPUTE_NofR

SUBROUTINE PRINT_RESULTS
  USE local, ONLY :  NofR, NofR2, nconf, maxr, ncenter
  USE parameters, ONLY : nhist
  IMPLICIT NONE
  REAL*8                     :: bin, SofZero, vol
  INTEGER                    :: ihist 
  
  OPEN(unit = 2,file = "NofR.dat")

  write(2,*) '#', 'R', 'NofR'

  bin =  maxr/float(nhist)

  DO ihist = 1, nhist

    SofZero = ( NofR2(ihist) - NofR(ihist)**2 ) / NofR(ihist)  
    !aver = float(gofr(ihist)) /  vol
    WRITE(2,*) bin*(float(ihist)-0.5), SofZero 

  END DO

  CLOSE(2)

END SUBROUTINE

SUBROUTINE SELECT_RANDOM_ORIGIN(origin,box)
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: box(3)
  REAL*8                     :: origin(3), r
  INTEGER                    :: i

  DO i=1,3
    CALL RANDOM_NUMBER(r)
    origin(i) = r*box(i)
  END DO

END SUBROUTINE SELECT_RANDOM_ORIGIN

SUBROUTINE AddONE(hist,ind,maxind)
  IMPLICIT NONE
  INTEGER, INTENT(in)        :: ind, maxind
  INTEGER*8                  :: hist(maxind)
  INTEGER                    :: i

  DO i=ind,maxind
    hist(i) = hist(i) + 1
  END DO

END SUBROUTINE AddONE

REAL*8 FUNCTION Distance(xyz1,xyz2)
  ! Distance between two points including pbc
  USE parameters, ONLY : box
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: xyz1(3),xyz2(3)
  REAL*8                     :: xyz(3)
  INTEGER                    :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = xyz1(i) - xyz2(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Distance = SQRT( SUM(xyz*xyz) )

END FUNCTION Distance
