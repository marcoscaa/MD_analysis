MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, nframes, nequil, stride
  INTEGER                            :: nhist, nk, ntype, nq
  INTEGER, ALLOCATABLE               :: atype(:), counter(:)
  DOUBLE PRECISION                             :: box(3), qmax
  DOUBLE PRECISION, ALLOCATABLE                :: pos(:,:), S(:)
  DOUBLE PRECISION, ALLOCATABLE                :: q_norm(:)
  DOUBLE PRECISION,PARAMETER                   :: Pi=3.1415
  DOUBLE PRECISION, ALLOCATABLE, DEVICE        :: pos_d(:,:), form_d(:,:), q_vector(:,:)
  DOUBLE PRECISION, ALLOCATABLE, DEVICE        :: Stmp(:), q_norm_d(:)
  INTEGER, ALLOCATABLE, DEVICE       :: atype_d(:)
END MODULE parameters

PROGRAM SofK
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  CALL Q_INDEX_TABLE
  
  DO frame = 1,nframes
    CALL READ_ATOM
    if (MOD(frame,stride)==0) THEN
      CALL COMPUTE_ATOMIC_FORM_FACTORS
      CALL CUDA_COMPUTE_SofK
    end if
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM SofK

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, pos, ntype, counter, stride,&
                         nhist, S, nframes, nequil, &
                         box, qmax, atype, pos_d, atype_d
  IMPLICIT NONE
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, ntype, nframes, nequil, stride, nhist, qmax
  READ(1,*) box !3 dimensions
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(S(nhist)); S=0
  ALLOCATE(counter(nhist)); counter=0

  !Device allocation
  allocate(pos_d(3,natoms))
  allocate(atype_d(natoms))

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
  INTEGER                    :: iat, ind
  DOUBLE PRECISION                     :: box_tmp(2)
 
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
    READ(1,*) ind, atype(ind), pos(:,ind)
    !pos(:,ind) = pos(:,ind) !* box !Non-reduced coordinates
  END DO

END SUBROUTINE READ_ATOM

SUBROUTINE Q_INDEX_TABLE
  USE parameters, ONLY : q_vector, q_norm, nq, Pi, qmax, box, nk, &
                         Stmp, q_norm_d, ntype, form_d
  IMPLICIT NONE
  INTEGER                    :: i,j,k,c
  DOUBLE PRECISION                     :: q, start, finish
  DOUBLE PRECISION, ALLOCATABLE        :: qtmp(:,:), formtmp(:,:)
  
  CALL cpu_time(start)

  nk=ceiling(qmax*maxval(box)/(2*Pi))
  nq=(2*nk+1)*(2*nk+1)*(2*nk+1)
  ALLOCATE(q_vector(3,nq)) !device
  ALLOCATE(qtmp(3,nq)) !host
  ALLOCATE(q_norm(nq)) !host
  ALLOCATE(q_norm_d(nq)) !device
  ALLOCATE(form_d(nq,ntype)) !device
 
  c=1
  DO i=-nk,nk
    DO j=-nk,nk
      DO k=-nk,nk
        q = 2.0*Pi*sqrt(float(i*i + j*j +k*k))
        qtmp(:,c) = 2.0*Pi*(/ float(i), float(j), float(k) /)
        q_norm(c) = q
        c = c + 1
      END DO
    END DO
  END DO

  if (c-1.ne.nq) then
    print *, "c ne q ", c, nq
    STOP
  end if

  q_vector = qtmp !Host to Device
  q_norm_d = q_norm ! H to D

  call cpu_time(finish)

  print *, 'Q_Table: ', finish-start

END SUBROUTINE Q_INDEX_TABLE

SUBROUTINE COMPUTE_ATOMIC_FORM_FACTORS
  USE parameters, ONLY : box, form_d, ntype, nq, q_norm_d
  USE cudafor
  IMPLICIT NONE
  DOUBLE PRECISION                     :: q, boxd 
  !DOUBLE PRECISION, DEVICE            :: Atomic_fq 
  INTEGER                              :: itype, iq, ig
  INTEGER                              :: get_natoms_type
  DOUBLE PRECISION, PARAMETER          :: fourpi=12.566370614359172
  DOUBLE PRECISION, DEVICE             :: a(4,2), b(4,2), c(2)
  LOGICAL, PARAMETER                   :: xray=.false.
  
  !Parameters of the atomic form factor
  !Parameters of the atomic form factor
  IF (xray) then
    !!Ti
    a(:,1) = (/ 9.7595, 7.3558, 1.6991, 1.9021 /)
    b(:,1) = (/ 7.8508, 0.5, 35.6338, 116.105 /)
    c(1)   = 1.2807
    !!O
    a(:,2) = (/ 3.0485, 2.2868, 1.5463, 0.867 /)
    b(:,2) = (/ 13.277, 5.7011, 0.3239, 32.9089 /)
    c(2)   = 0.2508
  ELSE
    !!Ti
    a(:,1) = (/ 3.5653, 2.8181, 1.8930, 0.4825 /)
    b(:,1) = (/ 81.9821, 19.0486, 3.5904, 0.3855 /)
    c(1)   = 0.0
    !!O
    a(:,2) = (/ 0.4548, 0.9173, 0.4719, 0.1384 /)
    b(:,2) = (/ 23.7803, 7.6220, 2.1440, 0.2959 /)
    c(2)   = 0.0
  ENDIF

  boxd = box(1)

  !$cuf kernel do(2)
  DO itype=1,ntype
    DO iq=1,nq
      q = q_norm_d(iq)/boxd !Only works for cubic cell
      !form_d(iq,itype) = Atomic_fq(q,itype) 
      form_d(iq,itype) = 0.0
      DO ig = 1,4
        form_d(iq,itype) = form_d(iq,itype) + &
                           a(ig,itype) * exp( -b(ig,itype) * (q/fourpi)**2 )
      END DO
      form_d(iq,itype) = form_d(iq,itype) + c(itype)
    END DO
  END DO

  !form_d(:,1)=22.0
  !form_d(:,2)=8.0

END SUBROUTINE COMPUTE_ATOMIC_FORM_FACTORS

SUBROUTINE CUDA_COMPUTE_SofK
  USE parameters, ONLY : natoms, pos, nhist, q_norm, qmax, &
                         pos_d, atype_d, box, S, atype, &
                         form_d, q_vector, nq, counter, Pi
  USE cudafor
  IMPLICIT NONE
  INTEGER                     :: iat, jat
  INTEGER                     :: qi, qind
  DOUBLE PRECISION                      :: f1,f2,normf
  DOUBLE PRECISION                      :: dp,start,finish 

  call cpu_time(start)

  pos_d=pos
  atype_d=atype

  !print *, 'Before starting'
  DO qi = 1,nq 
    qind=int(nhist*q_norm(qi)/qmax/box(1)) + 1
    if (qind.le.nhist) then
      f1=0.0;f2=0.0;normf=0.0
      !$cuf kernel do(1)<<<*,*>>>
      DO iat = 1, natoms
        dp = sum(q_vector(:,qi) * pos_d(:,iat))
        normf = normf + form_d(qi,atype_d(iat))**2
        !normf = 432*form_d(qi,1)**2 + 864*form_d(qi,2)**2
        f1 = f1 + cos(dp)*form_d(qi,atype_d(iat))
        f2 = f2 + sin(dp)*form_d(qi,atype_d(iat))
      END DO
      !normf = 432*form_d(qi,1)**2 + 864*form_d(qi,2)**2
      S(qind) = S(qind) + ( f1**2 + f2**2 ) / normf!**2
      counter(qind) = counter(qind) + 1
    end if
  END DO

  call cpu_time(finish)
  print *, 'S of K: ', finish-start
  !print *, 'Passed' 

END SUBROUTINE CUDA_COMPUTE_SofK

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, S, qmax, counter
  IMPLICIT NONE
  DOUBLE PRECISION                     :: bin!, aver, vol
  !DOUBLE PRECISION, PARAMETER          :: FourPioverThree = 4.18879
  INTEGER                    :: ihist
  
  OPEN(unit = 2,file = "sofk.dat")

  write(2,*) '#', 'K', 'SofK'

  S = S/real(counter)
  bin =  qmax/float(nhist)

  DO ihist = 1, nhist

    !vol = FourPioverThree * ( (ihist*bin)**3 - ((ihist-1)*bin)**3 )
    !aver = float(gofr(ihist)) / float(N*nframes) / vol
    WRITE(2,*) bin*(float(ihist)-0.5), S(ihist) 

  END DO

  CLOSE(2)

END SUBROUTINE

INTEGER FUNCTION get_natoms_type(atype,natoms,ind_rdf)
  !Return natoms with type ind_rdf
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: natoms, atype(natoms), ind_rdf
  INTEGER                    :: iat, N

  N = 0

  DO iat = 1, natoms

    IF ( ( atype(iat).eq.ind_rdf ) .or. &
       ( ind_rdf.eq.0 ) ) N = N + 1 

  END DO

  get_natoms_type = N

END FUNCTION get_natoms_type

DOUBLE PRECISION FUNCTION Dist(pos1,pos2,box)
  ! Distance between two points including pbc
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)         :: pos1(3),pos2(3),box(3)
  DOUBLE PRECISION                     :: xyz(3), SQRT
  INTEGER                    :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos1(i) - pos2(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = ( SUM(xyz*xyz) )**0.5

END FUNCTION Dist

DOUBLE PRECISION ATTRIBUTES(DEVICE) FUNCTION Atomic_Fq(atype,q)
  USE cudafor
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: atype
  DOUBLE PRECISION, INTENT(IN)         :: q
  DOUBLE PRECISION                     :: fourpi=12.566370614359172
  DOUBLE PRECISION                     :: fq, a(4,2), b(4,2), c(2)
  INTEGER                    :: ig

  !Parameters of the atomic form factor
  !Obtained from:
  !http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
  !!Ti
  a(:,1) = (/ 9.7595, 7.3558, 1.6991, 1.9021 /)
  b(:,1) = (/ 7.8508, 0.5, 35.6338, 116.105 /)
  c(1)   = 1.2807
  !!O
  a(:,2) = (/ 3.0485, 2.2868, 1.5463, 0.867 /)
  b(:,2) = (/ 13.277, 5.7011, 0.3239, 32.9089 /)
  c(2)   = 0.2508

  fq = 0.0
  DO ig = 1,4
    fq = fq + a(ig,atype) * exp( -b(ig,atype) * (q/fourpi)**2 )
  END DO

  fq = fq + c(atype)

  !Atomic_Fq = f

END FUNCTION Atomic_Fq
