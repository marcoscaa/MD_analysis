MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, nframes, nequil, stride
  INTEGER                            :: nhist, nk, ntype, nq
  INTEGER, ALLOCATABLE               :: atype(:), counter(:)
  REAL*8                             :: box(3), qmax
  REAL*8, ALLOCATABLE                :: pos(:,:), form_fact(:), S(:)
  REAL*8, ALLOCATABLE                :: q_norm(:)
  REAL*8,PARAMETER                   :: Pi=3.1415
  REAL*8, ALLOCATABLE, DEVICE        :: pos_d(:,:), form_d(:), q_vector(:,:)
  REAL*8, ALLOCATABLE, DEVICE        :: Stmp(:), q_norm_d(:)
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
    if (MOD(frame,stride)==0) CALL CUDA_COMPUTE_SofK
    !CALL MAKE_HISTOGRAM
    !if (MOD(frame,stride)==0) CALL CUDA_COMPUTE_SPHERICAL_SofK
  END DO

  CALL AVERAGE_SofK
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM SofK

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, pos, ntype, counter, stride,&
                         nhist, S, nframes, nequil, &
                         box, qmax, atype, form_fact, &
                         pos_d, atype_d, form_d
  IMPLICIT NONE
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, ntype, nframes, nequil, stride, nhist, qmax
  READ(1,*) box !3 dimensions
  ALLOCATE(form_fact(ntype))
  READ(1,*) form_fact
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(S(nhist)); S=0
  ALLOCATE(counter(nhist)); counter=0

  !Device allocation
  allocate(pos_d(3,natoms))
  allocate(form_d(ntype))
  allocate(atype_d(natoms))
  form_d=form_fact

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
    READ(1,*) ind, atype(ind), pos(:,ind)
    pos(:,ind) = pos(:,ind) !* box !Non-reduced coordinates
  END DO

END SUBROUTINE READ_ATOM

SUBROUTINE Q_INDEX_TABLE
  USE parameters, ONLY : q_vector, q_norm, nq, Pi, qmax, box, nk, &
                         Stmp, q_norm_d
  IMPLICIT NONE
  INTEGER                    :: i,j,k,c
  REAL*8                     :: q, start, finish
  REAL*8, ALLOCATABLE        :: qtmp(:,:)
  
  CALL cpu_time(start)

  nk=ceiling(qmax*maxval(box)/(2*Pi))
  nq=(2*nk+1)*(2*nk+1)*(nk+1)
  ALLOCATE(q_vector(3,nq)) !device
  ALLOCATE(qtmp(3,nq)) !host
  ALLOCATE(q_norm(nq)) !host
  ALLOCATE(q_norm_d(nq)) !device
 
  c=1
  DO i=-nk,nk
    DO j=-nk,nk
      DO k=0,nk
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

END SUBROUTINE

SUBROUTINE CUDA_COMPUTE_SofK
  USE parameters, ONLY : natoms, pos, nhist, q_norm, qmax, &
                         pos_d, atype_d, box, S, atype, &
                         form_d, q_vector, nq, counter, Pi
  USE cudafor
  IMPLICIT NONE
  INTEGER                     :: iat, jat
  INTEGER                     :: qi, qind
  REAL*8                      :: f1,f2
  REAL*8                      :: dp,start,finish 

  call cpu_time(start)

  pos_d=pos
  atype_d=atype

  !print *, 'Before starting'
  DO qi = 1,nq 
    qind=int(nhist*q_norm(qi)/qmax/box(1)) + 1
    if (qind.le.nhist) then
      f1=0.0;f2=0.0
      !$cuf kernel do(1)<<<*,*>>>
      DO iat = 1, natoms
        dp = sum(q_vector(:,qi) * pos_d(:,iat))
        f1 = f1 + cos(dp)*form_d(atype_d(iat))
        f2 = f2 + sin(dp)*form_d(atype_d(iat))
      END DO
      S(qind) = S(qind) + f1**2 + f2**2 
      counter(qind) = counter(qind) + 1
    end if
  END DO

  call cpu_time(finish)
  print *, 'S of K: ', finish-start
  !print *, 'Passed' 

END SUBROUTINE CUDA_COMPUTE_SofK

!SUBROUTINE MAKE_HISTOGRAM
!  !we are assuming the cell to be cubic
!  USE parameters, ONLY: S, Stmp, nq, nhist, qmax, box, counter, &
!                        q_norm
!  IMPLICIT NONE
!  INTEGER                    :: ihist, iq
!  REAL*8                     :: Stmp_h(nq), start, finish
!
!  call cpu_time(start)
!  print *, 'Make histogram, before memcopy'
!
!  Stmp_h = Stmp !device to host
!  print *, Stmp_h(nq), q_norm(nq)
!
!  print *, 'Make histogram, after memcopy'
!
!  DO iq = 1,nq
!    ihist = int(float(nhist)*q_norm(iq)/qmax/box(1)) + 1
!    if (ihist.le.nhist) then
!      S(ihist) = S(ihist) + Stmp_h(iq)
!      counter(ihist) = counter(ihist) + 1
!    end if
!  END DO
!
!  print *, Stmp_h(nq), counter(nhist)
!
!  call cpu_time(finish)
!  print *, 'Make Histogram: ', finish-start
!
!END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE CUDA_COMPUTE_SPHERICAL_SofK
  USE parameters, ONLY : natoms, pos_d, Stmp, atype_d, nq, form_d, &
                         q_norm_d, pos, atype
  USE cudafor
  IMPLICIT NONE
  INTEGER                     :: iat, jat,i
  INTEGER                     :: iq, c
  REAL*8                      :: q, r, S_d
  REAL*8, device              :: xyz(3)

  pos_d=pos
  atype_d=atype

  !print *, 'Before starting'

  !$cuf kernel do(3)<<<*,*>>>
  DO iq = 1,nq 
    DO iat = 1, natoms
      DO jat = 1,natoms
        DO i = 1,3
          xyz(i) = pos_d(i,iat) - pos_d(i,jat)
          xyz(i) = xyz(i) - nint( xyz(i) ) 
        END DO
        r = ( SUM(xyz*xyz) )**(0.5)
        q = q_norm_d(iq)
        if (iat.ne.jat) then
          Stmp(iq) = Stmp(iq) + &
              form_d(atype_d(iat))*form_d(atype_d(jat)) * sin(q*r)/(q*r)
        end if
      END DO
    END DO
  END DO

END SUBROUTINE CUDA_COMPUTE_SPHERICAL_SofK

SUBROUTINE AVERAGE_SofK
  USE parameters, ONLY : S, form_fact, atype, nframes, &
                         form_fact, natoms, ntype, counter
  IMPLICIT NONE
  INTEGER                    :: nattype(ntype), n1, n2
  INTEGER                    :: get_natoms_type
  REAL*8                     :: wijsum

  !DO n1 = 1,ntype
  !  nattype(n1) = get_natoms_type(atype,natoms,n1)
  !END DO

  wijsum=0
  DO n1 = 1,natoms
    DO n2 = 1,natoms
      wijsum = wijsum + form_fact(atype(n1))*form_fact(atype(n2))
      !wijsum = wijsum + form_fact(n1)*form_fact(n1)*dble(nattype(n1)*nattype(n1))
    END DO
  END DO

  S = S/wijsum/real(counter)


END SUBROUTINE AVERAGE_SofK

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, S, qmax
  IMPLICIT NONE
  REAL*8                     :: bin!, aver, vol
  !REAL*8, PARAMETER          :: FourPioverThree = 4.18879
  INTEGER                    :: ihist
  
  OPEN(unit = 2,file = "sofk.dat")

  write(2,*) '#', 'K', 'SofK'

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

REAL*8 FUNCTION Dist(pos1,pos2,box)
  ! Distance between two points including pbc
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos1(3),pos2(3),box(3)
  REAL*8                     :: xyz(3), SQRT
  INTEGER                    :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos1(i) - pos2(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = ( SUM(xyz*xyz) )**0.5

END FUNCTION Dist
