MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, nframes, nequil
  INTEGER                            :: nhist, nk, ntype, stride
  INTEGER, ALLOCATABLE               :: atype(:), counter(:)
  REAL*8                             :: box(3), qmax
  REAL*8, ALLOCATABLE                :: pos(:,:), form_fact(:), S(:)
  REAL*8,PARAMETER                   :: Pi=3.1415
END MODULE parameters

PROGRAM SofK
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM
    IF (MOD(frame,stride)==0) CALL COMPUTE_SofK_Faber_Ziman
  END DO

  !CALL AVERAGE_SofK
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM SofK

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, pos, ntype, counter, &
                         nhist, S, nframes, nequil, &
                         box, Pi, qmax, atype, nk, form_fact, &
                         stride
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

  nk=ceiling(qmax*maxval(box)/(2*Pi))

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(S(nhist)); S=0
  ALLOCATE(counter(nhist)); counter=0

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
    pos(:,ind) = pos(:,ind) * box !Non-reduced coordinates
  END DO

END SUBROUTINE READ_ATOM

SUBROUTINE COMPUTE_SofK
  USE parameters, ONLY : natoms, pos, ntype, Pi, qmax, &
                         nhist, S, atype, nk, form_fact, &
                         box, counter
  IMPLICIT NONE
  INTEGER                     :: iat
  INTEGER                     :: i,j,k
  INTEGER                     :: qind 
  REAL*8                      :: f(ntype,2), q(3), dp 

  DO i = 0,nk !-nk,nk
    DO j = 0,nk !-nk,nk
      DO k = 0,nk !-nk,nk
        q = (/ i,j,k /) * 2*Pi/box
        if ( norm2(q) < qmax ) THEN
          f=0.0
          DO iat = 1, natoms
            dp = dot_product(q,pos(:,iat))
            f(atype(iat),1) = f(atype(iat),1) + cos(dp)
            f(atype(iat),2) = f(atype(iat),2) + sin(dp)
          END DO
          qind=int(nhist*norm2(q)/qmax) + 1
          S(qind) = S(qind) + sum( form_fact*form_fact *( f(:,1)**2 + f(:,2)**2 ) )
          counter(qind) = counter(qind) + 1
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_SofK

SUBROUTINE COMPUTE_SofK_Faber_Ziman
  USE parameters, ONLY : natoms, pos, ntype, Pi, qmax, &
                         nhist, S, atype, nk, form_fact, &
                         box, counter
  IMPLICIT NONE
  INTEGER                     :: iat
  INTEGER                     :: i,j,k
  INTEGER                     :: qind 
  REAL*8                      :: f(ntype,2), q(3), dp 
  REAL*8                      :: newS

  DO i = -nk,nk !-nk,nk
    DO j = -nk,nk !-nk,nk
      DO k = 0,nk !-nk,nk
        q = (/ float(i),float(j),float(k) /) * 2.0*Pi/box
        !CALL COMPUTE_ATOMIC_FORM_FACTORS(norm2(q))
        if ( norm2(q) < qmax ) THEN
          f=0.0
          DO iat = 1, natoms
            dp = dot_product(q,pos(:,iat))
            f(atype(iat),1) = f(atype(iat),1) + cos(dp)
            f(atype(iat),2) = f(atype(iat),2) + sin(dp)
          END DO
          qind=int(nhist*norm2(q)/qmax) + 1
          !newS = sum( form_fact*form_fact *( f(:,1)**2 + f(:,2)**2 ) ) &
          newS = form_fact(1)*form_fact(1) *( f(1,1)**2 + f(1,2)**2 ) &
               + form_fact(2)*form_fact(2) *( f(2,1)**2 + f(2,2)**2 ) &
               + 2.d0 * form_fact(1) * form_fact(2) * ( f(1,1) * f(2,1) + f(1,2) * f(2,2) )
          newS = newS / dble(natoms) - (1./3.)*form_fact(1)**2 - (2.0/3.0)*form_fact(2)**2
          newS = newS / ((1./3.)*form_fact(1) + (2.0/3.0)*form_fact(2))**2 !/ dble(natoms)
          !newS = (f(1,1)+f(2,1))**2 + (f(1,2)+f(2,2))**2
          !newS = ( f(1,1)**2 + f(1,2)**2 ) + ( f(2,1)**2 + f(2,2)**2 ) &
          !     + 2.d0 * ( f(1,1) * f(2,1) + f(1,2) * f(2,2) )
          S(qind) = S(qind) + newS 
          counter(qind) = counter(qind) + 1
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_SofK_Faber_Ziman

SUBROUTINE COMPUTE_ATOMIC_FORM_FACTORS (q)
  USE parameters, ONLY : box, form_fact, ntype
  IMPLICIT NONE
  REAL*8, INTENT(IN)                   :: q 
  !DOUBLE PRECISION, DEVICE            :: Atomic_fq 
  INTEGER                              :: itype,ig
  REAL*8, PARAMETER                    :: fourpi=12.566370614359172
  REAL*8                               :: a(4,2), b(4,2), c(2)
  LOGICAL, PARAMETER                   :: xray=.false.
  
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

  DO itype=1,ntype
    form_fact(itype) = 0.0
    DO ig = 1,4
      form_fact(itype) = form_fact(itype) + &
                         a(ig,itype) * exp( -b(ig,itype) * (q/fourpi)**2 )
    END DO
    form_fact(itype) = form_fact(itype) + c(itype)
  END DO

END SUBROUTINE COMPUTE_ATOMIC_FORM_FACTORS

SUBROUTINE AVERAGE_SofK
  USE parameters, ONLY : S, form_fact, atype, nframes, &
                         form_fact, natoms, ntype, counter
  IMPLICIT NONE
  INTEGER                    :: nattype(ntype), n1, n2
  INTEGER                    :: get_natoms_type
  REAL*8                     :: wijsum

  DO n1 = 1,ntype
    nattype(n1) = get_natoms_type(atype,natoms,n1)
  END DO

  wijsum=0
  DO n1 = 1,ntype
    DO n2 = 1,ntype
      wijsum = wijsum + form_fact(n1)*form_fact(n2)*dble(nattype(n1)*nattype(n2))
    END DO
  END DO

  S = S/wijsum/real(counter)


END SUBROUTINE AVERAGE_SofK

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, S, qmax, counter
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
    WRITE(2,*) bin*(float(ihist)-0.5), S(ihist)/float(counter(ihist)) 

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

REAL*8 FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  REAL*8                     :: xyz(3)
  INTEGER                    :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos(i,ind1) - pos(i,ind2)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist
