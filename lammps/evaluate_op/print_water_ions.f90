!Computes the PT CV. This is not a general code!!!
!Only works if atype is ordered with O1 ... ON H1 ... H2N

MODULE OP
  REAL*8 :: lambda
  REAL*8, ALLOCATABLE :: charge(:)
  INTEGER,ALLOCATABLE :: min_ind(:)
END MODULE OP

PROGRAM Coord_number
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL COMPUTE_OP(frame)
    CALL PRINT_WATER_IONS 
  END DO

  CLOSE(1); CLOSE(2); CLOSE(3)
  CLOSE(4); CLOSE(5); CLOSE(6)

END PROGRAM Coord_number

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nwater, pos, atype, &
                         nframes, nequil
  USE OP, ONLY : lambda, min_ind, charge
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file, out_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)
  CALL getarg(3, out_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nwater, nframes, nequil, lambda
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(min_ind(2*nwater))
  ALLOCATE(charge(nwater))

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = out_file)
  OPEN(unit = 3,file = "pos-wi_O.dat")
  OPEN(unit = 4,file = "pos-wi_H.dat")
  OPEN(unit = 5,file = "pos-wi_O_env.dat")
  OPEN(unit = 6,file = "pos-wi_H_env.dat")

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_OP(step)
  !Compute center of mass along the z direction
  !Here we use the particle with index 1 as a reference
  USE parameters, ONLY : pos, nwater, natoms
  USE OP, ONLY: lambda, charge
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: step
  REAL*8                     :: compute_charge
  REAL*8                     :: dist_mat(natoms-nwater,nwater), sd
  REAL*8                     :: Dist, Total_charge
  REAL*8, PARAMETER          :: alpha=0.0001 !Smoothness of total charge
  INTEGER                    :: iw,jw

  CALL COMPUTE_DIST_MATRIX(dist_mat)

  CALL compute_charge_non_smooth(dist_mat,charge)

  sd = 0.d0

  DO iw = 1,nwater
!$omp parallel do reduction(+:sd)
    DO jw = iw+1,nwater
      sd = sd - charge(iw)*charge(jw)*Dist(iw,jw)
    END DO
!$end omp parallel do
  END DO

  Total_charge = sum(abs(charge))

  if (Total_charge>0) then
    WRITE(2,*) step, sd, Total_charge
  end if

END SUBROUTINE COMPUTE_OP

subroutine compute_charge_non_smooth(dist_mat, charge)
  use parameters, only: nwater, natoms
  use OP, only: min_ind
  implicit none
  integer                    :: iw, ih
  real*8                     :: charge(nwater)
  real*8                     :: mindist, dist_mat(natoms-nwater,nwater)

  min_ind=0
  charge = 0.d0
  do ih = nwater+1,natoms
    mindist=100.d0
    do iw = 1,nwater
      if (dist_mat(ih-nwater,iw)<mindist) then
        mindist=dist_mat(ih-nwater,iw)
        min_ind(ih-nwater)=iw
      endif
    end do
    charge(min_ind(ih-nwater)) = charge(min_ind(ih-nwater)) + 1.d0
  end do

  charge = charge - 2.d0

end subroutine compute_charge_non_smooth

SUBROUTINE COMPUTE_DIST_MATRIX(mat)
  USE parameters, only : natoms, nwater
  IMPLICIT NONE
  REAL*8                     :: mat(natoms-nwater,nwater)
  REAL*8                     :: Dist
  INTEGER                    :: iw, ih

  DO ih=nwater+1,natoms
    DO iw=1,nwater
      mat(ih-nwater,iw) = Dist(ih,iw)
    END DO
  END DO

END SUBROUTINE COMPUTE_DIST_MATRIX

SUBROUTINE PRINT_WATER_IONS
  USE parameters, only: natoms, nwater, pos
  USE op, only: min_ind, charge 
  IMPLICIT NONE
  INTEGER :: ih, iw, ind_h3o, c, n
  INTEGER, ALLOCATABLE :: ind_ion_pair(:)
  REAL*8 :: d(3), doh(3)
  REAL*8, ALLOCATABLE :: ion_pair(:,:)

  IF (.not.allocated(ion_pair)) then
    allocate(ion_pair(3,natoms))
    allocate(ind_ion_pair(natoms))
  END IF

  DO iw=1,nwater
    if(charge(iw)==1) then
      c=0
      ind_h3o=iw
      ion_pair(:,1) = (/ 0.0, 0.0, 0.0 /)
      ind_ion_pair(1) = iw
      DO ih=nwater+1,natoms
        if (min_ind(ih-nwater)==iw) then
          CALL DISTANCE_VECTOR( pos(:,iw), pos(:,ih), d)
          c=c+1
          ion_pair(:,c+1) = d
          ind_ion_pair(c+1) = ih
        end if
      END DO
      if (c/=3) print *, 'ERROR in H3O', c
    end if
  END DO

  DO iw=1,nwater
    if (charge(iw)==-1) then
      c=0
      CALL DISTANCE_VECTOR( pos(:,ind_h3o), pos(:,iw), doh)
      ion_pair(:,5) = doh
      ind_ion_pair(5)=iw
      DO ih=nwater+1,natoms
        if (min_ind(ih-nwater)==iw) then
          CALL DISTANCE_VECTOR( pos(:,iw), pos(:,ih), d)
          ion_pair(:,6) = d+doh
          ind_ion_pair(6) = ih
          c=c+1
        end if
      END DO
      if (c/=1) print *, 'ERROR in OH'
    end if
  END DO

  if (any(charge==1)) then

    CALL GET_ATOMS_AROUND_ION_PAIR( ind_ion_pair, ion_pair, n) 

    CALL PROJECT_TO_LOCAL_FRAME( ion_pair , n) 

    CALL WRITE_TO_FILE(ion_pair, ind_ion_pair, n)

  end if

END SUBROUTINE PRINT_WATER_IONS

SUBROUTINE WRITE_TO_FILE (ion_pair, ind, n)
  USE parameters, only : nwater
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: n, ind(n)
  REAL*8, INTENT(IN) :: ion_pair(3,n)
  INTEGER :: iw

  write(3,fmt="(4(F14.8,3X))") ion_pair(:,1), ion_pair(1,5)
  write(3,fmt="(4(F14.8,3X))") ion_pair(:,5), ion_pair(1,5)

  write(4,fmt="(4(F14.8,3X))") ion_pair(:,2), ion_pair(1,5)
  write(4,fmt="(4(F14.8,3X))") ion_pair(:,3), ion_pair(1,5)
  write(4,fmt="(4(F14.8,3X))") ion_pair(:,4), ion_pair(1,5)
  write(4,fmt="(4(F14.8,3X))") ion_pair(:,6), ion_pair(1,5)

  do iw=7,n
    IF (ind(iw).le.nwater) then
      write(5,fmt="(4(F14.8,3X))") ion_pair(:,iw), ion_pair(1,5)
    ELSE
      write(6,fmt="(4(F14.8,3X))") ion_pair(:,iw), ion_pair(1,5)
    END IF
  end do

END SUBROUTINE WRITE_TO_FILE

SUBROUTINE PROJECT_TO_LOCAL_FRAME (ion_pair, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL*8 :: ion_pair(3,n)
  REAL*8 :: rotmat(3,3), x(3), y(3), z(3)
  INTEGER :: i,j

  !X axis along O-O
  x = ion_pair(:,5)/norm2(ion_pair(:,5)) 
  !Y is perpendicular to the O(H3O)-O(oh)-H(oh) plane
  CALL CROSSPROD(ion_pair(:,5),ion_pair(:,6)-ion_pair(:,5),y)
  y = y/norm2(y)
  !Z is perpendicular to X and Y
  CALL CROSSPROD(x,y,z)
  z = z/norm2(z)
  
  DO i = 1,3
    rotmat(1,i)=x(i)
    rotmat(2,i)=y(i)
    rotmat(3,i)=z(i)
  END DO

  do i=1,n
    ion_pair(:,i) = MATMUL(rotmat,ion_pair(:,i))
  end do

END SUBROUTINE PROJECT_TO_LOCAL_FRAME

SUBROUTINE GET_ATOMS_AROUND_ION_PAIR(ind_ion_pair, ion_pair, n)
  USE parameters, only : pos, natoms
  IMPLICIT NONE
  INTEGER :: n, iat, c, ind_ion_pair(natoms)
  REAL*8 :: ion_pair(3,natoms), coord1(4), coord2(4), d
  REAL*8, PARAMETER :: rcut=3.0d0

  c=7

  DO iat=1,natoms
    IF(.not.any(ind_ion_pair(1:6)==iat)) THEN
      CALL DISTANCE_VECTOR(pos(:,ind_ion_pair(1)),pos(:,iat),coord1)
      CALL DISTANCE_VECTOR(pos(:,ind_ion_pair(5)),pos(:,iat),coord2)
      IF ((coord1(4)<rcut).or.(coord2(4)<rcut)) then
        ion_pair(:,c) = coord1(1:3)
        ind_ion_pair(c) = iat
        c=c+1
      END IF
    END IF
  END DO

  n=c-1

END SUBROUTINE GET_ATOMS_AROUND_ION_PAIR
