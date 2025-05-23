!Computes the number density of water ions at different Z positions

MODULE histogram 
  IMPLICIT NONE 
  INTEGER                                  :: nhist
  INTEGER, ALLOCATABLE                     :: hist(:,:)
  INTEGER                                  :: ntype, zdir
  LOGICAL                                  :: input_is_xyz=.false.
  REAL*8                                   :: zoffset
  REAL*8, parameter                        :: rcut_OTi=2.6
END MODULE histogram 

PROGRAM ZDF
  USE parameters
  USE histogram, only: input_is_xyz
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL(input_is_xyz)

  
  DO frame = 1,nframes
    !
    if (input_is_xyz) then
      CALL READ_XYZ
    else
      CALL READ_ATOM_REDUCED
    end if
    !
    IF ((frame==1).or.(MOD(frame,10*nskip)==0)) THEN
      CALL MAKE_NEIGHBOR_LIST(2,3,6.d0)
      CALL MAKE_NEIGHBOR_LIST(3,1,6.d0)
    END IF
    !
    IF (MOD(frame,nskip)==0) THEN
      !CALL SET_CENTER_OF_MASS_TO_ZERO
      CALL MAKE_HISTOGRAM
    END IF
    !
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : hist, ntype, nhist, input_is_xyz, zdir, zoffset
  USE parameters, ONLY : natoms, nframes, nequil, box,&
                         atype, pos, nskip, neighborlist, maxneighbor
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, ntype, nframes, nequil, nskip, nhist, zdir
  READ(1,*) zoffset

  maxneighbor=60
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(hist(5,nhist)); hist=0
  ALLOCATE(neighborlist(natoms,maxneighbor))

  !Read the index file
  IF (input_is_xyz) THEN
    READ(1,*) box 
    DO i = 1,natoms
      READ(1,fmt = "(I2)") atype(i)  
    END DO
  END IF

  CLOSE(1)
  OPEN(unit  =  1,file  =  pos_file)
 
  hist = 0

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  !Histogram the z coordinate of oxygen atoms. Their are
  !discriminated into 5 types:
  !1: bulk OH-, 2: bulk H3O+, 3: H2O, 4: bridging OH and 5: terminal OH
  USE histogram, ONLY  : hist, nhist, zdir, zoffset, rcut_OTi
  USE parameters, ONLY : pos, box, natoms, atype
  IMPLICIT NONE
  INTEGER          :: i, ind
  REAL*8           :: z,OTi_CN(natoms)
  INTEGER          :: Coordination_Number, OTiCN, typeind 
  INTEGER          :: charge(natoms), get_type_index

  CALL compute_charge_non_smooth(charge)
  
  DO i = 1,natoms

    if (atype(i)==3) then !oxygen atoms only

      OTiCN=Coordination_Number(i,1,rcut_OTi) ! 1 is the Ti type
      typeind = get_type_index(charge(i),OTiCN)
      
      IF (typeind>0) then
        !Applying PBC
        z = pos(zdir,i) + zoffset
        z = z  - nint(z/box(zdir))*box(zdir) 
        z = z + box(zdir)/2. ! for data analysis only

        !Assigning the index to the histogram
        ind = int( z * float(nhist) / box(zdir) ) + 1

        hist(typeind,ind) = hist(typeind,ind) + 1
      END IF
 
    end if

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram, ONLY  : nhist, hist, zdir
  USE parameters, ONLY : box, nframes, nskip
  IMPLICIT NONE
  INTEGER                    :: i
  REAL*8                     :: bin,dV

  OPEN(unit = 2,file = "ZDF_wi.dat")

  WRITE(2,*) '# OH- H3O+ H2O Bridging_OH Terminal_OH'
  bin = box(zdir)/float(nhist)
  dV = product(box)/box(zdir)*bin

    DO i = 1,nhist
      WRITE(2,fmt = "(F10.5, 5(3X, F12.7))"), bin/2.d0 + float(i-1)*bin - box(zdir)/2.d0, &
                           & float(hist(:,i))/(0.03316d0*dV*float(nframes)/float(nskip))
    END DO

  CLOSE(2)

END SUBROUTINE

subroutine compute_charge_non_smooth(charge)
  !Assign formal charge to oxygen atoms. Charge is defined
  !as the O-H coordination number minus 2.
  use parameters, only: natoms, atype, maxneighbor, neighborlist
  implicit none
  integer                    :: iw, ih, min_ind, charge(natoms)
  real*8                     :: Dist
  real*8                     :: mindist, d

  charge = 0
  do ih = 1,natoms
    !
    if (atype(ih)==2) then
      !
      mindist=100.d0
      min_ind=0
      do iw = 1,maxneighbor
        !
        if (neighborlist(ih,iw).eq.0) EXIT
        !
        d=Dist(ih,neighborlist(ih,iw))
        if (d<mindist) then
          mindist=d
          min_ind=neighborlist(ih,iw)
        endif
        !
      end do
      !
      IF (min_ind.ne.0) then
        charge(min_ind) = charge(min_ind) + 1
      ELSE
        WRITE(*,*) 'H atom not bound: ', ih
        WRITE(*,*) 'STOP!!!'
        STOP
      END IF
      !
    end if
    !
  end do
  !
  charge = charge - 2

end subroutine compute_charge_non_smooth

INTEGER FUNCTION get_type_index(charge,OTi_CN)
  !1: bulk OH-, 2: bulk H3O+, 3: H2O, 4: bridging OH and 5: terminal OH
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: charge,OTi_CN

  IF ((charge==0).and.(OTi_CN<2)) then
    get_type_index=3 
  ELSEIF (charge==1) then
    get_type_index=2 
  ELSEIF ((charge==-1).and.(OTi_CN==0)) then
    get_type_index=1
  ELSEIF ((charge==-1).and.(OTi_CN==2)) then
    get_type_index=4
  ELSEIF ((charge==-1).and.(OTi_CN==1)) then
    get_type_index=5
  ELSE
    get_type_index=-1
  END IF

END FUNCTION get_type_index
