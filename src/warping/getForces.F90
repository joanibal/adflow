
subroutine getForces(forces, npts, sps)
  use constants
  use communication, only : myid
  use blockPointers, only : BCData, nDom, nBocos, BCType
  use inputPhysics, only : forcesAsTractions
  use utils, only : setPointers, terminate, EChk
  use surfaceIntegrations, only : integrateSurfaces
  use surfaceFamilies, only : fullfamList
  use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
  use surfaceFamilies, only : familyExchange, BCFamExchange
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(inout) :: forces(3,npts)

  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: sss(3),v2(3),v1(3), qa, sepSensor, Cavitation
  real(kind=realType) :: sepSensorAvg(3)
  real(kind=realType) :: Fp(3), Fv(3), Mp(3), Mv(3), yplusmax, qf(3)
  real(kind=realType) :: localValues(nLocalValues)
  type(zipperMesh), pointer :: zipper
  type(familyexchange), pointer :: exch
  real(kind=realType), dimension(:), pointer :: localPtr
  real(kind=realType), dimension(nCostFunction) :: funcValues
  ! Make sure *all* forces are computed. Sectioning will be done
  ! else-where.

  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)
     localValues = zero
     call integrateSurfaces(localValues, fullFamList)
  end do domains

  if (forcesAsTractions) then
     ! Compute tractions if necessary
     call computeNodalTractions(sps)
  else
     call computeNodalForces(sps)
  end if

  ii = 0
  domains2: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1, nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           ! This is easy, just copy out F or T in continuous ordering.
           do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg, BCData(mm)%inEnd
                 ii = ii + 1
                 if (forcesAsTractions) then
                    Forces(:, ii) = bcData(mm)%Tp(i, j, :) + bcData(mm)%Tv(i, j, :)
                 else
                    Forces(:, ii) = bcData(mm)%F(i, j, :)
                 end if
              end do
           end do
        end if
     end do bocos
  end do domains2

  ! We know must consider additional forces that are required by the
  ! zipper mesh triangles on the root proc.

  ! Pointer for easier reading.
  zipper => zipperMeshes(iBCGroupWalls)
  exch => BCFamExchange(iBCGroupWalls, sps)
  ! No overset present or the zipper isn't allocated nothing to do:
  if (.not. oversetPresent .or. .not. zipper%allocated) then
     return
  end if

  if (.not. forcesAsTractions) then
     ! We have a zipper and regular forces are requested. This is not yet supported.
     call terminate('getForces', 'getForces() is not implmented for zipper meshes and '&
          &'forcesAsTractions=False')
  end if

  ! Loop over each dimension individually since we have a scalar
  ! scatter.
  dimLoop: do iDim=1,3

     call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Copy in the values we already have to the exchange.
     ii = size(LocalPtr)
     localPtr = forces(iDim, 1:ii)

     ! Restore the pointer
     call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Now scatter this to the zipper
     call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
          zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
          zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! The values we need are precisely what is in zipper%localVal
     call vecGetArrayF90(zipper%localVal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Just copy the received data into the forces array. Just on root proc:
     if (myid == 0) then
        forces(iDim, ii+1:ii+size(localPtr)) = localPtr
     end if

     call vecGetArrayF90(zipper%localVal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do dimLoop

end subroutine getForces

subroutine getForces_d(forces, forcesd, npts, sps)

  ! This routine performs the forward mode linearization getForces. It
  ! takes in perturbations defined on bcData(mm)%Fp, bcData(mm)%Fv and
  ! bcData(mm)%area and computes either the nodal forces or nodal
  ! tractions.
  use constants
  use communication, only : myid
  use blockPointers, only : nDom, nBocos, BCData, BCType, nBocos, BCDatad
  use inputPhysics, only : forcesAsTractions
  use surfaceFamilies, only: BCFamExchange, familyExchange
  use utils, only : setPointers, setPointers_d, EChk, terminate
  use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
  use surfaceFamilies, only : familyExchange, BCFamExchange
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(out), dimension(3, npts) :: forces, forcesd
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qa, qad, qf, qfd
  real(kind=realType), dimension(:), pointer :: localPtr, localPtrd
  type(zipperMesh), pointer :: zipper
  type(familyexchange), pointer :: exch

  if (forcesAsTractions) then
     call computeNodalTractions_d(sps)
  else
     call computeNodalForces_d(sps)
  end if

  ! Extract the values out into the output derivative array
  ii = 0
  domains2: do nn=1,nDom
     call setPointers_d(nn, 1_intType, sps)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1, nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           ! This is easy, just copy out F or T in continuous ordering.
           do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg, BCData(mm)%inEnd
                 ii = ii + 1
                 if (forcesAsTractions) then
                    Forcesd(:, ii) = bcDatad(mm)%Tp(i, j, :) + bcDatad(mm)%Tv(i, j, :)
                 else
                    Forcesd(:, ii) = bcDatad(mm)%F(i, j, :)
                 end if
              end do
           end do
        end if
     end do bocos
  end do domains2

  ! We know must consider additional forces that are required by the
  ! zipper mesh triangles on the root proc.

  ! Pointer for easier reading.
  zipper => zipperMeshes(iBCGroupWalls)
  exch => BCFamExchange(iBCGroupWalls, sps)
  ! No overset present or the zipper isn't allocated nothing to do:
  if (.not. oversetPresent .or. .not. zipper%allocated) then
     return
  end if

  if (.not. forcesAsTractions) then
     ! We have a zipper and regular forces are requested. This is not yet supported.
     call terminate('getForces', 'getForces() is not implmented for zipper meshes and '&
          &'forcesAsTractions=False')
  end if

  ! Loop over each dimension individually since we have a scalar
  ! scatter.
  dimLoop: do iDim=1,3

     call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Copy in the values we already have to the exchange.
     ii = size(LocalPtr)
     localPtr = forcesd(iDim, 1:ii)

     ! Restore the pointer
     call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Now scatter this to the zipper
     call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
          zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
          zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! The values we need are precisely what is in zipper%localVal
     call vecGetArrayF90(zipper%localVal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Just copy the received data into the forces array. Just on root proc:
     if (myid == 0) then
        forcesd(iDim, ii+1:ii+size(localPtr)) = localPtr
     end if

     call vecGetArrayF90(zipper%localVal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do dimLoop

end subroutine getForces_d

subroutine getForces_b(forcesd, npts, sps)

  ! This routine performs the reverse of getForces. It takes in
  ! forces_b and perfroms the reverse of the nodal averaging procedure
  ! in getForces to compute bcDatad(mm)%Fp, bcDatad(mm)%Fv and
  ! bcDatad(mm)%area.
  use constants
  use communication, only : myid
  use blockPointers, only : nDom, nBocos, BCData, BCType, nBocos, BCDatad
  use inputPhysics, only : forcesAsTractions
  use surfaceFamilies, only: BCFamExchange, familyExchange
  use communication
  use utils, only : EChk, setPointers, setPointers_d
  use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
  use surfaceFamilies, only : familyExchange, BCFamExchange
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(inout) :: forcesd(3, npts)
  integer(kind=intType) :: mm, nn, i, j, ii, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  type(zipperMesh), pointer :: zipper
  type(familyexchange), pointer :: exch
  real(kind=realType), dimension(:), pointer :: localPtr
  real(kind=realType), dimension(3, npts) :: forces

  ! Run nonlinear code to make sure that all intermediate values are
  ! updated.
  call getForces(forces, npts, sps)

  ! We now must consider additional forces that are required by the
  ! zipper mesh triangles on the root proc.

  ! Pointer for easier reading.
  zipper => zipperMeshes(iBCGroupWalls)
  exch => BCFamExchange(iBCGroupWalls, sps)
  ! No overset present or the zipper isn't allocated nothing to do:
  zipperReverse: if (oversetPresent .and. zipper%allocated) then

     ! Loop over each dimension individually since we have a scalar
     ! scatter.
     dimLoop: do iDim=1,3

        call vecGetArrayF90(zipper%localVal, localPtr, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ii = exch%nNodes
        ! Just copy the received data into the forces array. Just on root proc:
        if (myid == 0) then
           do i=1, size(localPtr)
              localPtr(i) = forcesd(iDim, ii+i)
              forcesd(iDim, ii+i) = zero
           end do
        end if

        call vecGetArrayF90(zipper%localVal, localPtr, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Zero the vector we are scatting into:
        call VecSet(exch%nodeValLocal, zero, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Scatter values from the root using the zipper scatter.
        call VecScatterBegin(zipper%scatter, zipper%localVal, &
             exch%nodeValLocal, ADD_VALUES, SCATTER_REVERSE, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call VecScatterEnd(zipper%scatter, zipper%localVal, &
             exch%nodeValLocal, ADD_VALUES, SCATTER_REVERSE, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Accumulate the scatted values onto forcesd
        ii = size(localPtr)
        forcesd(iDim, 1:ii) = forcesd(iDim, 1:ii) + localPtr

        ! Restore the pointer
        call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
        call EChk(ierr,__FILE__,__LINE__)

     end do dimLoop
  end if zipperReverse

  ! Set the incoming derivative values
  ii = 0
  domains2: do nn=1,nDom
     call setPointers_d(nn, 1_intType, sps)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1, nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           ! This is easy, just copy out F or T in continuous ordering.
           do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg, BCData(mm)%inEnd
                 ii = ii + 1
                 if (forcesAsTractions) then
                    bcDatad(mm)%Tp(i, j, :) = forcesd(:, ii)
                    bcDatad(mm)%Tv(i, j, :) = forcesd(:, ii)
                 else
                    bcDatad(mm)%F(i, j, :) = forcesd(:, ii)
                 end if
              end do
           end do
        end if
     end do bocos
  end do domains2

  if (.not. forcesAsTractions) then
     ! For forces, we can accumulate the nodal seeds on the Fp and Fv
     ! values. The area seed is zeroed.
     call computeNodalForces_b(sps)
  else
     call computeNodalTractions_b(sps)
  end if

end subroutine getForces_b


subroutine surfaceCellCenterToNode(exch)

  use constants
  use blockPointers, only : BCData, nDom, nBocos, BCType
  use surfaceFamilies, only : familyExchange
  use utils, only : setPointers, EChk
  use sorting, only : famInList
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif

  type(familyExchange) :: exch
  integer(kind=intType) ::  sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qv
  real(kind=realType), dimension(:), pointer :: localPtr

  ! We assume that normalization factor is already computed
  sps = exch%sps
  call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  localPtr = zero

  ! ii is the running counter through the pointer array.
  ii = 0
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        famInclude: if (famInList(BCData(mm)%famID, exch%famList)) then
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           do j=0,nj-2
              do i=0,ni-2
                 ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                 ! and always starts at one
                 qv = fourth * BCData(mm)%cellVal(i+1, j+1)
                 ind(1) = ii + (j  )*ni + i + 1
                 ind(2) = ii + (j  )*ni + i + 2
                 ind(3) = ii + (j+1)*ni + i + 2
                 ind(4) = ii + (j+1)*ni + i + 1
                 do jj=1,4
                    localPtr(ind(jj)) = localPtr(ind(jj)) + qv
                 end do
              end do
           end do
           ii = ii + ni*nj
        end if famInclude
     end do
  end do

  call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Globalize the current face based value
  call vecSet(exch%nodeValGlobal, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
       exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
       exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now divide by the weighting. We can do this with a vecpointwisemult
  call vecPointwiseMult(exch%nodeValGlobal, exch%nodeValGlobal, &
       exch%sumGlobal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Push back to the local values
  call VecScatterBegin(exch%scatter, exch%nodeValGlobal, &
       exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, exch%nodeValGlobal, &
  exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ii = 0
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        famInclude2: if (famInList(BCData(mm)%famID, exch%famList)) then
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           do j=1,nj
              do i=1,ni
                 ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                 ! and always starts at one
                 ii = ii + 1
                 BCData(mm)%nodeVal(i, j) = localPtr(ii)
              end do
           end do
        end if famInclude2
     end do
  end do

  call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine surfaceCellCenterToNode


subroutine surfaceCellCenterToNode_d(exch)

   use constants
   use blockPointers, only : BCData, BCDatad, nDom, nBocos, BCType
   use surfaceFamilies, only : familyExchange
   use utils, only : setPointers_d, EChk
   use sorting, only : famInList
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
   use petsc
   implicit none
#else
   implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif

   type(familyExchange) :: exch
   integer(kind=intType) ::  sps
   integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
   real(kind=realType) :: qv, qvd
   real(kind=realType), dimension(:), pointer :: nodeValLocPtr, nodeValLocPtrd
   Vec tmp


   sps = exch%sps

   call VecDuplicate(exch%sumGlobal, tmp, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! We assume that normalization factor is already computed
   call vecGetArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%nodeValLocal_d, nodeValLocPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   nodeValLocPtr = zero
   nodeValLocPtrd = zero

   ! ii is the running counter through the pointer array.
   ii = 0
   do nn=1, nDom
      call setPointers_d(nn, 1_intType, sps)
      do mm=1, nBocos
         famInclude: if (famInList(BCData(mm)%famID, exch%famList)) then
            iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
            jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
            ni = iEnd - iBeg + 1
            nj = jEnd - jBeg + 1
            do j=0,nj-2
               do i=0,ni-2
                  ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                  ! and always starts at one
                  qv = fourth * BCData(mm)%cellVal(i+1, j+1)
                  qvd = fourth * BCDatad(mm)%cellVal(i+1, j+1)
                  ind(1) = ii + (j  )*ni + i + 1
                  ind(2) = ii + (j  )*ni + i + 2
                  ind(3) = ii + (j+1)*ni + i + 2
                  ind(4) = ii + (j+1)*ni + i + 1
                  do jj=1,4
                     nodeValLocPtr(ind(jj)) = nodeValLocPtr(ind(jj)) + qv
                     nodeValLocPtrd(ind(jj)) = nodeValLocPtrd(ind(jj)) + qvd

                  end do
               end do
            end do
            ii = ii + ni*nj
         end if famInclude
      end do
   end do

   call vecRestoreArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   call vecRestoreArrayF90(exch%nodeValLocal_d, nodeValLocPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! Globalize the current face based value
   call vecSet(exch%nodeValGlobal, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
        exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
        exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)


     ! Globalize the current force derivative
   call vecSet(exch%nodeValGlobal_d, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterBegin(exch%scatter, exch%nodeValLocal_d, &
        exch%nodeValGlobal_d, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValLocal_d, &
        exch%nodeValGlobal_d, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

    ! The product rule here: (since we are multiplying)
    ! nodeValGlobal = nodeValGlobal * invArea
    ! nodeValGlobald = nodeValGlobald*invArea + nodeValGlobal*invAread

     ! First term:  nodeValGlobald = nodeValGlobald*invArea
   call vecPointwiseMult(exch%nodeValGlobal_d, exch%nodeValGlobal_d, &
   exch%sumGlobal, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! Second term:, tmp = nodeValGlobal*invAread
   call vecPointwiseMult(tmp, exch%nodeValGlobal, exch%sumGlobal_d, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Sum the second term into the first
   call VecAXPY(exch%nodeValGlobal_d, one, tmp, ierr)
   call EChk(ierr,__FILE__,__LINE__)




   ! Now divide by the weighting. We can do this with a vecpointwisemult
   call vecPointwiseMult(exch%nodeValGlobal, exch%nodeValGlobal, &
        exch%sumGlobal, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Push back to the local values
   call VecScatterBegin(exch%scatter, exch%nodeValGlobal, &
        exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValGlobal, &
   exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! Push back to the local derivative values
   call VecScatterBegin(exch%scatter, exch%nodeValGlobal_d, &
      exch%nodeValLocal_d, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValGlobal_d, &
      exch%nodeValLocal_d, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%nodeValLocal_d, nodeValLocPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ii = 0
   do nn=1, nDom
      call setPointers_d(nn, 1_intType, sps)
      do mm=1, nBocos
         famInclude2: if (famInList(BCData(mm)%famID, exch%famList)) then
            iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
            jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

            ni = iEnd - iBeg + 1
            nj = jEnd - jBeg + 1
            do j=1,nj
               do i=1,ni
                  ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                  ! and always starts at one
                  ii = ii + 1
                  BCData(mm)%nodeVal(i, j) = nodeValLocPtr(ii)
                  BCDatad(mm)%nodeVal(i, j) = nodeValLocPtrd(ii)
               end do
            end do
         end if famInclude2
      end do
   end do

   call vecRestoreArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecRestoreArrayF90(exch%nodeValLocal_d, nodeValLocPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)


  call VecDestroy(tmp, ierr)
  call EChk(ierr,__FILE__,__LINE__)

 end subroutine surfaceCellCenterToNode_d

 subroutine surfaceCellCenterToNode_b(exch)

  ! This routine performs the reverse of surfaceCellCenterToNode_d. It
  ! takes in BCDatad(mm)%nodeVal and perfroms the reverse of the
  ! nodal averaging procedure in getForces to compute BCDatad(mm)%cellVal
  ! and exch%sumGlobal_b.

   use constants
   use blockPointers, only : BCData, BCDatad, nDom, nBocos, BCType
   use surfaceFamilies, only : familyExchange
   use utils, only : setPointers, setPointers_b, EChk
   use sorting, only : famInList
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
   use petsc
   implicit none
#else
   implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif

   type(familyExchange) :: exch
   integer(kind=intType) ::  sps
   integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
   real(kind=realType) :: qv, qv_b
   real(kind=realType), dimension(:), pointer :: localPtr, nodeValLocPtr, nodeValLocPtr_b
   Vec tmp
   ! We assume that normalization factor is already computed

   call VecDuplicate(exch%sumGlobal, tmp, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   sps = exch%sps
   call vecGetArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)
   nodeValLocPtr = zero

   ! ii is the running counter through the pointer array.
   ii = 0
   do nn=1, nDom
      call setPointers(nn, 1_intType, sps)
      do mm=1, nBocos
         famInclude: if (famInList(BCData(mm)%famID, exch%famList)) then
            iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
            jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
            ni = iEnd - iBeg + 1
            nj = jEnd - jBeg + 1
            do j=0,nj-2
               do i=0,ni-2
                  ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                  ! and always starts at one
                  qv = fourth * BCData(mm)%cellVal(i+1, j+1)
                  ind(1) = ii + (j  )*ni + i + 1
                  ind(2) = ii + (j  )*ni + i + 2
                  ind(3) = ii + (j+1)*ni + i + 2
                  ind(4) = ii + (j+1)*ni + i + 1
                  do jj=1,4
                     nodeValLocPtr(ind(jj)) = nodeValLocPtr(ind(jj)) + qv
                  end do
               end do
            end do
            ii = ii + ni*nj
         end if famInclude
      end do
   end do

   call vecRestoreArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Globalize the current face based value
   call vecSet(exch%nodeValGlobal, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
        exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
        exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! ====================
   ! Do the reverse pass:
   ! ====================

   call vecGetArrayF90(exch%nodeValLocal_b, nodeValLocPtr_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ii = 0
   do nn=1, nDom
      call setPointers_b(nn, 1_intType, sps)
      do mm=1, nBocos
         famInclude2: if (famInList(BCData(mm)%famID, exch%famList)) then
            iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
            jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

            ni = iEnd - iBeg + 1
            nj = jEnd - jBeg + 1
            do j=1,nj
               do i=1,ni
                  ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                  ! and always starts at one
                  ii = ii + 1
                  nodeValLocPtr_b(ii) = BCDatad(mm)%nodeVal(i, j)
               end do
            end do
         end if famInclude2
      end do
   end do

   call vecRestoreArrayF90(exch%nodeValLocal_b, nodeValLocPtr_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! --- Globalize the local values ---
   call vecSet(exch%nodeValGlobal_b, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterBegin(exch%scatter, exch%nodeValLocal_b, &
        exch%nodeValGlobal_b, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValLocal_b, &
        exch%nodeValGlobal_b, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! ---  account for the weighting process ---
   ! in the forward pass
   ! exch%nodeValGlobal = exch%nodeValGlobal *  exch%sumGlobal
   ! becomes...
   ! exch%nodeValGlobal_b += exch%nodeValGlobal_b * exch%sumGlobal
   ! exch%sumGlobal_b += exch%nodeValGlobal* exch%nodeValGlobal_b
         ! (note non weighted value used ^)
   ! in the reverse pass

   call vecPointwiseMult(tmp, exch%nodeValGlobal, exch%nodeValGlobal_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   call vecPointwiseMult(exch%nodeValGlobal_b,  exch%sumGlobal, exch%nodeValGlobal_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Accumulate
   call vecAXPY(exch%sumGlobal_b, one, tmp, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   call vecGetArrayF90(exch%nodeValGlobal, localPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecRestoreArrayF90(exch%nodeValGlobal, localPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! --- Now communicate exch%nodeValGlobal_b back to the local patches ---

   call VecScatterBegin(exch%scatter, exch%nodeValGlobal_b, &
   exch%nodeValLocal_b, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValGlobal_b, &
   exch%nodeValLocal_b, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Copy the values into patches
   call vecGetArrayF90(exch%nodeValLocal_b, nodeValLocPtr_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! ii is the running counter through the pointer array.
   ii = 0
   do nn=1, nDom
      call setPointers_b(nn, 1_intType, sps)
      do mm=1, nBocos
         famInclude3: if (famInList(BCData(mm)%famID, exch%famList)) then
            iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
            jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
            ni = iEnd - iBeg + 1
            nj = jEnd - jBeg + 1
            do j=0,nj-2
               do i=0,ni-2


                  ind(1) = ii + (j  )*ni + i + 1
                  ind(2) = ii + (j  )*ni + i + 2
                  ind(3) = ii + (j+1)*ni + i + 2
                  ind(4) = ii + (j+1)*ni + i + 1
                  qv_b = zero
                  do jj=1,4
                     qv_b = qv_b + nodeValLocPtr_b(ind(jj))
                  end do
                  qv_b = qv_b*fourth

                  ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                  ! and always starts at one
                  BCDatad(mm)%cellVal(i+1, j+1) = &
                                 BCDatad(mm)%cellVal(i+1, j+1) + qv_b

               end do
            end do
            ii = ii + ni*nj
         end if famInclude3
      end do
   end do

   call vecRestoreArrayF90(exch%nodeValLocal_b, nodeValLocPtr_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecDestroy(tmp, ierr)
   call EChk(ierr,__FILE__,__LINE__)


end subroutine surfaceCellCenterToNode_b


subroutine computeWeighting(exch)

  use constants
  use blockPointers, only : BCData, nDom, nBocos, BCType
  use surfaceFamilies, only : familyExchange
  use utils, only : setPointers, EChk
  use sorting, only : famInList
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif
  type(familyExchange) :: exch
  integer(kind=intType) ::  sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qf, qa
  real(kind=realType), dimension(:), pointer :: nodeValLocPtr, sumGlobalPtr

  sps = exch%sps

  call vecGetArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  nodeValLocPtr = zero
  ! ii is the running counter through the pointer array.
  ii = 0
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        famInclude: if (famInList(BCData(mm)%famID, exch%famList)) then
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           do j=0,nj-2
              do i=0,ni-2

                 ! Scatter a quarter of the face value to each node:
                 ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                 ! and always starts at one
                 qa = fourth*BCData(mm)%cellVal(i+1, j+1)
                 ind(1) = ii + (j  )*ni + i + 1
                 ind(2) = ii + (j  )*ni + i + 2
                 ind(3) = ii + (j+1)*ni + i + 2
                 ind(4) = ii + (j+1)*ni + i + 1
                 do jj=1,4
                    nodeValLocPtr(ind(jj)) = nodeValLocPtr(ind(jj)) + qa
                 end do
              end do
           end do
           ii = ii + ni*nj
        end if famInclude
     end do
  end do

  call vecRestoreArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Globalize the face value
  call vecSet(exch%sumGlobal, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now compute the inverse of the weighting so that we can multiply
  ! instead of dividing. Note that we check dividing by zero and just
  ! set those to zero.

  call vecGetArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  do i=1, size(sumGlobalPtr)
     if (sumGlobalPtr(i) == zero) then
        sumGlobalPtr(i) = zero
     else
        sumGlobalPtr(i) = one/sumGlobalPtr(i)
     end if
  end do

  call vecRestoreArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine computeWeighting

subroutine computeWeighting_d(exch)

   use constants
   ! use blockPointers, only : BCData, nDom, nBocos, BCType
   use blockPointers, only : nDom, nBocos, BCData, BCType, nBocos, BCDatad

   use surfaceFamilies, only : familyExchange
   use utils, only : setPointers, setPointers_d, EChk
   use sorting, only : famInList
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
   use petsc
   implicit none
#else
   implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif
   type(familyExchange) :: exch
   integer(kind=intType) ::  sps
   integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
   real(kind=realType) :: qa, qad
   real(kind=realType), dimension(:), pointer :: nodeValLocPtr, nodeValLocPtrd
   real(kind=realType), dimension(:), pointer :: sumGlobalPtr, sumGlobalPtrd


   sps = exch%sps



   call vecGetArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%nodeValLocal_d, nodeValLocPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   nodeValLocPtrd = zero
   nodeValLocPtr = zero
   ! ii is the running counter through the pointer array.
   ii = 0
   do nn=1, nDom
      call setPointers_d(nn, 1_intType, sps)
      do mm=1, nBocos
         famInclude: if (famInList(BCData(mm)%famID, exch%famList)) then
            iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
            jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
            ni = iEnd - iBeg + 1
            nj = jEnd - jBeg + 1
            do j=0,nj-2
               do i=0,ni-2

                  ! Scatter a quarter of the face value to each node:
                  ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                  ! and always starts at one
                  qa = fourth*BCData(mm)%cellVal(i+1, j+1)
                  qad = fourth*BCDatad(mm)%cellVal(i+1, j+1)

                  ind(1) = ii + (j  )*ni + i + 1
                  ind(2) = ii + (j  )*ni + i + 2
                  ind(3) = ii + (j+1)*ni + i + 2
                  ind(4) = ii + (j+1)*ni + i + 1
                  do jj=1,4
                     nodeValLocPtr(ind(jj)) = nodeValLocPtr(ind(jj)) + qa
                     nodeValLocPtrd(ind(jj)) = nodeValLocPtrd(ind(jj)) + qad
                  end do
               end do
            end do
            ii = ii + ni*nj
         end if famInclude
      end do
   end do

   call vecRestoreArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecRestoreArrayF90(exch%nodeValLocal_d, nodeValLocPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Globalize the area
   call vecSet(exch%sumGlobal, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
        exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
        exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Globalize the area derivative
   call vecSet(exch%sumGlobal_d, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterBegin(exch%scatter, exch%nodeValLocal_d, &
        exch%sumGlobal_d, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValLocal_d, &
        exch%sumGlobal_d, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Now compute the inverse of the weighting so that we can multiply
   ! instead of dividing. Here we need the original value too:

   call vecGetArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%sumGlobal_d, sumGlobalPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   sumGlobalPtrd = -(sumGlobalPtrd/sumGlobalPtr**2)
   sumGlobalPtr = one/sumGlobalPtr

   call vecGetArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecRestoreArrayF90(exch%sumGlobal_d, sumGlobalPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

end subroutine computeWeighting_d


subroutine computeWeighting_b(exch)
   ! in: exch%sumGlobal_d
   ! out: BCDatad(mm)%cellVal

   use constants
   use blockPointers, only : BCData, BCDatad, nDom, nBocos, BCType
   use surfaceFamilies, only : familyExchange
   use utils, only : setPointers, setPointers_b, EChk
   use sorting, only : famInList
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
   use petsc
   implicit none
#else
   implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif
   type(familyExchange) :: exch
   integer(kind=intType) ::  sps
   integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
   real(kind=realType) :: qa, qa_b
   real(kind=realType), dimension(:), pointer :: nodeValLocPtr,   sumGlobalPtr,&
                                                 nodeValLocPtr_b, sumGlobalPtr_b

   sps = exch%sps

   call computeWeighting(exch)

   call vecGetArrayF90(exch%sumGlobal_b, sumGlobalPtr_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! Keep in mind sumGlobalPtr points to sumGlobal which already has
   ! been inversed so we just multiply.
   sumGlobalPtr_b = -sumGlobalPtr_b*sumGlobalPtr**2

   call vecRestoreArrayF90(exch%sumGlobal_b, sumGlobalPtr_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecRestoreArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)



   ! Push back to the local patches
   call VecScatterBegin(exch%scatter, exch%sumGlobal_b, &
         exch%nodeValLocal_b, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%sumGlobal_b, &
   exch%nodeValLocal_b, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%nodeValLocal_b, nodeValLocPtr_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! ii is the running counter through the pointer array.
   ii = 0
   do nn=1, nDom
      call setPointers_b(nn, 1_intType, sps)
      do mm=1, nBocos
         famInclude2: if (famInList(BCData(mm)%famID, exch%famList)) then
            iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
            jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
            ni = iEnd - iBeg + 1
            nj = jEnd - jBeg + 1
            do j=0,nj-2
               do i=0,ni-2

                  ! Note: No +iBeg, and +jBeg becuase cellVal is a pointer
                  ! and always starts at one
                  ind(1) = ii + (j  )*ni + i + 1
                  ind(2) = ii + (j  )*ni + i + 2
                  ind(3) = ii + (j+1)*ni + i + 2
                  ind(4) = ii + (j+1)*ni + i + 1
                  qa_b = zero
                  do jj=1,4
                     qa_b = qa_b + nodeValLocPtr_b(ind(jj))
                  end do
                  qa_b = fourth*qa_b
                  BCDatad(mm)%cellVal(i+1, j+1) = &
                                           BCDatad(mm)%cellVal(i+1, j+1) + qa_b
               end do
            end do
            ii = ii + ni*nj
         end if famInclude2
      end do
   end do


   call vecRestoreArrayF90(exch%nodeValLocal_b, nodeValLocPtr_b, ierr)
   call EChk(ierr,__FILE__,__LINE__)


 end subroutine computeWeighting_b


subroutine computeNodalTractions(sps)
  use constants
  use blockPointers, only : BCData, nDom, nBocos, BCType
  use surfaceFamilies, only : BCFamExchange, familyExchange
  use utils, only : setPointers, EChk, isWallType
  implicit none

  integer(kind=intType), intent(in) ::  sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qf, qa
  real(kind=realType), dimension(:), pointer :: localPtr
  type(familyExchange), pointer :: exch

  ! Set the pointer to the wall exchange:
  exch => BCfamExchange(iBCGroupWalls, sps)

  ! Set the weighting factors. In this case, area
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        bocoType1: if(isWallType(BCType(mm))) then
           BCData(mm)%cellVal => BCData(mm)%area(:, :)
        end if bocoType1
     end do
  end do

  call computeWeighting(exch)

  FpFvLoop: do iDim=1, 6
     ! ii is the running counter through the pointer array.
     ii = 0
     do nn=1, nDom
        call setPointers(nn, 1_intType, sps)
        do mm=1, nBocos
           bocoType2: if(isWallType(BCType(mm))) then
              if (iDim <= 3) then
                 BCData(mm)%cellVal => BCData(mm)%Fp(:, :, iDim)
                 BCData(mm)%nodeVal => BCData(mm)%Tp(:, :, iDim)
              else
                 BCData(mm)%cellVal => BCData(mm)%Fv(:, :, iDim-3)
                 BCData(mm)%nodeVal => BCData(mm)%Tv(:, :, iDim-3)
              end if
           end if bocoType2
        end do
     end do

     call surfaceCellCenterToNode(exch)

  end do FpFVLoop

end subroutine computeNodalTractions

subroutine computeNodalTractions_d(sps)

  ! Forward mode lineariation of nodal tractions

  use constants
  use blockPointers, only : nDom, nBocos, BCData, BCType, nBocos, BCDatad
  use inputPhysics, only : forcesAsTractions
  use surfaceFamilies, only: BCFamExchange, familyExchange
  use utils, only : setPointers, setPointers_d, EChk, isWallType
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif
  integer(kind=intType), intent(in) :: sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qa, qad, qf, qfd
  real(kind=realType), dimension(:), pointer :: localPtr, localPtrd
  type(familyExchange), pointer :: exch

  exch => BCFamExchange(iBCGroupWalls, sps)


  ! Set the weighting factors. In this case, area
   do nn=1, nDom
      call setPointers_d(nn, 1_intType, sps)
      do mm=1, nBocos
         bocoType1: if(isWallType(BCType(mm))) then
            BCData(mm)%cellVal => BCData(mm)%area(:, :)
            BCDatad(mm)%cellVal => BCDatad(mm)%area(:, :)
         end if bocoType1
      end do
   end do

   call computeWeighting_d(exch)

  FpFvLoop: do iDim=1, 6
   do nn=1, nDom
      call setPointers_d(nn, 1_intType, sps)
      do mm=1, nBocos
         bocoType2: if(isWallType(BCType(mm))) then
            if (iDim <= 3) then
               BCData(mm)%cellVal => BCData(mm)%Fp(:, :, iDim)
               BCData(mm)%nodeVal => BCData(mm)%Tp(:, :, iDim)

               BCDatad(mm)%cellVal => BCDatad(mm)%Fp(:, :, iDim)
               BCDatad(mm)%nodeVal => BCDatad(mm)%Tp(:, :, iDim)

            else
               BCData(mm)%cellVal => BCData(mm)%Fv(:, :, iDim-3)
               BCData(mm)%nodeVal => BCData(mm)%Tv(:, :, iDim-3)

               BCDatad(mm)%cellVal => BCDatad(mm)%Fv(:, :, iDim-3)
               BCDatad(mm)%nodeVal => BCDatad(mm)%Tv(:, :, iDim-3)
            end if
         end if bocoType2
      end do
   end do

   call surfaceCellCenterToNode_d(exch)

  end do FpFVLoop
end subroutine computeNodalTractions_d


subroutine computeNodalTractions_b(sps)

   ! This routine performs the reverse of computeNodalTractions. It
   ! takes in bcDatad%Tv and bcDatad%Tp and perfroms the reverse of the
   ! nodal averaging procedure in getForces to compute bcDatad(mm)%Fp,
   ! bcDatad(mm)%Fv and bcDatad(mm)%area.

   use constants
   use blockPointers, only : nDom, nBocos, BCData, BCType, nBocos, BCDatad
   use inputPhysics, only : forcesAsTractions
   use surfaceFamilies, only: BCFamExchange, familyExchange
   use communication
   use utils, only : EChk, setPointers, setPointers_d, setPointers_b, isWallType
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
   use petsc
   implicit none
#else
   implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif

   integer(kind=intType), intent(in) :: sps
   integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
   real(kind=realType) :: qf_b, qf, qa, qa_b
   real(kind=realType), dimension(:), pointer :: localPtr, localPtr_b, nodeValLocPtr_b, sumGlobalPtr_b
   type(familyExchange), pointer :: exch
   Vec tmp, T_b

   ! For better readibility
   exch => BCFamExchange(iBCGroupWalls, sps)


   call vecSet(exch%sumGlobal_b, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! For tractions it's (a lot) more difficult becuase we have to do
   ! the scatter/gather operation.

   ! ==================================
   !  Recompute the dual area
   ! ==================================


    ! Set the weighting factors. In this case, area
    do nn=1, nDom
       call setPointers(nn, 1_intType, sps)
       do mm=1, nBocos
          iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
          jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

          bocoType1: if(isWallType(BCType(mm))) then
             BCData(mm)%cellVal => BCData(mm)%area(:, :)
          end if bocoType1
       end do
    end do

    call computeWeighting(exch)

    FpFvLoop: do iDim=1, 6
    ! ii is the running counter through the pointer array.
    ii = 0
    do nn=1, nDom
       call setPointers_b(nn, 1_intType, sps)
       do mm=1, nBocos
          if(isWallType(BCType(mm))) then
             if (iDim <= 3) then
                BCData(mm)%cellVal => BCData(mm)%Fp(:, :, iDim)
                BCData(mm)%nodeVal => BCData(mm)%Tp(:, :, iDim)

                BCDatad(mm)%cellVal => BCDatad(mm)%Fp(:, :, iDim)
                BCDatad(mm)%nodeVal => BCDatad(mm)%Tp(:, :, iDim)
             else
                BCData(mm)%cellVal => BCData(mm)%Fv(:, :, iDim-3)
                BCData(mm)%nodeVal => BCData(mm)%Tv(:, :, iDim-3)

                BCDatad(mm)%cellVal => BCDatad(mm)%Fv(:, :, iDim-3)
                BCDatad(mm)%nodeVal => BCDatad(mm)%Tv(:, :, iDim-3)
             end if
          end if
       end do
    end do

    call surfaceCellCenterToNode_b(exch)

   end do FpFVLoop

   ! Finish the dual area sensitivity.
    do nn=1, nDom
       call setPointers_b(nn, 1_intType, sps)
       do mm=1, nBocos
          if(isWallType(BCType(mm))) then
            BCData(mm)%cellVal => BCData(mm)%area(:, :)
            BCDatad(mm)%cellVal => BCDatad(mm)%area(:, :)
         end if
       end do
    end do

    call computeWeighting_b(exch)

 end subroutine computeNodalTractions_b

subroutine computeNodalForces(sps)

  ! This subroutine averages the cell based forces and tractions to
  ! node based values. There is no need for communication since we are
  ! simplying summing a quarter of each value to each corner.

  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData
  use utils, only : setPointers
  implicit none

  integer(kind=intType), intent(in) ::  sps

  integer(kind=intType) :: mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: qf(3)

  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           BCData(mm)%F = zero
           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 qf = fourth*(BCData(mm)%Fp(i,j,:) + BCData(mm)%Fv(i,j,:))
                 BCData(mm)%F(i  , j,   :) = BCData(mm)%F(i  , j,   :) + qf
                 BCData(mm)%F(i-1, j,   :) = BCData(mm)%F(i-1, j  , :) + qf
                 BCData(mm)%F(i  , j-1, :) = BCData(mm)%F(i  , j-1, :) + qf
                 BCData(mm)%F(i-1, j-1, :) = BCData(mm)%F(i-1, j-1, :) + qf
              end do
           end do
        end if
     end do
  end do
end subroutine computeNodalForces

subroutine computeNodalForces_d(sps)

  ! Forward mode linearization of nodalForces

  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
  use utils, only : setPointers
  implicit none

  integer(kind=intType), intent(in) ::  sps

  integer(kind=intType) :: mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: qfd(3)

  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           BCDatad(mm)%F = zero
           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 qfd = fourth*(BCDatad(mm)%Fp(i,j,:) + BCDatad(mm)%Fv(i,j,:))
                 BCDatad(mm)%F(i  , j,   :) = BCDatad(mm)%F(i  , j,   :) + qfd
                 BCDatad(mm)%F(i-1, j,   :) = BCDatad(mm)%F(i-1, j  , :) + qfd
                 BCDatad(mm)%F(i  , j-1, :) = BCDatad(mm)%F(i  , j-1, :) + qfd
                 BCDatad(mm)%F(i-1, j-1, :) = BCDatad(mm)%F(i-1, j-1, :) + qfd
              end do
           end do
        end if
     end do
  end do
end subroutine computeNodalForces_d

subroutine computeNodalForces_b(sps)

  ! Reverse mode linearization of nodalForces

  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
  use utils, only : setPointers_d
  implicit none

  integer(kind=intType), intent(in) ::  sps

  integer(kind=intType) :: mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: qf_b(3)

  domains: do nn=1,nDom
     call setPointers_d(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           BCDatad(mm)%F = zero
           do j=jBeg, jEnd
              do i=iBeg, iEnd

                 qf_b = fourth*(BCDatad(mm)%F(i, j, :) + BCdatad(mm)%F(i-1, j, :) + &
                      BCDatad(mm)%F(i, j-1, :) + BCDatad(mm)%F(i-1, j-1, :))

                 ! Fp and Fv are face-based values
                 BCDatad(mm)%Fp(i, j, :) = BCDatad(mm)%Fp(i, j, :) + qf_b
                 BCDatad(mm)%Fv(i, j, :) = BCDatad(mm)%Fv(i, j, :) + qf_b
              end do
           end do
        end if
     end do
  end do domains
end subroutine computeNodalForces_b

subroutine getHeatFlux(hflux, npts, sps)
  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData
  use surfaceFamilies, only : BCFamExchange, familyExchange, &
       zeroCellVal, zeroNodeVal, fullfamList
  use surfaceIntegrations, only : integrateSurfaces

  use utils, only : setPointers
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(out) :: hflux(npts)
  real(kind=realType) :: localValues(nLocalValues)

  integer(kind=intType) :: mm, nn, i, j, ii
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  type(familyExchange), pointer :: exch

  ! Set the pointer to the wall exchange:
  exch => BCFamExchange(iBCGroupWalls, sps)

  domains: do nn=1,nDom
      call setPointers(nn, 1_intType, sps)
      localValues = zero
      call integrateSurfaces(localValues, fullFamList)
   end do domains


  ! Set the weighting factors. In this case, area
  ! the adiabatic wall cells have a weighting factor of zero
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos

        bocoType1: if (BCType(mm) == NSWallIsoThermal) then
           BCData(mm)%cellVal => BCData(mm)%area(:, :)
        else if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic) then
           BCData(mm)%cellVal => zeroCellVal
           BCData(mm)%nodeVal => zeroNodeVal
        end if bocoType1
     end do
  end do

  call computeWeighting(exch)


  ! compute the nodal values from the cell center values and the weighting

  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
        bocoType2: if (BCType(mm) == NSWallIsoThermal) then
           BCData(mm)%cellVal => BCData(mm)%cellHeatFlux(:, :)
           BCData(mm)%nodeVal => BCData(mm)%nodeHeatFlux(:, :)
        end if bocoType2
     end do
  end do

  call surfaceCellCenterToNode(exch)

  ! Now extract into the flat array:
  ii = 0
  do nn=1,nDom
     call setPointers(nn,1_intType,sps)

     ! Loop over the number of viscous boundary subfaces of this block.
     ! According to preprocessing/viscSubfaceInfo, visc bocos are numbered
     ! before other bocos. Therefore, mm_nViscBocos == mm_nBocos
     do mm=1,nBocos
        bocoType3: if (BCType(mm) == NSWallIsoThermal) then
           do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                 ii = ii + 1
                 hflux(ii) = BCData(mm)%nodeHeatFlux(i, j)
              end do
           end do

         ! Simply put in zeros for the other wall BCs
        else if (BCType(mm) == NSWallAdiabatic .or. BCType(mm) == EulerWall) then
           do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                 ii = ii + 1
                 hflux(ii) = zero
              end do
           end do
        end if bocoType3
     end do
  end do
end subroutine getHeatFlux

subroutine getHeatFlux_d(hflux, hfluxd, npts, sps)
   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
   use surfaceFamilies, only : BCFamExchange, familyExchange, &
        zeroCellVal, zeroNodeVal, fullfamList
   use utils, only : setPointers, setPointers_d, isWallType
   implicit none
   !
   !      Local variables.
   !
   integer(kind=intType), intent(in) :: npts, sps
   real(kind=realType), intent(out) :: hflux(npts),  hfluxd(npts)
   real(kind=realType) :: localValues(nLocalValues)

   integer(kind=intType) :: mm, nn, i, j, ii
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
   type(familyExchange), pointer :: exch

   ! Set the pointer to the wall exchange:
   exch => BCFamExchange(iBCGroupWalls, sps)


   ! compute BCData(mm)%cellHeatFlux for the cells
   ! do nn=1, nDom
   !  call setPointers_d(nn, 1_intType, sps)
   !  call heatFluxes_d()
   ! end do



   ! Set the weighting factors. In this case, area
   ! the adiabatic wall cells have a weighting factor of zero
   do nn=1, nDom
      call setPointers_d(nn, 1_intType, sps)
      do mm=1, nBocos

         write(*,*) mm, 'hfluxd ', minval(BCDatad(mm)%cellHeatFlux), maxval(BCDatad(mm)%cellHeatFlux)
         ! if (BCType(mm) == NSWallIsoThermal) then
         !    BCData(mm)%cellVal => BCData(mm)%area(:, :)
         !    BCDatad(mm)%cellVal => BCDatad(mm)%area(:, :)

         ! else if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic) then
         !    BCData(mm)%cellVal => zeroCellVal
         !    BCData(mm)%nodeVal => zeroNodeVal

         !    BCDatad(mm)%cellVal => zeroCellVal
         !    BCDatad(mm)%nodeVal => zeroNodeVal
         if (isWallType(BCType(mm))) then
            BCData(mm)%cellVal => BCData(mm)%area(:, :)
            BCDatad(mm)%cellVal => BCDatad(mm)%area(:, :)

         end if
      end do
   end do

   call computeWeighting_d(exch)


   ! compute the nodal values from the cell center values and the weighting

   do nn=1, nDom
      call setPointers_d(nn, 1_intType, sps)
      do mm=1, nBocos
         write(*,*) mm, 'hfluxd ', minval(BCDatad(mm)%cellHeatFlux), maxval(BCDatad(mm)%cellHeatFlux)

         ! if (BCType(mm) == NSWallIsoThermal) then
         !    BCData(mm)%cellVal => BCData(mm)%cellHeatFlux(:, :)
         !    BCData(mm)%nodeVal => BCData(mm)%nodeHeatFlux(:, :)

         !    BCDatad(mm)%cellVal => BCDatad(mm)%cellHeatFlux(:, :)
         !    BCDatad(mm)%nodeVal => BCDatad(mm)%nodeHeatFlux(:, :)

         if (isWallType(BCType(mm))) then
            BCData(mm)%cellVal => BCData(mm)%cellHeatFlux(:, :)
            BCData(mm)%nodeVal => BCData(mm)%nodeHeatFlux(:, :)

            BCDatad(mm)%cellVal => BCDatad(mm)%cellHeatFlux(:, :)
            BCDatad(mm)%nodeVal => BCDatad(mm)%nodeHeatFlux(:, :)

         end if
      end do
   end do

   call surfaceCellCenterToNode_d(exch)

   ! Now extract into the flat array:
   ii = 0
   do nn=1,nDom
      call setPointers_d(nn,1_intType,sps)

      ! Loop over the number of viscous boundary subfaces of this block.
      ! According to preprocessing/viscSubfaceInfo, visc bocos are numbered
      ! before other bocos. Therefore, mm_nViscBocos == mm_nBocos
      do mm=1,nBocos
         write(*,*) mm, 'hfluxd ', minval(BCDatad(mm)%cellHeatFlux), maxval(BCDatad(mm)%cellHeatFlux)

         ! if (BCType(mm) == NSWallIsoThermal) then
         !    do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
         !       do i=BCData(mm)%inBeg,BCData(mm)%inEnd
         !          ii = ii + 1
         !          hflux(ii) = BCData(mm)%nodeHeatFlux(i, j)
         !          hfluxd(ii) = BCDatad(mm)%nodeHeatFlux(i, j)
         !          write(*,*) ii, hflux(ii), hfluxd(ii)
         !       end do
         !    end do

         !  ! Simply put in zeros for the other wall BCs
         ! else if (BCType(mm) == NSWallAdiabatic .or. BCType(mm) == EulerWall) then
         !    do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
         !       do i=BCData(mm)%inBeg,BCData(mm)%inEnd
         !          ii = ii + 1
         !          hflux(ii) = zero
         !          hflux(ii) = zero
         !       end do
         !    end do
         ! end if

         if (isWallType(BCType(mm))) then
            do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
               do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                  ii = ii + 1
                  hflux(ii) = BCData(mm)%nodeHeatFlux(i, j)
                  hfluxd(ii) = BCDatad(mm)%nodeHeatFlux(i, j)
                  write(*,*) ii, hflux(ii), hfluxd(ii)
               end do
            end do
         end if
      end do
   end do
   write(*,*) 'hfluxd ', minval(hfluxd), maxval(hfluxd)
   write(*,*)
end subroutine getHeatFlux_d


subroutine getHeatFlux_b(hfluxd, npts, sps)
   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad, BCFaceID
   use surfaceFamilies, only : BCFamExchange, familyExchange, &
         zeroCellVal, zeroNodeVal, fullfamList
   use utils, only : setPointers, setPointers_b, setBCPointers

   implicit none
   !
   !      Local variables.
   !
   integer(kind=intType), intent(in) :: npts, sps
   real(kind=realType), intent(out) :: hfluxd(npts)
   real(kind=realType) :: localValues(nLocalValues)

   integer(kind=intType) :: mm, nn, i, j, ii
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
   type(familyExchange), pointer :: exch
   real(kind=realType), dimension(1, npts) :: heatFluxes

   ! Set the pointer to the wall exchange:
   exch => BCFamExchange(iBCGroupWalls, sps)


   call getHeatFlux(heatFluxes, npts, sps)
   ! ! compute BCData(mm)%cellHeatFlux for the cells
   ! do nn=1, nDom
   !    call setPointers(nn, 1_intType, sps)
   !    call heatFluxes()
   ! end do


   ! ! Set the weighting factors. In this case, area
   ! ! the adiabatic wall cells have a weighting factor of zero
   ! do nn=1, nDom
   !    call setPointers(nn, 1_intType, sps)
   !    do mm=1, nBocos

   !       if (BCType(mm) == NSWallIsoThermal) then
   !          BCData(mm)%cellVal => BCData(mm)%area(:, :)
   !       else if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic) then
   !          BCData(mm)%cellVal => zeroCellVal
   !          BCData(mm)%nodeVal => zeroNodeVal
   !       end if
   !    end do
   ! end do

   ! call computeWeighting(exch)


   ! ! compute the nodal values from the cell center values and the weighting

   ! do nn=1, nDom
   !    call setPointers(nn, 1_intType, sps)
   !    do mm=1, nBocos
   !       if (BCType(mm) == NSWallIsoThermal) then
   !          BCData(mm)%cellVal => BCData(mm)%cellHeatFlux(:, :)
   !          BCData(mm)%nodeVal => BCData(mm)%nodeHeatFlux(:, :)
   !       end if
   !    end do
   ! end do

   ! call surfaceCellCenterToNode(exch)

   ! ====================
   ! Do the reverse pass:
   ! ====================

   ! Now extract into the flat array:
   ii = 0
   do nn=1,nDom
      call setPointers_b(nn,1_intType,sps)

      ! Loop over the number of viscous boundary subfaces of this block.
      ! According to preprocessing/viscSubfaceInfo, visc bocos are numbered
      ! before other bocos. Therefore, mm_nViscBocos == mm_nBocos
      do mm=1,nBocos
         if (BCType(mm) == NSWallIsoThermal) then
            do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
               do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                  ii = ii + 1
                  BCDatad(mm)%nodeHeatFlux(i, j) = hfluxd(ii)
               end do
            end do

            ! Simply put in zeros for the other wall BCs
         end if
      end do
   end do


   do nn=1, nDom
      call setPointers_b(nn, 1_intType, sps)
      do mm=1, nBocos
         if (BCType(mm) == NSWallIsoThermal) then
            BCData(mm)%cellVal => BCData(mm)%cellHeatFlux(:, :)
            BCData(mm)%nodeVal => BCData(mm)%nodeHeatFlux(:, :)

            BCDatad(mm)%cellVal => BCDatad(mm)%cellHeatFlux(:, :)
            BCDatad(mm)%nodeVal => BCDatad(mm)%nodeHeatFlux(:, :)

         end if
      end do
   end do

   call surfaceCellCenterToNode_b(exch)


   do nn=1, nDom
      call setPointers_b(nn, 1_intType, sps)
      do mm=1, nBocos

         if (BCType(mm) == NSWallIsoThermal) then
            BCData(mm)%cellVal => BCData(mm)%area(:, :)
            BCDatad(mm)%cellVal => BCDatad(mm)%area(:, :)

         else if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic) then
            BCData(mm)%cellVal => zeroCellVal
            BCData(mm)%nodeVal => zeroNodeVal

            BCDatad(mm)%cellVal => zeroCellVal
            BCDatad(mm)%nodeVal => zeroNodeVal

         end if
      end do
   end do

   call computeWeighting_b(exch)

   ! bocos: do mm=1, nBocos

   !    ! Only do this on isoThermalWalls
   !    if (BCType(mm) == NSWallIsoThermal) then

   !       call setBCPointers(mm, .True.)


   !       ! Loop over the quadrilateral faces of the subface. Note that
   !       ! the nodal range of BCData must be used and not the cell
   !       ! range, because the latter may include the halo's in i and
   !       ! j-direction. The offset +1 is there, because inBeg and jnBeg
   !       ! refer to nodal ranges and not to cell ranges.
   !       !
   !       do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
   !          do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd

   !             ! Compute the normal heat flux on the face. Inward positive.

   !             ! cellHeatFlux = q vector dotted with the Area vector

   !             write(*,*) mm, i, jBCDatad(mm)%cellHeatFlux(i,j)

   !          enddo
   !       end do
   !    end if
   ! enddo bocos


end subroutine getHeatFlux_b



! subroutine heatFluxes
!    use constants
!    use blockPointers, only : BCData, nDom, nBocos, BCType, BCFaceID, viscSubFace
!    use BCPointers, only : ssi
!    use flowVarRefState, only : pRef, rhoRef
!    use utils, only : setBCPointers
!    implicit none
!    !
!    !      Local variables.
!    !
!    integer(kind=intType) :: i, j, ii, mm
!    real(kind=realType) :: fact, scaleDim, Q
!    real(kind=realType) :: qw, qA
!    logical :: heatedSubface

!    ! Set the actual scaling factor such that ACTUAL heat flux is computed
!    ! The factor is determined from stanton number
!    scaleDim = pRef*sqrt(pRef/rhoRef)

!    ! Loop over the boundary subfaces of this block.
!    bocos: do mm=1, nBocos

!       ! Only do this on isoThermalWalls
!       if (BCType(mm) == NSWallIsoThermal) then

!          call setBCPointers(mm, .True.)

!          select case (BCFaceID(mm))
!          case (iMin, jMin, kMin)
!             fact = -one
!          case (iMax, jMax, kMax)
!             fact = one
!          end select

!          ! Loop over the quadrilateral faces of the subface. Note that
!          ! the nodal range of BCData must be used and not the cell
!          ! range, because the latter may include the halo's in i and
!          ! j-direction. The offset +1 is there, because inBeg and jnBeg
!          ! refer to nodal ranges and not to cell ranges.
!          !
!          do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
!             do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd

!                ! Compute the normal heat flux on the face. Inward positive.

!                ! cellHeatFlux = q vector dotted with the Area vector

!                BCData(mm)%cellHeatFlux(i,j) = -fact*scaleDim* &
!                     ( viscSubface(mm)%q(i,j,1)*ssi(i,j,1) &
!                     + viscSubface(mm)%q(i,j,2)*ssi(i,j,2) &
!                     + viscSubface(mm)%q(i,j,3)*ssi(i,j,3))

!             enddo
!          end do
!       end if
!    enddo bocos

!  end subroutine heatFluxes


!  subroutine heatFluxes_d
!    ! forward mode differenitaion of heatfluexes
!    ! the inputs are ssid, viscSubFaced%q, and BCDatad(mm)%norm
!    ! the output is BCDatad(mm)%cellHeatFlux

!    use constants
!    use blockPointers, only : BCData, BCDatad, nDom, nBocos, BCType, BCFaceID, viscSubFace, viscSubFaced
!    use BCPointers, only : ssi, ssid
!    use flowVarRefState, only : pRef, rhoRef
!    use utils, only : setBCPointers_d
!    implicit none
!    !
!    !      Local variables.
!    !
!    integer(kind=intType) :: i, j, ii, mm
!    real(kind=realType) :: fact, scaleDim
!    real(kind=realType) :: qw, qA
!    logical :: heatedSubface

!    ! Set the actual scaling factor such that ACTUAL heat flux is computed
!    ! The factor is determined from stanton number
!    scaleDim = pRef*sqrt(pRef/rhoRef)

!    ! Loop over the boundary subfaces of this block.
!    bocos: do mm=1, nBocos

!       ! Only do this on isoThermalWalls
!       if (BCType(mm) == NSWallIsoThermal) then

!          call setBCPointers_d(mm, .True.)

!          select case (BCFaceID(mm))
!          case (iMin, jMin, kMin)
!             fact = -one
!          case (iMax, jMax, kMax)
!             fact = one
!          end select

!          ! Loop over the quadrilateral faces of the subface. Note that
!          ! the nodal range of BCData must be used and not the cell
!          ! range, because the latter may include the halo's in i and
!          ! j-direction. The offset +1 is there, because inBeg and jnBeg
!          ! refer to nodal ranges and not to cell ranges.
!          !
!          do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
!             do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd

!                ! Compute the normal heat flux on the face. Inward positive.
!                BCDatad(mm)%cellHeatFlux(i,j) = -fact*scaleDim* &
!                     ( viscSubfaced(mm)%q(i,j,1)*ssi(i,j,1) &
!                     + viscSubface(mm)%q(i,j,1)*ssid(i,j,1) &
!                     + viscSubfaced(mm)%q(i,j,2)*ssi(i,j,2) &
!                     + viscSubface(mm)%q(i,j,2)*ssid(i,j,2) &
!                     + viscSubfaced(mm)%q(i,j,3)*ssi(i,j,3) &
!                     + viscSubface(mm)%q(i,j,3)*ssid(i,j,3))
!                write(*,*) 'BCDatad(mm)%cellHeatFlux', i, j, BCDatad(mm)%cellHeatFlux(i,j)
!                write(*,*) 'viscSubFace', viscSubfaced(mm)%q(i,j,:)
!             enddo
!          end do
!       end if
!    enddo bocos

!  end subroutine heatFluxes_d


!  subroutine heatFluxes_b
!    ! forward mode differenitaion of heatfluexes
!    ! the inputs are ssid, viscSubFaced%q, and BCDatad(mm)%norm
!    ! the output is BCDatad(mm)%cellHeatFlux

!    use constants
!    use blockPointers, only : BCData, BCDatad, nDom, nBocos, BCType, BCFaceID, viscSubFace, viscSubFaced
!    use BCPointers, only : ssi, ssid
!    use flowVarRefState, only : pRef, rhoRef
!    use utils, only :  setBCPointers_d
!    implicit none
!    !
!    !      Local variables.
!    !
!    integer(kind=intType) :: i, j, ii, mm
!    real(kind=realType) :: fact, scaleDim
!    real(kind=realType) :: qw, qA
!    logical :: heatedSubface

!    ! Set the actual scaling factor such that ACTUAL heat flux is computed
!    ! The factor is determined from stanton number
!    scaleDim = pRef*sqrt(pRef/rhoRef)

!    ! Loop over the boundary subfaces of this block.
!    bocos: do mm=1, nBocos

!       ! Only do this on isoThermalWalls
!       if (BCType(mm) == NSWallIsoThermal) then

!          call setBCPointers_d(mm, .True.)

!          select case (BCFaceID(mm))
!          case (iMin, jMin, kMin)
!             fact = -one
!          case (iMax, jMax, kMax)
!             fact = one
!          end select

!          ! Loop over the quadrilateral faces of the subface. Note that
!          ! the nodal range of BCData must be used and not the cell
!          ! range, because the latter may include the halo's in i and
!          ! j-direction. The offset +1 is there, because inBeg and jnBeg
!          ! refer to nodal ranges and not to cell ranges.
!          !
!          do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
!             do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd

!                ! Compute the normal heat flux on the face. Inward positive.
!                ssid(i,j,1) = ssid(i,j,1) + viscSubface(mm)%q(i,j,1)*BCDatad(mm)%cellHeatFlux(i,j)
!                viscSubfaced(mm)%q(i,j,1) = viscSubfaced(mm)%q(i,j,1) + ssi(i,j,1)*BCDatad(mm)%cellHeatFlux(i,j)

!                ssid(i,j,2) = ssid(i,j,2) + viscSubface(mm)%q(i,j,2)*BCDatad(mm)%cellHeatFlux(i,j)
!                viscSubfaced(mm)%q(i,j,2) = viscSubfaced(mm)%q(i,j,2) + ssi(i,j,2)*BCDatad(mm)%cellHeatFlux(i,j)

!                ssid(i,j,3) = ssid(i,j,3) + viscSubface(mm)%q(i,j,3)*BCDatad(mm)%cellHeatFlux(i,j)
!                viscSubfaced(mm)%q(i,j,3) = viscSubfaced(mm)%q(i,j,3) + ssi(i,j,3)*BCDatad(mm)%cellHeatFlux(i,j)

!             enddo
!          end do
!       end if
!    enddo bocos

!  end subroutine heatFluxes_b


subroutine getTNSWall(tnsw, npts, sps)

  use constants
  use blockPointers, only : nDom, nBocos, BCData, BCType
  use flowVarRefState, only : TRef
  use utils, only : setPointers
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(out) :: tnsw(npts)

  ! Local Variables
  integer(kind=intType) :: mm, nn, i, j, ii
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  ii = 0
  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)
     ! Loop over the number of viscous boundary subfaces of this block.
     bocos: do mm=1,nBocos
        isoWall: if (BCType(mm) == NSWallIsoThermal) then
           jBeg = BCdata(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd
           do j=jBeg,jEnd
              do i=iBeg, iEnd
                 ii = ii + 1
                 tnsw(ii) = BCData(mm)%TNS_Wall(i,j)*Tref
              end do
           end do
        end if isoWall
     end do bocos

   end do domains

end subroutine getTNSWall


subroutine getWallTemperature(wallTempNodes, npts, sps)

   ! returns the wall nodal wall temperature for the walls in the mesh




   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData
   use BCPointers, only : ww1, ww2, pp1, pp2
   use flowVarRefState, only : TRef, RGas
   use surfaceFamilies, only : BCFamExchange, familyExchange, &
   zeroCellVal, zeroNodeVal
   use utils, only : setPointers, setBCPointers, isWallType
   implicit none



   ! Input Variables
   integer(kind=intType), intent(in) :: npts, sps
   real(kind=realType), intent(out) :: wallTempNodes(npts)

   ! Local Variables
   integer(kind=intType) :: mm, nn, i, j, ii, jj
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
   real(kind=realType) :: t1, t2, wallTempCells
   type(familyExchange), pointer :: exch

   exch => BCFamExchange(iBCGroupWalls, sps)

   ! compute the wall tempertaure for each BC that is a wall
   domains: do nn=1,nDom
      call setPointers(nn, 1_intType, sps)
      ! Loop over the number of viscous boundary subfaces of this block.
      bocos: do mm=1,nBocos
         wall: if(isWallType(BCType(mm))) then
            call setBCPointers(mm, .True.)

            do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
               do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd

                  ! The wall temperature is the average of the first off wall cell
                  ! and the first halo cell
                  t2 = pp2(i,j)/(RGas*ww2(i,j,irho))
                  t1 = pp1(i,j)/(RGas*ww1(i,j,irho))

                  BCData(mm)%cellTemperature(i,j) = (t2 + t1)/two*Tref *BCData(mm)%area(i,j)

               enddo
            end do
         end if wall
      end do bocos

   end do domains

   do nn=1, nDom
      call setPointers(nn, 1_intType, sps)
      do mm=1, nBocos
         iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
         jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

         bocoType1: if(isWallType(BCType(mm))) then
            BCData(mm)%cellVal => BCData(mm)%area(:, :)
         end if bocoType1
      end do
   end do
   call computeWeighting(exch)

   do nn=1, nDom
      call setPointers(nn, 1_intType, sps)
      do mm=1, nBocos
         wall2: if(isWallType(BCType(mm))) then
            BCData(mm)%cellVal => BCData(mm)%cellTemperature(:, :)
            BCData(mm)%nodeVal => BCData(mm)%nodeTemperature(:, :)
         end if wall2
      end do
   end do

 call surfaceCellCenterToNode(exch)

 ! Now extract into the flat array:
 ii = 0
 do nn=1,nDom
    call setPointers(nn,1_intType,sps)

    ! Loop over the number of viscous boundary subfaces of this block.
    ! According to preprocessing/viscSubfaceInfo, visc bocos are numbered
    ! before other bocos. Therefore, mm_nViscBocos == mm_nBocos
    do mm=1,nBocos
      wall3: if(isWallType(BCType(mm))) then
         do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
             do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                ii = ii + 1
                wallTempNodes(ii) = BCData(mm)%nodeTemperature(i, j)
             end do
         end do
      end if wall3
    end do
 end do
end subroutine getWallTemperature
