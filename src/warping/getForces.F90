
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
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
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

     call vecRestoreArrayF90(zipper%localVal, localPtr, ierr)
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
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
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
  use utils, only : EChk, setPointers, setPointers_d
  use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
  use surfaceFamilies, only : familyExchange, BCFamExchange
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
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

        call vecRestoreArrayF90(zipper%localVal, localPtr, ierr)
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
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none

  type(familyExchange) :: exch
  integer(kind=intType) ::  sps
  integer(kind=intType) :: mm, nn, i, j, ii, jj, iDim, ierr
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  real(kind=realType) :: qv
  real(kind=realType), dimension(:), pointer :: nodeValLocPtr

  ! We assume that normalization factor is already computed
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
                 BCData(mm)%nodeVal(i, j) = nodeValLocPtr(ii)
              end do
           end do
        end if famInclude2
     end do
  end do

  call vecRestoreArrayF90(exch%nodeValLocal, nodeValLocPtr, ierr)
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

   call vecGetArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
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


   call vecRestoreArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
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
   call vecSet(exch%nodeValGlobald, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterBegin(exch%scatter, exch%nodeValLocald, &
        exch%nodeValGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValLocald, &
        exch%nodeValGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

    ! The product rule here: (since we are multiplying)
    ! nodeValGlobal = nodeValGlobal * invArea
    ! nodeValGlobald = nodeValGlobald*invArea + nodeValGlobal*invAread

     ! First term:  nodeValGlobald = nodeValGlobald*invArea
   call vecPointwiseMult(exch%nodeValGlobald, exch%nodeValGlobald, &
   exch%sumGlobal, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! Second term:, tmp = nodeValGlobal*invAread
   call vecPointwiseMult(tmp, exch%nodeValGlobal, exch%sumGlobald, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Sum the second term into the first
   call VecAXPY(exch%nodeValGlobald, one, tmp, ierr)
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
   call VecScatterBegin(exch%scatter, exch%nodeValGlobald, &
      exch%nodeValLocald, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValGlobald, &
      exch%nodeValLocald, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
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

   call vecRestoreArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)


  call VecDestroy(tmp, ierr)
  call EChk(ierr,__FILE__,__LINE__)

 end subroutine surfaceCellCenterToNode_d

 subroutine surfaceCellCenterToNode_b(exch)

  ! This routine performs the reverse of surfaceCellCenterToNode_d. It
  ! takes in BCDatad(mm)%nodeVal and perfroms the reverse of the
  ! nodal averaging procedure in getForces to compute BCDatad(mm)%cellVal
  ! and exch%sumGlobald.

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
   real(kind=realType), dimension(:), pointer :: nodeValLocPtr, nodeValLocPtrd
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

   call vecGetArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
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
                  nodeValLocPtrd(ii) = BCDatad(mm)%nodeVal(i, j)

                  ! zero the input ad seed
                  BCDatad(mm)%nodeVal(i, j) = zero
               end do
            end do
         end if famInclude2
      end do
   end do

   call vecRestoreArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! --- Globalize the local values ---
   call vecSet(exch%nodeValGlobald, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterBegin(exch%scatter, exch%nodeValLocald, &
        exch%nodeValGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValLocald, &
        exch%nodeValGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! ---  account for the weighting process ---
   ! in the forward pass
   ! exch%nodeValGlobal = exch%nodeValGlobal *  exch%sumGlobal
   ! becomes...
   ! exch%nodeValGlobald += exch%nodeValGlobald * exch%sumGlobal
   ! exch%sumGlobald += exch%nodeValGlobal* exch%nodeValGlobald
         ! (note non weighted value used ^)
   ! in the reverse pass

   call vecPointwiseMult(tmp, exch%nodeValGlobal, exch%nodeValGlobald, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   call vecPointwiseMult(exch%nodeValGlobald,  exch%sumGlobal, exch%nodeValGlobald, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Accumulate
   call vecAXPY(exch%sumGlobald, one, tmp, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! --- Now communicate exch%nodeValGlobald back to the local patches ---

   call VecScatterBegin(exch%scatter, exch%nodeValGlobald, &
   exch%nodeValLocald, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValGlobald, &
   exch%nodeValLocald, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Copy the values into patches
   call vecGetArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
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
                     qv_b = qv_b + nodeValLocPtrd(ind(jj))
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

   call vecRestoreArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
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
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
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

   call vecGetArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
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

   call vecRestoreArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
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
   call vecSet(exch%sumGlobald, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterBegin(exch%scatter, exch%nodeValLocald, &
        exch%sumGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%nodeValLocald, &
        exch%sumGlobald, ADD_VALUES, SCATTER_FORWARD, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! Now compute the inverse of the weighting so that we can multiply
   ! instead of dividing. Here we need the original value too:

   call vecGetArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%sumGlobald, sumGlobalPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   sumGlobalPtrd = -(sumGlobalPtrd/sumGlobalPtr**2)
   sumGlobalPtr = one/sumGlobalPtr

   call vecGetArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecRestoreArrayF90(exch%sumGlobald, sumGlobalPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

end subroutine computeWeighting_d


subroutine computeWeighting_b(exch)
   ! in: exch%sumGlobald
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
   real(kind=realType) :: qa, qad
   real(kind=realType), dimension(:), pointer :: nodeValLocPtr,   sumGlobalPtr,&
                                                 nodeValLocPtrd, sumGlobalPtrd

   sps = exch%sps

   call computeWeighting(exch)

   call vecGetArrayF90(exch%sumGlobald, sumGlobalPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! Keep in mind sumGlobalPtr points to sumGlobal which already has
   ! been inversed so we just multiply.
   sumGlobalPtrd = -sumGlobalPtrd*sumGlobalPtr**2

   call vecRestoreArrayF90(exch%sumGlobald, sumGlobalPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecRestoreArrayF90(exch%sumGlobal, sumGlobalPtr, ierr)
   call EChk(ierr,__FILE__,__LINE__)



   ! Push back to the local patches
   call VecScatterBegin(exch%scatter, exch%sumGlobald, &
         exch%nodeValLocald, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call VecScatterEnd(exch%scatter, exch%sumGlobald, &
   exch%nodeValLocald, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   call vecGetArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
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
                  qad = zero
                  do jj=1,4
                     qad = qad + nodeValLocPtrd(ind(jj))
                  end do
                  BCDatad(mm)%cellVal(i+1, j+1) = &
                                  BCDatad(mm)%cellVal(i+1, j+1) + fourth*qad
               end do
            end do
            ii = ii + ni*nj
         end if famInclude2
      end do
   end do


   call vecRestoreArrayF90(exch%nodeValLocald, nodeValLocPtrd, ierr)
   call EChk(ierr,__FILE__,__LINE__)


   ! zero the iput AD Seed

   call vecSet(exch%sumGlobald, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)


end subroutine computeWeighting_b


subroutine computeNodalTractions(sps)
  use constants
  use blockPointers, only : BCData, nDom, nBocos, BCType
  use surfaceFamilies, only : BCFamExchange, familyExchange
  use utils, only : setPointers, EChk, isWallType
  implicit none

  integer(kind=intType), intent(in) ::  sps
  integer(kind=intType) :: mm, nn, iDim, ierr
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
  implicit none
  integer(kind=intType), intent(in) :: sps
  integer(kind=intType) :: mm, nn, iDim, ierr
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
#include <petsc/finclude/petsc.h>
   use petsc
   implicit none

   integer(kind=intType), intent(in) :: sps
   integer(kind=intType) :: mm, nn, iDim, ierr
   type(familyExchange), pointer :: exch

   ! For better readibility
   exch => BCFamExchange(iBCGroupWalls, sps)

   ! intialize weighting seed
   call vecSet(exch%sumGlobald, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! ==================================
   !  Recompute the dual area
   ! ==================================


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
  use utils, only : setPointers_d
  implicit none

  integer(kind=intType), intent(in) ::  sps

  integer(kind=intType) :: mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: qfd(3)

  do nn=1, nDom
     call setPointers_d(nn, 1_intType, sps)
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


  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
  use utils, only : setPointers_b
  implicit none

  integer(kind=intType), intent(in) ::  sps

  integer(kind=intType) :: mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: qf_b(3)

  domains: do nn=1,nDom
     call setPointers_b(nn, 1_intType, sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

            BCDatad(mm)%Fv = zero
            BCDatad(mm)%Fp = zero


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




subroutine getHeatXferRate(hxfer, npts, sps)
   ! computes and returns the nodal heat xfer rate from the cell center values.
   ! If heatxferratesAsfluxes, then the nodal heat fluxes are returned instead.
   ! This routine is modeled after the get forces routine above.

   ! THIS FUNCTION DOES NOT WORK WITH ZIPPER MESHES

   ! inputs
   !    BCDdata%area
   !             cell areas
   !    BCdata%cellHeatXferRate
   !             heat transfer rate defined at the cell centers

   ! outputs
   !    hxfer
   !       the heat transfer rate or the heat flux depending on `heatxferratesAsfluxes`


   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData
   use surfaceFamilies, only : BCFamExchange, familyExchange, &
        zeroCellVal, zeroNodeVal, fullfamList
   use surfaceIntegrations, only : integrateSurfaces
   use inputPhysics, only : heatxferratesAsfluxes

   use utils, only : setPointers, isWallType
   implicit none
   !
   !      Local variables.
   !
   integer(kind=intType), intent(in) :: npts, sps
   real(kind=realType), intent(out) :: hxfer(npts)
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


   if (heatxferratesAsfluxes) then
      call computeNodalHeatFlux(sps)
   else
      call computeNodalHeatXferRate(sps)
   end if

   ! Now extract into the flat array:
   ii = 0
   do nn=1,nDom
      call setPointers(nn,1_intType,sps)

      do mm=1,nBocos
         bocoType3: if(isWallType(BCType(mm))) then

             do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
                do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                  ii = ii + 1
                  if (heatxferratesAsfluxes) then
                     hxfer(ii) = BCData(mm)%nodeHeatFlux(i, j)
                  else
                     hxfer(ii) = BCData(mm)%nodeHeatXferRate(i, j)
                  end if
                end do
            end do
         end if bocoType3
      end do
   end do
end subroutine getHeatXferRate


subroutine getHeatXferRate_d(hxfer, hxferd, npts, sps)
   ! computes and returns the nodal heat xfer rate from the cell center values and the derivatives.
   ! If heatxferratesAsfluxes, then the nodal heat fluxes are returned instead.
   ! This routine is modeled after the get forces routine above.

   ! THIS FUNCTION DOES NOT WORK WITH ZIPPER MESHES

   ! inputs
   !    BCDdata%area
   !             cell areas
   !    BCdata%cellHeatXferRate
   !             heat transfer rate defined at the cell centers
   !    BCDdatad%area
   !             cell areas
   !    BCdatad%cellHeatXferRate
   !             heat transfer rate defined at the cell centers

   ! outputs
   !    hxfer
   !       the heat transfer rate or the heat flux depending on `heatxferratesAsfluxes`
   !    hxferd
   !       the heat transfer rate or the heat flux depending on `heatxferratesAsfluxes`


   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
   use surfaceFamilies, only : BCFamExchange, familyExchange, &
        zeroCellVal, zeroNodeVal, fullfamList
   use surfaceIntegrations, only : integrateSurfaces
   use inputPhysics, only : heatxferratesAsfluxes

   use utils, only : setPointers_d, isWallType
   implicit none
   !
   !      Local variables.
   !
   integer(kind=intType), intent(in) :: npts, sps
   real(kind=realType), intent(out) :: hxfer(npts)
   real(kind=realType), intent(out) :: hxferd(npts)


   integer(kind=intType) :: mm, nn, i, j, ii
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
   type(familyExchange), pointer :: exch


   ! Set the pointer to the wall exchange:
   exch => BCFamExchange(iBCGroupWalls, sps)

   ! this routine isn't part of the api and is always called after inegrateSurfaces_d


   if (heatxferratesAsfluxes) then
      call computeNodalHeatFlux_d(sps)
   else
      call computeNodalHeatXferRate_d(sps)
   end if

   ! Now extract into the flat array:
   ii = 0
   do nn=1,nDom
      call setPointers_d(nn,1_intType,sps)

      do mm=1,nBocos
         bocoType3: if(isWallType(BCType(mm))) then

             do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
                do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                  ii = ii + 1
                  if (heatxferratesAsfluxes) then
                     hxfer(ii) = BCData(mm)%nodeHeatFlux(i, j)
                     hxferd(ii) = BCDatad(mm)%nodeHeatFlux(i, j)
                  else
                     hxfer(ii) = BCData(mm)%nodeHeatXferRate(i, j)
                     hxferd(ii) = BCDatad(mm)%nodeHeatXferRate(i, j)
                  end if
                end do
            end do
         end if bocoType3
      end do
   end do
end subroutine getHeatXferRate_d


subroutine getHeatXferRate_b(hxferd, npts, sps)
   ! bwd mode of get getheatxferrate

   ! THIS FUNCTION DOES NOT WORK WITH ZIPPER MESHES

   ! inputs

   !    hxferd
   !       the heat transfer rate or the heat flux depending on `heatxferratesAsfluxes`

   ! outputs
   !    BCdatad%area
   !             heat transfer rate defined at the cell centers

   !    BCdatad%cellHeatXferRate
   !             heat transfer rate defined at the cell centers

   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad, BCFaceID
   use surfaceFamilies, only : BCFamExchange, familyExchange, &
         zeroCellVal, zeroNodeVal, fullfamList
   use utils, only : EChk, setPointers, setPointers_b, setBCPointers, isWallType
   use inputPhysics, only : heatxferratesAsfluxes

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
   real(kind=realType), intent(in) :: hxferd(npts)

   !
   !      Local variables.
   !

   integer(kind=intType) :: mm, nn, i, j, ii,  ierr
   type(familyExchange), pointer :: exch
   real(kind=realType) :: hflux(npts)

   ! Set the pointer to the wall exchange:
   exch => BCFamExchange(iBCGroupWalls, sps)

   ! Initalize the weighting seed
   call vecSet(exch%sumGlobald, zero, ierr)
   call EChk(ierr,__FILE__,__LINE__)

   ! ! ====================
   ! ! forward pass:
   ! ! ====================
   call getHeatXferRate(hflux, npts, sps)

      ! ====================
   ! Do the reverse pass:
   ! ====================

   ! Now extract the seeds:
   ii = 0
   do nn=1,nDom
      call setPointers_b(nn,1_intType,sps)

      do mm=1,nBocos
        bocoType3: if(isWallType(BCType(mm))) then
            do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
               do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                  ii = ii + 1
                  if (heatxferratesAsfluxes) then
                     BCDatad(mm)%nodeHeatFlux(i, j) = hxferd(ii)
                  else
                     BCDatad(mm)%nodeHeatXferRate(i, j) = hxferd(ii)
                  end if

               end do
            end do
        end if bocoType3
      end do
   end do


   if (heatxferratesAsfluxes) then
      call computeNodalHeatFlux_b(sps)
   else
      call computeNodalHeatXferRate_b(sps)
   end if

   ! this routine isn't part of the api and is always called *before* inegrateSurfaces_d


end subroutine getHeatXferRate_b




subroutine computeNodalHeatFlux(sps)
   ! inputs
   !    BCDdata%area
   !             cell areas
   !    BCdata%cellHeatXferRate
   !             heat transfer rate defined at the cell centers

   ! outputs
   !    BCdata%nodalHeatFlux
   !           The heat flux defined at the cell nodes


  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData
  use surfaceFamilies, only : BCFamExchange, familyExchange, &
       zeroCellVal, zeroNodeVal, fullfamList
  use surfaceIntegrations, only : integrateSurfaces

  use utils, only : setPointers, isWallType
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) ::  sps
  integer(kind=intType) :: mm, nn
  type(familyExchange), pointer :: exch

  ! Set the pointer to the wall exchange:
  exch => BCFamExchange(iBCGroupWalls, sps)

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

  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     do mm=1, nBocos
         bocoType2: if(isWallType(BCType(mm))) then

           BCData(mm)%cellVal => BCData(mm)%cellHeatXferRate(:, :)
           BCData(mm)%nodeVal => BCData(mm)%nodeHeatFlux(:, :)
        end if bocoType2

     end do
  end do

    ! compute the nodal values from the cell center values and the weighting
   call surfaceCellCenterToNode(exch)

end subroutine computeNodalHeatFlux


subroutine computeNodalHeatFlux_d(sps)
   ! forward mode differenitaion of computeNodalHeatFlux
   ! the inputs are BCDdatad%area and BCdatad%cellHeatXferRate
   ! the output is BCDdatad%nodeHeatflux

   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
   use surfaceFamilies, only : BCFamExchange, familyExchange, &
        zeroCellVal, zeroNodeVal, fullfamList
   use utils, only : setPointers, setPointers_d, isWallType
   implicit none

   integer(kind=intType), intent(in) ::  sps

   !
   !      Local variables.
   !
   integer(kind=intType) :: mm, nn, i, j, ii
   type(familyExchange), pointer :: exch

   ! Set the pointer to the wall exchange:
   exch => BCFamExchange(iBCGroupWalls, sps)


   ! Set the weighting factors. In this case, area
   ! the adiabatic wall cells have a weighting factor of zero
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


   ! compute the nodal values from the cell center values and the weighting

   do nn=1, nDom
      call setPointers_d(nn, 1_intType, sps)
      do mm=1, nBocos



         bocoType2: if(isWallType(BCType(mm))) then
            BCData(mm)%cellVal => BCData(mm)%cellHeatXferRate(:, :)
            BCData(mm)%nodeVal => BCData(mm)%nodeHeatFlux(:, :)

            BCDatad(mm)%cellVal => BCDatad(mm)%cellHeatXferRate(:, :)
            BCDatad(mm)%nodeVal => BCDatad(mm)%nodeHeatFlux(:, :)

        end if bocoType2
      end do
   end do

   call surfaceCellCenterToNode_d(exch)


end subroutine computeNodalHeatFlux_d


subroutine computeNodalHeatFlux_b(sps)
   ! in: bcdata%nodeheatflux
   ! out: bcdatad%area, bcdatad%cellHeatXferRate
   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad, BCFaceID
   use surfaceFamilies, only : BCFamExchange, familyExchange, &
         zeroCellVal, zeroNodeVal, fullfamList
   use utils, only : EChk, setPointers, setPointers_b, setBCPointers, isWallType

   integer(kind=intType), intent(in) ::  sps

   !
   !      Local variables.
   !

   integer(kind=intType) :: mm, nn, i, j, ii,  ierr
   type(familyExchange), pointer :: exch

   ! Set the pointer to the wall exchange:
   exch => BCFamExchange(iBCGroupWalls, sps)
   ! propagate the seeds from the nodes to the cell center and weighting factor
   do nn=1, nDom
      call setPointers_b(nn, 1_intType, sps)
      do mm=1, nBocos

        bocoType2: if(isWallType(BCType(mm))) then
            BCData(mm)%cellVal => BCData(mm)%cellHeatXferRate(:, :)
            BCData(mm)%nodeVal => BCData(mm)%nodeHeatFlux(:, :)

            BCDatad(mm)%cellVal => BCDatad(mm)%cellHeatXferRate(:, :)
            BCDatad(mm)%nodeVal => BCDatad(mm)%nodeHeatFlux(:, :)

        end if bocoType2


      end do
   end do

   call surfaceCellCenterToNode_b(exch)


   do nn=1, nDom
      call setPointers_b(nn, 1_intType, sps)
      do mm=1, nBocos

        bocoType1: if(isWallType(BCType(mm))) then
            BCData(mm)%cellVal => BCData(mm)%area(:, :)
            BCDatad(mm)%cellVal => BCDatad(mm)%area(:, :)
        end if bocoType1
      end do
   end do

   call computeWeighting_b(exch)

end subroutine computeNodalHeatFlux_b


subroutine computeNodalHeatXferRate(sps)
   ! computes the nodal heat transfer rate  (units of energy/time e.g. Watts)
   ! Based on the doc string for computeNodal forces, no scatter and gather is need

   ! inputs
   !    BCdata%cellHeatXferRate
   !             heat transfer rate defined at the cell centers

   ! outputs
   !    BCdata%nodeHeatXferRate
   !           The heat transfer rate defined at the cell nodes


   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData
   use utils, only : setPointers
   implicit none

   integer(kind=intType), intent(in) ::  sps

   integer(kind=intType) :: mm, nn, i, j
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
   real(kind=realType) :: qf

   do nn=1, nDom
      call setPointers(nn, 1_intType, sps)
      do mm=1, nBocos
         iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
         jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd

         if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
              BCType(mm) == NSWallIsothermal) then
            BCData(mm)%nodeHeatXferRate = zero
            do j=jBeg, jEnd
               do i=iBeg, iEnd
                  qf = fourth*BCData(mm)%cellHeatXferRate(i,j)
                  BCData(mm)%nodeHeatXferRate(i  , j  ) = BCData(mm)%nodeHeatXferRate(i  , j) + qf
                  BCData(mm)%nodeHeatXferRate(i-1, j  ) = BCData(mm)%nodeHeatXferRate(i-1, j) + qf
                  BCData(mm)%nodeHeatXferRate(i  , j-1) = BCData(mm)%nodeHeatXferRate(i  , j-1) + qf
                  BCData(mm)%nodeHeatXferRate(i-1, j-1) = BCData(mm)%nodeHeatXferRate(i-1, j-1) + qf
               end do
            end do
         end if
      end do
   end do
end subroutine computeNodalHeatXferRate

subroutine computeNodalHeatXferRate_d(sps)
   ! computes the nodal heat transfer rate  (units of energy/time e.g. Watts)
   ! Based on the doc string for computeNodal forces, no scatter and gather is need

   ! inputs
   !    BCdata%cellHeatXferRate
   !             heat transfer rate defined at the cell centers

   ! outputs
   !    BCdata%nodeHeatXferRate
   !           The heat transfer rate defined at the cell nodes


   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
   use utils, only : setPointers_d
   implicit none

   integer(kind=intType), intent(in) ::  sps

   integer(kind=intType) :: mm, nn, i, j
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
   real(kind=realType) :: qf, qfd

   do nn=1, nDom
      call setPointers_d(nn, 1_intType, sps)
      do mm=1, nBocos
         iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
         jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd

         if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
              BCType(mm) == NSWallIsothermal) then
            BCData(mm)%nodeHeatXferRate = zero
            BCDatad(mm)%nodeHeatXferRate = zero
            do j=jBeg, jEnd
               do i=iBeg, iEnd
                  qf = fourth*BCData(mm)%cellHeatXferRate(i,j)
                  qfd = fourth*BCDatad(mm)%cellHeatXferRate(i,j)

                  BCData(mm)%nodeHeatXferRate(i  , j  ) = BCData(mm)%nodeHeatXferRate(i  , j) + qf
                  BCData(mm)%nodeHeatXferRate(i-1, j  ) = BCData(mm)%nodeHeatXferRate(i-1, j) + qf
                  BCData(mm)%nodeHeatXferRate(i  , j-1) = BCData(mm)%nodeHeatXferRate(i  , j-1) + qf
                  BCData(mm)%nodeHeatXferRate(i-1, j-1) = BCData(mm)%nodeHeatXferRate(i-1, j-1) + qf

                  BCDatad(mm)%nodeHeatXferRate(i  , j  ) = BCDatad(mm)%nodeHeatXferRate(i  , j) + qfd
                  BCDatad(mm)%nodeHeatXferRate(i-1, j  ) = BCDatad(mm)%nodeHeatXferRate(i-1, j) + qfd
                  BCDatad(mm)%nodeHeatXferRate(i  , j-1) = BCDatad(mm)%nodeHeatXferRate(i  , j-1) + qfd
                  BCDatad(mm)%nodeHeatXferRate(i-1, j-1) = BCDatad(mm)%nodeHeatXferRate(i-1, j-1) + qfd
               end do
            end do
         end if
      end do
   end do
end subroutine computeNodalHeatXferRate_d

subroutine computeNodalHeatXferRate_b(sps)
   ! computes the nodal heat transfer rate  (units of energy/time e.g. Watts)
   ! Based on the doc string for computeNodal forces, no scatter and gather is need

   ! inputs
   !    BCdatad%nodeHeatXferRate
   !           The heat transfer rate defined at the cell nodes

   ! outputs
   !    BCdatad%cellHeatXferRate
   !             heat transfer rate defined at the cell centers



   use constants
   use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
   use utils, only : setPointers_b
   implicit none

   integer(kind=intType), intent(in) ::  sps

   integer(kind=intType) :: mm, nn, i, j
   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
   real(kind=realType) :: qf, qfd

   do nn=1, nDom
      call setPointers_b(nn, 1_intType, sps)
      do mm=1, nBocos
         iBeg = BCdata(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
         jBeg = BCdata(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd

         if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
              BCType(mm) == NSWallIsothermal) then

            BCDatad(mm)%cellHeatXferRate = zero

            do j=jBeg, jEnd
               do i=iBeg, iEnd
                  qfd = fourth*(BCDatad(mm)%nodeHeatXferRate(i  , j  ) + &
                                BCDatad(mm)%nodeHeatXferRate(i-1, j  ) + &
                                BCDatad(mm)%nodeHeatXferRate(i  , j-1) + &
                                BCDatad(mm)%nodeHeatXferRate(i-1, j-1))

                  BCDatad(mm)%cellHeatXferRate(i,j) = BCDatad(mm)%cellHeatXferRate(i,j) + qfd
               end do
            end do
         end if
      end do
   end do
end subroutine computeNodalHeatXferRate_b


! subroutine getHeatFluxCellCenter(hflux, npts, sps)
!    !TODO implement this subroutine

!   use constants
!   use blockPointers, only : nDom, nBocos, BCType, BCData
!   use surfaceFamilies, only : BCFamExchange, familyExchange, &
!        zeroCellVal, zeroNodeVal, fullfamList
!   use surfaceIntegrations, only : integrateSurfaces

!   use utils, only : setPointers
!   implicit none
!   !
!   !      Local variables.
!   !
!   integer(kind=intType), intent(in) :: npts, sps
!   real(kind=realType), intent(out) :: hflux(npts)
!   real(kind=realType) :: localValues(nLocalValues)

!   integer(kind=intType) :: mm, nn, i, j, ii
!   integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
!   type(familyExchange), pointer :: exch

!   ! Set the pointer to the wall exchange:
!   exch => BCFamExchange(iBCGroupWalls, sps)

!   domains: do nn=1,nDom
!       call setPointers(nn, 1_intType, sps)
!       localValues = zero
!       call integrateSurfaces(localValues, fullFamList)
!    end do domains


!   ! Now extract into the flat array:
!   ii = 0
!   do nn=1,nDom
!      call setPointers(nn,1_intType,sps)

!      ! Loop over the number of viscous boundary subfaces of this block.
!      ! According to preprocessing/viscSubfaceInfo, visc bocos are numbered
!      ! before other bocos. Therefore, mm_nViscBocos == mm_nBocos
!      do mm=1,nBocos
!         bocoType3: if (BCType(mm) == NSWallIsoThermal) then
!            do j=(BCData(mm)%jcBeg+1),(BCData(mm)%jcEnd-1)
!               do i=(BCData(mm)%icBeg+1),(BCData(mm)%icEnd-1)
!                  ii = ii + 1
!                  hflux(ii) = BCData(mm)%cellHeatXferRate(i, j)
!               end do
!            end do

!          ! Simply put in zeros for the other wall BCs
!         else if (BCType(mm) == NSWallAdiabatic .or. BCType(mm) == EulerWall) then
!            do j=(BCData(mm)%jcBeg+1),(BCData(mm)%jcEnd-1)
!               do i=(BCData(mm)%icBeg+1),(BCData(mm)%icEnd-1)
!                  ii = ii + 1
!                  hflux(ii) = zero
!               end do
!            end do
!         end if bocoType3
!      end do
!   end do

! end subroutine getHeatFluxCellCenter

