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

       call vecGetArrayF90(zipper%localVal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end do dimLoop

  end subroutine getForces

subroutine getHeatFlux(hflux, npts, sps)
  use constants
  use blockPointers, only : nDom, nBocos, BCType, BCData
  use surfaceFamilies, only : BCFamExchange, familyExchange, &
       zeroCellVal, zeroNodeVal
  use utils, only : setPointers
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts, sps
  real(kind=realType), intent(out) :: hflux(npts)

  integer(kind=intType) :: mm, nn, i, j, ii
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  type(familyExchange), pointer :: exch

  exch => BCFamExchange(iBCGroupWalls, sps)
  do nn=1, nDom
     call setPointers(nn, 1_intType, sps)
     call heatFluxes()

     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

        bocoType1: if (BCType(mm) == NSWallIsoThermal) then
           BCData(mm)%cellVal => BCData(mm)%area(:, :)
        else if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic) then
           BCData(mm)%cellVal => zeroCellVal
           BCData(mm)%nodeVal => zeroNodeVal
        end if bocoType1
     end do
  end do

  call computeWeighting(exch)

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

subroutine heatFluxes
  use constants
  use blockPointers, only : BCData, nDom, nBocos, BCType, BCFaceID, viscSubFace
  use BCPointers, only : ssi
  use flowVarRefState, only : pRef, rhoRef
  use utils, only : setPointers, setBCPointers
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, ii, mm
  real(kind=realType) :: fact, scaleDim, Q
  real(kind=realType) :: qw, qA
  logical :: heatedSubface

  ! Set the actual scaling factor such that ACTUAL heat flux is computed
  ! The factor is determined from stanton number
  scaleDim = pRef*sqrt(pRef/rhoRef)

  ! Loop over the boundary subfaces of this block.
  bocos: do mm=1, nBocos

     ! Only do this on isoThermalWalls
     if (BCType(mm) == NSWallIsoThermal) then

        ! Set a bunch of pointers depending on the face id to make
        ! a generic treatment possible. The routine setBcPointers
        ! is not used, because quite a few other ones are needed.
        call setBCPointers(mm, .True.)

        select case (BCFaceID(mm))
        case (iMin, jMin, kMin)
           fact = -one
        case (iMax, jMax, kMax)
           fact = one
        end select

        ! Loop over the quadrilateral faces of the subface. Note that
        ! the nodal range of BCData must be used and not the cell
        ! range, because the latter may include the halo's in i and
        ! j-direction. The offset +1 is there, because inBeg and jnBeg
        ! refer to nodal ranges and not to cell ranges.
        !
        do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
           do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd

              ! Compute the normal heat flux on the face. Inward positive.
              BCData(mm)%cellHeatFlux(i,j) = -fact*scaleDim* &
                   sqrt(ssi(i,j,1)**2 + ssi(i,j,2)**2 + ssi(i,j,3)**2) * &
                   ( viscSubface(mm)%q(i,j,1)*BCData(mm)%norm(i,j,1) &
                   + viscSubface(mm)%q(i,j,2)*BCData(mm)%norm(i,j,2) &
                   + viscSubface(mm)%q(i,j,3)*BCData(mm)%norm(i,j,3))

           enddo
        end do
     end if
  enddo bocos

end subroutine heatFluxes
