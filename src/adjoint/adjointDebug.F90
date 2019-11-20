! This is a special function that is sued to alloc derivative values
! in blockpointers_d for use with the AD code.
module adjointDebug

contains

#ifndef USE_COMPLEX

   subroutine computeMatrixFreeProductFwdFD(xvdot, extradot, wdot, bcDataValuesdot,&
      useSpatial, useState, famLists,&
      bcDataNames, bcDataValues, bcDataFamLists, bcVarsEmpty,&
      dwdot, funcsDot, fDot, hfdot, &
      costSize, fSize, nTime, h)

      ! This routine is used to debug master_d. It uses the forward seeds to set perturbations
      ! and then computes the value of the derivatives using forward finite diffenece


      use constants
      use adjointvars
      use blockPointers, only : nDom
      use communication, only : adflow_comm_world
      use inputTimeSpectral, only : nTimeIntervalsSpectral
      use inputPhysics, only :pointRefd, alphad, betad, equations, machCoefd, &
         machd, machGridd, rgasdimd
      use iteration, only : currentLevel, groundLevel
      use flowVarRefState, only : pInfDimd, rhoInfDimd, TinfDimd
      use blockPointers, only : nDom, il, jl, kl, wd, x, w, dw, dwd, nBocos, nViscBocos

      use adjointUtils, only : allocDerivativeValues, zeroADSeeds
      use masterRoutines, only : master
      use utils, only : isWallType, setPointers, setPointers_d, EChk
      use flowVarRefState, only : nw, nwf
      use wallDistanceData, only: xSurf, xSurfVec
      implicit none

      ! Input Variables
      real(kind=realType), dimension(:), intent(in) :: xvdot
      real(kind=realType), dimension(:), intent(in) :: extradot
      real(kind=realType), dimension(:), intent(in) :: wdot
      logical, intent(in) :: useSpatial, useState
      integer(kind=intType), dimension(:, :) :: famLists
      integer(kind=intType) :: costSize, fSize, nTime

      ! character, dimension(:, :), intent(in) :: bcDataNames
      ! real(kind=realType), dimension(:), intent(in) :: bcDataValues, bcDataValuesDot
      ! integer(kind=intType), dimension(:, :) :: bcDataFamLists
      character, dimension(:), intent(in) :: bcDataNames
      real(kind=realType), dimension(:,:), intent(inout) :: bcDataValues
      real(kind=realType), dimension(:,:), intent(in) :: bcDataValuesDot
      integer(kind=intType), dimension(:) :: bcDataFamLists
      logical, intent(in) :: BCVarsEmpty
      real(kind=realType), intent(in) :: h ! step size for Finite Difference

      ! Ouput Variables
      real(kind=realType), dimension(size(wdot)), intent(out) :: dwDot
      real(kind=realType), dimension(costSize, size(famLists,1)), intent(out) :: funcsDot
      real(kind=realType), dimension(3, fSize, nTime), intent(out) :: fDot
      real(kind=realType), dimension(1, fSize, nTime), intent(out) :: hfDot

      ! Working Variables
      integer(kind=intType) :: nn,sps, level
      integer(kind=intType) :: ierr, mm,i,j,k, l,  ii, jj, iRegion

      real(kind=realType), dimension(costSize, size(famLists,1)) :: funcs


      ! Input Arguments for master:
      real(kind=realType), dimension(costSize, size(famLists,1)) :: funcValues


         ! Working Variables
      real(kind=realType), dimension(:, :, :), allocatable :: forces
      real(kind=realType), dimension(:, :, :), allocatable :: heatfluxes

      fSize = size(fDot, 2)
      allocate(forces(3, fSize, nTimeIntervalsSpectral))

      fSize = size(hfDot, 2)
      allocate(heatfluxes(1, fSize, nTimeIntervalsSpectral))



      ! Need to trick the residual evalution to use coupled (mean flow and
      ! turbulent) together.
      level = 1
      currentLevel = level
      groundLevel = level

      ! Allocate the memory we need for derivatives if not done so
      ! already. Note this isn't deallocated until the adflow is
      ! destroyed.
      if (.not. derivVarsAllocated) then
         call allocDerivativeValues(level)
      end if

      ! Zero all AD seesd.
      do nn=1,nDom
         do sps=1,nTimeIntervalsSpectral
            call zeroADSeeds(nn,level, sps)
         end do
      end do

      ! Set the extra seeds now do the extra ones. Note that we are assuming the
      ! machNumber used for the coefficients follows the Mach number,
      ! not the grid mach number.
      alphad = extraDot(iAlpha)
      betad = extraDot(iBeta)
      machd = extraDot(iMach)
      machCoefd = extraDot(iMach)
      machGridd = extraDot(iMachGrid)
      PinfDimd = extraDot(iPressure)
      rhoinfDimd = extraDot(iDensity)
      tinfdimd = extraDot(iTemperature)
      pointrefd(1) = extraDot(iPointRefX)
      pointrefd(2) = extraDot(iPointRefY)
      pointrefd(3) = extraDot(iPointRefZ)
      rgasdimd = zero



      ! ----------------------------- Run Master ---------------------------------
      ! Run the super-dee-duper master rotuine
      if (bcVarsEmpty) then
         call master(useSpatial, famLists, funcValues, forces, heatfluxes)
      else
         call master(useSpatial, famLists, funcValues, forces, heatfluxes, &
                     bcDataNames, bcDataValues, bcDataFamLists)
      end if



      ! Copy out the residual derivative into the provided dwDot
      ii =0
      do nn=1, nDom
         do sps=1,nTimeIntervalsSpectral
            call setPointers_d(nn, 1, sps)
            do k=2, kl
               do j=2, jl
                  do i=2, il
                     do l=1, nw
                        ii = ii + 1
                        dwd(i,j,k,l) = -dw(i,j,k,l)
                     end do
                  end do
               end do
            end do
         end do
      end do

      fDot = -forces
      hfDot = -heatfluxes
      funcsDot = -funcValues

      !  --------------------- apply the perturbations ----------------------------
      ! Set the provided w and x seeds:
      ii = 0
      jj = 0
      domainLoop1: do nn=1,nDom
         spectalLoop1: do sps=1,nTimeIntervalsSpectral
            call setPointers(nn, 1, sps)
            do k=1, kl
               do j=1,jl
                  do i=1,il
                     do l=1,3
                        ii = ii + 1
                        x(i, j, k, l) = x(i, j, k, l) +  xvdot(ii)* h
                     end do
                  end do
               end do
            end do
            do k=2, kl
               do j=2,jl
                  do i=2,il
                     do l = 1, nw
                        jj = jj + 1
                        w(i, j, k, l) = w(i, j, k, l) +  wDot(jj)*h
                     end do
                  end do
               end do
            end do
         end do spectalLoop1
      end do domainLoop1

      write(*,*) '1 bcDataValues ', bcDataValues
      if (.not. bcVarsEmpty) then
         bcDataValues = bcDataValues + bcDataValuesDot*h
      endif 
      write(*,*) '2 bcDataValues ', bcDataValues


      ! ----------------------------- Run Master ---------------------------------
      ! Run the super-dee-duper master rotuine
      if (bcVarsEmpty) then
         call master(useSpatial, famLists, funcValues, forces, heatfluxes)
      else
         call master(useSpatial, famLists, funcValues, forces, heatfluxes, &
                     bcDataNames, bcDataValues, bcDataFamLists)
      end if


      ! Copy out the residual derivative into the provided dwDot and remove the
      ! perturbation
      ii = 0
      jj = 0
      do nn=1, nDom
         do sps=1,nTimeIntervalsSpectral
            call setPointers_d(nn, 1, sps)
            do k=1, kl
               do j=1,jl
                  do i=1,il
                     do l=1,3
                        ii = ii + 1
                        x(i, j, k, l) = x(i, j, k, l) -  xvdot(ii)* h
                     end do
                  end do
               end do
            end do
            do k=2, kl
               do j=2, jl
                  do i=2, il
                     do l=1, nw
                        jj = jj + 1
                        w(i, j, k, l) = w(i, j, k, l) -  wDot(jj)*h
                        dwd(i,j,k,l) = (dwd(i,j,k,l) + dw(i,j,k,l))/h
                        dwdot(jj) = dwd(i,j,k,l) ! copy values to output
                     end do
                  end do
               end do
            end do

         end do
      end do


      if (.not. bcVarsEmpty) then
         bcDataValues = bcDataValues - bcDataValuesDot*h
      endif 

      write(*,*) '1 heatflux', -hfDot
      write(*,*) '2 heatflux', heatfluxes
      fDot = (fDot + forces)/h
      hfDot = (hfDot + heatfluxes)/h
      funcsDot = (funcsDot + funcValues)/h
      write(*,*) 'FD fDot', minval(fDot), maxval(fDot)
      write(*,*) 'FD hfDot', minval(hfDot), maxval(hfDot)
      write(*,*) 'FD funcsDot', funcsDot

   end subroutine computeMatrixFreeProductFwdFD


   subroutine printADSeeds(nn, level, sps)

      use constants
      use block, only : flowDomsd, flowDoms
      use blockPointers
      use inputTimeSpectral
      use flowVarRefState
      use inputPhysics
      use BCPointers_b
      use communication
      use oversetData, only : oversetPresent
      use cgnsGrid, only : cgnsDoms, cgnsDomsd, cgnsNDom
      use actuatorRegionData, only : nActuatorRegions, actuatorRegionsd
      implicit none

      ! Input parameters
      integer(kind=intType) :: nn, level, sps

      ! Working parameters
      integer(kind=intType) :: mm, i, iDom
      integer(kind=intType) :: iBoco, iData, iDirichlet
      write(*,*) 'd2wall ', minval(flowDomsd(nn, level, sps)%d2wall), &
                        maxval(flowDomsd(nn, level, sps)%d2wall)
      write(*,*) 'x ', minval(flowDomsd(nn, level, sps)%x), &
                        maxval(flowDomsd(nn, level, sps)%x)
      write(*,*) 'si ', minval(flowDomsd(nn, level, sps)%si), &
                        maxval(flowDomsd(nn, level, sps)%si)
      write(*,*) 'sj ', minval(flowDomsd(nn, level, sps)%sj), &
                        maxval(flowDomsd(nn, level, sps)%sj)
      write(*,*) 'sk ', minval(flowDomsd(nn, level, sps)%sk), &
                        maxval(flowDomsd(nn, level, sps)%sk)
      write(*,*) 'vol ', minval(flowDomsd(nn, level, sps)%vol), &
                        maxval(flowDomsd(nn, level, sps)%vol)

      write(*,*) 's ', minval(flowDomsd(nn, level, sps)%s), &
                        maxval(flowDomsd(nn, level, sps)%s)
      write(*,*) 'sFaceI ', minval(flowDomsd(nn, level, sps)%sFaceI), &
                        maxval(flowDomsd(nn, level, sps)%sFaceI)
      write(*,*) 'sFaceJ ', minval(flowDomsd(nn, level, sps)%sFaceJ), &
                        maxval(flowDomsd(nn, level, sps)%sFaceJ)
      write(*,*) 'sFaceK ', minval(flowDomsd(nn, level, sps)%sFaceK), &
                        maxval(flowDomsd(nn, level, sps)%sFaceK)

      write(*,*) 'w ', minval(flowDomsd(nn, level, sps)%w), &
                        maxval(flowDomsd(nn, level, sps)%w)
      write(*,*) 'dw ', minval(flowDomsd(nn, level, sps)%dw), &
                        maxval(flowDomsd(nn, level, sps)%dw)
      write(*,*) 'fw ', minval(flowDomsd(nn, level, sps)%fw), &
                        maxval(flowDomsd(nn, level, sps)%fw)
      write(*,*) 'scratch ', minval(flowDomsd(nn, level, sps)%scratch), &
                        maxval(flowDomsd(nn, level, sps)%scratch)

      write(*,*) 'p ', minval(flowDomsd(nn, level, sps)%p), &
                        maxval(flowDomsd(nn, level, sps)%p)
      write(*,*) 'gamma ', minval(flowDomsd(nn, level, sps)%gamma), &
                        maxval(flowDomsd(nn, level, sps)%gamma)
      write(*,*) 'aa ', minval(flowDomsd(nn, level, sps)%aa), &
                        maxval(flowDomsd(nn, level, sps)%aa)

      write(*,*) 'rlv ', minval(flowDomsd(nn, level, sps)%rlv), &
                        maxval(flowDomsd(nn, level, sps)%rlv)
      write(*,*) 'rev ', minval(flowDomsd(nn, level, sps)%rev), &
                        maxval(flowDomsd(nn, level, sps)%rev)

      write(*,*) 'radI ', minval(flowDomsd(nn, level, sps)%radI), &
                        maxval(flowDomsd(nn, level, sps)%radI)
      write(*,*) 'radJ ', minval(flowDomsd(nn, level, sps)%radJ), &
                        maxval(flowDomsd(nn, level, sps)%radJ)
      write(*,*) 'radK ', minval(flowDomsd(nn, level, sps)%radK), &
                        maxval(flowDomsd(nn, level, sps)%radK)

      write(*,*) 'ux ', minval(flowDomsd(nn, level, sps)%ux), &
                        maxval(flowDomsd(nn, level, sps)%ux)
      write(*,*) 'uy ', minval(flowDomsd(nn, level, sps)%uy), &
                        maxval(flowDomsd(nn, level, sps)%uy)
      write(*,*) 'uz ', minval(flowDomsd(nn, level, sps)%uz), &
                        maxval(flowDomsd(nn, level, sps)%uz)
      write(*,*) 'vx ', minval(flowDomsd(nn, level, sps)%vx), &
                        maxval(flowDomsd(nn, level, sps)%vx)
      write(*,*) 'vy ', minval(flowDomsd(nn, level, sps)%vy), &
                        maxval(flowDomsd(nn, level, sps)%vy)
      write(*,*) 'vz ', minval(flowDomsd(nn, level, sps)%vz), &
                        maxval(flowDomsd(nn, level, sps)%vz)
      write(*,*) 'wx ', minval(flowDomsd(nn, level, sps)%wx), &
                        maxval(flowDomsd(nn, level, sps)%wx)
      write(*,*) 'wy ', minval(flowDomsd(nn, level, sps)%wy), &
                        maxval(flowDomsd(nn, level, sps)%wy)
      write(*,*) 'wz ', minval(flowDomsd(nn, level, sps)%wz), &
                        maxval(flowDomsd(nn, level, sps)%wz)
      write(*,*) 'qx ', minval(flowDomsd(nn, level, sps)%qx), &
                        maxval(flowDomsd(nn, level, sps)%qx)
      write(*,*) 'qy ', minval(flowDomsd(nn, level, sps)%qy), &
                        maxval(flowDomsd(nn, level, sps)%qy)
      write(*,*) 'qz ', minval(flowDomsd(nn, level, sps)%qz), &
                        maxval(flowDomsd(nn, level, sps)%qz)

      write(*,*) 'bmti1 ',minval(flowDomsd(nn, level, sps)%bmti1), &
                              maxval(flowDomsd(nn, level, sps)%bmti1)
      write(*,*) 'bmti2 ',minval(flowDomsd(nn, level, sps)%bmti2), &
                              maxval(flowDomsd(nn, level, sps)%bmti2)
      write(*,*) 'bmtj1 ',minval(flowDomsd(nn, level, sps)%bmtj1), &
                              maxval(flowDomsd(nn, level, sps)%bmtj1)
      write(*,*) 'bmtj2 ',minval(flowDomsd(nn, level, sps)%bmtj2), &
                              maxval(flowDomsd(nn, level, sps)%bmtj2)
      write(*,*) 'bmtk1 ',minval(flowDomsd(nn, level, sps)%bmtk1), &
                              maxval(flowDomsd(nn, level, sps)%bmtk1)
      write(*,*) 'bmtk2 ',minval(flowDomsd(nn, level, sps)%bmtk2), &
                              maxval(flowDomsd(nn, level, sps)%bmtk2)
      write(*,*) 'bvti1 ',minval(flowDomsd(nn, level, sps)%bvti1), &
                              maxval(flowDomsd(nn, level, sps)%bvti1)
      write(*,*) 'bvti2 ',minval(flowDomsd(nn, level, sps)%bvti2), &
                              maxval(flowDomsd(nn, level, sps)%bvti2)
      write(*,*) 'bvtj1 ',minval(flowDomsd(nn, level, sps)%bvtj1), &
                              maxval(flowDomsd(nn, level, sps)%bvtj1)
      write(*,*) 'bvtj2 ',minval(flowDomsd(nn, level, sps)%bvtj2), &
                              maxval(flowDomsd(nn, level, sps)%bvtj2)
      write(*,*) 'bvtk1 ',minval(flowDomsd(nn, level, sps)%bvtk1), &
                              maxval(flowDomsd(nn, level, sps)%bvtk1)
      write(*,*) 'bvtk2 ',minval(flowDomsd(nn, level, sps)%bvtk2), &
                              maxval(flowDomsd(nn, level, sps)%bvtk2)

      bocoLoop: do mm=1, flowDoms(nn, level, sps)%nBocos
         write(*,*) 'mm', mm, 'BCData(mm)%norm',minval(flowDomsd(nn, level, sps)%BCData(mm)%norm), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%norm)
         write(*,*) 'mm', mm, 'bcData(mm)%rface ',minval(flowDomsd(nn, level, sps)%bcData(mm)%rface), &
                                 maxval(flowDomsd(nn, level, sps)%bcData(mm)%rface)
         write(*,*) 'mm', mm, 'bcData(mm)%Fv ',minval(flowDomsd(nn, level, sps)%bcData(mm)%Fv), &
                                 maxval(flowDomsd(nn, level, sps)%bcData(mm)%Fv)
         write(*,*) 'mm', mm, 'bcData(mm)%Fp ',minval(flowDomsd(nn, level, sps)%bcData(mm)%Fp), &
                                 maxval(flowDomsd(nn, level, sps)%bcData(mm)%Fp)
         write(*,*) 'mm', mm, 'bcData(mm)%Tv ',minval(flowDomsd(nn, level, sps)%bcData(mm)%Tv), &
                                 maxval(flowDomsd(nn, level, sps)%bcData(mm)%Tv)
         write(*,*) 'mm', mm, 'bcData(mm)%Tp ',minval(flowDomsd(nn, level, sps)%bcData(mm)%Tp), &
                                 maxval(flowDomsd(nn, level, sps)%bcData(mm)%Tp)
         write(*,*) 'mm', mm, 'bcData(mm)%area ',minval(flowDomsd(nn, level, sps)%bcData(mm)%area), &
                                 maxval(flowDomsd(nn, level, sps)%bcData(mm)%area)
         write(*,*) 'mm', mm, 'BCData(mm)%uSlip ',minval(flowDomsd(nn, level, sps)%BCData(mm)%uSlip), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%uSlip)
         write(*,*) 'mm', mm, 'BCData(mm)%TNS_Wall ',minval(flowDomsd(nn, level, sps)%BCData(mm)%TNS_Wall), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%TNS_Wall)
         write(*,*) 'mm', mm, 'BCData(mm)%ptInlet ',minval(flowDomsd(nn, level, sps)%BCData(mm)%ptInlet), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%ptInlet)
         write(*,*) 'mm', mm, 'BCData(mm)%htInlet ',minval(flowDomsd(nn, level, sps)%BCData(mm)%htInlet), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%htInlet)
         write(*,*) 'mm', mm, 'BCData(mm)%ttInlet ',minval(flowDomsd(nn, level, sps)%BCData(mm)%ttInlet), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%ttInlet)
         write(*,*) 'mm', mm, 'BCData(mm)%turbInlet ',minval(flowDomsd(nn, level, sps)%BCData(mm)%turbInlet), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%turbInlet)
         write(*,*) 'mm', mm, 'BCData(mm)%ps ',minval(flowDomsd(nn, level, sps)%BCData(mm)%ps), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%ps)
         write(*,*) 'mm', mm, 'BCData(mm)%cellHeatFlux ',minval(flowDomsd(nn, level, sps)%BCData(mm)%cellHeatFlux), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%cellHeatFlux)
         write(*,*) 'mm', mm, 'BCData(mm)%nodeHeatFlux ',minval(flowDomsd(nn, level, sps)%BCData(mm)%nodeHeatFlux), &
                                 maxval(flowDomsd(nn, level, sps)%BCData(mm)%nodeHeatFlux)

      end do bocoLoop


      viscbocoLoop: do mm=1,flowDoms(nn, level, sps)%nViscBocos
         write(*,*) 'mm', mm, 'viscSubface(mm)%tau ',minval(flowDomsd(nn, level, sps)%viscSubface(mm)%tau), &
                                 maxval(flowDomsd(nn, level, sps)%viscSubface(mm)%tau)
         write(*,*) 'mm', mm, 'viscSubface(mm)%q ',minval(flowDomsd(nn, level, sps)%viscSubface(mm)%q), &
                                 maxval(flowDomsd(nn, level, sps)%viscSubface(mm)%q)
      end do viscbocoLoop

      ! For overset, the weights may be active in the comm structure. We
      ! need to zero them before we can accumulate.
      if (oversetPresent) then
         ! Pointers to the overset comms to make it easier to read
         sends: do i=1,commPatternOverset(level, sps)%nProcSend
            write(*,*) 'commPatternOverset(level, sps)%sendList(i)%interpd ',&
                        minval(commPatternOverset(level, sps)%sendList(i)%interpd), &
                        maxval(commPatternOverset(level, sps)%sendList(i)%interpd)
         end do sends
         write(*,*) 'internalOverset(level, sps)%donorInterpd ',minval(internalOverset(level, sps)%donorInterpd), &
                                 maxval(internalOverset(level, sps)%donorInterpd)
      end if

      write(*,*) 'alphad ',alphad
      write(*,*) 'betad ',betad
      write(*,*) 'machd ',machd
      write(*,*) 'machGridd ',machGridd
      write(*,*) 'machCoefd ',machCoefd
      write(*,*) 'pinfdimd ',pinfdimd
      write(*,*) 'tinfdimd ',tinfdimd
      write(*,*) 'rhoinfdimd ',rhoinfdimd
      write(*,*) 'rgasdimd ',rgasdimd
      write(*,*) 'pointrefd ',pointrefd
      write(*,*) 'prefd ',prefd
      write(*,*) 'rhoRefd ',rhoRefd
      write(*,*) 'Trefd ',Trefd
      write(*,*) 'murefd ',murefd
      write(*,*) 'urefd ',urefd
      write(*,*) 'hrefd ',hrefd
      write(*,*) 'timerefd ',timerefd
      write(*,*) 'pinfd ',pinfd
      write(*,*) 'pinfCorrd ',pinfCorrd
      write(*,*) 'rhoinfd ',rhoinfd
      write(*,*) 'uinfd ',uinfd
      write(*,*) 'rgasd ',rgasd
      write(*,*) 'muinfd ',muinfd
      write(*,*) 'gammainfd ',gammainfd
      write(*,*) 'winfd ',winfd
      write(*,*) 'veldirfreestreamd ',veldirfreestreamd
      write(*,*) 'liftdirectiond ',liftdirectiond
      write(*,*) 'dragdirectiond ',dragdirectiond

      ! Zero all the reverse seeds in the dirichlet input arrays
      do iDom=1, cgnsNDom
         do iBoco=1, cgnsDoms(iDom)%nBocos
            if (associated(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet)) then
               do iData=1, size(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet)
                  if (associated(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)) then
                     do iDirichlet = 1, size(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)
                        write(*,*) iDom, iBoco, iData, iDirichlet, 'dataArr(:) '&
                        ,cgnsDomsd(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays(iDirichlet)%dataArr(:)
                     end do
                  end if
               end do
            end if
         end do
      end do

      ! And the reverse seeds in the actuator zones
      do i=1, nActuatorRegions
         write(*,*) 'actuatorRegionsd(i)%F ',actuatorRegionsd(i)%F
         write(*,*) 'actuatorRegionsd(i)%T ',actuatorRegionsd(i)%T
      end do

   end subroutine printADSeeds


   subroutine computeDotProductTest(wdot, xvdot, bcdatavaluesdot, &
      fbar, hfbar, dwbar, funcsbar, &
      famlists, BCVarsEmpty, bcdatanames, bcdatavalues, bcdatafamlists, &
      costSize, fSize, nTime,  spatialSize, extraSize, stateSize)

      ! used to preform the dot product test on the master subroutines
      !
      !

      use constants
      use block, only : blockType, flowDomsd
      use communication, only : adflow_comm_world
      use blockPointers, only : nDom, dwd, il, jl, kl
      use inputTimeSpectral, only : nTimeIntervalsSpectral
      use inputPhysics, only : equations
      use iteration, only : currentLevel, groundLevel
      use flowVarRefState, only : nw, nwf
      use inputAdjoint, only : frozenTurbulence
      use ADjointPETSc, only : x_like, psi_like3
      use adjointvars, only : derivVarsAllocated
      use utils, only : setPointers_d, EChk
      use masterRoutines, only : master, master_b, master_d
      use blockPointers, only : nDom, nBocos, BCType, BCData, BCDatad
      use inputPhysics, only : pointRefd, alphad, betad, equations, machCoefd, &
      machd, machGridd, rgasdimd
      use flowVarRefState, only : pInfDimd, rhoInfDimd, TinfDimd
      use adjointUtils, only : allocDerivativeValues, zeroADSeeds
      implicit none

      real(kind=realtype), dimension(:), intent(in) :: wdot
      real(kind=realtype), dimension(:), intent(in) :: xvdot
      real(kind=realtype), dimension(:,:), intent(in) :: bcdatavaluesdot

      real(kind=realtype), dimension(:,:,:), intent(in) :: fbar
      real(kind=realtype), dimension(:,:,:), intent(in) :: hfbar
      real(kind=realType), dimension(:), intent(in) :: dwBar

      real(kind=realType), dimension(:, :), intent(in) :: funcsbar

      integer(kind=inttype), dimension(:,:) :: famlists
      ! character, dimension(:,:), intent(in) :: bcdatanames
      ! real(kind=realtype), dimension(:),intent(in) :: bcdatavalues
      ! integer(kind=inttype), dimension(:,:), intent(in) :: bcdatafamlists

      character, dimension(:), intent(in) :: bcdatanames
      real(kind=realtype), dimension(:,:),intent(in) :: bcdatavalues
      integer(kind=inttype), dimension(:), intent(in) :: bcdatafamlists
      
      logical, intent(in) :: BCVarsEmpty


      integer(kind=intType), intent(in) :: stateSize, extraSize, spatialSize
      integer(kind=intType), intent(in) :: costSize, fSize, nTime

      !--------------------- working variables ---------------------------------

      ! reverse mode output
      real(kind=realType), dimension(stateSize) :: wbar
      real(kind=realType), dimension(extraSize) :: extrabar
      real(kind=realType), dimension(spatialSize) :: xvbar
      real(kind=realType), dimension(size(bcDataValues,1), size(bcDataValues,2)) :: bcDataValuesbar
      real(kind=realType), dimension(size(funcsBar,1), size(funcsBar, 2)) :: funcValues


      ! forward mode output
      real(kind=realType), dimension(size(wdot)) :: dwDot
      real(kind=realType), dimension(costSize, size(famLists,1)) :: funcsDot
      real(kind=realType), dimension(3, fSize, nTime) :: fDot
      real(kind=realType), dimension(1, fSize, nTime) :: hfDot


      integer(kind=intType):: npts, sps, nn, level, mm, ierr, iim, nState
      real(kind=realType) :: hflux(10),  hflux_d(10), hflux_b(10)

      real(kind=realType) :: cellhf_b(4), cellhf_d(4)
      real(kind=realType) :: area_b(4), area_d(4)
      real(kind=realType) :: nodehf_b(10), nodehf_d(10)


      real(kind=realType), dimension(:), pointer ::  localPtr
      real(kind=realType), dimension(:), pointer :: nodeValLocPtr,   sumGlobalPtr,&
                                                    nodeValLocPtr_b, sumGlobalPtr_b

       ! Working Variables
      real(kind=realType), dimension(:, :, :), allocatable :: forces
      real(kind=realType), dimension(:, :, :), allocatable :: heatfluxes

      type(blockType) :: tmpDomsd(nDom, 1, nTimeIntervalsSpectral)



      real(kind=realType), dimension(:,:,:,:), allocatable  :: x, xtmp

      real(kind=realType), dimension(:,:,:,:,:), allocatable  :: Xold
      real(kind=realType), dimension(:,:,:,:), allocatable :: sI, sJ, sK
      real(kind=realType), dimension(:,:,:), allocatable :: tau, q
      real(kind=realType), dimension(:,:), allocatable :: cellHeatFlux, area, Fv, Fp


      allocate(forces(3, fSize, nTimeIntervalsSpectral))

      allocate(heatfluxes(1, fSize, nTimeIntervalsSpectral))

      if ( frozenTurbulence ) then
         nState = nwf
      else
         nState = nw
      endif

      sps = 1

      level = 1
      currentLevel = level
      groundLevel = level

         if (.not. derivVarsAllocated) then
            call allocDerivativeValues(level)
         end if

         ! Zero all AD seesd.
         do nn=1,nDom
            call zeroADSeeds(nn,level, sps)
         end do

         ! Set the extra seeds now do the extra ones. Note that we are assuming the
         ! machNumber used for the coefficients follows the Mach number,
         ! not the grid mach number.
         alphad = 0
         betad = 0
         machd = 0
         machCoefd = 0
         machGridd = 0
         PinfDimd = 0
         rhoinfDimd = 0
         tinfdimd = 0
         pointrefd(1) = 0
         pointrefd(2) = 0
         pointrefd(3) = 0
         rgasdimd = zero


         ! ----------------------------- Run Master ---------------------------------
         ! Run the super-dee-duper master rotuine
         ! if (bcVarsEmpty) then
         !    call master(.true., famLists, funcValues, forces, heatfluxes)
         ! else
         !    call master(.true., famLists, funcValues, forces, heatfluxes, &
         !                bcDataNames, bcDataValues, bcDataFamLists)
         ! end if


         !BWD
         ! Run the super-dee-duper master forward rotuine
         if (bcVarsEmpty) then
            call master_d(wDot, xVDot, fDot, hfDot, dwDot, famLists, funcValues, funcsDot)
         else
            call master_d(wDot, xVDot, fDot, hfDot, dwDot, &
               famLists, funcValues, funcsDot, bcDataNames, bcDataValues, bcDataValuesdot, bcDataFamLists)
         end if
         nn = 1
         level = 1
         sps = 1
         mm = 1

         ! tmpDomsd = flowdomsd

         ! allocate(sk, source = flowDomsd(nn, level, sps)%sk)
         ! allocate(q, source = flowDomsd(nn, level, sps)%viscSubface(mm)%q)
         ! allocate(area, source = flowDomsd(nn, level, sps)%BCData(mm)%area)
         ! allocate(cellHeatFlux, source = flowDomsd(nn, level, sps)%BCData(mm)%cellHeatFlux)


         do nn=1,nDom
            call zeroADSeeds(nn,level, sps)
         end do

         write(*,*) 'wbar',wbar
         write(*,*)
         write(*,*) 'xvbar',xvbar
         write(*,*)
         write(*,*) 'extraBar',extraBar
         write(*,*)
         write(*,*) 'fBar',fBar
         write(*,*)
         write(*,*) 'hfbar',hfbar
         write(*,*)
         write(*,*) 'dwbar',dwbar
         write(*,*)

         if (bcVarsEmpty) then
            call master_b(wbar, xvbar, extraBar, fBar, hfbar, dwbar, nState, famLists, &
               funcValues, funcsBar)
         else
            call master_b(wbar, xvbar, extraBar, fBar, hfbar, dwbar, nState, famLists, &
               funcValues, funcsBar, bcDataNames, bcDataValues, bcDataValuesbar, bcDataFamLists)
         end if

         write(*,*) 'wbar',wbar
         write(*,*)
         write(*,*) 'xvbar',xvbar
         write(*,*)
         write(*,*) 'extraBar',extraBar
         write(*,*)
         write(*,*) 'fBar',fBar
         write(*,*)
         write(*,*) 'hfbar',hfbar
         write(*,*)
         write(*,*) 'dwbar',dwbar
         write(*,*)

         ! ===================================================================
         ! check getHeatFluxes
         !===================================================================

         nn = 1
         mm = 1
         level = 1
         sps = 1

            write(*,*) '--------------------dot prod--------------------------------'

            ! write(*,*) sum(flowDomsd(nn, level, sps)%sk* tmpDomsd(nn, level, sps)%sk)
            ! write(*,*) sum(sk* flowDomsd(nn, level, sps)%sk)
            ! write(*,*) sum(q* flowDomsd(nn, level, sps)%viscSubface(mm)%q)
            ! write(*,*) sum(sk* flowDomsd(nn, level, sps)%sk) + sum(q* flowDomsd(nn, level, sps)%viscSubface(mm)%q)

            ! write(*,*) sum(cellHeatFlux* flowDomsd(nn, level, sps)%BCData(mm)%cellHeatFlux)
            ! write(*,*) sum(area* flowDomsd(nn, level, sps)%BCData(mm)%area)
            ! write(*,*) sum(cellHeatFlux* flowDomsd(nn, level, sps)%BCData(mm)%cellHeatFlux) + &
            ! sum(area* flowDomsd(nn, level, sps)%BCData(mm)%area)

            ! write(*,*) 'comparison'
            ! ! write(*,*) 'f', sum(fDot*fbar)
            ! write(*,*) 'hf', sum(hfDot*hfbar)
            ! write(*,*) 'hf'
            ! write(*,*) 'fwd', hfDot
            ! write(*,*) 'bwd', hfbar
            ! ! write(*,*) 'dw', sum(dwDot*dwbar)

            write(*,*) 'w'!, sum(wDot*wbar)
            write(*,*) 'fwd', wDot
            write(*,*) 'bwd', wbar

            write(*,*) 'xv'!, sum(xvDot*xvbar)
            write(*,*) 'fwd', xvDot
            write(*,*) 'bwd', xvbar


            write(*,*)  sum(wDot*wbar) + sum(xVDot*xvbar), sum(fDot*fbar) + sum(hfDot*hfbar) + sum(dwDot*dwbar)
            write(*,*)  sum(wDot*wbar) + sum(xVDot*xvbar)- (sum(fDot*fbar) + sum(hfDot*hfbar) + sum(dwDot*dwbar))
   end subroutine computeDotProductTest

#else
      
   subroutine computeMatrixFreeProductFwdCS(xvdot, extradot, wdot, bcDataValuesdot,&
      useSpatial, useState, famLists,&
      bcDataNames, bcDataValues, bcDataFamLists, bcVarsEmpty,&
      dwdot, funcsDot, fDot, hfdot, &
      costSize, fSize, nTime)

      ! This routine is used to debug master_d. It uses the forward seeds to set perturbations
      ! and then computes the value of the derivatives using forward finite diffenece


      use constants
      use adjointvars
      use blockPointers, only : nDom
      use communication, only : adflow_comm_world
      use inputTimeSpectral, only : nTimeIntervalsSpectral
      use inputPhysics, only :pointRef, alpha, beta, equations, machCoef, &
         mach, machGrid, rgasdim
      use iteration, only : currentLevel, groundLevel
      use flowVarRefState, only : pInfDim, rhoInfDim, TinfDim
      use blockPointers, only : nDom, il, jl, kl, wd, x, w, dw, dwd, nBocos, nViscBocos

      use adjointUtils, only : allocDerivativeValues, zeroADSeeds
      use masterRoutines, only : master
      use utils, only : isWallType, setPointers, setPointers_d, EChk
      use flowVarRefState, only : nw, nwf
      use wallDistanceData, only: xSurf, xSurfVec
      implicit none

      ! Input Variables
      complex(kind=realType), dimension(:), intent(in) :: xvdot
      complex(kind=realType), dimension(:), intent(in) :: extradot
      complex(kind=realType), dimension(:), intent(in) :: wdot
      logical, intent(in) :: useSpatial, useState
      integer(kind=intType), dimension(:, :) :: famLists
      integer(kind=intType) :: costSize, fSize, nTime

      character, dimension(:, :), intent(in) :: bcDataNames
      complex(kind=realType), dimension(:), intent(in) :: bcDataValues, bcDataValuesDot
      integer(kind=intType), dimension(:, :) :: bcDataFamLists
      logical, intent(in) :: BCVarsEmpty

      ! Ouput Variables
      complex(kind=realType), dimension(size(wdot)), intent(out) :: dwDot
      complex(kind=realType), dimension(costSize, size(famLists,1)), intent(out) :: funcsDot
      complex(kind=realType), dimension(3, fSize, nTime), intent(out) :: fDot
      complex(kind=realType), dimension(1, fSize, nTime), intent(out) :: hfDot

      ! Working Variables
      integer(kind=intType) :: nn,sps, level
      integer(kind=intType) :: ierr, mm,i,j,k, l,  ii, jj, iRegion

      complex(kind=realType), dimension(costSize, size(famLists,1)) :: funcs


      ! Input Arguments for master:
      complex(kind=realType), dimension(costSize, size(famLists,1)) :: funcValues


          ! Working Variables
      complex(kind=realType), dimension(:, :, :), allocatable :: forces
      complex(kind=realType), dimension(:, :, :), allocatable :: heatfluxes
      complex(kind=realType) :: h ! step size for Finite Difference

      h = cmplx(0, 1e-40)

      fSize = size(fDot, 2)
      allocate(forces(3, fSize, nTimeIntervalsSpectral))

      fSize = size(hfDot, 2)
      allocate(heatfluxes(1, fSize, nTimeIntervalsSpectral))



      ! Need to trick the residual evalution to use coupled (mean flow and
      ! turbulent) together.
      level = 1
      currentLevel = level
      groundLevel = level

      ! Allocate the memory we need for derivatives if not done so
      ! already. Note this isn't deallocated until the adflow is
      ! destroyed.
      if (.not. derivVarsAllocated) then
         call allocDerivativeValues(level)
      end if

      ! Zero all AD seesd.
      do nn=1,nDom
         do sps=1,nTimeIntervalsSpectral
            call zeroADSeeds(nn,level, sps)
         end do
      end do

      ! Set the extra seeds now do the extra ones. Note that we are assuming the
      ! machNumber used for the coefficients follows the Mach number,
      ! not the grid mach number.

      alpha = alpha + h*extraDot(iAlpha)
      beta = beta + h*extraDot(iBeta)
      mach = mach + h*extraDot(iMach)
      machCoef = machCoef + h*extraDot(iMach)
      machGrid = machGrid + h*extraDot(iMachGrid)
      PinfDim = PinfDim + h*extraDot(iPressure)
      rhoinfDim = rhoinfDim + h*extraDot(iDensity)
      tinfdim = tinfdim + h*extraDot(iTemperature)
      pointref(1) = pointref(1) + h*extraDot(iPointRefX)
      pointref(2) = pointref(2) + h*extraDot(iPointRefY)
      pointref(3) = pointref(3) + h*extraDot(iPointRefZ)
      rgasdim = rgasdim + h*zero


      !  --------------------- apply the perturbations ----------------------------
      ! Set the provided w and x seeds:
      ii = 0
      jj = 0
      domainLoop1: do nn=1,nDom
         spectalLoop1: do sps=1,nTimeIntervalsSpectral
            call setPointers(nn, 1, sps)
            do k=1, kl
               do j=1,jl
                  do i=1,il
                     do l=1,3
                        ii = ii + 1
                        x(i, j, k, l) = x(i, j, k, l) +  xvdot(ii)*h
                     end do
                  end do
               end do
            end do
            do k=2, kl
               do j=2,jl
                  do i=2,il
                     do l = 1, nw
                        jj = jj + 1
                        w(i, j, k, l) = w(i, j, k, l) +  wDot(jj)*h
                     end do
                  end do
               end do
            end do
         end do spectalLoop1
      end do domainLoop1



      ! ----------------------------- Run Master ---------------------------------
      ! Run the super-dee-duper master rotuine
      if (bcVarsEmpty) then
         call master(useSpatial, famLists, funcValues, forces, heatfluxes)
      else
         call master(useSpatial, famLists, funcValues, forces, heatfluxes, &
                     bcDataNames, bcDataValues, bcDataFamLists)
      end if


       ! Copy out the residual derivative into the provided dwDot and remove the
       ! perturbation
      ii = 0
      jj = 0
      do nn=1, nDom
         do sps=1,nTimeIntervalsSpectral
            call setPointers_d(nn, 1, sps)
            do k=1, kl
               do j=1,jl
                  do i=1,il
                     do l=1,3
                        ii = ii + 1
                        x(i, j, k, l) = x(i, j, k, l) -  xvdot(ii)* h
                     end do
                  end do
               end do
            end do
            do k=2, kl
               do j=2, jl
                  do i=2, il
                     do l=1, nw
                        jj = jj + 1
                        w(i, j, k, l) = w(i, j, k, l) -  wDot(jj)*h
                        dwd(i,j,k,l) = aimag(dw(i,j,k,l))/aimag(h)
                        dwdot(jj) = dwd(i,j,k,l) ! copy values to output
                     end do
                  end do
               end do
            end do

         end do
      end do

      fDot = aimag(forces)/aimag(h)
      hfDot = aimag(heatfluxes)/aimag(h)
      write(*,*) 'fDot', minval(real(fDot)), maxval(real(fDot))
      write(*,*) 'hfDot', minval(real(hfDot)), maxval(real(hfDot))

   end subroutine computeMatrixFreeProductFwdCS
 ! this isn't compliling anymore 

   ! subroutine printCSSeeds(nn, level, sps)

   !    use constants
   !    use block, only : flowDoms
   !    use blockPointers
   !    use inputTimeSpectral
   !    use flowVarRefState
   !    use inputPhysics
   !    use BCPointers

   !    use communication
   !    use oversetData, only : oversetPresent
   !    use cgnsGrid, only : cgnsDoms, cgnsNDom
   !    use actuatorRegionData, only : nActuatorRegions, actuatorRegions
   !    implicit none

   !    ! Input parameters
   !    integer(kind=intType) :: nn, level, sps

   !    ! Working parameters
   !    integer(kind=intType) :: mm, i, iDom
   !    integer(kind=intType) :: iBoco, iData, iDirichlet
   !    real(kind=realType) :: h=1e-40

   !    ! call setPointers(nn, level, sps)

   !    write(*,*) 'd2wall ', minval(imag(flowDoms(nn, level, sps)%d2wall)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%d2wall)/h)
   !    write(*,*) 'x ', minval(imag(flowDoms(nn, level, sps)%x)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%x)/h)
   !    write(*,*) 'si ', minval(imag(flowDoms(nn, level, sps)%si)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%si)/h)
   !    write(*,*) 'sj ', minval(imag(flowDoms(nn, level, sps)%sj)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%sj)/h)
   !    write(*,*) 'sk ', minval(imag(flowDoms(nn, level, sps)%sk)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%sk)/h)
   !    write(*,*) 'vol ', minval(imag(flowDoms(nn, level, sps)%vol)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%vol)/h)

   !    write(*,*) 's ', minval(imag(flowDoms(nn, level, sps)%s)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%s)/h)
   !    write(*,*) 'sFaceI ', minval(imag(flowDoms(nn, level, sps)%sFaceI)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%sFaceI)/h)
   !    write(*,*) 'sFaceJ ', minval(imag(flowDoms(nn, level, sps)%sFaceJ)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%sFaceJ)/h)
   !    write(*,*) 'sFaceK ', minval(imag(flowDoms(nn, level, sps)%sFaceK)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%sFaceK)/h)

   !    write(*,*) 'w ', minval(imag(flowDoms(nn, level, sps)%w)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%w)/h)
   !    write(*,*) 'dw ', minval(imag(flowDoms(nn, level, sps)%dw)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%dw)/h)
   !    write(*,*) 'fw ', minval(imag(flowDoms(nn, level, sps)%fw)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%fw)/h)
   !    write(*,*) 'scratch ', minval(imag(flowDoms(nn, level, sps)%scratch)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%scratch)/h)

   !    write(*,*) 'p ', minval(imag(flowDoms(nn, level, sps)%p)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%p)/h)
   !    write(*,*) 'gamma ', minval(imag(flowDoms(nn, level, sps)%gamma)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%gamma)/h)
   !    write(*,*) 'aa ', minval(imag(flowDoms(nn, level, sps)%aa)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%aa)/h)

   !    write(*,*) 'rlv ', minval(imag(flowDoms(nn, level, sps)%rlv)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%rlv)/h)
   !    write(*,*) 'rev ', minval(imag(flowDoms(nn, level, sps)%rev)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%rev)/h)

   !    write(*,*) 'radI ', minval(imag(flowDoms(nn, level, sps)%radI)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%radI)/h)
   !    write(*,*) 'radJ ', minval(imag(flowDoms(nn, level, sps)%radJ)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%radJ)/h)
   !    write(*,*) 'radK ', minval(imag(flowDoms(nn, level, sps)%radK)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%radK)/h)

   !    write(*,*) 'ux ', minval(imag(flowDoms(nn, level, sps)%ux)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%ux)/h)
   !    write(*,*) 'uy ', minval(imag(flowDoms(nn, level, sps)%uy)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%uy)/h)
   !    write(*,*) 'uz ', minval(imag(flowDoms(nn, level, sps)%uz)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%uz)/h)
   !    write(*,*) 'vx ', minval(imag(flowDoms(nn, level, sps)%vx)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%vx)/h)
   !    write(*,*) 'vy ', minval(imag(flowDoms(nn, level, sps)%vy)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%vy)/h)
   !    write(*,*) 'vz ', minval(imag(flowDoms(nn, level, sps)%vz)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%vz)/h)
   !    write(*,*) 'wx ', minval(imag(flowDoms(nn, level, sps)%wx)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%wx)/h)
   !    write(*,*) 'wy ', minval(imag(flowDoms(nn, level, sps)%wy)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%wy)/h)
   !    write(*,*) 'wz ', minval(imag(flowDoms(nn, level, sps)%wz)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%wz)/h)
   !    write(*,*) 'qx ', minval(imag(flowDoms(nn, level, sps)%qx)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%qx)/h)
   !    write(*,*) 'qy ', minval(imag(flowDoms(nn, level, sps)%qy)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%qy)/h)
   !    write(*,*) 'qz ', minval(imag(flowDoms(nn, level, sps)%qz)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%qz)/h)

   !    write(*,*) 'bmti1 ',minval(imag(flowDoms(nn, level, sps)%bmti1)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bmti1)/h)
   !    write(*,*) 'bmti2 ',minval(imag(flowDoms(nn, level, sps)%bmti2)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bmti2)/h)
   !    write(*,*) 'bmtj1 ',minval(imag(flowDoms(nn, level, sps)%bmtj1)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bmtj1)/h)
   !    write(*,*) 'bmtj2 ',minval(imag(flowDoms(nn, level, sps)%bmtj2)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bmtj2)/h)
   !    write(*,*) 'bmtk1 ',minval(imag(flowDoms(nn, level, sps)%bmtk1)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bmtk1)/h)
   !    write(*,*) 'bmtk2 ',minval(imag(flowDoms(nn, level, sps)%bmtk2)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bmtk2)/h)
   !    write(*,*) 'bvti1 ',minval(imag(flowDoms(nn, level, sps)%bvti1)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bvti1)/h)
   !    write(*,*) 'bvti2 ',minval(imag(flowDoms(nn, level, sps)%bvti2)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bvti2)/h)
   !    write(*,*) 'bvtj1 ',minval(imag(flowDoms(nn, level, sps)%bvtj1)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bvtj1)/h)
   !    write(*,*) 'bvtj2 ',minval(imag(flowDoms(nn, level, sps)%bvtj2)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bvtj2)/h)
   !    write(*,*) 'bvtk1 ',minval(imag(flowDoms(nn, level, sps)%bvtk1)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bvtk1)/h)
   !    write(*,*) 'bvtk2 ',minval(imag(flowDoms(nn, level, sps)%bvtk2)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bvtk2)/h)


   !    bocoLoop: do mm=1, flowDoms(nn, level, sps)%nBocos

   !       select case (flowDoms(nn, level, sps)%BCType(mm))

   !       case (NSWallAdiabatic)
   !          write(*,*) 'mm', mm, 'BCData(mm)%norm',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%norm)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%norm)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Fv ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fv)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fv)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Fp ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fp)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fp)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Tv ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tv)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tv)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Tp ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tp)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tp)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%area ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%area)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%area)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%uSlip ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%uSlip)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%uSlip)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%cellHeatFlux ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%cellHeatFlux)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%cellHeatFlux)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%nodeHeatFlux ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%nodeHeatFlux)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%nodeHeatFlux)/h)


   !          !=======================================================

   !       case (NSWallIsothermal)
   !          write(*,*) 'mm', mm, 'BCData(mm)%norm',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%norm)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%norm)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Fv ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fv)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fv)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Fp ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fp)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fp)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Tv ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tv)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tv)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Tp ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tp)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tp)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%area ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%area)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%area)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%uSlip ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%uSlip)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%uSlip)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%TNS_Wall ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%TNS_Wall)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%TNS_Wall)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%cellHeatFlux ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%cellHeatFlux)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%cellHeatFlux)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%nodeHeatFlux ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%nodeHeatFlux)/h), &
   !                      maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%nodeHeatFlux)/h)

   !          !=======================================================

   !       case (EulerWall)
   !          write(*,*) 'mm', mm, 'BCData(mm)%norm',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%norm)/h), &
   !                                  maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%norm)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Fv ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fv)/h), &
   !                                  maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fv)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Fp ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fp)/h), &
   !                                  maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Fp)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Tv ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tv)/h), &
   !                                  maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tv)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%Tp ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tp)/h), &
   !                                  maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%Tp)/h)
   !          write(*,*) 'mm', mm, 'bcData(mm)%area ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%area)/h), &
   !                                  maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%area)/h)



   !    !=======================================================

   !       case (farField)

   !          ! Just allocate the memory for the normal mesh
   !          ! velocity.

   !          write(*,*) 'mm', mm, 'bcData(mm)%rface ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%rface)/h), &
   !                                                   maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%rface)/h)

   !          !=======================================================

   !       case (symm, symmPolar)

   !          ! Allocate for symm as well. This is not necessary
   !          ! but we need it for the reverse AD.

   !          ! Modified by HDN

   !          write(*,*) 'mm', mm, 'bcData(mm)%rface ',minval(imag(flowDoms(nn, level, sps)%bcData(mm)%rface)/h), &
   !                                                   maxval(imag(flowDoms(nn, level, sps)%bcData(mm)%rface)/h)

   !          !=======================================================

   !       case (SupersonicInflow, DomainInterfaceAll)

   !          ! Supersonic inflow or a domain interface with
   !          ! all the data prescribed. Allocate the memory for
   !          ! the entire state vector to be prescribed.


   !          write(*,*) 'mm', mm, 'BCData(mm)%ps ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%ps)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%ps)/h)
   !          if(nt2 >= nt1) then
   !             write(*,*) 'mm', mm, 'BCData(mm)%turbInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%turbInlet)/h), &
   !                               maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%turbInlet)/h)
   !          endif

   !       case (SupersonicOutflow)
   !          ! No state is needed for this boco

   !       case (SubsonicInflow)

   !          ! Subsonic inflow. Allocate the memory for the
   !          ! variables needed. Note the there are two ways to
   !          ! specify boundary conditions for a subsonic inflow.

   !          write(*,*) 'mm', mm, 'BCData(mm)%ptInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%ptInlet)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%ptInlet)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%htInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%htInlet)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%htInlet)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%ttInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%ttInlet)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%ttInlet)/h)

   !          if(nt2 >= nt1) then
   !             write(*,*) 'mm', mm, 'BCData(mm)%turbInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%turbInlet)/h), &
   !                               maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%turbInlet)/h)

   !          endif

   !          !=======================================================

   !       case (SubsonicOutflow, MassBleedOutflow, &
   !             DomainInterfaceP)

   !          ! Subsonic outflow, outflow mass bleed or domain
   !          ! interface with prescribed pressure. Allocate the
   !          ! memory for the static pressure.

   !          write(*,*) 'mm', mm, 'BCData(mm)%ps ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%ps)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%ps)/h)


   !       case (DomainInterfaceRhoUVW)

   !          ! Domain interface with prescribed density and
   !          ! velocities, i.e. mass flow is prescribed. Allocate
   !          ! the memory for the variables needed.

   !          if(nt2 >= nt1) then
   !             write(*,*) 'mm', mm, 'BCData(mm)%turbInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%turbInlet)/h), &
   !                               maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%turbInlet)/h)

   !          endif


   !          !=======================================================

   !       case (DomainInterfaceTotal)

   !          ! Domain interface with prescribed total conditions.
   !          ! Allocate the memory for the variables needed.

   !          write(*,*) 'mm', mm, 'BCData(mm)%ptInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%ptInlet)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%ptInlet)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%htInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%htInlet)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%htInlet)/h)
   !          write(*,*) 'mm', mm, 'BCData(mm)%ttInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%ttInlet)/h), &
   !                            maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%ttInlet)/h)

   !          if(nt2 >= nt1) then
   !             write(*,*) 'mm', mm, 'BCData(mm)%turbInlet ',minval(imag(flowDoms(nn, level, sps)%BCData(mm)%turbInlet)/h), &
   !                               maxval(imag(flowDoms(nn, level, sps)%BCData(mm)%turbInlet)/h)

   !          endif

   !          !=======================================================

   !       case (domainInterfaceRho)

   !       end select

   !    enddo bocoLoop



   !    viscbocoLoop: do mm=1,flowDoms(nn, level, sps)%nViscBocos
   !       write(*,*) 'mm', mm, 'viscSubface(mm)%tau ',minval(imag(flowDoms(nn, level, sps)%viscSubface(mm)%tau)/h), &
   !                               maxval(imag(flowDoms(nn, level, sps)%viscSubface(mm)%tau)/h)
   !       write(*,*) 'mm', mm, 'viscSubface(mm)%q ',minval(imag(flowDoms(nn, level, sps)%viscSubface(mm)%q)/h), &
   !                               maxval(imag(flowDoms(nn, level, sps)%viscSubface(mm)%q)/h)
   !    end do viscbocoLoop

   !    ! For overset, the weights may be active in the comm structure. We
   !    ! need to zero them before we can accumulate.
   !    if (oversetPresent) then
   !       ! Pointers to the overset comms to make it easier to read
   !       sends: do i=1,commPatternOverset(level, sps)%nProcSend
   !          write(*,*) 'commPatternOverset(level, sps)%sendList(i)%interpd ',&
   !                      minval(imag(commPatternOverset(level, sps)%sendList(i)%interp)/h), &
   !                      maxval(imag(commPatternOverset(level, sps)%sendList(i)%interp)/h)
   !       end do sends
   !       write(*,*) 'internalOverset(level, sps)%donorInterpd ',minval(imag(internalOverset(level, sps)%donorInterp)/h), &
   !                               maxval(imag(internalOverset(level, sps)%donorInterp)/h)
   !    end if

   !    write(*,*) 'alphad ', imag(alpha)/h
   !    write(*,*) 'betad ', imag(beta)/h
   !    write(*,*) 'machd ', imag(mach)/h
   !    write(*,*) 'machGridd ', imag(machGrid)/h
   !    write(*,*) 'machCoefd ', imag(machCoef)/h
   !    write(*,*) 'pinfdimd ', imag(pinfdim)/h
   !    write(*,*) 'tinfdimd ', imag(tinfdim)/h
   !    write(*,*) 'rhoinfdimd ', imag(rhoinfdim)/h
   !    write(*,*) 'rgasdimd ', imag(rgasdim)/h
   !    write(*,*) 'pointrefd ', imag(pointref)/h
   !    write(*,*) 'prefd ', imag(pref)/h
   !    write(*,*) 'rhoRefd ', imag(rhoRef)/h
   !    write(*,*) 'Trefd ', imag(Tref)/h
   !    write(*,*) 'murefd ', imag(muref)/h
   !    write(*,*) 'urefd ', imag(uref)/h
   !    write(*,*) 'hrefd ', imag(href)/h
   !    write(*,*) 'timerefd ', imag(timeref)/h
   !    write(*,*) 'pinfd ', imag(pinf)/h
   !    write(*,*) 'pinfCorrd ', imag(pinfCorr)/h
   !    write(*,*) 'rhoinfd ', imag(rhoinf)/h
   !    write(*,*) 'uinfd ', imag(uinf)/h
   !    write(*,*) 'rgasd ', imag(rgas)/h
   !    write(*,*) 'muinfd ', imag(muinf)/h
   !    write(*,*) 'gammainfd ', imag(gammainf)/h
   !    write(*,*) 'winfd ', imag(winf)/h
   !    write(*,*) 'veldirfreestreamd ', imag(veldirfreestream)/h
   !    write(*,*) 'liftdirectiond ', imag(liftdirection)/h
   !    write(*,*) 'dragdirectiond ', imag(dragdirection)/h

   !    ! Zero all the reverse seeds in the dirichlet input arrays
   !    write(*,*) 'iDom, iBoco, iData, iDirichlet'
   !    do iDom=1, cgnsNDom
   !       do iBoco=1, cgnsDoms(iDom)%nBocos
   !          if (associated(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet)) then
   !             do iData=1, size(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet)
   !                if (associated(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)) then
   !                   do iDirichlet = 1, size(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)
   !                      write(*,*) iDom, iBoco, iData, iDirichlet, 'dataArr(:) '&
   !                      ,imag(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays(iDirichlet)%dataArr(:))/h
   !                   end do
   !                end if
   !             end do
   !          end if
   !       end do
   !    end do

   !    ! And the reverse seeds in the actuator zones
   !    do i=1, nActuatorRegions
   !       write(*,*) 'actuatorRegionsd(i)%F ',imag(actuatorRegions(i)%F)/h
   !       write(*,*) 'actuatorRegionsd(i)%T ',imag(actuatorRegions(i)%T)/h
   !    end do

   ! end subroutine printCSSeeds

#endif

end module adjointDebug
