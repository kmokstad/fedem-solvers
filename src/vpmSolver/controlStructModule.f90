!! $Id$
!!==============================================================================
!!
!!    F E D E M    T E C H N O L O G Y    A S
!!
!!    Copyright (C)
!!    1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008
!!    FEDEM Technology AS
!!    all rights reserved
!!
!!    This is UNPUBLISHED PROPRIETARY SOURCE CODE of FEDEM Technology AS;
!!    the contents of this file may not be disclosed to third parties,
!!    copied or duplicated in any form, in whole or in part, without
!!    the prior written permission of FEDEM Technology AS.
!!
!!==============================================================================

module ControlStructModule

  use MechanismTypeModule
  use ControlTypeModule

  implicit none

  type StructSensorGradType
     !!TODO,bh: Check if we need local dofstart for each node within sensorGrad_wrt_disp
     integer                    :: iCin            !! Which cIn this sensor is connected to
     !                                             !! j_index i.e. second (colum) index in Jacobi_CinToCout(:,:)
     !                                             !! i_index i.e. first (row) index in Jacobi_SensorToCIn(:,:)
     integer                    :: whichVreg       !! which VREG this is connected to (not necessarily equal to jIndex)
     integer                    :: lNode(2)        !! Local node number for the element matrecies
     integer                    :: dofStart(2)     !! Local dof start for this node
     integer                    :: nDofs(2)        !! number of dofs for this node
     type(CtrlPrm), pointer :: pCtrlPrm        !! pointer to the control parameter which has
     !                                             !! structural input
     type(SensorType),  pointer :: pSensor         !! The actual sensor pointer
     type(TriadType),   pointer :: pTriad1         !! non null if sensor on a triade disp, vel, or accel
     type(TriadType),   pointer :: pTriad2         !! non null if relative sensor on triads
     type(MasterSlaveJointType), pointer :: pJoint !! non null if sensor on a joint dof

     real(dp)                   :: sensorGrad_wrt_disp(12)!! 12 should be enough space
     !TODO, Missing which cIn this is connected to

     !! ?? How to handle pos, vel, or accel ???
  end type StructSensorGradType


  type ForceControlGradType
     integer                  :: iCout           !! Which cOut this force is connected to
     !                                           !! i_index (first)  index in Jacobi_CinToCout(:,:)
     !                                           !! j_index (second) index in Jacobi_ForceToCout(:,:)
     integer                  :: whichVreg       !! which VREG this is connected to (not necessarily equal to jIndex)
     integer                  :: lNode           !! Local node number for the element matrecies
     integer                  :: dofStart        !! Local dof start for this node
     integer                  :: nDofs           !! number of dofs for this node
     type(ForceType), pointer :: pForce          !! pointer to force with control out as input
     real(dp)                 :: forceGrad_wrt_cOut(6) !! Always dim 6
     !TODO: missing which cOut this is connected to
  end type ForceControlGradType

  type ControlStructType

     type(idType) :: id

     integer :: ctrlSysEigFlag

     integer :: samElNum
     integer :: nNods
     integer :: nDOFs
     integer, pointer :: samMNPC(:)
     integer, pointer :: local_MADOF(:)

     type(MechanismType),  pointer :: pMech
     type(ControlType),    pointer :: pControl

     !! CONTROL PERTURBATION

     integer          :: numVregIn       !! number of inputs to perturb, =size(whichVregIn)
     integer, pointer :: whichVregIn(:)  !! the vreg indexes that have a structural sensor

     integer          :: numVregOut      !! number of outputs to read perturbation from, =size(whichVregOut)
     integer, pointer :: whichVregOut(:) !! the vreg indexes which have a structural load vector

     !! SENSOR SIDE

     integer                             :: numStructCtrlParams
     type(StructSensorGradType), pointer :: structToControlSensors(:)

     !! FORCE SIDE

     integer                             :: numControlForces
     type(ForceControlGradType), pointer :: controlForces(:)


     !! The tangent matrecies
     !! These must hav a lookup table to internal dofs from both input side (sensor) and output
     !! (force) side where the input equal to output is properly handled (ie has the same dof).
     !! The dimension of the matrecies are "Number of unique input and output dofs"

     !! system equation: M*(d2x/dt2) + C*(dx/dt) + K*X + Q*int(x)dt = F
     real(dp), pointer :: ctrlProps(:,:,:) !! table of controller properties
     !                                     !!(no. of outputs from controller,
     !                                     !! no. of controller properties,
     !                                     !! no. of inputs to controller)
     !                                     !! ctrlProps(:,1,:) = Q
     !                                     !! ctrlProps(:,2,:) = K
     !                                     !! ctrlProps(:,3,:) = C
     !                                     !! ctrlProps(:,4,:) = M
     real(dp), pointer :: massMat(:,:)     !! mass matrix M
     real(dp), pointer :: dampMat(:,:)     !! damping matrix C
     real(dp), pointer :: stiffMat(:,:)    !! stiffness matrix K
     real(dp), pointer :: SSEEMat(:,:)     !! steady-state error elimination matrix Q

     real(dp), pointer :: Grad_CinWrtSensor(:,:) !! Dim (nCin,nDofs) Consists of rows of sensorGrad_wrt_disp
     real(dp), pointer :: Grad_CoutWrtCin(:,:)   !! Dim (nCout,nCin) What about mass, stiffness and damping

     !!
     real(dp), pointer :: Grad_ForceWrtCout(:,:) !! Dim (nDofs,nCout) Consists of columns of forceGrad_wrt_cOut

     integer :: SSEEMatIsNonSymmetric
     integer :: stiffMatIsNonSymmetric
     integer :: dampMatIsNonSymmetric
     integer :: massMatIsNonSymmetric
     integer :: anyMatIsNonSymmetric

  end type ControlStructType



contains

  subroutine NullifyControlStruct(ctrlStruct)
    use IdTypeModule, only : nullifyId

    type(controlStructType), intent(out) :: ctrlStruct

    call NullifyId(ctrlStruct%id)
    ctrlStruct%samElNum = 0
    ctrlStruct%nDOFs = 0
    nullify(ctrlStruct%pMech)
    nullify(ctrlStruct%pControl)

    nullify(ctrlStruct%massMat)
    nullify(ctrlStruct%dampMat)
    nullify(ctrlStruct%stiffMat)

  end subroutine NullifyControlStruct


  subroutine InitiateControlStruct(pCS)

    !!==============================================================================
    !! Initialize the control struct data type
    !!
    !!==============================================================================

    use ForceTypeModule
    use SensorTypeModule
    use ForceRoutinesModule
    use ControlRoutinesModule

    type(ControlStructType), intent(inout) :: pCS

    type(ForceType),            pointer :: pForce
    type(CtrlPrm),              pointer :: pCtrlPrm
    type(StructSensorGradType), pointer :: pS
    type(ForceControlGradType), pointer :: pF

    integer, allocatable :: lNode_from_SAM_node(:)
    integer, allocatable :: iCin_from_allVreg(:)
    integer, allocatable :: iCout_from_allVreg(:)
    integer, pointer     :: tmpIdx(:)

    integer :: i, n, nNodes, nDofs, ierr, iSize, iSAM, nElNodes, iStart, numVreg
    integer :: iCin, iCout, nCin, nCout, whichVreg

    ierr   = 0
    nNodes = 0
    nDofs  = 0

    !! Initialize the input side

    iSize = (size(pCS%pMech%triads) + size(pCS%pMech%sups) + size(pCS%pMech%joints)) * 10
    allocate( lNode_from_SAM_node(iSize))
    lNode_from_SAM_node = 0

    numVreg = size(pCS%pControl%vreg)

    allocate( iCin_from_AllVreg(numVreg) )
    iCin_from_AllVreg  = 0


    !! SENSOR side initialization
    !! Find the control input parameters that have structural sensors (triad dofs and joint dofs)

    nNodes = 0
    nDofs  = 0

    call getCtrlParamsWithStructSensors(pCS%pControl%input,pCS%numStructCtrlParams, tmpIdx)
    if (pCS%numStructCtrlParams <= 0 ) return

    allocate(pCS%structToControlSensors(pCS%numStructCtrlParams))

    do i = 1, pCS%numStructCtrlParams
       pS          => pCS%structToControlSensors(i)
       pS%pCtrlPrm => pCS%pControl%input(tmpIdx(i))

       pS%iCin = 0  !! i.e. not set yet
       pS%whichVreg = pS%pCtrlPrm%var
       iCin_from_AllVreg(pS%whichVreg) = 1 !! flag that this vreg is part of the cIn side

          !!TODO,bh: check this
          pS%pSensor => pS%pCtrlPrm%engine%args(1)%p


       if      (pS%pSensor%type == TRIAD_p) then
          pS%pTriad1 => pCS%pMech%triads(pS%pSensor%index(1))
          nullify(pS%pTriad2)
          nullify(pS%pJoint)
          lNode_from_SAM_node(pS%pTriad1%samNodNum) = 1

       else if (pS%pSensor%type == RELATIVE_TRIAD_p) then

          pS%pTriad1 => pCS%pMech%triads(pS%pSensor%index(1))
          pS%pTriad2 => pCS%pMech%triads(pS%pSensor%index(2))
          nullify(pS%pJoint)
          lNode_from_SAM_node(pS%pTriad1%samNodNum) = 1
          lNode_from_SAM_node(pS%pTriad2%samNodNum) = 1

       else if (pS%pSensor%type == JOINT_VARIABLE_p) then

          nullify(pS%pTriad1)
          nullify(pS%pTriad2)
          pS%pJoint => pCS%pMech%joints(pS%pSensor%index(1))
          lNode_from_SAM_node(pS%pJoint%samNodNum) = 1

       else
          !! Error, actually
       end if

    end do

    !! Count number of vreg with structural input
    pCS%numVregIn = 0
    do i = 1, numVreg
       if (iCin_from_AllVreg(i) > 0) then
          pCS%numVregIn = pCS%numVregIn + 1
          iCin_from_AllVreg(i) = pCS%numVregIn
       end if
    end do

    !! Set which vreg in (with structural input) are to be perturbed
    allocate(pCS%whichVregIn(pCS%numVregIn))

    !! Set which iCin this sensor is connected to
    do i = 1, pCS%numStructCtrlParams
       whichVreg                          = pCS%structToControlSensors(i)%whichVreg
       iCin                               = iCin_from_AllVreg(whichVreg)
       pCS%structToControlSensors(i)%iCin = iCin
       pCS%whichVregIn(iCin)              = whichVreg
    end do


    deallocate(iCin_from_AllVreg) !! done with this scratch table
    deallocate(tmpIdx)            !! done with this scratch

    !! FORCES side initialization
    !! Find the forces whose sensor is a controller variable

    allocate( iCout_from_AllVreg(numVreg) )
    iCout_from_AllVreg  = 0

    call getCtrlOutForces (pCS%pMech%forces, pCS%numControlForces, tmpIdx)
    if (pCS%numControlForces <= 0) return

    allocate(pCS%controlForces(pCS%numControlForces))

    do i = 1, pCS%numControlForces
       pF         => pCS%controlForces(i)
       pF%pForce  => pCS%pMech%forces(tmpIdx(i))

       pF%iCout = 0 !! i.e. not yet set
       pF%whichVreg = pF%pForce%engine%args(1)%p%dof
       iCout_from_AllVreg(pF%whichVreg) = 1  !! flag this as one to read from

       if      (associated(pF%pForce%triad)) then

          lNode_from_SAM_node(pF%pForce%triad%samNodNum) = 1

       else if (associated(pF%pForce%joint)) then

          lNode_from_SAM_node(pF%pForce%joint%samNodNum) = 1

       else
          !! Error, actually
       end if

    end do

    !! Count number of vreg which are to be read
    pCS%numVregOut = 0
    do i = 1, numVreg
       if (iCout_from_AllVreg(i) > 0) then
          pCS%numVregOut = pCS%numVregOut + 1
          iCout_from_AllVreg(i) = pCS%numVregOut
       end if
    end do

    !! Set which vreg in are to be read
    allocate(pCS%whichVregOut(pCS%numVregOut))

    do i = 1, pCS%numControlForces
       whichVreg                  = pCS%controlForces(i)%whichVreg
       iCout                      = iCout_from_AllVreg(whichVreg)
       pCS%controlForces(i)%iCout = iCout
       pCS%whichVregOut(iCout)    = whichVreg
    end do

    deallocate(iCout_from_AllVreg) !! done with this scratch table
    deallocate(tmpIdx)

    !! Count number of element nodes and set the local element node number

    nElNodes = 0
    do iSam = 1, size(lNode_from_SAM_node)   !! Loop over all sam node nums
       if ( lNode_from_SAM_node(iSam) > 0 ) then
          nELNodes = nElNodes + 1
          lNode_from_SAM_node(iSam) = nElNodes
       end if
    end do

    !! Allocate and initialize MNPC

    pCS%nNods = nElNodes

    allocate(pCS%samMNPC(nElNodes))
    allocate(pCS%local_MADOF(nElNodes+1))
    pCS%local_MADOF = 0

    do iSam = 1, size(lNode_from_SAM_node)   !! Loop over all sam node nums
       if ( lNode_from_SAM_node(iSam) > 0 ) then
          pCS%samMNPC( lNode_from_SAM_node(iSam) ) = iSam
       end if
    end do


    !! Set the local node association and dofStart for all the forces and sensors


    n = size(pCS%structToControlSensors)
    do i = 1,n
       pS => pCS%structToControlSensors(i)
       pS%lNode    = 0
       pS%dofStart = 0

       if ( associated(pS%pTriad1)) then
          pS%lNode(1)             = lNode_from_SAM_node(pS%pTriad1%samNodNum)
          pCS%local_MADOF(pS%lNode(1))= pS%pTriad1%nDOFs
          pS%nDOFs(1)             = pS%pTriad1%nDOFs
       end if

       if ( associated(pS%pTriad2)) then
          pS%lNode(2)             = lNode_from_SAM_node(pS%pTriad2%samNodNum)
          pCS%local_MADOF(pS%lNode(2))= pS%pTriad2%nDOFs
          pS%nDOFs(2)             = pS%pTriad2%nDOFs
       end if

       if ( associated(pS%pJoint))  then
          pS%lNode(1)             = lNode_from_SAM_node(pS%pJoint%samNodNum)
          pCS%local_MADOF(pS%lNode(2))= pS%pJoint%nJointDOFs
          pS%nDOFs(1)             = pS%pJoint%nJointDOFs
       end if

    end do

    n = size(pCS%controlForces)
    do i = 1,n
       pF => pCS%controlForces(i)
       pF%lNode    = 0
       pF%dofStart = 0

       if      (associated(pF%pForce%triad)) then
          pF%lNode             = lNode_from_SAM_node(pF%pForce%triad%samNodNum)
          pCS%local_MADOF(pF%lNode)= pF%pForce%triad%nDOFs
          pF%nDOFs             = pF%pForce%triad%nDOFs

       else if (associated(pF%pForce%joint)) then
          pF%lNode             = lNode_from_SAM_node(pF%pForce%joint%samNodNum)
          pCS%local_MADOF(pF%lNode)= pF%pForce%joint%nJointDOFs
          pF%nDOFs             = pF%pForce%joint%nJointDOFs

       end if

    end do

    !! Now accumulate the dofStart

    iStart = 1
    do i = 1, size(pCS%local_MADOF)
       nDofs = pCS%local_MADOF(i)
       pCS%local_MADOF(i) = iStart
       iStart = iStart + nDofs
    end do

    !! Number of total dofs for this element
    pCS%nDOFs = pCS%local_MADOF(pCS%nNods+1)-1

    nDofs = pCS%nDOFs
    nCin  = pCS%numVregIn
    nCout = pCS%numVregOut

    allocate(pCS%Grad_CinWrtSensor(pCS%numVregIn,pCS%nDofs))
    allocate(pCS%Grad_CoutWrtCin(pCS%numVregOut,pCS%numVregIn))
    allocate(pCS%Grad_ForceWrtCout(pCS%nDofs,pCS%numVregOut))
    allocate(pCS%ctrlProps(pCS%numVregOut,4,pCS%numVregIn))
    allocate(pCS%massMat(pCS%nDofs,pCS%nDofs))
    allocate(pCS%dampMat(pCS%nDofs,pCS%nDofs))
    allocate(pCS%stiffMat(pCS%nDofs,pCS%nDofs))
    allocate(pCS%SSEEMat(pCS%nDofs,pCS%nDofs))

    !! Also store the dofStart in the forces and sensors

    n = size(pCS%structToControlSensors)
    do i = 1,n
       pS => pCS%structToControlSensors(i)
       if (pS%lNode(1) > 0 ) pS%dofStart(1) = pCS%local_MADOF(pS%lNode(1))
       if (pS%lNode(2) > 0 ) pS%dofStart(2) = pCS%local_MADOF(pS%lNode(2))
    end do

    n = size(pCS%controlForces)
    do i = 1,n
       pF => pCS%controlForces(i)
       pF%dofStart = pCS%local_MADOF(pF%lNode)
    end do

    !! Clean up some scratch space
    deallocate( lNode_from_SAM_node )


  end subroutine InitiateControlStruct


  !!============================================================================
  !> @brief Finds all control input parameters coupled to structural DOFs.

  subroutine getCtrlParamsWithStructSensors (ctrlParams, numStructCtrlParams, &
       &                                     ctrlParamsWithStructSensors)

    use ControlTypeModule, only : CtrlPrm
    use ReportErrorModule, only : allocationError

    type(CtrlPrm), target, intent(in)  :: ctrlParams(:)
    integer,               intent(out) :: numStructCtrlParams
    integer,      pointer, intent(out) :: ctrlParamsWithStructSensors(:)

    !! Local variables
    integer :: i, iPrm

    !! --- Logic section ---

    numStructCtrlParams = 0
    do i = 1, size(ctrlParams)
       if (hasStructSensor(ctrlParams(i))) then
          numStructCtrlParams = numStructCtrlParams + 1
       end if
    end do

    allocate(ctrlParamsWithStructSensors(numStructCtrlParams),stat=iPrm)
    if (iPrm /= 0) then
       numStructCtrlParams = allocationError('getCtrlParamsWithStructSensors')
       return
    end if

    iPrm = 0
    do i = 1, size(ctrlParams)
       if (hasStructSensor(ctrlParams(i))) then
          iPrm = iPrm + 1
          ctrlParamsWithStructSensors(iPrm) = i
       end if
    end do

  contains

    !> @brief Checks if a control parameter is coupled to a structural DOF.
    logical function hasStructSensor (prm)
      use SensorTypeModule , only : SensorType, TRIAD_P, RELATIVE_TRIAD_p
      type(CtrlPrm), intent(in) :: prm
      integer :: iType
      type(SensorType), pointer :: pSensor => null()

      if (associated(prm%engine)) then
         !!TODO,bh: check this after the addition of the new argument array
         pSensor => prm%engine%args(1)%p
      end if

      if (associated(pSensor)) then
         iType = pSensor%type
         hasStructSensor = iType == TRIAD_P .or. iType == RELATIVE_TRIAD_p
      else
         hasStructSensor = .false.
      end if
    end function hasStructSensor

  end subroutine getCtrlParamsWithStructSensors


  !!============================================================================
  !> @brief Finds all control out forces.

  subroutine getCtrlOutForces (forces,numCtrlOutForces,ctrlOutForces)

    use ForceTypeModule  , only : ForceType
    use ReportErrorModule, only : allocationError

    type(ForceType), target, intent(in)  :: forces(:)
    integer,                 intent(out) :: numCtrlOutForces
    integer,        pointer, intent(out) :: ctrlOutForces(:)

    !! Local variables
    integer :: i, iCtrl

    !! --- Logic section ---

    numCtrlOutForces = 0
    do i = 1, size(forces)
       if (isCtrlOutForce(forces(i))) then
          numCtrlOutForces = numCtrlOutForces + 1
       end if
    end do
    if (numCtrlOutForces == 0) return

    allocate(ctrlOutForces(numCtrlOutForces),stat=iCtrl)
    if (iCtrl /= 0) then
       numCtrlOutForces = allocationError('getCtrlOutForces')
       return
    end if

    iCtrl = 0
    do i = 1, size(forces)
       if (isCtrlOutForce(forces(i))) then
          iCtrl = iCtrl + 1
          ctrlOutForces(iCtrl) = i
       end if
    end do

  contains

    !> @brief Checks if the source of a force is a control system variable.
    logical function isCtrlOutForce (force)
      use SensorTypeModule, only : CONTROL_p
      type(ForceType), intent(in) :: force
      isCtrlOutForce = .false.
      if (.not. associated(force%engine)) return
      if (.not. associated(force%engine%args)) return
      if (size(force%engine%args) < 1) return
      isCtrlOutForce = force%engine%args(1)%p%type == CONTROL_p
    end function isCtrlOutForce

  end subroutine getCtrlOutForces


  subroutine BuildStructControlJacobi(pCS,sys,sam)
    !!===========================================================================
    !!  Compute the gradient matrix of forces w.r.t. controller inputs
    !!===========================================================================

    use ForceRoutinesModule
    use EngineRoutinesModule
    use SystemTypeModule
    use SamModule

    implicit none

    type(ControlStructType), intent(inout) :: pCS
    type(SystemType), intent(inout) :: sys
    type(SamType), intent(in) :: sam

    integer  :: i, ii, j, iForce, iS, ierr
    integer  :: iNode, iCin, iCout, iStart, iEnd, nStep
    real(dp) :: eGrad
    real(dp), allocatable :: dScr(:,:)        !! Scratch for matrix multiply


    ierr = 0

    !! Compute the force gradients w.r.t control outputs

    pCS%Grad_ForceWrtCout = 0.0_dp

    do iForce = 1, size(pCS%controlForces)
       call calcExtForceGradient (pCS%controlForces(iForce)%pForce, &
            &                     pCS%controlForces(iForce)%forceGrad_wrt_cOut, ierr)

            !! Insert into Grad_ForceWrtCout
            iStart = pCS%controlForces(iForce)%dofStart
            iEnd   = iStart + pCS%controlForces(iForce)%nDofs - 1
            iCout  = pCS%controlForces(iForce)%iCout
            ii = 0
            do i = iStart, iEnd
               ii = ii + 1
               pCS%Grad_ForceWrtCout(i,iCout) = pCS%Grad_ForceWrtCout(i,iCout) &
               &                              + pCS%controlForces(iForce)%forceGrad_wrt_cOut(ii)
            end do
    end do

    !! Compute the control input gradients w.r.t displacement dofs

    pCS%Grad_CinWrtSensor = 0.0_dp

    do iS = 1, size(pCS%structToControlSensors)
       call SensorGradient_wrt_disp(pCS%structToControlSensors(iS)%pSensor, &
            &                       pCS%structToControlSensors(iS)%sensorGrad_wrt_disp, ierr)

       eGrad = EngineRate(pCS%structToControlSensors(iS)%pCtrlPrm%engine,ierr)

       do i = 1, size(pCS%structToControlSensors(iS)%sensorGrad_wrt_disp)
          pCS%structToControlSensors(iS)%sensorGrad_wrt_disp(i) = &
               & eGrad * pCS%structToControlSensors(iS)%sensorGrad_wrt_disp(i)
       end do

       do iNode = 1,2
          if ( pCS%structToControlSensors(iS)%lNode(iNode) <= 0 ) cycle

          !! Insert into Grad_CinWrtSensor
          iStart = pCS%structToControlSensors(iS)%dofStart(iNode)
          iEnd   = iStart + pCS%structToControlSensors(iS)%nDofs(iNode) - 1
          iCin   = pCS%structToControlSensors(iS)%iCin
          ii = (iNode-1)*6
          do i = iStart, iEnd
             ii = ii + 1
             pCS%Grad_CinWrtSensor(iCin,i) = pCS%Grad_CinWrtSensor(iCin,i) &  !!TODO,bh: Looka at the += form here
             &                             + pCS%structToControlSensors(iS)%sensorGrad_wrt_disp(ii)
          end do
       end do
    end do


    !! TODO, Magne: Add controller gradients here and build the full gradients
    pCS%ctrlProps = 0.0_dp

    select case (pCS%ctrlSysEigFlag)
    case (1) ! nPerturb = 3: P, I and D gains
       call EstimateControllerProperties01 (sys, pCS%pMech, pCS%pControl, sam%mpar, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               3, pCS%ctrlProps, ierr)

    case (2) ! nPerturb = 4
       call EstimateControllerProperties02 (sys, pCS%pMech, pCS%pControl, sam%mpar, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               4, pCS%ctrlProps, ierr)

    case (3) ! nPerturb = 6
       call EstimateControllerProperties03 (sys, pCS%pMech, pCS%pControl, sam%mpar, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               6, pCS%ctrlProps, ierr)

    case (4:8) ! nPerturb = 5
       nStep = pCS%ctrlSysEigFlag-3 ! nStep = 1...5
       call EstimateControllerProperties04 (sys, pCS%pControl, sam%mpar, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               nStep, pCS%ctrlProps, ierr)

    case (9:11) ! nPerturb = 5
       nStep = 10**(pCS%ctrlSysEigFlag-8) ! nStep = 10,100,1000
       call EstimateControllerProperties04 (sys, pCS%pControl, sam%mpar, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               nStep, pCS%ctrlProps, ierr)

    case (500) ! nPerturb = 1
       call EstimateControllerProperties500 (sys, pCS%pControl, sam%mpar, &
            &                                pCS%numVregIn, pCS%whichVregIn, &
            &                                pCS%numVregOut, pCS%whichVregOut, &
            &                                pCS%ctrlProps, ierr)

    case default ! Error
    end select


    allocate(dScr(pCS%nDofs,pCS%numVregIn))

    !! Steady-state error elimination matrix (Q)
    dScr         = 0.0_dp
    pCS%SSEEMat  = 0.0_dp
    dScr         = matmul(pCS%Grad_ForceWrtCout,pCS%ctrlProps(:,1,:))
    pCS%SSEEMat  = matmul(dScr,pCS%Grad_CinWrtSensor)
    pCS%SSEEMat  = pCS%SSEEMat*(-1.0_dp)

    !! Stiffness matrix (K)
    dScr         = 0.0_dp
    pCS%stiffMat = 0.0_dp
    dScr         = matmul(pCS%Grad_ForceWrtCout,pCS%ctrlProps(:,2,:))
    pCS%stiffMat = matmul(dScr,pCS%Grad_CinWrtSensor)
    pCS%stiffMat = pCS%stiffMat*(-1.0_dp)

    !! Damping matrix (C)
    dScr         = 0.0_dp
    pCS%dampMat  = 0.0_dp
    dScr         = matmul(pCS%Grad_ForceWrtCout,pCS%ctrlProps(:,3,:))
    pCS%dampMat  = matmul(dScr,pCS%Grad_CinWrtSensor)
    pCS%dampMat  = pCS%dampMat*(-1.0_dp)

    !! Mass matrix (M)
    dScr         = 0.0_dp
    pCS%massMat  = 0.0_dp
    dScr         = matmul(pCS%Grad_ForceWrtCout,pCS%ctrlProps(:,4,:))
    pCS%massMat  = matmul(dScr,pCS%Grad_CinWrtSensor)
    pCS%massMat  = pCS%massMat*(-1.0_dp)

    !! Check for symmetry
    pCS%SSEEMatIsNonSymmetric = 0
    pCS%stiffMatIsNonSymmetric = 0
    pCS%dampMatIsNonSymmetric = 0
    pCS%massMatIsNonSymmetric = 0
    pCS%anyMatIsNonSymmetric = 0

    outer: do i = 1,pCS%nDofs
       inner: do j = i,pCS%nDofs

          if (pCS%SSEEMatIsNonSymmetric == 0) then
             if (pCS%SSEEMat(i,j) /= pCS%SSEEMat(j,i)) then
             pCS%SSEEMatIsNonSymmetric = 1
             pCS%anyMatIsNonSymmetric = 1
             end if
          end if

          if (pCS%stiffMatIsNonSymmetric == 0) then
             if (pCS%stiffMat(i,j) /= pCS%stiffMat(j,i)) then
             pCS%stiffMatIsNonSymmetric = 1
             pCS%anyMatIsNonSymmetric = 1
             end if
          end if

          if (pCS%dampMatIsNonSymmetric == 0) then
             if (pCS%dampMat(i,j) /= pCS%dampMat(j,i)) then
             pCS%dampMatIsNonSymmetric = 1
             pCS%anyMatIsNonSymmetric = 1
             end if
          end if

          if (pCS%massMatIsNonSymmetric == 0) then
             if (pCS%massMat(i,j) /= pCS%massMat(j,i)) then
             pCS%massMatIsNonSymmetric = 1
             pCS%anyMatIsNonSymmetric = 1
             end if
          end if

          if (pCS%SSEEMatIsNonSymmetric == 1 .and. &
          &   pCS%stiffMatIsNonSymmetric == 1 .and. &
          &   pCS%dampMatIsNonSymmetric == 1 .and. &
          &   pCS%massMatIsNonSymmetric == 1) exit outer

       end do inner
    end do outer

!    !! Symmetrize
!    !! TODO,Magne: Create input option on this
!    do i = 1,pCS%nDofs
!       do j = i,pCS%nDofs
!            pCS%SSEEMat(i,j)  = (pCS%SSEEMat(i,j) + pCS%SSEEMat(j,i)) * 0.5_dp
!            pCS%SSEEMat(j,i)  =  pCS%SSEEMat(i,j)
!            pCS%stiffMat(i,j) = (pCS%stiffMat(i,j) + pCS%stiffMat(j,i)) * 0.5_dp
!            pCS%stiffMat(j,i) =  pCS%stiffMat(i,j)
!            pCS%dampMat(i,j)  = (pCS%dampMat(i,j) + pCS%dampMat(j,i)) * 0.5_dp
!            pCS%dampMat(j,i)  =  pCS%dampMat(i,j)
!            pCS%massMat(i,j)  = (pCS%massMat(i,j) + pCS%massMat(j,i)) * 0.5_dp
!            pCS%massMat(j,i)  =  pCS%massMat(i,j)
!       end do
!    end do

    deallocate(dScr)

  end subroutine BuildStructControlJacobi


  subroutine EstimateControllerProperties01(sys,mech,ctrl,msim,       &
    &                                   numVregIn, whichVregIn,   &
    &                                   numVregOut, whichVregOut, &
    &                                   nPerturb, ctrlProps, ierr)

    !!==========================================================================
    !! Purpose:
    !! Use a perturbation method, simular to the Matrix Stiffness Method /
    !! (Virtual) Displacement Method / Unit Load Method to find
    !! the equivalent mechanical properties of the controller. These controller properties
    !! will be added to existing matrices when conducting modal analysis / eigenvalue analysis.
    !! In this routine, the controller is limited to be of type PID.
    !! The equation for the system is:
    !! M*x'' + C*x' + K*x + Q*int(x)dt = F or M*(d2x/dt2) + C*(dx/dt) + K*X + Q*int(x)dt = F
    !! where M is mass, C is damping, K is stiffness and Q is steady state error elimination
    !! Controller input = y, controller output = u.
    !! The values of interest are:
    !! du/dy: Change du in output from controller with respect to
    !!        change dy in input to controller.
    !!        du/dy = proportional gain, Kp
    !! du/(int(dy)dt): Change du in output from controller with respect to
    !!                 change int(dy)dt in input to controller.
    !!                 du/(int(dy)dt) = integral gain, Ki
    !! du/(dy/dt): Change du in output from controller with respect to
    !!             change dy/dt in input to controller.
    !!             du/(dy/dt) = derivative gain, Kd
    !!
    !! Working order:
    !! 1) Do one initial perturbation on the controller with dy = 0 and dt /= 0.
    !!    This is to insure dy/dt = 0.
    !! 2) Get the initial values y0 and u0 for the controller.
    !! 3) Establish dy(j) and dt(j). j = number of perturbations.
    !!    For a PID-controller: j = 1...3.
    !! 4) Calculate d(int(dy(j)dt) and d(dy/dt).
    !! 5) Calculate y(j) and t(j).
    !! 6) Iterate the controller with these new values for the input y(j) and time t(j)
    !!    and save the reaction from the controller u(j) due to the change in the input.
    !! 7) Calculate du(j) based on u0 and u(j).
    !! 8) Calculate Kp, Ki and Kd.
    !! 9) Based on sensor type (position, velocity or acceleration), calculate Q, K, C and M.
    !!
    !! Programmer : Magne Bratland
    !! date/rev   : 25 Nov 2010 / 1.0
    !!==========================================================================

    use KindModule
    use ReportErrorModule
    use SystemTypeModule   ,   only : SystemType
    use MechanismTypeModule,   only : MechanismType
    use ControlTypeModule
    use ControlRoutinesModule, only : IterateControlSystem
    use SensorTypeModule
    use DenseMatrixModule, only : solveAxB

    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(in)    :: mech
    type(ControlType)  , intent(inout) :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: numVregIn         !! Number of vreg in to perturb
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: numVregOut        !! Number of vreg out to read variation from
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
    integer,             intent(in)    :: nPerturb          !! number of perturbations on controller
    real(dp),            intent(out)   :: ctrlProps(:,:,:)  !! table for storing controller properties
    !						  		                        !!(no. of outputs from controller,
    !				 				                        !! no. of controller properties,
    !		   			   			                        !! no. of inputs to controller)
    !                                                       !! ctrlProps(:,1,:) = Q
    !                                                       !! ctrlProps(:,2,:) = K
    !                                                       !! ctrlProps(:,3,:) = C
    !                                                       !! ctrlProps(:,4,:) = M

    integer            , intent(inout) :: ierr

    ! Local variables

    type(SensorType)   , pointer :: sensor
    type(ControlType)  , pointer :: ctrlCopy
    type(ControlType)  , pointer :: ctrlCopyCopy

    real(dp) :: orgTime                          !! original/initial value for the time
    real(dp) :: orgTimeStep                      !! original/initial value for the time step
    real(dp) :: t0                               !! quazi-initial value for the time
    real(dp) :: dt0                              !! quazi-initial time step
    real(dp) :: dt                               !! incremental time step
    real(dp) :: dy                               !! incremental controller input step
    real(dp) :: dintydt                          !! definite integral of y with respect to time t
    real(dp) :: ddydt                            !! d(dy/dt), 1st derivative of y
    !										     !! with respect to time t
    real(dp), allocatable :: y0(:)               !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)              !! u(y0), output from controller
    !										     !! when input to controller is y0
    real(dp), allocatable :: uy(:)               !! u(y), output from controller
    !										     !! when input to controller is y = y0+dy
    real(dp), allocatable :: du(:)               !! du = u(y)-u(y0)
    real(dp) :: dyMatrix(nPerturb,nPerturb)      !! matrix of perturbation parameters
    real(dp), allocatable :: duTable(:,:)        !! table for storing du = u(y)-u(y0)
    integer :: ctrlSysMode = 3                   !! ctrlSysMode = 3 = controller integration

    integer :: i, j, iInput

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    !! Reset counters
    i = 0   !! i used for numVregIn
    j = 0   !! j used for nPerturb
    iInput = 0

    !! Make copy of controller
    call AllocateCopyControlType(ctrl,ctrlCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return
    call AllocateCopyControlType(ctrlCopy,ctrlCopyCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return

    !! Reset variables
    t0      = 0.0_dp
    dt      = 0.0_dp
    dy      = 0.0_dp
    dintydt = 0.0_dp
    ddydt   = 0.0_dp

    !! Define length for arrays
    allocate(y0(numVregIn),uy0(numVregOut),uy(numVregOut),du(numVregOut))
    y0  = 0.0_dp
    uy0 = 0.0_dp
    uy  = 0.0_dp
    du  = 0.0_dp

    !! Define dimensions of tables
    allocate(duTable(nPerturb,numVregOut), STAT=ierr)

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Do one initial perturbation with dy = 0
    dt0 = orgTimeStep*1.0E-1_dp   !! TODO Magne: Change this value?
    t0 = orgTime+dt0

    do i = 1, numVregIn
       iInput = whichVregIn(i)
       sys%time = t0
       sys%timeStep = dt0
       dy = 0.0_dp
       call CopyControlType(ctrl,ctrlCopy)
       !! Initial perturbation
       call PerturbController(sys,ctrlCopy,msim,iInput,dt0,dy,numVregOut, &
            &                    whichVregOut,ctrlSysMode,uy0,ierr)
       !! Save the value of y0
       y0(i) = abs(ctrlCopy%vreg(iInput))

       !! Perturbation
       do j = 1, nPerturb
          !! Establish dy-matrix
          dt = orgTimeStep*1.0E-1_dp*j
          dy = dt
          dintydt = (y0(i)+0.5_dp*dy)*dt
          ddydt = dy/dt

          dyMatrix(j,1) = dintydt
          dyMatrix(j,2) = dy
          dyMatrix(j,3) = ddydt

          !! The perturbation sequence
          !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
          call CopyControlType(ctrlCopy,ctrlCopyCopy)
          sys%timeStep = dt
          uy  = 0.0_dp
          !! Perturb
          call PerturbController(sys,ctrlCopyCopy,msim,iInput,dt,dy,numVregOut, &
               &                    whichVregOut,ctrlSysMode,uy,ierr)
          if ( ierr /= 0 ) return
          !! Calculate du
          du(:) = uy(:)-uy0(:)
          !! Store du-results in a table
          duTable(j,:) = du(:)
       end do

       !! Calculate controller properties
       call solveAxB (dyMatrix,duTable,ierr)

       select case (ctrlCopy%input(iInput)%engine%args(1)%p%entity)
       case (POS_p) ! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = duTable(1,:)
          ctrlProps(:,2,i) = duTable(2,:)
          ctrlProps(:,3,i) = duTable(3,:)
       case (VEL_p) ! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,2,i) = duTable(1,:)
          ctrlProps(:,3,i) = duTable(2,:)
          ctrlProps(:,4,i) = duTable(3,:)
       case (ACC_p) ! find dintydt(j), dy(j)
          ctrlProps(:,3,i) = duTable(1,:)
          ctrlProps(:,4,i) = duTable(2,:)
       case default
          !! Error
       end select
    end do

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call DeallocateControlType(ctrlCopy,ierr)
    if ( ierr /= 0 ) return
    call DeallocateControlType(ctrlCopyCopy,ierr)
    if ( ierr /= 0 ) return

    deallocate(y0,uy0,uy,du,duTable)

  end subroutine EstimateControllerProperties01


  subroutine PerturbController(sys,ctrl,msim,iPerturb,dt,dy,numVregOut, &
       &                whichVregOut,ctrlSysMode,uy,ierr)

    !!==========================================================================
    !! Purpose:
    !! Perturb one of the inputs of the controller and get the
    !! reaction for all of the outputs of the controllers.
    !!
    !! Working order:
    !! 1) Change the input y for the controller from y0 to y = y0+dy,
    !!    where du is a small number. Change the time from t0 to t = t0+dt,
    !!    where dt is a small number.
    !! 2) Run this new input through the controller to get the reaction
    !!    (output) u(y) from the controller due to the change in the input.
    !!
    !! Programmer : Magne Bratland
    !! date/rev   : 18 Sept 2009 / 1.0
    !!==========================================================================


    use KindModule
    use ReportErrorModule
    use SystemTypeModule   , only : SystemType
    use ControlTypeModule  , only : ControlType
    use ControlRoutinesModule, only : IterateControlSystem
    use SensorTypeModule

    implicit none

    type(SystemType) , intent(inout) :: sys
    type(ControlType), intent(inout) :: ctrl
    integer          , intent(in)    :: msim(:)
    integer          , intent(in)    :: iPerturb          !! Which input to perturb
    real(dp)         , intent(in)    :: dt                !! incremental time step
    real(dp)         , intent(in)    :: dy                !! incremental step for use in du/dy
    integer          , intent(in)    :: numVregOut        !! number of outputs from the controller to read
    integer          , intent(in)    :: whichVregOut(:)   !! which outputs from the controller to read
    !								 					  !! which are outputs from the controller
    integer          , intent(in)    :: ctrlSysMode       !! ctrlSysMode = 3 = controller integration
    real(dp)		 , intent(out)   :: uy(:)             !! u(y) = u(y0+dy), output from controller
    !                                                     !! when input to controller is y = y0+dy

    integer          , intent(inout) :: ierr

    !! --- Local variables

    integer :: i, j

    logical :: bSetInput    !! switch to set if the controller should include the new value y=y+dy

    !! --- Logic section ---
    uy = 0.0_dp

    !! Change time step to dt
    sys%timeStep = dt

    !! Change input from y0 to y (y = y0+dy)
    i = ctrl%input(iPerturb)%var !!TODO,bh: Check this
    ctrl%vreg(i) = ctrl%vreg(i) + dy

    !! Iterate the controller with the new values for y (y = y0+dy) and t (t = t0+dt)
    !! to derive the value of u(y)
    bSetInput = .false.
    call IterateControlSystem (sys,ctrl,ctrlSysMode,msim,ierr,bSetInput)

    !! Save the value of u(y) in an array
    do i = 1, numVregOut
       uy(i) = ctrl%vreg(whichVregOut(i))
    end do

  end subroutine PerturbController


  subroutine addInControlStructMat (pCS,elMat,sysMat,sam,err)

    !!==========================================================================
    !! Add element matrix to system matrix
    !! date/rev   : 2009
    !!==========================================================================

    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType
    use AsmExtensionModule , only : csAddEM
    use reportErrorModule  , only : getErrorFile, reportError
    use reportErrorModule  , only : error_p, debugFileOnly_p

    type(ControlStructType), intent(inout) :: pCS
    real(dp)               , intent(in)    :: elMat(:,:)
    type(SysMatrixType)    , intent(inout) :: sysMat
    type(SamType)          , intent(in)    :: sam
    integer                , intent(out)   :: err

    !! Local variables
    integer           :: i, n, nDim, lpu

    !! --- Logic section ---

    lpu = getErrorFile()

    !! Add this matrix to the system matrix
    call csAddEM (sam,pCS%samElNum,pCS%nDofs,elMat,sysMat,err)
    if (err == 0) return

    call reportError (debugFileOnly_p,'AddInControlStructMat')

  end subroutine AddInControlStructMat


  subroutine EstimateControllerProperties02(sys,mech,ctrl,msim,       &
    &                                   numVregIn, whichVregIn,   &
    &                                   numVregOut, whichVregOut, &
    &                                   nPerturb, ctrlProps, ierr)

    !!==========================================================================
    !! Purpose:
    !! Use a perturbation method, simular to the Matrix Stiffness Method /
    !! (Virtual) Displacement Method / Unit Load Method to find
    !! the equivalent mechanical properties of the controller. These controller properties
    !! will be added to existing matrices when conducting modal analysis / eigenvalue analysis.
    !! There are a total of 4 values of interest:
    !! mass M, damping C, stiffness K and steady state error elimination Q.
    !! M*x'' + C*x' + K*x + Q*int(x)dt = F or M*(d2x/dt2) + C*(dx/dt) + K*X + Q*int(x)dt = F
    !! Controller input = y, controller output = u.
    !! The method is depending on sensor input (position x, velocity x' or acceleration x'')
    !! - If the sensor input is position:
    !!   du = Q*d(int(y)dt) + K*dy + C*d(dy/dt) + M*d(d2y/dt2)
    !! - If the sensor input is velocity:
    !!   du = Q*d(int(int(y)dt)dt) + K*d(int(y)dt)+ C*dy + M*d(dy/dt)
    !! - If the sensor input is acceleration:
    !!   du = Q*d(int(int(int(y)dt)dt)dt) + K*d(int(int(y)dt)dt) + C*d(int(y)dt)+ M*dy
    !!
    !! Working order:
    !! Per controller input, estimate all controller outputs:
    !! 1) Do one initial perturbation on the controller with dy = 0 and dt /= 0.
    !!    This is to insure dy0/dt = 0.
    !! 2) Get the initial values y0 and u0 for the controller.
    !! 3) Establish dy(j) and dt(j). j = number of perturbations. j = 1...4.
    !! 4) Depending on sensor input:
    !!    Position: Calculate d(int(y(j))dt), d(dy(j)/dt(j)) and d(d2y(j)/dt(j)2).
    !!    Velocity: Calculate d(int(int(y(j))dt)dt), d(int(y(j))dt) and d(dy(j)/dt(j))
    !!    Acceleration: Calculate d(int(int(int(y(j))dt)dt)dt), d(int(int(y(j))dt)dt)
    !!    and d(int(y(j))dt)
    !! 5) Calculate y(j) and t(j).
    !! 6) Iterate the controller with these new values for the input y(j) and time t(j)
    !!    and save the reaction from the controller u(j) due to the change in the input.
    !! 7) Calculate du(j) based on u0 and u(j).
    !! 8) Calculate Q, K, C and M.
    !!
    !! Programmer : Magne Bratland
    !! date/rev   : 01 June 2010 / 1.0
    !!==========================================================================

    use KindModule
    use ReportErrorModule
    use SystemTypeModule   ,   only : SystemType
    use MechanismTypeModule,   only : MechanismType
    use ControlTypeModule
    use ControlRoutinesModule, only : IterateControlSystem
    use SensorTypeModule
    use DenseMatrixModule, only : solveAxB

    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(in)    :: mech
    type(ControlType)  , intent(inout) :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: numVregIn         !! Number of vreg in to perturb
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: numVregOut        !! Number of vreg out to read variation from
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
    integer,             intent(in)    :: nPerturb          !! number of perturbations on controller
    real(dp),            intent(out)   :: ctrlProps(:,:,:)  !! table for storing controller properties
    !						  		                        !!(no. of outputs from controller,
    !				 				                        !! no. of controller properties,
    !		   			   			                        !! no. of inputs to controller)
    !                                                       !! ctrlProps(:,1,:) = Q
    !                                                       !! ctrlProps(:,2,:) = K
    !                                                       !! ctrlProps(:,3,:) = C
    !                                                       !! ctrlProps(:,4,:) = M

    integer            , intent(inout) :: ierr

    ! Local variables

    type(SensorType)   , pointer :: sensor
    type(ControlType)  , pointer :: ctrlCopy
    type(ControlType)  , pointer :: ctrlCopyCopy

    real(dp) :: orgTime                         !! original/initial value for the time
    real(dp) :: orgTimeStep                     !! original/initial value for the time step
    real(dp) :: t0                              !! quazi-initial value for the time
    real(dp) :: dt0                             !! quazi-initial time step
    real(dp) :: dy0                             !! quazi-initial input perturbation, dy0 = 0
    real(dp), allocatable :: y0(:)              !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)             !! u(y0), output from controller
    !										    !! when input to controller is y0
    real(dp), allocatable :: uy1(:)             !! u(y1)
    real(dp), allocatable :: uy(:)              !! u(y), output from controller
    !										    !! when input to controller is y = y0+dy
    real(dp), allocatable :: du(:)              !! du = u(y)-u(y0)
    real(dp), allocatable :: dt(:)              !! incremental time step
    real(dp), allocatable :: dy(:)              !! incremental controller input step
    real(dp), allocatable :: dy1(:)             !! incremental controller input step no. 1
    !                                           !! to be used when deriving d(d2y/dt2)
    real(dp), allocatable :: dy2(:)             !! incremental controller input step no. 2
    !                                           !! to be used when deriving d(d2y/dt2)
    real(dp), allocatable :: y1(:)              !! y1 = y0+dy1
    real(dp), allocatable :: ddydt(:)           !! d(dy/dt), 1st derivative of y
    !										    !! with respect to time t
    real(dp), allocatable :: dd2ydt2(:)         !! d(d2y/dt2), 2nd derivative of y
    !										    !! with respect to time t
    real(dp), allocatable :: dintydt(:)         !! definite integral of y with respect to time t
    real(dp), allocatable :: dintintydt(:)      !! definite double integral of y with respect to time t
    real(dp), allocatable :: dintintintydt(:)   !! definite triple integral of y with respect to time t
    real(dp) :: dyMatrix(nPerturb,nPerturb)     !! matrix of perturbation parameters
    real(dp), allocatable :: duTable(:,:)       !! table for storing du = u(y)-u(y0)
    integer :: ctrlSysMode = 3                  !! ctrlSysMode = 3 = controller integration

    integer :: i, j, iInput

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    !! Reset constants
    i = 0
    j = 0

    !! Make copy of controller
    call AllocateCopyControlType(ctrl,ctrlCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return
    call AllocateCopyControlType(ctrlCopy,ctrlCopyCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return

    !! Define length for arrays y0, uy0, uy1, uy and du
    allocate(y0(numVregIn),uy0(numVregOut),uy1(numVregOut),uy(numVregOut),du(numVregOut))
    y0 = 0.0_dp
    uy0 = 0.0_dp
    uy1 = 0.0_dp
    uy = 0.0_dp
    du = 0.0_dp

    !! Define dimensions for the du-table
    allocate(duTable(nPerturb,numVregOut), STAT=ierr)

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Do one time perturbation
    dt0 = sys%timeStep*1.0E-1_dp   !! TODO Magne: Change this value?
    t0 = sys%time+dt0
    sys%time = t0
    dy0 = 0.0_dp
    !! Start perturbation with time step dt0
    do i = 1, numVregIn
       !! Perturb
       iInput = whichVregIn(i)
       call PerturbController(sys,ctrlCopy,msim,iInput,dt0,dy0,numVregOut, &
            &                    whichVregOut,ctrlSysMode,uy,ierr)
    end do

    !! Save the value of y0 in an array
    do i = 1, numVregIn
       y0(i) = ctrlCopy%vreg(whichVregIn(i))
    end do

    !! Save the value of u(y0) in an array
    do i = 1, numVregOut
       uy0(i) = ctrlCopy%vreg(whichVregOut(i))
    end do

    !! Depending on sensor input, there are 3 various ways to perturb the system
    do i = 1, numVregIn
       sensor => ctrlCopy%input(i)%engine%args(1)%p
       if (sensor%entity == POS_p) then
          allocate(dt(nPerturb),dy1(nPerturb),dy2(nPerturb),y1(nPerturb),dintydt(nPerturb), &
               &   ddydt(nPerturb),dd2ydt2(nPerturb))

          do j = 1, nPerturb

             !! Establish dy-matrix
             dt(j) = orgTimeStep*1.0E-1_dp*j
             dy1(j) = dt(j)
             dy2(j) = dy1(j)*(-1)
             y1(j) = y0(i)+dy1(j)
             dintydt(j) = (y1(j)+(1.0_dp/2.0_dp)*dy2(j))*dt(j)
             ddydt(j) = (dy2(j)-dy1(j))/dt(j)
             dd2ydt2(j) = (dy2(j)-2*dy1(j))/dt(j)**2

             dyMatrix(j,1) = dintydt(j)
             dyMatrix(j,2) = dy2(j)
             dyMatrix(j,3) = ddydt(j)
             dyMatrix(j,4) = dd2ydt2(j)

             !! The perturbation sequence
             !! To derive d(d2y/dt2), the system has to be perturbed three times (two + initial)
             iInput = whichVregIn(i)
             !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
             call CopyControlType(ctrlCopy,ctrlCopyCopy)
             !! First perturbation
             call PerturbController(sys,ctrlCopyCopy,msim,iInput,dt(j),dy1(j),numVregOut, &
                  &                    whichVregOut,ctrlSysMode,uy,ierr)
             if ( ierr /= 0 ) return
             !! Save the value of u(y1) in an array
             uy1(:) = uy(:)
             !! Second perturbation
             call PerturbController(sys,ctrlCopyCopy,msim,iInput,dt(j),dy2(j),numVregOut, &
                  &                    whichVregOut,ctrlSysMode,uy,ierr)
             if ( ierr /= 0 ) return

             !! Calculate du
             du(:) = uy(:)-uy1(:)
             !! Store du-results in a table
             duTable(j,:) = du(:)
             !! Reset time
             sys%time = t0
          end do
          deallocate(dt,dy1,dy2,y1,dintydt,ddydt,dd2ydt2)

       else if (sensor%entity == VEL_p) then
          allocate(dt(nPerturb),dy(nPerturb),dintydt(nPerturb),dintintydt(nPerturb),ddydt(nPerturb))

          !! Establish dy-matrix
          do j = 1, nPerturb
             dt(j) = orgTimeStep*1.0E-1_dp*j
             dy(j) = dt(j)
             dintydt(j) = (y0(i)+(1.0_dp/2.0_dp)*dy(j))*dt(j)
             dintintydt(j) = ((1.0_dp/2.0_dp)*y0(i)+(1.0_dp/6.0_dp)*dy(j))*dt(j)**2
             ddydt(j) = dy(j)/dt(j)

             dyMatrix(j,1) = dintintydt(j)
             dyMatrix(j,2) = dintydt(j)
             dyMatrix(j,3) = dy(j)
             dyMatrix(j,4) = ddydt(j)

             !! The perturbation sequence
             iInput = whichVregIn(i)
             !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
             call CopyControlType(ctrlCopy,ctrlCopyCopy)
             !! Perturb
             call PerturbController(sys,ctrlCopyCopy,msim,iInput,dt(j),dy(j),numVregOut, &
                  &                    whichVregOut,ctrlSysMode,uy,ierr)
             if ( ierr /= 0 ) return

             !! Calculate du
             du(:) = uy(:)-uy0(:)
             !! Store du-results in a table
             duTable(j,:) = du(:)
             !! Reset time
             sys%time = t0
          end do
          deallocate(dt,dy,dintydt,dintintydt,ddydt)

       else if (sensor%entity == ACC_p) then
          allocate(dt(nPerturb),dy(nPerturb),dintydt(nPerturb),dintintydt(nPerturb), &
               &   dintintintydt(nPerturb))

          !! Establish dy-matrix
          do j = 1, nPerturb
             dt(j) = orgTimeStep*1.0E-1_dp*j
             dy(j) = dt(j)
             dintydt(j) = (y0(i)+(1.0_dp/2.0_dp)*dy(j))*dt(j)
             dintintydt(j) = ((1.0_dp/2.0_dp)*y0(i)+(1.0_dp/6.0_dp)*dy(j))*dt(j)**2
             dintintintydt(j) = ((1.0_dp/6.0_dp)*y0(i)+(1.0_dp/24.0_dp)*dy(j))*dt(j)**3

             dyMatrix(j,1) = dintintintydt(j)
             dyMatrix(j,2) = dintintydt(j)
             dyMatrix(j,3) = dintydt(j)
             dyMatrix(j,4) = dy(j)

             !! The perturbation sequence
             iInput = whichVregIn(i)
             !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
             call CopyControlType(ctrlCopy,ctrlCopyCopy)
             !! Perturb
             call PerturbController(sys,ctrlCopyCopy,msim,iInput,dt(j),dy(j),numVregOut, &
                  &                    whichVregOut,ctrlSysMode,uy,ierr)
             if ( ierr /= 0 ) return

             !! Calculate du
             du(:) = uy(:)-uy0(:)
             !! Store du-results in a table
             duTable(j,:) = du(:)
             !! Reset time
             sys%time = t0
          end do
          deallocate(dt,dy,dintydt,dintintydt,dintintintydt)

       else !! Error

       end if

       !! Calculate controller properties
       call solveAxB (dyMatrix,duTable,ierr)

       ctrlProps(:,:,i) = transpose(duTable)
    end do

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call DeallocateControlType(ctrlCopy,ierr)
    if ( ierr /= 0 ) return
    call DeallocateControlType(ctrlCopyCopy,ierr)
    if ( ierr /= 0 ) return

    deallocate(y0,uy0,uy1,uy,du,duTable)

  end subroutine EstimateControllerProperties02


subroutine EstimateControllerProperties03(sys,mech,ctrl,msim,       &
    &                                   numVregIn, whichVregIn,   &
    &                                   numVregOut, whichVregOut, &
    &                                   nPerturb, ctrlProps, ierr)

    !!==========================================================================
    !! Purpose:
    !! Use a perturbation method, simular to the Matrix Stiffness Method /
    !! (Virtual) Displacement Method / Unit Load Method to find
    !! the equivalent mechanical properties of the controller. These controller properties
    !! will be added to existing matrices when conducting modal analysis / eigenvalue analysis.
    !! There are a total of 4 values of interest:
    !! mass M, damping C, stiffness K and steady state error elimination Q.
    !! M*x'' + C*x' + K*x + Q*int(x)dt = F or M*(d2x/dt2) + C*(dx/dt) + K*X + Q*int(x)dt = F
    !! Controller input = y, controller output = u.
    !! The method is depending on sensor input (position x, velocity x' or acceleration x'')
    !! - If the sensor input is position:
    !!   du = Q*d(int(y)dt) + K*dy + C*d(dy/dt) + M*d(d2y/dt2)
    !! - If the sensor input is velocity:
    !!   du = Q*d(int(int(y)dt)dt) + K*d(int(y)dt)+ C*dy + M*d(dy/dt)
    !! - If the sensor input is acceleration:
    !!   du = Q*d(int(int(int(y)dt)dt)dt) + K*d(int(int(y)dt)dt) + C*d(int(y)dt)+ M*dy
    !! This is a just-in-case-algorithm; all possible parameters of interest will be derived:
    !! - d(int(int(int(y(j))dt)dt)dt)
    !! - d(int(int(y(j))dt)dt)
    !! - d(int(y(j))dt
    !! - dy(j)
    !! - d(dy(j)/dt(j))
    !! - d(d2y(j)/dt(j)2)
    !!
    !! Working order:
    !! Per controller input, estimate all controller outputs:
    !! 1) Do one initial perturbation on the controller with dy = 0 and dt /= 0.
    !!    This is to insure dy0/dt = 0.
    !! 2) Get the initial values y0 and u0 for the controller.
    !! 3) Establish dy(j) and dt(j). j = number of perturbations. j = 1...6.
    !! 4) Calculate d(int(int(int(y(j))dt)dt)dt), d(int(int(y(j))dt)dt), d(int(y(j))dt)
    !!    d(dy(j)/dt(j)) and d(d2y(j)/dt(j)2).
    !! 5) Calculate y(j) and t(j).
    !! 6) Iterate the controller with these new values for the input y(j) and time t(j)
    !!    and save the reaction from the controller u(j) due to the change in the input.
    !! 7) Calculate du(j) based on u0 and u(j).
    !! 8) Calculate Q, K, C and M.
    !!
    !! Programmer : Magne Bratland
    !! date/rev   : 01 June 2010 / 1.0
    !!==========================================================================

    use KindModule
    use ReportErrorModule
    use SystemTypeModule   ,   only : SystemType
    use MechanismTypeModule,   only : MechanismType
    use ControlTypeModule
    use ControlRoutinesModule, only : IterateControlSystem
    use SensorTypeModule
    use DenseMatrixModule, only : solveAxB

    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(in)    :: mech
    type(ControlType)  , intent(inout) :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: numVregIn         !! Number of vreg in to perturb
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: numVregOut        !! Number of vreg out to read variation from
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
    integer,             intent(in)    :: nPerturb          !! number of perturbations on controller
    real(dp),            intent(out)   :: ctrlProps(:,:,:)  !! table for storing controller properties
    !						  		                        !!(no. of outputs from controller,
    !				 				                        !! no. of controller properties,
    !		   			   			                        !! no. of inputs to controller)
    !                                                       !! ctrlProps(:,1,:) = Q
    !                                                       !! ctrlProps(:,2,:) = K
    !                                                       !! ctrlProps(:,3,:) = C
    !                                                       !! ctrlProps(:,4,:) = M

    integer            , intent(inout) :: ierr

    ! Local variables

    type(SensorType)   , pointer :: sensor
    type(ControlType)  , pointer :: ctrlCopy
    type(ControlType)  , pointer :: ctrlCopyCopy

    real(dp) :: orgTime                          !! original/initial value for the time
    real(dp) :: orgTimeStep                      !! original/initial value for the time step
    real(dp) :: t0                               !! quazi-initial value for the time
    real(dp) :: dt0                              !! quazi-initial time step
    real(dp) :: dy0                              !! quazi-initial input perturbation, dy0 = 0
    real(dp), allocatable :: y0(:)               !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)              !! u(y0), output from controller
    !										     !! when input to controller is y0
    real(dp), allocatable :: uy1(:)              !! u(y1)
    real(dp), allocatable :: uy(:)               !! u(y), output from controller
    !										     !! when input to controller is y = y0+dy
    real(dp), allocatable :: du(:)               !! du = u(y)-u(y0)
    real(dp), allocatable :: dt(:)               !! incremental time step
    real(dp), allocatable :: dy(:)               !! incremental controller input step
    real(dp), allocatable :: dy1(:)              !! incremental controller input step no. 1
    !                                            !! to be used when deriving d(d2y/dt2)
    real(dp), allocatable :: dy2(:)              !! incremental controller input step no. 2
    !                                            !! to be used when deriving d(d2y/dt2)
    real(dp), allocatable :: y1(:)               !! y1 = y0+dy1
    real(dp), allocatable :: ddydt(:)            !! d(dy/dt), 1st derivative of y
    !										     !! with respect to time t
    real(dp), allocatable :: dd2ydt2(:)          !! d(d2y/dt2), 2nd derivative of y
    !										     !! with respect to time t
    real(dp), allocatable :: dintydt(:)          !! definite integral of y with respect to time t
    real(dp), allocatable :: dintintydt(:)       !! definite double integral of y with respect to time t
    real(dp), allocatable :: dintintintydt(:)    !! definite triple integral of y with respect to time t
    real(dp) :: dyMatrix(nPerturb,nPerturb)      !! matrix of perturbation parameters
    real(dp), allocatable :: duTable(:,:)        !! table for storing du-results
    integer :: ctrlSysMode = 3                   !! ctrlSysMode = 3 = controller integration

    integer :: i, j, iInput

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    !! Reset constants
    i = 0
    j = 0

    !! Make copy of controller
    call AllocateCopyControlType(ctrl,ctrlCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return
    call AllocateCopyControlType(ctrlCopy,ctrlCopyCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return

    !! Define length for arrays y0, uy0, uy1, uy and du
    allocate(y0(numVregIn),uy0(numVregOut),uy1(numVregOut),uy(numVregOut),du(numVregOut))
    y0 = 0.0_dp
    uy0 = 0.0_dp
    uy1 = 0.0_dp
    uy = 0.0_dp
    du = 0.0_dp

    !! Define dimensions of tables
    allocate(duTable(nPerturb,numVregOut), STAT=ierr)

    allocate(dt(nPerturb),dy1(nPerturb),dy2(nPerturb),y1(nPerturb),          &
         &   dintydt(nPerturb),dintintydt(nPerturb),dintintintydt(nPerturb), &
         &   ddydt(nPerturb),dd2ydt2(nPerturb))

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Do one time perturbation
    dt0 = sys%timeStep*1.0E-1_dp   !! TODO Magne: Change this value?
    t0 = sys%time+dt0
    sys%time = t0
    dy0 = 0.0_dp
    !! Start perturbation with time step dt0
    do i = 1, numVregIn
       !! Perturb
       iInput = whichVregIn(i)
       call PerturbController(sys,ctrlCopy,msim,iInput,dt0,dy0,numVregOut, &
            &                    whichVregOut,ctrlSysMode,uy,ierr)
    end do

    !! Save the value of y0 in an array
    do i = 1, numVregIn
       y0(i) = ctrlCopy%vreg(whichVregIn(i))
    end do

    !! Save the value of u(y0) in an array
    do i = 1, numVregOut
       uy0(i) = ctrlCopy%vreg(whichVregOut(i))
    end do

    !! New way: do to many perturbations and exlude redundant results. Do then
    !! 6 three-step perturbations

    do i = 1, numVregIn
       do j = 1, nPerturb
          !! Establish dy-matrix
          dt(j) = orgTimeStep*1.0E-1_dp*j
          dy1(j) = dt(j)
          dy2(j) = dy1(j)*(-1)
          y1(j) = y0(i)+dy1(j)
          dintydt(j) = (y1(j)+(1.0_dp/2.0_dp)*dy2(j))*dt(j)
          dintintydt(j) = ((1.0_dp/2.0_dp)*y1(j)+(1.0_dp/6.0_dp)*dy2(j))*dt(j)**2
          dintintintydt(j) = ((1.0_dp/6.0_dp)*y1(j)+(1.0_dp/24.0_dp)*dy2(j))*dt(j)**3
          ddydt(j) = (dy2(j)-dy1(j))/dt(j)
          dd2ydt2(j) = (dy2(j)-2*dy1(j))/dt(j)**2

          dyMatrix(j,1) = dintintintydt(j)
          dyMatrix(j,2) = dintintydt(j)
          dyMatrix(j,3) = dintydt(j)
          dyMatrix(j,4) = dy2(j)
          dyMatrix(j,5) = ddydt(j)
          dyMatrix(j,6) = dd2ydt2(j)

          !! The perturbation sequence
          !! To derive d(d2y/dt2), the system has to be perturbed three times (two + initial)
          iInput = whichVregIn(i)
          !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
          call CopyControlType(ctrlCopy,ctrlCopyCopy)
          !! First perturbation
          call PerturbController(sys,ctrlCopyCopy,msim,iInput,dt(j),dy1(j),numVregOut, &
               &                    whichVregOut,ctrlSysMode,uy,ierr)
          if ( ierr /= 0 ) return
          !! Save the value of u(y1) in an array
          uy1(:) = uy(:)
          !! Second perturbation
          call PerturbController(sys,ctrlCopyCopy,msim,iInput,dt(j),dy2(j),numVregOut, &
               &                    whichVregOut,ctrlSysMode,uy,ierr)
          if ( ierr /= 0 ) return

          !! Calculate du
          du(:) = uy(:)-uy1(:)
          !! Store du-results in a table
          duTable(j,:) = du(:)
          !! Reset time
          sys%time = t0
       end do

       !! Calculate controller properties
       call solveAxB (dyMatrix,duTable,ierr)

       select case (ctrlCopy%input(i)%engine%args(1)%p%entity)
       case (POS_p) ! dintydt(j), dy(j), ddydt(j) and dd2ydt2(j)
          ctrlProps(:,1,i) = duTable(3,:)
          ctrlProps(:,2,i) = duTable(4,:)
          ctrlProps(:,3,i) = duTable(5,:)
          ctrlProps(:,4,i) = duTable(6,:)
       case (VEL_p) ! dintintydt(j), dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = duTable(2,:)
          ctrlProps(:,2,i) = duTable(3,:)
          ctrlProps(:,3,i) = duTable(4,:)
          ctrlProps(:,4,i) = duTable(5,:)
       case (ACC_p) ! dintintintydt(j), dintintydt(j), dintydt(j), dy(j)
          ctrlProps(:,1,i) = duTable(1,:)
          ctrlProps(:,2,i) = duTable(2,:)
          ctrlProps(:,3,i) = duTable(3,:)
          ctrlProps(:,4,i) = duTable(4,:)
       case default
          !! Error
       end select

    end do

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call DeallocateControlType(ctrlCopy,ierr)
    if ( ierr /= 0 ) return
    call DeallocateControlType(ctrlCopyCopy,ierr)
    if ( ierr /= 0 ) return

    deallocate(uy0,uy1,uy,duTable)
    deallocate(y0,du,dt,dy1,dy2,y1,dintydt,dintintydt, &
           &   dintintintydt,ddydt,dd2ydt2)

  end subroutine EstimateControllerProperties03


  subroutine EstimateControllerProperties04 (sys, ctrl, msim, &
       &                                     numVregIn, whichVregIn, &
       &                                     numVregOut, whichVregOut, &
       &                                     nStep, ctrlProps, ierr)

    !!==========================================================================
    !! Purpose:
    !! Use a perturbation method, simular to the Matrix Stiffness Method /
    !! (Virtual) Displacement Method / Unit Load Method to find
    !! the equivalent mechanical properties of the controller. These controller properties
    !! will be added to existing matrices when conducting modal analysis / eigenvalue analysis.
    !! In this routine, the controller is limited to be of type PID. This algorithm has additional
    !! steps in between perturbations to check for accuracy.
    !! The equation for the system is:
    !! M*x'' + C*x' + K*x + Q*int(x)dt = F or M*(d2x/dt2) + C*(dx/dt) + K*X + Q*int(x)dt = F
    !! where M is mass, C is damping, K is stiffness and Q is steady state error elimination
    !! Controller input = y, controller output = u.
    !! The values of interest are:
    !! du/dy: Change du in output from controller with respect to
    !!        change dy in input to controller.
    !!        du/dy = proportional gain, Kp
    !! du/(int(dy)dt): Change du in output from controller with respect to
    !!                 change int(dy)dt in input to controller.
    !!                 du/(int(dy)dt) = integral gain, Ki
    !! du/(dy/dt): Change du in output from controller with respect to
    !!             change dy/dt in input to controller.
    !!             du/(dy/dt) = derivative gain, Kd
    !!
    !! Working order:
    !! 1) Do one initial perturbation on the controller with dy = 0 and dt /= 0.
    !!    This is to insure dy/dt = 0.
    !! 2) Get the initial values y0 and u0 for the controller.
    !! 3) Establish dy(j) and dt(j). j = number of perturbations.
    !!    For a PID-controller: j = 1...3.
    !! 4) Calculate d(int(dy(j)dt) and d(dy/dt).
    !! 5) Calculate y(j) and t(j).
    !! 6) Iterate the controller with these new values for the input y(j) and time t(j)
    !!    and save the reaction from the controller u(j) due to the change in the input.
    !! 7) Calculate du(j) based on u0 and u(j).
    !! 8) Calculate Kp, Ki and Kd.
    !! 9) Based on sensor type (position, velocity or acceleration), calculate Q, K, C and M.
    !!
    !! Programmer : Magne Bratland
    !! date/rev   : 03 June 2010 / 1.0
    !!==========================================================================

    use KindModule
    use ReportErrorModule
    use SystemTypeModule   ,   only : SystemType
    use ControlTypeModule
    use ControlRoutinesModule, only : IterateControlSystem
    use SensorTypeModule
    use DenseMatrixModule, only : solveAxB

    type(SystemType)   , intent(inout) :: sys
    type(ControlType)  , intent(inout) :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: numVregIn         !! Number of vreg in to perturb
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: numVregOut        !! Number of vreg out to read variation from
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
    integer,             intent(in)    :: nStep             !! number of steps in between perturbations
    real(dp),            intent(out)   :: ctrlProps(:,:,:)  !! table for storing controller properties
    !						  		                        !!(no. of outputs from controller,
    !				 				                        !! no. of controller properties,
    !		   			   			                        !! no. of inputs to controller)
    !                                                       !! ctrlProps(:,1,:) = Q
    !                                                       !! ctrlProps(:,2,:) = K
    !                                                       !! ctrlProps(:,3,:) = C
    !                                                       !! ctrlProps(:,4,:) = M

    integer            , intent(inout) :: ierr

    ! Local variables

    integer, parameter :: nPerturb = 5 !< Number of controller perturbations

    type(SensorType)   , pointer :: sensor
    type(ControlType)  , pointer :: ctrlCopy
    type(ControlType)  , pointer :: ctrlCopyCopy

    real(dp) :: orgTime                          !! original/initial value for the time
    real(dp) :: orgTimeStep                      !! original/initial value for the time step
    real(dp) :: t0                               !! quazi-initial value for the time
    real(dp) :: dt0                              !! quazi-initial time step
    real(dp) :: dy0                              !! quazi-initial input perturbation, dy0 = 0
    real(dp) :: ddt                              !! incremental steps of dt
    real(dp) :: ddy                              !! incremental steps of dy
    real(dp), allocatable :: y0(:)               !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)              !! u(y0), output from controller
    !										     !! when input to controller is y0
    real(dp), allocatable :: uy(:)               !! u(y), output from controller
    !										     !! when input to controller is y = y0+dy
    real(dp), allocatable :: du(:)               !! du = u(y)-u(y0)
    real(dp), allocatable :: dt(:)               !! incremental time step
    real(dp), allocatable :: dy(:)               !! incremental controller input step
    real(dp), allocatable :: ddydt(:)            !! d(dy/dt), 1st derivative of y
    !										     !! with respect to time t
    real(dp), allocatable :: dintydt(:)          !! definite integral of y with respect to time t
    real(dp), allocatable :: dintintydt(:)       !! definite double integral of y with respect to time t
    real(dp), allocatable :: dintintintydt(:)    !! definite triple integral of y with respect to time t
    real(dp) :: dyMatrix(nPerturb,nPerturb)      !! matrix of perturbation parameters
    real(dp), allocatable :: duTable(:,:)        !! table for storing du = u(y)-u(y0)
    integer :: ctrlSysMode = 3                   !! ctrlSysMode = 3 = controller integration

    integer :: i, ii, j, iInput

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    !! Reset constants
    i = 0
    ii = 0
    j = 0

    !! Make copy of controller
    call AllocateCopyControlType(ctrl,ctrlCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return
    call AllocateCopyControlType(ctrlCopy,ctrlCopyCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return

    !! Define length for arrays y0, uy0, uy1, uy and du
    allocate(y0(numVregIn),uy0(numVregOut),uy(numVregOut),du(numVregOut))
    y0 = 0.0_dp
    uy0 = 0.0_dp
    uy = 0.0_dp
    du = 0.0_dp

    !! Define dimensions of tables
    allocate(duTable(nPerturb,numVregOut), STAT=ierr)

    allocate(dt(nPerturb),dy(nPerturb),dintydt(nPerturb),dintintydt(nPerturb),dintintintydt(nPerturb),ddydt(nPerturb))

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Do one initial time perturbation
    dt0 = sys%timeStep*1.0E-1_dp   !! TODO Magne: Change this value?
    t0 = sys%time+dt0
    sys%time = t0
    dy0 = 0.0_dp
    !! Start perturbation with time step dt0
    do i = 1, numVregIn
       !! Perturb
       iInput = whichVregIn(i)
       ddt = dt0/real(nStep,dp)
       do ii = 1, nStep
          call PerturbController(sys,ctrlCopy,msim,iInput,ddt,dy0,numVregOut, &
               &                    whichVregOut,ctrlSysMode,uy,ierr)
       end do
    end do

    !! Save the value of y0 in an array
    do i = 1, numVregIn
       y0(i) = ctrlCopy%vreg(whichVregIn(i))
    end do

    !! Save the value of u(y0) in an array
    do i = 1, numVregOut
       uy0(i) = ctrlCopy%vreg(whichVregOut(i))
    end do

    !! New way: do to many perturbations and exlude redundant results. Do then
    !! 6 three-step perturbations

    do i = 1, numVregIn
       do j = 1, nPerturb
          !! Establish dy-matrix
          dt(j)            = orgTimeStep*1.0E-1_dp*j
          dy(j)            = dt(j)
          dintydt(j)       = ((1.0_dp/1.0_dp)*y0(i)+(1.0_dp/2.0_dp)*dy(j))*dt(j)
          dintintydt(j)    = ((1.0_dp/2.0_dp)*y0(i)+(1.0_dp/6.0_dp)*dy(j))*dt(j)**2
          dintintintydt(j) = ((1.0_dp/6.0_dp)*y0(i)+(1.0_dp/24.0_dp)*dy(j))*dt(j)**3
          ddydt(j)         = dy(j)/dt(j)

          dyMatrix(j,1) = dintintintydt(j)
          dyMatrix(j,2) = dintintydt(j)
          dyMatrix(j,3) = dintydt(j)
          dyMatrix(j,4) = dy(j)
          dyMatrix(j,5) = ddydt(j)

          !! The perturbation sequence
          iInput = whichVregIn(i)
          !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
          call CopyControlType(ctrlCopy,ctrlCopyCopy)
          !! First perturbation
          ddt = dt(j)/real(nStep,dp)
          ddy = dy(j)/real(nStep,dp)
          do ii = 1, nStep
             call PerturbController(sys,ctrlCopyCopy,msim,iInput,ddt,ddy,numVregOut, &
                  &                    whichVregOut,ctrlSysMode,uy,ierr)
             if ( ierr /= 0 ) return
          end do

          !! Calculate du
          du(:) = uy(:)-uy0(:)
          !! Store du-results in a table
          duTable(j,:) = du(:)
          !! Reset time
          sys%time = orgTime+dt0
       end do

       !! Calculate controller properties
       call solveAxB (dyMatrix,duTable,ierr)

       select case (ctrlCopy%input(i)%engine%args(1)%p%entity)
       case (POS_p) ! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = duTable(3,:)
          ctrlProps(:,2,i) = duTable(4,:)
          ctrlProps(:,3,i) = duTable(5,:)
       case (VEL_p) ! find dintintydt(j), dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = duTable(2,:)
          ctrlProps(:,2,i) = duTable(3,:)
          ctrlProps(:,3,i) = duTable(4,:)
          ctrlProps(:,4,i) = duTable(5,:)
       case (ACC_p) ! find dintintintydt(j),dintintydt(j),dintydt(j), dy(j)
          ctrlProps(:,1,i) = duTable(1,:)
          ctrlProps(:,2,i) = duTable(2,:)
          ctrlProps(:,3,i) = duTable(3,:)
          ctrlProps(:,4,i) = duTable(4,:)
       case default
          !! Error
       end select

    end do

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call DeallocateControlType(ctrlCopy,ierr)
    if ( ierr /= 0 ) return
    call DeallocateControlType(ctrlCopyCopy,ierr)
    if ( ierr /= 0 ) return

    deallocate(uy0,uy,duTable)
    deallocate(y0,du,dt,dintydt,dintintydt,dintintintydt,ddydt)

  end subroutine EstimateControllerProperties04


  !!============================================================================
  !> @brief Perturbation method without initial perturbation.
  !> @details Perturbation w.r.t. a single integral.
  !> @a nPerturb = 1 in this subroutine, therefore no loop necessary.
  !>
  !> @author Magne Bratland
  !>
  !> @date 7 July 2010

  subroutine EstimateControllerProperties500 (sys, ctrl, msim, &
       &                                      numVregIn, whichVregIn, &
       &                                      numVregOut, whichVregOut, &
       &                                      ctrlProps, ierr)

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p

    implicit none

    type(SystemType) , intent(inout) :: sys
    type(ControlType), intent(inout) :: ctrl
    integer          , intent(in)    :: msim(:)
    integer          , intent(in)    :: numVregIn        !! Number of vreg in to perturb
    integer          , intent(in)    :: whichVregIn(:)   !! Which vreg in to perturb
    integer          , intent(in)    :: numVregOut       !! Number of vreg out to read variation from
    integer          , intent(in)    :: whichVregOut(:)  !! Which vreg out to read variation from
    real(dp)         , intent(out)   :: ctrlProps(:,:,:) !! table for storing controller properties
    integer          , intent(inout) :: ierr

    !! Local variables

    type(ControlType), pointer :: ctrlCopy
    type(ControlType), pointer :: ctrlCopyCopy

    real(dp) :: y0          ! initial values of the controller inputs
    real(dp), allocatable :: uy0(:) ! u(y0), output from controller when input to controller is y0
    real(dp), allocatable :: uy(:)  ! u(y), output from controller when input to controller is y = y0+dy
    real(dp) :: orgTime     ! original/initial value for the time
    real(dp) :: orgTimeStep ! original/initial value for the time step
    real(dp) :: dt          ! incremental time step
    real(dp) :: dy          ! incremental controller input step
    real(dp) :: ddydt       ! d(dy/dt), 1st derivative of y with respect to time t
    real(dp) :: dintydt     ! definite integral of y with respect to time t

    integer, parameter :: ctrlSysMode = 3 ! = controller integration

    integer :: i, iEntity

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    !! Make copy of controller
    call AllocateCopyControlType(ctrl,ctrlCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return
    call AllocateCopyControlType(ctrlCopy,ctrlCopyCopy,ierr) !! Allocate and copy
    if ( ierr /= 0 ) return

    allocate(uy0(numVregOut),uy(numVregOut), STAT=ierr)

    !! Store the initial time values
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Do one initial time perturbation
    sys%time = sys%time + sys%timeStep

    !! Save the value of u(y0) in an array
    do i = 1, numVregOut
       uy0(i) = ctrlCopy%vreg(whichVregOut(i))
    end do

    !! Perturbation
    do i = 1, numVregIn
       y0 = ctrlCopy%vreg(whichVregIn(i))

       !! Establish dy-matrix
       dt      = orgTimeStep
       dy      = dt
       ddydt   = dy/dt
       dintydt = (y0 + dy/2.0_dp)*dt

       !! The perturbation sequence
       !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
       call CopyControlType(ctrlCopy,ctrlCopyCopy)
       !! Perturb the controller
       call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy, &
            &                  numVregOut,whichVregOut,ctrlSysMode,uy,ierr)
       if ( ierr /= 0 ) return

       iEntity = ctrlCopy%input(i)%engine%args(1)%p%entity
       select case (iEntity)
       case (POS_p,VEL_p,ACC_p)
          !! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,iEntity,i) = (1.0_dp/dintydt) * (uy - uy0)
       case default ! Error
       end select

    end do

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call DeallocateControlType(ctrlCopy,ierr)
    if ( ierr /= 0 ) return
    call DeallocateControlType(ctrlCopyCopy,ierr)
    if ( ierr /= 0 ) return

    deallocate(uy0,uy)

  end subroutine EstimateControllerProperties500

end module ControlStructModule
