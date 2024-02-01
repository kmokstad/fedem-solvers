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

  use ControlTypeModule, only : CtrlPrm
  use SensorTypeModule , only : SensorType, IdType, dp
  use ForceTypeModule  , only : ForceType

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
     type(CtrlPrm)    , pointer :: pCtrlPrm        !! pointer to the control parameter which has
     !                                             !! structural input
     type(SensorType),  pointer :: pSensor         !! The actual sensor pointer

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
  end type ForceControlGradType

  type ControlStructType

     integer :: ctrlSysEigFlag

     integer :: samElNum
     integer :: nDOFs
     integer, pointer :: samMNPC(:)
     integer, pointer :: local_MADOF(:)

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

  end type ControlStructType


contains

  !!============================================================================
  !> @brief Initializes the control struct data type.
  !>

  subroutine InitiateControlStruct (pCS,inputs,triads,joints,forces, &
       &                            numVreg,numNod,ierr)

    use TriadTypeModule           , only : TriadType
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use SensorTypeModule          , only : TRIAD_p, RELATIVE_TRIAD_p
    use SensorTypeModule          , only : JOINT_VARIABLE_p
    use ForceTypeModule           , only : ForceType
    use ControlTypeModule         , only : getSensor
    use ReportErrorModule         , only : allocationError
    use ReportErrorModule         , only : reportError, debugFileOnly_p

    type(ControlStructType)   , intent(inout)      :: pCS
    type(CtrlPrm)             , intent(in), target :: inputs(:)
    type(TriadType)           , intent(in), target :: triads(:)
    type(MasterSlaveJointType), intent(in), target :: joints(:)
    type(ForceType)           , intent(in), target :: forces(:)
    integer                   , intent(in)         :: numVreg, numNod
    integer                   , intent(out)        :: ierr

    !! Local variables

    type(TriadType)           , pointer :: pTriad1, pTriad2
    type(MasterSlaveJointType), pointer :: pJoint
    type(StructSensorGradType), pointer :: pS
    type(ForceControlGradType), pointer :: pF

    integer, allocatable :: lNode_from_SAM_node(:)
    integer, allocatable :: iCin_from_allVreg(:)
    integer, allocatable :: iCout_from_allVreg(:)
    integer, pointer     :: tmpIdx(:)

    integer :: i, n, iCin, iCout, iSAM, nElNodes, iStart, whichVreg

    !! --- Logic section ---

    ierr = 0

    !! Find the control input parameters that have structural sensors
    !! (triad- and joint dofs only)

    call getCtrlParamsWithStructSensors (inputs,pCS%numStructCtrlParams,tmpIdx)
    if (pCS%numStructCtrlParams < 0) goto 915
    if (pCS%numStructCtrlParams == 0) return

    !! Initialize the input side

    allocate(lNode_from_SAM_node(numNod), iCin_from_AllVreg(numVreg), &
         &   pCS%structToControlSensors(pCS%numStructCtrlParams), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct')
       return
    end if

    lNode_from_SAM_node = 0
    iCin_from_AllVreg  = 0


    !! SENSOR side initialization

    do i = 1, pCS%numStructCtrlParams
       pS          => pCS%structToControlSensors(i)
       pS%pCtrlPrm => inputs(tmpIdx(i))
       pS%pSensor  => getSensor(pS%pCtrlPrm)

       pS%iCin = 0  !! i.e. not set yet
       pS%whichVreg = pS%pCtrlPrm%var
       iCin_from_AllVreg(pS%whichVreg) = 1 !! flag that this vreg is part of the cIn side

       select case (pS%pSensor%type)
       case (TRIAD_p)
          pTriad1 => triads(pS%pSensor%index)
          lNode_from_SAM_node(pTriad1%samNodNum) = 1

       case (RELATIVE_TRIAD_p)
          pTriad1 => triads(pS%pSensor%index)
          pTriad2 => triads(pS%pSensor%index2)
          lNode_from_SAM_node(pTriad1%samNodNum) = 1
          lNode_from_SAM_node(pTriad2%samNodNum) = 1

       case (JOINT_VARIABLE_p)
          pJoint => joints(pS%pSensor%index)
          lNode_from_SAM_node(pJoint%samNodNum) = 1

       end select

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
    allocate(pCS%whichVregIn(pCS%numVregIn), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct 2')
       return
    end if

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

    call getCtrlOutForces (forces, pCS%numControlForces, tmpIdx)
    if (pCS%numControlForces < 0) goto 915

    allocate(pCS%controlForces(pCS%numControlForces), &
         &   iCout_from_AllVreg(numVreg), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct 3')
       return
    end if

    iCout_from_AllVreg = 0

    do i = 1, pCS%numControlForces
       pF         => pCS%controlForces(i)
       pF%pForce  => forces(tmpIdx(i))

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
    allocate(pCS%whichVregOut(pCS%numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct 3')
       return
    end if

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

    allocate(pCS%samMNPC(nElNodes), pCS%local_MADOF(nElNodes+1), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct 4')
       return
    end if

    pCS%local_MADOF = 0

    do iSam = 1, size(lNode_from_SAM_node)   !! Loop over all sam node nums
       if ( lNode_from_SAM_node(iSam) > 0 ) then
          pCS%samMNPC( lNode_from_SAM_node(iSam) ) = iSam
       end if
    end do


    !! Set the local node association and dofStart for all the forces and sensors

    do i = 1, size(pCS%structToControlSensors)
       pS => pCS%structToControlSensors(i)
       pS%lNode    = 0
       pS%nDOFs    = 0
       pS%dofStart = 0

       select case (pS%pSensor%type)
       case (TRIAD_p)
          pTriad1 => triads(pS%pSensor%index)
          pS%lNode(1) = lNode_from_SAM_node(pTriad1%samNodNum)
          pS%nDOFs(1) = pTriad1%nDOFs

       case (RELATIVE_TRIAD_p)
          pTriad1 => triads(pS%pSensor%index)
          pTriad2 => triads(pS%pSensor%index2)
          pS%lNode(1) = lNode_from_SAM_node(pTriad1%samNodNum)
          pS%nDOFs(1) = pTriad1%nDOFs
          pS%lNode(2) = lNode_from_SAM_node(pTriad2%samNodNum)
          pS%nDOFs(2) = pTriad2%nDOFs

       case (JOINT_VARIABLE_p)
          pJoint => joints(pS%pSensor%index)
          pS%lNode(1) = lNode_from_SAM_node(pJoint%samNodNum)
          pS%nDOFs(1) = pJoint%nJointDOFs
       end select

       do n = 1, 2
          if (pS%lNode(n) > 0) pCS%local_MADOF(pS%lNode(n)) = pS%nDOFs(n)
       end do
    end do

    do i = 1, size(pCS%controlForces)
       pF => pCS%controlForces(i)

       if      (associated(pF%pForce%triad)) then
          pF%lNode = lNode_from_SAM_node(pF%pForce%triad%samNodNum)
          pF%nDOFs = pF%pForce%triad%nDOFs
       else if (associated(pF%pForce%joint)) then
          pF%lNode = lNode_from_SAM_node(pF%pForce%joint%samNodNum)
          pF%nDOFs = pF%pForce%joint%nJointDOFs
       else
          pF%lNode = 0
          pF%nDOFs = 0
       end if

       if (pF%lNode > 0) pCS%local_MADOF(pF%lNode) = pF%nDOFs
    end do

    !! Clean up some scratch space
    deallocate(lNode_from_SAM_node)

    !! Now accumulate the dofStart

    iStart = 1
    do i = 1, size(pCS%local_MADOF)
       n = pCS%local_MADOF(i)
       pCS%local_MADOF(i) = iStart
       iStart = iStart + n
    end do

    !! Number of total dofs for this element
    pCS%nDOFs = pCS%local_MADOF(nElNodes+1)-1

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

    allocate(pCS%Grad_CinWrtSensor(pCS%numVregIn,pCS%nDofs), &
         &   pCS%Grad_CoutWrtCin(pCS%numVregOut,pCS%numVregIn), &
         &   pCS%Grad_ForceWrtCout(pCS%nDofs,pCS%numVregOut), &
         &   pCS%ctrlProps(pCS%numVregOut,4,pCS%numVregIn), &
         &   pCS%massMat(pCS%nDofs,pCS%nDofs), &
         &   pCS%dampMat(pCS%nDofs,pCS%nDofs), &
         &   pCS%stiffMat(pCS%nDofs,pCS%nDofs), &
         &   pCS%SSEEMat(pCS%nDofs,pCS%nDofs), STAT=ierr)
    if (ierr /= 0) ierr = allocationError('InitiateControlStruct 5')

    return

915 ierr = -1
    call reportError (debugFileOnly_p,'InitiateControlStruct')

  end subroutine InitiateControlStruct


  !!============================================================================
  !> @brief Finds all control input parameters coupled to structural DOFs.

  subroutine getCtrlParamsWithStructSensors (ctrlParams, numStructCtrlParams, &
       &                                     ctrlParamsWithStructSensors)

    use ReportErrorModule, only : allocationError, reportError, warning_p

    type(CtrlPrm), target, intent(in)  :: ctrlParams(:)
    integer,               intent(out) :: numStructCtrlParams
    integer,      pointer, intent(out) :: ctrlParamsWithStructSensors(:)

    !! Local variables
    integer :: i, iPrm

    !! --- Logic section ---

    numStructCtrlParams = 0
    do i = 1, size(ctrlParams)
       if (hasStructSensor(ctrlParams(i),.false.)) then
          numStructCtrlParams = numStructCtrlParams + 1
       end if
    end do
    if (numStructCtrlParams == 0) then
       call reportError (warning_p,'No structural inputs in control system.', &
            'Coupled control system modal analysis is therefore switched off.')
       return
    end if

    allocate(ctrlParamsWithStructSensors(numStructCtrlParams),stat=iPrm)
    if (iPrm /= 0) then
       numStructCtrlParams = allocationError('getCtrlParamsWithStructSensors')
       return
    end if

    iPrm = 0
    do i = 1, size(ctrlParams)
       if (hasStructSensor(ctrlParams(i),.true.)) then
          iPrm = iPrm + 1
          ctrlParamsWithStructSensors(iPrm) = i
       end if
    end do

  contains

    !> @brief Checks if a control parameter is coupled to a structural DOF.
    logical function hasStructSensor (prm,notify)
      use ControlTypeModule, only : getSensor
      use SensorTypeModule , only : SensorType, sensorType_p
      use SensorTypeModule , only : TRIAD_P, RELATIVE_TRIAD_p
      use SensorTypeModule , only : JOINT_VARIABLE_p, TIME_p
      type(CtrlPrm), intent(in) :: prm
      logical      , intent(in) :: notify
      type(SensorType), pointer :: sensor
      sensor => getSensor(prm)
      if (associated(sensor)) then
         hasStructSensor = sensor%type == TRIAD_P .or. &
              &            sensor%type == RELATIVE_TRIAD_p .or. &
              &            sensor%type == JOINT_VARIABLE_p
         if (.not.hasStructSensor .and. sensor%type/=TIME_p .and. notify) then
            call reportError (warning_p,'Unsupported sensor type '// &
                 &            trim(sensorType_p(sensor%type))//' (ignored)', &
                 &            addString='getCtrlParamsWithStructSensors')
         end if
      else ! Logic error, should not happen
         hasStructSensor = .false.
      end if
    end function hasStructSensor

  end subroutine getCtrlParamsWithStructSensors


  !!============================================================================
  !> @brief Finds all control out forces.

  subroutine getCtrlOutForces (forces,numCtrlOutForces,ctrlOutForces)

    use ForceTypeModule  , only : ForceType
    use ReportErrorModule, only : allocationError, reportError, warning_p

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
    if (numCtrlOutForces == 0) then
       call reportError (warning_p,'No control output forces in the model')
       return
    end if

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


  !!============================================================================
  !> @brief Computes the gradient matrices of forces w.r.t. controller inputs.

  subroutine BuildStructControlJacobi (pCS,ctrl,sys,msim,ierr)

    use ControlTypeModule   , only : ControlType
    use SystemTypeModule    , only : SystemType
    use ForceRoutinesModule , only : calcExtTriadForceGradient
    use EngineRoutinesModule, only : SensorGradient
    use ReportErrorModule   , only : reportError, debugFileOnly_p

    implicit none

    type(ControlStructType), intent(inout) :: pCS
    type(ControlType)      , intent(inout) :: ctrl
    type(SystemType)       , intent(inout) :: sys
    integer                , intent(in)    :: msim(:)
    integer                , intent(out)   :: ierr

    !! Local variables
    integer  :: i, iNode, iCin, iCout, iStart, nStep
    real(dp) :: sGrad(12), fGrad(6)

    !! --- Logic section ---

    ierr = 0

    !! Compute the force gradients w.r.t control outputs

    pCS%Grad_ForceWrtCout = 0.0_dp

    do i = 1, size(pCS%controlForces)
       call calcExtTriadForceGradient (pCS%controlForces(i)%pForce, fGrad, ierr)
       if (ierr < 0) goto 915

       !! Insert into Grad_ForceWrtCout
       iStart = pCS%controlForces(i)%dofStart
       iCout  = pCS%controlForces(i)%iCout
       call DAXPY (pCS%controlForces(i)%nDofs, 1.0_dp, &
            &      fGrad(1), 1, pCS%Grad_ForceWrtCout(iStart,iCout), 1)
    end do

    !! Compute the control input gradients w.r.t displacement dofs

    pCS%Grad_CinWrtSensor = 0.0_dp

    do i = 1, size(pCS%structToControlSensors)
       call SensorGradient (pCS%structToControlSensors(i)%pCtrlPrm, sGrad, ierr)
       if (ierr < 0) goto 915

       do iNode = 1, 2
          if (pCS%structToControlSensors(i)%lNode(iNode) > 0) then
             !! Insert into Grad_CinWrtSensor
             iStart = pCS%structToControlSensors(i)%dofStart(iNode)
             iCin   = pCS%structToControlSensors(i)%iCin
             call DAXPY (pCS%structToControlSensors(i)%nDofs(iNode), 1.0_dp, &
                  &      sGrad(iNode*6-5), 1, &
                  &      pCS%Grad_CinWrtSensor(iCin,iStart), pCS%numVregIn)
          end if
       end do
    end do


    !! TODO, Magne: Add controller gradients here and build the full gradients
    pCS%ctrlProps = 0.0_dp

    select case (pCS%ctrlSysEigFlag)
    case (1) ! nPertub = 3: P, I and D gains
       call EstimateControllerProperties01 (sys, ctrl, msim, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               pCS%ctrlProps, ierr)

    case (2) ! nPerturb = 4
       call EstimateControllerProperties02 (sys, ctrl, msim, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               pCS%ctrlProps, ierr)

    case (3) ! nPerturb = 6
       call EstimateControllerProperties03 (sys, ctrl, msim, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               pCS%ctrlProps, ierr)

    case (4:8) ! nPerturb = 5
       nStep = pCS%ctrlSysEigFlag-3 ! nStep = 1...5
       call EstimateControllerProperties04 (sys, ctrl, msim, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               nStep, pCS%ctrlProps, ierr)

    case (9) ! nPerturb = 5
       nStep = 10
       call EstimateControllerProperties04 (sys, ctrl, msim, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               nStep, pCS%ctrlProps, ierr)

    case (10) ! nPerturb = 5
       nStep = 100
       call EstimateControllerProperties04 (sys, ctrl, msim ,&
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               nStep, pCS%ctrlProps, ierr)

    case (11) ! nPerturb = 5
       nStep = 1000
       call EstimateControllerProperties04 (sys, ctrl, msim, &
            &                               pCS%numVregIn, pCS%whichVregIn, &
            &                               pCS%numVregOut, pCS%whichVregOut, &
            &                               nStep, pCS%ctrlProps, ierr)

    case (500) ! nPerturb = 1
       call EstimateControllerProperties500 (sys, ctrl, msim, &
            &                                pCS%numVregIn, pCS%whichVregIn, &
            &                                pCS%numVregOut, pCS%whichVregOut, &
            &                                pCS%ctrlProps, ierr)

    case default ! Error
       if (pCS%ctrlSysEigFlag > 0) then
          ierr = -pCS%ctrlSysEigFlag
       else
          ierr = -999
       end if
    end select
    if (ierr < 0) goto 915

    !! Steady-state error elimination matrix (Q)
    call dMMM (pCS%SSEEMat,pCS%ctrlProps(:,1,:),ierr)
    if (ierr < 0) goto 915

    !! Stiffness matrix (K)
    call dMMM (pCS%stiffMat,pCS%ctrlProps(:,2,:),ierr)
    if (ierr < 0) goto 915

    !! Damping matrix (C)
    call dMMM (pCS%dampMat,pCS%ctrlProps(:,3,:),ierr)
    if (ierr < 0) goto 915

    !! Mass matrix (M)
    call dMMM (pCS%massMat,pCS%ctrlProps(:,4,:),ierr)
    if (ierr < 0) goto 915

    !! Check for symmetry
    pCS%SSEEMatIsNonSymmetric  = isNonSymmetric(pCS%SSEEMat)
    pCS%stiffMatIsNonSymmetric = isNonSymmetric(pCS%stiffMat)
    pCS%dampMatIsNonSymmetric  = isNonSymmetric(pCS%dampMat)
    pCS%massMatIsNonSymmetric  = isNonSymmetric(pCS%massMat)

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

    return

915 call reportError (debugFileOnly_p,'BuildStructControlJacobi')

  contains

    !> @brief Performs the matrix-matrix multiplication Cmat = A*ctrlProp*B.
    !> @details A = pCS%Grad_ForceWrtCout and B = pCS%Grad_CinWrtSensor.
    !> Using LAPACK::DGEMM('N','N',m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC)
    subroutine dmmm (Cmat,ctrlProp,ierr)
      use scratchArrayModule, only : getRealScratchArray
      real(dp), intent(out) :: Cmat(:,:)
      real(dp), intent(in)  :: ctrlProp(:,:)
      integer , intent(out) :: ierr
      real(dp), pointer     :: rWork(:)
      rWork => getRealScratchArray(pCS%nDofs*pCS%numVregIn,ierr)
      if (ierr == 0) then
         !! rWork(nDofs,nVregIn) = Grad_ForceWrtCout(nDofs,nVregOut)
         !!                      * ctrlProp(nVregOut,nVregIn)
         call DGEMM ('N','N', pCS%nDofs, pCS%numVregIn, pCS%numVregOut, &
              &       1.0_dp, pCS%Grad_ForceWrtCout(1,1), pCS%nDofs, &
              &               ctrlProp(1,1), pCS%numVregOut, &
              &       0.0_dp, rWork(1), pCS%nDofs)
         !! Cmat(nDofs,nDofs) = -rWork(nDofs,nVregIn)
         !!                   * Grad_CinWrtSensor(nVregIn,nDofs)
         call DGEMM ('N','N', pCS%nDofs,  pCS%nDofs, pCS%numVregIn,&
              &      -1.0_dp, rWork(1), pCS%nDofs, &
              &               pCS%Grad_CinWrtSensor(1,1), pCS%numVregIn, &
              &       0.0_dp, Cmat(1,1), pCS%nDofs)
      end if
    end subroutine dmmm

    !> @brief Checks if the matrix @b A has non-symmetric terms.
    function isNonSymmetric (A)
      real(dp), intent(in) :: A(:,:)
      integer :: i, j, isNonSymmetric
      isNonSymmetric = 0
      do i = 1, size(A,1)-1
         do j = i+1, size(A,2)
            if (abs(A(i,j)-A(j,i)) > 1.0e-15_dp) then
               isNonSymmetric = 1
               return
            end if
         end do
      end do
    end function isNonSymmetric

  end subroutine BuildStructControlJacobi


  subroutine EstimateControllerProperties01 (sys, ctrl, msim, &
       &                                     numVregIn, whichVregIn,   &
       &                                     numVregOut, whichVregOut, &
       &                                     ctrlProps, ierr)

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

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use DenseMatrixModule, only : solveAxB
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    implicit none

    type(SystemType)   , intent(inout) :: sys
    type(ControlType)  , intent(inout) :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: numVregIn         !! Number of vreg in to perturb
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: numVregOut        !! Number of vreg out to read variation from
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
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

    integer, parameter :: nPerturb = 3 !< Number of controller perturbations

    type(ControlType), pointer :: ctrlCopy => null()
    type(ControlType), pointer :: ctrlCopyCopy => null()

    real(dp) :: orgTime                          !! original/initial value for the time
    real(dp) :: orgTimeStep                      !! original/initial value for the time step
    real(dp) :: dt0                              !! quazi-initial time step
    real(dp) :: dt                               !! incremental time step
    real(dp) :: dy                               !! incremental controller input step
    real(dp) :: dintydt                          !! definite integral of y with respect to time
    real(dp), allocatable :: y0(:)               !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)              !! u(y0), output from controller when input is y0
    real(dp), allocatable :: uy(:)               !! u(y), output from controller when input is y = y0+dy
    real(dp) :: dyMatrix(nPerturb,nPerturb)      !! matrix of perturbation parameters
    real(dp), allocatable :: fullCtrlProps(:,:)  !! Matrix of unnecessarily many controller properties

    integer :: i, j, iInput

    !! --- Logic section ---

    allocate(y0(numVregIn),uy0(numVregOut),uy(numVregOut), &
         &   fullCtrlProps(nPerturb,numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties01')
       return
    end if

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    dt0 = orgTimeStep*0.1_dp   !! TODO Magne: Change this value?
    sys%time = orgTime + dt0

    do i = 1, numVregIn
       iInput = whichVregIn(i)

       call copyCtrl (ctrl,ctrlCopy,ierr)
       if (ierr < 0) goto 915

       dy = 0.0_dp ! Initial perturbation with dy = 0
       call PerturbController (sys,ctrlCopy,msim,iInput,dt0,dy, &
            &                  numVregOut,whichVregOut,uy0,ierr)
       if (ierr < 0) goto 915

       !! Save the value of y0
       y0(i) = abs(ctrlCopy%vreg(iInput))

       !! Perturbation
       do j = 1, nPerturb
          !! Establish dy-matrix
          dt = orgTimeStep*0.1_dp*j
          dy = dt
          dintydt = (y0(i) + 0.5_dp*dy)*dt

          dyMatrix(j,1) = dintydt
          dyMatrix(j,2) = dy
          dyMatrix(j,3) = dy/dt

          !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
          call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
          if (ierr < 0) goto 915

          !! Perturb
          call PerturbController (sys,ctrlCopyCopy,msim,iInput,dt,dy, &
               &                  numVregOut,whichVregOut,uy,ierr)
          if (ierr < 0) goto 915

          !! Calculate du and store in a table
          fullCtrlProps(j,:) = uy - uy0
       end do

       !! Calculate controller properties
       call solveAxB (dyMatrix,fullCtrlProps,ierr)
       if (ierr < 0) goto 915

       select case (ctrlCopy%input(iInput)%engine%args(1)%p%entity)
       case (POS_p) ! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = fullCtrlProps(1,:)
          ctrlProps(:,2,i) = fullCtrlProps(2,:)
          ctrlProps(:,3,i) = fullCtrlProps(3,:)
          ctrlProps(:,4,i) = 0.0_dp
       case (VEL_p) ! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = 0.0_dp
          ctrlProps(:,2,i) = fullCtrlProps(1,:)
          ctrlProps(:,3,i) = fullCtrlProps(2,:)
          ctrlProps(:,4,i) = fullCtrlProps(3,:)
       case (ACC_p) ! find dintydt(j), dy(j)
          ctrlProps(:,1,i) = 0.0_dp
          ctrlProps(:,2,i) = 0.0_dp
          ctrlProps(:,3,i) = fullCtrlProps(1,:)
          ctrlProps(:,4,i) = fullCtrlProps(2,:)
       case default
          !! Error
       end select
    end do

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    call deallocateCtrlCopy (ctrlCopyCopy)
    deallocate(y0,uy0,uy,fullCtrlProps)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties01')
    goto 900

  end subroutine EstimateControllerProperties01


  subroutine EstimateControllerProperties02 (sys, ctrl, msim, &
       &                                     numVregIn, whichVregIn,   &
       &                                     numVregOut, whichVregOut, &
       &                                     ctrlProps, ierr)

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

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use DenseMatrixModule, only : solveAxB
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    implicit none

    type(SystemType)   , intent(inout) :: sys
    type(ControlType)  , intent(inout) :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: numVregIn         !! Number of vreg in to perturb
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: numVregOut        !! Number of vreg out to read variation from
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
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

    integer, parameter :: nPerturb = 4 !< Number of controller perturbations

    type(ControlType), pointer :: ctrlCopy => null()
    type(ControlType), pointer :: ctrlCopyCopy => null()

    real(dp) :: orgTime                         !! original/initial value for the time
    real(dp) :: orgTimeStep                     !! original/initial value for the time step
    real(dp) :: dt0                             !! quazi-initial time step
    real(dp) :: dy0                             !! quazi-initial input perturbation, dy0 = 0
    real(dp) :: dt              !! incremental time step
    real(dp) :: dy              !! incremental controller input step
    real(dp) :: dy1             !! incremental controller input step no. 1
    !                                           !! to be used when deriving d(d2y/dt2)
    real(dp) :: dy2             !! incremental controller input step no. 2
    !                                           !! to be used when deriving d(d2y/dt2)
    real(dp) :: y1              !! y1 = y0+dy1
    real(dp) :: ddydt           !! d(dy/dt), 1st derivative of y with respect to time
    real(dp) :: dd2ydt2         !! d(d2y/dt2), 2nd derivative of y with respect to time
    real(dp) :: dintydt         !! definite integral of y with respect to time
    real(dp) :: dintintydt      !! definite double integral of y with respect to time
    real(dp) :: dintintintydt   !! definite triple integral of y with respect to time
    real(dp), allocatable :: y0(:)              !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)             !! u(y0), output from controller when input is y0
    real(dp), allocatable :: uy1(:)             !! u(y1)
    real(dp), allocatable :: uy(:)              !! u(y), output from controller when input is y = y0+dy
    real(dp) :: dyMatrix(nPerturb,nPerturb)     !! matrix of perturbation parameters
    real(dp), allocatable :: fullCtrlProps(:,:) !! Matrix of unnecessarily many controller properties

    integer :: i, j

    !! --- Logic section ---

    allocate(y0(numVregIn),uy0(numVregOut),uy1(numVregOut),uy(numVregOut), &
         &   fullCtrlProps(nPerturb,numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties02')
       return
    end if

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Make copy of the controller
    call copyCtrl (ctrl,ctrlCopy,ierr)
    if (ierr < 0) goto 915

    !! Do one time perturbation
    dt0 = sys%timeStep*0.1_dp   !! TODO Magne: Change this value?
    dy0 = 0.0_dp
    sys%time = sys%time + dt0

    !! Start perturbation with time step dt0
    do i = 1, numVregIn
       call PerturbController (sys,ctrlCopy,msim,whichVregIn(i),dt0,dy0, &
            &                  numVregOut,whichVregOut,uy,ierr)
       if (ierr < 0) goto 915
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
       select case (ctrlCopy%input(i)%engine%args(1)%p%entity)
       case (POS_p)
          do j = 1, nPerturb
             !! Establish dy-matrix
             dt  = orgTimeStep*0.1_dp*j
             dy1 = dt
             dy2 = -dy1
             y1  = y0(i) + dy1
             dintydt = (y1  + 0.5*dp*dy2)*dt
             ddydt   = (dy2 -        dy1)/dt
             dd2ydt2 = (dy2 - 2.0_dp*dy1)/(dt*dt)

             dyMatrix(j,1) = dintydt
             dyMatrix(j,2) = dy2
             dyMatrix(j,3) = ddydt
             dyMatrix(j,4) = dd2ydt2

             !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
             call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
             if (ierr < 0) goto 915

             !! To derive d(d2y/dt2), the system has to be perturbed three times (two + initial)

             !! First perturbation
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy1, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915

             !! Save the value of u(y1) in an array
             uy1 = uy

             !! Second perturbation
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy2, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915

             !! Calculate du and store in a table
             fullCtrlProps(j,:) = uy - uy1
          end do

       case (VEL_p)
          !! Establish dy-matrix
          do j = 1, nPerturb
             dt = orgTimeStep*0.1_dp*j
             dy = dt
             dintydt    = (y0(i)        + dy/2.0_dp)*dt
             dintintydt = (y0(i)/2.0_dp + dy/6.0_dp)*dt*dt

             dyMatrix(j,1) = dintintydt
             dyMatrix(j,2) = dintydt
             dyMatrix(j,3) = dy
             dyMatrix(j,4) = dy/dt

             !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
             call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
             if (ierr < 0) goto 915

             !! Perturb
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915

             !! Calculate du and store in a table
             fullCtrlProps(j,:) = uy - uy1
          end do

       case (ACC_p)
          do j = 1, nPerturb
             !! Establish dy-matrix
             dt = orgTimeStep*0.1_dp*j
             dy = dt
             dintydt       = (y0(i)        + dy/2.0_dp )*dt
             dintintydt    = (y0(i)/2.0_dp + dy/6.0_dp )*dt*dt
             dintintintydt = (y0(i)/6.0_dp + dy/24.0_dp)*dt*dt*dt

             dyMatrix(j,1) = dintintintydt
             dyMatrix(j,2) = dintintydt
             dyMatrix(j,3) = dintydt
             dyMatrix(j,4) = dy

             !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
             call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
             if (ierr < 0) goto 915

             !! Perturb
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915

             !! Calculate du and store in a table
             fullCtrlProps(j,:) = uy - uy1
          end do

       case default
          !! Error
          goto 915
       end select

       !! Calculate controller properties
       call solveAxB (dyMatrix,fullCtrlProps,ierr)
       if (ierr < 0) goto 915

       ctrlProps(:,:,i) = transpose(fullCtrlProps)
    end do

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    call deallocateCtrlCopy (ctrlCopyCopy)
    deallocate(y0,uy0,uy1,uy,fullCtrlProps)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties02')
    goto 900

  end subroutine EstimateControllerProperties02


  subroutine EstimateControllerProperties03 (sys, ctrl, msim, &
       &                                     numVregIn, whichVregIn, &
       &                                     numVregOut, whichVregOut, &
       &                                     ctrlProps, ierr)

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

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use DenseMatrixModule, only : solveAxB
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    implicit none

    type(SystemType)   , intent(inout) :: sys
    type(ControlType)  , intent(inout) :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: numVregIn         !! Number of vreg in to perturb
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: numVregOut        !! Number of vreg out to read variation from
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
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

    integer, parameter :: nPerturb = 6 !< Number of controller perturbations

    type(ControlType), pointer :: ctrlCopy => null()
    type(ControlType), pointer :: ctrlCopyCopy => null()

    real(dp) :: orgTime                          !! original/initial value for the time
    real(dp) :: orgTimeStep                      !! original/initial value for the time step
    real(dp) :: dt0                              !! quazi-initial time step
    real(dp) :: dy0                              !! quazi-initial input perturbation, dy0 = 0
    real(dp) :: dt               !! incremental time step
    real(dp) :: dy1              !! incremental controller input step no. 1
    !                                            !! to be used when deriving d(d2y/dt2)
    real(dp) :: dy2              !! incremental controller input step no. 2
    !                                            !! to be used when deriving d(d2y/dt2)
    real(dp) :: y1               !! y1 = y0+dy1
    real(dp) :: ddydt            !! d(dy/dt), 1st derivative of y
    !										     !! with respect to time t
    real(dp) :: dd2ydt2          !! d(d2y/dt2), 2nd derivative of y
    !										     !! with respect to time t
    real(dp) :: dintydt          !! definite integral of y with respect to time t
    real(dp) :: dintintydt       !! definite double integral of y with respect to time t
    real(dp) :: dintintintydt    !! definite triple integral of y with respect to time t
    real(dp), allocatable :: y0(:)               !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)              !! u(y0), output from controller when input is y0
    real(dp), allocatable :: uy1(:)              !! u(y1)
    real(dp), allocatable :: uy2(:)              !! u(y2), output from controller when input is y = y0+dy
    real(dp) :: dyMatrix(nPerturb,nPerturb)      !! matrix of perturbation parameters
    real(dp), allocatable :: fullCtrlProps(:,:)  !! Matrix of unnecessarily many controller properties

    integer :: i, j

    !! --- Logic section ---

    allocate(y0(numVregIn),uy0(numVregOut),uy1(numVregOut),uy2(numVregOut), &
         &   fullCtrlProps(nPerturb,numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties03')
       return
    end if

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Make copy of controller
    call copyCtrl (ctrl,ctrlCopy,ierr)
    if (ierr < 0) goto 915

    !! Do one time perturbation
    dt0 = sys%timeStep*0.1_dp   !! TODO Magne: Change this value?
    dy0 = 0.0_dp
    sys%time = sys%time + dt0

    !! Start perturbation with time step dt0
    do i = 1, numVregIn
       call PerturbController (sys,ctrlCopy,msim,whichVregIn(i),dt0,dy0, &
            &                  numVregOut,whichVregOut,uy0,ierr)
       if (ierr < 0) goto 915
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
          dt  = orgTimeStep*0.1_dp*j
          dy1 = dt
          dy2 = -dy1
          y1  = y0(i) + dy1
          dintydt       = (y1        + dy2/2.0_dp )*dt
          dintintydt    = (y1/2.0_dp + dy2/6.0_dp )*dt*dt
          dintintintydt = (y1/6.0_dp + dy2/24.0_dp)*dt*dt*dt
          ddydt   = (dy2 -        dy1)/dt
          dd2ydt2 = (dy2 - 2.0_dp*dy1)/(dt*dt)

          dyMatrix(j,1) = dintintintydt
          dyMatrix(j,2) = dintintydt
          dyMatrix(j,3) = dintydt
          dyMatrix(j,4) = dy2
          dyMatrix(j,5) = ddydt
          dyMatrix(j,6) = dd2ydt2

          !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
          call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
          if (ierr < 0) goto 915

          !! To derive d(d2y/dt2), the system has to be perturbed three times (two + initial)

          !! First perturbation
          call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy1, &
               &                  numVregOut,whichVregOut,uy1,ierr)
          if (ierr < 0) goto 915

          !! Second perturbation
          call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy2, &
               &                  numVregOut,whichVregOut,uy2,ierr)
          if (ierr < 0) goto 915

          !! Calculate du and store in a table
          fullCtrlProps(j,:) = uy2 - uy1
       end do

       !! Calculate controller properties
       call solveAxB (dyMatrix,fullCtrlProps,ierr)
       if (ierr < 0) goto 915

       select case (ctrlCopy%input(i)%engine%args(1)%p%entity)
       case (POS_p) ! dintydt, dy, ddydt and dd2ydt2
          ctrlProps(:,1,i) = fullCtrlProps(3,:)
          ctrlProps(:,2,i) = fullCtrlProps(4,:)
          ctrlProps(:,3,i) = fullCtrlProps(5,:)
          ctrlProps(:,4,i) = fullCtrlProps(6,:)
       case (VEL_p) ! dintintydt, dintydt, dy, ddydt
          ctrlProps(:,1,i) = fullCtrlProps(2,:)
          ctrlProps(:,2,i) = fullCtrlProps(3,:)
          ctrlProps(:,3,i) = fullCtrlProps(4,:)
          ctrlProps(:,4,i) = fullCtrlProps(5,:)
       case (ACC_p) ! dintintintydt, dintintydt, dintydt, dy
          ctrlProps(:,1,i) = fullCtrlProps(1,:)
          ctrlProps(:,2,i) = fullCtrlProps(2,:)
          ctrlProps(:,3,i) = fullCtrlProps(3,:)
          ctrlProps(:,4,i) = fullCtrlProps(4,:)
       case default
          !! Error
       end select

    end do

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    call deallocateCtrlCopy (ctrlCopyCopy)
    deallocate(y0,uy0,uy1,uy2,fullCtrlProps)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties03')
    goto 900

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

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use DenseMatrixModule, only : solveAxB
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    implicit none

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

    type(ControlType), pointer :: ctrlCopy => null()
    type(ControlType), pointer :: ctrlCopyCopy => null()

    real(dp) :: orgTime                          !! original/initial value for the time
    real(dp) :: orgTimeStep                      !! original/initial value for the time step
    real(dp) :: dt0                              !! quazi-initial time step
    real(dp) :: dy0                              !! quazi-initial input perturbation, dy0 = 0
    real(dp) :: ddt                              !! incremental steps of dt
    real(dp) :: ddy                              !! incremental steps of dy
    real(dp) :: dt               !! incremental time step
    real(dp) :: dy               !! incremental controller input step
    real(dp) :: dintydt          !! definite integral of y with respect to time
    real(dp) :: dintintydt       !! definite double integral of y with respect to time
    real(dp) :: dintintintydt    !! definite triple integral of y with respect to time
    real(dp), allocatable :: y0(:)               !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)              !! u(y0), output from controller when input is y0
    real(dp), allocatable :: uy(:)               !! u(y), output from controller when input is y = y0+dy
    real(dp) :: dyMatrix(nPerturb,nPerturb)      !! matrix of perturbation parameters
    real(dp), allocatable :: fullCtrlProps(:,:)  !! Matrix of unnecessarily many controller properties

    integer :: i, ii, j

    !! --- Logic section ---

    allocate(y0(numVregIn),uy0(numVregOut),uy(numVregOut), &
         &   fullCtrlProps(nPerturb,numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties04')
       return
    end if

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Make copy of controller
    call copyCtrl (ctrl,ctrlCopy,ierr)
    if (ierr < 0) goto 915

    !! Do one initial time perturbation
    dt0 = sys%timeStep*0.1_dp   !! TODO Magne: Change this value?
    dy0 = 0.0_dp
    sys%time = sys%time + dt0

    !! Start perturbation with time step dt0
    do i = 1, numVregIn
       !! Perturb
       ddt = dt0/real(nStep,dp)
       do ii = 1, nStep
          call PerturbController (sys,ctrlCopy,msim,whichVregIn(i),ddt,dy0, &
               &                  numVregOut,whichVregOut,uy,ierr)
          if (ierr < 0) goto 915
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
          dt            = orgTimeStep*0.1_dp*j
          dy            = dt
          dintydt       = (y0(i)/1.0_dp + dy/2.0_dp )*dt
          dintintydt    = (y0(i)/2.0_dp + dy/6.0_dp )*dt*dt
          dintintintydt = (y0(i)/6.0_dp + dy/24.0_dp)*dt*dt*dt

          dyMatrix(j,1) = dintintintydt
          dyMatrix(j,2) = dintintydt
          dyMatrix(j,3) = dintydt
          dyMatrix(j,4) = dy
          dyMatrix(j,5) = dy/dt

          !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
          call copyCtrl (ctrlCopy,ctrlCopyCopy)
          if (ierr < 0) goto 915

          !! First perturbation
          ddt = dt/real(nStep,dp)
          ddy = dy/real(nStep,dp)
          do ii = 1, nStep
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),ddt,ddy, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915
          end do

          !! Calculate and store du in a table
          fullCtrlProps(j,:) = uy - uy0
       end do

       !! Calculate controller properties
       call solveAxB (dyMatrix,fullCtrlProps,ierr)
       if (ierr < 0) goto 915

       select case (ctrlCopy%input(i)%engine%args(1)%p%entity)
       case (POS_p) ! find dintydt, dy, ddydt
          ctrlProps(:,1,i) = fullCtrlProps(3,:)
          ctrlProps(:,2,i) = fullCtrlProps(4,:)
          ctrlProps(:,3,i) = fullCtrlProps(5,:)
          ctrlProps(:,4,i) = 0.0_dp
       case (VEL_p) ! find dintintydt, dintydt, dy, ddydt
          ctrlProps(:,1,i) = fullCtrlProps(2,:)
          ctrlProps(:,2,i) = fullCtrlProps(3,:)
          ctrlProps(:,3,i) = fullCtrlProps(4,:)
          ctrlProps(:,4,i) = fullCtrlProps(5,:)
       case (ACC_p) ! find dintintintydt,dintintydt,dintydt, dy
          ctrlProps(:,1,i) = fullCtrlProps(1,:)
          ctrlProps(:,2,i) = fullCtrlProps(2,:)
          ctrlProps(:,3,i) = fullCtrlProps(3,:)
          ctrlProps(:,4,i) = fullCtrlProps(4,:)
       case default
          !! Error
       end select

    end do

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    call deallocateCtrlCopy (ctrlCopyCopy)
    deallocate(y0,uy0,uy,fullCtrlProps)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties04')
    goto 900

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
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    implicit none

    type(SystemType)   , intent(inout) :: sys
    type(ControlType)  , intent(inout) :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: numVregIn         !! Number of vreg in to perturb
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: numVregOut        !! Number of vreg out to read variation from
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
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

    type(ControlType), pointer :: ctrlCopy => null()

    real(dp) :: orgTime                          !! original/initial value for the time
    real(dp) :: orgTimeStep                      !! original/initial value for the time step
    real(dp) :: dt0                              !! quazi-initial time step
    real(dp) :: dy0                              !! quazi-initial input perturbation, dy0 = 0
    real(dp) :: dt               !! incremental time step
    real(dp) :: dy               !! incremental controller input step
    real(dp) :: dyMatrix         !! matrix of perturbation parameters (here just a 1x1 matrix)
    real(dp), allocatable :: y0(:)               !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)              !! u(y0), output from controller when input is y0
    real(dp), allocatable :: uy(:)               !! u(y), output from controller when input is y = y0+dy

    integer :: i

    !! --- Logic section ---

    allocate(y0(numVregIn),uy0(numVregOut),uy(numVregOut),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties500')
       return
    end if

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Do one initial time perturbation
    dt0 = sys%timeStep   !! TODO Magne: Change this value?
    dy0 = 0.0_dp
    sys%time = sys%time + dt0

    !! Save the value of y0 in an array
    do i = 1, numVregIn
       y0(i) = ctrl%vreg(whichVregIn(i))
    end do

    !! Save the value of u(y0) in an array
    do i = 1, numVregOut
       uy0(i) = ctrl%vreg(whichVregOut(i))
    end do

    !! Perturbation
    do i = 1, numVregIn
       !! Establish dy-matrix
       dt = orgTimeStep
       dy = dt

       dyMatrix = (y0(i)+dy/2.0_dp)*dt

       !! Reset current controller (ctrlCopy) to original state (ctrl)
       call copyCtrl (ctrl,ctrlCopy,ierr)
       if (ierr < 0) goto 915

       !! Perturb
       call PerturbController (sys,ctrlCopy,msim,whichVregIn(i),dt,dy, &
            &                  numVregOut,whichVregOut,uy,ierr)
       if (ierr < 0) goto 915

       !! Calculate controller properties
       uy = (1.0_dp/dyMatrix) * (uy - uy0)

       select case (ctrl%input(i)%engine%args(1)%p%entity)
       case (POS_p) ! find dintydt, dy, ddydt
          ctrlProps(:,1,i) = uy
          ctrlProps(:,2,i) = 0.0_dp
          ctrlProps(:,3,i) = 0.0_dp
          ctrlProps(:,4,i) = 0.0_dp
       case (VEL_p) ! find dintydt, dy, ddydt
          ctrlProps(:,1,i) = 0.0_dp
          ctrlProps(:,2,i) = uy
          ctrlProps(:,3,i) = 0.0_dp
          ctrlProps(:,4,i) = 0.0_dp
       case (ACC_p) ! find dintydt, dy
          ctrlProps(:,1,i) = 0.0_dp
          ctrlProps(:,2,i) = 0.0_dp
          ctrlProps(:,3,i) = uy
          ctrlProps(:,4,i) = 0.0_dp
       case default
          !! Error
       end select

    end do

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    deallocate(y0,uy0,uy)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties500')
    goto 900

  end subroutine EstimateControllerProperties500


  !!============================================================================
  !> @brief Perturbs one of the inputs of the control system.
  !>
  !> @param sys System level model data
  !> @param ctrl Control system data
  !> @param[in] msim Matrix of simulation parameters
  !> @param[in] iPert Which input to perturb
  !> @param[in] dt Incremental time step
  !> @param[in] dy Incremental step for use in du/dy
  !> @param[in] numVregOut Number of outputs from the controller to read
  !> @param[in] whichVregOut Which outputs from the controller to read
  !> @param[out] uy Perturbed output from controller, u(y) = u(y0+dy)
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine perturbs one of the control inputs and calculates
  !> the reaction in all of the control outputs.
  !>
  !> Working order:
  !> 1) Change the input y for the controller from y0 to y = y0+dy,
  !>    where dy is a small number. Change the time from t0 to t = t0+dt,
  !>    where dt is a small number.
  !> 2) Run this new input through the controller to get the reaction
  !>    (output) u(y) from the controller due to the change in the input.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Magne Bratland
  !>
  !> @date 18 Sept 2009

  subroutine PerturbController (sys,ctrl,msim,iPert,dt,dy, &
       &                        numVregOut,whichVregOut,uy,ierr)

    use SystemTypeModule     , only : SystemType, dp
    use ControlTypeModule    , only : ControlType
    use ControlRoutinesModule, only : IterateControlSystem

    implicit none

    type(SystemType) , intent(inout) :: sys
    type(ControlType), intent(inout) :: ctrl
    integer          , intent(in)    :: msim(:), iPert
    real(dp)         , intent(in)    :: dt, dy
    integer          , intent(in)    :: numVregOut, whichVregOut(:)
    real(dp)         , intent(out)   :: uy(:)
    integer          , intent(out)   :: ierr

    !! Local variables
    integer, parameter :: ctrlSysMode = 3 !< controller integration mode

    !! --- Logic section ---

    !! Change time step to dt
    sys%timeStep = dt

    !! Change input from y0 to y (y = y0+dy)
    ctrl%vreg(ctrl%input(iPert)%var) = ctrl%vreg(ctrl%input(iPert)%var) + dy

    !! Iterate the controller with the new values for y (y = y0+dy)
    !! and t (t = t0+dt), to derive the value of u(y)
    call IterateControlSystem (sys,ctrl,ctrlSysMode,msim,ierr,setInput=.false.)

    !! Save the value of u(y) in an array
    call DGATHR (numVregOut,whichVregOut(1),ctrl%vreg(1),uy(1),1)

  end subroutine PerturbController


  !!============================================================================
  !> @brief Deallocates a control system copy.
  !>
  !> @param ctrl Pointer to controltypemodule::controltype object to deallocate.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Jan 2024

  subroutine deallocateCtrlCopy (ctrl)

    use ControlTypeModule, only : ControlType, deallocateCtrl

    type(ControlType), pointer :: ctrl

    !! --- Logic section ---

    if (associated(ctrl)) then
       call deallocateCtrl (ctrl,.false.)
       deallocate(ctrl)
    end if

  end subroutine deallocateCtrlCopy


  subroutine setFixedControlDOFsOnEigValCalc(forces)

    !!==========================================================================
    !! Purpose:
    !! Set DOFs for triades, where a force, whose source is a control system, is acting,
    !! to fixed during eigenvalue calculations and free else.
    !!
    !! Working order:
    !! 1) Get all forces in an array
    !! 2) Search through that array and find all forces whose source (i.e. engine)
    !! 	  is a control system
    !! 3) For all forces whose source is a control system:
    !!  3.1) Find out the direction of the force
    !!  3.2) Set value for BCs for the triad, where the force is working, to 2
    !!
    !! Programmer : Magne Bratland
    !! date/rev   : 19 Feb 2009 / 1.0
    !!==========================================================================

    use ForceTypeModule   , only : ForceType
    use FunctionTypeModule, only : EngineType
    use SensorTypeModule  , only : SensorType, CONTROL_p
    use ReportErrorModule , only : reportError, note_p

    implicit none

    type(ForceType), intent(inout), target :: forces(:)

    type(EngineType), pointer :: enginePointer
    type(SensorType), pointer :: argSensorPointer

    integer :: i

    !! Search through array called forces and find all forces whose source is a control system
    do i = 1, size(forces)

       enginePointer => forces(i)%engine
       if (.not. associated(enginePointer)) cycle
       if (size(enginePointer%args) < 1) cycle

       argSensorPointer => enginePointer%args(1)%p
       if (.not. associated(argSensorPointer)) cycle

       if (.not. argSensorPointer%type == CONTROL_p) cycle

       if (associated(forces(i)%triad)) then
          !! Find out in which direction the force works
          select case (forces(i)%dof)
          case (-2) ! The force is a multi-dimensional moment
             forces(i)%triad%BC(4:6) = 2
             call reportError (note_p, &
                  & 'MAGNE: Actually setting fixed BC for control eig dof, all rotation')

          case (-1) ! The force is a multi-dimensional force
             forces(i)%triad%BC(1:3) = 2
             call reportError (note_p, &
                  & 'MAGNE: Actually setting fixed BC for control eig dof, all translation')

          case (1:6) ! The direction of the force is a pure one-dimensional force
             forces(i)%triad%BC(forces(i)%dof) = 2
             call reportError (note_p, &
                  & 'MAGNE: Actually setting fixed BC for control eig dof, single dof')
          case default
             call reportError (note_p,'MAGNE: Actually NOT setting fixed BC for control eig dof')
          end select
       else if (associated(forces(i)%joint)) then
          !! Set joint dof to zero
       else
          !! Error condition ?
       end if

    end do

  end subroutine setFixedControlDOFsOnEigValCalc

end module ControlStructModule
