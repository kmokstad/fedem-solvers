!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine LoadStrainRosetteFromCore (supId,iros,data,ndata,stat)
  !DEC$ ATTRIBUTES DLLEXPORT :: LoadStrainRosetteFromCore

  integer, intent(in)  :: supId, iros, ndata
  real   , intent(out) :: data(ndata)
  integer, intent(out) :: stat

  !! Dummy implementation, a real model-specific instance needs to be generated
  !! using the subroutine stressRecoveryModule::writeRosettes2Ftn
  if (supId < 0) then
     print *,' ** LoadStrainRosetteFromCore dummy: ',supId,iros,ndata
  end if
  data(1) = 0.0
  stat = 0

end subroutine LoadStrainRosetteFromCore
