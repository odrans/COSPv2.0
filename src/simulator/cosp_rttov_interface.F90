! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! May 2015 - D. Swales - Original version
! Apr 2015 - D. Swales - Modified for RTTOVv11.3
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_RTTOV_INTERFACE
  USE COSP_KINDS,       ONLY: wp
  USE MOD_COSP_CONFIG,  ONLY: RTTOV_MAX_CHANNELS,rttovDir
  use mod_cosp_rttov,   only: platform,satellite,sensor,nChannels,iChannel,  &
       opts, errorstatus_success, rttov_exit, coefs, &
       emis_atlas, brdf_atlas, atlas_type, dosolar, nthreads

  USE rttov_const, ONLY :     &
         surftype_sea,        &
         surftype_land,       &
         sensor_id_mw,        &
         sensor_id_po,        &
         inst_name,           &
         platform_name

  ! The rttov_emis_atlas_data type must be imported separately
  USE mod_rttov_emis_atlas, ONLY : &
       rttov_emis_atlas_data, &
       atlas_type_ir, atlas_type_mw
  
  ! The rttov_brdf_atlas_data type must be imported separately
  USE mod_rttov_brdf_atlas, ONLY : rttov_brdf_atlas_data

  
  
  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

! Use emissivity atlas
#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

! Use BRDF atlas
#include "rttov_setup_brdf_atlas.interface"
#include "rttov_get_brdf.interface"
#include "rttov_deallocate_brdf_atlas.interface"


  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE rttov_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type rttov_in
     integer,pointer :: &
          nPoints,      & ! Number of profiles to simulate
          nLevels,      & ! Number of levels
          nSubCols,     & ! Number of subcolumns
          month           ! Month (needed for surface emissivity calculation)
     real(wp),pointer :: &
          zenang,       & ! Satellite zenith angle
          co2,          & ! Carbon dioxide 
          ch4,          & ! Methane 
          n2o,          & ! n2o 
          co              ! Carbon monoxide
     real(wp),dimension(:),pointer :: &
          surfem          ! Surface emissivities for the channels
     real(wp),dimension(:),pointer :: &
          h_surf,       & ! Surface height
          u_surf,       & ! U component of surface wind
          v_surf,       & ! V component of surface wind
          t_skin,       & ! Surface skin temperature
          p_surf,       & ! Surface pressure
          t2m,          & ! 2 m Temperature
          q2m,          & ! 2 m Specific humidity
          lsmask,       & ! land-sea mask
          latitude,     & ! Latitude
          longitude,    & ! Longitude
          seaice          ! Sea-ice? 
     real(wp),dimension(:,:),pointer :: &
          p,            & ! Pressure @ model levels
          ph,           & ! Pressure @ model half levels
          t,            & ! Temperature 
          q,            & ! Specific humidity
          o3              ! Ozone
     
     ! These fields below are needed ONLY for the RTTOV all-sky brightness temperature
     real(wp),dimension(:,:),pointer :: &
          tca,          & ! Cloud fraction
          cldIce,       & ! Cloud ice
          cldLiq,       & ! Cloud liquid
          fl_rain,      & ! Precipitation flux (startiform+convective rain) (kg/m2/s)
          fl_snow         ! Precipitation flux (stratiform+convective snow)
  end type rttov_in
  
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(NchanIN,platformIN,satelliteIN,instrumentIN,channelsIN)
    integer,intent(in) :: & 
         NchanIN,     & ! Number of channels
         platformIN,  & ! Satellite platform
         satelliteIN, & ! Satellite
         instrumentIN   ! Instrument
    integer,intent(in),dimension(RTTOV_MAX_CHANNELS) :: &
         channelsIN     ! RTTOV channels

    ! Local variables
    character(len=256) :: coef_filename, cld_coef_filename, sat
    integer :: errorstatus, dosolar, imonth

    imonth=1
    dosolar = 0
    nthreads = 1
    nChannels  = NchanIN
    iChannel=channelsIN    

    sat="_"
    IF(satelliteIN.NE.0) THEN
       write(sat,*) satelliteIN
       sat="_"//trim(adjustl(sat))//"_"
    END IF
    COEF_FILENAME= "/pf/b/b380333/work/Tools/RTTOV/rttov121/rtcoef_rttov12/rttov7pred54L/rtcoef_"//trim(platform_name(platformIN))//trim(sat)//trim(inst_name(instrumentIN))//".dat"
    CLD_COEF_FILENAME= "/pf/b/b380333/work/Tools/RTTOV/rttov121/rtcoef_rttov12/cldaer_ir/sccldcoef_"//trim(platform_name(platformIN))//trim(sat)//trim(inst_name(instrumentIN))//".dat"    

    
    ! --------------------------------------------------------------------------
    ! 1. Initialise RTTOV options structure
    ! --------------------------------------------------------------------------
    IF (dosolar == 1) THEN
       opts % rt_ir % addsolar = .TRUE.           ! Include solar radiation
    ELSE
       opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
    ENDIF
    opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
    opts % interpolation % interp_mode = 1       ! Set interpolation method
!!    opts % interpolation % reg_limit_extrap = .TRUE.

    opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
    opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects
    opts % rt_ir % addclouds           = .TRUE.  ! Don't include cloud effects

    opts % rt_ir % ir_scatt_model      = 1 ! DOM
    opts % rt_ir % vis_scatt_model     = 1 ! DOM
    opts % rt_ir % dom_nstreams        = 8 

    opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
    opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
    opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
    opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
    opts % rt_ir % co_data             = .FALSE. !
    opts % rt_ir % so2_data            = .FALSE. !
    opts % rt_mw % clw_data            = .FALSE. 

    opts % config % verbose            = .FALSE.  ! Enable printing of warnings
    opts % config % do_checkinput      = .FALSE.

    ! --------------------------------------------------------------------------
    ! 2. Read coefficients
    ! --------------------------------------------------------------------------
    CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename, &
         file_sccld=cld_coef_filename)
    IF (errorstatus /= errorstatus_success) THEN
       WRITE(*,*) 'fatal error reading coefficients'
       CALL rttov_exit(errorstatus)
    ENDIF

    ! Ensure the options and coefficients are consistent
    CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
    IF (errorstatus /= errorstatus_success) THEN
       WRITE(*,*) 'error in rttov options'
       CALL rttov_exit(errorstatus)
    ENDIF


    ! --------------------------------------------------------------------------
    ! Initialize the emissivity and (if solar) the brdf atlases
    ! --------------------------------------------------------------------------

    IF (coefs%coef%id_sensor == sensor_id_mw .OR. &
         coefs%coef%id_sensor == sensor_id_po) THEN
       atlas_type = atlas_type_mw ! MW atlas
    ELSE
       atlas_type = atlas_type_ir ! IR atlas
    ENDIF
    CALL rttov_setup_emis_atlas(          &
         errorstatus,              &
         opts,                     &
         imonth,                   &
         atlas_type,               & ! Selects MW (1) or IR (2)
         emis_atlas,               &
         path = '/pf/b/b380333/work/Tools/RTTOV/rttov121/emis_data', & ! The default path to atlas data
         coefs = coefs) 


    IF (errorstatus /= errorstatus_success) THEN
       WRITE(*,*) 'error initialising emissivity atlas'
       CALL rttov_exit(errorstatus)
    ENDIF

    IF (opts % rt_ir % addsolar) THEN

       ! Initialise the RTTOV BRDF atlas
       CALL rttov_setup_brdf_atlas(        &
            errorstatus,            &
            opts,                   &
            imonth,                 &
            brdf_atlas,             &
            path='/pf/b/b380333/work/Tools/RTTOV/rttov121/brdf_data', &  ! The default path to atlas data
            coefs = coefs) ! If supplied the BRDF atlas is initialised for this sensor and
       ! this makes the atlas much faster to access
       IF (errorstatus /= errorstatus_success) THEN
          WRITE(*,*) 'error initialising BRDF atlas'
          CALL rttov_exit(errorstatus)
       ENDIF

    ENDIF

 
  END SUBROUTINE COSP_RTTOV_INIT
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
