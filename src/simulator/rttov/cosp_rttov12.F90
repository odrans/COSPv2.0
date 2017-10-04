! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2016, Regents of the University of Colorado
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
! March 2016 - M. Johnston - Original version
! April 2016 - D. Swales   - Modified for use in COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_cosp_rttov

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         surftype_sea,        &
         surftype_land,       &
         watertype_fresh_water, &
         watertype_ocean_water, &         
         sensor_id_mw,        &
         sensor_id_po, &
         wcl_id_stco, &
         wcl_id_stma


  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance,   &
         rttov_opt_param

  ! The rttov_emis_atlas_data type must be imported separately
  USE mod_rttov_emis_atlas, ONLY : &
        rttov_emis_atlas_data, &
        atlas_type_ir, atlas_type_mw

  ! The rttov_brdf_atlas_data type must be imported separately
  USE mod_rttov_brdf_atlas, ONLY : rttov_brdf_atlas_data
  
  USE rttov_unix_env, ONLY : rttov_exit  

  USE parkind1, ONLY : jpim, jprb, jplm
  
  use cosp_kinds,          only : wp
  use mod_cosp_config,     only : RTTOV_MAX_CHANNELS,N_HYDRO,rttovDir

  implicit none

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

  

  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output
  
  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances
  TYPE(rttov_opt_param)            :: cld_opt_param            ! Input cloud optical parameters
  
  TYPE(rttov_emis_atlas_data)      :: emis_atlas               ! Data structure for emissivity atlas
  TYPE(rttov_brdf_atlas_data)      :: brdf_atlas               ! Data structure for BRDF atlas
  
  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: atlas_type  
  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=11)  :: NameOfRoutine = 'example_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: prof_filename, cld_coef_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  REAL(KIND=jprb)    :: trans_out(10)
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch
  INTEGER(KIND=jpim) :: np, nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff, ichan
  INTEGER            :: ios


  !!! OLD STUFFS //
  ! Module parameters
  integer, parameter :: maxlim =  10000
  real(wp),parameter :: eps    =  0.622

  ! Initialization parameters
  integer :: &
       platform,   & ! RTTOV platform
       sensor,     & ! RTTOV instrument
       satellite
  
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
       iChannel      ! RTTOV channel numbers

  ! RTTOV setup and options (set during initialization)
  !!!! \\
  
contains

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE rttov_column
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine rttov_column(nPoints,nLevels_in,nSubCols,q,p,t,o3,ph,h_surf,u_surf,v_surf,     &
                          p_surf,t_skin,t2m,q2m,lsmask,lon,lat,seaice,co2,ch4,n2o,co,    &
                          zenang,lCleanup,                                               &
                          ! Outputs
                          Tb,error,                                                      &
                          surfem,imonth,tca,ciw,clw,rain,snow)
    ! Inputs
    integer,intent(in) :: &
         nPoints, & ! Number of gridpoints
         nLevels_in, & ! Number of vertical levels
         nSubCols   ! Number of subcolumns
    real(wp),intent(in) :: &
         co2,     & ! CO2 mixing ratio (kg/kg)
         ch4,     & ! CH4 mixing ratio (kg/kg)
         n2o,     & ! N2O mixing ratio (kg/kg)
         co,      & ! CO mixing ratio (kg/kg)
         zenang     ! Satellite zenith angle
    real(wp),dimension(nPoints),intent(in) :: &
         h_surf,  & ! Surface height (m)
         u_surf,  & ! Surface u-wind (m/s)
         v_surf,  & ! Surface v-wind (m/s)
         p_surf,  & ! Surface pressure (Pa)
         t_skin,  & ! Skin temperature (K)
         t2m,     & ! 2-meter temperature (K)
         q2m,     & ! 2-meter specific humidity (kg/kg)
         lsmask,  & ! Land/sea mask
         lon,     & ! Longitude (deg)
         lat,     & ! Latitude (deg)
         seaice     ! Seaice fraction (0-1)
    real(wp),dimension(nPoints,nLevels_in),intent(in) :: &
         q,       & ! Specific humidity (kg/kg)
         p,       & ! Pressure(Pa)
         t,       & ! Temperature (K)
         o3         ! Ozone
    real(wp),dimension(nPoints,nLevels_in+1),intent(in) :: &
         ph         ! Pressure @ half-levels (Pa)
    logical,intent(in) :: &
         lCleanup   ! Flag to determine whether to deallocate RTTOV types

    integer(KIND=jpim) :: &
         imonth  ! Month, needed to calculate surface emissivity

    real(wp),dimension(nChannels) :: &
         surfem     ! Surface emissivity for each RTTOV channel

    real(wp),dimension(nPoints,nLevels_in) :: &
         tca       ! Total column cloud amount (0-1)
    real(wp),dimension(nPoints,nLevels_in) :: &
         ciw,    & ! Cloud ice
         clw,    & ! Cloud liquid
         rain,   & ! Precipitation flux (kg/m2/s)
         snow      ! Precipitation flux (kg/m2/s)

    ! Outputs
    real(wp),dimension(nPoints,nChannels) :: &
         Tb        ! RTTOV brightness temperature.
    character(len=128) :: &
         error     ! Error messages (only populated if error encountered)
    
    ! Local variables
    integer :: &
         nloop,rmod,il,istart,istop,za,i,j,subcol,errorstatus,npts_it
    integer :: &
         alloc_status
    real(wp),dimension(nPoints) :: &
         sh_surf
    real(wp),dimension(nPoints,nLevels_in) :: &
         sh,totalice
    real(wp),dimension(nPoints,nChannels) :: &
         Tbs ! Subcolumn brightness temperature
    logical :: &
         use_totalice, mmr_snowrain, cfrac
    logical :: &
         lallSky, & ! Control for type of brightness temperature calculation
                    ! (False(default) => clear-sky brightness temperature, True => All-sky)
         lsfcEmis   ! Control for surface emissivity calculation (true => compute surface emissivity,
                    ! provided that the field "month" is available)

    Tb(:,:)    = 0._wp
    errorstatus = 0_jpim
    
    nprof=nPoints; nlevels=nlevels_in
    ALLOCATE(channel_list(nchannels)); channel_list= ichannel(1:nchannels)
    nthreads=1  
    
    ! Ensure input number of channels is not higher than number stored in coefficient file
    IF (nchannels > coefs % coef % fmv_chn) THEN
       nchannels = coefs % coef % fmv_chn
    ENDIF

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance, &
        init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  
  
  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------
  nch = 0_jpim
  DO j = 1, nprof
     DO jch = 1, nchannels
        nch = nch + 1_jpim
        chanprof(nch)%prof = j
        chanprof(nch)%chan = channel_list(jch)
     ENDDO
  ENDDO

  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  profiles(:) % gas_units = 1 !! tmp, check this, is necessary for water in kg/kg?
  
  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

     profiles(iprof) % p = p(iprof,:)*1E-2 
     profiles(iprof) % t = t(iprof,:)
     profiles(iprof) % q(:) = q(iprof,:) ! I left kg/kg, automatic choice in RRTOV?
     
     Profiles(iprof) % s2m % t = t2m(iprof)
     profiles(iprof) % s2m % q = q2m(iprof)
     profiles(iprof) % s2m % p = p_surf(iprof)*1E-2
     profiles(iprof) % s2m % u = u_surf(iprof)
     profiles(iprof) % s2m % v = v_surf(iprof)
     profiles(iprof) % s2m % wfetc = 100000 !! used typical value given in documentation
    
     profiles(iprof) % skin % t = t_skin(iprof)
     profiles(iprof) % skin % salinity = 0.0 !! tmp, use other typical value
     profiles(iprof) % skin % fastem = (/3.0, 5.0, 15.0, 0.1, 0.3/) !! tmp, typical for land, adjust
    
    if (lsmask(iprof) < 0.5) then
       profiles(iprof)%skin%surftype  = surftype_sea
    else
       profiles(iprof)%skin%surftype  = surftype_land
    endif
    profiles(iprof) %skin % watertype = watertype_fresh_water !! tmp, adapt this to truth, fresh more likely for ICON-DE simulations
    
    profiles(iprof)%elevation = h_surf(iprof) !! tmp, this is set to be 0, why? The elevation is an input of the COSP
    profiles(iprof) % latitude = lat(iprof)
    profiles(iprof) % longitude = lon(iprof)

    profiles(iprof) % zenangle = 0 !! tmp, not sure how to deal with that. Need inst coordinates.
    profiles(iprof) % azangle = 0 !! tmp, not sure how to deal with that. Need inst coordinates.
    profiles(iprof) % sunzenangle = 30 !! tmp, use libradtran function and lat,lon 
    profiles(iprof) % sunazangle = 0 !! tmp, use libradtran function and lat,lon
    
    profiles(iprof)%cfraction  =  0. !! tmp, adapt when clouds are in
    profiles(iprof)%ctp        =  500. !! tmp, adapt when clouds are in

    profiles(iprof) % mmr_cldaer = .TRUE.
    profiles(iprof) % ice_scheme = 2  !! Use baran
    
    profiles(iprof) % cfrac = tca(iprof,:)
    if (lsmask(iprof) < 0.5) then
       profiles(iprof) % cloud(wcl_id_stco,:) = clw(iprof,:)
    else
       profiles(iprof) % cloud(wcl_id_stma,:) = clw(iprof,:)
    end if
    profiles(iprof) % cloud(6,:) = ciw(iprof,:)

  ENDDO


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  ! In this example we have no values for input reflectances
  reflectance(:) % refl_in = 0._jprb

  ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
  ! (all channels in this case)
  calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)

  ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
  reflectance(:) % refl_cloud_top = 0._jprb

  
  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_direct(                &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance) ! inout input/output BRDFs per channel
  ELSE
    CALL rttov_parallel_direct(     &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance,&! inout input/output BRDFs per channel
            nthreads    = nthreads)    ! in    number of threads to use
  ENDIF

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --- Output the results --------------------------------------------------
  
  DO iprof = 1, nprof
     joff = (iprof-1_jpim) * nchannels
     ichan=1
     DO j=1+joff, nchannels+joff
        tb(iprof,ichan) = radiance % bt(j)
        ichan=ichan+1
     END DO
  END DO
  
  ! --- End of output section -----------------------------------------------


  
  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
     WRITE(*,*) 'mem dellocation error'
  ENDIF
  
  CALL rttov_alloc_direct( &
       errorstatus,             &
       0_jpim,                  &  ! 0 => deallocate
       nprof,                   &
       nchanprof,               &
       nlevels,                 &
       chanprof,                &
       opts,                    &
       profiles,                &
       coefs,                   &
       transmission,            &
       radiance,                &
       calcemis=calcemis,       &
       emissivity=emissivity,   &
       calcrefl=calcrefl,       &
       reflectance=reflectance)
  IF (errorstatus /= errorstatus_success) THEN
     WRITE(*,*) 'deallocation error for rttov_direct structures'
     CALL rttov_exit(errorstatus)
  ENDIF

  IF (lCleanup) THEN
     CALL rttov_dealloc_coefs(errorstatus, coefs)
     IF (errorstatus /= errorstatus_success) THEN
        WRITE(*,*) 'coefs deallocation error'
     ENDIF

     ! Deallocate emissivity atlas
     CALL rttov_deallocate_emis_atlas(emis_atlas)

     IF (opts % rt_ir % addsolar) THEN
        ! Deallocate BRDF atlas
        CALL rttov_deallocate_brdf_atlas(brdf_atlas)
     ENDIF

  END IF
     
end subroutine rttov_column

end module mod_cosp_rttov
