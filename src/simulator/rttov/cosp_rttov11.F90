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
         sensor_id_po


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

! Clouds 
#include "rttov_init_opt_param.interface"
#include "rttov_bpr_init.interface"
#include "rttov_bpr_calc.interface"
#include "rttov_bpr_dealloc.interface"
#include "rttov_legcoef_calc.interface"
  
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
  CHARACTER(LEN=256) :: prof_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: imonth  
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
  INTEGER(KIND=jpim) :: iprof, joff
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
                          ! Optional arguments for surface emissivity calculation.
                          surfem,month,                                                  &
                          ! Optional arguments for all-sky calculation.
                          tca,ciw,clw,rain,snow)
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

    ! Optional inputs (Needed for surface emissivity calculation)
    integer,optional :: &
         month      ! Month (needed to determine table to load)
    real(wp),dimension(nChannels),optional :: &
         surfem     ! Surface emissivity for each RTTOV channel

    ! Optional inputs (Needed for all-sky calculation)
    real(wp),dimension(nPoints,nLevels_in),optional :: &
         tca       ! Total column cloud amount (0-1)
    real(wp),dimension(nPoints,nSubCols,nLevels_in),optional :: &
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
    real(wp),dimension(nSubCols,nPoints,nChannels) :: &
         Tbs ! Subcolumn brightness temperature
    logical :: &
         use_totalice, mmr_snowrain, cfrac
    logical :: &
         lallSky, & ! Control for type of brightness temperature calculation
                    ! (False(default) => clear-sky brightness temperature, True => All-sky)
         lsfcEmis   ! Control for surface emissivity calculation (true => compute surface emissivity,
                    ! provided that the field "month" is available)

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

    !! Adapt the units //

    !! \\

    
    write(*,*) "Start with rttov"
    
      ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  COEF_FILENAME="/pf/b/b380333/work/Tools/RTTOV/rttov121/rtcoef_rttov12/rttov7pred54L/rtcoef_msg_4_seviri.dat" 
  PROF_FILENAME="/pf/b/b380333/work/Tools/RTTOV/rttov121/rttov_test/test_example.1/prof.dat"
  nprof=4
  nlevels=nlevels_in
  imonth=1  
  dosolar=0
  nchannels=2
  ALLOCATE(channel_list(nchannels))  
  channel_list=(/7,8/)
  nthreads=1  

  
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
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects
  opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects

  opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
  opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
  opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
  opts % rt_ir % co_data             = .FALSE. !
  opts % rt_ir % so2_data            = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !

  opts % config % verbose            = .TRUE.  ! Enable printing of warnings


  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

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

  ! Initialise the RTTOV emissivity atlas
  ! (this loads the default IR/MW atlases: use the atlas_id argument to select alternative atlases)
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
              coefs = coefs) ! This is mandatory for the CNRM MW atlas, ignored by TELSEM2;
                             ! if supplied for IR atlases they are initialised for this sensor
                             ! and this makes the atlas much faster to access
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

  !===============================================
  !========== Read profiles == start =============

  profiles(:) % gas_units = 1 !! tmp, check this, is necessary for water in kg/kg?
  
  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof
     
     profiles(iprof) % p = p(iprof,:)*1E-2 
     profiles(iprof) % t = t(iprof,:)
     profiles(iprof) % q(:) = q(iprof,:) ! I left kg/kg, automatic choice in RRTOV?
     
     profiles(iprof) % s2m % t = t2m(iprof)
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
    
  ENDDO

  !========== READ profiles == end =============
  !=============================================


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
     write(*,*) "run rttov now"
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

  ! Open output file where results are written
  OPEN(ioout, file='output_'//NameOfRoutine//'.dat', status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' Instrument ', inst_name(coefs % coef % id_inst)
  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' '
  CALL rttov_print_opts(opts, lu=ioout)

  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    nprint = 1 + INT((nchannels-1)/10)
    WRITE(ioout,*)' '
    WRITE(ioout,*)' Profile ', iprof

    CALL rttov_print_profile(profiles(iprof), lu=ioout)

    WRITE(ioout,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
    WRITE(ioout,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
    WRITE(ioout,222) (radiance % bt(j), j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SATELLITE REFLECTANCES (BRF):'
      WRITE(ioout,444) (radiance % refl(j), j = 1+joff, nchannels+joff)
    ENDIF
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED RADIANCES (mW/m2/sr/cm-1):'
    WRITE(ioout,222) (radiance % total(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED OVERCAST RADIANCES:'
    WRITE(ioout,222) (radiance % cloudy(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE TO SPACE TRANSMITTANCE:'
    WRITE(ioout,4444) (transmission % tau_total(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE EMISSIVITIES:'
    WRITE(ioout,444) (emissivity(j) % emis_out, j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SURFACE BRDF:'
      WRITE(ioout,444) (reflectance(j) % refl_out, j = 1+joff, nchannels+joff)
    ENDIF

    IF (nchannels <= 20) THEN
      DO np = 1, nprint
          WRITE(ioout,*)' '
          WRITE(ioout,*)'Level to space transmittances for channels'
          WRITE(ioout,1115) (chanprof(j+joff) % chan, &
                    j = 1+(np-1)*10, MIN(np*10, nchannels))
          DO ilev = 1, nlevels
            DO j = 1 + (np-1)*10, MIN(np*10, nchannels)
              ! Select transmittance based on channel type (VIS/NIR or IR)
              IF (coefs % coef % ss_val_chn(chanprof(j+joff) % chan) == 2) THEN
                trans_out(j - (np-1)*10) = transmission % tausun_levels_path1(ilev,j+joff)
              ELSE
                trans_out(j - (np-1)*10) = transmission % tau_levels(ilev,j+joff)
              ENDIF
            ENDDO
            WRITE(ioout,4445) ilev, trans_out(1:j-1-(np-1)*10)
          ENDDO
          WRITE(ioout,1115) (chanprof(j+joff) % chan, &
                    j = 1+(np-1)*10, MIN(np*10, nchannels))
      ENDDO
    ENDIF
  ENDDO

  ! Close output file
  CLOSE(ioout, iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error closing the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  ! --- End of output section -----------------------------------------------


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for rttov_direct
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

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
1115 FORMAT(3X,10I8)
222  FORMAT(1X,10F8.2)
444  FORMAT(1X,10F8.3)
4444 FORMAT(1X,10F8.4)
4445 FORMAT(1X,I2,10F8.4)
777  FORMAT(/,A,A9,I3)
    

    write(*,*) "Done with rttov!"
    stop

    
    ! ! Initialize some things
    ! totalice   = 0._wp
    ! Tbs(:,:,:) = 0._wp
    ! Tb(:,:)    = 0._wp
    ! error      = ''

    ! ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ! Setup for call to RTTOV
    ! ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ! First, check to see if we are doing an all-sky or clear-sky calculation brightness
    ! ! temperature
    ! lallSky = .false.
    ! if (present(tca) .and. present(clw) .and. present(ciw) .and. present(rain)           &
    !     .and. present(snow)) lallSky=.true.

    ! ! Check to see if we need to compute the surface emissivity (defualt is to compute
    ! ! surface emissivity using the atlas tables)
    ! lsfcEmis = .true.
    ! if (present(surfem)) lsfcEmis = .false.
    
    ! ! We also need the month for the emissivity atlas, so check...
    ! if (.not. present(month)) lsfcEmis = .false.

    ! if (lsfcEmis .eq. .false. .and. .not. present(surfem)) then
    !    error = 'ERROR (rttov_column): User did not provide surface emissivity and did not '//&
    !            'request the surface emissivity to be calculated!!!'
    !    return
    ! endif

    ! ! Convert specific humidity to ppmv
    ! sh       = ( q   / ( q   + eps * ( 1._wp - q ) ) ) * 1e6   
    ! sh_surf  = ( q2m / ( q2m + eps * ( 1._wp - q2m ) ) ) * 1e6

    ! ! Settings unique to all-sky call.
    ! use_totalice = .false.
    ! mmr_snowrain = .true.
    ! cfrac        = .true.
    ! opts_scatt%lusercfrac = cfrac

    ! ! RTTOV can handle only about 100 profiles at a time (fixme: check this with roger),
    ! ! so we are putting a loop of 100
    ! nloop  =  npoints / maxlim
    ! rmod   =  mod( npoints, maxlim )
    ! if( rmod .ne. 0 ) then
    !    nloop = nloop + 1
    ! endif

    ! ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ! Initialize emissivity atlas data for chosen sensor.
    ! ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! call rttov_setup_emis_atlas(errorstatus,opts,month,coef_rttov,path=trim(rttovDir)//"emis_data/")
    ! if (errorstatus /= errorstatus_success) then
    !    error = 'ERROR (rttov_column): Error reading emis atlas data!'
    !    return
    ! endif
    
    ! ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ! Some quality control prior to RTTOV call
    ! ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ! Ensure the options and coefficients are consistent
    ! if(opts_scatt%config%do_checkinput) then
    !    call rttov_user_options_checkinput(errorstatus, opts, coef_rttov)
    !    if (errorstatus /= errorstatus_success) then
    !       error =  'ERROR (rttov_column): Error when checking input data!'
    !       return
    !    endif
    ! endif
    
    ! ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ! Call to RTTOV
    ! ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ! Looping over maxlim number of profiles
    ! do il = 1, nloop
    !    istart  =  (il - 1) * maxlim + 1
    !    istop   =  min(il * maxlim, npoints)     
    !    if( ( il .eq. nloop ) .and. ( rmod .ne. 0 ) ) then
    !       npts_it = rmod
    !    else
    !       npts_it = maxlim
    !    endif
    !    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    ! Clear-sky brightness temperature
    !    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    !    if (.not. lallSky) then
    !       call rttov_multprof(nChannels,iChannel,surfem,npts_it,nLevels,platform,         &
    !                           satellite,sensor,opts,coef_rttov,zenang,                    &
    !                           p(istart:istop,:)/100._wp,t(istart:istop,:),                &
    !                           sh(istart:istop,:),(mdry/mo3)*o3(istart:istop,:)*1e6,       &
    !                           (mdry/mco2)*co2*1e6,(mdry/mch4)*ch4*1e6,(mdry/mn2o)*n2o*1e6,&
    !                           (mdry/mco)*co*1e6,h_surf(istart:istop),u_surf(istart:istop),&
    !                           v_surf(istart:istop),t_skin(istart:istop),                  &
    !                           p_surf(istart:istop)/100.,t2m(istart:istop),                &
    !                           sh_surf(istart:istop),lsmask(istart:istop),                 &
    !                           seaice(istart:istop),lat(istart:istop),lon(istart:istop),   &
    !                           Tb(istart:istop,:))  
    !    endif
    !    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    ! All-sky brightness temperature
    !    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    if (lallSky) then
    !       ! Loop over all subcolumns
    !       do subcol = 1, nSubCols	
    !          ! Call RTTOV
    !          call cosp_rttov_mwscatt(nChannels,iChannel,surfem,nPoints,nlevels,platform,  &
    !                                  satellite,sensor,opts,opts_scatt,coef_rttov,         &
    !                                  coef_scatt,zenang,p(istart:istop,:)/100._wp,         &
    !                                  ph(istart:istop,:)/100._wp,t(istart:istop, :),       &
    !                                  sh(istart:istop, :),                                 &
    !                                  (mdry/mo3)*o3(istart:istop,:)*1e6,                   &
    !                                  clw(istart:istop,subcol,:),                          &
    !                                  ciw(istart:istop,subcol,:),tca(istart:istop, :),     &
    !                                  totalice(istart:istop,:),snow(istart:istop,subcol,:),& 
    !                                  rain(istart:istop,subcol,:),(mdry/mco2)*co2*1e6,     &
    !                                  (mdry/mch4)*ch4*1e6,(mdry/mn2o)*n2o*1e6,             &
    !                                  (mdry/mco)*co*1e6,h_surf(istart:istop),              &
    !                                  u_surf(istart:istop),v_surf(istart:istop),           &
    !                                  t_skin(istart:istop), p_surf(istart:istop)/100.,     &
    !                                  t2m(istart:istop),sh_surf(istart:istop),             &
    !                                  lsmask(istart:istop),seaice(istart:istop),           &
    !                                  lat(istart:istop),lon(istart:istop), use_totalice,   &
    !                                  mmr_snowrain,cfrac,Tbs(subcol,istart:istop,:))
    !       enddo 
    !    endif 
    ! enddo

    ! ! For all-sky calculation we need to average together all of the cloudy subcolumns.
    ! if (lallSky) then
    !    do subcol = 1, nSubCols
    !       Tb = Tb + tbs(subcol,:,:)
    !    enddo
    !    Tb = Tb/nSubCols
    ! endif
    
    ! ! Free up space
    ! if (lCleanup) then
    !    call rttov_dealloc_coefs(errorstatus,coef_rttov)
    !    call rttov_deallocate_emis_atlas(coef_rttov)
    !    if (lallSky) call rttov_dealloc_scattcoeffs(coef_scatt)    
    ! endif
    
  end subroutine rttov_column

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE rttov_multprof
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   subroutine rttov_multprof( &
!        nch_in,     & ! number of channels
!        ichan_in,   & ! channel indices
!        surfem_in,  & ! surface emissivity values
!        prf_num_in, & ! number of profiles to simulate
!        nlevels_in, & ! number of pressure levels
!        plat_in,    & ! platform number
!        sat_in,     & ! satellite number
!        sens_in,    & ! instrument number
!        opts,       &
!        coef_rttov, &
!        zenang_in,  & ! zenith angle
!        p_in,       & ! pressure [hpa]
!        t_in,       & ! temperature [ k ]
!        q_in,       & ! specific humidity [ ppmv ]
!        o3_in,      & ! ozone vmr [ ppmv ]
!        co2_in,     & ! co2 vmr [ ppmv ] *this is a single value*
!        ch4_in,     & ! ch4 vmr [ ppmv ] *this is a single value*
!        n2o_in,     & ! n2o vmr [ ppmv ] *this is a single value*
!        co_in,      & ! co vmr [ ppmv ]  *this is a single value*
!        h_surf,     & ! surface height [ m ]
!        u_surf,     & ! u wind at 10 m  [ m/s ]
!        v_surf,     & ! v wind at 10 m [ m/s ]
!        t_skin,     & ! skin temperatre [ k ]
!        p_surf,     & ! surface pressure
!        t_surf,     & ! 1.5 m temperature [ k ]
!        q_surf,     & ! 1.5 m specific humidity [ ppmv ]
!        lsmask,     & ! land sea mask
!        seaice,     & ! seaice fraction
!        latitude,   & ! latitude [ deg north ]
!        longitude,  & ! longitude [ deg east ]
!        tbs         & ! brightness temperature [ k ] (output)
!        )

!     !------ input arguments. no rttov kinds should be used here -----------------
!     integer, intent(in)  :: nch_in             ! number of channels to be computed
!     integer, intent(in)  :: ichan_in(nch_in)   ! indices of selected channels
!     real(wp), intent(in)     :: surfem_in(nch_in)  ! surface emissivities for the channels
!     integer, intent(in)  :: prf_num_in
!     integer, intent(in)  :: nlevels_in
!     integer, intent(in)  :: plat_in  ! satellite platform
!     integer, intent(in)  :: sat_in   ! satellite number
!     integer, intent(in)  :: sens_in  ! satellite sensor
!     real(wp), intent(in)     :: zenang_in          ! satellite zenith angle

!     type(rttov_options)  :: opts
!     type(rttov_coefs)    :: coef_rttov

!     real(wp), intent(in)     :: p_in(prf_num_in, nlevels_in)  ! pressure profiles
!     real(wp), intent(in)     :: t_in(prf_num_in, nlevels_in)  ! temperature profiles
!     real(wp), intent(in)     :: q_in(prf_num_in, nlevels_in)  ! humidity profiles
!     real(wp), intent(in)     :: o3_in(prf_num_in, nlevels_in) ! ozone profiles

!     ! the following trace gases contain constant values
!     real(wp), intent(in) ::  co2_in ! carbon dioxide
!     real(wp), intent(in) ::  ch4_in ! methane
!     real(wp), intent(in) ::  n2o_in ! n2o
!     real(wp), intent(in) ::  co_in  ! carbon monoxide
!     real(wp), intent(in) ::  h_surf(prf_num_in)         ! surface height
!     real(wp), intent(in) ::  u_surf(prf_num_in)         ! u component of surface wind
!     real(wp), intent(in) ::  v_surf(prf_num_in)         ! v component of surface wind
!     real(wp), intent(in) ::  t_skin(prf_num_in)         ! surface skin temperature
!     real(wp), intent(in) ::  p_surf(prf_num_in)                  ! surface pressure
!     real(wp), intent(in) ::  t_surf(prf_num_in)         ! 1.5 m temperature
!     real(wp), intent(in) ::  q_surf(prf_num_in)         ! 1.5 m specific humidity
!     real(wp), intent(in) ::  lsmask(prf_num_in)         ! land-sea mask
!     real(wp), intent(in) ::  seaice(prf_num_in)         ! sea-ice fraction
!     real(wp), intent(in) ::  latitude(prf_num_in)       ! latitude
!     real(wp), intent(in) ::  longitude(prf_num_in)      ! longitude

!     real(wp), intent(inout) :: tbs(prf_num_in, nch_in)  ! tbs (in the right format)

!     !------ local variables. use only rttov kinds or derived types.
!     !       logical variables are declared with the same kind
!     !       as integers, as they are affected inthe same way by flags like -qintsize=8

!     !     type(rttov_options)                  :: opts           ! options structure
!     !     type(rttov_coefs),       allocatable :: coefs(:)       ! coefficients structure
!     type(rttov_chanprof),    allocatable :: chanprof(:)    ! input channel/profile list
!     type(profile_type),      allocatable :: profiles(:)    ! input profiles
!     logical,      allocatable :: calcemis(:)    ! flag to indicate calculation of emissivity within rttov
!     type(rttov_emissivity),  allocatable :: emissivity(:)  ! input/output surface emissivity
!     type(transmission_type)              :: transmission   ! output transmittances
!     type(radiance_type)                  :: radiance       ! output radiances

!     integer, allocatable :: instrument(:,:) ! instrument id (3 x n_instruments)
!     integer, allocatable :: nchan(:) ! number of channels per instrument
!     integer, allocatable :: ichan(:,:)   ! channel list per instrument

!     integer :: asw
!     integer :: mxchn
!     integer :: nrttovid     ! maximum number of instruments
!     integer :: no_id        ! instrument loop index
!     integer :: i, j, jch
!     integer :: nprof  ! number of calls to rttov
!     integer :: nch ! intermediate variable
!     integer :: errorstatus
!     integer :: ich, ich_temp, nchanprof, nchannels, chan
!     integer :: alloc_status(60)

!     real(wp),    allocatable :: input_emissivity (:)

!     character (len=14)  :: nameofroutine = 'rttov_multprof'

!     logical :: refrac, solrad, laerosl, lclouds, lsun, all_channels

!     ! local variables for input arguments that need type casting to avoid type-mismatch with
!     ! rttov kinds. this happens with some compiler flags (-qintsize=8).
!     integer  :: prof_num
!     integer  :: nlevels

!     ! --------------------------------------------------------------------------
!     ! 0. initialise cosp-specific things
!     ! --------------------------------------------------------------------------

!     ! type-casting of input arguments that need to be passed to rttov
!     prof_num = prf_num_in
!     nlevels  = nlevels_in
!     nprof = prof_num

!     ! currently we plan to calculate only 1 instrument per call
!     nrttovid  =  1
!     mxchn  =  nch_in

!     errorstatus     = 0
!     alloc_status(:) = 0

!     !     allocate(coefs(nrttovid), stat = alloc_status(1))

!     !     allocate(instrument(3, nrttovid), stat = alloc_status(4))

!     !maximum number of channels allowed for one instrument is mxchn
!     !    allocate(surfem(nch_in, nrttovid), stat = alloc_status(11))
!     allocate(ichan(nch_in, nrttovid), stat = alloc_status(12))
!     call rttov_error('ichan mem allocation error for profile array' , lalloc = .true.)


!     do no_id = 1, nrttovid
!        ichan(:, no_id)   = ichan_in
!     enddo

!     asw = 1 ! switch for allocation passed into rttov subroutines

!     ! allocate input profile arrays
!     allocate(profiles(nprof), stat = alloc_status(1))
!     call rttov_error('Profile mem allocation error' , lalloc = .true.)

!     call rttov_alloc_prof(     &
!          errorstatus,             &
!          nprof,                   &
!          profiles,                &
!          nlevels,                 &
!          opts,                    &
!          asw,                     &
!          coefs = coef_rttov,      &
!          init = .true.)
!     call rttov_error('Profile 2 mem allocation error' , lalloc = .true.)
!     ! --------------------------------------------------------------------------
!     ! 5. store profile data in profile type
!     ! --------------------------------------------------------------------------
!     do i = 1, nprof
!        profiles(i)%p(:) =  p_in(i, :)
!        profiles(i)%t(:) =  t_in(i, :)
!        profiles(i)%q(:) =  q_in(i, :)

!        where(profiles(i)%q(:) < 1e-4)
!           profiles(i)%q(:) = 1e-4
!        end where

!        profiles(i)%cfraction  =  0.
!        profiles(i)%ctp        =  500.

!        ! 2m parameters
!        profiles(i)%s2m%p  =  p_surf(i)
!        profiles(i)%s2m%t  =  t_surf(i)
!        profiles(i)%s2m%q  =  q_surf(i)
!        profiles(i)%s2m%u  =  u_surf(i) ! dar: hard-coded at 2ms-1?
!        profiles(i)%s2m%v  =  v_surf(i) ! dar: hard-coded at 2ms-1?
!        profiles(i)%s2m%wfetc  =  10000. ! dar: default?

!        ! skin variables for emissivity calculations
!        profiles(i)%skin%t          =  t_skin(i)

!        ! fastem coefficients - for mw calculations
!        profiles(i)%skin%fastem(1)  =  3.0
!        profiles(i)%skin%fastem(2)  =  5.0
!        profiles(i)%skin%fastem(3)  =  15.0
!        profiles(i)%skin%fastem(4)  =  0.1
!        profiles(i)%skin%fastem(5)  =  0.3

!        profiles(i)%zenangle      = zenang_in ! pass in from cosp

!        profiles(i)%azangle       = 0. ! hard-coded in rttov9 int

!        profiles(i)%latitude      = latitude(i)
!        profiles(i)%longitude     = longitude(i)
!        profiles(i)%elevation     = h_surf(i)

!        profiles(i)%sunzenangle   = 0. ! hard-coded in rttov9 int
!        profiles(i)%sunazangle    = 0. ! hard-coded in rttov9 int

!        ! surface type
!        ! land-sea mask indicates proportion of land in grid
!        if (lsmask(i) < 0.5) then
!           profiles(i)%skin%surftype  = surftype_sea
!        else
!           profiles(i)%skin%surftype  = surftype_land
!        endif
!        ! sea-ice fraction
!        if (seaice(i) >= 0.5) then
!           profiles(i)%skin%surftype  = surftype_seaice
!        endif

!        ! dar: hard-coded to 1 (=ocean water) in rttov 9 int
!        profiles(i)%skin%watertype = 1
!        profiles(i) %idg         = 0.
!        profiles(i) %ish         = 0.
!     enddo
!     ! end of 5.

!     ich_temp = 1
!     nchannels = nch_in
!     do no_id = 1, nrttovid

!        ! --------------------------------------------------------------------------
!        ! 3. build the list of profile/channel indices in chanprof
!        ! --------------------------------------------------------------------------

!        allocate(nchan(nprof))     ! number of channels per profile
!        nchan(:) = size(ichan(:,no_id))  ! = nch_in

!        ! size of chanprof array is total number of channels over all profiles
!        ! square in this case - here same channels done for all profiles
!        nchanprof = sum(nchan(:))

!        ! pack channels and input emissivity arrays
!        allocate(chanprof(nchanprof))
!        !      allocate(emis(nchanprof))
!        chanprof(:)%chan =0
       
!        nch = 0
!        do j = 1, nprof
!           do jch = 1, nchan(j)
!              nch = nch + 1
!              chanprof(nch)%prof = j
!              if(ichan(jch, no_id) < 1) then
!                 errorstatus = errorstatus_fatal
!                 call rttov_error('Sensor channel number must be 1 or greater' , lalloc = .true.)
!              else
!                 chanprof(nch)%chan = ichan(jch, no_id)
!              endif
!           enddo
!        enddo
!        ! end of 3.

!        ! allocate output radiance arrays
!        call rttov_alloc_rad( &
!             errorstatus,        &
!             nchanprof,          &
!             radiance,           &
!             nlevels - 1,   & ! nlayers
!             asw)
!        call rttov_error('allocation error for radiance arrays' , lalloc = .true.)

!        ! allocate transmittance structure
!        call rttov_alloc_transmission( &
!             errorstatus,             &
!             transmission,            &
!             nlevels - 1,        &
!             nchanprof,               &
!             asw,                     &
!             init=.true.)
!        call rttov_error('allocation error for transmission arrays' , lalloc = .true.)

!        ! allocate arrays for surface emissivity
!        allocate(calcemis(nchanprof), stat=alloc_status(1))
!        allocate(emissivity(nchanprof), stat=alloc_status(2))
!        call rttov_error('mem allocation error for emissivity arrays' , lalloc = .true.)

!        call rttov_get_emis(       &
!             & errorstatus, &
!             & opts,        &
!             & chanprof,    &
!             & profiles,    &
!             & coef_rttov,  &
!                                 !& resolution=resolution, &  ! *** MW atlas native
!                                 !                  resolution is 0.25 degree lat/lon; if you know better
!                                 !                  value for satellite footprint (larger than this) then
!                                 !                  you can specify it here
!             & emissivity=emissivity(:)%emis_in)
!        !            & emissivity(:)%emis_in)

!        call rttov_error('Get emissivity error' , lalloc = .true.)
!        calcemis(:) = .false.
!        ! calculate emissivity for missing and ocean location (fastem)
!        where (emissivity(:)%emis_in <= 0.0)
!           calcemis(:) = .true.
!        endwhere

!        call rttov_direct(         &
!             errorstatus,               &! out
!             chanprof,                  &
!             opts,                      &
!             profiles,                  &! in
!             coef_rttov,              &! in
!             transmission,              &! out
!             radiance,                  &
!             calcemis = calcemis,       &! in
!             emissivity = emissivity)    ! inout
!        call rttov_error('rttov_direct error', lalloc = .true.)

!        tbs(1:prof_num , ich_temp:ich_temp + size(ichan(:,no_id)) - 1) = &
!             transpose(reshape(radiance%bt(1:nchanprof), (/ size(ichan(:,no_id)), prof_num/) ))

!        ich_temp = ich_temp + size(ichan(:,no_id))

!        ! --------------------------------------------------------------------------
!        ! 8. deallocate all rttov arrays and structures
!        ! --------------------------------------------------------------------------
!        deallocate (nchan,       stat=alloc_status(3))
!        deallocate (chanprof,    stat=alloc_status(4))
!        deallocate (emissivity,  stat=alloc_status(5))
!        deallocate (calcemis,    stat=alloc_status(6))
!        call rttov_error('rttov array deallocation error', lalloc = .true.)

!        asw = 0 ! switch for deallocation passed into rttov subroutines

!        ! deallocate radiance arrays
!        call rttov_alloc_rad(errorstatus, nchannels, radiance, nlevels - 1, asw)
!        call rttov_error('radiance deallocation error', lalloc = .true.)

!        ! deallocate transmission arrays
!        call rttov_alloc_transmission(errorstatus, transmission, nlevels - 1, nchannels, asw)
!        call rttov_error('transmission deallocation error', lalloc = .true.)

!     enddo

!     ! deallocate profile arrays
!     call rttov_alloc_prof(errorstatus, nprof, profiles, nlevels, opts, asw)
!     call rttov_error('profile deallocation error', lalloc = .true.)

!     deallocate(profiles, stat=alloc_status(1))
!     call rttov_error('mem deallocation error for profile array', lalloc= .true.)

!   contains

!     subroutine rttov_error(msg, lalloc)
!       character(*) :: msg
!       logical  :: lalloc

!       if(lalloc) then
!          if (any(alloc_status /= 0)) then
!             write(*,*) msg
!             errorstatus = 1
!             call rttov_exit(errorstatus)
!          endif
!       else
!          if (errorstatus /= errorstatus_success) then
!             write(*,*) msg
!             call rttov_exit(errorstatus)
!          endif
!       endif
!     end subroutine rttov_error

!   end subroutine rttov_multprof
!   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   !----------------- subroutine cosp_rttov_mwscatt ---------------
!   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!   subroutine cosp_rttov_mwscatt(&
!        nch_in,     & ! number of channels
!        ichan_in,   & ! channel indices
!        surfem_in,  & ! surface emissivity values
!        prf_num_in, & ! number of profiles to simulate
!        nlevels_in, & ! number of pressure levels
!        plat_in,    & ! platform number
!        sat_in,     & ! satellite number
!        sens_in,    & ! instrument number
!        opts,       &
!        opts_scatt, &
!        coef_rttov, &
!        coef_scatt, &
!        zenang_in,  & ! zenith angle
!        p_in,       & ! pressure [hpa]
!        ph_in,      & ! pressure on half levels [hpa]
!        t_in,       & ! temperature [ k ]
!        q_in,       & ! specific humidity [ ppmv ]
!        o3_in,      & ! ozone vmr [ ppmv ]
!        clw_in,     & ! cloud water [0-1]
!        ciw_in,     & ! cloud ice [0-1]
!        cc_in,      & ! effective cloud fraction [0-1]
!        totalice_in,& ! total ice, except snow [kg/kg] or [kg/m2/s]
!        sp_in,      & ! solid precip with snow [kg/kg] or [kg/m2/s]
!        rain_in,    & ! total liquid water [kg/kg] or [kg/m2/s]
!        co2_in,     & ! co2 vmr [ ppmv ] *this is a single value*
!        ch4_in,     & ! ch4 vmr [ ppmv ] *this is a single value*
!        n2o_in,     & ! n2o vmr [ ppmv ] *this is a single value*
!        co_in,      & ! co vmr [ ppmv ]  *this is a single value*
!        h_surf,     & ! surface height [ m ]
!        u_surf,     & ! u wind at 10 m  [ m/s ]
!        v_surf,     & ! v wind at 10 m [ m/s ]
!        t_skin,     & ! skin temperatre [ k ]
!        p_surf,     & ! surface pressure
!        t_surf,     & ! 1.5 m temperature [ k ]
!        q_surf,     & ! 1.5 m specific humidity [ ppmv ]
!        lsmask,     & ! land sea mask
!        seaice,     & ! seaice fraction
!        latitude,   & ! latitude [ deg north ]
!        longitude,  & ! longitude [ deg east ]
!        use_totalice,& ! separate ice and snow, or total ice hydrometeor types
!        mmr_snowrain,& ! set units for snow and rain: if true units are kg/kg (the default)
!        cfrac,       & ! opts_scatt%lusercfrac=true., supply the effective cloud fraction
!        tbs          & ! brightness temperature [ k ] (output)
!        )





!     implicit none

!     !------ input arguments. no rttov kinds should be used here -----------------
!     integer, intent(in)  :: nch_in             ! number of channels to be computed
!     integer, intent(in)  :: ichan_in(nch_in)   ! indices of selected channels
!     real(wp), intent(in)     :: surfem_in(nch_in)  ! surface emissivities for the channels
!     integer, intent(in)  :: prf_num_in
!     integer, intent(in)  :: nlevels_in
!     integer, intent(in)  :: plat_in  ! satellite platform
!     integer, intent(in)  :: sat_in   ! satellite number
!     integer, intent(in)  :: sens_in  ! satellite sensor
!     real(wp), intent(in)     :: zenang_in          ! satellite zenith angle

!     type(rttov_options)       :: opts
!     type(rttov_options_scatt) :: opts_scatt
!     type(rttov_coefs)         :: coef_rttov
!     type(rttov_scatt_coef)    :: coef_scatt

!     real(wp), intent(in)     :: p_in(prf_num_in, nlevels_in)  ! pressure profiles
!     real(wp), intent(in)     :: t_in(prf_num_in, nlevels_in)  ! temperature profiles
!     real(wp), intent(in)     :: q_in(prf_num_in, nlevels_in)  ! humidity profiles
!     real(wp), intent(in)     :: o3_in(prf_num_in, nlevels_in) ! ozone profiles
!     real(wp), intent(in)     :: clw_in(prf_num_in, nlevels_in)
!     real(wp), intent(in)     :: ciw_in(prf_num_in, nlevels_in)
!     real(wp), intent(in)     :: cc_in(prf_num_in, nlevels_in)
!     real(wp), intent(in)     :: totalice_in(prf_num_in, nlevels_in)
!     real(wp), intent(in)     :: sp_in(prf_num_in, nlevels_in)
!     real(wp), intent(in)     :: rain_in(prf_num_in, nlevels_in)
!     real(wp), intent(in)     :: ph_in(prf_num_in, nlevels_in+1)

!     ! the following trace gases contain constant values
!     real(wp), intent(in) ::  co2_in ! carbon dioxide
!     real(wp), intent(in) ::  ch4_in ! methane
!     real(wp), intent(in) ::  n2o_in ! n2o
!     real(wp), intent(in) ::  co_in  ! carbon monoxide
!     real(wp), intent(in) ::  h_surf(prf_num_in)         ! surface height
!     real(wp), intent(in) ::  u_surf(prf_num_in)         ! u component of surface wind
!     real(wp), intent(in) ::  v_surf(prf_num_in)         ! v component of surface wind
!     real(wp), intent(in) ::  t_skin(prf_num_in)         ! surface skin temperature
!     real(wp), intent(in) ::  p_surf(prf_num_in)                  ! surface pressure
!     real(wp), intent(in) ::  t_surf(prf_num_in)         ! 1.5 m temperature
!     real(wp), intent(in) ::  q_surf(prf_num_in)         ! 1.5 m specific humidity
!     real(wp), intent(in) ::  lsmask(prf_num_in)         ! land-sea mask
!     real(wp), intent(in) ::  seaice(prf_num_in)         ! seaice fraction
!     real(wp), intent(in) ::  latitude(prf_num_in)       ! latitude
!     real(wp), intent(in) ::  longitude(prf_num_in)      ! longitude
!     logical, intent(in) :: cfrac, use_totalice, mmr_snowrain

!     real(wp), intent(inout) :: tbs(prf_num_in, nch_in)  ! tbs (in the right format)
!     !****************** local variables **********************************************
!     logical       , allocatable :: calcemis   (:)
!     type(rttov_emissivity)   , allocatable :: emissivity  (:)
!     integer       , allocatable :: frequencies (:)
!     type(rttov_chanprof)     , allocatable :: chanprof    (:) ! channel and profile indices
!     type(profile_type)       , allocatable :: profiles    (:)
!     type(profile_cloud_type) , allocatable :: cld_profiles(:)

!     integer, allocatable :: ichan(:,:)   ! channel list per instrument

!     integer         :: errorstatus
!     type (radiance_type)       :: radiance
!     ! 	type (rttov_options)       :: opts     ! defaults to everything optional switched off
!     ! 	type (rttov_options_scatt) :: opts_scatt
!     ! 	type (rttov_coefs)         :: coef_rttov
!     ! 	type (rttov_scatt_coef)    :: coef_scatt

!     ! 	integer, allocatable :: instrument (:,:)
!     integer :: j,k,asw
!     integer :: nchanxnprof, ninstruments
!     real(wp)     :: zenangle
!     character (len=256) :: outstring
!     integer :: alloc_status(60)

! #include "rttov_init_rad.interface"
! #include "rttov_scatt_setupindex.interface"
! #include "rttov_scatt.interface"
! #include "rttov_alloc_rad.interface"
! #include "rttov_alloc_prof.interface"
! #include "rttov_alloc_scatt_prof.interface"
! #include "rttov_get_emis.interface"
! #include "rttov_boundaryconditions.interface"

!     errorstatus     = 0
!     alloc_status(:) = 0
!     ninstruments = 1 ! number of sensors or platforms

!     allocate(ichan(nch_in, ninstruments), stat = alloc_status(3))

!     do j = 1, ninstruments
!        ichan(:, j) = ichan_in
!     enddo

!     nchanxnprof =  prf_num_in * nch_in            ! total channels to simulate * profiles

!     allocate (chanprof(nchanxnprof))
!     allocate (frequencies(nchanxnprof))
!     allocate (emissivity(nchanxnprof))
!     allocate (calcemis(nchanxnprof))
!     allocate (profiles(prf_num_in))
!     allocate (cld_profiles(prf_num_in))

!     ! request rttov / fastem to calculate surface emissivity
!     calcemis  = .true.
!     emissivity % emis_in = 0.0

!     ! setup indices
!     call rttov_scatt_setupindex ( 	&
! 	 & prf_num_in,      			& ! in
! 	 & nch_in,          			& ! in
! 	 & coef_rttov%coef, 			& ! in
! 	 & nchanxnprof,       			& ! in
! 	 & chanprof,        			& ! out
! 	 & frequencies)       		      ! out

!     ! allocate profiles (input) and radiance (output) structures
!     asw = 1
!     call rttov_alloc_prof( errorstatus,prf_num_in,profiles,nlevels_in,opts,asw, init = .true.)
!     call rttov_alloc_scatt_prof(prf_num_in,cld_profiles, nlevels_in, .false., 1, init = .true.)
!     call rttov_alloc_rad(errorstatus,nchanxnprof,radiance,nlevels_in-1,asw)

!     ! fill the profile structures with data
!     do j = 1, prf_num_in
!        profiles(j)%latitude    = latitude(j)
!        profiles(j)%longitude   = longitude(j)
!        profiles(j)%elevation   = h_surf(j)
!        profiles(j)%sunzenangle = 0.0 ! hard-coded in rttov9 int
!        profiles(j)%sunazangle  = 0.0 ! hard-coded in rttov9 int
!        profiles(j)%azangle     = 0.0
!        profiles(j)%zenangle    = zenang_in
!        profiles(j)%s2m%t       = t_surf(j)
!        profiles(j)%s2m%q       = q_surf(j)
!        profiles(j)%s2m%u       = u_surf(j)
!        profiles(j)%s2m%v       = v_surf(j)
!        profiles(j)%s2m%wfetc   = 10000.
!        profiles(j)%skin%t      = t_skin(j)
!        profiles(j)%skin%watertype  = 1 ! ocean water
!        if (lsmask(j) < 0.5) then
!           profiles(j)%skin%surftype  = surftype_sea
!        else
!           profiles(j)%skin%surftype  = surftype_land
!        endif
!        if (seaice(j) >= 0.5) then
!           profiles(j)%skin%surftype  = surftype_seaice
!        endif
!        profiles(j)%skin%fastem(1)  =  3.0
!        profiles(j)%skin%fastem(2)  =  5.0
!        profiles(j)%skin%fastem(3)  =  15.0
!        profiles(j)%skin%fastem(4)  =  0.1
!        profiles(j)%skin%fastem(5)  =  0.3
!        profiles(j)%cfraction       = 0.0
!        profiles(j)%ctp 	           = 500.0 ! not used but still required by rttov
!        profiles(j)%p(:)            = p_in(j,:)
!        profiles(j)%t(:)            = t_in(j,:)
!        profiles(j)%q(:)	           = q_in(j,:)
!        profiles(j)%idg             = 0.
!        profiles(j)%ish             = 0.
!        where(profiles(j)%q(:) < 1e-4)
!           profiles(j)%q(:) = 1e-4
!        end where
!        cld_profiles(j)%ph(:)   = ph_in(j,:)
!        cld_profiles(j)%cc(:)   = cc_in(j,:)
!        cld_profiles(j)%clw(:)  = clw_in(j,:)
!        cld_profiles(j)%ciw(:)  = ciw_in(j,:)
!        cld_profiles(j)%rain(:) = rain_in(j,:)
!        cld_profiles(j)%sp(:)   = sp_in(j,:)
!        profiles(j)%s2m%p       = cld_profiles(j)%ph(nlevels_in+1)
!     enddo

!     call rttov_get_emis(  &
!          & errorstatus, &
!          & opts,       &
!          & chanprof,   &
!          & profiles,   &
!          & coef_rttov,      &
!          !               & resolution=resolution, &  ! *** MW atlas native resolution is
!          !                 0.25 degree lat/lon; if you know better value for satellite
!          !                 footprint (larger than this) then you can specify it here
!          & emissivity=emissivity(:)%emis_in)
!     if (errorstatus /= errorstatus_success) then
!        write(*,*) 'In COSP_RTTOV11: Error RTTOV_GET_EMIS!'
!        call rttov_exit(errorstatus)
!     endif

!     calcemis(:) = .false.
!     where (emissivity(:)%emis_in <= 0.)
!        calcemis(:) = .true.
!     endwhere

!     call rttov_scatt (&
!          & errorstatus,         &! out
!          & opts_scatt,          &! in
!          & nlevels_in,          &! in
!          & chanprof,            &! in
!          & frequencies,         &! in
!          & profiles,            &! in
!          & cld_profiles,        &! in
!          & coef_rttov,          &! in
!          & coef_scatt,          &! in
!          & calcemis,           &! in
!          & emissivity,          &! in
!          & radiance)             ! out

!     if (errorstatus /= errorstatus_success) then
!        write(*,*) 'In COSP_RTTOV11: Error RTTOV_SCATT!'
!        call rttov_exit(errorstatus)
!     endif

!     !write(*,*) 'Checking emissivities: ', maxval(emissivity(:)%emis_out), \
!     !         minval(emissivity(:)%emis_out)
!     tbs(1:prf_num_in,1:1+size(ichan(:,1))-1) = &
!          transpose(reshape(radiance%bt(1:nchanxnprof),(/ size(ichan(:,1)),prf_num_in/) ))

!     ! deallocate all storage
!     asw = 0
!     ! 	call rttov_dealloc_coefs(errorstatus,coef_rttov)
!     ! 	call rttov_dealloc_scattcoeffs(coef_scatt)
!     call rttov_alloc_prof(errorstatus,prf_num_in,profiles,nlevels_in,opts,asw)
!     call rttov_alloc_scatt_prof(prf_num_in,cld_profiles,nlevels_in,.false.,asw)
!     call rttov_alloc_rad(errorstatus,nchanxnprof,radiance,nlevels_in-1,asw)
!     deallocate(ichan,chanprof,frequencies,emissivity,calcemis)  !instrument,
!     !***************************************************************************
!     !-------- end section --------
!     !***************************************************************************
!   end subroutine cosp_rttov_mwscatt
!   function construct_rttov_coeffilename(platform,satellite,instrument)
!     ! Inputs
!     integer,intent(in) :: platform,satellite,instrument
!     ! Outputs
!     character(len=256) :: construct_rttov_coeffilename
!     ! Local variables
!     character(len=256) :: coef_file
!     integer :: error

!     ! Initialize
!     error = 0
    
!     ! Platform
!     if (platform .eq. 1)  coef_file = 'rtcoef_noaa_'
!     if (platform .eq. 10) coef_file = 'rtcoef_metop_'
!     if (platform .eq. 11) coef_file = 'rtcoef_envisat_'
!     if (platform .ne. 1 .and. platform .ne. 10 .and. platform .ne. 11) then
!        error=error+1
!        write ( *,* ) 'Unsupported platform ID ',platform
!        return
!     endif

!     ! Satellite
!     if (satellite .lt. 10) then
!        coef_file = trim(coef_file) // char(satellite+48)
!     else if (satellite .lt. 100) then
!        coef_file = trim(coef_file) // char(int(satellite/10)+48)
!        coef_file = trim(coef_file) // char(satellite-int(satellite/10)*10+48)
!     else
!        error=error+1
!        write ( *,* ) 'Unsupported satellite number ',satellite
!        return
!     endif

!     ! Sensor
!     if (sensor .eq. 3)  coef_file = trim(coef_file) // '_amsua.dat'
!     if (sensor .eq. 5)  coef_file = trim(coef_file) // '_avhrr.dat'
!     if (sensor .eq. 49) coef_file = trim(coef_file) // '_mwr.dat'
!     if (sensor .ne. 3 .and. sensor .ne. 5 .and. sensor .ne. 49) then
!        error=error+1
!        write ( *,* ) 'Unsupported sensor number ', sensor
!        return
!     endif

!     if (error .eq. 0) construct_rttov_coeffilename=coef_file
    
!   end function construct_rttov_coeffilename
!   function construct_rttov_scatfilename(platform,satellite,instrument)
!     ! Inputs
!     integer,intent(in) :: platform,satellite,instrument
!     ! Outputs
!     character(len=256) :: construct_rttov_scatfilename
!     ! Local variables
!     character(len=256) :: coef_file
!     integer :: error

!     ! Initialize
!     error = 0
    
!     ! Platform
!     if (platform .eq. 1)  coef_file = 'sccldcoef_noaa_'
!     if (platform .eq. 10) coef_file = 'sccldcoef_metop_'
!     if (platform .eq. 11) coef_file = 'sccldcoef_envisat_'
!     if (platform .ne. 1 .and. platform .ne. 10 .and. platform .ne. 11) then
!        error=error+1
!        write ( *,* ) 'Unsupported platform ID ',platform
!        return
!     endif

!     ! Satellite
!     if (satellite .lt. 10) then
!        coef_file = trim(coef_file) // char(satellite+48)
!     else if (satellite .lt. 100) then
!        coef_file = trim(coef_file) // char(int(satellite/10)+48)
!        coef_file = trim(coef_file) // char(satellite-int(satellite/10)*10+48)
!     else
!        error=error+1
!        write ( *,* ) 'Unsupported satellite number ',satellite
!        return
!     endif

!     ! Sensor
!     if (sensor .eq. 3)  coef_file = trim(coef_file) // '_amsua.dat'
!     if (sensor .eq. 5)  coef_file = trim(coef_file) // '_avhrr.dat'
!     if (sensor .eq. 49) coef_file = trim(coef_file) // '_mwr.dat'
!     if (sensor .ne. 3 .and. sensor .ne. 5 .and. sensor .ne. 49) then
!        error=error+1
!        write ( *,* ) 'Unsupported sensor number ', sensor
!        return
!     endif

!     if (error .eq. 0) construct_rttov_scatfilename=coef_file
    
!   end function construct_rttov_scatfilename
  
end module mod_cosp_rttov
