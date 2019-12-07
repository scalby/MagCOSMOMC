    !Use CAMB
    module Calculator_CAMB
    use CosmologyTypes
    use CosmoTheory
    use CAMB, only : CAMB_GetResults, CAMB_GetAge, CAMBParams, CAMB_SetDefParams, &
        AccuracyBoost,  Cl_scalar, Cl_tensor, Cl_lensed, outNone, w_lam, wa_ppf,&
        CAMBParams_Set, MT, CAMBdata, NonLinear_Pk, Nonlinear_lens, Reionization_GetOptDepth, CAMB_GetZreFromTau, &
        CAMB_GetTransfers,CAMB_FreeCAMBdata,CAMB_InitCAMBdata, CAMB_TransfersToPowers, Transfer_SetForNonlinearLensing, &
        initial_adiabatic,initial_vector,initial_iso_baryon,initial_iso_CDM, initial_iso_neutrino, initial_iso_neutrino_vel, &
        HighAccuracyDefault, highL_unlensed_cl_template, ThermoDerivedParams, nthermo_derived, BackgroundOutputs, &
        Transfer_SortAndIndexRedshifts,  &
        Recombination_Name, reionization_name, power_name, threadnum, version, tensor_param_rpivot 
    use Errors !CAMB
    use settings
    use likelihood
    use Calculator_Cosmology
    use GeneralTypes
!-----
! MagCosmoMC
!   use, intrinsic :: ieee_arithmetic
    use constants
!-----
    implicit none
    private
!-----
! MagCosmoMC:
! adding definitions for magnetic modes.
integer :: MagneticMode ! 0 = no magnetic, 1 = compensated, 2 = passive
integer, parameter :: scalar_comp_dd = 1, scalar_comp_dp = 2, scalar_comp_pp = 3, &
		      vector_comp = 4, scalar_tensor_passive = 5, no_mag = 0

!adding other physical variables
real(mcp) :: Rg, Rv, lrat

!------

    Type, extends(TTheoryIntermediateCache) :: CAMBTransferCache
       Type (CAMBdata) :: Transfers !this can be used even for passive magnetic modes (scalar and tensor).
!-----------
! MAgCosmoMC: other CAMB data needed to store the different transfer functions
	Type (CAMBdata) :: VectorTransfers
	Type (CAMBdata) :: ScalarCompTransfersDD
	Type (CAMBdata) :: ScalarCompTransfersPP
	Type (CAMBdata) :: ScalarCompTransfersDP
	!Type (CAMBdata) :: TensorCompTransfers
!-----------
    contains
    procedure :: Clear => CAMBTransferCache_Clear
    end Type CAMBTransferCache

    Type, extends(TCosmologyCalculator) :: CAMB_Calculator
        integer :: ncalls = 0
        integer :: nerrors = 0
        logical :: CAMB_timing = .false.
        real(mcp) :: k_eta_max_scalar = -1._mcp
        logical :: accurate_BB =.false.
        type(CAMBParams)  CAMBP
        character(LEN=:), allocatable :: highL_theory_cl_template_file
        real(mcp), allocatable :: highL_lensedCL_template(:,:)
    contains
    !New
    procedure :: CMBToCAMB => CAMBCalc_CMBToCAMB
    procedure :: SetDerived => CAMBCalc_SetDerived
    procedure :: SetPowersFromCAMB => CAMBCalc_SetPowersFromCAMB
    procedure :: InitCAMB => CAMBCalc_InitCAMB
    procedure :: InitCAMBParams => CAMBCalc_InitCAMBParams
    procedure :: SetCAMBInitPower => CAMBCalc_SetCAMBInitPower
    procedure :: SetPkFromCAMB => CAMBCalc_SetPkFromCAMB
    procedure :: GetNLandRatios => CAMBCalc_GetNLandRatios
    !Overridden inherited
    procedure :: ReadParams => CAMBCalc_ReadParams
    procedure :: InitForLikelihoods => CAMBCalc_InitForLikelihoods
    procedure :: BAO_D_v => CAMBCalc_BAO_D_v
    procedure :: AngularDiameterDistance => CAMBCalc_AngularDiameterDistance
    procedure :: ComovingRadialDistance => CAMBCalc_ComovingRadialDistance
    procedure :: AngularDiameterDistance2 => CAMBCalc_AngularDiameterDistance2
    procedure :: LuminosityDistance => CAMBCalc_LuminosityDistance
    procedure :: Hofz => CAMBCalc_Hofz
    procedure :: CMBToTheta => CAMBCalc_CMBToTheta
    procedure :: GetNewPowerData => CAMBCalc_GetNewPowerData
    procedure :: GetNewTransferData => CAMBCalc_GetNewTransferData
    procedure :: GetCosmoTheoryForImportance => CAMBCalc_GetTheoryForImportance
    procedure :: GetZreFromTau  => CAMBCalc_GetZreFromTau
    procedure :: GetOpticalDepth => CAMBCalc_GetOpticalDepth
    procedure :: SetBackgroundTheoryData => CAMBCalc_SetBackgroundTheoryData
    procedure :: SetParamsForBackground => CAMBCalc_SetParamsForBackground
    procedure :: VersionTraceOutput => CAMBCalc_VersionTraceOutput
    procedure, private :: LoadFiducialHighLTemplate
    end type CAMB_Calculator


    !    integer, parameter :: ScalClOrder(5) = (/1,3,2,4,5/), TensClOrder(4) = (/1,4,2,3/)
    !Mapping of CAMB CL array ordering to TT , TE, EE, BB, phi, phiT


    public CAMB_Calculator
    contains
    !-----
! MagCosmoMC
   !subroutine CAMBCalc_CMBToCAMB(this,CMB,P)
    subroutine CAMBCalc_CMBToCAMB(this,CMB,P, ScComP, VecP)
    use LambdaGeneral
    use CAMBmain, only : ALens
    use constants, only : default_nnu
    use lensing, only : ALens_Fiducial
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    type(CAMBParams)  P
!-----
! MagCosmoMC: adding parameters for Vector modes.
type(CAMBParams), optional :: VecP
type(CAMBParams), optional :: ScComP

! when I'll add the negligible tensor compensated modes..
!type(CAMBParams), optional :: TenComP
!-----
    real(dl) neff_massive_standard

    P = this%CAMBP
    P%omegab = CMB%omb
    P%omegan = CMB%omnu
    P%omegac = CMB%omc
    P%omegav = CMB%omv
    P%H0 = CMB%H0
    P%Reion%redshift= CMB%zre
    P%Reion%delta_redshift = CMB%zre_delta
    w_lam = CMB%w
!    wa_ppf = CMB%wa 
!magcosmomc
    ALens = CMB%ALens
    ALens_Fiducial = CMB%ALensf
    P%InitialConditionVector(initial_iso_CDM) = &
        sign(sqrt(abs(CMB%iso_cdm_correlated) /(1-abs(CMB%iso_cdm_correlated))),CMB%iso_cdm_correlated)
    P%Num_Nu_Massive = 0
    P%Nu_mass_numbers = 0
    P%Num_Nu_Massless = CMB%nnu
    if (CMB%omnuh2>0) then
        P%Nu_mass_eigenstates=0
        if (CMB%omnuh2>CMB%omnuh2_sterile) then
            neff_massive_standard = CosmoSettings%num_massive_neutrinos*default_nnu/3
            P%Num_Nu_Massive = CosmoSettings%num_massive_neutrinos
            P%Nu_mass_eigenstates=P%Nu_mass_eigenstates+1
            if (CMB%nnu > neff_massive_standard) then
                P%Num_Nu_Massless = CMB%nnu - neff_massive_standard
            else
                P%Num_Nu_Massless = 0
                neff_massive_standard=CMB%nnu
            end if
            P%Nu_mass_numbers(P%Nu_mass_eigenstates) = CosmoSettings%num_massive_neutrinos
            P%Nu_mass_degeneracies(P%Nu_mass_eigenstates) = neff_massive_standard
            P%Nu_mass_fractions(P%Nu_mass_eigenstates) = (CMB%omnuh2-CMB%omnuh2_sterile)/CMB%omnuh2
        else
            neff_massive_standard=0
        end if
        if (CMB%omnuh2_sterile>0) then
            if (CMB%nnu<default_nnu) call MpiStop('nnu < 3.046 with massive sterile')
            P%Num_Nu_Massless = default_nnu - neff_massive_standard
            P%Num_Nu_Massive=P%Num_Nu_Massive+1
            P%Nu_mass_eigenstates=P%Nu_mass_eigenstates+1
            P%Nu_mass_numbers(P%Nu_mass_eigenstates) = 1
            P%Nu_mass_degeneracies(P%Nu_mass_eigenstates) = max(1d-6,CMB%nnu - default_nnu)
            P%Nu_mass_fractions(P%Nu_mass_eigenstates) = CMB%omnuh2_sterile/CMB%omnuh2
        end if
    end if

    P%YHe = CMB%YHe
#ifdef COSMOREC
    if (P%Recomb%fdm/=0._mcp) P%Recomb%runmode = 3
    P%Recomb%fdm = CMB%fdm * 1e-23_mcp
#else
    if (CMB%fdm/=0._mcp) call MpiStop('Compile with CosmoRec to use fdm')
#endif
    call this%SetCAMBInitPower(P,CMB,1)

!---------------------------
! MagCosmoMC: setting up rhe parameters for Scalar Compensated Modes
if(present(ScComP)) then
    ScComP = P
    ScComP%Scalar_initial_condition = 6
    ScComP%WantScalars = .true.
    ScComP%WantTensors = .false.
    ScComP%WantVectors = .false.
    ScComP%WantTransfer = .false.
    ScComP%OnlyTransfers = .true.
    ScComP%NonLinear = 0
    ScComP%MassiveNuMethod = 1
end if
!---------------------------
! MagCosmoMC: setting up the parameters for Vector Modes.
if(present(VecP)) then
    VecP = this%CAMBP
    VecP%omegab = CMB%omb
    VecP%omegan = 0.d0
! CAMB uses only massless neutrinos in vector modes. If there are massive neutrinos, give their density 
! to CDM. 
    VecP%omegac = CMB%omc + CMB%omnu
    VecP%omegav = CMB%omv
    VecP%H0 = CMB%H0
    VecP%Reion%redshift= CMB%zre
    VecP%Reion%delta_redshift = CMB%zre_delta
    VecP%Scalar_initial_condition = 1

    ! Working with massless neutrinos only for vector modes. I hope it will not create problems.
    VecP%Num_Nu_Massive = 0
    VecP%Nu_mass_numbers = 0
    VecP%Num_Nu_Massless = CMB%nnu

    VecP%WantScalars = .false.
    VecP%WantTensors = .false.
    VecP%WantVectors = .true.
    VecP%WantCls = .true.
    VecP%WantTransfer = .false.
    VecP%OnlyTransfers = .true.
    P%YHe = CMB%YHe
#ifdef COSMOREC
    if (VecP%Recomb%fdm/=0._mcp) VecP%Recomb%runmode = 3
    VecP%Recomb%fdm = CMB%fdm * 1e-23_mcp
#else
    if (CMB%fdm/=0._mcp) call MpiStop('Compile with CosmoRec to use fdm')
#endif
!---------------------------

!I need to add the part for tensor compensated... (in the future)
end if


    end subroutine CAMBCalc_CMBToCAMB

    subroutine CAMBCalc_SetDerived(this,Theory)
    class(CAMB_Calculator) :: this
    Class(TCosmoTheoryPredictions) Theory
    integer noutputs, i

    noutputs = size(BackgroundOutputs%z_outputs)
    Theory%numderived = nthermo_derived + noutputs*4
    if (Theory%numderived > max_derived_parameters) &
        call MpiStop('numderived > max_derived_parameters: increase in CosmologyTypes.f90')
    Theory%derived_parameters(1:nthermo_derived) = ThermoDerivedParams(1:nthermo_derived)
    do i=1, noutputs
        Theory%derived_parameters(nthermo_derived+(i-1)*4+1) = BackgroundOutputs%rs_by_D_v(i)
        Theory%derived_parameters(nthermo_derived+(i-1)*4+2) = BackgroundOutputs%H(i)*const_c/1e3_mcp
        Theory%derived_parameters(nthermo_derived+(i-1)*4+3) = BackgroundOutputs%DA(i)
        Theory%derived_parameters(nthermo_derived+(i-1)*4+4) = (1+BackgroundOutputs%z_outputs(i))* &
            BackgroundOutputs%DA(i) * BackgroundOutputs%H(i) !F_AP parameter
    end do
    end subroutine CAMBCalc_SetDerived

    subroutine CAMBCalc_SetParamsForBackground(this,CMB)
    use Camb, only: CAMBParams_Set
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    Type(CAMBParams)  P

    !set background parameters, but don't calculate thermal history
    call this%CMBToCAMB(CMB, P)
    call CAMBParams_Set(P)

    end subroutine CAMBCalc_SetParamsForBackground

    subroutine CAMBCalc_SetBackgroundTheoryData(this, CMB,Theory,error)
    use cambmain, only: initvars
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    Class(TCosmoTheoryPredictions) Theory
    integer error

    call InitVars !calculate thermal history, e.g. z_drag etc.
    if (global_error_flag/=0) then
        error=global_error_flag
        return
    end if
    call this%SetDerived(Theory)

    end subroutine CAMBCalc_SetBackgroundTheoryData

    subroutine CAMBCalc_GetNewTransferData(this, CMB,Info,Theory,error)
!-------
! MagCosmoMC: added this modules:
!use InitialPower 
use ModelParams, only : tens0, use_spline_template !I need to avoid using spline_template in magnetic modes.
use magnetic !ALEX
!--------
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error
    type(CAMBParams)  P
!------
! MagCosmoMC: I need at least one more structure for the MagCAMB parameters.
Type(CAMBParams) :: VecP
Type(CAMBParams) :: ScComP
! I might need the Tensor Compensated mode..
!------
    real(mcp) time


    allocate(CAMBTransferCache::Info)

    select type (Info)
    class is (CAMBTransferCache)
        call CAMB_InitCAMBdata(Info%Transfers)
	if(CosmoSettings%vector_compensated) call CAMB_InitCAMBdata(Info%VectorTransfers)
	if(CosmoSettings%scalar_compensated) then
		call CAMB_InitCAMBdata(Info%ScalarCompTransfersDD)
		call CAMB_InitCAMBdata(Info%ScalarCompTransfersPP)
		call CAMB_InitCAMBdata(Info%ScalarCompTransfersDP)
	end if
	!if(CosmoSettings%tensor_compensated) call CAMB_InitCAMBdata(Info%TensorCompTransfers)
	
        !call CAMB_InitCAMBdata(Info%ScalPassTrans)!ALEX
    if((.NOT.(CosmoSettings%scalar_compensated)) .AND. &
	   (.NOT.(CosmoSettings%vector_compensated)))then
		call this%CMBToCAMB(CMB, P)!, VecP, ScComP)
	else if (CosmoSettings%scalar_compensated .AND.&
	        (.NOT.(CosmoSettings%vector_compensated))) then
		call this%CMBToCAMB(CMB, P, ScComP=ScComP)
	else if ((.NOT.(CosmoSettings%scalar_compensated))&
		 .AND. CosmoSettings%vector_compensated) then
		call this%CMBToCAMB(CMB, P, VecP=VecP)
	else if (CosmoSettings%scalar_compensated .AND. CosmoSettings%vector_compensated) then
		call this%CMBToCAMB(CMB, P, ScComP, VecP)
	end if

!----------- Primary and Passive Modes -------------------
! MagCosmoMC: Transfer functions for primary and passive modes are the same.
!             Using the same for the two modes. 

   if(CosmoSettings%compute_tensors) tens0 = 1._dl

if (CosmoSettings%scalar_passive .OR. CosmoSettings%tensor_passive) then
! Magnetic fields must be initialized. 
   !call mag_set_alfven(CMB%magamp, CMB%magind)
   pib=0._dl
   delb=0._dl
    ! this is for tensors passive magnetic....

   if (CosmoSettings%tensor_passive)   tens0 = 1._dl
end if

        if (Feedback > 1) write(*,*) 'Calling CAMB: Primary and Passive Transfer functions.'

            Threadnum =num_threads

            if (this%CAMB_timing) time = TimerTime()	

            call CAMB_GetTransfers(P, Info%Transfers, error)

            if (this%CAMB_timing) call Timer('GetTransfers', time)

!**************** COMPENSATED MODES ************************
! compensated modes require different transfer functions.
!-----SCALAR: 3 types of transfer functions ----------------
if(CosmoSettings%scalar_compensated) then

    tens0 = 0._dl

! Delta-Delta
    delb = 1._dl
    pib  = 0._dl
!    if (Feedback>1)
    !write(*,*) "Scalar Compensated Modes: Delta-Delta"
    call CAMB_GetTransfers(ScComP, Info%ScalarCompTransfersDD, error)

! Pi-Pi
   delb = 0._dl
   pib  = 1._dl
!    if (Feedback>1) 
     !write(*,*) "Scalar Compensated Modes: Pi-Pi"
   call CAMB_GetTransfers(ScComP, Info%ScalarCompTransfersPP, error)

! Delta-Pi
   delb = 1._dl
   pib  = 1._dl
!    if (Feedback>1) 
   !write(*,*) "Scalar Compensated Modes: Delta-Pi"
   call CAMB_GetTransfers(ScComP, Info%ScalarCompTransfersDP, error)
end if

!------VECTOR:-------------------------------------------------
if(CosmoSettings%vector_compensated) then
   if (Feedback > 1) write(*,*) 'Magnetic Vector Transfer Functions.'
   !if (Feedback > 1)
   ! write(*,*) 'Vector Compensated Transfer Functions'
   !setting some stuff for vector modes
   tens0=0._dl
   vec_sig0 = 0._dl
   delb = 0._dl
   pib= 1._dl
   call CAMB_GetTransfers(VecP, Info%VectorTransfers, error)
end if
!-----------------------------------------------------------
    class default
        call MpiStop('CAMB_GetNewTransferData with wrong TTheoryIntermediateCache type')

    end select
    if (error==0) then
        call this%SetDerived(Theory)
    else
        if (stop_on_error) call MpiStop('CAMB error '//trim(global_error_message))
        if (Feedback > 0) write(*,*) 'CAMB returned error '//trim(global_error_message)
        this%nerrors=this%nerrors+1
    end if
    this%ncalls=this%ncalls+1
    if (mod(this%ncalls,100)==0 .and. LogFile%Opened()) then
        write(LogFile%unit,'("CAMB called ",I0," times; ",I0," errors")') this%ncalls, this%nerrors
    end if
    if (Feedback > 1) write(*,*) 'CAMB done'
! added just this write statement
!    write(*,*) 'CAMB done'
    end subroutine CAMBCalc_GetNewTransferData


    subroutine CAMBCalc_GetNewPowerData(this, CMB, Info, Theory, error)
!--------
! MagCosmoMC: two modules needed
use ModelParams
use magnetic
!I can compute all the Cls within this subroutine
    class(CAMB_Calculator) :: this
    class(CMBParams) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error,i,j
!------------
! MagCosmoMC:
! I need to store the CAMB params somewhere..
Type(CAMBParams) :: Pbuffer
!adding some variables for power spectrum
real(mcp) :: amp, camb_ind
!------------

    select type (Info)
    class is (CAMBTransferCache)
!------
! MAgCosmoMC: storing the CAMB params
Pbuffer = Info%Transfers%Params
!------

!------- PRIMARY --------

MagneticMode = no_mag
write(*,*) "no_mag",MagneticMode !yun check    
        call this%SetCAMBInitPower(Info%Transfers%Params,CMB,1)
if (Feedback>1) write(*,*) "Transfer to powers PRIMARY"
        call CAMB_TransfersToPowers(Info%Transfers)
if (Feedback>1)	write(*,*) "done."
        !this sets slow CAMB params correctly from value stored in Transfers
        if (global_error_flag/=0) then
            error=global_error_flag
            return
        end if
        !JD 08/13 added so we dont have to fill Cls unless using CMB
        if(CosmoSettings%use_CMB)then
            call this%SetPowersFromCAMB(CMB,Theory)
            do i=1, min(3,CosmoSettings%num_cls)
                if (CosmoSettings%cl_lmax(i,i)>0) then
                    if (any(Theory%cls(i,i)%Cl(:) < 0 )) then
                        error = 1
			            write(*,*) "error Theory%cls(i,i)%Cl(:) < 0 " !yun check 
                        call MpiStop('Calculator_CAMB: negative PRIMARY C_l (could edit to silent error here)')
                        return
                    end if
                end if
                do j= i, 1, -1
                    if (CosmoSettings%cl_lmax(i,j)>0) then
                        !if ( any(ieee_is_nan(Theory%cls(i,j)%Cl))) then !magcosmomc 
                        if ( any(isNan(Theory%cls(i,j)%Cl))) then
                            error=1
                            write(*,*) 'WARNING: NaN CL?', i, j
                            return
                        end if
                    end if
                end do
            end do
        end if

!*************************** MAGNETIC CLS ***************************
if(CosmoSettings%scalar_compensated .OR. CosmoSettings%scalar_passive .OR. &
   CosmoSettings%vector_compensated .OR. CosmoSettings%tensor_passive) then


	! use_spline_template set to false to avoid the interpolation error.
	use_spline_template = .false.

	!Setup of some common variables
        Rg=1._dl / (1._dl+7._dl/8._dl *(4._dl/11._dl)**(4._dl/3._dl) * (Info%Transfers%Params%Num_Nu_Massless + Info%Transfers%Params%Num_Nu_Massive))
        Rv=1._dl - Rg
    	lrat = CMB%maglrat * log(10.d0) ! convert log10 in ln
end if

if(CosmoSettings%do_helical) then
            write(*,*) "Doing helical PMFs!"         
            if(CosmoSettings%maximal_hel) then
                CMB%helical_ind = CMB%magind
                CMB%helical_amp = CMB%magamp*sqrt(GAMMA((CMB%helical_ind+4.d0)/2.d0)/GAMMA((CMB%magind+3.d0)/2.d0))
            write(*,*) "Computing maximal helical contributions: "
            end if
end if
        
!------------ MAGNETIC COMPENSATED MODES ----------
if(CosmoSettings%scalar_compensated) then  ! Scalar  ! Delta-Delta  
    !write(*,*) "------ Scal-Comp Delta-Delta - Delta-Pi, The density perturbed mode (delb = 1, pib = 0)
    MagneticMode= scalar_comp_dd
    write(*,*) "scalar_comp_dd",MagneticMode !yun check 
    !setting parameters
    call this%SetCAMBInitPower(Info%ScalarCompTransfersDD%Params,CMB,1)
    ! This stuff must be initialized properly since I've changed the procedure to compute the P.S.
    Info%ScalarCompTransfersDD%Params%InitPower%nn = 1
    Info%ScalarCompTransfersDD%Params%InitPower%rat(1) = 1._dl   
    Info%ScalarCompTransfersDD%Params%InitPower%n_run(1) = 0._dl
    Info%ScalarCompTransfersDD%Params%InitPower%nt_run(1) = 0._dl  !tensor 
    Info%ScalarCompTransfersDD%Params%InitPower%k_0_scalar = 2*pi
    Info%ScalarCompTransfersDD%Params%InitPower%k_0_tensor = 2*pi

!Setting up k_D
    Info%ScalarCompTransfersDD%Params%InitPower%kDissip = (5.5d0 * 10.d0**4)**(1.d0/(CMB%magind + 5.d0)) * (CMB%magamp)**(-2.d0/(CMB%magind+5.d0)) * &
        (2*pi)**((CMB%magind+3)/(CMB%magind+5))*(Info%ScalarCompTransfersDD%Params%H0/100)**(1/(CMB%magind+5))&
        *((Info%ScalarCompTransfersDD%Params%omegab*(Info%ScalarCompTransfersDD%Params%H0/100)**2)/0.022)**(1/(CMB%magind+5))


    if(CMB%magind .ge. -2) then !n >= -2
!write(*,*) "using FITTING Functions" 
        Info%ScalarCompTransfersDD%Params%InitPower%an(1) = CMB%magind
        Info%ScalarCompTransfersDD%Params%InitPower%CorrType = 1
        if (Feedback>1) write(*,*) "Using FITTING Functions"
        Info%ScalarCompTransfersDD%Params%InitPower%ScalarPowerAmp(1)=mag_amplitude(CMB%magind, CMB%magamp)
	    write(*,*) "Using FITTING Functions, Scal-Comp Delta-Delta - Delta-Pi", Info%ScalarCompTransfersDD%Params%InitPower%CorrType
    else !-3<n<-2
!	write(*,*) "using Interpolation TABLE"
        Info%ScalarCompTransfersDD%Params%InitPower%an(1) = 1._dl + 2._dl*(3 + CMB%magind)
        Info%ScalarCompTransfersDD%Params%InitPower%CorrType = 0
        Info%ScalarCompTransfersDD%Params%InitPower%ScalarPowerAmp(1)= mag_psamp(CMB%magind, CMB%magamp, 1) - mag_psamp(CMB%magind, CMB%magamp, 2)        
        !if (Feedback>1) write(*,*) "Using interpolation TABLE"
        write(*,*) "Using interpolation TABLE,Scal-Comp Delta-Delta - Delta-Pi"
    endif
 !------helical yun 
    if(CosmoSettings%do_helical) then
       if(CMB%helical_ind < -2.d0) then
         write(*,*) "Helical Using interpolation TABLE,Scal-Comp Delta-Delta - Delta-Pi"
         Info%ScalarCompTransfersDD%Params%InitPower%CorrType_hel = 0
         Info%ScalarCompTransfersDD%Params%InitPower%ScalarPowerAmp_hel(1)=mag_psamp_hel(CMB%helical_ind, CMB%helical_amp,6) -&
                                                       mag_psamp_hel(CMB%helical_ind, CMB%helical_amp,7)
         Info%ScalarCompTransfersDD%Params%InitPower%an_hel1(1) =  1.d0 + 2.d0*(CMB%helical_ind +3.d0)
       else  !Use first fitting formula
         Info%ScalarCompTransfersDD%Params%InitPower%an_hel1(1) = CMB%helical_ind
         Info%ScalarCompTransfersDD%Params%InitPower%ScalarPowerAmp_hel(1) = mag_amplitude_hel(CMB%helical_ind, CMB%helical_amp)
         Info%ScalarCompTransfersDD%Params%InitPower%CorrType_hel = 1
         write(*,*) "USING FITTING FORMULAS for helical P_ADD-P_ADPi" 
       end if                
    end if        
!-----      
    delb = 1._dl ! Set up perturbations
    pib = 0._dl

    Info%ScalarCompTransfersDD%Params%InitPower%ant(1) = 0.d0
    Info%ScalarCompTransfersDD%ClTransScal%NumSources = 2
	Info%ScalarCompTransfersDD%Params%DoLensing = .false.

    if (Feedback>1) write(*,*) "Transfer to powers  Delta Delta.."
	call CAMB_TransfersToPowers(Info%ScalarCompTransfersDD)
	if(Feedback>1) write(*,*) "..done."

	if (global_error_flag/=0) then
        error=global_error_flag
        return
    end if

	if(CosmoSettings%use_CMB)then
        call this%SetPowersFromCAMB(CMB,Theory) ! need to modify SetPowersFromCAMB...
        do i=1, min(3,CosmoSettings%num_cls)
            if (CosmoSettings%cl_lmax(i,i)>0) then
                if (any(Theory%MagScalCompClsDD(i,i)%Cl(:) < 0 )) then
                    !error = 1
		    !write(*,*) "CMB = ", CMB
                    !call MpiStop('Calculator_CAMB: negative ScalComp DD C_l (could edit to silent error here)')
                    !return
                end if
            end if
            do j= i, 1, -1
                if (CosmoSettings%cl_lmax(i,j)>0) then
                    if ( any(isNan(Theory%MagScalCompClsDD(i,j)%Cl))) then
                        error=1
                        write(*,*) 'WARNING: NaN CL?', i, j
                        return
                    end if
                end if
            end do
        end do
    end if
	

MagneticMode = scalar_comp_pp
   write(*,*) "scalar_comp_pp",MagneticMode !yun check 
! pp, The stress perturbed mode (delb = 0, pib = 1), write(*,*) "Pi-Pi - Delta-Pi"
 call this%SetCAMBInitPower(Info%ScalarCompTransfersPP%Params,CMB,1)
    Info%ScalarCompTransfersPP%Params%InitPower%nn = 1
    Info%ScalarCompTransfersPP%Params%InitPower%rat(1) = 1._dl  
    Info%ScalarCompTransfersPP%Params%InitPower%n_run(1) = 0._dl
    Info%ScalarCompTransfersPP%Params%InitPower%nt_run(1) = 0._dl!tensor
    Info%ScalarCompTransfersPP%Params%InitPower%k_0_scalar = 2*pi
    Info%ScalarCompTransfersPP%Params%InitPower%k_0_tensor = 2*pi
    Info%ScalarCompTransfersPP%Params%DoLensing = .false.
    Info%ScalarCompTransfersPP%ClTransScal%NumSources = 2

!Setting up k_D
        Info%ScalarCompTransfersPP%Params%InitPower%kDissip = (5.5d0 * 10.d0**4)**(1.d0/(CMB%magind + 5.d0)) * (CMB%magamp)**(-2.d0/(CMB%magind+5.d0)) * &
        (2*pi)**((CMB%magind+3)/(CMB%magind+5))*(Info%ScalarCompTransfersPP%Params%H0/100)**(1/(CMB%magind+5))&
        *((Info%ScalarCompTransfersPP%Params%omegab*(Info%ScalarCompTransfersPP%Params%H0/100)**2)/0.022)**(1/(CMB%magind+5))

    if(CMB%magind .ge. -2) then !n>=-2
        write(*,*) "Using FITTING Functions, P_PiPi-P_DPi"
        Info%ScalarCompTransfersPP%Params%InitPower%CorrType = 2
        if (Feedback>1) write(*,*) "Using FITTING Functions"
        Info%ScalarCompTransfersPP%Params%InitPower%an(1) =  CMB%magind
        Info%ScalarCompTransfersPP%Params%InitPower%ScalarPowerAmp(1)=mag_amplitude(CMB%magind, CMB%magamp)
    else
        write(*,*) "Using interpolation TABLE, P_PiPi-P_DPi"
        Info%ScalarCompTransfersPP%Params%InitPower%CorrType = 0
        if(Feedback>1) write(*,*) "Using interpolation TABLE"
        Info%ScalarCompTransfersPP%Params%InitPower%an(1) = 1._dl + 2._dl*(3 + CMB%magind)
        Info%ScalarCompTransfersPP%Params%InitPower%ScalarPowerAmp(1)= mag_psamp(CMB%magind, CMB%magamp, 3) - mag_psamp(CMB%magind, CMB%magamp, 2)
    endif
    delb = 0._dl ! Set up perturbations
    pib = 1._dl
!---helical ---yun
    if(CosmoSettings%do_helical) then
      if(CMB%helical_ind < -2.d0) then
           write(*,*) "Using interpolation TABLE for helical P_APiPi-P_ADPi part"
           Info%ScalarCompTransfersPP%Params%InitPower%ScalarPowerAmp_hel(1) = mag_psamp_hel(CMB%helical_ind, CMB%helical_amp,8) -&
                                                       mag_psamp_hel(CMB%helical_ind, CMB%helical_amp,7)
           Info%ScalarCompTransfersDD%Params%InitPower%CorrType_hel = 0
           Info%ScalarCompTransfersDD%Params%InitPower%an_hel1(1) =  1.d0 + 2.d0*(CMB%helical_ind +3.d0)
      else  !Use first fitting formula
           write(*,*) "USING FITTING FORMULAS for helical P_APiPi-P_ADPi part"
           Info%ScalarCompTransfersPP%Params%InitPower%CorrType_hel = 2
           Info%ScalarCompTransfersDD%Params%InitPower%an_hel1(1) = CMB%helical_ind
           Info%ScalarCompTransfersDD%Params%InitPower%ScalarPowerAmp_hel(1) = mag_amplitude_hel(CMB%helical_ind, CMB%helical_amp)
       end if               
     end if
!-------

	if(Feedback>1) write(*,*) "transfer to powers Pi-Pi"
        call CAMB_TransfersToPowers(Info%ScalarCompTransfersPP)
	    if(Feedback>1) write(*,*) "done."

        if (global_error_flag/=0) then
                error=global_error_flag
                return
        end if

        if(CosmoSettings%use_CMB)then
       	   call this%SetPowersFromCAMB(CMB,Theory) ! need to modify SetPowersFromCAMB...
           do i=1, min(3,CosmoSettings%num_cls)
              if (CosmoSettings%cl_lmax(i,i)>0) then
                  if (any(Theory%MagScalCompClsPP(i,i)%Cl(:) < 0 )) then
                    !error = 1
		    !write(*,*) "CMB = ", CMB
                    !call MpiStop('Calculator_CAMB: negative ScalComp PP C_l (could edit to silent error here)')
                    !return
                  end if
              end if
              do j= i, 1, -1
                if (CosmoSettings%cl_lmax(i,j)>0) then
                    if ( any(isNan(Theory%MagScalCompClsPP(i,j)%Cl))) then
                        error=1
                        write(*,*) 'WARNING: NaN CL?', i, j
                        return
                    end if
                end if
              end do
           end do
        end if

end if 

!dp, write(*,*) "------ Scal-Comp Delta-Pi -------"!! The combined perturbation (delb = 1, pib = 1)
MagneticMode = scalar_comp_dp
    write(*,*) "scalar_comp_dp",MagneticMode !yun check 
    !Set common stuff.
    call this%SetCAMBInitPower(Info%ScalarCompTransfersDP%Params,CMB,1)
    Info%ScalarCompTransfersDP%Params%InitPower%nn = 1
    Info%ScalarCompTransfersDP%Params%InitPower%rat(1) = 1._dl 
    Info%ScalarCompTransfersDP%Params%InitPower%n_run(1) = 0._dl
    Info%ScalarCompTransfersDP%Params%InitPower%nt_run(1) = 0._dl!tensor
    Info%ScalarCompTransfersDP%Params%InitPower%k_0_scalar = 2*pi
    Info%ScalarCompTransfersDP%Params%InitPower%k_0_tensor = 2*pi

!Setting up k_D
        Info%ScalarCompTransfersDP%Params%InitPower%kDissip = (5.5d0 * 10.d0**4)**(1.d0/(CMB%magind + 5.d0)) * (CMB%magamp)**(-2.d0/(CMB%magind+5.d0)) * &
        (2*pi)**((CMB%magind+3)/(CMB%magind+5))*(Info%ScalarCompTransfersDP%Params%H0/100)**(1/(CMB%magind+5))&
        *((Info%ScalarCompTransfersDP%Params%omegab*(Info%ScalarCompTransfersDP%Params%H0/100)**2)/0.022)**(1/(CMB%magind+5))

 
    if(CMB%magind .ge. -2) then !n>=-2
        write(*,*) "Using FITTING Functions, P_DP"
        Info%ScalarCompTransfersDP%Params%InitPower%an(1) =  CMB%magind
        Info%ScalarCompTransfersDP%Params%InitPower%CorrType = 3
        if(Feedback>1) write(*,*) "Using FITTING Functions"
        Info%ScalarCompTransfersDP%Params%InitPower%ScalarPowerAmp(1)= mag_amplitude(CMB%magind, CMB%magamp)* psconst(2)
    else
        write(*,*) "Using interpolation TABLE, P_DP"
        Info%ScalarCompTransfersDP%Params%InitPower%an(1) = 1._dl + 2._dl*(3 + CMB%magind)
        Info%ScalarCompTransfersDP%Params%InitPower%CorrType = 0
        if(Feedback>1) write(*,*) "Using interpolation TABLE"
        Info%ScalarCompTransfersDP%Params%InitPower%ScalarPowerAmp(1)= mag_psamp(CMB%magind, CMB%magamp, 2)
    endif
!----helical-----yun
    if(CosmoSettings%do_helical) then
      if(CMB%helical_ind < -2.d0) then
        write(*,*) "Using interpolation TABLE for helical P_ADP part"
        CP%InitPower%ScalarPowerAmp_hel(1) = mag_psamp_hel(CMB%helical_ind, CMB%helical_amp,7)
        Info%ScalarCompTransfersDD%Params%InitPower%CorrType_hel = 0
        Info%ScalarCompTransfersDD%Params%InitPower%an_hel1(1) =  1.d0 + 2.d0*(CMB%helical_ind +3.d0)             
      else
!Use first fitting formula
        write(*,*) "USING FITTING FORMULAS for helical P_ADP part"
        Info%ScalarCompTransfersDD%Params%InitPower%CorrType_hel = 3
        Info%ScalarCompTransfersDD%Params%InitPower%an_hel1(1) = CMB%helical_ind
        Info%ScalarCompTransfersDD%Params%InitPower%ScalarPowerAmp_hel(1) = mag_amplitude_hel(CMB%helical_ind, CMB%helical_amp)
      end if               
     end if
!------------   
    
delb = 1._dl ! Set up perturbations
pib = 1._dl

    Info%ScalarCompTransfersDP%Params%DoLensing = .false.
    Info%ScalarCompTransfersDP%ClTransScal%NumSources = 2

   if (Feedback>1) write(*,*) "transfer to powers Delta-Pi"
   call CAMB_TransfersToPowers(Info%ScalarCompTransfersDP)
   if (Feedback>1) write(*,*) "done."

        if (global_error_flag/=0) then
                error=global_error_flag
                return
        end if

    if(CosmoSettings%use_CMB)then
        call this%SetPowersFromCAMB(CMB,Theory) ! need to modify SetPowersFromCAMB...
        do i=1, min(3,CosmoSettings%num_cls)
            do j= i, 1, -1
                if (CosmoSettings%cl_lmax(i,j)>0) then
                    if ( any(isNan(Theory%MagScalCompClsDP(i,j)%Cl))) then
                        error=1
                        write(*,*) 'WARNING: NaN CL?', i, j
                        return
                    end if
                end if
            end do
        end do
    end if
!end scalar_compensated

if(CosmoSettings%vector_compensated) then
!VECTOR
MagneticMode = vector_comp
   write(*,*) "vector",MagneticMode !yun check 

    Rg=1._dl / (1._dl+7._dl/8._dl *(4._dl/11._dl)**(4._dl/3._dl) * (Info%Transfers%Params%Num_Nu_Massless + Info%Transfers%Params%Num_Nu_Massive))
    Rv=1._dl - Rg

    Info%VectorTransfers%Params%InitPower%nn = 1
    Info%VectorTransfers%Params%InitPower%rat(1) = 1._dl  
    Info%VectorTransfers%Params%InitPower%n_run(1) = 0._dl
    Info%VectorTransfers%Params%InitPower%nt_run(1) = 0._dl!tensor
    Info%VectorTransfers%Params%InitPower%k_0_scalar = 2*pi
    Info%VectorTransfers%Params%InitPower%k_0_tensor = 2*pi

!Setting up k_D
        Info%VectorTransfers%Params%InitPower%kDissip = (5.5d0 * 10.d0**4)**(1.d0/(CMB%magind + 5.d0)) * (CMB%magamp)**(-2.d0/(CMB%magind+5.d0)) * &
        (2*pi)**((CMB%magind+3)/(CMB%magind+5))*(Info%VEctorTransfers%Params%H0/100)**(1/(CMB%magind+5))&
        *((Info%VectorTransfers%Params%omegab*(Info%VectorTransfers%Params%H0/100)**2)/0.022)**(1/(CMB%magind+5))


    if(CMB%magind .ge. -2.d0) then
        !if (Feedbacklevel>1) 
        write(*,*) "USING FITTING FUNCTION Vector"
        Info%VectorTransfers%Params%InitPower%CorrType = 5 !VECTOR
        Info%VectorTransfers%Params%InitPower%ScalarPowerAmp(1)= mag_amplitude(CMB%magind, CMB%magamp)* psconst(4)
        Info%VectorTransfers%Params%InitPower%an(1) = CMB%magind
        !write(*,*) "P.S. Amplitude: ", amp
    else
        !if (Feedbacklevel>1) 
        write(*,*) "Using Table for integrals Vector"
        Info%VectorTransfers%Params%InitPower%CorrType = 0 ! use the approximated results.
        Info%VectorTransfers%Params%InitPower%an(1) = 1._dl + 2._dl*(3 + CMB%magind)
        Info%VectorTransfers%Params%InitPower%ScalarPowerAmp(1)= mag_psamp(CMB%magind, CMB%magamp, 4)
    end if
    
!-----helical----yun
  
    if(CosmoSettings%do_helical) then
       if(CMB%helical_ind < -2.d0) then !check
          write(*,*) "Using interpolation TABLE for helical P_Avpp part"
          Info%VectorTransfers%Params%InitPower%CorrType_hel = 0
          Info%VectorTransfers%Params%InitPower%an_hel1(1) = 1.d0 + 2.d0*(3.d0+CMB%helical_ind)
          Info%VectorTransfers%Params%InitPower%ScalarPowerAmp_hel(1) = mag_psamp_hel(CMB%helical_ind, CMB%helical_amp, 9)
       else
!Use first fitting formula
          write(*,*) "USING FITTING FORMULAS for helical P_Avpp part"
          Info%VectorTransfers%Params%InitPower%CorrType_hel = 5
          Info%VectorTransfers%Params%InitPower%an_hel1(1) = CMB%helical_ind
          Info%VectorTransfers%Params%InitPower%ScalarPowerAmp_hel(1) = mag_amplitude_hel(CMB%helical_ind, CMB%helical_amp)*psconst(9)

        end if                
     end if    
 !------------
    vec_sig0 = 0._dl
    delb = 0._dl ! Set up perturbations
    pib = 1._dl
    Info%VectorTransfers%Params%WantVectors = .true.
    Info%VectorTransfers%Params%DoLensing = .false.

    call CAMB_TransfersToPowers(Info%VectorTransfers)
   
    if (global_error_flag/=0) then
        error=global_error_flag
        return
    end if

     if(CosmoSettings%use_CMB)then
        call this%SetPowersFromCAMB(CMB,Theory) ! need to modify SetPowersFromCAMB...
        do i=1, min(3,CosmoSettings%num_cls)
            if (CosmoSettings%cl_lmax(i,i)>0) then
                if (any(Theory%MagVecCompCls(i,i)%Cl(:) < 0 )) then
                    !error = 1
		    !write(*,*) "CMB = ", CMB
                    !call MpiStop('Calculator_CAMB: negative VectComp C_l (could edit to silent error here)')
                    !return
                end if
            end if
            do j= i, 1, -1
                if (CosmoSettings%cl_lmax(i,j)>0) then
                    if ( any(isNan(Theory%MagVecCompCls(i,j)%Cl))) then
                        error=1
                        write(*,*) 'WARNING: NaN CL?', i, j
                        return
                    end if
                end if
            end do
        end do
    end if

end if !vector_compensated

!------------ PASSIVE -------------
if(CosmoSettings%scalar_passive .or. CosmoSettings%tensor_passive) then
    MagneticMode = scalar_tensor_passive
    write(*,*) "passive",MagneticMode !yun check 
      
    Rg=1._dl / (1._dl+7._dl/8._dl *(4._dl/11._dl)**(4._dl/3._dl) * (Info%Transfers%Params%Num_Nu_Massless + Info%Transfers%Params%Num_Nu_Massive))
    Rv=1._dl - Rg
        ! Set common power spectrum stuff
    Info%Transfers%Params%InitPower%nn = 1
    Info%Transfers%Params%InitPower%rat(1) = 1._dl  
    Info%Transfers%Params%InitPower%n_run(1) = 0._dl
    Info%Transfers%Params%InitPower%nt_run(1) = 0._dl!tensor
    Info%Transfers%Params%InitPower%k_0_scalar = 2*pi
    Info%Transfers%Params%InitPower%k_0_tensor = 2*pi
    lrat = CMB%maglrat * log(10.d0) ! convert log10 in ln

!Setting up k_D
	Info%Transfers%Params%InitPower%kDissip = (5.5d0 * 10.d0**4)**(1.d0/(CMB%magind + 5.d0)) * (CMB%magamp)**(-2.d0/(CMB%magind+5.d0)) * &
        (2*pi)**((CMB%magind+3)/(CMB%magind+5))*(Info%Transfers%Params%H0/100)**(1/(CMB%magind+5))&
        *((Info%Transfers%Params%omegab*(Info%Transfers%Params%H0/100)**2)/0.022)**(1/(CMB%magind+5))


    !Set them to .true. according to which mode we want to compute, check 
    Info%Transfers%Params%WantScalars = .false.
    Info%Transfers%Params%WantTensors = .false.

    if(CosmoSettings%scalar_passive) then

        Info%Transfers%Params%WantScalars = .true.
        Info%Transfers%Params%Scalar_initial_condition = 1
        
	    delb = 0._dl ! Set up perturbations
        pib = 0._dl

        if(CMB%magind.ge. -2.d0) then
            write(*,*) "Using FITTING FUNCTIONS, scalar passive"
            Info%Transfers%Params%InitPower%an(1) = CMB%magind
            Info%Transfers%Params%InitPower%CorrType = 4
            Info%Transfers%Params%InitPower%ScalarPowerAmp(1)=mag_amplitude(CMB%magind, CMB%magamp)*psconst(3)*&
                                                            (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2

        else
            write(*,*) "Using INTERPOLATION TABLE, scalar passive"
            Info%Transfers%Params%InitPower%an(1) = 1._dl + 2._dl*(3 + CMB%magind)
            Info%Transfers%Params%InitPower%ScalarPowerAmp(1)= mag_psamp(CMB%magind, CMB%magamp, 3) * (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
            Info%Transfers%Params%InitPower%CorrType = 0
        end if
!---helical---yun
        if(CosmoSettings%do_helical) then
          if(CMB%helical_ind < -2.d0) then
             write(*,*) "Using interpolation TABLE for helical passive P_App part"
             Info%Transfers%Params%InitPower%CorrType_hel = 0   
             Info%Transfers%Params%InitPower%an_hel1(1) = 1.d0 + 2.d0*(3.d0+CMB%helical_ind)
             Info%Transfers%Params%InitPower%ScalarPowerAmp_hel(1)=mag_psamp_hel(CMB%helical_ind, CMB%helical_amp, 8) &
                                                    * (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
          else
!Use first fitting formula
             write(*,*) "USING FITTING FORMULAS for helical passive P_App part"
             Info%Transfers%Params%InitPower%an_hel1(1) = CMB%helical_ind
             Info%Transfers%Params%InitPower%ScalarPowerAmp_hel(1)=mag_amplitude_hel(CMB%helical_ind, CMB%helical_amp)*psconst(8)*&
                                                    (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
             Info%Transfers%Params%InitPower%CorrType_hel = 4
           end if                
         end if
!-------------
    end if

    if(CosmoSettings%tensor_passive) then
        Info%Transfers%Params%WantTensors = .true.!YUN
        Info%Transfers%Params%InitPower%tensor_parameterization = 3
        delb = 0._dl ! Set up perturbations
        pib = 0._dl
        tens0 = 1._dl !yun
        !Choose between n>0 and n<0
        if(CMB%magind < -2.d0) then
            write(*,*) "USING INTERPOLATION TABLE, tensor passive"
            Info%Transfers%Params%InitPower%ant(1) = 2._dl*(3 + CMB%magind)
            Info%Transfers%Params%InitPower%TensorPowerAmp(1) = mag_psamp(CMB%magind, CMB%magamp, 5) * (6._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
            Info%Transfers%Params%InitPower%CorrType = 0
        else
            !Use first fitting formula
            write(*,*) "USING FITTING FORMULAS, tensor passive"
            Info%Transfers%Params%InitPower%ant(1) = CMB%magind
            Info%Transfers%Params%InitPower%TensorPowerAmp(1)=mag_amplitude(CMB%magind, CMB%magamp)*psconst(5)*(6._dl*Rg*(lrat+(5._dl/(8._dl*Rv)-1)))**2
            Info%Transfers%Params%InitPower%CorrType = 6
        end if
 !!!helical-----yun !tensor helical passive
        if(CosmoSettings%do_helical) then
           if(CMB%helical_ind < -2.d0) then
              Info%Transfers%Params%InitPower%CorrType_hel = 0
              Info%Transfers%Params%InitPower%TensorPowerAmp_hel(1) = mag_psamp_hel(CMB%helical_ind, CMB%helical_amp,10)&
                                                       *(6._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2 
              Info%Transfers%Params%InitPower%ant_hel1(1) = 2.d0*(CMB%helical_ind +3.d0)
              write(*,*) "Using interpolation TABLE for helical P_tpp"
           else 
              write(*,*) "Using fitting functions for helical P_tpp" 
              Info%Transfers%Params%InitPower%ant_hel1(1) = CMB%helical_ind
              Info%Transfers%Params%InitPower%TensorPowerAmp_hel(1) = mag_amplitude_hel(CMB%helical_ind, CMB%helical_amp)*psconst(10)&
                                                       *(6._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
              Info%Transfers%Params%InitPower%CorrType_hel = 6 
            end if
         end if       
!-------------------------
    end if

    Info%Transfers%Params%DoLensing = .false.

    call CAMB_TransfersToPowers(Info%Transfers)

    if (global_error_flag/=0) then
        error=global_error_flag
        return
    end if

    if(CosmoSettings%use_CMB)then
        call this%SetPowersFromCAMB(CMB,Theory) ! need to modify SetPowersFromCAMB...
        do i=1, min(3,CosmoSettings%num_cls)
            if (CosmoSettings%cl_lmax(i,i)>0) then
                if (any(Theory%MagPassCls(i,i)%Cl(:) < 0 )) then
                    !error = 1
		    !write(*,*) "CMB = ", CMB
                    !call MpiStop('Calculator_CAMB: negative PASSIVE C_l (could edit to silent error here)')
                    !return
                end if
            end if
            do j= i, 1, -1
                if (CosmoSettings%cl_lmax(i,j)>0) then
                    if ( any(isNan(Theory%MagPassCls(i,j)%Cl))) then
                        error=1
                        write(*,*) 'WARNING: NaN CL?', i, j
                        return
                    end if
                end if
            end do
        end do
    end if

end if !scalar-tensor passive
!---------------------------------------------------

!After calculating all the things, re-setting the parameters
if(CosmoSettings%scalar_passive .OR. CosmoSettings%scalar_compensated &
    .OR. CosmoSettings%vector_compensated .OR. CosmoSettings%tensor_passive) then
    delb = 0.d0
    pib = 0.d0

!    call mag_reset_alfven
end if

!----------- SUMMING UP ALL THE CLS -----------------

    if(CosmoSettings%use_CMB)then
! need to sum up all the Cls...
       do i=1, min(3,CosmoSettings%num_cls)
            if (CosmoSettings%cl_lmax(i,i)>0) then

                if(CosmoSettings%scalar_compensated) Theory%Cls(i,i)%CL = Theory%Cls(i,i)%CL + Theory%MagScalCompClsDD(i,i)%CL &
                                                    +Theory%MagScalCompClsDP(i,i)%CL + Theory%MagScalCompClsPP(i,i)%CL

                if(CosmoSettings%vector_compensated)  Theory%Cls(i,i)%CL = Theory%Cls(i,i)%CL + Theory%MagVecCompCls(i,i)%CL

                if(CosmoSettings%scalar_passive .OR. CosmoSettings%tensor_passive) Theory%Cls(i,i)%CL = Theory%Cls(i,i)%CL + Theory%MagPassCls(i,i)%CL

            end if
            do j= i, 1, -1
                if (CosmoSettings%cl_lmax(i,j)>0) then

                    if(CosmoSettings%scalar_compensated) Theory%Cls(i,j)%CL = Theory%Cls(i,j)%CL + Theory%MagScalCompClsDD(i,j)%CL &
                                                        +Theory%MagScalCompClsDP(i,j)%CL + Theory%MagScalCompClsPP(i,j)%CL

                    if(CosmoSettings%vector_compensated)  Theory%Cls(i,j)%CL = Theory%Cls(i,j)%CL + Theory%MagVecCompCls(i,j)%CL

                    if(CosmoSettings%scalar_passive .OR. CosmoSettings%tensor_passive) Theory%Cls(i,j)%CL = Theory%Cls(i,j)%CL + Theory%MagPassCls(i,j)%CL

                end if
            end do
        end do
    end if


! Reset the standard parameters:
Info%Transfers%Params = Pbuffer
!-----------------------------------------------------

!redshifts are in increasing order, so last index is redshift zero
        if (CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8) then
            Theory%sigma_8 = Info%Transfers%MTrans%sigma_8(size(Info%Transfers%MTrans%sigma_8,1),1)
        else
            Theory%sigma_8 = 0
        end if

        if (CosmoSettings%Use_LSS) then
            call this%SetPkFromCAMB(Info%Transfers%MTrans,Theory,error)
            if (error/=0) then
                write(*,*) 'WARNING: NaN PK?'
                return
            end if
        end if
    end select

    end subroutine CAMBCalc_GetNewPowerData

    subroutine CAMBCalc_GetTheoryForImportance(this, CMB, Theory, error)
    use ModelParams, only : ThreadNum
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    class(TCosmoTheoryPredictions) Theory
    class(CAMBTransferCache), allocatable :: Info
    integer error,i,j
    logical :: DoCls, DoPk
    type(CAMBParams)  P

    error = 0
    DoCls = this%ImportanceOptions%redo_cls
    DoPk = this%ImportanceOptions%redo_pk

    if (DoCls .or. DoPk) then
        allocate(Info)
        call CAMB_InitCAMBdata(Info%Transfers)
        call this%CMBToCAMB(CMB, P)
        Threadnum =num_threads

        P%WantCls = DoCls

        if (.not. DoPk .and. .not. (CosmoSettings%CMB_Lensing .and. &
            CosmoSettings%use_nonlinear_lensing)) P%WantTransfer = .false.

        if (this%CAMB_timing) call Timer()
        if (Feedback > 1) write (*,*) 'Calling CAMB'
        call CAMB_GetTransfers(P, Info%Transfers, error)
        if (Feedback > 1) write (*,*) 'CAMB Done'
        if (this%CAMB_timing) call Timer('CAMB_GetTransfers')

        if (error==0) then
            call this%SetCAMBInitPower(Info%Transfers%Params,CMB,1)
            call CAMB_TransfersToPowers(Info%Transfers)
            error=global_error_flag
        end if
    else
        call this%GetNewBackgroundData(CMB,Theory,error)
    end if

    if (DoCls .and. error==0) then
        call this%SetPowersFromCAMB(CMB,Theory)
        if (any(Theory%cls(1,1)%Cl(:) < 0 )) then
            error = 1
            call MpiStop('Calculator_CAMB: negative C_l (could edit to silent error here)')
        end if
        do i=1, min(3,CosmoSettings%num_cls)
            if(error/=0) exit
            do j= i, 1, -1
                if (CosmoSettings%cl_lmax(i,j)>0) then
                    if ( any(isNan(Theory%cls(i,j)%Cl))) then
                        error=1
                        write(*,*) 'WARNING: NaN CL?'
                        exit
                    end if
                end if
            end do
        end do
    end if

    if (DoPK .and. error==0) then
        Theory%sigma_8 = Info%Transfers%MTrans%sigma_8(size(Info%Transfers%MTrans%sigma_8,1),1)
        call this%SetPkFromCAMB(Info%Transfers%MTrans,Theory,error)
        if (error/=0) write(*,*) 'WARNING: NaN PK?'
    end if

    if (error==0) call this%SetDerived(Theory)

    if (DoCls .or. DoPk) call Info%Clear()

    end subroutine CAMBCalc_GetTheoryForImportance

    subroutine CAMBCalc_SetPowersFromCAMB(this,CMB,Theory)
    use constants
    use InitialPower
    use ModelData
    class(CAMB_Calculator) :: this
    class(CMBParams) :: CMB
    class(TCosmoTheoryPredictions), target :: Theory
    real(mcp), parameter :: cons =  (COBE_CMBTemp*1e6)**2
    integer l
    real(mcp) :: highL_norm = 0
    real(mcp) lens_recon_scale, rms
    integer i,j, lmx, lmaxCL
    integer, save, allocatable :: indicesS(:,:), indicesT(:,:)


    if (.not. allocated(indicesS)) then
        allocate(indicesS(3,3))
        allocate(indicesT(3,3))
        indicesS=0
        indicesS(1,1) =  C_Temp
        indicesS(2,2) =  C_E
        indicesT = indicesS
        indicesT(2,1) =  CT_Cross
        indicesT(3,3) =  CT_B
        indicesS(2,1) =  C_Cross
    end if


!-------- NO MAGNETIC magcosmomc----------
if(MagneticMode == no_mag) then

    lens_recon_scale = CMB%InitPower(Aphiphi_index)

    do i=1, min(3,CosmoSettings%num_cls)
        do j= i, 1, -1
            lmaxCL = CosmoSettings%cl_lmax(i,j)
            lmx = min(CosmoSettings%lmax_computed_cl, lmaxCL)
            if (lmx/=0) then
                associate( CL => Theory%Cls(i,j)%CL)
                    if (indicesT(i,j)==0) then
                        CL=0
                    else
                        if (CosmoSettings%CMB_Lensing) then
                            CL(2:lmx) = cons*Cl_lensed(2:lmx,1, indicesT(i,j))
                        else
                            if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_Scalar(2:lmx,1, indicesS(i,j))
                        end if
                        if (CosmoSettings%lmax_computed_cl < lmaxCL) then
                            if (highL_norm ==0) & !normally normalize off TT
                            & highL_norm = CL(lmx)/this%highL_lensedCL_template(lmx,indicesT(i,j))
                            CL(lmx+1:lmaxCL) =  highL_norm*this%highL_lensedCL_template(lmx+1:lmaxCL,indicesT(i,j))
                        end if
                        if (CosmoSettings%compute_tensors) then
                            lmx = min(lmx,CosmoSettings%lmax_tensor)
                            CL(2:lmx) =  CL(2:lmx) + cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                        end if
                    end if
                end associate
            end if
        end do
    end do

    if (CosmoSettings%use_lensing_potential) then
        !CMB lensing potential
        !in camb Cphi is l^4 C_l, we want [l(l+1)]^2Cphi/2pi
        lmx = min(CosmoSettings%lmax_computed_cl, CosmoSettings%cl_lmax(CL_Phi,CL_Phi))
        if (lmx/=0) then
            if (.not. CosmoSettings%CMB_lensing) call MpiStop('Must have lensing on to use lensing potential')
            associate(CL=> Theory%Cls(CL_Phi,CL_Phi)%CL)
                do l=2, lmx
                    CL(L) =  Cl_scalar(L,1, C_Phi)*(real(l+1)**2/l**2)/twopi * lens_recon_scale
                end do
                CL(lmx+1:)=0
            end associate
        end if
        lmx = min(CosmoSettings%lmax_computed_cl, CosmoSettings%cl_lmax(CL_Phi,CL_T))
        if (lmx/=0) then
            !lensing-temp
            do l=2, lmx
                Theory%Cls(CL_phi,CL_T)%CL = Cl_scalar(l,1, C_PhiTemp)/real(l)**3 * sqrt(lens_recon_scale)
            end do
        end if
    end if

    if (CosmoSettings%CMB_Lensing .and. this%CAMBP%max_l>=2000) then
        !Get RMS deflection angle in arcmin
        rms=0
        do L=2, 2000
            rms = rms +  Cl_scalar(L,1, C_Phi)*(real(l+1)**2/l**2)/twopi*(L+0.5_mcp)/(L*(L+1))
        end do
        Theory%Lensing_rms_deflect = sqrt(rms)*180/pi*60
    else
        Theory%Lensing_rms_deflect = 0
    end if

    if (CosmoSettings%compute_tensors) then
        Theory%tensor_ratio_02 = TensorPower(0.002d0,1)/ScalarPower(0.002d0,1)
        Theory%tensor_AT = TensorPower(CosmoSettings%tensor_pivot_k,1)
        Theory%tensor_ratio_BB = TensorPower(0.01d0,1)/ScalarPower(0.01d0,1)
        Theory%tensor_ratio_C10 = Cl_tensor(10, 1, 1)/Cl_scalar(10,1, 1)
    else
        Theory%tensor_ratio_02 = 0
        Theory%tensor_ratio_BB = 0
        Theory%tensor_ratio_C10 = 0
        Theory%tensor_AT = 0
    end if

!----------------------------
! COMPENSATED MODES
! SCALAR ! dd,  Delta-Delta - Delta-Pi, The density perturbed mode (delb = 1, pib = 0)
else if(MagneticMode == scalar_comp_dd) then
	
do i=1, min(3,CosmoSettings%num_cls)
        do j= i, 1, -1
            lmaxCL = CosmoSettings%cl_lmax(i,j)
            lmx = min(CosmoSettings%lmax_computed_cl, lmaxCL)
            if (lmx/=0) then
                associate( CL => Theory%MagScalCompClsDD(i,j)%CL)
                if (indicesT(i,j)==0) then
                    CL=0
                else
                    !if (CosmoSettings%CMB_Lensing) then
                    !    CL(2:lmx) = cons*Cl_lensed(2:lmx,1, indicesT(i,j))
          	    !else
                        if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_Scalar(2:lmx,1, indicesS(i,j))
                        !if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    !end if
                    !if (CosmoSettings%lmax_computed_cl < lmaxCL) then
                    !    if (highL_norm ==0) & !normally normalize off TT
                    !        & highL_norm = CL(lmx)/this%highL_lensedCL_template(lmx,indicesT(i,j))
                    !    CL(lmx+1:lmaxCL) =  highL_norm*this%highL_lensedCL_template(lmx+1:lmaxCL,indicesT(i,j))
                    !end if
                    ! Since I'm doing scalar passive, I don't need tensor now....
                    !if (CosmoSettings%compute_tensors) then
                    !    lmx = min(lmx,CosmoSettings%lmax_tensor)
                    !    CL(2:lmx) =  CL(2:lmx) + cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    !end if
    		end if
                end associate
            end if
        end do
    end do

! ! pp, The stress perturbed mode (delb = 0, pib = 1), write(*,*) "Pi-Pi - Delta-Pi"
else if(MagneticMode == scalar_comp_pp) then

do i=1, min(3,CosmoSettings%num_cls)
        do j= i, 1, -1
            lmaxCL = CosmoSettings%cl_lmax(i,j)
            lmx = min(CosmoSettings%lmax_computed_cl, lmaxCL)
            if (lmx/=0) then
                associate( CL => Theory%MagScalCompClsPP(i,j)%CL)
                if (indicesT(i,j)==0) then
                    CL=0
                else
                    !if (CosmoSettings%CMB_Lensing) then
                    !    CL(2:lmx) = cons*Cl_lensed(2:lmx,1, indicesT(i,j))
                    !else
                        if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_Scalar(2:lmx,1, indicesS(i,j))
                        !if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    !end if
                    !if (CosmoSettings%lmax_computed_cl < lmaxCL) then
                    !    if (highL_norm ==0) & !normally normalize off TT
                    !        & highL_norm = CL(lmx)/this%highL_lensedCL_template(lmx,indicesT(i,j))
                    !    CL(lmx+1:lmaxCL) =  highL_norm*this%highL_lensedCL_template(lmx+1:lmaxCL,indicesT(i,j))
                    !end if
                    ! Since I'm doing scalar passive, I don't need tensor now....
                    !if (CosmoSettings%compute_tensors) then
                    !    lmx = min(lmx,CosmoSettings%lmax_tensor)
                    !    CL(2:lmx) =  CL(2:lmx) + cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    !end if
                end if
          	end associate
            end if
        end do
    end do



! dp, Delta-Pi,The combined perturbation (delb = 1, pib = 1)
else if(MagneticMode == scalar_comp_dp) then

do i=1, min(3,CosmoSettings%num_cls)
        do j= i, 1, -1
            lmaxCL = CosmoSettings%cl_lmax(i,j)
            lmx = min(CosmoSettings%lmax_computed_cl, lmaxCL)
            if (lmx/=0) then
                associate( CL => Theory%MagScalCompClsDP(i,j)%CL)
                if (indicesT(i,j)==0) then
                    CL=0
                else
                    !if (CosmoSettings%CMB_Lensing) then
                    !    CL(2:lmx) = cons*Cl_lensed(2:lmx,1, indicesT(i,j))
                    !else
                        if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_Scalar(2:lmx,1, indicesS(i,j))
                        !if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    !end if
                    !if (CosmoSettings%lmax_computed_cl < lmaxCL) then
                    !    if (highL_norm ==0) & !normally normalize off TT
                    !        & highL_norm = CL(lmx)/this%highL_lensedCL_template(lmx,indicesT(i,j))
                    !    CL(lmx+1:lmaxCL) =  highL_norm*this%highL_lensedCL_template(lmx+1:lmaxCL,indicesT(i,j))
                    !end if
                    ! Since I'm doing scalar passive, I don't need tensor now....
                    !if (CosmoSettings%compute_tensors) then
                    !    lmx = min(lmx,CosmoSettings%lmax_tensor)
                    !    CL(2:lmx) =  CL(2:lmx) + cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    !end if
                end if
          	end associate
            end if
        end do
    end do



!VECTOR
else if(MagneticMode == vector_comp) then

do i=1, min(3,CosmoSettings%num_cls)
        do j= i, 1, -1
            lmaxCL = CosmoSettings%cl_lmax(i,j)
            lmx = min(CosmoSettings%lmax_computed_cl, lmaxCL)
            if (lmx/=0) then
                associate( CL => Theory%MagVecCompCls(i,j)%CL)
                if (indicesT(i,j)==0) then
                    CL=0
                else
                    !if (CosmoSettings%CMB_Lensing) then
                    !    CL(2:lmx) = cons*Cl_lensed(2:lmx,1, indicesT(i,j))
                    !else
                        !if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_vector(2:lmx,1, indicesS(i,j))
                CL(2:lmx) = cons*Cl_vector(2:lmx,1, indicesT(i,j))

			!if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    !end if
                    !if (CosmoSettings%lmax_computed_cl < lmaxCL) then
                    !    if (highL_norm ==0) & !normally normalize off TT
                    !        & highL_norm = CL(lmx)/this%highL_lensedCL_template(lmx,indicesT(i,j))
                    !    CL(lmx+1:lmaxCL) =  highL_norm*this%highL_lensedCL_template(lmx+1:lmaxCL,indicesT(i,j))
                    !end if
                    ! Since I'm doing scalar passive, I don't need tensor now....
                    !if (CosmoSettings%compute_tensors) then
                    !    lmx = min(lmx,CosmoSettings%lmax_tensor)
                    !    CL(2:lmx) =  CL(2:lmx) + cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    !end if

                end if
                end associate
            end if
        end do
    end do


!----------------------------
! PASSIVE tensor or scalar
else if (MagneticMode == scalar_tensor_passive) then
    do i=1, min(3,CosmoSettings%num_cls)
        do j= i, 1, -1
            lmaxCL = CosmoSettings%cl_lmax(i,j)
            lmx = min(CosmoSettings%lmax_computed_cl, lmaxCL)
            if (lmx/=0) then
                associate( CL => Theory%MagPassCls(i,j)%CL)
                if (indicesT(i,j)==0) then
                    CL=0
                else
                    !if (CosmoSettings%CMB_Lensing) then
                    !    CL(2:lmx) = cons*Cl_lensed(2:lmx,1, indicesT(i,j))
                    !else
                    if (CosmoSettings%scalar_passive) then
                        if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_Scalar(2:lmx,1, indicesS(i,j))
                    end if
			!if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    !end if
                    !if (CosmoSettings%lmax_computed_cl < lmaxCL) then
                    !    if (highL_norm ==0) & !normally normalize off TT
                    !        & highL_norm = CL(lmx)/this%highL_lensedCL_template(lmx,indicesT(i,j))
                    !    CL(lmx+1:lmaxCL) =  highL_norm*this%highL_lensedCL_template(lmx+1:lmaxCL,indicesT(i,j))
                    !end if
                    if (CosmoSettings%tensor_passive) then
                        lmx = min(lmx,CosmoSettings%lmax_tensor)
                        CL(2:lmx) =  CL(2:lmx) + cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                    end if

                end if
                end associate
            end if
        end do
    end do

end if !Magnetic Mode

    end subroutine CAMBCalc_SetPowersFromCAMB

    subroutine CAMBCalc_SetPkFromCAMB(this,M,Theory,error)
    use Transfer
    use camb, only : CP
    class(CAMB_Calculator) :: this
    class(TCosmoTheoryPredictions) Theory
    Type(MatterTransferData) M
    integer :: error
    real(mcp), allocatable :: k(:), z(:), PK(:,:)
    integer zix,nz,nk, nR
    real(mcp), allocatable :: NL_Ratios(:,:)
    real(mcp) :: dR, R, minR
    integer i

    !Free theory arrays as they may resize between samples
    call Theory%FreePK()

    nz=CP%Transfer%PK_num_redshifts

    if (.not. allocated(Theory%growth_z)) allocate(Theory%growth_z, Theory%sigma8_z)
    call Theory%growth_z%InitForSize(nz)
    call Theory%sigma8_z%InitForSize(nz)
    allocate(z(nz))

    do zix=1,nz
        z(zix) = CP%Transfer%PK_redshifts(nz-zix+1)
        Theory%sigma8_z%F(zix) = M%sigma_8(nz-zix+1,1)
        Theory%growth_z%F(zix) = M%sigma2_vdelta_8(nz-zix+1,1)/M%sigma_8(nz-zix+1,1)
    end do
    Theory%sigma8_z%X=z
    Theory%growth_z%X=z

    if (CosmoSettings%use_matterpower) then
        nk=M%num_q_trans
        nz=CP%Transfer%PK_num_redshifts
        allocate(PK(nk,nz))
        allocate(k(nk))

        k = log(M%TransferData(Transfer_kh,:,1))

        call Transfer_GetUnsplinedPower(M, PK,transfer_power_var,transfer_power_var)
        PK = Log(PK)
        if (any(isNan(PK))) then
            error = 1
            return
        end if
        allocate(Theory%MPK)
        call Theory%MPK%Init(k,z,PK)
    end if


    if (CosmoSettings%use_Weylpower) then
        call Transfer_GetUnsplinedPower(M, PK,transfer_Weyl,transfer_Weyl,hubble_units=.false.)
        PK = Log(PK)
        if (any(isNan(PK))) then
            error = 1
            return
        end if
        allocate(Theory%MPK_WEYL)
        call Theory%MPK_WEYL%Init(k,z,PK)
    end if

    if (CosmoSettings%use_SigmaR) then
        !Note R is in k*h units
        dR = log(1.2)/AccuracyLevel
        minR = 1/CP%Transfer%kmax
        nR = nint(log(150/minR)/dR) +1
        if (.not. allocated(Theory%Sigma_R)) allocate(Theory%Sigma_R)
        call Theory%Sigma_R%InitForSize(nR)
        do i=1, nR
            Theory%Sigma_R%X(i) = exp((i-1)*dR)*minR
        end do
        call Transfer_GetSigmaRArray(M, Theory%Sigma_R%X, Theory%Sigma_R%F, &
            var1 = transfer_nonu,var2=transfer_nonu)
    end if

    if(CosmoSettings%use_nonlinear)then
        call this%GetNLandRatios(M,Theory,NL_Ratios,error)
        if(error/=0) return
    end if

    end subroutine CAMBCalc_SetPkFromCAMB


    subroutine CAMBCalc_GetNLandRatios(this,M,Theory,Ratios,error)
    use Transfer
    class(CAMB_Calculator) :: this
    class(TCosmoTheoryPredictions) Theory
    Type(MatterTransferData) M
    real(mcp), allocatable, intent(out) :: Ratios(:,:)
    Type(MatterPowerData) :: CPK
    real(mcp), allocatable :: PK(:,:)
    integer error,zix,nz

    CPK%num_k = Theory%MPK%nx
    CPK%num_z = Theory%MPK%ny

    !Allocate Theory arrays
    allocate(Theory%NL_MPK)
    allocate(Ratios(CPK%num_k,CPK%num_z))

    !fill PK with Linear MPK
    allocate(PK(CPK%num_k,CPK%num_z))
    PK=Theory%MPK%z

    allocate(CPK%matpower(CPK%num_k,CPK%num_z))
    allocate(CPK%ddmat(CPK%num_k,CPK%num_z))
    allocate(CPK%nonlin_ratio(CPK%num_k,CPK%num_z))
    allocate(CPK%log_kh(CPK%num_k))
    allocate(CPK%redshifts(CPK%num_z))
    CPK%log_kh = Theory%MPK%x
    CPK%redshifts = Theory%MPK%y
    CPK%matpower = PK

    !need splines to get nonlinear ratios
    call MatterPowerdata_getsplines(CPK)
    call NonLinear_GetRatios(CPK)
    Ratios = CPK%nonlin_ratio
    call MatterPowerdata_Free(CPK)

    PK = PK+2*log(Ratios)
    if (any(isNan(PK))) then
        error = 1
        return
    end if
    call Theory%NL_MPK%Init(Theory%MPK%x,Theory%MPK%y,PK)

    if (allocated(Theory%MPK_WEYL)) then
        !Assume Weyl scales the same way under non-linear correction
        allocate(Theory%NL_MPK_WEYL)
        PK = Theory%MPK_WEYL%z + 2*log(Ratios)
        call Theory%NL_MPK_WEYL%Init(Theory%MPK%x,Theory%MPK%y,PK)
    end if

    end subroutine CAMBCalc_GetNLandRatios

    subroutine CAMBCalc_InitCAMB(this,CMB,error, DoReion)
    class(CAMB_Calculator) :: this
    class(CMBParams), intent(in) :: CMB
    logical, optional, intent(in) :: DoReion
    logical WantReion
    type(CAMBParams)  P
    integer error

    if (present(DoReion)) then
        WantReion = DoReion
    else
        WantReion = .true.
    end if

    call this%CMBToCAMB(CMB, P)
    call CAMBParams_Set(P,error,WantReion)

    end subroutine CAMBCalc_InitCAMB

    function CAMBCalc_GetOpticalDepth(this,CMB) result(GetOpticalDepth)
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    real(mcp) GetOpticalDepth
    type(CAMBParams)  P
    integer error

    call this%CMBToCAMB(CMB, P)
    call CAMBParams_Set(P,error)

    if (error/= 0) then
        GetOpticalDepth = -1
    else
        GetOpticalDepth = Reionization_GetOptDepth(P%Reion, P%ReionHist)
    end if
    end function CAMBCalc_GetOpticalDepth

    function CAMBCalc_GetZreFromTau(this,CMB, tau) result(GetZreFromTau)
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    real(mcp), intent(in) :: tau
    real(mcp) GetZreFromTau
    type(CAMBParams)  P

    call this%CMBToCAMB(CMB, P)
    GetZreFromTau = CAMB_GetZreFromTau(P,dble(tau))

    end function CAMBCalc_GetZreFromTau

    function CAMBCalc_CMBToTheta(this,CMB) result(CMBToTheta)
    use ModelParams
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    real(mcp) CMBToTheta
    integer error

    call this%InitCAMB(CMB,error,.false.)
    CMBToTheta = CosmomcTheta()

    end function CAMBCalc_CMBToTheta


    real(mcp) function CAMBCalc_BAO_D_v(this, z)
    use CAMB, only : BAO_D_v
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_BAO_D_v = BAO_D_v(z)

    end function CAMBCalc_BAO_D_v


    real(mcp) function CAMBCalc_AngularDiameterDistance(this, z)
    use CAMB, only : AngularDiameterDistance  !!angular diam distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_AngularDiameterDistance = AngularDiameterDistance(z)

    end function CAMBCalc_AngularDiameterDistance

    real(mcp) function CAMBCalc_ComovingRadialDistance(this, z)
    use CAMB, only : ComovingRadialDistance  !!comoving radial distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_ComovingRadialDistance = ComovingRadialDistance(z)

    end function CAMBCalc_ComovingRadialDistance

    real(mcp) function CAMBCalc_AngularDiameterDistance2(this, z1, z2)
    use CAMB, only : AngularDiameterDistance2  !!angular diam distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z1, z2

    CAMBCalc_AngularDiameterDistance2 = AngularDiameterDistance2(z1, z2)

    end function CAMBCalc_AngularDiameterDistance2

    real(mcp) function CAMBCalc_LuminosityDistance(this, z)
    use CAMB, only : LuminosityDistance  !! distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_LuminosityDistance = LuminosityDistance(z)

    end function CAMBCalc_LuminosityDistance

    real(mcp) function CAMBCalc_Hofz(this, z)
    use CAMB, only : Hofz  !!angular diam distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_Hofz = Hofz(z)

    end function CAMBCalc_Hofz


    subroutine CAMBCalc_InitCAMBParams(this,P)
    use lensing
    use ModelParams
    class(CAMB_Calculator) :: this
    type(CAMBParams)  P
    integer zix
    !JD Changed P%Transfer%redshifts and P%Transfer%num_redshifts to
    !P%Transfer%PK_redshifts and P%Transfer%PK_num_redshifts respectively
    !for nonlinear lensing of CMB + LSS compatibility
    Threadnum =num_threads
    w_lam = -1
    !wa_ppf = 0._dl
    call CAMB_SetDefParams(P)

    HighAccuracyDefault = .true.
    P%OutputNormalization = outNone

    !JD Modified to save computation time when only using MPK
    if(CosmoSettings%use_CMB) then
        P%WantScalars = .true.
        P%WantTensors = CosmoSettings%compute_tensors
        P%Max_l=CosmoSettings%lmax_computed_cl
        P%Max_eta_k=CosmoSettings%lmax_computed_cl*2
        P%Max_l_tensor=CosmoSettings%lmax_tensor
        P%Max_eta_k_tensor=CosmoSettings%lmax_tensor*5./2
    else
        P%WantCls = .false.
    end if

    P%WantTransfer = CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8
    P%Transfer%k_per_logint=0

    if (CosmoSettings%use_nonlinear) then
        P%NonLinear = NonLinear_pk
        P%Transfer%kmax = max(1.2_mcp,CosmoSettings%power_kmax)
    else
        P%Transfer%kmax = max(0.8_mcp,CosmoSettings%power_kmax)
    end if

    if (AccuracyLevel > 1 .or. HighAccuracyDefault) then
        if (CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8) then
            P%Transfer%high_precision=.true.
            P%Transfer%kmax=P%Transfer%kmax + 0.2
        end if
        AccuracyBoost = AccuracyLevel
        lAccuracyBoost = AccuracyLevel
        lSampleBoost = AccuracyLevel
        P%AccurateReionization = .true.
    end if

    P%AccurateBB = this%accurate_BB

    if (max_transfer_redshifts < CosmoSettings%num_power_redshifts) then
        stop 'Need to manually set max_transfer_redshifts larger in CAMB''s modules.f90'
    end if

    if (CosmoSettings%num_power_redshifts>1) then
        P%Transfer%PK_num_redshifts = CosmoSettings%num_power_redshifts
        do zix=1, CosmoSettings%num_power_redshifts
            !CAMB's ordering is from highest to lowest
            P%Transfer%PK_redshifts(zix) = CosmoSettings%power_redshifts(CosmoSettings%num_power_redshifts-zix+1)
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts(1) = 0
    end if

    P%Num_Nu_Massive = 3
    P%Num_Nu_Massless = 0.046
    P%InitPower%nn = 1
    P%AccuratePolarization = CosmoSettings%num_cls/=1
    P%Reion%use_optical_depth = .false.
    P%OnlyTransfers = .true.

    if (CosmoSettings%CMB_Lensing) then
        P%DoLensing = .true.
        P%Max_l = CosmoSettings%lmax_computed_cl +100 + 50 !+50 in case accuracyBoost>1 and so odd l spacing
        P%Max_eta_k = P%Max_l*2
    end if

    if (HighAccuracyDefault) then
        P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl*AccuracyLevel,P%Max_eta_k)
        if (CosmoSettings%CMB_Lensing .and. (CosmoSettings%use_lensing_potential .or. CosmoSettings%use_nonlinear_lensing)) &
            P%Max_eta_k = max(P%Max_eta_k, 14000*AccuracyLevel)
        !k_etamax=18000 give c_phi_phi accurate to sub-percent at L=1000, <4% at L=2000
        !k_etamax=10000 is just < 1% at L<=500
    end if
    if (this%k_eta_max_scalar>0) then
        P%Max_eta_k = this%k_eta_max_scalar
    end if
    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
    if (CosmoSettings%CMB_Lensing .and. CosmoSettings%use_nonlinear_lensing) then
        P%WantTransfer = .true.
        P%NonLinear = NonLinear_lens
        call Transfer_SetForNonlinearLensing(P%Transfer)
        if(CosmoSettings%use_nonlinear) P%NonLinear = NonLinear_both
    end if
    call Transfer_SortAndIndexRedshifts(P%Transfer)
    !End JD modifications
    lensing_includes_tensors = .false.

    P%Scalar_initial_condition = initial_vector
    P%InitialConditionVector = 0
    P%InitialConditionVector(initial_adiabatic) = -1

    BackgroundOutputs%z_outputs => CosmoSettings%z_outputs

    end subroutine CAMBCalc_InitCAMBParams

    !Mapping between array of power spectrum parameters and CAMB
    subroutine CAMBCalc_SetCAMBInitPower(this,P,CMB,ix)
    class(CAMB_Calculator) :: this
    type(CAMBParams)  P
    class(CMBParams) CMB
    integer, intent(in) :: ix

    if (Power_Name == 'power_tilt') then
        P%InitPower%k_0_scalar = CosmoSettings%pivot_k
        P%InitPower%k_0_tensor = CosmoSettings%tensor_pivot_k
        if (P%InitPower%k_0_tensor/=P%InitPower%k_0_scalar) P%InitPower%tensor_parameterization = tensor_param_rpivot
        P%InitPower%ScalarPowerAmp(ix) = cl_norm*CMB%InitPower(As_index)
        P%InitPower%rat(ix) = CMB%InitPower(amp_ratio_index)

        P%InitPower%an(ix) = CMB%InitPower(ns_index)

        if (P%InitPower%rat(ix)>0 .and. .not. CosmoSettings%compute_tensors) &
            call MpiStop('computing r>0 but compute_tensors=F')
        P%InitPower%n_run(ix) = CMB%InitPower(nrun_index)
        P%InitPower%n_runrun(ix) = CMB%InitPower(nrunrun_index)

        if (CosmoSettings%inflation_consistency) then
            if (CMB%InitPower(nt_index)/=0 .or. CMB%InitPower(ntrun_index)/=0) &
                & call MpiStop('Error: inflation_consistency but n_t not set to zero')
            ! first-order consistency relation
            !P%InitPower%ant(ix) = - CMB%InitPower(amp_ratio_index)/8
            !next order consistency relation
            P%InitPower%ant(ix) = - CMB%InitPower(amp_ratio_index)/8*(2-CMB%InitPower(ns_index) - CMB%InitPower(amp_ratio_index)/8)
            P%InitPower%nt_run(ix) = CMB%InitPower(amp_ratio_index)/8* &
                & (CMB%InitPower(amp_ratio_index)/8 + CMB%InitPower(ns_index) - 1)
            !note input n_T, nt run is ignored, so should be fixed
        else
            P%InitPower%ant(ix) = CMB%InitPower(nt_index)
            P%InitPower%nt_run(ix) = CMB%InitPower(ntrun_index)
        end if
    else
        stop 'CAMB_Calculator:Wrong initial power spectrum'
    end if

    end subroutine CAMBCalc_SetCAMBInitPower

    subroutine CAMBCalc_ReadParams(this,Ini)
    use NonLinear
    class(CAMB_Calculator) :: this
    class(TSettingIni) :: Ini

    call this%TCosmologyCalculator%ReadParams(Ini)
    this%calcName ='CAMB'

    this%CAMB_timing = Ini%Read_Logical('CAMB_timing',.false.)

    this%highL_theory_cl_template_file = Ini%ReadFilename('highL_theory_cl_template',DataDir,.true.)

    if (Ini%HasKey('highL_unlensed_cl_template')) then
        highL_unlensed_cl_template=  Ini%ReadFilename('highL_unlensed_cl_template')
    else
        highL_unlensed_cl_template = concat(LocalDir,'camb/',highL_unlensed_cl_template)
    end if

    this%k_eta_max_scalar = Ini%Read_Double('k_eta_max_scalar',-1._mcp)
    this%accurate_BB = Ini%Read_Logical('accurate_BB',.false.)

    halofit_version = Ini%Read_Int('halofit_version',halofit_default)

    end subroutine CAMBCalc_ReadParams


    subroutine CAMBCalc_InitForLikelihoods(this)
    !Called later after likelihoods etc loaded
    class(CAMB_Calculator) :: this

    if (CosmoSettings%use_CMB .and. CosmoSettings%lmax_computed_cl /= CosmoSettings%lmax) then
        if (CosmoSettings%compute_tensors .and. CosmoSettings%lmax_tensor > CosmoSettings%lmax_computed_cl) &
            & call MpiStop('lmax_tensor > lmax_computed_cl')
        call this%LoadFiducialHighLTemplate()
    end if

    call this%InitCAMBParams(this%CAMBP)

    if (Feedback > 0 .and. MPIRank==0) then
        if(CosmoSettings%use_CMB) write(*,*) 'max_eta_k         = ', real(this%CAMBP%Max_eta_k)
        write(*,*) 'transfer kmax     = ', real(this%CAMBP%Transfer%kmax)
    end if

    this%CAMBP%WantTensors = CosmoSettings%compute_tensors

    end subroutine CAMBCalc_InitForLikelihoods


    subroutine CAMBCalc_VersionTraceOutput(this, ReadValues)
    use GaugeInterface, only : Eqns_name
    class(CAMB_Calculator) :: this
    class(TNameValueList) :: ReadValues

    !Store for the record any useful info about version etc.
    call ReadValues%Add('Compiled_CAMB_version', version)
    call ReadValues%Add('Compiled_Recombination', Recombination_Name)
    call ReadValues%Add('Compiled_Equations', Eqns_name)
    call ReadValues%Add('Compiled_Reionization', Reionization_Name)
    call ReadValues%Add('Compiled_InitialPower', Power_Name)

    end subroutine CAMBCalc_VersionTraceOutput



    subroutine LoadFiducialHighLTemplate(this)
    class(CAMB_Calculator) :: this
    !This should be a lensed scalar CMB power spectrum, e.g. for including at very high L where foregrounds etc. dominate anyway
    integer L,  status
    real(mcp) array(4)
    Type(TTextFile) :: F

    allocate(this%highL_lensedCL_template(2:CosmoSettings%lmax, 4))
    call F%Open(this%highL_theory_cl_template_file)
    do
        read(F%unit,*, iostat=status) L , array
        if (status/=0 .or. L>CosmoSettings%lmax) exit
        if (L>=2) this%highL_lensedCL_template(L,:) = array
    end do
    call F%Close()

    if (this%highL_lensedCL_template(2,1) < 100) &
        call MpiStop('highL_theory_cl_template must be in muK^2')

    if (L<CosmoSettings%lmax) call MpiStop('highL_theory_cl_template does not go to lmax')

    end subroutine LoadFiducialHighLTemplate


    !!! CAMBTransferCache

    subroutine CAMBTransferCache_Clear(Info)
    class(CAMBTransferCache) Info

    call CAMB_FreeCAMBdata(Info%Transfers)
!---------
! MagCosmoMC: other CAMBdata
    call CAMB_FreeCAMBdata(Info%VectorTransfers)
    call CAMB_FreeCAMBdata(Info%ScalarCompTransfersDD)
    call CAMB_FreeCAMBdata(Info%ScalarCompTransfersPP)
    call CAMB_FreeCAMBdata(Info%ScalarCompTransfersDP)
!---------
    call Info%TTheoryIntermediateCache%Clear()

    end subroutine CAMBTransferCache_Clear

    end module Calculator_CAMB
