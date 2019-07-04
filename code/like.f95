module like
  use params
  use util
  use NeutrinoFuncs
  use WIMPFuncs
  use LabFuncs
  implicit none

contains

!====================================like.f95==================================!
! Modeule for doing all of the Likelihood analysis
!
! Contents:
!
! 1. Compute limits for input experiment
! CYGNUSLimit: Generates discovery limit (DL) for CYGNUS (He and F)
! GetLimits: Generates discovery limit for arbitrary experiment
! DiscoveryLimit: Mass scan for Median discovery limit at 3sigma
!
! 2. Likelihoods
! lnPF = sum over log(Poisson pdf)
! lnGF = sum over log(Gaussian pdf)
! llhood1 = signal+background likelihood
! llhood0 = background likelihood
!==============================================================================!



!==================================Generate limits=============================!
!---------------------------------------CYGNUS---------------------------------!
subroutine CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	! Limits from m_chi = m_min to m_max, with nm values
	! Cross section scan from sigma_min to sigma_max with ns values
	! DL = Discovery Limit
	double precision :: E_th_F,E_th_He
	double precision :: m_min,m_max,sigma_min,sigma_max,m_vals(nm),DL(nm),DLF(nm),DLHe(nm)
	integer :: i,nm,nf,ns
	character(len=*) :: filename

  	write(*,*) '----------------------------------------------------'
 	write(*,*) 'Start Limits:'
	write(*,*) 'Filename:',trim(filename)

	! CYGNUS is at Boubly
	lat = Boulby(1)
	long = Boulby(2)

	! Calculate exposure for specified TPC Volume x Time
	! 1000 m^3 of SF6 at 20 torr or He at 740 torr is 0.16 tons
  	if (searchmode) then
    !	 1000 m^3 of SF6 at 200 torr or He at 740 torr is 1.6 tons
	  	Exposure = 10.0*VolTime*(0.16/1000.0d0) ! Convert m^3-years into ton-years
   		nside = 0 ! switch off directionality in search mode
  	else
   	 ! 1000 m^3 of SF6 at 20 torr or He at 740 torr is 0.16 tons
   	 	Exposure = VolTime*(0.16/1000.0d0) ! Convert m^3-years into ton-years
  	end if

 	! Calculate Helium limits
  	if (searchmode) then
    	DLHe = 0.0d0 ! No helium used in search mode
  	else
    	nucleus = Helium
    	E_th = E_th_He
		Exposure = Exposure*(755.0/740.0)
  		call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)
		Exposure = Exposure/(755.0/740.0)
  	end if


	! Calculate Fluorine limits
  	nucleus = Fluorine
  	E_th = E_th_F
	Exposure = Exposure*(5.0/20.0)
  	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)


	! Save Data
	open(unit=123,file=trim(filename))
	write(123,*) m_vals
	write(123,*) DLF
	write(123,*) DLHe
	close(123)
	write(*,*) 'writing to: ',trim(filename)
  write(*,*) '----------------------------------------------------'
end subroutine

!---------------------------------Abitrary limit-------------------------------!
subroutine NwimpEvents(Nwimp,m_min,m_max,nm,m_vals,DL)
	! Limits from m_chi = m_min to m_max, with nm values
	! Cross section scan from sigma_min to sigma_max with ns values
	! DL = Discovery Limit
	double precision :: Nwimp,m_min,m_max,sigma_min,sigma_max,m_vals(nm),DL(nm)
	integer :: i,nm,nf,ns

    call GetNuFluxes ! Load Neutrinos
    call PreAllocate ! Allocate data size (readout dependent)
    call BackgroundRecoilDistribution ! Load Background	model
    call SHM ! Load halo model
    DL = 0.0
    m_vals = logspace(m_min,m_max,nm)
    do i = 1,nm
      m_chi = m_vals(i)
      call WIMPRecoilDistribution	! Call WIMP recoil distribution for each new mass
      DL(i) = (Nwimp*1.0)/sum(RD_wimp)
    end do
    call UnAllocate ! Reset
end subroutine

subroutine GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL)
	! Limits from m_chi = m_min to m_max, with nm values
	! Cross section scan from sigma_min to sigma_max with ns values
	! DL = Discovery Limit
	double precision :: m_min,m_max,sigma_min,sigma_max,m_vals(nm),DL(nm)
	integer :: i,nm,nf,ns
    write(*,*) 'Nucleus = ',nucleus,'Exposure = ',Exposure,'ton years'
    write(*,*) '[E_th,E_max] = ',E_th,E_max,'keV'
    write(*,*) '----------------------------------------------------'
    call GetNuFluxes ! Load Neutrinos
    call PreAllocate ! Allocate data size (readout dependent)
    call BackgroundRecoilDistribution ! Load Background	model
    call SHM ! Load halo model
    call DiscoveryLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL)
    call UnAllocate ! Reset
end subroutine

subroutine GetNuFloor_Simple(N_nu_events,m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL)
	double precision :: N_nu_events,N_tot_bg,EB8,m_min,m_max,sigma_min,sigma_max,m_vals(nm)
  double precision :: DL_low(nm),DL_high(nm),DL(nm)
	integer :: i,nm,nf,ns,si
  character(len=100) :: fn_end


  ! Set all directional stuff to something, not relevant for this file
	lat = Boulby(1)
	long = Boulby(2)
	energy_on = .true. ! energy info is currently turned on for best limits
	searchmode = .false.
	angres_on = .false.
	efficiency_on = .false.
	headtail_on = .false.
	energyres_on = .false.
	call LoadReadout(9,	fn_end)

  ! LOW NU FLOOR
  E_th = 1.0d-5
  Exposure = 1.0
  call GetNuFluxes ! Load Neutrinos
  call PreAllocate ! Allocate data size (readout dependent)
  call BackgroundRecoilDistribution ! Load Background	model
  N_tot_bg = 0.0d0
  do si = 1,n_bg
    N_tot_bg = N_tot_bg + sum(R_bg(si)*RD_bg(:,si))
  end do
  Exposure = N_nu_events/N_tot_bg ! Choose right exposure for N_nu events
  RD_bg = RD_bg*Exposure
  call SHM ! Load halo model
  call DiscoveryLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL_low)
  call UnAllocate ! Reset


  ! HIGH NU FLOOR
  E_th = MaxNuRecoilEnergy(6)+0.1d0 ! For there to be no solar nu
  Exposure = 1.0
  call GetNuFluxes ! Load Neutrinos
  call PreAllocate ! Allocate data size (readout dependent)
  call BackgroundRecoilDistribution ! Load Background	model
  N_tot_bg = 0.0d0
  do si = 1,n_bg
    N_tot_bg = N_tot_bg + sum(R_bg(si)*RD_bg(:,si))
  end do
  Exposure = N_nu_events/N_tot_bg
  RD_bg = RD_bg*Exposure
  call SHM ! Load halo model
  call DiscoveryLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL_high)
  call UnAllocate ! Reset

  where (DL_low.eq.0.0)
    Dl_low = 100.0
  end where

  where (DL_high.eq.0.0)
    DL_high = 100.0
  end where

  do i = 1,nm
      DL(i) = min(DL_high(i),DL_low(i))
  end do

end subroutine



!---------------------------------Abitrary limit-------------------------------!
subroutine GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL,Nsig,Nbg)
	! Limits from m_chi = m_min to m_max, with nm values
	! Cross section scan from sigma_min to sigma_max with ns values
	! DL = Discovery Limit
	double precision :: ex_min,ex_max,m,sigma_min,sigma_max,m_vals(1)
  double precision :: DL(n_ex),ex_vals(n_ex),Nsig(n_ex),Nbg(n_ex)
	integer :: i,nf,ns,n_ex,si
  write(*,*) 'Nucleus = ',nucleus,'Exposure = ',Exposure,'ton years'
  write(*,*) '[E_th,E_max] = ',E_th,E_max,'keV'
  write(*,*) '----------------------------------------------------'
  Exposure = 1.0
  call GetNuFluxes ! Load Neutrinos
  call PreAllocate ! Allocate data size (readout dependent)
  call BackgroundRecoilDistribution ! Load Background	model
  call SHM ! Load halo model
  ex_vals = logspace(ex_min,ex_max,n_ex)
  do i=1,n_ex
    Exposure = ex_vals(i)
    RD_bg = RD_bg*Exposure
    write(*,*) 'Exposure = ',Exposure
    call DiscoveryLimit(m,m,1,sigma_min,sigma_max,ns,	m_vals,DL(i))

    ! numbers of events
    Nsig(i) = sum(RD_wimp*DL(i))
    Nbg(i) = 0.0d0
    do si = 1,n_bg
      Nbg(i) = Nbg(i) + sum(R_bg(si)*RD_bg(:,si))
    end do

    RD_bg = RD_bg/Exposure
  end do
  call UnAllocate ! Reset
end subroutine



!---------------------------------Abitrary limit-------------------------------!
subroutine GetLimits_ExposureGradient(m_min,m_max,nm,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,DL_ex)
  integer :: i,j,k,k1,nm,nf,ns,n_ex,si,MAXFUNEVALS,IPRINT,NLOOP,IQUAD,ifault0,ii
	double precision :: m_min,m_max,ex_min,ex_max,m,sigma_min,sigma_max,m_vals(nm),D_prev,s_prev,grad_ex,grad_ex_prev
  double precision :: DL_ex(nm),DL(nm,n_ex),ex_vals(n_ex),Nsig(n_ex),Nbg(n_ex),sigma_p_vals(ns),N_tot_bg
  double precision,dimension(:),allocatable :: x_in0,x_in1,step0,N_exp,N_exp_bg
  double precision :: D01,L1,L0,SIMP,STOPCR0,var(2)
  write(*,*) 'Nucleus = ',nucleus,'Exposure = ',ex_min,ex_max
  write(*,*) '[E_th,E_max] = ',E_th,E_max,'keV'
  write(*,*) '----------------------------------------------------'
  MAXFUNEVALS = 10000 ! Maximum function evaluations
  IPRINT = -1 ! Print results from minimisaition
  NLOOP = 1 ! Number of iterations before looping
  IQUAD = 0 ! Can't remeber what this does
  SIMP = 0.1 ! Nor this
	STOPCR0 = 1.0d-10 ! Accuracy of max likelihood

  Exposure = 1.0
  call GetNuFluxes ! Load Neutrinos
  call PreAllocate ! Allocate data size (readout dependent)
  allocate(x_in0(n_bg))
  allocate(x_in1(n_bg+1))
  allocate(step0(n_bg))
  allocate(N_exp_bg(nTot_bins))
  allocate(N_exp(nTot_bins))
  call BackgroundRecoilDistribution ! Load Background	model
  call SHM ! Load halo model
  ex_vals = logspace(ex_min,ex_max,n_ex)
  m_vals = logspace(m_min,m_max,nm)
  sigma_p_vals = logspace(sigma_min,sigma_max,ns)
  DL = 0.0
  DL_ex = 0.0

  do i = 1,nm
    k1 = 1
    m_chi = m_vals(i)
    Exposure = 1.0
    call WIMPRecoilDistribution	! Call WIMP recoil distribution for each new mass
    if (sum(RD_wimp).gt.0.0) then
      grad_ex_prev = -1000.0
      do j=1,n_ex

        Exposure = ex_vals(j)
        RD_wimp = RD_wimp*Exposure
        RD_bg = RD_bg*Exposure
        N_exp_bg = 0.0d0
        do si = 1,n_bg
          N_exp_bg = N_exp_bg + R_bg(si)*RD_bg(:,si)
        end do
        N_tot_bg = sum(N_exp_bg)

        do k = 1,ns
          sigma_p = sigma_p_vals(k)
          !write(*,*) i,j,k,sigma_p,N_tot_bg,sum(RD_wimp*sigma_p)
          if (sum(RD_wimp*sigma_p).gt.1.0d0) then	! Generally need >0.5 events to see DM
            N_exp = N_exp_bg + RD_wimp*sigma_p
            N_obs = N_exp  ! Observed=Expected for Asimov data
            X_in1= (/log10(sigma_p),R_bg/)
            call llhood1(X_in1,L1) ! Asimov data maximises likelihood at correct value

            X_in0 = R_bg
            step0 = R_bg_err*R_bg
            call llhood0(X_in0,L0)
            call MINIM(X_in0,step0,n_bg,L0,MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,IQUAD,SIMP,VAR,llhood0,IFAULT0)

            D01 = -2.0*(L1-L0)
            if (D01.ge.9.0d0) then ! Median 3sigma detection -> D = 9
              ! Do interpolation to find discovery limit cross section
              DL(i,j) = 10.0d0**(interp1D((/D_prev,D01/),(/log10(s_prev),log10(sigma_p)/),2,9.0d0))
              write(*,*) 'm = ',m_chi,'|| Exposure = ',Exposure,'|| DL = ',DL(i,j)
              k1 = k-1
              exit
            end if
            s_prev = sigma_p ! Reset for interpolation
            D_prev = D01
          end if
        end do

        RD_bg = RD_bg/Exposure
        RD_wimp = RD_wimp/Exposure

        if (j.gt.1) then
          grad_ex = DL(i,j)-DL(i,j-1)
          write(*,*) DL(i,j),DL(i,j-1),grad_ex,grad_ex_prev,j
          if (grad_ex.lt.grad_ex_prev) then
            DL_ex(i) = DL(i,j)
            exit
          end if
          grad_ex_prev = grad_ex
        end if

      end do
    end if
  end do
  call UnAllocate ! Reset
end subroutine

subroutine GetLimits_MassExposure2(m_min,m_max,nm,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,DL,Nsig,Nbg)
  integer :: i,j,k,k1,nm,nf,ns,n_ex,si,MAXFUNEVALS,IPRINT,NLOOP,IQUAD,ifault0,ii
	double precision :: m_min,m_max,ex_min,ex_max,m,sigma_min,sigma_max,m_vals(nm),D_prev,ex_prev
  double precision :: DL(nm,ns),ex_vals(n_ex),Nsig(nm,ns),Nbg(nm,ns),sigma_p_vals(ns),N_tot_bg
  double precision,dimension(:),allocatable :: x_in0,x_in1,step0,N_exp,N_exp_bg
  double precision :: D01,L1,L0,SIMP,STOPCR0,var(2)
  write(*,*) 'Nucleus = ',nucleus,'Exposure = ',ex_min,ex_max
  write(*,*) '[E_th,E_max] = ',E_th,E_max,'keV'
  write(*,*) '----------------------------------------------------'
  MAXFUNEVALS = 10000 ! Maximum function evaluations
  IPRINT = -1 ! Print results from minimisaition
  NLOOP = 2 ! Number of iterations before looping
  IQUAD = 0 ! Can't remeber what this does
  SIMP = 0.1 ! Nor this
	STOPCR0 = 1.0d-6 ! Accuracy of max likelihood

  Exposure = 1.0
  call GetNuFluxes ! Load Neutrinos
  call PreAllocate ! Allocate data size (readout dependent)
  allocate(x_in0(n_bg))
  allocate(x_in1(n_bg+1))
  allocate(step0(n_bg))
  allocate(N_exp_bg(nTot_bins))
  allocate(N_exp(nTot_bins))
  call BackgroundRecoilDistribution ! Load Background	model
  call SHM ! Load halo model
  ex_vals = logspace(ex_min,ex_max,n_ex)
  m_vals = logspace(m_min,m_max,nm)
  sigma_p_vals = logspace(sigma_min,sigma_max,ns)
  DL = 0.0d0
  Nsig = 0.0d0
  Nbg = 0.0d0
  do i = 1,nm
    k1 = 1
    m_chi = m_vals(i)
    Exposure = 1.0
    call WIMPRecoilDistribution	! Call WIMP recoil distribution for each new mass

    if (sum(RD_wimp).gt.0.0) then
      do j=1,ns

        sigma_p = sigma_p_vals(ns+1-j)

        do k = k1,n_ex
          Exposure = ex_vals(k)
          RD_wimp = RD_wimp*Exposure
          RD_bg = RD_bg*Exposure
          N_exp_bg = 0.0d0
          do si = 1,n_bg
            N_exp_bg = N_exp_bg + R_bg(si)*RD_bg(:,si)
          end do
          N_tot_bg = sum(N_exp_bg)

          !write(*,*) i,j,k,sigma_p,N_tot_bg,sum(RD_wimp*sigma_p)
          if (sum(RD_wimp*sigma_p).gt.0.5d0) then	! Generally need >0.5 events to see DM
            N_exp = N_exp_bg + RD_wimp*sigma_p
            N_obs = N_exp  ! Observed=Expected for Asimov data
            X_in1= (/log10(sigma_p),R_bg/)
            X_in0 = R_bg
            step0 = R_bg_err*R_bg

            if (sum(RD_wimp*sigma_p).gt.1.0d3) then
             call llhood1_bigN(X_in1,L1) ! Asimov data maximises likelihood at correct value
             call llhood0_bigN(X_in0,L0)
             call MINIM(X_in0,step0,n_bg,L0,MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,IQUAD,SIMP,VAR,llhood0_bigN,IFAULT0)
            else
              call llhood1(X_in1,L1) ! Asimov data maximises likelihood at correct value
              call llhood0(X_in0,L0)
              call MINIM(X_in0,step0,n_bg,L0,MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,IQUAD,SIMP,VAR,llhood0,IFAULT0)
            end if

            D01 = -2.0*(L1-L0)
            if (D01.ge.9.0d0) then ! Median 3sigma detection -> D = 9
              ! Do interpolation to find discovery limit cross section
              if (k.eq.1) then
                DL(i,j) = ex_min
                k1 = 1
              else
                DL(i,j) = 10.0d0**(interp1D((/D_prev,D01/),(/log10(ex_prev),log10(Exposure)/),2,9.0d0))
                k1 = k-1
              end if
              Nsig(i,j) = sum(RD_wimp*sigma_p)
              Nbg(i,j) = N_tot_bg
              write(*,*) 'm = ',m_chi,'|| sigma = ',sigma_p,'|| DL = ',DL(i,j),'|| N_obs = ',sum(N_obs),L0,L1
              RD_bg = RD_bg/Exposure
              RD_wimp = RD_wimp/Exposure
              exit
            end if
            ex_prev = Exposure ! Reset for interpolation
            D_prev = D01
            !write(*,*) Exposure,D01
          end if
          RD_bg = RD_bg/Exposure
          RD_wimp = RD_wimp/Exposure
        end do

        if ((k-1).eq.n_ex) then
          DL(i,j:ns) = ex_max
          exit
        end if
      end do

      !open(unit=1000,file='test2.dat')
      !write(1000,*) DL(i,:)
      !write(1000,*) sigma_p_vals
      !close(1000)
      !stop
    end if
  end do
  call UnAllocate ! Reset
end subroutine










!=====================================DISCOVERY LIMIT==========================!
subroutine DiscoveryLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL)
	! * Calculates discovery limits (DL) from m_chi = m_min to m_max
	! * Parameters of likelihood are X = (sigma_p,R_bg)
	! * Assumes signal model RD_wimp and background model RD_bg
	! * RD_wimp = (Number of WIMP events in each bin)/sigma_p
	! * RD_bg = (Number of BG events each bin)/R_bg
	! * Parameters enter likelihood function to rescale RD_wimp and RD_bg
	! * Each mass scans cross sections from sigma_min to sigma_min
	! * Done this way because the minimisation is fastest for low event numbers
	! * Uses Asimov asymptotic result to get median limit for 3sigma detection
	integer :: nf,nm,i,ns,im,j,MAXFUNEVALS,IPRINT,NLOOP,IQUAD,ifault0,si,ii
	double precision :: sigma_p_vals(ns),DL(nm),m_vals(nm),m_min,m_max,sigma_min,sigma_max
	double precision :: x_in0(n_bg),x_in1(n_bg+1),step0(n_bg),N_exp(nTot_bins),N_exp_bg(nTot_bins)
	double precision :: D01,L1,L0,SIMP,STOPCR0,var(2),N_tot_bg,D_prev,s_prev

	! SET UP MINIMIMISATION
  MAXFUNEVALS = 10000 ! Maximum function evaluations
  IPRINT = -1 ! Print results from minimisaition
  NLOOP = 1 ! Number of iterations before looping
  IQUAD = 0 ! Can't remeber what this does
  SIMP = 0.1 ! Nor this
	STOPCR0 = 1.0d-6 ! Accuracy of max likelihood

	! Mass and cross section discretisation
  m_vals = logspace(m_min,m_max,nm)
	sigma_p_vals = logspace(sigma_min,sigma_max,ns)
	DL = 0.0d0

	! GENERATE BACKGROUND DATA
	N_exp_bg = 0.0d0
	do si = 1,n_bg
		N_exp_bg = N_exp_bg + R_bg(si)*RD_bg(:,si)
	end do
	N_tot_bg = sum(N_exp_bg)




	! MASS SCAN:
	do im = 1,nm
		m_chi = m_vals(im)
		call WIMPRecoilDistribution	! Call WIMP recoil distribution for each new mass
		! CROSS SECTION SCAN
		do j = 1,ns
			sigma_p = sigma_p_vals(j)
			if (sum(RD_wimp*sigma_p).gt.0.1d0) then	! Generally need >0.5 events to see DM



				N_exp = N_exp_bg + RD_wimp*sigma_p
				N_obs = N_exp  ! Observed=Expected for Asimov data

			  !------ Signal + Background ------!
				! X_in = (log10(sigma_p), background rates)
			  X_in1(1) = log10(sigma_p)
				X_in1(2:n_bg+1) = R_bg
			    call llhood1(X_in1,L1) ! Asimov data maximises likelihood at correct value

			  !------ Background only ------!
			  X_in0 = R_bg
			  step0 = R_bg_err*R_bg
				call llhood0(X_in0,L0)
			    call MINIM(X_in0,step0,n_bg,L0,&
				  MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,&
				  IQUAD,SIMP,VAR,llhood0,IFAULT0)

				! Test statistic
				D01 = -2.0*(L1-L0)
				if (D01.ge.9.0d0) then ! Median 3sigma detection -> D = 9
					! Do interpolation to find discovery limit cross section
					DL(im) = 10.0d0**(interp1D((/D_prev,D01/),(/log10(s_prev),log10(sigma_p)/),2,9.0d0))
					exit
				end if
				s_prev = sigma_p ! Reset for interpolation
				D_prev = D01
			end if
		end do
		write(*,*) 'CYG-DL:',im,'m = ',m_chi,'DL = ',DL(im),'Signal:',sum(RD_wimp*sigma_p),'BG:',N_tot_bg
		!stop
	end do
end subroutine






!===================================LIKELIHOODS================================!
!---------------------------------SIGNAL+BACKGROUND----------------------------!
 subroutine llhood1(X,  LL)
    ! input: Parameters X = (log10(sigma_p),background rates(1:n_bg))
	! output: -1*LogLikelihood = LL
	double precision :: X(n_bg+1),LL,N_exp0(nTot_bins),N_exp1(nTot_bins)
	integer :: i
	! Background events
	N_exp0 = 0.0d0 ! Expected number of events
	do i = 1,n_bg
		N_exp0 = N_exp0 + X(i+1)*RD_bg(:,i) ! Sum over backgrounds
	end do

	! Signal events
	N_exp1 = N_exp0 + RD_wimp*(10.0d0**X(1)) ! Add signal events sig


	! LL = log(Poiss. for N_obs events) + log(Gauss. for R_bg normalisations)
	LL = -1.0*(lnPF(nTot_bins,N_obs,N_exp1)+lnGF(X(2:n_bg+1),R_bg,R_bg_err*R_bg))

 end subroutine llhood1

!---------------------------------SIGNAL+BACKGROUND----------------------------!
subroutine llhood0(X,  LL)
    ! input: Parameters X = background rates(1:n_bg)
	! output: -1*LogLikelihood = LL
    double precision :: X(n_bg),LL,N_exp0(nTot_bins)
	integer :: i
	N_exp0 = 0.0d0 ! Expected number of events
	do i = 1,n_bg
		N_exp0 = N_exp0 + X(i)*RD_bg(:,i) ! Sum over backgrouds
	end do

	! LL = log(Poiss. for N_obs events) + log(Gauss. for R_bg normalisations)
	LL = -1.0*(lnPF(nTot_bins,N_obs,N_exp0)+lnGF(X,R_bg,R_bg_err*R_bg))
  end subroutine llhood0

!===================================USEFUL SUMS================================!
function lnPF(nbins,Nob,Nex)! SUM OF LOG(POISSON PDF)
	! Uses log-gamma function to generalise the factorial
	! to non-integer number of observed events (i.e. Asimov data)
	double precision :: lnPF,Nex(nbins),Nob(nbins)
	integer :: ii,nbins
	lnPF = 0.0d0
	do ii = 1,nbins
		lnPF = lnPF + Nob(ii)*log(Nex(ii)) - Nex(ii) - lgamma(Nob(ii)+1.0)
	end do
end function

function lnGF(Rob,Rex,Rer) ! SUM OF LOG(GAUSSIAN PDF)
	! for x=Rob, mu = Rex, sig = Rer
	double precision :: Rob(:),Rex(:),Rer(:),lnGF
	lnGF = sum(-1.0d0*log(Rer)-0.5d0*log(2.0d0*pi)&
	       	 -(Rob-Rex)**2.0d0/(2.0d0*Rer**2.0d0))
end function


!---------------------------------Large N limit----------------------------!
subroutine llhood1_bigN(X,  LL)
   ! input: Parameters X = (log10(sigma_p),background rates(1:n_bg))
 ! output: -1*LogLikelihood = LL
 double precision :: X(n_bg+1),LL,N_exp0(nTot_bins),N_exp1(nTot_bins)
 integer :: i
 ! Background events
 N_exp0 = 0.0d0 ! Expected number of events
 do i = 1,n_bg
   N_exp0 = N_exp0 + X(i+1)*RD_bg(:,i) ! Sum over backgrounds
 end do

 ! Signal events
 N_exp1 = N_exp0 + RD_wimp*(10.0d0**X(1)) ! Add signal events sig

 LL = 0.5*sum((N_obs-N_exp1)**2.0/N_exp1)
 LL = LL-lnGF(X(2:n_bg+1),R_bg,R_bg_err*R_bg)

end subroutine llhood1_bigN

subroutine llhood0_bigN(X,  LL)
   ! input: Parameters X = background rates(1:n_bg)
 ! output: -1*LogLikelihood = LL
   double precision :: X(n_bg),LL,N_exp0(nTot_bins)
 integer :: i
 N_exp0 = 0.0d0 ! Expected number of events
 do i = 1,n_bg
   N_exp0 = N_exp0 + X(i)*RD_bg(:,i) ! Sum over backgrouds
 end do

 LL = -1.0*(-0.5*sum((N_obs-N_exp0)**2.0/N_exp0)+lnGF(X,R_bg,R_bg_err*R_bg))
end subroutine llhood0_bigN


!---------------------------------Large N limit----------------------------!
subroutine llhood1_hugeN(X,  LL)
   ! input: Parameters X = (log10(sigma_p),background rates(1:n_bg))
 ! output: -1*LogLikelihood = LL
 double precision :: X(n_bg+1),LL,N_exp0(nTot_bins),N_exp1(nTot_bins)
 integer :: i
 ! Background events
 N_exp0 = 0.0d0 ! Expected number of events
 do i = 1,n_bg
   N_exp0 = N_exp0 + X(i+1)*RD_bg(:,i) ! Sum over backgrounds
 end do

 ! Signal events
 N_exp1 = N_exp0 + RD_wimp*(10.0d0**X(1)) ! Add signal events sig

 LL = 0.5*sum((N_obs-N_exp1)**2.0/N_exp1)
 LL = LL-lnGF(X(2:n_bg+1),R_bg,R_bg_err*R_bg)

end subroutine llhood1_hugeN

subroutine llhood0_hugeN(X,  LL)
   ! input: Parameters X = background rates(1:n_bg)
 ! output: -1*LogLikelihood = LL
   double precision :: X(n_bg),LL,N_exp0(nTot_bins)
 integer :: i
 N_exp0 = 0.0d0 ! Expected number of events
 do i = 1,n_bg
   N_exp0 = N_exp0 + 10.0**X(i)*RD_bg(:,i) ! Sum over backgrouds
 end do

 LL = -1.0*(-0.5*sum((N_obs-N_exp0)**2.0/N_exp0)+lnGF(10.0**X,R_bg,R_bg_err*R_bg))
end subroutine llhood0_hugeN


end module like
