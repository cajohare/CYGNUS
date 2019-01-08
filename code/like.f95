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
subroutine CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	! Limits from m_chi = m_min to m_max, with nm values
	! Cross section scan from sigma_min to sigma_max with ns values
	! DL = Discovery Limit
	double precision :: m_min,m_max,sigma_min,sigma_max,m_vals(nm),DLF(nm),DLHe(nm)
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
    ! 1000 m^3 of SF6 at 200 torr or He at 740 torr is 1.6 tons
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
    E_th = 1.8d0
  	E_max = 200.0d0
  	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)
  end if


	! Calculate Fluorine limits
  nucleus = Fluorine
  E_th = 3.0d0
	E_max = 200.0d0
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



end module like
