module like	
  use params
  use util  
  use NeutrinoFuncs
  use WIMPFuncs
  use LabFuncs
  implicit none

contains
	
	
	subroutine CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
		double precision :: m_min,m_max,sigma_min,sigma_max,m_vals(nm),DLF(nm),DLHe(nm)
		integer :: i,nm,nf,ns
		character(len=*) :: filename
		
	    write(*,*) '----------------------------------------------------'
	    write(*,*) 'Start CYGNUS Limits:'
		write(*,*) 'Experiment:',trim(filename)
		
		! CYGNUS is at Boubly
		lat = Boulby(1)
		long = Boulby(2)
		
		! Calculate exposure for specified TPC Volume x Time
		Exposure = VolTime*(0.16/1000.0d0) ! Convert m^3-years into ton-years
		
		! Calculate Fluorine limits
	    nucleus = Fluorine
	    E_th = 3.0d0
		E_max = 200.0d0
	    call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,DLF) 
   
   
   	 	! Calculate Helium limits
	    nucleus = Helium
	    E_th = 1.12d0
		E_max = 200.0d0
		call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,DLHe)
		
		
		! Save Data
		open(unit=123,file=trim(filename))
		write(123,*) m_vals
		write(123,*) DLF
		write(123,*) DLHe	
		close(123)	
		write(*,*) 'writing to: ',trim(filename)
	    write(*,*) '----------------------------------------------------'
	end subroutine


	subroutine GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,DL)
		double precision :: m_min,m_max,sigma_min,sigma_max,m_vals(nm),DL(nm)
		integer :: i,nm,nf,ns
	    write(*,*) 'Nucleus = ',nucleus,'Exposure = ',Exposure,'ton years'
	    write(*,*) '[E_th,E_max] = ',E_th,E_max,'keV'
	    write(*,*) '----------------------------------------------------'				
	    call GetNuFluxes ! Load Neutrinos	
	    call PreAllocate ! Allocate data size (readout dependent)		
	    call BackgroundRecoilDistribution ! Load Background	model				
	    call SHM ! Load halo model
		
	    call Sensitivity_Asimov(sigma_min,sigma_max,ns,m_min,m_max,nm,m_vals,DL) ! Discovery limit
		
	    call UnAllocate ! Reset
    end subroutine	

	
	!===================================USEFUL SUMS================================!
	function lnPF(nbins,Nob,Nex)! LOGARITHM OF POISSON PDF
		! Uses log-gamma function to generalise the factorial 
		! to non-integer number of observed events (i.e. Asimov data)
		double precision :: lnPF,Nex(nbins),Nob(nbins)
		integer :: ii,nbins
		lnPF = 0.0d0
		do ii = 1,nbins
			lnPF = lnPF + Nob(ii)*log(Nex(ii)) - Nex(ii) - lgamma(Nob(ii)+1.0)
		end do	
	end function
	
	function lnGF(Rob,Rex,Rer) ! LOGARITHM OF GAUSSIAN PDF
		! for x=Rob, mu = Rex, sig = Rer
		double precision :: Rob(:),Rex(:),Rer(:),lnGF
		lnGF = sum(-1.0d0*log(Rer)-0.5d0*log(2.0d0*pi)&
		        -(Rob-Rex)**2.0d0/(2.0d0*Rer**2.0d0))
	end function
	
	
	!===================================LIKELIHOODS================================!  
	 subroutine llhood1(X,  LL) ! LLHOOD FOR SHM+BACKGROUND HYPOTHESIS
	    double precision :: X(n_bg+1),LL,N_exp0(nTot_bins),N_exp1(nTot_bins)
		integer :: i
		! X = (sigma_p,backgrounds(1:n_bg))		
		! Background events	
		N_exp0 = 0.0d0
		do i = 1,n_bg
			N_exp0 = N_exp0 + X(i+1)*RD_bg(:,i)
		end do
	
		! Signal events	
		N_exp1 = N_exp0 + RD_wimp*(10.0d0**X(1))
		
		! logLike = log(Poisson factor) + log(Gaussian factor)
		LL = -1.0*(lnPF(nTot_bins,N_obs,N_exp1)+lnGF(X(2:n_bg+1),R_bg,R_bg_err*R_bg))
	 end subroutine llhood1
  
	 subroutine llhood0(X,  LL) ! LLHOOD FOR BACKGROUND ONLY HYPOTHESIS
	    double precision :: X(n_bg),LL,N_exp0(nTot_bins)
		integer :: i
		! X = backgrounds(1:n_bg)
		N_exp0 = 0.0d0
		do i = 1,n_bg
			N_exp0 = N_exp0 + X(i)*RD_bg(:,i)
		end do
		LL = -1.0*(lnPF(nTot_bins,N_obs,N_exp0)+lnGF(X,R_bg,R_bg_err*R_bg))
	  end subroutine llhood0
  
	!===============================SENSITIVITY LIMIT================================!  
	subroutine Sensitivity_Asimov(sigma_p_lower,sigma_p_upper,ns,m_lower,m_upper,nm,m_vals,DL)
		integer :: nf,nm,i,ns,im,j,MAXFUNEVALS,IPRINT,NLOOP,IQUAD,ifault0,si,ii
		double precision :: sigma_p_vals(ns),DL(nm),m_vals(nm),x_in0(n_bg),x_in1(n_bg+1),step0(n_bg),var(2),R_bg_obs(n_bg),N_exp(nTot_bins),N_exp_bg(nTot_bins),L1,L0
		double precision :: m_lower,m_upper,E_max_lim,sigma_p_upper,sigma_p_lower,D01,SIMP,STOPCR0,N_tot_bg,D_prev,s_prev	
	    ! SET UP MINIMIMISATION
	    MAXFUNEVALS = 10000
	    IPRINT = -1
	    NLOOP = 1
	    IQUAD = 0
	    SIMP = 0.1
		STOPCR0 = 1.0d-3
		
	    m_vals = logspace(m_lower,m_upper,nm)
		sigma_p_vals = logspace(sigma_p_lower,sigma_p_upper,ns)		
		DL = 0.0d0		
	
		! GENERATE BACKGROUND DATA
		R_bg_obs = R_bg
		N_exp_bg = 0.0d0
		do si = 1,n_bg
			N_exp_bg = N_exp_bg + R_bg_obs(si)*RD_bg(:,si)
		end do
		N_tot_bg = sum(N_exp_bg)
	
		! MASS SCAN:
		do im = 1,nm 
			m_chi = m_vals(im)
			! Recoil distributions must be called for each new mass
			call WIMPRecoilDistribution	
			do j = 1,ns
				sigma_p = sigma_p_vals(j)
				if (sum(RD_wimp*sigma_p).gt.0.5d0) then		
					N_exp = N_exp_bg + RD_wimp*sigma_p
					N_obs = N_exp  

				    ! MODEL 1: (SHM+BG)
				    X_in1(1) = log10(sigma_p)
					X_in1(2:n_bg+1) = R_bg_obs
				    call llhood1(X_in1,L1)	

				    ! MODEL 0: (BG)
				    X_in0 = R_bg_obs
				    step0 = R_bg_err*R_bg_obs
					call llhood0(X_in0,L0)
				    call MINIM(X_in0,step0,n_bg,L0,&
					  MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,&
					  IQUAD,SIMP,VAR,llhood0,IFAULT0)
					D01 = -2.0*(L1-L0)	
												
					if (D01.ge.9.0d0) then
						DL(im) = 10.0d0**(interp1D((/D_prev,D01/),(/log10(s_prev),log10(sigma_p)/),2,9.0d0))
						exit
					end if
					s_prev = sigma_p
					D_prev = D01
					!write(*,*) D01,sum(RD_wimp*sigma_p),X_in0
				end if				
			end do							
				!end if
			write(*,*) 'CYG-DL:',im,'m = ',m_chi,'DL = ',DL(im),'Signal:',sum(RD_wimp*sigma_p),'BG:',N_tot_bg
		end do
	end subroutine
	
  
 

end module like