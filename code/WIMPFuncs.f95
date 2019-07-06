module WIMPFuncs
  use params
  use util
  use LabFuncs
  implicit none

contains


!================================WIMPFuncs.f95==========================================!
! Just contains a bunch of boring but useful/essential stuff

! Contents:
! 1. Halo models
!    SHM: Loads in old SHM parameters
!    SHMpp: Loads in updated SHM paramers

! 2. Loading recoil distribution for WIMPs (RD_wimp)
!    WIMPRecoilDistribution: Loads RD_wimp to be used by likelihood
!    WIMPRD_Energy: Energy dependent recoil distribution
!    WIMPRD_3D: Direction dependent recoil distribution


! 3. Energy dependnt and direction dependent rates
!    WIMPRate_Direction: Directional recoil rate d^2R(t)/dEdO
!    WIMPRate_Energy: Non-directional recoil rate dR(t)/dE
!	 MaxWIMPEnergy: Maximum recoil energy for a given date
!=====================================================================================!



!=====================================================================================!
!------------------------------ OLD STANDARD HALO MODEL ------------------------------!
subroutine SHM
	integer :: bin
	v_LSR = 220.0d0
	sig_v = v_LSR/sqrt(2.0d0)
	v_esc = 544.0d0
	rho_0 = 0.3d0
	!if (nT_bins.gt.1) then
  	  	do bin = 1,nT_bins
		   v_lab_all(bin,:) = LabVelocity(T_bin_centers(bin))
    	end do
		!else
	!	 v_lab_all(1,:) = v_pec + (/0.0d0,v_LSR,0.0d0/)
	!end if
end subroutine

!----------------------------- NEW STANDARD HALO MODEL ---------------------------------!
subroutine SHMpp
	integer :: bin
	v_LSR = 235.0d0
	sig_v = v_LSR/sqrt(2.0d0)
	v_esc = 580.0d0
	rho_0 = 0.55d0
	if (nT_bins.gt.1) then
  	  	do bin = 1,nT_bins
		   v_lab_all(bin,:) = LabVelocity(T_bin_centers(bin))
    	end do
	else
		v_lab_all(1,:) = v_pec + (/0.0d0,v_LSR,0.0d0/)
	end if
end subroutine


!====================================================================================!
!---------- Load in full WIMP recoil distribution for arbitrary binning -------------!
subroutine WIMPRecoilDistribution
	double precision :: RD(nTot_bins_full),RD_red(nTot_bins)
	integer :: i,i1,i2,ii
	! Set a fiducial cross section just so numbers aren't crazy
  !m_chi = 2.0
	sigma_p = 1.0d-45
	! Choose whether directional or non-directional is used
	ii = 1
	do i = 1,nT_bins
		i1 = ii
		if (nside.gt.0) then
			i2 = i1+nE_bins*npix-1
			call WIMPRD_3D(RD(i1:i2),i)
		else
			i2 = i1+nE_bins-1
			call WIMPRD_Energy(RD(i1:i2),i)
		end if
		ii = i2+1
	end do

	! If directional, also smear by angular resolution
	if (angres_on) then
		if (nside.gt.0) then
			if (sum(RD).gt.0.0d0) then
				call SmearRD(RD)
			end if
		end if
	end if

	! If direction only then integrate over energies
  if (nside.gt.0) then
  	if (Energy_on) then
  		RD_wimp = RD
  	else
  		call IntegrateOverEnergies(RD,RD_red)
  		RD_wimp = RD_red
  	end if
  else
    RD_wimp = RD
  end if


	! Multiply whole thing by Exposure so RD = Num events/sigma_p
  !write(*,*) 'N_tot = ',sum(RD_wimp)
	RD_wimp = RD_wimp*Exposure/(nT_bins*sigma_p)

!  open(unit=1000,file='RD_wimp_1.dat')
!  i1 = 1
!  do i = 1,npix!
!    i2 = i1+nE_bins-1
!    write(1000,*) RD_wimp(i1:i2)
!    i1 = i2+1
  !end do
  !stop

end subroutine WIMPRecoilDistribution

!-------------------- Directional recoil distribution --------------------------------!
subroutine WIMPRD_3D(RD,tbin)
	double precision :: fE_r1,fE_r2,E_r1(3),E_r2(3)
	double precision :: dpix,RD(nE_bins*npix),E_r,wid,R,R_s
	double precision,dimension(nbins_full) :: E_full,dRdE_full,dRdE_full_s,f_s,eff,eff_HT,sig_E
	integer :: i,j,k,ii,tbin,ia
	wid = (log10(E_upper)-log10(E_lower))/(1.0d0*nbins_full-1.0d0)
	do i = 1,nbins_full
	   E_full(i) = E_lower*10.0d0**((i-1)*wid)
	end do
	
	RD = 0.0
	dpix = 4*pi/(npix*1.0d0)
	eff = Efficiency(E_full,nbins_full)
	sig_E = EnergyResolution(E_full,nbins_full)
	eff_HT = HeadTailEfficiency(E_full,nbins_full)

	ii = 1
	do i = 1,nT_bins
	    do k = 1,npix
		
			! Correct for energy resolution
			do ia = 1,nbins_full
				dRdE_full(ia) = eff(ia)*(eff_HT(ia)*WIMPRate_Direction(E_full(ia)*x_pix(k,:),tbin) + &
 		   			(1.0-eff_HT(1))*WIMPRate_Direction(-1.0d0*E_full(ia)*x_pix(k,:),tbin))
			end do
			R = sum((E_full(2:nbins_full)-E_full(1:nbins_full-1))*0.5d0*(dRdE_full(1:nbins_full-1) + dRdE_full(2:nbins_full)))
			
			if ((sum(sig_E).gt.0.0).and.(R.gt.0.0)) then
				do ia = 1,nbins_full
				   E_r = E_full(ia)
				   f_s = 1.0d0/(sqrt(2.0d0*pi)*sig_E)&
				        *exp(-0.5d0*((E_r-E_full)**2.0/sig_E**2.0d0))
				   dRdE_full_s(ia) = sum((E_full(2:nbins_full)-E_full(1:nbins_full-1))&
				        *0.5d0*(dRdE_full(1:nbins_full-1)*f_s(1:nbins_full-1) + dRdE_full(2:nbins_full)*f_s(2:nbins_full)))
				end do
				R_s = sum((E_full(2:nbins_full)-E_full(1:nbins_full-1))*0.5d0*(dRdE_full_s(1:nbins_full-1) + dRdE_full_s(2:nbins_full)))
				dRdE_full_s = dRdE_full_s*R/R_s
			else
				dRdE_full_s = dRdE_full
			end if
	
			call interp1(E_full,dRdE_full_s,nbins_full,E_bin_edges(1),fE_r1)
			do j = 1,nE_bins
				call interp1(E_full,dRdE_full_s,nbins_full,E_bin_edges(j+1),fE_r2)
				RD(ii) = (dpix/(1.0d0*nT_bins))*(E_bin_edges(j+1)-E_bin_edges(j))*(fE_r1 + fE_r2)/2.0
				fE_r1 = fE_r2
				ii = ii+1
			end do
		  end do
	end do
	
end subroutine WIMPRD_3D

!-------------------- Non-directional recoil distribution-------------------------------------!
subroutine WIMPRD_Energy(RD,tbin)
	double precision,dimension(nbins_full) :: E_full,dRdE_full,dRdE_full_s,f_s,eff_full
	double precision :: fE_r1,fE_r2,E_r1,E_r2,RD(nE_bins),E_r,eff(nE_bins+1),sig_E(nbins_full),wid,R,R_s
	integer :: i,j,ii,nbins,ia,tbin
	wid = (log10(E_upper)-log10(E_lower))/(1.0d0*nbins_full-1.0d0)
	do i = 1,nbins_full
	   E_full(i) = E_lower*10.0d0**((i-1)*wid)
	end do
	
	
	eff = Efficiency(E_bin_edges,nE_bins+1)
	if (energyres_on) then
		eff_full = Efficiency(E_full,nbins_full)
		sig_E = EnergyResolution(E_full,nbins_full)

	    do ia = 1,nbins_full
			dRdE_full(ia) = eff_full(ia)*WIMPRate_Energy(E_full(ia),tbin)
	    end do
		R = sum((E_full(2:nbins_full)-E_full(1:nbins_full-1))*0.5d0*(dRdE_full(1:nbins_full-1) + dRdE_full(2:nbins_full)))
		if (R.gt.0.0) then
		    do ia = 1,nbins_full
		       E_r = E_full(ia)
		       f_s = 1.0d0/(sqrt(2.0d0*pi)*sig_E)&
		            *exp(-0.5d0*((E_r-E_full)**2.0/sig_E**2.0d0))
		       dRdE_full_s(ia) = sum((E_full(2:nbins_full)-E_full(1:nbins_full-1))&
		            *0.5d0*(dRdE_full(1:nbins_full-1)*f_s(1:nbins_full-1) + dRdE_full(2:nbins_full)*f_s(2:nbins_full)))
		    end do
			R_s = sum((E_full(2:nbins_full)-E_full(1:nbins_full-1))*0.5d0*(dRdE_full_s(1:nbins_full-1) + dRdE_full_s(2:nbins_full)))
			dRdE_full_s = dRdE_full_s*R/R_s
		else
			dRdE_full_s = 0.0
		end if
		
		! Compute actual RD
	    RD = 0.0d0
		E_r1 = E_bin_edges(1)
		!fE_r1 = WIMPRate_Energy(E_r1,sig,v_lab)
	    call interp1(E_full,dRdE_full_s,nbins_full,E_r1,   fE_r1)
		ii = 1
		do j = 1,nE_bins
	  		E_r2 = E_bin_edges(j+1)
		  	!fE_r2 = WIMPRate_Energy(E_r2,sig,v_lab)
			!write(*,*) j,E_r2
	        call interp1(E_full,dRdE_full_s,nbins_full,E_r2,   fE_r2)
		  	RD(ii) = (E_r2-E_r1)*(fE_r1 + fE_r2)/2.0d0
		  	E_r1 = E_r2
		  	fE_r1 = fE_r2
		  	ii = ii+1
		end do

	else
		! Compute actual RD
	    RD = 0.0d0
		E_r1 = E_bin_edges(1)
		fE_r1 = WIMPRate_Energy(E_r1,tbin)
		ii = 1
		do j = 1,nE_bins
	  		E_r2 = E_bin_edges(j+1)
		  	fE_r2 = WIMPRate_Energy(E_r2,tbin)
		  	RD(ii) = (E_r2-E_r1)*(fE_r1*eff(j) + fE_r2*eff(j+1))/2.0d0
		  	E_r1 = E_r2
		  	fE_r1 = fE_r2
		  	ii = ii+1
		end do
	end if

end subroutine WIMPRD_Energy





!=====================================================================================!
!-------------------- Energy-Time dependent recoil rate-------------------------------!
function WIMPRate_Energy(E_r,tbin) result(dRdE)
	double precision :: E_r,dRdE,m_p,m_N,m_chi_kg,mu_p,c_cm
	double precision :: x,y,z,v_e,N_esc,gvmin,m_N_keV,mu_N,v_min,v_0,v_lab(3)
	integer :: A,tbin

	! relevant constants
	A = sum(nucleus) ! mass number of nucleus
	m_p = 0.9315*1e6 ! mass of proton
	m_N = m_p*A ! mass of nucleus
	m_chi_kg = m_chi*1e6*1.783d-33 ! mass of dm in kg
	mu_p = 1.0d6*m_chi*m_p/(1.0d6*m_chi + m_p) ! reduced proton mass
	c_cm = 3.0d8*100.0 ! speed of light in cm/s
	m_N_keV = A*0.9315*1.0d6 ! nucleus mass in keV
	mu_N = 1.0d6*m_chi*m_N_keV/(1.0d6*m_chi + m_N_keV) ! reduced nucleus mass
	v_min = (sqrt(2.0d0*m_N_keV*E_r)/(2.0d0*mu_N))*3.0d8/1000.0d0 ! vmin in km/s


	! Compute g(vmin) mean inverse speed halo integral
	N_esc = erf(v_esc/(sqrt(2.0d0)*sig_v))&
	   -sqrt(2.0d0/pi)*(v_esc/sig_v)*exp(-v_esc**2.0d0/(2.0d0*sig_v**2.0d0))
	v_0 = sqrt(2.0d0)*sig_v
	x = v_min/v_0
	v_e = sqrt(sum(v_lab_all(tbin,:)**2.0d0))
	y = v_e/v_0
	z = v_esc/v_0
	if (x.lt.abs(y-z)) then
		if (z.gt.y) then
			gvmin = (1.0d0/(2*N_esc*v_0*y))*(erf(x+y)-erf(x-y)&
			     -(4.0d0/sqrt(pi))*y*exp(-z**2.0d0))
		else
			gvmin = (1.0d0/(v_0*y))
		end if
	else if ((abs(y-z).lt.x).and.(x.lt.(y+z))) then
		gvmin = (1.0d0/(2.0d0*N_esc*v_0*y))*(erf(z)-erf(x-y)&
				-(2.0d0/sqrt(pi))*(y+z-x)*exp(-z**2.0d0))
	else
		gvmin = 0.0
	end if
	gvmin = gvmin/(1000.0d0*100.0d0) ! convert to cm^-1 s

	! compute rate
	dRdE = (c_cm**2.0)*((rho_0*1.0d6*A**2.0*sigma_p)/(2*m_chi_kg*mu_p**2.0))*gvmin
	dRdE = dRdE*3600*24*365*1000.0d0 ! convert to per ton-year
	dRdE = dRdE*FormFactorHelm(E_r,A)**2.0d0 ! apply form factor

end function WIMPRate_Energy


!-------------------- Direction dependent recoil rate--------------------------------------!
function WIMPRate_Direction(E,tbin) result(dRdEdO)
	double precision :: Erot(3),E(3),x(3),E_r,dRdEdO,m_p,m_N,m_chi_kg,mu_p,c_cm,fhat
	double precision :: v_min,N_esc,vlabdotq,m_N_keV,mu_N,v_lab(3)
	integer :: A,tbin

	! Calculate relavant directions
	v_lab = v_lab_all(tbin,:)	! lab velocity
	E_r = sqrt(sum(E**2.0)) ! Recoil energy
	x = (/E(1)/E_r,E(2)/E_r,E(3)/E_r/) ! Recoil direction
	vlabdotq = sum(x*v_lab) ! Lab-recoil projection


	! Relevant constants
	A = sum(nucleus)
	m_p = 0.9315*1e6
	m_chi_kg = m_chi*1.0d6*1.783d-33
	mu_p = 1.0d6*m_chi*m_p/(1.0d6*m_chi + m_p)
	c_cm = 3.0d8*100.0
	m_N = m_p*A
	m_N_keV = A*0.9315*1.0d6
	mu_N = 1.0d6*m_chi*m_N_keV/(1.0d6*m_chi + m_N_keV)
	v_min = (sqrt(2.0d0*m_N_keV*E_r)/(2.0d0*mu_N))*3.0d8/1000.0d0

	! Compute Radon transform fhat
	N_esc = erf(v_esc/(sqrt(2.0)*sig_v))&
	   -sqrt(2.0/pi)*(v_esc/sig_v)*exp(-v_esc**2.0/(2*sig_v**2.0))
	if ((v_min+vlabdotq).le.(v_esc)) then
		fhat = (1/(N_esc*sqrt(2*pi*sig_v**2.0)))*&
	          (exp(-(v_min+vlabdotq)**2.0/(2*sig_v**2.0))-&
	          exp(-v_esc**2.0/(2*sig_v**2.0)))
	else
		fhat = 0.0d0
	end if
	fhat = fhat/(1000.0d0*100.0d0) ! convert to cm^-1 s

	! Compute rate
	dRdEdO = (c_cm**2.0)*((rho_0*1.0d6*A**2.0*sigma_p)/(4*pi*m_chi_kg*mu_p**2.0))*fhat
	dRdEdO = dRdEdO*3600*24*365*1000.0d0 ! convert to per ton-year
	dRdEdO = dRdEdO*FormFactorHelm(E_r,A)**2.0d0 ! correct for form factor

end function WIMPRate_Direction


!=================================Max Energy==================================!
function MaxWIMPEnergy(A,v_lab) result(E_max_lim)
	double precision :: E_max_lim,mu_N,m_N,v_lab(3)
	integer :: A
	m_N = 0.9315*1.0d6*A
	mu_N = 1.0d6*m_N*m_chi/(1.0d6*m_chi+m_N)
	E_max_lim = 2.0d0*mu_N**2.0d0*((v_esc+sqrt(sum(v_lab**2.0)))*1000.0/3.0d8)**2.0d0/m_N
end function MaxWIMPEnergy

end module WIMPFuncs
