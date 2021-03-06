module NeutrinoFuncs
  use params
  use LabFuncs
  use util
  implicit none

contains

!================================NeutrinoFuncs.f95==========================================!
!
! Contents:
! 1. BackgroundRecoilDistribution: Generate Background model RD_bg = (Num of bg events)/Rate
!
! 2. Neutrino event rates and fluxes
! GetNuFluxes = Determines which neutrino backgrounds are needed and loads them in
!
!===========================================================================================!



!===========================================================================================!
subroutine BackgroundRecoilDistribution
	double precision :: RD(nTot_bins_full,n_bg),RD_red(nTot_bins)
	integer :: i,i1,i2,ii,j,s,k

	RD_bg = 0.0d0
	RD = 0.0d0
	RD_red = 0.0d0

	! Load all Neutrino RDs
	call NeutrinoRD(nTot_bins_full,RD)

	! Get total rates and rescale RD_bg by them
	do s = 1,n_bg
		R_bg(s) = sum(RD(:,s))
		RD(:,s) = RD(:,s)/R_bg(s)
	end do

	if (sum(RD).gt.1.0d10) then
		WRITE(*,*) 'one'
		stop
	end if

	! If directional, also smear by angular resolution
	if (nside.gt.0) then
		if (angres_on) then
			do s = 1,2
				call SmearRD(RD(:,s))
			end do
		end if
	end if

	if (sum(RD).gt.1.0d10) then
		WRITE(*,*) 'TWO'
		ii = 1
		do i = 1,npix
			do j = 1,nE_bins
				write(*,*) i,j,RD(ii,1:2)
				ii = ii+1
			end do
		end do
		stop
	end if

	! If direction only then integrate over energies
  	if (Energy_on) then
  		RD_bg = RD
  	else
  		do s = 1,n_bg
  			call IntegrateOverEnergies(RD(:,s),RD_red)
  		    RD_bg(:,s) = RD_red
  		end do
  	end if


	! Multiply whole thing by Exposure so RD = Num events/R_bg
	RD_bg = RD_bg*Exposure

!  open(unit=2000,file='RD_nu_1.txt')
!  do i = 1,nTot_bins
!    write(2000,*) RD_bg(i,:)
  !end do

end subroutine BackgroundRecoilDistribution
!---------------------------------------------------------------------------------------------!






!===============================Neutrino Recoil distributions=================================!
subroutine NeutrinoRD(n1,RD) ! Generates an RD for all neutrinos
	double precision :: RD(n1,n_bg),fE_r1,fE_r2,E_r1,E_r2,fEmin_r1,fEmin_r2,dpix,E_r
	double precision :: dRdE(nE_bins),R_tot,Flux_t,wid,R,R_s
	integer :: i,j,si,n1,i1,i2,ii,k,n_nu,ia
	double precision,dimension(nbins_full) :: E_full,dRdE_full,dRdE_full_s,f_s,eff,eff_HT,sig_E
	wid = (log10(E_upper)-log10(E_lower))/(1.0d0*nbins_full-1.0d0)
	do i = 1,nbins_full
	   E_full(i) = E_lower*10.0d0**((i-1)*wid)
	end do
	
	RD = 0.0D0
	n_nu = n_bg

	eff = Efficiency(E_full,nbins_full)
	sig_E = EnergyResolution(E_full,nbins_full)

	!---------------------------  NON-DIRECTIONAL ----------------------!
	if (nside.eq.0) then
		do si = 1,n_nu			
			! Correct for energy resolution
		    do ia = 1,nbins_full
				dRdE_full(ia) = eff(ia)*NeutrinoRecoilEnergySpectrum(E_full(ia),E_nu_all(:,si),Flux_all(:,si))
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
			
			E_r1 = E_bin_edges(1)
		    call interp1(E_full,dRdE_full_s,nbins_full,E_r1,   fE_r1)
			ii = 1
			do j = 1,nE_bins
		  		E_r2 = E_bin_edges(j+1)
		        call interp1(E_full,dRdE_full_s,nbins_full,E_r2,   fE_r2)
			  	dRdE(ii) = (E_r2-E_r1)*(fE_r1 + fE_r2)/2.0d0
			  	E_r1 = E_r2
			  	fE_r1 = fE_r2
			  	ii = ii+1
			end do		   

			! Correct for annual modulation
			ii = 1
			do i = 1,nT_bins
				i1 = ii
				i2 = i1+nE_bins-1
				if (si.lt.n_nu-1) then
					! Solar neutrinos:
					! Flux is scaled by 1/EarthSunDistance^2 but since Flux is already averaged
					! We need to also divide by Integral(1/R^2) over one year
					! Integral_inv_EarthSun_sq is defined in params.f95
					RD(i1:i2,si) = (dRdE/Integral_inv_EarthSun_sq)&
							*(1.0d0/EarthSunDistance(T_bin_centers(i))**2.0d0)
				else
					! dsnb and atm neutrinos:
					! Currently am not incorporating any modulation since the
					! event rate is nonexistent anyway
					RD(i1:i2,si) = dRdE
				end if
				ii = i2+1
			end do
		RD(:,si) = RD(:,si)/(1.0d0*nT_bins)
	    end do
	end if





	!---------------------------  DIRECTIONAL ---------------------------!
	if (nside.gt.0) then
		! Load Head-tail efficiency
		eff_HT = HeadTailEfficiency(E_full,nbins_full)

	    ! Solar neutrinos (set each pixel to energy only spectrum/npixels*ntimes)
	    dpix = 4*pi/(npix*1.0d0)
	    do si = 1,n_nu-2
  			ii = 1
  			do i = 1,nT_bins
  			    do k = 1,npix
					
					! Correct for energy resolution
					do ia = 1,nbins_full
						dRdE_full(ia) = eff(ia)*(eff_HT(ia)*NeutrinoRecoilSpectrum_Solar(E_full(ia)*x_pix(k,:),i,E_nu_all(:,si),Flux_all(:,si)) + &
  		 		   			(1.0-eff_HT(1))*NeutrinoRecoilSpectrum_Solar(-1.0d0*E_full(ia)*x_pix(k,:),i,E_nu_all(:,si),Flux_all(:,si)))
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
						RD(ii,si) = (dpix/(1.0d0*nT_bins))*(E_bin_edges(j+1)-E_bin_edges(j))*(fE_r1 + fE_r2)/2.0
						fE_r1 = fE_r2
						ii = ii+1
					end do
  				  end do
  			end do
		  end do
		  
		  
		  

	    ! Isotropic neutrinos (set each pixel to energy only spectrum/npixels*ntimes)
	    do si = (n_nu-1),(n_nu) ! last two background are always isotropic (DSNB+Atm)
		    do ia = 1,nbins_full
				dRdE_full(ia) = eff(ia)*NeutrinoRecoilEnergySpectrum(E_full(ia),E_nu_all(:,si),Flux_all(:,si))
		    end do
			R = sum((E_full(2:nbins_full)-E_full(1:nbins_full-1))*0.5d0*(dRdE_full(1:nbins_full-1) + dRdE_full(2:nbins_full)))
			
			if ((sum(sig_E).gt.0.0).and.(sum(dRdE_full).gt.0.0)) then
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
			
			E_r1 = E_bin_edges(1)
		    call interp1(E_full,dRdE_full_s,nbins_full,E_r1,   fE_r1)
			ii = 1
			do j = 1,nE_bins
		  		E_r2 = E_bin_edges(j+1)
		        call interp1(E_full,dRdE_full_s,nbins_full,E_r2,   fE_r2)
			  	dRdE(ii) = (E_r2-E_r1)*(fE_r1 + fE_r2)/2.0d0
			  	E_r1 = E_r2
			  	fE_r1 = fE_r2
			  	ii = ii+1
			end do	
			
		   ii = 1
			do i = 1,nT_bins
				do k = 1,npix
					i1 = ii
		 		   	i2 = ii+nE_bins-1
		  		 	RD(i1:i2,si) = dRdE/(nT_bins*npix*1.0d0)
		  			ii = i2+1
		  		end do
			end do
	    end do
	end if

end subroutine NeutrinoRD
!---------------------------------------------------------------------------------------------!







!========================================Neutrino data========================================!
subroutine GetNuFluxes
	! Reads each neutrino flux data file
	! Each flux file has 1000 rows apart from monochromatic ones
	! the energies are stored in E_nu_all, fluxes in Flux_all
	integer :: i,nvals,ii,j
	integer :: sel(11),s
	double precision :: E_r_max

	! CURRENTLY USING BARCELONA FLUXES
	NuFlux = NuFlux_B17GS98
	NuUnc = NuUnc_B17GS98
	nvals = 1000 ! All neutrino backgrounds saved with 1000 entries


	! ORDER OF NEUTRINOS
	! 1. pp
	! 2. pep
	! 3. hep
	! 4. 7Be [343 keV]
	! 5. 7Be [843 keV]
	! 6. 8B
	! 7. 13 N
	! 8. 15 O
	! 9. 17 F
	! 10. DSNB
	! 11. Atm

	! Monochromatic neutrinos (2, 4, 5) have a negative value for E_nu which is
	! Used to tell the rate formula to use monochromatic result

	! Figure out which backgrounds give recoils above E_th
	sel = 0
	n_bg = 0
	do i = 1,11
		E_r_max = MaxNuRecoilEnergy(i) ! Max recoil energy for neutrino
		if (E_r_max.gt.E_th) then
			sel(i) = i ! Set sel if neutrino needs to be loaded
			n_bg = n_bg+1
		end if
		write(*,*) E_r_max,E_th,sel(i),i
	end do


	! Load in all backgrounds with sel != 0
	allocate(E_nu_all(1000,n_bg))
	allocate(Flux_all(1000,n_bg))
	allocate(R_bg_err(n_bg))
	allocate(R_bg(n_bg))
	E_nu_all = 0.0d0
	Flux_all = 0.0d0
	ii = 1
	do i = 1,11
		s = sel(i)
		if (s.gt.0) then
		    if (s.eq.1) then
		       open(unit=10,file='../neutrinos/pp-1000.txt')
		    else if (s.eq.2) then
		       open(unit=10,file='../neutrinos/pep-1000.txt')
		    else if (s.eq.3) then
		       open(unit=10,file='../neutrinos/hep-1000.txt')
		    else if (s.eq.4) then
		       open(unit=10,file='../neutrinos/7Be1-1000.txt')
		    else if (s.eq.5) then
		       open(unit=10,file='../neutrinos/7Be2-1000.txt')
		    else if (s.eq.6) then
		       open(unit=10,file='../neutrinos/8B-1000.txt')
		    else if (s.eq.7) then
		       open(unit=10,file='../neutrinos/13N-1000.txt')
		    else if (s.eq.8) then
		       open(unit=10,file='../neutrinos/15O-1000.txt')
		    else if (s.eq.9) then
		       open(unit=10,file='../neutrinos/17F-1000.txt')
		    else if (s.eq.10) then
		       open(unit=10,file='../neutrinos/DSNB-1000.txt')
		    else if (s.eq.11) then
		       open(unit=10,file='../neutrinos/Atm-1000.txt')
		    end if
		    do j = 1,nvals
		      read(10,*) E_nu_all(j,ii),Flux_all(j,ii)
		    end do
		    close(10)
		    Flux_all(:,ii) = NuFlux(i)*Flux_all(:,ii) ! Select rate normalisations
			R_bg_err(ii) = NuUnc(i) ! Select rate normalisation uncertainties
			ii = ii+1
		end if
	end do

end subroutine GetNuFluxes

 !---------------------------------------------------------------------------------------------!
function MaxNuRecoilEnergy(s) result(E_r_max) ! Max recoil energy for neutrino number s
	integer :: s
	double precision :: E_r_max,m_N_keV
	m_N_keV = 0.93141941*(nucleus(1)+nucleus(2))*1.0e6
	E_r_max = 2*m_N_keV*(1000.0*NuMaxEnergy(s))**2.0&
	   /(m_N_keV+1000*NuMaxEnergy(s))**2.0
end function

!===================================nu spectra=================================!
function NeutrinoRecoilEnergySpectrum(E_r,E_nu,Flux) result(dRdE)
	integer :: N,Z,nn,i
	double precision :: Q_W,E_r,dRdE,m_N_keV,m_N_GeV
	double precision, dimension(1000) :: E_nu,Flux,diff_sigma,dR
	nn = 1000
	N = nucleus(1) ! number of neutrons
	Z = nucleus(2) ! number of protons
	Q_W = N-(1-4.0d0*sintheta_Wsq)*Z ! weak nuclear hypercharge
	m_N_GeV = 0.93141941d0*(N+Z) ! nucleus mass in GeV
	m_N_keV = m_N_GeV*1.0d6 ! nucleus mass in keV

	! differential cross section
	diff_sigma = (G_F_GeV**2.0d0/(4.0d0*pi))*(Q_W**2.0d0)*&
	     m_N_GeV*(1.0d0-(m_N_keV*E_r)/(2.0d0*(E_nu*1000.0d0)**2.0d0))*&
	     (0.197d-13)**2.0d0*(1.0d-6)*1000.0d0/(1.0d0*N+1.0d0*Z)*(N_A)*&
	     FormFactorHelm(E_r,N+Z)**2.0d0

	! diff_sigma goes negative for kinematically forbidden energies so just set them 0 and sum
	where (diff_sigma.lt.0.0d0)
	   diff_sigma=0.0d0
	end where

	! Event rate is integral(cross section x flux)
	if (Flux(2).gt.0.0) then
		dR = diff_sigma*Flux
	 	dRdE = sum(0.5d0*(E_nu(2:nn)-E_nu(1:nn-1))*(dR(2:nn)+dR(1:nn-1)))
	else
		dRdE = diff_sigma(1)*Flux(1)*E_nu(1) ! for monochromatic nu's
	end if

	! Convert into /ton/year/keV
	dRdE = dRdE*(365.0d0*3600.0d0*24d0)*(1000.0)
end function NeutrinoRecoilEnergySpectrum

!=================================Directional spectrum=========================================!
function NeutrinoRecoilSpectrum_Solar(E3,tbin,E_nu,Flux_in) result(dR)
	double precision :: E3(3),t,Q_W,E_r,x_r(3),x_sun(3),dRdEdO(1),F0,F1,E0,E1,F_value
	double precision :: E_nu_min,Eps,costh,diff_sigma,aa,bb,costh_r,dR,E_r_max
	double precision :: m_N_keV,m_N_GeV
	integer :: i,N,Z,nn,ii,tbin
	double precision, dimension(1000) :: E_nu,Flux,Flux_in
	nn = 1000

	E_r = sqrt(sum(E3**2.0)) ! Recoil energy
	x_r = E3/E_r ! Recoil direction

	x_sun = 1.0*SolarDirection(T_bin_centers(tbin)) ! Direction of sun at tbin
  Flux = (Flux_in/Integral_inv_EarthSun_sq)*(1.0/EarthSunDistance(t)**2.0) ! Modulated flux

	dR = 0.0d0
	N = nucleus(1)
	Z = nucleus(2)
	Q_W = N-(1-4*sintheta_Wsq)*Z ! Weak nuclear hypercharge
	m_N_GeV = 0.93141941*(N+Z)
	m_N_keV = m_N_GeV*1.0e6
	m_N_keV = 0.9315*1.0d6*(N+Z) ! keV

  if (Flux(2).gt.0.0) then
	   E_r_max = 2*m_N_keV*(1000.0*E_nu(nn))**2.0/(m_N_keV+1000.0*E_nu(nn))**2.0
  else
     E_r_max = 2*m_N_keV*(1000.0*E_nu(1))**2.0/(m_N_keV+1000.0*E_nu(1))**2.0
  end if

	costh = sum(-(x_sun*x_r)) ! Angle between sun and recoil
	Eps = 0.0
	if (E_r.le.E_r_max) then ! recoil energy is permitted
	  E_nu_min = sqrt(m_N_keV*E_r/2.0)	! Min nu energy to create
	  if (Flux(2).gt.0.0d0)	 then
	      ! CHROMATIC NEUTRINOS
	      if (costh.gt.(E_nu_min/m_N_keV)) then
		  Eps = 1.0/(costh/E_nu_min - 1.0/m_N_keV)
		  if (Eps.le.1000.0*E_nu(nn)) then

		    diff_sigma = (G_F_GeV**2.0/(4*pi))*(Q_W**2.0)*&
			  m_N_GeV*(1.0-(m_N_keV*E_r)/(2.0*Eps**2.0))*&
			  (0.197e-13)**2.0*(1.0e-6)*1000/(N+Z)*(N_A)*&
			  FormFactorHelm(E_r,sum(nucleus))**2.0

		    ! INTERP
		    F_value = 0.0
		    if (Eps.gt.(1000.3650*E_nu(1))) then
			do i = 2,nn
			  if ((1000.0*E_nu(i)).gt.Eps) then
			      F0 = Flux(i-1)
			      F1 = Flux(i)
			      E0 = 1000.0*E_nu(i-1)
			      E1 = 1000.0*E_nu(i)
			      F_value = (F0+(F1-F0)*(Eps-E0)/(E1-E0))
			      exit

			  end if
			end do
		    else
			F_value = Flux(1)
		    end if
		    dRdEdO = diff_sigma*(F_value/1000)*(Eps**2.0)/(E_nu_min)
		  else
		    dRdEdO = 0.0
		  end if
	      else
		  dRdEdO = 0.0
	      end if
	  else
		  ! MONOCHROMATIC NEUTRINOS
	      costh_r=((E_nu(1)*1000.0d0+m_N_keV)/(1000.0d0*E_nu(1)))&
	           *sqrt(E_r/(2*m_N_keV))
	      if (abs(costh-costh_r)<0.01) then
	         Eps = 1000*E_nu(1)
	         diff_sigma = (G_F_GeV**2.0/(4*pi))*(Q_W**2.0)*&
	              m_N_GeV*(1.0-(m_N_keV*E_r)/(2.0*Eps**2.0))*&
	              (0.197e-13)**2.0*(1.0e-6)*1000/(N+Z)*(N_A)*&
	              FormFactorHelm(E_r,N+Z)**2.0
	         dRdEdO = diff_sigma*(Flux(1)/1000.0)*E_nu(1)*1000.0
	      else
	         dRdEdO = 0.0
	      end if
	  end if
	else
	   dRdEdO = 0.0
	end if

	dRdEdO = dRdEdO*(365.0d0*3600.0d0*24d0)*(1000)/(2.0*pi)
	dR = dR + dRdEdO(1)

end function NeutrinoRecoilSpectrum_Solar

end module NeutrinoFuncs
