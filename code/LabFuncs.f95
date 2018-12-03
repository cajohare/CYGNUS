module LabFuncs
	use params
	use util

	implicit none

contains


!================================LabFuncs.f95==================================!
! Contains all functions needed to compute lab dependent and detector dependent effects

! Contents:

! 1. Nuclear Form Factors
! FormFactorHelm: Helm form factor, currently the only one I can be bothered with

! 2. Detector performance
! EnergyResolution: charge detection resolution sig_E [keV]
! Efficiency: charge detection efficiency eff = (0 -> 1)
! AngularResolution: energy dependent angular resolution curve in radians
! HeadTailEfficiency: prob of correct sense recognition eff_HT = (0 -> 1)
! LoadReadout: loads readout detector performance data files

! 3. Detector directional integrals
! IntegrateOverEnergies: Converts an RD into direction-only data
! SmearRD: Applies the preloaded angular resolution curve sig_gamma to the RD
! Smear: The function for actually smearing the angular distribution
! GaussianKernel: The function for gaussian smoothing on a sphere

! 4. Lab velocity
! LabVelocitySimple: Simple lab velocity in Galactic system without Earth rotation
! LabVelocity: Full lab velocity in (N,W,Z) with Earth rotation
! JulianDay: JulianDay at dd-mm-yyyy hh:hh

! 5. Solar direction
! EarthSunDistance: Distance between Earth and Sun as a function of time
! SolarDirection: Direction of the sun at a given time

! 6. Co-ordinate transformations
! gal2eqt: Galactic system to equatorial system
! gal2lab: Galactic system to lab system
! Rgal2lab: Full galactic system to lab system
! lab2gal: Lab system to galactic system
!==============================================================================!






!====================================Form Factors==============================!
!----------------------------------------HELM----------------------------------!
function FormFactorHelm(E_r,A)
	integer :: A
	double precision :: E_r,FormFactorHelm,q,c1,s,R_1
	q = sqrt(2*A*931.5*1000*E_r)*1.0d-12/1.97d-7

	c1 = 1.23d0*A**(1.0d0/3.0d0)-0.6d0
	s = 0.9d0
	R_1 = sqrt(c1**2 + (7.0d0/3.0d0)*pi**2.0d0*(0.52d0**2.0d0) - 5*s**2.0)

	FormFactorHelm = (3*(sin(q*R_1) - q*R_1*cos(q*R_1))*exp(-q*q*s*s/2.0)/(q*R_1)**3)
	if (q.eq.0.0d0) then
	FormFactorHelm = 1.0d0
	end if
end function FormFactorHelm
!------------------------------------------------------------------------------!




!==============================Detector Performanc=============================!
!----------------------------Load Detector Resolution curve--------------------!
function EnergyResolution(E_r,ni) result(sig_E)
	integer :: i,nuc,ni
	double precision :: E_r(ni),sig_E(ni)
	if (energyres_on) then
		if (nucleus(1).eq.Fluorine(1)) then
			nuc = 1
		elseif (nucleus(1).eq.Helium(1)) then
			nuc = 2
		end if
		do i = 1,ni
			sig_E(i) = E_vals(i)*interp1D(E_vals,energyres_data(:,nuc),1000,E_r(i))
		end do
	else
		sig_E = 1.0d0
	end if
end function

!---------------------------------Load  Efficiency curve-----------------------!
function Efficiency(E_r,ni) result(eff)
	integer :: i,nuc,ni
	double precision :: E_r(ni),eff(ni),s1,Ec,EC2
	if (efficiency_on) then
		if (nucleus(1).eq.Fluorine(1)) then
			nuc = 1
		elseif (nucleus(1).eq.Helium(1)) then
			nuc = 2
		end if
		do i = 1,ni
			eff(i) = interp1D(E_vals,efficiency_data(:,nuc),1000,E_r(i))
		end do
	else
		eff = 1.0d0
	end if
end function

!----------------------------Load Angular resolution curve---------------------!
function AngularResolution(E_r,ni) result(sig_gamma)
	integer :: i,nuc,ni
	double precision :: E_r(ni),sig_gamma(ni)
	if (angres_on) then
		if (nucleus(1).eq.Fluorine(1)) then
			nuc = 1
		elseif (nucleus(1).eq.Helium(1)) then
			nuc = 2
		end if
		do i = 1,ni
			sig_gamma(i) = interp1D(E_vals,angres_data(:,nuc),1000,E_r(i))*(pi/180.0d0)
			if (sig_gamma(i)>1.0) then
					sig_gamma(i) = 0.99
			end if
		end do
	else
		sig_gamma = 0.0d0
	end if
end function

!----------------------------Load Head-tail efficiency curve-------------------!
function HeadTailEfficiency(E_r,ni) result(eff_HT)
	integer :: i,ni,nuc
	double precision :: E_r(ni),eff_HT(ni)
	if (headtail_on) then
		if (nucleus(1).eq.Fluorine(1)) then
			nuc = 1
		elseif (nucleus(1).eq.Helium(1)) then
			nuc = 2
		end if
		do i = 1,ni
			eff_HT(i) = interp1D(E_vals,headtail_data(:,nuc),1000,E_r(i))
		end do
	else
		eff_HT = 1.0d0
	end if
end function


!----------------------------Load Everything-----------------------------------!
subroutine LoadReadout(ro, ReadoutName)
	integer :: ro,i
	character(len=100) :: ReadoutName
	double precision :: Ei

	! Set readout name
	if (ro.eq.0) then
		ReadoutName = 'Ideal'
	elseif (ro.eq.1) then
		ReadoutName = 'Pixel'
	elseif (ro.eq.2) then
		ReadoutName = 'Predrift'
	elseif (ro.eq.3) then
		ReadoutName = 'Postdrift'
	elseif (ro.eq.4) then
		ReadoutName = 'Strip'
	elseif (ro.eq.5) then
		ReadoutName = 'Optical'
	elseif (ro.eq.6) then
		ReadoutName = 'Wire'
	elseif (ro.eq.7) then
		ReadoutName = 'Pad'
	elseif (ro.eq.8) then
		ReadoutName = 'Planar'
	elseif (ro.eq.9) then
		ReadoutName = 'Nondirectional'
	end if

	! Load all the data files needed
	open(unit=1300,file='../readouts/energyres/'//trim(ReadoutName)//'-EnergyRes.txt')
	open(unit=1301,file='../readouts/efficiency/'//trim(ReadoutName)//'-Efficiency.txt')
	open(unit=1302,file='../readouts/angres/'//trim(ReadoutName)//'-AngRes.txt')
	open(unit=1303,file='../readouts/headtail/'//trim(ReadoutName)//'-HeadTail.txt')

	! Allocate all data to the right arrays
	do i = 1,1000
		read(1300,*) E_vals(i),energyres_data(i,:)
		read(1301,*) Ei,efficiency_data(i,:)
		read(1302,*) Ei,angres_data(i,:)
		read(1303,*) Ei,headtail_data(i,:)
	end do

	! Remeber to close the file or you will let "him" out
	close(1300)
	close(1301)
	close(1302)
	close(1303)

end subroutine
!------------------------------------------------------------------------------!





!======Directional performance integrals needed to modify recoil distribitions=!
subroutine IntegrateOverEnergies(RD,RD_red)
	double precision :: RD(nTot_bins_full),RD_red(nTot_bins)
	integer :: i1,i2,i,k
	i1 = 1
	do i = 1,nT_bins
		do k = 1,npix
			i2 = i1+nE_bins-1
			RD_red = sum(RD(i1:i2))
			i1 = i2+1
		end do
	end do
end subroutine

!-------------------------Smear the full Energy-Time-Direction RD--------------!
subroutine SmearRD(RD)
	integer :: i,ii,k,ibin
	double precision :: sig_gamma(nE_bins),RD(nTot_bins_full)
	double precision :: RDpix(npix),RD_smeared(nTot_bins_full)
	! weird loop structures just getting
	! round the stupid way I set up the binning order
	sig_gamma = AngularResolution(E_bin_centers,nE_bins)

	RD_smeared = RD
	do ibin = 1,nE_bins
		do i = 1,nT_bins
			do k = 1,npix
				ii = (i-1)*npix*nE_bins + (k-1)*nE_bins + ibin
				RDpix(k) = RD(ii)
			end do
		end do

		if (sum(RDpix).gt.0.0d0) then
			call Smear(RDpix,sig_gamma(ibin))

			do i = 1,nT_bins
				do k = 1,npix
					ii = (i-1)*npix*nE_bins + (k-1)*nE_bins + ibin
					!write(*,*) ibin,i,k,ii,RD(ii),RDpix(k)
					RD_smeared(ii) = RDpix(k)
					!write(*,*) k,RD(ii),RD_smeared(ii),sig_gamma(ibin)
				end do
			end do
		end if
	end do
	RD = RD_smeared
end subroutine

!-------------------------Smear just one set of pixels-------------------------!
subroutine Smear(RDpix_in,sig_gammai)
	double precision :: RDpix(npix),RD_K(npix),RDpix_in(npix)
	double precision :: sig_gammai,Rtot,x0(3),Rpixtot
	integer :: k,k2
	! Smears pixelised angular recoil distribution by angular res. sig_gamma
	! assumes that angular distribution is at a single energy
	Rpixtot = sum(RDpix_in) ! used to normalise the smearing
	RDpix = RDpix_in
	do k = 1,npix
		x0 = x_pix(k,:)
		RD_K = 0.0d0
		do k2 = 1,npix
			RD_K(k2) = RDpix_in(k2)*GaussianKernelCosth(x_pix(k2,:),x0,sig_gammai)
		end do
		RD_K(k) = RDpix_in(k)
		RDpix(k) = sum(RD_K)
	end do
	RDpix = RDpix/sum(RDpix)
	RDpix = RDpix*Rpixtot
	! update RDpix_in
	RDpix_in = RDpix
end subroutine

!------------------Gaussian (theta) smoothing kernal on a sphere---------------!
function GaussianKernel(x,x0,sig_gammai) result(K)
	double precision :: x(3),x0(3),K,gamma,sig_gammai
	gamma = (x(1)*x0(1) + x(2)*x0(2) + x(3)*x0(3))
	if (gamma.gt.1.0d0) then
		gamma = 1.0d0
		end if
	if (gamma.lt.-1.0d0) then
		gamma = -1.0d0
		end if
	gamma = acos(gamma)
	K = exp(-gamma*gamma/(2.0d0*sig_gammai*sig_gammai))
end function
!------------------------------------------------------------------------------!


!------------------Gaussian (costh) smoothing kernal on a sphere---------------!
function GaussianKernelCosth(x,x0,sig_gammai) result(K)
	double precision :: x(3),x0(3),K,cosgamma,cossig,sig_gammai
	cosgamma = (x(1)*x0(1) + x(2)*x0(2) + x(3)*x0(3))
	cossig = 1.0-cos(asin(sig_gammai/1.0))
	! sig_gammai can't be more than 1 radian or this will break
	!if (gamma.gt.1.0d0) then
		!gamma = 1.0d0
		!end if
	!if (gamma.lt.-1.0d0) then
		!gamma = -1.0d0
		!end if
	K = exp(-(1.0-cosgamma)**2.0/(2.0d0*cossig*cossig))
end function
!------------------------------------------------------------------------------!


!=================================Lab Velocity=================================!
!------------------------------------------------------------------------------!
function LabVelocitySimple(day) result(vlab)
	! Simple version in Galactic system without Earth ro.
	double precision :: day,e1(3),e2(3),t1,w,vlab(3)
	w = 2*pi/365.0d0
	t1 = 79.0d0
	e1 = (/0.9941d0,0.1088d0,0.0042d0/)
	e2 = (/-0.0504d0,0.4946d0,-0.8677d0/)
	vlab = vv_earthrev*(cos(w*(day-t1))*e1 + sin(w*(day-t1))*e2) ! Earth rev.
	vlab = vlab + v_pec ! Add peculiar velocity
	vlab(2) = vlab(2) + v_LSR ! add LSR velocity
end function LabVelocitySimple

!------------------------------------------------------------------------------!
function LabVelocity(day) ! More complex version in Lab frame system
	double precision :: day,JD,LabVelocity(3),UT,MJD,T_0,t_GAST,t_lab
	double precision :: vv_galrot,v_galrot(3),v_solar(3),v_earthrot(3)
	double precision :: e,lambda_0,L,g,lambda_sun,beta(3),lambda_i(3),v_earthrev(3)

	! Convert day into phase of Earth rotation t_lab
	JD = day+Jan1
	UT = 24*(JD+0.5-floor(JD+0.5)) ! Universal time
	MJD = JD - 2400000.5 ! Modified Julian Day
	T_0 = (floor(MJD)-55197.5)/36525.0
	t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
	t_lab = t_GAST + long/15
	t_lab = 15*t_lab ! Lab time in degrees

	! Galactic (LSR) Rotation
	vv_galrot = v_LSR
	v_galrot = (/0.0d0,vv_galrot,0.0d0/)
	call gal2lab(v_galrot,t_lab) ! transform to lab co-ords

	! Peculiar solar Motion
	v_solar = v_pec
	call gal2lab(v_solar,t_lab) ! transform to lab co-ords

	! Earth's revolution (first calculate in galactic frame then transform)
	e = eccentricity
	lambda_0 = orb_long_ecliptic
	L = 281.0298 + 36000.77*T_0 + 0.04107*UT
	g = 357.9258 + 35999.05*T_0 + 0.04107*UT
	lambda_sun = L + (1.915 - 0.0048*T_0)*sin(g*pi/180.0)&
				+ 0.020*sin(2*g*pi/180.0)
	beta = lat_ecl_gal
	lambda_i = long_ecl_gal
	v_earthrev = vv_earthrev*(1-e*sin(pi/180.0*(lambda_sun-lambda_0)))*&
				(cos(beta*pi/180.0)*sin(pi/180.0*(lambda_sun-lambda_i)))
	call gal2lab(v_earthrev,t_lab) ! transform to lab co-ords

	! Earth's rotation
	v_earthrot = 0.465102*cos(lat*pi/180)*(/0.0,-1.0,0.0/) ! already in lab co-ords

	! Total
	LabVelocity = v_earthrot+v_earthrev+v_solar+v_galrot

end function LabVelocity

!------------------------------------------------------------------------------!
function JulianDay(month,day,year,hour) ! calculates Julian day from input date
	integer :: month,day,year,year_r,month_r
	double precision :: hour,JulianDay
	year_r = year+4800-floor((14-month)/12.0)
	month_r = month+12*floor((14-month)/12.0)-3
	JulianDay = day + floor((153*month_r+2)/5.0) + 365*year_r &
				+ floor(year_r/4.0) -floor(year_r/100.0) &
				+ floor(year_r/400.0) - 32045 + (hour-12.0)/24.0
end function JulianDay
!------------------------------------------------------------------------------!





!===================================Solar direction============================!
!------------------------------------------------------------------------------!
function EarthSunDistance(day) result(r_es) ! Earth-sun distance at time = day
	double precision :: day,r_es,g,D
	double precision :: JD
	JD = day+Jan1
	D = JD-2451545.0
	g = 357.529 + 0.98560028*D
	g = g*pi/180
	r_es = 1.00014 - 0.01671*cos(g) - 0.00014*cos(2*g)
	r_es = r_es*AstronomicalUnit
end function  EarthSunDistance

!------------------------------------------------------------------------------!
function SolarDirection(day) ! Solar direction in lab coords at time = day
	double precision :: day,n,Omega,L,g,ll,ep,ra,dec,gmst,lmst,JD
	double precision :: omega_hour,theta_z,parallax,gamma,longr,latr
	double precision :: UT,MJD,T_0,t_GAST,t_lab,x_sun(3)
	double precision :: SolarDirection(3)

	! Compute Solar vector in equatorial system
	JD = day+Jan1
	n = JD - 2451545.0
	Omega = 2.1429-0.0010394594*n
	L = 4.8950630 + 0.017202791698*n
	g = 6.2400600 + 0.0172019699*n
	ll = L+0.03341607*sin(g) + 0.00034894*sin(2*g) &
	- 0.0001134 - 0.0000203*sin(Omega)
	ep = 0.4090928 - 6.214e-9*n + 0.0000396*cos(Omega)
	ra = atan2((cos(ep)*sin(ll)),cos(ll)) ! Right ascension of Sun
	dec = asin(sin(ep)*sin(ll)) ! Declination of sun
	x_sun(1) = cos(dec)*cos(ra)
	x_sun(2) = cos(dec)*sin(ra)
	x_sun(3) = sin(dec)

	! Lab time conversion
	UT = 24*(JD+0.5-floor(JD+0.5))
	MJD = JD - 2400000.5
	T_0 = (floor(MJD)-55197.5)/36525.0
	t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
	t_lab = t_GAST + long/15
	t_lab = 15*t_lab ! DEGREES

	! Convert vector from equatorial system into lab system
	call eqt2lab(x_sun,t_lab)
	SolarDirection = x_sun
end function SolarDirection
!------------------------------------------------------------------------------!
!==============================================================================!







!=================================Co-ord transformations=======================!
!------------------------------------------------------------------------------!
subroutine eqt2lab(v,t_lab) ! Equatorial (x_e,y_e,z_e) to Laboratory (N,W,Z)
	double precision:: v(3),t_lab,t,vp(3),latr
	t = t_lab*pi/180.0
	latr = lat*pi/180.0
	vp = v
	v(1) = -cos(t)*sin(latr)*vp(1) - sin(t)*sin(latr)*vp(2) + cos(latr)*vp(3)
	v(2) = sin(t)*vp(1) - cos(t)*vp(2)
	v(3) = cos(t)*cos(latr)*vp(1) + cos(latr)*sin(t)*vp(2) + sin(latr)*vp(3)
end subroutine eqt2lab

!------------------------------------------------------------------------------!
subroutine gal2eqt(v) ! Galactic (x_g,y_g,z_g) to Equatorial (x_e,y_e,z_e)
	double precision :: v(3),vp(3)
	vp = v
	v(1) = -0.06699*vp(1) + 0.4927*vp(2) - 0.8676*vp(3)
	v(2) = -0.8728*vp(1) -0.4503*vp(2) -0.1884*vp(3)
	v(3) = -0.4835*vp(1) + 0.7446*vp(2) + 0.4602*vp(3)
end subroutine gal2eqt

!------------------------------------------------------------------------------!
subroutine gal2lab(v,t_lab) ! Galactic (x_g,y_g,z_g) to Laboratory (N,W,Z)
	double precision :: v(3),t_lab
	call gal2eqt(v)
	call eqt2lab(v,t_lab)
end subroutine gal2lab

!------------------------------------------------------------------------------!
function Rgal2lab(v,day) result(vrot)
	! Full Galactic (x_g,y_g,z_g) to Laboratory (N,W,Z)
	double precision :: vrot(3),v(3),t_lab,JD,day
	double precision :: UT,MJD,T_0,t_GAST
	JD = Jan1 + day
	! Lab time conversion
	UT = 24*(JD+0.5-floor(JD+0.5))
	MJD = JD - 2400000.5
	T_0 = (floor(MJD)-55197.5)/36525.0
	t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
	t_lab = t_GAST + long/15
	t_lab = 15*t_lab ! DEGREES

	vrot = v
	call gal2eqt(vrot)
	call eqt2lab(vrot,t_lab)
end function

!------------------------------------Inverse of above--------------------------!
subroutine lab2gal(v,day)
	double precision :: day,JD,v(3) ! input
	double precision :: UT,MJD,T_0,t_GAST,t_lab,vp(3),latr
	double precision :: M(3,3),l2e(3,3),e2g(3,3),e2l(3,3),g2e(3,3)
	JD = day+Jan1

	! Lab time conversion
	UT = 24*(JD+0.5-floor(JD+0.5)) ! Universal time
	MJD = JD - 2400000.5 ! Modified Julian Day
	T_0 = (floor(MJD)-55197.5)/36525.0
	t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
	t_lab = t_GAST + long/15
	t_lab = 15*t_lab ! Lab time in degrees

	t_lab = t_lab*pi/180.0
	latr = lat*pi/180.0

	g2e(1,:) = (/-0.06699d0, 0.4927d0, -0.8676d0/)
	g2e(2,:) = (/-0.8728d0,-0.4503d0,-0.1884d0/)
	g2e(3,:) = (/-0.4835d0,0.7446d0,0.4602d0/)
	call inverse(g2e,e2g,3)

	e2l(1,:) = (/-sin(latr)*cos(t_lab), -sin(latr)*sin(t_lab), cos(latr)/)
	e2l(2,:) = (/sin(t_lab), -cos(t_lab),0.0d0/)
	e2l(3,:) = (/cos(latr)*cos(t_lab), cos(latr)*sin(t_lab), sin(latr)/)
	call inverse(e2l,l2e,3)
	M =  matmul(e2g,l2e)
	vp = v
	v(1) = M(1,1)*vp(1) + M(1,2)*vp(2) + M(1,3)*vp(3)
	v(2) = M(2,1)*vp(1) + M(2,2)*vp(2) + M(2,3)*vp(3)
	v(3) = M(3,1)*vp(1) + M(3,2)*vp(2) + M(3,3)*vp(3)
end subroutine lab2gal
!------------------------------------------------------------------------------!




end module LabFuncs
