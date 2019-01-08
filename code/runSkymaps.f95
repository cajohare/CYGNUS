program runSkymaps
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	use pix_tools

	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax,i1,i2,i,pix
	double precision :: m_min,m_max,sigma_min,sigma_max,Vol,Time,x_pix_db(3)
	double precision,dimension(:),allocatable :: m_vals,DLHe,DLF
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)


	!open(unit=123,file='../pixels/xpix32.txt')
	!nside = 32
 	!npix = 12*nside**2
 	!do pix = 0,npix-1
  	!call pix2vec_ring(nside,pix,  x_pix_db)
		!		write(123,*) x_pix_db
 	!end do
	!close(123)
	!stop

	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!
	  nE_bins = 40 ! Number of energy bins		  		!
	  nT_bins = 1	 ! Number of time bins					  !
	  nside = 32  ! Order of pixelation (2,4 or 8)		!
	!----- Mass range for limits -------------------!
	  nm = 200 !  Number of mass points							!
	  m_min = 0.5d0 ! Min mass											!
	  m_max = 10000.0d0 ! Max mass										!
	!----- Cross section range for limits ----------!
	  ns = 500 ! resolution of cs scan							!
	  sigma_min = 1.0d-50 ! min. expected cs				!
	  sigma_max = 1.0d-41 ! max expected cs					!
	!-----------------------------------------------!
	!-----------------------------------------------!
	allocate(m_vals(nm))
	allocate(DLF(nm))
	allocate(DLHe(nm))

	VolTime = 1000.0*3.0d0
	Exposure = VolTime*(0.16/1000.0d0) ! Convert m^3-years into ton-years
	lat = Boulby(1)
	long = Boulby(2)
	E_max = 20.0d0
	E_th = 3.0d0
	nucleus = Fluorine
	readout = 1
	energy_on = .true. ! energy info is currently turned on for best limits
	energyres_on = .false.
	searchmode = .false.
	efficiency_on = .true.

	call GetNuFluxes ! Load Neutrinos
	call PreAllocate ! Allocate data size (readout dependent)
	call SHM ! Load halo model
	!v_lab_all(1,:) = -v_pec - (/0.0d0,v_LSR,0.0d0/)
	m_chi = 10.0

	write(*,*) T_bin_centers(1)+Jan1
	write(*,*) T_bin_centers(1)
	write(*,*) LabVelocity(T_bin_centers(1))
	write(*,*) SolarDirection(T_bin_centers(1))

	angres_on = .false.
	headtail_on = .false.
	call LoadReadout(readout,	fn_end)
	call BackgroundRecoilDistribution ! Load Background	model
	call WIMPRecoilDistribution	! Call WIMP recoil distribution for each new mass
	open(unit=1000,file='RD_wimp_0deg.dat')
	open(unit=2000,file='RD_nu_0deg.dat')
	i1 = 1
	do i = 1,npix!
    i2 = i1+nE_bins-1
    write(1000,*) RD_wimp(i1:i2)
		write(2000,*) RD_bg(i1:i2,1)
    i1 = i2+1
 	end do
	close(1000)


	angres_on = .true.
	headtail_on = .false.
	call LoadReadout(readout,	fn_end)
	call BackgroundRecoilDistribution ! Load Background	model
	call WIMPRecoilDistribution	! Call WIMP recoil distribution for each new mass
	open(unit=1000,file='RD_wimp_30deg.dat')
	open(unit=2000,file='RD_nu_30deg.dat')
	i1 = 1
	do i = 1,npix!
    i2 = i1+nE_bins-1
    write(1000,*) RD_wimp(i1:i2)
		write(2000,*) RD_bg(i1:i2,1)
    i1 = i2+1
 	end do
	close(1000)


	angres_on = .true.
	headtail_on = .true.
	call LoadReadout(readout,	fn_end)
	call BackgroundRecoilDistribution ! Load Background	model
	call WIMPRecoilDistribution	! Call WIMP recoil distribution for each new mass
	open(unit=1000,file='RD_wimp_30deg_HT.dat')
	open(unit=2000,file='RD_nu_30deg_HT.dat')
	i1 = 1
	do i = 1,npix!
    i2 = i1+nE_bins-1
    write(1000,*) RD_wimp(i1:i2)
		write(2000,*) RD_bg(i1:i2,1)
    i1 = i2+1
 	end do
	close(1000)


end program runSkymaps
