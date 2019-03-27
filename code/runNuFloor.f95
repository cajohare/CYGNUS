program runNuFloor
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax,i1,i2,i,n_ex,n_Eth
	double precision :: m_min,m_max,sigma_min,sigma_max,Vol,Time,ex_min,ex_max
	double precision,dimension(:),allocatable :: m_vals
	double precision,dimension(:),allocatable :: DL
	double precision :: Eth_max,Eth_min,EB8
	character(len=100) :: filename1,filename2,filename3
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)

	!-----Set binning ------------------------------!
	nE_bins = 200 ! Number of energy bins		  		!
	nT_bins = 1	 ! Number of time bins					  !
	nside = 0  ! Order of pixelation (2,4 or 8)		!
	!---------

	!---------Set precision
	nm = 1000
	m_min = 0.1
	m_max = 1.0d4
	ns = 500 ! resolution of cs scan							!
	sigma_min = 1.0d-50 ! min. expected cs				!
	sigma_max = 1.0d-38 ! max expected cs					!
	allocate(m_vals(nm))
	allocate(DL(nm))
  !---------

	nucleus = Xenon131
	call GetNuFloor_Simple(1.0d4,m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL)
	filename1 = '../data/WIMPLimits/SI/nuFloorSimple-Xe.txt'
	open(unit=1000,file=trim(filename1))
	write(1000,*) m_vals
	write(1000,*) DL

	nucleus = Fluorine
	call GetNuFloor_Simple(1.0d4,m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL)
	filename1 = '../data/WIMPLimits/SI/nuFloorSimple-F.txt'
	open(unit=1000,file=trim(filename1))
	write(1000,*) m_vals
	write(1000,*) DL

	nucleus = Helium
	call GetNuFloor_Simple(1.0d4,m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL)
	filename1 = '../data/WIMPLimits/SI/nuFloorSimple-He.txt'
	open(unit=1000,file=trim(filename1))
	write(1000,*) m_vals
	write(1000,*) DL

end program runNuFloor
