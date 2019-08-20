program runCYGNUS_thresholds
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax,n_ex
	double precision :: m_min,m_max,sigma_min,sigma_max,ex_min,ex_max,Vol,Time,m,m1,m2,m3
	double precision,dimension(:),allocatable :: m_vals,DLHe,DLF
	double precision :: E_th_F_edisc,E_th_He_edisc,E_th_F_1e,E_th_He_1e,E_th_He_dir,E_th_F_dir
	double precision,dimension(:),allocatable :: ex_vals,DL1_Ex,DL2_Ex,DL3_Ex,Nbg1,Nsig1,Nbg2,Nsig2,Nbg3,Nsig3,DLHe_1wimp,DLF_1wimp
	logical :: saveNwimp
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)


! 	write(*,*) '====================================================='
! 	write(*,*) 'Select run (0-10):'
! 	write(*,*) '0 --- Idealised detector (3D)'
! 	write(*,*) '1 --- Pre drift grid readout (3D)'
! 	write(*,*) '2 --- Post drift grid readout (3D)'
! 	write(*,*) '3 --- Gem readout'
! 	write(*,*) '4 --- Wire readout'
! 	write(*,*) '5 --- Pad readout'
! 	write(*,*) '6 --- Strip readout'
!	write(*,*) '7 --- Pixel readout'
! 	write(*,*) '8 --- Nondirectional'
! 	write(*,*) '9 --- RUN ALL READOUTS'
! 	write(*,*) '10 --- RUN ALL PIXEL READOUTS'

	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!
	  nE_bins = 20 ! Number of energy bins		  		
	  nT_bins = 1	 ! Number of time bins					  
	  nside = 4  ! Order of pixelation (2,4 or 8)		
	!----- Mass range for limits -------------------!
	  nm = 200 !  Number of mass points							
	  m_min = 0.1d0 ! Min mass											
	  m_max = 10000.0d0 ! Max mass										
	!----- Cross section range for limits ----------!
	  ns = 300 ! resolution of cs scan							
	  sigma_min = 1.0d-50 ! min. expected cs				
	  sigma_max = 1.0d-38 ! max expected cs					

	Time = 6.0
	VolTime = 1000.0*Time
	lat = Boulby(1)
	long = Boulby(2)

	E_th_F_1e = 0.26
	E_th_He_1e = 0.26
	E_th_F_dir = 8.0
	E_th_He_dir = 8.0
	
	energy_on = .true. ! energy info is currently turned on for best limits
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .true.
	energyres_on = .true.
	searchmode = .false.
	saveNwimp = .true.
	
	
	fn_end = 'Predrift'
	readout = 1 ! 
	VolTime = 10.0*Time
	call LoadReadout(readout,	fn_end)
	call CYGNUSLimit(0.25d0,0.25d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS10-'//trim(fn_end)//'_250eV.txt',saveNwimp)
	call LoadReadout(readout,	fn_end)
	call CYGNUSLimit(8.0d0,8.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS10-'//trim(fn_end)//'_8keV.txt',saveNwimp)
	
	
	fn_end = 'Predrift'
	readout = 1 ! 
	VolTime = 10.0d3*Time
	call LoadReadout(readout,	fn_end)
	call CYGNUSLimit(0.25d0,0.25d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_250eV.txt',saveNwimp)
	call LoadReadout(readout,	fn_end)
	call CYGNUSLimit(8.0d0,8.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_8keV.txt',saveNwimp)
	
	
	fn_end = 'Predrift'
	readout = 1 ! 
	VolTime = 10.0d5*Time
	call LoadReadout(readout,	fn_end)
	call CYGNUSLimit(0.25d0,0.25d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS100k-'//trim(fn_end)//'_250eV.txt',saveNwimp)
	call LoadReadout(readout,	fn_end)
	call CYGNUSLimit(8.0d0,8.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS100k-'//trim(fn_end)//'_8keV.txt',saveNwimp)
	
	
	
!
! 	fn_end = 'Predrift'
! 	readout = 1 !
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(0.25d0,0.25d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_250eV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(1.0d0,1.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_1keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(2.0d0,2.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_2keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(3.0d0,3.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_3keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(4.0d0,4.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_4keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(5.0d0,5.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_5keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(6.0d0,6.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_6keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(7.0d0,7.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_7keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(8.0d0,8.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_8keV.txt',saveNwimp)


	
! 	fn_end = 'Postdrift'
! 	readout = 2 !
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(0.25d0,0.25d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_250eV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(1.0d0,1.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_1keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(2.0d0,2.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_2keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(3.0d0,3.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_3keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(4.0d0,4.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_4keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(5.0d0,5.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_5keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(6.0d0,6.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_6keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(7.0d0,7.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_7keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(8.0d0,8.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_8keV.txt',saveNwimp)
!
!
!
! 	fn_end = 'Strip'
! 	readout = 6 !
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(0.25d0,0.25d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_250eV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(1.0d0,1.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_1keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(2.0d0,2.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_2keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(3.0d0,3.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_3keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(4.0d0,4.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_4keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(5.0d0,5.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_5keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(6.0d0,6.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_6keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(7.0d0,7.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_7keV.txt',saveNwimp)
! 	call LoadReadout(readout,	fn_end)
! 	call CYGNUSLimit(8.0d0,8.0d0,m_min,m_max,nm,sigma_min,sigma_max,ns,'../data/CYGNUS1000-'//trim(fn_end)//'_8keV.txt',saveNwimp)



end program runCYGNUS_thresholds
