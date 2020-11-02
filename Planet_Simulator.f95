program Planet_Orbit
   ! Kelly Lewis
   ! Planet Orbit Simulation
    implicit none
    integer :: j, i, n
	double precision :: Planet(4)
	real*8 :: kinetic, potiential, fP(4) , totalenergy, h, total_time, KE, PE, time
	real*8 :: k1(4), k2(4), k3(4), k4(4) 
	

	! ask for initial variables, read directly into array
	print*, "Initial x-position: "
	read*, Planet(1)
	print*, "Initial x-velocity: "
	read*, Planet(2)
	print*, "Initial y-position: "
	read*,Planet(3)
	print*, "Initial y-velocity: "
	read*, Planet(4) 
	print*, "Time increment: "
	read*, h
	print*, "Total Time: "
	read*, total_time
	
	PE = potiential(Planet(1), Planet(3))
	KE = kinetic(Planet(2) , Planet(4))
	totalenergy = kinetic(Planet(2) , Planet(4))  + potiential(Planet(1), Planet(3))
	
	! print initial energies
	print*, "Initial Kinetic Energy: ", KE
	print*, "Initial Potential Energy: ", PE
	print*, "Initial Total Energy: " , totalenergy
	! blank line
	print*, "    "
	
	! open file
	open(45, file= "planet_bonus.txt")
	
	! number of iterations
	n = total_time/h
	
	! Simulation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Forward Euler 
	do j = 1, n
		!call fPrime(Planet, fP)
		call fPrime(Planet, k1)
		call fPrime(Planet, k2)
		call fPrime(Planet, k3)
		call fPrime(Planet, k4)
		
		do i = 1, 4
			!Planet(i) = Planet(i) + h * fP(i)
			Planet(i) = Planet(i) + (h/6)*(k1(i)+2*k2(i)+2*k3(i)+k4(i))
		enddo 
		
		PE = potiential(Planet(1), Planet(3))
		KE = kinetic(Planet(2) , Planet(4))
		totalenergy = kinetic(Planet(2) , Planet(4))  + potiential(Planet(1), Planet(3))
		time = time + h
		
		! write to the file
		write(45,*) Planet(1), Planet(3), PE, KE , totalenergy, time
	enddo
	
	
	! print final energies
	print*, "Final Kinetic Energy: ", KE 
	print*, "Final Potential Energy: ", PE
	print*, "Final Total Energy: ", totalenergy

	! close file 
	close(45) 

end program Planet_Orbit

! f prime subfunction 
subroutine fPrime(Planet, fP)
	implicit none
	real*8 :: Planet(4), fP(4), r
	r = sqrt(Planet(1)**2 + Planet(3)**2) 
	fP(1) = Planet(2) 
	fP(2) = - Planet(1) / r**3
	fP(3) = Planet(4) 
	fP(4) = - Planet(3) / r**3 

end subroutine fPrime


function kinetic(vx, vy) result(k) 
	implicit none 
	real*8 :: k, vx, vy
	k = 0.5*(vx**2 + vy**2) 

end function

function potiential(x, y) result(p) 
	implicit none 
	real*8 :: p , x, y
	p = -1.0 / sqrt(x**2 + y**2)

end function
