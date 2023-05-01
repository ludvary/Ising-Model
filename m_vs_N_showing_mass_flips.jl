# After the lattice reaches equilibrium state, all the spins flip from -1 to +1 and vice versa together in a mass at equal interval. The frequency of this mass spin flip depends on the temperature and the box size
using Printf
function Ising()
	mcs_steps, M, T, step_size, L, index = 1000000, 0, 2.1, 10, 5, 0 #initializing values (for mass flips to happen, the temperature should be near the critical temp which for this system is about 2.2)

	lattice = rand([1, -1], L, L) # random lattice config

	# initializing the arrays
	time_arr = zeros(Int64, div(mcs_steps ,step_size)) 
	M_arr = zeros(Float64, div(mcs_steps ,step_size))
	M = sum(lattice)

	# printing initial magnetization per spin
	println("mag is $(@sprintf("%.30f", M/(L*L)))") 

	function MCS(lattice, L, T, M)
		for m in 1:mcs_steps
			i, j = rand(1:L), rand(1:L) # choosing a spin at random

			# finding the intial value of energy of the randomly selected spin
			eInit = -(lattice[i, j] *(lattice[i, mod1(j+1, L)] + lattice[i, mod1(j-1, L)]+ lattice[mod1(i+1, L), j]+ lattice[mod1(i-1, L), j] ))

			# trail flip
			lattice[i, j] = -lattice[i, j]

			# finding the energy after trial flip
			efinal = -((lattice[i, j]) *(lattice[i, mod1(j+1, L)] + lattice[i, mod1(j-1, L)]+ lattice[mod1(i+1, L), j]+ lattice[mod1(i-1, L), j] ))

			# finding the difference
			delE = efinal - eInit
			
			# applying Metropolis-Hastings
			if delE < 0
				M += 2*lattice[i, j] # find instantaneous magnetization and accept the flip
			else
				prob = exp(-delE/T)
				random = rand()
				if random < prob
					M += 2*lattice[i, j] # find instantaneous magnetization and accept the flip
				else
					lattice[i, j] = -lattice[i, j] # reject the flip
				end
			end
			if m%step_size == 0 # not collecting data in every monte carlo step but collecting it every 10 steps
				index += 1
				time_arr[index] = m 
				M_arr[index] = M/(L*L)
			end
		end
		return time_arr, M_arr, lattice, M
	end
	
	time_arr, M_arr, lattice, M = MCS(lattice, L, T, M) # calling the MCS() on lattice

	# writing to the dat file
	file = open("spin_flips_in_mass.dat", "w")
	for a in 1:length(time_arr)
		@printf(file, "%d\t%.30f\n", time_arr[a], M_arr[a])
	end
	close(file)
end
println(@elapsed Ising())
