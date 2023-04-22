using Printf
function Ising()

    mcs_steps::Int64, L::Int64, T::Float64 = 100000, 25, 1
    J::Float64, M::Float64, E::Float64 = 1, 0, 0

    # starting from a random lattice config
    LATTICE = rand([1, -1], L, L)

    # this finds Magnetization of the lattice (not magnetization per spin, just the magnetization)
    function find_magnetization(lattice)
       return sum(lattice) 
    end

    # this finds Energy of the whole lattice (again not per spin)
    function find_energy(lattice)
        E = 0
        for v in 1:L
		    for b in 1:L
	    	    E -= J * (lattice[v, b] * (lattice[v, mod1(b+1, L)] + lattice[v, mod1(b-1, L)] + lattice[mod1(v+1, L), b] + lattice[mod1(v-1, L), b]))
    		end
    	end
        return E
    end
    # this applies Metropolis-Hastingsnwe collect E and M
    function do_simlulation(lattice, T)

        E = find_energy(lattice) # finding the initial energy and magnetization of the lattice
        M = find_magnetization(lattice)
        E_arr = zeros(Float64, mcs_steps) # initializing the arrrays
        M_arr = zeros(Float64, mcs_steps)
        Time_arr = zeros(Float64, mcs_steps)
        for t in 1:mcs_steps
            p, q = rand(1:L), rand(1:L) # choosing random lattice point

            del_E = 2*J * (lattice[p, q] * (lattice[p, mod1(q+1, L)] + lattice[p, mod1(q-1, L)] + lattice[mod1(p+1, L), q] + lattice[mod1(p-1, L), q]))# finding the energy change before and after the randomly choosen spin is flipped. (we could've done this by finding E_initial and flipping the spin and finding E_final and then subtracting both but this way is computationally cheaper)

            # apply Metropolis-Hastings
            if del_E < 0
                lattice[p, q] = -lattice[p, q] # spin flip is accepted
                E += del_E                     # finding instantaneous energy and magnetization
                M += 2*lattice[p, q]
            
            else
                prob = exp(-del_E/T)
                random = rand()
                
                if random < prob
                    lattice[p, q] = -lattice[p, q] # spin flip is accepted
                    E += del_E                     # finding the instantaneous energy and magnetization
                    M += 2*lattice[p, q]
                end
            end
            E_arr[t] = E/(L*L) # storing the values of E and M inside arrays 
            M_arr[t] = M/(L*L)
            Time_arr[t] = t
        end
        return E_arr, M_arr, Time_arr
    end

    E_arr, M_arr, Time_arr = do_simlulation(LATTICE, T)


# writing to dat file
    file = open("lattice_eqm.dat", "w")
    for m in 1:length(Time_arr)   
        @printf(file, "%d\t%.30f\t%.30f\n", Time_arr[m], E_arr[m], M_arr[m])
    end
    close(file)
end
Ising()
