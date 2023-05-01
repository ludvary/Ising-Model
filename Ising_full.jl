# ISING MODEL SIMULATION USING MONTE CARLO. HERE WE PLOT <E per spin> vs <M per spin> against temperature.
using Statistics
using Printf
function main()

    mcs_steps::Int64, sim_steps::Int64, L::Int64 = 5000, 1000, 25
    T_init::Int64, T_final::Int64, T_step::Float64, J::Float64, M::Float64, E::Float64 = 1, 10, 0.001, 1, 0, 0

    # starting from a random lattice config
    LATTICE = rand([1, -1], L, L)
    
    # finding the number of loops over temperature
    T_iters = trunc(Int64, (T_final - T_init)/T_step)
    
    # initializing the arrays for temperature, Average Energy per spin and Average Magnetization per spin
    AVG_E = zeros(Float64, T_iters)
    AVG_M = zeros(Float64, T_iters)
    T_ARR = zeros(Float64, T_iters)

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


    # this applies the metropolis hastings algorithm to the lattice at every temperature so that it comes to equillibrium so that we can start collecting data
    function do_equillibrium(lattice, Temparature)
        E = find_energy(lattice)
        for t in 1:mcs_steps
            p, q = rand(1:L), rand(1:L) # choosing a random lattice point

            del_E = 2* J * (lattice[p, q] * (lattice[p, mod1(q+1, L)] + lattice[p, mod1(q-1, L)] + lattice[mod1(p+1, L), q] + lattice[mod1(p-1, L), q])) # finding the energy change before and after the randomly choosen spin is flipped.

            # now we employ Metropolis-Hastings
            if del_E < 0
                lattice[p, q] = -lattice[p, q] # spin flip is accepted
                E += del_E                     # instantaneous values of E and M are calculated
                M += 2*lattice[p, q]           # we could've just called find_energy() and find_magnetization() here but that would've been overhead
            
            else
                prob = exp(-del_E/Temparature)
                random = rand()
                
                if random < prob
                    lattice[p, q] = -lattice[p, q] # spin flip is accepted
                    E += del_E                     # instantaneous values of E and M are calculated
                    M += 2*lattice[p, q]
                end
            end
        end

       return lattice # we don't want to collect any data when this function is called, hence it only returns the lattice which has been equilibrated
    end


    # this applies Metropolis-Hastings at equillibrium and this time we collect E and M for each temperature
    function do_simlulation(lattice, T)

        E = find_energy(lattice) # finding the initial energy and magnetization of the equilibrated lattice
        M = find_magnetization(lattice)
        E_arr = zeros(Float64, sim_steps) # initializing the arrrays. New arrays will be created at every temperature
        M_arr = zeros(Float64, sim_steps)

        for t in 1:sim_steps
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
            E_arr[t] = E/(L*L) # storing the values of E and M inside arrays. the averages of these arrays will be calculated and stored each time do_simlulation() is called and then the arrays will be reintialized
            M_arr[t] = M/(L*L)
        end
        return E_arr, M_arr
    end



    function find_average(array)
        return mean(array) 
    end
    


    for i in 1:T_iters 
        T = T_init + i*T_step # updating T
        LATTICE = do_equillibrium(LATTICE, T) # bringing the lattice at equilibrium at T and storing the lattice config in LATTICE
        E_arr, M_arr = do_simlulation(LATTICE, T) # input that equilibrated lattice and collect data

        AVG_E[i] = find_average(E_arr) # store it here for a particular T
        AVG_M[i] = find_average(M_arr)
        T_ARR[i] = T
    end # repeat for all T

    # writing to that dat file 
    file = open("Ising_full.dat", "w")
    for j in 1:T_iters
        @printf(file, "%.15f\t%.15f\t%.15f\n", T_ARR[j], AVG_E[j], AVG_M[j])
    end

end
main()