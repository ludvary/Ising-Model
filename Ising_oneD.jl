using Statistics
using Printf

function main()

    L :: Int16 = 1000
    J :: Float64 = 1
    eqm_steps :: Int64 = 10000
    sim_steps :: Int64 = 5000
    T_init :: Float64 = 1
    T_final :: Float64 = 10
    T_step :: Float64 = 0.001
    no_of_T :: Int64 = (T_final - T_init)/T_step

    Var_E = zeros(Float64, no_of_T)
    Var_M = zeros(Float64, no_of_T)
    Temperature_arr = zeros(Float64, no_of_T)


    LATTICE = rand([1, -1], L)

    
    function find_energy(lattice)
        E = 0
        for i in 1:L
            E -= J * (lattice[i]) * (lattice[mod1(i+1, L)] + lattice[mod1(i-1, L)])
        end
        return E
    end


    function find_magnetization(lattice)
        return sum(lattice)
    end

    function find_var(array)
        return var(array)
    end
    
    function do_eqm(lattice, temperature, J_ising)

        E = find_energy(lattice)
        M = find_magnetization(lattice)

        for k in 1:eqm_steps
            p = rand(1:L)

            del_E = 2 * J_ising * (lattice[p]) * (lattice[mod1(p+1, L)] + lattice[mod1(p-1, L)])

            if del_E < 0
                lattice[p] = -lattice[p]
                E += del_E
                M += 2 * lattice[p]

            else
                prob = exp(-del_E/temperature)
                rand_no = rand() 
                if rand_no < prob lattice[p] = -lattice[p]
                    E += del_E
                    M += 2 * lattice[p]
                end
            end
        end
        return lattice
    end


    function do_sim(lattice, temperature, J_ising)

        E = find_energy(lattice)
        M = find_magnetization(lattice)
        E_arr = zeros(Float64, L)
        M_arr = zeros(Float64, L)

        for k in 1:eqm_steps
            p = rand(1:L)

            del_E = 2 * J_ising * (lattice[p]) * (lattice[mod1(p+1, L)] + lattice[mod1(p-1, L)])

            if del_E < 0
                lattice[p] = -lattice[p]
                E += del_E
                M += 2 * lattice[p]

            else
                prob = exp(-del_E/temperature)
                rand_no = rand()
                
                if rand_no < prob
                    lattice[p] = -lattice[p]
                    E += del_E
                    M += 2 * lattice[p]
                end
            end
            E_arr[k] = E/L
            M_arr[k] = M/L
            
            return E_arr, M_arr
        end
    end 

    for i in 1:no_of_T
        T = T_init + i*T_step
        LATTICE = do_eqm(LATTICE, T, J)
        E_arr, M_arr = do_sim(LATTICE, T, J)
        Var_E[i] = find_var(E_arr)
        Var_M[i] = find_var(M_arr)
        Temperature_arr[i] = T
    end

    file = open("oneD.dat", "w")
    for i in 1:no_of_T
        @printf(file, "%.15f\t%.15f\t%.15f\n", Temperature_arr[i], Var_E[i], Var_M[i])
    end
end
main()
