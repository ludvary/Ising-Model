using Printf

function main()

    L :: Int16 = 1000
    J :: Float64 = 1
    eqm_steps :: Int64 = 1000000
    T :: Float64 = 1

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

    
    function do_eqm(lattice, temperature, J_ising)

        E = find_energy(lattice)
        M = find_magnetization(lattice)
        E_arr = zeros(Float64, eqm_steps)
        M_arr = zeros(Float64, eqm_steps)
        t_arr = zeros(Int64, eqm_steps)

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
            E_arr[k] = E/L
            M_arr[k] = M/L
            t_arr[k] = k
        end
        return E_arr, M_arr, t_arr
    end

    E_arr, M_arr, t_arr = do_eqm(LATTICE, T, J)

    file = open("oneD_eqm.dat", "w")
    for i in 1:eqm_steps
        @printf(file, "%.15f\t%.15f\t%.15f\n", t_arr[i], E_arr[i], M_arr[i])
    end
end
main()
