function leapfrog_integration(system_func, R, initial_values, t_start, t_end, step_size )

    #Construct a time vector based upon t_start, t_end and step_size
    N_steps = floor(Int,(t_end - t_start)/step_size)
    t_vector = LinRange(t_start,t_end,N_steps)

    #Create a matrix to hold the solved values of the ODE(s)
    #Note the that the no. of rows correspond to the no. ODEs
    x_solved = zeros(size(initial_values,1),size(initial_values,2),size(t_vector,1))

    #Set the initial values as the first values of x and transfer them to x_solved
    x_vect = copy(initial_values)
    x_solved[:,:,1] = x_vect

    #for loop to implement Euler's method
    for i in 2:N_steps
        dx = system_func(x_vect,R)
        dt = (t_vector[i] - t_vector[i-1])
        x_vect += dx*dt
        x_solved[:,:,i] = x_vect

    end    

    return_matrix = x_solved
    
    return return_matrix 

    

end

@info "Leapfrog is now defined"
