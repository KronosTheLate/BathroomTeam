function euler_integration(system_matrix, initial_values, t_start, t_end, step_size )
    #This function does numerical integration using Euler's method in order to
    #solve system of Ordinary Differential Equations (ODEs)

    #Arguments:
    # system_matrix (A):    a matrix describing the system of eqs. , i.e. A*x=dx/dt
    # initial_values:       initial conditions for the system
    # t_start:              start time of integration
    # t_end:                end time of integration
    # step_size:            the step size of integration (larger -> faster but error)

    #Returns:
    # x_solved:             A matrix consisting of solved x values for different values in time
    # t_vector:             time values used to evaluate the ODE(s)


    #Construct a time vector based upon t_start, t_end and step_size
    N_steps = floor(Int,(t_end - t_start)/step_size)
    t_vector = LinRange(t_start,t_end,N_steps)

    #Create a matrix to hold the solved values of the ODE(s)
    #Note the that the no. of rows correspond to the no. ODEs
    x_solved = zeros(size(initial_values,1),size(t_vector,1))

    #Set the initial values as the first values of x and transfer them to x_solved
    x_vect = initial_values
    x_solved[:,1] = x_vect

    #for loop to implement Euler's method
    for i in 2:N_steps
        dx = system_matrix*x_vect
        x_vect += dx * (t_vector[i] - t_vector[i-1])
        x_solved[:,i] = x_vect

    end    

    return_matrix = [x_solved; transpose(collect(t_vector))]
    
    return return_matrix 

    

end

@info "Euler is now defined"