function reference_function(t,γ)
    #Reference_function computes the analytical solution to the 2nd order ODE
    #for initial conditions y(0) = 1 and y'(0) = 0
    #γ is also used as an input, as the reference function depends on whether the
    #system is critically damped. γ can be thought of as the normalized damping 
    #coefficient

    #NOTE: See problems_03_solution.pdf for derivation

    if(abs(γ - 2) < 1e-7) #True for (near) critically damped cases
        solution =  @.[(1 + t)*exp(-t) -t*exp(-t)]
    else 
        lambda_1 = (-γ + sqrt(γ^2-(4+0im))) / 2
        lambda_2 = (-γ - sqrt(γ^2-(4+0im))) / 2
        a = - lambda_2 / (lambda_1 - lambda_2)
        b = 1 - a
        solution = @.[a*exp(lambda_1*t) + b*exp(lambda_2*t);]    #y'(t)
    end
    return solution
end

@info "Reference function is now defined"