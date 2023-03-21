function sys_func(x_vect)
    #system function describing the kepler system as given 
    #in problem set 2

    #Structure of x_vect: [r_x r_y p_x p_y]
    drv = [0; 0; 0; 0;]

    drv[1] = x_vect[3]
    drv[2] = x_vect[4]
    drv[3] = -x_vect[1]/norm(x_vect[1:2])^3
    drv[4] = -x_vect[2]/norm(x_vect[1:2])^3

    return drv
end

@info "system function defined"