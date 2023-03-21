function kepler_eom(init_cond,R)
    #init_cond should have the following shape:

    # row 1:3           = position vector
    # row 4:6           = momentum vector
    # No. of colums (N) = Number of particles

    #get particle count N
    N = size(init_cond,2)

    #Momentum and position vectors
    drv = copy(init_cond)

    

    for i in 1:N
        r           = init_cond[1:3,i]
        p           = init_cond[4:6,i]
        #global dr  = drv[1:3,i]
        #global dp  = drv[4:6,i]


        drv[4:6,i] = -norm(r-R)^2 .*(r-R)   #Momentum
        drv[1:3,i] = p                      #Position
       
    end

    return drv
end