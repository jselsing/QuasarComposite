## Lisa Winter -- Feb. 19, 2010

def poly(y, c): # input y value and array with constant
    a = 0.0
    
    for i, item in enumerate(c):
        a += c[i]*y**i
    return a

def cardelli_reddening(wave, flux, ebv, R_V=3.1): #Compute the reddening curve value at this wavelength
    #### Assumed wave - angstroms and in the optical
    
    for i, mywave in enumerate(wave):
    
        mywave = mywave*0.0001 #convert angstrom wave to micron
    
        x = 1./mywave

        if (x >= 0.3) and (x < 1.1): #IR
            y = x
            a = 0.574*x**(1.61)
            b = -0.527*x**(1.61)

        elif (x >= 1.1) and (x < 3.3): #optical to NIR
            y = x - 1.82
        
            #a = 1 + (0.17699*y) - (0.50447*y**2) - (0.02427*y**3)
            #a = a + (0.72085*y**4) + (0.01979*y**5) - (0.77530*y**6) + (0.32999*y**7)
        
            #b = (1.41338*y) + (2.28305*y**2) + (1.07233*y**3)
            #b = b - (5.38434*y**4) - (0.62251*y**5) + (5.30260*y**6) - (2.09002*y**7)
            ## New values from O'Donnell 1994
            c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137, -1.718,   -0.827,    1.647, -0.505 ]
            c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985, 11.102,    5.491,  -10.805,  3.347 ]
    
            a = poly(y, c1)
            b = poly(y, c2)
        
        elif (x >= 3.3) and (x < 8): #mid-UV (1250A - 3333A)
            y = x
            F_a = 0
            F_b = 0
        
            if y > 5.9:
                y1 = y - 5.9
                F_a = -0.04473*y1**2 - 0.009779 * y1**3
                F_b = 0.2130 * y1**2  +  0.1207 * y1**3
        
            a = 1.752 - 0.316*y - (0.104 / ( (y-4.67)**2 + 0.341 )) + F_a
            b = -3.090 + 1.825*y + (1.206 / ( (y-4.62)**2 + 0.263 )) + F_b
    
        elif (x >= 8) and (x <= 11): #Far-UV (909A - 1250A)
            y = x-8.
        
            c1 = [ -1.073, -0.628,  0.137, -0.070 ]
            c2 = [ 13.670,  4.257, -0.420,  0.374 ]
        
            a = poly(y, c1)
            b = poly(y, c2)
    
        A_V = R_V * ebv
        A_lambda = A_V*(a + b/R_V)
    
        flux[i] = flux[i] * 10.**(0.4*A_lambda)
    
    return flux