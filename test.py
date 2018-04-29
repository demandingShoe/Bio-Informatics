
#To import python libraries

import scipy.integrate as spi
import numpy as np
import pylab as pl

#Initialising the variables (which are used in the equations) to fixed values which have already been studied 

beta=1.63e-7
gamma=float(raw_input("Enter gamma value: ")) #Will have to modify this value to improve our model


I0=846
S0=4292419
INP = (S0, I0, 0.0)

#Writing down the equations

def diff_eqs(INP,t):
    '''The main set of equations'''
    Y=np.zeros((3)) #Y=[0,0,0]
    V=INP
    Y[0] = -(beta * V[0] * V[1])
    Y[1] = (beta * V[0] * V[1]) - (gamma * V[1])
    Y[2] = gamma * V[1]
    return Y   # For odient in scipy to solve the differential equation



t_start = float(raw_input("Enter the starting time: "))
ND=float(raw_input("Enter the ending time: "))
TS=float(raw_input("Enter the step interval: "))
t_end = ND
t_inc = TS

#to build an array using numpy

t_range = np.arange(t_start, t_end+t_inc, t_inc)

#callable function,array(initial conditions),sequence of time points for which to solve for y

RES = spi.odeint(diff_eqs,INP,t_range)

print RES

#Ploting and labelling

pl.plot(RES[:,0], '-bs', label='Susceptibles') 
pl.plot(RES[:,2], '-g^', label='Recovereds')  
pl.plot(RES[:,1], '-ro', label='Infectious')
pl.legend(loc=0)
pl.title('SIR epidemic without births or deaths')
pl.xlabel('Time')
pl.ylabel('Susceptibles, Recovereds, and Infectious')
pl.savefig('2.1-SIR-high.png', dpi=900)
pl.show()

pl.figure()
pl.plot(RES[:,1], '-ro',label='Infectious')
pl.legend(loc=0)
pl.title("Only infectious ones")
pl.xlabel('Time')
pl.ylabel('Infectious')
pl.ylim([0,5000])
pl.show()
        
