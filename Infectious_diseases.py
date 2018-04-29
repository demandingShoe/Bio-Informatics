# Module for simulating infectious disease dynamics with the SIR (Susceptible-Infected-Recovered) model

import numpy as np, scipy.integrate, scipy.optimize, sys

class SIRsystem:

    def __init__(self, beta, gamma, S, I, R):
        
        self.beta=beta
        self.gamma=gamma
        self.S=S
        self.I=I
        self.R=R
        self.t=0
        self.N=S+I+R
        self.trajectory=np.array([[self.S,self.I,self.R]],dtype=float)/self.N
        self.times=None
        

    def reset(self, S, I, R, t=0.):
       
        self.t=t
        self.S=S
        self.I=I
        self.R=R
        self.trajectory=np.array([[self.S,self.I,self.R]],dtype=float)/self.N
        
        
        
    def get_total_infected(self):
       return self.I + self.R

class DeterministicSIRsystem (SIRsystem):

   
    def dydt(self, y, t):
        
        s,i,r=y
        
        dsdt=-self.beta * self.S * self.I
        didt=(self.beta * self.S * self.I) - (self.gamma * self.I)
        drdt= self.gamma * self.I

        return np.array([dsdt,didt,drdt])
        
        

    def run(self, T, dt=None):
        
        y0=np.array([[self.S,self.I,self.R]],dtype=float)/self.N

        if dt is None:
            self.times=np.linspace(0.,T,int(T+1),endpoint=True)
        else:
            self.times=np.arange(0,T,dt)
        self.trajectory=scipy.integrate.odeint(self.dydt,y0,self.times)

    
def SimulateDeterministicOutbreakSize(N, R0_range=np.arange(0.,5.,0.1)):
    
    gamma = 1.0
    beta = 0.0
    Nf = float(N)
    dsir = DeterministicSIRsystem(beta, gamma, (N-1), 1, 0)
    sizes = []
    for R0 in R0_range:
        beta = R0 * gamma
        dsir.beta = beta
        dsir.reset((N-1),1,0,0.)
        dsir.run(100.)
        dsir.S, dsir.I, dsir.R = list(map(int, N*dsir.trajectory[-1]))
        R_inf = dsir.get_total_infected()
        sizes.append((R0, R_inf))
    return np.array(sizes)

def CalculateDeterministicOutbreakSize(N, R0_range=np.arange(0.,5.,0.1)):
    
    func = lambda R_inf, R0: R_inf - (1.-np.exp(-R0*R_inf))
    sizes = []
    for R0 in R0_range:
        R_inf = scipy.optimize.fsolve(func, 0.5, args=(R0,))[0]
        sizes.append((R0, R_inf))
    return np.array(sizes)


def FractionLargeOutbreaks(osd, N):
   
    Nthresh = 0.1*N
    return (1.*np.sum(osd<Nthresh))/len(osd)

def yesno():
    response = input('    Continue? (y/n) ')
    if len(response)==0:        # [CR] returns true
        return True
    elif response[0] == 'n' or response[0] == 'N':
        return False
    else:                       # Default
        return True
    
def demo():
    import pylab
    N = 1000
    print("SIR demo")
    print("Deterministic SIR dynamics")
    pylab.figure(1)
    pylab.clf()
    dsir = DeterministicSIRsystem(1.5, 1.0, N-1, 1, 0)
    #dsir.run(30,0.1)
    pylab.plot(dsir.times, dsir.trajectory[:,0], 'b-', label='S')
    pylab.plot(dsir.times, dsir.trajectory[:,1], 'r-', label='I')
    pylab.plot(dsir.times, dsir.trajectory[:,2], 'g-', label='R')
    pylab.legend(loc='upper right')
    if not yesno(): return
    print("Deterministic outbreak size")
    R0_range = np.arange(0.,5.,0.1)
    simulated_sizes = SimulateDeterministicOutbreakSize(N=N, R0_range=R0_range)
    theoretical_sizes = CalculateDeterministicOutbreakSize(N=N, R0_range=R0_range)
    pylab.figure(2)
    pylab.clf()
    pylab.plot(R0_range, N*theoretical_sizes[:,1], 'b-', label='theory')
    pylab.plot(R0_range, simulated_sizes[:,1], 'bo', label='simulations')
    pylab.xlabel('R0')
    pylab.ylabel('total outbreak size')
    if not yesno(): return
    print("Stochastic SIR dynamics")
    pylab.figure(3)
    pylab.clf()
    pylab.plot(dsir.times, N*dsir.trajectory[:,1], 'r-', linewidth=2)
    for n in range(20):
        tfinal = ssir.run(20)
        pylab.plot(ssir.times, ssir.trajectory[:,1], 'b-')
    if not yesno(): return


demo()
