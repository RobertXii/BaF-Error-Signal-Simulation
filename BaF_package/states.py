import numpy as np
from spin_params import S,LAMBDA,SIGMA,OMEGA
class SigmaLevel:
    def __init__(self,G,N,F1,F,mF=None):
        self.G = G
        self.N = N
        self.F1 = F1
        self.F = F
        self.mF = mF
        
    def __repr__(self):
        return f"|G = {self.G}, N = {self.N}, F1 = {self.F1}, F = {self.F}, mF = {self.mF}>"
    def __str__(self):
        return f"|G = {self.G}, N = {self.N}, F1 = {self.F1}, F = {self.F}, mF = {self.mF}>"
        
    """
    def GtoJ(self):
        
    """


class PiLevelParity:
    def __init__(self,parity,J,F1,F,mF=None):
        self.J = J
        self.F1 = F1
        self.F = F
        self.mF = mF
        self.parity = parity
          
        
    def __str__(self):
        if self.parity == 1:
            return f"|J = {self.J}+, F1 = {self.F1}, F = {self.F}, mF = {self.mF}>"
        elif self.parity == -1:
            return f"|J = {self.J}-, F1 = {self.F1}, F = {self.F}, mF = {self.mF}>"
        
    def __repr__(self):
        if self.parity == 1:
            return f"|J = {self.J}+, F1 = {self.F1}, F = {self.F}, mF = {self.mF}>"
        elif self.parity == -1:
            return f"|J = {self.J}-, F1 = {self.F1}, F = {self.F}, mF = {self.mF}>"
        
    def parity_to_omega(self):
        state1 = PiLevelOmega(LAMBDA,SIGMA,OMEGA,self)
        state2 = PiLevelOmega(-LAMBDA,-SIGMA,-OMEGA,self)
        return Superposition([1/np.sqrt(2),1/np.sqrt(2)*self.parity*(-1)**(self.J-S)],[state1,state2]) 
        
        

class PiLevelOmega:
    def __init__(self,Lambda,Sigma,Omega,parity_state:PiLevelParity):
        self.Lambda = Lambda
        self.Sigma = Sigma
        self.Omega = Omega
        self.parity_state = parity_state
    
    def __repr__(self):
        return f"|LAMBDA = {self.Lambda}, SIGMA = {self.Sigma}, OMEGA = {self.Omega} "+ \
                f"; J = {self.parity_state.J}, F1 = {self.parity_state.F1},F = {self.parity_state.F}, mF = {self.parity_state.mF}>"
    def __str__(self):
        return f"|LAMBDA = {self.Lambda}, SIGMA = {self.Sigma}, OMEGA = {self.Omega} "+ \
                f"; J = {self.parity_state.J}, F1 = {self.parity_state.F1},F = {self.parity_state.F}, mF = {self.parity_state.mF}>"
        
        
    
        
class Superposition:
    """Coefficients is a list of coeffieicients and states is a list of the SigmaLevel or PiLevel objects"""
    def __init__(self,amplitude:list,states:list):
        self.amplitude = amplitude
        self.states = states
    
    def __str__(self):
        val=""
        for i,amp_val in enumerate(self.amplitude):
            if i == len(self.amplitude)-1 and np.abs(amp_val)>1e-4:
                val += str(np.round(amp_val,4)) + ' ' + str(self.states[i])
            else:
                val += str(np.round(amp_val,4)) + ' ' + str(self.states[i]) + ' + \n'
        return val
    
    def __repr__(self):
        val=""
        for i,amp_val in enumerate(self.amplitude):
            if i == len(self.amplitude)-1 and np.abs(amp_val)>1e-4:
                val += str(np.round(amp_val,4)) + ' ' + str(self.states[i])
            else:
                val += str(np.round(amp_val,4)) + ' ' + str(self.states[i]) + ' + \n'
        return val
            
            
            
        
        
