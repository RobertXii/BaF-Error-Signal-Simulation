import numpy as np
import scipy
import warnings
from spin_params import *
from molecular_parameters import T00
from hamiltonian import H0_sigma,HZeeman_sigma,H0_pi_parity_basis,HZeeman_pi_parity_basis,H_int
from states import SigmaLevel, PiLevelParity, Superposition
#from numba import jit
from joblib import Parallel, delayed
import multiprocessing
import time

class System():
    """ What does the system class contain!"""
    def __init__(self,N_sigma,J_pi,B_field = [0,0,0],ignore_mF = False):
        """
        J_pi is a string or a list of string where each element is of the form 'J-' or 'J+' where 
        J is the total angular momentum and the sign following is the parity
        """ 
        if type(N_sigma) == list:
            self.N_list = N_sigma
        else:
            self.N_list = [N_sigma]
        
        if type(J_pi) == list:
            self.J_list = J_pi
        else:
            self.J_list = [J_pi]
        
        self.B_field = B_field
        self.sigma_states = []
        self.pi_states = []
        self.F_plus_sigma_all = []
        self.F_plus_pi_all = []
        self.generate_sigma_states(ignore_mF = ignore_mF)
        self.generate_pi_states(ignore_mF = ignore_mF)
        

        self.sigma_Hamiltonian = self.SigmaHamiltonian(self.sigma_states,self.F_plus_sigma_all,self.B_field)
        self.pi_Hamiltonian = self.PiHamiltonian(self.pi_states,self.F_plus_pi_all,self.B_field)
        
        self.interaction_Hamiltonian = None
        self.branching_ratios = None
        #self.amu = amu
        
       
        
    class SigmaHamiltonian():
        def __init__(self,states,F_plus,B_field):
            self.bare = []
            self.Zeeman = self.Zeeman(states,F_plus) ####3
            self.states = states
            self.diagonalized_states = []
            self.diagonalized_Hamiltonian = None
            self.B_field = B_field
            #self.species = species
        
        class Zeeman():
            def __init__(self,states,F_plus):
                self.Nlevels = len(states)
                self.states = states
                #print(self.states)
                self.F_plus,self.F_minus = self.create_ladder_operators(F_plus)
                self.F_x = 1/2*(self.F_plus+self.F_minus)
                self.F_y = -1j/2*(self.F_plus-self.F_minus)
                self.X,self.Y,self.Z =[],[],[]
                #self.generate_Zeeman()

            def create_ladder_operators(self,F_plus):
                A_plus = np.zeros((self.Nlevels,self.Nlevels))
                current_index = 0
                for i,sub_F_plus in enumerate(F_plus):
                    m = np.shape(sub_F_plus)[0] #sub_F_plus is a numpy array of dim 2F+1 x 2F+1
                    A_plus[current_index:current_index+m,current_index:current_index+m] = sub_F_plus
                    current_index += m
                A_minus = np.transpose(A_plus)
                return A_plus,A_minus

            def generate_Zeeman(self):
                #check if the states have been generated with the mF levels
                if self.states[0].mF == None:
                    raise ValueError("Error. Cannot generate Zeeman Hamiltonian. States generated without magnetic sublevels")
                HZ = np.zeros((self.Nlevels,self.Nlevels))
                for row in range(self.Nlevels):
                    for col in range(row+1):
                        HZ[row,col] = HZeeman_sigma(self.states[row],self.states[col])
                        if row != col:
                            HZ[col,row] = np.conjugate(HZ[row,col]) #########################conjugated
                self.Z = HZ
                
                #generate the X Hamiltonian
                # X is obtained by rotating the Z Hamiltonian by pi/2 about the y axis
                Ux = scipy.linalg.expm(-1j*self.F_y*np.pi/2)
                Uxd = scipy.linalg.expm(1j*self.F_y*np.pi/2)
                self.X = Ux@self.Z@Uxd
                #print(self.X)
                #same for Y
                #Y is obtained by rotating the Z Hamiltonian by -pi/2 about the x axis
                Uy = scipy.linalg.expm(1j*self.F_x*np.pi/2)
                Uyd = scipy.linalg.expm(-1j*self.F_x*np.pi/2)
                self.Y = Uy@self.Z@Uyd
                
                
        
        def generate_bare(self):
            num = len(self.states)
            print(num)
            H0 = np.zeros((num,num),dtype=np.complex_)
            for row in range(num):
                for col in range(row+1):
                    H0[row,col] = H0_sigma(self.states[row],self.states[col])
                    if row != col:
                        H0[col,row] = H0[row,col]
            self.bare = H0
            
        
            
        def diagonalize(self):
            if len(self.bare) == 0 and len(self.Zeeman.Z) == 0:
                raise ValueError("Error. Hamiltonian not generated.")
            if len(self.Zeeman.Z) == 0 and len(self.bare) != 0:
                H_temp = self.bare
            if len(self.Zeeman.Z) != 0 and len(self.bare) != 0:
                H_temp = self.bare + self.B_field[0]*self.Zeeman.X+ \
                                     self.B_field[1]*self.Zeeman.Y+ \
                                     self.B_field[2]*self.Zeeman.Z
            #print(H_temp)
            w,v = np.linalg.eigh(H_temp) #columns of v are the eigenvectors
            #print(w)
            sort_idx = np.argsort(w)
            #print(v)
            w,v = w[sort_idx],v[:,sort_idx]
            #print(v)
            self.diagonalized_Hamiltonian = np.conjugate(v.transpose())@H_temp@v
            w=np.round(w,4)
            v=np.round(v,4)
            
            #print(w)
            #to represent the diagonalized states
            for i in range(len(w)):
                v_temp = np.round(v[:,i],2)###############################################################
                #v_temp = v[:,i]
                non_zero_idx = np.nonzero(v_temp)[0]
                amp = []
                st = []
                for idx in non_zero_idx:
                    amp.append(v_temp[idx])
                    st.append(self.states[idx])
                self.diagonalized_states.append(Superposition(amp,st))
        
        """
        def delete_states(self,delete_list: list):
            nonlocal branching_ratios
            nonlocal interaction_Hamiltonian
            #delete_list is the index of all the states to be deleted
            for i in delete_list:
                #check if there is any branching from the state (in case of ground levels)
                # and any excitation to the state (in case of excited levels)
                branching_to_current = branching_ratios[i,:]
                excitation_to_current = interaction_Hamiltonian[i,:] #could chosse the columns too
                branching_to_current = True if (np.sum(branching_ratios[i,:])>1e-4) else False #branching ratio matrix has ground levels along the rows
                doesBranch = True in (x>1e-4 for x in branching_to_current)
                doesExcite = True in (np.abs(x)>1e-4 for x in excitation_to_current)
                if doesBranch:
                    warnings.warn('Non-zero branching to this level index :'+str(i))
                if doesExcite:
                    warnings.warn('Non-zero excitation from/to this level index :'+str(i))
                branching_ratios = np.delete(branching_ratios,i,axis = 0) #axis implying row
                self.diagonalized_Hamiltonian = np.delete(self.diagonalized_Hamiltonian,i,axis = 0)
                self.diagonalized_Hamiltonian = np.delete(self.diagonalized_Hamiltonian,i,axis = 1)
                self.diagonalized_states.pop(i)
         """
                
                
        
    
    
    
    def generate_sigma_states(self,ignore_mF = False):
        """Use ignore_mF = True to generatte states without mFs"""
        for N in self.N_list:
            for G in np.arange(np.abs(I1-S),I1+S+1):
                for F1 in np.arange(np.abs(N-G),N+G+1):
                    for F in np.arange(np.abs(F1-I2),F1+I2+1):
                        if ignore_mF:
                            self.sigma_states.append(SigmaLevel(G,N,F1,F))
                        else:
                            #construct the rotation matrix here for 
                            self.F_plus_sigma_all.append(self.create_F_plus(F))
                            for mF in np.arange(-F,F+1):
                                self.sigma_states.append(SigmaLevel(G,N,F1,F,mF))
    
    
    def create_F_plus(self,F):
        MF = np.arange(-F,F+1)
        F_plus = np.zeros((len(MF),len(MF)))
        for i,mF in enumerate(MF):
            try:
                F_plus[i+1,i] = np.sqrt(F*(F+1)-mF*(mF+1))
            except:
                continue
        return F_plus
    
    
    def generate_pi_states(self,ignore_mF = False):
        """Use ignore_mF = True to generatte states without mFs"""
        for J_str in self.J_list:
            if J_str[-1] == '+':
                parity = 1
            elif J_str[-1] == '-':
                parity = -1
            J = eval(J_str[:-1])
            for F1 in np.arange(np.abs(J-I1),J+I1+1):
                for F in np.arange(np.abs(F1-I2),F1+I2+1):
                    if ignore_mF:
                        self.pi_states.append(PiLevelParity(parity,J,F1,F))
                    else:
                        #construct the rotation matrix here for 
                        self.F_plus_pi_all.append(self.create_F_plus(F))
                        for mF in np.arange(-F,F+1):
                            self.pi_states.append(PiLevelParity(parity,J,F1,F,mF))
                        

        
    class PiHamiltonian():
        def __init__(self,states,F_plus,B_field):
            self.bare = []
            #print("F_plus", F_plus)
            self.Zeeman = self.Zeeman(states,F_plus)
            self.states = states
            self.diagonalized_states = []
            self.diagonalized_Hamiltonian = None
            self.B_field = B_field
        
        class Zeeman():
            def __init__(self,states,F_plus):
                self.Nlevels = len(states)
                self.states = states
                #print(self.states)
                self.F_plus,self.F_minus = self.create_ladder_operators(F_plus)
                self.F_x = 1/2*(self.F_plus+self.F_minus)
                self.F_y = -1j/2*(self.F_plus-self.F_minus)
                self.X,self.Y,self.Z =[],[],[]
                #self.generate_Zeeman()
                
                
            
            def create_ladder_operators(self,F_plus):
                A_plus = np.zeros((self.Nlevels,self.Nlevels))
                current_index = 0
                for i,sub_F_plus in enumerate(F_plus):
                    m = np.shape(sub_F_plus)[0] #sub_F_plus is a numpy array of dim 2F+1 x 2F+1
                    A_plus[current_index:current_index+m,current_index:current_index+m] = sub_F_plus
                    current_index += m
                A_minus = np.transpose(A_plus)
                return A_plus,A_minus
            
            
            def generate_Zeeman(self):
                #check if the states have been generated with the mF levels
                if self.states[0].mF == None:
                    raise ValueError("Error. Cannot generate Zeeman Hamiltonian. States generated without magnetic sublevels")
                
                #generate the Z hamiltonian first
                HZ = np.zeros((self.Nlevels,self.Nlevels))
                for row in range(self.Nlevels):
                    for col in range(row+1):
                        HZ[row,col] = HZeeman_pi_parity_basis(self.states[row],self.states[col])
                        if row != col:
                            HZ[col,row] = np.conjugate(HZ[row,col]) ##############################################conjugated!
                self.Z = HZ
                #generate the X Hamiltonian
                # X is obtained by rotating the Z Hamiltonian by pi/2 about the y axis
                Ux = scipy.linalg.expm(-1j*self.F_y*np.pi/2)
                self.X = Ux@self.Z@np.transpose(np.conjugate(Ux))
                
                #generate the Y Hamiltonian
                #Y is obtained by rotating the Z Hamiltonian by -pi/2 about the x axis
                Uy = scipy.linalg.expm(1j*self.F_x*np.pi/2)
                self.Y = Uy@self.Z@np.transpose(np.conjugate(Uy))
                
        
        def generate_bare(self):
            num = len(self.states)
            #check if the states have been generated with the mF levels
           
            H0 = np.zeros((num,num),dtype=np.complex_)
            for row in range(num):
                for col in range(row+1):
                    H0[row,col] = H0_pi_parity_basis(self.states[row],self.states[col])
                    if row != col:
                        H0[col,row] = H0[row,col]
            self.bare = H0  + T00*np.identity(np.shape(H0)[0])   
        
            
        def diagonalize(self):
            if len(self.bare) == 0 and len(self.Zeeman.Z) == 0:
                raise ValueError("Error. Hamiltonian not generated.")
            if len(self.Zeeman.Z) == 0 and len(self.bare) != 0:
                H_temp = self.bare
            if len(self.Zeeman.Z) != 0 and len(self.bare) != 0:
                H_temp = self.bare + self.B_field[0]*self.Zeeman.X+ \
                                     self.B_field[1]*self.Zeeman.Y+ \
                                     self.B_field[2]*self.Zeeman.Z
            
            w,v = np.linalg.eigh(H_temp) #columns of v are the eigenvectors
            sort_idx = np.argsort(w)
            w,v = w[sort_idx],v[:,sort_idx]
            self.diagonalized_Hamiltonian = np.conjugate(v.transpose())@H_temp@v
            w=np.round(w,4)
            v=np.round(v,4)
            
            
            #to represent the diagonalized states
            for i in range(len(w)):
                v_temp = np.round(v[:,i],2) #############################################################
                #v_temp = v[:,i]
                non_zero_idx = np.nonzero(v_temp)[0]
                amp = []
                st = []
                for idx in non_zero_idx:
                    amp.append(v_temp[idx])
                    st.append(self.states[idx])
                self.diagonalized_states.append(Superposition(amp,st))
    
    
    
    @staticmethod
    def generate_interaction_matrix(state1:list,state2:list,pol):
        num_1 = len(state1)
        num_2 = len(state2)
        Hmat = np.zeros((num_1,num_2),dtype=np.complex_)
        """
        for m in range(num_1):
            for n in range(num_2):
                Hmat[m,n] = H_int(state1[m],state2[n],pol)
        
        """
        Htemp = Parallel(n_jobs = -1)(delayed(H_int)(state1[m],state2[n],pol) for m in range(num_1) for n in range(num_2))
        #Htemp = Parallel(n_jobs = -1)(delayed(lambda m,n :(print(f"({m},{n})"),H_int(state1[m],state2[n],pol))[1])(m,n) for m in range(num_1) for n in range(num_2))
        #print(f"Pol : {pol}")
        for m in range(num_1):
            #print(m,end=",")
            for n in range(num_2):
                Hmat[m,n] = Htemp[0]
                Htemp.pop(0)
        
        return Hmat
    
    
    def generate_interaction_Hamiltonian(self,state1:list,state2:list,pol = 0):
        num_1 = len(state1)
        num_2 = len(state2)
        Hint = np.zeros((num_1+num_2,num_1+num_2),dtype=np.complex_)
        
        Htemp = Parallel(n_jobs = -1)(delayed(H_int)(state1[m],state2[n],pol) for m in range(num_1) for n in range(num_2))
        for m in range(num_1):
            for n in range(num_2):
                Hint[m,num_1+n] = Htemp[0]
                Hint[num_1+n,m] = np.conj(Hint[m,num_1+n])
                Htemp.pop(0)
        #Hint = Hint + np.conj(Hint.T)
        self.interaction_Hamiltonian = Hint
        
    def generate_branching_ratios(self,ground_state:list,excited_state:list):
        start = time.time()
        Trans_z = np.abs(self.generate_interaction_matrix(ground_state,excited_state,pol=0))**2
        stop_pi = time.time()
        print(f"Pi branching took : {stop_pi-start} sec")
        Trans_sigma_plus = np.abs(self.generate_interaction_matrix(ground_state,excited_state,pol=1))**2 #pol=-1 is sigma plus because of the convecntion. the matrix element represents transition from state2 to state1 whcih are normally the excited and ground states
        stop_sigmaminus = time.time()
        print(f"Sigma- branching took : {stop_sigmaminus-stop_pi} sec")
        Trans_sigma_minus = np.abs(self.generate_interaction_matrix(ground_state,excited_state,pol=-1))**2
        stop_sigmaplus = time.time()
        print(f"Sigma+ branching took : {stop_sigmaplus-stop_sigmaminus} sec")
        Trans_tot = Trans_z+Trans_sigma_plus+Trans_sigma_minus
        sum_over_ground = np.sum(Trans_tot,axis=0) #sum over the rows
        rows,cols = np.shape(Trans_tot)
        for i in range(cols):
            Trans_tot[:,i] /= sum_over_ground[i]
        self.branching_ratios = Trans_tot
        
 
