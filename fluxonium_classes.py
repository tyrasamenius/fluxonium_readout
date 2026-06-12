import numpy as np
from qutip import *
import matplotlib.pyplot as plt

class Device:
    e = 1.602*10**(-19)
    h = 6.626*10**(-34)
    hbar = h/(2*np.pi)
    phi_0 = h/(2*e)
    unit_fF_inv_to_GHz = e**2 /h * 10**6
    unit_nH_inv_to_GHz = 1/h * (phi_0/(2*np.pi))**2

    def get_h(self, *args, **kwargs): ...

    def plot_transitions(self, *args, **kwargs): ...

    

class Resonator(Device):
    def __init__(self, L_r, C_r, fock_dim):
        self.L_r = L_r 
        self.C_r = C_r
        self.E_C_r = 1/(2*C_r) * self.unit_fF_inv_to_GHz
        self.E_L_r = 1/L_r  * self.unit_nH_inv_to_GHz
        self.freq = 1/np.sqrt(L_r*C_r)/(2*np.pi)
        self.fock_dim = fock_dim
        self.n_zpf = np.sqrt(self.hbar*self.C_r/(2*self.L_r))

        self.E_C_eff = None
        self.eff_freq = None
    
    def get_h(self, freq):
        a = destroy(self.fock_dim)
        H = freq*(a.dag()*a)
        return H
    
    def get_h_bare(self):
        return self.get_h(self.freq)
    
    def get_h_eff(self):
        if self.E_C_eff is None:
            raise ValueError("E_C_eff not set")
        eff_freq = np.sqrt(self.E_C_eff * self.E_L_r)
        return self.get_h(eff_freq)
    
    def get_charge_op(self):
        a = destroy(self.fock_dim)
        n_op = 1j * self.n_zpf * (a.dag() + a)
        return n_op
    
    def get_freq(self):
        return self.freq

    def get_n_zpf(self): 
        return self.n_zpf
    
    def set_E_C_eff(self, E_C_eff):
        self.E_C_eff = E_C_eff

class Fluxonium(Device):

    def __init__(self, E_J, #GHz 
                 C_J, # fF
                   C, # fF
                   E_L, # GHz
                   resonators = None, 
                   N_phi = 301, 
                   nbr_periods = 6):
        self.N_phi = N_phi
        self.nbr_periods = nbr_periods
        self.phi_vals = np.linspace(-nbr_periods*np.pi, nbr_periods*np.pi, N_phi)
        self.E_J = E_J
        self.C_J = C_J
        self.C = C
        self.E_C = 1/(2*(C_J+C)) * self.unit_fF_inv_to_GHz
        self.E_L = E_L
        self.resonators = resonators
        self.E_C_eff = None

    def get_h(self, E_C, phi_e):
        N_phi = self.N_phi
        phi_vals = self.phi_vals
        nbr_periods = self.nbr_periods

        n_up = Qobj(np.diag(np.ones(N_phi - 1), 1))
        n_down = Qobj(np.diag(np.ones(N_phi - 1), -1))
        n_diag = Qobj(np.diag(-2*np.ones(N_phi)))
        delta_phi = nbr_periods*2*np.pi/(N_phi-1)
        d_dphase_op = (n_up+n_down+n_diag)/delta_phi**2
        phase_op = Qobj(np.diag(phi_vals))
        
        Cterm = -4*E_C*d_dphase_op
        Jterm = -self.E_J*phase_op.cosm()
        Lterm = 1/2*self.E_L*(phase_op-phi_e)**2
        H = Cterm + Jterm + Lterm
        return H
    
    def get_h_bare(self, phi_e):
        return self.get_h(self.E_C, phi_e)
    
    def get_h_eff(self, phi_e):
        if self.E_C_eff is None:
            raise ValueError("E_C_eff not set")
        return self.get_h(self.E_C_eff, phi_e)

    def get_charge_op(self): ...

    def plot_wavefuncs(self, nbr_wavfuncs, phi_e):
        H = self.get_h_bare(phi_e)
        eigenenergies, eigenstates = H.eigenstates()
        fig, axes = plt.subplots(1, nbr_wavfuncs, figsize=(12, 3))

        for wavfunc_idx in range(nbr_wavfuncs):
            ax = axes[wavfunc_idx]
            psi = eigenstates[wavfunc_idx].full().flatten()
            ax.plot(self.phi_vals, psi.real)
            ax.set_xlabel('$\phi$')
            ax.set_title(fr"Phase distribution, $|{wavfunc_idx}\rangle$")

        plt.tight_layout()
        plt.show()
    
    def plot_transistions_over_flux(self, phi_es, N_levels):
        N_flux = len(phi_es)
        energy_over_flux_mat = np.empty((N_flux, N_levels))

        for phi_e_idx, phi_e in enumerate(phi_es):
            H = self.get_h_bare(phi_e)
            eigenenergies, eigenstates = H.eigenstates()
            for i in range(1, N_levels):
                energy_over_flux_mat[phi_e_idx, i] = eigenenergies[i] - eigenenergies[0]


        for i in range(1, N_levels):
            j = N_levels - i
            plt.plot(phi_es, energy_over_flux_mat[:, j], label = rf'E_{j}')
        plt.xlabel('$\Phi_e/\Phi_0$')
        plt.ylabel('$E-E_0$ (GHz)')
        plt.legend()
    
    def add_resonator_coupling(self, resonators, capacitance_matrices): 
        
        if len(resonators) != len(capacitance_matrices):
            print('Number of resonators and couplings need to match')
            return 
        self.resonators = resonators
        res_couplings = []
        for res, cap_mat in zip(resonators, capacitance_matrices):
            inv_cap_mat = np.linalg.inv(cap_mat)
            E_Cr = inv_cap_mat[1,1]/2 * self.unit_fF_inv_to_GHz
            E_C = inv_cap_mat[0,0]/2 * self.unit_fF_inv_to_GHz
            res.set_E_C_eff(E_Cr)
            self.set_E_C_eff(E_C)
            g = 2*(inv_cap_mat[0,1] + inv_cap_mat[1,0]) * self.unit_fF_inv_to_GHz

            res_couplings.append(g * tensor(self.get_charge_op(), res.get_charge_op()))
        self.res_couplings = res_couplings

    def find_chi_SW(): ...

    def find_chi_full(self, res_ind, phi_es):
        if not hasattr(self, "res_couplings"):
            raise RuntimeError("Call add_resonator_coupling first")

        res = self.resonators[res_ind]
        a = destroy(res.fock_dim)
        a_dag = tensor(qeye(self.N_phi), a.dag())

        chis = []

        for phi_ext in phi_es:
            H_0 = tensor( self.get_h_eff(phi_e = phi_ext), res.get_h_eff())
            H_tot = H_0 + self.res_couplings[res_ind]
            energies, states = H_tot.eigenstates()

            ground_state = states[0]
            state_q1_r0 = states[1] # assumes qubit frequency is below resonator
            f_01 = energies[1]-energies[0]

            identified_state_ids = self.identify_state_index([a_dag*ground_state, a_dag*state_q1_r0], states)
            state_id_q0_r1 = identified_state_ids[0]
            state_id_q1_r1 = identified_state_ids[1]
            f01_plus_chi = energies[state_id_q1_r1]-energies[state_id_q0_r1]
            
            chi = f01_plus_chi- f_01
            chis.append(chi)
        plt.plot(phi_es, chis)
        plt.xlabel('\phi_e')
        plt.ylabel('\chi/2\pi')
        return chis

    def get_charge_op(self):
        N_phi = self.N_phi
        delta_phi = self.nbr_periods*np.pi/(N_phi-1)
        n_up = Qobj(np.diag(np.ones(N_phi - 1), 1))
        n_down = Qobj(np.diag(-1*np.ones(N_phi - 1), -1))
        n_op = -1j/(2*delta_phi)*(n_up+n_down)
        return n_op
    
    def set_E_C_eff(self, E_C_eff):
        self.E_C_eff = E_C_eff

    def identify_state_index(self, states_to_identify, eigenstates):

        indices = []
        overlaps = np.empty((len(states_to_identify), len(eigenstates)))

        for unknown_id in range(len(states_to_identify)):

            for state_ind in range(len(eigenstates)):
                eigstate = eigenstates[state_ind]
                overlap = np.abs(eigstate.overlap(states_to_identify[unknown_id])) 
                overlaps[unknown_id, state_ind] = overlap
            eigstate_ind = np.argmax(overlaps[unknown_id,:])
            indices.append(eigstate_ind)
        return indices
    
class Double_Junc_Fluxonium(Device):
    def __init__(self, E_J, E_C, E_L, resonators = None, N_phi = 301, nbr_periods = 6):
        self.N_phi = N_phi
        self.nbr_periods = nbr_periods
        self.phi_vals = np.linspace(-nbr_periods*np.pi, nbr_periods*np.pi, N_phi)
        self.E_J = E_J
        self.E_C = E_C
        self.E_L = E_L
        self.resonators = resonators
    



