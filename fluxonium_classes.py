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
        self.freq = np.sqrt(8*self.E_C_r*self.E_L_r) 
        self.fock_dim = fock_dim
        self.n_zpf = 1/2*(self.E_L_r/(2*self.E_C_r))**(1/4)

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
        eff_freq = np.sqrt(8*self.E_C_eff * self.E_L_r)
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
        self.E_C_effs = None

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
    
    def get_h_effs(self, phi_e):

        if self.E_C_effs is None:
            raise ValueError("E_C_eff not set")
        
        H_effs = []
        for E_C_eff in self.E_C_effs:

            h = self.get_h(E_C_eff, phi_e)
            H_effs.append(h)

        return H_effs

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
        res_coupling_strengths = []

        fluxonium_E_C_effs = []

        for res, cap_mat in zip(resonators, capacitance_matrices):
            inv_cap_mat = np.linalg.inv(cap_mat)
            E_Cr = inv_cap_mat[1,1]/2 * self.unit_fF_inv_to_GHz
            E_C = inv_cap_mat[0,0]/2 * self.unit_fF_inv_to_GHz
            g = 2*(inv_cap_mat[0,1] + inv_cap_mat[1,0]) * self.unit_fF_inv_to_GHz

            res.set_E_C_eff(E_Cr)
            fluxonium_E_C_effs.append(E_C)
            res_couplings.append(g * tensor(self.get_charge_op(), res.get_charge_op()))

            res_coupling_strengths.append(g)
        self.set_E_C_effs(fluxonium_E_C_effs)
        self.res_couplings = res_couplings
        self.res_coupling_strengths = res_coupling_strengths

    def find_chi_SW(self, eigenergies, eigstates, n_op, g, fr, nr_zpf, number_states):

        chi = 0

        for state_ind in range(number_states):

            chi = chi + get_chi_ij(0, state_ind, eigenergies, eigstates, n_op, g, fr, nr_zpf)-get_chi_ij(1, state_ind, eigenergies, eigstates, n_op, g, fr, nr_zpf)

        return chi
    
    def find_chi_over_flux_SW(self, res_ind, number_states, phi_es):

        g = self.res_coupling_strengths[res_ind]
        res = self.resonators[res_ind]
        fr = res.get_freq()
        nr_zpf = res.get_n_zpf()
        n_op = self.get_charge_op()

        chis = []

        for phi_ext in phi_es:
            h_effs = self.get_h_effs(phi_ext)
            h_eff = h_effs[res_ind]
            eigenergies, eigstates = h_eff.eigenstates()
    
            chi = self.find_chi_SW(eigenergies, eigstates, n_op, g, fr, nr_zpf, number_states)
            chis.append(chi)
        E_C_effs = self.E_C_effs
        plt.plot(phi_es, chis)
        plt.xlabel(r'$\phi_e$')
        plt.ylabel(r'$\chi/2\pi$ (MHz)')

        return chis

    def find_chi_full(self, res_ind, phi_es):
        if not hasattr(self, "res_couplings"):
            raise RuntimeError("Call add_resonator_coupling first")

        res = self.resonators[res_ind]
        a = destroy(res.fock_dim)
        a_dag = tensor(qeye(self.N_phi), a.dag())

        chis = []
        fqs = []
        frs = []
        for phi_ext in phi_es:
            h_effs = self.get_h_effs(phi_ext)
            h_eff = h_effs[res_ind]
            H_0 = tensor( h_eff, qeye(res.fock_dim) ) + tensor( qeye(self.N_phi), res.get_h_eff() )
            H_tot = H_0 + self.res_couplings[res_ind]
            energies, states = H_tot.eigenstates()

            ground_state = states[0]
            state_q1_r0 = states[1] # assumes qubit frequency is below resonator
            f_01 = energies[1]-energies[0]
            fqs.append(f_01)

            identified_state_ids = self.identify_state_index([(a_dag * ground_state).unit(), (a_dag * state_q1_r0).unit()], states)
            state_id_q0_r1 = identified_state_ids[0]
            state_id_q1_r1 = identified_state_ids[1]
            f01_plus_chi = energies[state_id_q1_r1]-energies[state_id_q0_r1]

            fr = energies[state_id_q0_r1] - energies[0]
            frs.append(fr)
            
            chi = (f01_plus_chi - f_01)*10**3
            chis.append(chi)
        plt.plot(phi_es, chis)
        plt.xlabel(r'$\phi_e$')
        plt.ylabel(r'$\chi/2\pi$ (MHz)')
        return chis, frs, fqs

    def get_charge_op(self):
        N_phi = self.N_phi
        delta_phi = self.nbr_periods*np.pi/(N_phi-1)
        n_up = Qobj(np.diag(np.ones(N_phi - 1), 1))
        n_down = Qobj(np.diag(-1*np.ones(N_phi - 1), -1))
        n_op = -1j/(2*delta_phi)*(n_up+n_down)
        return n_op
    
    def set_E_C_effs(self, E_C_effs):
        self.E_C_effs = E_C_effs

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



def get_g_ij(state_i, state_j, eigstates,  n_op, g, nr_zpf):

    g_ij = nr_zpf * eigstates[state_i].overlap(2*np.pi*(g*n_op) * eigstates[state_j]) # radians

    return g_ij

def get_w_ij(i, j, eigenergies):

    w_ij = 2*np.pi*(eigenergies[i] - eigenergies[j])

    return w_ij

def get_chi_ij(i, j, eigenergies, eigstates, n_op, g, fr, nr_zpf):

    w_ij = get_w_ij(i, j, eigenergies)
    g_ij = get_g_ij(i, j, eigstates, n_op, g, nr_zpf)

    w_res = 2*np.pi*fr
    chi_ij = np.abs(g_ij)**2*(1/(w_ij-w_res)+1/(w_ij+w_res))*10**3/(2*np.pi) # MHz

    return chi_ij



    



