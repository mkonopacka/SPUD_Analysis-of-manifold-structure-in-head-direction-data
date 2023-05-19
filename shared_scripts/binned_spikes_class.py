''' 
November 8th 2018
Interface to the binned spike counts/rate estimates. 
'''
# General imports
import numpy as np
import sys, os

# Import shared scripts
gen_fn_dir = os.path.abspath('.') + '/shared_scripts'
sys.path.append(gen_fn_dir)
from general_file_fns import load_pickle_file, save_pickle_file

# Load general parameters
gen_params = load_pickle_file('./general_params/general_params.pkl')

# Define spike_counts class
class spike_counts:
    '''Interface to the estimates of the rates/spike counts, allowing for 
    several operations that either pull out data from particular intervals 
    during a state or concatenate across intervals, and that filter data by
    brain region. Can use either binned spike counts or rate estimates. For historical reasons
    these are saved in similar structures but with slight differences.

    Args:
        session_id(str): id of used session
        params(dict): dict containing parameters like sigma, dt, binwidth
        count_type(str): 'rate' or 'binned'
        anat_region(str): 'all' or a specific area ("ADn", "PoS"...)

    Attributes:
        spike_data(dict): dict containing the spike counts/rate estimates loaded from gen_params['spike_counts_dir'] file (for all 62 cells)
        cell_ids(list): list of cell ids in the session filtered by anatomical region

    '''
    def __init__(self, session_id: str, params: dict, count_type: str = 'rate', anat_region: str = 'all'):
        # Create spike_data_path based on chosen count_type
        if count_type == 'binned':
            bw = params['binwidth']
            spike_data_path = f"{gen_params['spike_counts_dir']}{bw * 1000:.0f}ms/{session_id}_{bw * 1000:.0f}ms.pkl"
        
        elif count_type == 'rate':
            sigma = params['sigma']
            dt = params['dt']
            if not (dt * 100) % (0.05 * 100) == 0:
                ValueError('Please specify a dt that is a multiple of 50ms using create_general_params.py script.')
            else:
                spike_data_path = f"{gen_params['kernel_rates_dir']}{sigma * 1000:.0f}ms_sigma/{session_id}.pkl"
                self.dt = dt
        
        # Load spike_data
        print(f"Loading spike counts from {spike_data_path}...")
        self.spike_data = np.load(spike_data_path, allow_pickle=True)

        # Save parameters
        self.cell_ids = sorted(self.spike_data['cells'])
        self.session_id = self.spike_data['session']
        self.count_type = count_type
        self.anat_region = anat_region
        self.params = params

        # Filter by anatomical region if necessary
        if anat_region is not 'all':
            area_info = load_pickle_file('./general_params/area_shank_info.pkl')
            relevant_shanks = area_info[session_id][anat_region]
            relevant_cells = [id for id in self.cell_ids if id[0] in relevant_shanks] # TODO shouldn't it be id[1]? There are no cells with shank 10
            self.cell_ids = relevant_cells

    def single_intvl_dict(self, state, interval_num=0):
        '''Get spike counts/rates for a single interval as a dict of the form
        {cell_id : counts}. Also returns the appropriately smoothed measured 
        head angle for that interval, and the end-points of the interval 
        (since we're calling this by interval number rather than interval
        end points)
        '''
        interval = sorted(self.spike_data[state].keys())[interval_num]
        
        spike_counts = {}
        if self.count_type == 'binned':
            for cell in self.cell_ids:
                spike_counts[cell] = self.spike_data[state][interval]['num_spikes'][cell].copy()
            interval_angles = np.array(self.spike_data[state][interval]['avg_angles'])
        elif self.count_type == 'rate':
            jump = np.round(self.dt / 0.05).astype(int)
            for cell in self.cell_ids:
                spike_counts[cell] = self.spike_data[state][interval]['rates'][cell][::jump].copy()
            # Check types
            interval_angles = np.array(self.spike_data[state][interval]['angles'][::jump])
            
        return interval, spike_counts, interval_angles

    def all_intvl_dict(self, state):
        '''Concatenate data across all intervals for the given state (in temporal order)
        and return a dict of similar format to single_intvl_dict along with the 
        correspondingly concatenated angles.
        '''
        spike_counts = {cell : np.array([]) for cell in self.cell_ids}
        all_angles = np.array([])
        interval_list = sorted(self.spike_data[state].keys())
        
        for i, interval in enumerate(interval_list):
            # Get the appropriat single_intvl_dict
            tmp_interval, curr_interval_dict, curr_angles = self.single_intvl_dict(
                state, interval_num=i)

            for cell in self.cell_ids:
                spike_counts[cell] = np.append(spike_counts[cell], 
                    curr_interval_dict[cell])
            
            all_angles = np.append(all_angles, curr_angles)
        return spike_counts, all_angles

    def get_spike_matrix(self, state, interval_num = 'all') -> np.ndarray:
        if interval_num != 'all':
            # single interval
            interval_bounds, spike_counts, angles = self.single_intvl_dict(state, interval_num)
        else:
            # concatenate
            spike_counts, angles = self.all_intvl_dict(state)
        
        count_matrix = np.array([spike_counts[cell] for cell in self.cell_ids]).T
        
        if interval_num != 'all':
            return interval_bounds, count_matrix.copy(), angles
        else:
            return count_matrix.copy(), angles


