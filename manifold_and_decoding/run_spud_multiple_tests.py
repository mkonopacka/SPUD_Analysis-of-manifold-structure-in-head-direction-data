'''March 21st 2019
Fit the wake manifold and do decoding. This version runs multiple cross-validated
tries and plots the distribution of decoding errors.
'''
# General imports
import numpy as np
import numpy.linalg as la
import sys, os 
import time, datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Import shared modules
gen_fn_dir = os.path.abspath('.') + '/shared_scripts'
sys.path.append(gen_fn_dir)
import general_file_fns as gff
from binned_spikes_class import spike_counts
from dim_red_fns import run_dim_red
import manifold_fit_and_decode_fns as mff

# Set random seed
sd = int((time.time()%1)*(2**31))
np.random.seed(sd)
print(f"Running run_spud_multiple_tests.py with random seed {sd} -----------------------------------")

# Load and print general params and create directory to load dimensionality reduction results if needed
gen_params = gff.load_pickle_file('./general_params/general_params.pkl')
print(f"General params used for this session:\n{gen_params}")

# Get current date and create directory to save results
curr_date = datetime.datetime.now().strftime('%Y_%m_%d')+'_'
dir_to_save = gff.return_dir(gen_params['results_dir'] + curr_date + '_curve_fits/')

# Set up command line or default parameters
cmd_line = False
run_dim_red_here = True # run dimensionality reduction, if False loads from `dim_red_dir`
dim_red_dir = gen_params['results_dir'] + 'dim_red/'
plot_projected_points = True
plot_path = f"./figures/"

if cmd_line:
    session = sys.argv[1]
    fit_dim = int(sys.argv[2])
    nKnots = int(sys.argv[3])
    knot_order = sys.argv[4]
    penalty_type = sys.argv[5]
    nTests = int(sys.argv[6])
    train_frac = float(sys.argv[7])
# TODO describe the parameters
else:
    session = 'Mouse28-140313'
    fit_dim = 3
    nKnots = 15
    knot_order = 'wt_per_len'
    penalty_type = 'mult_len'
    nTests = 10
    train_frac = 0.8

area = 'ADn'    # Anatomical region to use
state = 'Wake'  # Sleep stage to use
dt_kernel = 0.1 # Kernel time step(?)
sigma = 0.1     # Kernel width
method = 'iso'  # Dimensionality reduction method, `iso` stands for isomap
n_neighbors = 5 # Number of neighbors for isomap

# Print parameters
print(('Session: %s, fit dim: %d, nKnots: %d, knot_order: %s, penalty: %s, nTests: %d, train_frac: %.2f'%(
    session, fit_dim, nKnots, knot_order, penalty_type, nTests, train_frac)))

# Run dimensionality reduction if needed
if run_dim_red_here:
    print('Running initial dimensionality reduction...')
    rate_params = {'dt' : dt_kernel, 'sigma' : sigma}
    dim_red_params = {'n_neighbors' : n_neighbors, 'target_dim' : fit_dim}
    desired_nSamples = 15000
    # Get counts from spike_counts interface
    session_rates = spike_counts(session, rate_params, count_type='rate', anat_region=area) # for ADn there are 21 cell ids
    counts, tmp_angles = session_rates.get_spike_matrix(state)
    selected_counts = counts[:desired_nSamples]
    # Run dimensionality reduction and save projected points for later
    projection = run_dim_red(selected_counts, params = dim_red_params, method = method)
    np.save(dim_red_dir + f"{session}_{state}_{method}_{dim_red_params['target_dim']}_projection_{curr_date}.pkl", projection)
    assert projection.shape == (desired_nSamples, fit_dim)
    # If plot_projected_points, plot projection in 3d and save to file
    if plot_projected_points and fit_dim == 3:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        sc = ax.scatter(projection[:, 0], projection[:, 1], projection[:, 2], alpha=0.5)
        ax.set_title('Projection of %s %s data onto first two isomap dimensions'%(session, state))
        plt.colorbar(sc)
        plot_filename = plot_path + f"{session}_{state}_{method}_{dim_red_params['target_dim']}_projection_{curr_date}.png"
        plt.savefig(plot_filename)

    # Create embeddings dictionary
    embeddings = {state : projection, 'meas_angles' : tmp_angles[:desired_nSamples]}
    embeddings_fname = 'not_saved'
# If not running dimensionality reduction, load embeddings from file with specified pattern
else:
    file_pattern = '%s_%s_kern_%dms_sigma_%dms_binsep_%s_embeddings_%s_%ddims_%dneighbors_*.pkl'%(
        session, area, sigma * 1000, dt_kernel * 1000, state, method, fit_dim, n_neighbors)
    embeddings, embeddings_fname = gff.load_file_from_pattern(dim_red_dir+file_pattern)

curr_mani = embeddings[state]
nPoints = len(curr_mani)
nTrain = np.round(train_frac * nPoints).astype(int)

# Use measured angles to set origin and direction of coordinate increase
ref_angles = embeddings['meas_angles']
fit_params = {'dalpha' : 0.005, 'knot_order' : knot_order, 'penalty_type' : penalty_type, 'nKnots' : nKnots}

results = {}
tic = time.time()
k = (session, fit_dim, nKnots, knot_order, penalty_type, train_frac)
print('Fitting manifold')
for curr_sample in range(nTests):
    curr_fit_params = dict(fit_params)
        
    train_idx = np.random.choice(nPoints, size=nTrain, replace=False)
    test_idx = np.array([idx for idx in range(nPoints) if idx not in train_idx])
    data_to_fit = curr_mani[train_idx].copy()
    data_to_decode = curr_mani[test_idx].copy()

    curr_fit_result = mff.fit_manifold(data_to_fit, curr_fit_params)
    dec_angle, mse = mff.decode_from_passed_fit(data_to_decode, curr_fit_result['tt'][:-1], 
        curr_fit_result['curve'][:-1], ref_angles[test_idx])
    if k in results:
        results[k].append([mse, curr_fit_result['fit_err'], 
            np.array(curr_fit_result['final_knots'])])
    else:
        results[k] = [[mse, curr_fit_result['fit_err'], 
            np.array(curr_fit_result['final_knots'])]]
print('Time ', time.time()-tic)

to_save = {'fit_results' : results, 'session' : session, 'area' : area, 'state' : state, 
    'embeddings_file' : embeddings_fname} 
gff.save_pickle_file(to_save, dir_to_save + '%s_%s_dim%d_trainfrac%.2f_decode_errors_sd%d.pkl'%(
    session, state, fit_dim, train_frac, sd))

rmse_to_plot = np.sqrt([x[0] for x in results[k]])

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
if nTests<20:
    ax.scatter(np.ones_like(rmse_to_plot), rmse_to_plot)
    ax.set_xticks([1])
    # ax.set_xticklabels(['Samples'])
else:
    vp = ax.violinplot([rmse_to_plot], positions=[1], points=100,
        widths=0.75, showmeans=True, showextrema=False, showmedians=True)
ax.set_xlim([0,2])
ax.set_ylabel('Root mean squared error (rad)')
ax.set_ylim([0,1.8])
plt.show()


