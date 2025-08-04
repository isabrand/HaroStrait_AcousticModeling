import arlpy.uwapm as pm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import warnings

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', titlesize=20)
plt.rc('axes', labelsize=18)
plt.rc('legend', fontsize=20)

np.set_printoptions(suppress=True)

# ================================================TABLE OF CONTENTS================================================
# 1. Environmental Set Up Functions
#    - test_my_file()
#    - set_up_environment(field_depth, field_range, center_freq, freqs_checked, c_sea, surface_reflection, bottom_reflection, num_receivers, depth_start, d_spacing, range_start, r_spacing)
#    - set_up_source_array(num_receivers, num_depths, num_ranges, test_d_min, test_d_max, test_r_min, test_r_max)
#    - calculate_geometry(source_pos, receiver_pos, source_shape, field_depth)
#    - monte_carlo_setup
# 2. Beamforming Calculation Functions
#    - p_field_via_math
#    - p_field_bellhop
#    - get_first_few_arrivals
#    - add_noise_to_field
#    - calculate_K_matrices
#    - get_replicate
#    - calculate_beamformer
# 3. Plotting Functions
#    - plot_image_source_geometry
#    - plot_beamformer
#    - plot_beamformer_crosssections
#    - plot_monte_carloesque
# ================================================ENVIRONMENTAL SET UP FUNCTIONS================================================

def test_my_file():
    print("This is a test function to check if the file is being imported correctly.")
    return "File imported successfully!"

def set_up_environment(field_depth, field_range, center_freq, freqs_checked, c_sea, surface_coef, bottom_coef, num_receivers, depth_start, d_spacing, range_start, r_spacing):
    # field depth and field ranfe are in meters, for scale, the mooring to shore is ~800m
    # c_sea, speed of sound in seawater, is in m/s
    wavelengths = c_sea / freqs_checked              # (m), wavelength
    k = 2 * np.pi / wavelengths             # (rad/m), wavenumber

    print(f"Wavelengths: {wavelengths}\nWavenumbers: {k}")

    print(f"Field depth: {field_depth} m")
    print(f"Field range: {field_range} m")
    print(f"Center frequency: {center_freq} Hz")
    print(f"Sound speed average: {c_sea} m/s")

    max_d = c_sea/(2*center_freq)          # (m), max spacing to avoid spatial aliasing
    receiver_pos = np.array([[range_start + i *r_spacing, depth_start + i * d_spacing] for i in range(num_receivers)])
    ff_start = (2*(int(num_receivers*d_spacing))^2)/(c_sea / center_freq)  # fraunhofer distance (m)

    print(f"Maximum receiver spacing to avoid spatial aliasing for center frequency {center_freq} Hz: {max_d} m")
    print(f"Fraunhofer distance for center frequency {center_freq} Hz: {ff_start} m")
    print(f"There are {num_receivers} receivers in the array.")
    #print(f"The receivers start with {depth_start} m depth and {range_start} m range, with a spacing of {d_spacing} m depth and {r_spacing} m range per hydrophone.")
    print(f"Receiver positions: {receiver_pos}")

    return receiver_pos, c_sea, wavelengths, k, max_d, ff_start


def set_up_source_array(num_receivers, num_depths, num_ranges, test_d_min, test_d_max, test_r_min, test_r_max):
    # set up source array, which is a grid of sources in the field that we want to test

    d_grid, r_grid = np.meshgrid(np.linspace(test_d_min, test_d_max, num_depths), np.linspace(test_r_min, test_r_max, num_ranges), indexing='ij')
    source_pos = np.stack((r_grid, d_grid), axis=-1)
    source_shape = source_pos.shape
    source_shape = (source_shape[0], source_shape[1], num_receivers)

    # CHANGE HERE IF YOU WANT A RANDOM SOURCE OR FIXED SOURCE AFTER SETTING UP THE SOURCE ARRAY
    #source_pos[0, 0, :] = [random.randint(10, field_range), random.randint(1, 80)]
    #source_pos[0, 0, :] = [1100, 25]

    print(f"\nThe test source positions for this run are: \n{source_pos} m")

    return source_pos, source_shape, num_depths, num_ranges


def calculate_geometry(source_pos, receiver_pos, source_shape, field_depth):
    r_d, r_s, r_b = np.zeros(source_shape, dtype=float), np.zeros(source_shape, dtype=float), np.zeros(source_shape, dtype=float)
    ba_d, ba_s, ba_b = np.zeros(source_shape, dtype=float), np.zeros(source_shape, dtype=float), np.zeros(source_shape, dtype=float)
    ba_d_deg, ba_s_deg, ba_b_deg = np.zeros(source_shape, dtype=float), np.zeros(source_shape, dtype=float), np.zeros(source_shape, dtype=float)

    source_pos_surf, source_pos_bot = source_pos.copy(), source_pos.copy()
    source_pos_surf[:, :, 1] *= -1
    source_pos_bot[:, :, 1] = field_depth + (field_depth - source_pos_bot[:, :, 1])

    source_pos_dir_flat = source_pos.reshape(-1, 2)
    source_pos_surf_flat = source_pos_surf.reshape(-1, 2)
    source_pos_bot_flat = source_pos_bot.reshape(-1, 2)


    for d, s, b in zip(source_pos_dir_flat, source_pos_surf_flat, source_pos_bot_flat):
        row_idx, col_idx = np.where((source_pos == d).all(axis=2))

        # Calculate the distance from source to EACH receiver for a direct path, surface reflection, and bottom reflection
        r_d[row_idx, col_idx] = np.sqrt((receiver_pos[:, 0] - d[0])**2 + (d[1] - receiver_pos[:, 1])**2)
        r_s[row_idx, col_idx] = np.sqrt((receiver_pos[:, 0] - s[0])**2 + (s[1] - receiver_pos[:, 1])**2)
        r_b[row_idx, col_idx] = np.sqrt((receiver_pos[:, 0] - b[0])**2 + (b[1] - receiver_pos[:, 1])**2)

        # the bearing angle θi​ is the azimuth from EACH receiver to the source, measured clockwise from the x-axis
        ba_d[row_idx, col_idx] = np.arctan2((receiver_pos[:, 1] - d[1]), (d[0] - receiver_pos[:, 0]))
        ba_s[row_idx, col_idx] = np.arctan2((receiver_pos[:, 1] - s[1]), (s[0] - receiver_pos[:, 0]))
        ba_b[row_idx, col_idx] = np.arctan2((receiver_pos[:, 1] - b[1]), (b[0] - receiver_pos[:, 0]))

        ba_d_deg[row_idx, col_idx], ba_s_deg[row_idx, col_idx], ba_b_deg[row_idx, col_idx] = np.degrees(ba_d[row_idx, col_idx]), np.degrees(ba_s[row_idx, col_idx]), np.degrees(ba_b[row_idx, col_idx])

    all_r = [r_d, r_s, r_b]
    all_ba_deg = [ba_d_deg, ba_s_deg, ba_b_deg]
    all_source_pos_flat = [source_pos_dir_flat, source_pos_surf_flat, source_pos_bot_flat]

    return all_r, all_ba_deg, all_source_pos_flat

    
def monte_carlo_setup(receiver_pos, field_range, field_depth, num_zones, num_test_points):
    zone_width = field_range / num_zones

    zones = [(i * zone_width, (i + 1) * zone_width, 0, 150) for i in range(num_zones)]

    fig, ax = plt.subplots(figsize=(20, 20))

    ax.scatter(receiver_pos[:, 0], receiver_pos[:, 1], color='black', marker="o", s=10, label='Hydrophones')

    ax.fill_between(x=[-10, field_range], y1=-240, y2=0, color='aliceblue', alpha=0.3)    
    ax.fill_between(x=[-10, field_range], y1=0, y2=240, color='lightsteelblue', alpha=0.3)
    ax.fill_between(x=[-10, field_range], y1=240, y2=480, color='burlywood', alpha=0.3)

    ax.set_xlim(-10, field_range)
    ax.set_ylim(-50, (field_depth)+50)
    ax.axhline(y=240, color='brown', linestyle='-', linewidth=0.5)
    ax.axhline(y=0, color='blue', linestyle='-', linewidth=0.5)
    ax.axhline(y=zones[0][3], color='gray', linestyle='--', linewidth=0.5)

    for zone in zones:
        ax.axvline(x=zone[1], ymin=0.41, ymax=0.85, color='gray', linestyle='--', linewidth=0.5)

    ax.invert_yaxis()

    monte_carlo_points = [
        np.column_stack((
            np.random.uniform(zone[0], zone[1], num_test_points),
            np.random.uniform(zone[2], zone[3], num_test_points)
        )) for zone in zones
    ]

    # Combine all zones into one array if needed
    zone_colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']

    for mc, color in zip(monte_carlo_points, zone_colors):
        for pt in tqdm(mc.reshape(-1, 2), desc="Monte Carlo Simulation", unit="point"):
            ax.scatter(pt[0], pt[1], color=color, marker='o', s=10)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="lower right", fontsize=10)

    ax.set_aspect('equal')
    ax.set_xlabel("Range (m)", fontsize=14)
    ax.set_ylabel("Depth (m)", fontsize=14)
    ax.set_title("Monte Carlo Field Breakdown", fontsize=16)
    plt.show()

    return monte_carlo_points


# ================================================BEAMFORMING CALCULATION FUNCTIONS================================================

def p_field_via_math(receiver_pos, pt_pos, k, freqs_checked, r, ba_deg, c_sea, surface_reflection, bottom_reflection, get_rep, p_dir_only, add_noise, snr_db):
    P_per_freqs = np.zeros((len(freqs_checked), len(receiver_pos)), dtype=complex)
    r_direct, r_surface, r_bottom = r[0], r[1], r[2]
    ba_d_deg, ba_s_deg, ba_b_deg = ba_deg[0], ba_deg[1], ba_deg[2]
    arrivals_math = pd.DataFrame(columns=['toa', 'arr_amp', 'ang_of_arr', 'surf_bounce', 'bot_bounce'])

    for k_ind, k_use in enumerate(k):
        Pdirect = (np.exp(1j * k_use * r_direct)) / r_direct       # THE DIRECT PATH! ALWAYS INCORPORATE AT BARE MINIMUM
        Psurface = ((surface_reflection)*(np.exp(1j * k_use * r_surface))) / r_surface
        Pbottom = ((bottom_reflection)*(np.exp(1j * k_use * r_bottom))) / r_bottom
        #print(f"Pdirect is {Pdirect}\nPsurface is {Psurface}\nPbottom is {Pbottom}")
        if p_dir_only:
            P = Pdirect
            if k_ind==0 and not get_rep:
                arrivals_math['toa'] = [r_direct[0]/c_sea]
                arrivals_math['ang_of_arr'] = [ba_d_deg[0]]
                arrivals_math['arr_amp'] = [[complex(round(x.real, 8), round(x.imag, 8)) for x in Pdirect][0]]
                arrivals_math['surf_bounce'] = [0]
                arrivals_math['bot_bounce'] = [0]
        else:
            P = Pdirect + Psurface + Pbottom
            if k_ind == 0 and not get_rep: 
                arrivals_math['toa'] = [r_direct[0]/c_sea, r_surface[0]/c_sea, r_bottom[0]/c_sea]
                arrivals_math['ang_of_arr'] = [ba_d_deg[0], ba_s_deg[0], ba_b_deg[0]]
                arrivals_math['arr_amp'] = [[complex(round(x.real, 8), round(x.imag, 8)) for x in Pdirect][0], [complex(round(x.real, 8), round(x.imag, 8)) for x in Psurface][0], [complex(round(x.real, 8), round(x.imag, 8)) for x in Pbottom][0]]
                arrivals_math['surf_bounce'] = [0, 1, 0]
                arrivals_math['bot_bounce'] = [0, 0, 1]

        if add_noise:
            print("Adding noise to P_signals before calculating K...")
            P = add_noise_to_field(P, snr_db)

        P_per_freqs[k_ind] = P

    return P_per_freqs, arrivals_math


def p_field_bellhop(source_pos, receiver_pos, ssp, bathy, bottom_reflection, freqs_checked, add_noise, snr_db, first_few, num):
    P_vec = np.zeros((len(freqs_checked), len(receiver_pos)), dtype=complex)
    arrivals_bellhop = pd.DataFrame(columns=['toa', 'arr_amp', 'ang_of_arr', 'surf_bounce', 'bot_bounce'])

    for freq_ind, freq in enumerate(freqs_checked):
        for r_i, rec in enumerate(receiver_pos):
            env = pm.create_env2d(
                depth=bathy,
                soundspeed=ssp,
                bottom_soundspeed=1480,
                bottom_density=1200,
                bottom_absorption=1-bottom_reflection,
                frequency=freq,
                tx_depth=source_pos[1],
                rx_depth=rec[1],
                rx_range=source_pos[0]
            )
            env['nbeams'] = 6001

            #pm.plot_env(env, width=900)    
            #pm.print_env(env)
            #rays = pm.compute_eigenrays(env)
            
            arrivals = pm.compute_arrivals(env)
            arrivals = arrivals[['time_of_arrival', 'arrival_amplitude', 'angle_of_arrival', 'surface_bounces', 'bottom_bounces']]
            if first_few:
                arrivals = get_first_few_arrivals(arrivals, num)
                arrivals = arrivals.astype({'time_of_arrival': 'float64', 'arrival_amplitude': 'complex', 'angle_of_arrival': 'float64', 'surface_bounces': 'float64', 'bottom_bounces': 'float64'})
                arrivals_bellhop["toa"], arrivals_bellhop["arr_amp"], arrivals_bellhop["ang_of_arr"], arrivals_bellhop["surf_bounce"], arrivals_bellhop["bot_bounce"] = arrivals['time_of_arrival'], arrivals['arrival_amplitude'], arrivals['angle_of_arrival'], arrivals['surface_bounces'], arrivals['bottom_bounces']
            
            #pm.plot_rays(rays, env=env, width=900)
            #pm.plot_arrivals(arrivals, dB=True, width=900)
            P_vec[freq_ind, r_i] = arrivals['arrival_amplitude'].sum()

        if add_noise:
            print("Adding noise to P_signals before calculating K...")
            P_vec = add_noise_to_field(P_vec, snr_db)

    return P_vec, arrivals_bellhop


def get_first_few_arrivals(arrivals, count):
    surf_paths = [0, 0, 1, 1]
    bot_paths = [0, 1, 0, 1]
    arrivals = arrivals.sort_values(by='time_of_arrival')
    dirs, surfs, bots, surfbots, botsurfs = pd.DataFrame(columns=arrivals.columns), pd.DataFrame(columns=arrivals.columns), pd.DataFrame(columns=arrivals.columns), pd.DataFrame(columns=arrivals.columns), pd.DataFrame(columns=arrivals.columns)
    for _, row in arrivals.iterrows():
        if row['surface_bounces'] == surf_paths[0] and row['bottom_bounces'] == bot_paths[0]:
            dirs.loc[len(dirs)] = row
        elif row['surface_bounces'] == surf_paths[1] and row['bottom_bounces'] == bot_paths[1]:
            surfs.loc[len(surfs)] = row
        elif row['surface_bounces'] == surf_paths[2] and row['bottom_bounces'] == bot_paths[2]:
            bots.loc[len(bots)] = row
        elif row['surface_bounces'] == surf_paths[3] and row['bottom_bounces'] == bot_paths[3]:
            if row['angle_of_arrival'] > 0:
                botsurfs.loc[len(botsurfs)] = row
            elif row['angle_of_arrival'] < 0:
                surfbots.loc[len(surfbots)] = row

    subset = pd.DataFrame(columns=arrivals.columns)
    def get_strongest_arrival(arrivals_per_path):
        if len(arrivals_per_path) == 0: return None
        elif len(arrivals_per_path) > 2: return arrivals_per_path.iloc[np.argmax(arrivals_per_path['arrival_amplitude'])]
        else: return arrivals_per_path.iloc[0]

    for arrivals_per_path in [dirs, surfs, bots, surfbots, botsurfs]:
        strongest = get_strongest_arrival(arrivals_per_path)
        if strongest is not None:
            subset.loc[len(subset)] = strongest

    subset = subset.sort_values(by='time_of_arrival').head(count)

    return subset


# ADD NOISE, set SNR higher for noise to be less significant
def add_noise_to_field(field, snr_db):
    flattened = field.reshape(-1, field.shape[-1]) 
    noisy_field = flattened.copy()

    for ind, sig in enumerate(flattened):
        signal_power = np.mean(np.abs(sig)**2)
        noise_power = signal_power / (10 ** (snr_db / 10))

        noise_real = np.random.normal(0, np.sqrt(noise_power/2), len(sig))        # complex Gaussian noise
        noise_imag = np.random.normal(0, np.sqrt(noise_power/2), len(sig))        # complex Gaussian noise
        
        noise = noise_real + 1j * noise_imag
        noisy_field[ind] = sig + noise

    return noisy_field.reshape(field.shape)


def calculate_K_matrices(P, num_frequencies, row_idx, col_idx):
    # Calculate K as the cross spectral density matrix of P_signals

    K_freq_avg = np.zeros((P.shape[3], P.shape[3]), dtype=complex)
    for freq_idx in range(num_frequencies):
        sig = P[row_idx, col_idx, freq_idx]
        K_freq_avg += np.outer(sig, np.conj(sig))       # add together the K matrices for each frequency snapshot

    K_freq_avg /= num_frequencies                       

    trace_K = np.trace(K_freq_avg)
    if np.abs(trace_K) > 1e-10:
        K_freq_avg = K_freq_avg / trace_K          # Normalize matrix K, don't divide by zero/very small numbers

    return K_freq_avg


def get_replicate(field_depth, search_ranges, search_depths, receiver_pos, freqs_checked, k, c_sea, surface_reflection, bottom_reflection, num_receivers, get_rep=True):
    search_field_greens = np.zeros((len(search_ranges), len(search_depths), len(freqs_checked), len(receiver_pos)), dtype=complex)
    search_field_multi = np.zeros((len(search_ranges), len(search_depths), len(freqs_checked), len(receiver_pos)), dtype=complex)

    # form the replica vector for each freq. over the search field
    for r_ind, r in enumerate(search_ranges):
        for d_ind, d in enumerate(search_depths):
            search_source = np.array([[[r, d]]])
            r_paths_w, ba_deg_paths_w, all_source_pos_flat_w = calculate_geometry(search_source, receiver_pos, source_shape=(1, 1, num_receivers), field_depth=field_depth)
            # solve using only the direct path and e^ikr -> GREENS FUNCTION
            r_w = [r_paths_w[0], r_paths_w[1], r_paths_w[2]]
            ba_deg_w = [ba_deg_paths_w[0], ba_deg_paths_w[1], ba_deg_paths_w[2]]
            search_field_greens[r_ind, d_ind], arrivals_math = p_field_via_math(receiver_pos, search_source, k, freqs_checked, r_w, ba_deg_w, c_sea, surface_reflection, bottom_reflection, get_rep, p_dir_only=True, add_noise=False, snr_db=0)
            # solve using the pressure field calculation mathematically incorporating surface and bottom reflections
            search_field_multi[r_ind, d_ind], arrivals_math = p_field_via_math(receiver_pos, search_source, k, freqs_checked, r_w, ba_deg_w, c_sea, surface_reflection, bottom_reflection, get_rep, p_dir_only=False, add_noise=False, snr_db=0)
    return search_field_greens, search_field_multi, search_depths, search_ranges


def calculate_beamformer(method, pt_pos_all, source_pos, search_field, P_signals, freqs_checked, K_mat, search_depths, search_ranges):
    B = np.zeros((source_pos.shape[0], source_pos.shape[1], len(search_depths), len(search_ranges)), dtype=complex)

    for pt_pos in pt_pos_all:

        row_idx, col_idx = np.where((source_pos == pt_pos).all(axis=2))
        pt_B = np.zeros((len(search_depths), len(search_ranges)), dtype=complex)
        K = K_mat[row_idx, col_idx]

        for r_ind, r in enumerate(search_ranges):
            for d_ind, d in enumerate(search_depths):
                B_freqs = []
                for f_ind, freq in enumerate(freqs_checked):
                    sig = P_signals[row_idx, col_idx, f_ind].reshape(-1, 1)
                    w = search_field[r_ind, d_ind, f_ind]

                    # Bartlett Beamformer
                    if method == 'bartlett':
                        numerator = np.abs(np.vdot(sig, w))**2                          # correlates the signal with the replicate and sees how well they match
                        denominator = np.vdot(w, w)
                        B_freqs.append((numerator / denominator))

                    # Classic Beamformer    -> B = w^H * K * w
                    elif method == 'classic':
                        w = w.reshape(-1, 1)
                        numerator = np.real((w.conj().T @ K @ w).item())       # represents the power of the system, correlation of w and signal, taking into account the covariance matrix
                        denominator = np.vdot(w.flatten(), w.flatten())
                        B_freqs.append((numerator / denominator))
                    
                    elif method == 'mvdr':
                        w = w.reshape(-1, 1)
                        K_inv = np.linalg.pinv(K)                   # Regularized inverse
                        B_freqs.append(1 / (w.conj().T @ K_inv @ w).real)

                    else:
                        raise ValueError("Unknown beamforming method.")

                pt_B[d_ind, r_ind] = np.mean(B_freqs, dtype=complex)
        B[row_idx, col_idx] = pt_B
    return B


# ================================================PLOTTING FUNCTIONS================================================

def plot_image_source_geometry(receiver_pos, pt_pos_all, pt_pos_surf_all, pt_pos_bot_all, field_range, field_depth, num_receivers):
    plt.figure(figsize=(12, 8))

    for i in range(num_receivers):
        hp = receiver_pos[i, :]  # Current receiver position
        
        # Plot all receivers first (background)
        plt.scatter(hp[0], hp[1], color='black', marker='o', label='Receiver', s=25)
        
        # Plot ALL source points and paths for this receiver
        for j in range(len(pt_pos_all)):
            pt = pt_pos_all[j]
            pt_surf = pt_pos_surf_all[j]
            pt_bot = pt_pos_bot_all[j]
            
            # Direct path
            plt.scatter(pt[0], pt[1], color='red', marker='o', s=20, label='Source - Direct Path')
            plt.plot([pt[0], hp[0]], [pt[1], hp[1]], 'k--', linewidth=0.5)

            # Surface reflected path
            plt.scatter(pt_surf[0], pt_surf[1], color='cadetblue', marker='o', s=20, label='Source Projection - Surface Reflection')
            plt.plot([pt_surf[0], hp[0]], [pt_surf[1], hp[1]], color='gray', linestyle='--', linewidth=0.5)

            # Bottom reflected path
            plt.scatter(pt_bot[0], pt_bot[1], color='olive', marker='o', s=20, label='Source Projection - Bottom Reflection')
            plt.plot([pt_bot[0], hp[0]], [pt_bot[1], hp[1]], color='gray', linestyle='--', linewidth=0.5)


        plt.fill_between(x=[-10, field_range], y1=-240, y2=0, color='aliceblue', alpha=0.3)    
        plt.fill_between(x=[-10, field_range], y1=0, y2=240, color='lightsteelblue', alpha=0.3)
        plt.fill_between(x=[-10, field_range], y1=240, y2=480, color='burlywood', alpha=0.3)
        plt.xlim(-10, field_range)
        plt.ylim(-240, (field_depth)+240)
        plt.axhline(y=0, color='blue', linestyle='-', linewidth=0.3)
        plt.axhline(y=240, color='brown', linestyle='-', linewidth=0.3)
        plt.gca().invert_yaxis()
        plt.gca().set_aspect('equal')
        plt.gca().set_xlabel("Range (m)", fontsize=12)
        plt.gca().set_ylabel("Depth (m)", fontsize=12)
        plt.gca().set_title(f"Calculated paths of arrival for each test source to each hydrophone", fontsize=18)

        return


def plot_beamformer(pt_pos_all, source_pos, B, search_ranges, search_depths, receiver_pos, num_ranges, num_depths,title):
    """
    Plot the beamformer results for each test source.
    
    Parameters:
    - pt_pos_all: Array of points tested
    - source_pos: Source positions
    - B: Beamformer results
    - search_ranges: Ranges for the search field.
    - search_depths: Depths for the search field.
    - receiver_pos: Receiver positions.
    """

    fig, axs = plt.subplots(num_depths, num_ranges, figsize=(6*num_ranges, 3*num_depths), sharex=True, sharey=True)
    if len(pt_pos_all) == 1: axs = np.array([axs])
    axs = axs.ravel()
    rec_mean = np.array(np.mean(receiver_pos, axis=0))

    for count, (pt_pos, ax) in enumerate(zip(pt_pos_all, axs)):
        row_idx, col_idx = np.where((source_pos == pt_pos).all(axis=2))
        pt_B = np.real(B[row_idx[0], col_idx[0]])

        im = ax.imshow(
            pt_B,
            extent=[0, search_ranges[-1], search_depths[-1], 0],
            aspect='auto',
            cmap='viridis',
        )

        peak_index = np.unravel_index(np.argmax(pt_B), pt_B.shape)
        peak_depth, peak_range = search_depths[peak_index[0]], search_ranges[peak_index[1]]

        range_err = peak_range - pt_pos[0]
        depth_err = peak_depth - pt_pos[1]
        dist_err = np.sqrt(range_err**2 + depth_err**2)
        print(f"Range error: {range_err:.2f} m, Depth error: {depth_err:.2f} m, Dist error: {dist_err:.2f} m")
        if np.abs(range_err) <= 150 and np.abs(depth_err) <= 20: outline_color = 'green'
        else: outline_color = 'red'

        v1 = pt_pos - rec_mean
        v2 = np.array([peak_range, peak_depth]) - rec_mean
        cos_theta = (np.dot(v1, v2)) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
        angle_deg = np.degrees(angle_rad)
        print(f"Angle between source and peak: {angle_deg.item():.2f} degrees")

        if angle_deg < 1: line_col = 'green'
        else: line_col = 'red'

        ax.plot([np.mean(receiver_pos[:,0]), pt_pos[0]], [np.mean(receiver_pos[:,1]), pt_pos[1]], color='black', linestyle='--', linewidth=1, label='Source to Receiver Path')
        ax.axline(xy1=(np.mean(receiver_pos[:,0]), np.mean(receiver_pos[:,1])), xy2=(peak_range, peak_depth), color=line_col, linestyle='--', linewidth=1)

        ax.scatter(pt_pos[0], pt_pos[1], color='red', marker='o', label=f"Test Source ({pt_pos[0]} m, {pt_pos[1]} m)", s=70)
        ax.set_title(f'Test at {pt_pos[0]:0.0f}m, {pt_pos[1]:0.0f}m --- Est at {peak_range:0.0f}m, {peak_depth:0.0f}m', fontsize=10)
        ax.scatter(peak_range, peak_depth, color='red', marker='x', label=f'Estimated Peak ({peak_range:0.0f} m, {peak_depth:0.0f} m)', s=120)
        ax.scatter(*receiver_pos.T, s=100, color='white', marker='.', label="Receivers")

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.suptitle(title, fontsize=20)

    return peak_depth, peak_range, dist_err


def plot_beamformer_crosssections(pt_pos_all, source_pos, B, search_ranges, search_depths, num_depths, num_ranges, peak_range, peak_depth, dist_err):
    fig, ax = plt.subplots(len(pt_pos_all), 2, figsize=(10, 8), sharey=True)

    for count, pt_pos in enumerate(pt_pos_all):
        range_err = peak_range - pt_pos[0]
        depth_err = peak_depth - pt_pos[1]
        dist_err = np.sqrt(range_err**2 + depth_err**2)
        print(f"Range error: {range_err:.2f} m, Depth error: {depth_err:.2f} m, Dist error: {dist_err:.2f} m")

        row_idx, col_idx = np.where((source_pos == pt_pos).all(axis=2))
        pt_B = B[row_idx[0], col_idx[0]]

        # Find closest indices for cross-sections
        cs_depth = np.argmin(np.abs(search_depths - pt_pos[1]))
        cs_range = np.argmin(np.abs(search_ranges - pt_pos[0]))
        
        # Get peak location
        peak_index = np.unravel_index(np.argmax(pt_B), pt_B.shape)
        peak_depth = search_depths[peak_index[0]]
        peak_range = search_ranges[peak_index[1]]
                
        # Range cross-section
        ax[count, 0].plot(search_ranges, pt_B[cs_depth, :])
        if count == 0:
            ax[count, 0].axvline(x=pt_pos[0], color='red', linestyle='--', label='Test Source')
            ax[count, 0].axvline(x=peak_range, color='green', linestyle=':', label='Estimated Location', linewidth=4)
        else:
            ax[count, 0].axvline(x=pt_pos[0], color='red', linestyle='--')
            ax[count, 0].axvline(x=peak_range, color='green', linestyle=':', linewidth=4)
        ax[count, 0].set_xlabel('Range (m)', fontsize=12)
        ax[count, 0].set_ylabel('Amp', fontsize=12)
        ax[count, 0].set_title(f'Range values at {pt_pos[1]:.1f} m of depth', fontsize=15)
        ax[count, 0].grid(True)

        # Depth cross-section
        ax[count, 1].plot(search_depths, pt_B[:, cs_range])
        ax[count, 1].axvline(x=pt_pos[1], color='red', linestyle='--'), ax[count, 1].axvline(x=peak_depth, color='green', linestyle=':', linewidth=4)
        ax[count, 1].set_xlabel('Depth (m)', fontsize=12), ax[count, 1].set_ylabel('Amp', fontsize=12)
        ax[count, 1].set_title(f'Depth values at {pt_pos[0]:.1f} m of range', fontsize=15)
        ax[count, 1].grid(True)

    fig.suptitle('Beampower Cross-Sections', fontsize=20)
    fig.legend(loc='lower center', fontsize=10)
    fig.tight_layout()

    plt.show()


def plot_monte_carloesque(pt_pos_all, source_pos, B, search_ranges, search_depths, receiver_pos, num_ranges, num_depths,title):
    # Plot the beamformer results for a bunch of tests all on one plot. Don't need to plot the actual ambiguity surface, just success or fail.

    fig, axs = plt.subplots(1, 1, figsize=(6*num_ranges, 5*num_depths), sharex=True, sharey=True)
    if len(pt_pos_all) == 1:    axs = np.array([axs])
    axs = axs.ravel()

    for count, (pt_pos, ax) in enumerate(zip(pt_pos_all, axs)):
        row_idx, col_idx = np.where((source_pos == pt_pos).all(axis=2))
        pt_B = np.real(B[row_idx[0], col_idx[0]])

        im = ax.imshow(
            pt_B,
            extent=[0, search_ranges[-1], search_depths[-1], 0],
            aspect='auto',
            cmap='viridis',
        )

        peak_index = np.unravel_index(np.argmax(pt_B), pt_B.shape)
        peak_depth, peak_range = search_depths[peak_index[0]], search_ranges[peak_index[1]]

        range_err = peak_range - pt_pos[0]
        depth_err = peak_depth - pt_pos[1]
        dist_err = np.sqrt(range_err**2 + depth_err**2)
        #print(f"Range error: {range_err:.2f} m, Depth error: {depth_err:.2f} m, Dist error: {dist_err:.2f} m")
        if np.abs(range_err) <= 150 and np.abs(depth_err) <= 20: outline_color = 'green'
        else: outline_color = 'red'
        for spine in ax.spines.values():
            spine.set_color(outline_color)
            spine.set_linewidth(3)
        
        ax.scatter(pt_pos[0], pt_pos[1], color='red', marker='o', label=f"Test Source ({pt_pos[0]} m, {pt_pos[1]} m)", s=70)
        ax.set_title(f'Test at {pt_pos[0]:0.0f}m, {pt_pos[1]:0.0f}m --- Est at {peak_range:0.0f}m, {peak_depth:0.0f}m', fontsize=12)
        ax.scatter(peak_range, peak_depth, color='red', marker='x', label=f'Estimated Peak ({peak_range:0.0f} m, {peak_depth:0.0f} m)', s=120)
        ax.scatter(*receiver_pos.T, s=100, color='white', marker='.', label="Receivers")
        fig.suptitle(title, fontsize=20)

    return peak_depth, peak_range, dist_err

