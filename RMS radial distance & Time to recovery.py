import pandas as pd
import numpy as np


df = pd.read_csv("1.csv")

print("CSV loaded successfully.")
print(df.head())






def compute_rms_radial_distance(x, y):
    """
    Computes RMS Radial Distance from mean CoP position.

    RMS-RD = sqrt( mean( (x - mean_x)^2 + (y - mean_y)^2 ) )
    """

    mean_x = np.mean(x)
    mean_y = np.mean(y)

    r = np.sqrt((x - mean_x)**2 + (y - mean_y)**2)
    rms_rd = np.sqrt(np.mean(r**2))

    return rms_rd, r







def compute_excursion_events(r, dt=0.01, threshold=None):
    """
    Computes:
        - Number of excursion events
        - Time-to-recovery list (seconds)
        - Threshold used

    Parameters:
        r : radial distance time series
        dt : sampling interval
        threshold : optional; if None uses mean + 2*std

    Returns:
        num_events : count
        ttr_list   : list of durations
        threshold  : threshold used
    """

    if threshold is None:
        threshold = np.mean(r) + 2*np.std(r)

    above = r > threshold
    num_events = 0
    ttr_list = []
    i = 0
    N = len(r)

    while i < N:
        if above[i]:  
            num_events += 1
            start = i

            # move until r recovers below threshold
            while i < N and above[i]:
                i += 1
            
            end = i
            ttr = (end - start) * dt
            ttr_list.append(ttr)
        else:
            i += 1

    return num_events, ttr_list, threshold





pairs = [('X1','Z1'), ('X3','Z3'), ('X5','Z5'), ('X7','Z7')]
postural_results = {}

dt = 1  # adjust if your sampling rate differs

for x_col, z_col in pairs:
    x = df[x_col].values
    y = df[z_col].values

    # ---- RMS Radial Distance ----
    rms_rd, r_series = compute_rms_radial_distance(x, y)

    # ---- Excursion Events & TTR ----
    num_events, ttr_list, used_threshold = compute_excursion_events(r_series, dt)

    postural_results[f"{x_col}_{z_col}"] = {
        "RMS_radial_distance": rms_rd,
        "num_excursions": num_events,
        "mean_TTR": np.mean(ttr_list) if len(ttr_list)>0 else 0,
        "max_TTR": np.max(ttr_list) if len(ttr_list)>0 else 0,
        "threshold": used_threshold
    }


print("\n===== POSTURAL STABILITY RESULTS =====")
for trial, values in postural_results.items():
    print(f"\nTrajectory {trial}:")
    print(f"RMS Radial Distance : {values['RMS_radial_distance']:.6f}")
    print(f"Excursion Events     : {values['num_excursions']}")
    print(f"Mean TTR (s)         : {values['mean_TTR']:.4f}")
    print(f"Max TTR (s)          : {values['max_TTR']:.4f}")
    print(f"Threshold (radius)   : {values['threshold']:.6f}")
