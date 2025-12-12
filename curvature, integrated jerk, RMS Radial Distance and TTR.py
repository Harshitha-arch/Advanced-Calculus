import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

df = pd.read_csv("1.csv")
print("CSV loaded successfully!")
print(df.head())

# List of CoP trajectories to analyze
pairs = [('X1','Z1'), ('X3','Z3'), ('X5','Z5'), ('X7','Z7')]
dt = 0.01

def compute_curvature(x, y, dt=0.01, smooth_window=21, polyorder=3):

    x_s = savgol_filter(x, smooth_window, polyorder)
    y_s = savgol_filter(y, smooth_window, polyorder)

    dx = np.gradient(x_s, dt)
    dy = np.gradient(y_s, dt)

    ddx = np.gradient(dx, dt)
    ddy = np.gradient(dy, dt)

    numerator = np.abs(dx * ddy - dy * ddx)
    denominator = (dx**2 + dy**2)**1.5
    denominator[denominator == 0] = 1e-8

    curvature = numerator / denominator
    return curvature

def compute_integrated_jerk(x, y, dt=0.01, smooth_window=21, polyorder=3):

    x_s = savgol_filter(x, smooth_window, polyorder)
    y_s = savgol_filter(y, smooth_window, polyorder)

    dx = np.gradient(x_s, dt)
    dy = np.gradient(y_s, dt)

    ddx = np.gradient(dx, dt)
    ddy = np.gradient(dy, dt)

    dddx = np.gradient(ddx, dt)
    dddy = np.gradient(ddy, dt)

    jerk_mag = np.sqrt(dddx**2 + dddy**2)
    integrated_jerk = np.sum(jerk_mag) * dt

    return integrated_jerk, jerk_mag


def compute_rms_radial_distance(x, y):

    mean_x = np.mean(x)
    mean_y = np.mean(y)

    r = np.sqrt((x - mean_x)**2 + (y - mean_y)**2)
    rms_rd = np.sqrt(np.mean(r**2))

    return rms_rd, r


def compute_excursion_events(r, dt=0.01, threshold=None):

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
            while i < N and above[i]:
                i += 1
            end = i
            ttr_list.append((end - start) * dt)
        else:
            i += 1

    return num_events, ttr_list, threshold


results = {}


for x_col, z_col in pairs:

    print(f"\n\nANALYZING {x_col}_{z_col}")

    x = df[x_col].values
    y = df[z_col].values

    # ---- Curvature ----
    curvature = compute_curvature(x, y, dt)
    mean_curv = np.mean(curvature)
    max_curv = np.max(curvature)

    # ---- Jerk ----
    integrated_jerk, jerk_ts = compute_integrated_jerk(x, y, dt)
    mean_jerk = np.mean(jerk_ts)
    max_jerk = np.max(jerk_ts)

    # ---- RMS Radial Distance ----
    rms_rd, r_series = compute_rms_radial_distance(x, y)

    # ---- Excursion Events ----
    num_events, ttr_list, threshold = compute_excursion_events(r_series, dt)
    mean_ttr = np.mean(ttr_list) if len(ttr_list) else 0
    max_ttr = np.max(ttr_list) if len(ttr_list) else 0

    # Save results
    results[f"{x_col}_{z_col}"] = {
        "mean_curvature": mean_curv,
        "max_curvature": max_curv,
        "integrated_jerk": integrated_jerk,
        "mean_jerk": mean_jerk,
        "max_jerk": max_jerk,
        "RMS_RD": rms_rd,
        "num_excursions": num_events,
        "mean_TTR": mean_ttr,
        "max_TTR": max_ttr,
        "threshold": threshold
    }

    t = np.arange(len(x)) * dt

    # --- 1. Trajectory Plot ---
    plt.figure(figsize=(5,5))
    plt.plot(x, y, label="CoP Path")
    plt.scatter(np.mean(x), np.mean(y), color='red', label='Mean Position')
    plt.title(f"Trajectory Path: {x_col}_{z_col}")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend()
    plt.grid()
    plt.show()

    # --- 2. Curvature Time Series ---
    plt.figure(figsize=(10,4))
    plt.plot(t, curvature)
    plt.title(f"Curvature Time Series: {x_col}_{z_col}")
    plt.xlabel("Time (s)")
    plt.ylabel("Curvature")
    plt.grid()
    plt.show()

    # --- 3. Jerk Time Series ---
    plt.figure(figsize=(10,4))
    plt.plot(t, jerk_ts)
    plt.title(f"Jerk Time Series: {x_col}_{z_col}")
    plt.xlabel("Time (s)")
    plt.ylabel("Jerk Magnitude")
    plt.grid()
    plt.show()

    # --- 4. Radial Distance + Threshold ---
    plt.figure(figsize=(10,4))
    plt.plot(t, r_series, label="Radial Distance")
    plt.axhline(threshold, color='red', linestyle='--', label="Threshold")
    plt.title(f"Radial Distance & Excursions: {x_col}_{z_col}")
    plt.xlabel("Time (s)")
    plt.ylabel("Radial Distance")
    plt.legend()
    plt.grid()
    plt.show()


print("\n\nFINAL RESULTS\n")

for key, val in results.items():
    print(f"\n---- {key} ----")
    for metric, number in val.items():
        print(f"{metric:20s} : {number:.6f}")
