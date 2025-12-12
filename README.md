# Center of Pressure (CoP) Stability Analysis

This repository contains Python code for analyzing postural stability using Center of Pressure (CoP) trajectories.  
The script computes four key balance metrics from raw CoP data collected from force platforms.



## Included Analysis Metrics

### 1. Curvature
Measures how sharply the CoP path bends over time.  
- Computed from smoothed CoP derivatives  
- Output: mean curvature, max curvature  
- Visualization: Curvature Time Series


### 2. Integrated Jerk
Represents the smoothness of postural control.  
- Jerk = third derivative of CoP displacement  
- Higher jerk = more abrupt corrective movements  
- Output: integrated jerk (scalar), jerk time series  
- Visualization: Jerk Time Series Plot



### 3. RMS Radial Distance (RMS-RD)
Quantifies sway magnitude relative to the mean CoP position.  
- Computed as RMS distance from mean (X,Y)  
- Output: RMS radial distance  
- Visualization: Radial Distance plot



### 4. Time to Recovery (TTR) from Excursions
Detects instability events where CoP exceeds a statistical threshold.  
- Threshold = mean(r) + 2·std(r)  
- Computes number of excursions & duration until return  
- Output: num excursions, mean/max TTR  
- Visualization: Radial Distance with Threshold



## Generated Plots
The script outputs 4 plots for each CoP pair (`X1–Z1`, `X3–Z3`, `X5–Z5`, `X7–Z7`):

- Trajectory Path (X vs Y)  
- Curvature Time Series  
- Jerk Time Series  
- Radial Distance + Threshold (Excursions)

---

## Input Data
Place your input CoP file as:
Your CSV must contain columns:  
`X1, Z1, X3, Z3, X5, Z5, X7, Z7`



## How to Run

1. Install required libraries:
2. Run the script:
3. Output includes:

- Visual plots  
- Final numerical metrics printed in the terminal



## File Explanation

| File | Description |
|------|-------------|
| `analysis.py` | Main script performing all computations and visualizations |
| `1.csv` | Input CoP trajectory data |
| results (terminal) | Summary of curvature, jerk, RMS-RD, and excursion metrics |

---

## Summary
This project provides a full analysis pipeline for studying postural stability using Center of Pressure movement.  
It computes and visualizes four major stability metrics commonly used in biomechanics and rehabilitation research.
