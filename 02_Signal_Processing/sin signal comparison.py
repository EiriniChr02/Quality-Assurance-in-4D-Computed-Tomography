import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, butter, filtfilt

# 1. Loading Data
df = pd.read_csv('waveform_export_SIN.csv')
time = df['time_s'].values
raw_signal = df['Ant.-Post. [millimeter]'].values

# 2. 
p_high = np.percentile(raw_signal, 99)
p_low  = np.percentile(raw_signal, 1)
midpoint = (p_high + p_low) / 2
centered = raw_signal - midpoint

# 3. Low-pass Butterworth 2Hz
fs = 1 / (time[1] - time[0]) 
b, a = butter(4, 2.0 / (fs/2), btype='low')
filtered = filtfilt(b, a, centered)

# 4. Parameters
nominal_amplitude = 10.0
nominal_frequency = 0.25

best_phase, best_rms = 0, np.inf
# Brute force
for phase in np.linspace(-np.pi, np.pi, 1000):
    theo = nominal_amplitude * np.sin(2 * np.pi * nominal_frequency * time + phase)
    rms = np.sqrt(np.mean((filtered - theo)**2))
    if rms < best_rms:
        best_rms, best_phase = rms, phase

theoretical = nominal_amplitude * np.sin(2 * np.pi * nominal_frequency * time + best_phase)

# 5. Errors Calculation
error_sig = filtered - theoretical
measured_amplitude = (np.percentile(filtered, 99) - np.percentile(filtered, 1)) / 2
error_amp_pct = (abs(measured_amplitude - nominal_amplitude) / nominal_amplitude) * 100

# 6. Figures
fig, axes = plt.subplots(3, 1, figsize=(14, 12))

# Subplot 1: Raw vs Filtered
axes[0].plot(time, centered, alpha=0.3, color='gray', label='Raw Data')
axes[0].plot(time, filtered, color='#1f77b4', label='Filtered (Butterworth 2Hz)')
axes[0].set_title('Signal Preprocessing (Raw vs Filtered)')
axes[0].legend()
axes[0].grid(True)

# Subplot 2: Comparison with Theoretical
axes[1].plot(time, filtered, color='#1f77b4', label='Measured')
axes[1].plot(time, theoretical, color='red', linestyle='--', label=f'Theoretical (Phase={best_phase:.2f} rad)')
axes[1].set_title(f'Empirical vs Theoretical Trajectory: {error_amp_pct:.2f}%')
axes[1].legend()
axes[1].grid(True)

# Subplot 3: Residuals
axes[2].fill_between(time, error_sig, 0, color='green', alpha=0.2)
axes[2].plot(time, error_sig, color='green', linewidth=0.8)
axes[2].set_title(f'Point-by-Point Residual Error (RMS = {best_rms:.2f} mm)')
axes[2].set_xlabel('Time (s)')
plt.tight_layout()
plt.show()

# 7.
peaks_err, _ = find_peaks(filtered, prominence=2.0, distance=int(fs*3))
errors_per_cycle = [np.sqrt(np.mean(error_sig[peaks_err[i]:peaks_err[i+1]]**2)) for i in range(len(peaks_err)-1)]
# Check (like POLIZZI) 
peak_times = time[peaks_err] 

measured_periods = np.diff(peak_times)
mean_period = np.mean(measured_periods)
std_period = np.std(measured_periods)

print(f"\nPeriod Check (TOLERANCE ±0.2s)")
print(f"Nominal Period: 4.00 s")
print(f"Measured Period:  {mean_period:.3f} s (± {std_period:.3f} s)")
print(f"\nSin Statistics (SIN):")
print(f"Nominal Amplitude: {nominal_amplitude} mm")
print(f"Measured Amplitude:  {measured_amplitude:.2f} mm")
print(f"Errors per cycle (RMS):   {np.mean(errors_per_cycle):.2f} mm")