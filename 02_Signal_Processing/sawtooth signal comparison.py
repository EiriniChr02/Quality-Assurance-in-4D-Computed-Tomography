import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, butter, filtfilt
from scipy import signal as scipy_signal

df = pd.read_csv('waveform_export_SAW.csv')
time       = df['time_s'].values
raw_signal = df['Ant.-Post. [millimeter]'].values

midpoint = (np.max(raw_signal) + np.min(raw_signal)) / 2
centered = -(raw_signal - midpoint)

fs   = 1 / (time[1] - time[0])
b, a = butter(4, 2.0 / (fs/2), btype='low')
filtered = filtfilt(b, a, centered)

# Parameters
nominal_amplitude = 10.0
nominal_frequency = 0.25
from scipy.signal import correlate

# BRUTE FORCE:
best_phase, best_width, best_rms = 0, 0.85, np.inf

for width in np.linspace(0.70, 0.95, 25):
    for phase in np.linspace(0, 2*np.pi, 500):
        saw_raw      = scipy_signal.sawtooth(2*np.pi*nominal_frequency*time + phase, width=width)
        saw_mid      = (np.max(saw_raw) + np.min(saw_raw)) / 2
        saw_centered = saw_raw - saw_mid

        pos_max = np.max(saw_centered)
        neg_min = np.min(saw_centered)
        theo = np.where(
            saw_centered >= 0,
            nominal_amplitude * saw_centered / pos_max,
            nominal_amplitude * saw_centered / abs(neg_min)
        )

        rms = np.sqrt(np.mean((filtered - theo)**2))
        if rms < best_rms:
            best_rms, best_phase, best_width = rms, phase, width

print(f"Best width={best_width:.3f}")
print(f"Best phase={best_phase:.4f} rad")
print(f"Best RMS={best_rms:.2f} mm")

# Fianl Theoretical
saw_raw      = scipy_signal.sawtooth(2*np.pi*nominal_frequency*time + best_phase, width=best_width)
saw_mid      = (np.max(saw_raw) + np.min(saw_raw)) / 2
saw_centered = saw_raw - saw_mid
pos_max      = np.max(saw_centered)
neg_min      = np.min(saw_centered)
theoretical  = np.where(
    saw_centered >= 0,
    nominal_amplitude * saw_centered / pos_max,
    nominal_amplitude * saw_centered / abs(neg_min)
)
# Measurements
measured_amplitude = (np.max(filtered) - np.min(filtered)) / 2
error_amp  = (abs(measured_amplitude - nominal_amplitude) / nominal_amplitude) * 100
error_sig  = filtered - theoretical
rms_error  = np.sqrt(np.mean(error_sig**2))
# Plots
fig, axes = plt.subplots(3, 1, figsize=(14, 12))

axes[0].plot(time, centered, alpha=0.3, color='gray', linewidth=0.8, label='Raw')
axes[0].plot(time, filtered, color='#1f77b4', linewidth=1.5, label='Filtered (2Hz)')
axes[0].axhline( nominal_amplitude, color='orange', linewidth=0.8, linestyle='--', label=f'±{nominal_amplitude}mm')
axes[0].axhline(-nominal_amplitude, color='orange', linewidth=0.8, linestyle='--')
axes[0].axhline(0, color='gray', linewidth=0.8, linestyle=':')
axes[0].set_title('Nominal Signal: Raw vs Filtered')
axes[0].set_ylabel('Amplitude (mm)')
axes[0].set_ylim(-16, 16)
axes[0].legend()
axes[0].grid(True)

axes[1].plot(time, filtered,    color='#1f77b4', linewidth=1.5, label='Nominal (filtered)')
axes[1].plot(time, theoretical, color='red', linewidth=2.0, linestyle='--',
             label=f'Theoretical Sawtooth (A=±{nominal_amplitude}mm, T=4s)')
axes[1].axhline( nominal_amplitude, color='orange', linewidth=0.8, linestyle='--')
axes[1].axhline(-nominal_amplitude, color='orange', linewidth=0.8, linestyle='--')
axes[1].axhline(0, color='gray', linewidth=0.8, linestyle=':')
axes[1].set_title('Empirical vs Theoretical Trajectory')
axes[1].set_ylabel('Amplitude (mm)')
axes[1].set_ylim(-16, 16)
axes[1].legend()
axes[1].grid(True)

axes[2].plot(time, error_sig, color='green', linewidth=1.0, label='Σφάλμα')
axes[2].axhline(0, color='gray', linewidth=0.8, linestyle=':')
axes[2].fill_between(time, error_sig, 0, alpha=0.2, color='green')
axes[2].set_title(f'Error (Nominal − Theoretical)  |  RMS={rms_error:.2f} mm')
axes[2].set_xlabel('Time (s)')
axes[2].set_ylabel('Difference (mm)')
axes[2].grid(True)
axes[2].legend()

plt.tight_layout()
plt.show()

print(f"\nNominal Amplitude:  {nominal_amplitude:.2f} mm")
print(f"Measured Amplitude:   {measured_amplitude:.2f} mm")
print(f"Amplitude Error:     {error_amp:.2f}%")
print(f"RMS Error:      {rms_error:.2f} mm")

peaks_err, _ = find_peaks(filtered, prominence=2.0, distance=int(fs*3))

errors_per_cycle = []
for i in range(len(peaks_err) - 1):
    start = peaks_err[i]
    end   = peaks_err[i+1]
    rms_i = np.sqrt(np.mean(error_sig[start:end]**2))
    errors_per_cycle.append(rms_i)

cycle_numbers = np.arange(1, len(errors_per_cycle) + 1)

fig, axes = plt.subplots(2, 1, figsize=(12, 8))

# Plot RMS/cycle
axes[0].plot(cycle_numbers, errors_per_cycle,
             marker='o', color='steelblue', linewidth=1.5, markersize=5)
axes[0].axhline(np.mean(errors_per_cycle), color='red', linestyle='--',
                linewidth=1.5, label=f'Μέσος RMS = {np.mean(errors_per_cycle):.2f} mm')
axes[0].axhline(np.mean(errors_per_cycle) + np.std(errors_per_cycle),
                color='orange', linestyle=':', label=f'±1σ = {np.std(errors_per_cycle):.2f} mm')
axes[0].axhline(np.mean(errors_per_cycle) - np.std(errors_per_cycle),
                color='orange', linestyle=':')
axes[0].set_title('RMS error per breathing cycle')
axes[0].set_xlabel('Beathing Cycle (#)')
axes[0].set_ylabel('RMS Error (mm)')
axes[0].legend()
axes[0].grid(True)

# Histogram
axes[1].hist(errors_per_cycle, bins=10, color='steelblue',
             edgecolor='white', alpha=0.8)
axes[1].axvline(np.mean(errors_per_cycle), color='red', linestyle='--',
                linewidth=1.5, label=f'Mean = {np.mean(errors_per_cycle):.2f} mm')
axes[1].axvline(np.max(errors_per_cycle), color='orange', linestyle=':',
                linewidth=1.5, label=f'Max = {np.max(errors_per_cycle):.2f} mm')
axes[1].set_title('RMS Error/cycle')
axes[1].set_xlabel('RMS Error (mm)')
axes[1].set_ylabel('Number of Cycles')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()

print(f"\nAnalysis per cycle ({len(errors_per_cycle)} κύκλοι):")
print(f"Mean RMS:   {np.mean(errors_per_cycle):.2f} mm")
print(f"Std RMS:     {np.std(errors_per_cycle):.2f} mm")
print(f"Min RMS:     {np.min(errors_per_cycle):.2f} mm")
print(f"Max RMS:     {np.max(errors_per_cycle):.2f} mm")

# Check of real period (like POLIZZI)
peak_times = time[peaks_err] 
measured_periods = np.diff(peak_times)
mean_period = np.mean(measured_periods)
std_period = np.std(measured_periods)

print(f"\nPeriod Check(TOLERANCE ±0.2s) ---")
print(f"Nominal Period: 4.00 s")
print(f"Measured Period:  {mean_period:.3f} s (± {std_period:.3f} s)")