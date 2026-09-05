import pydicom
import numpy as np
import pandas as pd
from pathlib import Path 

# ---------------------------
# CONFIG
# ---------------------------
PROJECT_DIR = Path(__file__).parent.resolve()

input_filename = "SIN.dcm"
dicom_path = PROJECT_DIR / input_filename

output_csv = PROJECT_DIR / "waveform_export_SIN.csv"
output_csv_raw = PROJECT_DIR / "waveform_export_SIN_raw.csv"

if not dicom_path.exists():
    raise FileNotFoundError(f"File not found: {dicom_path}")

print(f"Work File: {PROJECT_DIR}")
print(f"Initial File: {dicom_path.name}")


import pydicom
import numpy as np
import pandas as pd
from pathlib import Path

# ---------------------------
# CONFIG
# ---------------------------
dicom_path = r"C:\Users\Name\Downloads\file.csv"
output_csv = r"C:\Users\Name\Downloads\file.csv"
output_csv_raw = r"C:\Users\name\Downloads\file.csv"

# ---------------------------
# LOAD DICOM
# ---------------------------
ds = pydicom.dcmread(dicom_path)

if "WaveformSequence" not in ds:
    raise ValueError("No Waveform Sequence found in this DICOM file.")

wf = ds.WaveformSequence[0]

n_channels = int(wf.NumberOfWaveformChannels)
n_samples = int(wf.NumberOfWaveformSamples)
fs = float(wf.SamplingFrequency)
sample_interp = wf.WaveformSampleInterpretation
bits_allocated = int(wf.WaveformBitsAllocated)

print(f"Channels: {n_channels}")
print(f"Samples: {n_samples}")
print(f"Sampling frequency: {fs} Hz")
print(f"Sample interpretation: {sample_interp}")
print(f"Bits allocated: {bits_allocated}")

# ---------------------------
# MAP DICOM waveform datatype
# ---------------------------
dtype_map = {
    "SB": np.int8,
    "UB": np.uint8,
    "MB": np.uint8,   # not common
    "AB": np.uint8,   # not common
    "SS": np.int16,
    "US": np.uint16,
    "SL": np.int32,
    "UL": np.uint32,
}

if sample_interp not in dtype_map:
    raise ValueError(f"Unsupported WaveformSampleInterpretation: {sample_interp}")

dtype = dtype_map[sample_interp]

raw = np.frombuffer(wf.WaveformData, dtype="<" + np.dtype(dtype).str[1:])

expected_count = n_samples * n_channels
if raw.size != expected_count:
    raise ValueError(
        f"Unexpected waveform length: got {raw.size}, expected {expected_count}"
    )

# Interleaved layout: [s0c0, s0c1, ..., s0c7, s1c0, ...]
raw = raw.reshape(n_samples, n_channels)

# ---------------------------
# EXTRACT CHANNEL METADATA
# ---------------------------
channel_defs = wf.ChannelDefinitionSequence

channel_labels = []
sensitivities = []
correction_factors = []
baselines = []
units = []

for ch in channel_defs:
    label = str(ch.get("ChannelLabel", f"Ch{len(channel_labels)+1}"))
    sens = float(ch.get("ChannelSensitivity", 1.0))
    corr = float(ch.get("ChannelSensitivityCorrectionFactor", 1.0))
    base = float(ch.get("ChannelBaseline", 0.0))

    unit = ""
    if "ChannelSensitivityUnitsSequence" in ch and len(ch.ChannelSensitivityUnitsSequence) > 0:
        unit = str(ch.ChannelSensitivityUnitsSequence[0].get("CodeMeaning", ""))

    channel_labels.append(label)
    sensitivities.append(sens)
    correction_factors.append(corr)
    baselines.append(base)
    units.append(unit)

sensitivities = np.array(sensitivities, dtype=float)
correction_factors = np.array(correction_factors, dtype=float)
baselines = np.array(baselines, dtype=float)


# DICOM waveform formula:
# physical = raw * sensitivity * correction + baseline

scaled = raw.astype(np.float64) * sensitivities * correction_factors + baselines

# Time axis
time_s = np.arange(n_samples) / fs

# ---------------------------
# EXPORT RAW CSV
# ---------------------------
raw_columns = ["time_s"] + [f"{lbl}_raw" for lbl in channel_labels]
raw_df = pd.DataFrame(np.column_stack([time_s, raw]), columns=raw_columns)
raw_df.to_csv(output_csv_raw, index=False, encoding="utf-8")

# ---------------------------
# EXPORT SCALED CSV
# ---------------------------
scaled_colnames = []
for lbl, unit in zip(channel_labels, units):
    if unit:
        scaled_colnames.append(f"{lbl} [{unit}]")
    else:
        scaled_colnames.append(lbl)

scaled_df = pd.DataFrame(np.column_stack([time_s, scaled]), columns=["time_s"] + scaled_colnames)
scaled_df.to_csv(output_csv, index=False, encoding="utf-8")

print(f"Saved raw waveform CSV to: {output_csv_raw}")
print(f"Saved scaled waveform CSV to: {output_csv}") 