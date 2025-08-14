# Multi-Channel Simulated Oscilloscope (Tkinter + Matplotlib)
3 Channel Simulated Oscilloscope Python

A desktop GUI oscilloscope simulation built with CustomTkinter and Matplotlib.
Create up to 3 channels, choose waveform types (Sine, Square, Sawtooth, Triangle, White Noise), set frequency, amplitude, duty cycle, phase, and offsets, then visualize signals with an oscilloscope-style grid. Live measurements (Freq, Vpp, Vrms, Vmax, Vmin, Duty Cycle) update per selected channel. Data can be exported to CSV.

**Features:**
- Up to 3 channels with per-channel settings & tabs
- Waveforms: Sine, Square (with Duty Cycle), Sawtooth, Triangle, White Noise
- Time/Div and Volts/Div controls with oscilloscope-style graticule
- Live measurements per selected channel: Frequency, Vpp, Vrms, Vmax, Vmin, Duty Cycle
- Start / Pause / Reset workflow
- CSV export: time & voltage columns for each channel
- Modern dark theme using CustomTkinter

**Requirements**
- Python 3.8+
- OS: Windows / macOS / Linux
- GUI backend: Tkinter (comes with most Python installers; on Linux you might need to install it)

**Python packges**
- pip install customtkinter matplotlib numpy scipy

**How to Use**
1) Add/Remove Channel
- Click Add Channel to create a channel tab (max 3).
- Use Remove Channel to delete the currently selected tab.
2) Set Oscilloscope Controls
- Time/Div (ms): horizontal scale across 10 divisions
- Volts/Div (V): vertical scale across 8 divisions
3) Configure Channel (per tab)
- Waveform Type: Sine, Square, Sawtooth, Triangle, White Noise
- Signal Freq (Hz)
- Signal Amplitude (V)
- Vertical Offset (V)
- Phase Offset (deg)
- Duty Cycle (%): shown only for Square
- Click Update Signal to regenerate waveforms immediately.
4) Run / Pause / Reset
- Start Oscilloscope begins periodic updates (~500 ms).
- Pause Oscilloscope stops updates but keeps plots.
- Reset clears channels, plots, and measurements.
5)Export to CSV
- Click Save to CSV (enabled when paused or stopped).
- File contains two columns per channel: Channel N Time (s), Channel N Voltage (V).
