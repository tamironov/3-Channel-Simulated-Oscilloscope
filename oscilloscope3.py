import customtkinter as ctk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import threading
import numpy as np
from scipy.signal import sawtooth, square
import csv
from datetime import datetime

# Set a modern dark theme for customtkinter
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

class OscilloscopeApp(ctk.CTk):
    """
    An advanced Python application for a simulated multi-channel oscilloscope.
    This version features a fully optimized and responsive GUI layout using grid.
    The layout has been improved to ensure all controls and data fit on the screen.
    """
    def __init__(self):
        super().__init__()

        # --- Window Configuration ---
        self.title("Multi-Channel Simulated Oscilloscope")
        self.geometry("1600x900")
        # Configure the main window grid: Control panel on the left, graph on the right
        self.grid_columnconfigure(0, weight=0)  # Control frame has fixed width
        self.grid_columnconfigure(1, weight=1)  # Graph frame takes all remaining space
        self.grid_rowconfigure(0, weight=1)

        # --- Frames ---
        # Main control frame on the left side
        self.control_frame = ctk.CTkFrame(self, width=400, corner_radius=10)
        self.control_frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
        self.control_frame.grid_columnconfigure(0, weight=1)
        self.control_frame.grid_rowconfigure(3, weight=1) # Allow the tabview to expand

        # Main graph and measurement frame on the right side
        self.graph_frame = ctk.CTkFrame(self, corner_radius=10)
        self.graph_frame.grid(row=0, column=1, padx=(0, 20), pady=20, sticky="nsew")
        self.graph_frame.grid_columnconfigure(0, weight=1)
        self.graph_frame.grid_rowconfigure(0, weight=1)  # Graph takes most of the space
        self.graph_frame.grid_rowconfigure(1, weight=0)  # Measurements section does not expand

        # --- Data containers and threading controls ---
        self.running = False
        self.channels = {}  # Dictionary to hold signal data for each channel
        self.channel_count = 0
        self.selected_channel_id = None
        self.channel_colors = ['#FFFF00', '#00FFFF', '#FF00FF'] # Bright, distinct colors

        self.data_lock = threading.Lock()
        self.after_id = None

        self.setup_ui()

    def setup_ui(self):
        """
        Builds all the UI widgets and lays them out using a grid system.
        The layout is now more robust and responsive.
        """
        # --- Control Frame Header and Separator ---
        ctk.CTkLabel(self.control_frame, text="Simulated Oscilloscope", font=ctk.CTkFont(size=20, weight="bold")).grid(row=0, column=0, padx=20, pady=(20, 5), sticky="ew")
        
        # --- Channel Management Section ---
        channel_mgt_section = ctk.CTkFrame(self.control_frame, fg_color="transparent")
        channel_mgt_section.grid(row=1, column=0, padx=20, pady=(10, 5), sticky="ew")
        channel_mgt_section.grid_columnconfigure(0, weight=1)
        channel_mgt_section.grid_columnconfigure(1, weight=1)
        
        self.add_channel_button = ctk.CTkButton(channel_mgt_section, text="Add Channel", command=self.add_channel)
        self.add_channel_button.grid(row=0, column=0, padx=5, pady=5, sticky="ew")
        
        self.remove_channel_button = ctk.CTkButton(channel_mgt_section, text="Remove Channel", command=self.remove_channel, state="disabled", fg_color="#E74C3C", hover_color="#C0392B")
        self.remove_channel_button.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        # --- Oscilloscope Controls Section ---
        scope_controls_section = ctk.CTkFrame(self.control_frame, fg_color="transparent")
        scope_controls_section.grid(row=2, column=0, padx=20, pady=(10, 5), sticky="ew")
        scope_controls_section.grid_columnconfigure(0, weight=1)
        scope_controls_section.grid_columnconfigure(1, weight=1)
        
        ctk.CTkLabel(scope_controls_section, text="Oscilloscope Controls", font=ctk.CTkFont(size=16, weight="bold")).grid(row=0, column=0, columnspan=2, padx=5, pady=(0, 10))
        
        ctk.CTkLabel(scope_controls_section, text="Time/Div (ms):", anchor="w").grid(row=1, column=0, padx=5, pady=(5, 2), sticky="w")
        self.time_div_entry = ctk.CTkEntry(scope_controls_section, placeholder_text="e.g., 2")
        self.time_div_entry.insert(0, "2")
        self.time_div_entry.grid(row=1, column=1, padx=5, pady=(5, 2), sticky="ew")
        
        ctk.CTkLabel(scope_controls_section, text="Volts/Div (V):", anchor="w").grid(row=2, column=0, padx=5, pady=(5, 2), sticky="w")
        self.volts_div_entry = ctk.CTkEntry(scope_controls_section, placeholder_text="e.g., 1")
        self.volts_div_entry.insert(0, "1")
        self.volts_div_entry.grid(row=2, column=1, padx=5, pady=(5, 2), sticky="ew")

        # --- Signal Generator Controls (Tab View for Channels) ---
        # This will hold the settings for each individual channel
        self.signal_tabview = ctk.CTkTabview(self.control_frame, command=self.on_tab_change)
        self.signal_tabview.grid(row=3, column=0, padx=20, pady=10, sticky="nsew")

        # --- Main Action Buttons and Status ---
        action_frame = ctk.CTkFrame(self.control_frame, fg_color="transparent")
        action_frame.grid(row=4, column=0, padx=20, pady=(5, 20), sticky="ew")
        action_frame.grid_columnconfigure(0, weight=1)
        action_frame.grid_columnconfigure(1, weight=1)
        
        self.start_continue_button = ctk.CTkButton(action_frame, text="Start Oscilloscope", command=self.handle_start_pause, font=ctk.CTkFont(size=14, weight="bold"), height=40)
        self.start_continue_button.grid(row=0, column=0, columnspan=2, padx=5, pady=(0, 5), sticky="ew")
        
        self.reset_button = ctk.CTkButton(action_frame, text="Reset", command=self.reset_test, state="disabled", fg_color="#3a7e9e", hover_color="#306985")
        self.reset_button.grid(row=1, column=0, padx=5, pady=5, sticky="ew")

        self.export_button = ctk.CTkButton(action_frame, text="Save to CSV", command=self.export_csv_as, state="disabled")
        self.export_button.grid(row=1, column=1, padx=5, pady=5, sticky="ew")

        self.status_label = ctk.CTkLabel(action_frame, text="Status: Idle", font=ctk.CTkFont(size=14, weight="bold"), anchor="center")
        self.status_label.grid(row=2, column=0, columnspan=2, pady=(5, 0), sticky="ew")
        
        # --- Graph Frame Widgets ---
        # A small, fixed-size figure for the graph
        self.fig, self.ax = plt.subplots(figsize=(6, 4), dpi=100)
        self.configure_plot_style()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.graph_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        # A new frame below the graph to display live measurements
        self.measurement_display_frame = ctk.CTkFrame(self.graph_frame, corner_radius=8, fg_color="#363636")
        self.measurement_display_frame.grid(row=1, column=0, padx=20, pady=(0, 20), sticky="ew")
        self.measurement_display_frame.grid_columnconfigure(0, weight=1)
        self.measurement_display_frame.grid_columnconfigure(1, weight=1)
        self.measurement_display_frame.grid_columnconfigure(2, weight=1)
        
        # Measurement labels now in a grid for better alignment
        ctk.CTkLabel(self.measurement_display_frame, text="Measurements (Selected Channel)", font=ctk.CTkFont(size=16, weight="bold")).grid(row=0, column=0, columnspan=3, pady=(10, 10), sticky="ew")
        
        self.freq_label = ctk.CTkLabel(self.measurement_display_frame, text="Frequency: --- Hz", font=ctk.CTkFont(size=14), anchor="w")
        self.vpp_label = ctk.CTkLabel(self.measurement_display_frame, text="V_pp: --- V", font=ctk.CTkFont(size=14), anchor="w")
        self.vrms_label = ctk.CTkLabel(self.measurement_display_frame, text="V_rms: --- V", font=ctk.CTkFont(size=14), anchor="w")
        self.vmax_label = ctk.CTkLabel(self.measurement_display_frame, text="V_max: --- V", font=ctk.CTkFont(size=14), anchor="w")
        self.vmin_label = ctk.CTkLabel(self.measurement_display_frame, text="V_min: --- V", font=ctk.CTkFont(size=14), anchor="w")
        self.duty_cycle_label_display = ctk.CTkLabel(self.measurement_display_frame, text="Duty Cycle: --- %", font=ctk.CTkFont(size=14), anchor="w")
        
        self.freq_label.grid(row=1, column=0, padx=(20, 5), pady=2, sticky="w")
        self.vpp_label.grid(row=2, column=0, padx=(20, 5), pady=2, sticky="w")
        self.vrms_label.grid(row=3, column=0, padx=(20, 5), pady=2, sticky="w")
        self.vmax_label.grid(row=1, column=1, padx=(5, 20), pady=2, sticky="w")
        self.vmin_label.grid(row=2, column=1, padx=(5, 20), pady=2, sticky="w")
        self.duty_cycle_label_display.grid(row=3, column=1, padx=(5, 20), pady=2, sticky="w")
        
        self.add_channel() # Start with one channel by default

    def configure_plot_style(self):
        """Sets the visual style for the matplotlib plot."""
        self.fig.set_facecolor("#2b2b2b")
        self.ax.set_facecolor("#000000") # Oscilloscope-like black background
        
        self.ax.set_title("Simulated Waveform", color='white', fontname="Consolas")
        self.ax.set_xlabel("Time (s)", color='white', fontname="Consolas")
        self.ax.set_ylabel("Voltage (V)", color='white', fontname="Consolas")
        
        self.ax.tick_params(axis='x', colors='green') # Oscilloscope-like grid colors
        self.ax.tick_params(axis='y', colors='green')
        self.ax.spines['left'].set_color('green')
        self.ax.spines['bottom'].set_color('green')
        self.ax.spines['right'].set_color('green')
        self.ax.spines['top'].set_color('green')
        
        # Add finer grid lines to emulate an oscilloscope's graticule
        self.ax.grid(True, which='major', linestyle='-', alpha=0.5, color='green')
        self.ax.grid(True, which='minor', linestyle=':', alpha=0.2, color='green')
        self.ax.minorticks_on()

    def add_channel(self):
        """Adds a new signal channel with default settings."""
        if self.channel_count >= 3:
            messagebox.showinfo("Limit Reached", "You can only have up to 3 channels.")
            return

        self.channel_count += 1
        channel_id = self.channel_count
        
        # Default channel parameters
        self.channels[channel_id] = {
            'waveform_type': ctk.StringVar(value="Sine"),
            'gen_freq': ctk.StringVar(value="100"),
            'gen_amp': ctk.StringVar(value="5"),
            'duty_cycle': ctk.StringVar(value="50"),
            'offset': ctk.StringVar(value="0"),
            'phase': ctk.StringVar(value="0"),
            'time_data': np.array([]),
            'signal_data': np.array([]),
            'color': self.channel_colors[len(self.channels) - 1]
        }
        
        # Create a new tab for the channel and load its controls
        tab_name = f"Channel {channel_id}"
        self.signal_tabview.add(tab_name)
        self.load_channel_controls(channel_id)
        
        # Select the new tab
        self.signal_tabview.set(tab_name)
        self.selected_channel_id = channel_id
        
        if self.running:
            self.generate_signals()
        
        self.remove_channel_button.configure(state="normal")

    def remove_channel(self):
        """Removes the currently selected channel."""
        if self.selected_channel_id:
            tab_name = f"Channel {self.selected_channel_id}"
            
            # Remove from data and GUI
            del self.channels[self.selected_channel_id]
            self.signal_tabview.delete(tab_name)
            
            self.selected_channel_id = None
            if self.channels:
                # Find the next available channel and select it
                next_channel_id = sorted(self.channels.keys())[0]
                self.signal_tabview.set(f"Channel {next_channel_id}")
                self.selected_channel_id = next_channel_id
            else:
                self.remove_channel_button.configure(state="disabled")
            
            self.channel_count -= 1
            if self.channel_count == 0:
                self.reset_data_and_ui()

            self.update_plots()
            self.calculate_and_display_measurements()
        
        if self.running and self.channels:
            self.generate_signals()

    def on_tab_change(self, tab_name):
        """Callback for when a tab is changed. Updates the selected channel and measurements."""
        if tab_name:
            self.selected_channel_id = int(tab_name.split()[-1])
            self.calculate_and_display_measurements()

    def load_channel_controls(self, channel_id):
        """Populates a new tab with the control widgets for a given channel."""
        tab_frame = self.signal_tabview.tab(f"Channel {channel_id}")
        tab_frame.grid_columnconfigure(1, weight=1)
        channel = self.channels[channel_id]
        
        # Waveform Type
        ctk.CTkLabel(tab_frame, text="Waveform Type:", anchor="w").grid(row=0, column=0, padx=5, pady=2, sticky="w")
        option_menu = ctk.CTkOptionMenu(tab_frame, variable=channel['waveform_type'], values=["Sine", "Square", "Sawtooth", "Triangle", "White Noise"], command=lambda v: self.on_waveform_change(v, channel_id))
        option_menu.grid(row=0, column=1, padx=5, pady=2, sticky="ew")

        # Frequency
        ctk.CTkLabel(tab_frame, text="Signal Freq (Hz):", anchor="w").grid(row=1, column=0, padx=5, pady=2, sticky="w")
        ctk.CTkEntry(tab_frame, textvariable=channel['gen_freq']).grid(row=1, column=1, padx=5, pady=2, sticky="ew")

        # Amplitude
        ctk.CTkLabel(tab_frame, text="Signal Amplitude (V):", anchor="w").grid(row=2, column=0, padx=5, pady=2, sticky="w")
        ctk.CTkEntry(tab_frame, textvariable=channel['gen_amp']).grid(row=2, column=1, padx=5, pady=2, sticky="ew")

        # Vertical Offset
        ctk.CTkLabel(tab_frame, text="Vertical Offset (V):", anchor="w").grid(row=3, column=0, padx=5, pady=2, sticky="w")
        ctk.CTkEntry(tab_frame, textvariable=channel['offset']).grid(row=3, column=1, padx=5, pady=2, sticky="ew")
        
        # Phase Offset
        ctk.CTkLabel(tab_frame, text="Phase Offset (deg):", anchor="w").grid(row=4, column=0, padx=5, pady=2, sticky="w")
        ctk.CTkEntry(tab_frame, textvariable=channel['phase']).grid(row=4, column=1, padx=5, pady=2, sticky="ew")

        # Duty Cycle (Square wave only)
        channel['duty_cycle_label_gen'] = ctk.CTkLabel(tab_frame, text="Duty Cycle (%):", anchor="w")
        channel['duty_cycle_entry_gen'] = ctk.CTkEntry(tab_frame, textvariable=channel['duty_cycle'])
        
        # Update button for manual signal regeneration
        self.update_button = ctk.CTkButton(tab_frame, text="Update Signal", command=self.generate_signals)
        self.update_button.grid(row=6, column=0, columnspan=2, padx=5, pady=(10, 5), sticky="ew")

        self.on_waveform_change(channel['waveform_type'].get(), channel_id)

    def on_waveform_change(self, value, channel_id):
        """Shows/hides the duty cycle entry based on the selected waveform."""
        channel = self.channels[channel_id]
        if value == "Square":
            channel['duty_cycle_label_gen'].grid(row=5, column=0, padx=5, pady=2, sticky="w")
            channel['duty_cycle_entry_gen'].grid(row=5, column=1, padx=5, pady=2, sticky="ew")
            self.update_button.grid(row=6, column=0, columnspan=2, padx=5, pady=(10, 5), sticky="ew")
        else:
            channel['duty_cycle_label_gen'].grid_remove()
            channel['duty_cycle_entry_gen'].grid_remove()
            self.update_button.grid(row=5, column=0, columnspan=2, padx=5, pady=(10, 5), sticky="ew")

    def generate_signals(self):
        """
        Simulates signals for all active channels based on their parameters.
        """
        try:
            time_div = float(self.time_div_entry.get())
            x_divs = 10 
            total_time_s = time_div * x_divs / 1000.0
            
            with self.data_lock:
                for ch_id, ch_data in self.channels.items():
                    gen_freq = float(ch_data['gen_freq'].get())
                    gen_amp = float(ch_data['gen_amp'].get())
                    waveform_type = ch_data['waveform_type'].get()
                    offset = float(ch_data['offset'].get())
                    phase_deg = float(ch_data['phase'].get())
                    
                    duty_cycle = 0.5
                    if waveform_type == "Square":
                        duty_cycle = float(ch_data['duty_cycle'].get()) / 100.0
                        if not 0 < duty_cycle <= 1:
                            raise ValueError(f"Duty cycle for Channel {ch_id} must be between 1 and 100.")

                    sample_rate = max(gen_freq * 1000, 5000)
                    num_samples = int(total_time_s * sample_rate)
                    t = np.linspace(0, total_time_s, num_samples, endpoint=False)
                    phase_rad = np.deg2rad(phase_deg)
                    
                    if waveform_type == "Sine":
                        signal = gen_amp * np.sin(2 * np.pi * gen_freq * t + phase_rad) + offset
                    elif waveform_type == "Square":
                        signal = gen_amp * square(2 * np.pi * gen_freq * t + phase_rad, duty=duty_cycle) + offset
                    elif waveform_type == "Sawtooth":
                        signal = gen_amp * sawtooth(2 * np.pi * gen_freq * t + phase_rad) + offset
                    elif waveform_type == "Triangle":
                        signal = gen_amp * sawtooth(2 * np.pi * gen_freq * t + phase_rad, width=0.5) + offset
                    elif waveform_type == "White Noise":
                        signal = np.random.normal(0, gen_amp * 0.5, num_samples) + offset
                    
                    # Add some minor Gaussian noise to all signals
                    signal += np.random.normal(0, gen_amp * 0.01, num_samples)
                    
                    ch_data['time_data'] = t
                    ch_data['signal_data'] = signal
            
            self.status_label.configure(text="Status: Signals updated.")
            self.update_plots()
            self.calculate_and_display_measurements()
            
        except ValueError as e:
            messagebox.showerror("Input Error", f"Please enter valid numeric values: {e}")

    def handle_start_pause(self):
        """Handles the logic for the start/pause button."""
        if self.running:
            self.pause_test()
        else:
            self.start_test()

    def start_test(self):
        """Starts the oscilloscope measurement loop."""
        if not self.channels:
            messagebox.showerror("No Channels", "Please add at least one channel to start the oscilloscope.")
            return

        try:
            float(self.time_div_entry.get())
            float(self.volts_div_entry.get())
            
            self.running = True
            self.start_continue_button.configure(text="Pause Oscilloscope", fg_color="#E74C3C", hover_color="#C0392B")
            self.reset_button.configure(state="normal")
            self.export_button.configure(state="disabled")
            
            self.status_label.configure(text="Status: Oscilloscope running.")
            
            self.generate_signals()
            self.update_gui()

        except ValueError as e:
            messagebox.showerror("Input Error", f"Please enter valid numeric values for the oscilloscope controls: {e}")
            self.status_label.configure(text="Status: Idle (Error)")

    def pause_test(self):
        """Pauses the measurement and updates the UI."""
        self.running = False
        
        if self.after_id:
            self.after_cancel(self.after_id)
            self.after_id = None
        
        self.status_label.configure(text="Status: Paused")
        self.start_continue_button.configure(text="Continue Oscilloscope", fg_color="#1f6aa5", hover_color="#144876")
        self.reset_button.configure(state="normal")
        self.export_button.configure(state="normal")

    def reset_test(self):
        """Resets all data and UI elements to their initial state."""
        self.running = False
        
        if self.after_id:
            self.after_cancel(self.after_id)
            self.after_id = None
        
        self.reset_data_and_ui()
        self.status_label.configure(text="Status: Idle")

        self.start_continue_button.configure(text="Start Oscilloscope", fg_color="#1f6aa5", hover_color="#144876")
        self.reset_button.configure(state="disabled")
        self.export_button.configure(state="disabled")

    def reset_data_and_ui(self):
        """Clears data arrays and resets the graph and labels."""
        self.channels = {}
        self.channel_count = 0
        self.selected_channel_id = None
        
        # Clear all tabs from the tabview
        for tab_name in self.signal_tabview.tabnames:
            self.signal_tabview.delete(tab_name)
        
        # Reset measurement labels
        self.freq_label.configure(text="Frequency: --- Hz", text_color="gray")
        self.vpp_label.configure(text="V_pp: --- V", text_color="gray")
        self.vrms_label.configure(text="V_rms: --- V", text_color="gray")
        self.vmax_label.configure(text="V_max: --- V", text_color="gray")
        self.vmin_label.configure(text="V_min: --- V", text_color="gray")
        self.duty_cycle_label_display.grid_remove() # Use grid_remove instead of pack_forget
        
        self.ax.clear()
        self.configure_plot_style()
        self.canvas.draw()
        
    def update_gui(self):
        """
        The main GUI update loop, which is called periodically.
        It updates the plots and live value labels.
        """
        if self.running:
            self.generate_signals()
            self.after_id = self.after(500, self.update_gui)
        else:
            self.status_label.configure(text="Status: Paused")
            if self.channels:
                self.export_button.configure(state="normal")
                self.reset_button.configure(state="normal")
            self.after_id = None

    def update_plots(self):
        """Clears and redraws the oscilloscope plot for all channels."""
        self.ax.clear()
        
        try:
            time_div_ms = float(self.time_div_entry.get())
            volts_div = float(self.volts_div_entry.get())
            
            x_divs = 10
            y_divs = 8
            
            x_lim = (0, time_div_ms * x_divs / 1000.0)
            y_lim = (-volts_div * y_divs/2, volts_div * y_divs/2)
            
            self.ax.set_xlim(x_lim)
            self.ax.set_ylim(y_lim)

            self.ax.axhline(0, color='gray', linewidth=0.8, linestyle='--')
            
        except ValueError:
            pass
        
        # Plot all channels
        for ch_id, ch_data in self.channels.items():
            if ch_data['signal_data'].size > 0:
                self.ax.plot(ch_data['time_data'], ch_data['signal_data'], color=ch_data['color'], linewidth=2, label=f'Channel {ch_id}')

        # The legend is now placed outside the plot to avoid overlapping with data.
        self.ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.15), frameon=False, prop={'size': 12}, labelcolor='white')
        
        self.configure_plot_style()
        self.fig.tight_layout()
        self.canvas.draw()
    
    def calculate_and_display_measurements(self):
        """Calculates and updates the live measurement labels for the selected channel."""
        if not self.selected_channel_id or not self.channels[self.selected_channel_id]['signal_data'].size > 0:
            self.freq_label.configure(text="Frequency: --- Hz", text_color="gray")
            self.vpp_label.configure(text="V_pp: --- V", text_color="gray")
            self.vrms_label.configure(text="V_rms: --- V", text_color="gray")
            self.vmax_label.configure(text="V_max: --- V", text_color="gray")
            self.vmin_label.configure(text="V_min: --- V", text_color="gray")
            self.duty_cycle_label_display.grid_remove()
            return

        channel_data = self.channels[self.selected_channel_id]
        signal = channel_data['signal_data']
        color = channel_data['color']
        
        # Frequency
        try:
            gen_freq = float(channel_data['gen_freq'].get())
            self.freq_label.configure(text=f"Frequency: {gen_freq:.2f} Hz", text_color=color)
        except ValueError:
            self.freq_label.configure(text="Frequency: --- Hz", text_color="gray")

        # Peak-to-Peak Voltage
        v_pp = np.max(signal) - np.min(signal)
        self.vpp_label.configure(text=f"V_pp: {v_pp:.2f} V", text_color=color)
        
        # RMS Voltage
        v_rms = np.sqrt(np.mean(signal**2))
        self.vrms_label.configure(text=f"V_rms: {v_rms:.2f} V", text_color=color)
        
        # Min/Max Voltage
        v_max = np.max(signal)
        v_min = np.min(signal)
        self.vmax_label.configure(text=f"V_max: {v_max:.2f} V", text_color=color)
        self.vmin_label.configure(text=f"V_min: {v_min:.2f} V", text_color=color)
        
        # Duty Cycle (for square waves only)
        if channel_data['waveform_type'].get() == "Square":
            self.duty_cycle_label_display.grid(row=3, column=1, padx=(5, 20), pady=2, sticky="w")
            try:
                duty_cycle = float(channel_data['duty_cycle'].get())
                self.duty_cycle_label_display.configure(text=f"Duty Cycle: {duty_cycle:.2f} %", text_color=color)
            except ValueError:
                self.duty_cycle_label_display.configure(text="Duty Cycle: --- %", text_color="gray")
        else:
            self.duty_cycle_label_display.grid_remove()

    def export_csv_as(self):
        """Saves the current data from all channels to a CSV file."""
        if not self.channels:
            self.status_label.configure(text="Status: No data to export.")
            return

        filetypes = [('CSV files', '*.csv'), ('All files', '*.*')]
        filename = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=filetypes,
                                                initialfile=f"oscilloscope_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
        if not filename:
            self.status_label.configure(text="Status: Save cancelled.")
            return

        try:
            with open(filename, mode='w', newline='') as file:
                writer = csv.writer(file)
                
                # Write header row
                header = []
                for ch_id in sorted(self.channels.keys()):
                    header.append(f"Channel {ch_id} Time (s)")
                    header.append(f"Channel {ch_id} Voltage (V)")
                writer.writerow(header)

                # Get all signal data and find the maximum length
                all_signals_data = {ch_id: self.channels[ch_id]['signal_data'] for ch_id in self.channels}
                all_time_data = {ch_id: self.channels[ch_id]['time_data'] for ch_id in self.channels}
                max_len = max(len(s) for s in all_signals_data.values()) if all_signals_data else 0

                # Write data rows
                for i in range(max_len):
                    row = []
                    for ch_id in sorted(self.channels.keys()):
                        time_val = f"{all_time_data[ch_id][i]:.6f}" if i < len(all_time_data[ch_id]) else ""
                        volt_val = f"{all_signals_data[ch_id][i]:.4f}" if i < len(all_signals_data[ch_id]) else ""
                        row.extend([time_val, volt_val])
                    writer.writerow(row)

            self.status_label.configure(text=f"Status: Data successfully saved to {filename}")
        except Exception as e:
            messagebox.showerror("CSV Export Error", f"An error occurred while saving the file: {e}")
            self.status_label.configure(text="Status: CSV Export Error")

if __name__ == "__main__":
    app = OscilloscopeApp()
    app.mainloop()
