"""
MultiGaussFit v1.0
Signal preprocessing and multi-Gaussian deconvolution GUI.

Author: Juan Emanuel López Cervantes
https://github.com/1JELC1/MultiGaussFit
"""


# Imports

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import pandas as pd
import os
import sys
import traceback
import threading

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.optimize import differential_evolution
from scipy.fft import fft, ifft, fftfreq
from scipy.signal import savgol_filter


# Theme colors

BG_DARK       = "#1e1e1e"
BG_SIDEBAR    = "#2b2b2b"
BG_WIDGET     = "#3c3c3c"
BG_ENTRY      = "#444444"
FG_TEXT       = "#e0e0e0"
FG_DIM        = "#888888"
ACCENT        = "#4CAF50"
ACCENT_HOVER  = "#66BB6A"
ACCENT_BLUE   = "#42A5F5"
ACCENT_RED    = "#ef5350"
ACCENT_ORANGE = "#FFA726"
BORDER        = "#555555"
PLOT_BG       = "#1e1e2e"
PLOT_FG       = "#cccccc"
PLOT_GRID     = "#444466"

GAUSS_COLORS = [
    '#AB47BC', '#26A69A', '#FFCA28', '#EF5350', '#5C6BC0',
    '#8D6E63', '#78909C', '#66BB6A', '#FF7043', '#29B6F6',
]


# Signal processing functions

def gaussian(x, a, mu, sigma):
    """Evaluate a single Gaussian"""
    with np.errstate(divide='ignore', invalid='ignore'):
        return a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))


def multi_gaussian(x, params, n_centers, num_gaussians):
    """Sum of n_centers Gaussians (1 or 2 per center)"""
    y_fit = np.zeros_like(x, dtype=float)
    block = 3 if num_gaussians == 1 else 5
    for i in range(n_centers):
        p = params[i * block:(i + 1) * block]
        if num_gaussians == 1:
            y_fit += gaussian(x, p[0], p[1], p[2])
        else:
            y_fit += gaussian(x, p[0], p[2], p[3]) + gaussian(x, p[1], p[2], p[4])
    return y_fit


def filter_noise_ft(x, y, n_frequencies):
    """Filter signal keeping only the top n_frequencies Fourier components"""
    n = len(x)
    fft_y = fft(y)
    indices = np.argsort(np.abs(fft_y))[::-1]
    mask = np.zeros(n, dtype=bool)
    mask[indices[:n_frequencies]] = True
    return np.real(ifft(fft_y * mask))


def filter_noise_sg(x, y, window_length=11, polyorder=3):
    """Filter signal using Savitzky-Golay"""
    if window_length % 2 == 0:
        window_length += 1
    if window_length > len(y):
        window_length = len(y) - 1 if len(y) % 2 == 0 else len(y)
    if polyorder >= window_length:
        polyorder = window_length - 1
    return savgol_filter(y, window_length, polyorder)


def linear_baseline_correction(x, y):
    """Remove linear baseline by connecting minima in left and right halves"""
    N = len(y)
    i1 = np.argmin(y[:N // 2])
    i2 = np.argmin(y[N // 2:]) + N // 2
    x1, y1 = x[i1], y[i1]
    x2, y2 = x[i2], y[i2]
    if x1 == x2:
        return y - np.min(y)
    m = (y2 - y1) / (x2 - x1)
    c = y1 - m * x1
    y_corr = y - (m * x + c)
    y_corr -= np.min(y_corr)
    return y_corr


def shift_minimum(y):
    """Shift signal so that minimum is zero"""
    return y - np.min(y)


def remove_negatives(y):
    """Clip negative values to zero"""
    return np.clip(y, a_min=0, a_max=None)


def peak_summary(block_params):
    """Extract amplitude, center, width from a parameter block"""
    if len(block_params) == 3:
        return block_params[0], block_params[1], block_params[2]
    elif len(block_params) == 5:
        a1, a2, mu, s1, s2 = block_params
        return a1 + a2, mu, max(s1, s2)
    raise ValueError(f"Invalid parameter block length: {len(block_params)}")


def auto_detect_centers(x, y, intervals):
    """Find x position of the max y value within each interval"""
    centers = []
    for xl, xr in intervals:
        mask = (x >= xl) & (x <= xr)
        if np.any(mask):
            centers.append(x[mask][np.argmax(y[mask])])
        else:
            centers.append((xl + xr) / 2)
    return centers


def optimize_gaussians(x, y, centers, intervals, num_gaussians=1,
                       max_iter=1000, tolerance=1e-8):
    """Optimize Gaussian parameters via differential evolution"""
    n_centers = len(centers)
    ymax = float(np.max(y))
    bounds = []
    for xl, xr in intervals:
        width = xr - xl
        if num_gaussians == 1:
            bounds.extend([(0, ymax), (xl, xr), (0, width)])
        else:
            bounds.extend([
                (0, ymax), (0, ymax),
                (xl, xr),
                (0, width), (0, width),
            ])

    def objective(params):
        y_fit = multi_gaussian(x, params, n_centers, num_gaussians)
        err = 0.0
        for xl, xr in intervals:
            m = (x >= xl) & (x <= xr)
            err += np.sum((y[m] - y_fit[m]) ** 2)
        return err

    result = differential_evolution(
        func=objective, bounds=bounds,
        maxiter=max_iter, tol=tolerance,
        mutation=(0.5, 1.5), recombination=0.7,
        updating='deferred', workers=1,
    )
    return result.x, result.fun


def is_number_str(s):
    """Check if a string represents a number"""
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


# GUI application

class MultiGaussFitApp(tk.Tk):
    """Main application window"""

    def __init__(self):
        super().__init__()
        self.title("MultiGaussFit v1.0")
        self.geometry("1400x850")
        self.minsize(1100, 650)
        self.configure(bg=BG_DARK)

        # Start maximized on Windows
        try:
            self.state('zoomed')
        except tk.TclError:
            pass
            
        try:
            import ctypes
            myappid = 'jelc.multigaussfit.gui.1.0'
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
        except Exception:
            pass

        try:
            img = tk.PhotoImage(file=self._resource_path("logo.png"))
            self.iconphoto(True, img)
        except Exception:
            pass


        # Application state
        self.filepath = None
        self.df = None
        self.pairs = []           # [(x, y, name), ...]
        self.current_idx = None
        self.x = None
        self.y = None
        self.y_original = None
        self.undo_stack = []
        self.pipeline = []        # preprocessing steps for batch replay

        # Deconvolution state
        self.intervals = []       # [[xl, xr], ...]
        self.deconv_params = None
        self.deconv_n_gauss = None
        self.all_results = {}     # {name: (params, n_gauss, x, y, y_orig)}

        # Click mode
        self.click_mode = False
        self.first_click_x = None

        # Batch processing flag
        self._batch_running = False

        # Tk variables
        self.struct_var = tk.StringVar(value='XYYY')
        self.ngauss_var = tk.StringVar(value='1')
        self.is_light_mode = False

        # Build UI
        self._configure_styles()
        self._build_status_bar()
        self._build_sidebar()
        self._build_plot_area()

    def _resource_path(self, relative_path):
        """ Get absolute path to resource"""
        try:
            base_path = sys._MEIPASS
        except Exception:
            base_path = os.path.abspath(".")
        return os.path.join(base_path, relative_path)


    # STYLES AND THEMES

    def _toggle_theme(self):
        self.is_light_mode = getattr(self, 'is_light_mode', False)
        self.is_light_mode = not self.is_light_mode
        self._apply_theme()

    def _apply_theme(self):
        global BG_DARK, BG_SIDEBAR, BG_WIDGET, BG_ENTRY, FG_TEXT, FG_DIM, PLOT_BG, PLOT_FG, PLOT_GRID, BORDER
        if self.is_light_mode:
            BG_DARK       = "#f5f5f5"
            BG_SIDEBAR    = "#e0e0e0"
            BG_WIDGET     = "#c8c8c8"
            BG_ENTRY      = "#ffffff"
            FG_TEXT       = "#1e1e1e"
            FG_DIM        = "#666666"
            PLOT_BG       = "#ffffff"
            PLOT_FG       = "#1e1e1e"
            PLOT_GRID     = "#e0e0e0"
            BORDER        = "#aaaaaa"
        else:
            BG_DARK       = "#1e1e1e"
            BG_SIDEBAR    = "#2b2b2b"
            BG_WIDGET     = "#3c3c3c"
            BG_ENTRY      = "#444444"
            FG_TEXT       = "#e0e0e0"
            FG_DIM        = "#888888"
            PLOT_BG       = "#1e1e2e"
            PLOT_FG       = "#cccccc"
            PLOT_GRID     = "#444466"
            BORDER        = "#555555"

        self.configure(bg=BG_DARK)
        if hasattr(self, 'container'):
            self.container.configure(bg=BG_SIDEBAR)
            self._sb_canvas.configure(bg=BG_SIDEBAR)
            self.sb_inner.configure(style='TFrame') # Forces ttk frame to update if needed
        if hasattr(self, 'status_bar'):
            self.status_bar.configure(bg='#dddddd' if self.is_light_mode else '#1a1a1a', fg=FG_DIM)
        if hasattr(self, 'signal_list'):
            self.signal_list.configure(bg=BG_ENTRY, fg=FG_TEXT, highlightbackground=BORDER)
        if hasattr(self, 'results_text'):
            self.results_text.configure(bg=BG_ENTRY, fg=FG_TEXT, highlightbackground=BORDER, highlightcolor=BORDER)
        if hasattr(self, 'plot_frame'):
            self.plot_frame.configure(bg=BG_DARK)
            self.toolbar_frame.configure(bg=BG_DARK)
        if hasattr(self, 'lbl_filepath'):
            self.lbl_filepath.configure(foreground=FG_DIM)

        self._configure_styles()
        if hasattr(self, 'interval_frame'):
            self._refresh_intervals_ui()
        if hasattr(self, 'ax'):
            self._style_axes()
            if not self._batch_running:
                self._redraw()
            else:
                self.canvas_mpl.draw()

    def _configure_styles(self):
        s = ttk.Style()
        s.theme_use('clam')

        s.configure('.', background=BG_SIDEBAR, foreground=FG_TEXT,
                    font=('Segoe UI', 10))
        s.configure('TFrame', background=BG_SIDEBAR)
        s.configure('TLabel', background=BG_SIDEBAR, foreground=FG_TEXT)

        s.configure('TLabelframe', background=BG_SIDEBAR, foreground=ACCENT,
                    font=('Segoe UI', 10, 'bold'), borderwidth=1, relief='groove')
        s.configure('TLabelframe.Label', background=BG_SIDEBAR, foreground=ACCENT,
                    font=('Segoe UI', 10, 'bold'))

        # Standard button
        s.configure('TButton', background=BG_WIDGET, foreground=FG_TEXT,
                    borderwidth=0, focusthickness=0, padding=(8, 5))
        s.map('TButton',
              background=[('active', '#505050'), ('disabled', '#333333')],
              foreground=[('disabled', FG_DIM)])

        # Green accent button
        s.configure('Accent.TButton', background=ACCENT, foreground='white',
                    font=('Segoe UI', 10, 'bold'), padding=(10, 7))
        s.map('Accent.TButton',
              background=[('active', ACCENT_HOVER), ('disabled', '#2E7D32')],
              foreground=[('disabled', '#aaaaaa')])

        # Orange click-mode button
        s.configure('ClickMode.TButton', background=ACCENT_ORANGE, foreground='white',
                    font=('Segoe UI', 10, 'bold'), padding=(8, 5))
        s.map('ClickMode.TButton', background=[('active', '#FFB74D')])

        # Radio buttons
        s.configure('TRadiobutton', background=BG_SIDEBAR, foreground=FG_TEXT)
        s.map('TRadiobutton', background=[('active', BG_SIDEBAR)])

        # Combobox
        s.configure('TCombobox', fieldbackground=BG_ENTRY,
                    background=BG_WIDGET, foreground=FG_TEXT, arrowcolor=FG_TEXT)
        s.map('TCombobox',
              fieldbackground=[('readonly', BG_ENTRY)],
              foreground=[('readonly', FG_TEXT)])

        # Separator
        s.configure('TSeparator', background=BORDER)


    # status bar

    def _build_status_bar(self):
        self.status_bar = tk.Label(
            self, text="  Welcome to MultiGaussFit v1.0  —  Open a CSV file to begin",
            bg='#1a1a1a', fg=FG_DIM, anchor=tk.W, padx=12,
            font=('Segoe UI', 9))
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def set_status(self, text):
        self.status_bar.configure(text=f"  {text}")
        self.update_idletasks()


    # sidebar

    def _build_sidebar(self):
        self.container = tk.Frame(self, bg=BG_SIDEBAR, width=290)
        self.container.pack(side=tk.LEFT, fill=tk.Y)
        self.container.pack_propagate(False)

        # Scrollable canvas for sidebar content
        self._sb_canvas = tk.Canvas(self.container, bg=BG_SIDEBAR,
                                    highlightthickness=0, width=290)
        self.sb_inner = ttk.Frame(self._sb_canvas)
        self.sb_inner.bind('<Configure>',
                      lambda e: self._sb_canvas.configure(
                          scrollregion=self._sb_canvas.bbox('all')))
        self._sb_canvas.create_window((0, 0), window=self.sb_inner,
                                      anchor='nw', width=280)
        self._sb_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Mousewheel scrolling (only when hovering sidebar)
        self._sb_canvas.bind('<Enter>', self._on_sidebar_enter)
        self._sb_canvas.bind('<Leave>', self._on_sidebar_leave)

        PADX = 6

        # settings
        sec = ttk.LabelFrame(self.sb_inner, text="  ⚙️  Settings  ", padding=8)
        sec.pack(fill=tk.X, padx=PADX, pady=(8, 3))
        ttk.Button(sec, text="🌓 Toggle Light/Dark Mode",
                   command=self._toggle_theme).pack(fill=tk.X, pady=2)

        # file
        sec = ttk.LabelFrame(self.sb_inner, text="  📂  File  ", padding=8)
        sec.pack(fill=tk.X, padx=PADX, pady=3)

        ttk.Button(sec, text="Open CSV…",
                   command=self._open_file).pack(fill=tk.X, pady=(0, 4))

        self.lbl_filepath = ttk.Label(sec, text="No file loaded",
                                       foreground=FG_DIM, wraplength=245,
                                       font=('Segoe UI', 8))
        self.lbl_filepath.pack(fill=tk.X, pady=(0, 6))

        row = ttk.Frame(sec)
        row.pack(fill=tk.X)
        ttk.Label(row, text="Structure:").pack(side=tk.LEFT)
        ttk.Radiobutton(row, text="XYYY", variable=self.struct_var,
                        value='XYYY').pack(side=tk.LEFT, padx=(10, 4))
        ttk.Radiobutton(row, text="XYXY", variable=self.struct_var,
                        value='XYXY').pack(side=tk.LEFT)

        ttk.Button(sec, text="Load Data",
                   command=self._load_data).pack(fill=tk.X, pady=(6, 0))

        # signals
        sec = ttk.LabelFrame(self.sb_inner, text="  📊  Signals  ", padding=8)
        sec.pack(fill=tk.X, padx=PADX, pady=3)

        self.signal_list = tk.Listbox(
            sec, bg=BG_ENTRY, fg=FG_TEXT, selectbackground=ACCENT,
            selectforeground='white', font=('Segoe UI', 10), height=4,
            borderwidth=0, highlightthickness=1,
            highlightcolor=ACCENT, highlightbackground=BORDER,
            activestyle='none')
        self.signal_list.pack(fill=tk.X)
        self.signal_list.bind('<<ListboxSelect>>', self._on_signal_select)
        # Isolate listbox scroll from sidebar scroll
        self.signal_list.bind('<Enter>', self._on_listbox_enter)
        self.signal_list.bind('<Leave>', self._on_listbox_leave)

        # preprocessing
        sec = ttk.LabelFrame(self.sb_inner, text="  🔧  Preprocessing  ", padding=8)
        sec.pack(fill=tk.X, padx=PADX, pady=3)

        for text, cmd in [
            ("Shift Minimum",       self._preproc_shift),
            ("Baseline Correction", self._preproc_baseline),
            ("Remove Negatives",    self._preproc_clip_neg),
        ]:
            ttk.Button(sec, text=text, command=cmd).pack(fill=tk.X, pady=2)

        ttk.Separator(sec).pack(fill=tk.X, pady=5)

        for text, cmd in [
            ("Savitzky-Golay Filter",    self._dialog_sg),
            ("Fourier Transform Filter", self._dialog_ft),
        ]:
            ttk.Button(sec, text=text, command=cmd).pack(fill=tk.X, pady=2)

        ttk.Separator(sec).pack(fill=tk.X, pady=5)

        row = ttk.Frame(sec)
        row.pack(fill=tk.X)
        ttk.Button(row, text="↩ Undo",
                   command=self._undo).pack(
            side=tk.LEFT, expand=True, fill=tk.X, padx=(0, 2))
        ttk.Button(row, text="Reset All",
                   command=self._reset_all).pack(
            side=tk.RIGHT, expand=True, fill=tk.X, padx=(2, 0))

        # deconvolution
        sec = ttk.LabelFrame(self.sb_inner, text="  📈  Deconvolution  ", padding=8)
        sec.pack(fill=tk.X, padx=PADX, pady=(3, 8))

        row = ttk.Frame(sec)
        row.pack(fill=tk.X, pady=(0, 6))
        ttk.Label(row, text="Gaussians / peak:").pack(side=tk.LEFT)
        ttk.Combobox(row, textvariable=self.ngauss_var,
                     values=['1', '2'], width=3,
                     state='readonly').pack(side=tk.RIGHT)

        self.btn_add_interval = ttk.Button(
            sec, text="+ Add Interval  (click on plot)",
            command=self._toggle_click_mode)
        self.btn_add_interval.pack(fill=tk.X, pady=(0, 4))

        self.interval_frame = ttk.Frame(sec)
        self.interval_frame.pack(fill=tk.X, pady=(0, 4))

        ttk.Button(sec, text="Clear All Intervals",
                   command=self._clear_intervals).pack(fill=tk.X, pady=(0, 8))

        self.btn_run = ttk.Button(
            sec, text="▶  Run Deconvolution",
            style='Accent.TButton', command=self._run_deconvolution)
        self.btn_run.pack(fill=tk.X, pady=(0, 4))

        self.btn_batch = ttk.Button(
            sec, text="▶  Deconvolve All (same config)",
            command=self._run_batch_deconvolution)
        self.btn_batch.pack(fill=tk.X, pady=(0, 6))

        # Results display
        self.results_text = tk.Text(
            sec, bg=BG_ENTRY, fg=FG_TEXT, font=('Consolas', 9),
            height=5, borderwidth=0, highlightthickness=1,
            highlightcolor=BORDER, highlightbackground=BORDER,
            state='disabled', wrap='word')
        self.results_text.pack(fill=tk.X, pady=(0, 6))

        ttk.Button(sec, text="💾  Save Results…",
                   command=self._save_results).pack(fill=tk.X)

        # Isolate results_text scroll from sidebar scroll
        self.results_text.bind('<Enter>', self._on_results_enter)
        self.results_text.bind('<Leave>', self._on_results_leave)

    def _on_sidebar_enter(self, event):
        self.bind_all('<MouseWheel>', self._sb_scroll)

    def _on_sidebar_leave(self, event):
        self.unbind_all('<MouseWheel>')

    def _on_listbox_enter(self, event):
        self.bind_all('<MouseWheel>', self._listbox_scroll)

    def _on_listbox_leave(self, event):
        self.bind_all('<MouseWheel>', self._sb_scroll)

    def _on_results_enter(self, event):
        self.bind_all('<MouseWheel>', self._results_scroll)

    def _on_results_leave(self, event):
        self.bind_all('<MouseWheel>', self._sb_scroll)

    def _sb_scroll(self, event):
        self._sb_canvas.yview_scroll(int(-1 * (event.delta / 120)), 'units')

    def _listbox_scroll(self, event):
        self.signal_list.yview_scroll(int(-1 * (event.delta / 120)), 'units')
        return 'break'

    def _results_scroll(self, event):
        self.results_text.yview_scroll(int(-1 * (event.delta / 120)), 'units')
        return 'break'

    # plot area

    def _build_plot_area(self):
        self.plot_frame = tk.Frame(self, bg=BG_DARK)
        self.plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.fig = Figure(dpi=100, facecolor=BG_DARK)
        self.ax = self.fig.add_subplot(111)

        self.canvas_mpl = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas_mpl.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Navigation toolbar (zoom, pan, save)
        self.toolbar_frame = tk.Frame(self.plot_frame, bg=BG_DARK)
        self.toolbar_frame.pack(fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas_mpl, self.toolbar_frame)
        self.toolbar.update()

        # Connect click handler
        self.canvas_mpl.mpl_connect('button_press_event', self._on_plot_click)

        self._style_axes()
        self.ax.set_title('Load a CSV file to begin', fontsize=12, pad=10)
        self.canvas_mpl.draw()

    def _style_axes(self):
        ax = self.ax
        self.fig.patch.set_facecolor(BG_DARK)
        ax.set_facecolor(PLOT_BG)
        ax.tick_params(colors=PLOT_FG, which='both')
        ax.xaxis.label.set_color(PLOT_FG)
        ax.yaxis.label.set_color(PLOT_FG)
        ax.title.set_color(FG_TEXT)
        for spine in ax.spines.values():
            spine.set_color(BORDER)
        ax.grid(True, alpha=0.15, color=PLOT_GRID)
        ax.set_xlabel('Position')
        ax.set_ylabel('Intensity')

    # file operations

    def _open_file(self):
        fp = filedialog.askopenfilename(
            title="Select CSV file",
            filetypes=[("CSV files", "*.csv"),
                       ("Text files", "*.txt"),
                       ("All files", "*.*")])
        if fp:
            self.filepath = fp
            self.lbl_filepath.configure(text=os.path.basename(fp),
                                         foreground=FG_TEXT)
            self.set_status(f"File selected: {os.path.basename(fp)}  —  "
                           f"Choose structure and click 'Load Data'")

    def _load_data(self):
        if not self.filepath:
            messagebox.showwarning("No File", "Please open a CSV file first.")
            return
        try:
            df = pd.read_csv(self.filepath, encoding='latin-1')
            all_num = all(is_number_str(c) for c in df.columns)
            if all_num:
                df = pd.read_csv(self.filepath, encoding='latin-1', header=None)
            df = df.dropna(axis=1, how='all')
            if df.shape[1] < 2:
                messagebox.showerror("Error",
                                     "At least 2 columns with data are required.")
                return
            self.df = df
            self._parse_pairs()
            self._populate_signals()
            self.set_status(
                f"Loaded {len(self.pairs)} signal(s) from "
                f"{os.path.basename(self.filepath)}")
        except Exception as e:
            messagebox.showerror("Load Error", str(e))

    def _parse_pairs(self):
        df = self.df
        cols = list(df.columns)

        # Auto-rename purely-numeric headers
        num_count = sum(is_number_str(c) for c in cols)
        if num_count > len(cols) / 2:
            n = len(cols)
            if self.struct_var.get() == 'XYYY':
                df.columns = (['Position'] +
                              [f"Signal {i}" for i in range(1, n)])
            else:
                p = n // 2
                nc = []
                for i in range(1, p + 1):
                    nc += [f"Position {i}", f"Signal {i}"]
                df.columns = nc[:n]

        self.pairs = []
        if self.struct_var.get() == 'XYYY':
            x = df.iloc[:, 0].astype(float).values
            for j in range(1, df.shape[1]):
                y = df.iloc[:, j].astype(float).values
                self.pairs.append((x.copy(), y.copy(), str(df.columns[j])))
        else:
            for i in range(df.shape[1] // 2):
                try:
                    xi = df.iloc[:, 2 * i].astype(float).values
                    yi = df.iloc[:, 2 * i + 1].astype(float).values
                    self.pairs.append(
                        (xi.copy(), yi.copy(), str(df.columns[2 * i + 1])))
                except Exception:
                    pass

    def _populate_signals(self):
        self.signal_list.delete(0, tk.END)
        for _, _, name in self.pairs:
            self.signal_list.insert(tk.END, name)
        if self.pairs:
            self.signal_list.selection_set(0)
            self._on_signal_select(None)


    # signal selection

    def _on_signal_select(self, event):
        sel = self.signal_list.curselection()
        if not sel:
            return
        idx = sel[0]
        if idx == self.current_idx:
            return

        self.current_idx = idx
        x, y, name = self.pairs[idx]
        self.x = x.copy()
        self.y = y.copy()
        self.y_original = y.copy()
        self.undo_stack.clear()
        self.pipeline.clear()
        self.intervals.clear()
        self.deconv_params = None
        self.deconv_n_gauss = None
        self.click_mode = False
        self.first_click_x = None
        self.btn_add_interval.configure(
            text="+ Add Interval  (click on plot)", style='TButton')
        self._refresh_intervals_ui()
        self._clear_results_text()
        self._redraw()
        self.set_status(f"Signal: {name}")

    # plot drawing

    def _redraw(self):
        self.ax.clear()
        self._style_axes()

        if self.x is None:
            self.ax.set_title('Load a CSV file to begin', fontsize=12, pad=10)
            self.canvas_mpl.draw()
            return

        name = self.pairs[self.current_idx][2]

        # Original signal (faded background when preprocessed)
        has_preproc = len(self.undo_stack) > 0
        if has_preproc:
            self.ax.plot(self.x, self.y_original, color='#888888',
                         linewidth=1.0, label='Original', alpha=0.35,
                         zorder=2, linestyle='-')

        # Current (preprocessed) signal
        sig_label = 'Preprocessed' if has_preproc else 'Signal'
        self.ax.plot(self.x, self.y, color=ACCENT_BLUE, linewidth=1.4,
                     label=sig_label, alpha=0.85, zorder=3)

        # Intervals (yellow shading)
        for i, (xl, xr) in enumerate(self.intervals):
            self.ax.axvspan(xl, xr, alpha=0.10, color='#FFD54F', zorder=1)
            self.ax.axvline(xl, color='#FFD54F', lw=0.7, ls='--', alpha=0.5)
            self.ax.axvline(xr, color='#FFD54F', lw=0.7, ls='--', alpha=0.5)
            # Number label at top
            ylim = self.ax.get_ylim()
            self.ax.text((xl + xr) / 2, ylim[1] * 0.97, str(i + 1),
                        ha='center', va='top', fontsize=9, color='#FFD54F',
                        fontweight='bold', alpha=0.7, zorder=11)

        # Deconvolution results
        title_extra = ''
        if self.deconv_params is not None:
            ng = self.deconv_n_gauss
            nc = len(self.intervals)
            block = 3 if ng == 1 else 5

            y_fit = multi_gaussian(self.x, self.deconv_params, nc, ng)
            sse = np.sum((self.y - y_fit) ** 2)

            # Total fit line
            self.ax.plot(self.x, y_fit, color='#FF7043', lw=2,
                        ls='--', label='Total fit', zorder=5)

            # Individual Gaussians
            for i in range(nc):
                ci = i % len(GAUSS_COLORS)
                if ng == 1:
                    a, mu, sig = self.deconv_params[3*i:3*i+3]
                    yg = gaussian(self.x, a, mu, sig)
                    self.ax.plot(self.x, yg, color=GAUSS_COLORS[ci], lw=1.2,
                                ls=':', label=f'G{i+1}', zorder=4)
                    self.ax.fill_between(self.x, 0, yg,
                                         color=GAUSS_COLORS[ci], alpha=0.07)
                else:
                    a1, a2, mu, s1, s2 = self.deconv_params[5*i:5*i+5]
                    y1 = gaussian(self.x, a1, mu, s1)
                    y2 = gaussian(self.x, a2, mu, s2)
                    self.ax.plot(self.x, y1, color=GAUSS_COLORS[ci], lw=1.2,
                                ls=':', label=f'P{i+1}-G1', zorder=4)
                    self.ax.plot(self.x, y2, color=GAUSS_COLORS[ci], lw=1.2,
                                ls='-.', label=f'P{i+1}-G2', zorder=4)
                    self.ax.fill_between(self.x, 0, y1 + y2,
                                         color=GAUSS_COLORS[ci], alpha=0.05)

            title_extra = f'  |  SSE = {sse:.2e}'
            self.set_status(
                f"SSE = {sse:.2e}  |  "
                f"{nc} peak(s) fitted with {ng} Gaussian(s) each")

        # Temporary click line
        if self.first_click_x is not None:
            self.ax.axvline(self.first_click_x, color=ACCENT_ORANGE,
                           lw=1.5, ls='--', zorder=10)

        # Title and legend
        self.ax.set_title(f'{name}{title_extra}', fontsize=11, pad=10)
        self.ax.legend(loc='upper right', fontsize=8, framealpha=0.7,
                       facecolor=BG_SIDEBAR, edgecolor=BORDER,
                       labelcolor=FG_TEXT)
        self.fig.tight_layout()
        self.canvas_mpl.draw()

    # preprocessing

    def _push_undo(self):
        # save current y for undo
        self.undo_stack.append(self.y.copy())
        self.deconv_params = None
        self.deconv_n_gauss = None
        self._clear_results_text()
        self.all_results.clear()

    def _preproc_shift(self):
        if self.y is None:
            return
        self._push_undo()
        self.y = shift_minimum(self.y)
        self.pipeline.append(('shift', {}))
        self._redraw()
        self.set_status("Applied: Shift Minimum")

    def _preproc_baseline(self):
        if self.y is None:
            return
        self._push_undo()
        self.y = linear_baseline_correction(self.x, self.y)
        self.pipeline.append(('baseline', {}))
        self._redraw()
        self.set_status("Applied: Baseline Correction")

    def _preproc_clip_neg(self):
        if self.y is None:
            return
        self._push_undo()
        self.y = remove_negatives(self.y)
        self.pipeline.append(('clip_neg', {}))
        self._redraw()
        self.set_status("Applied: Remove Negatives")

    # Filter dialogs

    def _dialog_sg(self):
        if self.y is None:
            return
        self._show_param_dialog(
            title="Savitzky-Golay Filter",
            fields=[("Window Length (odd, 5–25):", "11"),
                    ("Polynomial Order (2–4):", "3")],
            callback=self._apply_sg)

    def _apply_sg(self, values):
        wl, po = int(values[0]), int(values[1])
        if wl < 3 or po < 1:
            messagebox.showwarning("Invalid", "Window ≥ 3 and order ≥ 1.")
            return False
        self._push_undo()
        self.y = filter_noise_sg(self.x, self.y, wl, po)
        if np.any(~np.isfinite(self.y)):
            self.y = np.nan_to_num(self.y, nan=0.0, posinf=0.0, neginf=0.0)
        self.pipeline.append(('sg', {'wl': wl, 'po': po}))
        self._redraw()
        self.set_status(f"Applied: SG Filter (window={wl}, order={po})")
        return True

    def _dialog_ft(self):
        if self.y is None:
            return
        self._show_param_dialog(
            title="Fourier Transform Filter",
            fields=[("Top frequencies to keep:", "10")],
            callback=self._apply_ft)

    def _apply_ft(self, values):
        nf = int(values[0])
        if nf <= 0:
            messagebox.showwarning("Invalid", "Must be a positive integer.")
            return False
        self._push_undo()
        self.y = filter_noise_ft(self.x, self.y, nf)
        if np.any(~np.isfinite(self.y)):
            self.y = np.nan_to_num(self.y, nan=0.0, posinf=0.0, neginf=0.0)
        self.pipeline.append(('ft', {'nf': nf}))
        self._redraw()
        self.set_status(f"Applied: FT Filter (n_freq={nf})")
        return True

    def _show_param_dialog(self, title, fields, callback):
        # show a small param dialog
        dlg = tk.Toplevel(self)
        dlg.title(title)
        dlg.configure(bg=BG_SIDEBAR)
        dlg.resizable(False, False)
        dlg.transient(self)
        dlg.grab_set()

        content = ttk.Frame(dlg, padding=20)
        content.pack(fill=tk.BOTH, expand=True)

        entries = []
        for label_text, default in fields:
            ttk.Label(content, text=label_text).pack(anchor=tk.W, pady=(0, 2))
            var = tk.StringVar(value=default)
            e = tk.Entry(content, textvariable=var, bg=BG_ENTRY, fg=FG_TEXT,
                         insertbackground=FG_TEXT, font=('Segoe UI', 10),
                         borderwidth=0, highlightthickness=1,
                         highlightcolor=ACCENT, highlightbackground=BORDER)
            e.pack(fill=tk.X, pady=(0, 10), ipady=4)
            entries.append(var)

        def on_apply():
            try:
                vals = [v.get() for v in entries]
                if callback(vals):
                    dlg.destroy()
            except ValueError:
                messagebox.showwarning("Invalid Input",
                                       "Please enter valid numbers.",
                                       parent=dlg)

        ttk.Button(content, text="Apply", style='Accent.TButton',
                   command=on_apply).pack(fill=tk.X, pady=(4, 0))

        # Center the dialog on the parent window
        dlg.update_idletasks()
        w = dlg.winfo_reqwidth()
        h = dlg.winfo_reqheight()
        x = self.winfo_x() + (self.winfo_width() - w) // 2
        y = self.winfo_y() + (self.winfo_height() - h) // 2
        dlg.geometry(f"+{x}+{y}")

    # Undo / Reset

    def _undo(self):
        if not self.undo_stack:
            self.set_status("Nothing to undo")
            return
        self.y = self.undo_stack.pop()
        if self.pipeline:
            self.pipeline.pop()
        self.deconv_params = None
        self.deconv_n_gauss = None
        self.all_results.clear()
        self._clear_results_text()
        self._redraw()
        self.set_status(
            f"Undo applied  ({len(self.undo_stack)} step(s) remaining)")

    def _reset_all(self):
        if self.y_original is None:
            return
        self.y = self.y_original.copy()
        self.undo_stack.clear()
        self.pipeline.clear()
        self.deconv_params = None
        self.deconv_n_gauss = None
        self.all_results.clear()
        self._clear_results_text()
        self._redraw()
        self.set_status("Reset to original signal")

    # interval management

    def _toggle_click_mode(self):
        if self.x is None:
            messagebox.showinfo("No Data", "Load a signal first.")
            return

        self.click_mode = not self.click_mode
        if self.click_mode:
            self.btn_add_interval.configure(
                text="🎯  Click 2 points on the plot…",
                style='ClickMode.TButton')
            self.set_status(
                "CLICK MODE  —  click two points on the plot "
                "to define a deconvolution interval")
        else:
            self.btn_add_interval.configure(
                text="+ Add Interval  (click on plot)", style='TButton')
            self.first_click_x = None
            self._redraw()
            self.set_status("Click mode cancelled")

    def _on_plot_click(self, event):
        # Ignore clicks outside axes, in toolbar mode, or when click mode is off
        if (not self.click_mode
                or event.inaxes != self.ax
                or event.xdata is None):
            return

        # Ignore if toolbar is in zoom/pan mode
        if self.toolbar.mode:
            return

        if self.first_click_x is None:
            # First click: draw temporary vertical line
            self.first_click_x = event.xdata
            self._redraw()
            self.set_status(
                f"First point: {event.xdata:.2f}  —  "
                f"now click the second point")
        else:
            # Second click: create interval
            x1, x2 = sorted([self.first_click_x, event.xdata])
            self.intervals.append([x1, x2])
            self.first_click_x = None
            self.click_mode = False
            self.btn_add_interval.configure(
                text="+ Add Interval  (click on plot)", style='TButton')
            self.deconv_params = None
            self.deconv_n_gauss = None
            self._clear_results_text()
            self._refresh_intervals_ui()
            self._redraw()
            self.set_status(f"Interval added: [{x1:.2f}, {x2:.2f}]")

    def _refresh_intervals_ui(self):
        for w in self.interval_frame.winfo_children():
            w.destroy()
        for i, (xl, xr) in enumerate(self.intervals):
            row = tk.Frame(self.interval_frame, bg=BG_SIDEBAR)
            row.pack(fill=tk.X, pady=1)
            tk.Label(row, text=f"  {i+1}:  [{xl:.1f}  –  {xr:.1f}]",
                     bg=BG_SIDEBAR, fg=FG_TEXT, font=('Segoe UI', 9),
                     anchor=tk.W).pack(side=tk.LEFT, fill=tk.X, expand=True)
            tk.Button(row, text="×", fg=ACCENT_RED, bg=BG_SIDEBAR,
                      font=('Segoe UI', 10, 'bold'), borderwidth=0,
                      activebackground=BG_SIDEBAR, activeforeground='#ff1744',
                      cursor='hand2',
                      command=lambda idx=i: self._remove_interval(idx)
                      ).pack(side=tk.RIGHT)

    def _remove_interval(self, idx):
        if 0 <= idx < len(self.intervals):
            self.intervals.pop(idx)
            self.deconv_params = None
            self.deconv_n_gauss = None
            self._clear_results_text()
            self._refresh_intervals_ui()
            self._redraw()

    def _clear_intervals(self):
        self.intervals.clear()
        self.deconv_params = None
        self.deconv_n_gauss = None
        self._clear_results_text()
        self._refresh_intervals_ui()
        self._redraw()
        self.set_status("All intervals cleared")

    # deconvolution

    def _run_deconvolution(self):
        if self.x is None:
            messagebox.showinfo("No Data", "Load a signal first.")
            return
        if not self.intervals:
            messagebox.showinfo(
                "No Intervals",
                "Add at least one interval by clicking on the plot.")
            return

        n_gauss = int(self.ngauss_var.get())
        self.btn_run.configure(state='disabled', text="⏳  Running…")
        self.set_status(
            "Running differential evolution… this may take a few seconds")

        # Freeze current state for thread safety
        x = self.x.copy()
        y = self.y.copy()
        intervals = [iv[:] for iv in self.intervals]

        thread = threading.Thread(
            target=self._deconv_worker,
            args=(x, y, intervals, n_gauss),
            daemon=True)
        thread.start()

    def _deconv_worker(self, x, y, intervals, n_gauss):
        try:
            centers = auto_detect_centers(x, y, intervals)
            params, fun = optimize_gaussians(
                x, y, centers, intervals, num_gaussians=n_gauss)
            self.after(0, self._deconv_done, params, n_gauss)
        except Exception as e:
            self.after(0, self._deconv_fail, str(e))

    def _deconv_done(self, params, n_gauss):
        self.deconv_params = params
        self.deconv_n_gauss = n_gauss
        # Store in all_results
        sig_name = self.pairs[self.current_idx][2]
        self.all_results[sig_name] = (
            params, n_gauss, self.x.copy(), self.y.copy(),
            self.y_original.copy())
        self.btn_run.configure(state='normal',
                               text="▶  Run Deconvolution")
        self._show_results()
        self._redraw()

    def _deconv_fail(self, error_msg):
        self.btn_run.configure(state='normal',
                               text="▶  Run Deconvolution")
        messagebox.showerror("Deconvolution Failed", error_msg)
        self.set_status(
            "Deconvolution failed — check intervals and try again")

    def _show_results(self):
        # fill results panel with peak info
        if self.deconv_params is None:
            return
        ng = self.deconv_n_gauss
        nc = len(self.intervals)
        block = 3 if ng == 1 else 5
        y_fit = multi_gaussian(self.x, self.deconv_params, nc, ng)
        sse = np.sum((self.y - y_fit) ** 2)

        lines = [f"SSE = {sse:.2e}", ""]
        for j in range(nc):
            bp = self.deconv_params[j * block:j * block + block]
            A, mu, w = peak_summary(bp)
            if ng == 1:
                lines.append(
                    f"Peak {j+1}:  A={A:.3f}  μ={mu:.3f}  σ={w:.3f}")
            else:
                a1, a2, mu_val, s1, s2 = bp
                lines.append(
                    f"Peak {j+1}:  μ={mu_val:.3f}  (total A={A:.3f})")
                lines.append(
                    f"   G1: A={a1:.3f}  σ={s1:.3f}")
                lines.append(
                    f"   G2: A={a2:.3f}  σ={s2:.3f}")

        self.results_text.configure(state='normal')
        self.results_text.delete('1.0', tk.END)
        self.results_text.insert('1.0', '\n'.join(lines))
        self.results_text.configure(state='disabled')

    def _clear_results_text(self):
        self.results_text.configure(state='normal')
        self.results_text.delete('1.0', tk.END)
        self.results_text.configure(state='disabled')

    # batch deconvolution

    def _run_batch_deconvolution(self):
        # apply same config to all signals
        if self.deconv_params is None:
            messagebox.showinfo(
                "No Configuration",
                "First deconvolve the current signal, then use this button\n"
                "to apply the same config to all remaining signals.")
            return
        if len(self.pairs) < 2:
            messagebox.showinfo("Single Signal", "Only one signal loaded.")
            return

        # Freeze current config
        pipeline = list(self.pipeline)
        intervals = [iv[:] for iv in self.intervals]
        n_gauss = self.deconv_n_gauss

        # Find remaining signal indices
        remaining = [i for i in range(len(self.pairs)) if i != self.current_idx]
        if not remaining:
            self.set_status("No remaining signals to process.")
            return

        self._batch_running = True
        self.btn_batch.configure(state='disabled', text="⏳  Processing…")
        self.btn_run.configure(state='disabled')

        thread = threading.Thread(
            target=self._batch_worker,
            args=(remaining, pipeline, intervals, n_gauss),
            daemon=True)
        thread.start()

    def _batch_worker(self, remaining, pipeline, intervals, n_gauss):
        # process each signal in turn
        for idx in remaining:
            x_raw, y_raw, name = self.pairs[idx]
            x = x_raw.copy()
            y = y_raw.copy()
            y_orig = y.copy()

            # Apply preprocessing pipeline
            for step_name, step_params in pipeline:
                if step_name == 'shift':
                    y = shift_minimum(y)
                elif step_name == 'baseline':
                    y = linear_baseline_correction(x, y)
                elif step_name == 'clip_neg':
                    y = remove_negatives(y)
                elif step_name == 'sg':
                    y = filter_noise_sg(x, y, step_params['wl'],
                                        step_params['po'])
                elif step_name == 'ft':
                    y = filter_noise_ft(x, y, step_params['nf'])
                if np.any(~np.isfinite(y)):
                    y = np.nan_to_num(y, nan=0.0, posinf=0.0, neginf=0.0)

            # Show preprocessing result on main thread
            self.after(0, self._batch_show_preproc, idx, x, y, y_orig, name)
            import time
            time.sleep(1.5)  # Pause so user can see the preprocessing

            # Deconvolve
            centers = auto_detect_centers(x, y, intervals)
            try:
                params, fun = optimize_gaussians(
                    x, y, centers, intervals, num_gaussians=n_gauss)
            except Exception as e:
                self.after(0, self.set_status,
                           f"Failed on {name}: {e}")
                continue

            # Store result
            self.all_results[name] = (params, n_gauss, x, y, y_orig)

            # Show deconvolution result on main thread
            self.after(0, self._batch_show_deconv,
                       idx, x, y, y_orig, name, params, n_gauss, intervals)
            time.sleep(2.0)  # Pause so user can see the fit

        self.after(0, self._batch_done)

    def _batch_show_preproc(self, idx, x, y, y_orig, name):
        # show preprocessing for a batch signal
        self.ax.clear()
        self._style_axes()
        self.ax.plot(x, y_orig, color='#888888', lw=1.0,
                     label='Original', alpha=0.35, zorder=2)
        self.ax.plot(x, y, color=ACCENT_BLUE, lw=1.4,
                     label='Preprocessed', alpha=0.85, zorder=3)
        self.ax.set_title(f'Preprocessing — {name}', fontsize=11, pad=10)
        self.ax.legend(loc='upper right', fontsize=8, framealpha=0.7,
                       facecolor=BG_SIDEBAR, edgecolor=BORDER,
                       labelcolor=FG_TEXT)
        self.fig.tight_layout()
        self.canvas_mpl.draw()
        self.signal_list.selection_clear(0, tk.END)
        self.signal_list.selection_set(idx)
        self.set_status(f"Preprocessing: {name}")

    def _batch_show_deconv(self, idx, x, y, y_orig, name,
                           params, n_gauss, intervals):
        # show deconv result for a batch signal
        nc = len(intervals)
        block = 3 if n_gauss == 1 else 5
        self.ax.clear()
        self._style_axes()

        self.ax.plot(x, y_orig, color='#888888', lw=1.0,
                     label='Original', alpha=0.25, zorder=2)
        self.ax.plot(x, y, color=ACCENT_BLUE, lw=1.4,
                     label='Preprocessed', alpha=0.85, zorder=3)

        y_fit = multi_gaussian(x, params, nc, n_gauss)
        sse = np.sum((y - y_fit) ** 2)

        self.ax.plot(x, y_fit, color='#FF7043', lw=2,
                    ls='--', label='Total fit', zorder=5)

        for i in range(nc):
            ci = i % len(GAUSS_COLORS)
            if n_gauss == 1:
                a, mu, sig = params[3*i:3*i+3]
                yg = gaussian(x, a, mu, sig)
                self.ax.plot(x, yg, color=GAUSS_COLORS[ci], lw=1.2,
                            ls=':', label=f'G{i+1}', zorder=4)
                self.ax.fill_between(x, 0, yg,
                                     color=GAUSS_COLORS[ci], alpha=0.07)
            else:
                a1, a2, mu, s1, s2 = params[5*i:5*i+5]
                y1 = gaussian(x, a1, mu, s1)
                y2 = gaussian(x, a2, mu, s2)
                self.ax.plot(x, y1, color=GAUSS_COLORS[ci], lw=1.2,
                            ls=':', label=f'P{i+1}-G1', zorder=4)
                self.ax.plot(x, y2, color=GAUSS_COLORS[ci], lw=1.2,
                            ls='-.', label=f'P{i+1}-G2', zorder=4)
                self.ax.fill_between(x, 0, y1 + y2,
                                     color=GAUSS_COLORS[ci], alpha=0.05)

        for xl, xr in intervals:
            self.ax.axvspan(xl, xr, alpha=0.10, color='#FFD54F', zorder=1)

        self.ax.set_title(
            f'{name}  |  SSE = {sse:.2e}',
            fontsize=11, pad=10)
        self.ax.legend(loc='upper right', fontsize=8, framealpha=0.7,
                       facecolor=BG_SIDEBAR, edgecolor=BORDER,
                       labelcolor=FG_TEXT)
        self.fig.tight_layout()
        self.canvas_mpl.draw()
        self.set_status(
            f"{name}  —  SSE = {sse:.2e}")

        # Update results text with this signal's parameters
        lines = [f"— {name} —", f"SSE = {sse:.2e}", ""]
        for j in range(nc):
            bp = params[j * block:j * block + block]
            A, mu_val, w = peak_summary(bp)
            if n_gauss == 1:
                lines.append(
                    f"Peak {j+1}:  A={A:.3f}  μ={mu_val:.3f}  σ={w:.3f}")
            else:
                a1, a2, mu_v, s1, s2 = bp
                lines.append(
                    f"Peak {j+1}:  μ={mu_v:.3f}  (total A={A:.3f})")
                lines.append(
                    f"   G1: A={a1:.3f}  σ={s1:.3f}")
                lines.append(
                    f"   G2: A={a2:.3f}  σ={s2:.3f}")
        self.results_text.configure(state='normal')
        self.results_text.delete('1.0', tk.END)
        self.results_text.insert('1.0', '\n'.join(lines))
        self.results_text.configure(state='disabled')

    def _batch_done(self):
        self._batch_running = False
        self.btn_batch.configure(state='normal',
                                  text="▶  Deconvolve All (same config)")
        self.btn_run.configure(state='normal')
        n = len(self.all_results)
        # restore current signal so _redraw shows the right plot
        if self.current_idx is not None and self.current_idx < len(self.pairs):
            x_raw, y_raw, name = self.pairs[self.current_idx]
            self.x = x_raw.copy()
            self.y_original = y_raw.copy()
            # if this signal was batch-processed, use the preprocessed y
            if name in self.all_results:
                _, _, _, y_proc, _ = self.all_results[name]
                self.y = y_proc.copy()
                res = self.all_results[name]
                self.deconv_params = res[0]
                self.deconv_n_gauss = res[1]
            else:
                self.y = y_raw.copy()
            self._redraw()
        self.set_status(
            f"Batch complete — {n} signal(s) deconvolved. "
            f"Click 'Save Results' to export all.")

    # save results

    def _save_results(self):
        has_preproc = len(self.pipeline) > 0
        has_deconv = (self.deconv_params is not None
                      or len(self.all_results) > 0)

        if not has_preproc and not has_deconv:
            messagebox.showinfo(
                "Nothing to Save",
                "Apply preprocessing or run deconvolution first.")
            return

        initial_dir = (os.path.dirname(self.filepath)
                       if self.filepath else None)
        base_name = (os.path.splitext(os.path.basename(self.filepath))[0]
                     if self.filepath else "output")
        saved_files = []

        # Save preprocessed CSV (only if preprocessing was applied)
        if has_preproc:
            preproc_path = os.path.join(
                initial_dir or '.', f"{base_name}_preprocessed.csv")
            fp_pre = filedialog.asksaveasfilename(
                title="Save Preprocessed Signal",
                defaultextension=".csv",
                filetypes=[("CSV files", "*.csv")],
                initialfile=os.path.basename(preproc_path),
                initialdir=initial_dir)
            if fp_pre:
                try:
                    df_pre = pd.DataFrame({
                        'Position': self.x,
                        f"{self.pairs[self.current_idx][2]}_preprocessed": self.y
                    })
                    # If batch was run, add all preprocessed signals
                    for name, (_, _, x, y, _) in self.all_results.items():
                        if name != self.pairs[self.current_idx][2]:
                            df_pre[f"{name}_preprocessed"] = y
                    df_pre.to_csv(fp_pre, index=False)
                    saved_files.append(os.path.basename(fp_pre))
                except Exception as e:
                    messagebox.showerror("Save Error", str(e))

        # Save deconvolution CSV (only if deconv was done)
        if has_deconv:
            deconv_path = os.path.join(
                initial_dir or '.', f"{base_name}_deconv.csv")
            fp_dec = filedialog.asksaveasfilename(
                title="Save Deconvolution Results",
                defaultextension=".csv",
                filetypes=[("CSV files", "*.csv")],
                initialfile=os.path.basename(deconv_path),
                initialdir=initial_dir)
            if fp_dec:
                try:
                    rows = []
                    # Collect results – current signal + any batch results
                    results_to_save = dict(self.all_results)
                    sig_name = self.pairs[self.current_idx][2]
                    if sig_name not in results_to_save and self.deconv_params is not None:
                        results_to_save[sig_name] = (
                            self.deconv_params, self.deconv_n_gauss,
                            self.x.copy(), self.y.copy(),
                            self.y_original.copy())

                    for name, (params, ng, x, y, _) in results_to_save.items():
                        block = 3 if ng == 1 else 5
                        nc = len(self.intervals)
                        y_fit = multi_gaussian(x, params, nc, ng)
                        sse = float(np.sum((y - y_fit) ** 2))
                        for j in range(nc):
                            bp = params[j * block:j * block + block]
                            A, mu, w = peak_summary(bp)
                            row = {
                                'Signal': name,
                                'Peak': j + 1,
                                'Amplitude (A)': round(float(A), 6),
                                'Center (µ)': round(float(mu), 6),
                                'Width (σ)': round(float(w), 6),
                                'SSE': sse,
                            }
                            if ng == 2:
                                row.update({
                                    'a1': float(bp[0]),
                                    'a2': float(bp[1]),
                                    'sigma1': float(bp[3]),
                                    'sigma2': float(bp[4]),
                                })
                            rows.append(row)

                    pd.DataFrame(rows).to_csv(fp_dec, index=False, encoding='utf-8-sig')
                    saved_files.append(os.path.basename(fp_dec))
                except Exception as e:
                    messagebox.showerror("Save Error", str(e))

        if saved_files:
            self.set_status(f"Saved: {', '.join(saved_files)}")
        else:
            self.set_status("Save cancelled")


# Entry point

if __name__ == "__main__":
    try:
        app = MultiGaussFitApp()
        app.mainloop()
    except Exception:
        traceback.print_exc()
        input("Press Enter to close...")
