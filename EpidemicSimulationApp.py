import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import os
from dataclasses import dataclass
from typing import List, Tuple
from PIL import Image


@dataclass
class SimulationParameters:
    """Класс для хранения параметров симуляции эпидемии"""
    grid_size: int = 100
    initial_infected: int = 5
    infection_probability: float = 0.5
    recovery_time: int = 14
    mortality_rate: float = 0.02
    initial_vaccinated_count: int = 0
    social_distancing: float = 0.0
    vaccinated_infection_probability: float = 0.0


class CellularAutomaton:
    """Класс клеточного автомата для моделирования эпидемии"""

    SUSCEPTIBLE = 0
    INFECTED = 1
    RECOVERED = 2
    DEAD = 3
    VACCINATED = 4
    NEUTRAL = 5

    def __init__(self, params: SimulationParameters):
        self.params = params
        self.grid_size = params.grid_size
        self.grid = np.zeros((self.grid_size, self.grid_size), dtype=int)
        self.infection_time = np.zeros((self.grid_size, self.grid_size), dtype=int)
        self.statistics = {
            'susceptible': [], 'infected': [], 'recovered': [],
            'dead': [], 'vaccinated': [], 'neutral': []
        }
        self.iteration = 0
        self.grid_history: List[np.ndarray] = []
        self._initialize_grid()

    def _initialize_grid(self):
        self.grid.fill(self.SUSCEPTIBLE)
        self.infection_time.fill(0)
        num_total_cells = self.grid_size * self.grid_size

        if self.grid_size <= 0:
            self._update_statistics()
            if not self.grid_history:
                self.grid_history.append(self.grid.copy())
            else:
                self.grid_history[0] = self.grid.copy()
            return

        available_for_vaccination_indices = np.where(self.grid.flatten() == self.SUSCEPTIBLE)[0]
        current_initial_vaccinated = min(self.params.initial_vaccinated_count, len(available_for_vaccination_indices))
        if current_initial_vaccinated < 0: current_initial_vaccinated = 0

        if current_initial_vaccinated > 0:
            chosen_vacc_indices = np.random.choice(
                available_for_vaccination_indices,
                current_initial_vaccinated,
                replace=False
            )
            for flat_idx in chosen_vacc_indices:
                i, j = np.unravel_index(flat_idx, (self.grid_size, self.grid_size))
                self.grid[i, j] = self.VACCINATED

        available_indices_for_infection = np.where(self.grid.flatten() == self.SUSCEPTIBLE)[0]
        current_initial_infected = min(self.params.initial_infected, len(available_indices_for_infection))
        if current_initial_infected < 0: current_initial_infected = 0

        if current_initial_infected > 0:
            infected_indices_chosen = np.random.choice(
                available_indices_for_infection, current_initial_infected, replace=False
            )
            for flat_idx in infected_indices_chosen:
                i, j = np.unravel_index(flat_idx, (self.grid_size, self.grid_size))
                self.grid[i, j] = self.INFECTED
                self.infection_time[i, j] = self.iteration

        self._update_statistics()
        if not self.grid_history:
            self.grid_history.append(self.grid.copy())
        else:
            self.grid_history[0] = self.grid.copy()
            self.grid_history = self.grid_history[:1]

    def step(self):
        self.iteration += 1
        new_grid = self.grid.copy()
        for i in range(self.grid_size):
            for j in range(self.grid_size):
                if self.grid[i, j] == self.NEUTRAL:
                    continue
                if self.grid[i, j] == self.INFECTED:
                    if self.iteration - self.infection_time[i, j] >= self.params.recovery_time:
                        if np.random.random() < self.params.mortality_rate:
                            new_grid[i, j] = self.DEAD
                        else:
                            new_grid[i, j] = self.RECOVERED
                elif self.grid[i, j] == self.SUSCEPTIBLE or self.grid[i, j] == self.VACCINATED:
                    infected_neighbors_count = self._count_infected_neighbors(i, j)
                    if infected_neighbors_count > 0:

                        if self.grid[i, j] == self.VACCINATED:

                            base_prob_for_event = self.params.vaccinated_infection_probability
                            num_events_to_check = 1
                        else:
                            base_prob_for_event = self.params.infection_probability
                            num_events_to_check = infected_neighbors_count

                        # Calculate the effective probability for each event, considering social distancing.
                        effective_event_prob = base_prob_for_event * (1 - self.params.social_distancing)

                        # Perform infection attempts based on the number of events.
                        # For vaccinated, this loop runs once.
                        # For susceptible, this loop runs up to infected_neighbors_count times.
                        for _ in range(num_events_to_check):
                            if np.random.random() < effective_event_prob:
                                new_grid[i, j] = self.INFECTED
                                self.infection_time[i, j] = self.iteration
                                break
        self.grid = new_grid
        self._update_statistics()
        self.grid_history.append(self.grid.copy())
        return self.grid

    def _count_infected_neighbors(self, i, j):
        count = 0
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                if di == 0 and dj == 0: continue
                ni, nj = (i + di) % self.grid_size, (j + dj) % self.grid_size
                if self.grid[ni, nj] == self.INFECTED:
                    count += 1
        return count

    def _update_statistics(self):
        counts = {state: 0 for state in range(6)}
        if self.grid_size > 0:
            for row in self.grid:
                for cell_state in row:
                    counts[cell_state] += 1
        if not self.statistics['susceptible']:
            for key in self.statistics: self.statistics[key] = []

        self.statistics['susceptible'].append(counts[self.SUSCEPTIBLE])
        self.statistics['infected'].append(counts[self.INFECTED])
        self.statistics['recovered'].append(counts[self.RECOVERED])
        self.statistics['dead'].append(counts[self.DEAD])
        self.statistics['vaccinated'].append(counts[self.VACCINATED])
        self.statistics['neutral'].append(counts[self.NEUTRAL])

    def get_statistics(self):
        return self.statistics

    def get_grid(self):
        return self.grid

    def _grid_to_colored_grid(self, grid_state: np.ndarray) -> np.ndarray:
        if self.grid_size <= 0:
            return np.zeros((0, 0, 4), dtype=np.float32)

        colored_grid = np.zeros((self.grid_size, self.grid_size, 4), dtype=np.float32)
        colors = {
            self.SUSCEPTIBLE: [0.0, 0.0, 0.0, 0.0],
            self.INFECTED: [0.9, 0.1, 0.1, 0.75],
            self.RECOVERED: [0.1, 0.7, 0.1, 0.75],
            self.DEAD: [0.1, 0.1, 0.1, 0.85],
            self.VACCINATED: [0.1, 0.1, 0.9, 0.75],
            self.NEUTRAL: [0.5, 0.5, 0.5, 0.75]
        }
        for state, color in colors.items():
            mask = grid_state == state
            colored_grid[mask] = color
        return colored_grid

    def get_colored_grid(self):
        return self._grid_to_colored_grid(self.grid)

    def update_current_snapshot_stats_and_history(self):
        current_counts = {s: 0 for s in range(6)}
        if self.grid_size > 0:
            for r_idx in range(self.grid_size):
                for c_idx in range(self.grid_size):
                    current_counts[self.grid[r_idx, c_idx]] += 1

        if not self.grid_history or not self.statistics['susceptible']:
            self.grid_history = [self.grid.copy()]
            for key in self.statistics: self.statistics[key] = []
            self.statistics['susceptible'].append(current_counts[self.SUSCEPTIBLE])
            self.statistics['infected'].append(current_counts[self.INFECTED])
            self.statistics['recovered'].append(current_counts[self.RECOVERED])
            self.statistics['dead'].append(current_counts[self.DEAD])
            self.statistics['vaccinated'].append(current_counts[self.VACCINATED])
            self.statistics['neutral'].append(current_counts[self.NEUTRAL])
        else:
            self.grid_history[-1] = self.grid.copy()
            self.statistics['susceptible'][-1] = current_counts[self.SUSCEPTIBLE]
            self.statistics['infected'][-1] = current_counts[self.INFECTED]
            self.statistics['recovered'][-1] = current_counts[self.RECOVERED]
            self.statistics['dead'][-1] = current_counts[self.DEAD]
            self.statistics['vaccinated'][-1] = current_counts[self.VACCINATED]
            self.statistics['neutral'][-1] = current_counts[self.NEUTRAL]


class EpidemicSimulationApp:
    def __init__(self, root_window):
        self.root = root_window
        self.root.title("Моделирование эпидемии - Клеточные автоматы")
        self.root.geometry("1350x900")

        self.params = SimulationParameters()
        self.automaton = None
        self.animation_obj = None
        self.is_playing = False
        self.speed = 200
        self.stats_annot = None
        self.last_annot_info = None
        self.hover_connection_id = None
        self.stats_ax_background = None  # For blitting optimization

        self.background_image_pil = None
        self.background_image_path = None
        self.processed_background_image_mpl = None
        self.background_image_alpha_var = tk.DoubleVar(value=0.5)
        self.background_image_visible = tk.BooleanVar(value=True)
        self.background_image_filename_var = tk.StringVar(value="Файл не выбран")
        self._background_changed_flag = False

        self.grid_lines_visible = tk.BooleanVar(value=True)

        self.drawing_mode_var = tk.StringVar(value="none")
        self.brush_radius_var = tk.IntVar(value=0)
        self.is_drawing_on_grid = False
        self.grid_press_cid = None
        self.grid_motion_cid = None
        self.grid_release_cid = None
        self.drawing_tools_frame = None
        self.grid_image_artist = None
        self.last_drawn_cells = set()
        self.vaccinated_infection_prob_var = tk.DoubleVar(value=self.params.vaccinated_infection_probability)

        self._create_ui()
        self._init_simulation()
        self.root.protocol("WM_DELETE_WINDOW", self._on_closing)

    def _on_closing(self):
        if self.is_playing and self.animation_obj and hasattr(self.animation_obj,
                                                              'event_source') and self.animation_obj.event_source:
            self.animation_obj.event_source.stop()
        self.is_playing = False
        if self.hover_connection_id and self.stats_canvas:
            try:
                self.stats_canvas.mpl_disconnect(self.hover_connection_id)
            except:
                pass
        if self.grid_press_cid and self.grid_canvas:
            try:
                self.grid_canvas.mpl_disconnect(self.grid_press_cid)
            except:
                pass
        if self.grid_motion_cid and self.grid_canvas:
            try:
                self.grid_canvas.mpl_disconnect(self.grid_motion_cid)
            except:
                pass
        if self.grid_release_cid and self.grid_canvas:
            try:
                self.grid_canvas.mpl_disconnect(self.grid_release_cid)
            except:
                pass
        plt.close('all')
        self.root.quit()
        self.root.destroy()

    def _create_scale_with_label(self, parent, text, row, variable, from_, to_, format_str="{:.2f}", column_offset=0):
        ttk.Label(parent, text=text).grid(row=row, column=0 + column_offset, padx=5, pady=5, sticky="w")
        scale = ttk.Scale(parent, from_=from_, to=to_, variable=variable, orient="horizontal")
        scale.grid(row=row, column=1 + column_offset, padx=5, pady=5, sticky="ew")
        value_label = ttk.Label(parent, text=format_str.format(variable.get()))
        value_label.grid(row=row, column=2 + column_offset, padx=5, pady=5, sticky="w")
        scale.config(
            command=lambda v, lbl=value_label, var=variable, fmt=format_str: lbl.config(text=fmt.format(var.get())))
        return scale, value_label

    def _create_ui(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Файл", menu=file_menu)

        save_menu = tk.Menu(file_menu, tearoff=0)
        file_menu.add_cascade(label="Сохранить...", menu=save_menu)
        save_menu.add_command(label="Общая статистика (PNG)", command=self._save_overall_stats)
        save_menu.add_command(label="Текущие статистика и итерация (PNG)", command=self._save_individual_plots)
        save_menu.add_command(label="Анимация итераций (GIF)", command=self._save_gif_animation)
        save_menu.add_command(label="Изображения всех итераций (PNG)", command=self._save_all_iterations_images)
        file_menu.add_separator()
        file_menu.add_command(label="Выход", command=self._on_closing)

        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Справка", menu=help_menu)
        help_menu.add_command(label="О программе", command=self._show_about_info)

        left_panel = ttk.Frame(self.root, padding="10")
        left_panel.grid(row=0, column=0, sticky="nswe")
        self.root.grid_columnconfigure(0, weight=0)
        self.root.grid_columnconfigure(1, weight=1)
        self.root.grid_rowconfigure(0, weight=1)

        params_frame = ttk.LabelFrame(left_panel, text="Параметры симуляции")
        params_frame.pack(fill=tk.X, pady=5)

        # Row 0: Grid Size
        ttk.Label(params_frame, text="Размер сетки:").grid(row=0, column=0, padx=5, pady=2, sticky="w")
        self.grid_size_var = tk.IntVar(value=self.params.grid_size)
        ttk.Spinbox(params_frame, from_=10, to=200, textvariable=self.grid_size_var, width=7).grid(row=0, column=1,
                                                                                                   columnspan=3, padx=5,
                                                                                                   pady=2,
                                                                                                   sticky="ew")

        # Row 1: Initial Infected and Initial Vaccinated
        ttk.Label(params_frame, text="Начальное значение Инфицированных:").grid(row=1, column=0, padx=5, pady=2, sticky="w")
        self.initial_infected_var = tk.IntVar(value=self.params.initial_infected)
        ttk.Spinbox(params_frame, from_=0, to=40000, textvariable=self.initial_infected_var, width=5).grid(row=1,
                                                                                                           # width changed
                                                                                                           column=1,
                                                                                                           padx=5,
                                                                                                           pady=2,
                                                                                                           sticky="ew")

        ttk.Label(params_frame, text="Вакцинированных:").grid(row=1, column=2, padx=5, pady=2, sticky="w")
        self.initial_vaccinated_var = tk.IntVar(value=self.params.initial_vaccinated_count)
        self.initial_vaccinated_spinbox = ttk.Spinbox(params_frame, from_=0, to=40000,
                                                      textvariable=self.initial_vaccinated_var,
                                                      width=5)
        self.initial_vaccinated_spinbox.grid(row=1, column=3,
                                             padx=5, pady=2, sticky="ew")

        # Row 2: Infection Probability
        self.infection_prob_var = tk.DoubleVar(value=self.params.infection_probability)
        self._create_scale_with_label(params_frame, "Вероятность инфицирования:", 2, self.infection_prob_var, 0.0, 1.0,
                                      column_offset=0)

        # Row 3: Vaccinated Infection Probability
        self.vaccinated_infection_prob_var = tk.DoubleVar(value=self.params.vaccinated_infection_probability)
        self._create_scale_with_label(params_frame, "Вероятность заражения вакцинированных:", 3, self.vaccinated_infection_prob_var,
                                      0.0, 1.0, column_offset=0)

        # Row 4: Recovery Time
        ttk.Label(params_frame, text="Время выздоровления (дни):").grid(row=4, column=0, padx=5, pady=2, sticky="w")
        self.recovery_time_var = tk.IntVar(value=self.params.recovery_time)
        ttk.Spinbox(params_frame, from_=1, to=100, textvariable=self.recovery_time_var, width=7).grid(row=4, column=1,
                                                                                                      columnspan=3,
                                                                                                      padx=5,
                                                                                                      pady=2,
                                                                                                      sticky="ew")

        # Row 5: Mortality Rate
        self.mortality_rate_var = tk.DoubleVar(value=self.params.mortality_rate)
        self._create_scale_with_label(params_frame, "Уровень смертности:", 5, self.mortality_rate_var, 0.0, 1.0,
                                      column_offset=0)

        # Initial Vaccinated related widgets were moved to row 1

        # Row 6: Social Distancing (previously row 7)
        self.social_distancing_var = tk.DoubleVar(value=self.params.social_distancing)
        self._create_scale_with_label(params_frame, "Социальное дистанцирование:", 6, self.social_distancing_var, 0.0,
                                      1.0, column_offset=0)

        # Row 7: Grid Lines Checkbox (previously row 8)
        self.grid_lines_checkbox = ttk.Checkbutton(params_frame, text="Показывать линии сетки",
                                                   variable=self.grid_lines_visible,
                                                   command=self._toggle_grid_lines_visibility)
        self.grid_lines_checkbox.grid(row=7, column=0, columnspan=4, padx=5, pady=2,
                                      sticky="w")

        # Row 8: Initialize/Update Button (previously row 9)
        ttk.Button(params_frame, text="Инициализировать/Обновить", command=self._update_params).grid(row=8, column=0,
                                                                                                     columnspan=4,
                                                                                                     padx=5, pady=10)
        params_frame.columnconfigure(1, weight=1)
        params_frame.columnconfigure(3, weight=1)

        self.drawing_tools_frame = ttk.LabelFrame(left_panel, text="Инструменты рисования")
        self.drawing_tools_frame.pack(fill=tk.X, pady=5)
        drawing_modes = [("Не рисовать", "none"), ("Инфицированные (ЛКМ)", "infected"),
                         ("Вакцинированные (ЛКМ)", "vaccinated"), ("Нейтральные (ЛКМ)", "neutral"),
                         ("Стереть (ЛКМ)", "erase")]
        for text, mode in drawing_modes:
            ttk.Radiobutton(self.drawing_tools_frame, text=text, variable=self.drawing_mode_var, value=mode).pack(
                anchor=tk.W, padx=5)
        brush_radius_frame = ttk.Frame(self.drawing_tools_frame)
        brush_radius_frame.pack(anchor=tk.W, padx=5, pady=2)
        ttk.Label(brush_radius_frame, text="Радиус кисти (0-одна клетка):").pack(side=tk.LEFT)
        ttk.Spinbox(brush_radius_frame, from_=0, to=20, textvariable=self.brush_radius_var, width=5).pack(side=tk.LEFT,
                                                                                                          padx=5)

        bg_map_frame = ttk.LabelFrame(left_panel, text="Карта фона")
        bg_map_frame.pack(fill=tk.X, pady=5)
        ttk.Button(bg_map_frame, text="Загрузить карту", command=self._load_background_image).grid(row=0, column=0,
                                                                                                   columnspan=3, padx=5,
                                                                                                   pady=5, sticky="ew")
        ttk.Label(bg_map_frame, text="Файл:").grid(row=1, column=0, padx=5, pady=2, sticky="w")
        self.bg_filename_label = ttk.Label(bg_map_frame, textvariable=self.background_image_filename_var,
                                           wraplength=200)
        self.bg_filename_label.grid(row=1, column=1, columnspan=2, padx=5, pady=2, sticky="ew")
        ttk.Label(bg_map_frame, text="Непрозрачность:").grid(row=2, column=0, padx=5, pady=2, sticky="w")
        self.opacity_scale = ttk.Scale(bg_map_frame, from_=0.0, to=1.0, variable=self.background_image_alpha_var,
                                       orient="horizontal", command=self._on_opacity_change)
        self.opacity_scale.grid(row=2, column=1, padx=5, pady=2, sticky="ew")
        self.opacity_value_label = ttk.Label(bg_map_frame, text=f"{self.background_image_alpha_var.get():.2f}")
        self.opacity_value_label.grid(row=2, column=2, padx=5, pady=2, sticky="w")
        self.opacity_scale.config(state=tk.DISABLED)
        self.toggle_bg_visibility_button = ttk.Button(bg_map_frame, text="Скрыть карту",
                                                      command=self._toggle_background_visibility)
        self.toggle_bg_visibility_button.grid(row=3, column=0, columnspan=3, padx=5, pady=5, sticky="ew")
        self.toggle_bg_visibility_button.config(state=tk.DISABLED)
        self.remove_bg_button = ttk.Button(bg_map_frame, text="Удалить карту", command=self._remove_background_image)
        self.remove_bg_button.grid(row=4, column=0, columnspan=3, padx=5, pady=5, sticky="ew")
        self.remove_bg_button.config(state=tk.DISABLED)
        bg_map_frame.columnconfigure(1, weight=1)

        stats_frame = ttk.LabelFrame(left_panel, text="Статистика")
        stats_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        self.stats_fig, self.stats_ax = plt.subplots()
        self.stats_canvas = FigureCanvasTkAgg(self.stats_fig, master=stats_frame)
        self.stats_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        right_panel = ttk.Frame(self.root)
        right_panel.grid(row=0, column=1, sticky="nswe", padx=10, pady=10)
        right_panel.grid_rowconfigure(0, weight=1)
        right_panel.grid_columnconfigure(0, weight=1)

        grid_frame_outer = ttk.LabelFrame(right_panel, text="Сетка симуляции")
        grid_frame_outer.grid(row=0, column=0, sticky="nswe")
        grid_frame_outer.grid_rowconfigure(0, weight=1)
        grid_frame_outer.grid_columnconfigure(0, weight=1)

        self.grid_fig, self.grid_ax = plt.subplots()
        self.grid_canvas = FigureCanvasTkAgg(self.grid_fig, master=grid_frame_outer)
        self.grid_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        playback_control_frame = ttk.Frame(right_panel)
        playback_control_frame.grid(row=1, column=0, sticky="ew", pady=10)
        playback_control_frame.columnconfigure(0, weight=1)
        playback_control_frame.columnconfigure(1, weight=1)
        playback_control_frame.columnconfigure(2, weight=1)

        buttons_subframe_main = ttk.Frame(playback_control_frame)
        buttons_subframe_main.pack(pady=5)
        ttk.Button(buttons_subframe_main, text="Шаг", command=self._step_forward).pack(side=tk.LEFT, padx=5)
        self.play_button = ttk.Button(buttons_subframe_main, text="Воспроизвести", command=self._toggle_play)
        self.play_button.pack(side=tk.LEFT, padx=5)
        ttk.Button(buttons_subframe_main, text="Сброс", command=self._reset_simulation).pack(side=tk.LEFT, padx=5)

        speed_control_frame_main = ttk.Frame(playback_control_frame)
        speed_control_frame_main.pack(pady=5)
        ttk.Label(speed_control_frame_main, text="Скорость:").pack(side=tk.LEFT, padx=5)
        self.speed_var = tk.IntVar(value=self.speed)
        speeds = [("Медленно", 400), ("Средне", 200), ("Быстро", 50)]
        for text, val in speeds:
            ttk.Radiobutton(speed_control_frame_main, text=text, value=val, variable=self.speed_var,
                            command=lambda v=val: self._set_speed(v)).pack(side=tk.LEFT, padx=2)

    def _load_background_image(self):
        filepath = filedialog.askopenfilename(
            title="Выберите изображение для фона",
            filetypes=[("Image files", "*.png *.jpg *.jpeg *.bmp *.gif *.tiff"), ("All files", "*.*")]
        )
        if not filepath: return
        try:
            self.background_image_pil = Image.open(filepath)
            self.background_image_path = filepath
            self.background_image_filename_var.set(os.path.basename(filepath))
            self.processed_background_image_mpl = None  # Force re-processing/re-rendering
            self.opacity_scale.config(state=tk.NORMAL)
            self.toggle_bg_visibility_button.config(state=tk.NORMAL)
            self.remove_bg_button.config(state=tk.NORMAL)

            # Ensure image is visible on new load and button text is correct
            if not self.background_image_visible.get():
                self.background_image_visible.set(True)
            self.toggle_bg_visibility_button.config(text="Скрыть карту")

            self._background_changed_flag = True
            self._update_grid_visualization()

        except Exception as e:
            messagebox.showerror("Ошибка загрузки изображения", f"Не удалось загрузить файл: {e}")
            self._reset_background_image_state()
            self._background_changed_flag = True
            self._update_grid_visualization()

    def _show_about_info(self):
        about_window = tk.Toplevel(self.root)
        about_window.title("О программе")
        # Установка предпочтительного размера, можно настроить
        about_window.geometry("580x280")
        about_window.resizable(False, False)  # Запрет изменения размера окна
        # Связать с главным окном (поверх него, минимизируется/восстанавливается с ним)
        about_window.transient(self.root)
        # Сделать окно модальным - блокирует взаимодействие с другими окнами приложения, пока это открыто
        about_window.grab_set()


        # Даем Tkinter время обработать очередь событий и обновить геометрию окон
        self.root.update_idletasks()
        about_window.update_idletasks()

        parent_x = self.root.winfo_x()
        parent_y = self.root.winfo_y()
        parent_width = self.root.winfo_width()
        parent_height = self.root.winfo_height()

        dialog_width = about_window.winfo_width()
        dialog_height = about_window.winfo_height()

        # Если winfo_width/height вернули 1 (окно еще не отображено полностью),
        # используем reqwidth/reqheight как запасной вариант для расчета размеров.
        if dialog_width <= 1: dialog_width = about_window.winfo_reqwidth()
        if dialog_height <= 1: dialog_height = about_window.winfo_reqheight()

        # Расчет координат X и Y для центрирования
        x_pos = parent_x + (parent_width // 2) - (dialog_width // 2)
        y_pos = parent_y + (parent_height // 2) - (dialog_height // 2)

        # Применение рассчитанных координат для позиционирования окна
        about_window.geometry(f"+{x_pos}+{y_pos}")

        # Основной фрейм для содержимого с внутренними отступами
        main_content_frame = ttk.Frame(about_window, padding="20")
        main_content_frame.pack(expand=True, fill=tk.BOTH)  # Растягивается по всему доступному пространству

        # Заголовок программы (берется из заголовка главного окна для консистентности)
        program_name_label = ttk.Label(main_content_frame, text=self.root.title(), font=("Arial", 16, "bold"))
        program_name_label.pack(pady=(0, 15))  # Отступ снизу для отделения от следующего текста

        # Текст с информацией о программе, версии и разработчике
        about_text_content = (
            "Программа для моделирования распространения эпидемии\n"
            "с использованием концепции клеточных автоматов.\n\n"
            "Версия: 1.0\n\n"
            "Разработчик: skail\n"
            "Email: iliushechkind@yandex.ru"
        )
        info_label = ttk.Label(main_content_frame, text=about_text_content, justify=tk.CENTER)
        info_label.pack(pady=10)  # Отступы сверху и снизу

        # Горизонтальный разделитель
        separator = ttk.Separator(main_content_frame, orient=tk.HORIZONTAL)
        separator.pack(fill=tk.X, pady=10)  # Растягивается по ширине, с отступами

        # Кнопка "OK" для закрытия диалогового окна
        ok_button = ttk.Button(main_content_frame, text="OK", command=about_window.destroy, width=12)
        ok_button.pack(pady=(5, 0))  # Небольшой отступ сверху перед кнопкой

        # Установить фокус на кнопку "OK", чтобы нажатие Enter ее активировало
        ok_button.focus_set()

        # Обработка закрытия окна по стандартному "крестику" в заголовке окна
        about_window.protocol("WM_DELETE_WINDOW", about_window.destroy)

        # Ожидание закрытия этого модального окна перед тем, как вернуть управление главному окну
        about_window.wait_window()

    def _remove_background_image(self):
        if not self.background_image_pil: return
        self._reset_background_image_state()
        self._background_changed_flag = True
        self._update_grid_visualization()

    def _reset_background_image_state(self):
        self.background_image_pil = None
        self.background_image_path = None
        self.background_image_filename_var.set("Файл не выбран")
        self.processed_background_image_mpl = None
        self.opacity_scale.config(state=tk.DISABLED)
        self.toggle_bg_visibility_button.config(state=tk.DISABLED, text="Скрыть карту")
        self.remove_bg_button.config(state=tk.DISABLED)
        self.background_image_visible.set(True)

    def _on_opacity_change(self, _=None):
        self.opacity_value_label.config(text=f"{self.background_image_alpha_var.get():.2f}")
        if self.background_image_pil:
            self._background_changed_flag = True
            self._update_grid_visualization()
            # Flag is reset inside _update_grid_visualization

    def _toggle_background_visibility(self):
        if not self.background_image_pil: return
        current_visibility = self.background_image_visible.get()
        self.background_image_visible.set(not current_visibility)
        self.toggle_bg_visibility_button.config(
            text="Скрыть карту" if self.background_image_visible.get() else "Показать карту")
        self._background_changed_flag = True
        self._update_grid_visualization()
        # Flag is reset inside _update_grid_visualization

    def _update_params(self):
        old_grid_size = self.params.grid_size
        self.params = SimulationParameters(
            grid_size=self.grid_size_var.get(),
            initial_infected=self.initial_infected_var.get(),
            infection_probability=self.infection_prob_var.get(),
            recovery_time=self.recovery_time_var.get(),
            mortality_rate=self.mortality_rate_var.get(),
            initial_vaccinated_count=self.initial_vaccinated_var.get(),
            social_distancing=self.social_distancing_var.get(),
            vaccinated_infection_probability=self.vaccinated_infection_prob_var.get()
        )
        if old_grid_size != self.params.grid_size:
            self.processed_background_image_mpl = None
            self.grid_image_artist = None
        self._init_simulation()

    def _init_simulation(self):
        if self.animation_obj is not None and hasattr(self.animation_obj,
                                                      'event_source') and self.animation_obj.event_source:
            self.animation_obj.event_source.stop()
        self.is_playing = False
        self.play_button.config(text="Воспроизвести")
        self.drawing_mode_var.set("none")
        self.grid_image_artist = None

        if self.drawing_tools_frame:
            for widget in self.drawing_tools_frame.winfo_children():
                if isinstance(widget, (ttk.Radiobutton, ttk.Spinbox, ttk.LabelFrame)):
                    try:
                        widget.config(state=tk.NORMAL)
                        if isinstance(widget, ttk.LabelFrame):
                            for sub_widget in widget.winfo_children():
                                if isinstance(sub_widget, ttk.Spinbox):
                                    sub_widget.config(state=tk.NORMAL)
                    except tk.TclError:
                        pass

        self.automaton = CellularAutomaton(self.params)
        for key in self.automaton.statistics: self.automaton.statistics[key] = []
        self.automaton.grid_history = []
        self.automaton._initialize_grid()

        self._update_grid_visualization()
        self._update_stats_visualization()

        if self.automaton and self.automaton.grid_size > 0:
            current_grid = self.automaton.get_grid()
            self.initial_infected_var.set(np.sum(current_grid == CellularAutomaton.INFECTED))
            self.initial_vaccinated_var.set(np.sum(current_grid == CellularAutomaton.VACCINATED))
        elif self.automaton and self.automaton.grid_size == 0:
            self.initial_infected_var.set(0)
            self.initial_vaccinated_var.set(0)

        self.animation_obj = animation.FuncAnimation(self.grid_fig, self._animate, interval=self.speed, blit=False,
                                                     save_count=50)
        if hasattr(self.animation_obj, 'event_source') and self.animation_obj.event_source:
            self.animation_obj.event_source.stop()

        if self.stats_canvas:
            if self.hover_connection_id: self.stats_canvas.mpl_disconnect(self.hover_connection_id)
            self.hover_connection_id = self.stats_canvas.mpl_connect('motion_notify_event', self._on_stats_hover)

        if self.grid_canvas:
            if self.grid_press_cid: self.grid_canvas.mpl_disconnect(self.grid_press_cid)
            if self.grid_motion_cid: self.grid_canvas.mpl_disconnect(self.grid_motion_cid)
            if self.grid_release_cid: self.grid_canvas.mpl_disconnect(self.grid_release_cid)
            self.grid_press_cid = self.grid_canvas.mpl_connect('button_press_event', self._on_grid_press)
            self.grid_motion_cid = self.grid_canvas.mpl_connect('motion_notify_event', self._on_grid_motion)
            self.grid_release_cid = self.grid_canvas.mpl_connect('button_release_event', self._on_grid_release)

        # Захватываем фон для статистики для оптимизации (blitting)
        if self.stats_canvas and self.stats_fig.canvas.supports_blit:
            self.stats_fig.canvas.draw()  # Начальная отрисовка
            self.stats_ax_background = self.stats_fig.canvas.copy_from_bbox(self.stats_ax.bbox)

    def _animate(self, i):
        if self.is_playing and self.automaton:
            if self.automaton.grid_size <= 0:
                self._toggle_play(force_pause=True)
                messagebox.showinfo("Симуляция остановлена", "Размер сетки 0, симуляция не может продолжаться.")
                return [self.grid_image_artist] if self.grid_image_artist else []  # Возвращаем список артистов

            self.automaton.step()

            # Обновляем заголовок сетки с новой итерацией
            if self.automaton:
                self.grid_ax.set_title(f"Итерация: {self.automaton.iteration}")

            if self.grid_image_artist:
                self.grid_image_artist.set_data(self.automaton.get_colored_grid())
                self.grid_canvas.draw_idle()  # Обновляем холст
            else:
                self._update_grid_visualization()  # Полная перерисовка, если артиста нет

            self._update_stats_visualization()
            stats = self.automaton.get_statistics()
            if self.automaton.iteration > 0 and stats['infected'] and stats['infected'][-1] == 0:
                was_infected_previously = any(count > 0 for count in stats['infected'][:-1]) if len(
                    stats['infected']) > 1 else False
                if len(stats['infected']) == 1 and self.params.initial_infected > 0: was_infected_previously = True
                if was_infected_previously:
                    self._toggle_play(force_pause=True)
                    messagebox.showinfo("Симуляция завершена", "Эпидемия закончилась (нет инфицированных).")

        return [self.grid_image_artist] if self.grid_image_artist else []

    def _paint_at_event_location(self, event, is_motion=False):
        if event.xdata is None or event.ydata is None or not self.automaton: return False
        if self.automaton.grid_size <= 0: return False

        center_col = int(round(event.xdata))
        center_row = int(round(event.ydata))

        current_brush_center_cell = (center_row, center_col)
        if is_motion and current_brush_center_cell == getattr(self, '_last_painted_center_cell', None):
            return False
        self._last_painted_center_cell = current_brush_center_cell

        if not (0 <= center_row < self.automaton.grid_size and 0 <= center_col < self.automaton.grid_size):
            return False

        mode = self.drawing_mode_var.get()
        target_state = -1
        if mode == "infected":
            target_state = CellularAutomaton.INFECTED
        elif mode == "vaccinated":
            target_state = CellularAutomaton.VACCINATED
        elif mode == "neutral":
            target_state = CellularAutomaton.NEUTRAL
        elif mode == "erase":
            target_state = CellularAutomaton.SUSCEPTIBLE
        else:
            return False

        brush_radius = self.brush_radius_var.get()
        changed_anything = False

        if brush_radius == 0:  # Single point drawing
            original_state_at_center = self.automaton.grid[center_row, center_col]

            can_paint = True
            if target_state == CellularAutomaton.VACCINATED:
                if original_state_at_center != CellularAutomaton.SUSCEPTIBLE:
                    can_paint = False
            elif original_state_at_center == CellularAutomaton.NEUTRAL and \
                    target_state != CellularAutomaton.SUSCEPTIBLE:
                can_paint = False # Do not paint over NEUTRAL cells unless erasing

            if can_paint and original_state_at_center != target_state:
                self.automaton.grid[center_row, center_col] = target_state
                changed_anything = True
                if target_state == CellularAutomaton.INFECTED:
                    self.automaton.infection_time[center_row, center_col] = self.automaton.iteration
                elif original_state_at_center == CellularAutomaton.INFECTED and \
                        target_state == CellularAutomaton.SUSCEPTIBLE:
                    self.automaton.infection_time[center_row, center_col] = 0
        else:  # Brush drawing (radius > 0)
            brush_radius_sq = brush_radius ** 2
            min_r = max(0, center_row - brush_radius)
            max_r = min(self.automaton.grid_size - 1, center_row + brush_radius)
            min_c = max(0, center_col - brush_radius)
            max_c = min(self.automaton.grid_size - 1, center_col + brush_radius)

            for r_idx in range(min_r, max_r + 1):
                for c_idx in range(min_c, max_c + 1):
                    if (r_idx - center_row) ** 2 + (c_idx - center_col) ** 2 <= brush_radius_sq:
                        original_state_at_cell = self.automaton.grid[r_idx, c_idx]

                        can_paint_cell = True
                        if target_state == CellularAutomaton.VACCINATED:
                            if original_state_at_cell != CellularAutomaton.SUSCEPTIBLE:
                                can_paint_cell = False
                        elif original_state_at_cell == CellularAutomaton.NEUTRAL and \
                                target_state != CellularAutomaton.SUSCEPTIBLE:
                            can_paint_cell = False # Do not paint over NEUTRAL cells unless erasing

                        if can_paint_cell and original_state_at_cell != target_state:
                            self.automaton.grid[r_idx, c_idx] = target_state
                            changed_anything = True
                            if target_state == CellularAutomaton.INFECTED:
                                self.automaton.infection_time[r_idx, c_idx] = self.automaton.iteration
                            elif original_state_at_cell == CellularAutomaton.INFECTED and \
                                    target_state == CellularAutomaton.SUSCEPTIBLE:
                                self.automaton.infection_time[r_idx, c_idx] = 0

        if changed_anything and self.grid_image_artist:
            self.grid_image_artist.set_data(self.automaton.get_colored_grid())
            if is_motion and hasattr(self.grid_canvas, 'supports_blit') and self.grid_canvas.supports_blit:
                self.grid_ax.draw_artist(self.grid_image_artist)
                self.grid_canvas.blit(self.grid_ax.bbox)
            else:
                self.grid_canvas.draw_idle()

        return changed_anything

    def _on_grid_press(self, event):
        if event.inaxes != self.grid_ax or not self.automaton or event.button != 1: return
        if self.automaton.grid_size <= 0 or self.is_playing: return
        self.is_drawing_on_grid = True
        self._last_painted_center_cell = None  # Сбрасываем для нового нажатия
        self.last_drawn_cells.clear()

        if self._paint_at_event_location(event, is_motion=False):  # is_motion=False для одиночного клика
            pass  # Обновление теперь внутри _paint_at_event_location (draw_idle)

    def _on_grid_motion(self, event):
        if not self.is_drawing_on_grid or event.inaxes != self.grid_ax or \
                not self.automaton or self.automaton.grid_size <= 0 or self.is_playing:
            self._last_painted_center_cell = None  # Сбрасываем, если вышли за пределы или не рисуем
            return


        self._paint_at_event_location(event, is_motion=True)

    def _on_grid_release(self, event):
        if event.button != 1 or not self.is_drawing_on_grid: return
        self.is_drawing_on_grid = False
        self._last_painted_center_cell = None  # Сброс после отпускания

        # Финализация нужна, если что-то было изменено во время всего жеста рисования (press + motion)
        # Проверяем, были ли какие-либо изменения в сетке (можно добавить флаг в _paint_at_event_location)
        # или просто вызываем финализацию, она проверит статистику.
        if self.automaton and not self.is_playing:
            self._finalize_manual_grid_change()  # Финализация всегда при отпускании, если рисовали
        self.last_drawn_cells.clear()  # Очищаем в любом случае

    def _finalize_manual_grid_change(self):
        if not self.automaton: return
        self.automaton.update_current_snapshot_stats_and_history()
        current_grid = self.automaton.get_grid()
        num_infected, num_vaccinated = 0, 0
        if self.automaton.grid_size > 0:
            num_infected = np.sum(current_grid == CellularAutomaton.INFECTED)
            num_vaccinated = np.sum(current_grid == CellularAutomaton.VACCINATED)
        self.initial_infected_var.set(num_infected)
        self.initial_vaccinated_var.set(num_vaccinated)
        if self.automaton.iteration == 0:
            self.params.initial_infected = num_infected
            self.params.initial_vaccinated_count = num_vaccinated

        self._update_grid_visualization()  # Полная перерисовка после финализации
        self._update_stats_visualization()

    def _update_grid_visualization(self):
        # Определите, требуется ли полная перерисовка из-за изменения фона.
        force_full_redraw_for_background = hasattr(self, '_background_changed_flag') and self._background_changed_flag

        # Полная перерисовка при изменении галочки, если артист еще не создан, или если фон изменился
        needs_full_redraw = (
            self.grid_image_artist is None or
            (hasattr(self, '_grid_lines_visibility_changed_flag') and self._grid_lines_visibility_changed_flag) or
            force_full_redraw_for_background
        )

        if needs_full_redraw:
            self.grid_ax.clear()
            self.grid_image_artist = None # Reset cell artist

            # Сброс флагов, вызвавших полную перерисовку
            if hasattr(self, '_grid_lines_visibility_changed_flag') and self._grid_lines_visibility_changed_flag:
                self._grid_lines_visibility_changed_flag = False
            if force_full_redraw_for_background: # Сбросить фоновый флаг, если он был триггером
                 self._background_changed_flag = False

            current_grid_size = self.params.grid_size if self.params else 0
            common_extent = (-0.5, current_grid_size - 0.5, current_grid_size - 0.5, -0.5)

            if self.background_image_pil and self.background_image_visible.get() and current_grid_size > 0:
                if self.processed_background_image_mpl is None or \
                        self.processed_background_image_mpl.shape[0] != current_grid_size or \
                        self.processed_background_image_mpl.shape[1] != current_grid_size:
                    try:
                        resized_pil_img = self.background_image_pil.resize((current_grid_size, current_grid_size),
                                                                           Image.Resampling.LANCZOS)
                        if resized_pil_img.mode not in ['RGB', 'RGBA']: resized_pil_img = resized_pil_img.convert(
                            'RGBA')
                        img_array = np.array(resized_pil_img)
                        if img_array.dtype in [np.float32, np.float64] and img_array.max() > 1.0:
                            img_array = img_array / 255.0
                        elif img_array.dtype != np.uint8 and not (img_array.dtype in [np.float32, np.float64]):
                            if img_array.max() <= 255 and img_array.min() >= 0: img_array = img_array.astype(np.uint8)
                        self.processed_background_image_mpl = img_array
                    except Exception as e:
                        print(f"Ошибка обработки фонового изображения: {e}");
                        self.processed_background_image_mpl = None
                if self.processed_background_image_mpl is not None:
                    self.grid_ax.imshow(self.processed_background_image_mpl,
                                        alpha=self.background_image_alpha_var.get(), extent=common_extent,
                                        aspect='equal', zorder=0)

            if self.automaton:
                self.grid_image_artist = self.grid_ax.imshow(self.automaton.get_colored_grid(), extent=common_extent,
                                                             aspect='equal', zorder=1, interpolation='nearest')
                self.grid_ax.set_title(f"Итерация: {self.automaton.iteration}")
            else:
                self.grid_ax.set_title("Инициализация...")

            if self.grid_lines_visible.get() and current_grid_size > 0:
                self.grid_ax.set_xticks(np.arange(-.5, current_grid_size, 1), minor=True)
                self.grid_ax.set_yticks(np.arange(-.5, current_grid_size, 1), minor=True)
                self.grid_ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5, zorder=2)
            else:
                self.grid_ax.grid(False)

            self.grid_ax.set_xticks([]);
            self.grid_ax.set_yticks([])

            if self.background_image_pil and self.background_image_visible.get() and self.background_image_path and current_grid_size > 0:
                filename_text = os.path.basename(self.background_image_path)
                max_len = max(10, current_grid_size // 7 if current_grid_size > 0 else 10)
                if len(filename_text) > max_len: filename_text = filename_text[:max_len - 3] + "..."
                self.grid_ax.text(0.02, 0.98, filename_text, transform=self.grid_ax.transAxes, color='white',
                                  fontsize=8, alpha=0.9,
                                  bbox=dict(facecolor='black', alpha=0.5, boxstyle='round,pad=0.2'), zorder=3,
                                  verticalalignment='top')

            for spine in self.grid_ax.spines.values(): spine.set_visible(current_grid_size > 0)
            self.grid_ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                                     labelbottom=False, labeltop=False, labelleft=False, labelright=False)

            self.grid_canvas.draw_idle()

        elif self.grid_image_artist and self.automaton:  # Оптимизированный путь (только обновление данных клеток)
            self.grid_image_artist.set_data(self.automaton.get_colored_grid())
            self.grid_ax.set_title(f"Итерация: {self.automaton.iteration}")
            self.grid_canvas.draw_idle()

    def _toggle_grid_lines_visibility(self):  # Новая функция для чекбокса
        self._grid_lines_visibility_changed_flag = True  # Устанавливаем флаг
        self._update_grid_visualization()

    def _update_stats_visualization(self):
        self.stats_ax.clear()
        if not self.automaton: self.stats_canvas.draw_idle(); return
        stats = self.automaton.get_statistics()
        iterations = range(len(stats.get('susceptible', [])))

        labels_counts = {
            'susceptible': ('Восприимчивые', 'y-'), 'infected': ('Инфицированные', 'r-'),
            'recovered': ('Выздоровевшие', 'g-'), 'dead': ('Умершие', 'k-'),
            'vaccinated': ('Вакцинированные', 'b-'), 'neutral': ('Нейтральные', 'gray')
        }
        lines = []
        for key, (label_text, style) in labels_counts.items():
            count_str, data_plot = " (0)", []
            if key in stats and stats.get(key):
                count_str = f" ({stats[key][-1]})" if stats[key] else " (0)"  # Handle empty list for key
                data_plot = stats[key]

            if len(data_plot) < len(iterations):
                data_plot.extend([data_plot[-1] if data_plot else 0] * (len(iterations) - len(data_plot)))
            elif len(data_plot) > len(iterations) and iterations:
                data_plot = data_plot[:len(iterations)]

            line_data = (iterations, data_plot) if iterations else ([], [])
            line, = self.stats_ax.plot(*line_data, style, label=f'{label_text}{count_str}')
            lines.append(line)

        self.stats_ax.set_xlabel('Итерация');
        self.stats_ax.set_ylabel('Количество')
        if lines: self.stats_ax.legend(handles=lines, loc='upper center', bbox_to_anchor=(0.5, -0.20), fancybox=True,
                                       shadow=True, ncol=3, fontsize='small')
        self.stats_fig.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.30)
        self.stats_ax.grid(True);

        if self.stats_canvas.supports_blit and hasattr(self, 'stats_ax_background') and self.stats_ax_background:
            self.stats_fig.canvas.draw()  # Полная перерисовка для обновления фона
            self.stats_ax_background = self.stats_fig.canvas.copy_from_bbox(self.stats_ax.bbox)
        else:
            self.stats_canvas.draw_idle()

    def _on_stats_hover(self, event):
        if self.is_playing:
            # Если симуляция воспроизводится, отключаем логику аннотаций.
            # Если аннотация осталась от предыдущего состояния паузы,
            # просто очищаем ссылку на нее. Анимация сама очистит холст.
            if self.stats_annot is not None:
                # Мы не будем вызывать .remove() или .set_visible(False) здесь,
                # так как объект аннотации может быть уже в "невалидном" состоянии
                # из-за ax.clear() в цикле анимации. Просто обнуляем ссылку.
                self.stats_annot = None
                self.last_annot_info = None
            return  # Выходим, не обрабатывая далее событие наведения

        ax = self.stats_ax
        # Используем флаг, чтобы избежать нескольких последовательных перерисовок, если возможно,
        # хотя в данной структуре с break и return, draw_idle() вызывается по месту.
        # Для строгого следования шаблону пользователя, draw_idle() вызывается сразу.

        if event.inaxes == ax:  # Курсор внутри осей
            found_target_for_annot = False
            for line in ax.get_lines():
                if not line.get_visible() or not line.get_xdata().size > 0:  # Проверка на видимость и наличие данных
                    continue

                contains, P_info = line.contains(event)
                if contains and P_info['ind'].size > 0:  # Убедимся, что P_info['ind'] не пуст
                    x_data_line, y_data_line = line.get_data()
                    min_dist_sq_to_event = float('inf')
                    actual_idx_on_line = -1
                    for i_in_ind in P_info['ind']:
                        if i_in_ind < len(x_data_line):  # Проверка границ
                            dist_sq = (x_data_line[i_in_ind] - event.xdata) ** 2
                            if dist_sq < min_dist_sq_to_event:
                                min_dist_sq_to_event = dist_sq
                                actual_idx_on_line = i_in_ind

                    if actual_idx_on_line == -1:  # Если по какой-то причине не нашли, пропускаем
                        continue

                    idx = actual_idx_on_line

                    # Если аннотация уже для этой точки, ничего не делаем
                    if self.stats_annot and self.last_annot_info == (line, idx):
                        return  # Ничего не изменилось, перерисовка не нужна

                    # Удаляем старую аннотацию, если она была
                    if self.stats_annot:
                        try:
                            self.stats_annot.remove()
                        except Exception:
                            # Игнорируем, если художник уже удален или произошла другая ошибка
                            pass
                        self.stats_annot = None  # Убедимся, что ссылка очищена

                    point_x, point_y = x_data_line[idx], y_data_line[idx]

                    label_base = line.get_label()
                    if '(' in label_base: label_base = label_base.split(' (')[0]

                    annotation_text = f"{label_base}\nИтерация: {int(point_x)}\nЗначение: {int(point_y)}"

                    self.stats_annot = ax.annotate(
                        annotation_text, xy=(point_x, point_y),
                        xytext=(10, -20), textcoords="offset points",  # Используем фиксированное смещение
                        bbox=dict(boxstyle="round,pad=0.3", fc="lemonchiffon", alpha=0.85, lw=0.5),
                        fontsize='small'
                    )
                    self.last_annot_info = (line, idx)
                    self.stats_canvas.draw_idle()  # Перерисовываем холст
                    found_target_for_annot = True
                    break  # Нашли линию и обработали, выходим из цикла по линиям

            if not found_target_for_annot and self.stats_annot:
                # Курсор все еще в осях, но не на линии (например, отодвинули от линии),
                # а аннотация осталась от предыдущего положения на линии. Удаляем её.
                try:
                    self.stats_annot.remove()
                except Exception:
                    pass
                self.stats_annot = None
                self.last_annot_info = None
                self.stats_canvas.draw_idle()  # Перерисовываем холст

        else:  # Курсор вне осей
            if self.stats_annot:
                # Аннотация есть, но курсор ушел с графика, удаляем аннотацию.
                try:
                    self.stats_annot.remove()
                except Exception:
                    pass
                self.stats_annot = None
                self.last_annot_info = None
                self.stats_canvas.draw_idle()  # Перерисовываем холст

    def _step_forward(self):
        if not self.automaton or self.automaton.grid_size <= 0:
            messagebox.showwarning("Предупреждение", "Размер сетки 0 или симуляция не инициализирована.");
            return
        current_stats = self.automaton.get_statistics()
        epidemic_ended = self.automaton.iteration > 0 and current_stats.get('infected') and current_stats[
            'infected'] and current_stats['infected'][-1] == 0 and any(c > 0 for c in current_stats['infected'][:-1])

        num_infected_on_grid = np.sum(self.automaton.grid == CellularAutomaton.INFECTED)
        if num_infected_on_grid == 0 and not epidemic_ended:
            messagebox.showwarning("Предупреждение",
                                   "Нет инфицированных для продолжения. Нарисуйте или инициализируйте инфицированных.")
            return
        if epidemic_ended:
            messagebox.showwarning("Предупреждение", "Эпидемия уже завершилась.")
            return

        if self.is_playing: self._toggle_play(force_pause=True)
        self.automaton.step()
        self._update_grid_visualization();
        self._update_stats_visualization()
        stats_after = self.automaton.get_statistics()
        if self.automaton.iteration > 0 and stats_after.get('infected') and stats_after['infected'] and \
                stats_after['infected'][-1] == 0:
            was_infected_previously = any(c > 0 for c in stats_after['infected'][:-1]) if len(
                stats_after['infected']) > 1 else (self.params.initial_infected > 0)
            if was_infected_previously:
                if hasattr(self.animation_obj,
                           'event_source') and self.animation_obj.event_source: self.animation_obj.event_source.stop()
                self.is_playing = False;
                self.play_button.config(text="Воспроизвести")
                messagebox.showinfo("Симуляция завершена", "Эпидемия закончилась (нет инфицированных).")

    def _toggle_play(self, force_pause=False):
        if not self.automaton: return
        if self.automaton.grid_size <= 0 and not force_pause:
            messagebox.showwarning("Предупреждение", "Размер сетки 0, симуляция не может быть запущена.");
            return

        if not force_pause and not self.is_playing:
            current_stats = self.automaton.get_statistics()
            no_infected_on_grid = np.sum(self.automaton.grid == CellularAutomaton.INFECTED) == 0
            epidemic_ended = self.automaton.iteration > 0 and current_stats.get('infected') and current_stats[
                'infected'] and current_stats['infected'][-1] == 0 and any(
                c > 0 for c in current_stats['infected'][:-1])
            if no_infected_on_grid and not epidemic_ended:
                messagebox.showwarning("Предупреждение",
                                       "Нет инфицированных для начала. Нарисуйте или инициализируйте.");
                return
            if epidemic_ended:
                messagebox.showinfo("Информация", "Эпидемия уже завершилась. Сбросьте для нового запуска.");
                return

        self.is_playing = not self.is_playing if not force_pause else False
        self.play_button.config(text="Пауза" if self.is_playing else "Воспроизвести")
        if self.drawing_tools_frame:
            for widget in self.drawing_tools_frame.winfo_children():
                if isinstance(widget, (ttk.Radiobutton, ttk.Spinbox, ttk.LabelFrame)):
                    try:
                        widget.config(state=tk.DISABLED if self.is_playing else tk.NORMAL)
                        if isinstance(widget, ttk.LabelFrame):
                            for sub_widget in widget.winfo_children():
                                if isinstance(sub_widget, ttk.Spinbox): sub_widget.config(
                                    state=tk.DISABLED if self.is_playing else tk.NORMAL)
                    except tk.TclError:
                        pass
        if hasattr(self.animation_obj, 'event_source') and self.animation_obj.event_source:
            if self.is_playing:
                self.animation_obj.event_source.start()
            else:
                self.animation_obj.event_source.stop()

    def _reset_simulation(self):
        self.params.initial_infected = self.initial_infected_var.get()
        self.params.initial_vaccinated_count = self.initial_vaccinated_var.get()
        self.params.grid_size = self.grid_size_var.get()
        self.params.vaccinated_infection_probability = self.vaccinated_infection_prob_var.get()
        self._init_simulation()

    def _set_speed(self, speed_value):
        self.speed = speed_value
        if self.animation_obj and hasattr(self.animation_obj, 'event_source') and self.animation_obj.event_source:
            is_playing_orig = self.is_playing
            if is_playing_orig: self.animation_obj.event_source.stop()
            self.animation_obj.event_source.interval = self.speed
            if is_playing_orig: self.animation_obj.event_source.start()

    def _save_overall_stats(self):
        if not self.automaton or not self.automaton.statistics.get('susceptible'): messagebox.showinfo("Информация",
                                                                                                       "Нет данных для сохранения."); return
        filepath = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")],
                                                title="Сохранить общую статистику",
                                                initialfile="overall_epidemic_statistics.png")
        if not filepath: return
        fig_save, ax_save = plt.subplots(figsize=(10, 7));
        stats = self.automaton.get_statistics()
        iterations = range(len(stats.get('susceptible', [])))
        if not iterations: plt.close(fig_save); messagebox.showinfo("Информация",
                                                                    "Нет итераций для сохранения статистики."); return
        labels_data = {'susceptible': ('Восприимчивые', 'y-'), 'infected': ('Инфицированные', 'r-'),
                       'recovered': ('Выздоровевшие', 'g-'), 'dead': ('Умершие', 'k-'),
                       'vaccinated': ('Вакцинированные', 'b-'), 'neutral': ('Нейтральные', 'gray')}
        lines_save = []
        for key, (label, style) in labels_data.items():
            count_str, data_plot = " (0)", []
            if key in stats and stats.get(key): count_str, data_plot = f" ({stats[key][-1]})" if stats[key] else " (0)", \
            stats[key]
            if len(data_plot) < len(iterations):
                data_plot.extend([data_plot[-1] if data_plot else 0] * (len(iterations) - len(data_plot)))
            elif len(data_plot) > len(iterations):
                data_plot = data_plot[:len(iterations)]
            line, = ax_save.plot(iterations, data_plot, style, linewidth=2, label=f"{label}{count_str}");
            lines_save.append(line)
        ax_save.set_title('Динамика эпидемии (Общая статистика)', fontsize=16);
        ax_save.set_xlabel('Итерация', fontsize=14);
        ax_save.set_ylabel('Количество', fontsize=14)
        ax_save.legend(handles=lines_save, loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True,
                       ncol=3, fontsize=10)
        fig_save.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.25);
        ax_save.grid(True)
        try:
            fig_save.savefig(filepath, dpi=300); messagebox.showinfo("Сохранение",
                                                                     f"Общая статистика сохранена в {filepath}")
        except Exception as e:
            messagebox.showerror("Ошибка сохранения", f"Не удалось сохранить файл: {e}")
        finally:
            plt.close(fig_save)

    def _save_individual_plots(self):
        if not self.automaton: messagebox.showinfo("Информация", "Симуляция не инициализирована."); return
        if self.automaton.grid_size <= 0: messagebox.showinfo("Информация",
                                                              "Сетка пуста, нет графиков для сохранения."); return
        save_dir = filedialog.askdirectory(title="Выберите папку для сохранения текущих графиков")
        if not save_dir: return
        try:
            self.grid_fig.savefig(os.path.join(save_dir, f"current_grid_iter_{self.automaton.iteration}.png"), dpi=300,
                                  bbox_inches='tight')
            if self.automaton.statistics.get('susceptible'): self.stats_fig.savefig(
                os.path.join(save_dir, f"current_stats_iter_{self.automaton.iteration}.png"), dpi=300,
                bbox_inches='tight')
            messagebox.showinfo("Сохранение", f"Текущие графики сохранены в {save_dir}")
        except Exception as e:
            messagebox.showerror("Ошибка сохранения", f"Не удалось сохранить файлы: {e}")

    def _save_gif_animation(self):
        if not self.automaton or not self.automaton.grid_history: messagebox.showinfo("Информация",
                                                                                      "Нет данных для GIF."); return
        if self.automaton.grid_size <= 0: messagebox.showinfo("Информация",
                                                              "Сетка пуста, нет анимации для сохранения."); return
        filepath = filedialog.asksaveasfilename(defaultextension=".gif", filetypes=[("GIF files", "*.gif")],
                                                title="Сохранить GIF анимацию", initialfile="epidemic_iterations.gif")
        if not filepath: return
        history = self.automaton.grid_history
        if not history: messagebox.showinfo("Информация", "Нет кадров для GIF."); return
        fig_gif, ax_gif = plt.subplots(figsize=self.grid_fig.get_size_inches());
        ax_gif.axis('off');
        ims = []
        bg_img_for_gif, bg_alpha_for_gif = None, self.background_image_alpha_var.get()
        current_grid_size = self.automaton.grid_size
        common_extent = (-0.5, current_grid_size - 0.5, current_grid_size - 0.5, -0.5)
        draw_grid_lines_for_gif = self.grid_lines_visible.get()

        if self.background_image_pil and self.background_image_visible.get() and current_grid_size > 0:
            try:
                if self.processed_background_image_mpl is not None and self.processed_background_image_mpl.shape[
                    0] == current_grid_size and self.processed_background_image_mpl.shape[1] == current_grid_size:
                    bg_img_for_gif = self.processed_background_image_mpl
                else:
                    resized_pil = self.background_image_pil.resize((current_grid_size, current_grid_size),
                                                                   Image.Resampling.LANCZOS)
                    if resized_pil.mode not in ['RGB', 'RGBA']: resized_pil = resized_pil.convert('RGBA')
                    bg_img_for_gif = np.array(resized_pil)
                    if bg_img_for_gif.dtype in [np.float32, np.float64] and bg_img_for_gif.max() > 1.0:
                        bg_img_for_gif /= 255.0
                    elif bg_img_for_gif.dtype != np.uint8 and not (bg_img_for_gif.dtype in [np.float32, np.float64]):
                        if bg_img_for_gif.max() <= 255 and bg_img_for_gif.min() >= 0: bg_img_for_gif = bg_img_for_gif.astype(
                            np.uint8)
            except Exception as e:
                print(f"Ошибка подготовки фона для GIF: {e}"); bg_img_for_gif = None

        filename_on_gif_text = None
        if bg_img_for_gif is not None and self.background_image_path and current_grid_size > 0:
            base_fn = os.path.basename(self.background_image_path)
            max_len = max(10, current_grid_size // 7 if current_grid_size > 0 else 10)
            if len(base_fn) > max_len: base_fn = base_fn[:max_len - 3] + "..."
            filename_on_gif_text = base_fn

        for i, grid_snapshot in enumerate(history):
            if grid_snapshot.shape[0] != current_grid_size or grid_snapshot.shape[1] != current_grid_size: continue
            frame_artists = []
            if bg_img_for_gif is not None: frame_artists.append(
                ax_gif.imshow(bg_img_for_gif, alpha=bg_alpha_for_gif, extent=common_extent, aspect='equal', zorder=0,
                              animated=True))
            frame_artists.append(
                ax_gif.imshow(self.automaton._grid_to_colored_grid(grid_snapshot), extent=common_extent, aspect='equal',
                              zorder=1, animated=True, interpolation='nearest'))
            if draw_grid_lines_for_gif and current_grid_size > 0:
                xticks_m, yticks_m = np.arange(-.5, current_grid_size, 1), np.arange(-.5, current_grid_size, 1)
                for x_v in xticks_m: frame_artists.append(
                    ax_gif.plot([x_v, x_v], [yticks_m[0], yticks_m[-1]], "gray", ls='-', lw=0.5, zorder=2,
                                animated=True)[0])
                for y_v in yticks_m: frame_artists.append(
                    ax_gif.plot([xticks_m[0], xticks_m[-1]], [y_v, y_v], "gray", ls='-', lw=0.5, zorder=2,
                                animated=True)[0])
            frame_artists.append(
                ax_gif.text(0.5, 1.01, f"Итерация: {i}", ha="center", va="bottom", transform=ax_gif.transAxes,
                            fontsize="small", animated=True))
            if filename_on_gif_text: frame_artists.append(
                ax_gif.text(0.02, 0.98, filename_on_gif_text, transform=ax_gif.transAxes, color='w', fontsize=8,
                            alpha=0.9, bbox={'fc': 'k', 'alpha': 0.5, 'boxstyle': 'round,pad=0.2'}, zorder=3, va='top',
                            animated=True))
            ax_gif.set_xlim(common_extent[0], common_extent[1]);
            ax_gif.set_ylim(common_extent[2], common_extent[3])
            ims.append(frame_artists)

        if not ims: messagebox.showinfo("Информация", "Не удалось создать кадры GIF."); plt.close(fig_gif); return
        ani_gif = animation.ArtistAnimation(fig_gif, ims, interval=self.speed, blit=True, repeat_delay=1000)
        progress_bar = None
        try:
            progress_bar = ttk.Progressbar(self.root, mode='determinate', maximum=len(ims));
            progress_bar.place(relx=0.5, rely=0.5, anchor=tk.CENTER, width=300);
            self.root.update_idletasks()

            def upd_prog(cur, tot):
                if progress_bar and progress_bar.winfo_exists(): progress_bar.config(
                    value=cur + 1); self.root.update_idletasks()

            ani_gif.save(filepath, writer='pillow', fps=max(1, min(30, int(1000 / self.speed))),
                         progress_callback=upd_prog)
            messagebox.showinfo("Сохранение", f"GIF анимация сохранена в {filepath}")
        except Exception as e:
            messagebox.showerror("Ошибка сохранения GIF", f"Не удалось сохранить GIF: {e}")
        finally:
            if progress_bar and progress_bar.winfo_exists(): progress_bar.destroy()
            plt.close(fig_gif)

    def _save_all_iterations_images(self):
        if not self.automaton or not self.automaton.grid_history: messagebox.showinfo("Информация",
                                                                                      "Нет истории итераций."); return
        if self.automaton.grid_size <= 0: messagebox.showinfo("Информация", "Сетка пуста, нет изображений."); return
        save_dir = filedialog.askdirectory(title="Выберите папку для сохранения изображений итераций")
        if not save_dir: return
        stats, history_proc = self.automaton.get_statistics(), []
        is_over = self.automaton.iteration > 0 and stats.get('infected') and stats['infected'] and stats['infected'][
            -1] == 0 and any(s > 0 for s in stats['infected'][:-1])
        num_hist = len(self.automaton.grid_history)
        if not self.is_playing:
            history_proc = self.automaton.grid_history if is_over else self.automaton.grid_history[
                                                                       :min(self.automaton.iteration + 1, num_hist)]
            if not history_proc: messagebox.showinfo("Информация", "Нет итераций для сохранения."); return
        else:
            messagebox.showinfo("Информация", "Остановите симуляцию для сохранения изображений."); return

        fig_iter, ax_iter = plt.subplots(figsize=self.grid_fig.get_size_inches());
        ax_iter.axis('off')
        bg_img_for_iter, bg_alpha_for_iter = None, self.background_image_alpha_var.get()
        current_grid_size = self.automaton.grid_size
        common_extent = (-0.5, current_grid_size - 0.5, current_grid_size - 0.5, -0.5)
        draw_grid_lines_for_iter = self.grid_lines_visible.get()

        if self.background_image_pil and self.background_image_visible.get() and current_grid_size > 0:
            try:
                if self.processed_background_image_mpl is not None and self.processed_background_image_mpl.shape[
                    0] == current_grid_size and self.processed_background_image_mpl.shape[1] == current_grid_size:
                    bg_img_for_iter = self.processed_background_image_mpl
                else:
                    resized_pil = self.background_image_pil.resize((current_grid_size, current_grid_size),
                                                                   Image.Resampling.LANCZOS)
                    if resized_pil.mode not in ['RGB', 'RGBA']: resized_pil = resized_pil.convert('RGBA')
                    bg_img_for_iter = np.array(resized_pil)
                    if bg_img_for_iter.dtype in [np.float32, np.float64] and bg_img_for_iter.max() > 1.0:
                        bg_img_for_iter /= 255.0
                    elif bg_img_for_iter.dtype != np.uint8 and not (bg_img_for_iter.dtype in [np.float32, np.float64]):
                        if bg_img_for_iter.max() <= 255 and bg_img_for_iter.min() >= 0: bg_img_for_iter = bg_img_for_iter.astype(
                            np.uint8)
            except Exception as e:
                print(f"Ошибка подготовки фона для изображений: {e}"); bg_img_for_iter = None

        filename_on_iter_text = None
        if bg_img_for_iter is not None and self.background_image_path and current_grid_size > 0:
            base_fn = os.path.basename(self.background_image_path)
            max_len = max(10, current_grid_size // 7 if current_grid_size > 0 else 10)
            if len(base_fn) > max_len: base_fn = base_fn[:max_len - 3] + "..."
            filename_on_iter_text = base_fn

        progress_bar = None
        try:
            if not history_proc: messagebox.showinfo("Информация", "Нет кадров для сохранения."); plt.close(
                fig_iter); return
            progress_bar = ttk.Progressbar(self.root, mode='determinate', maximum=len(history_proc));
            progress_bar.place(relx=0.5, rely=0.5, anchor=tk.CENTER, width=300);
            self.root.update_idletasks()
            num_saved = 0
            for i, grid_snap in enumerate(history_proc):
                if grid_snap.shape[0] != current_grid_size or grid_snap.shape[1] != current_grid_size:
                    if progress_bar and progress_bar.winfo_exists(): progress_bar[
                        'value'] = i + 1; self.root.update_idletasks()
                    continue
                ax_iter.clear()
                if bg_img_for_iter is not None: ax_iter.imshow(bg_img_for_iter, alpha=bg_alpha_for_iter,
                                                               extent=common_extent, aspect='equal', zorder=0)
                ax_iter.imshow(self.automaton._grid_to_colored_grid(grid_snap), extent=common_extent, aspect='equal',
                               zorder=1, interpolation='nearest')
                ax_iter.set_title(f"Итерация: {i}", fontsize="small")
                if draw_grid_lines_for_iter and current_grid_size > 0:
                    ax_iter.set_xticks(np.arange(-.5, current_grid_size, 1), minor=True);
                    ax_iter.set_yticks(np.arange(-.5, current_grid_size, 1), minor=True)
                    ax_iter.grid(which="minor", color="gray", linestyle='-', linewidth=0.5, zorder=2)
                    ax_iter.set_xticks([]);
                    ax_iter.set_yticks([])
                else:
                    ax_iter.grid(False); ax_iter.set_xticks([]); ax_iter.set_yticks([])
                if filename_on_iter_text: ax_iter.text(0.02, 0.98, filename_on_iter_text, transform=ax_iter.transAxes,
                                                       color='w', fontsize=8, alpha=0.9,
                                                       bbox={'fc': 'k', 'alpha': 0.5, 'boxstyle': 'round,pad=0.2'},
                                                       zorder=3, va='top')
                ax_iter.axis('off');
                ax_iter.set_xlim(common_extent[0], common_extent[1]);
                ax_iter.set_ylim(common_extent[2], common_extent[3])
                fig_iter.savefig(os.path.join(save_dir, f"iteration_{i}.png"), dpi=150, bbox_inches='tight');
                num_saved += 1
                if progress_bar and progress_bar.winfo_exists(): progress_bar[
                    'value'] = i + 1; self.root.update_idletasks()
            messagebox.showinfo("Сохранение завершено", f"{num_saved} изображений сохранено в {save_dir}")
        except Exception as e:
            messagebox.showerror("Ошибка сохранения", f"Ошибка при сохранении изображений: {e}")
        finally:
            if progress_bar and progress_bar.winfo_exists(): progress_bar.destroy()
            plt.close(fig_iter)


if __name__ == "__main__":
    root = tk.Tk()
    app = EpidemicSimulationApp(root)
    root.mainloop()