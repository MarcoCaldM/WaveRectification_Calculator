from functions_IN_IRN_2 import get_Beta, get_Theta, get_Current_Parameters, get_Voltage_Parameters, plot_gamma, plot_IN, plot_IRN,plot_volts
import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import traceback

def on_closing():
    if messagebox.askokcancel("Salir", "¿Quieres cerrar la aplicación?"):
        plt.close('all')
        root.destroy()
        
# Parámetros constantes
T = 2 * np.pi
f = 60
w = 2 * np.pi * f
v_ak = 0.7
n_div_x = 8
n_div_y = 10
fig_width = 10
fig_height = 8
flag_R = 0
canvas = None   # Para la gráfica

# Alpha de 0 a 180 grados (utilizado en cálculos y gráficas)
alpha_vals = np.linspace(0, 180, 180)

def update_fields():
    txt_selected = mode_combo.get()
    mode = op_CB_1.get(txt_selected)
    
    color_label = "black"

    if mode == 1:
        alpha_entry.config(state='normal')
        theta_entry.config(state='disabled')
        beta_entry.config(state='disabled')
        R_entry.config(state='normal')
        L_entry.config(state='normal')
        Vmax_entry.config(state='normal')
        Z_entry.config(state='disabled')
        
        # Actualizar color de los labels
        alpha_label.config(fg=color_label)
        theta_label.config(fg="#b7b7b7")
        beta_label.config(fg="#b7b7b7")  
        R_label.config(fg=color_label)  
        L_label.config(fg=color_label)  
        Vmax_label.config(fg=color_label)  
        Z_label.config(fg="#b7b7b7")  

    elif mode == 2:
        alpha_entry.config(state='normal')
        theta_entry.config(state='normal')
        beta_entry.config(state='disabled')
        R_entry.config(state='disabled')
        L_entry.config(state='disabled')
        Vmax_entry.config(state='normal')
        Z_entry.config(state='normal')
        
        # Actualizar color de los labels
        alpha_label.config(fg=color_label)
        theta_label.config(fg=color_label)
        beta_label.config(fg="#b7b7b7")  
        R_label.config(fg="#b7b7b7")  
        L_label.config(fg="#b7b7b7")  
        Vmax_label.config(fg=color_label)  
        Z_label.config(fg=color_label) 
        
    elif mode == 3:
        alpha_entry.config(state='normal')
        theta_entry.config(state='disabled')
        beta_entry.config(state='normal')
        R_entry.config(state='normal')
        L_entry.config(state='disabled')
        Vmax_entry.config(state='normal')
        Z_entry.config(state='disabled')
        
        # Actualizar color de los labels
        alpha_label.config(fg=color_label)
        theta_label.config(fg="#b7b7b7")
        beta_label.config(fg=color_label)  
        R_label.config(fg=color_label)  
        L_label.config(fg="#b7b7b7")  
        Vmax_label.config(fg=color_label)  
        Z_label.config(fg="#b7b7b7") 

def change_Graph():
    global canvas   # Necesaria
    
    # Elige el tipo de gráfica
    txt_selected = grafica_combo.get()
    grafica_tipo = op_CB_2.get(txt_selected)
    
    if canvas is not None:  # Si ya existe un canvas previo, destruirlo
            canvas.get_tk_widget().destroy() # Eliminar el widget de la gráfica
            canvas = None  # Limpiar la referencia al canvas
    
    if grafica_tipo != 0:
        # Si la gráfica no está ya mostrada, mostrarla
        if canvas is None:
            # Crea la figura para la gráfica
            fig = plt.figure(figsize=(fig_width, fig_height))
            canvas = FigureCanvasTkAgg(fig, master=root)
            canvas.get_tk_widget().grid(row=0, column=2, rowspan="50", padx=20, pady=20)
        
        # Limpia la figura antes de graficar la nueva
        canvas.figure.clf()
    
        if grafica_tipo == 1:   # Gráfica de i(t), v(t), v_ak(t)
            plot_volts(deg, i_t, v_t, v_ak_t, alpha, beta_prom, gamma_prom, Vmax)
        
        elif grafica_tipo == 2:  # Gráfica de Gamma vs Alpha
            plot_gamma(alpha_vals, gamma_vals)
        
        elif grafica_tipo == 3:  # Gráfica de IN vs Alpha
            plot_IN(alpha_vals, i_n_vals)
        
        elif grafica_tipo == 4:  # Gráfica de IRN vs Alpha
            plot_IRN(alpha_vals, i_rn_vals)
        
        canvas.draw()

def calculate_values():
    global gamma_vals, i_n_vals, i_rn_vals, deg, i_t, v_t, v_ak_t, alpha, beta_prom, gamma_prom, Vmax, v_ak
    try:
        plt.close('all')
        error_label.config(text="")
        sign_type = rectification_type.get()
        txt_selected = mode_combo.get()
        mode = op_CB_1.get(txt_selected)
        global flag_R 
        
        if mode == 1:
            alpha = float(alpha_entry.get())
            R = float(R_entry.get())
            L = float(L_entry.get())/1000
            Vmax = float(Vmax_entry.get())
            Z = np.sqrt(pow(R,2)+pow(w*L,2))
            i_b = Vmax/Z
            
            match L, R:
                case (0, _):  
                    flag_R = 1
                    R_L = 0
                    XL_R = 0
                    theta = 0
                case (_, 0):  
                    flag_R = 0
                    R_L = 0
                    XL_R = 0
                    theta = 90
                case (_, _):  
                    flag_R = 0
                    R_L = R / L
                    XL_R = (w * L) / R
                    theta = np.rad2deg(np.arctan(XL_R))

        elif mode == 2:
            alpha = float(alpha_entry.get())
            theta = float(theta_entry.get())
            Vmax = float(Vmax_entry.get())
            Z = float(Z_entry.get())
            R = 0
            i_b = Vmax/Z
            if theta == 0:
                R_L = float('inf') 
            else:
                R_L = w / np.tan(np.deg2rad(theta))
                
        elif mode == 3:
            alpha = float(alpha_entry.get())
            beta_prom = float(beta_entry.get())
            R = float(R_entry.get())
            Vmax = float(Vmax_entry.get())
            theta = get_Theta(alpha, beta_prom)
            gamma_prom = beta_prom - alpha
            L = (R*np.tan(np.deg2rad(theta))) / (2*np.pi*f)
            Z = np.sqrt(pow(R,2)+pow(w*L,2))
            i_b = Vmax/Z
            R_L = R/L
        
        if mode == 1 or mode == 2: 
            beta = 0
            gamma = 0
            
            # Se calcula Beta y Gamma en una función aparte
            beta_prom, gamma_prom = get_Beta(alpha, theta)
        
        i_n = 0
        i_rn = 0
        deg = np.arange(0,360,0.01) if sign_type == 1 else np.arange(0,450,0.01)
            
        # Se calcula io_avg, io_rms, i_n e i_rn en una función aparte
        io_avg, io_rms, i_n, i_rn, i_t = get_Current_Parameters(alpha, theta, beta_prom, deg, i_b, sign_type)
        vo_avg, vo_rms, v_t, v_ak_t = get_Voltage_Parameters(alpha, beta_prom, deg, Vmax, sign_type)
        
        # Actualiza los labels de los resultados
        L_label_res.config(text=f"{(L*1000):.2f} mH") if mode == 3 else L_label_res.config(text="...")
        theta_label_res.config(text=f"{theta:.2f} grados")
        beta_label_res.config(text=f"{beta_prom:.2f} grados")
        gamma_label_res.config(text=f"{gamma_prom:.2f} grados")
        in_label_res.config(text=f"{i_n:.4f}")
        irn_label_res.config(text=f"{i_rn:.4f}")
        io_avg_label_res.config(text=f"{io_avg:.4f} A") if io_avg > 1 else io_avg_label_res.config(text=f"{io_avg*1000:.2f} mA") 
        io_rms_label_res.config(text=f"{io_rms:.4f} A") if io_rms > 1 else io_rms_label_res.config(text=f"{io_rms*1000:.2f} mA")
        vo_avg_label_res.config(text=f"{vo_avg:.2f} V") if vo_avg > 1 else vo_avg_label_res.config(text=f"{vo_avg*1000:.2f} mV")
        vo_rms_label_res.config(text=f"{vo_rms:.2f} V") if vo_rms > 1 else vo_rms_label_res.config(text=f"{vo_rms*1000:.2f} mV")
        
        # Calcular gamma, IN e IRN para alpha de 0 a 180 grados
        gamma_vals = []
        i_n_vals = []
        i_rn_vals = []
        
        for alpha_i in alpha_vals:
            # Se calcula Beta y Gamma en una función aparte
            beta, gamma = get_Beta(alpha_i, theta)
            # Se calcula io_avg, io_rms, i_n e i_rn en una función aparte
            _, _, i_n, i_rn, _ = get_Current_Parameters(alpha_i, theta, beta, deg, i_b, sign_type)
            gamma_vals.append(gamma)
            i_n_vals.append(i_n)
            i_rn_vals.append(i_rn)
        
        grafica_combo.grid()  # Al pulsar "Calcular" muestra el ComboBox y le label de gráficas
        change_Graph()
        
    except Exception as e:
        error_message = traceback.format_exc()
        #error_label.config(text=f"Error: {e}\nDetalles:\n{error_message}")


# ------------------------------------------------------------
# Funciones y objetos de la interfaz grafica

# Crear ventana principal
root = tk.Tk()
root.title("Calculadora de Beta, Gamma, IN e IRN")
root.resizable(True, True)

screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

match [screen_width, screen_height]:
    case [1920, 1080]:
        fig_width, fig_height = [14.5, 10.5]
    case [1366, 768]:
        fig_width, fig_height = [8.5, 6.5]
    case [1600, 900]:
        fig_width, fig_height = [11, 8]
    case [2560, 1440]:
        fig_width, fig_height = [18, 12]
    case [3840, 2160]: 
        fig_width, fig_height = [20, 14]
    case _:
        fig_width, fig_height = [10, 7] 

root.configure(bg="#f0f0f0")  
estilo = ttk.Style()
estilo.configure("TLabel", background="#f0f0f0", foreground="#404040")  
estilo.configure("TButton", background="#d0d0d0", foreground="#404040")  
estilo.configure("TCombobox", background="#d0d0d0", foreground="#404040")
estilo.configure("TEntry", background="#d0d0d0", foreground="#404040")  

# Variables
mode_var = tk.StringVar(value="1")
grafica_var = tk.StringVar(value="0")
alpha_var = tk.DoubleVar()
theta_var = tk.DoubleVar()
beta_var = tk.DoubleVar()
R_var = tk.DoubleVar()
L_var = tk.DoubleVar()
Vmax_var = tk.DoubleVar()
Z_var = tk.DoubleVar()
rectification_type = tk.IntVar()

# Etiquetas y campos de entrada con mayor tamaño de fuente
font_settings = ('Century Gothic', 14)

# Labels
alpha_label = tk.Label(root, text="Alpha (grados):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
alpha_label.grid(row=4, column=0, padx=10, pady=2)
alpha_entry = ttk.Entry(root, textvariable=alpha_var, font=font_settings)
alpha_entry.grid(row=4, column=1, padx=10, pady=2)

theta_label = tk.Label(root, text="Theta (grados):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
theta_label.grid(row=5, column=0, padx=10, pady=2)
theta_entry = ttk.Entry(root, textvariable=theta_var, font=font_settings)
theta_entry.grid(row=5, column=1, padx=10, pady=2)

beta_label = tk.Label(root, text="β (grados):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
beta_label.grid(row=6, column=0, padx=10, pady=2)
beta_entry = ttk.Entry(root, textvariable= beta_var, font=font_settings)
beta_entry.grid(row=6, column=1, padx=10, pady=2)

R_label = tk.Label(root, text="Resistencia (Ω):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
R_label.grid(row=7, column=0, padx=10, pady=2)
R_entry = ttk.Entry(root, textvariable=R_var, font=font_settings)
R_entry.grid(row=7, column=1, padx=10, pady=2)

L_label = tk.Label(root, text="Inductancia (mH):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
L_label.grid(row=8, column=0, padx=10, pady=2)
L_entry = ttk.Entry(root, textvariable=L_var, font=font_settings)
L_entry.grid(row=8, column=1, padx=10, pady=2)

Vmax_label = tk.Label(root, text="Vmax (V):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
Vmax_label.grid(row=9, column=0, padx=10, pady=2)
Vmax_entry = ttk.Entry(root, textvariable= Vmax_var, font=font_settings)
Vmax_entry.grid(row=9, column=1, padx=10, pady=1)

Z_label = tk.Label(root, text="Z (Ω):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
Z_label.grid(row=10, column=0, padx=10, pady=1)
Z_entry = ttk.Entry(root, textvariable= Z_var, font=font_settings)
Z_entry.grid(row=10, column=1, padx=10, pady=1)

# Botón para calcular
calc_btn = ttk.Button(root, text="Calcular", command=calculate_values)
calc_btn.grid(row=11, column=0, columnspan=2, padx=10, pady=10)

# Etiqueta para mostrar resultados
# Creación de estilos individuales para cada etiqueta
style = ttk.Style() 

# Estilo para las etiquetas de resultados
style.configure("Modern.TLabel", 
                background="#f4f2b9",  
                foreground="#000000",  
                font=("Century Gothic", 12, "bold"), 
                padding=(10, 2),
                #relief="solid",  
                borderwidth=1)

L_label_txt =  ttk.Label(root, text="L:", style="Modern.TLabel")
L_label_txt.grid(row=13, column=0, padx=10, sticky="we")

theta_label_txt = ttk.Label(root, text="Theta:", style="Modern.TLabel")
theta_label_txt.grid(row=14, column=0, padx=10, sticky="we")

beta_label_txt = ttk.Label(root, text="Beta:", style="Modern.TLabel")
beta_label_txt.grid(row=15, column=0, padx=10, sticky="we")

gamma_label_txt = ttk.Label(root, text="Gamma:", style="Modern.TLabel")
gamma_label_txt.grid(row=16, column=0, padx=10, sticky="we")

in_label_txt = ttk.Label(root, text="In:", style="Modern.TLabel")
in_label_txt.grid(row=17, column=0, padx=10, sticky="we")

irn_label_txt = ttk.Label(root, text="Irn:", style="Modern.TLabel")
irn_label_txt.grid(row=18, column=0, padx=10, sticky="we")

io_avg_label_txt = ttk.Label(root, text="Io Avg:", style="Modern.TLabel")
io_avg_label_txt.grid(row=19, column=0, padx=10, sticky="we")

io_rms_label_txt = ttk.Label(root, text="Io RMS:", style="Modern.TLabel")
io_rms_label_txt.grid(row=20, column=0, padx=10, sticky="we")

vo_avg_label_txt = ttk.Label(root, text="Vo Avg:", style="Modern.TLabel")
vo_avg_label_txt.grid(row=21, column=0, padx=10, sticky="we")

vo_rms_label_txt = ttk.Label(root, text="Vo RMS:", style="Modern.TLabel")
vo_rms_label_txt.grid(row=22, column=0, padx=10, sticky="we")

#--------
L_label_res = ttk.Label(root, text="", style="Modern.TLabel")
L_label_res.grid(row=13, column=1, padx=10, sticky="we")

theta_label_res = ttk.Label(root, text="", style="Modern.TLabel")
theta_label_res.grid(row=14, column=1, padx=10, sticky="we")

beta_label_res = ttk.Label(root, text="", style="Modern.TLabel")
beta_label_res.grid(row=15, column=1, padx=10, sticky="we")

gamma_label_res = ttk.Label(root, text="", style="Modern.TLabel")
gamma_label_res.grid(row=16, column=1, padx=10, sticky="we")

in_label_res = ttk.Label(root, text="", style="Modern.TLabel")
in_label_res.grid(row=17, column=1, padx=10, sticky="we")

irn_label_res = ttk.Label(root, text="", style="Modern.TLabel")
irn_label_res.grid(row=18, column=1, padx=10, sticky="we")

io_avg_label_res = ttk.Label(root, text="", style="Modern.TLabel")
io_avg_label_res.grid(row=19, column=1, padx=10, sticky="we")

io_rms_label_res = ttk.Label(root, text="", style="Modern.TLabel")
io_rms_label_res.grid(row=20, column=1, padx=10, sticky="we")

vo_avg_label_res = ttk.Label(root, text="", style="Modern.TLabel")
vo_avg_label_res.grid(row=21, column=1, padx=10, sticky="we")

vo_rms_label_res = ttk.Label(root, text="", style="Modern.TLabel")
vo_rms_label_res.grid(row=22, column=1, padx=10, sticky="we")

error_label = ttk.Label(root, text="")
error_label.grid(row=23, column=0, padx=10, sticky="w")

autor1_label = tk.Label(root, text="Autores:", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
autor1_label.grid(row=24, column=0, padx=10, sticky="w")

autor2_label = tk.Label(root, text="Marco A. Calderón M.", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
autor2_label.grid(row=24, column=1, padx=10, sticky="w")

autor2_label = tk.Label(root, text="Alejandro Morales H.", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
autor2_label.grid(row=25, column=1, padx=10, sticky="w")
# ComboBox del modo
op_CB_1 = {
    "R y L" : 1,
    "θ y Z" : 2,
    "β y R" : 3
}

tk.Label(root, text="Modo:", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 16)).grid(row=2, column=0, padx=10, pady=10)

mode_combo = ttk.Combobox(root, textvariable=mode_var, values=list(op_CB_1.keys()), state="readonly", font=font_settings, text='')
mode_combo.set("R y L")
mode_combo.grid(row=2, column=1, padx=10, pady=10)
mode_combo.bind("<<ComboboxSelected>>", lambda e: update_fields())


# ComboBox del tipo de gráfica
op_CB_2 = {
    "Ninguna" : 0,
    "V-I" : 1,
    "α vs ɣ" : 2,
    "IN" : 3,
    "IRN" : 4,
}

tk.Label(root, text="Gráfica:", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 16)).grid(row=12, column=0, padx=10, pady=10)
grafica_combo = ttk.Combobox(root, textvariable=grafica_var, values=list(op_CB_2.keys()), state="readonly", font=font_settings)
grafica_combo.grid(row=12, column=1, padx=10, pady=10)
grafica_combo.set("Ninguna")
grafica_combo.grid_remove()  # Esconde el combo box al inicio
grafica_combo.bind("<<ComboboxSelected>>", lambda e: change_Graph())

#Radio button para la selección del tipo de rectificación
radio_halfwave_singlephase = tk.Radiobutton(root, text="Media onda monofásico", variable=rectification_type, value=1)
radio_halfwave_singlephase.grid(row=0, column=0, padx=5, pady=5)
radio_halfwave_singlephase.select()

radio_fullwave_singlephase = tk.Radiobutton(root, text="Onda completa monofásico", variable=rectification_type, value=2)
radio_fullwave_singlephase.grid(row=1, column=0, padx=5, pady=5)

radio_halfwave_threephase = tk.Radiobutton(root, text="Media onda trifásico", variable=rectification_type, value=3)
radio_halfwave_threephase.grid(row=0, column=1, padx=5, pady=5)

radio_fullwave_threephase = tk.Radiobutton(root, text="Onda completa trifásico", variable=rectification_type, value=4)
radio_fullwave_threephase.grid(row=1, column=1, padx=5, pady=5)
# Iniciar con el modo predeterminado
update_fields()

root.protocol("WM_DELETE_WINDOW", on_closing)

# Iniciar la interfaz
root.mainloop()