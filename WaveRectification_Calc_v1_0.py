import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
import mplcursors
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from scipy.optimize import newton

def on_closing():
    if messagebox.askokcancel("Salir", "¿Quieres cerrar la aplicación?"):
        plt.close('all')
        root.destroy()
        
# Parámetros constantes
T = 2 * np.pi
f = 60
w = 2 * np.pi * f
n_div_x = 8
n_div_y = 10
fig_width = 10
fig_height = 8
flag_R = 0
canvas = None   # Para la gráfica

# Alpha de 0 a 180 grados (utilizado en cálculos y gráficas)
alpha_vals = np.linspace(0, 180, 180)

def actualizar_campos():
    txt_selected = mode_combo.get()
    mode = op_CB_1.get(txt_selected)
    
    color_label = "black"

    if mode == 1:
        alpha_entry.config(state='normal')
        theta_entry.config(state='disabled')
        R_entry.config(state='normal')
        L_entry.config(state='normal')
        Vmax_entry.config(state='normal')
        Z_entry.config(state='disabled')

    elif mode == 2:
        alpha_entry.config(state='normal')
        theta_entry.config(state='normal')
        R_entry.config(state='disabled')
        L_entry.config(state='disabled')
        Vmax_entry.config(state='normal')
        Z_entry.config(state='normal')

    # Actualizar color de los labels
    alpha_label.config(fg=color_label)
    theta_label.config(fg=color_label if mode != 1 else "#b7b7b7")  
    R_label.config(fg=color_label if mode != 2 else "#b7b7b7")  
    L_label.config(fg=color_label if mode != 2 else "#b7b7b7")  
    Vmax_label.config(fg=color_label)  
    Z_label.config(fg=color_label if mode == 2 else "#b7b7b7")  

def calcular_valores():
    global gamma_vals, i_n_vals, i_rn_vals, deg, i_t, v_t, v_ak_t, alpha, beta_prom, gamma_prom, Vmax
    try:
        plt.close('all')
        error_label.config(text="")
        txt_selected = mode_combo.get()
        mode = op_CB_1.get(txt_selected)
        global flag_R 
        
        if mode == 1:
            alpha = float(alpha_entry.get())
            R = float(R_entry.get())
            L = float(L_entry.get())
            Vmax = float(Vmax_entry.get())
            Z = np.sqrt(pow(R,2)+pow(w*L,2))
            i_b = Vmax/Z
            if L == 0:
                flag_R = 1
            else:
                flag_R = 0
            if R == 0:
                R_L = 0
                theta = 90
            else:
                if flag_R == 0:
                    R_L = R / L
                    XL_R = (w * L) / R
                    theta = np.rad2deg(np.arctan(XL_R))
                else:
                    R_L = 0
                    XL_R = 0
                    theta = 0                    

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
        
        beta = 0
        gamma = 0
        i_n = 0
        i_rn = 0
        deg = np.arange(0,360,0.01)
        
        # Se calcula Beta y Gamma en una función aparte
        beta_prom, gamma_prom = calcular_Beta(alpha, theta)
        
        # Se calcula io_avg, io_rms, i_n e i_rn en una función aparte
        io_avg, io_rms, i_n, i_rn = calcular_IN_IRN(alpha, theta, beta_prom, deg, Vmax, Z, i_b)
        
        # Se calculan los vectores para las gráficas de voltajes y corrientes
        v_t, v_ak_t, i_t = calcular_It_Vt_Vakt(deg, alpha, beta_prom, theta, Vmax, Z, i_b, R_L, R, flag_R)
        vo_avg = (Vmax/(2*np.pi)) * (np.cos(np.deg2rad(alpha)) - np.cos(np.deg2rad(beta_prom)))
        vo_rms = Vmax * np.sqrt(1/(4*np.pi) * (np.deg2rad(beta_prom) - np.deg2rad(alpha) + (np.sin(np.deg2rad(2*alpha)))/2 - (np.sin(np.deg2rad(2*beta_prom)))/2))
        
        # Actualiza los labels de los resultados
        theta_label_res.config(text=f"{theta:.2f} grados")
        beta_label_res.config(text=f"{beta_prom:.2f} grados")
        gamma_label_res.config(text=f"{gamma_prom:.2f} grados")
        in_label_res.config(text=f"{i_n:.4f} A")
        irn_label_res.config(text=f"{i_rn:.4f} A")
        io_avg_label_res.config(text=f"{io_avg:.4f} A")
        io_rms_label_res.config(text=f"{io_rms:.4f} A")
        vo_avg_label_res.config(text=f"{vo_avg:.2f} V")
        vo_rms_label_res.config(text=f"{vo_rms:.2f} V")
        
        # Calcular gamma, IN e IRN para alpha de 0 a 180 grados
        gamma_vals = []
        i_n_vals = []
        i_rn_vals = []
        
        for alpha_i in alpha_vals:
            # Se calcula Beta y Gamma en una función aparte
            beta, gamma = calcular_Beta(alpha_i, theta)
            # Se calcula io_avg, io_rms, i_n e i_rn en una función aparte
            _, _, i_n, i_rn = calcular_IN_IRN(alpha_i, theta, beta, deg, Vmax, Z, i_b)
            gamma_vals.append(gamma)
            i_n_vals.append(i_n)
            i_rn_vals.append(i_rn)
        
        grafica_combo.grid()  # Al pulsar "Calcular" muestra el ComboBox y le label de gráficas
        cambiar_Grafica()
        
    except Exception as e:
        error_label.config(text=f"Error: {e}")

import numpy as np

def calcular_Beta(alpha, theta):
    if theta == 0:
        R_L = float('inf')
    else:
        R_L = w / np.tan(np.deg2rad(theta))

    theta_rad = np.deg2rad(theta)
    alpha_rad = np.deg2rad(alpha)

    # Generar valores de beta
    beta_vals = np.arange(np.pi - np.deg2rad(alpha + 1), (2 * np.pi) - np.deg2rad(alpha - 1), 0.001)

    # Calcular los términos en un solo paso
    term1 = np.sin(beta_vals - theta_rad)
    term2 = np.sin(alpha_rad - theta_rad) * np.exp((-R_L) * (beta_vals - alpha_rad) / w)

    # Encontrar índices donde la diferencia es menor que la tolerancia
    tolerance = 1e-3
    beta_indices = np.where(np.abs(term1 - term2) < tolerance)[0]

    beta_rad = beta_vals[beta_indices]
    beta_deg = []
    gamma = []

    if len(beta_rad) > 0:
        if alpha == 0:
            beta_deg_temp = np.rad2deg(beta_rad[0])
            beta_deg.append(beta_deg_temp)
            gamma.append(beta_deg_temp)
        else:
            for beta in beta_rad:
                beta_deg_temp = np.rad2deg(beta)
                if (theta < 10 and alpha <= 10) or (theta == 90 and 1 <= alpha <= 5):
                    beta_deg.append(beta_deg_temp)
                    gamma.append(beta_deg_temp - alpha)
                    break
                elif beta_deg_temp - alpha > 1:
                    beta_deg.append(beta_deg_temp)
                    gamma.append(beta_deg_temp - alpha)

    beta_prom = np.mean(beta_deg) if beta_deg else 0
    gamma_prom = np.mean(gamma) if gamma else 0

    return beta_prom, gamma_prom

def calcular_IN_IRN(alpha, theta, beta, deg, Vmax, Z, i_b):
    # Inicializar la lista de corriente con ceros
    i_t_temp = np.zeros(len(deg))

    if theta == 0:
        R_L = float('inf')
    else:
        R_L = w / np.tan(np.deg2rad(theta))

    # Conversión de grados a tiempo en segundos
    t = deg * (0.0083333 / 180)

    alpha_rad = np.deg2rad(alpha)
    theta_rad = np.deg2rad(theta)

    # Calcular i_t_temp usando vectorización
    limit_indices = (deg >= alpha) & (deg <= beta)
    if R_L == float('inf'):
        i_t_temp[limit_indices] = (Vmax / Z) * np.sin(w * t[limit_indices])
    else:
        i_t_temp[limit_indices] = (
            i_b * (np.sin(w * t[limit_indices] - theta_rad) - 
                    np.sin(alpha_rad - theta_rad) * np.exp((-R_L) * (w * t[limit_indices] - alpha_rad) / w))
        )

    # Índices para alpha y beta
    alpha_index = np.searchsorted(deg, alpha)
    beta_index = np.searchsorted(deg, beta)

    # Calcular io_avg usando integración numérica
    if alpha_index + 1 < beta_index:
        x = np.deg2rad(deg[alpha_index + 1:beta_index] - deg[alpha_index:beta_index - 1]) / (2 * np.pi)
        io_avg = np.sum(x * (i_t_temp[alpha_index:beta_index - 1] + i_t_temp[alpha_index + 1:beta_index]) / 2)
    else:
        io_avg = 0  # o manejar de otra forma

    # Calcular io_rms
    if alpha_index + 1 < beta_index:
        y_rms = (i_t_temp[alpha_index:beta_index - 1] ** 2 + i_t_temp[alpha_index + 1:beta_index] ** 2) / 2
        io_rms = np.sqrt(np.sum(x * y_rms))
    else:
        io_rms = 0  # O manejar de otra forma

    # Calcular corrientes normalizadas
    i_n = io_avg / i_b
    i_rn = io_rms / i_b

    return io_avg, io_rms, i_n, i_rn

def calcular_It_Vt_Vakt(deg, alpha, beta, theta, Vmax, Z, i_b, R_L, R, flag_R):
    i_t = []
    v_t = []
    v_ak = 1
    v_ak_t = []
    
    alpha_rad = np.deg2rad(alpha)
    theta_rad = np.deg2rad(theta)
    
    for i in range (len(deg)):
        t = deg[i] * (0.0083333/180)
        if deg[i] < alpha or deg[i] > beta:
            i_t.append(0)
            v_t.append(0)
            v_ak_t.append(Vmax * np.sin(w*t))
        else:
            if flag_R == 0:
                if R_L == float('inf'):
                    i_t.append((Vmax/Z) *  np.sin(w*t))
                else:
                    i_t.append((i_b) * (np.sin(w*t - theta_rad) - np.sin(alpha_rad - theta_rad) * np.exp((-R_L) * (w*t - alpha_rad) / w)))
            else:
                i_t.append((Vmax/R) *  np.sin(w*t))
            if deg[i] < 180:
                v_t.append((Vmax-v_ak) * np.sin(w*t))
            else:
                v_t.append((Vmax * np.sin(w*t)))
            v_ak_t.append(v_ak)
            
    return v_t, v_ak_t, i_t



# ------------------------------------------------------------
# Funciones para las gráficas
def cambiar_Grafica():
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
            graficar_volts(deg, i_t, v_t, v_ak_t, alpha, beta_prom, gamma_prom, Vmax)
        
        elif grafica_tipo == 2:  # Gráfica de Gamma vs Alpha
            graficar_gamma(alpha_vals, gamma_vals)
        
        elif grafica_tipo == 3:  # Gráfica de IN vs Alpha
            graficar_IN(alpha_vals, i_n_vals)
        
        elif grafica_tipo == 4:  # Gráfica de IRN vs Alpha
            graficar_IRN(alpha_vals, i_rn_vals)
        
        canvas.draw()

def graficar_volts(deg, i_t, v_t, v_ak_t, alpha, beta, gamma, Vmax):
    # Grafica de Corriente de la carga
    plt.subplot(3, 1, 1)  
    plt.plot(deg, i_t, label='i_t (A)', color='blue')
    plt.xlabel('Degrees (deg)')
    plt.ylabel('Amplitude (A)')
    plt.title('i(t) vs Degrees')
    plt.legend()

    plt.xlim(0, 360)  
    plt.ylim(-max(i_t) - (0.01 * max(i_t)), max(i_t) + (0.01 * max(i_t)))  

    x_step = 360 / n_div_x
    y_step = (2 * max(i_t)) / n_div_y

    plt.xticks(np.arange(0, round(max(deg)) + 1, 45))  
    plt.yticks(np.arange(-max(i_t), max(i_t) + y_step , y_step))    
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray') 

    cursor = mplcursors.cursor(hover=True)
    cursor.connect("add", lambda sel: sel.annotation.set_text(f"({sel.target[0]:.2f}, {sel.target[1]:.2f})"))

    graficar_lineas(alpha, beta, gamma, i_t)
    
    # Grafica de Voltaje de la carga
    plt.subplot(3, 1, 2) 
    plt.plot(deg, v_t, label='v_t (V)')
    plt.xlabel('Degrees (deg)')
    plt.ylabel('Amplitude (V)')
    plt.title('v(t) vs Degrees')
    plt.legend()

    plt.xlim(0, 360)  
    plt.ylim(-Vmax - (0.01 * Vmax), Vmax + (0.01 * Vmax))  

    x_step = 360 / n_div_x
    y_step = (2 * Vmax) / n_div_y

    plt.xticks(np.arange(0, round(max(deg)) + 1, x_step))  
    plt.yticks(np.arange(-Vmax , Vmax + y_step, y_step))    
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray') 

    cursor2 = mplcursors.cursor(hover=True)
    cursor2.connect("add", lambda sel: sel.annotation.set_text(f"({sel.target[0]:.2f}, {sel.target[1]:.2f})"))
    
    graficar_lineas(alpha, beta, gamma, v_t)
    
    # Grafica de Voltaje entre terminales
    plt.subplot(3, 1, 3) 
    plt.plot(deg, v_ak_t, label='v_ak_t (V)', color='orange')
    plt.xlabel('Degrees (deg)')
    plt.ylabel('Amplitude (V)')
    plt.title('v_ak(t) vs Degrees')
    plt.legend()

    plt.xlim(0, 360)  
    plt.ylim(-Vmax - (0.01 * Vmax), Vmax + (0.01 * Vmax))  

    plt.xticks(np.arange(0, round(max(deg)) + 1, x_step))  
    plt.yticks(np.arange(-Vmax , Vmax + y_step, y_step))  
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')  

    cursor3 = mplcursors.cursor(hover=True)
    cursor3.connect("add", lambda sel: sel.annotation.set_text(f"({sel.target[0]:.2f}, {sel.target[1]:.2f})"))
    
    graficar_lineas(alpha, beta, gamma, v_ak_t)

    plt.tight_layout()
    
def graficar_lineas(alpha, beta, gamma, graf_var):
    plt.scatter(alpha, 0, color='#fa3c59', s=50)  
    plt.text(alpha, -0.1 * max(graf_var), f'α ({alpha:.2f}°)', fontsize=12, ha='center', va='top', color='#fa3c59') 
    
    plt.scatter(beta, 0, color='#5cb849', s=50)  
    plt.text(beta, -0.1 * max(graf_var), f'β ({beta:.2f}°)', fontsize=12, ha='center', va='top', color='#5cb849') 
    
    plt.annotate('', xy=(beta, -0.6 * max(graf_var)), xytext=(alpha, -0.6 * max(graf_var)), arrowprops=dict(arrowstyle='->', color='#b4ae25', lw=1))
    plt.annotate('', xy=(alpha, -0.6 * max(graf_var)), xytext=(beta, -0.6 * max(graf_var)), arrowprops=dict(arrowstyle='->', color='#b4ae25', lw=1))
    plt.text((alpha + beta)/2, -0.7*max(graf_var), f'ɣ ({gamma:.2f}°)', fontsize=12, ha='right', va='top', color='#b4ae25')    

def graficar_gamma(alpha_vals, gamma_vals):
    plt.plot(alpha_vals, gamma_vals, label='Alpha vs Gamma', color='blue')
    plt.title('Gráfica de Alpha vs Gamma')
    plt.xlabel('Alpha (grados)')
    plt.ylabel('Gamma (grados)')
    plt.legend()

    plt.xticks(np.arange(0, 180, 10))
    plt.yticks(np.arange(0, 360, 10))
    plt.xlim(0, 181)
    plt.ylim(0, 361) 
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray') 

    cursor = mplcursors.cursor(hover=True)
    cursor.connect("add", lambda sel: sel.annotation.set_text(f"({sel.target[0]:.2f}, {sel.target[1]:.2f})"))

def graficar_IN(alpha_vals, i_n_vals):
    plt.plot(alpha_vals, i_n_vals, label='Alpha vs IN', color='blue')
    plt.title('Gráfica de Alpha vs IN')
    plt.xlabel('Alpha (grados)')
    plt.ylabel('IN')
    plt.legend()

    plt.xticks(np.arange(0, 180, 10))
    plt.yticks(np.arange(0, 1, 0.025))
    plt.xlim(0, 181)
    plt.ylim(0, 1.01) 
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray') 

    cursor = mplcursors.cursor(hover=True)
    cursor.connect("add", lambda sel: sel.annotation.set_text(f"({sel.target[0]:.2f}, {sel.target[1]:.4f})"))

def graficar_IRN(alpha_vals, i_rn_vals):
    plt.plot(alpha_vals, i_rn_vals, label='Alpha vs IRN', color='blue')
    plt.title('Gráfica de Alpha vs IRN')
    plt.xlabel('Alpha (grados)')
    plt.ylabel('IRN')
    plt.legend()

    plt.xticks(np.arange(0, 180, 10))
    plt.yticks(np.arange(0, 1.26, 0.03))
    plt.xlim(0, 181)
    plt.ylim(0, 1.26) 
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray') 

    cursor = mplcursors.cursor(hover=True)
    cursor.connect("add", lambda sel: sel.annotation.set_text(f"({sel.target[0]:.2f}, {sel.target[1]:.4f})"))
    
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
        fig_width, fig_height = [15, 10.5]
    case [1280, 720]:
        fig_width, fig_height = [6, 5]
    case [1366, 768]:
        fig_width, fig_height = [9, 7]
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
R_var = tk.DoubleVar()
L_var = tk.DoubleVar()
Vmax_var = tk.DoubleVar()
Z_var = tk.DoubleVar()

# Etiquetas y campos de entrada con mayor tamaño de fuente
font_settings = ('Century Gothic', 14)

# Labels
alpha_label = tk.Label(root, text="Alpha (grados):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
alpha_label.grid(row=1, column=0, padx=10, pady=10)
alpha_entry = ttk.Entry(root, textvariable=alpha_var, font=font_settings)
alpha_entry.grid(row=1, column=1, padx=10, pady=10)

theta_label = tk.Label(root, text="Theta (grados):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
theta_label.grid(row=2, column=0, padx=10, pady=10)
theta_entry = ttk.Entry(root, textvariable=theta_var, font=font_settings)
theta_entry.grid(row=2, column=1, padx=10, pady=10)

R_label = tk.Label(root, text="Resistencia (R):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
R_label.grid(row=3, column=0, padx=10, pady=10)
R_entry = ttk.Entry(root, textvariable=R_var, font=font_settings)
R_entry.grid(row=3, column=1, padx=10, pady=10)

L_label = tk.Label(root, text="Inductancia (L):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
L_label.grid(row=4, column=0, padx=10, pady=10)
L_entry = ttk.Entry(root, textvariable=L_var, font=font_settings)
L_entry.grid(row=4, column=1, padx=10, pady=10)

Vmax_label = tk.Label(root, text="Vmax (V):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
Vmax_label.grid(row=5, column=0, padx=10, pady=10)
Vmax_entry = ttk.Entry(root, textvariable= Vmax_var, font=font_settings)
Vmax_entry.grid(row=5, column=1, padx=10, pady=10)

Z_label = tk.Label(root, text="Z (R):", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 12))
Z_label.grid(row=6, column=0, padx=10, pady=10)
Z_entry = ttk.Entry(root, textvariable= Z_var, font=font_settings)
Z_entry.grid(row=6, column=1, padx=10, pady=10)

# Botón para calcular
calc_btn = ttk.Button(root, text="Calcular", command=calcular_valores)
calc_btn.grid(row=7, column=0, columnspan=2, padx=10, pady=10)

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

theta_label_txt = ttk.Label(root, text="Theta:", style="Modern.TLabel")
theta_label_txt.grid(row=9, column=0, padx=10, sticky="we")

beta_label_txt = ttk.Label(root, text="Beta:", style="Modern.TLabel")
beta_label_txt.grid(row=10, column=0, padx=10, sticky="we")

gamma_label_txt = ttk.Label(root, text="Gamma:", style="Modern.TLabel")
gamma_label_txt.grid(row=11, column=0, padx=10, sticky="we")

in_label_txt = ttk.Label(root, text="In:", style="Modern.TLabel")
in_label_txt.grid(row=12, column=0, padx=10, sticky="we")

irn_label_txt = ttk.Label(root, text="Irn:", style="Modern.TLabel")
irn_label_txt.grid(row=13, column=0, padx=10, sticky="we")

io_avg_label_txt = ttk.Label(root, text="Io Avg:", style="Modern.TLabel")
io_avg_label_txt.grid(row=14, column=0, padx=10, sticky="we")

io_rms_label_txt = ttk.Label(root, text="Io RMS:", style="Modern.TLabel")
io_rms_label_txt.grid(row=15, column=0, padx=10, sticky="we")

vo_avg_label_txt = ttk.Label(root, text="Vo Avg:", style="Modern.TLabel")
vo_avg_label_txt.grid(row=16, column=0, padx=10, sticky="we")

vo_rms_label_txt = ttk.Label(root, text="Vo RMS:", style="Modern.TLabel")
vo_rms_label_txt.grid(row=17, column=0, padx=10, sticky="we")

#--------

theta_label_res = ttk.Label(root, text="", style="Modern.TLabel")
theta_label_res.grid(row=9, column=1, padx=10, sticky="we")

beta_label_res = ttk.Label(root, text="", style="Modern.TLabel")
beta_label_res.grid(row=10, column=1, padx=10, sticky="we")

gamma_label_res = ttk.Label(root, text="", style="Modern.TLabel")
gamma_label_res.grid(row=11, column=1, padx=10, sticky="we")

in_label_res = ttk.Label(root, text="", style="Modern.TLabel")
in_label_res.grid(row=12, column=1, padx=10, sticky="we")

irn_label_res = ttk.Label(root, text="", style="Modern.TLabel")
irn_label_res.grid(row=13, column=1, padx=10, sticky="we")

io_avg_label_res = ttk.Label(root, text="", style="Modern.TLabel")
io_avg_label_res.grid(row=14, column=1, padx=10, sticky="we")

io_rms_label_res = ttk.Label(root, text="", style="Modern.TLabel")
io_rms_label_res.grid(row=15, column=1, padx=10, sticky="we")

vo_avg_label_res = ttk.Label(root, text="", style="Modern.TLabel")
vo_avg_label_res.grid(row=16, column=1, padx=10, sticky="we")

vo_rms_label_res = ttk.Label(root, text="", style="Modern.TLabel")
vo_rms_label_res.grid(row=17, column=1, padx=10, sticky="we")

error_label = ttk.Label(root, text="")
error_label.grid(row=18, column=0, padx=10, sticky="w")

# ComboBox del modo
op_CB_1 = {
    "R y L" : 1,
    "θ y Z" : 2
}

tk.Label(root, text="Modo:", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 16)).grid(row=0, column=0, padx=10, pady=10)

mode_combo = ttk.Combobox(root, textvariable=mode_var, values=list(op_CB_1.keys()), state="readonly", font=font_settings, text='')
mode_combo.set("R y L")
mode_combo.grid(row=0, column=1, padx=10, pady=10)
mode_combo.bind("<<ComboboxSelected>>", lambda e: actualizar_campos())


# ComboBox del tipo de gráfica
op_CB_2 = {
    "Ninguna" : 0,
    "V-I" : 1,
    "α vs ɣ" : 2,
    "IN" : 3,
    "IRN" : 4,
}

tk.Label(root, text="Gráfica:", bg="#f0f0f0", fg="#404040", font=("Century Gothic", 16)).grid(row=8, column=0, padx=10, pady=10)
grafica_combo = ttk.Combobox(root, textvariable=grafica_var, values=list(op_CB_2.keys()), state="readonly", font=font_settings)
grafica_combo.grid(row=8, column=1, padx=10, pady=10)
grafica_combo.set("Ninguna")
grafica_combo.grid_remove()  # Esconde el combo box al inicio
grafica_combo.bind("<<ComboboxSelected>>", lambda e: cambiar_Grafica())

# Iniciar con el modo predeterminado
actualizar_campos()

root.protocol("WM_DELETE_WINDOW", on_closing)

# Iniciar la interfaz
root.mainloop()