import numpy as np
import matplotlib.pyplot as plt
import mplcursors

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

def get_Beta(alpha, theta):
    global flag_R
    
    if flag_R == 1 or theta == 0:
        return 180, 180 - alpha
    
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

def get_Current_Parameters(alpha, theta, beta, deg, i_b, sign_type):
    # Inicializar la lista de corriente con ceros
    i_t_temp = np.zeros(len(deg))
    T = 2*np.pi if sign_type == 1 else np.pi
    
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
    limit_indices_2 = (deg >= alpha + 180) & (deg <= beta + 180) 
    
    if R_L == float('inf'):
        i_t_temp[limit_indices] = i_b * np.sin(w * t[limit_indices])
        if sign_type == 2:
            i_t_temp[limit_indices_2] = i_b * np.sin(w * t[limit_indices_2])
    else:
        i_t_temp[limit_indices] = (
            i_b * (np.sin(w * t[limit_indices] - theta_rad) - np.sin(alpha_rad - theta_rad) * np.exp((-R_L) * (w * t[limit_indices] - alpha_rad) / w))
        )
        if sign_type == 2:
            i_t_temp[limit_indices_2] = (
                i_b * (np.sin(w * t[limit_indices_2] - theta_rad) - np.sin(alpha_rad - theta_rad) * np.exp((-R_L) * (w * t[limit_indices_2] - alpha_rad) / w))
            )
        
        
    # Índices para alpha y beta
    alpha_index = np.searchsorted(deg, alpha)
    beta_index = np.searchsorted(deg, beta)
    alpha_pi_index = np.searchsorted(deg, alpha + 180)

    # Calcular io_avg e io_rms usando integración numérica
    if alpha_index + 1 < beta_index:
        if beta < 180 + alpha:
            x = np.deg2rad(deg[alpha_index + 1:beta_index] - deg[alpha_index:beta_index - 1]) / T 
            io_avg = np.sum(x * (i_t_temp[alpha_index:beta_index - 1] + i_t_temp[alpha_index + 1:beta_index]) / 2)    
            y_rms = (i_t_temp[alpha_index:beta_index - 1] ** 2 + i_t_temp[alpha_index + 1:beta_index] ** 2) / 2
            io_rms = np.sqrt(np.sum(x * y_rms))
        else:
            x = np.deg2rad(deg[alpha_index + 1:alpha_pi_index] - deg[alpha_index:alpha_pi_index - 1]) / T
            io_avg = np.sum(x * (i_t_temp[alpha_index:alpha_pi_index - 1] + i_t_temp[alpha_index + 1:alpha_pi_index]) / 2)    
            y_rms = (i_t_temp[alpha_index:alpha_pi_index - 1] ** 2 + i_t_temp[alpha_index + 1:alpha_pi_index] ** 2) / 2
            io_rms = np.sqrt(np.sum(x * y_rms))
    else:
        io_avg = 0
        io_rms = 0
        

    # Calcular corrientes normalizadas
    i_n = io_avg / i_b
    i_rn = io_rms / i_b

    return io_avg, io_rms, i_n, i_rn, i_t_temp

def get_Voltage_Parameters(alpha, beta, deg, Vmax, sign_type):
    # Inicializar la lista de voltajes con ceros
    v_t_temp = np.zeros(len(deg))
    v_ak_temp = np.zeros(len(deg))
    
    # Conversión de grados a tiempo en segundos
    t = deg * (0.0083333 / 180)
    T = 2*np.pi if sign_type == 1 else np.pi
    
    # Calcular v_t_temp usando vectorización
    limit_indices = (deg >= alpha) & (deg <= beta)
    limit_indices_2 = (deg >= alpha + 180) & (deg <= beta + 180) if beta < 180 + alpha else (deg >= alpha + 180) & (deg <= 360 + alpha)
    limit_indices_ak = ~limit_indices
    limit_indices_ak_2 = (deg >= 360 + alpha) if beta < 180 + alpha else (deg >= 360 + alpha)
    
    v_t_temp[limit_indices] = (Vmax-v_ak)*np.sin(w * t[limit_indices])
    
    if sign_type == 2:
        v_t_temp[limit_indices_2] = (Vmax-v_ak)*-np.sin(w * t[limit_indices_2]) 
    
    v_ak_temp[limit_indices_ak] = Vmax * np.sin(w * t[limit_indices_ak])
    v_ak_temp[limit_indices] = v_ak 
    
    if sign_type == 2:
        v_ak_temp[limit_indices_2] = 2 * Vmax * np.sin(w * t[limit_indices_2])
        v_ak_temp[limit_indices_ak_2] = v_ak
    
    # Índices para alpha y beta
    alpha_index = np.searchsorted(deg, alpha)
    beta_index = np.searchsorted(deg, beta)
    alpha_pi_index = np.searchsorted(deg, alpha + 180)
    
    # Calcular io_avg e io_rms usando integración numérica
    if alpha_index + 1 < beta_index:
        if beta < 180 + alpha:
            x = np.deg2rad(deg[alpha_index + 1:beta_index] - deg[alpha_index:beta_index - 1]) / T
            vo_avg = np.sum(x* (v_t_temp[alpha_index:beta_index - 1] + v_t_temp[alpha_index + 1:beta_index]) / 2)
            y_rms_v = (v_t_temp[alpha_index:beta_index - 1] ** 2 + v_t_temp[alpha_index + 1:beta_index] ** 2) / 2
            vo_rms = np.sqrt(np.sum(x * y_rms_v))
        else:
            x = np.deg2rad(deg[alpha_index + 1:alpha_pi_index] - deg[alpha_index:alpha_pi_index - 1]) / T
            vo_avg = np.sum(x* (v_t_temp[alpha_index:alpha_pi_index - 1] + v_t_temp[alpha_index + 1:alpha_pi_index]) / 2)
            y_rms_v = (v_t_temp[alpha_index:alpha_pi_index - 1] ** 2 + v_t_temp[alpha_index + 1:alpha_pi_index] ** 2) / 2
            vo_rms = np.sqrt(np.sum(x * y_rms_v))
    else:
        vo_avg = 0
        vo_rms = 0
        
    return vo_avg, vo_rms, v_t_temp, v_ak_temp

# ------------------------------------------------------------
# Funciones para las gráficas
def plot_volts(deg, i_t, v_t, v_ak_t, alpha, beta, gamma, Vmax):
    min_v_t = min(v_t) if abs(min(v_t)) > max(v_t) else -max(v_t)
    max_v_t = max(v_t) if max(v_t) > abs(min(v_t)) else abs(min(v_t))
    min_v_ak = min(v_ak_t) if abs(min(v_ak_t)) > max(v_ak_t) else -max(v_ak_t)
    max_v_ak = max(v_ak_t) if max(v_ak_t) > abs(min(v_ak_t)) else abs(min(v_ak_t))
    min_i_t = min(i_t) if abs(min(i_t)) > max(i_t) else -max(i_t)
    max_i_t = max(i_t) if max(i_t) > abs(min(i_t)) else abs(min(i_t))
    
    # Grafica de Corriente de la carga
    plt.subplot(3, 1, 1)  
    plt.plot(deg, i_t, label='i_t (A)', color='blue')
    plt.xlabel('Degrees (deg)')
    plt.ylabel('Amplitude (A)')
    plt.title('i(t) vs Degrees')
    plt.legend()

    plt.xlim(0, 360)  
    plt.ylim(min_i_t, max_i_t)  

    x_step = 360 / n_div_x
    y_step = (2 * max_i_t) / n_div_y

    plt.xticks(np.arange(0, round(max(deg)) + 1, 45))  
    plt.yticks(np.arange(min_i_t, max_i_t + y_step , y_step))    
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray') 

    cursor = mplcursors.cursor(hover=True)
    cursor.connect("add", lambda sel: sel.annotation.set_text(f"({sel.target[0]:.2f}, {sel.target[1]:.2f})"))

    plot_lines(alpha, beta, gamma, i_t)
    
    # Grafica de Voltaje de la carga
    plt.subplot(3, 1, 2) 
    plt.plot(deg, v_t, label='v_t (V)')
    plt.xlabel('Degrees (deg)')
    plt.ylabel('Amplitude (V)')
    plt.title('v(t) vs Degrees')
    plt.legend()
    
    plt.xlim(0, 360)  
    plt.ylim(min_v_t, max_v_t)  

    x_step = 360 / n_div_x
    y_step = (2 * max_v_t) / n_div_y

    plt.xticks(np.arange(0, round(max(deg)) + 1, x_step))  
    plt.yticks(np.arange(min_v_t , max_v_t + y_step, y_step))    
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray') 

    cursor2 = mplcursors.cursor(hover=True)
    cursor2.connect("add", lambda sel: sel.annotation.set_text(f"({sel.target[0]:.2f}, {sel.target[1]:.2f})"))
    
    plot_lines(alpha, beta, gamma, v_t)
    
    # Grafica de Voltaje entre terminales
    plt.subplot(3, 1, 3) 
    plt.plot(deg, v_ak_t, label='v_ak_t (V)', color='orange')
    plt.xlabel('Degrees (deg)')
    plt.ylabel('Amplitude (V)')
    plt.title('v_ak(t) vs Degrees')
    plt.legend()

    plt.xlim(0, 360)  
    plt.ylim(min_v_ak, max_v_ak)  
    
    x_step = 360 / n_div_x
    y_step = (2 * max_v_ak) / n_div_y
    
    plt.xticks(np.arange(0, round(max(deg)) + 1, x_step))  
    plt.yticks(np.arange(min_v_ak , max_v_ak + y_step, y_step))  
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')  

    cursor3 = mplcursors.cursor(hover=True)
    cursor3.connect("add", lambda sel: sel.annotation.set_text(f"({sel.target[0]:.2f}, {sel.target[1]:.2f})"))
    
    plot_lines(alpha, beta, gamma, v_ak_t)

    plt.tight_layout()
    
def plot_lines(alpha, beta, gamma, graf_var):
    plt.scatter(alpha, 0, color='#fa3c59', s=50)  
    plt.text(alpha, -0.1 * max(graf_var), f'α ({alpha:.2f}°)', fontsize=12, ha='center', va='top', color='#fa3c59') 
    
    plt.scatter(beta, 0, color='#5cb849', s=50)  
    plt.text(beta, -0.1 * max(graf_var), f'β ({beta:.2f}°)', fontsize=12, ha='center', va='top', color='#5cb849') 
    
    plt.annotate('', xy=(beta, -0.6 * max(graf_var)), xytext=(alpha, -0.6 * max(graf_var)), arrowprops=dict(arrowstyle='->', color='#b4ae25', lw=1))
    plt.annotate('', xy=(alpha, -0.6 * max(graf_var)), xytext=(beta, -0.6 * max(graf_var)), arrowprops=dict(arrowstyle='->', color='#b4ae25', lw=1))
    plt.text((alpha + beta)/2, -0.7*max(graf_var), f'ɣ ({gamma:.2f}°)', fontsize=12, ha='right', va='top', color='#b4ae25')    

def plot_gamma(alpha_vals, gamma_vals):
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

def plot_IN(alpha_vals, i_n_vals):
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

def plot_IRN(alpha_vals, i_rn_vals):
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