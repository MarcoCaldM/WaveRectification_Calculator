import numpy as np
import matplotlib.pyplot as plt
import mplcursors  # Importar la librería para la interactividad
from scipy.optimize import fsolve
import math 

# Parámetros constantes
f = 60
w = 2 * np.pi * f
Vmax = 100
Z = 100
i_b = Vmax / Z

# Datos de ejemplo para deg e i_t
deg = np.arange(0, 361, 1)  # Asumiendo que deg tiene valores de 0 a 360
i_t = np.random.uniform(0, 10, size=len(deg))  # Valores de corriente de ejemplo

# Función para calcular beta y gamma
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

# Función de cálculo de i_n e i_{rn}
import numpy as np

def calcular_in_irn(alpha, theta, beta):
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
    x = np.deg2rad(deg[alpha_index + 1:beta_index] - deg[alpha_index:beta_index - 1]) / (2 * np.pi)
    io_avg = np.sum(x * (i_t_temp[alpha_index:beta_index - 1] + i_t_temp[alpha_index + 1:beta_index]) / 2)

    # Calcular io_rms
    y_rms = (i_t_temp[alpha_index:beta_index - 1] ** 2 + i_t_temp[alpha_index + 1:beta_index] ** 2) / 2
    io_rms = np.sqrt(np.sum(x * y_rms))

    # Calcular corrientes normalizadas
    i_n = io_avg / i_b
    i_rn = io_rms / i_b

    return i_n, i_rn


# Valores de alpha y theta
alphas = np.arange(0, 180, 1)  
thetas = np.arange(0, 95, 5)  
#thetas = [89.5,89.9]
#thetas = 10, 50

# Almacenar resultados para las gráficas
gamma_values = {theta: [] for theta in thetas}
IN_values = {theta: [] for theta in thetas}
IRN_values = {theta: [] for theta in thetas}

for theta in thetas:
    for alpha in alphas:
        beta, gamma = calcular_Beta(alpha, theta)
        i_n, i_rn = calcular_in_irn(alpha, theta, beta)  # Pasar beta a la función
        gamma_values[theta].append(gamma)
        IN_values[theta].append(i_n)
        IRN_values[theta].append(i_rn)

# Crear gráficas separadas
fig1, ax1 = plt.subplots(figsize=(16, 10))
fig1.subplots_adjust(left=0.04, right=0.99, top=0.97, bottom=0.05)
for theta in thetas:
    ax1.plot(alphas, gamma_values[theta], label=f'θ = {theta}°')
ax1.set_title('Alpha vs Gamma', loc='center')
ax1.annotate('Autores: Marco A. Calderón M. y Alejandro Morales H.', 
             xy=(1, 1), xycoords='axes fraction', 
             fontsize=12, ha='right', va='bottom')
ax1.set_xlabel('Alpha (degrees)')
ax1.set_ylabel('Gamma (degrees)')
ax1.legend()
ax1.grid(True)
#mplcursors.cursor(hover=True) 
ax1.tick_params(axis='both', labelsize=8)
ax1.set_xticks(np.arange(0, 181, 5))  
ax1.set_yticks(np.arange(0, 361, 5))
ax1.set_xlim(0,180)
ax1.set_ylim(0,360)
ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
#plt.tight_layout()
plt.show()

fig2, ax2 = plt.subplots(figsize=(10, 6))
fig2.subplots_adjust(left=0.04, right=0.99, top=0.97, bottom=0.05)
for theta in thetas:
    ax2.plot(alphas, IN_values[theta], label=f'θ = {theta}°')
ax2.set_title('Alpha vs IN', loc='center')
ax2.annotate('Autores: Marco A. Calderón M. y Alejandro Morales H.', 
             xy=(1, 1), xycoords='axes fraction', 
             fontsize=12, ha='right', va='bottom')
ax2.set_xlabel('Alpha (degrees)')
ax2.set_ylabel('IN (A)')
ax2.legend()
ax2.grid(True)
ax2.yaxis.set_major_formatter(plt.FormatStrFormatter('%.4f'))
ax2.tick_params(axis='both', labelsize=8)
ax2.set_xticks(np.arange(0, 181, 5))  
ax2.set_yticks(np.arange(0, 1.0125, 0.0125))
ax2.set_xlim(0,180)
ax2.set_ylim(0,1)
ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
#mplcursors.cursor(hover=True)  # Añadir interactividad a la segunda gráfica
#plt.tight_layout()
plt.show()

fig3, ax3 = plt.subplots(figsize=(10, 6))
fig3.subplots_adjust(left=0.04, right=0.99, top=0.97, bottom=0.05)
for theta in thetas:
    ax3.plot(alphas, IRN_values[theta], label=f'θ = {theta}°')
ax3.set_title('Alpha vs IRN', loc='center')
ax3.annotate('Autores: Marco A. Calderón M. y Alejandro Morales H.', 
             xy=(1, 1), xycoords='axes fraction', 
             fontsize=12, ha='right', va='bottom')
ax3.set_xlabel('Alpha (degrees)')
ax3.set_ylabel('IRN (A)')
ax3.legend()
ax3.grid(True)
ax3.yaxis.set_major_formatter(plt.FormatStrFormatter('%.3f'))
ax3.tick_params(axis='both', labelsize=8)
ax3.set_xticks(np.arange(0, 181, 5))  
ax3.set_yticks(np.arange(0, 1.26, 0.015))
ax3.set_xlim(0,180)
ax3.set_ylim(0,1.245)
ax3.grid(True, which='both', linestyle='--', linewidth=0.5)
#mplcursors.cursor(hover=True)  # Añadir interactividad a la tercera gráfica
#plt.tight_layout()
plt.show()