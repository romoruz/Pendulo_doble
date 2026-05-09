import numpy as np
import matplotlib.pyplot as plt

# =============================================================
# 1. MATEMÁTICA ESTRICTA (Exactamente tu código original)
# =============================================================

def derivadas(estado, p):
    theta1, theta2, omega1, omega2 = estado
    m1, m2, l1, l2, g = p['m1'], p['m2'], p['l1'], p['l2'], p['g']
    delta = theta1 - theta2
    cd, sd = np.cos(delta), np.sin(delta)

    M00 = (m1 + m2) * l1
    M01 = m2 * l2 * cd
    M10 = m2 * l1 * cd
    M11 = m2 * l2
    detM = M00 * M11 - M01 * M10

    C0 = m2 * l2 * omega2 * omega2 * sd
    C1 = -m2 * l1 * omega1 * omega1 * sd
    G0 = (m1 + m2) * g * np.sin(theta1)
    G1 = m2 * g * np.sin(theta2)

    rhs0 = -(C0 + G0)
    rhs1 = -(C1 + G1)
    alpha1 = (M11 * rhs0 - M01 * rhs1) / detM
    alpha2 = (-M10 * rhs0 + M00 * rhs1) / detM

    return np.array([omega1, omega2, alpha1, alpha2])

def rk4_modificado(estado, dt, p):
    k1 = derivadas(estado, p)
    k2 = derivadas(estado + dt * (k1 / 3.0), p)
    k3 = derivadas(estado + dt * (-k1 / 3.0 + k2), p)
    k4 = derivadas(estado + dt * (k1 - k2 + k3), p)
    return estado + (dt / 8.0) * (k1 + 3.0 * k2 + 3.0 * k3 + k4)

def integrar(estado_inicial, t_max, dt, p):
    n = int(t_max / dt)
    traj = np.empty((n + 1, 4))
    traj[0] = estado_inicial
    for i in range(n):
        traj[i + 1] = rk4_modificado(traj[i], dt, p)
    return traj

def seccion_poincare(traj):
    th1_wrap = np.arctan2(np.sin(traj[:, 0]), np.cos(traj[:, 0]))
    pts = []
    for i in range(len(th1_wrap) - 1):
        a, b = th1_wrap[i], th1_wrap[i + 1]
        if a == 0.0:
            if traj[i, 2] < 0:
                pts.append((traj[i, 1], traj[i, 3]))
            continue
        if a * b < 0 and abs(a - b) < np.pi:
            if traj[i, 2] < 0:
                t = a / (a - b)
                th2 = traj[i, 1] + t * (traj[i + 1, 1] - traj[i, 1])
                w2 = traj[i, 3] + t * (traj[i + 1, 3] - traj[i, 3])
                th2 = np.arctan2(np.sin(th2), np.cos(th2))
                pts.append((th2, w2))
    return np.array(pts) if pts else np.empty((0, 2))

# =============================================================
# 2. GENERADOR DE GRÁFICAS CORREGIDO PARA LATEX
# =============================================================

def procesar_angulos(theta, omega):
    """ Envuelve ángulos en [-pi, pi] y rompe las líneas falsas de proyección """
    th_wrap = np.arctan2(np.sin(theta), np.cos(theta))
    w = omega.copy()
    
    # Detectar dónde el ángulo brinca de pi a -pi (o viceversa)
    saltos = np.where(np.abs(np.diff(th_wrap)) > np.pi)[0]
    
    # Insertar NaNs para que matplotlib no una esos puntos con una línea
    th_wrap[saltos] = np.nan
    w[saltos] = np.nan
    return th_wrap, w

def generar_graficas(estado_inicial, t_max, dt, params, archivo, titulo_principal):
    print(f"Calculando {titulo_principal}...")
    traj = integrar(estado_inicial, t_max, dt, params)
    
    # Procesamiento topológico correcto (evita líneas falsas que cruzan la gráfica)
    th1, w1 = procesar_angulos(traj[:, 0], traj[:, 2])
    th2, w2 = procesar_angulos(traj[:, 1], traj[:, 3])
    
    poincare_pts = seccion_poincare(traj)

    # Configuración de figura
    plt.rcParams['font.family'] = 'serif'
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), dpi=300)
    fig.suptitle(titulo_principal, fontsize=14, fontweight='bold', y=1.02)

    lw_linea = 0.3
    s_poincare = 1.5

    # Puntos de equilibrio (Centros y Sillas)
    centros = [(0, 0)]
    sillas = [(-np.pi, 0), (np.pi, 0)]

    # --- Fase 1 ---
    axes[0].plot(th1, w1, '-', lw=lw_linea, color='steelblue')
    axes[0].set_title('Plano de Fase 1')
    axes[0].set_xlabel(r'$\theta_1$')
    axes[0].set_ylabel(r'$\omega_1$')

    # --- Fase 2 ---
    axes[1].plot(th2, w2, '-', lw=lw_linea, color='darkorange')
    axes[1].set_title('Plano de Fase 2')
    axes[1].set_xlabel(r'$\theta_2$')
    axes[1].set_ylabel(r'$\omega_2$')

    # Añadir equilibrios en ambos planos de fase
    for ax in axes[:2]:
        for c in centros: ax.plot(c[0], c[1], 'go', markersize=6, markeredgecolor='black', zorder=5)
        for s in sillas:  ax.plot(s[0], s[1], 'rx', markersize=6, markeredgewidth=1.5, zorder=5)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-np.pi - 0.2, np.pi + 0.2)

    # --- Poincaré ---
    axes[2].set_title('Sección de Poincaré')
    axes[2].set_xlabel(r'$\theta_2$')
    axes[2].set_ylabel(r'$\omega_2$')
    axes[2].grid(True, alpha=0.3)
    
    if len(poincare_pts) > 0:
        axes[2].scatter(poincare_pts[:, 0], poincare_pts[:, 1], s=s_poincare, c='purple', alpha=0.8, edgecolors='none')
    
    # Centro proyectado en Poincaré
    axes[2].plot(0, 0, 'go', markersize=6, markeredgecolor='black', zorder=5)
    axes[2].set_xlim(-np.pi - 0.2, np.pi + 0.2)

    plt.tight_layout()
    plt.savefig(archivo, format='png', bbox_inches='tight', dpi=300)
    plt.close()
    print(f"-> Guardado: {archivo}\n")


if __name__ == '__main__':
    # Usando tus parámetros
    p = {'m1': 1.0, 'm2': 1.0, 'l1': 1.0, 'l2': 1.0, 'g': 9.81}
    dt = 0.005

    # 1. Régimen Regular
    # t_max a 150s es suficiente para mostrar toros regulares sin sobreescribir demasiado
    generar_graficas(
        estado_inicial=np.array([0.20, 0.25, 0.0, 0.0]), 
        t_max=150.0, dt=dt, params=p, 
        archivo='regimen_regular.png', 
        titulo_principal='Régimen Regular'
    )

    # 2. Régimen Caótico
    # t_max alto para que el polvo fractal en Poincaré sea evidente
    generar_graficas(
        estado_inicial=np.array([np.pi - 0.05, np.pi - 0.05, 0.0, 0.0]), 
        t_max=600.0, dt=dt, params=p, 
        archivo='regimen_caotico.png', 
        titulo_principal='Régimen Caótico'
    )