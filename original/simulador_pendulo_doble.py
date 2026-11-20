"""
Simulador del Péndulo Doble
===========================
Basado en el informe "Análisis Dinámico Riguroso y Simulación
Computacional del Péndulo Doble" (Rodrigo Moreno Cruz, ITAM).

Implementa el campo vectorial f(x) derivado por el formalismo lagrangiano
en la sección 1, integrado vía Runge-Kutta de orden 4 (clásico o
modificado en su variante 3/8 de Kutta).

Estado del sistema:  x = [theta1, theta2, omega1, omega2]  ∈  T² × R²

Disposición de la ventana:
    Parte superior:
        - Izquierda: panel de parámetros (sliders, botones).
        - Derecha:   animación del péndulo doble.
    Parte inferior:
        - Retrato de fase (theta1, omega1) con equilibrios.
        - Retrato de fase (theta2, omega2) con equilibrios.
        - Sección de Poincaré  Σ = { theta1 = 0, omega1 < 0 }.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.animation import FuncAnimation
import time


# =============================================================
#   1. CAMPO VECTORIAL Y ENERGÍA  (Sección 1.2 - 1.3 del informe)
# =============================================================

def derivadas(estado, p):
    """
    f(x) = dx/dt  para el péndulo doble.

    Parte de la formulación matricial
        M(theta) · alpha + C(theta, omega) + G(theta) = 0
    y despeja alpha = -M^{-1} (C + G).
    """
    theta1, theta2, omega1, omega2 = estado
    m1, m2, l1, l2, g = p['m1'], p['m2'], p['l1'], p['l2'], p['g']

    delta = theta1 - theta2
    cd, sd = np.cos(delta), np.sin(delta)

    # Matriz de inercia (sin dividir entre l_i, según ecuación 10 del informe)
    M00 = (m1 + m2) * l1
    M01 = m2 * l2 * cd
    M10 = m2 * l1 * cd
    M11 = m2 * l2
    detM = M00 * M11 - M01 * M10  # > 0 por la sección 2.1

    # Vectores C (centrífugo) y G (gravitatorio)
    C0 = m2 * l2 * omega2 * omega2 * sd
    C1 = -m2 * l1 * omega1 * omega1 * sd
    G0 = (m1 + m2) * g * np.sin(theta1)
    G1 = m2 * g * np.sin(theta2)

    rhs0 = -(C0 + G0)
    rhs1 = -(C1 + G1)
    alpha1 = (M11 * rhs0 - M01 * rhs1) / detM
    alpha2 = (-M10 * rhs0 + M00 * rhs1) / detM

    return np.array([omega1, omega2, alpha1, alpha2])


def energia(estado, p):
    """Integral primera E(x) = T + V (sección 3.1)."""
    theta1, theta2, omega1, omega2 = estado
    m1, m2, l1, l2, g = p['m1'], p['m2'], p['l1'], p['l2'], p['g']
    delta = theta1 - theta2
    T = (0.5 * (m1 + m2) * l1 * l1 * omega1 * omega1
         + 0.5 * m2 * l2 * l2 * omega2 * omega2
         + m2 * l1 * l2 * omega1 * omega2 * np.cos(delta))
    V = ((m1 + m2) * g * l1 * (1 - np.cos(theta1))
         + m2 * g * l2 * (1 - np.cos(theta2)))
    return T + V


def energias_criticas(p):
    """Umbrales de bifurcación E_crit_i de la sección 2.4."""
    m1, m2, l1, l2, g = p['m1'], p['m2'], p['l1'], p['l2'], p['g']
    Ec1 = 2 * m2 * g * l2                       # masa inferior arriba
    Ec2 = 2 * (m1 + m2) * g * l1                # masa superior arriba
    Ec3 = 2 * (m1 + m2) * g * l1 + 2 * m2 * g * l2
    return Ec1, Ec2, Ec3


# =============================================================
#   2. INTEGRADORES NUMÉRICOS
# =============================================================

def rk4_clasico(estado, dt, p):
    """Runge-Kutta clásico de orden 4."""
    k1 = derivadas(estado, p)
    k2 = derivadas(estado + 0.5 * dt * k1, p)
    k3 = derivadas(estado + 0.5 * dt * k2, p)
    k4 = derivadas(estado + dt * k3, p)
    return estado + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def rk4_modificado(estado, dt, p):
    """
    Runge-Kutta de orden 4 modificado (regla 3/8 de Kutta).
    Misma precisión que RK4 clásico pero con coeficientes que reducen
    el sesgo direccional en sistemas oscilatorios.
    """
    k1 = derivadas(estado, p)
    k2 = derivadas(estado + dt * (k1 / 3.0), p)
    k3 = derivadas(estado + dt * (-k1 / 3.0 + k2), p)
    k4 = derivadas(estado + dt * (k1 - k2 + k3), p)
    return estado + (dt / 8.0) * (k1 + 3.0 * k2 + 3.0 * k3 + k4)


def integrar(estado_inicial, t_max, dt, p, metodo='rk4'):
    paso = rk4_modificado if metodo == 'rk4_mod' else rk4_clasico
    n = int(t_max / dt)
    traj = np.empty((n + 1, 4))
    traj[0] = estado_inicial
    for i in range(n):
        traj[i + 1] = paso(traj[i], dt, p)
    return np.arange(n + 1) * dt, traj


# =============================================================
#   3. SECCIÓN DE POINCARÉ   Σ = { theta1 = 0,  omega1 < 0 }
# =============================================================

def seccion_poincare(traj):
    """Detecta cruces theta1 ≡ 0 (mod 2π) con omega1 < 0 e interpola."""
    th1_wrap = np.arctan2(np.sin(traj[:, 0]), np.cos(traj[:, 0]))
    pts = []
    for i in range(len(th1_wrap) - 1):
        a, b = th1_wrap[i], th1_wrap[i + 1]
        if a == 0.0:
            if traj[i, 2] < 0:
                pts.append((traj[i, 1], traj[i, 3]))
            continue
        if a * b < 0 and abs(a - b) < np.pi:    # cruce real (sin wrap)
            if traj[i, 2] < 0:                  # omega1 < 0
                t = a / (a - b)
                th2 = traj[i, 1] + t * (traj[i + 1, 1] - traj[i, 1])
                w2 = traj[i, 3] + t * (traj[i + 1, 3] - traj[i, 3])
                # Wrap theta2 a [-π, π] para visualización
                th2 = np.arctan2(np.sin(th2), np.cos(th2))
                pts.append((th2, w2))
    return np.array(pts) if pts else np.empty((0, 2))


# =============================================================
#   4. INTERFAZ GRÁFICA
# =============================================================

class SimuladorPenduloDoble:
    DEF = dict(m1=1.0, m2=1.0, l1=1.0, l2=1.0, g=9.81,
               theta1_0=np.pi / 2, theta2_0=np.pi / 2,
               omega1_0=0.0, omega2_0=0.0,
               t_max=30.0, dt=0.005)

    def __init__(self):
        self.params = {k: self.DEF[k] for k in ('m1', 'm2', 'l1', 'l2', 'g')}
        self.estado_inicial = np.array(
            [self.DEF['theta1_0'], self.DEF['theta2_0'],
             self.DEF['omega1_0'], self.DEF['omega2_0']])
        self.t_max = self.DEF['t_max']
        self.dt = self.DEF['dt']
        self.metodo = 'rk4'

        self.fig = plt.figure(figsize=(16, 9))
        try:
            self.fig.canvas.manager.set_window_title(
                'Simulador del Péndulo Doble')
        except Exception:
            pass

        # Layout: panel izquierdo más delgado, área de info a la derecha
        gs_top = gridspec.GridSpec(1, 2, figure=self.fig,
                                   left=0.22, right=0.98,
                                   top=0.95, bottom=0.55,
                                   wspace=0.15,
                                   width_ratios=[3.5, 1.2])
        gs_bot = gridspec.GridSpec(1, 3, figure=self.fig,
                                   left=0.22, right=0.98,
                                   top=0.45, bottom=0.07,
                                   wspace=0.30)

        self.ax_pend = self.fig.add_subplot(gs_top[0, 0])
        self.ax_info = self.fig.add_subplot(gs_top[0, 1])
        self.ax_info.axis('off')
        
        self.ax_f1 = self.fig.add_subplot(gs_bot[0, 0])
        self.ax_f2 = self.fig.add_subplot(gs_bot[0, 1])
        self.ax_p = self.fig.add_subplot(gs_bot[0, 2])

        self._init_pendulo()
        self._init_fases()
        self._init_poincare()
        self._init_panel()

        self._recalcular()
        self._iniciar_animacion()

    # ---- subplots --------------------------------------------------

    def _init_pendulo(self):
        ax = self.ax_pend
        ax.set_aspect('equal')
        ax.set_title('Movimiento del Péndulo Doble', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.scatter([0], [0], c='black', zorder=5, s=30)
        self.line_pend, = ax.plot([], [], 'o-', lw=2.2,
                                  color='steelblue', markersize=9)
        self.trace, = ax.plot([], [], '-', lw=0.8, alpha=0.55,
                              color='crimson')
        self.text_info = self.ax_info.text(
            0.0, 0.95, '', transform=self.ax_info.transAxes,
            va='top', ha='left', fontsize=9, family='monospace',
            bbox=dict(facecolor='white', alpha=0.85, edgecolor='gray'))

    def _init_fases(self):
        for ax, sub in [(self.ax_f1, '1'), (self.ax_f2, '2')]:
            ax.set_xlabel(rf'$\theta_{sub}$ [rad]')
            ax.set_ylabel(rf'$\omega_{sub}$ [rad/s]')
            ax.set_title(rf'Retrato de Fase ($\theta_{sub}, \omega_{sub}$)')
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color='gray', lw=0.5)
            ax.axvline(0, color='gray', lw=0.5)
            # Equilibrios proyectados sobre cada plano: (0,0) centro,
            # (±π, 0) puntos silla, según sección 2.3.
            ax.plot([0], [0], 'o', color='forestgreen',
                    markersize=9, label='Centro estable')
            ax.plot([-np.pi, np.pi], [0, 0], 'x', color='red',
                    markersize=11, mew=2, label='Silla')
            ax.legend(loc='upper right', fontsize=8, framealpha=0.85)
        self.line_f1, = self.ax_f1.plot([], [], '-', lw=0.7,
                                        color='steelblue')
        self.line_f2, = self.ax_f2.plot([], [], '-', lw=0.7,
                                        color='darkorange')

    def _init_poincare(self):
        ax = self.ax_p
        ax.set_xlabel(r'$\theta_2$ [rad]')
        ax.set_ylabel(r'$\omega_2$ [rad/s]')
        ax.set_title(r'Sección de Poincaré ($\theta_1=0,\;\omega_1<0$)')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='gray', lw=0.5)
        ax.axvline(0, color='gray', lw=0.5)
        self.scatter_p = ax.scatter([], [], s=4, c='purple', alpha=0.65)

    # ---- panel de control -----------------------------------------

    def _init_panel(self):
        x0, w = 0.03, 0.14
        h = 0.015
        y0 = 0.95
        dy = 0.032

        self.fig.text(x0, y0 + 0.015, 'Parámetros del Sistema',
                      fontsize=11, fontweight='bold')

        defs = [
            ('m1', r'$m_1$ [kg]',     0.1, 5.0,  self.DEF['m1'],   '%.2f'),
            ('m2', r'$m_2$ [kg]',     0.1, 5.0,  self.DEF['m2'],   '%.2f'),
            ('l1', r'$l_1$ [m]',      0.1, 2.0,  self.DEF['l1'],   '%.2f'),
            ('l2', r'$l_2$ [m]',      0.1, 2.0,  self.DEF['l2'],   '%.2f'),
            ('g',  r'$g$ [m/s$^2$]',  1.0, 25.0, self.DEF['g'],    '%.2f'),
            ('theta1_0', r'$\theta_1^0$ [rad]',
             -np.pi, np.pi, self.DEF['theta1_0'], '%.3f'),
            ('theta2_0', r'$\theta_2^0$ [rad]',
             -np.pi, np.pi, self.DEF['theta2_0'], '%.3f'),
            ('omega1_0', r'$\omega_1^0$ [rad/s]',
             -10.0, 10.0, self.DEF['omega1_0'], '%.2f'),
            ('omega2_0', r'$\omega_2^0$ [rad/s]',
             -10.0, 10.0, self.DEF['omega2_0'], '%.2f'),
            ('t_max', r'$t_{max}$ [s]',
             5.0, 200.0, self.DEF['t_max'], '%.1f'),
            ('dt',    r'$dt$ [s]',
             0.001, 0.02, self.DEF['dt'], '%.4f'),
        ]

        self.sliders = {}
        for i, (k, lbl, vmin, vmax, vinit, fmt) in enumerate(defs):
            ax_s = self.fig.add_axes([x0, y0 - i * dy, w, h])
            self.sliders[k] = Slider(ax_s, lbl, vmin, vmax,
                                     valinit=vinit, valfmt=fmt)

        # Botones
        y_btn = y0 - len(defs) * dy - 0.005
        ax_b1 = self.fig.add_axes([x0, y_btn, 0.065, 0.035])
        ax_b2 = self.fig.add_axes([x0 + 0.075, y_btn, 0.065, 0.035])
        self.btn_apply = Button(ax_b1, 'Recalcular',
                                color='lightsteelblue',
                                hovercolor='steelblue')
        self.btn_reset = Button(ax_b2, 'Reset',
                                color='lavenderblush',
                                hovercolor='lightcoral')
        self.btn_apply.on_clicked(self._on_apply)
        self.btn_reset.on_clicked(self._on_reset)

        # Selector de método de integración
        y_rad = y_btn - 0.12
        ax_r = self.fig.add_axes([x0, y_rad, 0.14, 0.09])
        ax_r.set_title('Integrador', fontsize=9, loc='left')
        self.radio = RadioButtons(ax_r, ('RK4 clásico',
                                         'RK4 modificado (3/8)'),
                                  active=0)
        self.radio.on_clicked(self._on_metodo)

        # Selector de "preset" para condiciones iniciales conocidas
        y_pre = y_rad - 0.17
        ax_pre = self.fig.add_axes([x0, y_pre, 0.14, 0.13])
        ax_pre.set_title('Régimen sugerido', fontsize=9, loc='left')
        self.radio_preset = RadioButtons(
            ax_pre, ('Personalizado',
                     'Pequeñas oscilaciones',
                     'Cuasi-periódico',
                     'Caótico'), active=0)
        self.radio_preset.on_clicked(self._on_preset)

    # ---- callbacks ------------------------------------------------

    def _on_metodo(self, label):
        self.metodo = 'rk4_mod' if 'modificado' in label else 'rk4'

    def _on_apply(self, _):
        self._lee_sliders()
        self._recalcular()

    def _on_reset(self, _):
        for s in self.sliders.values():
            s.reset()
        self.radio_preset.set_active(0)
        self._on_apply(None)

    def _on_preset(self, label):
        if label == 'Personalizado':
            return
        if label == 'Pequeñas oscilaciones':
            ic = (0.20, 0.25, 0.0, 0.0)
        elif label == 'Cuasi-periódico':
            ic = (0.80, 0.60, 0.0, 0.0)
        else:  # Caótico
            ic = (np.pi - 0.05, np.pi - 0.05, 0.0, 0.0)
        self.sliders['theta1_0'].set_val(ic[0])
        self.sliders['theta2_0'].set_val(ic[1])
        self.sliders['omega1_0'].set_val(ic[2])
        self.sliders['omega2_0'].set_val(ic[3])
        self._on_apply(None)

    def _lee_sliders(self):
        for k in ('m1', 'm2', 'l1', 'l2', 'g'):
            self.params[k] = self.sliders[k].val
        self.estado_inicial = np.array(
            [self.sliders['theta1_0'].val, self.sliders['theta2_0'].val,
             self.sliders['omega1_0'].val, self.sliders['omega2_0'].val])
        self.t_max = self.sliders['t_max'].val
        self.dt = max(self.sliders['dt'].val, 1e-4)

    # ---- núcleo de simulación -------------------------------------

    def _recalcular(self):
        self.estado = self.estado_inicial.copy()
        self.t = 0.0
        self.t0_real = time.time()
        
        self.E0 = energia(self.estado, self.params)
        self.Ec1, self.Ec2, self.Ec3 = energias_criticas(self.params)

        self.hist_th1 = []
        self.hist_th2 = []
        self.hist_w1 = []
        self.hist_w2 = []
        self.hist_x2 = []
        self.hist_y2 = []
        self.pts_poincare = []

        l1, l2 = self.params['l1'], self.params['l2']
        L = l1 + l2
        self.ax_pend.set_xlim(-1.20 * L, 1.20 * L)
        self.ax_pend.set_ylim(-1.20 * L, 1.20 * L)

        self.line_f1.set_data([], [])
        self.line_f2.set_data([], [])
        self.scatter_p.set_offsets(np.empty((0, 2)))
        self.ax_p.set_title(r'Sección de Poincaré ($\theta_1=0,\;\omega_1<0$)')

        for ax in (self.ax_f1, self.ax_f2):
            ax.set_xlim(-np.pi - 0.2, np.pi + 0.2)
            ax.set_ylim(-10, 10)

        self.ax_p.set_xlim(-np.pi - 0.2, np.pi + 0.2)
        self.ax_p.set_ylim(-10, 10)

        self.fig.canvas.draw_idle()

    # ---- animación -------------------------------------------------

    def _iniciar_animacion(self):
        def update(_):
            target_t = time.time() - self.t0_real
            paso = rk4_modificado if self.metodo == 'rk4_mod' else rk4_clasico
            max_steps = 200
            steps = 0
            
            max_hist = int(self.t_max / self.dt)
            
            while self.t < target_t and steps < max_steps:
                prev_th1 = self.estado[0]
                prev_w1 = self.estado[2]
                prev_th2 = self.estado[1]
                prev_w2 = self.estado[3]
                
                self.estado = paso(self.estado, self.dt, self.params)
                self.t += self.dt
                steps += 1
                
                th1w = np.arctan2(np.sin(self.estado[0]), np.cos(self.estado[0]))
                th2w = np.arctan2(np.sin(self.estado[1]), np.cos(self.estado[1]))
                w1 = self.estado[2]
                w2 = self.estado[3]
                
                self.hist_th1.append(th1w)
                self.hist_th2.append(th2w)
                self.hist_w1.append(w1)
                self.hist_w2.append(w2)
                
                l1, l2 = self.params['l1'], self.params['l2']
                x1 = l1 * np.sin(self.estado[0])
                y1 = -l1 * np.cos(self.estado[0])
                x2 = x1 + l2 * np.sin(self.estado[1])
                y2 = y1 - l2 * np.cos(self.estado[1])
                self.hist_x2.append(x2)
                self.hist_y2.append(y2)
                
                a = np.arctan2(np.sin(prev_th1), np.cos(prev_th1))
                b = th1w
                if a == 0.0 and prev_w1 < 0:
                    self.pts_poincare.append((th2w, w2))
                elif a * b < 0 and abs(a - b) < np.pi:
                    if prev_w1 < 0:
                        tt = a / (a - b)
                        interp_th2 = prev_th2 + tt * (self.estado[1] - prev_th2)
                        interp_w2 = prev_w2 + tt * (w2 - prev_w2)
                        interp_th2w = np.arctan2(np.sin(interp_th2), np.cos(interp_th2))
                        self.pts_poincare.append((interp_th2w, interp_w2))
                
                if len(self.hist_th1) > max_hist:
                    self.hist_th1.pop(0)
                    self.hist_th2.pop(0)
                    self.hist_w1.pop(0)
                    self.hist_w2.pop(0)
                    self.hist_x2.pop(0)
                    self.hist_y2.pop(0)

            if not self.hist_th1:
                return ()

            x1 = l1 * np.sin(self.estado[0])
            y1 = -l1 * np.cos(self.estado[0])
            x2 = self.hist_x2[-1]
            y2 = self.hist_y2[-1]
            
            self.line_pend.set_data([0, x1, x2], [0, y1, y2])
            
            tail = max(0, len(self.hist_x2) - 600)
            self.trace.set_data(self.hist_x2[tail:], self.hist_y2[tail:])

            self.line_f1.set_data(self.hist_th1, self.hist_w1)
            self.line_f2.set_data(self.hist_th2, self.hist_w2)

            for ax, w in [(self.ax_f1, self.estado[2]), (self.ax_f2, self.estado[3])]:
                ymin, ymax = ax.get_ylim()
                if w < ymin + 1 or w > ymax - 1:
                    ax.set_ylim(min(ymin, w - 3), max(ymax, w + 3))

            if self.pts_poincare:
                pts = np.array(self.pts_poincare)
                self.scatter_p.set_offsets(pts)
                
                wmin, wmax = float(pts[:, 1].min()), float(pts[:, 1].max())
                ymin, ymax = self.ax_p.get_ylim()
                if wmin < ymin + 1 or wmax > ymax - 1:
                    margin = max(5.0, 0.15 * (wmax - wmin))
                    self.ax_p.set_ylim(wmin - margin, wmax + margin)
                
                self.ax_p.set_title(
                    rf'Sección de Poincaré ($\theta_1=0,\;\omega_1<0$)'
                    f' — {len(pts)} cruces')

            E = energia(self.estado, self.params)
            dE = E - self.E0

            if E < self.Ec1:
                regimen = 'Toros KAM (regular)'
            elif E < self.Ec2:
                regimen = 'rotaciones masa 2'
            elif E < self.Ec3:
                regimen = 'rotaciones masa 1'
            else:
                regimen = 'caos pleno'

            self.text_info.set_text(
                f't       = {self.t:7.2f} s\n'
                f'E       = {E:9.4f} J\n'
                f'ΔE      = {dE:+8.2e}\n'
                f'E_c1    = {self.Ec1:9.4f}\n'
                f'E_c2    = {self.Ec2:9.4f}\n'
                f'E_c3    = {self.Ec3:9.4f}\n'
                f'régimen = {regimen}\n'
                f'método  = {self.metodo}'
            )

            return self.line_pend, self.trace, self.text_info, self.line_f1, self.line_f2, self.scatter_p

        self.anim = FuncAnimation(self.fig, update, interval=20,
                                  blit=False, cache_frame_data=False)


def main():
    SimuladorPenduloDoble()
    plt.show()


if __name__ == '__main__':
    main()
