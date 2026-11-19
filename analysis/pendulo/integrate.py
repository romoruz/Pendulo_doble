"""Integradores numéricos para el péndulo doble.

Incluye:
- RK4 variante 3/8 de Kutta (escalar y vectorizado en forma (N, 4))
- Punto medio implícito en coordenadas canónicas (theta, p) (simpléctico, orden 2)
"""

from typing import Tuple, Optional
import numpy as np
from .config import ParametrosPendulo, PARAMS_DEFECTO
from .model import (
    campo_derivadas,
    campo_hamiltoniano,
    x_a_canonico,
    canonico_a_x,
)


def paso_rk4_38(
    x: np.ndarray,
    dt: float,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> np.ndarray:
    """Avanza un paso de tiempo dt usando RK4 variante 3/8 de Kutta.
    
    Acepta x con forma (4,) o (N, 4).
    k1 = f(x)
    k2 = f(x + dt/3 * k1)
    k3 = f(x - dt/3 * k1 + dt * k2)
    k4 = f(x + dt * k1 - dt * k2 + dt * k3)
    x_next = x + dt/8 * (k1 + 3*k2 + 3*k3 + k4)
    """
    k1 = campo_derivadas(x, params)
    k2 = campo_derivadas(x + (dt / 3.0) * k1, params)
    k3 = campo_derivadas(x - (dt / 3.0) * k1 + dt * k2, params)
    k4 = campo_derivadas(x + dt * k1 - dt * k2 + dt * k3, params)

    x_siguiente = x + (dt / 8.0) * (k1 + 3.0 * k2 + 3.0 * k3 + k4)
    return x_siguiente


def paso_rk4_38_hamiltoniano(
    z: np.ndarray,
    dt: float,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> np.ndarray:
    """Avanza un paso dt del sistema hamiltoniano en (theta, p) con RK4 3/8."""
    k1 = campo_hamiltoniano(z, params)
    k2 = campo_hamiltoniano(z + (dt / 3.0) * k1, params)
    k3 = campo_hamiltoniano(z - (dt / 3.0) * k1 + dt * k2, params)
    k4 = campo_hamiltoniano(z + dt * k1 - dt * k2 + dt * k3, params)

    z_siguiente = z + (dt / 8.0) * (k1 + 3.0 * k2 + 3.0 * k3 + k4)
    return z_siguiente


def integrar_rk4(
    x0: np.ndarray,
    t_max: float,
    dt: float,
    params: ParametrosPendulo = PARAMS_DEFECTO,
    t0: float = 0.0,
    guardar_todos: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Integra el sistema con RK4 3/8 desde t0 hasta t_max.
    
    Parámetros:
        x0: Estado inicial, forma (4,) o (N, 4).
        t_max: Tiempo final.
        dt: Paso de tiempo fijo.
        params: Parámetros del péndulo.
        guardar_todos: Si es False, solo retorna los tiempos y el último estado.
        
    Retorna:
        t: Arreglo 1D de tiempos.
        trayectoria: Arreglo de estados con forma (n_pasos + 1, ...) si guardar_todos=True,
                     o (1, ...) con el estado final si guardar_todos=False.
    """
    n_pasos = int(np.round((t_max - t0) / dt))
    t = t0 + np.arange(n_pasos + 1) * dt

    x = np.array(x0, dtype=np.float64, copy=True)

    if not guardar_todos:
        for _ in range(n_pasos):
            x = paso_rk4_38(x, dt, params)
        return t, x

    # Dimensiones para preasignar
    if x.ndim == 1:
        trayectoria = np.empty((n_pasos + 1, 4), dtype=np.float64)
    else:
        trayectoria = np.empty((n_pasos + 1, x.shape[0], 4), dtype=np.float64)

    trayectoria[0] = x
    for i in range(1, n_pasos + 1):
        x = paso_rk4_38(x, dt, params)
        trayectoria[i] = x

    return t, trayectoria


def integrar_rk4_hamiltoniano(
    x0: np.ndarray,
    t_max: float,
    dt: float,
    params: ParametrosPendulo = PARAMS_DEFECTO,
    t0: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Integra en variables canónicas (theta, p) usando RK4 3/8 y convierte a (theta, omega)."""
    n_pasos = int(np.round((t_max - t0) / dt))
    t = t0 + np.arange(n_pasos + 1) * dt

    z = x_a_canonico(np.array(x0, dtype=np.float64, copy=True), params)
    trayectoria = np.empty((n_pasos + 1, 4), dtype=np.float64)
    trayectoria[0] = canonico_a_x(z, params)

    for i in range(1, n_pasos + 1):
        z = paso_rk4_38_hamiltoniano(z, dt, params)
        trayectoria[i] = canonico_a_x(z, params)

    return t, trayectoria


# ---------------------------------------------------------------------------
# Punto Medio Implícito en Coordenadas Canónicas
# ---------------------------------------------------------------------------

def paso_punto_medio_implicito(
    z_n: np.ndarray,
    dt: float,
    params: ParametrosPendulo = PARAMS_DEFECTO,
    tol: float = 1e-13,
    max_iter: int = 100,
) -> np.ndarray:
    """Avanza un paso dt con el método del punto medio implícito canónico.
    
    z_(n+1) = z_n + dt * J_H((z_n + z_(n+1)) / 2)
    Resuelve por iteración de punto fijo:
    z^(k+1) = z_n + dt * J_H((z_n + z^(k)) / 2)
    """
    # Estimación inicial: paso de Euler canónico explícito
    z_sig = z_n + dt * campo_hamiltoniano(z_n, params)

    for iteracion in range(max_iter):
        z_medio = 0.5 * (z_n + z_sig)
        f_medio = campo_hamiltoniano(z_medio, params)
        z_nuevo = z_n + dt * f_medio

        error = np.max(np.abs(z_nuevo - z_sig))
        z_sig = z_nuevo

        if error < tol:
            return z_sig

    raise RuntimeError(
        f"Punto medio implícito no convergió tras {max_iter} iteraciones (error = {error:.2e} > tol = {tol:.2e})."
    )


def integrar_punto_medio(
    x0: np.ndarray,
    t_max: float,
    dt: float,
    params: ParametrosPendulo = PARAMS_DEFECTO,
    t0: float = 0.0,
    tol: float = 1e-13,
    max_iter: int = 100,
    guardar_cada: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Integra el péndulo doble con el integrador simpléctico del punto medio implícito.
    
    Parámetros:
        x0: Estado inicial en (theta, omega), forma (4,).
        t_max: Tiempo total.
        dt: Paso temporal.
        guardar_cada: Submuestreo para guardar puntos (útil en simulaciones largas de T=2000s).
        
    Retorna:
        t: Arreglo de tiempos guardados.
        trayectoria: Arreglo de estados (theta, omega) de forma (n_guardados, 4).
    """
    n_pasos = int(np.round((t_max - t0) / dt))
    z = x_a_canonico(np.array(x0, dtype=np.float64, copy=True), params)

    n_guardados = n_pasos // guardar_cada + 1
    t_out = np.empty(n_guardados, dtype=np.float64)
    trayectoria_out = np.empty((n_guardados, 4), dtype=np.float64)

    t_actual = t0
    idx_guardado = 0
    t_out[idx_guardado] = t_actual
    trayectoria_out[idx_guardado] = canonico_a_x(z, params)

    for i in range(1, n_pasos + 1):
        z = paso_punto_medio_implicito(z, dt, params, tol=tol, max_iter=max_iter)
        t_actual = t0 + i * dt
        if i % guardar_cada == 0:
            idx_guardado += 1
            t_out[idx_guardado] = t_actual
            trayectoria_out[idx_guardado] = canonico_a_x(z, params)

    return t_out, trayectoria_out
