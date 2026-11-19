"""Cálculo del Exponente de Lyapunov de Tiempo Finito (FTLE) y Mapas de Caos.

Implementa el método canónico de Benettin de dos trayectorias renormalizadas:
- Perturbación inicial d0 = 1e-8 en theta1
- Renormalización de la separación cada 10 pasos
- lambda_T = (1 / T) * sum(log(d_k / d0))
"""

from typing import Union, Tuple, Optional
import numpy as np
from .config import ParametrosPendulo, PARAMS_DEFECTO
from .integrate import paso_rk4_38


def calcular_ftle(
    x0: np.ndarray,
    t_max: float = 40.0,
    dt: float = 0.01,
    params: ParametrosPendulo = PARAMS_DEFECTO,
    d0: float = 1e-8,
    intervalo_renorm: int = 10,
) -> Union[float, np.ndarray]:
    """Calcula el FTLE (exponente de Lyapunov de tiempo finito) por el método de Benettin.
    
    Acepta x0 como vector de 1 condición inicial (4,) o un arreglo de condiciones iniciales (N, 4).
    
    Parámetros:
        x0: Estado o lote de estados iniciales (theta1, theta2, omega1, omega2).
        t_max: Tiempo total de integración (por defecto 40.0 s).
        dt: Paso temporal fijo para RK4 3/8 (por defecto 0.01 s).
        params: Parámetros del péndulo doble.
        d0: Magnitud de la perturbación inicial en theta1 (1e-8).
        intervalo_renorm: Cada cuántos pasos de integración se renormaliza (10 pasos).
        
    Retorna:
        lambda_T: Valor escalar o arreglo de forma (N,) con el FTLE en s^-1.
    """
    x0_arr = np.asarray(x0, dtype=np.float64)
    es_escalar = (x0_arr.ndim == 1)

    if es_escalar:
        x_base = x0_arr.copy()
        x_pert = x0_arr.copy()
        x_pert[0] += d0
        acum_log = 0.0
    else:
        x_base = x0_arr.copy()
        x_pert = x0_arr.copy()
        x_pert[:, 0] += d0
        acum_log = np.zeros(x0_arr.shape[0], dtype=np.float64)

    n_pasos = int(np.round(t_max / dt))
    tiempo_total = n_pasos * dt

    for paso in range(1, n_pasos + 1):
        x_base = paso_rk4_38(x_base, dt, params)
        x_pert = paso_rk4_38(x_pert, dt, params)

        if paso % intervalo_renorm == 0 or paso == n_pasos:
            delta = x_pert - x_base
            dist = np.linalg.norm(delta, axis=-1)

            # Evitar división por cero si dist fuera extremadamente pequeña
            dist_segura = np.maximum(dist, 1e-15)
            acum_log += np.log(dist_segura / d0)

            # Renormalizar la trayectoria perturbada a distancia d0 en la dirección de delta
            if es_escalar:
                x_pert = x_base + (d0 / dist_segura) * delta
            else:
                factor = (d0 / dist_segura)[..., np.newaxis]
                x_pert = x_base + factor * delta

    lambda_t = acum_log / tiempo_total
    return float(lambda_t) if es_escalar else lambda_t


def generar_malla_mapa_caos(
    n_puntos: int = 200,
    rango: Tuple[float, float] = (-np.pi, np.pi),
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Genera la malla 2D regular de condiciones iniciales en el plano (theta1(0), theta2(0)) desde el reposo.
    
    Retorna: (theta1_grid, theta2_grid, x0_batch)
        theta1_grid: Matriz 2D de tamaño (n_puntos, n_puntos)
        theta2_grid: Matriz 2D de tamaño (n_puntos, n_puntos)
        x0_batch: Arreglo de forma (n_puntos * n_puntos, 4) con velocidades cero.
    """
    th1 = np.linspace(rango[0], rango[1], n_puntos)
    th2 = np.linspace(rango[0], rango[1], n_puntos)
    th1_grid, th2_grid = np.meshgrid(th1, th2)

    n_total = n_puntos * n_puntos
    x0_batch = np.zeros((n_total, 4), dtype=np.float64)
    x0_batch[:, 0] = th1_grid.ravel()
    x0_batch[:, 1] = th2_grid.ravel()
    # omega1 y omega2 quedan en 0

    return th1_grid, th2_grid, x0_batch


def _ftle_worker(args):
    """Función auxiliar para paralelización de lotes de FTLE."""
    lote, t_max, dt, params = args
    return calcular_ftle(lote, t_max=t_max, dt=dt, params=params)


def calcular_mapa_ftle_malla(
    n_puntos: int = 200,
    t_max: float = 40.0,
    dt: float = 0.01,
    params: ParametrosPendulo = PARAMS_DEFECTO,
    tamano_lote: int = 2500,
    n_jobs: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calcula el mapa de FTLE sobre una malla (n_puntos x n_puntos) procesando por lotes.
    
    Aprovecha paralelismo multinúcleo en Linux para calcular la rejilla en pocos segundos.
    
    Retorna: (th1_grid, th2_grid, ftle_grid)
    """
    import os
    import multiprocessing as mp

    th1_grid, th2_grid, x0_batch = generar_malla_mapa_caos(n_puntos=n_puntos)
    n_total = x0_batch.shape[0]

    lotes = [
        (x0_batch[i:min(i + tamano_lote, n_total)], t_max, dt, params)
        for i in range(0, n_total, tamano_lote)
    ]

    if n_jobs is None or n_jobs <= 0:
        n_jobs = min(8, os.cpu_count() or 1)

    if n_jobs > 1 and hasattr(mp, "get_context"):
        try:
            ctx = mp.get_context("fork")
            with ctx.Pool(processes=n_jobs) as pool:
                resultados_lotes = pool.map(_ftle_worker, lotes)
            ftle_vec = np.concatenate(resultados_lotes)
        except Exception:
            # Fallback seguro secuencial
            resultados_lotes = [_ftle_worker(lote_args) for lote_args in lotes]
            ftle_vec = np.concatenate(resultados_lotes)
    else:
        resultados_lotes = [_ftle_worker(lote_args) for lote_args in lotes]
        ftle_vec = np.concatenate(resultados_lotes)

    ftle_grid = ftle_vec.reshape(n_puntos, n_puntos)
    return th1_grid, th2_grid, ftle_grid

