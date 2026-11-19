"""Detección y conteo riguroso de vueltas (pasos por arriba) del péndulo doble.

Definición matemática oficial:
Llevando el ángulo theta sin envolver:
k(t) = floor((theta + pi) / (2*pi))
Cada cambio de k entre pasos sucesivos representa un cruce por la vertical superior:
- Delta_k = k(t_{n+1}) - k(t_n)
- Conteo sin signo: sum(|Delta_k|)
- Conteo neto con signo: sum(Delta_k)
"""

from typing import Tuple, Union
import numpy as np


def calcular_indice_giro(theta: Union[float, np.ndarray]) -> Union[int, np.ndarray]:
    """Calcula el índice de sector k = floor((theta + pi) / (2*pi))."""
    return np.floor((theta + np.pi) / (2.0 * np.pi)).astype(int)


def contar_vueltas_trayectoria(
    theta_serie: np.ndarray,
    axis: int = 0,
) -> Tuple[Union[int, np.ndarray], Union[int, np.ndarray]]:
    """Calcula los pasos por arriba a lo largo de una trayectoria temporal de ángulos.
    
    Parámetros:
        theta_serie: Arreglo de ángulos continuos (sin envolver). El eje temporal
                     se especifica con 'axis' (por defecto 0).
        axis: Eje correspondiente al tiempo.
        
    Retorna:
        (conteo_sin_signo, conteo_con_signo):
            - conteo_sin_signo: Total absoluto de pasos por arriba.
            - conteo_con_signo: Conteo neto con signo (giros netos).
    """
    k = calcular_indice_giro(theta_serie)
    dk = np.diff(k, axis=axis)

    conteo_sin_signo = np.sum(np.abs(dk), axis=axis)
    conteo_con_signo = np.sum(dk, axis=axis)

    if np.ndim(conteo_sin_signo) == 0:
        return int(conteo_sin_signo), int(conteo_con_signo)
    return conteo_sin_signo, conteo_con_signo


def pasos_por_arriba_estado(
    x_traj: np.ndarray,
    masa: int = 2,
) -> Tuple[Union[int, np.ndarray], Union[int, np.ndarray]]:
    """Calcula las vueltas para la masa indicada (1 o 2) en una trayectoria de estados.
    
    x_traj: Arreglo con estados x = (theta1, theta2, omega1, omega2).
            La dimensión de estado (4) es la última.
    """
    col = 0 if masa == 1 else 1
    theta = x_traj[..., col]
    return contar_vueltas_trayectoria(theta, axis=0)
