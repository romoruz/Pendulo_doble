"""Modelo físico matemático del péndulo doble conservativo.

Fuente única de verdad para ecuaciones de movimiento (Lagrangiana y Hamiltoniana),
energías, matrices de inercia, linealización, umbrales y equilibrios.
"""

from typing import Dict, Tuple, Union
import numpy as np
from .config import ParametrosPendulo, PARAMS_DEFECTO


def posiciones_cartesianas(
    theta1: Union[float, np.ndarray],
    theta2: Union[float, np.ndarray],
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray],
           Union[float, np.ndarray], Union[float, np.ndarray]]:
    """Calcula las posiciones cartesianas (x1, y1) y (x2, y2).
    
    Ángulos absolutos desde la vertical hacia abajo, antihorario positivo:
    x1 = l1 * sin(theta1),  y1 = -l1 * cos(theta1)
    x2 = x1 + l2 * sin(theta2),  y2 = y1 - l2 * cos(theta2)
    """
    x1 = params.l1 * np.sin(theta1)
    y1 = -params.l1 * np.cos(theta1)
    x2 = x1 + params.l2 * np.sin(theta2)
    y2 = y1 - params.l2 * np.cos(theta2)
    return x1, y1, x2, y2


def matriz_inercia(
    theta1: Union[float, np.ndarray],
    theta2: Union[float, np.ndarray],
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray],
           Union[float, np.ndarray], Union[float, np.ndarray],
           Union[float, np.ndarray]]:
    """Elementos de M(theta) y su determinante.
    
    M11 = (m1 + m2) * l1^2
    M12 = M21 = m2 * l1 * l2 * cos(theta1 - theta2)
    M22 = m2 * l2^2
    det(M) = m2 * l1^2 * l2^2 * (m1 + m2 * sin^2(theta1 - theta2)) > 0
    
    Retorna: (M11, M12, M21, M22, det_M)
    """
    delta = theta1 - theta2
    cos_delta = np.cos(delta)
    sin_delta = np.sin(delta)

    m1, m2 = params.m1, params.m2
    l1, l2 = params.l1, params.l2

    m11 = (m1 + m2) * (l1 ** 2)
    m12 = m2 * l1 * l2 * cos_delta
    m21 = m12
    m22 = m2 * (l2 ** 2)

    det_m = m2 * (l1 ** 2) * (l2 ** 2) * (m1 + m2 * (sin_delta ** 2))
    return m11, m12, m21, m22, det_m


def energia_cinetica(
    x: np.ndarray,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> Union[float, np.ndarray]:
    """Energía cinética T = 1/2(m1+m2)l1^2 omega1^2 + 1/2 m2 l2^2 omega2^2 + m2 l1 l2 omega1 omega2 cos(Delta)."""
    theta1 = x[..., 0]
    theta2 = x[..., 1]
    omega1 = x[..., 2]
    omega2 = x[..., 3]
    delta = theta1 - theta2

    m1, m2 = params.m1, params.m2
    l1, l2 = params.l1, params.l2

    t = (
        0.5 * (m1 + m2) * (l1 ** 2) * (omega1 ** 2)
        + 0.5 * m2 * (l2 ** 2) * (omega2 ** 2)
        + m2 * l1 * l2 * omega1 * omega2 * np.cos(delta)
    )
    return t


def energia_potencial(
    x: np.ndarray,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> Union[float, np.ndarray]:
    """Energía potencial con referencia V(0,0) = 0:
    V = (m1 + m2) * g * l1 * (1 - cos(theta1)) + m2 * g * l2 * (1 - cos(theta2)).
    """
    theta1 = x[..., 0]
    theta2 = x[..., 1]
    m1, m2 = params.m1, params.m2
    l1, l2 = params.l1, params.l2
    g = params.g

    v = (m1 + m2) * g * l1 * (1.0 - np.cos(theta1)) + m2 * g * l2 * (1.0 - np.cos(theta2))
    return v


def energia_total(
    x: np.ndarray,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> Union[float, np.ndarray]:
    """Energía total conservada E = T + V."""
    return energia_cinetica(x, params) + energia_potencial(x, params)


def campo_derivadas(
    x: np.ndarray,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> np.ndarray:
    """Campo vectorial dx/dt = (omega1, omega2, alpha1, alpha2) resuelto con Cramer.
    
    Acepta tanto vector 1D (4,) como arreglos vectorizados (N, 4).
    """
    theta1 = x[..., 0]
    theta2 = x[..., 1]
    omega1 = x[..., 2]
    omega2 = x[..., 3]
    delta = theta1 - theta2

    m1, m2 = params.m1, params.m2
    l1, l2 = params.l1, params.l2
    g = params.g

    sin_delta = np.sin(delta)
    m11, m12, m21, m22, det_m = matriz_inercia(theta1, theta2, params)

    f1 = -m2 * l1 * l2 * sin_delta * (omega2 ** 2) - (m1 + m2) * g * l1 * np.sin(theta1)
    f2 = +m2 * l1 * l2 * sin_delta * (omega1 ** 2) - m2 * g * l2 * np.sin(theta2)

    # Regla de Cramer para M * alpha = f
    alpha1 = (f1 * m22 - f2 * m12) / det_m
    alpha2 = (m11 * f2 - m21 * f1) / det_m

    dxdt = np.empty_like(x)
    dxdt[..., 0] = omega1
    dxdt[..., 1] = omega2
    dxdt[..., 2] = alpha1
    dxdt[..., 3] = alpha2
    return dxdt


# ---------------------------------------------------------------------------
# Forma Hamiltoniana Canónica (theta, p)
# ---------------------------------------------------------------------------

def x_a_canonico(
    x: np.ndarray,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> np.ndarray:
    """Convierte el estado (theta, omega) a coordenadas canónicas (theta, p).
    
    p = M(theta) * omega
    """
    theta1 = x[..., 0]
    theta2 = x[..., 1]
    omega1 = x[..., 2]
    omega2 = x[..., 3]

    m11, m12, m21, m22, _ = matriz_inercia(theta1, theta2, params)
    p1 = m11 * omega1 + m12 * omega2
    p2 = m21 * omega1 + m22 * omega2

    z = np.empty_like(x)
    z[..., 0] = theta1
    z[..., 1] = theta2
    z[..., 2] = p1
    z[..., 3] = p2
    return z


def canonico_a_x(
    z: np.ndarray,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> np.ndarray:
    """Convierte coordenadas canónicas (theta, p) a estado (theta, omega).
    
    omega = M(theta)^(-1) * p
    """
    theta1 = z[..., 0]
    theta2 = z[..., 1]
    p1 = z[..., 2]
    p2 = z[..., 3]

    m11, m12, m21, m22, det_m = matriz_inercia(theta1, theta2, params)
    omega1 = (m22 * p1 - m12 * p2) / det_m
    omega2 = (-m21 * p1 + m11 * p2) / det_m

    x = np.empty_like(z)
    x[..., 0] = theta1
    x[..., 1] = theta2
    x[..., 2] = omega1
    x[..., 3] = omega2
    return x


def hamiltoniano(
    z: np.ndarray,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> Union[float, np.ndarray]:
    """Calcula H(theta, p) = 1/2 p^T M(theta)^(-1) p + V(theta)."""
    x = canonico_a_x(z, params)
    return energia_total(x, params)


def campo_hamiltoniano(
    z: np.ndarray,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> np.ndarray:
    """Ecuaciones canónicas de Hamilton para dz/dt = (theta_dot, p_dot).
    
    theta_dot = omega = M(theta)^(-1) * p
    p1_dot = -m2 * l1 * l2 * omega1 * omega2 * sin(Delta) - (m1 + m2) * g * l1 * sin(theta1)
    p2_dot = +m2 * l1 * l2 * omega1 * omega2 * sin(Delta) - m2 * g * l2 * sin(theta2)
    """
    theta1 = z[..., 0]
    theta2 = z[..., 1]
    p1 = z[..., 2]
    p2 = z[..., 3]
    delta = theta1 - theta2

    m1, m2 = params.m1, params.m2
    l1, l2 = params.l1, params.l2
    g = params.g

    m11, m12, m21, m22, det_m = matriz_inercia(theta1, theta2, params)
    omega1 = (m22 * p1 - m12 * p2) / det_m
    omega2 = (-m21 * p1 + m11 * p2) / det_m

    sin_delta = np.sin(delta)
    p1_dot = -m2 * l1 * l2 * omega1 * omega2 * sin_delta - (m1 + m2) * g * l1 * np.sin(theta1)
    p2_dot = +m2 * l1 * l2 * omega1 * omega2 * sin_delta - m2 * g * l2 * np.sin(theta2)

    dzdt = np.empty_like(z)
    dzdt[..., 0] = omega1
    dzdt[..., 1] = omega2
    dzdt[..., 2] = p1_dot
    dzdt[..., 3] = p2_dot
    return dzdt


# ---------------------------------------------------------------------------
# Umbrales de Energía
# ---------------------------------------------------------------------------

def umbrales_energia(params: ParametrosPendulo = PARAMS_DEFECTO) -> Dict[str, float]:
    """Calcula los umbrales de energía correspondientes a los puntos silla.
    
    Ec_m2 = V(0, pi) = 2 * m2 * g * l2
    Ec_m1 = V(pi, 0) = 2 * (m1 + m2) * g * l1
    Ec_ambas = V(pi, pi) = Ec_m1 + Ec_m2
    
    Nota: Cuál umbral es menor depende de si m2 * l2 < (m1 + m2) * l1.
    """
    m1, m2 = params.m1, params.m2
    l1, l2 = params.l1, params.l2
    g = params.g

    ec_m2 = 2.0 * m2 * g * l2
    ec_m1 = 2.0 * (m1 + m2) * g * l1
    ec_ambas = ec_m1 + ec_m2

    return {
        "Ec_m2": float(ec_m2),
        "Ec_m1": float(ec_m1),
        "Ec_ambas": float(ec_ambas),
        "primer_umbral": "Ec_m2" if ec_m2 < ec_m1 else "Ec_m1",
        "m2_l2": float(m2 * l2),
        "m_total_l1": float((m1 + m2) * l1),
    }


# ---------------------------------------------------------------------------
# Linealización en el origen y Modos Normales
# ---------------------------------------------------------------------------

def matriz_linealizacion(params: ParametrosPendulo = PARAMS_DEFECTO) -> Tuple[np.ndarray, np.ndarray]:
    """Calcula la matriz A del sistema linealizado x_dot = A * x en el origen.
    
    A = [[0, I], [-M0^(-1) * K, 0]]
    M0 = [[(m1+m2)*l1^2, m2*l1*l2], [m2*l1*l2, m2*l2^2]]
    K = diag((m1+m2)*g*l1, m2*g*l2)
    
    Retorna: (A, frecuencias_modos)
    """
    m1, m2 = params.m1, params.m2
    l1, l2 = params.l1, params.l2
    g = params.g

    m0 = np.array([
        [(m1 + m2) * (l1 ** 2), m2 * l1 * l2],
        [m2 * l1 * l2, m2 * (l2 ** 2)],
    ], dtype=np.float64)

    k = np.array([
        [(m1 + m2) * g * l1, 0.0],
        [0.0, m2 * g * l2],
    ], dtype=np.float64)

    m0_inv = np.linalg.inv(m0)
    omega_sq_mat = m0_inv @ k

    # Las frecuencias al cuadrado son los eigenvalores de M0^(-1) * K
    eigvals = np.linalg.eigvals(omega_sq_mat)
    # Ordenar de menor a mayor frecuencia
    eigvals = np.sort(np.real(eigvals))
    frecuencias = np.sqrt(eigvals)

    a = np.zeros((4, 4), dtype=np.float64)
    a[0:2, 2:4] = np.eye(2)
    a[2:4, 0:2] = -omega_sq_mat

    return a, frecuencias


# ---------------------------------------------------------------------------
# Equilibrios y Estabilidad
# ---------------------------------------------------------------------------

EQUILIBRIOS_NOMBRES = {
    "origen": (0.0, 0.0, 0.0, 0.0),
    "silla_pi_0": (np.pi, 0.0, 0.0, 0.0),
    "silla_0_pi": (0.0, np.pi, 0.0, 0.0),
    "silla_pi_pi": (np.pi, np.pi, 0.0, 0.0),
}


def jacobiano_numerico(
    x: np.ndarray,
    params: ParametrosPendulo = PARAMS_DEFECTO,
    eps: float = 1e-7,
) -> np.ndarray:
    """Calcula el Jacobiano df/dx por diferencias centrales en un estado x."""
    x = np.asarray(x, dtype=np.float64)
    j = np.zeros((4, 4), dtype=np.float64)
    for i in range(4):
        dx = np.zeros(4, dtype=np.float64)
        dx[i] = eps
        f_plus = campo_derivadas(x + dx, params)
        f_minus = campo_derivadas(x - dx, params)
        j[:, i] = (f_plus - f_minus) / (2.0 * eps)
    return j


def analizar_equilibrios(
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> Dict[str, Dict[str, Union[np.ndarray, str, float]]]:
    """Evalúa el campo en los 4 equilibrios y calcula sus eigenvalores."""
    resultados = {}
    for nombre, coord in EQUILIBRIOS_NOMBRES.items():
        x_eq = np.array(coord, dtype=np.float64)
        f_val = campo_derivadas(x_eq, params)
        j = jacobiano_numerico(x_eq, params)
        eigs = np.linalg.eigvals(j)

        # Clasificación según parte real e imaginaria
        partes_reales = np.real(eigs)
        partes_imag = np.imag(eigs)

        tol = 1e-5
        pares_reales = int(np.sum(np.abs(partes_imag) < tol) // 2)
        pares_imag = int(np.sum(np.abs(partes_reales) < tol) // 2)

        if nombre == "origen":
            tipo = "centro"
        elif nombre == "silla_pi_pi":
            tipo = "silla_doble_real"
        else:
            tipo = "silla_mixta"

        resultados[nombre] = {
            "coord": x_eq,
            "f_residual": float(np.linalg.norm(f_val)),
            "eigenvalores": eigs,
            "partes_reales": partes_reales,
            "partes_imaginarias": partes_imag,
            "tipo": tipo,
            "pares_reales": pares_reales,
            "pares_imaginarios": pares_imag,
        }
    return resultados
