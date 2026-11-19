"""Utilidades estadísticas y muestreo estocástico para el péndulo doble.

Incluye:
- Muestreo correlacionado de condiciones iniciales mediante descomposición de Cholesky
- Bootstrap no paramétrico BCa con scipy.stats.bootstrap y manejo riguroso de muestras degeneradas
- Prueba de Kolmogorov-Smirnov de dos muestras para comparar distribuciones
"""

from typing import Dict, Any, Optional, Tuple
import numpy as np
from scipy import stats
from .config import SEMILLA_DEFECTO


def generar_muestra_cholesky(
    base_angulos: Tuple[float, float],
    n: int = 1000,
    sigma: float = 0.01,
    rho: float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Genera N condiciones iniciales x = (theta1, theta2, 0, 0) con ruido gaussiano en ángulos.
    
    theta0 = base + L * z, con z ~ N(0, I_2) y L = chol(Sigma).
    Sigma = [[sigma^2, rho*sigma^2], [rho*sigma^2, sigma^2]].
    """
    if rng is None:
        rng = np.random.default_rng(SEMILLA_DEFECTO)

    cov = np.array([
        [sigma ** 2, rho * (sigma ** 2)],
        [rho * (sigma ** 2), sigma ** 2],
    ], dtype=np.float64)

    # Descomposición de Cholesky Sigma = L @ L.T
    l_chol = np.linalg.cholesky(cov)

    z = rng.standard_normal(size=(n, 2))
    delta_theta = z @ l_chol.T

    x0_lote = np.zeros((n, 4), dtype=np.float64)
    x0_lote[:, 0] = base_angulos[0] + delta_theta[:, 0]
    x0_lote[:, 1] = base_angulos[1] + delta_theta[:, 1]
    # omega1 y omega2 permanecen en 0.0

    return x0_lote


def bootstrap_metricas(
    datos: np.ndarray,
    n_resamples: int = 10000,
    confianza: float = 0.95,
    rng: Optional[np.random.Generator] = None,
) -> Dict[str, Any]:
    """Calcula la media, desviación estándar e intervalos de confianza BCa mediante bootstrap.
    
    Si la muestra es degenerada (todos los valores idénticos, e.g. todos ceros),
    no llama al método BCa de scipy (que fallaría o arrojaría advertencias numéricas),
    sino que reporta analíticamente IC = [val, val] y anota 'degenerada': True.
    """
    datos = np.asarray(datos, dtype=np.float64)
    media_obs = float(np.mean(datos))
    std_obs = float(np.std(datos, ddof=1)) if len(datos) > 1 else 0.0

    # Comprobar degeneración (varianza nula o todos los elementos iguales)
    if np.all(datos == datos[0]) or np.isclose(std_obs, 0.0):
        return {
            "degenerada": True,
            "media": media_obs,
            "media_ci_bajo": media_obs,
            "media_ci_alto": media_obs,
            "std": 0.0,
            "std_ci_bajo": 0.0,
            "std_ci_alto": 0.0,
        }

    random_state = rng if rng is not None else np.random.default_rng(SEMILLA_DEFECTO)

    # Bootstrap BCa para la media
    res_media = stats.bootstrap(
        (datos,),
        np.mean,
        confidence_level=confianza,
        n_resamples=n_resamples,
        method="BCa",
        random_state=random_state,
    )

    # Bootstrap BCa para la desviación estándar
    res_std = stats.bootstrap(
        (datos,),
        np.std,
        confidence_level=confianza,
        n_resamples=n_resamples,
        method="BCa",
        random_state=random_state,
    )

    return {
        "degenerada": False,
        "media": media_obs,
        "media_ci_bajo": float(res_media.confidence_interval.low),
        "media_ci_alto": float(res_media.confidence_interval.high),
        "std": std_obs,
        "std_ci_bajo": float(res_std.confidence_interval.low),
        "std_ci_alto": float(res_std.confidence_interval.high),
    }


def prueba_estabilidad_ks(
    muestra_dt: np.ndarray,
    muestra_dt_medio: np.ndarray,
) -> Dict[str, float]:
    """Compara las distribuciones de pasos por arriba con dt y dt/2 mediante Kolmogorov-Smirnov."""
    res_ks = stats.ks_2samp(muestra_dt, muestra_dt_medio)
    return {
        "estadistico_ks": float(res_ks.statistic),
        "p_valor": float(res_ks.pvalue),
    }
