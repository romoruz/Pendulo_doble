"""Prueba de consistencia entre las estadísticas de JavaScript (docs/physics.js) y Python.

Verifica que:
1. El generador PRNG mulberry32 en Python y JS produzca secuencias idénticas bit a bit.
2. El bootstrap percentil para media y desviación estándar coincida con error < 1e-12.
3. El intervalo de puntuación de Wilson coincida exacto con la fórmula matemática (< 1e-12).
"""

import json
import math
import shutil
import subprocess
from pathlib import Path
import numpy as np
import pytest

DIR_RAIZ = Path(__file__).resolve().parent.parent.parent


def to_int32(x: int) -> int:
    x = x & 0xFFFFFFFF
    return x - 0x100000000 if x >= 0x80000000 else x


def imul(a: int, b: int) -> int:
    return to_int32(a * b)


def mulberry32_py(seed: int):
    s = to_int32(seed)

    def rng():
        nonlocal s
        s = to_int32(s + 0x6D2B79F5)
        t = s
        t_u15 = (t & 0xFFFFFFFF) >> 15
        t = imul(t ^ t_u15, t | 1)
        t_u7 = (t & 0xFFFFFFFF) >> 7
        t = to_int32(t ^ to_int32(t + imul(t ^ t_u7, t | 61)))
        t_u14 = (t & 0xFFFFFFFF) >> 14
        return ((t ^ t_u14) & 0xFFFFFFFF) / 4294967296.0

    return rng


def percentil_lineal_py(arr: np.ndarray, q: float) -> float:
    n = len(arr)
    if n == 0:
        return 0.0
    if n == 1:
        return float(arr[0])
    idx = (n - 1) * q
    lo = math.floor(idx)
    hi = math.ceil(idx)
    weight = idx - lo
    return float((1.0 - weight) * arr[lo] + weight * arr[hi])


def bootstrap_percentil_py(datos: np.ndarray, B: int, seed: int, confianza: float = 0.95):
    datos = np.asarray(datos, dtype=np.float64)
    N = len(datos)
    rng = mulberry32_py(seed)

    media_obs = float(np.mean(datos))
    std_obs = float(np.std(datos, ddof=1)) if N > 1 else 0.0

    medias_b = np.empty(B, dtype=np.float64)
    stds_b = np.empty(B, dtype=np.float64)
    b_sample = np.empty(N, dtype=np.float64)

    for b in range(B):
        s_b = 0.0
        for j in range(N):
            idx = math.floor(rng() * N)
            val = datos[idx]
            b_sample[j] = val
            s_b += val
        m_b = s_b / N
        medias_b[b] = m_b

        s_sq_b = 0.0
        for j in range(N):
            diff = b_sample[j] - m_b
            s_sq_b += diff * diff
        stds_b[b] = math.sqrt(s_sq_b / (N - 1)) if N > 1 else 0.0

    medias_b.sort()
    stds_b.sort()

    alpha = 1.0 - confianza
    q_bajo = alpha * 0.5
    q_alto = 1.0 - q_bajo

    return {
        "media": media_obs,
        "media_ci_bajo": percentil_lineal_py(medias_b, q_bajo),
        "media_ci_alto": percentil_lineal_py(medias_b, q_alto),
        "std": std_obs,
        "std_ci_bajo": percentil_lineal_py(stds_b, q_bajo),
        "std_ci_alto": percentil_lineal_py(stds_b, q_alto),
    }


def wilson_interval_py(k: int, n: int, z: float = 1.959963984540054):
    if n <= 0:
        return {"p": 0.0, "bajo": 0.0, "alto": 0.0}
    p = k / n
    z2 = z * z
    denom = n + z2
    p_tilde = (k + 0.5 * z2) / denom
    se = (z / denom) * math.sqrt((k * (n - k) / n) + (z2 * 0.25))
    return {
        "p": p,
        "bajo": max(0.0, p_tilde - se),
        "alto": min(1.0, p_tilde + se),
    }


def test_consistencia_estadisticas_web():
    """Compara docs/physics.js contra la implementación idéntica de referencia en Python."""
    node_bin = shutil.which("node")
    if not node_bin:
        pytest.skip("Node.js no está disponible; se salta test_consistencia_estadisticas_web.")

    # Generamos una muestra fija representativa de vueltas de la masa 2
    datos_muestra = [0] * 320 + [1] * 120 + [2] * 50 + [3] * 10
    B = 2000
    semilla = 1234567

    # Resultados en Python
    py_boot = bootstrap_percentil_py(datos_muestra, B, seed=semilla)
    k_exitos = sum(1 for x in datos_muestra if x > 0)
    n_total = len(datos_muestra)
    py_wilson = wilson_interval_py(k_exitos, n_total)

    # Ejecución en Node.js cargando docs/physics.js
    script_node = f"""
    const PenduloFisica = require('./docs/physics');
    const datos = {json.dumps(datos_muestra)};
    const B = {B};
    const semilla = {semilla};
    const rng = PenduloFisica.mulberry32(semilla);
    const bootRes = PenduloFisica.bootstrapPercentil(datos, B, rng);
    const wilsonRes = PenduloFisica.intervaloWilson({k_exitos}, {n_total});
    console.log(JSON.stringify({{ boot: bootRes, wilson: wilsonRes }}));
    """

    res = subprocess.run(
        [node_bin, "-e", script_node],
        cwd=str(DIR_RAIZ),
        capture_output=True,
        text=True,
        check=True,
    )

    js_out = json.loads(res.stdout)
    js_boot = js_out["boot"]
    js_wilson = js_out["wilson"]

    # Comparación de bootstrap percentil (< 1e-12)
    for clave in ["media", "media_ci_bajo", "media_ci_alto", "std", "std_ci_bajo", "std_ci_alto"]:
        diff = abs(py_boot[clave] - js_boot[clave])
        assert diff < 1e-12, f"Diferencia en {clave}: Python={py_boot[clave]}, JS={js_boot[clave]}, diff={diff}"

    # Comparación de intervalo Wilson (< 1e-12)
    for clave in ["p", "bajo", "alto"]:
        diff = abs(py_wilson[clave] - js_wilson[clave])
        assert diff < 1e-12, f"Diferencia en Wilson {clave}: Python={py_wilson[clave]}, JS={js_wilson[clave]}, diff={diff}"
