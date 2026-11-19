"""Pruebas unitarias para dinámica, linealización, vueltas y FTLE (Fase 1: Tests 6, 8, 10)."""

import numpy as np
import pytest
from scipy.linalg import expm

from analysis.pendulo.config import (
    PARAMS_DEFECTO,
)
from analysis.pendulo.model import (
    matriz_linealizacion,
    energia_total,
    umbrales_energia,
)
from analysis.pendulo.integrate import (
    integrar_rk4,
)
from analysis.pendulo.flips import (
    pasos_por_arriba_estado,
)
from analysis.pendulo.chaos import (
    calcular_ftle,
)


def test_6_linealizacion_pequena_amplitud():
    """Test 6: Con amplitud 1e-3 rad, la diferencia relativa entre trayectoria no lineal y expm(A t) x0 es < 1e-3 en [0, 10]."""
    params = PARAMS_DEFECTO
    a_mat, frecs = matriz_linealizacion(params)

    x0 = np.array([1e-3, 1e-3, 0.0, 0.0], dtype=np.float64)
    dt = 1e-3
    t_max = 10.0

    t_vals, traj_no_lineal = integrar_rk4(x0, t_max=t_max, dt=dt, params=params)

    # Evaluar la solución lineal expm(A * t) * x0 en cada tiempo (submuestreado cada 100 pasos para velocidad)
    paso_eval = 100
    indices_eval = np.arange(0, len(t_vals), paso_eval)

    norm_x0 = np.linalg.norm(x0)
    for idx in indices_eval:
        t = t_vals[idx]
        sol_lin = expm(a_mat * t) @ x0
        sol_no_lin = traj_no_lineal[idx]
        dif_rel = np.linalg.norm(sol_no_lin - sol_lin) / norm_x0
        assert dif_rel < 1e-3, f"A t={t:.2f}s, la diferencia relativa {dif_rel:.3e} excede 1e-3"


def test_8_cero_vueltas_bajo_umbral():
    """Test 8: Si E0 < Ec_m2, el conteo de pasos por arriba de la masa 2 es estrictamente 0 hasta T = 200 s."""
    params = PARAMS_DEFECTO
    umbrales = umbrales_energia(params)
    ec_m2 = umbrales["Ec_m2"]

    # Condición de baja energía partiendo del reposo
    x0 = np.array([0.6, 0.6, 0.0, 0.0], dtype=np.float64)
    e0 = energia_total(x0, params)
    assert e0 < ec_m2, f"E0 ({e0:.2f}) no es menor que Ec_m2 ({ec_m2:.2f})"

    dt = 0.005
    t_max = 200.0

    _, traj = integrar_rk4(x0, t_max=t_max, dt=dt, params=params)
    vueltas_sin_signo, vueltas_con_signo = pasos_por_arriba_estado(traj, masa=2)

    assert vueltas_sin_signo == 0, f"Se detectaron {vueltas_sin_signo} vueltas a pesar de E0 < Ec_m2"
    assert vueltas_con_signo == 0, f"El conteo neto de vueltas no es 0: {vueltas_con_signo}"


def test_10_ftle_regular_vs_caotico():
    """Test 10: Cerca del origen (0.05 rad) lambda_T < 0.2 con T=40s; con theta1=theta2=2.8 rad, lambda_T > 0.5."""
    params = PARAMS_DEFECTO
    t_max = 40.0
    dt = 0.01

    # Condición cuasi-regular cerca del origen
    x0_regular = np.array([0.05, 0.05, 0.0, 0.0], dtype=np.float64)
    lambda_regular = calcular_ftle(x0_regular, t_max=t_max, dt=dt, params=params)
    assert lambda_regular < 0.2, f"FTLE cuasi-regular {lambda_regular:.3f} no es < 0.2"

    # Condición fuertemente caótica
    x0_caotico = np.array([2.8, 2.8, 0.0, 0.0], dtype=np.float64)
    lambda_caotico = calcular_ftle(x0_caotico, t_max=t_max, dt=dt, params=params)
    assert lambda_caotico > 0.5, f"FTLE caótico {lambda_caotico:.3f} no es > 0.5"
