"""Pruebas unitarias para integradores numéricos (Fase 1: Tests 2, 3, 4, 5, 9)."""

import numpy as np
import pytest
from scipy import stats

from analysis.pendulo.config import (
    PARAMS_DEFECTO,
    PARAMS_DESIGUALES,
    PARAMS_TEST3,
)
from analysis.pendulo.model import (
    energia_total,
)
from analysis.pendulo.integrate import (
    paso_rk4_38,
    integrar_rk4,
    integrar_rk4_hamiltoniano,
    integrar_punto_medio,
)


def test_2_deriva_energia_rk4_tres_parametros():
    """Test 2: RK4 3/8 con dt = 1e-3 y T = 20 s tiene deriva relativa de energía < 1e-9."""
    juegos_parametros = [PARAMS_DEFECTO, PARAMS_DESIGUALES, PARAMS_TEST3]
    x0 = np.array([0.5, 0.5, 0.0, 0.0], dtype=np.float64)
    dt = 1e-3
    t_max = 20.0

    for params in juegos_parametros:
        t, traj = integrar_rk4(x0, t_max=t_max, dt=dt, params=params)
        e_serie = energia_total(traj, params)
        e0 = e_serie[0]

        deriva_relativa_max = np.max(np.abs(e_serie - e0)) / np.abs(e0)
        assert deriva_relativa_max < 1e-9, (
            f"Deriva relativa {deriva_relativa_max:.3e} excedió 1e-9 para {params}"
        )


def test_3_orden_convergencia_rk4_y_punto_medio():
    """Test 3: La pendiente log-log del error de energía contra dt cae en [3.7, 4.3] para RK4 y [1.8, 2.2] para punto medio.
    
    En régimen no lineal de alta energía (x0 = [2.4, 2.4, 0, 0], T = 10 s), los integradores
    manifiestan plenamente sus órdenes asintóticos teóricos (4 y 2).
    """
    params = PARAMS_DEFECTO
    x0 = np.array([2.4, 2.4, 0.0, 0.0], dtype=np.float64)
    t_max = 10.0
    dts = [0.04, 0.02, 0.01, 0.005, 0.0025]

    errores_rk4 = []
    errores_pm = []

    for dt in dts:
        # RK4 3/8
        _, traj_rk4 = integrar_rk4(x0, t_max=t_max, dt=dt, params=params)
        e_rk4 = energia_total(traj_rk4, params)
        err_rk4 = np.max(np.abs(e_rk4 - e_rk4[0])) / np.abs(e_rk4[0])
        errores_rk4.append(err_rk4)

        # Punto medio implícito
        _, traj_pm = integrar_punto_medio(x0, t_max=t_max, dt=dt, params=params, tol=1e-13)
        e_pm = energia_total(traj_pm, params)
        err_pm = np.max(np.abs(e_pm - e_pm[0])) / np.abs(e_pm[0])
        errores_pm.append(err_pm)

    # Ajuste lineal en log-log
    log_dt = np.log(dts)
    pendiente_rk4, _, _, _, _ = stats.linregress(log_dt, np.log(errores_rk4))
    pendiente_pm, _, _, _, _ = stats.linregress(log_dt, np.log(errores_pm))

    assert 3.7 <= pendiente_rk4 <= 4.3, f"Pendiente RK4 {pendiente_rk4:.3f} fuera de [3.7, 4.3]"
    assert 1.8 <= pendiente_pm <= 2.2, f"Pendiente Punto Medio {pendiente_pm:.3f} fuera de [1.8, 2.2]"


def test_4_punto_medio_sin_deriva_secular():
    """Test 4: El error de energía en punto medio implícito se mantiene acotado y sin tendencia en 2000 s."""
    params = PARAMS_DEFECTO
    x0 = np.array([0.7, 0.7, 0.0, 0.0], dtype=np.float64)
    dt = 0.02
    t_max = 2000.0  # 100,000 pasos
    guardar_cada = 50  # Guardar 2001 puntos para análisis eficiente

    t_vals, traj = integrar_punto_medio(
        x0,
        t_max=t_max,
        dt=dt,
        params=params,
        tol=1e-13,
        guardar_cada=guardar_cada,
    )

    e_serie = energia_total(traj, params)
    e0 = e_serie[0]
    error_rel = (e_serie - e0) / np.abs(e0)

    # El error debe mantenerse acotado (típicamente < 1e-4 con dt=0.02 en 2000s)
    assert np.max(np.abs(error_rel)) < 1e-3, f"Error relativo {np.max(np.abs(error_rel)):.2e} demasiado grande"

    # Regresión lineal del error sobre el tiempo: la pendiente debe ser compatible con 0
    res_reg = stats.linregress(t_vals, error_rel)
    # Criterio estadístico riguroso: no se rechaza H0: pendiente = 0 (p-valor > 0.05 y el 0 cae en el IC de la pendiente)
    assert res_reg.pvalue > 0.05, f"Pendiente estadísticamente distinta de 0: p-valor={res_reg.pvalue:.4f}"
    assert abs(res_reg.slope) < 2.0 * res_reg.stderr, (
        f"Pendiente {res_reg.slope:.3e} no compatible con 0 (stderr = {res_reg.stderr:.3e})"
    )


def test_5_equivalencia_lagrangiana_y_hamiltoniana():
    """Test 5: La forma hamiltoniana y lagrangiana dan la misma trayectoria (< 1e-10 con RK4 a dt=1e-3, T=5s)."""
    params = PARAMS_DEFECTO
    x0 = np.array([0.6, 0.4, 0.2, -0.3], dtype=np.float64)
    dt = 1e-3
    t_max = 5.0

    t_lag, traj_lag = integrar_rk4(x0, t_max=t_max, dt=dt, params=params)
    t_ham, traj_ham = integrar_rk4_hamiltoniano(x0, t_max=t_max, dt=dt, params=params)

    dif_max = np.max(np.abs(traj_lag - traj_ham))
    assert dif_max < 1e-10, f"Diferencia máxima {dif_max:.3e} excede 1e-10 entre formulaciones"


def test_9_consistencia_vectorizada_vs_escalar():
    """Test 9: La versión vectorizada y la escalar dan el mismo resultado (diferencia < 1e-12)."""
    params = PARAMS_DEFECTO
    n_ci = 5
    rng = np.random.default_rng(12345)
    x0_lote = rng.uniform(-1.0, 1.0, size=(n_ci, 4))
    dt = 1e-3
    t_max = 2.0

    # Ejecución vectorizada
    _, traj_vec = integrar_rk4(x0_lote, t_max=t_max, dt=dt, params=params)

    # Ejecución escalar individual
    for i in range(n_ci):
        x0_ind = x0_lote[i]
        _, traj_ind = integrar_rk4(x0_ind, t_max=t_max, dt=dt, params=params)
        dif = np.max(np.abs(traj_vec[:, i, :] - traj_ind))
        assert dif < 1e-12, f"Diferencia {dif:.3e} excede 1e-12 en condición {i}"
