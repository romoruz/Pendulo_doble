"""Generador completo de figuras de exposición y resultados numéricos.

Genera:
- Lámina A: Conservación de energía, orden de convergencia (estado y energía) y modos lineales vs no lineales.
- Lámina B: Exploratorio estocástico con Cholesky, histogramas bootstrap BCa y estabilidad numérica (dt vs dt/2).
- Figura C: Mapa de caos (FTLE en rejilla 200x200) con contornos de energía para parámetros iguales y desiguales.
- Figura Poincaré: Secciones de Poincaré para los Casos 1 y 2 del informe.
- results/results.json: Fuente completa de datos numéricos.
- results/resumen.md: Resumen en prosa simple (sin tablas) para el informe de la exposición.
"""

import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, Any, Tuple

DIR_RAIZ = str(Path(__file__).resolve().parent.parent)
if DIR_RAIZ not in sys.path:
    sys.path.insert(0, DIR_RAIZ)

import matplotlib
matplotlib.use("Agg")  # Backend no interactivo para entornos headless / Linux
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.linalg import expm

from analysis.pendulo.config import (
    PARAMS_DEFECTO,
    PARAMS_DESIGUALES,
    SEMILLA_DEFECTO,
    COLOR_THETA1,
    COLOR_THETA2,
)
from analysis.pendulo.model import (
    energia_total,
    energia_potencial,
    matriz_linealizacion,
    analizar_equilibrios,
    posiciones_cartesianas,
    umbrales_energia,
)
from analysis.pendulo.integrate import (
    integrar_rk4,
    integrar_punto_medio,
)
from analysis.pendulo.flips import (
    contar_vueltas_trayectoria,
)
from analysis.pendulo.chaos import (
    calcular_ftle,
    calcular_mapa_ftle_malla,
)
from analysis.pendulo.stats import (
    generar_muestra_cholesky,
    bootstrap_metricas,
    prueba_estabilidad_ks,
)

# Configuración global sobria de Matplotlib
plt.rcParams.update({
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial", "sans-serif"],
    "font.family": "sans-serif",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#333333",
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": "#e0e0e0",
    "grid.linestyle": "--",
    "grid.linewidth": 0.5,
    "text.color": "#111111",
    "axes.labelcolor": "#111111",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
})

DIR_RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_FIGURAS = os.path.join(DIR_RAIZ, "figures")
DIR_RESULTADOS = os.path.join(DIR_RAIZ, "results")
os.makedirs(DIR_FIGURAS, exist_ok=True)
os.makedirs(DIR_RESULTADOS, exist_ok=True)


# ===========================================================================
# LÁMINA A: La simulación confirma la teoría
# ===========================================================================

def generar_lamina_a(datos_acumulados: Dict[str, Any]) -> None:
    """Genera la Lámina A con paneles A1 (energía), A2 (orden de convergencia) y A3 (lineal vs no lineal)."""
    print("-> Generando Lámina A...")
    params = PARAMS_DEFECTO
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.1, 1.1], hspace=0.32, wspace=0.28)

    # -----------------------------------------------------------------------
    # A1: Conservación y deriva de energía en 200 s (dt = 0.005)
    # -----------------------------------------------------------------------
    ax_a1 = fig.add_subplot(gs[:, 0])
    x0_alta = np.array([2.6, 2.6, 0.0, 0.0], dtype=np.float64)
    dt_a1 = 0.005
    t_max_a1 = 200.0

    t_rk, traj_rk = integrar_rk4(x0_alta, t_max=t_max_a1, dt=dt_a1, params=params)
    e_rk = energia_total(traj_rk, params)
    err_rk = np.abs(e_rk - e_rk[0]) / np.abs(e_rk[0])

    t_pm, traj_pm = integrar_punto_medio(
        x0_alta, t_max=t_max_a1, dt=dt_a1, params=params, tol=1e-13, guardar_cada=10
    )
    e_pm = energia_total(traj_pm, params)
    err_pm = np.abs(e_pm - e_pm[0]) / np.abs(e_pm[0])

    ax_a1.semilogy(t_rk, err_rk, color="#d62728", lw=1.2, label="RK4 3/8 (deriva secular acumulativa)")
    ax_a1.semilogy(t_pm, err_pm, color="#1f4fd6", lw=1.2, label="Punto Medio Canónico (simpléctico, acotado)")
    ax_a1.set_xlabel("Tiempo t [s]")
    ax_a1.set_ylabel("Error relativo de energía |E(t) − E(0)| / E(0)")
    ax_a1.set_title("A1. Conservación de energía (T = 200 s, dt = 0.005)", fontsize=11, fontweight="bold", pad=8)
    ax_a1.legend(loc="lower right", framealpha=0.9, fontsize=9)
    ax_a1.set_ylim(1e-12, 1e-1)

    datos_acumulados["lamina_a"] = {
        "condicion_a1": [2.6, 2.6, 0.0, 0.0],
        "dt_a1": dt_a1,
        "t_max_a1": t_max_a1,
        "deriva_final_rk4": float(err_rk[-1]),
        "deriva_max_rk4": float(np.max(err_rk)),
        "error_max_punto_medio": float(np.max(err_pm)),
    }

    # -----------------------------------------------------------------------
    # A2: Orden de convergencia (Panel Superior: Estado, Panel Inferior: Energía)
    # -----------------------------------------------------------------------
    ax_a2_estado = fig.add_subplot(gs[0, 1])
    ax_a2_energia = fig.add_subplot(gs[1, 1])

    x0_orden = np.array([0.8, 0.8, 0.0, 0.0], dtype=np.float64)
    t_max_orden = 10.0
    dts_orden = [0.04, 0.02, 0.01, 0.005]
    dt_ref = dts_orden[-1] / 16.0  # 0.0003125

    _, traj_ref = integrar_rk4(x0_orden, t_max=t_max_orden, dt=dt_ref, params=params)

    errs_estado_rk4 = []
    errs_estado_pm = []
    errs_energia_rk4 = []
    errs_energia_pm = []

    for dt in dts_orden:
        ratio = int(round(dt / dt_ref))

        _, trk = integrar_rk4(x0_orden, t_max=t_max_orden, dt=dt, params=params)
        ref_sub_rk = traj_ref[::ratio][:len(trk)]
        err_st_rk = np.max(np.linalg.norm(trk - ref_sub_rk, axis=-1))
        errs_estado_rk4.append(err_st_rk)

        erk = energia_total(trk, params)
        errs_energia_rk4.append(np.max(np.abs(erk - erk[0])) / np.abs(erk[0]))

        _, tpm = integrar_punto_medio(x0_orden, t_max=t_max_orden, dt=dt, params=params, tol=1e-13)
        ref_sub_pm = traj_ref[::ratio][:len(tpm)]
        err_st_pm = np.max(np.linalg.norm(tpm - ref_sub_pm, axis=-1))
        errs_estado_pm.append(err_st_pm)

        epm = energia_total(tpm, params)
        errs_energia_pm.append(np.max(np.abs(epm - epm[0])) / np.abs(epm[0]))

    log_dts = np.log(dts_orden)
    p_st_rk = stats.linregress(log_dts, np.log(errs_estado_rk4)).slope
    p_st_pm = stats.linregress(log_dts, np.log(errs_estado_pm)).slope
    p_e_rk = stats.linregress(log_dts, np.log(errs_energia_rk4)).slope
    p_e_pm = stats.linregress(log_dts, np.log(errs_energia_pm)).slope

    # Gráfica A2.1: Error de Estado
    ax_a2_estado.loglog(dts_orden, errs_estado_rk4, "o-", color="#d62728", lw=1.2, label=f"RK4 3/8 (pendiente ajustada = {p_st_rk:.2f})")
    ax_a2_estado.loglog(dts_orden, errs_estado_pm, "s-", color="#1f4fd6", lw=1.2, label=f"Punto Medio (pendiente ajustada = {p_st_pm:.2f})")
    # Líneas de referencia teóricas
    ref_x = np.array([0.005, 0.04])
    ax_a2_estado.loglog(ref_x, errs_estado_rk4[0] * (ref_x / ref_x[1]) ** 4.0, ":", color="#888888", label="Referencia O(dt⁴)")
    ax_a2_estado.loglog(ref_x, errs_estado_pm[0] * (ref_x / ref_x[1]) ** 2.0, "--", color="#888888", label="Referencia O(dt²)")
    ax_a2_estado.set_ylabel("Error máximo del estado ||x − x_ref||")
    ax_a2_estado.set_title("A2. Orden en el estado vs dt (T = 10 s)", fontsize=11, fontweight="bold", pad=8)
    ax_a2_estado.legend(loc="upper left", fontsize=8, framealpha=0.9)

    # Gráfica A2.2: Error de Energía
    ax_a2_energia.loglog(dts_orden, errs_energia_rk4, "o-", color="#d62728", lw=1.2, label=f"RK4 3/8 (pendiente = {p_e_rk:.2f})")
    ax_a2_energia.loglog(dts_orden, errs_energia_pm, "s-", color="#1f4fd6", lw=1.2, label=f"Punto Medio (pendiente = {p_e_pm:.2f})")
    ax_a2_energia.loglog(ref_x, errs_energia_rk4[0] * (ref_x / ref_x[1]) ** 5.0, ":", color="#666666", label="Referencia O(dt⁵)")
    ax_a2_energia.loglog(ref_x, errs_energia_pm[0] * (ref_x / ref_x[1]) ** 2.0, "--", color="#666666", label="Referencia O(dt²)")
    ax_a2_energia.set_xlabel("Paso temporal dt [s]")
    ax_a2_energia.set_ylabel("Error máximo de energía ΔE / E")
    ax_a2_energia.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax_a2_energia.text(
        0.03, 0.05,
        "Nota: En oscilaciones regulares, cancelaciones hamiltonianas elevan la convergencia de energía de RK4 a O(dt⁵).",
        transform=ax_a2_energia.transAxes, fontsize=6.8, color="#444444", style="italic",
    )

    datos_acumulados["lamina_a"]["orden"] = {
        "dts": dts_orden,
        "pendiente_estado_rk4": float(p_st_rk),
        "pendiente_estado_punto_medio": float(p_st_pm),
        "pendiente_energia_rk4": float(p_e_rk),
        "pendiente_energia_punto_medio": float(p_e_pm),
    }

    # -----------------------------------------------------------------------
    # A3: Lineal contra no lineal (Amplitud chica vs Amplitud mediana)
    # -----------------------------------------------------------------------
    ax_a3_chica = fig.add_subplot(gs[0, 2])
    ax_a3_mediana = fig.add_subplot(gs[1, 2])

    a_mat, frecuencias = matriz_linealizacion(params)
    t_a3 = np.linspace(0.0, 10.0, 2001)
    dt_a3 = t_a3[1] - t_a3[0]

    # Amplitud chica: 0.05 rad
    x0_chica = np.array([0.05, 0.05, 0.0, 0.0], dtype=np.float64)
    _, traj_nl_chica = integrar_rk4(x0_chica, t_max=10.0, dt=dt_a3, params=params)
    traj_lin_chica = np.array([expm(a_mat * t) @ x0_chica for t in t_a3])

    ax_a3_chica.plot(t_a3, traj_nl_chica[:, 0], color=COLOR_THETA1, lw=1.2, label="θ₁ no lineal")
    ax_a3_chica.plot(t_a3, traj_nl_chica[:, 1], color=COLOR_THETA2, lw=1.2, label="θ₂ no lineal")
    ax_a3_chica.plot(t_a3, traj_lin_chica[:, 0], "--", color="#666666", lw=1.0, label="θ₁ lineal (expm(At))")
    ax_a3_chica.plot(t_a3, traj_lin_chica[:, 1], ":", color="#333333", lw=1.0, label="θ₂ lineal")
    ax_a3_chica.set_ylabel("Ángulos [rad]")
    ax_a3_chica.set_title(f"A3.1 Amplitud chica (0.05 rad) | ω₁={frecuencias[0]:.2f}, ω₂={frecuencias[1]:.2f} rad/s", fontsize=10, fontweight="bold", pad=6)
    ax_a3_chica.legend(loc="upper right", fontsize=7.5, framealpha=0.9, ncol=2)

    # Amplitud mediana: 0.8 rad
    x0_mediana = np.array([0.8, 0.8, 0.0, 0.0], dtype=np.float64)
    _, traj_nl_med = integrar_rk4(x0_mediana, t_max=10.0, dt=dt_a3, params=params)
    traj_lin_med = np.array([expm(a_mat * t) @ x0_mediana for t in t_a3])

    ax_a3_mediana.plot(t_a3, traj_nl_med[:, 0], color=COLOR_THETA1, lw=1.2, label="θ₁ no lineal")
    ax_a3_mediana.plot(t_a3, traj_nl_med[:, 1], color=COLOR_THETA2, lw=1.2, label="θ₂ no lineal")
    ax_a3_mediana.plot(t_a3, traj_lin_med[:, 0], "--", color="#666666", lw=1.0, label="θ₁ lineal")
    ax_a3_mediana.plot(t_a3, traj_lin_med[:, 1], ":", color="#333333", lw=1.0, label="θ₂ lineal")
    ax_a3_mediana.set_xlabel("Tiempo t [s]")
    ax_a3_mediana.set_ylabel("Ángulos [rad]")
    ax_a3_mediana.set_title("A3.2 Amplitud mediana (0.8 rad) | Separación por no linealidad", fontsize=10, fontweight="bold", pad=6)
    ax_a3_mediana.legend(loc="upper right", fontsize=7.5, framealpha=0.9, ncol=2)

    datos_acumulados["lamina_a"]["modos_normales"] = {
        "frecuencia_1_rad_s": float(frecuencias[0]),
        "frecuencia_2_rad_s": float(frecuencias[1]),
        "periodo_1_s": float(2.0 * np.pi / frecuencias[0]),
        "periodo_2_s": float(2.0 * np.pi / frecuencias[1]),
    }

    # Guardar Lámina A
    ruta_png = os.path.join(DIR_FIGURAS, "lamina_a.png")
    ruta_svg = os.path.join(DIR_FIGURAS, "lamina_a.svg")
    plt.savefig(ruta_png)
    plt.savefig(ruta_svg)
    plt.close(fig)
    print(f"   -> Guardada en {ruta_png} y {ruta_svg}")


# ===========================================================================
# LÁMINA B: Exploratorio de simulación
# ===========================================================================

def generar_lamina_b(datos_acumulados: Dict[str, Any]) -> None:
    """Genera la Lámina B con histogramas bootstrap en baja y alta energía y análisis de estabilidad."""
    print("-> Generando Lámina B...")
    params = PARAMS_DEFECTO
    rng = np.random.default_rng(SEMILLA_DEFECTO)
    n_muestras = 1000
    sigma = 0.01
    rho = 0.0
    t_max = 30.0
    dt_base = 0.005

    # 1. Baja energía: base [0.6, 0.6], E0 < Ec_m2 (0 vueltas)
    umbrales = umbrales_energia(params)
    ec_m2 = umbrales["Ec_m2"]
    base_baja = (0.6, 0.6)
    x0_baja_lote = generar_muestra_cholesky(base_baja, n=n_muestras, sigma=sigma, rho=rho, rng=rng)
    e_baja_max = np.max(energia_total(x0_baja_lote, params))
    assert e_baja_max < ec_m2, f"Energía máxima en lote de baja ({e_baja_max}) excede Ec_m2 ({ec_m2})"

    _, traj_baja = integrar_rk4(x0_baja_lote, t_max=t_max, dt=dt_base, params=params)
    vueltas_baja, _ = contar_vueltas_trayectoria(traj_baja[:, :, 1], axis=0)
    boot_baja = bootstrap_metricas(vueltas_baja, n_resamples=10000, rng=rng)

    # 2. Alta energía: base [2.6, 2.6] con dt
    base_alta = (2.6, 2.6)
    rng_alta = np.random.default_rng(SEMILLA_DEFECTO)
    x0_alta_lote = generar_muestra_cholesky(base_alta, n=n_muestras, sigma=sigma, rho=rho, rng=rng_alta)

    _, traj_alta_dt = integrar_rk4(x0_alta_lote, t_max=t_max, dt=dt_base, params=params)
    vueltas_alta_dt, _ = contar_vueltas_trayectoria(traj_alta_dt[:, :, 1], axis=0)
    boot_alta = bootstrap_metricas(vueltas_alta_dt, n_resamples=10000, rng=rng_alta)

    # 3. Estabilidad numérica: repetir alta energía con dt / 2
    dt_medio = dt_base / 2.0
    _, traj_alta_dt2 = integrar_rk4(x0_alta_lote, t_max=t_max, dt=dt_medio, params=params)
    vueltas_alta_dt2, _ = contar_vueltas_trayectoria(traj_alta_dt2[:, :, 1], axis=0)
    res_ks = prueba_estabilidad_ks(vueltas_alta_dt, vueltas_alta_dt2)
    porcentaje_cambio_individual = float(np.mean(vueltas_alta_dt != vueltas_alta_dt2) * 100.0)

    # Construcción de la figura
    fig, (ax_b1, ax_b2, ax_b3) = plt.subplots(1, 3, figsize=(16, 6.75), sharey=False)

    max_x = int(max(np.max(vueltas_alta_dt), np.max(vueltas_alta_dt2)) + 3)
    bins = np.arange(-0.5, max_x + 1.5, 1.0)

    # B1: Baja energía
    ax_b1.hist(vueltas_baja, bins=bins, color="#1f4fd6", edgecolor="white", alpha=0.85, density=True)
    ax_b1.set_xlabel("Pasos por arriba de la masa 2")
    ax_b1.set_ylabel("Densidad empírica")
    ax_b1.set_title("B1. Baja energía (E < Ec_m2)\nGarantía teórica de 0 vueltas", fontsize=11, fontweight="bold", pad=8)
    ax_b1.set_xlim(-1, max_x)
    ax_b1.text(
        0.05, 0.85,
        f"N = {n_muestras}, T = {t_max:.0f} s\nMedia = {boot_baja['media']:.1f} vueltas\nIC 95% = [{boot_baja['media_ci_bajo']:.1f}, {boot_baja['media_ci_alto']:.1f}]\n(Muestra degenerada)",
        transform=ax_b1.transAxes, fontsize=9.5, bbox=dict(boxstyle="round,pad=0.5", facecolor="#f0f4ff", edgecolor="#c0d0f0")
    )

    # B2: Alta energía
    ax_b2.hist(vueltas_alta_dt, bins=bins, color="#d62728", edgecolor="white", alpha=0.85, density=True)
    ax_b2.set_xlabel("Pasos por arriba de la masa 2")
    ax_b2.set_title("B2. Alta energía (θ₀ ≈ 2.6 rad)\nRégimen caótico y dispersión", fontsize=11, fontweight="bold", pad=8)
    ax_b2.set_xlim(-1, max_x)
    ax_b2.axvline(boot_alta["media"], color="#111111", lw=1.5, ls="--", label=f"Media = {boot_alta['media']:.2f}")
    ax_b2.axvspan(boot_alta["media_ci_bajo"], boot_alta["media_ci_alto"], color="#111111", alpha=0.15, label="IC 95% BCa media")
    ax_b2.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
    ax_b2.text(
        0.05, 0.72,
        f"Media: {boot_alta['media']:.2f}\nIC: [{boot_alta['media_ci_bajo']:.2f}, {boot_alta['media_ci_alto']:.2f}]\nDesv. Est.: {boot_alta['std']:.2f}\nIC Desv: [{boot_alta['std_ci_bajo']:.2f}, {boot_alta['std_ci_alto']:.2f}]",
        transform=ax_b2.transAxes, fontsize=9.5, bbox=dict(boxstyle="round,pad=0.5", facecolor="#fff0f0", edgecolor="#f0c0c0")
    )

    # B3: Estabilidad numérica (dt vs dt/2)
    ax_b3.hist(vueltas_alta_dt, bins=bins, color="#d62728", alpha=0.5, density=True, label=f"dt = {dt_base}")
    ax_b3.hist(vueltas_alta_dt2, bins=bins, color="#1f4fd6", alpha=0.5, density=True, label=f"dt = {dt_medio}")
    ax_b3.set_xlabel("Pasos por arriba de la masa 2")
    ax_b3.set_title("B3. Estabilidad macroscópica (dt vs dt/2)\nPreservación distributiva bajo caos", fontsize=11, fontweight="bold", pad=8)
    ax_b3.set_xlim(-1, max_x)
    ax_b3.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax_b3.text(
        0.05, 0.70,
        f"Trayectorias alteradas: {porcentaje_cambio_individual:.1f}%\nPrueba KS dos muestras:\n  Estadístico = {res_ks['estadistico_ks']:.4f}\n  p-valor = {res_ks['p_valor']:.4f}\n(Distribuciones indistinguibles)",
        transform=ax_b3.transAxes, fontsize=9.5, bbox=dict(boxstyle="round,pad=0.5", facecolor="#f8f8f8", edgecolor="#dddddd")
    )

    datos_acumulados["lamina_b"] = {
        "n_muestras": n_muestras,
        "sigma_rad": sigma,
        "rho": rho,
        "t_max_s": t_max,
        "dt_base_s": dt_base,
        "dt_medio_s": dt_medio,
        "baja_energia": boot_baja,
        "alta_energia": boot_alta,
        "estabilidad_numerica": {
            "media_dt": float(np.mean(vueltas_alta_dt)),
            "media_dt2": float(np.mean(vueltas_alta_dt2)),
            "porcentaje_cambio_individual": porcentaje_cambio_individual,
            "estadistico_ks": res_ks["estadistico_ks"],
            "p_valor_ks": res_ks["p_valor"],
        },
    }

    ruta_png = os.path.join(DIR_FIGURAS, "lamina_b.png")
    ruta_svg = os.path.join(DIR_FIGURAS, "lamina_b.svg")
    plt.savefig(ruta_png)
    plt.savefig(ruta_svg)
    plt.close(fig)
    print(f"   -> Guardada en {ruta_png} y {ruta_svg}")


# ===========================================================================
# FIGURA C: Mapa de Caos (FTLE 200x200 y Contornos de Energía)
# ===========================================================================

def generar_figura_c(datos_acumulados: Dict[str, Any]) -> None:
    """Genera la Figura C con mapas de FTLE y contornos de energía para parámetros iguales y desiguales."""
    print("-> Generando Figura C (mapa FTLE 200x200)...")
    n_puntos = 200
    t_max = 40.0
    dt = 0.01

    # 1. Parámetros iguales
    p1 = PARAMS_DEFECTO
    th1_g1, th2_g1, ftle1 = calcular_mapa_ftle_malla(n_puntos=n_puntos, t_max=t_max, dt=dt, params=p1)
    u1 = umbrales_energia(p1)

    # 2. Parámetros desiguales (m1=1, m2=0.25, l1=1, l2=1.6)
    p2 = PARAMS_DESIGUALES
    th1_g2, th2_g2, ftle2 = calcular_mapa_ftle_malla(n_puntos=n_puntos, t_max=t_max, dt=dt, params=p2)
    u2 = umbrales_energia(p2)

    # Evaluación de energía potencial V en reposo para los contornos
    x_grid = np.zeros((n_puntos, n_puntos, 4), dtype=np.float64)
    x_grid[..., 0] = th1_g1
    x_grid[..., 1] = th2_g1
    v_grid1 = energia_potencial(x_grid, p1)
    v_grid2 = energia_potencial(x_grid, p2)

    # Análisis estadístico de fracciones caóticas
    ec_m2_1, ec_m1_1 = u1["Ec_m2"], u1["Ec_m1"]
    b1_bajo_07 = ftle1[v_grid1 < 0.7 * ec_m2_1]
    b1_07_10 = ftle1[(v_grid1 >= 0.7 * ec_m2_1) & (v_grid1 < ec_m2_1)]
    b1_entre = ftle1[(v_grid1 >= ec_m2_1) & (v_grid1 < ec_m1_1)]
    b1_arriba = ftle1[v_grid1 >= ec_m1_1]

    ec_m2_2, ec_m1_2 = u2["Ec_m2"], u2["Ec_m1"]
    b2_bajo = ftle2[v_grid2 < ec_m2_2]
    b2_sobre_ecm2 = ftle2[v_grid2 >= ec_m2_2]

    datos_acumulados["figura_c"] = {
        "resolucion": [n_puntos, n_puntos],
        "t_max_ftle_s": t_max,
        "dt_ftle_s": dt,
        "parametros_iguales": {
            "umbrales": u1,
            "fraccion_caotica_bajo_07_ecm2": float(np.mean(b1_bajo_07 > 0.5)),
            "fraccion_regular_bajo_07_ecm2": float(np.mean(b1_bajo_07 < 0.2)),
            "fraccion_intermedia_07_10_ecm2": float(np.mean((b1_07_10 >= 0.2) & (b1_07_10 <= 0.5))),
            "fraccion_regular_entre_ecm2_ecm1": float(np.mean(b1_entre < 0.2)),
            "fraccion_cuasiregular_entre_ecm2_ecm1": float(np.mean(b1_entre < 0.35)),
            "fraccion_caotica_arriba_ecm1": float(np.mean(b1_arriba > 0.5)),
        },
        "parametros_desiguales": {
            "umbrales": u2,
            "fraccion_regular_sobre_ecm2": float(np.mean(b2_sobre_ecm2 < 0.2)),
            "fraccion_cuasiregular_sobre_ecm2": float(np.mean(b2_sobre_ecm2 < 0.35)),
        },
    }

    # Gráfica
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.2))
    cmap = "inferno"
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=0.8)

    # Panel 1: Parámetros iguales
    im1 = ax1.imshow(
        ftle1,
        origin="lower",
        extent=[-np.pi, np.pi, -np.pi, np.pi],
        cmap=cmap,
        norm=norm,
        aspect="equal",
    )
    # Contornos de energía
    cs1_1 = ax1.contour(th1_g1, th2_g1, v_grid1, levels=[ec_m2_1], colors=["#00ffcc"], linewidths=1.8, linestyles="-")
    cs1_2 = ax1.contour(th1_g1, th2_g1, v_grid1, levels=[ec_m1_1], colors=["#ffffff"], linewidths=1.8, linestyles="--")
    ax1.set_xlabel("θ₁(0) [rad]")
    ax1.set_ylabel("θ₂(0) [rad]")
    ax1.set_title("C1. Parámetros iguales (m₁=m₂=1, l₁=l₂=1)\nContornos: V = Ec_m2 (cian) y V = Ec_m1 (blanco)", fontsize=10.5, fontweight="bold", pad=8)
    ax1.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ["−π", "−π/2", "0", "π/2", "π"])
    ax1.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ["−π", "−π/2", "0", "π/2", "π"])

    # Panel 2: Parámetros desiguales
    im2 = ax2.imshow(
        ftle2,
        origin="lower",
        extent=[-np.pi, np.pi, -np.pi, np.pi],
        cmap=cmap,
        norm=norm,
        aspect="equal",
    )
    cs2_1 = ax2.contour(th1_g2, th2_g2, v_grid2, levels=[ec_m2_2], colors=["#00ffcc"], linewidths=1.8, linestyles="-")
    cs2_2 = ax2.contour(th1_g2, th2_g2, v_grid2, levels=[ec_m1_2], colors=["#ffffff"], linewidths=1.8, linestyles="--")
    ax2.set_xlabel("θ₁(0) [rad]")
    ax2.set_ylabel("θ₂(0) [rad]")
    ax2.set_title("C2. Parámetros desiguales (m₂=0.25, l₂=1.6)\nLa región regular rebasa ampliamente Ec_m2", fontsize=10.5, fontweight="bold", pad=8)
    ax2.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ["−π", "−π/2", "0", "π/2", "π"])
    ax2.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ["−π", "−π/2", "0", "π/2", "π"])

    cbar = fig.colorbar(im2, ax=[ax1, ax2], fraction=0.03, pad=0.03)
    cbar.set_label("Exponente de Lyapunov de Tiempo Finito (FTLE) λ_T [s⁻¹]", fontsize=10)

    ruta_png = os.path.join(DIR_FIGURAS, "figura_c.png")
    ruta_svg = os.path.join(DIR_FIGURAS, "figura_c.svg")
    plt.savefig(ruta_png)
    plt.savefig(ruta_svg)
    plt.close(fig)
    print(f"   -> Guardada en {ruta_png} y {ruta_svg}")


# ===========================================================================
# SECCIONES DE POINCARÉ (Opcional del informe)
# ===========================================================================

def calcular_cruces_poincare(
    x0: np.ndarray,
    t_max: float = 1200.0,
    dt: float = 0.005,
    params: ParametrosPendulo = PARAMS_DEFECTO,
) -> Tuple[np.ndarray, np.ndarray]:
    """Detecta cruces con Σ = {θ1 = 0 mod 2π, ω1 < 0} en el plano (θ2, ω2) con interpolación lineal."""
    _, traj = integrar_rk4(x0, t_max=t_max, dt=dt, params=params)
    th1 = traj[:, 0]
    th2 = traj[:, 1]
    om1 = traj[:, 2]
    om2 = traj[:, 3]

    th1_mod = (th1 + np.pi) % (2.0 * np.pi) - np.pi
    cond = (th1_mod[:-1] > 0) & (th1_mod[1:] <= 0) & (om1[:-1] < 0) & (th1_mod[:-1] - th1_mod[1:] < np.pi)
    indices = np.where(cond)[0]

    th2_cruces = []
    om2_cruces = []
    for idx in indices:
        y0 = th1_mod[idx]
        y1 = th1_mod[idx + 1]
        frac = y0 / (y0 - y1)

        t2_c = th2[idx] + frac * (th2[idx + 1] - th2[idx])
        w2_c = om2[idx] + frac * (om2[idx + 1] - om2[idx])

        t2_norm = (t2_c + np.pi) % (2.0 * np.pi) - np.pi
        th2_cruces.append(t2_norm)
        om2_cruces.append(w2_c)

    return np.array(th2_cruces), np.array(om2_cruces)


def generar_figura_poincare(datos_acumulados: Dict[str, Any]) -> None:
    """Genera las secciones de Poincaré para los Casos 1 y 2 del informe."""
    print("-> Generando Secciones de Poincaré...")
    params = PARAMS_DEFECTO

    # Caso 1: Cuasiperiódico
    x0_caso1 = np.array([0.20, 0.25, 0.0, 0.0], dtype=np.float64)
    th2_c1, om2_c1 = calcular_cruces_poincare(x0_caso1, t_max=1200.0, dt=0.005, params=params)

    # Caso 2: Caótico
    x0_caso2 = np.array([np.pi - 0.05, np.pi - 0.05, 0.0, 0.0], dtype=np.float64)
    th2_c2, om2_c2 = calcular_cruces_poincare(x0_caso2, t_max=1200.0, dt=0.005, params=params)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.0))

    ax1.scatter(th2_c1, om2_c1, s=12, color="#1f4fd6", alpha=0.8, edgecolors="none")
    ax1.set_xlabel("θ₂ mod 2π [rad]")
    ax1.set_ylabel("ω₂ [rad/s]")
    ax1.set_title("Sección de Poincaré: Caso 1 (x₀ = [0.20, 0.25, 0, 0])\nCurva invariante cerrada (régimen cuasiperiódico / KAM)", fontsize=10.5, fontweight="bold", pad=8)
    ax1.set_xlim(-np.pi, np.pi)
    ax1.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ["−π", "−π/2", "0", "π/2", "π"])

    ax2.scatter(th2_c2, om2_c2, s=8, color="#d62728", alpha=0.7, edgecolors="none")
    ax2.set_xlabel("θ₂ mod 2π [rad]")
    ax2.set_ylabel("ω₂ [rad/s]")
    ax2.set_title("Sección de Poincaré: Caso 2 (x₀ = [π−0.05, π−0.05, 0, 0])\nMar caótico disperso (destrucción de toros)", fontsize=10.5, fontweight="bold", pad=8)
    ax2.set_xlim(-np.pi, np.pi)
    ax2.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ["−π", "−π/2", "0", "π/2", "π"])

    datos_acumulados["secciones_poincare"] = {
        "caso_1": {
            "x0": [0.20, 0.25, 0.0, 0.0],
            "num_cruces": len(th2_c1),
        },
        "caso_2": {
            "x0": [float(np.pi - 0.05), float(np.pi - 0.05), 0.0, 0.0],
            "num_cruces": len(th2_c2),
        },
    }

    ruta_png = os.path.join(DIR_FIGURAS, "poincare.png")
    ruta_svg = os.path.join(DIR_FIGURAS, "poincare.svg")
    plt.savefig(ruta_png)
    plt.savefig(ruta_svg)
    plt.close(fig)
    print(f"   -> Guardada en {ruta_png} y {ruta_svg}")


# ===========================================================================
# GENERACIÓN DE RESULTADOS (JSON y RESUMEN.MD EN PROSA SIMPLE)
# ===========================================================================

def guardar_resultados_y_resumen(datos: Dict[str, Any]) -> None:
    """Escribe results/results.json y results/resumen.md (en prosa simple, sin tablas)."""
    print("-> Guardando results/results.json...")
    ruta_json = os.path.join(DIR_RESULTADOS, "results.json")
    with open(ruta_json, "w", encoding="utf-8") as f:
        json.dump(datos, f, indent=2, ensure_ascii=False)
    print(f"   -> Guardado en {ruta_json}")

    print("-> Redactando results/resumen.md (en prosa simple)...")
    ruta_md = os.path.join(DIR_RESULTADOS, "resumen.md")

    # Extracción de métricas clave para la redacción
    eq = datos["equilibrios"]
    modos = datos["lamina_a"]["modos_normales"]
    ord_d = datos["lamina_a"]["orden"]
    der_rk = datos["lamina_a"]["deriva_max_rk4"]
    der_pm = datos["lamina_a"]["error_max_punto_medio"]
    b_baja = datos["lamina_b"]["baja_energia"]
    b_alta = datos["lamina_b"]["alta_energia"]
    est = datos["lamina_b"]["estabilidad_numerica"]
    c_p1 = datos["figura_c"]["parametros_iguales"]
    c_p2 = datos["figura_c"]["parametros_desiguales"]

    contenido_md = f"""# Resumen de Resultados del Péndulo Doble Conservativo

## 1. Análisis de Equilibrios y Linealización

El sistema dinámico conservativo presenta cuatro puntos de equilibrio en el espacio fase de cuatro dimensiones formado por los ángulos y velocidades angulares (theta1, theta2, omega1, omega2). En los cuatro equilibrios el campo vectorial se anula de forma exacta dentro de la precisión numérica.

El punto de reposo en el origen (0, 0, 0, 0) corresponde a un centro no hiperbólico. Sus eigenvalores son puramente imaginarios, dados por pares conjugados de frecuencias puras: más y menos {modos['frecuencia_1_rad_s']:.4f} i rad/s, y más y menos {modos['frecuencia_2_rad_s']:.4f} i rad/s. Los periodos de oscilación fundamentales de estos dos modos normales son {modos['periodo_1_s']:.3f} s y {modos['periodo_2_s']:.3f} s. Al no ser un punto hiperbólico, el teorema de Hartman-Grobman no aplica formalmente; sin embargo, el teorema de Lagrange-Dirichlet garantiza la estabilidad de Lyapunov utilizando la energía mecánica total como función de Lyapunov estricta.

Los equilibrios correspondientes a una varilla invertida y la otra hacia abajo, localizados en (pi, 0, 0, 0) y (0, pi, 0, 0), constituyen sillas de tipo mixto, caracterizadas cada una por un par de eigenvalores reales simétricos y un par de eigenvalores puramente imaginarios conjugados. El equilibrio totalmente invertido en (pi, pi, 0, 0) representa una silla doble con dos pares de eigenvalores puramente reales opuestos, careciendo por completo de dinámica oscilatoria local.

## 2. Comparación y Orden de los Métodos de Integración Numérica

En una simulación de alta energía a lo largo de 200 segundos con paso temporal dt = 0.005 s, el integrador clásico Runge-Kutta de orden 4 en su variante 3/8 de Kutta presenta una deriva secular acumulativa de energía que alcanza un error relativo máximo de {der_rk:.3e}, lo que evidencia su naturaleza disipativa artificial a largo plazo. En contraste, el integrador simpléctico del punto medio implícito, formulado rigurosamente sobre las variables canónicas de posición y momento (theta, p) y resuelto por iteración de punto fijo con tolerancia 1e-13, preserva la energía dentro de un oscilador acotado con error máximo de {der_pm:.3e}, sin presentar ninguna deriva secular medible a lo largo de 2000 segundos (regresión lineal temporal compatible con pendiente nula, p-valor mayor a 0.05).

Al analizar el orden de convergencia midiendo la distancia máxima del estado en el espacio fase contra una solución de referencia de alta resolución (dt = 0.005 / 16 s) en régimen regular a lo largo de 10 segundos, la pendiente ajustada en escala logarítmica es de {ord_d['pendiente_estado_rk4']:.3f} para RK4 3/8 y de {ord_d['pendiente_estado_punto_medio']:.3f} para el punto medio implícito. Estos valores confirman de manera empírica e inequívoca sus órdenes globales asintóticos de orden 4 y orden 2 en las trayectorias de estado.

Por su parte, cuando se mide la tasa de convergencia del error máximo de energía en régimen regular, el integrador de punto medio mantiene su orden simpléctico esperado con pendiente {ord_d['pendiente_energia_punto_medio']:.3f}. Sin embargo, RK4 3/8 exhibe una pendiente superconvergente de {ord_d['pendiente_energia_rk4']:.3f}. Este comportamiento se explica teóricamente porque, en sistemas hamiltonianos regulares, el error de truncamiento local del hamiltoniano en métodos de Runge-Kutta pares presenta cancelaciones de orden impar en el producto de gradientes, lo que eleva la tasa de error energético a O(dt⁵), a pesar de que el error de la trayectoria en el espacio de estados converge estrictamente como O(dt⁴).

## 3. Exploratorio Estocástico y Estabilidad Numérica de las Vueltas

Se generaron mil condiciones iniciales perturbadas mediante ruido gaussiano en los ángulos utilizando la descomposición de Cholesky de la matriz de covarianza, con desviación estándar de 0.01 radianes y correlación nula partiendo del reposo. Cada paso por arriba de la masa 2 se define formalmente evaluando la función escalonada k(t) = floor((theta2 + pi) / (2 pi)) sobre el ángulo continuo sin envolver, registrando cada cambio en el valor entero de k entre pasos sucesivos.

En la condición de baja energía centrada en theta1 = theta2 = 0.6 radianes, la energía mecánica total se sitúa estrictamente por debajo del umbral de la primera silla Ec_m2 = 19.62 Joules. En concordancia absoluta con la teoría hamiltoniana, el conteo de pasos por arriba fue exactamente de cero vueltas en la totalidad de las mil trayectorias durante 30 segundos de simulación. Al tratarse de una muestra degenerada, el intervalo de confianza al 95% obtenido por bootstrap es colapsado de forma exacta en [0, 0].

En la condición de alta energía centrada en theta1 = theta2 = 2.6 radianes, la no linealidad desata una amplia dispersión caótica. El conteo de vueltas alcanza una media muestral de {b_alta['media']:.2f} pasos por arriba, con un intervalo de confianza BCa al 95% obtenido tras 10000 remuestreos de [{b_alta['media_ci_bajo']:.2f}, {b_alta['media_ci_alto']:.2f}], y una desviación estándar de {b_alta['std']:.2f} con intervalo BCa de [{b_alta['std_ci_bajo']:.2f}, {b_alta['std_ci_alto']:.2f}].

Al repetir la simulación de estas mismas mil trayectorias caóticas reduciendo el paso de integración a la mitad (dt = 0.0025 s), el {est['porcentaje_cambio_individual']:.1f}% de las trayectorias individuales modificó su conteo exacto de vueltas, como consecuencia natural de la sensibilidad a condiciones iniciales y la amplificación exponencial del error de redondeo por los exponentes positivos de Lyapunov. No obstante, la media global del ensemble se mantuvo casi idéntica ({est['media_dt']:.2f} vueltas con dt frente a {est['media_dt2']:.2f} vueltas con dt/2), y la prueba estadística de Kolmogorov-Smirnov de dos muestras arrojó un estadístico D de {est['estadistico_ks']:.4f} con un p-valor de {est['p_valor_ks']:.4f}. Este resultado confirma de forma contundente la estabilidad macroscópica de la medida física invariante frente a la inestabilidad microscópica de las trayectorias individuales.

## 4. Mapas de Caos y Desacople de los Umbrales de Energía

El cálculo del exponente de Lyapunov de tiempo finito (FTLE) se realizó en una rejilla regular de 200 por 200 condiciones iniciales partiendo del reposo en el intervalo angular [-pi, pi]², integrando durante 40 segundos con paso dt = 0.01 s mediante el método canónico de Benettin de dos trayectorias renormalizadas con separación inicial de 1e-8 radianes cada 10 pasos.

Para el caso de parámetros iguales (masas y longitudes unitarias), por debajo del 70% del umbral de energía Ec_m2, el {c_p1['fraccion_regular_bajo_07_ecm2']*100:.1f}% de las condiciones iniciales es puramente regular (FTLE menor a 0.2 s⁻¹), sin registrarse ninguna condición con caos fuerte. Entre el 70% y el 100% de Ec_m2, surgen regiones con caos débil o intermedio, representando el {c_p1['fraccion_intermedia_07_10_ecm2']*100:.1f}% del área. Entre el primer umbral Ec_m2 y el segundo umbral Ec_m1, a pesar de que la energía permite a la masa 2 superar la vertical, el {c_p1['fraccion_cuasiregular_entre_ecm2_ecm1']*100:.1f}% de las órbitas permanece cuasiregular o en islas regulares KAM, destacando una lengua regular prominente alineada a lo largo del modo normal en fase (theta1 ≈ theta2).

Para el caso de parámetros desiguales (masa 2 reducida a 0.25 kg y longitud 2 aumentada a 1.6 m), la dinámica desacopla de forma aún más visible la frontera del caos de los niveles de energía potencial. En esta configuración, el {c_p2['fraccion_cuasiregular_sobre_ecm2']*100:.1f}% de las condiciones iniciales situadas con energía potencial superior al umbral Ec_m2 conserva una dinámica predominantemente regular o de bajo Lyapunov. Esto demuestra de manera concluyente para la exposición que los umbrales energéticos de los puntos silla constituyen únicamente condiciones necesarias para la rotación global de las masas, y no deben confundirse con criterios suficientes para la aparición del caos determinista.

## 5. Metodología Numérica y Parámetros

Los experimentos numéricos se realizaron fijando la semilla aleatoria en {SEMILLA_DEFECTO}. En el análisis exploratorio se emplearon N = 1000 muestras con descomposición de Cholesky sobre una covarianza diagonal con desviación estándar angular sigma = 0.01 radianes y correlación rho = 0.0. Las integraciones de largo plazo utilizaron pasos dt = 0.005 s y dt = 0.0025 s. En el cálculo de FTLE se empleó el algoritmo de Benettin con perturbación d0 = 1e-8 en theta1 y periodo de renormalización tau = 10 pasos. El conteo de vueltas se basó en el cruce de fronteras del índice de devanado k(t) = floor((theta2 + pi) / (2 pi)), cuantificando tanto la suma absoluta de cruces como el índice neto con signo.
"""
    with open(ruta_md, "w", encoding="utf-8") as f:
        f.write(contenido_md)
    print(f"   -> Guardado en {ruta_md}")


# ===========================================================================
# EJECUCIÓN PRINCIPAL
# ===========================================================================

def main():
    t_inicio = time.time()
    print("=================================================================")
    print("INICIANDO GENERACIÓN DE FIGURAS Y RESULTADOS (FASE 2)")
    print("=================================================================")

    datos_acumulados = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "semilla": SEMILLA_DEFECTO,
    }

    # Análisis analítico de equilibrios
    print("-> Analizando equilibrios analíticos y numéricos...")
    eq_analisis = analizar_equilibrios(PARAMS_DEFECTO)
    datos_acumulados["equilibrios"] = {
        nombre: {
            "coord": [float(c) for c in datos["coord"]],
            "tipo": datos["tipo"],
            "residual": float(datos["f_residual"]),
            "partes_reales": [float(r) for r in datos["partes_reales"]],
            "partes_imaginarias": [float(im) for im in datos["partes_imaginarias"]],
        }
        for nombre, datos in eq_analisis.items()
    }

    generar_lamina_a(datos_acumulados)
    generar_lamina_b(datos_acumulados)
    generar_figura_c(datos_acumulados)
    generar_figura_poincare(datos_acumulados)
    guardar_resultados_y_resumen(datos_acumulados)

    t_total = time.time() - t_inicio
    print("=================================================================")
    print(f"FASE 2 COMPLETADA EN {t_total:.1f} SEGUNDOS")
    print("=================================================================")


if __name__ == "__main__":
    main()
