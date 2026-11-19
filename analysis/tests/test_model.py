"""Pruebas unitarias para el modelo físico del péndulo doble (Fase 1: Tests 1 y 7)."""

import numpy as np
import pytest

from analysis.pendulo.config import (
    ParametrosPendulo,
    PARAMS_DEFECTO,
    PARAMS_DESIGUALES,
    PARAMS_TEST3,
)
from analysis.pendulo.model import (
    matriz_inercia,
    analizar_equilibrios,
    umbrales_energia,
)


def test_1_m_simetrica_y_det_positivo():
    """Test 1: M es simétrica y det M > 0 en una malla de ángulos con varios juegos de parámetros."""
    juegos_parametros = [
        PARAMS_DEFECTO,
        PARAMS_DESIGUALES,
        PARAMS_TEST3,  # m1=3, m2=0.5, l1=0.7, l2=1.4
    ]

    th1_vals = np.linspace(-np.pi, np.pi, 25)
    th2_vals = np.linspace(-np.pi, np.pi, 25)
    th1_grid, th2_grid = np.meshgrid(th1_vals, th2_vals)

    for params in juegos_parametros:
        m11, m12, m21, m22, det_m = matriz_inercia(th1_grid, th2_grid, params)

        # 1. Simetría: M12 == M21 en toda la malla
        assert np.allclose(m12, m21, atol=1e-15), f"Falla de simetría en M para {params}"

        # 2. Positividad estricta de det M
        # det M >= m2 * l1^2 * l2^2 * m1 > 0
        min_det_teorico = params.m2 * (params.l1 ** 2) * (params.l2 ** 2) * params.m1
        assert np.all(det_m >= min_det_teorico - 1e-12), f"det M menor al mínimo teórico para {params}"
        assert np.all(det_m > 0), f"det M no es estrictamente positivo para {params}"


def test_7_equilibrios_y_estructura_eigenvalores():
    """Test 7: f se anula en los cuatro equilibrios y los eigenvalores tienen la estructura esperada."""
    juegos_parametros = [PARAMS_DEFECTO, PARAMS_DESIGUALES]

    for params in juegos_parametros:
        res = analizar_equilibrios(params)

        for nombre_eq, datos in res.items():
            # Residual del campo debe ser 0 en equilibrio
            assert datos["f_residual"] < 1e-12, f"f no se anula en {nombre_eq} para {params}"

            eigs = datos["eigenvalores"]
            partes_reales = np.real(eigs)
            partes_imag = np.imag(eigs)

            tol = 1e-4

            if nombre_eq == "origen":
                # Centro: 4 eigenvalores puramente imaginarios (2 pares conjugados +- i*omega)
                assert np.all(np.abs(partes_reales) < tol), f"Origen tiene parte real no nula: {eigs}"
                assert np.sum(np.abs(partes_imag) > tol) == 4, f"Origen no tiene 4 frecuencias puras: {eigs}"

            elif nombre_eq in ("silla_pi_0", "silla_0_pi"):
                # Silla con un par real (+- lambda) y un par imaginario (+- i*omega)
                reales_no_cero = np.sum(np.abs(partes_reales) > tol)
                imag_no_cero = np.sum(np.abs(partes_imag) > tol)
                assert reales_no_cero == 2, f"{nombre_eq} debe tener 2 eigs reales, tiene: {eigs}"
                assert imag_no_cero == 2, f"{nombre_eq} debe tener 2 eigs imaginarios, tiene: {eigs}"

            elif nombre_eq == "silla_pi_pi":
                # Silla con dos pares reales (+- lambda1, +- lambda2)
                assert np.all(np.abs(partes_imag) < tol), f"silla_pi_pi debe ser puramente real: {eigs}"
                assert np.sum(np.abs(partes_reales) > tol) == 4, f"silla_pi_pi no tiene 4 eigs reales: {eigs}"


def test_umbrales_energia_independencia_orden():
    """Verifica que el orden de los umbrales no se asume fijo y depende de m2*l2 vs (m1+m2)*l1."""
    # Caso 1: m2*l2 < (m1+m2)*l1 -> Ec_m2 < Ec_m1
    u1 = umbrales_energia(ParametrosPendulo(m1=1.0, m2=0.5, l1=1.0, l2=1.0))
    assert u1["Ec_m2"] < u1["Ec_m1"]
    assert u1["primer_umbral"] == "Ec_m2"

    # Caso 2: m2*l2 > (m1+m2)*l1 -> Ec_m1 < Ec_m2
    u2 = umbrales_energia(ParametrosPendulo(m1=1.0, m2=10.0, l1=1.0, l2=2.0))
    assert u2["Ec_m1"] < u2["Ec_m2"]
    assert u2["primer_umbral"] == "Ec_m1"
