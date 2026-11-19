"""Exportador de datos, presets y referencias para el simulador web.

Genera:
- docs/presets.js: Presets de condiciones iniciales y parámetros físicos utilizables en la UI.
- docs/reference.json: Tres simulaciones completas de alta precisión de RK4 3/8 (2000 pasos c/u)
  utilizadas por docs/check_physics.js para verificar que la física en JavaScript coincida
  bit a bit (error < 10^-9) con la implementación canónica en Python.
"""

import json
import os
import sys
from pathlib import Path
import numpy as np

DIR_RAIZ = str(Path(__file__).resolve().parent.parent)
if DIR_RAIZ not in sys.path:
    sys.path.insert(0, DIR_RAIZ)

from analysis.pendulo.config import (
    PARAMS_DEFECTO,
    PARAMS_DESIGUALES,
    PRESETS_CI,
    COLOR_THETA1,
    COLOR_THETA2,
)
from analysis.pendulo.model import (
    matriz_linealizacion,
    analizar_equilibrios,
    umbrales_energia,
)
from analysis.pendulo.integrate import integrar_rk4

DIR_DOCS = os.path.join(DIR_RAIZ, "docs")
os.makedirs(DIR_DOCS, exist_ok=True)


def generar_presets_js() -> None:
    """Genera docs/presets.js con compatibilidad universal (Browser y Node.js)."""
    print("-> Generando docs/presets.js...")
    params = PARAMS_DEFECTO
    _, frecs = matriz_linealizacion(params)
    u_def = umbrales_energia(params)
    u_des = umbrales_energia(PARAMS_DESIGUALES)

    # Calculamos autofunciones / razones de amplitud de modos normales
    # Para m1=m2=1, l1=l2=1: razón modo 1 es +sqrt(2), razón modo 2 es -sqrt(2)
    presets_completos = {
        "cerca_origen": {
            "nombre": PRESETS_CI["cerca_origen"]["nombre"],
            "descripcion": PRESETS_CI["cerca_origen"]["descripcion"],
            "params": {"m1": 1.0, "m2": 1.0, "l1": 1.0, "l2": 1.0, "g": 9.81},
            "x0": PRESETS_CI["cerca_origen"]["x0"],
        },
        "baja_energia": {
            "nombre": PRESETS_CI["baja_energia"]["nombre"],
            "descripcion": PRESETS_CI["baja_energia"]["descripcion"],
            "params": {"m1": 1.0, "m2": 1.0, "l1": 1.0, "l2": 1.0, "g": 9.81},
            "x0": PRESETS_CI["baja_energia"]["x0"],
        },
        "alta_energia": {
            "nombre": PRESETS_CI["alta_energia"]["nombre"],
            "descripcion": PRESETS_CI["alta_energia"]["descripcion"],
            "params": {"m1": 1.0, "m2": 1.0, "l1": 1.0, "l2": 1.0, "g": 9.81},
            "x0": PRESETS_CI["alta_energia"]["x0"],
        },
        "caso_1_informe": {
            "nombre": PRESETS_CI["caso_1_informe"]["nombre"],
            "descripcion": PRESETS_CI["caso_1_informe"]["descripcion"],
            "params": {"m1": 1.0, "m2": 1.0, "l1": 1.0, "l2": 1.0, "g": 9.81},
            "x0": PRESETS_CI["caso_1_informe"]["x0"],
        },
        "caso_2_informe": {
            "nombre": PRESETS_CI["caso_2_informe"]["nombre"],
            "descripcion": PRESETS_CI["caso_2_informe"]["descripcion"],
            "params": {"m1": 1.0, "m2": 1.0, "l1": 1.0, "l2": 1.0, "g": 9.81},
            "x0": PRESETS_CI["caso_2_informe"]["x0"],
        },
        "modo_1": {
            "nombre": "Modo Normal 1 (En fase)",
            "descripcion": f"Oscilación armónica simétrica de baja frecuencia (ω₁ = {frecs[0]:.2f} rad/s)",
            "params": {"m1": 1.0, "m2": 1.0, "l1": 1.0, "l2": 1.0, "g": 9.81},
            "x0": [0.08, float(0.08 * np.sqrt(2)), 0.0, 0.0],
        },
        "modo_2": {
            "nombre": "Modo Normal 2 (Contrafase)",
            "descripcion": f"Oscilación armónica antisimétrica de alta frecuencia (ω₂ = {frecs[1]:.2f} rad/s)",
            "params": {"m1": 1.0, "m2": 1.0, "l1": 1.0, "l2": 1.0, "g": 9.81},
            "x0": [0.08, float(-0.08 * np.sqrt(2)), 0.0, 0.0],
        },
    }

    contenido_js = f"""/**
 * Presets y configuraciones físicas para el simulador web del péndulo doble.
 * Generado automáticamente por analysis/export_web.py.
 */
(function (root, factory) {{
    if (typeof module === 'object' && module.exports) {{
        module.exports = factory();
    }} else {{
        root.PRESETS_PENDULO = factory();
    }}
}}(typeof self !== 'undefined' ? self : this, function () {{
    return {json.dumps(presets_completos, indent=2, ensure_ascii=False)};
}}));
"""
    ruta_presets = os.path.join(DIR_DOCS, "presets.js")
    with open(ruta_presets, "w", encoding="utf-8") as f:
        f.write(contenido_js)
    print(f"   -> Guardado en {ruta_presets}")


def generar_reference_json() -> None:
    """Genera docs/reference.json con 3 juegos de parámetros × 2000 pasos de RK4 3/8 (dt=0.005)."""
    print("-> Generando docs/reference.json (3 juegos de parámetros × 2000 pasos)...")
    from analysis.pendulo.config import PARAMS_TEST3
    dt = 0.005
    n_pasos = 2000
    t_max = dt * n_pasos

    casos = [
        {
            "id": "params_defecto",
            "nombre": "Juego 1: Parámetros Defecto (m1=1, m2=1, l1=1, l2=1)",
            "params_obj": PARAMS_DEFECTO,
            "params": {"m1": 1.0, "m2": 1.0, "l1": 1.0, "l2": 1.0, "g": 9.81},
            "x0": [0.8, 0.8, 0.0, 0.0],
        },
        {
            "id": "params_desiguales",
            "nombre": "Juego 2: Parámetros Desiguales (m1=1, m2=0.25, l1=1, l2=1.6)",
            "params_obj": PARAMS_DESIGUALES,
            "params": {"m1": 1.0, "m2": 0.25, "l1": 1.0, "l2": 1.6, "g": 9.81},
            "x0": [0.6, 0.6, 0.0, 0.0],
        },
        {
            "id": "params_test3",
            "nombre": "Juego 3: Parámetros Test 3 (m1=3, m2=0.5, l1=0.7, l2=1.4)",
            "params_obj": PARAMS_TEST3,
            "params": {"m1": 3.0, "m2": 0.5, "l1": 0.7, "l2": 1.4, "g": 9.81},
            "x0": [0.5, -0.5, 0.0, 0.0],
        },
    ]

    datos_ref = {
        "descripcion": "Soluciones numéricas de referencia generadas con Python para validación de physics.js",
        "dt": dt,
        "n_pasos": n_pasos,
        "t_max": t_max,
        "casos": [],
    }

    for caso in casos:
        p_obj = caso["params_obj"]
        x0_arr = np.array(caso["x0"], dtype=np.float64)
        t_arr, traj_arr = integrar_rk4(x0_arr, t_max=t_max, dt=dt, params=p_obj)

        datos_ref["casos"].append({
            "id": caso["id"],
            "nombre": caso["nombre"],
            "params": caso["params"],
            "x0": caso["x0"],
            "trayectoria": traj_arr.tolist(),
        })

    ruta_ref = os.path.join(DIR_DOCS, "reference.json")
    with open(ruta_ref, "w", encoding="utf-8") as f:
        json.dump(datos_ref, f, indent=None)  # Compacto para velocidad
    print(f"   -> Guardado en {ruta_ref} ({os.path.getsize(ruta_ref) / 1024:.1f} KB)")


if __name__ == "__main__":
    generar_presets_js()
    generar_reference_json()
    print("Exportación web completada exitosamente.")
