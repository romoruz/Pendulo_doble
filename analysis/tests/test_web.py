"""Prueba de consistencia entre la simulación web (JavaScript / Node.js) y Python."""

import shutil
import subprocess
import pytest
from pathlib import Path

DIR_RAIZ = Path(__file__).resolve().parent.parent.parent


def test_consistencia_javascript_python():
    """Ejecuta docs/check_physics.js en Node.js y verifica que el error sea < 1e-9.
    
    Si node no está instalado en el sistema, la prueba se salta elegantemente.
    """
    node_bin = shutil.which("node")
    if not node_bin:
        pytest.skip("Node.js no está disponible en el entorno; se salta la prueba de consistencia JS.")

    script_check = str(DIR_RAIZ / "docs" / "check_physics.js")
    res = subprocess.run(
        [node_bin, script_check],
        cwd=str(DIR_RAIZ),
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, f"Fallo en la prueba de consistencia Node.js vs Python:\n{res.stdout}\n{res.stderr}"
    assert "RESULTADO: TODAS LAS PRUEBAS PASARON EXITOSAMENTE" in res.stdout


def test_orientacion_pantalla_debajo_pivote():
    """Verifica que con theta1 = theta2 = 0, ambas masas queden estrictamente debajo del pivote en pantalla.
    
    En física (eje Y positivo hacia arriba):
        y1 = -l1 * cos(theta1)
        y2 = y1 - l2 * cos(theta2)
    En pantalla / canvas (eje Y positivo hacia abajo):
        y_pantalla = y_pivote + escala * l * cos(theta) = y_pivote - y_cart * escala
    Con theta1 = theta2 = 0:
        y_pantalla_1 = y_pivote + escala * l1 > y_pivote
        y_pantalla_2 = y_pantalla_1 + escala * l2 > y_pantalla_1
    """
    from analysis.pendulo.model import posiciones_cartesianas
    from analysis.pendulo.config import PARAMS_DEFECTO

    x1, y1, x2, y2 = posiciones_cartesianas(0.0, 0.0, PARAMS_DEFECTO)

    y_pivote = 200.0
    escala = 80.0

    y_pantalla_1 = y_pivote - y1 * escala
    y_pantalla_2 = y_pivote - y2 * escala

    assert y_pantalla_1 > y_pivote, f"Masa 1 ({y_pantalla_1}) debe quedar debajo del pivote ({y_pivote})"
    assert y_pantalla_2 > y_pantalla_1, f"Masa 2 ({y_pantalla_2}) debe quedar debajo de masa 1 ({y_pantalla_1})"

