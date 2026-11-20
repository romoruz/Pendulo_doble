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


def test_reiniciar_simulacion():
    """Verifica el comportamiento del botón Reiniciar y los requisitos de UI:
    
    1. Regresa a t = 0 con la condición inicial vigente.
    2. Limpia la estela (estelaM2 = []).
    3. Limpia la gráfica theta(t) (historialTiempo = []).
    4. Limpia el contador de vueltas (vueltasM1 = 0, vueltasM2 = 0).
    5. Queda en pausa (corriendo = false, botón en 'Reanudar').
    6. Sin sliders de velocidad inicial en el DOM (w1=0, w2=0 siempre).
    7. Exactamente 4 presets: pequenas_oscilaciones, modo_normal_en_fase, caotico, casi_invertido.
    8. Pestaña 'Sensibilidad' y panel 'Sensibilidad a condiciones iniciales'.
    """
    html_path = DIR_RAIZ / "docs" / "index.html"
    app_js_path = DIR_RAIZ / "docs" / "app.js"
    check_js_path = DIR_RAIZ / "docs" / "check_physics.js"

    html_content = html_path.read_text(encoding="utf-8")
    app_content = app_js_path.read_text(encoding="utf-8")
    check_content = check_js_path.read_text(encoding="utf-8")

    # 1. Sin sliders de velocidad inicial
    assert 'id="slider-w1"' not in html_content, "slider-w1 no debe existir en index.html"
    assert 'id="slider-w2"' not in html_content, "slider-w2 no debe existir en index.html"

    # 2. Presets exactos
    assert 'data-preset="pequenas_oscilaciones"' in html_content
    assert 'data-preset="modo_normal_en_fase"' in html_content
    assert 'data-preset="caotico"' in html_content
    assert 'data-preset="casi_invertido"' in html_content

    # 3. Pestaña Sensibilidad y títulos
    assert 'data-tab="gemelos" id="tab-gemelos">Sensibilidad</button>' in html_content
    assert "Sensibilidad a condiciones iniciales" in html_content

    # 4. Monte Carlo texto σ = 0.57°
    assert "Cada péndulo empieza en el ángulo nominal más un error normal de σ = 0.57°, como un error de medición" in html_content

    # 5. Código de reiniciarSimulacion en app.js
    assert "function reiniciarSimulacion()" in app_content
    assert "estado.corriendo = false;" in app_content
    assert "dom.btnPlayPause.textContent = 'Reanudar';" in app_content
    assert "estado.t = 0.0;" in app_content
    assert "estado.estelaM2 = [];" in app_content
    assert "estado.historialTiempo = [];" in app_content
    assert "estado.vueltasM1 = 0;" in app_content
    assert "estado.vueltasM2 = 0;" in app_content

    # 6. Verificación en Node.js mediante check_physics.js
    node_bin = shutil.which("node")
    if node_bin:
        res = subprocess.run(
            [node_bin, str(check_js_path)],
            cwd=str(DIR_RAIZ),
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0
        assert "[✓] Botón Reiniciar verificado" in res.stdout

