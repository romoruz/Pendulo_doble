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
