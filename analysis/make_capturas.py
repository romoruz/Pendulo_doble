"""Script para generar todas las capturas de pantalla requeridas para el simulador web.

Fase 3c:
- Las 4 pestañas a 1920x1080 (proyector) y 1180x820 (iPad horizontal), tema oscuro.
- Péndulo a 1920x1080 en tema claro.
- Compatibilidad con nombres anteriores de proyector/ipad/laptop.
"""

import os
import shutil
import subprocess
import sys
import time
from http.server import HTTPServer, SimpleHTTPRequestHandler
import threading

DIR_RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_DOCS = os.path.join(DIR_RAIZ, "docs")
DIR_CAPTURAS = os.path.join(DIR_DOCS, "capturas")
DIR_PROFILE = "/tmp/ff_headless_profile"

os.makedirs(DIR_CAPTURAS, exist_ok=True)
os.makedirs(DIR_PROFILE, exist_ok=True)

PUERTO = 8098


class HandlerDocs(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=DIR_DOCS, **kwargs)

    def log_message(self, format, *args):
        # Silenciar logs HTTP
        pass


def arrancar_servidor():
    servidor = HTTPServer(("127.0.0.1", PUERTO), HandlerDocs)
    hilo = threading.Thread(target=servidor.serve_forever, daemon=True)
    hilo.start()
    return servidor


def capturar_pantalla(ancho: int, alto: int, tab: str, tema: str, ruta_salida: str):
    url = f"http://127.0.0.1:{PUERTO}/?tab={tab}&theme={tema}&capture=1"
    cmd = [
        "firefox",
        "--headless",
        "--no-remote",
        "--profile",
        DIR_PROFILE,
        "--screenshot",
        ruta_salida,
        f"--window-size={ancho},{alto}",
        url,
    ]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Error al capturar {ruta_salida}: {res.stderr}")
        return False
    tam_kb = os.path.getsize(ruta_salida) / 1024
    print(f"  [OK] {os.path.basename(ruta_salida)} ({ancho}x{alto}, {tab}, {tema}) -> {tam_kb:.1f} KB")
    return True


def main():
    print("Iniciando captura automatizada de vistas del simulador...")
    servidor = arrancar_servidor()
    time.sleep(1)

    capturas = [
        # 1920x1080 Oscuro (4 pestañas)
        (1920, 1080, "pendulo", "dark", "pendulo_1080_oscuro.png"),
        (1920, 1080, "gemelos", "dark", "gemelos_1080_oscuro.png"),
        (1920, 1080, "montecarlo", "dark", "montecarlo_1080_oscuro.png"),
        (1920, 1080, "mapacaos", "dark", "mapacaos_1080_oscuro.png"),

        # 1180x820 Oscuro (4 pestañas)
        (1180, 820, "pendulo", "dark", "pendulo_1180_oscuro.png"),
        (1180, 820, "gemelos", "dark", "gemelos_1180_oscuro.png"),
        (1180, 820, "montecarlo", "dark", "montecarlo_1180_oscuro.png"),
        (1180, 820, "mapacaos", "dark", "mapacaos_1180_oscuro.png"),

        # 1920x1080 Claro (Péndulo)
        (1920, 1080, "pendulo", "light", "pendulo_1080_claro.png"),

        # Formatos estándar para compatibilidad
        (1920, 1080, "pendulo", "dark", "proyector_1080_oscuro.png"),
        (1920, 1080, "pendulo", "light", "proyector_1080_claro.png"),
        (1180, 820, "pendulo", "dark", "ipad_horizontal_oscuro.png"),
        (1180, 820, "pendulo", "light", "ipad_horizontal_claro.png"),
        (1366, 768, "pendulo", "dark", "laptop_1366_oscuro.png"),
        (1366, 768, "pendulo", "light", "laptop_1366_claro.png"),
    ]

    for ancho, alto, tab, tema, nombre in capturas:
        ruta = os.path.join(DIR_CAPTURAS, nombre)
        capturar_pantalla(ancho, alto, tab, tema, ruta)

    # Limpiar archivos de prueba si existen
    for archivo in os.listdir(DIR_CAPTURAS):
        if archivo.startswith("test_"):
            try:
                os.remove(os.path.join(DIR_CAPTURAS, archivo))
            except OSError:
                pass

    servidor.shutdown()
    print("Todas las capturas se generaron correctamente.")


if __name__ == "__main__":
    main()
