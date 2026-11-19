# Péndulo Doble Conservativo

Simulación, análisis dinámico y visualización interactiva del péndulo doble sin fricción.
Proyecto para Sistemas Dinámicos II, Departamento de Matemáticas Aplicadas, ITAM.

## Instalación y Configuración del Entorno

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Ejecución de Pruebas

```bash
pytest analysis/tests/
```

## Generación de Figuras y Resultados

```bash
python analysis/make_figures.py
```

## Exportación para Simulador Web

```bash
python analysis/export_web.py
```
