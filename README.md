# Péndulo doble

Este repositorio contiene la implementación del análisis dinámico y la simulación numérica del péndulo doble, con énfasis en el estudio de sistemas hamiltonianos y la transición hacia el caos determinista. Proyecto desarrollado para la asignatura de Sistemas Dinámicos II en el ITAM.

## Descripción técnica

El proyecto aborda el sistema desde una perspectiva analítica y computacional:

1. Modelado Físico: Derivación de las ecuaciones de Euler-Lagrange y formulación del campo vectorial en el espacio de fases tetradimensional.
2. Métodos Numéricos: Implementación de algoritmos Runge-Kutta de cuarto orden (clásico y variante 3/8 de Kutta) para la integración de las EDOs.
3. Análisis Topológico: Visualización de retratos de fase en la variedad cilíndrica del espacio de fases.
4. Reducción Dimensional: Construcción de Secciones de Poincaré para identificar la ruptura de toros KAM y la emergencia de enredos homoclínicos.

## Requisitos

- Python 3.x
- NumPy
- Matplotlib

## Estructura del repositorio

- `simulator.py`: Aplicación interactiva con interfaz gráfica para la exploración de parámetros en tiempo real.
- `retratos_fase.py`: Script para la generación de gráficas de alta resolución y análisis de estabilidad.
- `main.pdf`: Documento técnico con la deducción analítica completa y discusión de resultados.

## Instrucciones de uso

Para ejecutar el simulador interactivo:
```bash
python simulator.py
