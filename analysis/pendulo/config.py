"""Configuración central y presets para el péndulo doble.

Fuente única de verdad para parámetros físicos por defecto, presets de simulación,
semillas aleatorias y estilos visuales.
"""

from dataclasses import dataclass
import numpy as np


GRAVEDAD_DEFECTO = 9.81  # m/s^2
SEMILLA_DEFECTO = 20261119

# Colores institucionales / requeridos
COLOR_THETA1 = "#1f4fd6"  # Azul
COLOR_THETA2 = "#d62728"  # Rojo


@dataclass(frozen=True)
class ParametrosPendulo:
    """Parámetros físicos del péndulo doble."""
    m1: float = 1.0  # kg
    m2: float = 1.0  # kg
    l1: float = 1.0  # m
    l2: float = 1.0  # m
    g: float = GRAVEDAD_DEFECTO

    @property
    def mu(self) -> float:
        """Cociente de masas m2 / m1."""
        return self.m2 / self.m1

    @property
    def lambda_(self) -> float:
        """Cociente de longitudes l2 / l1."""
        return self.l2 / self.l1

    @property
    def tiempo_caracteristico(self) -> float:
        """Unidad de tiempo adimensional: sqrt(l1 / g)."""
        return np.sqrt(self.l1 / self.g)

    @property
    def energia_caracteristica(self) -> float:
        """Unidad de energía adimensional: (m1 + m2) * g * l1."""
        return (self.m1 + self.m2) * self.g * self.l1


# Parámetros por defecto para figuras y comparaciones
PARAMS_DEFECTO = ParametrosPendulo()

# Parámetros secundarios (Figura C y pruebas)
PARAMS_DESIGUALES = ParametrosPendulo(m1=1.0, m2=0.25, l1=1.0, l2=1.6)
PARAMS_TEST3 = ParametrosPendulo(m1=3.0, m2=0.5, l1=0.7, l2=1.4)

# Presets de condiciones iniciales x = [theta1, theta2, omega1, omega2]
PRESETS_CI = {
    "cerca_origen": {
        "nombre": "Cerca del origen",
        "descripcion": "Amplitud pequeña en régimen lineal",
        "x0": [0.05, 0.05, 0.0, 0.0],
    },
    "baja_energia": {
        "nombre": "Baja energía",
        "descripcion": "E < Ec_m2 (garantizado 0 vueltas)",
        "x0": [0.6, 0.6, 0.0, 0.0],
    },
    "alta_energia": {
        "nombre": "Alta energía",
        "descripcion": "E > Ec_m2 con comportamiento caótico y vueltas frecuentes",
        "x0": [2.6, 2.6, 0.0, 0.0],
    },
    "pequenas_oscilaciones": {
        "nombre": "Pequeñas oscilaciones",
        "descripcion": "Régimen lineal cercano (0.20, 0.25)",
        "x0": [0.20, 0.25, 0.0, 0.0],
    },
    "cuasi_periodico": {
        "nombre": "Cuasi-periódico",
        "descripcion": "Trayectoria regular en toro invariante (0.80, 0.60)",
        "x0": [0.80, 0.60, 0.0, 0.0],
    },
    "caotico": {
        "nombre": "Caótico",
        "descripcion": "Alta sensibilidad cerca de la cima invertida (π-0.05, π-0.05)",
        "x0": [float(np.pi - 0.05), float(np.pi - 0.05), 0.0, 0.0],
    },
    "caso_1_informe": {
        "nombre": "Caso 1 del informe",
        "descripcion": "Condición cuasiperiódica regular",
        "x0": [0.20, 0.25, 0.0, 0.0],
    },
    "caso_2_informe": {
        "nombre": "Caso 2 del informe",
        "descripcion": "Condición caótica cerca de la inversión vertical",
        "x0": [float(np.pi - 0.05), float(np.pi - 0.05), 0.0, 0.0],
    },
}
