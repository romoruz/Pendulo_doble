/**
 * Presets y configuraciones físicas para el simulador web del péndulo doble.
 * Generado automáticamente por analysis/export_web.py.
 */
(function (root, factory) {
    if (typeof module === 'object' && module.exports) {
        module.exports = factory();
    } else {
        root.PRESETS_PENDULO = factory();
    }
}(typeof self !== 'undefined' ? self : this, function () {
    return {
  "cerca_origen": {
    "nombre": "Cerca del origen",
    "descripcion": "Amplitud pequeña en régimen lineal",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      0.05,
      0.05,
      0.0,
      0.0
    ]
  },
  "baja_energia": {
    "nombre": "Baja energía",
    "descripcion": "E < Ec_m2 (garantizado 0 vueltas)",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      0.6,
      0.6,
      0.0,
      0.0
    ]
  },
  "alta_energia": {
    "nombre": "Alta energía",
    "descripcion": "E > Ec_m2 con comportamiento caótico y vueltas frecuentes",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      2.6,
      2.6,
      0.0,
      0.0
    ]
  },
  "pequenas_oscilaciones": {
    "nombre": "Pequeñas oscilaciones",
    "descripcion": "Régimen lineal cercano (0.20, 0.25)",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      0.2,
      0.25,
      0.0,
      0.0
    ]
  },
  "cuasi_periodico": {
    "nombre": "Cuasi-periódico",
    "descripcion": "Trayectoria regular en toro invariante (0.80, 0.60)",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      0.8,
      0.6,
      0.0,
      0.0
    ]
  },
  "caotico": {
    "nombre": "Caótico",
    "descripcion": "Alta sensibilidad cerca de la cima invertida (π-0.05, π-0.05)",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      3.0915926535897933,
      3.0915926535897933,
      0.0,
      0.0
    ]
  },
  "caso_1_informe": {
    "nombre": "Caso 1 del informe",
    "descripcion": "Condición cuasiperiódica regular",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      0.2,
      0.25,
      0.0,
      0.0
    ]
  },
  "caso_2_informe": {
    "nombre": "Caso 2 del informe",
    "descripcion": "Condición caótica cerca de la inversión vertical",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      3.0915926535897933,
      3.0915926535897933,
      0.0,
      0.0
    ]
  },
  "modo_1": {
    "nombre": "Modo Normal 1 (En fase)",
    "descripcion": "Oscilación armónica simétrica de baja frecuencia (ω₁ = 2.40 rad/s)",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      0.08,
      0.11313708498984762,
      0.0,
      0.0
    ]
  },
  "modo_2": {
    "nombre": "Modo Normal 2 (Contrafase)",
    "descripcion": "Oscilación armónica antisimétrica de alta frecuencia (ω₂ = 5.79 rad/s)",
    "params": {
      "m1": 1.0,
      "m2": 1.0,
      "l1": 1.0,
      "l2": 1.0,
      "g": 9.81
    },
    "x0": [
      0.08,
      -0.11313708498984762,
      0.0,
      0.0
    ]
  }
};
}));
