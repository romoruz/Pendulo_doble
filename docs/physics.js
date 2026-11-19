/**
 * physics.js - Motor de física y simulación del Péndulo Doble Conservativo.
 * Implementación canónica compartida para el navegador y Node.js.
 * 
 * Implementa:
 * - Parámetros por defecto y evaluación de energías (cinética, potencial, total).
 * - Campo vectorial lagrangiano resuelto vía regla de Cramer.
 * - Integrador Runge-Kutta 4 (variante 3/8 de Kutta) idéntico a Python.
 * - Conteo exacto de vueltas continuo k(t) = floor((theta2 + pi) / (2*pi)).
 * - Umbrales de energía potencial crítica (Ec_m2, Ec_m1, Ec_ambas).
 * - Matriz de linealización A y propagador lineal para modos normales.
 */

(function (root, factory) {
    if (typeof module === 'object' && module.exports) {
        module.exports = factory();
    } else {
        root.PenduloFisica = factory();
    }
}(typeof self !== 'undefined' ? self : this, function () {
    'use strict';

    var GRAVEDAD_DEFECTO = 9.81;

    function crearParametros(m1, m2, l1, l2, g) {
        return {
            m1: (m1 !== undefined) ? Number(m1) : 1.0,
            m2: (m2 !== undefined) ? Number(m2) : 1.0,
            l1: (l1 !== undefined) ? Number(l1) : 1.0,
            l2: (l2 !== undefined) ? Number(l2) : 1.0,
            g: (g !== undefined) ? Number(g) : GRAVEDAD_DEFECTO
        };
    }

    /**
     * Campo vectorial dx/dt = [omega1, omega2, alpha1, alpha2]
     * x = [theta1, theta2, omega1, omega2]
     */
    function campoDerivadas(x, p) {
        var th1 = x[0];
        var th2 = x[1];
        var w1 = x[2];
        var w2 = x[3];

        var m1 = p.m1;
        var m2 = p.m2;
        var l1 = p.l1;
        var l2 = p.l2;
        var g = p.g;

        var delta = th1 - th2;
        var cosDelta = Math.cos(delta);
        var sinDelta = Math.sin(delta);

        var m11 = (m1 + m2) * l1 * l1;
        var m12 = m2 * l1 * l2 * cosDelta;
        var m21 = m12;
        var m22 = m2 * l2 * l2;

        var detM = m2 * (l1 * l1) * (l2 * l2) * (m1 + m2 * sinDelta * sinDelta);

        var f1 = -m2 * l1 * l2 * sinDelta * (w2 * w2) - (m1 + m2) * g * l1 * Math.sin(th1);
        var f2 = m2 * l1 * l2 * sinDelta * (w1 * w1) - m2 * g * l2 * Math.sin(th2);

        var alpha1 = (f1 * m22 - f2 * m12) / detM;
        var alpha2 = (m11 * f2 - m21 * f1) / detM;

        return [w1, w2, alpha1, alpha2];
    }

    /**
     * Paso de integración RK4 3/8 de Kutta
     * x_next = x + dt/8 * (k1 + 3*k2 + 3*k3 + k4)
     */
    function pasoRK4(x, dt, p) {
        var dt_3 = dt / 3.0;
        var dt_8 = dt / 8.0;

        // k1
        var k1 = campoDerivadas(x, p);

        // k2 = f(x + dt/3 * k1)
        var x2 = [
            x[0] + dt_3 * k1[0],
            x[1] + dt_3 * k1[1],
            x[2] + dt_3 * k1[2],
            x[3] + dt_3 * k1[3]
        ];
        var k2 = campoDerivadas(x2, p);

        // k3 = f(x - dt/3 * k1 + dt * k2)
        var x3 = [
            x[0] - dt_3 * k1[0] + dt * k2[0],
            x[1] - dt_3 * k1[1] + dt * k2[1],
            x[2] - dt_3 * k1[2] + dt * k2[2],
            x[3] - dt_3 * k1[3] + dt * k2[3]
        ];
        var k3 = campoDerivadas(x3, p);

        // k4 = f(x + dt * k1 - dt * k2 + dt * k3)
        var x4 = [
            x[0] + dt * k1[0] - dt * k2[0] + dt * k3[0],
            x[1] + dt * k1[1] - dt * k2[1] + dt * k3[1],
            x[2] + dt * k1[2] - dt * k2[2] + dt * k3[2],
            x[3] + dt * k1[3] - dt * k2[3] + dt * k3[3]
        ];
        var k4 = campoDerivadas(x4, p);

        return [
            x[0] + dt_8 * (k1[0] + 3.0 * k2[0] + 3.0 * k3[0] + k4[0]),
            x[1] + dt_8 * (k1[1] + 3.0 * k2[1] + 3.0 * k3[1] + k4[1]),
            x[2] + dt_8 * (k1[2] + 3.0 * k2[2] + 3.0 * k3[2] + k4[2]),
            x[3] + dt_8 * (k1[3] + 3.0 * k2[3] + 3.0 * k3[3] + k4[3])
        ];
    }

    /**
     * Paso de integración RK4 clásico (regla de Simpson 1/6, 1/3, 1/3, 1/6)
     * x_next = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
     */
    function pasoRK4Clasico(x, dt, p) {
        var dt_2 = dt * 0.5;
        var dt_6 = dt / 6.0;

        var k1 = campoDerivadas(x, p);

        var x2 = [
            x[0] + dt_2 * k1[0],
            x[1] + dt_2 * k1[1],
            x[2] + dt_2 * k1[2],
            x[3] + dt_2 * k1[3]
        ];
        var k2 = campoDerivadas(x2, p);

        var x3 = [
            x[0] + dt_2 * k2[0],
            x[1] + dt_2 * k2[1],
            x[2] + dt_2 * k2[2],
            x[3] + dt_2 * k2[3]
        ];
        var k3 = campoDerivadas(x3, p);

        var x4 = [
            x[0] + dt * k3[0],
            x[1] + dt * k3[1],
            x[2] + dt * k3[2],
            x[3] + dt * k3[3]
        ];
        var k4 = campoDerivadas(x4, p);

        return [
            x[0] + dt_6 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
            x[1] + dt_6 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
            x[2] + dt_6 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
            x[3] + dt_6 * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3])
        ];
    }

    /**
     * Energía cinética T
     */
    function energiaCinetica(x, p) {
        var th1 = x[0], th2 = x[1], w1 = x[2], w2 = x[3];
        var delta = th1 - th2;
        return 0.5 * (p.m1 + p.m2) * (p.l1 * p.l1) * (w1 * w1) +
               0.5 * p.m2 * (p.l2 * p.l2) * (w2 * w2) +
               p.m2 * p.l1 * p.l2 * w1 * w2 * Math.cos(delta);
    }

    /**
     * Energía potencial V con referencia V(0,0) = 0
     */
    function energiaPotencial(x, p) {
        var th1 = x[0], th2 = x[1];
        return (p.m1 + p.m2) * p.g * p.l1 * (1.0 - Math.cos(th1)) +
               p.m2 * p.g * p.l2 * (1.0 - Math.cos(th2));
    }

    /**
     * Energía total E = T + V
     */
    function energiaTotal(x, p) {
        return energiaCinetica(x, p) + energiaPotencial(x, p);
    }

    /**
     * Umbrales de energía potencial crítica
     */
    function umbralesEnergia(p) {
        var ec_m2 = 2.0 * p.m2 * p.g * p.l2;
        var ec_m1 = 2.0 * (p.m1 + p.m2) * p.g * p.l1;
        var ec_ambas = ec_m1 + ec_m2;
        return {
            Ec_m2: ec_m2,
            Ec_m1: ec_m1,
            Ec_ambas: ec_ambas
        };
    }

    /**
     * Índice de vuelta continua para la masa 2
     * kappa(t) = floor((theta2 + pi) / (2*pi))
     */
    function indiceVuelta(theta2) {
        return Math.floor((theta2 + Math.PI) / (2.0 * Math.PI));
    }

    /**
     * Posiciones cartesianas de los dos extremos (x1, y1) y (x2, y2)
     * En coordenadas físicas (y hacia arriba, x hacia la derecha).
     */
    function posicionesCartesianas(x, p) {
        var th1 = x[0], th2 = x[1];
        var x1 = p.l1 * Math.sin(th1);
        var y1 = -p.l1 * Math.cos(th1);
        var x2 = x1 + p.l2 * Math.sin(th2);
        var y2 = y1 - p.l2 * Math.cos(th2);
        return { x1: x1, y1: y1, x2: x2, y2: y2 };
    }

    /**
     * Matriz Jacobiana de linealización alrededor del origen (0,0,0,0)
     * dx/dt = A * x
     */
    function linealizacion(p) {
        var mu = p.m2 / p.m1;
        var l1 = p.l1;
        var l2 = p.l2;
        var g = p.g;

        // A = [[0, 0, 1, 0],
        //      [0, 0, 0, 1],
        //      [-(1+mu)*g/l1, mu*g/l1, 0, 0],
        //      [(1+mu)*g/l2, -(1+mu)*g/l2, 0, 0]]
        var a20 = -(1.0 + mu) * g / l1;
        var a21 = mu * g / l1;
        var a30 = (1.0 + mu) * g / l2;
        var a31 = -(1.0 + mu) * g / l2;

        // Frecuencias propias resolviendo det(A - lambda*I) = 0 -> lambda^2 = -omega^2
        // omega^4 - Tr*omega^2 + Det = 0
        var b = -(a20 + a31);
        var c = a20 * a31 - a21 * a30;
        var disc = Math.sqrt(Math.max(0.0, b * b - 4.0 * c));
        var w1_sq = 0.5 * (b - disc);
        var w2_sq = 0.5 * (b + disc);

        var omega1 = Math.sqrt(Math.max(0.0, w1_sq));
        var omega2 = Math.sqrt(Math.max(0.0, w2_sq));

        return {
            A: [
                [0, 0, 1, 0],
                [0, 0, 0, 1],
                [a20, a21, 0, 0],
                [a30, a31, 0, 0]
            ],
            omega1: omega1,
            omega2: omega2
        };
    }

    /**
     * Paso temporal lineal para el "fantasma lineal"
     * dx/dt = A * x resuelto con RK4 (consistente)
     */
    function pasoLinealRK4(x, dt, linInfo) {
        var A = linInfo.A;
        function fLin(v) {
            return [
                v[2],
                v[3],
                A[2][0] * v[0] + A[2][1] * v[1],
                A[3][0] * v[0] + A[3][1] * v[1]
            ];
        }

        var dt_3 = dt / 3.0;
        var dt_8 = dt / 8.0;

        var k1 = fLin(x);
        var x2 = [x[0] + dt_3 * k1[0], x[1] + dt_3 * k1[1], x[2] + dt_3 * k1[2], x[3] + dt_3 * k1[3]];
        var k2 = fLin(x2);
        var x3 = [x[0] - dt_3 * k1[0] + dt * k2[0], x[1] - dt_3 * k1[1] + dt * k2[1], x[2] - dt_3 * k1[2] + dt * k2[2], x[3] - dt_3 * k1[3] + dt * k2[3]];
        var k3 = fLin(x3);
        var x4 = [x[0] + dt * k1[0] - dt * k2[0] + dt * k3[0], x[1] + dt * k1[1] - dt * k2[1] + dt * k3[1], x[2] + dt * k1[2] - dt * k2[2] + dt * k3[2], x[3] + dt * k1[3] - dt * k2[3] + dt * k3[3]];
        var k4 = fLin(x4);

        return [
            x[0] + dt_8 * (k1[0] + 3.0 * k2[0] + 3.0 * k3[0] + k4[0]),
            x[1] + dt_8 * (k1[1] + 3.0 * k2[1] + 3.0 * k3[1] + k4[1]),
            x[2] + dt_8 * (k1[2] + 3.0 * k2[2] + 3.0 * k3[2] + k4[2]),
            x[3] + dt_8 * (k1[3] + 3.0 * k2[3] + 3.0 * k3[3] + k4[3])
        ];
    }

    /**
     * Generador pseudoaleatorio Mulberry32 de 32 bits determinista y sembrable.
     * Devuelve función rng() que entrega flotantes uniformes en [0, 1).
     */
    function mulberry32(seed) {
        var s = seed | 0;
        return function () {
            s = (s + 0x6D2B79F5) | 0;
            var t = Math.imul(s ^ (s >>> 15), 1 | s);
            t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
            return ((t ^ (t >>> 14)) >>> 0) / 4294967296.0;
        };
    }

    /**
     * Percentil por interpolación lineal (método 'linear' idéntico a numpy.percentile).
     * arr debe estar ordenado ascendentemente. q in [0, 1].
     */
    function percentilLineal(arr, q) {
        var n = arr.length;
        if (n === 0) return 0;
        if (n === 1) return arr[0];
        var idx = (n - 1) * q;
        var lo = Math.floor(idx);
        var hi = Math.ceil(idx);
        var weight = idx - lo;
        return (1.0 - weight) * arr[lo] + weight * arr[hi];
    }

    /**
     * Bootstrap no paramétrico percentil para media y desviación estándar muestral.
     * datos: array numérico de tamaño N.
     * B: número de remuestreos (por defecto 2000).
     * rng: generador aleatorio uniforme [0, 1). Si es null/undefined, usa mulberry32(12345).
     * confianza: nivel de confianza (por defecto 0.95).
     */
    function bootstrapPercentil(datos, B, rng, confianza) {
        if (!datos || datos.length === 0) {
            return {
                media: 0, media_ci_bajo: 0, media_ci_alto: 0,
                std: 0, std_ci_bajo: 0, std_ci_alto: 0
            };
        }
        var N = datos.length;
        B = B || 2000;
        confianza = (confianza !== undefined) ? confianza : 0.95;
        if (!rng) {
            rng = mulberry32(12345);
        } else if (typeof rng === 'number') {
            rng = mulberry32(rng);
        }

        var sum = 0;
        for (var i = 0; i < N; i++) sum += datos[i];
        var mediaObs = sum / N;

        var sumSq = 0;
        for (var i = 0; i < N; i++) {
            var diff = datos[i] - mediaObs;
            sumSq += diff * diff;
        }
        var stdObs = (N > 1) ? Math.sqrt(sumSq / (N - 1)) : 0;

        if (stdObs === 0) {
            return {
                media: mediaObs, media_ci_bajo: mediaObs, media_ci_alto: mediaObs,
                std: 0, std_ci_bajo: 0, std_ci_alto: 0
            };
        }

        var mediasB = new Float64Array(B);
        var stdsB = new Float64Array(B);
        var bSample = new Float64Array(N);

        for (var b = 0; b < B; b++) {
            var sB = 0;
            for (var j = 0; j < N; j++) {
                var idx = Math.floor(rng() * N);
                var val = datos[idx];
                bSample[j] = val;
                sB += val;
            }
            var mB = sB / N;
            mediasB[b] = mB;

            var sSqB = 0;
            for (var j = 0; j < N; j++) {
                var d = bSample[j] - mB;
                sSqB += d * d;
            }
            stdsB[b] = (N > 1) ? Math.sqrt(sSqB / (N - 1)) : 0;
        }

        mediasB.sort();
        stdsB.sort();

        var alpha = 1.0 - confianza;
        var qBajo = alpha * 0.5;
        var qAlto = 1.0 - qBajo;

        return {
            media: mediaObs,
            media_ci_bajo: percentilLineal(mediasB, qBajo),
            media_ci_alto: percentilLineal(mediasB, qAlto),
            std: stdObs,
            std_ci_bajo: percentilLineal(stdsB, qBajo),
            std_ci_alto: percentilLineal(stdsB, qAlto)
        };
    }

    /**
     * Intervalo de puntuación de Wilson para proporción binaria k/n sin corrección de continuidad.
     * k: éxitos
     * n: total de ensayos
     * z: cuantil normal estándar (por defecto 1.959963984540054 para 95% de confianza)
     */
    function intervaloWilson(k, n, z) {
        if (n <= 0) {
            return { p: 0, bajo: 0, alto: 0 };
        }
        if (z === undefined || z === null) {
            z = 1.959963984540054;
        }
        var p = k / n;
        var z2 = z * z;
        var denom = n + z2;
        var pTilde = (k + 0.5 * z2) / denom;
        var se = (z / denom) * Math.sqrt((k * (n - k) / n) + (z2 * 0.25));
        var bajo = Math.max(0.0, pTilde - se);
        var alto = Math.min(1.0, pTilde + se);
        return {
            p: p,
            bajo: bajo,
            alto: alto
        };
    }

    return {
        crearParametros: crearParametros,
        campoDerivadas: campoDerivadas,
        pasoRK4: pasoRK4,
        pasoRK4Clasico: pasoRK4Clasico,
        energiaCinetica: energiaCinetica,
        energiaPotencial: energiaPotencial,
        energiaTotal: energiaTotal,
        umbralesEnergia: umbralesEnergia,
        indiceVuelta: indiceVuelta,
        posicionesCartesianas: posicionesCartesianas,
        linealizacion: linealizacion,
        pasoLinealRK4: pasoLinealRK4,
        mulberry32: mulberry32,
        percentilLineal: percentilLineal,
        bootstrapPercentil: bootstrapPercentil,
        intervaloWilson: intervaloWilson
    };
}));
