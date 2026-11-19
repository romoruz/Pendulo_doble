/**
 * app.js - Controlador Principal del Simulador del Péndulo Doble Conservativo.
 * Laboratorio de Sistemas Dinámicos II | ITAM
 * 
 * Implementa:
 * - Tema oscuro por defecto y alternancia a tema claro con paleta estricta.
 * - Integración en tiempo real (RK4 3/8 y RK4 Clásico) con conservación de energía.
 * - Rejilla 2x2 de gráficas en vivo (Ángulos, Retrato de fase con cortes, Energía/Log, Panel contextual).
 * - Retratos de fase envueltos en (-pi, pi] cortando saltos (|Δθ| > pi) con moveTo.
 * - Mapa de caos adimensional FTLE con ejes, barra de color, contornos y retícula interactiva.
 * - Laboratorio de Monte Carlo con Web Worker (Blob), bootstrap percentil (B=2000), Wilson y comparativa A/B.
 * - Soporte táctil integral (iPad Retina) y modo presentación.
 */

(function () {
    'use strict';

    // =========================================================================
    // 1. CONFIGURACIÓN Y PALETAS DE COLOR
    // =========================================================================
    var TEMAS = {
        dark: {
            bg: '#0b0f17',
            panel: '#121826',
            panelSub: '#172033',
            border: '#1f293d',
            text: '#e6e9ef',
            textMuted: '#8892b0',
            grid: '#1f293d',
            th1: '#4cc9f0',        // Cian
            th2: '#f72585',        // Magenta
            energy: '#ffb703',     // Ámbar
            kinetic: '#2ec4b6',    // Verde Menta
            potential: '#9b5de5',  // Violeta
            traceSlow: [76, 201, 240],   // Cian
            traceMid: [255, 183, 3],     // Ámbar
            traceFast: [247, 37, 133],   // Magenta caliente
            contourM2: '#4cc9f0',
            contourM1: '#ffb703'
        },
        light: {
            bg: '#f4f6fa',
            panel: '#ffffff',
            panelSub: '#f8fafc',
            border: '#cbd5e1',
            text: '#0f172a',
            textMuted: '#475569',
            grid: '#e2e8f0',
            th1: '#0284c7',        // Cian oscuro
            th2: '#d946ef',        // Magenta
            energy: '#b45309',     // Ámbar
            kinetic: '#0d9488',    // Menta oscuro
            potential: '#7c3aed',  // Violeta
            traceSlow: [2, 132, 199],
            traceMid: [180, 83, 9],
            traceFast: [217, 70, 239],
            contourM2: '#0284c7',
            contourM1: '#111827'
        }
    };

    var temaActual = 'dark';
    function colorTema(nombre) {
        return TEMAS[temaActual][nombre];
    }

    // =========================================================================
    // 2. ESTADO GLOBAL DE LA SIMULACIÓN
    // =========================================================================
    var estado = {
        // Modo activo: 'normal', 'gemelo', 'lineal', 'nube', 'montecarlo'
        modo: 'normal',
        corriendo: true,
        velocidad: 1.0,
        mostrarEstela: true,
        algoritmo: 'rk4_38', // 'rk4_38' o 'rk4_clasico'
        subfaseRetrato: 'both', // 'both', 'p1', 'p2'

        // Parámetros físicos
        params: PenduloFisica.crearParametros(1.0, 1.0, 1.0, 1.0, 9.81),

        // Estado dinámico canónico x = [theta1, theta2, omega1, omega2]
        x: [Math.PI / 4, Math.PI / 4, 0.0, 0.0],
        x0: [Math.PI / 4, Math.PI / 4, 0.0, 0.0],
        t: 0.0,
        dtBase: 0.005,

        // Energía inicial de referencia
        e0: 0.0,
        deltaERel: 0.0,

        // Vueltas continuas y pulsos
        theta2Continuo: Math.PI / 4,
        vueltasM2: 0,
        vueltasM1: 0,
        theta1Continuo: Math.PI / 4,

        // Historial de estela térmica de la masa 2 (máx 1500)
        estelaM2: [],

        // Historial temporal para gráficas (últimos 20s)
        historialTiempo: [],
        ventanaHistorialS: 20.0,

        // Estado gemelo (Modo Gemelo)
        xGemelo: null,
        historialGemelo: [], // { t, d, logD }
        satGemelo: false,

        // Estado lineal (Modo Lineal)
        xLineal: null,
        historialLineal: [], // { t, th2Nonlin, th2Lin, diff }

        // Estado Nube (Modo Nube)
        nubePendulos: [],
        nubeFlips: [],

        // Arrastre interactivo
        arrastrandoMasa: 0, // 0: ninguna, 1: masa 1, 2: masa 2

        // Mapa FTLE
        ftleData: null,
        ftleWorker: null,
        ftleCalculando: false
    };

    // =========================================================================
    // 3. ELEMENTOS DOM Y CANVAS
    // =========================================================================
    var dom = {
        body: document.body,
        simCanvas: document.getElementById('sim-canvas'),
        canvasWrapper: document.getElementById('canvas-wrapper'),
        linearWarning: document.getElementById('linear-range-warning'),
        
        // HUD
        valTime: document.getElementById('val-time'),
        valEnergyTotal: document.getElementById('val-energy-total'),
        valDeltaE: document.getElementById('val-delta-e'),
        valFlipsM2: document.getElementById('val-flips-m2'),

        // Barra de Energía
        fillKinetic: document.getElementById('energy-fill-kinetic'),
        fillPotential: document.getElementById('energy-fill-potential'),
        markerM2: document.getElementById('marker-ec-m2'),
        markerM1: document.getElementById('marker-ec-m1'),
        markerAmbas: document.getElementById('marker-ec-ambas'),
        valEnergyK: document.getElementById('val-energy-k'),
        valEnergyV: document.getElementById('val-energy-v'),
        valEnergyE: document.getElementById('val-energy-e'),
        lblEcM2Val: document.getElementById('lbl-ec-m2-val'),
        lblEcM1Val: document.getElementById('lbl-ec-m1-val'),
        lblEcAmbasVal: document.getElementById('lbl-ec-ambas-val'),

        // Gráficas en vivo (2x2)
        chartAnglesCanvas: document.getElementById('chart-angles'),
        chartPhaseCanvas: document.getElementById('chart-phase'),
        chartEnergyCanvas: document.getElementById('chart-energy'),
        contextualTitle: document.getElementById('contextual-panel-title'),
        contextualMeta: document.getElementById('contextual-meta-text'),
        
        // Vistas contextuales
        ftleContainer: document.getElementById('ftle-map-container'),
        ftleCanvas: document.getElementById('ftle-canvas'),
        ftleTh1: document.getElementById('ftle-th1'),
        ftleTh2: document.getElementById('ftle-th2'),
        ftleVal: document.getElementById('ftle-val'),
        gemeloContainer: document.getElementById('gemelo-container'),
        gemeloCanvas: document.getElementById('gemelo-canvas'),
        linealContainer: document.getElementById('lineal-container'),
        linealCanvas: document.getElementById('lineal-canvas'),
        nubeContainer: document.getElementById('nube-container'),
        nubeCanvas: document.getElementById('nube-canvas'),

        // Botones y Transporte
        btnPlayPause: document.getElementById('btn-play-pause'),
        btnStep: document.getElementById('btn-step'),
        btnReset: document.getElementById('btn-reset'),
        selectSpeed: document.getElementById('select-speed'),
        checkTrace: document.getElementById('check-trace'),
        btnThemeToggle: document.getElementById('btn-theme-toggle'),
        themeIcon: document.getElementById('theme-icon'),
        themeText: document.getElementById('theme-text'),
        btnDrawerToggle: document.getElementById('btn-drawer-toggle'),
        btnDrawerClose: document.getElementById('btn-drawer-close'),
        drawerPanel: document.getElementById('drawer-panel'),
        btnPresentation: document.getElementById('btn-presentation'),

        // Sliders de parámetros
        sliderG: document.getElementById('slider-g'),
        lblG: document.getElementById('lbl-g'),
        sliderL1: document.getElementById('slider-l1'),
        lblL1: document.getElementById('lbl-l1'),
        sliderL2: document.getElementById('slider-l2'),
        lblL2: document.getElementById('lbl-l2'),
        sliderM1: document.getElementById('slider-m1'),
        lblM1: document.getElementById('lbl-m1'),
        sliderM2: document.getElementById('slider-m2'),
        lblM2: document.getElementById('lbl-m2'),
        lblMu: document.getElementById('lbl-mu'),
        lblLambda: document.getElementById('lbl-lambda'),

        // Sliders CI
        sliderTh1: document.getElementById('slider-th1'),
        lblTh1: document.getElementById('lbl-th1'),
        sliderTh2: document.getElementById('slider-th2'),
        lblTh2: document.getElementById('lbl-th2'),
        sliderW1: document.getElementById('slider-w1'),
        lblW1: document.getElementById('lbl-w1'),
        sliderW2: document.getElementById('slider-w2'),
        lblW2: document.getElementById('lbl-w2'),

        // Integradores
        btnIntegRK4: document.getElementById('btn-integ-rk4'),
        btnIntegClasico: document.getElementById('btn-integ-clasico'),

        // Monte Carlo
        montecarloPanel: document.getElementById('montecarlo-panel'),
        mcBaseLabel: document.getElementById('mc-base-label'),
        btnMcUseCurrent: document.getElementById('btn-mc-use-current'),
        mcSliderN: document.getElementById('mc-slider-n'),
        mcLblN: document.getElementById('mc-lbl-n'),
        mcSliderSigma: document.getElementById('mc-slider-sigma'),
        mcLblSigma: document.getElementById('mc-lbl-sigma'),
        mcSliderRho: document.getElementById('mc-slider-rho'),
        mcLblRho: document.getElementById('mc-lbl-rho'),
        mcSliderT: document.getElementById('mc-slider-t'),
        mcLblT: document.getElementById('mc-lbl-t'),
        mcSeedInput: document.getElementById('mc-seed-input'),
        btnMcRun: document.getElementById('btn-mc-run'),
        btnMcCompareLow: document.getElementById('btn-mc-compare-low'),
        btnMcReplayCloud: document.getElementById('btn-mc-replay-cloud'),
        mcProgressFill: document.getElementById('mc-progress-fill'),
        mcProgressText: document.getElementById('mc-progress-text'),
        mcValMean: document.getElementById('mc-val-mean'),
        mcCiMean: document.getElementById('mc-ci-mean'),
        mcValStd: document.getElementById('mc-val-std'),
        mcCiStd: document.getElementById('mc-ci-std'),
        mcValPflip: document.getElementById('mc-val-pflip'),
        mcCiWilson: document.getElementById('mc-ci-wilson'),
        mcHistogramCanvas: document.getElementById('mc-histogram-canvas')
    };

    // =========================================================================
    // 4. INICIALIZACIÓN DE CONTEXTOS CANVAS (CON DEVICE PIXEL RATIO)
    // =========================================================================
    function ajustarCanvas(canvas) {
        if (!canvas) return null;
        var dpr = window.devicePixelRatio || 1;
        var rect = canvas.getBoundingClientRect();
        var ancho = Math.floor(rect.width * dpr);
        var alto = Math.floor(rect.height * dpr);

        if (canvas.width !== ancho || canvas.height !== alto) {
            canvas.width = ancho;
            canvas.height = alto;
        }
        var ctx = canvas.getContext('2d');
        ctx.resetTransform();
        ctx.scale(dpr, dpr);
        return { ctx: ctx, width: rect.width, height: rect.height, dpr: dpr };
    }

    // =========================================================================
    // 5. GESTIÓN DEL TEMA (OSCURO / CLARO)
    // =========================================================================
    function alternarTema() {
        if (temaActual === 'dark') {
            temaActual = 'light';
            dom.body.classList.remove('theme-dark');
            dom.body.classList.add('theme-light');
            dom.themeIcon.textContent = '☀️';
            dom.themeText.textContent = 'Claro';
        } else {
            temaActual = 'dark';
            dom.body.classList.remove('theme-light');
            dom.body.classList.add('theme-dark');
            dom.themeIcon.textContent = '🌙';
            dom.themeText.textContent = 'Oscuro';
        }
        actualizarUmbralesEnergiaUI();
    }

    // =========================================================================
    // 6. MOTOR DE FÍSICA Y REINICIO
    // =========================================================================
    function calcularPasoFisica(xActual, dt, params) {
        if (estado.algoritmo === 'rk4_clasico') {
            return PenduloFisica.pasoRK4Clasico(xActual, dt, params);
        }
        return PenduloFisica.pasoRK4(xActual, dt, params);
    }

    function reiniciarSimulacion() {
        estado.t = 0.0;
        estado.x = [estado.x0[0], estado.x0[1], estado.x0[2], estado.x0[3]];
        estado.theta1Continuo = estado.x0[0];
        estado.theta2Continuo = estado.x0[1];
        estado.vueltasM1 = 0;
        estado.vueltasM2 = 0;
        estado.estelaM2 = [];
        estado.historialTiempo = [];

        estado.e0 = PenduloFisica.energiaTotal(estado.x, estado.params);
        estado.deltaERel = 0.0;

        // Modo Gemelo
        if (estado.modo === 'gemelo') {
            var pert = 1e-7;
            estado.xGemelo = [estado.x0[0] + pert, estado.x0[1], estado.x0[2], estado.x0[3]];
            estado.historialGemelo = [];
            estado.satGemelo = false;
        }

        // Modo Lineal
        if (estado.modo === 'lineal') {
            estado.xLineal = [estado.x0[0], estado.x0[1], estado.x0[2], estado.x0[3]];
            estado.historialLineal = [];
        }

        // Modo Nube
        if (estado.modo === 'nube') {
            inicializarNube(200, 0.01);
        }

        actualizarUmbralesEnergiaUI();
        actualizarLabelsParametros();
    }

    function inicializarNube(n, sigma) {
        var rng = PenduloFisica.mulberry32(20261119);
        estado.nubePendulos = [];
        estado.nubeFlips = new Int32Array(n);

        for (var i = 0; i < n; i++) {
            // Box-Muller simple
            var u1 = Math.max(1e-12, rng());
            var u2 = rng();
            var z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
            var z2 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);

            var th1 = estado.x0[0] + sigma * z1;
            var th2 = estado.x0[1] + sigma * z2;
            estado.nubePendulos.push({
                x: [th1, th2, 0.0, 0.0],
                th2Cont: th2,
                flips: 0
            });
        }
    }

    function actualizarUmbralesEnergiaUI() {
        var u = PenduloFisica.umbralesEnergia(estado.params);
        var ecAmbas = u.Ec_ambas;

        dom.lblEcM2Val.textContent = u.Ec_m2.toFixed(2);
        dom.lblEcM1Val.textContent = u.Ec_m1.toFixed(2);
        dom.lblEcAmbasVal.textContent = ecAmbas.toFixed(2);

        // Posiciones relativas en la barra de energía
        var pM2 = Math.min(100, Math.max(0, (u.Ec_m2 / ecAmbas) * 100));
        var pM1 = Math.min(100, Math.max(0, (u.Ec_m1 / ecAmbas) * 100));
        var pAmbas = 100;

        dom.markerM2.style.left = pM2 + '%';
        dom.markerM1.style.left = pM1 + '%';
        dom.markerAmbas.style.left = pAmbas + '%';
    }

    function actualizarLabelsParametros() {
        dom.lblG.textContent = estado.params.g.toFixed(2) + ' m/s²';
        dom.lblL1.textContent = estado.params.l1.toFixed(2) + ' m';
        dom.lblL2.textContent = estado.params.l2.toFixed(2) + ' m';
        dom.lblM1.textContent = estado.params.m1.toFixed(2) + ' kg';
        dom.lblM2.textContent = estado.params.m2.toFixed(2) + ' kg';

        var mu = estado.params.m2 / estado.params.m1;
        var lambda = estado.params.l2 / estado.params.l1;
        dom.lblMu.textContent = mu.toFixed(2);
        dom.lblLambda.textContent = lambda.toFixed(2);

        dom.lblTh1.textContent = Math.round(estado.x0[0] * 180 / Math.PI) + '°';
        dom.lblTh2.textContent = Math.round(estado.x0[1] * 180 / Math.PI) + '°';
        dom.lblW1.textContent = estado.x0[2].toFixed(1) + ' rad/s';
        dom.lblW2.textContent = estado.x0[3].toFixed(1) + ' rad/s';

        // Actualizar valores de sliders
        dom.sliderG.value = estado.params.g;
        dom.sliderL1.value = estado.params.l1;
        dom.sliderL2.value = estado.params.l2;
        dom.sliderM1.value = estado.params.m1;
        dom.sliderM2.value = estado.params.m2;
        dom.sliderTh1.value = Math.round(estado.x0[0] * 180 / Math.PI);
        dom.sliderTh2.value = Math.round(estado.x0[1] * 180 / Math.PI);
        dom.sliderW1.value = estado.x0[2];
        dom.sliderW2.value = estado.x0[3];
    }

    // =========================================================================
    // 7. PASO DE SIMULACIÓN Y ACTUALIZACIÓN DE ESTADO
    // =========================================================================
    function pasoSimulacion(dtPaso) {
        var xPrev = estado.x;
        var xNext = calcularPasoFisica(xPrev, dtPaso, estado.params);

        // Detección de vueltas mediante ángulo acumulado continuo
        var dTh2 = xNext[1] - xPrev[1];
        // Desenrollar diferencias angulares mayores a pi
        dTh2 = ((dTh2 + Math.PI) % (2.0 * Math.PI)) - Math.PI;
        estado.theta2Continuo += dTh2;
        var nuevasVueltasM2 = Math.abs(PenduloFisica.indiceVuelta(estado.theta2Continuo));

        if (nuevasVueltasM2 !== estado.vueltasM2) {
            estado.vueltasM2 = nuevasVueltasM2;
            animarPulsoVueltas();
        }

        var dTh1 = xNext[0] - xPrev[0];
        dTh1 = ((dTh1 + Math.PI) % (2.0 * Math.PI)) - Math.PI;
        estado.theta1Continuo += dTh1;
        estado.vueltasM1 = Math.abs(PenduloFisica.indiceVuelta(estado.theta1Continuo));

        estado.x = xNext;
        estado.t += dtPaso;

        // Cálculos de energía
        var eK = PenduloFisica.energiaCinetica(estado.x, estado.params);
        var eV = PenduloFisica.energiaPotencial(estado.x, estado.params);
        var eTotal = eK + eV;
        if (Math.abs(estado.e0) > 1e-12) {
            estado.deltaERel = Math.abs(eTotal - estado.e0) / Math.abs(estado.e0);
        } else {
            estado.deltaERel = Math.abs(eTotal - estado.e0);
        }

        // Registrar estela de la masa 2
        var pos = PenduloFisica.posicionesCartesianas(estado.x, estado.params);
        // Rapidez puntual de la masa 2
        var vx2 = estado.params.l1 * Math.cos(estado.x[0]) * estado.x[2] +
                  estado.params.l2 * Math.cos(estado.x[1]) * estado.x[3];
        var vy2 = estado.params.l1 * Math.sin(estado.x[0]) * estado.x[2] +
                  estado.params.l2 * Math.sin(estado.x[1]) * estado.x[3];
        var rapidez = Math.sqrt(vx2 * vx2 + vy2 * vy2);

        estado.estelaM2.push({
            x: pos.x2,
            y: pos.y2,
            v: rapidez,
            t: estado.t
        });
        if (estado.estelaM2.length > 1500) {
            estado.estelaM2.shift();
        }

        // Historial para gráficas en vivo (últimos 20 s)
        estado.historialTiempo.push({
            t: estado.t,
            th1: estado.x[0],
            th2: estado.x[1],
            w1: estado.x[2],
            w2: estado.x[3],
            T: eK,
            V: eV,
            E: eTotal,
            errLog: Math.max(-14, Math.log10(Math.max(1e-14, estado.deltaERel)))
        });

        // Poda de ventana deslizante
        var tMin = estado.t - estado.ventanaHistorialS;
        while (estado.historialTiempo.length > 1 && estado.historialTiempo[0].t < tMin) {
            estado.historialTiempo.shift();
        }

        // Actualizaciones de modos especiales
        if (estado.modo === 'gemelo' && estado.xGemelo) {
            estado.xGemelo = calcularPasoFisica(estado.xGemelo, dtPaso, estado.params);
            var d1 = estado.x[0] - estado.xGemelo[0];
            var d2 = estado.x[1] - estado.xGemelo[1];
            var dw1 = estado.x[2] - estado.xGemelo[2];
            var dw2 = estado.x[3] - estado.xGemelo[3];
            var sep = Math.sqrt(d1 * d1 + d2 * d2 + dw1 * dw1 + dw2 * dw2);
            var logSep = Math.log10(Math.max(1e-15, sep));
            estado.historialGemelo.push({ t: estado.t, d: sep, logD: logSep });
            if (sep >= 0.5) {
                estado.satGemelo = true;
            }
            if (estado.historialGemelo.length > 2000) estado.historialGemelo.shift();
        }

        if (estado.modo === 'lineal' && estado.xLineal) {
            var linInfo = PenduloFisica.linealizacion(estado.params);
            estado.xLineal = PenduloFisica.pasoLinealRK4(estado.xLineal, dtPaso, linInfo);
            var diffTh2 = Math.abs(estado.x[1] - estado.xLineal[1]);
            estado.historialLineal.push({
                t: estado.t,
                th2Nonlin: estado.x[1],
                th2Lin: estado.xLineal[1],
                diff: diffTh2
            });
            if (estado.historialLineal.length > 2000) estado.historialLineal.shift();
        }

        if (estado.modo === 'nube') {
            for (var k = 0; k < estado.nubePendulos.length; k++) {
                var p = estado.nubePendulos[k];
                var prevTh2 = p.x[1];
                p.x = calcularPasoFisica(p.x, dtPaso, estado.params);
                var delta = ((p.x[1] - prevTh2 + Math.PI) % (2.0 * Math.PI)) - Math.PI;
                p.th2Cont += delta;
                p.flips = Math.abs(PenduloFisica.indiceVuelta(p.th2Cont));
                estado.nubeFlips[k] = p.flips;
            }
        }
    }

    function animarPulsoVueltas() {
        dom.valFlipsM2.classList.remove('pulse');
        void dom.valFlipsM2.offsetWidth; // Trigger reflow
        dom.valFlipsM2.classList.add('pulse');
    }

    // =========================================================================
    // 8. RENDERIZADO DEL ESCENARIO DEL PÉNDULO (#sim-canvas)
    // =========================================================================
    function dibujarEscenario() {
        var info = ajustarCanvas(dom.simCanvas);
        if (!info) return;
        var ctx = info.ctx;
        var w = info.width;
        var h = info.height;

        ctx.clearRect(0, 0, w, h);

        // Centro del pivote y escala física
        var cx = w * 0.5;
        var cy = h * 0.42;
        var lTotal = estado.params.l1 + estado.params.l2;
        var maxRadioPx = Math.min(w, h) * 0.40;
        var escala = maxRadioPx / lTotal;

        var l1Px = estado.params.l1 * escala;
        var l2Px = estado.params.l2 * escala;

        var th1 = estado.x[0];
        var th2 = estado.x[1];

        var x1 = cx + l1Px * Math.sin(th1);
        var y1 = cy + l1Px * Math.cos(th1);
        var x2 = x1 + l2Px * Math.sin(th2);
        var y2 = y1 + l2Px * Math.cos(th2);

        // 1. Estela térmica de la masa 2
        if (estado.mostrarEstela && estado.estelaM2.length > 1) {
            ctx.save();
            if (temaActual === 'dark') {
                ctx.globalCompositeOperation = 'lighter';
            }
            var nPuntos = estado.estelaM2.length;
            for (var i = 1; i < nPuntos; i++) {
                var pA = estado.estelaM2[i - 1];
                var pB = estado.estelaM2[i];

                var pxA = cx + pA.x * escala;
                var pyA = cy - pA.y * escala; // y física hacia arriba
                var pxB = cx + pB.x * escala;
                var pyB = cy - pB.y * escala;

                var fraccEdad = i / nPuntos;
                var vNorm = Math.min(1.0, pB.v / 12.0); // Normalizado hasta 12 m/s

                // Interpolación térmica: frío (cian) -> templado (ámbar) -> caliente (magenta)
                var r, g, b;
                if (vNorm < 0.5) {
                    var f = vNorm * 2.0;
                    r = Math.round(TEMAS[temaActual].traceSlow[0] * (1 - f) + TEMAS[temaActual].traceMid[0] * f);
                    g = Math.round(TEMAS[temaActual].traceSlow[1] * (1 - f) + TEMAS[temaActual].traceMid[1] * f);
                    b = Math.round(TEMAS[temaActual].traceSlow[2] * (1 - f) + TEMAS[temaActual].traceMid[2] * f);
                } else {
                    var f2 = (vNorm - 0.5) * 2.0;
                    r = Math.round(TEMAS[temaActual].traceMid[0] * (1 - f2) + TEMAS[temaActual].traceFast[0] * f2);
                    g = Math.round(TEMAS[temaActual].traceMid[1] * (1 - f2) + TEMAS[temaActual].traceFast[1] * f2);
                    b = Math.round(TEMAS[temaActual].traceMid[2] * (1 - f2) + TEMAS[temaActual].traceFast[2] * f2);
                }

                var alpha = fraccEdad * (temaActual === 'dark' ? 0.75 : 0.60);
                ctx.strokeStyle = 'rgba(' + r + ',' + g + ',' + b + ',' + alpha + ')';
                ctx.lineWidth = 1.2 + 2.0 * fraccEdad;
                ctx.beginPath();
                ctx.moveTo(pxA, pyA);
                ctx.lineTo(pxB, pyB);
                ctx.stroke();
            }
            ctx.restore();
        }

        // 2. Líneas verticales punteadas de referencia
        ctx.save();
        ctx.setLineDash([4, 4]);
        ctx.strokeStyle = colorTema('textMuted');
        ctx.lineWidth = 1;
        ctx.globalAlpha = 0.5;

        // Vertical en pivote
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(cx, cy + l1Px + 20);
        ctx.stroke();

        // Vertical en masa 1
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x1, y1 + l2Px + 20);
        ctx.stroke();
        ctx.restore();

        // 3. Cuñas de ángulo sombreadas (~45 px de radio)
        var radioCuña = Math.min(48, Math.min(l1Px, l2Px) * 0.7);

        // Cuña θ1
        ctx.save();
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        // Ángulo de la vertical abajo es +PI/2
        var angVert = Math.PI * 0.5;
        var angVarilla1 = angVert + th1;
        ctx.arc(cx, cy, radioCuña, angVert, angVarilla1, th1 < 0);
        ctx.closePath();
        ctx.fillStyle = temaActual === 'dark' ? 'rgba(76, 201, 240, 0.22)' : 'rgba(2, 132, 199, 0.22)';
        ctx.fill();
        ctx.strokeStyle = colorTema('th1');
        ctx.lineWidth = 1.5;
        ctx.stroke();

        // Cuña θ2
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        var angVarilla2 = angVert + th2;
        ctx.arc(x1, y1, radioCuña, angVert, angVarilla2, th2 < 0);
        ctx.closePath();
        ctx.fillStyle = temaActual === 'dark' ? 'rgba(247, 37, 133, 0.22)' : 'rgba(217, 70, 239, 0.22)';
        ctx.fill();
        ctx.strokeStyle = colorTema('th2');
        ctx.lineWidth = 1.5;
        ctx.stroke();
        ctx.restore();

        // 4. Modo Lineal: Péndulo Fantasma
        if (estado.modo === 'lineal' && estado.xLineal) {
            var th1L = estado.xLineal[0];
            var th2L = estado.xLineal[1];
            var x1L = cx + l1Px * Math.sin(th1L);
            var y1L = cy + l1Px * Math.cos(th1L);
            var x2L = x1L + l2Px * Math.sin(th2L);
            var y2L = y1L + l2Px * Math.cos(th2L);

            ctx.save();
            ctx.setLineDash([5, 5]);
            ctx.strokeStyle = colorTema('textMuted');
            ctx.lineWidth = 2;
            ctx.globalAlpha = 0.55;
            ctx.beginPath();
            ctx.moveTo(cx, cy);
            ctx.lineTo(x1L, y1L);
            ctx.lineTo(x2L, y2L);
            ctx.stroke();

            ctx.fillStyle = colorTema('textMuted');
            ctx.beginPath();
            ctx.arc(x1L, y1L, 6, 0, 2 * Math.PI);
            ctx.arc(x2L, y2L, 6, 0, 2 * Math.PI);
            ctx.fill();
            ctx.restore();

            // Advertencia de rango lineal (|θ| > 0.3 rad)
            if (Math.abs(th1) > 0.3 || Math.abs(th2) > 0.3) {
                dom.linearWarning.style.display = 'block';
            } else {
                dom.linearWarning.style.display = 'none';
            }
        } else {
            dom.linearWarning.style.display = 'none';
        }

        // 5. Modo Gemelo: Péndulo Gemelo Faint
        if (estado.modo === 'gemelo' && estado.xGemelo) {
            var th1G = estado.xGemelo[0];
            var th2G = estado.xGemelo[1];
            var x1G = cx + l1Px * Math.sin(th1G);
            var y1G = cy + l1Px * Math.cos(th1G);
            var x2G = x1G + l2Px * Math.sin(th2G);
            var y2G = y1G + l2Px * Math.cos(th2G);

            ctx.save();
            ctx.strokeStyle = colorTema('energy');
            ctx.lineWidth = 1.5;
            ctx.globalAlpha = 0.6;
            ctx.beginPath();
            ctx.moveTo(cx, cy);
            ctx.lineTo(x1G, y1G);
            ctx.lineTo(x2G, y2G);
            ctx.stroke();

            ctx.fillStyle = colorTema('energy');
            ctx.beginPath();
            ctx.arc(x2G, y2G, 6, 0, 2 * Math.PI);
            ctx.fill();
            ctx.restore();
        }

        // 6. Modo Nube: Enjambre de masas 2
        if (estado.modo === 'nube' && estado.nubePendulos.length > 0) {
            ctx.save();
            if (temaActual === 'dark') ctx.globalCompositeOperation = 'lighter';
            for (var k = 0; k < estado.nubePendulos.length; k++) {
                var p = estado.nubePendulos[k];
                var px1 = cx + l1Px * Math.sin(p.x[0]);
                var py1 = cy + l1Px * Math.cos(p.x[0]);
                var px2 = px1 + l2Px * Math.sin(p.x[1]);
                var py2 = py1 + l2Px * Math.cos(p.x[1]);

                var colDot = p.flips > 0 ? colorTema('th2') : colorTema('th1');
                ctx.fillStyle = colDot;
                ctx.globalAlpha = 0.45;
                ctx.beginPath();
                ctx.arc(px2, py2, 3, 0, 2 * Math.PI);
                ctx.fill();
            }
            ctx.restore();
        }

        // 7. Varillas Físicas Principales
        ctx.save();
        ctx.lineCap = 'round';
        ctx.lineJoin = 'round';

        // Varilla 1
        ctx.strokeStyle = colorTema('th1');
        ctx.lineWidth = 4.5;
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(x1, y1);
        ctx.stroke();

        // Varilla 2
        ctx.strokeStyle = colorTema('th2');
        ctx.lineWidth = 4.0;
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();

        // Pivote central
        ctx.fillStyle = colorTema('text');
        ctx.beginPath();
        ctx.arc(cx, cy, 5.5, 0, 2 * Math.PI);
        ctx.fill();
        ctx.restore();

        // 8. Masas con radio cbrt(m), halo y gradiente
        var rBase = 15;
        var r1 = Math.max(8, rBase * Math.cbrt(estado.params.m1));
        var r2 = Math.max(8, rBase * Math.cbrt(estado.params.m2));

        // Dibujar Masa 1
        dibujarMasaConHalo(ctx, x1, y1, r1, colorTema('th1'), 'm₁: ' + estado.params.m1.toFixed(1) + ' kg');

        // Dibujar Masa 2
        dibujarMasaConHalo(ctx, x2, y2, r2, colorTema('th2'), 'm₂: ' + estado.params.m2.toFixed(1) + ' kg');
    }

    function dibujarMasaConHalo(ctx, x, y, radio, colorHex, etiqueta) {
        ctx.save();

        // Halo difuso suave
        var gradHalo = ctx.createRadialGradient(x, y, radio * 0.5, x, y, radio * 1.8);
        gradHalo.addColorStop(0, colorHex);
        gradHalo.addColorStop(1, 'transparent');
        ctx.fillStyle = gradHalo;
        ctx.globalAlpha = 0.35;
        ctx.beginPath();
        ctx.arc(x, y, radio * 1.8, 0, 2 * Math.PI);
        ctx.fill();

        // Disco principal con gradiente esférico
        ctx.globalAlpha = 1.0;
        var gradBob = ctx.createRadialGradient(x - radio * 0.3, y - radio * 0.3, radio * 0.1, x, y, radio);
        gradBob.addColorStop(0, '#ffffff');
        gradBob.addColorStop(0.3, colorHex);
        gradBob.addColorStop(1, '#000000');
        ctx.fillStyle = gradBob;
        ctx.beginPath();
        ctx.arc(x, y, radio, 0, 2 * Math.PI);
        ctx.fill();

        ctx.strokeStyle = colorHex;
        ctx.lineWidth = 1.5;
        ctx.stroke();

        // Etiqueta al lado del disco con alto contraste
        ctx.font = '600 11px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
        var textW = ctx.measureText(etiqueta).width;
        var textX = x + radio + 8;
        var textY = y + 4;

        // Pastilla de fondo para legibilidad
        ctx.fillStyle = temaActual === 'dark' ? 'rgba(18, 24, 38, 0.88)' : 'rgba(255, 255, 255, 0.88)';
        ctx.fillRect(textX - 4, textY - 11, textW + 8, 15);
        ctx.strokeStyle = colorTema('border');
        ctx.lineWidth = 1;
        ctx.strokeRect(textX - 4, textY - 11, textW + 8, 15);

        ctx.fillStyle = colorTema('text');
        ctx.fillText(etiqueta, textX, textY);

        ctx.restore();
    }

    // =========================================================================
    // 9. RENDERIZADO DE LAS 4 GRÁFICAS EN VIVO (REJILLA 2x2)
    // =========================================================================

    // Gráfica 1: θ1(t) y θ2(t) en ventana rodante de 20s
    function dibujarGraficaAngulos() {
        var info = ajustarCanvas(dom.chartAnglesCanvas);
        if (!info) return;
        var ctx = info.ctx, w = info.width, h = info.height;

        ctx.clearRect(0, 0, w, h);
        dibujarMarcoYRejilla(ctx, w, h, 'θ (rad)');

        var hist = estado.historialTiempo;
        if (hist.length < 2) return;

        var tMax = estado.t;
        var tMin = Math.max(0, tMax - estado.ventanaHistorialS);
        var dtSpan = Math.max(0.1, tMax - tMin);

        var escalaY = (h - 24) / (2.0 * Math.PI); // de -pi a +pi
        var yCentro = h * 0.5;

        // Curva θ1 (Cian)
        ctx.save();
        ctx.strokeStyle = colorTema('th1');
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        for (var i = 0; i < hist.length; i++) {
            var pt = hist[i];
            var xPx = 36 + ((pt.t - tMin) / dtSpan) * (w - 44);
            // Envolver a [-pi, pi]
            var th1Wrap = ((pt.th1 + Math.PI) % (2 * Math.PI)) - Math.PI;
            var yPx = yCentro - th1Wrap * escalaY;
            if (i === 0) ctx.moveTo(xPx, yPx);
            else ctx.lineTo(xPx, yPx);
        }
        ctx.stroke();

        // Curva θ2 (Magenta)
        ctx.strokeStyle = colorTema('th2');
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        for (var j = 0; j < hist.length; j++) {
            var pt2 = hist[j];
            var xPx2 = 36 + ((pt2.t - tMin) / dtSpan) * (w - 44);
            var th2Wrap = ((pt2.th2 + Math.PI) % (2 * Math.PI)) - Math.PI;
            var yPx2 = yCentro - th2Wrap * escalaY;
            if (j === 0) ctx.moveTo(xPx2, yPx2);
            else ctx.lineTo(xPx2, yPx2);
        }
        ctx.stroke();
        ctx.restore();
    }

    // Gráfica 2: Retratos de Fase (θ, ω) con envoltura y cortes en |Δθ| > pi
    function dibujarGraficaFase() {
        var info = ajustarCanvas(dom.chartPhaseCanvas);
        if (!info) return;
        var ctx = info.ctx, w = info.width, h = info.height;

        ctx.clearRect(0, 0, w, h);
        dibujarMarcoFase(ctx, w, h);

        var hist = estado.historialTiempo;
        if (hist.length < 2) return;

        var sub = estado.subfaseRetrato;
        var cx = 34 + (w - 42) * 0.5;
        var cy = h * 0.5;
        var escalaX = (w - 46) * 0.5 / Math.PI; // -pi a +pi
        var maxOmega = 12.0;
        var escalaY = (h - 20) * 0.5 / maxOmega;

        // Péndulo 1
        if (sub === 'both' || sub === 'p1') {
            ctx.save();
            ctx.strokeStyle = colorTema('th1');
            ctx.lineWidth = 1.4;
            ctx.beginPath();
            var prevTh1 = null;

            for (var i = 0; i < hist.length; i++) {
                var pt = hist[i];
                var thWrap = ((pt.th1 + Math.PI) % (2 * Math.PI)) - Math.PI;
                var px = cx + thWrap * escalaX;
                var py = cy - Math.max(-maxOmega, Math.min(maxOmega, pt.w1)) * escalaY;

                if (prevTh1 === null || Math.abs(thWrap - prevTh1) > Math.PI) {
                    ctx.moveTo(px, py); // CORTE VITAL: Evita raya horizontal
                } else {
                    ctx.lineTo(px, py);
                }
                prevTh1 = thWrap;
            }
            ctx.stroke();
            ctx.restore();
        }

        // Péndulo 2
        if (sub === 'both' || sub === 'p2') {
            ctx.save();
            ctx.strokeStyle = colorTema('th2');
            ctx.lineWidth = 1.4;
            ctx.beginPath();
            var prevTh2 = null;

            for (var j = 0; j < hist.length; j++) {
                var pt2 = hist[j];
                var th2Wrap = ((pt2.th2 + Math.PI) % (2 * Math.PI)) - Math.PI;
                var px2 = cx + th2Wrap * escalaX;
                var py2 = cy - Math.max(-maxOmega, Math.min(maxOmega, pt2.w2)) * escalaY;

                if (prevTh2 === null || Math.abs(th2Wrap - prevTh2) > Math.PI) {
                    ctx.moveTo(px2, py2); // CORTE VITAL
                } else {
                    ctx.lineTo(px2, py2);
                }
                prevTh2 = th2Wrap;
            }
            ctx.stroke();
            ctx.restore();
        }
    }

    // Gráfica 3: Energía apilada (T menta, V violeta, E ámbar) + Tira Logarítmica
    function dibujarGraficaEnergia() {
        var info = ajustarCanvas(dom.chartEnergyCanvas);
        if (!info) return;
        var ctx = info.ctx, w = info.width, h = info.height;

        ctx.clearRect(0, 0, w, h);

        var hLog = 28;
        var hEnergia = h - hLog - 6;

        var hist = estado.historialTiempo;
        var tMax = estado.t;
        var tMin = Math.max(0, tMax - estado.ventanaHistorialS);
        var dtSpan = Math.max(0.1, tMax - tMin);

        // Marco área de energía
        ctx.save();
        ctx.strokeStyle = colorTema('border');
        ctx.lineWidth = 1;
        ctx.strokeRect(34, 4, w - 42, hEnergia - 6);

        var maxE = Math.max(10.0, estado.e0 * 1.25);
        var escalaY = (hEnergia - 12) / maxE;
        var yBase = hEnergia - 2;

        if (hist.length >= 2) {
            // Área de T (cinética - menta)
            ctx.beginPath();
            ctx.moveTo(34, yBase);
            for (var i = 0; i < hist.length; i++) {
                var pt = hist[i];
                var px = 34 + ((pt.t - tMin) / dtSpan) * (w - 42);
                var pyT = yBase - pt.T * escalaY;
                ctx.lineTo(px, pyT);
            }
            ctx.lineTo(34 + (w - 42), yBase);
            ctx.closePath();
            ctx.fillStyle = colorTema('kinetic');
            ctx.globalAlpha = 0.65;
            ctx.fill();

            // Área de V apilada encima de T (violeta)
            ctx.beginPath();
            ctx.globalAlpha = 0.55;
            ctx.fillStyle = colorTema('potential');
            for (var j = 0; j < hist.length; j++) {
                var pt2 = hist[j];
                var px2 = 34 + ((pt2.t - tMin) / dtSpan) * (w - 42);
                var pyT2 = yBase - pt2.T * escalaY;
                var pyTotal = yBase - (pt2.T + pt2.V) * escalaY;
                if (j === 0) ctx.moveTo(px2, pyT2);
                else ctx.lineTo(px2, pyTotal);
            }
            for (var k = hist.length - 1; k >= 0; k--) {
                var ptK = hist[k];
                var pxK = 34 + ((ptK.t - tMin) / dtSpan) * (w - 42);
                var pyTK = yBase - ptK.T * escalaY;
                ctx.lineTo(pxK, pyTK);
            }
            ctx.closePath();
            ctx.fill();

            // Línea de energía total E (ámbar)
            ctx.globalAlpha = 1.0;
            ctx.strokeStyle = colorTema('energy');
            ctx.lineWidth = 2.0;
            ctx.beginPath();
            for (var m = 0; m < hist.length; m++) {
                var ptM = hist[m];
                var pxM = 34 + ((ptM.t - tMin) / dtSpan) * (w - 42);
                var pyE = yBase - ptM.E * escalaY;
                if (m === 0) ctx.moveTo(pxM, pyE);
                else ctx.lineTo(pxM, pyE);
            }
            ctx.stroke();
        }

        // Tira logarítmica inferior |ΔE|/E0
        var yTopLog = h - hLog;
        ctx.fillStyle = colorTema('panelSub');
        ctx.fillRect(34, yTopLog, w - 42, hLog - 2);
        ctx.strokeStyle = colorTema('border');
        ctx.strokeRect(34, yTopLog, w - 42, hLog - 2);

        ctx.font = '500 8.5px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
        ctx.fillStyle = colorTema('textMuted');
        ctx.fillText('log₁₀|ΔE|', 2, yTopLog + 12);

        if (hist.length >= 2) {
            ctx.strokeStyle = colorTema('th1');
            ctx.lineWidth = 1.2;
            ctx.beginPath();
            for (var n = 0; n < hist.length; n++) {
                var ptN = hist[n];
                var pxN = 34 + ((ptN.t - tMin) / dtSpan) * (w - 42);
                // Rango logarítmico de -14 a 0
                var logNorm = (ptN.errLog - (-14)) / 14.0;
                var pyLog = (yTopLog + hLog - 3) - Math.max(0, Math.min(1, logNorm)) * (hLog - 4);
                if (n === 0) ctx.moveTo(pxN, pyLog);
                else ctx.lineTo(pxN, pyLog);
            }
            ctx.stroke();
        }
        ctx.restore();
    }

    // Gráfica 4: Panel Contextual dinámico
    function dibujarGraficaContextual() {
        if (estado.modo === 'gemelo') {
            dibujarGraficaGemelo();
        } else if (estado.modo === 'lineal') {
            dibujarGraficaLineal();
        } else if (estado.modo === 'nube') {
            dibujarGraficaNube();
        } else {
            // Modo Normal: el mapa FTLE ya se dibuja en su propio canvas
            // Solo redibujamos la retícula si hace falta
            dibujarReticulaFtle();
        }
    }

    function dibujarGraficaGemelo() {
        var info = ajustarCanvas(dom.gemeloCanvas);
        if (!info) return;
        var ctx = info.ctx, w = info.width, h = info.height;
        ctx.clearRect(0, 0, w, h);

        dibujarMarcoYRejilla(ctx, w, h, 'log₁₀ d(t)');
        var hist = estado.historialGemelo;
        if (hist.length < 2) return;

        var tMin = hist[0].t;
        var tMax = hist[hist.length - 1].t;
        var dtSpan = Math.max(0.1, tMax - tMin);

        var logMin = -15.0;
        var logMax = 0.5;
        var escalaY = (h - 26) / (logMax - logMin);

        ctx.save();
        ctx.strokeStyle = colorTema('energy');
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        for (var i = 0; i < hist.length; i++) {
            var px = 34 + ((hist[i].t - tMin) / dtSpan) * (w - 42);
            var py = (h - 16) - (hist[i].logD - logMin) * escalaY;
            if (i === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();

        // Estimación de lambda mediante regresión lineal en la fase creciente (d < 0.2)
        var puntosAjuste = [];
        for (var k = 0; k < hist.length; k++) {
            if (hist[k].d < 0.2) {
                puntosAjuste.push(hist[k]);
            }
        }

        var lambdaEst = 0;
        if (puntosAjuste.length >= 10) {
            var sumT = 0, sumL = 0, sumTL = 0, sumTT = 0;
            var N = puntosAjuste.length;
            for (var m = 0; m < N; m++) {
                sumT += puntosAjuste[m].t;
                sumL += puntosAjuste[m].logD;
                sumTL += puntosAjuste[m].t * puntosAjuste[m].logD;
                sumTT += puntosAjuste[m].t * puntosAjuste[m].t;
            }
            var pendLog10 = (N * sumTL - sumT * sumL) / (N * sumTT - sumT * sumT);
            // d(t) ~ d0 * e^(lambda*t) -> log10(d) = pend*t -> lambda = pend * ln(10)
            lambdaEst = pendLog10 * Math.LN10;
        }

        var tCaract = Math.sqrt(estado.params.l1 / estado.params.g);
        var lambdaAdim = lambdaEst * tCaract;

        ctx.font = '600 11px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
        ctx.fillStyle = colorTema('text');
        var textoLambda = 'λ ≈ ' + Math.max(0, lambdaEst).toFixed(2) + ' s⁻¹ | λ·√(l₁/g) ≈ ' + Math.max(0, lambdaAdim).toFixed(2);
        if (estado.satGemelo) {
            textoLambda += ' (Saturado)';
        }
        ctx.fillText(textoLambda, 40, 20);
        ctx.restore();
    }

    function dibujarGraficaLineal() {
        var info = ajustarCanvas(dom.linealCanvas);
        if (!info) return;
        var ctx = info.ctx, w = info.width, h = info.height;
        ctx.clearRect(0, 0, w, h);

        dibujarMarcoYRejilla(ctx, w, h, 'θ₂ (rad)');
        var hist = estado.historialLineal;
        if (hist.length < 2) return;

        var tMax = estado.t;
        var tMin = Math.max(0, tMax - estado.ventanaHistorialS);
        var dtSpan = Math.max(0.1, tMax - tMin);

        var escalaY = (h - 26) / (2.0 * Math.PI);
        var cy = h * 0.5;

        // θ2 No lineal (Magenta)
        ctx.save();
        ctx.strokeStyle = colorTema('th2');
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        for (var i = 0; i < hist.length; i++) {
            var px = 34 + ((hist[i].t - tMin) / dtSpan) * (w - 42);
            var py = cy - hist[i].th2Nonlin * escalaY;
            if (i === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();

        // θ2 Lineal Fantasma (Punteado)
        ctx.setLineDash([4, 4]);
        ctx.strokeStyle = colorTema('textMuted');
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        for (var j = 0; j < hist.length; j++) {
            var px2 = 34 + ((hist[j].t - tMin) / dtSpan) * (w - 42);
            var py2 = cy - hist[j].th2Lin * escalaY;
            if (j === 0) ctx.moveTo(px2, py2);
            else ctx.lineTo(px2, py2);
        }
        ctx.stroke();

        ctx.font = '600 11px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
        ctx.fillStyle = colorTema('text');
        ctx.fillText('Residual: |θ₂ - θ₂,lin| = ' + (hist[hist.length - 1].diff).toFixed(4) + ' rad', 40, 20);
        ctx.restore();
    }

    function dibujarGraficaNube() {
        var info = ajustarCanvas(dom.nubeCanvas);
        if (!info) return;
        var ctx = info.ctx, w = info.width, h = info.height;
        ctx.clearRect(0, 0, w, h);

        ctx.font = '600 11px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
        ctx.fillStyle = colorTema('text');
        ctx.fillText('Distribución en vivo de vueltas (N = 200)', 40, 18);

        var flips = estado.nubeFlips;
        var maxF = 0;
        for (var i = 0; i < flips.length; i++) if (flips[i] > maxF) maxF = flips[i];
        maxF = Math.max(5, maxF);

        var conteo = new Int32Array(maxF + 1);
        for (var j = 0; j < flips.length; j++) conteo[flips[j]]++;

        var maxConteo = 1;
        for (var k = 0; k <= maxF; k++) if (conteo[k] > maxConteo) maxConteo = conteo[k];

        var bw = (w - 50) / (maxF + 1);
        ctx.save();
        for (var b = 0; b <= maxF; b++) {
            var bh = (conteo[b] / maxConteo) * (h - 45);
            var bx = 40 + b * bw;
            var by = (h - 20) - bh;

            ctx.fillStyle = b === 0 ? colorTema('th1') : colorTema('th2');
            ctx.fillRect(bx, by, Math.max(2, bw - 3), bh);

            ctx.font = '9px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
            ctx.fillStyle = colorTema('textMuted');
            ctx.fillText(b, bx + bw * 0.2, h - 6);
        }
        ctx.restore();
    }

    // Auxiliares de marcos de gráficas
    function dibujarMarcoYRejilla(ctx, w, h, etiquetaY) {
        ctx.save();
        ctx.strokeStyle = colorTema('border');
        ctx.lineWidth = 1;
        ctx.strokeRect(34, 4, w - 42, h - 22);

        // Línea central horizontal
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(34, h * 0.5);
        ctx.lineTo(w - 8, h * 0.5);
        ctx.stroke();

        ctx.font = '500 9px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
        ctx.fillStyle = colorTema('textMuted');
        ctx.fillText(etiquetaY, 2, 12);
        ctx.fillText('0', 20, h * 0.5 + 3);
        ctx.fillText('t (s)', w - 30, h - 4);
        ctx.restore();
    }

    function dibujarMarcoFase(ctx, w, h) {
        ctx.save();
        ctx.strokeStyle = colorTema('border');
        ctx.lineWidth = 1;
        ctx.strokeRect(34, 4, w - 42, h - 22);

        var cx = 34 + (w - 42) * 0.5;
        var cy = h * 0.5;

        // Ejes ortogonales
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(34, cy);
        ctx.lineTo(w - 8, cy);
        ctx.moveTo(cx, 4);
        ctx.lineTo(cx, h - 18);
        ctx.stroke();

        ctx.font = '500 9px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
        ctx.fillStyle = colorTema('textMuted');
        ctx.fillText('-π', 36, cy + 12);
        ctx.fillText('+π', w - 24, cy + 12);
        ctx.fillText('ω', cx + 4, 14);
        ctx.fillText('θ', w - 18, cy - 4);
        ctx.restore();
    }

    // =========================================================================
    // 10. MAPA DE CAOS ADIMENSIONAL (FTLE) CON WEB WORKER (BLOB)
    // =========================================================================
    var CODIGO_WORKER_FTLE = [
        "self.onmessage = function (e) {",
        "    var d = e.data;",
        "    var res = d.res;",
        "    var mu = d.mu;",
        "    var lambda = d.lambda;",
        "    var tMax = d.tMax;",
        "    var dt = d.dt;",
        "    var p = { m1: 1.0, m2: mu, l1: 1.0, l2: lambda, g: 1.0 };",
        "    var ftleGrid = new Float32Array(res * res);",
        "    var d0 = 1e-7;",
        "    for (var r = 0; r < res; r++) {",
        "        var th2 = -Math.PI + (2.0 * Math.PI * r) / (res - 1);",
        "        for (var c = 0; c < res; c++) {",
        "            var th1 = -Math.PI + (2.0 * Math.PI * c) / (res - 1);",
        "            var xA = [th1, th2, 0.0, 0.0];",
        "            var xB = [th1 + d0, th2, 0.0, 0.0];",
        "            var sumLog = 0.0;",
        "            var nPasos = Math.floor(tMax / dt);",
        "            for (var paso = 0; paso < nPasos; paso++) {",
        "                xA = pasoRK4Worker(xA, dt, p);",
        "                xB = pasoRK4Worker(xB, dt, p);",
        "                if ((paso % 8) === 0) {",
        "                    var dTh1 = xB[0] - xA[0];",
        "                    var dTh2 = xB[1] - xA[1];",
        "                    var dist = Math.sqrt(dTh1*dTh1 + dTh2*dTh2);",
        "                    if (dist > 1e-15) {",
        "                        sumLog += Math.log(dist / d0);",
        "                        xB[0] = xA[0] + (dTh1 / dist) * d0;",
        "                        xB[1] = xA[1] + (dTh2 / dist) * d0;",
        "                    }",
        "                }",
        "            }",
        "            ftleGrid[r * res + c] = Math.max(0.0, sumLog / tMax);",
        "        }",
        "        self.postMessage({ progreso: (r + 1) / res });",
        "    }",
        "    self.postMessage({ completado: true, grid: ftleGrid });",
        "};",
        "function campoDerivWorker(x, p) {",
        "    var th1 = x[0], th2 = x[1], w1 = x[2], w2 = x[3];",
        "    var delta = th1 - th2;",
        "    var cosD = Math.cos(delta), sinD = Math.sin(delta);",
        "    var m11 = (p.m1 + p.m2) * p.l1 * p.l1;",
        "    var m12 = p.m2 * p.l1 * p.l2 * cosD;",
        "    var m22 = p.m2 * p.l2 * p.l2;",
        "    var detM = p.m2 * (p.l1*p.l1) * (p.l2*p.l2) * (p.m1 + p.m2 * sinD * sinD);",
        "    var f1 = -p.m2 * p.l1 * p.l2 * sinD * (w2*w2) - (p.m1 + p.m2) * p.g * p.l1 * Math.sin(th1);",
        "    var f2 = p.m2 * p.l1 * p.l2 * sinD * (w1*w1) - p.m2 * p.g * p.l2 * Math.sin(th2);",
        "    return [w1, w2, (f1*m22 - f2*m12)/detM, (m11*f2 - m12*f1)/detM];",
        "}",
        "function pasoRK4Worker(x, dt, p) {",
        "    var k1 = campoDerivWorker(x, p);",
        "    var dt3 = dt/3.0, dt8 = dt/8.0;",
        "    var x2 = [x[0]+dt3*k1[0], x[1]+dt3*k1[1], x[2]+dt3*k1[2], x[3]+dt3*k1[3]];",
        "    var k2 = campoDerivWorker(x2, p);",
        "    var x3 = [x[0]-dt3*k1[0]+dt*k2[0], x[1]-dt3*k1[1]+dt*k2[1], x[2]-dt3*k1[2]+dt*k2[2], x[3]-dt3*k1[3]+dt*k2[3]];",
        "    var k3 = campoDerivWorker(x3, p);",
        "    var x4 = [x[0]+dt*k1[0]-dt*k2[0]+dt*k3[0], x[1]+dt*k1[1]-dt*k2[1]+dt*k3[1], x[2]+dt*k1[2]-dt*k2[2]+dt*k3[2], x[3]+dt*k1[3]-dt*k2[3]+dt*k3[3]];",
        "    var k4 = campoDerivWorker(x4, p);",
        "    return [",
        "        x[0]+dt8*(k1[0]+3*k2[0]+3*k3[0]+k4[0]),",
        "        x[1]+dt8*(k1[1]+3*k2[1]+3*k3[1]+k4[1]),",
        "        x[2]+dt8*(k1[2]+3*k2[2]+3*k3[2]+k4[2]),",
        "        x[3]+dt8*(k1[3]+3*k2[3]+3*k3[3]+k4[3])",
        "    ];",
        "}"
    ].join('\n');

    function iniciarCalculoFtle() {
        if (estado.ftleCalculando) return;
        estado.ftleCalculando = true;

        var res = (window.innerWidth < 1200) ? 96 : 128;
        var mu = estado.params.m2 / estado.params.m1;
        var lambda = estado.params.l2 / estado.params.l1;

        // Vista previa rápida síncrona inmediata para carga instantánea
        if (!estado.ftleData) {
            var resP = 40;
            var gridP = new Float32Array(resP * resP);
            var pP = { m1: 1.0, m2: mu, l1: 1.0, l2: lambda, g: 1.0 };
            var dtP = 0.05;
            var tMaxP = 10.0;
            var nP = Math.floor(tMaxP / dtP);
            for (var rP = 0; rP < resP; rP++) {
                var th2P = -Math.PI + (2.0 * Math.PI * rP) / (resP - 1);
                for (var cP = 0; cP < resP; cP++) {
                    var th1P = -Math.PI + (2.0 * Math.PI * cP) / (resP - 1);
                    var xAP = [th1P, th2P, 0.0, 0.0];
                    var xBP = [th1P + 1e-7, th2P, 0.0, 0.0];
                    var sLogP = 0.0;
                    for (var sP = 0; sP < nP; sP++) {
                        xAP = PenduloFisica.pasoRK4(xAP, dtP, pP);
                        xBP = PenduloFisica.pasoRK4(xBP, dtP, pP);
                        if ((sP % 8) === 0) {
                            var d1 = xBP[0] - xAP[0];
                            var d2 = xBP[1] - xAP[1];
                            var distP = Math.hypot(d1, d2);
                            if (distP > 1e-15) {
                                sLogP += Math.log(distP / 1e-7);
                                xBP[0] = xAP[0] + (d1 / distP) * 1e-7;
                                xBP[1] = xAP[1] + (d2 / distP) * 1e-7;
                            }
                        }
                    }
                    gridP[rP * resP + cP] = Math.max(0.0, sLogP / tMaxP);
                }
            }
            estado.ftleData = { res: resP, grid: gridP, mu: mu, lambda: lambda };
            renderizarMapaFtle();
        }

        var blob = new Blob([CODIGO_WORKER_FTLE], { type: 'application/javascript' });
        var worker = new Worker(URL.createObjectURL(blob));

        worker.onmessage = function (ev) {
            if (ev.data.completado) {
                estado.ftleData = {
                    res: res,
                    grid: ev.data.grid,
                    mu: mu,
                    lambda: lambda
                };
                estado.ftleCalculando = false;
                worker.terminate();
                renderizarMapaFtle();
            }
        };

        worker.postMessage({
            res: res,
            mu: mu,
            lambda: lambda,
            tMax: 20.0,
            dt: 0.02
        });
    }

    function renderizarMapaFtle() {
        var info = ajustarCanvas(dom.ftleCanvas);
        if (!info || !estado.ftleData) return;
        var ctx = info.ctx, w = info.width, h = info.height;
        ctx.clearRect(0, 0, w, h);

        var data = estado.ftleData;
        var res = data.res;
        var grid = data.grid;

        var padL = 28, padB = 20, padR = 30, padT = 10;
        var mapW = w - padL - padR;
        var mapH = h - padT - padB;

        // Renderizar mapa de color en un canvas auxiliar
        var offCanvas = document.createElement('canvas');
        offCanvas.width = res;
        offCanvas.height = res;
        var offCtx = offCanvas.getContext('2d');
        var imgData = offCtx.createImageData(res, res);

        var maxVal = 0.55;
        for (var r = 0; r < res; r++) {
            for (var c = 0; c < res; c++) {
                var val = grid[r * res + c];
                var t = Math.min(1.0, val / maxVal);

                // Colormap Magma / Inferno
                var rgb = colormapInferno(t);
                var idx = (r * res + c) * 4;
                imgData.data[idx] = rgb[0];
                imgData.data[idx + 1] = rgb[1];
                imgData.data[idx + 2] = rgb[2];
                imgData.data[idx + 3] = 255;
            }
        }
        offCtx.putImageData(imgData, 0, 0);

        // Dibujar mapa escalado
        ctx.drawImage(offCanvas, padL, padT, mapW, mapH);

        // Contornos de energía Ec_m2 y Ec_m1
        dibujarContornosEnergiaFtle(ctx, padL, padT, mapW, mapH);

        // Ejes y marcas
        ctx.save();
        ctx.strokeStyle = colorTema('border');
        ctx.lineWidth = 1;
        ctx.strokeRect(padL, padT, mapW, mapH);

        ctx.font = '500 9px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
        ctx.fillStyle = colorTema('textMuted');
        ctx.fillText('-π', padL, h - 8);
        ctx.fillText('0', padL + mapW * 0.5 - 3, h - 8);
        ctx.fillText('+π', padL + mapW - 12, h - 8);
        ctx.fillText('θ₁', padL + mapW * 0.5 - 5, h - 1);

        ctx.fillText('+π', 2, padT + 8);
        ctx.fillText('0', 12, padT + mapH * 0.5 + 3);
        ctx.fillText('-π', 2, padT + mapH - 2);
        ctx.fillText('θ₂', 2, padT + mapH * 0.5 - 8);

        // Barra de color lateral
        var cbX = padL + mapW + 8;
        var cbW = 10;
        var cbH = mapH;
        for (var y = 0; y < cbH; y++) {
            var frac = 1.0 - (y / cbH);
            var col = colormapInferno(frac);
            ctx.fillStyle = 'rgb(' + col[0] + ',' + col[1] + ',' + col[2] + ')';
            ctx.fillRect(cbX, padT + y, cbW, 1);
        }
        ctx.strokeStyle = colorTema('border');
        ctx.strokeRect(cbX, padT, cbW, cbH);
        ctx.fillText('λ·t₀', cbX - 2, padT - 3);
        ctx.fillText('0.5', cbX + cbW + 2, padT + 8);
        ctx.fillText('0', cbX + cbW + 2, padT + cbH);
        ctx.restore();

        dibujarReticulaFtle();
    }

    function dibujarContornosEnergiaFtle(ctx, padL, padT, mapW, mapH) {
        // V(th1, th2) = (m1+m2)*g*l1*(1 - cos(th1)) + m2*g*l2*(1 - cos(th2))
        // En unidades adimensionales: V_adim(th1, th2) = (1 + mu)*(1 - cos th1) + mu*lambda*(1 - cos th2)
        // Ec_m2 = 2 * mu * lambda
        // Ec_m1 = 2 * (1 + mu)
        var mu = estado.params.m2 / estado.params.m1;
        var lambda = estado.params.l2 / estado.params.l1;
        var ecM2 = 2.0 * mu * lambda;
        var ecM1 = 2.0 * (1.0 + mu);

        ctx.save();
        ctx.lineWidth = 1.8;

        // Contorno Ec_m2 (Cian)
        ctx.strokeStyle = colorTema('contourM2');
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        for (var i = 0; i <= 80; i++) {
            var th1 = -Math.PI + (2.0 * Math.PI * i) / 80;
            var term1 = (1.0 + mu) * (1.0 - Math.cos(th1));
            var rem = ecM2 - term1;
            if (rem >= 0 && rem <= 2.0 * mu * lambda) {
                var cosTh2 = 1.0 - rem / (mu * lambda);
                cosTh2 = Math.max(-1.0, Math.min(1.0, cosTh2));
                var th2A = Math.acos(cosTh2);
                var px = padL + ((th1 + Math.PI) / (2 * Math.PI)) * mapW;
                var pyA = padT + ((th2A + Math.PI) / (2 * Math.PI)) * mapH;
                var pyB = padT + ((-th2A + Math.PI) / (2 * Math.PI)) * mapH;
                ctx.rect(px, pyA, 1.5, 1.5);
                ctx.rect(px, pyB, 1.5, 1.5);
            }
        }
        ctx.stroke();

        // Contorno Ec_m1 (Ámbar / Negro)
        ctx.strokeStyle = colorTema('contourM1');
        ctx.setLineDash([2, 3]);
        ctx.beginPath();
        for (var j = 0; j <= 80; j++) {
            var th1B = -Math.PI + (2.0 * Math.PI * j) / 80;
            var term1B = (1.0 + mu) * (1.0 - Math.cos(th1B));
            var remB = ecM1 - term1B;
            if (remB >= 0 && remB <= 2.0 * mu * lambda) {
                var cosTh2B = 1.0 - remB / (mu * lambda);
                cosTh2B = Math.max(-1.0, Math.min(1.0, cosTh2B));
                var th2BA = Math.acos(cosTh2B);
                var pxB = padL + ((th1B + Math.PI) / (2 * Math.PI)) * mapW;
                var pyBA = padT + ((th2BA + Math.PI) / (2 * Math.PI)) * mapH;
                var pyBB = padT + ((-th2BA + Math.PI) / (2 * Math.PI)) * mapH;
                ctx.rect(pxB, pyBA, 1.5, 1.5);
                ctx.rect(pxB, pyBB, 1.5, 1.5);
            }
        }
        ctx.stroke();
        ctx.restore();
    }

    function dibujarReticulaFtle() {
        if (!estado.ftleData || estado.modo !== 'normal') return;
        var info = ajustarCanvas(dom.ftleCanvas);
        if (!info) return;
        var ctx = info.ctx, w = info.width, h = info.height;

        var padL = 28, padB = 20, padR = 30, padT = 10;
        var mapW = w - padL - padR;
        var mapH = h - padT - padB;

        // Posición actual x0
        var th1Wrap = ((estado.x0[0] + Math.PI) % (2 * Math.PI)) - Math.PI;
        var th2Wrap = ((estado.x0[1] + Math.PI) % (2 * Math.PI)) - Math.PI;

        var px = padL + ((th1Wrap + Math.PI) / (2 * Math.PI)) * mapW;
        var py = padT + ((th2Wrap + Math.PI) / (2 * Math.PI)) * mapH;

        ctx.save();
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 2.0;
        ctx.shadowColor = 'rgba(0,0,0,0.8)';
        ctx.shadowBlur = 4;

        // Cruz con círculo
        ctx.beginPath();
        ctx.arc(px, py, 5, 0, 2 * Math.PI);
        ctx.moveTo(px - 9, py);
        ctx.lineTo(px + 9, py);
        ctx.moveTo(px, py - 9);
        ctx.lineTo(px, py + 9);
        ctx.stroke();
        ctx.restore();
    }

    function colormapInferno(t) {
        // Rampa suave estilo Inferno
        t = Math.max(0.0, Math.min(1.0, t));
        var r = Math.round(255 * Math.min(1.0, 1.5 * Math.pow(t, 0.8)));
        var g = Math.round(255 * Math.pow(t, 2.0));
        var b = Math.round(255 * (t < 0.6 ? 1.4 * Math.sin(t * Math.PI) : 0.4 * (1.0 - t)));
        return [Math.min(255, r), Math.min(255, g), Math.min(255, b)];
    }

    // =========================================================================
    // 11. LABORATORIO DE MONTE CARLO Y REPRODUCCIÓN ESTOCÁSTICA
    // =========================================================================
    var CODIGO_WORKER_MC = [
        "self.onmessage = function (e) {",
        "    var d = e.data;",
        "    var N = d.N, sigma = d.sigma, rho = d.rho, T = d.T, dt = 0.005, baseTh1 = d.baseTh1, baseTh2 = d.baseTh2;",
        "    var seed = d.seed;",
        "    var p = d.params;",
        "    var rng = mulberry32Worker(seed);",
        "    var flipsArray = new Int32Array(N);",
        "    var nPasos = Math.floor(T / dt);",
        "    var l11 = sigma;",
        "    var l21 = rho * sigma;",
        "    var l22 = sigma * Math.sqrt(Math.max(0.0, 1.0 - rho*rho));",
        "    for (var i = 0; i < N; i++) {",
        "        var u1 = Math.max(1e-12, rng()), u2 = rng();",
        "        var z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);",
        "        var z2 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);",
        "        var th1 = baseTh1 + l11 * z1;",
        "        var th2 = baseTh2 + l21 * z1 + l22 * z2;",
        "        var x = [th1, th2, 0.0, 0.0];",
        "        var th2Cont = th2;",
        "        for (var step = 0; step < nPasos; step++) {",
        "            var prevTh2 = x[1];",
        "            x = pasoRK4Worker(x, dt, p);",
        "            var dTh2 = ((x[1] - prevTh2 + Math.PI) % (2.0 * Math.PI)) - Math.PI;",
        "            th2Cont += dTh2;",
        "        }",
        "        flipsArray[i] = Math.abs(Math.floor((th2Cont + Math.PI) / (2.0 * Math.PI)));",
        "        if ((i % 50) === 0) self.postMessage({ progreso: (i + 1) / N });",
        "    }",
        "    self.postMessage({ completado: true, flips: Array.from(flipsArray) });",
        "};",
        "function mulberry32Worker(seed) {",
        "    var s = seed | 0;",
        "    return function () {",
        "        s = (s + 0x6D2B79F5) | 0;",
        "        var t = Math.imul(s ^ (s >>> 15), 1 | s);",
        "        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;",
        "        return ((t ^ (t >>> 14)) >>> 0) / 4294967296.0;",
        "    };",
        "}",
        "function campoDerivWorker(x, p) {",
        "    var th1 = x[0], th2 = x[1], w1 = x[2], w2 = x[3];",
        "    var delta = th1 - th2;",
        "    var cosD = Math.cos(delta), sinD = Math.sin(delta);",
        "    var m11 = (p.m1 + p.m2) * p.l1 * p.l1;",
        "    var m12 = p.m2 * p.l1 * p.l2 * cosD;",
        "    var m22 = p.m2 * p.l2 * p.l2;",
        "    var detM = p.m2 * (p.l1*p.l1) * (p.l2*p.l2) * (p.m1 + p.m2 * sinD * sinD);",
        "    var f1 = -p.m2 * p.l1 * p.l2 * sinD * (w2*w2) - (p.m1 + p.m2) * p.g * p.l1 * Math.sin(th1);",
        "    var f2 = p.m2 * p.l1 * p.l2 * sinD * (w1*w1) - p.m2 * p.g * p.l2 * Math.sin(th2);",
        "    return [w1, w2, (f1*m22 - f2*m12)/detM, (m11*f2 - m12*f1)/detM];",
        "}",
        "function pasoRK4Worker(x, dt, p) {",
        "    var k1 = campoDerivWorker(x, p);",
        "    var dt3 = dt/3.0, dt8 = dt/8.0;",
        "    var x2 = [x[0]+dt3*k1[0], x[1]+dt3*k1[1], x[2]+dt3*k1[2], x[3]+dt3*k1[3]];",
        "    var k2 = campoDerivWorker(x2, p);",
        "    var x3 = [x[0]-dt3*k1[0]+dt*k2[0], x[1]-dt3*k1[1]+dt*k2[1], x[2]-dt3*k1[2]+dt*k2[2], x[3]-dt3*k1[3]+dt*k2[3]];",
        "    var k3 = campoDerivWorker(x3, p);",
        "    var x4 = [x[0]+dt*k1[0]-dt*k2[0]+dt*k3[0], x[1]+dt*k1[1]-dt*k2[1]+dt*k3[1], x[2]+dt*k1[2]-dt*k2[2]+dt*k3[2], x[3]+dt*k1[3]-dt*k2[3]+dt*k3[3]];",
        "    var k4 = campoDerivWorker(x4, p);",
        "    return [",
        "        x[0]+dt8*(k1[0]+3*k2[0]+3*k3[0]+k4[0]),",
        "        x[1]+dt8*(k1[1]+3*k2[1]+3*k3[1]+k4[1]),",
        "        x[2]+dt8*(k1[2]+3*k2[2]+3*k3[2]+k4[2]),",
        "        x[3]+dt8*(k1[3]+3*k2[3]+3*k3[3]+k4[3])",
        "    ];",
        "}"
    ].join('\n');

    var estadoMc = {
        corriendo: false,
        flipsActuales: null,
        flipsBaja: null,
        modoComparacion: false
    };

    function ejecutarMonteCarlo(esBajaEnergia) {
        if (estadoMc.corriendo) return;
        estadoMc.corriendo = true;
        dom.mcProgressFill.style.width = '0%';
        dom.mcProgressText.textContent = 'Simulando trayectorias...';

        var N = parseInt(dom.mcSliderN.value, 10);
        var sigma = parseFloat(dom.mcSliderSigma.value);
        var rho = parseFloat(dom.mcSliderRho.value);
        var T = parseFloat(dom.mcSliderT.value);
        var seed = parseInt(dom.mcSeedInput.value, 10) || 20261119;

        var baseTh1 = esBajaEnergia ? 0.6 : estado.x0[0];
        var baseTh2 = esBajaEnergia ? 0.6 : estado.x0[1];

        var blob = new Blob([CODIGO_WORKER_MC], { type: 'application/javascript' });
        var worker = new Worker(URL.createObjectURL(blob));

        worker.onmessage = function (ev) {
            if (ev.data.progreso !== undefined) {
                var p = Math.round(ev.data.progreso * 100);
                dom.mcProgressFill.style.width = p + '%';
                dom.mcProgressText.textContent = 'Calculando ensemble (' + p + '%)...';
            } else if (ev.data.completado) {
                var flips = ev.data.flips;
                worker.terminate();
                estadoMc.corriendo = false;
                dom.mcProgressFill.style.width = '100%';
                dom.mcProgressText.textContent = 'Simulación finalizada';

                if (esBajaEnergia) {
                    estadoMc.flipsBaja = flips;
                    estadoMc.modoComparacion = true;
                } else {
                    estadoMc.flipsActuales = flips;
                }

                // Cálculo estadístico
                procesarEstadisticasMonteCarlo(flips, seed);
                renderizarHistogramaMonteCarlo();
                dom.btnMcReplayCloud.style.display = 'inline-flex';
            }
        };

        worker.postMessage({
            N: N,
            sigma: sigma,
            rho: rho,
            T: T,
            seed: seed,
            baseTh1: baseTh1,
            baseTh2: baseTh2,
            params: estado.params
        });
    }

    function procesarEstadisticasMonteCarlo(flips, seed) {
        var rng = PenduloFisica.mulberry32(seed + 999);
        var bootRes = PenduloFisica.bootstrapPercentil(flips, 2000, rng, 0.95);

        var kExitos = 0;
        for (var i = 0; i < flips.length; i++) {
            if (flips[i] > 0) kExitos++;
        }
        var wilsonRes = PenduloFisica.intervaloWilson(kExitos, flips.length);

        dom.mcValMean.textContent = bootRes.media.toFixed(2);
        dom.mcCiMean.textContent = 'IC 95%: [' + bootRes.media_ci_bajo.toFixed(2) + ', ' + bootRes.media_ci_alto.toFixed(2) + ']';

        dom.mcValStd.textContent = bootRes.std.toFixed(2);
        dom.mcCiStd.textContent = 'IC 95%: [' + bootRes.std_ci_bajo.toFixed(2) + ', ' + bootRes.std_ci_alto.toFixed(2) + ']';

        dom.mcValPflip.textContent = (wilsonRes.p * 100).toFixed(1) + '%';
        dom.mcCiWilson.textContent = 'Wilson 95%: [' + (wilsonRes.bajo * 100).toFixed(1) + '%, ' + (wilsonRes.alto * 100).toFixed(1) + '%]';
    }

    function renderizarHistogramaMonteCarlo() {
        var info = ajustarCanvas(dom.mcHistogramCanvas);
        if (!info) return;
        var ctx = info.ctx, w = info.width, h = info.height;
        ctx.clearRect(0, 0, w, h);

        var fA = estadoMc.flipsActuales;
        var fB = estadoMc.modoComparacion ? estadoMc.flipsBaja : null;
        if (!fA) return;

        var maxVueltas = 0;
        for (var i = 0; i < fA.length; i++) if (fA[i] > maxVueltas) maxVueltas = fA[i];
        if (fB) {
            for (var j = 0; j < fB.length; j++) if (fB[j] > maxVueltas) maxVueltas = fB[j];
        }
        maxVueltas = Math.max(8, maxVueltas);

        var conteoA = new Int32Array(maxVueltas + 1);
        for (var k = 0; k < fA.length; k++) conteoA[fA[k]]++;

        var conteoB = fB ? new Int32Array(maxVueltas + 1) : null;
        if (fB) {
            for (var m = 0; m < fB.length; m++) conteoB[fB[m]]++;
        }

        var maxFreq = 1;
        for (var n = 0; n <= maxVueltas; n++) {
            if (conteoA[n] > maxFreq) maxFreq = conteoA[n];
            if (conteoB && conteoB[n] > maxFreq) maxFreq = conteoB[n];
        }

        var padL = 40, padB = 25, padT = 25, padR = 20;
        var plotW = w - padL - padR;
        var plotH = h - padT - padB;
        var bw = plotW / (maxVueltas + 1);

        // Marco
        ctx.save();
        ctx.strokeStyle = colorTema('border');
        ctx.lineWidth = 1;
        ctx.strokeRect(padL, padT, plotW, plotH);

        // Ejes
        ctx.font = '600 10px ' + getComputedStyle(dom.body).getPropertyValue('--font-mono');
        ctx.fillStyle = colorTema('textMuted');
        ctx.fillText('Frecuencia', 2, padT - 8);
        ctx.fillText('Número de Vueltas de la Masa 2', padL + plotW * 0.35, h - 6);

        // Barras
        for (var v = 0; v <= maxVueltas; v++) {
            var bx = padL + v * bw;

            if (conteoB) {
                // Modo A/B Comparativo: dos barras lado a lado
                var wSub = (bw - 4) * 0.48;
                // Barra A (Condición Actual - Magenta)
                var bhA = (conteoA[v] / maxFreq) * plotH;
                ctx.fillStyle = colorTema('th2');
                ctx.fillRect(bx + 2, padT + plotH - bhA, wSub, bhA);

                // Barra B (Baja Energía - Cian)
                var bhB = (conteoB[v] / maxFreq) * plotH;
                ctx.fillStyle = colorTema('th1');
                ctx.fillRect(bx + 2 + wSub, padT + plotH - bhB, wSub, bhB);
            } else {
                // Solo condición actual
                var bh = (conteoA[v] / maxFreq) * plotH;
                ctx.fillStyle = v === 0 ? colorTema('th1') : colorTema('th2');
                ctx.fillRect(bx + 2, padT + plotH - bh, Math.max(2, bw - 4), bh);
            }

            // Etiquetas eje X
            if (maxVueltas <= 20 || (v % 2 === 0)) {
                ctx.fillStyle = colorTema('textMuted');
                ctx.fillText(v, bx + bw * 0.25, padT + plotH + 14);
            }
        }

        // Leyenda
        ctx.font = '600 11px ' + getComputedStyle(dom.body).getPropertyValue('--font-sans');
        if (conteoB) {
            ctx.fillStyle = colorTema('th2');
            ctx.fillText('■ Condición Actual', padL + 10, padT - 8);
            ctx.fillStyle = colorTema('th1');
            ctx.fillText('■ Baja Energía (0.6, 0.6 rad)', padL + 150, padT - 8);
        } else {
            ctx.fillStyle = colorTema('text');
            ctx.fillText('Histograma de Pasos por Arriba (N = ' + fA.length + ')', padL + 10, padT - 8);
        }
        ctx.restore();
    }

    // =========================================================================
    // 12. EVENTOS DE INTERFAZ Y CONTROLES
    // =========================================================================
    function vincularEventos() {
        // Alternancia de tema
        dom.btnThemeToggle.addEventListener('click', alternarTema);

        // Cajón plegable de controles
        function toggleDrawer() {
            dom.drawerPanel.classList.toggle('open');
        }
        dom.btnDrawerToggle.addEventListener('click', toggleDrawer);
        dom.btnDrawerClose.addEventListener('click', toggleDrawer);

        // Teclas rápidas
        window.addEventListener('keydown', function (e) {
            if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
            if (e.key === 'c' || e.key === 'C') {
                toggleDrawer();
            } else if (e.key === ' ') {
                e.preventDefault();
                togglePausa();
            } else if (e.key === 'r' || e.key === 'R') {
                reiniciarSimulacion();
            } else if (e.key === 's' || e.key === 'S') {
                estado.corriendo = false;
                dom.btnPlayPause.textContent = 'Reanudar';
                pasoSimulacion(estado.dtBase);
            } else if (e.key === 'p' || e.key === 'P') {
                alternarPresentacion();
            }
        });

        // Modo presentación
        function alternarPresentacion() {
            dom.body.classList.toggle('mode-presentation');
            if (dom.body.classList.contains('mode-presentation')) {
                dom.btnPresentation.style.background = colorTema('energy');
            } else {
                dom.btnPresentation.style.background = '';
            }
        }
        dom.btnPresentation.addEventListener('click', alternarPresentacion);

        // Botón Play / Pausa
        function togglePausa() {
            estado.corriendo = !estado.corriendo;
            dom.btnPlayPause.textContent = estado.corriendo ? 'Pausa' : 'Reanudar';
            if (estado.corriendo) {
                dom.btnPlayPause.classList.add('btn-primary');
            } else {
                dom.btnPlayPause.classList.remove('btn-primary');
            }
        }
        dom.btnPlayPause.addEventListener('click', togglePausa);

        // Botón Paso
        dom.btnStep.addEventListener('click', function () {
            estado.corriendo = false;
            dom.btnPlayPause.textContent = 'Reanudar';
            dom.btnPlayPause.classList.remove('btn-primary');
            pasoSimulacion(estado.dtBase);
        });

        // Botón Reiniciar
        dom.btnReset.addEventListener('click', reiniciarSimulacion);

        // Selector de velocidad
        dom.selectSpeed.addEventListener('change', function () {
            estado.velocidad = parseFloat(this.value);
        });

        // Checkbox estela
        dom.checkTrace.addEventListener('change', function () {
            estado.mostrarEstela = this.checked;
            if (!this.checked) estado.estelaM2 = [];
        });

        // Selector de integrador
        dom.btnIntegRK4.addEventListener('click', function () {
            estado.algoritmo = 'rk4_38';
            dom.btnIntegRK4.classList.add('active');
            dom.btnIntegClasico.classList.remove('active');
        });
        dom.btnIntegClasico.addEventListener('click', function () {
            estado.algoritmo = 'rk4_clasico';
            dom.btnIntegClasico.classList.add('active');
            dom.btnIntegRK4.classList.remove('active');
        });

        // Botones de gravedad rápida
        var btnGravs = document.querySelectorAll('.btn-gravity-preset');
        btnGravs.forEach(function (btn) {
            btn.addEventListener('click', function () {
                var gVal = parseFloat(this.getAttribute('data-g'));
                estado.params.g = gVal;
                dom.sliderG.value = gVal;
                dom.lblG.textContent = gVal.toFixed(2) + ' m/s²';
                actualizarUmbralesEnergiaUI();
            });
        });

        // Slider de gravedad
        dom.sliderG.addEventListener('input', function () {
            var val = parseFloat(this.value);
            estado.params.g = val;
            dom.lblG.textContent = val.toFixed(2) + ' m/s²';
            actualizarUmbralesEnergiaUI();
        });

        // Sliders de geometría y masa
        function setupSlider(slider, lbl, prop, fmt, unit) {
            slider.addEventListener('input', function () {
                var v = parseFloat(this.value);
                estado.params[prop] = v;
                lbl.textContent = v.toFixed(fmt) + (unit ? ' ' + unit : '');
                actualizarLabelsParametros();
                actualizarUmbralesEnergiaUI();
            });
        }
        setupSlider(dom.sliderL1, dom.lblL1, 'l1', 2, 'm');
        setupSlider(dom.sliderL2, dom.lblL2, 'l2', 2, 'm');
        setupSlider(dom.sliderM1, dom.lblM1, 'm1', 2, 'kg');
        setupSlider(dom.sliderM2, dom.lblM2, 'm2', 2, 'kg');

        // Sliders de Condiciones Iniciales
        dom.sliderTh1.addEventListener('input', function () {
            var deg = parseInt(this.value, 10);
            estado.x0[0] = deg * Math.PI / 180;
            dom.lblTh1.textContent = deg + '°';
            reiniciarSimulacion();
        });
        dom.sliderTh2.addEventListener('input', function () {
            var deg2 = parseInt(this.value, 10);
            estado.x0[1] = deg2 * Math.PI / 180;
            dom.lblTh2.textContent = deg2 + '°';
            reiniciarSimulacion();
        });
        dom.sliderW1.addEventListener('input', function () {
            var w1 = parseFloat(this.value);
            estado.x0[2] = w1;
            dom.lblW1.textContent = w1.toFixed(1) + ' rad/s';
            reiniciarSimulacion();
        });
        dom.sliderW2.addEventListener('input', function () {
            var w2 = parseFloat(this.value);
            estado.x0[3] = w2;
            dom.lblW2.textContent = w2.toFixed(1) + ' rad/s';
            reiniciarSimulacion();
        });

        // Selector de Pestañas de Modo
        var tabBtns = document.querySelectorAll('.mode-tab-btn');
        tabBtns.forEach(function (btn) {
            btn.addEventListener('click', function () {
                tabBtns.forEach(function (b) { b.classList.remove('active'); });
                this.classList.add('active');
                cambiarModo(this.getAttribute('data-mode'));
            });
        });

        // Pestañas del Retrato de Fase
        var phaseBtns = document.querySelectorAll('.phase-tab-btn');
        phaseBtns.forEach(function (btn) {
            btn.addEventListener('click', function () {
                phaseBtns.forEach(function (b) { b.classList.remove('active'); });
                this.classList.add('active');
                estado.subfaseRetrato = this.getAttribute('data-phase');
            });
        });

        // Botones de Presets
        var presetBtns = document.querySelectorAll('.btn-preset-card');
        presetBtns.forEach(function (btn) {
            btn.addEventListener('click', function () {
                var pId = this.getAttribute('data-preset');
                cargarPreset(pId);
            });
        });

        // Inspección interactiva del mapa FTLE (hover y click)
        dom.ftleCanvas.addEventListener('mousemove', inspeccionarFtle);
        dom.ftleCanvas.addEventListener('touchmove', function (e) {
            if (e.touches.length > 0) inspeccionarFtle(e.touches[0]);
        }, { passive: true });

        dom.ftleCanvas.addEventListener('click', fijarCondicionFtle);
        dom.ftleCanvas.addEventListener('touchend', function (e) {
            if (e.changedTouches.length > 0) fijarCondicionFtle(e.changedTouches[0]);
        });

        // Arrastre interactivo en canvas del péndulo (cuando está en pausa)
        dom.simCanvas.addEventListener('pointerdown', iniciarArrastreMasa);
        dom.simCanvas.addEventListener('pointermove', moverArrastreMasa);
        dom.simCanvas.addEventListener('pointerup', terminarArrastreMasa);
        dom.simCanvas.addEventListener('pointercancel', terminarArrastreMasa);

        // Monte Carlo UI
        dom.btnMcUseCurrent.addEventListener('click', function () {
            var d1 = Math.round(estado.x0[0] * 180 / Math.PI);
            var d2 = Math.round(estado.x0[1] * 180 / Math.PI);
            dom.mcBaseLabel.textContent = 'Actual (' + d1 + '°, ' + d2 + '°)';
        });
        dom.mcSliderN.addEventListener('input', function () {
            dom.mcLblN.textContent = this.value;
        });
        dom.mcSliderSigma.addEventListener('input', function () {
            dom.mcLblSigma.textContent = parseFloat(this.value).toFixed(3) + ' rad';
        });
        dom.mcSliderRho.addEventListener('input', function () {
            dom.mcLblRho.textContent = parseFloat(this.value).toFixed(2);
        });
        dom.mcSliderT.addEventListener('input', function () {
            dom.mcLblT.textContent = this.value + ' s';
        });
        dom.btnMcRun.addEventListener('click', function () {
            estadoMc.modoComparacion = false;
            ejecutarMonteCarlo(false);
        });
        dom.btnMcCompareLow.addEventListener('click', function () {
            ejecutarMonteCarlo(true);
        });
        dom.btnMcReplayCloud.addEventListener('click', function () {
            cambiarModo('nube');
            document.getElementById('tab-mode-nube').click();
        });
    }

    function cambiarModo(nuevoModo) {
        estado.modo = nuevoModo;

        // Ocultar todos los contenedores contextuales
        dom.ftleContainer.style.display = 'none';
        dom.gemeloContainer.style.display = 'none';
        dom.linealContainer.style.display = 'none';
        dom.nubeContainer.style.display = 'none';
        dom.montecarloPanel.style.display = 'none';

        if (nuevoModo === 'normal') {
            dom.ftleContainer.style.display = 'flex';
            dom.contextualTitle.textContent = 'Mapa de Caos Adimensional (FTLE)';
            dom.contextualMeta.textContent = 'τ = t√(g/l₁)';
        } else if (nuevoModo === 'gemelo') {
            dom.gemeloContainer.style.display = 'block';
            dom.contextualTitle.textContent = 'Separación Exponencial Gemelo';
            dom.contextualMeta.textContent = 'd(t) en log₁₀';
        } else if (nuevoModo === 'lineal') {
            dom.linealContainer.style.display = 'block';
            dom.contextualTitle.textContent = 'Lineal vs No Lineal (θ₂)';
            dom.contextualMeta.textContent = 'Modos Normales';
        } else if (nuevoModo === 'nube') {
            dom.nubeContainer.style.display = 'block';
            dom.contextualTitle.textContent = 'Distribución en Vivo Nube';
            dom.contextualMeta.textContent = 'Ensemble N = 200';
        } else if (nuevoModo === 'montecarlo') {
            dom.montecarloPanel.style.display = 'flex';
        }

        reiniciarSimulacion();
    }

    function cargarPreset(pId) {
        if (!window.PRESETS_PENDULO || !window.PRESETS_PENDULO[pId]) return;
        var p = window.PRESETS_PENDULO[pId];

        estado.x0 = [p.x0[0], p.x0[1], p.x0[2], p.x0[3]];
        estado.params.m1 = p.params.m1;
        estado.params.m2 = p.params.m2;
        estado.params.l1 = p.params.l1;
        estado.params.l2 = p.params.l2;
        estado.params.g = p.params.g;

        actualizarLabelsParametros();
        actualizarUmbralesEnergiaUI();
        reiniciarSimulacion();
    }

    // Inspección interactiva del mapa FTLE
    function inspeccionarFtle(e) {
        if (!estado.ftleData || dom.ftleContainer.style.display === 'none') return;
        var rect = dom.ftleCanvas.getBoundingClientRect();
        var mx = e.clientX - rect.left;
        var my = e.clientY - rect.top;

        var padL = 28, padB = 20, padR = 30, padT = 10;
        var mapW = rect.width - padL - padR;
        var mapH = rect.height - padT - padB;

        if (mx < padL || mx > padL + mapW || my < padT || my > padT + mapH) return;

        var fracX = (mx - padL) / mapW;
        var fracY = (my - padT) / mapH;

        var th1 = -Math.PI + fracX * (2.0 * Math.PI);
        var th2 = -Math.PI + fracY * (2.0 * Math.PI);

        var res = estado.ftleData.res;
        var c = Math.max(0, Math.min(res - 1, Math.floor(fracX * res)));
        var r = Math.max(0, Math.min(res - 1, Math.floor(fracY * res)));
        var ftleVal = estado.ftleData.grid[r * res + c];

        dom.ftleTh1.textContent = (th1 * 180 / Math.PI).toFixed(1) + '°';
        dom.ftleTh2.textContent = (th2 * 180 / Math.PI).toFixed(1) + '°';
        dom.ftleVal.textContent = ftleVal.toFixed(3);
    }

    function fijarCondicionFtle(e) {
        if (!estado.ftleData || dom.ftleContainer.style.display === 'none') return;
        var rect = dom.ftleCanvas.getBoundingClientRect();
        var mx = e.clientX - rect.left;
        var my = e.clientY - rect.top;

        var padL = 28, padB = 20, padR = 30, padT = 10;
        var mapW = rect.width - padL - padR;
        var mapH = rect.height - padT - padB;

        if (mx < padL || mx > padL + mapW || my < padT || my > padT + mapH) return;

        var fracX = (mx - padL) / mapW;
        var fracY = (my - padT) / mapH;

        var th1 = -Math.PI + fracX * (2.0 * Math.PI);
        var th2 = -Math.PI + fracY * (2.0 * Math.PI);

        estado.x0[0] = th1;
        estado.x0[1] = th2;
        estado.x0[2] = 0.0;
        estado.x0[3] = 0.0;

        actualizarLabelsParametros();
        reiniciarSimulacion();
    }

    // Arrastre interactivo de masas
    function iniciarArrastreMasa(e) {
        if (estado.corriendo) return; // Solo cuando está en pausa
        var rect = dom.simCanvas.getBoundingClientRect();
        var mx = e.clientX - rect.left;
        var my = e.clientY - rect.top;

        var cx = rect.width * 0.5;
        var cy = rect.height * 0.42;
        var lTotal = estado.params.l1 + estado.params.l2;
        var escala = (Math.min(rect.width, rect.height) * 0.40) / lTotal;

        var l1Px = estado.params.l1 * escala;
        var l2Px = estado.params.l2 * escala;

        var x1 = cx + l1Px * Math.sin(estado.x[0]);
        var y1 = cy + l1Px * Math.cos(estado.x[0]);
        var x2 = x1 + l2Px * Math.sin(estado.x[1]);
        var y2 = y1 + l2Px * Math.cos(estado.x[1]);

        var dist1 = Math.hypot(mx - x1, my - y1);
        var dist2 = Math.hypot(mx - x2, my - y2);

        if (dist2 < 30) {
            estado.arrastrandoMasa = 2;
        } else if (dist1 < 30) {
            estado.arrastrandoMasa = 1;
        }
    }

    function moverArrastreMasa(e) {
        if (!estado.arrastrandoMasa) return;
        var rect = dom.simCanvas.getBoundingClientRect();
        var mx = e.clientX - rect.left;
        var my = e.clientY - rect.top;

        var cx = rect.width * 0.5;
        var cy = rect.height * 0.42;

        if (estado.arrastrandoMasa === 1) {
            var ang1 = Math.atan2(mx - cx, my - cy);
            estado.x[0] = ang1;
            estado.x0[0] = ang1;
            estado.x[2] = 0.0;
            estado.x0[2] = 0.0;
        } else if (estado.arrastrandoMasa === 2) {
            var lTotal = estado.params.l1 + estado.params.l2;
            var escala = (Math.min(rect.width, rect.height) * 0.40) / lTotal;
            var l1Px = estado.params.l1 * escala;
            var x1 = cx + l1Px * Math.sin(estado.x[0]);
            var y1 = cy + l1Px * Math.cos(estado.x[0]);

            var ang2 = Math.atan2(mx - x1, my - y1);
            estado.x[1] = ang2;
            estado.x0[1] = ang2;
            estado.x[3] = 0.0;
            estado.x0[3] = 0.0;
        }

        actualizarLabelsParametros();
        estado.e0 = PenduloFisica.energiaTotal(estado.x, estado.params);
    }

    function terminarArrastreMasa() {
        estado.arrastrandoMasa = 0;
    }

    // =========================================================================
    // 13. BUCLE PRINCIPAL DE ANIMACIÓN A 60 FPS
    // =========================================================================
    var ultimoTimestamp = performance.now();

    function buclePrincipal(timestamp) {
        var dtReal = (timestamp - ultimoTimestamp) / 1000.0;
        ultimoTimestamp = timestamp;
        dtReal = Math.min(0.1, dtReal); // Limitar delta t para evitar saltos

        if (estado.corriendo) {
            var tiempoASimular = dtReal * estado.velocidad;
            var nSubpasos = Math.max(1, Math.round(tiempoASimular / estado.dtBase));
            var dtSub = tiempoASimular / nSubpasos;

            for (var s = 0; s < nSubpasos; s++) {
                pasoSimulacion(dtSub);
            }
        }

        // Actualizar HUD Cuantitativo
        dom.valTime.textContent = estado.t.toFixed(2) + ' s';
        var eK = PenduloFisica.energiaCinetica(estado.x, estado.params);
        var eV = PenduloFisica.energiaPotencial(estado.x, estado.params);
        var eTotal = eK + eV;

        dom.valEnergyTotal.textContent = eTotal.toFixed(2) + ' J';
        dom.valDeltaE.textContent = estado.deltaERel.toExponential(2);
        dom.valFlipsM2.textContent = estado.vueltasM2;

        dom.valEnergyK.textContent = eK.toFixed(2);
        dom.valEnergyV.textContent = eV.toFixed(2);
        dom.valEnergyE.textContent = eTotal.toFixed(2);

        // Actualizar barra de energía
        var u = PenduloFisica.umbralesEnergia(estado.params);
        var eEscala = Math.max(u.Ec_ambas, eTotal);
        var pctK = Math.min(100, Math.max(0, (eK / eEscala) * 100));
        var pctV = Math.min(100, Math.max(0, (eV / eEscala) * 100));
        dom.fillKinetic.style.width = pctK + '%';
        dom.fillPotential.style.width = pctV + '%';

        // Renderizado
        dibujarEscenario();
        dibujarGraficaAngulos();
        dibujarGraficaFase();
        dibujarGraficaEnergia();
        dibujarGraficaContextual();

        requestAnimationFrame(buclePrincipal);
    }

    // =========================================================================
    // 14. ARRANQUE
    // =========================================================================
    window.addEventListener('DOMContentLoaded', function () {
        vincularEventos();
        actualizarLabelsParametros();
        actualizarUmbralesEnergiaUI();
        if (window.location.search.indexOf('theme=light') !== -1) {
            alternarTema();
        }
        reiniciarSimulacion();
        iniciarCalculoFtle();
        if (window.location.search.indexOf('capture=1') !== -1) {
            for (var k = 0; k < 220; k++) {
                pasoSimulacion(0.01);
            }
        }
        requestAnimationFrame(buclePrincipal);
    });

})();
