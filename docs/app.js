/**
 * app.js - Controlador Principal del Simulador del Péndulo Doble Conservativo.
 * Sistemas Dinámicos II | ITAM - Fase 3c: Simplificación y Claridad Visual
 *
 * Características:
 * - Barra lateral izquierda fija (300 px) con scroll independiente.
 * - 4 pestañas directas: Péndulo, Gemelos (caos), Monte Carlo, Mapa de caos.
 * - Estilo visual de divulgación ("Mates Mike"): péndulo grande, ángulos en números grandes, cuñas sombreadas.
 * - Sin emojis en toda la interfaz (advertencia lineal en texto plano).
 * - Termómetro de energía horizontal con umbrales en palabras y el axioma "Poder dar la vuelta no implica caos."
 * - Gemelos con ajuste de exponente de Lyapunov dimensional (s⁻¹) y adimensional (λ·√(l₁/g)).
 * - Monte Carlo dual side-by-side (Baja energía vs Alta energía) con animación sincronizada de N=300 y estadísticos bootstrap/Wilson.
 * - Mapa de caos precargado (96×96) instantáneo, nunca negro, con una sola cruz interactiva.
 */

(function () {
    'use strict';

    // =========================================================================
    // 1. PALETAS DE COLOR Y TEMAS
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
            rod: '#d0d7de',
            rodGhost: 'rgba(150, 160, 180, 0.35)',
            th1: '#4cc9f0',        // Cian
            th2: '#f72585',        // Magenta
            energy: '#ffb703',     // Ámbar
            kinetic: '#2ec4b6',    // Verde Menta
            potential: '#9b5de5',  // Violeta
            contourM2: '#4cc9f0',
            contourM1: '#ffb703',
            cross: '#ffffff',
            traceSlow: [76, 201, 240],
            traceMid: [255, 183, 3],
            traceFast: [247, 37, 133]
        },
        light: {
            bg: '#f4f6fa',
            panel: '#ffffff',
            panelSub: '#f8fafc',
            border: '#cbd5e1',
            text: '#0f172a',
            textMuted: '#475569',
            grid: '#e2e8f0',
            rod: '#1e293b',
            rodGhost: 'rgba(100, 116, 139, 0.35)',
            th1: '#0284c7',        // Cian oscuro
            th2: '#d946ef',        // Magenta
            energy: '#b45309',     // Ámbar
            kinetic: '#0d9488',    // Menta oscuro
            potential: '#7c3aed',  // Violeta
            contourM2: '#0284c7',
            contourM1: '#111827',
            cross: '#000000',
            traceSlow: [2, 132, 199],
            traceMid: [180, 83, 9],
            traceFast: [217, 70, 239]
        }
    };

    var temaActual = 'dark';
    function color(nombre) {
        return TEMAS[temaActual][nombre] || '#ffffff';
    }

    function canvasFont(sizePx, weight, family) {
        var scale = (typeof document !== 'undefined' && document.documentElement && document.documentElement.classList.contains('presentation-mode')) ? 1.25 : 1.0;
        var finalSize = Math.max(9, Math.round(sizePx * scale));
        var wt = weight ? (weight + ' ') : '';
        var fam = family || '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
        return wt + finalSize + 'px ' + fam;
    }

    // =========================================================================
    // 2. ESTADO GLOBAL DE LA APLICACIÓN
    // =========================================================================
    var estado = {
        tabActiva: 'pendulo', // 'pendulo', 'gemelos', 'montecarlo', 'mapacaos'
        corriendo: true,
        velocidad: 1.0,
        mostrarEstela: true,
        algoritmo: 'rk4_38',
        modoLinealActivo: false,
        vistaEspacioFase: false,

        // Parámetros físicos
        params: PenduloFisica.crearParametros(1.0, 1.0, 1.0, 1.0, 9.81),

        // Estado dinámico canónico x = [theta1, theta2, omega1, omega2]
        x: [Math.PI / 4, Math.PI / 4, 0.0, 0.0],
        x0: [Math.PI / 4, Math.PI / 4, 0.0, 0.0],
        t: 0.0,
        dtBase: 0.005,

        // Energía
        e0: 0.0,
        deltaERel: 0.0,

        // Conteo continuo de vueltas
        theta2Continuo: Math.PI / 4,
        vueltasM2: 0,
        theta1Continuo: Math.PI / 4,
        vueltasM1: 0,

        // Historiales
        estelaM2: [],
        historialTiempo: [], // { t, th1, th2, w1, w2, th1Lin, th2Lin }
        ventanaS: 20.0,

        // Estado Lineal
        xLineal: [Math.PI / 4, Math.PI / 4, 0.0, 0.0],

        // Estado Gemelos
        deltaGemelos: 0.001,
        xGemelo: null,
        historialGemelos: [], // { t, d, logD }
        lambdaFit: 0.0,
        lambdaFitAdim: 0.0,

        // Estado Monte Carlo Dual
        mc: {
            corriendo: false,
            t: 0.0,
            tMax: 30.0,
            dt: 0.01,
            n: 300,
            sigma: 0.010,
            targetDerecho: 'alta', // 'alta' o 'actual'
            left: {
                x0: [0.6, 0.6, 0.0, 0.0],
                pendulos: [], // array de [th1, th2, w1, w2]
                flips: [],
                th2Continuo: []
            },
            right: {
                x0: [2.6, 2.6, 0.0, 0.0],
                pendulos: [],
                flips: [],
                th2Continuo: []
            },
            statsLeft: null,
            statsRight: null
        },

        // Mapa FTLE
        ftle: {
            gridObj: null,
            imgBitmap: null,
            offCanvas: null,
            inspectTh1: 0.0,
            inspectTh2: 0.0,
            inspectVal: 0.0,
            calculando: false,
            worker: null
        },

        // Arrastre interactivo en el péndulo
        arrastre: {
            activo: false,
            masa: 0 // 1 o 2
        }
    };

    // =========================================================================
    // 3. REFERENCIAS A ELEMENTOS DEL DOM
    // =========================================================================
    var dom = {};

    function initDomReferences() {
        dom.body = document.body;

        // Pestañas
        dom.tabBtns = document.querySelectorAll('.mode-tab-btn');
        dom.tabViews = {
            pendulo: document.getElementById('view-pendulo'),
            gemelos: document.getElementById('view-gemelos'),
            montecarlo: document.getElementById('view-montecarlo'),
            mapacaos: document.getElementById('view-mapacaos')
        };
        dom.controlsTab = {
            pendulo: document.getElementById('controls-tab-pendulo'),
            gemelos: document.getElementById('controls-tab-gemelos'),
            montecarlo: document.getElementById('controls-tab-montecarlo'),
            mapacaos: document.getElementById('controls-tab-mapacaos')
        };

        // Cabecera
        dom.btnThemeToggle = document.getElementById('btn-theme-toggle');
        dom.themeIcon = document.getElementById('theme-icon');
        dom.themeText = document.getElementById('theme-text');
        dom.btnToggleSidebar = document.getElementById('btn-toggle-sidebar');
        dom.btnPresentation = document.getElementById('btn-presentation');
        dom.sidebar = document.getElementById('sidebar-controls');

        // Transporte
        dom.btnPlayPause = document.getElementById('btn-play-pause');
        dom.btnReset = document.getElementById('btn-reset');
        dom.btnStep = document.getElementById('btn-step');

        // Presets
        dom.presetBtns = document.querySelectorAll('.btn-preset');

        // Sliders CI
        dom.sliderTh1 = document.getElementById('slider-th1');
        dom.sliderTh2 = document.getElementById('slider-th2');
        dom.lblTh1 = document.getElementById('lbl-th1');
        dom.lblTh2 = document.getElementById('lbl-th2');

        // Sliders Parámetros
        dom.sliderM1 = document.getElementById('slider-m1');
        dom.sliderM2 = document.getElementById('slider-m2');
        dom.sliderL1 = document.getElementById('slider-l1');
        dom.sliderL2 = document.getElementById('slider-l2');
        dom.lblM1 = document.getElementById('lbl-m1');
        dom.lblM2 = document.getElementById('lbl-m2');
        dom.lblL1 = document.getElementById('lbl-l1');
        dom.lblL2 = document.getElementById('lbl-l2');
        dom.lblMu = document.getElementById('lbl-mu');
        dom.lblLambda = document.getElementById('lbl-lambda');

        // Gravedad
        dom.gravityBtns = document.querySelectorAll('.btn-gravity');
        dom.sliderG = document.getElementById('slider-g');
        dom.lblG = document.getElementById('lbl-g');

        // Péndulo
        dom.canvasPenduloMain = document.getElementById('canvas-pendulo-main');
        dom.canvasPenduloChart = document.getElementById('canvas-pendulo-chart');
        dom.bigAngleTh1 = document.getElementById('big-angle-th1');
        dom.bigAngleTh2 = document.getElementById('big-angle-th2');
        dom.bigTime = document.getElementById('big-time');
        dom.linearWarning = document.getElementById('linear-range-warning');
        dom.checkLinearMode = document.getElementById('check-linear-mode');
        dom.btnTogglePhaseChart = document.getElementById('btn-toggle-phase-chart');
        dom.chartPenduloTitle = document.getElementById('chart-pendulo-title');
        dom.chartPenduloCaption = document.getElementById('chart-pendulo-caption');

        // Termómetro
        dom.thermometerFill = document.getElementById('thermometer-fill');
        dom.markM2 = document.getElementById('mark-m2');
        dom.markM1 = document.getElementById('mark-m1');
        dom.markAmbas = document.getElementById('mark-ambas');
        dom.markLblM2 = document.getElementById('mark-lbl-m2');
        dom.markLblM1 = document.getElementById('mark-lbl-m1');
        dom.markLblAmbas = document.getElementById('mark-lbl-ambas');
        dom.readoutTotalE = document.getElementById('readout-total-e');
        dom.readoutFlips = document.getElementById('readout-flips');

        // Gemelos
        dom.canvasGemelosStage = document.getElementById('canvas-gemelos-stage');
        dom.canvasGemelosSeparation = document.getElementById('canvas-gemelos-separation');
        dom.sliderDeltaGemelos = document.getElementById('slider-delta-gemelos');
        dom.lblDeltaGemelos = document.getElementById('lbl-delta-gemelos');
        dom.btnGemelosRegular = document.getElementById('btn-gemelos-regular');
        dom.btnGemelosCaotico = document.getElementById('btn-gemelos-caotico');
        dom.lblLambdaFit = document.getElementById('lbl-lambda-fit');
        dom.lblLambdaAdim = document.getElementById('lbl-lambda-adim');
        dom.gemelosBadge = document.getElementById('gemelos-status-badge');

        // Monte Carlo
        dom.canvasMcStageLeft = document.getElementById('canvas-mc-stage-left');
        dom.canvasMcHistLeft = document.getElementById('canvas-mc-hist-left');
        dom.canvasMcStageRight = document.getElementById('canvas-mc-stage-right');
        dom.canvasMcHistRight = document.getElementById('canvas-mc-hist-right');
        dom.sliderMcN = document.getElementById('slider-mc-n');
        dom.sliderMcSigma = document.getElementById('slider-mc-sigma');
        dom.sliderMcT = document.getElementById('slider-mc-t');
        dom.lblMcN = document.getElementById('lbl-mc-n');
        dom.lblMcSigma = document.getElementById('lbl-mc-sigma');
        dom.lblMcT = document.getElementById('lbl-mc-t');
        dom.selectMcRightTarget = document.getElementById('select-mc-right-target');
        dom.btnRunMontecarlo = document.getElementById('btn-run-montecarlo');
        dom.mcRightTitle = document.getElementById('mc-right-title');
        dom.mcRightSub = document.getElementById('mc-right-sub');

        dom.mcLeftMean = document.getElementById('mc-left-mean');
        dom.mcLeftCiMean = document.getElementById('mc-left-ci-mean');
        dom.mcLeftStd = document.getElementById('mc-left-std');
        dom.mcLeftCiStd = document.getElementById('mc-left-ci-std');
        dom.mcLeftPflip = document.getElementById('mc-left-pflip');
        dom.mcLeftCiWilson = document.getElementById('mc-left-ci-wilson');

        dom.mcRightMean = document.getElementById('mc-right-mean');
        dom.mcRightCiMean = document.getElementById('mc-right-ci-mean');
        dom.mcRightStd = document.getElementById('mc-right-std');
        dom.mcRightCiStd = document.getElementById('mc-right-ci-std');
        dom.mcRightPflip = document.getElementById('mc-right-pflip');
        dom.mcRightCiWilson = document.getElementById('mc-right-ci-wilson');

        // Mapa de Caos
        dom.canvasMapaCaos = document.getElementById('canvas-mapacaos-main');
        dom.mapacaosBox = document.getElementById('mapacaos-box');
        dom.mapaInspectTh1 = document.getElementById('mapa-inspect-th1');
        dom.mapaInspectTh2 = document.getElementById('mapa-inspect-th2');
        dom.mapaInspectLambda = document.getElementById('mapa-inspect-lambda');
        dom.btnVerEnPendulo = document.getElementById('btn-ver-en-pendulo');

        // Avanzado
        dom.selectIntegrator = document.getElementById('select-integrator');
        dom.selectSpeed = document.getElementById('select-speed');
        dom.checkTrace = document.getElementById('check-trace');
        dom.valErrorDeltaE = document.getElementById('val-error-delta-e');
    }

    // =========================================================================
    // 4. INICIALIZACIÓN Y CONFIGURACIÓN FÍSICA
    // =========================================================================
    function reiniciarSimulacion() {
        estado.corriendo = false;
        dom.btnPlayPause.textContent = 'Reanudar';
        dom.btnPlayPause.classList.remove('btn-primary');
        inicializarSimulacion();
        renderizarVistaActual();
    }

    function inicializarSimulacion() {
        estado.params = PenduloFisica.crearParametros(
            parseFloat(dom.sliderM1.value),
            parseFloat(dom.sliderM2.value),
            parseFloat(dom.sliderL1.value),
            parseFloat(dom.sliderL2.value),
            parseFloat(dom.sliderG.value)
        );
        estado.linInfo = PenduloFisica.linealizacion(estado.params);

        var th1 = parseFloat(dom.sliderTh1.value) * (Math.PI / 180.0);
        var th2 = parseFloat(dom.sliderTh2.value) * (Math.PI / 180.0);
        var w1 = 0.0;
        var w2 = 0.0;

        estado.x0 = [th1, th2, w1, w2];
        estado.x = [th1, th2, w1, w2];
        estado.xLineal = [th1, th2, w1, w2];
        estado.t = 0.0;
        estado.theta1Continuo = th1;
        estado.theta2Continuo = th2;
        estado.vueltasM1 = 0;
        estado.vueltasM2 = 0;

        estado.e0 = PenduloFisica.energiaTotal(estado.x, estado.params);
        estado.deltaERel = 0.0;

        estado.estelaM2 = [];
        estado.historialTiempo = [];

        // Inicializar Gemelo
        var dVal = parseFloat(dom.sliderDeltaGemelos.value);
        estado.deltaGemelos = dVal;
        estado.xGemelo = [th1 + dVal, th2 + dVal, w1, w2];
        estado.historialGemelos = [];
        estado.lambdaFit = 0.0;
        estado.lambdaFitAdim = 0.0;

        registrarHistorialPunto();
        actualizarTermometroEnergia();
        actualizarIndicadoresHUD();
        actualizarParametrosLabels();

        // Si la pestaña actual es Monte Carlo, reiniciar también Monte Carlo
        if (estado.tabActiva === 'montecarlo') {
            reiniciarMonteCarlo();
        }

        // Si el mapa ya tiene retícula, redibujar la cruz
        redibujarMapaCaos();
    }

    function registrarHistorialPunto() {
        var x = estado.x;
        var p = estado.params;
        var xLin = estado.xLineal;

        // Posiciones cartesianas masa 2 para estela
        var pos = PenduloFisica.posicionesCartesianas(x, p);
        var v2Sq = x[3] * x[3] * (p.l2 * p.l2);
        estado.estelaM2.push({ x: pos.x2, y: pos.y2, v2Sq: v2Sq });
        if (estado.estelaM2.length > 1200) {
            estado.estelaM2.shift();
        }

        // Historial tiempo (ventana rodante de 20s)
        estado.historialTiempo.push({
            t: estado.t,
            th1: x[0],
            th2: x[1],
            w1: x[2],
            w2: x[3],
            th1Lin: xLin[0],
            th2Lin: xLin[1]
        });

        var tMin = estado.t - estado.ventanaS;
        while (estado.historialTiempo.length > 2 && estado.historialTiempo[0].t < tMin) {
            estado.historialTiempo.shift();
        }

        // Historial Gemelos
        if (estado.xGemelo) {
            var diff0 = estado.x[0] - estado.xGemelo[0];
            var diff1 = estado.x[1] - estado.xGemelo[1];
            var diff2 = estado.x[2] - estado.xGemelo[2];
            var diff3 = estado.x[3] - estado.xGemelo[3];
            var d = Math.sqrt(diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3);
            if (d < 1e-12) d = 1e-12;
            var logD = Math.log10(d);

            estado.historialGemelos.push({
                t: estado.t,
                d: d,
                logD: logD
            });

            if (estado.historialGemelos.length > 3000) {
                estado.historialGemelos.shift();
            }

            calcularAjusteLyapunov();
        }
    }

    function calcularAjusteLyapunov() {
        var h = estado.historialGemelos;
        if (h.length < 50) return;

        // Buscamos la fase inicial de crecimiento exponencial antes de saturar (d < 2.0)
        var puntos = [];
        var t0 = h[0].t;
        for (var i = 0; i < h.length; i++) {
            if (h[i].d > 2.5) break; // Saturó a escala del atractor
            puntos.push({ t: h[i].t - t0, logE: Math.log(h[i].d) });
        }

        if (puntos.length < 30) {
            estado.lambdaFit = 0.0;
            estado.lambdaFitAdim = 0.0;
            return;
        }

        // Regresión lineal simple: log(d) = a + lambda * t
        var sumT = 0, sumL = 0, sumT2 = 0, sumTL = 0;
        var n = puntos.length;
        for (var j = 0; j < n; j++) {
            sumT += puntos[j].t;
            sumL += puntos[j].logE;
            sumT2 += puntos[j].t * puntos[j].t;
            sumTL += puntos[j].t * puntos[j].logE;
        }
        var denom = (n * sumT2 - sumT * sumT);
        if (Math.abs(denom) > 1e-12) {
            var lambda = (n * sumTL - sumT * sumL) / denom;
            if (lambda < 0) lambda = 0;
            estado.lambdaFit = lambda;
            var escalaTiempo = Math.sqrt(estado.params.l1 / estado.params.g);
            estado.lambdaFitAdim = lambda * Math.sqrt(estado.params.l1 / estado.params.g);
        }

        // Actualizar badges
        dom.lblLambdaFit.textContent = 'λ ≈ ' + estado.lambdaFit.toFixed(2) + ' s⁻¹';
        dom.lblLambdaAdim.textContent = '(adimensional: λ·√(l₁/g) ≈ ' + estado.lambdaFitAdim.toFixed(2) + ')';

        if (estado.lambdaFit > 0.35) {
            dom.gemelosBadge.textContent = 'Caótico (separación exponencial)';
            dom.gemelosBadge.style.background = 'rgba(247, 37, 133, 0.2)';
            dom.gemelosBadge.style.color = '#f72585';
            dom.gemelosBadge.style.borderColor = 'rgba(247, 37, 133, 0.4)';
        } else {
            dom.gemelosBadge.textContent = 'Regular (separación lenta/acotada)';
            dom.gemelosBadge.style.background = 'rgba(76, 201, 240, 0.2)';
            dom.gemelosBadge.style.color = '#4cc9f0';
            dom.gemelosBadge.style.borderColor = 'rgba(76, 201, 240, 0.4)';
        }
    }

    // =========================================================================
    // 5. PASO DE SIMULACIÓN NUMÉRICA
    // =========================================================================
    function pasoSimulacion(dtEfectivo) {
        var pasoFunc = (estado.algoritmo === 'rk4_clasico')
            ? PenduloFisica.pasoRK4Clasico
            : PenduloFisica.pasoRK4;

        // Integrar péndulo principal
        var xAnt = estado.x;
        estado.x = pasoFunc(estado.x, dtEfectivo, estado.params);
        estado.t += dtEfectivo;

        // Vueltas continuas y conteo acumulativo de cruces (ROADMAP.md y flips.py)
        var dTh2 = estado.x[1] - xAnt[1];
        var th2Ant = estado.theta2Continuo;
        estado.theta2Continuo += dTh2;
        estado.vueltasM2 = PenduloFisica.contarCruces(th2Ant, estado.theta2Continuo, estado.vueltasM2);

        var dTh1 = estado.x[0] - xAnt[0];
        var th1Ant = estado.theta1Continuo;
        estado.theta1Continuo += dTh1;
        estado.vueltasM1 = PenduloFisica.contarCruces(th1Ant, estado.theta1Continuo, estado.vueltasM1);

        // Integrar modelo lineal si está activo o en historial
        if (!estado.linInfo) {
            estado.linInfo = PenduloFisica.linealizacion(estado.params);
        }
        estado.xLineal = PenduloFisica.pasoLinealRK4(estado.xLineal, dtEfectivo, estado.linInfo);

        // Integrar gemelo
        if (estado.xGemelo) {
            estado.xGemelo = pasoFunc(estado.xGemelo, dtEfectivo, estado.params);
        }

        // Conservación de energía
        var eActual = PenduloFisica.energiaTotal(estado.x, estado.params);
        if (Math.abs(estado.e0) > 1e-12) {
            estado.deltaERel = Math.abs(eActual - estado.e0) / Math.abs(estado.e0);
        } else {
            estado.deltaERel = Math.abs(eActual - estado.e0);
        }
    }

    // =========================================================================
    // 6. MONTE CARLO DUAL LADO A LADO
    // =========================================================================
    function reiniciarMonteCarlo() {
        var n = parseInt(dom.sliderMcN.value, 10);
        var sigma = parseFloat(dom.sliderMcSigma.value);
        var tMax = parseFloat(dom.sliderMcT.value);
        var targetDerecho = dom.selectMcRightTarget.value;

        estado.mc.n = n;
        estado.mc.sigma = sigma;
        estado.mc.tMax = tMax;
        estado.mc.t = 0.0;
        estado.mc.targetDerecho = targetDerecho;
        estado.mc.corriendo = true;

        // Panel Izquierdo: Siempre Baja Energía (0.6, 0.6)
        estado.mc.left.x0 = [0.6, 0.6, 0.0, 0.0];
        estado.mc.left.pendulos = [];
        estado.mc.left.puntosIniciales = [];
        estado.mc.left.flips = new Array(n).fill(0);
        estado.mc.left.th2Continuo = new Array(n).fill(0.6);

        // Panel Derecho: Alta Energía (2.6, 2.6) o Condición Actual
        if (targetDerecho === 'actual') {
            estado.mc.right.x0 = [estado.x0[0], estado.x0[1], 0.0, 0.0];
            dom.mcRightTitle.textContent = 'Condición actual (' + (estado.x0[0] * 180 / Math.PI).toFixed(0) + '°, ' + (estado.x0[1] * 180 / Math.PI).toFixed(0) + '°)';
            dom.mcRightSub.textContent = 'E = ' + estado.e0.toFixed(1) + ' J';
        } else {
            estado.mc.right.x0 = [2.6, 2.6, 0.0, 0.0];
            dom.mcRightTitle.textContent = 'Alta energía (2.6 rad, 2.6 rad)';
            dom.mcRightSub.textContent = 'E > Ec_m2 (régimen caótico)';
        }
        estado.mc.right.pendulos = [];
        estado.mc.right.puntosIniciales = [];
        estado.mc.right.flips = new Array(n).fill(0);
        estado.mc.right.th2Continuo = new Array(n).fill(estado.mc.right.x0[1]);

        // Generar partículas con mulberry32 determinista
        var rngL = PenduloFisica.mulberry32(12345);
        for (var i = 0; i < n; i++) {
            // Box-Muller para distribución normal N(0, sigma)
            var u1 = Math.max(1e-12, rngL());
            var u2 = rngL();
            var z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
            var z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);

            var th1L = estado.mc.left.x0[0] + sigma * z0;
            var th2L = estado.mc.left.x0[1] + sigma * z1;
            estado.mc.left.pendulos.push([th1L, th2L, 0.0, 0.0]);
            estado.mc.left.puntosIniciales.push([th1L * 180 / Math.PI, th2L * 180 / Math.PI]);
            estado.mc.left.th2Continuo[i] = th2L;
        }

        var rngR = PenduloFisica.mulberry32(67890);
        for (var j = 0; j < n; j++) {
            var u1r = Math.max(1e-12, rngR());
            var u2r = rngR();
            var z0r = Math.sqrt(-2.0 * Math.log(u1r)) * Math.cos(2.0 * Math.PI * u2r);
            var z1r = Math.sqrt(-2.0 * Math.log(u1r)) * Math.sin(2.0 * Math.PI * u2r);

            var th1R = estado.mc.right.x0[0] + sigma * z0r;
            var th2R = estado.mc.right.x0[1] + sigma * z1r;
            estado.mc.right.pendulos.push([th1R, th2R, 0.0, 0.0]);
            estado.mc.right.puntosIniciales.push([th1R * 180 / Math.PI, th2R * 180 / Math.PI]);
            estado.mc.right.th2Continuo[j] = th2R;
        }

        actualizarEstadisticasMonteCarlo();
    }

    function pasoMonteCarlo(dtTotal) {
        if (!estado.mc.corriendo) return;

        var subDt = 0.005;
        var subPasos = Math.max(1, Math.min(20, Math.round(dtTotal / subDt)));
        var dtStep = dtTotal / subPasos;
        var n = estado.mc.n;
        var p = estado.params;
        var pasoFunc = PenduloFisica.pasoRK4;

        for (var s = 0; s < subPasos; s++) {
            // Avanzar izquierda
            for (var i = 0; i < n; i++) {
                var xAntL = estado.mc.left.pendulos[i];
                var xNextL = pasoFunc(xAntL, dtStep, p);
                estado.mc.left.pendulos[i] = xNextL;

                var thAntL = estado.mc.left.th2Continuo[i];
                var thNextL = thAntL + (xNextL[1] - xAntL[1]);
                estado.mc.left.flips[i] = PenduloFisica.contarCruces(thAntL, thNextL, estado.mc.left.flips[i]);
                estado.mc.left.th2Continuo[i] = thNextL;
            }

            // Avanzar derecha
            for (var j = 0; j < n; j++) {
                var xAntR = estado.mc.right.pendulos[j];
                var xNextR = pasoFunc(xAntR, dtStep, p);
                estado.mc.right.pendulos[j] = xNextR;

                var thAntR = estado.mc.right.th2Continuo[j];
                var thNextR = thAntR + (xNextR[1] - xAntR[1]);
                estado.mc.right.flips[j] = PenduloFisica.contarCruces(thAntR, thNextR, estado.mc.right.flips[j]);
                estado.mc.right.th2Continuo[j] = thNextR;
            }

            estado.mc.t += dtStep;
            if (estado.mc.t >= estado.mc.tMax) {
                estado.mc.corriendo = false;
                break;
            }
        }

        actualizarEstadisticasMonteCarlo();
    }

    function actualizarEstadisticasMonteCarlo() {
        var n = estado.mc.n;
        var flipsL = estado.mc.left.flips;
        var flipsR = estado.mc.right.flips;

        // Estadísticos Panel Izquierdo
        var bootL = PenduloFisica.bootstrapPercentil(flipsL, 1000, 42, 0.95);
        var kFlipsL = 0;
        for (var i = 0; i < n; i++) {
            if (flipsL[i] > 0) kFlipsL++;
        }
        var wilsonL = PenduloFisica.intervaloWilson(kFlipsL, n);

        dom.mcLeftMean.textContent = bootL.media.toFixed(2);
        dom.mcLeftCiMean.textContent = 'IC 95%: [' + bootL.media_ci_bajo.toFixed(2) + ', ' + bootL.media_ci_alto.toFixed(2) + ']';
        dom.mcLeftStd.textContent = bootL.std.toFixed(2);
        dom.mcLeftCiStd.textContent = 'IC 95%: [' + bootL.std_ci_bajo.toFixed(2) + ', ' + bootL.std_ci_alto.toFixed(2) + ']';
        dom.mcLeftPflip.textContent = (wilsonL.p * 100).toFixed(1) + '%';
        dom.mcLeftCiWilson.textContent = 'Wilson: [' + (wilsonL.bajo * 100).toFixed(1) + '%, ' + (wilsonL.alto * 100).toFixed(1) + '%]';

        // Estadísticos Panel Derecho
        var bootR = PenduloFisica.bootstrapPercentil(flipsR, 1000, 43, 0.95);
        var kFlipsR = 0;
        for (var j = 0; j < n; j++) {
            if (flipsR[j] > 0) kFlipsR++;
        }
        var wilsonR = PenduloFisica.intervaloWilson(kFlipsR, n);

        dom.mcRightMean.textContent = bootR.media.toFixed(2);
        dom.mcRightCiMean.textContent = 'IC 95%: [' + bootR.media_ci_bajo.toFixed(2) + ', ' + bootR.media_ci_alto.toFixed(2) + ']';
        dom.mcRightStd.textContent = bootR.std.toFixed(2);
        dom.mcRightCiStd.textContent = 'IC 95%: [' + bootR.std_ci_bajo.toFixed(2) + ', ' + bootR.std_ci_alto.toFixed(2) + ']';
        dom.mcRightPflip.textContent = (wilsonR.p * 100).toFixed(1) + '%';
        dom.mcRightCiWilson.textContent = 'Wilson: [' + (wilsonR.bajo * 100).toFixed(1) + '%, ' + (wilsonR.alto * 100).toFixed(1) + '%]';

        estado.mc.statsLeft = { boot: bootL, wilson: wilsonL };
        estado.mc.statsRight = { boot: bootR, wilson: wilsonR };
    }

    // =========================================================================
    // 7. MAPA DE CAOS (FTLE ADIMENSIONAL)
    // =========================================================================
    function inicializarMapaCaos() {
        if (window.FTLE_DEFAULT_GRID) {
            estado.ftle.gridObj = window.FTLE_DEFAULT_GRID;
            generarBitmapMapaCaos();
        }
    }

    var INFERNO_FLAT = new Uint8Array([0, 0, 4, 1, 0, 5, 1, 1, 6, 1, 1, 8, 2, 1, 10, 2, 2, 12, 2, 2, 14, 3, 2, 16, 4, 3, 18, 4, 3, 20, 5, 4, 23, 6, 4, 25, 7, 5, 27, 8, 5, 29, 9, 6, 31, 10, 7, 34, 11, 7, 36, 12, 8, 38, 13, 8, 41, 14, 9, 43, 16, 9, 45, 17, 10, 48, 18, 10, 50, 20, 11, 52, 21, 11, 55, 22, 11, 57, 24, 12, 60, 25, 12, 62, 27, 12, 65, 28, 12, 67, 30, 12, 69, 31, 12, 72, 33, 12, 74, 35, 12, 76, 36, 12, 79, 38, 12, 81, 40, 11, 83, 41, 11, 85, 43, 11, 87, 45, 11, 89, 47, 10, 91, 49, 10, 92, 50, 10, 94, 52, 10, 95, 54, 9, 97, 56, 9, 98, 57, 9, 99, 59, 9, 100, 61, 9, 101, 62, 9, 102, 64, 10, 103, 66, 10, 104, 68, 10, 104, 69, 10, 105, 71, 11, 106, 73, 11, 106, 74, 12, 107, 76, 12, 107, 77, 13, 108, 79, 13, 108, 81, 14, 108, 82, 14, 109, 84, 15, 109, 85, 15, 109, 87, 16, 110, 89, 16, 110, 90, 17, 110, 92, 18, 110, 93, 18, 110, 95, 19, 110, 97, 19, 110, 98, 20, 110, 100, 21, 110, 101, 21, 110, 103, 22, 110, 105, 22, 110, 106, 23, 110, 108, 24, 110, 109, 24, 110, 111, 25, 110, 113, 25, 110, 114, 26, 110, 116, 26, 110, 117, 27, 110, 119, 28, 109, 120, 28, 109, 122, 29, 109, 124, 29, 109, 125, 30, 109, 127, 30, 108, 128, 31, 108, 130, 32, 108, 132, 32, 107, 133, 33, 107, 135, 33, 107, 136, 34, 106, 138, 34, 106, 140, 35, 105, 141, 35, 105, 143, 36, 105, 144, 37, 104, 146, 37, 104, 147, 38, 103, 149, 38, 103, 151, 39, 102, 152, 39, 102, 154, 40, 101, 155, 41, 100, 157, 41, 100, 159, 42, 99, 160, 42, 99, 162, 43, 98, 163, 44, 97, 165, 44, 96, 166, 45, 96, 168, 46, 95, 169, 46, 94, 171, 47, 94, 173, 48, 93, 174, 48, 92, 176, 49, 91, 177, 50, 90, 179, 50, 90, 180, 51, 89, 182, 52, 88, 183, 53, 87, 185, 53, 86, 186, 54, 85, 188, 55, 84, 189, 56, 83, 191, 57, 82, 192, 58, 81, 193, 58, 80, 195, 59, 79, 196, 60, 78, 198, 61, 77, 199, 62, 76, 200, 63, 75, 202, 64, 74, 203, 65, 73, 204, 66, 72, 206, 67, 71, 207, 68, 70, 208, 69, 69, 210, 70, 68, 211, 71, 67, 212, 72, 66, 213, 74, 65, 215, 75, 63, 216, 76, 62, 217, 77, 61, 218, 78, 60, 219, 80, 59, 221, 81, 58, 222, 82, 56, 223, 83, 55, 224, 85, 54, 225, 86, 53, 226, 87, 52, 227, 89, 51, 228, 90, 49, 229, 92, 48, 230, 93, 47, 231, 94, 46, 232, 96, 45, 233, 97, 43, 234, 99, 42, 235, 100, 41, 235, 102, 40, 236, 103, 38, 237, 105, 37, 238, 106, 36, 239, 108, 35, 239, 110, 33, 240, 111, 32, 241, 113, 31, 241, 115, 29, 242, 116, 28, 243, 118, 27, 243, 120, 25, 244, 121, 24, 245, 123, 23, 245, 125, 21, 246, 126, 20, 246, 128, 19, 247, 130, 18, 247, 132, 16, 248, 133, 15, 248, 135, 14, 248, 137, 12, 249, 139, 11, 249, 140, 10, 249, 142, 9, 250, 144, 8, 250, 146, 7, 250, 148, 7, 251, 150, 6, 251, 151, 6, 251, 153, 6, 251, 155, 6, 251, 157, 7, 252, 159, 7, 252, 161, 8, 252, 163, 9, 252, 165, 10, 252, 166, 12, 252, 168, 13, 252, 170, 15, 252, 172, 17, 252, 174, 18, 252, 176, 20, 252, 178, 22, 252, 180, 24, 251, 182, 26, 251, 184, 29, 251, 186, 31, 251, 188, 33, 251, 190, 35, 250, 192, 38, 250, 194, 40, 250, 196, 42, 250, 198, 45, 249, 199, 47, 249, 201, 50, 249, 203, 53, 248, 205, 55, 248, 207, 58, 247, 209, 61, 247, 211, 64, 246, 213, 67, 246, 215, 70, 245, 217, 73, 245, 219, 76, 244, 221, 79, 244, 223, 83, 244, 225, 86, 243, 227, 90, 243, 229, 93, 242, 230, 97, 242, 232, 101, 242, 234, 105, 241, 236, 109, 241, 237, 113, 241, 239, 117, 241, 241, 121, 242, 242, 125, 242, 244, 130, 243, 245, 134, 243, 246, 138, 244, 248, 142, 245, 249, 146, 246, 250, 150, 248, 251, 154, 249, 252, 157, 250, 253, 161, 252, 255, 164]);

    function infernoColormap(t) {
        // Colormap exacto inferno de matplotlib (256 colores)
        t = Math.max(0.0, Math.min(1.0, t));
        var idx = Math.floor(t * 255);
        if (idx > 255) idx = 255;
        var off = idx * 3;
        return [INFERNO_FLAT[off], INFERNO_FLAT[off + 1], INFERNO_FLAT[off + 2]];
    }

    function generarBitmapMapaCaos() {
        if (!estado.ftle.gridObj) return;
        var g = estado.ftle.gridObj;
        var res = g.res;
        var data = g.grid;

        var offCanvas = document.createElement('canvas');
        offCanvas.width = res;
        offCanvas.height = res;
        var ctx = offCanvas.getContext('2d');
        var imgData = ctx.createImageData(res, res);

        // Encontrar max para normalizar (saturar en 3.5 para contraste óptimo)
        var maxVal = 3.5;

        for (var r = 0; r < res; r++) {
            for (var c = 0; c < res; c++) {
                // El grid está en orden row (theta1 de -pi a pi) y col (theta2 de -pi a pi)
                // Para que el eje Y apunte hacia arriba (theta2 en Y), invertimos r
                var idxData = r * res + c;
                var val = data[idxData];
                var t = Math.min(1.0, val / maxVal);
                var rgb = infernoColormap(t);

                // En canvas Y=0 es arriba, así que Y=(res - 1 - r) o r directo
                var pixelIdx = ((res - 1 - r) * res + c) * 4;
                imgData.data[pixelIdx + 0] = rgb[0];
                imgData.data[pixelIdx + 1] = rgb[1];
                imgData.data[pixelIdx + 2] = rgb[2];
                imgData.data[pixelIdx + 3] = 255;
            }
        }

        ctx.putImageData(imgData, 0, 0);
        estado.ftle.offCanvas = offCanvas;
        redibujarMapaCaos();
    }

    function redibujarMapaCaos() {
        var canvas = dom.canvasMapaCaos;
        if (!canvas) return;
        var rect = canvas.getBoundingClientRect();
        if (rect.width <= 0 || rect.height <= 0) return;

        var dpr = window.devicePixelRatio || 1;
        var w = Math.floor(rect.width * dpr);
        var h = Math.floor(rect.height * dpr);
        if (canvas.width !== w || canvas.height !== h) {
            canvas.width = w;
            canvas.height = h;
        }

        var ctx = canvas.getContext('2d');
        ctx.save();
        ctx.scale(dpr, dpr);

        var cw = rect.width;
        var ch = rect.height;

        // Limpiar el lienzo SIEMPRE
        ctx.fillStyle = color('bg');
        ctx.fillRect(0, 0, cw, ch);

        // Márgenes para ejes y colorbar
        var padLeft = 46;
        var padBottom = 38;
        var padTop = 16;
        var padRight = 50;
        var plotW = cw - padLeft - padRight;
        var plotH = ch - padTop - padBottom;

        if (plotW <= 10 || plotH <= 10) {
            ctx.restore();
            return;
        }

        // 1. Dibujar el mapa si está listo
        if (estado.ftle.offCanvas) {
            ctx.imageSmoothingEnabled = true;
            ctx.drawImage(estado.ftle.offCanvas, padLeft, padTop, plotW, plotH);
        } else {
            ctx.fillStyle = color('panelSub');
            ctx.fillRect(padLeft, padTop, plotW, plotH);
            ctx.fillStyle = color('textMuted');
            ctx.font = canvasFont(13);
            ctx.textAlign = 'center';
            ctx.fillText('Calculando mapa FTLE...', padLeft + plotW / 2, padTop + plotH / 2);
        }

        // 2. Dibujar curvas de contorno analíticas (Ec_m2 y Ec_m1)
        dibujarContornosEnergia(ctx, padLeft, padTop, plotW, plotH);

        // 3. Ejes y etiquetas (-π a π)
        ctx.strokeStyle = color('border');
        ctx.lineWidth = 1.5;
        ctx.strokeRect(padLeft, padTop, plotW, plotH);

        ctx.fillStyle = color('textMuted');
        ctx.font = canvasFont(11);
        ctx.textAlign = 'center';

        var ticks = [
            { val: -Math.PI, txt: '-π' },
            { val: -Math.PI / 2, txt: '-π/2' },
            { val: 0, txt: '0' },
            { val: Math.PI / 2, txt: 'π/2' },
            { val: Math.PI, txt: 'π' }
        ];

        // Eje X: θ₁
        for (var i = 0; i < ticks.length; i++) {
            var tx = padLeft + ((ticks[i].val + Math.PI) / (2 * Math.PI)) * plotW;
            ctx.beginPath();
            ctx.moveTo(tx, padTop + plotH);
            ctx.lineTo(tx, padTop + plotH + 5);
            ctx.stroke();
            ctx.fillText(ticks[i].txt, tx, padTop + plotH + 16);
        }
        ctx.fillStyle = color('text');
        ctx.font = canvasFont(12);
        ctx.fillText('Ángulo inicial θ₁(0) [rad]', padLeft + plotW / 2, padTop + plotH + 32);

        // Eje Y: θ₂
        ctx.fillStyle = color('textMuted');
        ctx.font = canvasFont(11);
        ctx.textAlign = 'right';
        for (var j = 0; j < ticks.length; j++) {
            var ty = padTop + plotH - ((ticks[j].val + Math.PI) / (2 * Math.PI)) * plotH;
            ctx.beginPath();
            ctx.moveTo(padLeft, ty);
            ctx.lineTo(padLeft - 5, ty);
            ctx.stroke();
            ctx.fillText(ticks[j].txt, padLeft - 8, ty + 4);
        }
        ctx.save();
        ctx.translate(padLeft - 32, padTop + plotH / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillStyle = color('text');
        ctx.font = canvasFont(12);
        ctx.textAlign = 'center';
        ctx.fillText('Ángulo inicial θ₂(0) [rad]', 0, 0);
        ctx.restore();

        // 4. Colorbar vertical a la derecha
        var cbX = padLeft + plotW + 12;
        var cbW = 12;
        var cbH = plotH;
        for (var py = 0; py < cbH; py++) {
            var normVal = 1.0 - (py / cbH);
            var rgb = infernoColormap(normVal);
            ctx.fillStyle = 'rgb(' + rgb[0] + ',' + rgb[1] + ',' + rgb[2] + ')';
            ctx.fillRect(cbX, padTop + py, cbW, 1);
        }
        ctx.strokeStyle = color('border');
        ctx.lineWidth = 1;
        ctx.strokeRect(cbX, padTop, cbW, cbH);

        ctx.fillStyle = color('textMuted');
        ctx.font = canvasFont(10);
        ctx.textAlign = 'left';
        ctx.fillText('≥3.5', cbX + cbW + 4, padTop + 8);
        ctx.fillText('1.7', cbX + cbW + 4, padTop + cbH / 2 + 4);
        ctx.fillText('0.0', cbX + cbW + 4, padTop + cbH);

        // 5. EXACTAMENTE UNA SOLA CRUZ en la condición inicial actual (x0[0], x0[1])
        var th1Act = estado.x0[0];
        var th2Act = estado.x0[1];
        // Mapear a [-pi, pi]
        th1Act = Math.atan2(Math.sin(th1Act), Math.cos(th1Act));
        th2Act = Math.atan2(Math.sin(th2Act), Math.cos(th2Act));

        var crossX = padLeft + ((th1Act + Math.PI) / (2 * Math.PI)) * plotW;
        var crossY = padTop + plotH - ((th2Act + Math.PI) / (2 * Math.PI)) * plotH;

        if (crossX >= padLeft && crossX <= padLeft + plotW && crossY >= padTop && crossY <= padTop + plotH) {
            ctx.lineWidth = 3;
            ctx.strokeStyle = '#000000';
            ctx.beginPath();
            ctx.moveTo(crossX - 9, crossY); ctx.lineTo(crossX + 9, crossY);
            ctx.moveTo(crossX, crossY - 9); ctx.lineTo(crossX, crossY + 9);
            ctx.stroke();

            ctx.lineWidth = 1.8;
            ctx.strokeStyle = '#ffffff';
            ctx.beginPath();
            ctx.moveTo(crossX - 8, crossY); ctx.lineTo(crossX + 8, crossY);
            ctx.moveTo(crossX, crossY - 8); ctx.lineTo(crossX, crossY + 8);
            ctx.stroke();

            // Círculo concéntrico sutil
            ctx.beginPath();
            ctx.arc(crossX, crossY, 5, 0, 2 * Math.PI);
            ctx.strokeStyle = '#000000';
            ctx.lineWidth = 2.5;
            ctx.stroke();
            ctx.strokeStyle = '#ffff00';
            ctx.lineWidth = 1.5;
            ctx.stroke();
        }

        // 6. Leyenda interna dentro del mapa
        var legW = 155;
        var legH = 46;
        var legX = padLeft + plotW - legW - 8;
        var legY = padTop + 8;
        ctx.fillStyle = (temaActual === 'dark') ? 'rgba(18, 24, 38, 0.88)' : 'rgba(255, 255, 255, 0.88)';
        ctx.strokeStyle = color('border');
        ctx.lineWidth = 1;
        ctx.fillRect(legX, legY, legW, legH);
        ctx.strokeRect(legX, legY, legW, legH);

        // Ec,m2 (Cian sólido 2px)
        ctx.strokeStyle = color('contourM2');
        ctx.lineWidth = 2.0;
        ctx.beginPath();
        ctx.moveTo(legX + 8, legY + 14);
        ctx.lineTo(legX + 30, legY + 14);
        ctx.stroke();

        ctx.fillStyle = color('text');
        ctx.font = canvasFont(11);
        ctx.textAlign = 'left';
        ctx.fillText('E = Ec,m2 (vuelta m₂)', legX + 36, legY + 18);

        // Ec,m1 (Dashed 1.6px)
        ctx.strokeStyle = color('contourM1');
        ctx.lineWidth = 1.6;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(legX + 8, legY + 32);
        ctx.lineTo(legX + 30, legY + 32);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillText('E = Ec,m1 (vuelta m₁)', legX + 36, legY + 36);

        ctx.restore();
    }

    function dibujarContornosEnergia(ctx, x0, y0, w, h) {
        var p = estado.params;
        var umb = PenduloFisica.umbralesEnergia(p);
        var ec2 = umb.Ec_m2;
        var ec1 = umb.Ec_m1;

        // V(th1, th2) = -(m1+m2)g*l1*cos(th1) - m2*g*l2*cos(th2)
        // V(th1, th2) = Ec => curva analítica
        // cos(th2) = (Ec + (m1+m2)*g*l1*cos(th1)) / (-m2*g*l2)
        function dibujarCurvaUmbral(ec, estiloColor, esDashed, ancho) {
            ctx.save();
            ctx.strokeStyle = estiloColor;
            ctx.lineWidth = ancho || 1.6;
            if (esDashed) ctx.setLineDash([4, 4]);

            var denom = -p.m2 * p.g * p.l2;
            var numBase = ec;
            var numTh1Factor = (p.m1 + p.m2) * p.g * p.l1;

            var nSteps = 180;
            // Rama positiva y negativa de acos
            for (var sign = -1; sign <= 1; sign += 2) {
                ctx.beginPath();
                var iniciado = false;
                for (var s = 0; s <= nSteps; s++) {
                    var th1 = -Math.PI + (s / nSteps) * 2 * Math.PI;
                    var cosTh2 = (numBase + numTh1Factor * Math.cos(th1)) / denom;
                    if (cosTh2 >= -1.0 && cosTh2 <= 1.0) {
                        var th2 = sign * Math.acos(cosTh2);
                        var px = x0 + ((th1 + Math.PI) / (2 * Math.PI)) * w;
                        var py = y0 + h - ((th2 + Math.PI) / (2 * Math.PI)) * h;
                        if (!iniciado) {
                            ctx.moveTo(px, py);
                            iniciado = true;
                        } else {
                            ctx.lineTo(px, py);
                        }
                    } else {
                        iniciado = false;
                    }
                }
                ctx.stroke();
            }
            ctx.restore();
        }

        // Contorno Ec_m2 (Cian sólido 2.0 px)
        dibujarCurvaUmbral(ec2, color('contourM2'), false, 2.0);
        // Contorno Ec_m1 (Ámbar punteado 1.6 px)
        dibujarCurvaUmbral(ec1, color('contourM1'), true, 1.6);
    }

    // =========================================================================
    // 8. RENDERIZADO VISUAL DEL PÉNDULO ("MATES MIKE")
    // =========================================================================
    function dibujarEscenarioPendulo() {
        var canvas = dom.canvasPenduloMain;
        if (!canvas) return;
        var rect = canvas.getBoundingClientRect();
        if (rect.width <= 0 || rect.height <= 0) return;

        var dpr = window.devicePixelRatio || 1;
        var w = Math.floor(rect.width * dpr);
        var h = Math.floor(rect.height * dpr);
        if (canvas.width !== w || canvas.height !== h) {
            canvas.width = w;
            canvas.height = h;
        }

        var ctx = canvas.getContext('2d');
        ctx.save();
        ctx.scale(dpr, dpr);

        var cw = rect.width;
        var ch = rect.height;

        // Fondo
        ctx.fillStyle = color('panel');
        ctx.fillRect(0, 0, cw, ch);

        var origenX = cw * 0.5;
        var origenY = ch * 0.30;
        var escalaLongitud = Math.min(cw, ch) * 0.28 / Math.max(0.1, (estado.params.l1 + estado.params.l2));

        var xCart = PenduloFisica.posicionesCartesianas(estado.x, estado.params);
        var p1x = origenX + xCart.x1 * escalaLongitud;
        var p1y = origenY - xCart.y1 * escalaLongitud;
        var p2x = origenX + xCart.x2 * escalaLongitud;
        var p2y = origenY - xCart.y2 * escalaLongitud;

        // 1. Estela térmica de la masa 2
        if (estado.mostrarEstela && estado.estelaM2.length > 1) {
            ctx.save();
            ctx.lineWidth = 2.5;
            ctx.lineCap = 'round';
            var trailLen = estado.estelaM2.length;
            for (var k = 1; k < trailLen; k++) {
                var ptPrev = estado.estelaM2[k - 1];
                var ptCurr = estado.estelaM2[k];
                var alpha = (k / trailLen) * 0.75;

                // Color según velocidad
                var speedRatio = Math.min(1.0, Math.sqrt(ptCurr.v2Sq) / 10.0);
                var rgb;
                if (speedRatio < 0.5) {
                    rgb = color('traceSlow');
                } else if (speedRatio < 0.8) {
                    rgb = color('traceMid');
                } else {
                    rgb = color('traceFast');
                }

                ctx.strokeStyle = 'rgba(' + rgb[0] + ',' + rgb[1] + ',' + rgb[2] + ',' + alpha + ')';
                ctx.beginPath();
                ctx.moveTo(origenX + ptPrev.x * escalaLongitud, origenY - ptPrev.y * escalaLongitud);
                ctx.lineTo(origenX + ptCurr.x * escalaLongitud, origenY - ptCurr.y * escalaLongitud);
                ctx.stroke();
            }
            ctx.restore();
        }

        // 2. Modelo Lineal Fantasma (si está activo el checkbox)
        if (estado.modoLinealActivo) {
            var xCartLin = PenduloFisica.posicionesCartesianas(estado.xLineal, estado.params);
            var lin1x = origenX + xCartLin.x1 * escalaLongitud;
            var lin1y = origenY - xCartLin.y1 * escalaLongitud;
            var lin2x = origenX + xCartLin.x2 * escalaLongitud;
            var lin2y = origenY - xCartLin.y2 * escalaLongitud;

            ctx.save();
            ctx.setLineDash([5, 5]);
            ctx.strokeStyle = color('rodGhost');
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.moveTo(origenX, origenY);
            ctx.lineTo(lin1x, lin1y);
            ctx.lineTo(lin2x, lin2y);
            ctx.stroke();

            // Bob fantasma 1
            ctx.beginPath();
            ctx.arc(lin1x, lin1y, 9, 0, 2 * Math.PI);
            ctx.fillStyle = 'rgba(76, 201, 240, 0.3)';
            ctx.fill();

            // Bob fantasma 2
            ctx.beginPath();
            ctx.arc(lin2x, lin2y, 9, 0, 2 * Math.PI);
            ctx.fillStyle = 'rgba(247, 37, 133, 0.3)';
            ctx.fill();
            ctx.restore();
        }

        // 3. Cuñas sombreadas de ángulo (radio 45px desde vertical punteada)
        dibujarCunaAngulo(ctx, origenX, origenY, estado.x[0], 45, color('th1'));
        dibujarCunaAngulo(ctx, p1x, p1y, estado.x[1], 45, color('th2'));

        // 4. Varillas rígidas de masa despreciable
        ctx.strokeStyle = color('rod');
        ctx.lineWidth = 4.5;
        ctx.lineCap = 'round';
        ctx.beginPath();
        ctx.moveTo(origenX, origenY);
        ctx.lineTo(p1x, p1y);
        ctx.stroke();

        ctx.beginPath();
        ctx.moveTo(p1x, p1y);
        ctx.lineTo(p2x, p2y);
        ctx.stroke();

        // Pivote fijo
        ctx.fillStyle = color('text');
        ctx.beginPath();
        ctx.arc(origenX, origenY, 5, 0, 2 * Math.PI);
        ctx.fill();

        // 5. Masas puntuales (m1 en Cian, m2 en Magenta) con etiquetas laterales
        // Masa 1
        var radioM1 = 13 + Math.min(6, estado.params.m1 * 2);
        ctx.fillStyle = color('th1');
        ctx.beginPath();
        ctx.arc(p1x, p1y, radioM1, 0, 2 * Math.PI);
        ctx.fill();
        ctx.strokeStyle = color('bg');
        ctx.lineWidth = 2;
        ctx.stroke();

        // Etiqueta m1 al lado derecho
        ctx.fillStyle = color('text');
        ctx.font = canvasFont(13, 'bold');
        ctx.textAlign = 'left';
        ctx.fillText('m₁ (' + estado.params.m1.toFixed(1) + ' kg)', p1x + radioM1 + 6, p1y + 4);

        // Masa 2
        var radioM2 = 13 + Math.min(6, estado.params.m2 * 2);
        ctx.fillStyle = color('th2');
        ctx.beginPath();
        ctx.arc(p2x, p2y, radioM2, 0, 2 * Math.PI);
        ctx.fill();
        ctx.strokeStyle = color('bg');
        ctx.lineWidth = 2;
        ctx.stroke();

        // Etiqueta m2 al lado derecho
        ctx.fillStyle = color('text');
        ctx.font = canvasFont(13, 'bold');
        ctx.textAlign = 'left';
        ctx.fillText('m₂ (' + estado.params.m2.toFixed(1) + ' kg)', p2x + radioM2 + 6, p2y + 4);

        ctx.restore();
    }

    function dibujarCunaAngulo(ctx, cx, cy, angulo, radio, estiloColor) {
        ctx.save();
        // Línea vertical punteada de referencia (hacia abajo)
        ctx.strokeStyle = color('textMuted');
        ctx.lineWidth = 1.2;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(cx, cy + radio + 10);
        ctx.stroke();
        ctx.setLineDash([]);

        // Cuña sombreada
        var angNorm = Math.atan2(Math.sin(angulo), Math.cos(angulo));
        var anguloInicio = Math.PI * 0.5; // Vertical hacia abajo en canvas (+Y)
        var anguloFin = anguloInicio - angNorm; // Hacia la varilla
        var sentidoAntihorario = angNorm > 0;

        ctx.fillStyle = estiloColor.replace(')', ', 0.22)').replace('rgb', 'rgba');
        ctx.strokeStyle = estiloColor;
        ctx.lineWidth = 1.8;

        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.arc(cx, cy, radio, anguloInicio, anguloFin, sentidoAntihorario);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();

        ctx.restore();
    }

    // =========================================================================
    // 9. RENDERIZADO DE GRÁFICAS DE LA COLUMNA DERECHA EN PÉNDULO
    // =========================================================================
    function dibujarGraficaDerechaPendulo() {
        var canvas = dom.canvasPenduloChart;
        if (!canvas) return;
        var rect = canvas.getBoundingClientRect();
        if (rect.width <= 0 || rect.height <= 0) return;

        var dpr = window.devicePixelRatio || 1;
        var w = Math.floor(rect.width * dpr);
        var h = Math.floor(rect.height * dpr);
        if (canvas.width !== w || canvas.height !== h) {
            canvas.width = w;
            canvas.height = h;
        }

        var ctx = canvas.getContext('2d');
        ctx.save();
        ctx.scale(dpr, dpr);

        var cw = rect.width;
        var ch = rect.height;

        ctx.fillStyle = color('panel');
        ctx.fillRect(0, 0, cw, ch);

        var padLeft = 46;
        var padRight = 16;
        var padTop = 16;
        var padBottom = 30;
        var plotW = cw - padLeft - padRight;
        var plotH = ch - padTop - padBottom;

        if (plotW <= 10 || plotH <= 10) {
            ctx.restore();
            return;
        }

        var historial = estado.historialTiempo;
        if (historial.length < 2) {
            ctx.restore();
            return;
        }

        if (!estado.vistaEspacioFase) {
            // MODO 1: ÁNGULOS EN EL TIEMPO θ₁(t) y θ₂(t)
            var tMax = estado.t;
            var tMin = Math.max(0.0, tMax - estado.ventanaS);
            if (tMax - tMin < 1.0) tMax = tMin + 1.0;

            var yMin = -Math.PI;
            var yMax = Math.PI;

            // Cuadrícula y ejes
            ctx.strokeStyle = color('grid');
            ctx.lineWidth = 1;
            ctx.beginPath();
            var yZero = padTop + plotH * 0.5;
            ctx.moveTo(padLeft, yZero);
            ctx.lineTo(padLeft + plotW, yZero);
            ctx.stroke();

            ctx.fillStyle = color('textMuted');
            ctx.font = canvasFont(11);
            ctx.textAlign = 'right';
            ctx.fillText('π', padLeft - 6, padTop + 10);
            ctx.fillText('0', padLeft - 6, yZero + 4);
            ctx.fillText('-π', padLeft - 6, padTop + plotH);

            ctx.textAlign = 'center';
            ctx.fillText((tMin).toFixed(0) + 's', padLeft, padTop + plotH + 18);
            ctx.fillText((tMax).toFixed(0) + 's', padLeft + plotW, padTop + plotH + 18);

            // Curvas θ₁(t) y θ₂(t)
            dibujarSerieTemporal(ctx, historial, 'th1', tMin, tMax, yMin, yMax, padLeft, padTop, plotW, plotH, color('th1'), false);
            dibujarSerieTemporal(ctx, historial, 'th2', tMin, tMax, yMin, yMax, padLeft, padTop, plotW, plotH, color('th2'), false);

            // Si modo lineal está activo, dibujar curvas punteadas
            if (estado.modoLinealActivo) {
                dibujarSerieTemporal(ctx, historial, 'th1Lin', tMin, tMax, yMin, yMax, padLeft, padTop, plotW, plotH, color('th1'), true);
                dibujarSerieTemporal(ctx, historial, 'th2Lin', tMin, tMax, yMin, yMax, padLeft, padTop, plotW, plotH, color('th2'), true);
            }
        } else {
            // MODO 2: ESPACIO FASE (θ, ω) envolviendo en (-π, π] y cortando saltos
            ctx.strokeStyle = color('grid');
            ctx.lineWidth = 1;
            ctx.strokeRect(padLeft, padTop, plotW, plotH);

            var cx = padLeft + plotW * 0.5;
            var cy = padTop + plotH * 0.5;
            ctx.beginPath();
            ctx.moveTo(padLeft, cy); ctx.lineTo(padLeft + plotW, cy);
            ctx.moveTo(cx, padTop); ctx.lineTo(cx, padTop + plotH);
            ctx.stroke();

            ctx.fillStyle = color('textMuted');
            ctx.font = canvasFont(11);
            ctx.textAlign = 'center';
            ctx.fillText('θ ∈ (-π, π]', cx, padTop + plotH + 18);

            dibujarTrayectoriaFaseCorte(ctx, historial, 'th1', 'w1', padLeft, padTop, plotW, plotH, color('th1'));
            dibujarTrayectoriaFaseCorte(ctx, historial, 'th2', 'w2', padLeft, padTop, plotW, plotH, color('th2'));
        }

        ctx.restore();
    }

    function dibujarSerieTemporal(ctx, data, campo, tMin, tMax, yMin, yMax, x0, y0, w, h, estiloColor, esDashed) {
        ctx.save();
        ctx.strokeStyle = estiloColor;
        ctx.lineWidth = esDashed ? 1.5 : 2.2;
        if (esDashed) ctx.setLineDash([4, 4]);

        ctx.beginPath();
        var iniciado = false;
        var prevVal = 0;

        for (var i = 0; i < data.length; i++) {
            var pt = data[i];
            if (pt.t < tMin) continue;

            // Envolver a [-pi, pi]
            var val = pt[campo];
            var valWrapped = Math.atan2(Math.sin(val), Math.cos(val));

            var px = x0 + ((pt.t - tMin) / (tMax - tMin)) * w;
            var py = y0 + h - ((valWrapped - yMin) / (yMax - yMin)) * h;

            // Cortar saltos cuando cruza frontera ±pi
            if (!iniciado) {
                ctx.moveTo(px, py);
                iniciado = true;
            } else {
                if (Math.abs(valWrapped - prevVal) > Math.PI) {
                    ctx.moveTo(px, py);
                } else {
                    ctx.lineTo(px, py);
                }
            }
            prevVal = valWrapped;
        }
        ctx.stroke();
        ctx.restore();
    }

    function dibujarTrayectoriaFaseCorte(ctx, data, campoTh, campoW, x0, y0, w, h, estiloColor) {
        ctx.save();
        ctx.strokeStyle = estiloColor;
        ctx.lineWidth = 1.8;

        var wMax = 12.0; // Rango típico de velocidad angular [-12, 12]

        ctx.beginPath();
        var iniciado = false;
        var prevTh = 0;

        for (var i = 0; i < data.length; i++) {
            var pt = data[i];
            var thWrapped = Math.atan2(Math.sin(pt[campoTh]), Math.cos(pt[campoTh]));
            var vel = Math.max(-wMax, Math.min(wMax, pt[campoW]));

            var px = x0 + ((thWrapped + Math.PI) / (2 * Math.PI)) * w;
            var py = y0 + h * 0.5 - (vel / (2 * wMax)) * h;

            if (!iniciado) {
                ctx.moveTo(px, py);
                iniciado = true;
            } else {
                if (Math.abs(thWrapped - prevTh) > Math.PI) {
                    ctx.moveTo(px, py);
                } else {
                    ctx.lineTo(px, py);
                }
            }
            prevTh = thWrapped;
        }
        ctx.stroke();
        ctx.restore();
    }

    // =========================================================================
    // 10. RENDERIZADO DE LA PESTAÑA GEMELOS
    // =========================================================================
    function dibujarPestañaGemelos() {
        // 1. Escenario de péndulos gemelos superpuestos
        var canvasStage = dom.canvasGemelosStage;
        if (canvasStage) {
            var rectS = canvasStage.getBoundingClientRect();
            var dprS = window.devicePixelRatio || 1;
            canvasStage.width = Math.floor(rectS.width * dprS);
            canvasStage.height = Math.floor(rectS.height * dprS);
            var ctxS = canvasStage.getContext('2d');
            ctxS.save();
            ctxS.scale(dprS, dprS);

            var cw = rectS.width;
            var ch = rectS.height;
            ctxS.fillStyle = color('panel');
            ctxS.fillRect(0, 0, cw, ch);

            var ox = cw * 0.5;
            var oy = ch * 0.30;
            var esc = Math.min(cw, ch) * 0.28 / Math.max(0.1, (estado.params.l1 + estado.params.l2));

            var colBase = color('th1');     // Cian (#4cc9f0 / #0284c7)
            var colGem = color('energy');    // Ámbar (#ffb703 / #b45309)

            // 1. Péndulo Base (Cian completo: varillas y ambas masas)
            var cBase = PenduloFisica.posicionesCartesianas(estado.x, estado.params);
            var b1x = ox + cBase.x1 * esc, b1y = oy - cBase.y1 * esc;
            var b2x = ox + cBase.x2 * esc, b2y = oy - cBase.y2 * esc;

            ctxS.strokeStyle = colBase;
            ctxS.lineWidth = 3.5;
            ctxS.lineCap = 'round';
            ctxS.beginPath();
            ctxS.moveTo(ox, oy); ctxS.lineTo(b1x, b1y); ctxS.lineTo(b2x, b2y);
            ctxS.stroke();

            ctxS.fillStyle = colBase;
            ctxS.beginPath(); ctxS.arc(b1x, b1y, 9, 0, 2 * Math.PI); ctxS.fill();
            ctxS.strokeStyle = color('bg'); ctxS.lineWidth = 1.5; ctxS.stroke();
            ctxS.beginPath(); ctxS.arc(b2x, b2y, 10, 0, 2 * Math.PI); ctxS.fill();
            ctxS.stroke();

            // 2. Péndulo Gemelo (Ámbar completo: varillas y ambas masas)
            if (estado.xGemelo) {
                var cGem = PenduloFisica.posicionesCartesianas(estado.xGemelo, estado.params);
                var g1x = ox + cGem.x1 * esc, g1y = oy - cGem.y1 * esc;
                var g2x = ox + cGem.x2 * esc, g2y = oy - cGem.y2 * esc;

                ctxS.strokeStyle = colGem;
                ctxS.lineWidth = 3.0;
                ctxS.lineCap = 'round';
                ctxS.beginPath();
                ctxS.moveTo(ox, oy); ctxS.lineTo(g1x, g1y); ctxS.lineTo(g2x, g2y);
                ctxS.stroke();

                ctxS.fillStyle = colGem;
                ctxS.beginPath(); ctxS.arc(g1x, g1y, 9, 0, 2 * Math.PI); ctxS.fill();
                ctxS.strokeStyle = color('bg'); ctxS.lineWidth = 1.5; ctxS.stroke();
                ctxS.beginPath(); ctxS.arc(g2x, g2y, 10, 0, 2 * Math.PI); ctxS.fill();
                ctxS.stroke();
            }

            // Pivote pequeño (r = 3 px)
            ctxS.fillStyle = color('text');
            ctxS.beginPath();
            ctxS.arc(ox, oy, 3, 0, 2 * Math.PI);
            ctxS.fill();

            // Leyenda de identificación en escenario
            var legW = 195;
            var legH = 46;
            var legX = 12;
            var legY = 12;
            ctxS.fillStyle = (temaActual === 'dark') ? 'rgba(18, 24, 38, 0.88)' : 'rgba(255, 255, 255, 0.88)';
            ctxS.strokeStyle = color('border');
            ctxS.lineWidth = 1;
            ctxS.fillRect(legX, legY, legW, legH);
            ctxS.strokeRect(legX, legY, legW, legH);

            // Item 1: Péndulo base (Cian)
            ctxS.fillStyle = colBase;
            ctxS.beginPath(); ctxS.arc(legX + 16, legY + 15, 6, 0, 2 * Math.PI); ctxS.fill();
            ctxS.fillStyle = color('text');
            ctxS.font = canvasFont(11);
            ctxS.textAlign = 'left';
            ctxS.fillText('Péndulo 1 (base)', legX + 28, legY + 19);

            // Item 2: Péndulo gemelo (Ámbar)
            ctxS.fillStyle = colGem;
            ctxS.beginPath(); ctxS.arc(legX + 16, legY + 31, 6, 0, 2 * Math.PI); ctxS.fill();
            ctxS.fillStyle = color('text');
            ctxS.font = canvasFont(11);
            ctxS.fillText('Péndulo 2 (+0.001 rad)', legX + 28, legY + 35);

            ctxS.restore();
        }

        // 2. Gráfica de Separación log10(d) vs t
        var canvasSep = dom.canvasGemelosSeparation;
        if (canvasSep) {
            var rectC = canvasSep.getBoundingClientRect();
            var dprC = window.devicePixelRatio || 1;
            canvasSep.width = Math.floor(rectC.width * dprC);
            canvasSep.height = Math.floor(rectC.height * dprC);
            var ctxC = canvasSep.getContext('2d');
            ctxC.save();
            ctxC.scale(dprC, dprC);

            var pw = rectC.width;
            var ph = rectC.height;
            ctxC.fillStyle = color('panel');
            ctxC.fillRect(0, 0, pw, ph);

            var padL = 46, padR = 16, padT = 16, padB = 30;
            var gw = pw - padL - padR;
            var gh = ph - padT - padB;

            if (gw > 10 && gh > 10) {
                var h = estado.historialGemelos;
                var logMin = -4.0; // d = 10^-4
                var logMax = 1.0;  // d = 10^1 = 10
                var tMax = Math.max(10.0, estado.t);

                // Cuadrícula y ejes
                ctxC.strokeStyle = color('grid');
                ctxC.lineWidth = 1;
                ctxC.strokeRect(padL, padT, gw, gh);

                ctxC.fillStyle = color('textMuted');
                ctxC.font = canvasFont(11);
                ctxC.textAlign = 'right';
                ctxC.fillText('10⁰', padL - 6, padT + gh * 0.2);
                ctxC.fillText('10⁻²', padL - 6, padT + gh * 0.6);
                ctxC.fillText('10⁻⁴', padL - 6, padT + gh);

                ctxC.textAlign = 'center';
                ctxC.fillText('0s', padL, padT + gh + 18);
                ctxC.fillText(tMax.toFixed(0) + 's', padL + gw, padT + gh + 18);

                // Curva log10(d)
                if (h.length > 1) {
                    ctxC.strokeStyle = color('energy');
                    ctxC.lineWidth = 2.2;
                    ctxC.beginPath();
                    for (var k = 0; k < h.length; k++) {
                        var px = padL + (h[k].t / tMax) * gw;
                        var py = padT + gh - ((h[k].logD - logMin) / (logMax - logMin)) * gh;
                        py = Math.max(padT, Math.min(padT + gh, py));
                        if (k === 0) ctxC.moveTo(px, py);
                        else ctxC.lineTo(px, py);
                    }
                    ctxC.stroke();
                }
            }
            ctxC.restore();
        }
    }

    // =========================================================================
    // 11. RENDERIZADO DE LA PESTAÑA MONTE CARLO
    // =========================================================================
    function dibujarPestañaMonteCarlo() {
        // Panel Izquierdo: Stage + Histograma
        dibujarMiniEnsemble(
            dom.canvasMcStageLeft,
            estado.mc.left.pendulos,
            color('th1'),
            estado.mc.left.puntosIniciales,
            [0.6 * 180 / Math.PI, 0.6 * 180 / Math.PI]
        );
        dibujarHistogramaFlips(dom.canvasMcHistLeft, estado.mc.left.flips, color('th1'));

        // Panel Derecho: Stage + Histograma
        dibujarMiniEnsemble(
            dom.canvasMcStageRight,
            estado.mc.right.pendulos,
            color('th2'),
            estado.mc.right.puntosIniciales,
            [estado.mc.right.x0[0] * 180 / Math.PI, estado.mc.right.x0[1] * 180 / Math.PI]
        );
        dibujarHistogramaFlips(dom.canvasMcHistRight, estado.mc.right.flips, color('th2'));
    }

    function dibujarMiniEnsemble(canvas, pendulos, estiloColor, puntosIniciales, centroNominalDeg) {
        if (!canvas) return;
        var rect = canvas.getBoundingClientRect();
        if (rect.width <= 0 || rect.height <= 0) return;

        var dpr = window.devicePixelRatio || 1;
        canvas.width = Math.floor(rect.width * dpr);
        canvas.height = Math.floor(rect.height * dpr);
        var ctx = canvas.getContext('2d');
        ctx.save();
        ctx.scale(dpr, dpr);

        var cw = rect.width;
        var ch = rect.height;
        ctx.fillStyle = color('panelSub');
        ctx.fillRect(0, 0, cw, ch);

        var ox = cw * 0.5;
        var oy = ch * 0.30;
        var esc = Math.min(cw, ch) * 0.28 / Math.max(0.1, (estado.params.l1 + estado.params.l2));

        var n = pendulos.length;
        // Dibujar varillas muy tenues
        ctx.strokeStyle = 'rgba(150, 160, 180, 0.08)';
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (var i = 0; i < n; i += 2) {
            var c = PenduloFisica.posicionesCartesianas(pendulos[i], estado.params);
            ctx.moveTo(ox, oy);
            ctx.lineTo(ox + c.x1 * esc, oy - c.y1 * esc);
            ctx.lineTo(ox + c.x2 * esc, oy - c.y2 * esc);
        }
        ctx.stroke();

        // Dibujar puntos brillantes de la masa 2
        ctx.fillStyle = estiloColor;
        for (var j = 0; j < n; j++) {
            var c2 = PenduloFisica.posicionesCartesianas(pendulos[j], estado.params);
            var px = ox + c2.x2 * esc;
            var py = oy - c2.y2 * esc;
            ctx.beginPath();
            ctx.arc(px, py, 2.5, 0, 2 * Math.PI);
            ctx.fill();
        }

        // Pivote pequeño
        ctx.fillStyle = color('text');
        ctx.beginPath();
        ctx.arc(ox, oy, 2.5, 0, 2 * Math.PI);
        ctx.fill();

        // Recuadro "Condiciones iniciales" con los N puntos en grados
        if (puntosIniciales && puntosIniciales.length > 0 && centroNominalDeg) {
            var boxW = Math.min(150, cw * 0.36);
            var boxH = Math.min(110, ch * 0.55);
            var boxX = 10;
            var boxY = 10;

            ctx.fillStyle = (temaActual === 'dark') ? 'rgba(18, 24, 38, 0.88)' : 'rgba(255, 255, 255, 0.88)';
            ctx.strokeStyle = color('border');
            ctx.lineWidth = 1;
            ctx.fillRect(boxX, boxY, boxW, boxH);
            ctx.strokeRect(boxX, boxY, boxW, boxH);

            ctx.fillStyle = color('text');
            ctx.font = canvasFont(10, 'bold');
            ctx.textAlign = 'left';
            ctx.fillText('Condiciones iniciales', boxX + 6, boxY + 14);

            var sigDeg = (estado.mc.sigma * 180.0 / Math.PI).toFixed(2);
            ctx.fillStyle = color('textMuted');
            ctx.font = canvasFont(9);
            ctx.fillText('error σ = ' + sigDeg + '°', boxX + 6, boxY + 26);

            var spPadX = 6;
            var spTop = boxY + 32;
            var spW = boxW - 12;
            var spH = boxH - 38;

            ctx.strokeStyle = color('border');
            ctx.strokeRect(boxX + spPadX, spTop, spW, spH);

            var midX = boxX + spPadX + spW * 0.5;
            var midY = spTop + spH * 0.5;
            ctx.strokeStyle = 'rgba(150, 160, 180, 0.3)';
            ctx.beginPath();
            ctx.moveTo(midX, spTop); ctx.lineTo(midX, spTop + spH);
            ctx.moveTo(boxX + spPadX, midY); ctx.lineTo(boxX + spPadX + spW, midY);
            ctx.stroke();

            var spanDeg = Math.max(1.5, parseFloat(sigDeg) * 2.8);
            var cTh1 = centroNominalDeg[0];
            var cTh2 = centroNominalDeg[1];

            ctx.fillStyle = estiloColor;
            for (var pIdx = 0; pIdx < puntosIniciales.length; pIdx++) {
                var pt = puntosIniciales[pIdx];
                var dTh1 = pt[0] - cTh1;
                var dTh2 = pt[1] - cTh2;
                var ppx = midX + (dTh1 / spanDeg) * (spW * 0.46);
                var ppy = midY - (dTh2 / spanDeg) * (spH * 0.46);
                if (ppx >= boxX + spPadX && ppx <= boxX + spPadX + spW && ppy >= spTop && ppy <= spTop + spH) {
                    ctx.beginPath();
                    ctx.arc(ppx, ppy, 1.8, 0, 2 * Math.PI);
                    ctx.fill();
                }
            }
        }

        ctx.restore();
    }

    function dibujarHistogramaFlips(canvas, flips, estiloColor) {
        if (!canvas) return;
        var rect = canvas.getBoundingClientRect();
        if (rect.width <= 0 || rect.height <= 0) return;

        var dpr = window.devicePixelRatio || 1;
        canvas.width = Math.floor(rect.width * dpr);
        canvas.height = Math.floor(rect.height * dpr);
        var ctx = canvas.getContext('2d');
        ctx.save();
        ctx.scale(dpr, dpr);

        var cw = rect.width;
        var ch = rect.height;
        ctx.fillStyle = color('panelSub');
        ctx.fillRect(0, 0, cw, ch);

        var padL = 50, padR = 14, padT = 16, padB = 42;
        var gw = cw - padL - padR;
        var gh = ch - padT - padB;

        if (gw <= 10 || gh <= 10) {
            ctx.restore();
            return;
        }

        // Determinar max flips entre ambos paneles para escala compartida
        var maxFlipObs = 0;
        var n = flips.length;
        for (var i = 0; i < n; i++) {
            if (flips[i] > maxFlipObs) maxFlipObs = flips[i];
        }
        var maxBins = Math.max(6, maxFlipObs + 1);

        // Contar frecuencias
        var cuentas = new Array(maxBins).fill(0);
        for (var j = 0; j < n; j++) {
            var f = flips[j];
            if (f < maxBins) cuentas[f]++;
        }

        var maxCuenta = Math.max(1, Math.max.apply(null, cuentas));

        // Ejes
        ctx.strokeStyle = color('border');
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(padL, padT);
        ctx.lineTo(padL, padT + gh);
        ctx.lineTo(padL + gw, padT + gh);
        ctx.stroke();

        // Etiqueta vertical "Número de péndulos"
        ctx.save();
        ctx.translate(16, padT + gh * 0.5);
        ctx.rotate(-Math.PI * 0.5);
        ctx.fillStyle = color('text');
        ctx.font = canvasFont(11);
        ctx.textAlign = 'center';
        ctx.fillText('Número de péndulos', 0, 0);
        ctx.restore();

        // Barras
        var barW = gw / maxBins;
        ctx.fillStyle = estiloColor;
        ctx.strokeStyle = color('bg');
        ctx.lineWidth = 1;

        for (var b = 0; b < maxBins; b++) {
            var bH = (cuentas[b] / maxCuenta) * (gh - 10);
            var bx = padL + b * barW;
            var by = padT + gh - bH;

            ctx.fillRect(bx + 2, by, barW - 4, bH);
            ctx.strokeRect(bx + 2, by, barW - 4, bH);

            // Etiquetas eje X
            if (maxBins <= 10 || b % 2 === 0) {
                ctx.fillStyle = color('textMuted');
                ctx.font = canvasFont(10);
                ctx.textAlign = 'center';
                ctx.fillText(b.toString(), bx + barW * 0.5, padT + gh + 14);
                ctx.fillStyle = estiloColor;
            }
        }

        // Título del eje X
        ctx.fillStyle = color('text');
        ctx.font = canvasFont(11);
        ctx.textAlign = 'center';
        ctx.fillText('Vueltas de la masa 2 en T = ' + Math.round(estado.mc.tMax) + ' s', padL + gw * 0.5, padT + gh + 34);

        ctx.restore();
    }

    // =========================================================================
    // 12. ACTUALIZACIÓN DEL TERMÓMETRO DE ENERGÍA Y HUD
    // =========================================================================
    function actualizarTermometroEnergia() {
        var p = estado.params;
        var umb = PenduloFisica.umbralesEnergia(p);
        var ec2 = umb.Ec_m2;
        var ec1 = umb.Ec_m1;
        var ecAmbas = umb.Ec_ambas;

        // Actualizar textos de los marcadores con sus valores en Joules
        dom.markLblM2.textContent = 'm₂ (' + ec2.toFixed(1) + ' J)';
        dom.markLblM1.textContent = 'm₁ (' + ec1.toFixed(1) + ' J)';
        dom.markLblAmbas.textContent = 'ambas (' + ecAmbas.toFixed(1) + ' J)';

        // Posiciones relativas en la barra (escala de 0 a Ec_ambas * 1.15)
        var maxEscala = ecAmbas * 1.15;
        var posM2 = Math.min(100, Math.max(0, (ec2 / maxEscala) * 100));
        var posM1 = Math.min(100, Math.max(0, (ec1 / maxEscala) * 100));
        var posAmbas = Math.min(100, Math.max(0, (ecAmbas / maxEscala) * 100));

        dom.markM2.style.left = posM2 + '%';
        dom.markM1.style.left = posM1 + '%';
        dom.markAmbas.style.left = posAmbas + '%';

        // Llenado actual
        var eActual = PenduloFisica.energiaTotal(estado.x, estado.params);
        var posActual = Math.min(100, Math.max(0, (eActual / maxEscala) * 100));
        dom.thermometerFill.style.width = posActual + '%';

        dom.readoutTotalE.textContent = eActual.toFixed(2);
        dom.readoutFlips.textContent = estado.vueltasM2.toString();
    }

    function actualizarIndicadoresHUD() {
        // Números grandes de ángulos
        var th1 = estado.x[0];
        var th2 = estado.x[1];
        var th1Deg = th1 * (180.0 / Math.PI);
        var th2Deg = th2 * (180.0 / Math.PI);

        dom.bigAngleTh1.textContent = th1.toFixed(2) + ' rad (' + th1Deg.toFixed(0) + '°)';
        dom.bigAngleTh2.textContent = th2.toFixed(2) + ' rad (' + th2Deg.toFixed(0) + '°)';
        dom.bigTime.textContent = estado.t.toFixed(2) + ' s';

        // Advertencia lineal en texto plano (NO emojis)
        if (estado.modoLinealActivo) {
            if (Math.abs(th1) > 0.3 || Math.abs(th2) > 0.3) {
                dom.linearWarning.style.display = 'block';
            } else {
                dom.linearWarning.style.display = 'none';
            }
        } else {
            dom.linearWarning.style.display = 'none';
        }

        // Error relativo de energía
        dom.valErrorDeltaE.textContent = estado.deltaERel.toExponential(2);
    }

    function actualizarParametrosLabels() {
        dom.lblTh1.textContent = parseFloat(dom.sliderTh1.value).toFixed(0) + '°';
        dom.lblTh2.textContent = parseFloat(dom.sliderTh2.value).toFixed(0) + '°';

        dom.lblM1.textContent = parseFloat(dom.sliderM1.value).toFixed(2) + ' kg';
        dom.lblM2.textContent = parseFloat(dom.sliderM2.value).toFixed(2) + ' kg';
        dom.lblL1.textContent = parseFloat(dom.sliderL1.value).toFixed(2) + ' m';
        dom.lblL2.textContent = parseFloat(dom.sliderL2.value).toFixed(2) + ' m';

        var m1 = parseFloat(dom.sliderM1.value);
        var m2 = parseFloat(dom.sliderM2.value);
        var l1 = parseFloat(dom.sliderL1.value);
        var l2 = parseFloat(dom.sliderL2.value);
        dom.lblMu.textContent = (m2 / m1).toFixed(2);
        dom.lblLambda.textContent = (l2 / l1).toFixed(2);

        dom.lblG.textContent = parseFloat(dom.sliderG.value).toFixed(2) + ' m/s²';
        dom.lblDeltaGemelos.textContent = parseFloat(dom.sliderDeltaGemelos.value).toFixed(4) + ' rad';

        dom.lblMcN.textContent = dom.sliderMcN.value;
        var sigDeg = (parseFloat(dom.sliderMcSigma.value) * 180.0 / Math.PI).toFixed(2);
        dom.lblMcSigma.textContent = sigDeg + '°';
        dom.lblMcT.textContent = dom.sliderMcT.value + ' s';
    }

    function aplicarPreset(pKey) {
        if (pKey === 'pequenas_oscilaciones') {
            if (dom.sliderTh1) dom.sliderTh1.value = (0.20 * 180 / Math.PI).toFixed(1);
            if (dom.sliderTh2) dom.sliderTh2.value = (0.25 * 180 / Math.PI).toFixed(1);
        } else if (pKey === 'modo_normal_en_fase') {
            var m1 = parseFloat(dom.sliderM1.value);
            var m2 = parseFloat(dom.sliderM2.value);
            var l1 = parseFloat(dom.sliderL1.value);
            var l2 = parseFloat(dom.sliderL2.value);
            var g = parseFloat(dom.sliderG.value);
            var mu = m2 / m1;
            var a20 = -(1.0 + mu) * g / l1;
            var a21 = mu * g / l1;
            var a30 = (1.0 + mu) * g / l2;
            var a31 = -(1.0 + mu) * g / l2;
            var b = -(a20 + a31);
            var c = a20 * a31 - a21 * a30;
            var disc = Math.sqrt(Math.max(0.0, b * b - 4.0 * c));
            var w1_sq = 0.5 * (b - disc);
            var ratio = Math.abs(a21) > 1e-9 ? (-(a20 + w1_sq) / a21) : Math.SQRT2;
            var th1 = 0.1;
            var th2 = ratio * th1;
            if (dom.sliderTh1) dom.sliderTh1.value = (th1 * 180 / Math.PI).toFixed(1);
            if (dom.sliderTh2) dom.sliderTh2.value = (th2 * 180 / Math.PI).toFixed(1);
        } else if (pKey === 'caotico') {
            if (dom.sliderTh1) dom.sliderTh1.value = (2.6 * 180 / Math.PI).toFixed(1);
            if (dom.sliderTh2) dom.sliderTh2.value = (2.6 * 180 / Math.PI).toFixed(1);
        } else if (pKey === 'casi_invertido') {
            var th = Math.PI - 0.05;
            if (dom.sliderTh1) dom.sliderTh1.value = (th * 180 / Math.PI).toFixed(1);
            if (dom.sliderTh2) dom.sliderTh2.value = (th * 180 / Math.PI).toFixed(1);
        } else if (window.PRESETS_PENDULO && window.PRESETS_PENDULO[pKey]) {
            var preset = window.PRESETS_PENDULO[pKey];
            if (dom.sliderTh1) dom.sliderTh1.value = (preset.x0[0] * 180 / Math.PI).toFixed(1);
            if (dom.sliderTh2) dom.sliderTh2.value = (preset.x0[1] * 180 / Math.PI).toFixed(1);
            if (preset.params) {
                if (dom.sliderM1) dom.sliderM1.value = preset.params.m1;
                if (dom.sliderM2) dom.sliderM2.value = preset.params.m2;
                if (dom.sliderL1) dom.sliderL1.value = preset.params.l1;
                if (dom.sliderL2) dom.sliderL2.value = preset.params.l2;
                if (dom.sliderG) dom.sliderG.value = preset.params.g;
            }
        }

        actualizarParametrosLabels();
        reiniciarSimulacion();
    }

    // =========================================================================
    // 13. NAVEGACIÓN Y CAMBIO DE PESTAÑAS
    // =========================================================================
    function cambiarPestaña(nuevaTab) {
        if (!dom.tabViews[nuevaTab]) return;

        estado.tabActiva = nuevaTab;

        // Botones nav
        dom.tabBtns.forEach(function (btn) {
            if (btn.getAttribute('data-tab') === nuevaTab) {
                btn.classList.add('active');
            } else {
                btn.classList.remove('active');
            }
        });

        // Vistas de contenido
        Object.keys(dom.tabViews).forEach(function (tabKey) {
            if (tabKey === nuevaTab) {
                dom.tabViews[tabKey].style.display = 'block';
                dom.tabViews[tabKey].classList.add('active');
            } else {
                dom.tabViews[tabKey].style.display = 'none';
                dom.tabViews[tabKey].classList.remove('active');
            }
        });

        // Controles contextuales en barra lateral
        Object.keys(dom.controlsTab).forEach(function (ctrlKey) {
            if (ctrlKey === nuevaTab) {
                dom.controlsTab[ctrlKey].style.display = 'block';
            } else {
                dom.controlsTab[ctrlKey].style.display = 'none';
            }
        });

        if (nuevaTab === 'gemelos') {
            aplicarPreset('caotico');
        } else if (nuevaTab === 'montecarlo') {
            reiniciarMonteCarlo();
        } else if (nuevaTab === 'mapacaos') {
            redibujarMapaCaos();
        }
    }

    // =========================================================================
    // 14. BUCLE PRINCIPAL DE ANIMACIÓN Y RENDERIZADO
    // =========================================================================
    var ultimoTiempoMs = 0;

    function renderizarVistaActual() {
        actualizarIndicadoresHUD();
        actualizarTermometroEnergia();

        if (estado.tabActiva === 'pendulo') {
            dibujarEscenarioPendulo();
            dibujarGraficaDerechaPendulo();
        } else if (estado.tabActiva === 'gemelos') {
            dibujarPestañaGemelos();
        } else if (estado.tabActiva === 'montecarlo') {
            dibujarPestañaMonteCarlo();
        } else if (estado.tabActiva === 'mapacaos') {
            redibujarMapaCaos();
        }
    }

    function bucleAnimacion(tiempoMs) {
        requestAnimationFrame(bucleAnimacion);

        if (!ultimoTiempoMs) ultimoTiempoMs = tiempoMs;
        var deltaMs = tiempoMs - ultimoTiempoMs;
        ultimoTiempoMs = tiempoMs;

        if (deltaMs > 100) deltaMs = 100; // Limitar salto largo

        if (estado.corriendo) {
            // Pasos de integración sub-frame para estabilidad (dt = 0.005 s)
            var nSubPasos = Math.round((deltaMs / 1000.0) * estado.velocidad / estado.dtBase);
            if (nSubPasos < 1) nSubPasos = 1;
            if (nSubPasos > 20) nSubPasos = 20;

            for (var sp = 0; sp < nSubPasos; sp++) {
                pasoSimulacion(estado.dtBase);
            }
            registrarHistorialPunto();

            // Si Monte Carlo está activo, avanzar ensemble
            if (estado.tabActiva === 'montecarlo' && estado.mc.corriendo) {
                var dtMc = (deltaMs / 1000.0) * estado.velocidad;
                pasoMonteCarlo(dtMc);
            }
        }

        renderizarVistaActual();
    }

    // =========================================================================
    // 15. CONFIGURACIÓN DE EVENTOS E INTERACCIÓN
    // =========================================================================
    function vincularEventos() {
        // Pestañas
        dom.tabBtns.forEach(function (btn) {
            btn.addEventListener('click', function () {
                cambiarPestaña(btn.getAttribute('data-tab'));
            });
        });

        // Tema Claro / Oscuro
        dom.btnThemeToggle.addEventListener('click', function () {
            if (temaActual === 'dark') {
                temaActual = 'light';
                dom.body.classList.remove('theme-dark');
                dom.body.classList.add('theme-light');
                dom.themeIcon.textContent = '';
                dom.themeText.textContent = 'Claro';
            } else {
                temaActual = 'dark';
                dom.body.classList.remove('theme-light');
                dom.body.classList.add('theme-dark');
                dom.themeIcon.textContent = '';
                dom.themeText.textContent = 'Oscuro';
            }
            redibujarMapaCaos();
        });

        // Alternar barra lateral
        dom.btnToggleSidebar.addEventListener('click', function () {
            dom.sidebar.classList.toggle('hidden');
        });

        // Presentación
        dom.btnPresentation.addEventListener('click', function () {
            document.documentElement.classList.toggle('presentation-mode');
            document.body.classList.toggle('mode-presentation');
            setTimeout(renderizarVistaActual, 30);
        });

        // Transporte
        dom.btnPlayPause.addEventListener('click', function () {
            estado.corriendo = !estado.corriendo;
            dom.btnPlayPause.textContent = estado.corriendo ? 'Pausa' : 'Reanudar';
            if (estado.corriendo) {
                dom.btnPlayPause.classList.add('btn-primary');
            } else {
                dom.btnPlayPause.classList.remove('btn-primary');
            }
        });

        dom.btnReset.addEventListener('click', function () {
            reiniciarSimulacion();
        });

        dom.btnStep.addEventListener('click', function () {
            estado.corriendo = false;
            dom.btnPlayPause.textContent = 'Reanudar';
            dom.btnPlayPause.classList.remove('btn-primary');
            pasoSimulacion(estado.dtBase * 2);
            registrarHistorialPunto();
            actualizarIndicadoresHUD();
            actualizarTermometroEnergia();
        });

        // Presets
        dom.presetBtns.forEach(function (btn) {
            btn.addEventListener('click', function () {
                var pKey = btn.getAttribute('data-preset');
                aplicarPreset(pKey);
            });
        });

        // Sliders CI
        function onSliderCIChange() {
            actualizarParametrosLabels();
            inicializarSimulacion();
        }
        dom.sliderTh1.addEventListener('input', onSliderCIChange);
        dom.sliderTh2.addEventListener('input', onSliderCIChange);

        // Sliders Parámetros
        function onSliderParamChange() {
            actualizarParametrosLabels();
            inicializarSimulacion();
        }
        dom.sliderM1.addEventListener('input', onSliderParamChange);
        dom.sliderM2.addEventListener('input', onSliderParamChange);
        dom.sliderL1.addEventListener('input', onSliderParamChange);
        dom.sliderL2.addEventListener('input', onSliderParamChange);

        // Gravedad
        dom.gravityBtns.forEach(function (btn) {
            btn.addEventListener('click', function () {
                var gVal = parseFloat(btn.getAttribute('data-g'));
                dom.sliderG.value = gVal;
                onSliderParamChange();
            });
        });
        dom.sliderG.addEventListener('input', onSliderParamChange);

        // Checkbox Lineal
        dom.checkLinearMode.addEventListener('change', function () {
            estado.modoLinealActivo = dom.checkLinearMode.checked;
        });

        // Alternar Retrato de Fase
        dom.btnTogglePhaseChart.addEventListener('click', function () {
            estado.vistaEspacioFase = !estado.vistaEspacioFase;
            if (estado.vistaEspacioFase) {
                dom.btnTogglePhaseChart.textContent = 'Ver ángulos θ(t)';
                dom.chartPenduloTitle.textContent = 'Espacio Fase (θ, ω)';
                dom.chartPenduloCaption.textContent = 'Plano de fases envuelto en (-π, π]. Las curvas cerradas son oscilaciones periódicas.';
            } else {
                dom.btnTogglePhaseChart.textContent = 'Ver espacio fase (θ, ω)';
                dom.chartPenduloTitle.textContent = 'Ángulos en el tiempo';
                dom.chartPenduloCaption.textContent = 'Ángulos en el tiempo. Cerca del reposo se ven como suma de dos oscilaciones (modos normales).';
            }
        });

        // Gemelos: Separación δ y presets
        dom.sliderDeltaGemelos.addEventListener('input', function () {
            actualizarParametrosLabels();
            inicializarSimulacion();
        });

        dom.btnGemelosRegular.addEventListener('click', function () {
            dom.sliderTh1.value = (0.20 * 180 / Math.PI).toFixed(0);
            dom.sliderTh2.value = (0.25 * 180 / Math.PI).toFixed(0);
            actualizarParametrosLabels();
            reiniciarSimulacion();
        });

        dom.btnGemelosCaotico.addEventListener('click', function () {
            dom.sliderTh1.value = (2.6 * 180 / Math.PI).toFixed(0);
            dom.sliderTh2.value = (2.6 * 180 / Math.PI).toFixed(0);
            actualizarParametrosLabels();
            reiniciarSimulacion();
        });

        // Monte Carlo: Controles
        dom.sliderMcN.addEventListener('input', actualizarParametrosLabels);
        dom.sliderMcSigma.addEventListener('input', actualizarParametrosLabels);
        dom.sliderMcT.addEventListener('input', actualizarParametrosLabels);
        dom.selectMcRightTarget.addEventListener('change', reiniciarMonteCarlo);
        dom.btnRunMontecarlo.addEventListener('click', reiniciarMonteCarlo);

        // Avanzado
        dom.selectIntegrator.addEventListener('change', function () {
            estado.algoritmo = dom.selectIntegrator.value;
        });
        dom.selectSpeed.addEventListener('change', function () {
            estado.velocidad = parseFloat(dom.selectSpeed.value);
        });
        dom.checkTrace.addEventListener('change', function () {
            estado.mostrarEstela = dom.checkTrace.checked;
            if (!estado.mostrarEstela) estado.estelaM2 = [];
        });

        // Mapa de caos: Hover y Click
        vincularEventosMapaCaos();

        // Botón "Ver este péndulo"
        dom.btnVerEnPendulo.addEventListener('click', function () {
            cambiarPestaña('pendulo');
        });

        // Arrastre interactivo en canvas del péndulo
        vincularArrastrePendulo();
    }

    function vincularEventosMapaCaos() {
        var canvas = dom.canvasMapaCaos;
        if (!canvas) return;

        function obtenerCoordsMapa(evt) {
            var rect = canvas.getBoundingClientRect();
            var mx = evt.clientX - rect.left;
            var my = evt.clientY - rect.top;

            var padLeft = 46;
            var padBottom = 38;
            var padTop = 16;
            var padRight = 50;
            var plotW = rect.width - padLeft - padRight;
            var plotH = rect.height - padTop - padBottom;

            if (mx < padLeft || mx > padLeft + plotW || my < padTop || my > padTop + plotH) {
                return null;
            }

            var th1 = -Math.PI + ((mx - padLeft) / plotW) * (2 * Math.PI);
            var th2 = Math.PI - ((my - padTop) / plotH) * (2 * Math.PI);

            return { th1: th1, th2: th2 };
        }

        canvas.addEventListener('mousemove', function (evt) {
            var coords = obtenerCoordsMapa(evt);
            if (!coords) return;

            dom.mapaInspectTh1.textContent = (coords.th1 * 180 / Math.PI).toFixed(1) + '°';
            dom.mapaInspectTh2.textContent = (coords.th2 * 180 / Math.PI).toFixed(1) + '°';

            // Buscar valor en grid
            if (estado.ftle.gridObj) {
                var res = estado.ftle.gridObj.res;
                var col = Math.floor(((coords.th1 + Math.PI) / (2 * Math.PI)) * res);
                var row = Math.floor(((coords.th2 + Math.PI) / (2 * Math.PI)) * res);
                col = Math.max(0, Math.min(res - 1, col));
                row = Math.max(0, Math.min(res - 1, row));
                var val = estado.ftle.gridObj.grid[row * res + col];
                dom.mapaInspectLambda.textContent = (val !== undefined) ? val.toFixed(2) : '--';
            }
        });

        canvas.addEventListener('click', function (evt) {
            var coords = obtenerCoordsMapa(evt);
            if (!coords) return;

            dom.sliderTh1.value = (coords.th1 * 180 / Math.PI).toFixed(0);
            dom.sliderTh2.value = (coords.th2 * 180 / Math.PI).toFixed(0);

            actualizarParametrosLabels();
            inicializarSimulacion();
        });
    }

    function vincularArrastrePendulo() {
        var canvas = dom.canvasPenduloMain;
        if (!canvas) return;

        function getMousePos(evt) {
            var rect = canvas.getBoundingClientRect();
            var clientX = evt.clientX || (evt.touches && evt.touches[0].clientX);
            var clientY = evt.clientY || (evt.touches && evt.touches[0].clientY);
            return {
                x: clientX - rect.left,
                y: clientY - rect.top,
                cw: rect.width,
                ch: rect.height
            };
        }

        function onPointerDown(evt) {
            var pos = getMousePos(evt);
            var ox = pos.cw * 0.5;
            var oy = pos.ch * 0.30;
            var esc = Math.min(pos.cw, pos.ch) * 0.28 / Math.max(0.1, (estado.params.l1 + estado.params.l2));

            var c = PenduloFisica.posicionesCartesianas(estado.x, estado.params);
            var p1x = ox + c.x1 * esc, p1y = oy - c.y1 * esc;
            var p2x = ox + c.x2 * esc, p2y = oy - c.y2 * esc;

            var d1 = Math.hypot(pos.x - p1x, pos.y - p1y);
            var d2 = Math.hypot(pos.x - p2x, pos.y - p2y);

            if (d2 < 24) {
                estado.arrastre.activo = true;
                estado.arrastre.masa = 2;
                estado.corriendo = false;
                dom.btnPlayPause.textContent = 'Reanudar';
                dom.btnPlayPause.classList.remove('btn-primary');
            } else if (d1 < 24) {
                estado.arrastre.activo = true;
                estado.arrastre.masa = 1;
                estado.corriendo = false;
                dom.btnPlayPause.textContent = 'Reanudar';
                dom.btnPlayPause.classList.remove('btn-primary');
            }
        }

        function onPointerMove(evt) {
            if (!estado.arrastre.activo) return;
            var pos = getMousePos(evt);
            var ox = pos.cw * 0.5;
            var oy = pos.ch * 0.30;
            var esc = Math.min(pos.cw, pos.ch) * 0.28 / Math.max(0.1, (estado.params.l1 + estado.params.l2));

            var c = PenduloFisica.posicionesCartesianas(estado.x, estado.params);
            var p1x = ox + c.x1 * esc, p1y = oy - c.y1 * esc;

            if (estado.arrastre.masa === 1) {
                var dx1 = pos.x - ox;
                var dy1 = pos.y - oy;
                var newTh1 = Math.atan2(dx1, dy1);
                dom.sliderTh1.value = (newTh1 * 180 / Math.PI).toFixed(0);
                actualizarParametrosLabels();
                inicializarSimulacion();
            } else if (estado.arrastre.masa === 2) {
                var dx2 = pos.x - p1x;
                var dy2 = pos.y - p1y;
                var newTh2 = Math.atan2(dx2, dy2);
                dom.sliderTh2.value = (newTh2 * 180 / Math.PI).toFixed(0);
                actualizarParametrosLabels();
                inicializarSimulacion();
            }
        }

        function onPointerUp() {
            estado.arrastre.activo = false;
            estado.arrastre.masa = 0;
        }

        canvas.addEventListener('mousedown', onPointerDown);
        window.addEventListener('mousemove', onPointerMove);
        window.addEventListener('mouseup', onPointerUp);

        canvas.addEventListener('touchstart', onPointerDown, { passive: true });
        window.addEventListener('touchmove', onPointerMove, { passive: true });
        window.addEventListener('touchend', onPointerUp);
    }

    // =========================================================================
    // 16. MODO CAPTURA Y CONFIGURACIÓN POR URL (?tab=...&capture=1)
    // =========================================================================
    function aplicarParametrosUrl() {
        var params = new URLSearchParams(window.location.search);

        var tabParam = params.get('tab');
        if (tabParam && ['pendulo', 'gemelos', 'montecarlo', 'mapacaos'].indexOf(tabParam) !== -1) {
            cambiarPestaña(tabParam);
        }

        var themeParam = params.get('theme');
        if (themeParam === 'light') {
            temaActual = 'light';
            dom.body.classList.remove('theme-dark');
            dom.body.classList.add('theme-light');
            dom.themeIcon.textContent = '';
            dom.themeText.textContent = 'Claro';
        }

        var captureParam = params.get('capture');
        if (captureParam === '1') {
            // Calentar simulaciones para que las gráficas y estadísticas estén completas
            if (estado.tabActiva === 'montecarlo') {
                reiniciarMonteCarlo();
                for (var mcStep = 0; mcStep < 300; mcStep++) {
                    pasoMonteCarlo(0.05);
                }
            } else if (estado.tabActiva === 'gemelos') {
                for (var gStep = 0; gStep < 1400; gStep++) {
                    pasoSimulacion(estado.dtBase);
                    if (gStep % 3 === 0) registrarHistorialPunto();
                }
            } else {
                for (var step = 0; step < 800; step++) {
                    pasoSimulacion(estado.dtBase);
                    if (step % 4 === 0) registrarHistorialPunto();
                }
            }
            renderizarVistaActual();
        }
    }

    // =========================================================================
    // 17. ARRANQUE AL CARGAR EL DOM
    // =========================================================================
    function arrancarApp() {
        initDomReferences();
        inicializarSimulacion();
        inicializarMapaCaos();
        vincularEventos();
        aplicarParametrosUrl();
        renderizarVistaActual();

        // Iniciar loop continuo
        requestAnimationFrame(bucleAnimacion);

        // Exponer para verificación automatizada y pruebas
        if (typeof window !== 'undefined') {
            window.__simuladorApp = {
                estado: estado,
                dom: dom,
                reiniciarSimulacion: reiniciarSimulacion,
                inicializarSimulacion: inicializarSimulacion,
                aplicarPreset: aplicarPreset,
                cambiarPestaña: cambiarPestaña
            };
        }
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', arrancarApp);
    } else {
        arrancarApp();
    }

})();
