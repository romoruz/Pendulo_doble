/**
 * app.js - Controlador principal del Simulador de Péndulo Doble Conservativo.
 * 
 * Cumple con todas las especificaciones de diseño:
 * - Sin CDNs ni dependencias externas; funciona en file:// y GitHub Pages.
 * - Soporte táctil con Pointer Events (áreas >= 44px) para Safari iPad en horizontal.
 * - Integrador RK4 3/8 con dt fijo de 0.001 s y acumulador de tiempo real.
 * - Renderizado en Canvas con escala fija y radio proporcional a la raíz cúbica de m.
 * - 4 Modos: Normal, Gemelo (separación log10), Lineal (fantasma A) y Nube (Box-Muller con histograma).
 * - Barra de energía vertical con umbrales en orden físico real y badges de régimen.
 * - Mapa de caos FTLE progresivo generado en Web Worker vía Blob (compatible con file://).
 */

(function () {
    'use strict';

    // =========================================================================
    // CONSTANTES Y CONFIGURACIÓN INICIAL
    // =========================================================================
    const DT = 0.001; // Paso temporal fijo de física en segundos
    const L_MAX_SUM = 4.0; // l1_max + l2_max (2.0 + 2.0 m) para escala fija
    const COLOR_TH1 = '#1f4fd6'; // Azul institucional
    const COLOR_TH2 = '#d62728'; // Rojo institucional
    const COLOR_GHOST = '#777777';
    const COLOR_TWIN = '#e65100';

    // Estado físico del péndulo principal
    let params = PenduloFisica.crearParametros(1.0, 1.0, 1.0, 1.0, 9.81);
    let state = [Math.PI / 4, Math.PI / 4, 0.0, 0.0]; // [th1, th2, w1, w2]
    let stateInitial = [...state];
    let simTime = 0.0;
    let initialEnergy = PenduloFisica.energiaTotal(state, params);
    let flipsM1 = 0;
    let flipsM2 = 0;
    let k1Prev = PenduloFisica.indiceVuelta(state[0]);
    let k2Prev = PenduloFisica.indiceVuelta(state[1]);

    // Estado del bucle de animación
    let isPlaying = true;
    let speedMultiplier = 1.0;
    let showTrace = true;
    let tracePoints = []; // [{x, y, alpha}]
    const MAX_TRACE_POINTS = 350;

    // Modo activo ('normal', 'gemelo', 'lineal', 'nube')
    let currentMode = 'normal';

    // Estado del modo Gemelo
    let twinDelta = 1e-3; // rad
    let twinState = [state[0] + twinDelta, state[1], state[2], state[3]];
    let twinSeparationHistory = []; // [{t, logSep}]

    // Estado del modo Lineal
    let linInfo = PenduloFisica.linealizacion(params);
    let linearState = [...state];

    // Estado del modo Nube
    let cloudN = 200;
    let cloudSigma = 0.01;
    let cloudStates = []; // Array of [th1, th2, w1, w2]
    let cloudFlips = []; // Array of flips count per pendulum

    // Arrastre interactivo
    let draggingBob = 0; // 0: ninguno, 1: masa 1, 2: masa 2

    // =========================================================================
    // ELEMENTOS DEL DOM
    // =========================================================================
    const canvas = document.getElementById('sim-canvas');
    const ctx = canvas.getContext('2d');
    const canvasWrapper = document.getElementById('canvas-wrapper');

    const overlayPanel = document.getElementById('overlay-panel');
    const overlayTitle = document.getElementById('overlay-title');
    const overlayCanvas = document.getElementById('overlay-canvas');
    const overlayCtx = overlayCanvas.getContext('2d');

    const chaosCanvas = document.getElementById('chaos-canvas');
    const chaosCtx = chaosCanvas.getContext('2d');
    const chaosMarker = document.getElementById('chaos-target-marker');

    // Botones de transporte
    const btnPlayPause = document.getElementById('btn-play-pause');
    const btnStep = document.getElementById('btn-step');
    const btnReset = document.getElementById('btn-reset');
    const selectSpeed = document.getElementById('select-speed');
    const checkTrace = document.getElementById('check-trace');
    const btnPresentation = document.getElementById('btn-presentation');
    const dragHint = document.getElementById('drag-hint');

    // Métricas y datos de energía
    const labelDeltaE = document.getElementById('label-delta-e');
    const valEnergyTotal = document.getElementById('val-energy-total');
    const valEnergyKinetic = document.getElementById('val-energy-kinetic');
    const valEnergyPotential = document.getElementById('val-energy-potential');
    const valTime = document.getElementById('val-time');
    const valFlipsM2 = document.getElementById('val-flips-m2');
    const valFlipsM1 = document.getElementById('val-flips-m1');
    const regimeBadge = document.getElementById('regime-badge');

    // Barras de energía y umbrales
    const energyBarKinetic = document.getElementById('energy-bar-kinetic');
    const energyBarPotential = document.getElementById('energy-bar-potential');
    const lineEcM2 = document.getElementById('line-ec-m2');
    const lineEcM1 = document.getElementById('line-ec-m1');
    const lineEcAmbas = document.getElementById('line-ec-ambas');
    const lblEcM2 = document.getElementById('lbl-ec-m2');
    const lblEcM1 = document.getElementById('lbl-ec-m1');
    const lblEcAmbas = document.getElementById('lbl-ec-ambas');

    // Sliders
    const sliderL1 = document.getElementById('slider-l1');
    const sliderL2 = document.getElementById('slider-l2');
    const sliderM1 = document.getElementById('slider-m1');
    const sliderM2 = document.getElementById('slider-m2');
    const sliderTh1 = document.getElementById('slider-th1');
    const sliderTh2 = document.getElementById('slider-th2');
    const sliderW1 = document.getElementById('slider-w1');
    const sliderW2 = document.getElementById('slider-w2');

    const lblL1 = document.getElementById('lbl-l1');
    const lblL2 = document.getElementById('lbl-l2');
    const lblM1 = document.getElementById('lbl-m1');
    const lblM2 = document.getElementById('lbl-m2');
    const lblTh1 = document.getElementById('lbl-th1');
    const lblTh2 = document.getElementById('lbl-th2');
    const lblW1 = document.getElementById('lbl-w1');
    const lblW2 = document.getElementById('lbl-w2');
    const lblMu = document.getElementById('lbl-mu');
    const lblLambda = document.getElementById('lbl-lambda');

    const modeControlsSection = document.getElementById('mode-controls-section');
    const modeControlsTitle = document.getElementById('mode-controls-title');
    const modeControlsContent = document.getElementById('mode-controls-content');

    // =========================================================================
    // DIMENSIONAMIENTO DEL CANVAS
    // =========================================================================
    let canvasWidth = 0;
    let canvasHeight = 0;
    let pivotX = 0;
    let pivotY = 0;
    let scalePxPerM = 100;

    function resizeCanvas() {
        const rect = canvasWrapper.getBoundingClientRect();
        const dpr = window.devicePixelRatio || 1;
        canvasWidth = rect.width;
        canvasHeight = rect.height;

        canvas.width = Math.floor(canvasWidth * dpr);
        canvas.height = Math.floor(canvasHeight * dpr);
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

        pivotX = canvasWidth / 2;
        pivotY = canvasHeight * 0.28;

        // Escala fija en px/m calculada con la longitud máxima posible (L_MAX_SUM = 4.0 m)
        // Deja márgenes confortables arriba y abajo para visualización completa
        const usableHeight = canvasHeight * 0.65;
        scalePxPerM = Math.max(30, usableHeight / L_MAX_SUM);
    }
    window.addEventListener('resize', resizeCanvas);

    // =========================================================================
    // REINICIO Y ACTUALIZACIÓN DE ESTADO
    // =========================================================================
    function actualizarRatios() {
        const mu = params.m2 / params.m1;
        const lam = params.l2 / params.l1;
        lblMu.textContent = mu.toFixed(2);
        lblLambda.textContent = lam.toFixed(2);
    }

    function inicializarNube() {
        cloudStates = [];
        cloudFlips = [];
        // Box-Muller para ruido normal
        for (let i = 0; i < cloudN; i++) {
            const u1 = Math.max(1e-12, Math.random());
            const u2 = Math.random();
            const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
            const z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);

            cloudStates.push([
                stateInitial[0] + cloudSigma * z0,
                stateInitial[1] + cloudSigma * z1,
                stateInitial[2],
                stateInitial[3]
            ]);
            cloudFlips.push(0);
        }
    }

    function resetSimulation(keepInitial = true) {
        if (!keepInitial) {
            stateInitial = [
                parseFloat(sliderTh1.value) * Math.PI / 180.0,
                parseFloat(sliderTh2.value) * Math.PI / 180.0,
                parseFloat(sliderW1.value),
                parseFloat(sliderW2.value)
            ];
        }
        state = [...stateInitial];
        simTime = 0.0;
        initialEnergy = PenduloFisica.energiaTotal(state, params);
        flipsM1 = 0;
        flipsM2 = 0;
        k1Prev = PenduloFisica.indiceVuelta(state[0]);
        k2Prev = PenduloFisica.indiceVuelta(state[1]);
        tracePoints = [];

        // Modos secundarios
        twinState = [state[0] + twinDelta, state[1], state[2], state[3]];
        twinSeparationHistory = [];

        linInfo = PenduloFisica.linealizacion(params);
        linearState = [...state];

        if (currentMode === 'nube') {
            inicializarNube();
        }

        actualizarRatios();
        actualizarPosicionMarcadorCaos();
        actualizarPanelesMetricas();
    }

    // =========================================================================
    // SUB-PASO DE FÍSICA (RK4 3/8)
    // =========================================================================
    function stepPhysics() {
        // 1. Péndulo principal
        state = PenduloFisica.pasoRK4(state, DT, params);
        simTime += DT;

        // Conteo exacto de pasos por arriba
        const k1Now = PenduloFisica.indiceVuelta(state[0]);
        const k2Now = PenduloFisica.indiceVuelta(state[1]);
        flipsM1 += Math.abs(k1Now - k1Prev);
        flipsM2 += Math.abs(k2Now - k2Prev);
        k1Prev = k1Now;
        k2Prev = k2Now;

        // 2. Modos secundarios
        if (currentMode === 'gemelo') {
            twinState = PenduloFisica.pasoRK4(twinState, DT, params);
            // Distancia en el espacio fase ||x1 - x2||
            const dth1 = state[0] - twinState[0];
            const dth2 = state[1] - twinState[1];
            const dw1 = state[2] - twinState[2];
            const dw2 = state[3] - twinState[3];
            const dist = Math.sqrt(dth1 * dth1 + dth2 * dth2 + dw1 * dw1 + dw2 * dw2);
            const logDist = Math.log10(Math.max(1e-15, dist));
            if (twinSeparationHistory.length === 0 || simTime - twinSeparationHistory[twinSeparationHistory.length - 1].t >= 0.02) {
                twinSeparationHistory.push({ t: simTime, logSep: logDist });
                if (twinSeparationHistory.length > 500) twinSeparationHistory.shift();
            }
        } else if (currentMode === 'lineal') {
            linearState = PenduloFisica.pasoLinealRK4(linearState, DT, linInfo);
        } else if (currentMode === 'nube') {
            for (let i = 0; i < cloudStates.length; i++) {
                const sPrev = cloudStates[i];
                const sNext = PenduloFisica.pasoRK4(sPrev, DT, params);
                const kPrev = PenduloFisica.indiceVuelta(sPrev[1]);
                const kNext = PenduloFisica.indiceVuelta(sNext[1]);
                cloudFlips[i] += Math.abs(kNext - kPrev);
                cloudStates[i] = sNext;
            }
        }

        // Historial de la estela de la masa 2
        if (showTrace) {
            const pos = PenduloFisica.posicionesCartesianas(state, params);
            tracePoints.push({
                x: pivotX + pos.x2 * scalePxPerM,
                y: pivotY - pos.y2 * scalePxPerM,
                time: simTime
            });
            if (tracePoints.length > MAX_TRACE_POINTS) {
                tracePoints.shift();
            }
        }
    }

    // =========================================================================
    // ACTUALIZACIÓN DE INTERFAZ Y MÉTRICAS
    // =========================================================================
    function actualizarPanelesMetricas() {
        const E = PenduloFisica.energiaTotal(state, params);
        const T = PenduloFisica.energiaCinetica(state, params);
        const V = PenduloFisica.energiaPotencial(state, params);
        const deltaE = Math.abs(E - initialEnergy) / Math.max(1.0, Math.abs(initialEnergy));

        valEnergyTotal.textContent = E.toFixed(2) + ' J';
        valEnergyKinetic.textContent = T.toFixed(2) + ' J';
        valEnergyPotential.textContent = V.toFixed(2) + ' J';
        valTime.textContent = simTime.toFixed(2) + ' s';
        valFlipsM2.textContent = flipsM2.toString();
        valFlipsM1.textContent = flipsM1.toString();
        labelDeltaE.textContent = `|ΔE|/E₀: ${deltaE.toExponential(2)}`;

        // Umbrales reales según parámetros
        const u = PenduloFisica.umbralesEnergia(params);
        const ecM2 = u.Ec_m2;
        const ecM1 = u.Ec_m1;
        const ecAmbas = u.Ec_ambas;

        // Escala vertical de la barra de energía
        const maxDisplayE = Math.max(ecAmbas * 1.15, E * 1.1, 10.0);
        const pctKinetic = Math.min(100, Math.max(0, (T / maxDisplayE) * 100));
        const pctPotential = Math.min(100, Math.max(0, (V / maxDisplayE) * 100));

        energyBarKinetic.style.height = pctKinetic + '%';
        energyBarPotential.style.bottom = pctKinetic + '%';
        energyBarPotential.style.height = pctPotential + '%';

        // Posiciones físicas de las líneas de umbral
        const pctLineM2 = (ecM2 / maxDisplayE) * 100;
        const pctLineM1 = (ecM1 / maxDisplayE) * 100;
        const pctLineAmbas = (ecAmbas / maxDisplayE) * 100;

        lineEcM2.style.bottom = pctLineM2 + '%';
        lineEcM1.style.bottom = pctLineM1 + '%';
        lineEcAmbas.style.bottom = pctLineAmbas + '%';

        lblEcM2.textContent = `Ec_m2: ${ecM2.toFixed(1)} J`;
        lblEcM1.textContent = `Ec_m1: ${ecM1.toFixed(1)} J`;
        lblEcAmbas.textContent = `Ec_ambas: ${ecAmbas.toFixed(1)} J`;

        // Clasificación rigurosa del régimen energético
        const primerUmbral = Math.min(ecM2, ecM1);
        const segundoUmbral = Math.max(ecM2, ecM1);

        if (E < primerUmbral) {
            regimeBadge.style.color = '#006622';
            regimeBadge.style.borderColor = '#b3e6c2';
            regimeBadge.style.background = '#f0fff4';
            regimeBadge.textContent = 'E < Ec_min: Ninguna masa puede dar la vuelta completa (garantizado teóricamente por conservación de energía).';
        } else if (E < segundoUmbral) {
            regimeBadge.style.color = '#996600';
            regimeBadge.style.borderColor = '#ffe099';
            regimeBadge.style.background = '#fffdf0';
            const masaQueGira = (ecM2 < ecM1) ? 'La masa 2 (inferior)' : 'La masa 1 (superior)';
            regimeBadge.innerHTML = `<strong>${primerUmbral.toFixed(1)} J &le; E < ${segundoUmbral.toFixed(1)} J</strong>: ${masaQueGira} tiene energía suficiente para superar la vertical.<br><span style="font-size: 10.5px; color: #666;">Nota: Poder dar la vuelta no implica caos necesariamente.</span>`;
        } else {
            regimeBadge.style.color = '#990000';
            regimeBadge.style.borderColor = '#f5c6cb';
            regimeBadge.style.background = '#fff5f5';
            regimeBadge.innerHTML = `<strong>E &ge; ${segundoUmbral.toFixed(1)} J</strong>: Ambas masas poseen energía suficiente para rotar completamente sobre sus pivotes.`;
        }
    }

    // =========================================================================
    // DIBUJO EN CANVAS PRINCIPAL
    // =========================================================================
    function radioMasa(m) {
        // Radio proporcional a la raíz cúbica de m (acotado entre 8px y 25px)
        return Math.max(8.0, Math.min(25.0, 9.0 * Math.cbrt(m)));
    }

    function renderPendulo(st, p, colorVarillas, alphaVarillas, colorBob1, colorBob2, drawLabels = true) {
        const pos = PenduloFisica.posicionesCartesianas(st, p);
        const x1 = pivotX + pos.x1 * scalePxPerM;
        const y1 = pivotY - pos.y1 * scalePxPerM;
        const x2 = pivotX + pos.x2 * scalePxPerM;
        const y2 = pivotY - pos.y2 * scalePxPerM;

        ctx.save();
        ctx.globalAlpha = alphaVarillas;

        // Varillas
        ctx.lineWidth = 2.5;
        ctx.strokeStyle = colorVarillas;
        ctx.beginPath();
        ctx.moveTo(pivotX, pivotY);
        ctx.lineTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();

        // Masa 1
        const r1 = radioMasa(p.m1);
        ctx.fillStyle = colorBob1;
        ctx.beginPath();
        ctx.arc(x1, y1, r1, 0, 2.0 * Math.PI);
        ctx.fill();
        ctx.lineWidth = 1.2;
        ctx.strokeStyle = '#111111';
        ctx.stroke();

        // Masa 2
        const r2 = radioMasa(p.m2);
        ctx.fillStyle = colorBob2;
        ctx.beginPath();
        ctx.arc(x2, y2, r2, 0, 2.0 * Math.PI);
        ctx.fill();
        ctx.lineWidth = 1.2;
        ctx.strokeStyle = '#111111';
        ctx.stroke();

        // Etiquetas
        if (drawLabels) {
            ctx.fillStyle = '#111111';
            ctx.font = 'bold 11px var(--font-stack)';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText(`${p.m1.toFixed(1)} kg`, x1, y1);
            ctx.fillText(`${p.m2.toFixed(1)} kg`, x2, y2);
        }

        ctx.restore();
        return { x1, y1, x2, y2 };
    }

    function render() {
        ctx.clearRect(0, 0, canvasWidth, canvasHeight);

        // Guía del pivote
        ctx.fillStyle = '#111111';
        ctx.beginPath();
        ctx.arc(pivotX, pivotY, 4.5, 0, 2 * Math.PI);
        ctx.fill();

        // 1. Estela de la masa 2
        if (showTrace && tracePoints.length > 1) {
            ctx.lineWidth = 1.5;
            for (let i = 1; i < tracePoints.length; i++) {
                const alpha = (i / tracePoints.length) * 0.75;
                ctx.strokeStyle = `rgba(214, 39, 40, ${alpha})`;
                ctx.beginPath();
                ctx.moveTo(tracePoints[i - 1].x, tracePoints[i - 1].y);
                ctx.lineTo(tracePoints[i].x, tracePoints[i].y);
                ctx.stroke();
            }
        }

        // 2. Modos secundarios (debajo del péndulo principal)
        if (currentMode === 'gemelo') {
            renderPendulo(twinState, params, COLOR_TWIN, 0.45, COLOR_TWIN, COLOR_TWIN, false);
        } else if (currentMode === 'lineal') {
            renderPendulo(linearState, params, COLOR_GHOST, 0.4, COLOR_GHOST, COLOR_GHOST, false);
        } else if (currentMode === 'nube') {
            // Dibujar masa 2 de cada miembro de la nube
            ctx.fillStyle = 'rgba(214, 39, 40, 0.35)';
            for (let i = 0; i < cloudStates.length; i++) {
                const cPos = PenduloFisica.posicionesCartesianas(cloudStates[i], params);
                const cx2 = pivotX + cPos.x2 * scalePxPerM;
                const cy2 = pivotY - cPos.y2 * scalePxPerM;
                ctx.beginPath();
                ctx.arc(cx2, cy2, 3.0, 0, 2 * Math.PI);
                ctx.fill();
            }
        }

        // 3. Péndulo principal
        const coords = renderPendulo(state, params, '#111111', 1.0, COLOR_TH1, COLOR_TH2, true);

        // 4. Líneas guía verticales punteadas y arcos de ángulos
        const len1 = params.l1 * scalePxPerM;
        const len2 = params.l2 * scalePxPerM;

        ctx.save();
        ctx.setLineDash([3, 4]);
        ctx.strokeStyle = '#888888';
        ctx.lineWidth = 1.0;

        // Línea vertical por el pivote
        ctx.beginPath();
        ctx.moveTo(pivotX, pivotY);
        ctx.lineTo(pivotX, pivotY + len1 * 1.15);
        ctx.stroke();

        // Línea vertical por la masa 1
        ctx.beginPath();
        ctx.moveTo(coords.x1, coords.y1);
        ctx.lineTo(coords.x1, coords.y1 + len2 * 1.15);
        ctx.stroke();
        ctx.restore();

        // Arcos de ángulo
        // θ1 (azul)
        const arcR1 = Math.min(36, len1 * 0.45);
        ctx.save();
        ctx.strokeStyle = COLOR_TH1;
        ctx.lineWidth = 2.0;
        ctx.beginPath();
        const startAng1 = Math.PI / 2; // Vertical hacia abajo
        const endAng1 = Math.PI / 2 + state[0];
        ctx.arc(pivotX, pivotY, arcR1, startAng1, endAng1, state[0] < 0);
        ctx.stroke();

        // Etiqueta θ1
        ctx.fillStyle = COLOR_TH1;
        ctx.font = 'bold 11px var(--font-stack)';
        const midAng1 = Math.PI / 2 + state[0] / 2;
        ctx.fillText('θ₁', pivotX + (arcR1 + 10) * Math.cos(midAng1), pivotY + (arcR1 + 10) * Math.sin(midAng1));

        // θ2 (rojo)
        const arcR2 = Math.min(30, len2 * 0.45);
        ctx.strokeStyle = COLOR_TH2;
        ctx.lineWidth = 2.0;
        ctx.beginPath();
        const startAng2 = Math.PI / 2;
        const endAng2 = Math.PI / 2 + state[1];
        ctx.arc(coords.x1, coords.y1, arcR2, startAng2, endAng2, state[1] < 0);
        ctx.stroke();

        // Etiqueta θ2
        ctx.fillStyle = COLOR_TH2;
        const midAng2 = Math.PI / 2 + state[1] / 2;
        ctx.fillText('θ₂', coords.x1 + (arcR2 + 10) * Math.cos(midAng2), coords.y1 + (arcR2 + 10) * Math.sin(midAng2));
        ctx.restore();

        // Renderizado del panel de diagnóstico superpuesto
        renderOverlay();
    }

    // =========================================================================
    // RENDERIZADO DEL PANEL SUPERPUESTO (GEMELO / NUBE)
    // =========================================================================
    function renderOverlay() {
        if (currentMode === 'gemelo') {
            overlayPanel.classList.add('visible');
            overlayTitle.textContent = `Sensibilidad Caótica (δ₀ = ${twinDelta.toExponential(0)})`;

            overlayCtx.clearRect(0, 0, 250, 120);
            overlayCtx.fillStyle = '#555555';
            overlayCtx.font = '9px var(--font-mono)';
            overlayCtx.fillText('log10 ||Δx||', 6, 12);
            overlayCtx.fillText('t [s]', 220, 115);

            if (twinSeparationHistory.length > 1) {
                overlayCtx.strokeStyle = COLOR_TWIN;
                overlayCtx.lineWidth = 1.5;
                overlayCtx.beginPath();
                const minLog = -8;
                const maxLog = 1;
                const maxT = Math.max(10.0, simTime);

                for (let i = 0; i < twinSeparationHistory.length; i++) {
                    const pt = twinSeparationHistory[i];
                    const px = 20 + (pt.t / maxT) * 220;
                    const py = 100 - ((pt.logSep - minLog) / (maxLog - minLog)) * 85;
                    if (i === 0) overlayCtx.moveTo(px, py);
                    else overlayCtx.lineTo(px, py);
                }
                overlayCtx.stroke();
            }
        } else if (currentMode === 'nube') {
            overlayPanel.classList.add('visible');
            overlayTitle.textContent = `Histograma Vueltas Nube (N = ${cloudN})`;

            overlayCtx.clearRect(0, 0, 250, 120);
            // Conteo de frecuencias de vueltas
            const bins = {};
            let maxCount = 1;
            let maxFlip = 0;
            for (let f of cloudFlips) {
                bins[f] = (bins[f] || 0) + 1;
                if (bins[f] > maxCount) maxCount = bins[f];
                if (f > maxFlip) maxFlip = f;
            }

            const rangeX = Math.max(10, maxFlip + 2);
            const barW = Math.max(3, 200 / rangeX);

            overlayCtx.fillStyle = COLOR_TH2;
            for (let f = 0; f <= rangeX; f++) {
                const count = bins[f] || 0;
                const h = (count / maxCount) * 85;
                const bx = 20 + f * barW;
                const by = 105 - h;
                overlayCtx.fillRect(bx, by, barW - 1, h);
            }

            overlayCtx.fillStyle = '#444444';
            overlayCtx.font = '9px var(--font-mono)';
            overlayCtx.fillText(`Vueltas m₂ (Max: ${maxFlip})`, 20, 118);
        } else {
            overlayPanel.classList.remove('visible');
        }
    }

    // =========================================================================
    // MAPA DE CAOS (FTLE EN WEB WORKER VIA BLOB)
    // =========================================================================
    let workerFTLE = null;

    function iniciarWorkerFTLE() {
        // Código fuente del Web Worker autocontenido (compatible con file://)
        const workerSource = `
            self.onmessage = function (e) {
                const { res, params, tMax, dt, pass } = e.data;
                const m1 = params.m1, m2 = params.m2, l1 = params.l1, l2 = params.l2, g = params.g;
                const d0 = 1e-7;
                const renormSteps = 10;
                const dt_3 = dt / 3.0, dt_8 = dt / 8.0;

                function f(th1, th2, w1, w2) {
                    const delta = th1 - th2;
                    const sinD = Math.sin(delta), cosD = Math.cos(delta);
                    const m11 = (m1 + m2) * l1 * l1;
                    const m12 = m2 * l1 * l2 * cosD;
                    const m22 = m2 * l2 * l2;
                    const det = m2 * (l1 * l1) * (l2 * l2) * (m1 + m2 * sinD * sinD);
                    const f1 = -m2 * l1 * l2 * sinD * (w2 * w2) - (m1 + m2) * g * l1 * Math.sin(th1);
                    const f2 = m2 * l1 * l2 * sinD * (w1 * w1) - m2 * g * l2 * Math.sin(th2);
                    return [w1, w2, (f1 * m22 - f2 * m12) / det, (m11 * f2 - m12 * f1) / det];
                }

                function rk(s) {
                    const k1 = f(s[0], s[1], s[2], s[3]);
                    const k2 = f(s[0]+dt_3*k1[0], s[1]+dt_3*k1[1], s[2]+dt_3*k1[2], s[3]+dt_3*k1[3]);
                    const k3 = f(s[0]-dt_3*k1[0]+dt*k2[0], s[1]-dt_3*k1[1]+dt*k2[1], s[2]-dt_3*k1[2]+dt*k2[2], s[3]-dt_3*k1[3]+dt*k2[3]);
                    const k4 = f(s[0]+dt*k1[0]-dt*k2[0]+dt*k3[0], s[1]+dt*k1[1]-dt*k2[1]+dt*k3[1], s[2]+dt*k1[2]-dt*k2[2]+dt*k3[2], s[3]+dt*k1[3]-dt*k2[3]+dt*k3[3]);
                    return [
                        s[0] + dt_8*(k1[0]+3*k2[0]+3*k3[0]+k4[0]),
                        s[1] + dt_8*(k1[1]+3*k2[1]+3*k3[1]+k4[1]),
                        s[2] + dt_8*(k1[2]+3*k2[2]+3*k3[2]+k4[2]),
                        s[3] + dt_8*(k1[3]+3*k2[3]+3*k3[3]+k4[3])
                    ];
                }

                const nSteps = Math.floor(tMax / dt);
                const nRenorms = Math.floor(nSteps / renormSteps);
                const grid = new Float32Array(res * res);

                for (let j = 0; j < res; j++) {
                    const th2 = -Math.PI + (2 * Math.PI * j) / (res - 1);
                    for (let i = 0; i < res; i++) {
                        const th1 = -Math.PI + (2 * Math.PI * i) / (res - 1);
                        let sBase = [th1, th2, 0.0, 0.0];
                        let sPert = [th1 + d0, th2, 0.0, 0.0];
                        let sumLog = 0.0;

                        for (let r = 0; r < nRenorms; r++) {
                            for (let step = 0; step < renormSteps; step++) {
                                sBase = rk(sBase);
                                sPert = rk(sPert);
                            }
                            const d0_x = sPert[0] - sBase[0];
                            const d1_x = sPert[1] - sBase[1];
                            const d2_x = sPert[2] - sBase[2];
                            const d3_x = sPert[3] - sBase[3];
                            const dist = Math.sqrt(d0_x*d0_x + d1_x*d1_x + d2_x*d2_x + d3_x*d3_x);
                            sumLog += Math.log(dist / d0);
                            const factor = d0 / dist;
                            sPert = [
                                sBase[0] + d0_x * factor,
                                sBase[1] + d1_x * factor,
                                sBase[2] + d2_x * factor,
                                sBase[3] + d3_x * factor
                            ];
                        }
                        grid[j * res + i] = sumLog / tMax;
                    }
                }
                self.postMessage({ grid, res, pass });
            };
        `;

        const blob = new Blob([workerSource], { type: 'application/javascript' });
        const blobUrl = URL.createObjectURL(blob);
        workerFTLE = new Worker(blobUrl);

        workerFTLE.onmessage = function (e) {
            const { grid, res, pass } = e.data;
            dibujarMapaFTLE(grid, res);
            // Resolución progresiva: 32 -> 64 -> 96
            if (pass === 1) {
                workerFTLE.postMessage({ res: 64, params, tMax: 15.0, dt: 0.015, pass: 2 });
            } else if (pass === 2) {
                workerFTLE.postMessage({ res: 96, params, tMax: 20.0, dt: 0.01, pass: 3 });
            }
        };
    }

    function dispararCalculoMapaCaos() {
        if (!workerFTLE) {
            iniciarWorkerFTLE();
        }
        // Inicia primer pase rápido de 32x32
        workerFTLE.postMessage({ res: 32, params, tMax: 10.0, dt: 0.02, pass: 1 });
    }

    // Mapa de colores inferno simplificado y nítido
    function colorInferno(v) {
        // v en [0, 1]
        const t = Math.max(0.0, Math.min(1.0, v));
        const r = Math.floor(Math.min(255, t < 0.5 ? 2 * t * 180 : 180 + 2 * (t - 0.5) * 75));
        const g = Math.floor(Math.min(255, t < 0.6 ? Math.pow(t, 2) * 160 : 60 + (t - 0.6) * 450));
        const b = Math.floor(Math.min(255, t < 0.3 ? t * 300 : Math.max(0, 100 - (t - 0.3) * 150)));
        return [r, g, b];
    }

    function dibujarMapaFTLE(grid, res) {
        const w = chaosCanvas.width;
        const h = chaosCanvas.height;
        const imgData = chaosCtx.createImageData(w, h);
        const data = imgData.data;

        const maxFTLE = 0.8;

        for (let py = 0; py < h; py++) {
            const j = Math.floor(((h - 1 - py) / (h - 1)) * (res - 1));
            for (let px = 0; px < w; px++) {
                const i = Math.floor((px / (w - 1)) * (res - 1));
                const val = grid[j * res + i];
                const normVal = Math.max(0.0, Math.min(1.0, val / maxFTLE));
                const [r, g, b] = colorInferno(normVal);

                const pIndex = (py * w + px) * 4;
                data[pIndex] = r;
                data[pIndex + 1] = g;
                data[pIndex + 2] = b;
                data[pIndex + 3] = 255;
            }
        }
        chaosCtx.putImageData(imgData, 0, 0);

        // Contornos analíticos de energía potencial Ec_m2 (cian) y Ec_m1 (negro/blanco)
        const u = PenduloFisica.umbralesEnergia(params);
        chaosCtx.save();
        chaosCtx.lineWidth = 1.6;

        // V = (m1+m2)*g*l1*(1-cos(th1)) + m2*g*l2*(1-cos(th2))
        function drawContour(targetV, strokeColor, isDashed) {
            chaosCtx.strokeStyle = strokeColor;
            chaosCtx.setLineDash(isDashed ? [3, 3] : []);
            chaosCtx.beginPath();
            let started = false;

            // Recorrer el borde de la curva implícita
            for (let px = 0; px < w; px++) {
                const th1 = -Math.PI + (2 * Math.PI * px) / (w - 1);
                const v1 = (params.m1 + params.m2) * params.g * params.l1 * (1.0 - Math.cos(th1));
                const remV = targetV - v1;
                const cosTh2 = 1.0 - remV / (params.m2 * params.g * params.l2);

                if (cosTh2 >= -1.0 && cosTh2 <= 1.0) {
                    const th2_pos = Math.acos(cosTh2);
                    const py_pos = (1.0 - (th2_pos + Math.PI) / (2 * Math.PI)) * (h - 1);
                    const py_neg = (1.0 - (-th2_pos + Math.PI) / (2 * Math.PI)) * (h - 1);

                    // Puntos del contorno
                    chaosCtx.fillRect(px, py_pos, 1.2, 1.2);
                    chaosCtx.fillRect(px, py_neg, 1.2, 1.2);
                }
            }
        }

        drawContour(u.Ec_m2, '#00ffcc', false); // Ec_m2 en cian sólido
        drawContour(u.Ec_m1, '#ffffff', true);  // Ec_m1 en blanco punteado
        chaosCtx.restore();
    }

    function actualizarPosicionMarcadorCaos() {
        // Mapea stateInitial[0] y stateInitial[1] a coordenadas en chaos-map-canvas-wrapper
        const normX = (stateInitial[0] + Math.PI) / (2.0 * Math.PI);
        const normY = 1.0 - (stateInitial[1] + Math.PI) / (2.0 * Math.PI);
        chaosMarker.style.left = (normX * 100) + '%';
        chaosMarker.style.top = (normY * 100) + '%';
    }

    // Interacción al hacer click / touch en el mapa de caos
    chaosCanvas.parentElement.addEventListener('pointerdown', (e) => {
        const rect = chaosCanvas.getBoundingClientRect();
        const px = e.clientX - rect.left;
        const py = e.clientY - rect.top;
        const normX = Math.max(0.0, Math.min(1.0, px / rect.width));
        const normY = Math.max(0.0, Math.min(1.0, py / rect.height));

        const newTh1 = -Math.PI + normX * 2.0 * Math.PI;
        const newTh2 = Math.PI - normY * 2.0 * Math.PI;

        sliderTh1.value = Math.round(newTh1 * 180 / Math.PI);
        sliderTh2.value = Math.round(newTh2 * 180 / Math.PI);
        sliderW1.value = 0;
        sliderW2.value = 0;

        lblTh1.textContent = sliderTh1.value + '°';
        lblTh2.textContent = sliderTh2.value + '°';
        lblW1.textContent = '0.0 rad/s';
        lblW2.textContent = '0.0 rad/s';

        resetSimulation(false);
    });

    // =========================================================================
    // ARRASTRE DIRECTO DE LAS MASAS CON POINTER EVENTS
    // =========================================================================
    canvasWrapper.addEventListener('pointerdown', (e) => {
        if (isPlaying) return; // Solo interactivo en pausa

        const rect = canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        const pos = PenduloFisica.posicionesCartesianas(state, params);
        const x1 = pivotX + pos.x1 * scalePxPerM;
        const y1 = pivotY - pos.y1 * scalePxPerM;
        const x2 = pivotX + pos.x2 * scalePxPerM;
        const y2 = pivotY - pos.y2 * scalePxPerM;

        const d2 = Math.hypot(mx - x2, my - y2);
        const d1 = Math.hypot(mx - x1, my - y1);

        if (d2 <= Math.max(30, radioMasa(params.m2) + 15)) {
            draggingBob = 2;
            canvasWrapper.setPointerCapture(e.pointerId);
        } else if (d1 <= Math.max(30, radioMasa(params.m1) + 15)) {
            draggingBob = 1;
            canvasWrapper.setPointerCapture(e.pointerId);
        }
    });

    canvasWrapper.addEventListener('pointermove', (e) => {
        if (!draggingBob) return;
        const rect = canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        if (draggingBob === 1) {
            // Ajustar theta1
            const dx = mx - pivotX;
            const dy = my - pivotY;
            const angle = Math.atan2(dx, dy); // Vertical hacia abajo es 0
            state[0] = angle;
            state[2] = 0.0;
            sliderTh1.value = Math.round(angle * 180 / Math.PI);
            lblTh1.textContent = sliderTh1.value + '°';
        } else if (draggingBob === 2) {
            // Ajustar theta2 relativo a masa 1
            const pos = PenduloFisica.posicionesCartesianas(state, params);
            const x1 = pivotX + pos.x1 * scalePxPerM;
            const y1 = pivotY - pos.y1 * scalePxPerM;
            const dx = mx - x1;
            const dy = my - y1;
            const angle = Math.atan2(dx, dy);
            state[1] = angle;
            state[3] = 0.0;
            sliderTh2.value = Math.round(angle * 180 / Math.PI);
            lblTh2.textContent = sliderTh2.value + '°';
        }

        stateInitial = [...state];
        initialEnergy = PenduloFisica.energiaTotal(state, params);
        actualizarPanelesMetricas();
        actualizarPosicionMarcadorCaos();
        render();
    });

    function stopDragging(e) {
        if (draggingBob) {
            draggingBob = 0;
            try { canvasWrapper.releasePointerCapture(e.pointerId); } catch (_) {}
        }
    }
    canvasWrapper.addEventListener('pointerup', stopDragging);
    canvasWrapper.addEventListener('pointercancel', stopDragging);

    // =========================================================================
    // CONTROLES DE TRANSPORTE Y MODOS
    // =========================================================================
    btnPlayPause.addEventListener('click', () => {
        isPlaying = !isPlaying;
        btnPlayPause.textContent = isPlaying ? 'Pausa' : 'Iniciar';
        dragHint.style.display = isPlaying ? 'none' : 'block';
    });

    btnStep.addEventListener('click', () => {
        if (isPlaying) {
            isPlaying = false;
            btnPlayPause.textContent = 'Iniciar';
            dragHint.style.display = 'block';
        }
        for (let i = 0; i < 20; i++) {
            stepPhysics();
        }
        actualizarPanelesMetricas();
        render();
    });

    btnReset.addEventListener('click', () => {
        resetSimulation(true);
    });

    selectSpeed.addEventListener('change', () => {
        speedMultiplier = parseFloat(selectSpeed.value);
    });

    checkTrace.addEventListener('change', () => {
        showTrace = checkTrace.checked;
        if (!showTrace) tracePoints = [];
    });

    btnPresentation.addEventListener('click', () => {
        document.body.classList.toggle('presentation-mode');
        const isPres = document.body.classList.contains('presentation-mode');
        btnPresentation.textContent = isPres ? 'Salir Presentación' : 'Presentación';
        setTimeout(resizeCanvas, 50);
    });

    // Pestañas de Modos
    document.querySelectorAll('.mode-tab-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.mode-tab-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            currentMode = btn.getAttribute('data-mode');
            configurarControlesModo();
            resetSimulation(true);
        });
    });

    function configurarControlesModo() {
        if (currentMode === 'normal') {
            modeControlsSection.style.display = 'none';
        } else if (currentMode === 'gemelo') {
            modeControlsSection.style.display = 'block';
            modeControlsTitle.textContent = 'Controles del Modo Gemelo';
            modeControlsContent.innerHTML = `
                <div class="control-item">
                    <div class="control-label-row">
                        <span>Separación inicial &delta;<sub>0</sub>:</span>
                        <span id="lbl-twin-delta">${twinDelta.toExponential(0)} rad</span>
                    </div>
                    <input type="range" id="slider-twin-delta" class="control-slider" min="-5" max="-2" step="0.5" value="-3">
                </div>
            `;
            const slDelta = document.getElementById('slider-twin-delta');
            slDelta.addEventListener('input', () => {
                twinDelta = Math.pow(10, parseFloat(slDelta.value));
                document.getElementById('lbl-twin-delta').textContent = twinDelta.toExponential(0) + ' rad';
                resetSimulation(true);
            });
        } else if (currentMode === 'lineal') {
            modeControlsSection.style.display = 'block';
            modeControlsTitle.textContent = 'Modo Lineal (Fantasma At)';
            modeControlsContent.innerHTML = `
                <p style="font-size: 11.5px; color: var(--text-muted); line-height: 1.4;">
                    El péndulo gris integra la linealización <strong>dx/dt = A x</strong>.<br>
                    Para ángulos chicos (&le; 5°) ambos se mueven idénticos en modos normales (&omega;<sub>1</sub> = ${linInfo.omega1.toFixed(2)}, &omega;<sub>2</sub> = ${linInfo.omega2.toFixed(2)} rad/s).
                </p>
            `;
        } else if (currentMode === 'nube') {
            modeControlsSection.style.display = 'block';
            modeControlsTitle.textContent = 'Controles del Modo Nube';
            modeControlsContent.innerHTML = `
                <div class="control-item">
                    <div class="control-label-row">
                        <span>Péndulos N:</span>
                        <span id="lbl-cloud-n">${cloudN}</span>
                    </div>
                    <input type="range" id="slider-cloud-n" class="control-slider" min="50" max="500" step="25" value="${cloudN}">
                </div>
                <div class="control-item" style="margin-top: 8px;">
                    <div class="control-label-row">
                        <span>Dispersión &sigma;:</span>
                        <span id="lbl-cloud-sigma">${cloudSigma.toExponential(1)} rad</span>
                    </div>
                    <input type="range" id="slider-cloud-sigma" class="control-slider" min="-4" max="-1.3" step="0.1" value="-2">
                </div>
            `;
            const slN = document.getElementById('slider-cloud-n');
            const slSig = document.getElementById('slider-cloud-sigma');

            slN.addEventListener('input', () => {
                cloudN = parseInt(slN.value, 10);
                document.getElementById('lbl-cloud-n').textContent = cloudN;
                inicializarNube();
            });
            slSig.addEventListener('input', () => {
                cloudSigma = Math.pow(10, parseFloat(slSig.value));
                document.getElementById('lbl-cloud-sigma').textContent = cloudSigma.toExponential(1) + ' rad';
                inicializarNube();
            });
        }
    }

    // =========================================================================
    // SLIDERS DE PARÁMETROS Y CONDICIONES INICIALES
    // =========================================================================
    function vincularSlider(slider, label, formateador, alCambiar) {
        slider.addEventListener('input', () => {
            const val = parseFloat(slider.value);
            label.textContent = formateador(val);
            alCambiar(val);
            resetSimulation(false);
        });
    }

    vincularSlider(sliderL1, lblL1, v => v.toFixed(2) + ' m', v => {
        params.l1 = v;
        linInfo = PenduloFisica.linealizacion(params);
        dispararCalculoMapaCaos();
    });
    vincularSlider(sliderL2, lblL2, v => v.toFixed(2) + ' m', v => {
        params.l2 = v;
        linInfo = PenduloFisica.linealizacion(params);
        dispararCalculoMapaCaos();
    });
    vincularSlider(sliderM1, lblM1, v => v.toFixed(1) + ' kg', v => {
        params.m1 = v;
        linInfo = PenduloFisica.linealizacion(params);
        dispararCalculoMapaCaos();
    });
    vincularSlider(sliderM2, lblM2, v => v.toFixed(1) + ' kg', v => {
        params.m2 = v;
        linInfo = PenduloFisica.linealizacion(params);
        dispararCalculoMapaCaos();
    });

    vincularSlider(sliderTh1, lblTh1, v => Math.round(v) + '°', () => {});
    vincularSlider(sliderTh2, lblTh2, v => Math.round(v) + '°', () => {});
    vincularSlider(sliderW1, lblW1, v => v.toFixed(1) + ' rad/s', () => {});
    vincularSlider(sliderW2, lblW2, v => v.toFixed(1) + ' rad/s', () => {});

    document.getElementById('btn-recalc-map').addEventListener('click', dispararCalculoMapaCaos);

    // =========================================================================
    // PRESETS
    // =========================================================================
    document.querySelectorAll('.btn-preset').forEach(btn => {
        btn.addEventListener('click', () => {
            const key = btn.getAttribute('data-preset');
            const p = window.PRESETS_PENDULO && window.PRESETS_PENDULO[key];
            if (!p) return;

            // Parámetros físicos
            params = PenduloFisica.crearParametros(p.params.m1, p.params.m2, p.params.l1, p.params.l2, p.params.g);
            sliderL1.value = params.l1;
            sliderL2.value = params.l2;
            sliderM1.value = params.m1;
            sliderM2.value = params.m2;
            lblL1.textContent = params.l1.toFixed(2) + ' m';
            lblL2.textContent = params.l2.toFixed(2) + ' m';
            lblM1.textContent = params.m1.toFixed(1) + ' kg';
            lblM2.textContent = params.m2.toFixed(1) + ' kg';

            // Ángulos y velocidades
            const deg1 = Math.round(p.x0[0] * 180 / Math.PI);
            const deg2 = Math.round(p.x0[1] * 180 / Math.PI);
            sliderTh1.value = deg1;
            sliderTh2.value = deg2;
            sliderW1.value = p.x0[2];
            sliderW2.value = p.x0[3];
            lblTh1.textContent = deg1 + '°';
            lblTh2.textContent = deg2 + '°';
            lblW1.textContent = p.x0[2].toFixed(1) + ' rad/s';
            lblW2.textContent = p.x0[3].toFixed(1) + ' rad/s';

            stateInitial = [...p.x0];
            linInfo = PenduloFisica.linealizacion(params);
            dispararCalculoMapaCaos();
            resetSimulation(true);
        });
    });

    // =========================================================================
    // BUCLE PRINCIPAL DE ANIMACIÓN
    // =========================================================================
    let lastTime = performance.now();
    let accumulator = 0.0;

    function animLoop(now) {
        let frameElapsed = (now - lastTime) / 1000.0;
        lastTime = now;
        if (frameElapsed > 0.1) frameElapsed = 0.1; // Limitar lag de ventana inactiva

        if (isPlaying) {
            accumulator += frameElapsed * speedMultiplier;
            let subSteps = 0;
            const MAX_STEPS = 120;
            while (accumulator >= DT && subSteps < MAX_STEPS) {
                stepPhysics();
                accumulator -= DT;
                subSteps++;
            }
            actualizarPanelesMetricas();
        }

        render();
        requestAnimationFrame(animLoop);
    }

    // Inicialización al cargar la ventana
    window.addEventListener('DOMContentLoaded', () => {
        resizeCanvas();
        resetSimulation(true);
        dispararCalculoMapaCaos();
        dragHint.style.display = 'none';
        requestAnimationFrame(animLoop);
    });

})();
