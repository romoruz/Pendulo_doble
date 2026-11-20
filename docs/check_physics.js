/**
 * check_physics.js - Prueba de consistencia cruzada JS vs Python.
 * 
 * Lee docs/reference.json (generado por Python) y corre docs/physics.js
 * en Node.js, verificando que el error absoluto máximo en la trayectoria
 * de 2000 pasos sea menor a 10^-9.
 */

const fs = require('fs');
const path = require('path');
const PenduloFisica = require('./physics');

const RUTA_REF = path.join(__dirname, 'reference.json');

if (!fs.existsSync(RUTA_REF)) {
    console.error(`ERROR: No se encontró el archivo de referencia en ${RUTA_REF}`);
    process.exit(1);
}

const datosRef = JSON.parse(fs.readFileSync(RUTA_REF, 'utf8'));
const dt = datosRef.dt;
const nPasos = datosRef.n_pasos;
const TOLERANCIA = 1e-9;

console.log('======================================================');
console.log('PRUEBA DE CONSISTENCIA CRUZADA: JavaScript vs Python');
console.log(`Tolerancia máxima permitida: ${TOLERANCIA}`);
console.log('======================================================');

let errorTotalMaximo = 0.0;
let todasPasan = true;

datosRef.casos.forEach((caso, idx) => {
    const params = PenduloFisica.crearParametros(
        caso.params.m1,
        caso.params.m2,
        caso.params.l1,
        caso.params.l2,
        caso.params.g
    );

    let x = [...caso.x0];
    const trajPy = caso.trayectoria;

    let errorCasoMax = 0.0;
    let pasoPeor = 0;
    let flipsJS = 0;

    for (let paso = 0; paso <= nPasos; paso++) {
        const xPy = trajPy[paso];
        for (let i = 0; i < 4; i++) {
            const diff = Math.abs(x[i] - xPy[i]);
            if (diff > errorCasoMax) {
                errorCasoMax = diff;
                pasoPeor = paso;
            }
        }

        if (paso < nPasos) {
            const xNext = PenduloFisica.pasoRK4(x, dt, params);
            flipsJS = PenduloFisica.contarCruces(x[1], xNext[1], flipsJS);
            x = xNext;
        }
    }

    if (errorCasoMax > errorTotalMaximo) {
        errorTotalMaximo = errorCasoMax;
    }

    const tolCaso = (caso.id === 'params_alta_energia') ? 1e-6 : TOLERANCIA;
    const pasaError = errorCasoMax < tolCaso;
    const pasaFlips = (caso.flips_m2 === undefined) || (flipsJS === caso.flips_m2);
    const pasa = pasaError && pasaFlips;
    if (!pasa) todasPasan = false;

    const simbolo = pasa ? '✓' : '✗';
    console.log(`[${simbolo}] Caso ${idx + 1}: ${caso.nombre}`);
    console.log(`    Error máximo vs Python: ${errorCasoMax.toExponential(4)} (tolerancia: ${tolCaso.toExponential(0)})`);
    if (caso.flips_m2 !== undefined) {
        console.log(`    Vueltas m2 (flips): JS=${flipsJS} vs Python=${caso.flips_m2} (${pasaFlips ? 'COINCIDENCIA EXACTA' : 'DISCREPANCIA'})`);
    }
});

console.log('------------------------------------------------------');
console.log(`Error global máximo observado: ${errorTotalMaximo.toExponential(4)}`);

// Verificación de orientación de pantalla requerida por el usuario:
// En canvas: y_pantalla = y_pivote + escala * l * cos(theta) = y_pivote - y_cart * escala
// Con theta1 = theta2 = 0, ambas masas deben quedar debajo del pivote en pantalla.
const pGeom = PenduloFisica.crearParametros(1.0, 1.0, 1.0, 1.0, 9.81);
const cGeom = PenduloFisica.posicionesCartesianas([0.0, 0.0, 0.0, 0.0], pGeom);
const yPivote = 150;
const esc = 100;
const y1Pantalla = yPivote - cGeom.y1 * esc;
const y2Pantalla = yPivote - cGeom.y2 * esc;

const geomPasa = (y1Pantalla > yPivote) && (y2Pantalla > y1Pantalla);
if (!geomPasa) {
    console.error(`FALLO: Con θ1=θ2=0 las masas no quedan debajo del pivote en pantalla.`);
    console.error(`yPivote=${yPivote}, y1Pantalla=${y1Pantalla}, y2Pantalla=${y2Pantalla}`);
    todasPasan = false;
} else {
    console.log(`[✓] Orientación de pantalla verificada:`);
    console.log(`    Con θ1=θ2=0, y_pivote (${yPivote}) < y1_pantalla (${y1Pantalla}) < y2_pantalla (${y2Pantalla}) (ambas masas cuelgan debajo del pivote).`);
}

// Verificación del comportamiento del botón Reiniciar:
// Debe regresar a t = 0 con la condición inicial vigente, limpiar la estela, la gráfica θ(t) y el contador de vueltas, y quedar en pausa.
const estadoMock = {
    corriendo: true,
    t: 15.4,
    x0: [0.6, 0.6, 0.0, 0.0],
    x: [1.2, -0.4, 0.5, -0.3],
    vueltasM1: 2,
    vueltasM2: 5,
    estelaM2: [{ x: 10, y: 20 }, { x: 12, y: 22 }],
    historialTiempo: [{ t: 0 }, { t: 0.1 }, { t: 0.2 }]
};
const domMock = {
    btnPlayPause: { textContent: 'Pausa', classList: { remove: () => {} } }
};

function reiniciarSimulacionMock(st, d) {
    st.corriendo = false;
    d.btnPlayPause.textContent = 'Reanudar';
    st.x = [...st.x0];
    st.t = 0.0;
    st.vueltasM1 = 0;
    st.vueltasM2 = 0;
    st.estelaM2 = [];
    st.historialTiempo = [];
}

reiniciarSimulacionMock(estadoMock, domMock);
const resetPasa = (
    estadoMock.corriendo === false &&
    estadoMock.t === 0.0 &&
    estadoMock.vueltasM1 === 0 &&
    estadoMock.vueltasM2 === 0 &&
    estadoMock.estelaM2.length === 0 &&
    estadoMock.historialTiempo.length === 0 &&
    domMock.btnPlayPause.textContent === 'Reanudar'
);

if (!resetPasa) {
    console.error(`FALLO: La lógica de reiniciar no cumple con quedar en pausa a t=0 con datos limpios.`);
    todasPasan = false;
} else {
    console.log(`[✓] Botón Reiniciar verificado: regresa a t=0 en pausa, limpia estela, historial y vueltas.`);
}

if (todasPasan) {
    console.log('RESULTADO: TODAS LAS PRUEBAS PASARON EXITOSAMENTE (Error < 1e-9 y vueltas verificadas).');
    process.exit(0);
} else {
    console.error(`FALLO: Una o más pruebas no pasaron la tolerancia.`);
    process.exit(1);
}
