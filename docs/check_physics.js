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
            x = PenduloFisica.pasoRK4(x, dt, params);
        }
    }

    if (errorCasoMax > errorTotalMaximo) {
        errorTotalMaximo = errorCasoMax;
    }

    const pasa = errorCasoMax < TOLERANCIA;
    if (!pasa) todasPasan = false;

    const simbolo = pasa ? '✓' : '✗';
    console.log(`[${simbolo}] Caso ${idx + 1}: ${caso.nombre}`);
    console.log(`    Error máximo vs Python: ${errorCasoMax.toExponential(4)} (en paso ${pasoPeor})`);
});

console.log('------------------------------------------------------');
console.log(`Error global máximo observado: ${errorTotalMaximo.toExponential(4)}`);

if (todasPasan) {
    console.log('RESULTADO: TODAS LAS PRUEBAS PASARON EXITOSAMENTE (Error < 1e-9).');
    process.exit(0);
} else {
    console.error(`FALLO: El error global (${errorTotalMaximo.toExponential(4)}) excede la tolerancia (${TOLERANCIA}).`);
    process.exit(1);
}
