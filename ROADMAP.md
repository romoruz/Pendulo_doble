# Péndulo doble conservativo: roadmap para Claude Code

Lee este archivo completo antes de escribir código. Luego dame un plan breve (máximo 15 líneas) y pregúntame solo lo que de verdad sea ambiguo. Trabaja por fases: al terminar cada una, corre las pruebas, haz commit y detente para que yo revise antes de seguir.

## 1. Contexto

Somos un equipo de Matemáticas Aplicadas del ITAM (Sistemas Dinámicos II). Tenemos que exponer y entregar un informe breve sobre el péndulo doble **sin fricción**. El profesor no pidió código: esto es un plus, y tienes libertad de diseño dentro de las restricciones de este documento.

El informe base del equipo tiene errores conocidos. La sección 7 lista los que ya identificamos. El código debe seguir **este** documento, no el informe; si encuentras otra discrepancia, avísame.

Mi computadora es Linux. Nada en el flujo debe depender de macOS ni de software de pago.

### Entregables

1. **Simulador web interactivo**: una página estática en HTML/JS publicada en GitHub Pages. Se tiene que poder proyectar desde un iPad o desde cualquier computadora sin instalar nada.
2. **Análisis en Python** que genera las figuras de la expo:
   - **Lámina A, "La simulación confirma la teoría"**:
     - La energía se conserva. Se compara RK4 contra un integrador simpléctico y se valida el orden de cada uno.
     - Trayectoria no lineal contra la linealizada cerca del origen (modos normales).
   - **Lámina B, "Exploratorio de simulación"**: 1000 condiciones iniciales con ruido normal pequeño alrededor de un punto de baja energía y otro de alta. Se grafica el histograma de cuántas veces la masa 2 pasa por arriba en un tiempo fijo T, con IC bootstrap.
   - **Figura C, mapa de caos**: el plano (θ1(0), θ2(0)) partiendo del reposo, coloreado por el exponente de Lyapunov de tiempo finito (FTLE), con los contornos de energía de los puntos silla encima. Muestra que la frontera del caos no coincide con los umbrales de energía.

## 2. Modelo físico (fuente única de verdad)

No cambies nada de esta sección sin preguntarme.

- Masas puntuales `m1`, `m2` y varillas rígidas sin masa `l1`, `l2`, no necesariamente iguales. `g = 9.81`. Unidades SI. **Sin fricción.**
- Ángulos **absolutos, medidos desde la vertical hacia abajo**, positivos en sentido antihorario.
  - `x1 = l1 sin θ1`, `y1 = −l1 cos θ1`
  - `x2 = x1 + l2 sin θ2`, `y2 = y1 − l2 cos θ2`
- Estado: `x = (θ1, θ2, ω1, ω2)`. Notación: `Δ = θ1 − θ2`.
- Energía cinética: `T = ½(m1+m2) l1² ω1² + ½ m2 l2² ω2² + m2 l1 l2 ω1 ω2 cos Δ`
- Energía potencial, con referencia `V(0,0) = 0` como en el informe: `V = (m1+m2) g l1 (1 − cos θ1) + m2 g l2 (1 − cos θ2)`. Energía total: `E = T + V`, que se conserva.
- Ecuaciones de movimiento `M(θ) ω̇ = f`:
  - `M11 = (m1+m2) l1²`, `M12 = M21 = m2 l1 l2 cos Δ`, `M22 = m2 l2²`
  - `det M = m2 l1² l2² (m1 + m2 sin²Δ) > 0`
  - `f1 = −m2 l1 l2 sin Δ · ω2² − (m1+m2) g l1 sin θ1`
  - `f2 = +m2 l1 l2 sin Δ · ω1² − m2 g l2 sin θ2`
  - Resuelve el sistema 2×2 con Cramer.
- Forma hamiltoniana, necesaria para el integrador simpléctico. Momentos canónicos `p = M(θ) ω`, `H(θ, p) = ½ pᵀ M(θ)⁻¹ p + V(θ)`:
  - `θ̇ = M(θ)⁻¹ p` (es decir, θ̇ = ω)
  - `ṗ1 = −m2 l1 l2 ω1 ω2 sin Δ − (m1+m2) g l1 sin θ1`
  - `ṗ2 = +m2 l1 l2 ω1 ω2 sin Δ − m2 g l2 sin θ2`
  - Aquí ω se calcula como `M(θ)⁻¹ p`.
- **Equilibrios**, en el orden de estado (θ1, θ2, ω1, ω2):
  - `(0,0,0,0)` es un centro: eigenvalores puramente imaginarios, estable por Lagrange-Dirichlet con E como función de Lyapunov.
  - `(π,0,0,0)` y `(0,π,0,0)` son sillas con un par real y un par imaginario.
  - `(π,π,0,0)` tiene dos pares reales.
  - Verifica los eigenvalores numéricamente y repórtalos.
- **Linealización en el origen**: `ẋ = A x`, con `A = [[0, I], [−M0⁻¹K, 0]]`, donde
  - `M0 = [[(m1+m2) l1², m2 l1 l2], [m2 l1 l2, m2 l2²]]`
  - `K = diag((m1+m2) g l1, m2 g l2)`
  - Las frecuencias de los modos normales son las raíces de los eigenvalores de `M0⁻¹K`.
  - Ojo: el origen **no es hiperbólico**, así que Hartman-Grobman no aplica. La comparación lineal contra no lineal se presenta como aproximación de amplitud pequeña.
- **Umbrales de energía**, que son los valores de V en las sillas:
  - `Ec_m2 = V(0,π) = 2 m2 g l2`: si `E < Ec_m2`, la masa 2 nunca puede pasar por arriba. Como E se conserva, esto vale para todo t.
  - `Ec_m1 = V(π,0) = 2 (m1+m2) g l1`: igual, para la masa 1.
  - `Ec_ambas = V(π,π) = Ec_m1 + Ec_m2`.
  - **Cuál umbral es el primero depende de los parámetros**: el de la masa 2 va primero solo si `m2 l2 < (m1+m2) l1`. No asumas el orden en ningún lado.
  - Los umbrales son condiciones **necesarias** para dar vueltas, **no** criterios de caos.
- **Escalamiento**: con unidad de tiempo `√(l1/g)` y unidad de energía `(m1+m2) g l1`, la dinámica solo depende de `μ = m2/m1` y `λ = l2/l1`. Úsalo para explicar en la página que multiplicar todas las masas (o todas las longitudes) por la misma constante no cambia la forma del movimiento.
- **Definición de "vuelta" (paso por arriba) de la masa 2**: lleva θ2 sin envolver. Sea `k(t) = floor((θ2 + π) / (2π))`. Cada cambio de `k` entre un paso y el siguiente cuenta como un paso por arriba. Guarda el conteo sin signo y el neto con signo, y documenta esta definición en el README y en el resumen de resultados.

## 3. Integradores

- **RK4, variante 3/8 de Kutta** (la que menciona el informe), de paso fijo. Mismo algoritmo en Python y en JS. En Python, vectorízalo sobre N condiciones iniciales con arrays de forma `(N, 4)`.
- **Punto medio implícito en variables canónicas (θ, p)**: es simpléctico y de orden 2. Resuelve la ecuación implícita por iteración de punto fijo con tolerancia 1e-13 y tope de iteraciones; si no converge, levanta un error.
  - Importante: aplicado en (θ, ω) **no** es simpléctico. Debe ser en (θ, p).
  - Solo se usa en Python para la lámina A.
- **FTLE**: método de dos trayectorias renormalizadas (Benettin), con separación inicial d0 = 1e-8 en θ1 y renormalización cada 10 pasos. El resultado es `λ_T = (1/T) Σ log(d_k/d0)`.

## 4. Estructura del repo

```
pendulo-doble/
  README.md
  requirements.txt          # numpy, scipy, matplotlib, pytest. Nada más.
  analysis/
    pendulo/
      __init__.py
      config.py             # parámetros por defecto y presets (fuente única)
      model.py              # M, f, energía, forma hamiltoniana, A, umbrales, equilibrios
      integrate.py          # RK4 3/8 escalar y vectorizado, punto medio implícito
      flips.py              # conteo de pasos por arriba
      chaos.py              # FTLE y mapa de caos
      stats.py              # bootstrap y utilidades
    make_figures.py         # un solo comando genera todas las figuras y resultados
    export_web.py           # genera docs/presets.js y docs/reference.json
    tests/
  figures/                  # PNG y SVG generados (se versionan)
  results/
    results.json
    resumen.md
  docs/                     # raíz de GitHub Pages
    index.html
    style.css
    physics.js
    app.js
    presets.js              # generado desde config.py
    reference.json          # generado desde Python para la prueba de consistencia
    check_physics.js        # prueba en Node: JS vs Python
    .nojekyll
```

## 5. Fases

### Fase 0: preparación

- `git init`, `.gitignore` (`.venv`, `__pycache__`, `.pytest_cache`), `requirements.txt`, esqueleto de carpetas.
- Instrucciones en el README: `python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt`.

### Fase 1: núcleo en Python y pruebas

Pruebas con pytest; todas deben pasar:

1. M es simétrica y `det M > 0` en una malla de ángulos, con varios juegos de parámetros (incluye desiguales, por ejemplo m1=3, m2=0.5, l1=0.7, l2=1.4).
2. RK4 3/8 con dt = 1e-3 y T = 20 s: deriva relativa de energía < 1e-9, para tres juegos de parámetros.
3. Orden: la pendiente log-log del error máximo de energía contra dt cae en [3.7, 4.3] para RK4 y en [1.8, 2.2] para punto medio.
4. Punto medio implícito: el error de energía se mantiene acotado y sin tendencia en 2000 s. Prueba de que no hay deriva secular: la pendiente de una regresión del error sobre t es compatible con 0.
5. La forma hamiltoniana y la lagrangiana dan la misma trayectoria (diferencia < 1e-10 con RK4 a dt = 1e-3, T = 5 s).
6. Linealización: con amplitud 1e-3 rad, la diferencia relativa entre la trayectoria no lineal y `expm(A t) x0` es < 1e-3 en t ∈ [0, 10].
7. Equilibrios: f se anula en los cuatro, y los eigenvalores del Jacobiano numérico tienen la estructura de la sección 2.
8. Si `E0 < Ec_m2`, el conteo de pasos por arriba es 0 hasta T = 200 s.
9. La versión vectorizada y la escalar dan el mismo resultado (diferencia < 1e-12).
10. FTLE: cerca del origen (amplitud 0.05 rad) λ_T < 0.2 con T = 40 s; con θ1 = θ2 = 2.8 rad desde el reposo, λ_T > 0.5. Parámetros iguales a 1.

### Fase 2: figuras y resultados en Python

En `config.py` define el juego de parámetros de las figuras. Por defecto: m1 = m2 = 1 kg, l1 = l2 = 1 m, como en el informe. Semilla fija: `np.random.default_rng(20261119)`.

**Lámina A**

- **A1 Energía.** Error relativo de energía contra t en 200 s, con RK4 3/8 y punto medio implícito al mismo dt (0.005, como en el informe). Escala log en y. Se debe ver que RK4 deriva y que el punto medio oscila acotado. Usa una condición de energía alta.
- **A2 Orden.** Error máximo de energía contra dt en log-log, para ambos integradores, con dt ∈ {0.04, 0.02, 0.01, 0.005, 0.0025} y T = 10 s. Agrega rectas de referencia de pendiente 4 y 2 y anota las pendientes ajustadas. Si algún dt toca el piso de redondeo, quítalo y dilo en el resumen.
- **A3 Lineal contra no lineal.** θ1(t) y θ2(t) del sistema no lineal contra los de `expm(A t) x0`, en dos amplitudes: una chica (≈0.05 rad), donde casi coinciden, y una mediana (≈0.8 rad), donde se separan. Anota las frecuencias de los modos normales.

**Lámina B**

- Dos puntos base, partiendo del reposo:
  - **Baja energía**: verifica en el código que `E0 < Ec_m2`, con lo que la teoría garantiza 0 vueltas.
  - **Alta energía**: por ejemplo θ1 = θ2 ≈ 2.6 rad. Ajústalo para que, con T = 30 s, la media de pasos por arriba sea ≥ 3 y el histograma quede claramente abierto.
- N = 1000. El ruido va solo en los ángulos: `θ0 = base + L z`, con `z ~ N(0, I₂)` y `L = chol(Σ)`, donde Σ es 2×2 con σ = 0.01 rad y correlación ρ (por defecto 0, configurable; la idea es mostrar Cholesky).
- **B1 Histogramas.** Los de baja y alta energía lado a lado, con el mismo eje x.
- **Bootstrap** con `scipy.stats.bootstrap`, B = 10000, método BCa, para la media y la desviación estándar del conteo. Si la muestra es degenerada (todo ceros), no llames a BCa: reporta IC = [0, 0] y anótalo como degenerado.
- **Estabilidad numérica.** En alta energía, repite con dt/2. Los conteos individuales van a cambiar, porque el caos amplifica el error numérico, pero la distribución no debe cambiar. Compáralas con KS de dos muestras y con los IC. Documéntalo en el resumen: es punto de discusión de la expo.

**Figura C, mapa de caos**

- Rejilla de 200×200 en (θ1(0), θ2(0)) ∈ [−π, π]², partiendo del reposo. Color: FTLE con T = 40 s, dt = 0.01, escala de color fija.
- Encima, los contornos `V = Ec_m2` y `V = Ec_m1`, cada uno con su etiqueta.
- Dos paneles: parámetros iguales, y m1 = 1, m2 = 0.25, l1 = 1, l2 = 1.6.
- En `results.json`, la fracción de condiciones con λ_T > 0.5 por banda de energía (en unidades de `Ec_m2`).
- Ya verificamos con un cálculo rápido lo que debería salir; si no te sale, investiga y repórtalo:
  - Con parámetros iguales, debajo de `Ec_m2` casi todo es regular, pero desde ~0.7·Ec_m2 aparecen zonas con FTLE intermedio (0.3–0.5).
  - Entre `Ec_m2` y `Ec_m1`, alrededor de 20% sigue siendo regular, incluida una "lengua" regular cerca del modo en fase.
  - En el segundo juego de parámetros, la región regular rebasa claramente el contorno de `Ec_m2`.

**Opcional: secciones de Poincaré**

Reproduce los Casos 1 y 2 del informe con Σ = {θ1 = 0 mod 2π, ω1 < 0}, en el plano (θ2, ω2):
- Caso 1: x0 = [0.20, 0.25, 0, 0].
- Caso 2: x0 = [π−0.05, π−0.05, 0, 0].
- Detecta los cruces con interpolación (el truco de Hénon o una interpolación cúbica entre pasos).

**Estilo de las figuras**

- Fondo blanco y sobrio. θ1 en azul (`#1f4fd6`) y θ2 en rojo (`#d62728`), consistente en todas las figuras.
- Ejes en español con unidades. Fuente sans. Nada decorativo y sin títulos rimbombantes: que no parezca hecho por una IA.
- Tamaño 16:9 para láminas (12 × 6.75 in, 200 dpi). Exporta PNG y SVG a `figures/`.
- Todos los números van a `results/results.json`.

**`results/resumen.md`**

Prosa simple con los números para el informe: eigenvalores de los cuatro equilibrios, frecuencias de los modos normales, pendientes de orden, magnitud de la deriva de cada integrador, estadísticas de los histogramas con sus IC, estabilidad con dt/2, fracciones caóticas por banda de energía, y la metodología (semilla, N, σ, ρ, dt, T, definición de vuelta, definición de FTLE).

### Fase 3: simulador web

**Restricciones técnicas**

- Página estática, sin build, sin dependencias, sin CDN y sin ES modules: tiene que funcionar también abriéndola como `file://`.
- `physics.js` es un script clásico que expone `window.Pendulo` y además hace `module.exports` cuando corre en Node.
- Para el Web Worker del mapa de caos, crea el worker desde un Blob armado con el código fuente (`Function.prototype.toString`). Un `new Worker('archivo.js')` falla en `file://`.
- Objetivo principal: Safari de iPad en horizontal. También laptop y proyector a 1920×1080.
- Usa Pointer Events para el touch. Nada que dependa de hover. Áreas táctiles de al menos 44 px.
- Física con RK4 3/8, dt fijo de 1e-3 s y acumulador de tiempo real, con tope de subpasos por cuadro.

**Dibujo**

- Pivote arriba al centro. La escala en px/m es fija, calculada con la longitud máxima posible (l1 + l2 con los sliders al tope), para que alargar una varilla se vea más larga de verdad. Recalcula al cambiar el tamaño de la ventana.
- Varillas como líneas negras. Masas como círculos con radio proporcional a la raíz cúbica de m (con mínimo y máximo en px) y su etiqueta ("1 kg", "14 kg").
- Líneas verticales punteadas por el pivote y por la masa 1, con arcos de ángulo: θ1 en azul y θ2 en rojo.
- Estela de la masa 2 que se desvanece, activable.

**Controles**, cada slider con su valor numérico visible:

- `l1`, `l2` ∈ [0.2, 2.0] m, paso 0.05
- `m1`, `m2` ∈ [0.1, 20] kg, paso 0.1
- `θ1(0)`, `θ2(0)` ∈ [−180°, 180°]
- `ω1(0)`, `ω2(0)` ∈ [−10, 10] rad/s, por defecto 0
- Iniciar/Pausa, Reiniciar, velocidad (0.25×, 0.5×, 1×, 2×)
- Con la simulación en pausa, se pueden arrastrar las masas con el dedo o el mouse para fijar la condición inicial.
- Cambiar cualquier parámetro reinicia la simulación desde la condición inicial vigente.
- Junto a los sliders, mostrar `m2/m1` y `l2/l1`, con una nota corta: "la forma del movimiento solo depende de estos dos cocientes".

**Panel de energía y umbrales**, en vivo:

- E, error relativo de energía, t, pasos por arriba de cada masa.
- Una barra vertical de energía con marcas en `Ec_m2`, `Ec_m1` y `Ec_ambas`, en el orden real según los parámetros, y la E actual.
- Texto según la región:
  - "Ninguna masa puede dar la vuelta (garantizado por energía)"
  - "La masa 2 ya puede dar la vuelta"
  - etc.
  - Aclara que poder dar la vuelta no implica caos.

**Mapa de caos**

- Panel cuadrado con el plano (θ1(0), θ2(0)) desde el reposo, coloreado por FTLE para los parámetros actuales.
- Se calcula en el Web Worker con resolución progresiva (32×32 → 64×64 → 96×96), T = 20 s y dt = 0.01. Se recalcula al cambiar masas o longitudes.
- Encima, los contornos `V = Ec_m2` y `V = Ec_m1`.
- La condición inicial actual aparece como un punto. Tocar el mapa fija esa condición inicial (con ω = 0).
- Rótulo claro: "indicador aproximado: exponente de Lyapunov de tiempo finito".

**Modos**

1. **Normal.**
2. **Gemelo**: un segundo péndulo semitransparente que parte con θ1 + δ (δ configurable, por defecto 1e-3 rad). Mini gráfica en vivo de `log10 |separación en el espacio de estados|` contra t; en régimen caótico la pendiente inicial se ve recta.
3. **Lineal contra no lineal**: un péndulo fantasma gris integra `ẋ = A x`.
4. **Nube**: N péndulos (por defecto 200, slider de 50 a 500) con ángulos `base + L z`.
   - Box-Muller para z, slider de σ en escala log de 1e-4 a 5e-2 rad, y ρ opcional.
   - Se dibuja la masa 2 de cada uno como punto, con varillas muy tenues.
   - Abajo, un histograma en vivo de pasos por arriba por péndulo.

**Presets**

- Botones "Cerca del origen", "Baja energía", "Alta energía", "Caso 1 del informe" y "Caso 2 del informe", con los mismos puntos de Python. Se leen de `presets.js`, que genera `export_web.py` desde `config.py`.
- Botón "Modo presentación" que oculta los controles y agranda el canvas. En iPad la Fullscreen API es limitada; no dependas de ella.

**Estilo**

- Que parezca herramienta de laboratorio, no landing page.
- Fondo blanco, texto negro, fuente del sistema, sin gradientes, sin sombras, sin emojis.
- Interfaz en español.

**Prueba de consistencia**

- `export_web.py` escribe en `reference.json` trayectorias de Python: 3 juegos de parámetros × 2000 pasos.
- `node docs/check_physics.js` las compara con el RK4 de JS. Diferencia máxima absoluta < 1e-9.
- Agrega una prueba de pytest que corra esto si `node` está disponible; si no, que la salte.

### Fase 4: publicación

- README sencillo, sin decoración:
  - Qué es el proyecto y el link al simulador.
  - Cómo regenerar figuras (`python analysis/make_figures.py`) y cómo correr las pruebas (`pytest`).
  - La definición de vuelta y la de FTLE.
- Instrucciones para GitHub Pages: Settings → Pages → Deploy from a branch → `main` → `/docs`.
- Si `gh` está instalado, propón los comandos para crear el repo y hacer push, pero **no hagas push sin preguntarme**.
- Sección "Respaldo sin internet" en el README:
  - Grabar el simulador en video (en GNOME: Ctrl+Shift+Alt+R, o con OBS).
  - Llevar las figuras en PDF.
  - El simulador también abre como archivo local en una laptop.

## 6. Reglas generales

- No cambies el modelo físico ni las definiciones sin preguntar.
- Si un criterio de aceptación no se cumple, explícame por qué. No relajes el umbral en silencio.
- Commits pequeños, al menos uno por fase.
- Código legible, con comentarios en español donde aporten. Los nombres de variables pueden ir en inglés.
- No agregues dependencias fuera de numpy, scipy, matplotlib y pytest.
- Antes de dar por terminada una fase, abre o renderiza lo que produjiste (figuras, página) y revísalo visualmente.

## 7. Errores conocidos del informe base (no los copies al código)

1. **Sección 2.2: los equilibrios tienen los ángulos en las posiciones de las velocidades.** Con x = (θ1, θ2, ω1, ω2), los correctos son (0,0,0,0), (π,0,0,0), (0,π,0,0) y (π,π,0,0).
2. **"tr(J) = 0 implica preservación de volumen (Liouville)"** es incorrecto por dos razones:
   - Que la traza se anule en los equilibrios no dice nada del flujo en general.
   - En (θ, ω) la divergencia del campo no es cero. Liouville vale en las variables canónicas (θ, p).
3. **"Con E < Ec1 todo es KAM y al rebasarlo se desata el caos."** Los umbrales son condiciones necesarias para dar vueltas, no la frontera del caos. Puede haber caos débil debajo del umbral y zonas regulares arriba (ver figura C).
4. **"Ec1 es el primer umbral."** Solo es cierto si `m2 l2 < (m1+m2) l1`. Con otros parámetros, el umbral de la masa 1 es menor.
5. **"Con E como función de Lyapunov, las trayectorias quedan en curvas de nivel cerradas."** En 4D los conjuntos de nivel cerca del mínimo son 3-esferas. La conclusión de estabilidad sí es correcta.
6. **La ecuación (13) da el determinante de la M simétrica**, pero la M de la ecuación (10) está escrita dividida entre l1 y l2, y su determinante es `m2 l1 l2 (m1 + m2 sin²Δ)`. Las k_ij de (24)–(25) sí son consistentes con la ecuación (10).
<ADDITIONAL_METADATA>
The current local time is: 2026-11-19T08:39:03-06:00.
</ADDITIONAL_METADATA>
<USER_SETTINGS_CHANGE>
The user changed setting `Model Selection` from None to Gemini 3.8 Flash (High). No need to comment on this change if the user doesn't ask about it. If reporting what model you are, please use a human readable name instead of the exact string.
</USER_SETTINGS_CHANGE>