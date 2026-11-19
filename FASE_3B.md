# Fase 3b: rediseño visual del simulador

Este documento complementa a ROADMAP.md, que sigue siendo la referencia del modelo físico; aquí no cambia nada de la física. Trabaja en `docs/`. Al terminar, corre las pruebas (incluida `node docs/check_physics.js`), haz commit y detente para que yo revise. No hagas push.

## 1. Diagnóstico de la versión actual

Revisé una captura de pantalla de la versión actual. Hay que corregir esto:

1. **No hay gráficas en vivo.** La mitad del lienzo está vacía y todo lo cuantitativo son números en una tabla. El simulador tiene que *mostrar* la dinámica.
2. **La parte de simulación y estadística no se ve.** Solo existe el modo Nube, con un histograma chico. Falta una pestaña de Monte Carlo con intervalos de confianza (ver sección 4).
3. **Arcos de ángulo rotos.** El arco de θ2 se dibuja como un círculo completo alrededor de la masa 1, y el de θ1 es diminuto. Deben ser cuñas sombreadas de unos 45 px de radio, desde la vertical punteada hasta la varilla.
4. **Etiquetas de masa ilegibles.** "3.1 kg" se dibuja encima del disco. La etiqueta va al lado del disco, con contraste.
5. **Barra de energía.**
   - Las etiquetas están cortadas ("as: 144.3 J", "1: 124.7 J").
   - La barra es muy angosta.
   - Sigue usando rojo y azul para T y V, y esos colores ya son de θ1 y θ2.
6. **Mapa de caos.**
   - No tiene ejes ni barra de color.
   - Los contornos que promete el texto no se ven.
   - El punto de la condición inicial apenas se distingue.

## 2. Qué reciclar del simulador original del equipo (repo `romoruz/Pendulo_doble`)

El archivo `simulador_pendulo_doble.py` tiene ecuaciones correctas: verifiqué `derivadas()` contra nuestro modelo y la diferencia es de 1e-15. Hay que portar a la web estas ideas:

- **Retratos de fase en vivo (θ1, ω1) y (θ2, ω2)**, con el ángulo envuelto a (−π, π] y la línea cortada en cada salto. Es el truco de `procesar_angulos()` con NaN; en canvas equivale a `moveTo` en vez de `lineTo` cuando |Δθ| > π. Ojo: la versión animada original no corta las líneas y deja rayas horizontales. No repitas ese error.
- **Selector de integrador**: RK4 clásico o RK4 3/8. Muestra en vivo el error de energía de cada uno.
- **Slider de g** con botones Tierra (9.81), Luna (1.62) y Júpiter (24.79). Además de divertido, es una demostración: g solo cambia la escala de tiempo, así que el mapa de caos no debe cambiar si se calcula en tiempo adimensional (ver sección 5).
- **Presets** "Pequeñas oscilaciones" (0.20, 0.25), "Cuasi-periódico" (0.80, 0.60) y "Caótico" (π−0.05, π−0.05), además de los que ya hay.

**No reciclar:**
- Secciones de Poincaré: no van en el proyecto.
- La etiqueta de "régimen" por umbral de energía: es falsa, los umbrales no son frontera de caos.
- El comentario de que la variante 3/8 "reduce el sesgo direccional": no tiene sustento.
- Las "sillas" marcadas en los planos proyectados (θi, ωi). Si se marcan, hay que aclarar que son proyecciones de puntos de R⁴.

## 3. Estética: modo oscuro y más vida

- **Tema oscuro por defecto** (fondo `#0b0f17`, paneles `#121826`, texto `#e6e9ef`) con un botón para cambiar a claro. El modo claro se queda: en un salón iluminado, el proyector suele lavar los temas oscuros.
  - Define todos los colores con variables CSS, y los colores del canvas desde un objeto de tema.
- **Paleta:** θ1 y masa 1 en cian (`#4cc9f0`), θ2 y masa 2 en magenta (`#f72585`), energía en ámbar (`#ffb703`), cinética en verde menta (`#2ec4b6`), potencial en violeta (`#9b5de5`). Úsala consistentemente en todo el simulador.
  - En modo claro, usa versiones más oscuras de los mismos tonos.
- **Estela de la masa 2:**
  - Brillo con `globalCompositeOperation = 'lighter'` y desvanecimiento por edad.
  - Color según la rapidez: frío cuando va lenta, caliente cuando va rápida.
  - Tope de unos 1500 puntos.
- **Masas:** disco con gradiente radial y halo suave, con radio proporcional a la raíz cúbica de m.
- **Estado actual:** números grandes en tipografía monoespaciada (t, E, error de energía, vueltas). Cuando cambia el conteo de vueltas, el número pulsa brevemente.
- **Animaciones:** transiciones de 150 a 250 ms en pestañas y paneles, sin exagerar. Se respeta `prefers-reduced-motion`.
- **Canvas:** consciente del devicePixelRatio, para que se vea nítido en el iPad.
- **Prohibido:** emojis y gradientes de fondo tipo landing page. El brillo va en los datos, no en la decoración.

## 4. Distribución y gráficas en vivo

En horizontal, el péndulo ocupa ~55% del ancho a la izquierda y a la derecha va una rejilla de 2×2 con gráficas en vivo. Los controles van en un cajón plegable (tecla `C` o botón), cerrado por defecto en el modo presentación.

Las gráficas son canvas propios, sin librerías:

1. **θ1(t) y θ2(t):** desplazamiento continuo con ventana de 20 s.
2. **Retratos de fase:** (θ1, ω1) y (θ2, ω2) superpuestos con su color, o con un selector para alternarlos. Envueltos y con cortes, como se describe en la sección 2.
3. **Energía:** T y V como áreas apiladas, E como línea encima, y un renglón con |ΔE|/E₀ en escala log.
4. **Panel según el modo:**
   - **Normal:** el mapa de caos (sección 5).
   - **Gemelo:** log10 de la separación contra t, con un ajuste lineal en vivo sobre la parte creciente. Muestra la pendiente como "λ estimado ≈ … 1/s" y avisa cuando la separación se satura.
   - **Lineal:** θ2 no lineal contra lineal, superpuestos, con su diferencia.
   - **Nube:** histograma en vivo de vueltas por péndulo.

**Pestaña nueva "Monte Carlo".** Es la parte estadística, la sección 6 de la expo:

- **Entradas:** punto base (el actual o un preset), N (100–2000), σ, ρ y T.
- **Botón "Correr".** El cálculo va en un Web Worker con barra de progreso.
- **Resultados:**
  - Histograma de vueltas de la masa 2.
  - Media y desviación estándar con IC 95% por bootstrap percentil (B = 2000).
  - P(al menos una vuelta) con intervalo de Wilson.
- **Comparación A/B:** botón "Comparar con baja energía" que corre el preset de baja energía y pone los dos histogramas lado a lado con el mismo eje. Es la versión interactiva de la lámina B.
- **Animación de la nube:** al terminar, se reproduce la nube de esos N péndulos (solo la masa 2, como puntos brillantes coloreados por su número final de vueltas).
- **Reproducibilidad:** semilla fija y visible, editable, con un generador simple sembrable (mulberry32 o similar). Documenta el generador en el README.
- **Prueba:** agrega una en Node que verifique que el bootstrap y Wilson de JS coinciden con Python en un arreglo fijo.

## 5. Mapa de caos, versión legible

- **Ejes con marcas:** θ1(0) y θ2(0) de −π a π.
- **Barra de color** con unidades.
- **Contornos** `V = Ec_m2` y `V = Ec_m1`, gruesos y visibles en ambos temas.
- **Cruz** en la condición inicial actual.
- **Lectura al pasar el cursor o tocar:** muestra θ1, θ2 y el FTLE del pixel. Al soltar, fija la condición inicial.
- **Resolución progresiva:** hasta 128×128 en laptop y 96×96 en el iPad, detectado por ancho y `hardwareConcurrency`.
- **Mapa de color tipo magma o inferno.** Funciona bien en tema oscuro.
- **Tiempo adimensional.** Calcula el FTLE con T en unidades de √(l1/g) y reporta λ adimensional. Así, al cambiar g, el mapa queda idéntico, y se demuestra en vivo que solo importan m2/m1 y l2/l1. Agrega una línea de texto que lo diga.

## 6. Criterios de aceptación

- Sigue funcionando en `file://`, sin dependencias ni CDN.
- `node docs/check_physics.js` y `pytest` pasan.
- A 60 fps en una laptop normal en los modos Normal y Gemelo, y al menos a 30 fps en Nube con N = 200.
- Nada se encima ni se corta a 1920×1080, a 1366×768 ni a 1180×820 (iPad horizontal). Toma capturas de pantalla en esos tres tamaños con el modo de dispositivo del navegador y en ambos temas, y guárdalas en `docs/capturas/` para que yo las revise.
- Los colores de θ1, θ2, T, V y E son los mismos en el péndulo, en todas las gráficas y en el texto.
<ADDITIONAL_METADATA>
The current local time is: 2026-11-19T15:40:40-06:00.
</ADDITIONAL_METADATA>