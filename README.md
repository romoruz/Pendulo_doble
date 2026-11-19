# Péndulo Doble Conservativo

Simulación, análisis dinámico y visualización interactiva del péndulo doble sin fricción.  
Proyecto de Sistemas Dinámicos II, Departamento de Matemáticas Aplicadas, ITAM.

## Entregables Principales

1. **Simulador Web Interactivo**: Página estática en HTML/JS lista para publicarse en GitHub Pages (o ejecutarse de manera local sin servidor ni conexión a internet).
2. **Análisis Numérico y Teórico en Python**: Código para la derivación de equilibrios, linealización, órdenes de convergencia, análisis bootstrap y mapas de caos.
3. **Láminas de Exposición**:
   - `figures/lamina_a.png` (y `.svg`): Conservación de energía, orden de convergencia en el espacio de estados y en energía, y comparación de trayectorias lineales vs no lineales.
   - `figures/lamina_b.png` (y `.svg`): Exploratorio de 1000 condiciones iniciales con ruido correlacionado (Cholesky), histogramas de vueltas con intervalos bootstrap BCa y análisis de estabilidad numérica macroscópica ($dt$ vs $dt/2$).
   - `figures/figura_c.png` (y `.svg`): Mapas de exponente de Lyapunov de tiempo finito (FTLE) en malla de 200×200 con contornos de energía de los puntos silla ($E_{c,m2}$ y $E_{c,m1}$), demostrando que la frontera del caos no coincide con los umbrales de energía.
   - `figures/figura_d.png` (y `.svg`): Paisaje de energía potencial $V(\theta_1, \theta_2)$ en el toro plano $[-\pi, \pi]^2$ con isolíneas, los 4 equilibrios y los umbrales de energía de las sillas ($E_{c,m2}, E_{c,m1}, E_{c,\text{ambas}}$).

---

## Instalación y Configuración del Entorno

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

## Ejecución de Pruebas Unitarias

El proyecto cuenta con una suite completa de pruebas unitarias y de consistencia validadas con `pytest`:

```bash
pytest -v
```

---

## Generación de Figuras y Resultados

Un solo comando reproduce todas las figuras en alta resolución (PNG y SVG), los datos brutos en `results/results.json` y el informe en prosa simple en `results/resumen.md`:

```bash
python analysis/make_figures.py
```

---

## Definiciones Metodológicas Clave

### 1. Definición de "Vuelta" (Paso por Arriba) de la Masa 2
Llevando el ángulo $\theta_2(t)$ de manera continua sin envolver (sin restringir artificialmente a $[-\pi, \pi]$):
$$\kappa(t) = \left\lfloor \frac{\theta_2(t) + \pi}{2\pi} \right\rfloor$$
Cada cambio en el valor entero de $\kappa$ entre pasos temporales sucesivos representa un paso por la vertical superior ($\theta_2 = \pi \pmod{2\pi}$):
- **Conteo sin signo (vueltas totales)**: $\sum_{n} |\kappa(t_{n+1}) - \kappa(t_n)|$
- **Conteo neto con signo (giros netos)**: $\sum_{n} (\kappa(t_{n+1}) - \kappa(t_n)) = \kappa(t_{\text{final}}) - \kappa(t_{\text{inicial}})$

### 2. Exponente de Lyapunov de Tiempo Finito (FTLE)
Se calcula mediante el método canónico de dos trayectorias renormalizadas de Benettin:
- Condición base $x_0$ y condición perturbada $x_{0,\text{pert}} = x_0 + [d_0, 0, 0, 0]^T$ con separación inicial $d_0 = 10^{-8}$ rad en $\theta_1$.
- Cada $\tau = 10$ pasos de integración con RK4 3/8, se calcula la distancia euclidiana $d_k = \|x_{\text{pert}} - x\|_2$, se acumula $\ln(d_k / d_0)$ y se renormaliza la trayectoria perturbada en la dirección de la separación a longitud $d_0$.
- El exponente de Lyapunov de tiempo finito acumulado en el horizonte $T$ es:
$$\lambda_T = \frac{1}{T} \sum_{k=1}^{K} \ln\left(\frac{d_k}{d_0}\right)$$

### 3. Orden de Convergencia: Estado vs Energía
- **Error en el Estado**: Medido contra una solución de referencia de alta resolución ($dt_{\text{ref}} = \min(dt) / 16$) en régimen regular. RK4 3/8 presenta una pendiente asintótica en escala log-log de $\approx 4.0$, confirmando su orden 4 global. El Punto Medio Implícito Canónico presenta pendiente $\approx 2.0$, confirmando su orden 2.
- **Error en la Energía**: Para el oscilador armónico lineal ($z = \omega\, dt$), el polinomio de estabilidad de RK4 cumple $|R(iz)|^2 = 1 - z^6/72 + O(z^8)$. Por ende, el error de energía por paso es $O(dt^6)$ y la acumulación en $T/dt$ pasos es $O(dt^5)$. En régimen no lineal regular esto constituye una heurística que verificamos numéricamente con una pendiente empírica de $\approx 4.93$.

---

## Guía de Exposición (Mapeo Conceptual para Sistemas Dinámicos II)

Orden sugerido de presentación de láminas y conceptos clave a defender:

1. **Lámina A (A1, A2, A3) — Fundamentos Numéricos y Linealización**:
   - **A1 (Conservación de energía a corto vs largo plazo)**: Defiende la diferencia entre métodos estándar y simplécticos. En horizontes cortos ($T \le 40$ s), RK4 es más preciso por su alto orden; a largo plazo ($T > 400$ s), la deriva secular de RK4 domina y el Punto Medio Implícito gana al mantener la energía acotada en una variedad simpléctica cercana (teorema de backward error analysis).
   - **A2 (Órdenes de convergencia estado y energía)**: Valida la implementación rigurosa: orden 4 en estado para RK4, orden 2 para Punto Medio, y la superconvergencia $O(dt^5)$ en energía de RK4.
   - **A3 (Linealización y modos normales)**: Explica la aproximación armónica alrededor del equilibrio estable $(0,0,0,0)$. Valida que para amplitudes pequeñas ($0.05$ rad) la solución lineal describe fielmente el sistema desacoplado en modos normales ($\omega_1, \omega_2$), mientras que a amplitudes moderadas ($0.8$ rad) aparecen armónicos superiores y corrimiento de frecuencias. Se discute por qué Hartman-Grobman no aplica en el origen (centro no hiperbólico).

2. **Figura D — Paisaje de Energía Potencial y Puntos de Equilibrio**:
   - Muestra la topología del espacio de fases sobre el toro plano $[-\pi, \pi]^2$ y las curvas de nivel de $V(\theta_1, \theta_2)$.
   - Identifica los cuatro equilibrios: mínimo $(0,0)$, dos sillas $(\pi,0)$ y $(0,\pi)$ con sus energías críticas $E_{c,m2}$ y $E_{c,m1}$, y el máximo $(\pi,\pi)$ con $E_{c,\text{ambas}}$.
   - Introduce la noción de variedades invariantes y separatrices asociadas a los puntos de silla que organizan el flujo global.

3. **Figura C — Transición al Caos y Mapas de FTLE**:
   - Muestra cómo la frontera del caos no coincide con los umbrales de energía de los puntos silla ($E_{c,m2}$ o $E_{c,m1}$).
   - Destaca la presencia de lenguas regulares (en fase $\theta_1 \approx \theta_2$) que sobreviven muy por encima de $E_{c,m2}$, y zonas de caos que emergen por debajo de umbrales en masas desiguales.
   - Relaciona el exponente de Lyapunov con la divergencia exponencial y la pérdida de previsibilidad determinista.

4. **Lámina B — Estadística Colectiva y Sensibilidad Macroscópica**:
   - Compara ensembles de 1000 trayectorias: a baja energía el número de vueltas es idénticamente cero (IC degenerado $[0,0]$); a alta energía el histograma exhibe asimetría y variabilidad caótica, cuantificada con bootstrap BCa.
   - **Estabilidad macroscópica ($dt$ vs $dt/2$)**: Muestra la paradoja del caos: aunque el $87.5\%$ de las trayectorias individuales divergen microscópicamente al refinar el paso temporal, la distribución estadística macroscópica se preserva (prueba Kolmogorov-Smirnov no rechaza con $p = 0.888$), demostrando propiedades ergódicas del atractor en el espacio de fases.

---

## Respaldo sin Internet para la Exposición

1. **Simulador Web Local**: El archivo `docs/index.html` es una página web estática sin dependencias externas, sin CDNs ni llamadas de red. Puede abrirse directamente en cualquier navegador arrastrándolo o mediante la URI `file:///.../docs/index.html`.
2. **Figuras Vectoriales y en Alta Definición**: Las carpetas `figures/` contienen tanto archivos PNG a 200 DPI como gráficos vectoriales SVG listos para proyectar sin pérdida de calidad.
3. **Grabación de Pantalla**: En entornos Linux GNOME es posible iniciar una grabación local con `Ctrl+Shift+Alt+R`, o bien con herramientas como OBS Studio.
