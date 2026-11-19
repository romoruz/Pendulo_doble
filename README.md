# Péndulo Doble Conservativo

Simulación, análisis dinámico y visualización interactiva del péndulo doble sin fricción.  
Proyecto de Sistemas Dinámicos II, Departamento de Matemáticas Aplicadas, ITAM.

## Entregables Principales

1. **Simulador Web Interactivo**: Página estática en HTML/JS lista para publicarse en GitHub Pages (o ejecutarse de manera local sin servidor ni conexión a internet).
2. **Análisis Numérico y Teórico en Python**: Código para la derivación de equilibrios, linealización, órdenes de convergencia, análisis bootstrap y mapas de caos.
3. **Láminas de Exposición**:
   - `figures/lamina_a.png` (y `.svg`): Conservación de energía, orden de convergencia en el espacio de estados y en energía, y comparación de trayectorias lineales vs no lineales. Exportada también en paneles individuales listos para proyección: `figures/a1.png` / `.svg`, `figures/a2.png` / `.svg`, y `figures/a3.png` / `.svg`.
   - `figures/lamina_b.png` (y `.svg`): Exploratorio de 1000 condiciones iniciales con ruido correlacionado (Cholesky), histogramas de vueltas con intervalos bootstrap BCa y análisis de estabilidad numérica macroscópica ($dt$ vs $dt/2$).
   - `figures/figura_c.png` (y `.svg`): Mapas de exponente de Lyapunov de tiempo finito (FTLE) en malla de 200×200 con contornos de energía de los puntos silla ($E_{c,m2}$ en cian y $E_{c,m1}$ en negro punteado), demostrando que la frontera del caos no coincide con los umbrales de energía.
   - `figures/figura_d.png` (y `.svg`): Paisaje de energía potencial $V(\theta_1, \theta_2)$ en el toro plano $[-\pi, \pi]^2$ con isolíneas, clasificación rigurosa de los 4 equilibrios por eigenvalores y los umbrales de energía de las sillas ($E_{c,m2}, E_{c,m1}, E_{c,\text{ambas}}$).

---

## Instalación y Configuración del Entorno

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

## Ejecución de Pruebas Unitarias

El proyecto cuenta con una suite completa de pruebas unitarias y de consistencia validadas con `pytest` y Node.js:

```bash
# Suite completa de pruebas (Python y consistencia web)
pytest -v

# Validación directa de física en JavaScript contra Python (error < 1e-9)
node docs/check_physics.js
```

---

## Generación de Figuras y Resultados

Un solo comando reproduce todas las figuras en alta resolución (PNG y SVG), los datos brutos en `results/results.json` y el informe en prosa simple en `results/resumen.md`:

```bash
python analysis/make_figures.py
```

Para actualizar los presets y archivos de referencia del simulador web:

```bash
python analysis/export_web.py
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

El orden de la exposición sigue estrictamente el índice del informe del equipo:

1. **Sección 4 del informe — Formulación hamiltoniana y espacio fase (Figura D)**:
   - **Figura D (Paisaje de potencial y clasificación de equilibrios)**: Muestra la topología del toro plano $[-\pi, \pi]^2$ y las curvas de nivel de $V(\theta_1, \theta_2)$. Defiende la clasificación rigurosa según los eigenvalores del Jacobiano: $(0,0)$ centro con eigenvalores imaginarios puros; $(0,\pm\pi)$ y $(\pm\pi,0)$ sillas-centro no hiperbólicas (un par imaginario puro y un par real simétrico); y $(\pm\pi,\pm\pi)$ silla doble estrictamente hiperbólica. Identifica los niveles de energía crítica de las separatrices: $E_{c,m2}$, $E_{c,m1}$ y $E_{c,\text{ambas}}$.

2. **Sección 5 del informe — Linealización y modos normales (A3 y Figura D)**:
   - **A3 (Comparación lineal vs no lineal)**: Explica la dinámica local alrededor del reposo $(0,0,0,0)$ mediante el sistema desacoplado en modos normales ($\omega_1, \omega_2$). Para pequeñas amplitudes ($0.05$ rad) la solución lineal $\dot{x}=Ax$ describe fielmente el sistema; para amplitudes moderadas ($0.8$ rad) surgen armónicos superiores y desfase por no linealidad.
   - **Conexión teórica con Figura D**: Explica por qué el Teorema de Hartman-Grobman no aplica en el origen (centro no hiperbólico), requiriendo el uso de la energía mecánica total como función de Lyapunov estricta (Teorema de Lagrange-Dirichlet) para demostrar estabilidad.

3. **Sección 6 del informe — Métodos numéricos, convergencia y estadística de vueltas (A1, A2 y Lámina B)**:
   - **A1 (Conservación de energía a largo plazo)**: Compara RK4 3/8 contra el integrador simpléctico de punto medio implícito. En horizontes moderados ($T \le 40$ s con $dt = 0.005$ s), RK4 es superior con error relativo $< 7 \times 10^{-6}$. A largo plazo ($T = 1500$ s), la deriva creciente de RK4 cruza y supera en $t \approx 884$ s la envolvente del error del punto medio, demostrando la superioridad simpléctica para integraciones prolongadas.
   - **A2 (Órdenes de convergencia estado y energía)**: Demuestra el orden 4 global en estado para RK4 y orden 2 para punto medio. En energía, valida la superconvergencia $O(dt^5)$ de RK4 en oscilaciones armónicas regulares derivada analíticamente de $|R(iz)|^2 = 1 - z^6/72 + O(z^8)$ (pendiente observada 4.93).
   - **Lámina B (Exploratorio estocástico y estabilidad macroscópica)**: Con 1000 condiciones iniciales perturbadas vía Cholesky, B1 demuestra 0 vueltas analíticamente garantizadas bajo $E < E_{c,m2}$; B2 muestra el régimen caótico con intervalos bootstrap BCa; y B3 demuestra la estabilidad macroscópica del atractor: aunque el $87.5\%$ de las trayectorias individuales divergen al reducir el paso a $dt/2$, la prueba de Kolmogorov-Smirnov no rechaza la igualdad distributiva ($p = 0.89$), con la observación de que KS es conservador para datos discretos.

4. **Sección 7 del informe — Transición al caos y mapas de FTLE (Figura C)**:
   - **Figura C (Desacople de umbrales energéticos y frontera del caos)**: Analiza los mapas de FTLE en malla 200×200 con contornos de energía de las sillas ($E_{c,m2}$ en cian y $E_{c,m1}$ en negro punteado) para masas iguales y desiguales ($m_2=0.25$).
   - Defiende que los umbrales de energía potencial son necesarios pero no suficientes para el caos: lenguas regulares en modo en fase ($\theta_1 \approx \theta_2$) sobreviven muy por encima de $E_{c,m2}$, y el caos puede emerger antes de $E_{c,m2}$ bajo asimetría de masas.

---

## Respaldo sin Internet para la Exposición

1. **Simulador Web Local**: El archivo `docs/index.html` es una página web estática sin dependencias externas, sin CDNs ni llamadas de red. Puede abrirse directamente en cualquier navegador arrastrándolo o mediante la URI `file:///.../docs/index.html`.
2. **Figuras Vectoriales y en Alta Definición**: Las carpetas `figures/` contienen tanto archivos PNG a 200 DPI como gráficos vectoriales SVG listos para proyectar sin pérdida de calidad.
3. **Grabación de Pantalla**: En entornos Linux GNOME es posible iniciar una grabación local con `Ctrl+Shift+Alt+R`, o bien con herramientas como OBS Studio.
