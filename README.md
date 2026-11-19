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
   - `figures/poincare.png` (y `.svg`): Secciones de Poincaré en $\Sigma = \{\theta_1 = 0 \pmod{2\pi}, \omega_1 < 0\}$ (curva invariante cuasiperiódica vs mar caótico).

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
- **Error en la Energía**: Para el Punto Medio Implícito, el error de energía conserva su tasa de segundo orden ($O(dt^2)$). Sin embargo, para RK4 3/8 en oscilaciones regulares, las cancelaciones en el producto interno de gradientes de un sistema hamiltoniano provocan que el error en la energía mecánica converja a una tasa superconvergente de $O(dt^5)$ (pendiente ajustada $\approx 4.93$).

---

## Respaldo sin Internet para la Exposición

1. **Simulador Web Local**: El archivo `docs/index.html` es una página web estática sin dependencias externas, sin CDNs ni llamadas de red. Puede abrirse directamente en cualquier navegador arrastrándolo o mediante la URI `file:///.../docs/index.html`.
2. **Figuras Vectoriales y en Alta Definición**: Las carpetas `figures/` contienen tanto archivos PNG a 200 DPI como gráficos vectoriales SVG listos para proyectar sin pérdida de calidad.
3. **Grabación de Pantalla**: En entornos Linux GNOME es posible iniciar una grabación local con `Ctrl+Shift+Alt+R`, o bien con herramientas como OBS Studio.
