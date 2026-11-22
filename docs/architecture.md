# Arquitectura

## Capas y sus responsabilidades

```
┌─────────────────────────────────────────────────────────┐
│ simulation/Problem                                       │
│   amarra mesh + material + element_formulation + BCs     │
│   corre la historia de carga, produce (strain, stress)   │
└───────────────────────┬───────────────────────────────────┘
                         │ usa
        ┌────────────────┼────────────────┐
        ▼                ▼                ▼
┌───────────────┐ ┌──────────────┐ ┌─────────────────────┐
│ assembly/      │ │ solvers/      │ │ boundary_conditions/ │
│ Assembler      │ │ NewtonRaphson │ │ DirichletBC          │
│ (malla global  │ │ ArcLength     │ │ NeumannBC            │
│  -> residuo y  │ │ LinearSolver  │ │ (datos puros)        │
│  tangente       │ │ (no saben     │ └─────────────────────┘
│  globales)      │ │  nada de FEM) │
└───────┬────────┘ └──────────────┘
        │ usa
        ▼
┌───────────────┐        ┌─────────────────────┐
│ elements/      │──usa──▶│ materials/           │
│ FiniteElement  │        │ MaterialModel         │
│ (Quad4, ...)   │        │ (OgdenHyperfoam,      │
│                │        │  NeoHookean, ...)     │
└───────┬────────┘        └─────────────────────┘
        │ usa
        ▼
┌───────────────┐        ┌─────────────────────┐
│ mesh/          │        │ kinematics/          │
│ Mesh, Node,    │        │ F, C, E, invariantes  │
│ Connectivity   │        │ (funciones puras)     │
└───────────────┘        └─────────────────────┘
```

## Reglas de dependencia (quién puede incluir a quién)

- `mesh/` y `kinematics/` no dependen de nada más arriba: son geometría y
  álgebra puras.
- `materials/` depende solo de `core/` (tipos). No sabe que existen
  mallas, elementos ni solvers.
- `elements/` depende de `materials/` (a través de la interfaz
  `MaterialModel`, nunca de una clase concreta como `OgdenHyperfoam`
  directamente) y de `kinematics/`.
- `assembly/` depende de `mesh/` y `elements/`. Es el único lugar que
  traduce "nodo local de un elemento" a "índice global en el vector de
  grados de libertad".
- `solvers/` **no depende de nada de FEM**. `NewtonRaphsonSolver` y
  `ArcLengthSolver` operan sobre `std::function` abstractas
  (`residual(x)`, `tangent(x)`). Esto es deliberado: permite
  implementarlos y probarlos ahora mismo, sin esperar a que el resto del
  proyecto esté listo, y permite reusarlos en cualquier problema no
  lineal futuro que no tenga nada que ver con espumas.
- `simulation/Problem` es la única clase que conoce todas las piezas a la
  vez.

## Por qué un `MaterialModel` abstracto

El objetivo explícito del proyecto es poder comparar distintos modelos
constitutivos (Ogden hyperfoam fenomenológico, un modelo micromecánico
basado en una celda unitaria, o más adelante un sustituto entrenado con
datos) bajo la misma malla, el mismo elemento y el mismo solver. Por eso
`elements/` y `assembly/` solo conocen `MaterialModel` (la interfaz en
`materials/material_base.hpp`), nunca una clase concreta.

## `tangent_moduli()` es el tangente LAGRANGIANO (dS/dE), no uno espacial

Decisión tomada durante la implementación del Paso 1, tras un error real
que vale la pena documentar en vez de esconder:

La primera versión de `tangent_moduli()` devolvía un tangente *espacial*
(relacionando una tasa de esfuerzo de Cauchy con la tasa de deformación
D), con la forma cerrada clásica c_ijkl = (λ/J)δδ + G(δδ+δδ). Al
validarla con un tangente numérico (diferencias finitas directas sobre
`cauchy_stress()`), **no coincidían**, y la discrepancia no era ruido
numérico sino un error conceptual real: una diferencia finita simple
sobre `cauchy_stress()` mide la tasa *simple* del esfuerzo de Cauchy,
que difiere de la tasa *objetiva* (de Truesdell) -- la que esa fórmula
cerrada en realidad relaciona con D -- en términos convectivos extra
(D·σ + σ·D - tr(D)·σ) que no se anulan salvo que σ=0.

En vez de corregir esa diferencia finita para que incorporara la
corrección de Truesdell (viable, pero más propenso a errores de signo),
se redefinió `tangent_moduli()` para devolver el tangente **lagrangiano**
D = dS/dE (2do tensor de Piola-Kirchhoff respecto a la deformación de
Green-Lagrange, ambos en la configuración de referencia). Esta cantidad
no tiene ninguna ambigüedad de objetividad: S y E son invariantes ante
una rotación superpuesta por construcción, así que la diferencia finita
directa (perturbando C, reconstruyendo una F representativa vía raíz
cuadrada simétrica de C, sin rotación) es correcta sin ningún término de
corrección. Esto también encaja mejor con el plan para `Quad4Element`
(Paso 2): una formulación Lagrangiana Total clásica (K = B_Lᵀ D B_L +
rigidez geométrica a partir de S), que evita por completo la necesidad
de una tasa objetiva en la implementación del elemento.

`NeoHookean::tangent_moduli()` ahora tiene la forma cerrada lagrangiana
correcta (derivada en el propio código, verificada a mano reduciéndola
al caso C=I, donde debe coincidir con el tensor de elasticidad lineal
isótropo estándar), cruzada contra `numerical_tangent_moduli()` en una
prueba unitaria -- y ambas coinciden. `OgdenHyperfoam::tangent_moduli()`
sigue usando el tangente numérico por defecto (ningún tangente analítico
lagrangiano se ha derivado todavía para Ogden; es el siguiente candidato
natural una vez que `Quad4Element` necesite uno más rápido que
diferencias finitas).

## Decisiones que quedaron pendientes a propósito

- **Elementos de orden superior / integración reducida**: no están en el
  esqueleto. Añadir cuando el Quad4 esté validado.
- **Contacto / autocontacto** (relevante si se quiere capturar el
  colapso de celdas a nivel micromecánico): no está contemplado todavía;
  sería un módulo nuevo bajo `elements/` o una capa aparte.
- **Paralelismo**: el ensamblaje y la factorización son, por ahora,
  secuenciales. `Eigen::SparseLU` es un solver directo de un solo hilo;
  si las mallas crecen, la ruta natural es un solver iterativo con
  precondicionador (o enlazar PETSc), sin cambiar la interfaz de
  `LinearSolver`.
