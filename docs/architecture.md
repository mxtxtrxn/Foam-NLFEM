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

## Por qué un tangente numérico es una opción legítima al principio

Derivar a mano el tangente consistente de una energía tipo Ogden es
tedioso y fácil de equivocar. Antes de confiar en una derivación
analítica dentro de Newton-Raphson, vale la pena:

1. Implementar y validar `cauchy_stress()` solo (más fácil de verificar).
2. Obtener un tangente por diferencias finitas sobre `cauchy_stress()`
   (más lento, pero correcto por construcción si el esfuerzo lo es).
3. Usar ese tangente numérico como *referencia* para validar, con una
   prueba unitaria, el tangente analítico una vez esté escrito.

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
