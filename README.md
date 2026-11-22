# foam-nlfem

Código de investigación en C++ para analizar y simular el comportamiento
mecánico no lineal de espumas poliméricas (hiperelasticidad compresible,
grandes deformaciones, inestabilidades tipo *snap-through* en la meseta
de colapso) mediante el método de elementos finitos, construido con el
mínimo de dependencias externas ("programación pura" + Eigen para álgebra
lineal).

## Estado actual

**Pasos 1, 2 y 3 de la hoja de ruta están completos.** Hay un pipeline
de elementos finitos no lineales de punta a punta (malla → elemento →
ensamblaje → Newton-Raphson con condiciones de frontera) que reproduce
exactamente el resultado de un ensayo de un solo punto material, más un
solver de longitud de arco para cuando haga falta trazar inestabilidades
que el control de desplazamiento puro no pueda seguir.

- **Implementada y validada**: todo lo de los Pasos 1 y 2 (ver
  `docs/roadmap.md` para el detalle), más -- del Paso 3 --
  **`ArcLengthSolver`** (Riks/Wempner, interfaz genérica vía ∂R/∂λ por
  diferencias finitas, reducción automática de paso ante un fallo).
  Validado contra un problema con solución cerrada conocida (curva S
  clásica de snap-through), no contra el modelo de espuma: se verificó
  primero que el material ilustrativo usado hasta ahora **no tiene
  ablandamiento real** (la rigidez baja mucho en la meseta pero nunca
  cruza a negativo), así que el Newton-Raphson con control de
  desplazamiento del Paso 2 ya le basta -- ver la nota completa en
  `docs/roadmap.md`, Paso 3.
- **Declarado con `TODO` explícito**: `write_vtk()` (Paso 4, opcional).

**Nota importante sobre `tangent_moduli()`**: es el tangente
LAGRANGIANO dS/dE (2do Piola-Kirchhoff respecto a Green-Lagrange), no
un tangente espacial sobre el esfuerzo de Cauchy. La primera versión
intentó esto último y tenía un error real (confundía la tasa simple con
la tasa objetiva de Truesdell); ver `docs/architecture.md` para el
detalle completo. `Quad4Element` (Paso 2) ya está formulado como
Lagrangiano Total para ser consistente con esta definición -- vale la
pena leer esa nota antes de modificar su código.

Ver `docs/roadmap.md` para el detalle de lo hecho y el orden
recomendado para el resto.

## Dependencias

- CMake ≥ 3.16
- Compilador C++17 (probado con GCC 13)
- Eigen3 ≥ 3.3
  - Ubuntu/Debian: `sudo apt install libeigen3-dev`
  - macOS: `brew install eigen`
  - o vía vcpkg/conan si prefieres un gestor de paquetes
  - si no está instalado, CMake usa la copia vendorizada en
    `external/eigen` (Eigen 3.4.0, MPL-2.0)
- doctest (framework de pruebas) — **ya viene incluido** como header
  único en `external/doctest/doctest.h` (v2.4.11, vendorizado), no hay
  que instalarlo aparte.

## Compilar y correr

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j

# correr las pruebas unitarias
ctest --output-on-failure

# correr el ejecutable de ejemplo: valida OgdenHyperfoam con un ensayo
# de compresión uniaxial libre (dos variantes, beta=0 y beta=0.75),
# imprime módulo tangente inicial / esfuerzo a 10%-40%-70% de
# deformación, y escribe foam_beta0_curve.csv / foam_beta_pos_curve.csv
./apps/single_element_test

# correr el ensayo de compresión sobre una malla 2D real (Paso 2):
# compara, escalón a escalón, contra el resultado de un solo punto
# material (Paso 1) -- deberían coincidir casi exactamente
./apps/mesh_compression_test

# demo de ArcLengthSolver (Paso 3) sobre la curva S clásica de
# snap-through: muestra que el control de carga puro falla justo en el
# punto límite, y que la longitud de arco traza la curva completa
./apps/snap_through_demo

# graficar una curva (opcional, requiere matplotlib)
python3 ../../scripts/plot_stress_strain.py apps/foam_beta0_curve.csv
```

## Estructura del repositorio

```
include/foam/        # headers públicos, organizados por capa
  core/               # alias de Eigen (Vector3, Matrix3, SparseMatrix...)
  mesh/               # Node, ElementConnectivity, Mesh, generadores de malla
  kinematics/         # gradiente de deformación, medidas de deformación
  materials/          # MaterialModel (interfaz), OgdenHyperfoam, NeoHookean
  material_point/     # ensayos a nivel de un solo punto material (sin malla):
                       # uniaxial_test.hpp, usado para validar el Paso 1
  elements/           # FiniteElement (interfaz), Quad4Element
  assembly/           # Assembler: malla + elemento + material -> sistema global
  boundary_conditions/# DirichletBC, NeumannBC (datos puros)
  solvers/            # LinearSolver, NewtonRaphsonSolver, ArcLengthSolver
  io/                 # write_csv, write_vtk
  simulation/         # Problem: amarra todo lo anterior

src/                  # implementaciones (.cpp), misma organización que include/
apps/                 # ejecutables de punta a punta
tests/                # pruebas unitarias (doctest)
external/doctest/     # doctest vendorizado (single header)
data/experimental/    # curvas experimentales publicadas, para validar contra ellas
data/meshes/          # mallas de entrada (cuando se use una malla real)
scripts/              # utilidades en Python (graficar resultados, etc.) -- opcional
docs/                 # architecture.md, roadmap.md
```

**Por qué esta separación:** la malla (`mesh`) no sabe nada de elementos
finitos ni de materiales; un `FiniteElement` no sabe nada de
`MaterialModel` salvo a través de la interfaz abstracta
`materials/material_base.hpp`; el `Assembler` es el único lugar que
conoce la numeración global de grados de libertad; y
`NewtonRaphsonSolver`/`LinearSolver` no saben nada de FEM -- operan sobre
funciones `residual(x)`/`tangent(x)` abstractas. Esto permite implementar
y probar cada capa por separado, y cambiar una pieza (p. ej. el modelo de
material, o pasar de Newton-Raphson a longitud de arco) sin tocar las
demás.

## Licencia

MIT (ver `LICENSE`). El header vendorizado en
`external/doctest/doctest.h` mantiene su propia licencia MIT
(Copyright (c) 2016-2023 Viktor Kirilov).
