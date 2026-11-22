# foam-nlfem

Código de investigación en C++ para analizar y simular el comportamiento
mecánico no lineal de espumas poliméricas (hiperelasticidad compresible,
grandes deformaciones, inestabilidades tipo *snap-through* en la meseta
de colapso) mediante el método de elementos finitos, construido con el
mínimo de dependencias externas ("programación pura" + Eigen para álgebra
lineal).

## Estado actual

**Pasos 1 y 2 de la hoja de ruta están completos**: hay un pipeline de
elementos finitos no lineales de punta a punta (malla → elemento →
ensamblaje → Newton-Raphson con condiciones de frontera) que reproduce
exactamente el resultado de un ensayo de un solo punto material.

- **Implementada y validada**: todo lo del Paso 1 (ver `docs/roadmap.md`
  para el detalle), más -- del Paso 2 -- `compute_deformation_gradient()`
  (2D con embebido de deformación plana, y 3D),
  `second_piola_kirchhoff_stress()` (nueva utilidad pública),
  `make_structured_grid_2d()`, **`Quad4Element`** (formulación
  Lagrangiana Total completa: fuerza interna + rigidez material +
  rigidez geométrica), **`Assembler`** (ensamblaje estándar vía
  *scatter-add*), y **`Problem::run_compression_test()`** (eliminación
  de grados de libertad Dirichlet, pasos de carga con *continuation*).
- **Declarada con `TODO` explícito**: `ArcLengthSolver` (Paso 3),
  `write_vtk()` (Paso 4, opcional).

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
