# foam-nlfem

Código de investigación en C++ para analizar y simular el comportamiento
mecánico no lineal de espumas poliméricas (hiperelasticidad compresible,
grandes deformaciones, inestabilidades tipo *snap-through* en la meseta
de colapso) mediante el método de elementos finitos, construido con el
mínimo de dependencias externas ("programación pura" + Eigen para álgebra
lineal).

## Estado actual

Este es un **esqueleto compilable**, no una implementación terminada.
Cada pieza está en uno de dos estados:

- **Implementada de verdad**: `Mesh`/`Node`/`ElementConnectivity`,
  `make_single_hex_element()`, `principal_invariants()`,
  `jacobian_determinant()`, `green_lagrange_strain()`,
  `NeoHookean::cauchy_stress()`, `LinearSolver` (wrapper sobre
  `Eigen::SparseLU`), `NewtonRaphsonSolver` (iteración genérica,
  desacoplada de la física vía `std::function`), `write_csv()`.
- **Declarada con `TODO` explícito** (lanzan `std::logic_error` si se
  llaman): `OgdenHyperfoam::cauchy_stress/tangent_moduli`,
  `NeoHookean::tangent_moduli`, `compute_deformation_gradient()`,
  `Quad4Element`, `Assembler`, `ArcLengthSolver`,
  `Problem::run_compression_test()`, `write_vtk()`,
  `make_structured_grid_2d()`.

Ver `docs/roadmap.md` para el orden recomendado en el que cerrar esos
`TODO`.

## Dependencias

- CMake ≥ 3.16
- Compilador C++17 (probado con GCC 13)
- Eigen3 ≥ 3.3. Si no hay una instalación de sistema, CMake usa
  automáticamente la copia header-only vendorizada en `external/eigen/`
  (Eigen 3.4.0, licencia MPL-2.0). No se necesita `apt-get`, vcpkg ni conan.
- doctest (framework de pruebas) — **ya viene incluido** como header
  único en `external/doctest/doctest.h` (v2.4.11, vendorizado, MIT), no hay
  que instalarlo aparte.

## Compilar y correr

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j

# correr las pruebas unitarias
ctest --output-on-failure

# correr el ejecutable de ejemplo (sanity check del solver, aún no la
# física de la espuma -- ver docs/roadmap.md)
./apps/single_element_test
```

## Estructura del repositorio

```
include/foam/        # headers públicos, organizados por capa
  core/               # alias de Eigen (Vector3, Matrix3, SparseMatrix...)
  mesh/               # Node, ElementConnectivity, Mesh, generadores de malla
  kinematics/         # gradiente de deformación, medidas de deformación
  materials/          # MaterialModel (interfaz), OgdenHyperfoam, NeoHookean
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

MIT (ver `LICENSE`). El código es propiedad del autor y no está ligado a
ninguna institución ni empleador. Los componentes vendorizados conservan
sus propias licencias: doctest es MIT y Eigen es MPL-2.0 (ver los avisos
en `LICENSE` y en `external/eigen/`).
