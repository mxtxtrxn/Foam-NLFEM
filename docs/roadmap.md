# Hoja de ruta

Orden recomendado para cerrar los `TODO` del esqueleto. Cada paso debe
quedar validado (con pruebas unitarias o contra datos publicados) antes
de pasar al siguiente -- el objetivo es nunca tener que depurar la física
del material y el solver de elementos finitos al mismo tiempo.

## Paso 1 -- Modelo constitutivo sobre un solo elemento ✅ (completado)

- [x] Implementar `OgdenHyperfoam::cauchy_stress()`: descomposición
      espectral de B = FFᵀ (autovectores = direcciones principales
      espaciales, autovalores = λₖ²) y la fórmula cerrada estándar de
      Ogden hyperfoam, σₖ = Σᵢ (2μᵢ/(αᵢJ))(λₖ^αᵢ - J^(-αᵢβᵢ)).
- [x] Validar de forma independiente (no solo contra sí misma):
      `tests/test_materials.cpp` cruza el esfuerzo cerrado contra una
      diferencia finita de `strain_energy()` (energía, no esfuerzo),
      verifica objetividad (σ' = R σ Rᵀ bajo F → RF), y confirma el caso
      especial β=0 (esfuerzo lateral exactamente cero sin estiramiento
      lateral -- la propiedad de "Poisson efectivo cero").
- [x] Implementar `OgdenHyperfoam::tangent_moduli()` -- vía
      `numerical_tangent_moduli()` (diferencias finitas), no una forma
      cerrada todavía. **Importante**: este tangente es el LAGRANGIANO
      dS/dE, no uno espacial -- ver `docs/architecture.md` para el
      porqué (un error real de la primera versión, no solo una
      simplificación de alcance).
- [x] En el camino: se construyó `foam/material_point/uniaxial_test.hpp`
      -- un ensayo de compresión uniaxial *libre* (superficie lateral sin
      confinar, no deformación uniaxial confinada) a nivel de un solo
      punto material, resolviendo el estiramiento lateral de equilibrio
      con el propio `NewtonRaphsonSolver` genérico. Esto generaliza
      "validar sobre `make_single_hex_element()`" a algo reusable y
      con pruebas propias (`tests/test_uniaxial.cpp`), sin esperar a
      `Assembler`.
- [x] Reemplazar `apps/single_element_test.cpp` (antes un *sanity check*
      del solver con x²=2) por el ensayo real: dos espumas de Ogden
      (β=0 y β=0.75) barridas hasta 70% de deformación de ingeniería,
      con reporte de módulo tangente inicial, esfuerzo en 10%/40%/máxima
      deformación, y estiramiento lateral -- confirma visualmente las
      tres regiones (codo elástico, meseta, inicio de densificación).

**Pendiente opcional dentro de este mismo paso** (no bloquea el Paso 2):
derivar una forma cerrada para `OgdenHyperfoam::tangent_moduli()`
(lagrangiana, dS/dE) que reemplace el tangente numérico -- más rápido,
pero más fácil de equivocar; el numérico seguirá sirviendo como
referencia de validación una vez exista.

**Parámetros usados**: μ=80 kPa, α=13, β∈{0, 0.75} -- *ilustrativos*,
elegidos para producir la forma cualitativa correcta, no ajustados a un
material real publicado. Ajustar contra datos experimentales reales
sigue siendo el Paso 4.

## Paso 2 -- Malla 2D y solver no lineal completo ✅ (completado)

Formulación Lagrangiana Total, consistente con `tangent_moduli()` =
dS/dE del Paso 1: fuerza interna vía S (2do Piola-Kirchhoff) y el
operador deformación-desplazamiento no lineal B_L, rigidez = B_Lᵀ D B_L
(parte material) + rigidez geométrica construida directamente de S. Sin
necesidad de ninguna tasa de esfuerzo objetiva dentro del elemento.

- [x] `compute_deformation_gradient()`: generalizado para dim=2 (con
      embebido de deformación plana, F33=1) y dim=3.
- [x] `second_piola_kirchhoff_stress()`: nueva utilidad pública en
      `materials/material_base.hpp` (antes interna, usada solo por el
      tangente numérico); ahora también la usa `Quad4Element`.
- [x] `make_structured_grid_2d()`: malla estructurada de Quad4,
      orientación CCW consistente con las funciones de forma.
- [x] `Quad4Element::internal_force()` / `tangent_stiffness()`:
      cuadratura de Gauss 2x2, B_L derivado y verificado a mano (ver
      comentarios en `src/elements/quad4.cpp`), rigidez geométrica
      derivada por separado.
- [x] `Assembler::assemble_residual()` / `assemble_tangent()`: ensamblaje
      estándar vía *scatter-add* / triplets.
- [x] `Problem::run_compression_test()`: eliminación de grados de
      libertad restringidos (Dirichlet), pasos de carga con
      *continuation* (la solución convergida de un paso arranca el
      siguiente), esfuerzo nominal vía fuerza de reacción en el borde
      cargado. Incluye `uniaxial_compression_bcs()`, un helper para
      montar las condiciones de frontera de un ensayo de compresión
      libre sobre una malla estructurada.

**Validaciones** (cada una más fuerte que la anterior):
- Consistencia tangente↔fuerza por diferencias finitas a nivel de
  elemento (`tests/test_quad4.cpp`) -- el análogo, para el elemento, del
  chequeo que atrapó el bug del tangente en el Paso 1.
- *Patch test* clásico: una deformación homogénea aplicada en toda la
  malla da fuerza interna ≈0 en todo nodo interior
  (`tests/test_assembler.cpp`).
- **Cross-check de punta a punta**: un solo elemento Quad4, con los
  desplazamientos nodales exactos de la solución de compresión libre del
  Paso 1, da fuerza lateral ≈0 (ya que β=0 implica S22=0) y la fuerza
  axial total en el borde cargado coincide con el esfuerzo nominal del
  Paso 1 hasta precisión de máquina (`tests/test_quad4.cpp`).
- **El mismo cross-check, pero con una malla completa** (3x6 elementos,
  `Assembler`+`Problem` de punta a punta, no un solo elemento a mano):
  coincide con el resultado del punto material del Paso 1 en cada
  escalón de la curva, con diferencia ~0.00% (`tests/test_problem.cpp`,
  `apps/mesh_compression_test.cpp`).

**Pendiente, movido fuera de este paso**: la prueba de parche ya está
hecha (ver arriba); validar contra datos experimentales reales de una
espuma sigue siendo el Paso 4, ya que los parámetros usados aquí siguen
siendo ilustrativos (ver nota del Paso 1).

## Paso 3 -- Inestabilidades y meseta de colapso ✅ (completado, con una advertencia honesta)

- [x] **Antes de implementar nada**: se verificó numéricamente que el
      material de Ogden hyperfoam usado hasta ahora (μ=80 kPa, α=13,
      β=0) **no tiene ablandamiento real** -- la rigidez tangente local
      baja mucho en la meseta pero nunca cruza a negativo (mínimo
      verificado ≈ -4.5e6, siempre del signo esperado, en todo el rango
      de deformación barrido). Es decir: con estos parámetros, el
      Newton-Raphson con control de desplazamiento del Paso 2 ya basta;
      no hay snap-through real que forzar.
- [x] `ArcLengthSolver` implementado de todas formas, con una decisión
      de diseño importante: **se rediseñó la interfaz** respecto a la
      declaración original del Paso 1 (que asumía control de carga puro,
      R(u,λ)=F_int(u)-λ·F_ref). En vez de eso, toma R(u,λ) como una
      función abstracta y calcula ∂R/∂λ internamente por diferencias
      finitas -- mismo espíritu que `numerical_tangent_moduli()` -- lo
      cual funciona igual de bien para control de carga que para control
      de desplazamiento (que es lo que usa `Problem`), sin que el solver
      necesite saber cuál de los dos es.
- [x] Método: Riks/Wempner de "plano normal" (constraint lineal, no la
      esférica exacta de Crisfield con selección de raíz cuadrática --
      más simple de implementar correctamente, con el mismo resultado en
      la práctica). Con reducción automática de paso ante un fallo
      (`max_retries_per_step`), una mejora de robustez estándar.
- [x] **Validado contra un problema con solución cerrada conocida**
      (no contra el modelo de espuma, que no tiene inestabilidad que
      exhibir -- ver arriba): R(u,λ)=u³-3u-λ=0, cuya forma cerrada
      λ(u)=u³-3u es la curva S clásica de snap-through, con puntos
      límite en u=-1 (λ=2) y u=+1 (λ=-2). El solver traza la curva
      completa, incluida la rama inestable entre ambos puntos límite,
      coincidiendo con la forma cerrada hasta la tolerancia de Newton en
      cada paso. Como contraste, `apps/snap_through_demo` muestra que un
      Newton-Raphson con control de carga puro (el `NewtonRaphsonSolver`
      del Paso 1, incrementando λ de forma monótona) efectivamente falla
      al llegar al punto límite (λ≈2).

**Pendiente, honestamente fuera de alcance por ahora**: no se aplicó
`ArcLengthSolver` al modelo de espuma real, porque los parámetros
actuales no lo necesitan (ver arriba). Para una demostración foam-a-foam
de esto haría falta o bien una ley de material con ablandamiento real, o
bien una malla fina donde aparezca pandeo/localización geométrica --
cualquiera de las dos es un paso natural pero no trivial, no forzado
aquí para evitar un ejemplo artificioso.

## Paso 4 (opcional, en función de hacia dónde quieras llevar el proyecto)

- [ ] Dependencia de la velocidad de deformación / viscoelasticidad.
- [ ] Modelo micromecánico (celda unitaria tipo Kelvin o Voronoi,
      vigas representando los puntales) como alternativa o
      complemento al modelo fenomenológico de Ogden.
- [ ] `write_vtk()` para inspeccionar la malla deformada en ParaView.
- [ ] Ajuste de parámetros del material contra datos experimentales
      (mínimos cuadrados no lineales sobre μᵢ, αᵢ, βᵢ).
