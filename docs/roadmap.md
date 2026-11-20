# Hoja de ruta

Orden recomendado para cerrar los `TODO` del esqueleto. Cada paso debe
quedar validado (con pruebas unitarias o contra datos publicados) antes
de pasar al siguiente -- el objetivo es nunca tener que depurar la física
del material y el solver de elementos finitos al mismo tiempo.

## Paso 1 -- Modelo constitutivo sobre un solo elemento

- [ ] Implementar `OgdenHyperfoam::cauchy_stress()`:
      descomponer F (o C = FᵀF) para obtener los estiramientos
      principales λ₁, λ₂, λ₃ y las direcciones principales; evaluar
      ∂W/∂λₖ por cada término de Ogden; reconstruir el tensor de
      esfuerzos en la base cartesiana.
- [ ] Validar esa función contra una curva de compresión uniaxial
      conocida (analítica o de `data/experimental/`) usando
      `make_single_hex_element()` con una deformación homogénea
      prescrita directamente (sin pasar todavía por `Assembler`).
- [ ] Implementar `OgdenHyperfoam::tangent_moduli()` y verificarlo con
      diferencias finitas sobre `cauchy_stress()` (ver
      `docs/architecture.md`, sección sobre el tangente numérico).
- [ ] Reemplazar `apps/single_element_test.cpp` (hoy un *sanity check*
      del solver con x²=2) por el test real de compresión homogénea.

## Paso 2 -- Malla 2D/axisimétrica y solver no lineal completo

- [ ] Implementar `compute_deformation_gradient()`.
- [ ] Implementar `Quad4Element::internal_force()` /
      `tangent_stiffness()` (cuadratura de Gauss 2x2).
- [ ] Implementar `make_structured_grid_2d()`.
- [ ] Implementar `Assembler::assemble_residual()` /
      `assemble_tangent()`.
- [ ] Implementar `Problem::run_compression_test()`, conectando
      `Assembler` + `NewtonRaphsonSolver` (éste ya está implementado) con
      un esquema de pasos de carga/desplazamiento.
- [ ] Prueba de parche (*patch test*) estándar: un campo de deformación
      uniforme debe reproducirse exactamente en la malla.
- [ ] Validar la curva tensión-deformación resultante contra datos
      experimentales reales de una espuma (PU, EPS, etc.) en
      `data/experimental/`.

## Paso 3 -- Inestabilidades y meseta de colapso

- [ ] Implementar `ArcLengthSolver` (Riks) una vez que Newton-Raphson
      puro empiece a fallar de verdad en la meseta (no antes -- es más
      fácil depurar longitud de arco cuando ya se sabe que todo lo demás
      funciona).
- [ ] Comparar la curva con y sin longitud de arco para confirmar que
      ambas coinciden donde Newton-Raphson puro sí converge.

## Paso 4 (opcional, en función de hacia dónde quieras llevar el proyecto)

- [ ] Dependencia de la velocidad de deformación / viscoelasticidad.
- [ ] Modelo micromecánico (celda unitaria tipo Kelvin o Voronoi,
      vigas representando los puntales) como alternativa o
      complemento al modelo fenomenológico de Ogden.
- [ ] `write_vtk()` para inspeccionar la malla deformada en ParaView.
- [ ] Ajuste de parámetros del material contra datos experimentales
      (mínimos cuadrados no lineales sobre μᵢ, αᵢ, βᵢ).
