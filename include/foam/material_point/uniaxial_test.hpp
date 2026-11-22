#pragma once

#include <vector>

#include "foam/materials/material_base.hpp"

namespace foam {

// Simulates an unconfined ("free") uniaxial compression/tension test at
// a SINGLE material point -- no mesh, no Assembler, just a material
// model and F = diag(axial_stretch, lateral_stretch, lateral_stretch).
// This is what a real compression-test specimen does: the axial
// stretch is prescribed by the test machine, but the lateral surfaces
// are free, so the lateral stretch is whatever makes the lateral
// (transverse) Cauchy stress vanish -- NOT fixed at 1, which would be
// uniaxial STRAIN (a confined, different, and much easier to get wrong
// test case to validate against real compression data).
//
// This is roadmap step 1's validation harness: it deliberately stays
// below the FEM machinery (no Assembler/Quad4Element needed) while
// still exercising the generic NewtonRaphsonSolver, same as it will be
// reused for later.
struct UniaxialPoint {
  Real axial_stretch;
  Real lateral_stretch;
  Real engineering_strain;      // 1 - axial_stretch (positive in compression)
  Real cauchy_stress_axial;     // sigma_11
  Real nominal_stress_axial;    // P_11 (1st Piola-Kirchhoff), what a real
                                 // test machine reports as force/original area
  bool converged;
};

// Solves for the lateral stretch at a single prescribed axial stretch.
UniaxialPoint solve_uniaxial_point(const MaterialModel& material, Real axial_stretch,
                                    Real initial_lateral_guess = 1.0);

// Sweeps axial_stretch from axial_stretch_min to axial_stretch_max over
// num_steps points, warm-starting each lateral-stretch solve from the
// previous step's converged value (continuation), which is both faster
// and more robust than guessing 1.0 at every step once the material is
// well into the densification region.
std::vector<UniaxialPoint> uniaxial_compression_curve(const MaterialModel& material,
                                                        Real axial_stretch_min,
                                                        Real axial_stretch_max, int num_steps);

}  // namespace foam
