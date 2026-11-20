#pragma once

#include <string>

#include "foam/core/types.hpp"
#include "foam/mesh/mesh.hpp"

namespace foam {

// Exports the mesh plus a nodal displacement field (and, later, element
// stresses) to legacy VTK (.vtk) format for visualization in ParaView.
// Only worth building once there's an actual 2D/3D mesh to look at --
// deferred until roadmap step 2 is producing real deformed geometry.
void write_vtk(const std::string& path, const Mesh& mesh, const VectorX& nodal_displacements);

}  // namespace foam
