#include "foam/io/vtk_writer.hpp"

#include <stdexcept>

namespace foam {

void write_vtk(const std::string& /*path*/, const Mesh& /*mesh*/,
               const VectorX& /*nodal_displacements*/) {
  // TODO(roadmap step 2): write POINTS/CELLS/POINT_DATA sections in
  // legacy VTK ASCII format. Keep it minimal at first (displacement
  // vectors only); add CELL_DATA for stress once elements report it.
  throw std::logic_error("write_vtk: not implemented yet");
}

}  // namespace foam
