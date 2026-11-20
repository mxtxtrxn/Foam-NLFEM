#pragma once

#include <string>
#include <utility>
#include <vector>

#include "foam/core/types.hpp"

namespace foam {

// Writes a simple two-column CSV (e.g. engineering strain vs. nominal
// stress). This is the first "result" most of this project will
// produce -- the compression stress-strain curve to validate the
// material model against published/experimental data -- so it's worth
// having working from day one, independent of how far the FEM side is.
void write_csv(const std::string& path,
               const std::vector<std::pair<Real, Real>>& xy_data,
               const std::string& x_label = "strain",
               const std::string& y_label = "stress");

}  // namespace foam
