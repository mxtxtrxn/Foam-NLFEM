#include "foam/io/csv_writer.hpp"

#include <fstream>
#include <stdexcept>

namespace foam {

void write_csv(const std::string& path, const std::vector<std::pair<Real, Real>>& xy_data,
               const std::string& x_label, const std::string& y_label) {
  std::ofstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("write_csv: could not open '" + path + "' for writing");
  }
  file << x_label << "," << y_label << "\n";
  for (const auto& [x, y] : xy_data) {
    file << x << "," << y << "\n";
  }
}

}  // namespace foam
