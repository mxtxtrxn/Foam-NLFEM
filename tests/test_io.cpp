#include <cstdio>
#include <fstream>
#include <string>

#include "doctest.h"
#include "foam/io/csv_writer.hpp"

TEST_CASE("write_csv: writes header and rows in the expected format") {
  const std::string path = "test_write_csv_output.csv";
  foam::write_csv(path, {{0.0, 0.0}, {0.1, 1.5}}, "strain", "stress");

  std::ifstream file(path);
  REQUIRE(file.is_open());

  std::string line;
  std::getline(file, line);
  CHECK(line == "strain,stress");
  std::getline(file, line);
  CHECK(line == "0,0");
  std::getline(file, line);
  CHECK(line == "0.1,1.5");

  file.close();
  std::remove(path.c_str());
}
