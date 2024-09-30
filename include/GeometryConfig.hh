#pragma once
#include <yaml-cpp/yaml.h>

#include <string>
#include <unordered_map>

#include "G4VisAttributes.hh"

class GeometryConfig {
public:
  static void LoadMaterials(const char *path);
  static void LoadVolumes(const char *path);

private:
  YAML::Node node_;
  static std::unordered_map<std::string, G4VisAttributes> fMaterialVisAttributes;

  explicit GeometryConfig(const char *path);
  void ProcessMaterials();
  void ProcessVolumes();
};
