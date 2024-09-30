#include "GeometryConfig.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Color.hh"
#include "G4Element.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

using namespace std;

namespace {

G4LogicalVolume *CreateBoxVolume(const string &name, G4double hx, G4double hy, G4double hz, G4Material *material)
{
  G4cout << " * Box " << name << ": " << hx * 2 / mm << ", " << hy * 2 / mm << ", " << hz * 2 / mm << " (mm)" << G4endl;
  auto box = new G4Box(name, hx, hy, hz);
  return new G4LogicalVolume(box, material, name);
}

pair<bool, G4double> ParsePhysicsVariable(const string &variable)
{
  G4double value;
  G4String unit;
  istringstream iss(variable);

  if(!(iss >> value)) {
    G4cerr << "ERROR: failed to parse physics variable: " << variable << G4endl;
    exit(EXIT_FAILURE);
  }
  if(iss >> unit) {
    if(unit == "%") return { true, value / 100 };
    value *= G4UnitDefinition::GetValueOf(unit);
  }
  return { false, value };
}

G4double ParseAbsolutePhysicsVariable(const string &variable)
{
  auto [relevant, value] = ParsePhysicsVariable(variable);
  if(relevant) {
    G4cerr << "ERROR: expect absolute value: " << variable << G4endl;
    exit(EXIT_FAILURE);
  }
  return value;
}

G4Material *ParseMaterial(const string &name)
{
  G4Material *material = G4Material::GetMaterial(name, false);
  if(!material) {
    G4cerr << "ERROR: unknown material: " << name << G4endl;
    exit(EXIT_FAILURE);
  }
  return material;
}

G4Element *ParseElement(const string &name)
{
  G4Element *element = G4Element::GetElement(name, false);
  if(!element) {
    G4cerr << "ERROR: unknown element: " << name << G4endl;
    exit(EXIT_FAILURE);
  }
  return element;
}

G4LogicalVolume *ProcessBox(const string &name, YAML::Node node)
{
  G4double x = ParseAbsolutePhysicsVariable(node["x"].as<string>());
  G4double y = ParseAbsolutePhysicsVariable(node["y"].as<string>());
  G4double z = ParseAbsolutePhysicsVariable(node["z"].as<string>());
  G4Material *material = ParseMaterial(node["material"].as<string>());
  return CreateBoxVolume(name, x / 2, y / 2, z / 2, material);
}

vector<G4LogicalVolume *> ProcessStackComponents(
    YAML::Node node, G4double &hx_s, G4double &hy_s, G4double &hz_s, G4double &hx_m, G4double &hy_m, G4double &hz_m)
{
  hx_s = hy_s = hz_s = 0.0;
  hx_m = hy_m = hz_m = 0.0;
  vector<G4LogicalVolume *> children;
  children.reserve(node["components"].size());
  for(YAML::Node child_name : node["components"]) {
    G4LogicalVolume *child = G4LogicalVolumeStore::GetInstance()->GetVolume(child_name.as<string>());
    if(!child) {
      G4cerr << "ERROR: unknown logical volume: " << child_name.as<string>() << G4endl;
      exit(EXIT_FAILURE);
    }
    G4Box *box = dynamic_cast<G4Box *>(child->GetSolid());
    if(!box) {
      G4cerr << "ERROR: expect box component" << G4endl;
      exit(EXIT_FAILURE);
    }
    children.push_back(child);
    hx_s += box->GetXHalfLength();
    hy_s += box->GetYHalfLength();
    hz_s += box->GetZHalfLength();
    hx_m = max(hx_m, box->GetXHalfLength());
    hy_m = max(hy_m, box->GetYHalfLength());
    hz_m = max(hz_m, box->GetZHalfLength());
  }
  return children;
}

void ProcessStackSize(const string &name, G4double &current, G4double target)
{
  if(current <= target) {
    current = target;
    return;
  }
  G4cerr << "ERROR: " << name << " size " << target * 2 << " smaller than needed " << current * 2 << G4endl;
  exit(EXIT_FAILURE);
}

void ProcessStackSize(const string &name, YAML::Node node, G4double &hx, G4double &hy, G4double &hz)
{
  if(node["padding"]) {
    auto [relevant, padding] = ParsePhysicsVariable(node["padding"].as<string>());
    if(relevant) {
      hx *= 1 + padding;
      hy *= 1 + padding;
      hz *= 1 + padding;
    } else {
      hx += padding;
      hy += padding;
      hz += padding;
    }
  } else {
    if(node["x"]) { ProcessStackSize(name + ":x", hx, ParseAbsolutePhysicsVariable(node["x"].as<string>()) / 2); }
    if(node["y"]) { ProcessStackSize(name + ":y", hy, ParseAbsolutePhysicsVariable(node["y"].as<string>()) / 2); }
    if(node["z"]) { ProcessStackSize(name + ":z", hz, ParseAbsolutePhysicsVariable(node["z"].as<string>()) / 2); }
  }
}

void ProcessOffset(
    const string &name, YAML::Node node, G4double hx, G4double hy, G4double hz, G4double &x, G4double &y, G4double &z)
{
  if(!node["offset"]) { return; }
  if(node["offset"].size() != 3) {
    G4cerr << "ERROR: " << name << ": offset not of length 3" << G4endl;
    exit(EXIT_FAILURE);
  }
  G4double *h[3] = { &hx, &hy, &hz };
  G4double *p[3] = { &x, &y, &z };
  for(size_t i = 0; i < 3; ++i) {
    auto [relevant, offset] = ParsePhysicsVariable(node["offset"][i].as<string>());
    *p[i] += relevant ? *h[i] * offset : offset;
  }
}

G4LogicalVolume *ProcessBottomUp(const string &name, YAML::Node node)
{
  G4double hx_s, hy_s, hz_s, hx_m, hy_m, hz_m;
  auto children = ProcessStackComponents(node, hx_s, hy_s, hz_s, hx_m, hy_m, hz_m);
  G4double hx = hx_m, hy = hy_m, hz = hz_s;
  size_t duplicate = 1;
  if(node["duplicate"]) { duplicate = node["duplicate"].as<size_t>(); }
  hz *= duplicate;
  G4double x = 0.0, y = 0.0, z = -hz;
  ProcessStackSize(name, node, hx, hy, hz);
  ProcessOffset(name, node, hx, hy, hz, x, y, z);
  G4Material *material = ParseMaterial(node["material"].as<string>());

  auto logical = CreateBoxVolume(name, hx, hy, hz, material);
  size_t i = 0;
  for(size_t d = 0; d < duplicate; ++d)
    for(G4LogicalVolume *child : children) {
      G4double ht = ((G4Box *)child->GetSolid())->GetZHalfLength();
      string child_name = name + "_" + to_string(i++) + ":" + child->GetName();
      new G4PVPlacement(0, { x, y, z + ht }, child, child_name, logical, false, 0, true);
      z += ht * 2;
    }
  return logical;
}

G4LogicalVolume *ProcessLeftRight(const string &name, YAML::Node node)
{
  G4double hx_s, hy_s, hz_s, hx_m, hy_m, hz_m;
  auto children = ProcessStackComponents(node, hx_s, hy_s, hz_s, hx_m, hy_m, hz_m);
  G4double hx = hx_s, hy = hy_m, hz = hz_m;
  size_t duplicate = 1;
  if(node["duplicate"]) { duplicate = node["duplicate"].as<size_t>(); }
  hx *= duplicate;
  G4double x = -hx, y = 0.0, z = 0.0;
  ProcessStackSize(name, node, hx, hy, hz);
  ProcessOffset(name, node, hx, hy, hz, x, y, z);
  G4Material *material = ParseMaterial(node["material"].as<string>());

  auto logical = CreateBoxVolume(name, hx, hy, hz, material);
  size_t i = 0;
  for(size_t d = 0; d < duplicate; ++d)
    for(G4LogicalVolume *child : children) {
      G4double ht = ((G4Box *)child->GetSolid())->GetXHalfLength();
      string child_name = name + "_" + to_string(i++) + ":" + child->GetName();
      new G4PVPlacement(0, { x + ht, y, z }, child, child_name, logical, false, 0, true);
      x += ht * 2;
    }
  return logical;
}

void ProcessRotation(const string &name, G4RotationMatrix *rotation, const string &axis, G4double degree)
{
  if(degree != (int)degree || (int)degree % 90) {
    G4cerr << "ERROR: " << name << ": Rotation degree must be multiple of 90: " << degree << G4endl;
    exit(EXIT_FAILURE);
  }
  if(axis == "x") {
    rotation->rotateX(degree * CLHEP::deg);
  } else if(axis == "y") {
    rotation->rotateY(degree * CLHEP::deg);
  } else if(axis == "z") {
    rotation->rotateZ(degree * CLHEP::deg);
  } else {
    G4cerr << "ERROR: " << name << ": Unknown rotation axis: " << axis << G4endl;
    exit(EXIT_FAILURE);
  }
}

G4LogicalVolume *ProcessRotation(const string &name, YAML::Node node)
{
  G4LogicalVolume *child = G4LogicalVolumeStore::GetInstance()->GetVolume(node["components"][0].as<string>());
  G4Box *box = dynamic_cast<G4Box *>(child->GetSolid());
  if(!box) {
    G4cerr << "ERROR: expect box component" << G4endl;
    exit(EXIT_FAILURE);
  }
  auto rotation = new G4RotationMatrix();
  for(size_t i = 1; i < node["components"].size(); ++i) {
    YAML::Node item = node["components"][i];
    G4double degree = ParseAbsolutePhysicsVariable(item[1].as<string>()) / CLHEP::deg;
    ProcessRotation(name, rotation, item[0].as<string>(), degree);
  }
  G4ThreeVector v(box->GetXHalfLength(), box->GetYHalfLength(), box->GetZHalfLength());
  v = *rotation * v;
  auto logical = CreateBoxVolume(name, fabs(v.x()), fabs(v.y()), fabs(v.z()), child->GetMaterial());
  node["material"] = (string)logical->GetMaterial()->GetName();
  string child_name = name + "_0:" + child->GetName();
  new G4PVPlacement(rotation, { 0, 0, 0 }, child, child_name, logical, false, 0, true);  // rotation owned by us
  G4AutoDelete::Register(rotation);
  return logical;
}

void ProcessVisAttributes(YAML::Node node, G4VisAttributes &attr)
{
  if(node["hidden"]) {
    attr.SetVisibility(false);
    return;
  }
  if(node["color"]) {
    G4Color color;
    if(!G4Color::GetColour(node["color"].as<string>(), color)) {
      G4cerr << "ERROR: unknown color: " << node["color"].as<string>() << G4endl;
      exit(EXIT_FAILURE);
    }
    if(node["alpha"]) { color.SetAlpha(node["alpha"].as<G4double>()); }
    attr.SetColor(color);
  }
  if(node["line_style"]) { attr.SetLineStyle((G4VisAttributes::LineStyle)node["line_style"].as<size_t>()); }
  if(node["line_width"]) { attr.SetLineWidth(ParseAbsolutePhysicsVariable(node["line_width"].as<string>())); }
  if(node["color"] && !node["line_style"] && !node["line_width"]) {
    attr.SetLineWidth(0);
    attr.SetForceSolid();
  }
}

void ProcessNist(const string &name, YAML::Node node)
{
  G4NistManager *manager = G4NistManager::Instance();
  if(node["type"] && node["type"].as<string>() == "element") {
    if(!manager->FindOrBuildElement(name)) {
      G4cerr << "ERROR: unknown NIST element: " << name << G4endl;
      exit(EXIT_FAILURE);
    }
  } else {
    if(!manager->FindOrBuildMaterial(name)) {
      G4cerr << "ERROR: unknown NIST material: " << name << G4endl;
      exit(EXIT_FAILURE);
    }
  }
}

void AddElement(G4Material *material, G4Element *child, YAML::Node node)
{
  G4cout << " * " << material->GetName() << " <- element " << child->GetName() << ", " << node.as<string>() << G4endl;
  try {
    material->AddElement(child, node.as<G4int>());
  } catch(const YAML::BadConversion &) {
    material->AddElement(child, node.as<G4double>());
  }
}

void AddMaterial(G4Material *material, G4Material *child, YAML::Node node)
{
  G4cout << " * " << material->GetName() << " <- material " << child->GetName() << ", " << node.as<string>() << G4endl;
  material->AddMaterial(child, node.as<G4double>());
}

void AddComponent(G4Material *material, YAML::Node node)
{
  auto it = node.begin();
  string type = "material";
  if(node.size() == 3) {
    type = it++->as<string>();
  } else if(node.size() != 2) {
    G4cerr << "ERROR: expect material component of size 2: " << node << G4endl;
    exit(EXIT_FAILURE);
  }
  if(type != "material" && type != "element") {
    G4cerr << "ERROR: unknown material type: " << type << G4endl;
    exit(EXIT_FAILURE);
  }
  string name = it++->as<string>();
  if(type == "material") {
    G4Material *child = ParseMaterial(name);
    AddMaterial(material, child, *it);
  } else {
    G4Element *child = ParseElement(name);
    AddElement(material, child, *it);
  }
}

G4Material *ProcessMaterial(const string &name, YAML::Node node)
{
  if(node["alias"]) {
    G4cout << " * " << node["alias"].as<string>() << " -> " << name << G4endl;
    G4Material *child = ParseMaterial(node["alias"].as<string>());
    G4Material *material = new G4Material(name, child->GetDensity(), 1);
    G4AutoDelete::Register(material);  // owned by us
    material->AddMaterial(child, 1.0);
    return material;
  }
  G4double density = ParseAbsolutePhysicsVariable(node["density"].as<string>());
  G4cout << " * " << name << " <- density " << density / (g / cm3) << " g/cm3" << G4endl;
  G4Material *material = new G4Material(name, density, node["components"].size());
  G4AutoDelete::Register(material);  // owned by us
  for(YAML::Node component : node["components"]) { AddComponent(material, component); }
  return material;
}

G4Element *ProcessElement(const string &name, YAML::Node)
{
  G4cerr << "ERROR: element mixture unimplemented: " << name << G4endl;
  exit(EXIT_FAILURE);
}

}  // namespace

unordered_map<string, G4VisAttributes> GeometryConfig::fMaterialVisAttributes;

GeometryConfig::GeometryConfig(const char *path)
{
  ifstream file(path);
  node_ = YAML::LoadAll(file);
}

void GeometryConfig::LoadVolumes(const char *path)
{
  G4cout << "Loading volumes from " << path << G4endl;
  GeometryConfig(path).ProcessVolumes();
}

void GeometryConfig::LoadMaterials(const char *path)
{
  G4cout << "Loading materials from " << path << G4endl;
  GeometryConfig(path).ProcessMaterials();
}

void GeometryConfig::ProcessVolumes()
{
  for(const pair<YAML::Node, YAML::Node> &pair : node_[1]) {
    auto &[name, node] = pair;
    if(node["alternative"]) {
      if(G4LogicalVolumeStore::GetInstance()->GetVolume(name.as<string>(), false)) { continue; }
    }
    G4cout << "Building volume " << name << G4endl;
    G4LogicalVolume *logical;
    if(node["solid"].as<string>() == "box") {
      logical = ProcessBox(name.as<string>(), node);
    } else if(node["solid"].as<string>() == "bottom_up") {
      logical = ProcessBottomUp(name.as<string>(), node);
    } else if(node["solid"].as<string>() == "left_right") {
      logical = ProcessLeftRight(name.as<string>(), node);
    } else if(node["solid"].as<string>() == "rotation") {
      logical = ProcessRotation(name.as<string>(), node);
    } else {
      G4cerr << "ERROR: unknown solid type: " << node["solid"].as<string>() << G4endl;
      exit(EXIT_FAILURE);
    }
    G4VisAttributes attr = fMaterialVisAttributes[node["material"].as<string>()];
    ProcessVisAttributes(node, attr);
    logical->SetVisAttributes(attr);
  }
}

void GeometryConfig::ProcessMaterials()
{
  for(const pair<YAML::Node, YAML::Node> &pair : node_[1]) {
    auto &[name, node] = pair;
    string type = "material";
    if(node["type"]) { type = node["type"].as<string>(); }
    if(type != "material" && type != "element") {
      G4cerr << "ERROR: unknown material type: " << type << G4endl;
      exit(EXIT_FAILURE);
    }
    if(node["alternative"]) {
      if(type == "material") {
        if(G4Material::GetMaterial(name.as<string>(), false)) continue;
      } else {
        if(G4Element::GetElement(name.as<string>(), false)) continue;
      }
    }
    G4cout << "Building " << type << " " << name << G4endl;
    if(node["from"]) {
      if(node["from"].as<string>() == "nist") {
        ProcessNist(name.as<string>(), node);
      } else {
        G4cerr << "ERROR: unknown material source: " << node["from"].as<string>() << G4endl;
        exit(EXIT_FAILURE);
      }
    } else if(type == "material") {
      ProcessMaterial(name.as<string>(), node);
    } else {
      ProcessElement(name.as<string>(), node);
    }
    G4VisAttributes attr;
    ProcessVisAttributes(node, attr);
    if(type == "material") { fMaterialVisAttributes.insert_or_assign(name.as<string>(), attr); }
  }
}
