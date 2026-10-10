
class Reader:
    def __init__(self, filename):
        self.filename = filename
        
        self._material = {}
        self._nodes = {}
        self._elements = []
        self._degrees_of_freedom = {}
        self._forces = {}
        self._boundary_conditions = {}
        
        self._parse_file()
        
    def _parse_file(self):
        section = None
        
        with open(self.filename, "r", encoding="utf-8") as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                if line.startswith("*"):
                    section = line.removeprefix("*")
                    print(section)
                    continue
            
                data = [item.strip() for item in line.split(";")]
                
                if section == "MATERIAL":
                    self._material["E"] = float(data[0])
                    self._material["A"] = float(data[1])
                    
                elif section == "NODES":
                    node_id = int(data[0]);
                    x, y = float(data[1]), float(data[2])
                    self._nodes[node_id] = [x, y]
                    
                elif section == "ELEMENTS":
                    element_id = int(data[0])
                    node_init, node_final = int(data[1]), int(data[2])
                    self._elements.append((element_id, node_init, node_final))
                
                elif section == "DEGREES OF FREEDOM":
                    node_id = int(data[0])
                    fix_x = int(data[1])
                    fix_y = int(data[2])
                    self._degrees_of_freedom[node_id] = (fix_x, fix_y)
                    
                elif section == "FORCES":
                    node_id = int(data[0])
                    Fx, Fy = float(data[1]), float(data[2])
                    self._forces[(node_id, "x")] = Fx
                    self._forces[(node_id, "y")] = Fy
                    
                elif section == "BOUNDARY CONDITIONS":
                    node_id = int(data[0])
                    ux = float(data[1])
                    uy = float(data[2])
                    self._boundary_conditions[node_id] = (ux, uy)
                        
    def get_material(self):
        return self._material
    
    def get_nodes(self):
        return self._nodes
    
    def get_elements(self):
        return self._elements
    
    def get_degrees_of_freedom(self):
        return self._degrees_of_freedom
    
    def get_forces(self):
        return self._forces
    
    def get_boundary_conditions(self):
        return self._boundary_conditions

    
reader = Reader("data.txt")
print(reader.get_material())
print(reader.get_nodes())
print(reader.get_elements())
print(reader.get_degrees_of_freedom())
print(reader.get_forces())
print(reader.get_boundary_conditions())