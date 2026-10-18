class Reader:
    def __init__(self, filename):
        self.filename = filename

        self._material = {}
        self._nodes = {}
        self._elements = []
        self._degrees_of_freedom = {}
        self._forces = []
        self._boundary_conditions = []
        
        self._parse_file()
        
    def _parse_file(self):
        section = None
        
        with open(self.filename, "r", encoding="utf-8") as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                if line.startswith("*"):
                    section = line.removeprefix("*").strip().upper().replace(" ", "_")
                    continue
            
                data = [item.strip() for item in line.split(";")]
                
                if section == "MATERIAL":
                    self._material["E"] = float(data[0])
                    if len(data) > 1 and data[1] != "":
                        self._material["A"] = float(data[1])
                    
                elif section == "NODES":
                    node_id = int(data[0])
                    x, y = float(data[1]), float(data[2])
                    self._nodes[node_id] = [x, y]
                    
                elif section == "ELEMENTS":
                    element_id = int(data[0])
                    node_init, node_final = int(data[1]), int(data[2])
                    elem_E = float(data[3]) if len(data) > 3 and data[3] != "" else None
                    elem_A = float(data[4]) if len(data) > 4 and data[4] != "" else None
                    self._elements.append((element_id, node_init, node_final, elem_E, elem_A))
                    
                elif section == "FORCES":
                    node_id = int(data[0])
                    # Format: node; direction ('x' or 'y'); value
                    if data[1].lower() in ("x", "y"):
                        direction = data[1].lower()
                        force_val = float(data[2])
                        self._forces.append((node_id, direction, force_val))
                    # Format: node; Fx; Fy
                    else:
                        fx = float(data[1])
                        fy = float(data[2])
                        self._forces.append((node_id, "x", fx))
                        self._forces.append((node_id, "y", fy))
                    
                elif section == "BOUNDARY_CONDITIONS":
                    node_id = int(data[0])
                    # Format: node; direction ('x' or 'y'); value
                    if data[1].lower() in ("x", "y"):
                        direction = data[1].lower()
                        displacement = float(data[2])
                        self._boundary_conditions.append((node_id, direction, displacement))
                    # Format: node; ux; uy
                    elif len(data) >= 3:
                        ux = float(data[1])
                        uy = float(data[2])
                        self._boundary_conditions.append((node_id, "x", ux))
                        self._boundary_conditions.append((node_id, "y", uy))
                        
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

    def load_model(self, model):
        for node_id, (x, y) in self.get_nodes().items():
            model.add_node(node_id, x, y)

        mat = self.get_material()
        default_E = mat.get("E")
        default_A = mat.get("A")

        for elem in self.get_elements():
            elem_id, node_i, node_j = elem[0], elem[1], elem[2]
            elem_E = elem[3] if len(elem) > 3 and elem[3] is not None else default_E
            elem_A = elem[4] if len(elem) > 4 and elem[4] is not None else default_A
            model.add_element(elem_id, node_i, node_j, elem_E, elem_A)

        for node_id, direction, displacement in self.get_boundary_conditions():
            model.add_boundary_conditions(node_id, direction, displacement)

        for node_id, direction, force in self.get_forces():
            model.add_forces(node_id, direction, force)

        return model


if __name__ == "__main__":
    reader = Reader("data.txt")
    print("Material:", reader.get_material())
    print("Nós:", reader.get_nodes())
    print("Elementos:", reader.get_elements())
    print("Graus de Liberdade:", reader.get_degrees_of_freedom())
    print("Forças:", reader.get_forces())
    print("Condições de Fronteira:", reader.get_boundary_conditions())