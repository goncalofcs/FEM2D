import numpy as np

class FEM2D:
    def __init__ (self):
        self.nodes = {}
        self.elements = []
        self. forces = {}
        self.boundary_conditions = {}

    def add_node(self, node, x, y):
        self.nodes[node] = [x, y]

    def add_element(self, element_id, node_i, node_j, E, A, L = None):
        element = Element2D(element_id, node_i, node_j, E, A, L)
        element.get_angle_length(self.nodes)
        self.elements.append(element)

    def get_dof(self, node, direction):
        if direction == 'x':
            return 2 * (node - 1)
        else:
            return 2 * (node - 1) + 1


