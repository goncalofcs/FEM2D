import numpy as np

class FEM2D:
    def __init__ (self):
        self.nodes = {}
        self.elements_list = []
        self. forces = {}
        self.boundary_conditions = {}

    def add_node(self, node, x, y):
        self.nodes[node] = [x, y]

    def add_element(self, element_id, node_i, node_j, E, A):
        element = Element2D(element_id, node_i, node_j, E, A)
        self.elements_list.append(element)

    def get_dof(self, node, direction):
        if direction == 'x':
            return 2 * node
        else:
            return 2 * node + 1

    def add_forces(self, node, direction, force):
        dof = self.get_dof(node, direction)
        self.forces[dof] = force

    def add_boundary_conditions(self, node, direction, u):
        dof = self.get_dof(node, direction)
        self.boundary_conditions[dof] = u

    def global_stiffness_matrix(self):
        num_nodes = len(self.nodes)
        num_dofs = num_nodes * 2 
        k_global = np.zeros((num_dofs, num_dofs))

        for element in self.elements_list:
            dof1 = self.get_dof(element.node_i,'x')
            dof2 = self.get_dof(element.node_j,'y')
            dof3 = self.get_dof(element.node_i,'x')
            dof4 = self.get_dof(element.node_j,'y')

            dofs_list = [dof1, dof2, dof3, dof4]

            k_e = element.stiffness_matrix()

            for i in range(4):
                for j in range(4):
                    k_global[dofs_list[i], dofs_list[j]] += k_e[i,j]

        return k_global



