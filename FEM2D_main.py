import numpy as np

np.set_printoptions(precision=8, suppress=True)

class Element2D:
    def __init__(self, element_id, node_i, node_j, E, A):
        self.element_id = element_id
        self.node_i = node_i
        self.node_j = node_j
        self.E = E
        self.A = A
        
        
    def calculate_angle_length(self, nodes):
        xi, yi = nodes[self.node_i]
        xj, yj = nodes[self.node_j]
        dx = xj - xi
        dy = yj - yi
        
        theta = np.arctan2(dy, dx) #rad
        #better than arctan because it solves for dx=0 --> theta = pi/2
        length = np.sqrt(dx**2 + dy**2)
        return theta, length
    
    def calculate_stress(self, u, nodes): # u list
        # node i has degrees of freedom 2*i, 2*i+1
        ui_x = u[2 * self.node_i]
        ui_y = u[2 * self.node_i + 1]
        uj_x = u[2 * self.node_j]
        uj_y = u[2 * self.node_j + 1]
        
        angle, L = self.calculate_angle_length(nodes)
        c = np.cos(angle)
        s = np.sin(angle)
        
        extension = ( (uj_x - ui_x)*c + (uj_y - ui_y)*s ) / L
        
        sigma = self.E * extension
        return sigma
    
    def stiffness_matrix(self, nodes):
        angle, L = self.calculate_angle_length(nodes=nodes)
        c = np.cos(angle)
        s = np.sin(angle)        
        k = self.E * self.A / L
        K_i = k * np.array([[  c**2,   c*s,  -(c**2), -(c*s)], 
                            [  c*s,    s**2, -(c*s), -(s**2)], 
                            [-(c**2), -(c*s), c**2, c*s], 
                            [-(c*s), -(s**2), c*s, s**2]])
        return K_i


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
            dof2 = self.get_dof(element.node_i,'y')
            dof3 = self.get_dof(element.node_j,'x')
            dof4 = self.get_dof(element.node_j,'y')

            dofs_list = [dof1, dof2, dof3, dof4]

            k_i = element.stiffness_matrix(self.nodes)

            for i in range(4):
                for j in range(4):
                    k_global[dofs_list[i], dofs_list[j]] += k_i[i,j]

        return k_global

    def solve(self):
            num_nodes = len(self.nodes)
            num_dofs = 2 * num_nodes
            K = self.global_stiffness_matrix()
            
            F = np.zeros(num_dofs)
            for dof, force in self.forces.items():
                F[dof] = force
                            
            fixed_index = list(self.boundary_conditions.keys())
            free_index = [index for index in range(num_dofs) if index not in fixed_index]
            u_fixed = np.array(list(self.boundary_conditions.values()))
            
            K_free = K[np.ix_(free_index, free_index)]
            K_fixed = K[np.ix_(free_index, fixed_index)]
            
            # Resolver K_free * u_free = F_free - K_fixed * u_fixed
            u_free = np.linalg.solve(K_free, F[free_index] - K_fixed @ u_fixed)
            
            u = np.zeros(num_dofs)
            u[fixed_index] = u_fixed
            u[free_index] = u_free
            
            R = K @ u - F
            return u, R

model = FEM2D()

model.add_node(0, 1.35, 1.1) 
model.add_node(1, 2.25, 1.1) 
model.add_node(2, 2.7, 0)
model.add_node(3, 1.8, 0)
model.add_node(4, 0.45, 1.1)
model.add_node(5, 0.9, 0)
model.add_node(6, 0, 0)

model.add_element(0, 0, 1, 70e9, 560e-6)
model.add_element(1, 1, 2, 70e9, 560e-6)
model.add_element(2, 3, 2, 70e9, 560e-6)
model.add_element(3, 3, 1, 70e9, 560e-6)
model.add_element(4, 4, 0, 70e9, 560e-6)
model.add_element(5, 4, 5, 70e9, 560e-6)
model.add_element(6, 5, 3, 70e9, 560e-6)
model.add_element(7, 0, 3, 70e9, 560e-6)
model.add_element(8, 5, 0, 70e9, 560e-6)
model.add_element(9, 6, 5, 70e9, 560e-6)
model.add_element(10, 6, 4, 70e9, 560e-6)

model.add_forces(0, 'y', -3000)
model.add_forces(1, 'x', 3000)
model.add_forces(4, 'y', -5000)
model.add_forces(5, 'y', -5000)

model.add_boundary_conditions(6, 'x', 0)
model.add_boundary_conditions(6, 'y', 0)
model.add_boundary_conditions(2, 'y', 0)

#print(model.global_stiffness_matrix())

u,R = model.solve()
print(u)

'''
print(R)

tensao = []
for i in range(3):
    sigma = model.elements_list[i].calculate_stress(u,model.nodes)
    tensao.append(float(sigma))
print(tensao)

#elements_sig_rounded = np.round(elements_sig / 1e6, 2)
'''