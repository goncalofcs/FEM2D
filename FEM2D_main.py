import numpy as np
import matplotlib.pyplot as plt  
from reader import Reader
from visualize import Visualize

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
            
            # Solve K_free * u_free = F_free - K_fixed * u_fixed
            u_free = np.linalg.solve(K_free, F[free_index] - K_fixed @ u_fixed)
            
            u = np.zeros(num_dofs)
            u[fixed_index] = u_fixed
            u[free_index] = u_free
            
            R = K @ u - F
            return u, R


model = FEM2D()
reader = Reader(filename="data.txt")

model = reader.load_model(model)


#print(model.global_stiffness_matrix())

u,R = model.solve()

tensao = []
for i in range(len(model.elements_list)):
    sigma = model.elements_list[i].calculate_stress(u,model.nodes)
    tensao.append(float(sigma))
    
nodes_def = {}
for n in range(len(model.nodes)):
    xi, yi = model.nodes[n]
    xf = xi + u[2*n]
    yf = yi + u[2*n+1]
    nodes_def[n] = xf, yf
    
# results presentation
print("RESULTADOS")
print("")
print("Coordenadas da configuração original")
for n in range(len(model.nodes)):
    x = model.nodes[n][0]
    y = model.nodes[n][1]
    print(f"Nó {n + 1}: ({x},{y}) m")
print("")
print("----//----")
print("")
print("Coordenadas da configuração deformada")
for n in range(len(model.nodes)):
    x,y = np.round(nodes_def[n],6)
    print(f"Nó {n+1}: ({x},{y}) m")
print("")
print("----//----")
print("")
print("Deslocamentos nodais")
for n in range(len(model.nodes)):
    disp_x = np.round(u[2*n] * 1e3, 3)
    disp_y = np.round(u[2*n+1] * 1e3, 3)
    print(f"Nó {n+1}: ux = {disp_x} mm, uy = {disp_y} mm")
print("")
print("----//----")
print("")
print("Forças de reação (R)")
for n in range(len(model.nodes)):
    reaction_x=np.round(R[2*n], 1)
    reaction_y=np.round(R[2*n+1], 1)
    print(f"Nó {n+1}: ({reaction_x}, {reaction_y}) N")
print("")
print("----//----")
print("")
print("Tensão (S11)")
for n in range(len(model.elements_list)):
    stress=np.round(tensao[n]/1e6, 1)
    element = model.elements_list[n]
    print(f"Elemento {n+1} ({element.node_i} --> {element.node_j}): {stress} MPa")
print("")   
print("")
print("============//=============")
print("============//=============")
print("")

# Second case: preserve the original forces and supports, and impose
# an additional vertical displacement at node D (node ID 3).
model_2 = FEM2D()
model_2.nodes = {node_id: coordinates.copy() for node_id, coordinates in model.nodes.items()}
model_2.elements_list = model.elements_list.copy()
model_2.forces = model.forces.copy()
model_2.boundary_conditions = model.boundary_conditions.copy()
model_2.add_boundary_conditions(3, 'y', -0.002)
u_2,R_2 = model_2.solve()

print("Força para o deslocamento imposto uD_y = -2 mm")
print("O nó D é o nó 4")
# in calculus node 3
FD2 = np.round(R_2[2*3+1],2)
print(f"Força em D = {FD2} N")


# GRAPHIC VISUALIZATION
visualizer = Visualize(model, scale=150)

fig, subplots = plt.subplots(1, 2, figsize=(18, 6))

visualizer.plot_configuration(subplots[0], u, "Estrutura Deformada", show_legend=False)
visualizer.plot_configuration(subplots[1], u_2, "Estrutura Deformada com uD_y = -2 mm", show_legend=True)

fig.tight_layout()
plt.show()
