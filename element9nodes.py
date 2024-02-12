import matplotlib.pyplot as plt
import numpy as np
import sympy as sym
import dill 
dill.settings['recurse'] = True

############################### ELEMENT GEOMETRY ######################################

EDGE_SIZE = 0.1/5   # to overwrite
ELEMENT_CONSTRUCTION = np.array([[0, 0], [0.5, 0], [1.0, 0], [1.0, 0.5], [1.0, 1.0], [0.5, 1.0], [0.0, 1.0], [0.0, 0.5], [0.5, 0.5]])
ELEMENT_CONSTRUCTION_LOCAL = (ELEMENT_CONSTRUCTION - 0.5) * EDGE_SIZE

######################## SYMBOLIC CALCULATIONS #########################################

s, t = sym.Symbol('s'), sym.Symbol('t')                     # local coordinates
a, b = sym.Symbol('a'), sym.Symbol('b')                     # element edges sizes
kx,ky = sym.Symbol('k_x'), sym.Symbol('k_y')                # conductivity coefficients
p,q = sym.Symbol('p'),sym.Symbol('q')                       # Neuman's boundary coefficients
alpha,beta = sym.Symbol('\\alpha'), sym.Symbol('\\beta')    # Neuman's boundary coefficients for edge

# each node of element
CONSTRUCTION_SYMBOLIC = sym.Matrix([[-a, -b], [0, -b], [a, -b], [a, 0], [a, b], [0, b], [-a, b], [-a, 0], [0, 0]])

# interpolation function
interp_components = sym.Matrix([[1, s, t, s * t, s ** 2, t ** 2, s ** 2 * t, t ** 2 * s, s ** 2 * t ** 2]])

# interpolation
Vandermonde = sym.Matrix(np.ones(shape=(CONSTRUCTION_SYMBOLIC.shape[0], (interp_components.shape[1]))))
for r in range(CONSTRUCTION_SYMBOLIC.shape[0]):
    for c in range(len(interp_components)):
        Vandermonde[r, c] = interp_components[0, c].subs(s, CONSTRUCTION_SYMBOLIC[r, 0]).subs(t, CONSTRUCTION_SYMBOLIC[r, 1])
function_sym_N = interp_components @ Vandermonde.inv()   # [N1 N2 N3 N4 N5 N6 N7 N8 N9] - shape function

# helpful matrices
function_sym_B = sym.Matrix([sym.diff(function_sym_N,s), sym.diff(function_sym_N,t)])
matrix_C = np.array([[kx,0],[0,ky]])

# k_k       - conductivity
# k_p, r_q  - Neuman boundaries for element inside
# k_a, r_b  - Neuman boundaries for element's edge

# symbolic matrices
mat_sym_k_k = sym.integrate(  sym.integrate(  function_sym_B.T@matrix_C@function_sym_B , (t,-b,b) ) , (s,-a,a) )

mat_sym_k_a_inside = -1*alpha*function_sym_N.T@function_sym_N
mat_sym_k_a =              [ sym.integrate( (mat_sym_k_a_inside).subs(t,-b) , (s,-a, a) ),
                                sym.integrate( (mat_sym_k_a_inside).subs(s, a) , (t,-b, b) ),
                                sym.integrate( (mat_sym_k_a_inside).subs(t, b) , (s,-a, a) ),
                                sym.integrate( (mat_sym_k_a_inside).subs(s,-a) , (t,-b, b) ) ]  # for each side

mat_sym_r_b_inside = beta*function_sym_N.T
mat_sym_r_b =              [ sym.integrate( (mat_sym_r_b_inside).subs(t,-b) , (s,-a, a) ),
                                sym.integrate( (mat_sym_r_b_inside).subs(s, a) , (t,-b, b) ),
                                sym.integrate( (mat_sym_r_b_inside).subs(t, b) , (s,-a, a) ),
                                sym.integrate( (mat_sym_r_b_inside).subs(s,-a) , (t,-b, b) ) ]  # for each side

mat_sym_k_k=mat_sym_k_k.subs({a: EDGE_SIZE / 2, b: EDGE_SIZE / 2})
# mat_sym_k_p=mat_sym_k_p.subs({a:BOK/2, b:BOK/2})
# mat_sym_r_q=mat_sym_r_q.subs({a:BOK/2, b:BOK/2})

mat_sym_k_a= [item.subs({a: EDGE_SIZE / 2, b: EDGE_SIZE / 2}) for item in mat_sym_k_a]
mat_sym_r_b= [item.subs({a: EDGE_SIZE / 2, b: EDGE_SIZE / 2}) for item in mat_sym_r_b]

####################### NUMERIC FUNCTIONS #############################################

mat_fun_interp_comp = sym.lambdify([s, t], interp_components)

mat_fun_k_k = sym.lambdify([kx,ky],mat_sym_k_k)
# mat_fun_k_p = sym.lambdify([p],mat_sym_k_p)
# mat_fun_r_q = sym.lambdify([q],mat_sym_r_q)

mat_fun_k_a = [sym.lambdify([alpha],item) for item in mat_sym_k_a]
mat_fun_r_b = [sym.lambdify([beta],item) for item in mat_sym_r_b]


################### GLOBAL MATRICES INITIALIZATION ####################################

NODES_COORDS = np.array([0])
global_MAT_k_k = np.array([0])
global_MAT_k_a = np.array([0])
global_MAT_r_b = np.array([0])
u = np.array([0])               # solution
#################### CLASS ELEMENT ###################################

class Element:
    def __init__(self,id,element_nodes_id):
        self.id = id                                                    # element id
        self.nodes = element_nodes_id                                   # inside nodes id in order
        self.directions = [self.nodes[item] for item in [1,3,5,7]]      # down, right, up, left wall
    def agr_conduct(self,k):
        """heat conductivity agregation"""
        idx = np.ix_(self.nodes,self.nodes)
        global_MAT_k_k[idx] += mat_fun_k_k(k[0],k[1])
    def agr_heat_flow(self,nr_bok, q0):
        """heatflow agregation"""
        #print('q0:',q0)
        idx = np.ix_(self.nodes)
        global_MAT_r_b[idx] += mat_fun_r_b[nr_bok](-1*q0/1)
    def agr_convection(self,nr_bok, h, Too):    # nr_bok = element.directions[dir]
        """convection agregation"""
        #print('h:',h,' Too:',Too)
        idx = np.ix_(self.nodes,self.nodes)     # nodes - elementy węzła
        global_MAT_k_a[idx] += mat_fun_k_a[nr_bok](-h)
        idx = np.ix_(self.nodes)
        global_MAT_r_b[idx] += mat_fun_r_b[nr_bok](+1*h*(Too+273.15))
    def agr_T(self,nr_bok, T):
        """temperature agregation"""
        #print('T:',T)
        self.agr_convection(nr_bok, 100000000, T)
    def plot_by_color(self,minT,maxT):
        # val = mat_fun_interp_comp(EDGE_SIZE / 2, EDGE_SIZE / 2) @ u[self.nodes]
        for ntp in [[0,1,8,7,0],[1,2,3,8,1],[7,8,5,6,7],[3,4,5,8,3]]:
            val = np.average(u[self.nodes[ntp]])
            val = 0 + (val-minT)/(maxT-minT)*(1-0)
            plt.fill(NODES_COORDS[self.nodes[ntp],0],NODES_COORDS[self.nodes[ntp],1],color=cmap(val))
        return val
