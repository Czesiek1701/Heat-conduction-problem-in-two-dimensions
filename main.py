import numpy as np
from matplotlib import pyplot as plt
import element9nodes as e9n
from matplotlib import cm
import matplotlib as mpl

############################### SET UP'Y ######################################

np.set_printoptions(linewidth=320)

FILE_INPUT_GEO = "input/geometria_pipe.txt"
FILE_INPUT_BC = "input/warunki_brzegowe_pipe.txt"
FILE_ELEMENTS = "robocze/elementy.txt"
FILE_NODES = "robocze/wezly.txt"

GENERATE_NEW_GRID = True
PLOT_GRID = False
SHOW_RESULTS_COLOR = True

VSV = 1e-8

val_kx,val_ky = 55,55

############################### FUNCTIONS ###########################################

def create_initial_points(file_from, file_to="robocze/wp.txt"):
    """infill points in rectangles from file"""
    m_rectangles = np.loadtxt(file_from)
    with open(file_to,'w') as wfile:
        initial_grid=np.array([[m_rectangles[0, 0], m_rectangles[0, 1]]]).astype('float')
        for prostokat in m_rectangles:
            xs = prostokat[0]
            ys = prostokat[1]
            xe = prostokat[2]
            ye = prostokat[3]
            temp_grid = np.round(np.mgrid[xs:(xe + VSV):e9n.EDGE_SIZE, ys:(ye + VSV):e9n.EDGE_SIZE].T.reshape(-1, 2), 6)
            initial_grid = np.concatenate([initial_grid, temp_grid], axis=0)
            initial_grid = np.unique(initial_grid, axis=0)
    initial_grid = np.unique(initial_grid, axis=0)
    np.savetxt(file_to, initial_grid, fmt='%.6f')

def create_initial_elements(file_punkty="robocze/wp.txt", file_elementy="robocze/er.txt"):
    """stack points to initial elements"""
    m_points = np.round(np.loadtxt(file_punkty), 6)
    l_points = m_points.tolist()
    ELEMENTS = []
    el_id=0
    for p_id, row in enumerate(l_points):
        all_in=True
        for addit in [[e9n.EDGE_SIZE, 0], [e9n.EDGE_SIZE, e9n.EDGE_SIZE], [0, e9n.EDGE_SIZE]]:
            if [round(row[0]+addit[0],6),round(row[1]+addit[1],6)] not in l_points:
                all_in=False
        if all_in:
            ELEMENTS.append([p_id])
            for addit in [[e9n.EDGE_SIZE, 0], [e9n.EDGE_SIZE, e9n.EDGE_SIZE], [0, e9n.EDGE_SIZE]]:
                ELEMENTS[-1].append(l_points.index([round(row[0] + addit[0], 6), round(row[1] + addit[1], 6)]))
            el_id+=1
    np.savetxt(file_elementy,np.array(ELEMENTS).astype('int'),fmt='%i')

def create_nodes(file_output_wezly, file_punkty="robocze/wp.txt", file_elementy="robocze/er.txt"):
    """create global nodes based on local in elements"""
    m_points = np.loadtxt(file_punkty)
    m_nodes = np.array([[0, 0]])
    m_init_elements = np.loadtxt(file_elementy).astype('int')
    for element in m_init_elements:
        for i in range(element.shape[0]):
            p_start = (m_points[element[i]])
            p_end = (m_points[element[(i + 1) % 4]])
            p_between = (p_start + p_end)/2
            m_nodes = np.concatenate([m_nodes, np.array([p_between.tolist()])], axis=0)
        p_center = (m_points[element[0]] + m_points[element[2]]) / 2
        m_nodes = np.concatenate([m_nodes, np.array([p_center.tolist()])], axis=0)
    m_nodes = np.unique(m_nodes[1:], axis=0)
    np.savetxt(file_output_wezly, np.concatenate([m_points, m_nodes], axis=0), fmt='%.8f')

def group_elements(file_wezly, output_file_elementy, file_elementy_init="robocze/er.txt"):
    """group nodes to target elements"""
    m_nodes = np.round(np.loadtxt(file_wezly),6)
    l_nodes = m_nodes.tolist()
    for node in l_nodes:
        print(node)
    m_init_elements = np.loadtxt(file_elementy_init).astype('int')
    ELEMENTS = []
    for init_element in m_init_elements:
        start_point = m_nodes[init_element[0]]
        ELEMENTS.append([])
        for addit in e9n.ELEMENT_CONSTRUCTION * e9n.EDGE_SIZE:
            act_p = np.round((start_point+addit),6)
            point_id = l_nodes.index([act_p[0], act_p[1]])
            ELEMENTS[-1].append(point_id)
    np.savetxt(output_file_elementy,ELEMENTS,fmt='%i')

def nodes_in_rect(candidates, p1, p2):
    """Return list of nodes iniside rectangle"""
    pxmin = min(p1[0], p2[0])
    pxmax = max(p1[0], p2[0])
    pymin = min(p1[1], p2[1])
    pymax = max(p1[1], p2[1])
    list_in=[]
    for id,candi in enumerate(candidates):
        if pxmin <= candi[0] <= pxmax and pymin <= candi[1] <= pymax:
            list_in.append(id)
    return list_in

def show_results(ELEMENTS, temperature):

    fig = plt.figure(1)
    plt.axis('equal')
    plt.title("Temperature $^oC$")
    plt.grid()

    # colors = plt.cm.jet(np.linspace(np.min(temperature),np.max(temperature),256))
    e9n.cmap = cm.get_cmap('RdYlBu_r')
    minT = np.min(temperature)
    maxT = np.max(temperature)

    for element in ELEMENTS:
        element.plot_by_color(minT,maxT)
    v = np.linspace(minT, maxT, 8, endpoint=True)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.05])
    cb = mpl.colorbar.ColorbarBase(ax, orientation='horizontal',
                                   cmap = e9n.cmap,
                                   norm = mpl.colors.Normalize(vmin=minT, vmax=maxT),
                                   ticks = v,
                                   format = mpl.ticker.FormatStrFormatter('%.1f'))
    plt.show()

################################ MAIN ################################################

def main():

    # generating new grid
    if GENERATE_NEW_GRID:
        create_initial_points(FILE_INPUT_GEO)
        create_initial_elements()
        create_nodes(FILE_NODES)
        group_elements(FILE_NODES, FILE_ELEMENTS)

    # read nodes coords and elements (stack nodes)
    ELEMENT_ID_NODES = np.loadtxt("robocze/elementy.txt").astype(int)
    e9n.NODES_COORDS = np.loadtxt("robocze/wezly.txt")

    # plot grid
    if True:
        fig = plt.figure(1)
        if PLOT_GRID:
            plt.plot(e9n.NODES_COORDS[:, 0], e9n.NODES_COORDS[:, 1], 'r.')
        for element in ELEMENT_ID_NODES:
            # print(element)
            # plt.plot(WEZLY[element[8],0],WEZLY[element[8],1],'*')
            x,y=[],[]
            for p in element[:-1]:
                x.append(e9n.NODES_COORDS[p, 0])
                y.append(e9n.NODES_COORDS[p, 1])
            x.append(e9n.NODES_COORDS[element[0], 0])
            y.append(e9n.NODES_COORDS[element[0], 1])
            plt.plot(x,y,'c-',linewidth=.1)
        if PLOT_GRID:
            for id,node in enumerate(e9n.NODES_COORDS):
                plt.text(node[0] - e9n.EDGE_SIZE / 8, node[1] - e9n.EDGE_SIZE / 8, str(id))
        plt.axis('equal')
        if PLOT_GRID:
            plt.show()

    # initialize objects class Element
    ELEMENTS = []
    for id, elid in enumerate(ELEMENT_ID_NODES):
        ELEMENTS.append(e9n.Element(id, elid))

    # empty global matrices
    e9n.global_MAT_k_k = np.zeros(shape=(e9n.NODES_COORDS.shape[0], e9n.NODES_COORDS.shape[0]))
    e9n.global_MAT_k_a = np.zeros(shape=(e9n.NODES_COORDS.shape[0], e9n.NODES_COORDS.shape[0]))
    e9n.global_MAT_r_b = np.zeros(shape=(e9n.NODES_COORDS.shape[0], 1))

    # agregation of conductivity
    for element in ELEMENTS:
        element.agr_conduct([val_kx,val_ky])

    # boundary conditions and agregation
    with open(FILE_INPUT_BC, 'r') as rfile:
        for line in rfile:
            sline = line.split(' ')
            p1 = [float(item) for item in sline[0:2]]   # first point of rect
            p2 = [float(item) for item in sline[2:4]]   # second point of rect
            plt.plot([p1[0],p2[0]],[p1[1],p2[1]],'r-')  # plot where boundary is set
            nodes_inside = nodes_in_rect(e9n.NODES_COORDS, p1, p2)  # get nodes inside rect
            for element in ELEMENTS:
                for dir in [0,1,2,3]:   # iterate on directions (walls of element)
                    if element.directions[dir] in nodes_inside:
                        # print(element.id, element.directions[dir])
                        if sline[4] == 'hf':    # heatflow
                            element.agr_heat_flow(dir,float(sline[5]))
                        if sline[4] == 'cnv':   # convection
                            element.agr_convection(dir,float(sline[5]),float(sline[6]))
                        if sline[4] == 'T':     # temperature
                            element.agr_T(dir,float(sline[5]))

    # calculate distribution of temperature
    e9n.u = np.linalg.inv(e9n.global_MAT_k_k + e9n.global_MAT_k_a) @ e9n.global_MAT_r_b - 273.15

    # print temperature of nodes on plot
    if PLOT_GRID:
        for i in range(0, e9n.NODES_COORDS.shape[0]):
            plt.text(e9n.NODES_COORDS[i, 0], e9n.NODES_COORDS[i, 1],
                     '{:.1f}'.format(e9n.u[i, 0]), fontsize = 8 )

    # plt.axis('equal')
    # plt.title('Temperature $^oC$')

    if SHOW_RESULTS_COLOR:
        show_results(ELEMENTS , e9n.u)
    else:
        plt.show()

    return(np.min(e9n.u))


############################ MAIN FUNCTION CALL ################################################

if __name__ == "__main__":
    main()
