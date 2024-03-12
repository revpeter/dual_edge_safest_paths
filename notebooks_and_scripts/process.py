import networkx as nx
import numpy as np
import copy

def get_srlg_cfp(srlgs, srlg_probs, scf):
    if scf in srlgs:
        cfp = srlg_probs[srlgs.index(scf)]
    else:
        cfp = 0
    return cfp


def get_st_graph(g_base, srlgs, srlg_probs, g_dual, pos_dual, s, t):
    g_dual_st = copy.deepcopy(g_dual)
    pos_dual_st = copy.deepcopy(pos_dual)

    nx.set_node_attributes(g_dual_st, "red", "c")
    nx.set_edge_attributes(g_dual_st, "red", "c")

    edges_for_s = [i for i in g_dual_st.nodes() if s in i[:2]]
    edges_for_t = [i for i in g_dual_st.nodes() if t in i[:2]]

    g_dual_st.add_node("Start", c = "blue")
    for pe in edges_for_s:
        g_dual_st.add_edge("Start", pe, c = "blue")

    g_dual_st.add_node("Target", c = "green")
    for pe in edges_for_t:
        g_dual_st.add_edge("Target", pe, c = "green")


    pos_dual_st.update({"Start":(g_base.nodes[s]["Longitude"], g_base.nodes[s]["Latitude"])})
    pos_dual_st.update({"Target":(g_base.nodes[t]["Longitude"], g_base.nodes[t]["Latitude"])})

    ### Calculate cost of the dual edges ###
    edge_cost_dict = {}
    for f, g,_ in g_dual_st.edges(keys=True):
        if f in ["Start", "Target"]:
            e_cfp = get_srlg_cfp(srlgs, srlg_probs, {g})
            cost = 0.5*e_cfp
        elif g in ["Start", "Target"]:
            e_cfp = get_srlg_cfp(srlgs, srlg_probs, {f})
            cost = 0.5*e_cfp
        else:
            f_cfp = get_srlg_cfp(srlgs, srlg_probs, {f})
            g_cfp = get_srlg_cfp(srlgs, srlg_probs, {g})
            fg_cfp = get_srlg_cfp(srlgs, srlg_probs, {f,g})
            
            if (f_cfp == 0) or (g_cfp == 0):
                cost = max([f_cfp, g_cfp])
            else:
                cost = (0.5*f_cfp) + (0.5*g_cfp) - fg_cfp

        edge_cost_dict.update({(f,g,0):cost})

    nx.set_edge_attributes(g_dual_st, edge_cost_dict, "cost")

    return(g_dual_st, pos_dual_st)


def get_fuggetlen_graph(g_base, srlgs, srlg_probs):
    g_fuggetlen = copy.deepcopy(g_base)

    edge_cost_dict = {}
    for edge in g_base.edges(keys=True):
        p_e = srlg_probs[srlgs.index({edge})]
        edge_cost_dict.update({edge:-np.log(1-p_e)})
    
    nx.set_edge_attributes(g_fuggetlen, edge_cost_dict, "cost")

    return g_fuggetlen