from backend import *
import pandas as pd
import numpy as np
import networkx as nx
import copy
import matplotlib.pyplot as plt
import itertools


class Pipeline:

    def __init__(self, network, filter_th = False, th_mod = 1, equal_reward = False) -> None:
        self.equal_reward = equal_reward
        
        ### Init network srlgs ###
        self.srlgs, srlg_probs = get_SRLGs(f'../PSRLGs/{network}.xml')
        self.srlg_probs = np.asarray(srlg_probs)
        #self.srlg_lengths = np.asarray([len(srlg) for srlg in self.srlgs])

        ### Init the original graph ####
        self.g_base = nx.read_gml(f'../networks/{network}.gml', label="id")
        self.g_base_pos = {node:(self.g_base.nodes[node]['Longitude'], self.g_base.nodes[node]['Latitude']) for node in self.g_base.nodes}
        if network == "italy_995":
            self.g_base.remove_edges_from([(19, 21, 1)])

        ### Calculate threshold and modify srlg list if needed ###
        if filter_th:
            srlg_index_list_above_thcut = self.calculate_edges_above_th(th_mod)
            self.srlgs = [self.srlgs[i] for i in srlg_index_list_above_thcut]
            self.srlg_probs = [self.srlg_probs[i] for i in srlg_index_list_above_thcut]

        ### Create the edge dual graph ###
        self.g_dual = nx.line_graph(self.g_base)
        self.g_dual_pos = {node:(
                            (self.g_base.nodes[node[0]]["Longitude"] + self.g_base.nodes[node[1]]["Longitude"]) /2, 
                            (self.g_base.nodes[node[0]]["Latitude"] + self.g_base.nodes[node[1]]["Latitude"]) /2
                        )
                    for node in self.g_dual.nodes()
                    }

        # self.g_thcut = copy.deepcopy(self.g_base)
        # for i in self.g_base.edges(keys = True):
        #     if i not in self.edges_above_th:
        #         self.g_thcut.remove_edge(i[0], i[1])
        # self.g_thcut_dual = nx.line_graph(self.g_thcut)

    def calculate_edges_above_th(self, th_mod):
        srlg_desc = {"idx":[], "edgeNum":[], "prob":[], "connected":[], "diameter":[]}
        for i,val in enumerate(self.srlgs):
            srlg_desc["idx"].append(i)
            srlg_desc["edgeNum"].append(len(val))
            srlg_desc["prob"].append(self.srlg_probs[i])

            subG = nx.Graph([s[:2] for s in self.srlgs[i]])
            connected = nx.is_connected(subG)

            srlg_desc["connected"].append(connected)
            if connected:
                srlg_desc["diameter"].append(nx.diameter(subG))
            else:
                cc = nx.connected_components(subG)
                S = [subG.subgraph(c).copy() for c in cc]
                biggest_cc_dia = max([nx.diameter(i) for i in S])
                srlg_desc["diameter"].append(1000)

        self.srlgDf = pd.DataFrame(srlg_desc)
        self.highest_bad = self.srlgDf.query("diameter > 2 | connected == False").sort_values("prob", ascending=False).iloc[0]["prob"] * th_mod
        
        #self.highest_bad_above_srlg = sum(self.srlgDf.query("diameter <= 2 & connected == True")["prob"] > self.highest_bad) / len(self.srlgDf.query("diameter <= 2 & connected == True")["prob"])

        srglIndexList = self.srlgDf.query(f"diameter <= 2 & connected == True & prob > {self.highest_bad}")["idx"].to_list()
        
        self.edges_above_th = set()
        for i in srglIndexList:
            if len(self.srlgs[i]) <= 2:
                self.edges_above_th = self.edges_above_th.union(self.srlgs[i])
            else:
                continue
        
        return srglIndexList

    def get_srlg_cfp(self, scf):
        if scf in self.srlgs:
            cfp = self.srlg_probs[self.srlgs.index(scf)]
        else:
            cfp = 0
        return cfp

    def match_edge_with_srlgs(self, path_nodes):
        path_edges = []
        one_edge_srlgs = [list(srlg)[0] for srlg in self.srlgs if len(srlg) == 1]
        for ind,val in enumerate(path_nodes[:-1]):
            if (val,path_nodes[ind+1],0) in one_edge_srlgs:
                path_edges.append((val,path_nodes[ind+1],0))
            elif (path_nodes[ind+1],val,0) in one_edge_srlgs:
                path_edges.append((path_nodes[ind+1],val,0))
            else:
                path_edges.append((val,path_nodes[ind+1],0))
    
        return path_edges
    
    def sort_edge_nodes(self, list_of_edges):
        list_of_edges_mod = copy.deepcopy(list_of_edges)
        for i, edge in enumerate(list_of_edges):
            if edge[0] > edge[1]:
                nodes_sorted = sorted(edge[:2])
                nodes_sorted.append(0)
                list_of_edges_mod[i] = tuple(nodes_sorted)
            else:
                continue
        return list_of_edges_mod
        

    def edge_dual_availability(self, edge_list):
        ed_avb = 0
        for edge in edge_list:
            ed_avb += self.get_srlg_cfp({edge})
        
        for ind, edge in enumerate(edge_list[:-1]):
            ed_avb -= self.get_srlg_cfp({edge, edge_list[ind+1]})
        
        return ed_avb
    
    def edge_dual_availability_metaData(self, edge_list):
        eda_metaData = {"pos":[], "neg":[]}
        for edge in edge_list:
            eda_metaData["pos"].append(({edge}, self.get_srlg_cfp({edge})))
        
        for ind, edge in enumerate(edge_list[:-1]):
            eda_metaData["neg"].append( ({edge, edge_list[ind+1]}, self.get_srlg_cfp({edge, edge_list[ind+1]})) )
        
        return eda_metaData

    def create_st_graph(self, s, t):
        g_dual_st = copy.deepcopy(self.g_dual)
        dual_st_pos = copy.deepcopy(self.g_dual_pos)

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


        dual_st_pos.update({"Start":(self.g_base.nodes[s]["Longitude"], self.g_base.nodes[s]["Latitude"])})
        dual_st_pos.update({"Target":(self.g_base.nodes[t]["Longitude"], self.g_base.nodes[t]["Latitude"])})

        ### Calculate cost of the dual edges ###
        edge_cost_dict = {}
        for f, g, _ in g_dual_st.edges(keys=True):
            if f in ["Start", "Target"]:
                e_cfp = self.get_srlg_cfp({g})
                cost = 0.5*e_cfp
            elif g in ["Start", "Target"]:
                e_cfp = self.get_srlg_cfp({f})
                cost = 0.5*e_cfp
            else:
                f_cfp = self.get_srlg_cfp({f})
                g_cfp = self.get_srlg_cfp({g})
                fg_cfp = self.get_srlg_cfp({f,g})
                
                if (f_cfp == 0) | (g_cfp == 0):
                    if self.equal_reward:
                        cost = 0
                    else:
                        cost = max([f_cfp, g_cfp]) * 0.5
                else:
                    cost = (0.5*f_cfp) + (0.5*g_cfp) - fg_cfp

            edge_cost_dict.update({(f,g,0):cost})

        nx.set_edge_attributes(g_dual_st, edge_cost_dict, "cost")

        return(g_dual_st, dual_st_pos)


    def create_fuggetlen_graph(self):
        g_fuggetlen = copy.deepcopy(self.g_base)

        edge_cost_dict = {}
        for edge in self.g_base.edges(keys=True):
            p_e = self.get_srlg_cfp({edge})
            edge_cost_dict.update({edge:-np.log(1-p_e)})
        
        nx.set_edge_attributes(g_fuggetlen, edge_cost_dict, "cost")

        return g_fuggetlen

    def calculate_all_path(self):
        resDict = {"simID":[], "start":[], "target":[],
            "len_sht":[],"avb_sht":[], "ed_avb_sht":[], 
            "len_djk":[], "avb_djk":[], "ed_avb_djk":[], 
            "len_fugg":[], "avb_fugg":[], "ed_avb_fugg":[]}
        resDictPaths = {}
        resDictPathCfps = {}
        resDictPathAvbMeta = {}
        simID = 0

        for s,t in tqdm(list(itertools.combinations(self.g_base.nodes, 2))):
            ### Dual graph with our cost function ###
            g_dual_st, _ = self.create_st_graph(s,t)
            djk_path = nx.dijkstra_path(g_dual_st, "Start", "Target", weight="cost")
            djk_path = self.sort_edge_nodes(djk_path)
            djk_avb = np.prod([1- self.get_srlg_cfp({i}) for i in djk_path[1:-1]])
            djk_ed_avb = self.edge_dual_availability(djk_path[1:-1])

            ### Original graph with Fuggetlen cost function ###
            g_fuggetlen = self.create_fuggetlen_graph()
            djk_path_fuggetlen_nodes = nx.dijkstra_path(g_fuggetlen, s, t, "cost")
            djk_path_fuggetlen_edges = self.match_edge_with_srlgs(djk_path_fuggetlen_nodes)
            djk_path_fuggetlen_edges = self.sort_edge_nodes(djk_path_fuggetlen_edges)
            fuggetlen_avb = np.prod([1 - self.get_srlg_cfp({i}) for i in djk_path_fuggetlen_edges])
            fuggetlen_ed_avb = self.edge_dual_availability(djk_path_fuggetlen_edges)

            ### Original graph with Shortest path ###
            sht_path_nodes = nx.shortest_path(self.g_base, s,t)
            sht_path_edges = self.match_edge_with_srlgs(sht_path_nodes)
            sht_path_edges = self.sort_edge_nodes(sht_path_edges)
            sht_avb = np.prod([1 - self.get_srlg_cfp({i}) for i in sht_path_edges])
            sht_ed_avb = self.edge_dual_availability(sht_path_edges)

            ### Collect results ###
            for val, key in zip([simID, s, t, len(sht_path_edges), sht_avb, sht_ed_avb, len(djk_path)-2, djk_avb, djk_ed_avb, len(djk_path_fuggetlen_edges), fuggetlen_avb, fuggetlen_ed_avb], resDict.keys()):
                resDict[key].append(val)

            resDictPaths.update({simID:{"djk_path":djk_path[1:-1], "sht_path":sht_path_edges, "fuggetlen_path":djk_path_fuggetlen_edges}})
            resDictPathCfps.update({simID:{"djk_path":[self.get_srlg_cfp({i}) for i in djk_path[1:-1]],
                                           "sht_path":[self.get_srlg_cfp({i}) for i in sht_path_edges],
                                           "fuggetlen_path":[self.get_srlg_cfp({i}) for i in djk_path_fuggetlen_edges]
                                           }
                                    })

            resDictPathAvbMeta.update({simID:{"djk_path":self.edge_dual_availability_metaData(djk_path[1:-1]),
                                           "sht_path": self.edge_dual_availability_metaData(sht_path_edges),
                                           "fuggetlen_path":self.edge_dual_availability_metaData(djk_path_fuggetlen_edges)
                                           }
                                    })
            
            simID += 1
        
        

        resDf = pd.DataFrame(resDict)
        resDf["avb_diff"] = resDf["avb_djk"] - resDf["avb_fugg"] # djk > fugg == (+) | djk < fugg == (-)
        resDf["ed_avb_diff"] = resDf["ed_avb_djk"] - resDf["ed_avb_fugg"]
        resDf["ed_avb_pct_diff"] = (resDf["ed_avb_djk"] - resDf["ed_avb_fugg"]) / resDf["ed_avb_fugg"] * 100 # Ha (-) a duális cost-ja ennyi %-al kevesebb mint a fugg cost
        resDf["pct_zero_djk"] = resDf["simID"].apply(lambda x: resDictPathCfps[x]["djk_path"].count(0) / len(resDictPathCfps[x]["djk_path"]))
        resDf["pct_zero_fugg"] = resDf["simID"].apply(lambda x: resDictPathCfps[x]["fuggetlen_path"].count(0) / len(resDictPathCfps[x]["fuggetlen_path"]))
        resDf["pct_zero_sht"] = resDf["simID"].apply(lambda x: resDictPathCfps[x]["sht_path"].count(0) / len(resDictPathCfps[x]["sht_path"]))
        resDf["first_unmatch"] = resDf["simID"].apply(lambda x: self.first_unmatch(x, resDictPaths))
        resDf["len_diff"] = resDf["len_djk"] - resDf["len_fugg"] # Ha (+) a dualis megoldás annyival hosszabb
        resDf["pct_zero_diff"] = resDf["pct_zero_djk"] - resDf["pct_zero_fugg"]
        resDf["ed_avb_class"] = resDf["ed_avb_diff"].apply(lambda x: "dual" if x < 0 else "fugg" if x > 0 else "eq")

        return(resDf, resDictPaths, resDictPathCfps, resDictPathAvbMeta)

    def calculate_one_path(self,s,t):
        self.g_dual_st, self.pos_dual_st = self.create_st_graph(s,t)
        _djk_path = nx.dijkstra_path(self.g_dual_st, "Start", "Target", weight="cost")
        self.djk_path = self.sort_edge_nodes(_djk_path)

        self.g_fuggetlen = self.create_fuggetlen_graph()
        self.djk_path_fuggetlen_nodes = nx.dijkstra_path(self.g_fuggetlen, s, t, "cost")
        _djk_path_fuggetlen_edges = self.match_edge_with_srlgs(self.djk_path_fuggetlen_nodes)
        self.djk_path_fuggetlen_edges = self.sort_edge_nodes(_djk_path_fuggetlen_edges)

        self.sht_path_nodes = nx.shortest_path(self.g_base, s,t)
        _sht_path_edges = self.match_edge_with_srlgs(self.sht_path_nodes)
        self.sht_path_edges = self.sort_edge_nodes(_sht_path_edges)

    def first_unmatch(self,sid, resDictPaths):
        dual_list = resDictPaths[sid]["djk_path"]
        fugg_list = resDictPaths[sid]["fuggetlen_path"]
        if dual_list != fugg_list:
            for pos in range(min([len(dual_list), len(fugg_list)])):
                if (pos+1 == len(dual_list)) or (pos+1 == len(fugg_list)):
                    return pos
                elif dual_list[pos] == fugg_list[pos]:
                    continue
                else:
                    return pos
        else:
            return np.nan

    def get_srlg_info_table(self, miniNetworkEdges):
        srlg_desc = {"idx":[], "edgeNum":[], "prob":[], "connected":[], "diameter":[], "usable":[]}

        for idx, srlg in enumerate(self.srlgs):
            if sum([i in miniNetworkEdges for i in srlg]) == len(srlg):
                srlg_desc["idx"].append(idx)
                srlg_desc["edgeNum"].append(len(srlg))
                srlg_desc["prob"].append(self.srlg_probs[idx])

                subG = nx.Graph([s[:2] for s in self.srlgs[idx]])
                connected = nx.is_connected(subG)

                srlg_desc["connected"].append(connected)
                if connected:
                    srlg_desc["diameter"].append(nx.diameter(subG))
                else:
                    cc = nx.connected_components(subG)
                    S = [subG.subgraph(c).copy() for c in cc]
                    srlg_desc["diameter"].append(1000)

                if (srlg_desc["connected"][-1] == True) & (srlg_desc["diameter"][-1] <= 2):
                    srlg_desc["usable"].append(True)
                else:
                    srlg_desc["usable"].append(False)
        
        srlg_info_df = pd.DataFrame(srlg_desc)

        if hasattr(self, "highest_bad"):
            srlg_info_df["above_th"] = srlg_info_df["prob"].apply(lambda x: True if x > self.highest_bad else False)

        return srlg_info_df

    ################## PLOTS #####################
        
    def plot_base_and_dual(self,fgSize = (10,9)):
        options = {
            'node_size': 300,
            'node_color': 'gray',
            'width': [2+1*self.g_base.edges[edge]['onspine'] for edge in self.g_base.edges],
        }

        fig,ax = plt.subplots(figsize=fgSize)
        nx.draw_networkx(self.g_base, self.g_base_pos, ax, **options)
        nx.draw_networkx(self.g_dual, self.g_dual_pos, with_labels=False, node_size = 50, edge_color = "red", node_color = "red", alpha = 0.5)
        plt.show()

    def plot_dual_with_st_nodes(self,fgSize = (10,9)):
        node_color_list = [j["c"] for i,j in self.g_dual_st.nodes(data=True)]
        edge_color_list = [k["c"] for i,j,k in self.g_dual_st.edges(data=True)]  
        
        plt.figure(figsize=fgSize)
        nx.draw_networkx(self.g_dual_st, self.pos_dual_st, with_labels=False, node_size = 50, edge_color = edge_color_list, node_color = node_color_list)
        plt.show()


    def plot_dual_with_path(self, fgSize = (10,9)):
        node_color_list = [j["c"] for i,j in self.g_dual_st.nodes(data=True)]
        edge_color_list = [k["c"] for i,j,k in self.g_dual_st.edges(data=True)]  
    
        plt.figure(figsize=fgSize)
        nx.draw_networkx(self.g_dual_st, self.pos_dual_st, with_labels=False, node_size = 50, edge_color = edge_color_list, node_color = node_color_list)

        for path_node in self.djk_path:
            plt.scatter(self.pos_dual_st[path_node][0], self.pos_dual_st[path_node][1], facecolors='none', edgecolors='k', s = 200, linewidths=2)

        plt.show()

    def plot_base_with_path(self,s,t,fgSize=(10,9), export_to = False):
        options = {
            'node_size': 200,
            'node_color': 'gray',
            "width": [3 if edge in self.djk_path else 1 for edge in self.g_base.edges],
            "edge_color": ["orange" if edge in self.djk_path else "black" for edge in self.g_base.edges]
        }
        options2 = copy.deepcopy(options)
        options2["edge_color"] = ["purple" if edge in self.djk_path_fuggetlen_edges else "black" for edge in self.g_base.edges]
        options2["width"] = [3 if edge in self.djk_path_fuggetlen_edges else 1 for edge in self.g_base.edges]

        options3 = copy.deepcopy(options)
        options3["edge_color"] = ["pink" if edge in self.sht_path_edges else "black" for edge in self.g_base.edges]
        options3["width"] = [3 if edge in self.sht_path_edges else 1 for edge in self.g_base.edges]

        fig, ax = plt.subplots(ncols = 3, figsize=fgSize)
        plt.subplots_adjust(wspace=0.01)
        nx.draw_networkx(self.g_base, self.g_base_pos, **options, ax = ax[0])
        nx.draw_networkx(self.g_base, self.g_base_pos, **options2, ax = ax[1])
        nx.draw_networkx(self.g_base, self.g_base_pos, **options3, ax = ax[2])

        ax[0].set_title("Edge-dual method")
        ax[1].set_title("Független method")
        ax[2].set_title("Shortest path")
        fig.suptitle(f"Route {s} --> {t}")
        
        if export_to:
            plt.savefig(export_to + f"route_{s}_{t}.jpeg", dpi = 150)
            plt.close()
        else:
            plt.show()

############## Out of Class functions ##################

def print_results(resDf):
    print("Független féle availability alapján")
    print("Független a jobb",sum(resDf["avb_diff"] < 0))
    print("Egyelnő",sum(resDf["avb_diff"] == 0))
    print("Dualis a jobb",sum(resDf["avb_diff"] > 0))

    print("\n")

    print("Duális féle cost függvény alapján")
    fugg_mean_avb = f'{np.mean(resDf[resDf["ed_avb_diff"] > 0]["ed_avb_fugg"]):.2e}'
    fugg_mean_len = round(resDf[resDf["ed_avb_diff"] > 0]["len_fugg"].mean(), 2)
    print(f"Független a jobb {sum(resDf['ed_avb_diff'] > 0)} esetben | Átlag cost: {fugg_mean_avb} | Átlag hossz: {fugg_mean_len}")
    
    equal_mean_avb = round(np.mean(resDf[resDf["ed_avb_diff"] == 0]["ed_avb_djk"]),4)
    equal_djk_mean_len = round(resDf[resDf["ed_avb_diff"] == 0]["len_djk"].mean(), 2)
    equal_fugg_mean_len = round(resDf[resDf["ed_avb_diff"] == 0]["len_fugg"].mean(), 2)     
    print(f"Egyelnő {sum(resDf['ed_avb_diff'] == 0)} esetben | Átlag cost: {equal_mean_avb} | Duális mean hossz: {equal_djk_mean_len} | Független mean hossz: {equal_fugg_mean_len}",)
    
    dual_mean_avb = round(np.mean(resDf[resDf["ed_avb_diff"] < 0]["ed_avb_djk"]),4)
    dual_mean_len = round(resDf[resDf["ed_avb_diff"] < 0]["len_djk"].mean(), 2)
    print(f"Dualis a jobb {sum(resDf['ed_avb_diff'] < 0)} esetben | Átlag cost: {dual_mean_avb} | Átlag hossz: {dual_mean_len}")

def gather_results(network, resDf, resDictPathCfps):
    better_fugg_mean_avb = np.mean(resDf[resDf["ed_avb_diff"] > 0]["ed_avb_fugg"])
    fugg_mean_avb = resDf["ed_avb_fugg"].mean()
    fugg_mean_len = round(resDf[resDf["ed_avb_diff"] > 0]["len_fugg"].mean(), 2)

    equal_mean_avb = np.mean(resDf[resDf["ed_avb_diff"] == 0]["ed_avb_djk"])
    equal_djk_mean_len = round(resDf[resDf["ed_avb_diff"] == 0]["len_djk"].mean(), 2)
    equal_fugg_mean_len = round(resDf[resDf["ed_avb_diff"] == 0]["len_fugg"].mean(), 2)
    
    better_dual_mean_avb = np.mean(resDf[resDf["ed_avb_diff"] < 0]["ed_avb_djk"])
    dual_mean_avb = resDf["ed_avb_djk"].mean()
    dual_mean_len = round(resDf[resDf["ed_avb_diff"] < 0]["len_djk"].mean(), 2)    

    mean_pct_zero_djk = resDf[resDf["ed_avb_diff"] != 0]["simID"].apply(lambda x: resDictPathCfps[x]["djk_path"].count(0) / len(resDictPathCfps[x]["djk_path"])).mean()
    mean_pct_zero_fugg = resDf[resDf["ed_avb_diff"] != 0]["simID"].apply(lambda x: resDictPathCfps[x]["fuggetlen_path"].count(0) / len(resDictPathCfps[x]["fuggetlen_path"])).mean()
    
    d = {"network":network,
         "fugg_jobb":sum(resDf['ed_avb_diff'] > 0), "eq":sum(resDf['ed_avb_diff'] == 0), "dual_jobb": sum(resDf['ed_avb_diff'] < 0),
         "fugg_meanLen":fugg_mean_len, "eq_f_meanLen":equal_fugg_mean_len, "eq_d_meanLen":equal_djk_mean_len, "dual_meanLen":dual_mean_len,
         "fuggBet_meanEdAvb": better_fugg_mean_avb, "eq_meanEdAvb":equal_mean_avb, "dualBet_meanEdAvb":better_dual_mean_avb,
         "fugg_meanZpct":mean_pct_zero_fugg, "dual_meanZpct":mean_pct_zero_djk,
         "fugg_meanEdAvb":fugg_mean_avb, "dual_meanEdAvb":dual_mean_avb}
    
    return d

def plot_paths(g_base, g_base_pos, djk_path, djk_path_fuggetlen_edges,s,t,fgSize=(10,9)):
    options = {
        'node_size': 200,
        'node_color': 'gray',
        "width": [3 if edge in djk_path else 1 for edge in g_base.edges],
        "edge_color": ["orange" if edge in djk_path else "black" for edge in g_base.edges]
    }
    options2 = copy.deepcopy(options)
    options2["edge_color"] = ["purple" if edge in djk_path_fuggetlen_edges else "black" for edge in g_base.edges]
    options2["width"] = [3 if edge in djk_path_fuggetlen_edges else 1 for edge in g_base.edges]

    fig, ax = plt.subplots(ncols = 2, figsize=fgSize)
    plt.subplots_adjust(wspace=0.01)
    nx.draw_networkx(g_base, g_base_pos, **options, ax = ax[0])
    nx.draw_networkx(g_base, g_base_pos, **options2, ax = ax[1])

    ax[0].set_title("Edge-dual method")
    ax[1].set_title("Független method")

    fig.suptitle(f"Route {s} --> {t}")
    plt.show()