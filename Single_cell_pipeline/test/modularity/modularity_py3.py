# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 17:43:55 2019
Test only
@author: huangzhihao
"""

def inputdata(filePath):
    f = open(filePath,'r')
    vector_dict = {}
    edge_dict = {}
    for line in f.readlines():
        lines = line.strip().split()
        for i in range(2):
            if lines[i] not in vector_dict:
                #put the vertex into the vector_dict
                vector_dict[lines[i]] = True
                #put the edges into the edge_dict; initialize as 1
                edge_list = []
                edge_list.append(lines[1-i]+":"+"1")
                edge_dict[lines[i]] = edge_list
            else:
                edge_list = edge_dict[lines[i]]
                edge_list.append(lines[1-i]+":"+"1")
                edge_dict[lines[i]] = edge_list
    return(vector_dict, edge_dict)

def modularity(vector_dict, edge_dict):
    Q = 0.0
    # m represents the total wight
    m = 0
    for i in edge_dict.keys():
        edge_list = edge_dict[i]
        for j in range(len(edge_list)):
            l = edge_list[j].strip().split(":")
            m += float(l[1].strip())
    # cal community of every vertex
    #find member in every community
    community_dict = {}
    for i in vector_dict.keys():
        if vector_dict[i] not in community_dict:
            community_list = []
        else:
            community_list = community_dict[vector_dict[i]]

        community_list.append(i)
        community_dict[vector_dict[i]] = community_list
    for i in community_dict.keys():
        sum_in = 0.0
        sum_tot = 0.0
        #vector num
        vector_list = community_dict[i]
        #print "vector_list : ", vector_list
        #two loop cal inner link
        if len(vector_list) == 1:
            tmp_list = edge_dict[vector_list[0]]
            tmp_dict = {}
            for link_mem in tmp_list:
                l = link_mem.strip().split(":")
                tmp_dict[l[0]] = l[1]
            #判断该顶点是否在连接的顶点所在的社区里
            if vector_list[0] in tmp_dict:
                sum_in = float(tmp_dict[vector_list[0]])
            else:
                sum_in = 0.0
        else:
            for j in range(0,len(vector_list)):
                link_list = edge_dict[vector_list[j]]
                tmp_dict = {}
                for link_mem in link_list:
                    l = link_mem.strip().split(":")
                    #split the vertex and weight
                    tmp_dict[l[0]] = l[1]
                for k in range(0, len(vector_list)):
                    if vector_list[k] in tmp_dict:
                        sum_in += float(tmp_dict[vector_list[k]])

        #cal degree
        for vec in vector_list:
            link_list = edge_dict[vec]
            for i in link_list:
                l = i.strip().split(":")
                sum_tot += float(l[1])        
        Q += ((sum_in / m) - (sum_tot/m)*(sum_tot/m))
    return(Q)

def chage_community(vector_dict, edge_dict, Q):
    vector_tmp_dict = {}
    for key in vector_dict:
        vector_tmp_dict[key] = vector_dict[key]
    #for every vertex chose it's neighbor
    for key in vector_tmp_dict.keys():
        neighbor_vector_list = edge_dict[key]
        for vec in neighbor_vector_list:
            ori_com = vector_tmp_dict[key]
            vec_v = vec.strip().split(":")
            #vec_v[0] is one of the neighbors
            #compare the list_member with ori_com
            #if the first value of key(for example 64:1) has a higher Q, 64 temporarily belong to 1.  
            if ori_com != vector_tmp_dict[vec_v[0]]:
                vector_tmp_dict[key] = vector_tmp_dict[vec_v[0]]
                Q_new = modularity(vector_tmp_dict, edge_dict)
                print(vector_tmp_dict)
                print(Q_new)
                if (Q_new - Q) > 0:
                    Q = Q_new
                else:
                    vector_tmp_dict[key] = ori_com
    return(vector_tmp_dict, Q)

def modify_community(vector_dict):
    #modify the community
    community_dict = {}
    community_num = 0
    for community_values in vector_dict.values():
        if community_values not in community_dict:
            community_dict[community_values] = community_num
            community_num += 1
    for key in vector_dict.keys():
        vector_dict[key] = community_dict[vector_dict[key]]
    return(community_num)

def louvain_test(vector_dict, edge_dict):
    for i in vector_dict.keys():
        vector_dict[i] = i
    Q = modularity(vector_dict, edge_dict)  
    #2. for every vertex, chose the community
    Q_new = 0.0
    while (Q_new != Q):
        Q_new = Q
        vector_dict, Q = chage_community(vector_dict, edge_dict, Q)
    community_num = modify_community(vector_dict)
    print("Best Modularity = ", Q)
    print('cell numbers: ',len(edge_dict.keys()))
    print('community numbers: ',community_num)

if __name__ == "__main__":    
    vector_dict, edge_dict=inputdata("E:\project\hierarchical cluster\SNN_edge_file_noedge.txt")
    louvain_test(vector_dict, edge_dict)
