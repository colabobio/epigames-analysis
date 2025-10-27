from datetime import timedelta
import pandas as pd
import numpy as np
import networkx as nx

# Data parsing functions

# It calculates the weights in minutes
def get_contact_list(user_index, events, infections, def_contact_time, uuid_to_id=None, data_format=2, print_data_warnings=False):
    mili_to_seconds = 60 * 1000
    
    contacts = events[events["type"] == "contact"]

    node0 = contacts.user_id.values
    node1 = contacts.peer_id.values
    length = contacts.contact_length.values

    clist = {}
    for id0, id1, l01 in zip(node0, node1, length):
        n0 = user_index[id0]
        n1 = -1
        if 0 < data_format:
            if id1 in user_index:
                n1 = user_index[id1]
            elif print_data_warnings:
                print("Cannot find peer", id1)
        else:
            # Old format (sims before 2022): p2p id is in the infection column
            if id1 in uuid_to_id:
                n1 = user_index[uuid_to_id[id1]]
            elif print_data_warnings:
                print("Cannot find peer", id1)
    
        if -1 < n1:
            if n0 < n1:
                p01 = (n0, n1)
            else:
                p01 = (n1, n0)
            if p01 in clist:
                c = clist[p01]
            else: 
                c = 0
            
            clist[p01] = c + round(l01 / mili_to_seconds)

    for p in clist:
        clist[p] /= 2
    
    # Adding contacts from transmissions if they are not registered as contacts already
    for (n0, n1) in infections:
        if n0 < n1:
            p01 = (n0, n1)
        else:
            p01 = (n1, n0)
        if not p01 in clist:
            if 0 < def_contact_time:
                clist[p01] = def_contact_time
            if print_data_warnings: print("Cannot find contact between", n0, "and", n1)            

    return clist

def get_raw_contact_list(user_index, events, uuid_to_id=None, data_format=2, print_data_warnings=False):
    contacts = events[events["type"] == "contact"]

    node0 = contacts.user_id.values
    node1 = contacts.peer_id.values
    time = contacts.time.values    
    length = contacts.contact_length.values

    clist = []
    for id0, id1, t, l in zip(node0, node1, time, length):
        n0 = user_index[id0]
        n1 = -1
        if 0 < data_format:
            if id1 in user_index:
                n1 = user_index[id1]
            elif print_data_warnings:
                print("Cannot find peer", id1)
        else:
            # Old format (sims before 2022): p2p id is in the infection column
            if id1 in uuid_to_id:
                n1 = user_index[uuid_to_id[id1]]
            elif print_data_warnings:
                print("Cannot find peer", id1)
    
        if -1 < n1:
            clist += [(n0, n1, t, l)]

    return clist

def get_infection_list(user_index, events, discard_reinfections, time_delta_sec, uuid_to_id=None, data_format=2, print_data_warnings=False):
    infections = events[(events["type"] == "infection")]
    
    ilist = []
    itimes = {}
    infected = infections.user_id.values
    peers = infections.inf.values
    timestamp = infections.time.values
    for id1, peer0, ts in zip(infected, peers, timestamp):
        n1 = user_index[id1]
            
        if "PEER" in peer0:
            if 0 < data_format:
                # New schema
                id0 = int(peer0[peer0.index("[") + 1:peer0.index(":")])
                if id0 in user_index:
                    n0 = user_index[id0]                    
                    add_infection = True
                    for e in ilist:
                        if e[1] == n1:
                            if discard_reinfections:
                                add_infection = False
                                break
                            pid0 = index_user[e[0]]
                            ts0 = itimes[(pid0, id1)]
                            if abs(ts - ts0) <= time_delta_sec:
                                add_infection = False
                                if print_data_warnings:
                                    if pid0 == id0:
                                        print("Duplicated infection:", id1, "was already infected by", id0, "in the last", time_step_min, "minutes")
                                    else:
                                        print("Multiple infection:", id1, "is being infected by", id0, "but was already infected by", pid0, "in the last", time_step_min, "minutes")
                                break    

                    if add_infection: 
                        ilist += [(n0, n1)]
                        itimes[(id0, id1)] = ts
                elif print_data_warnings:
                    print("Cannot find peer", id0)                    
            else:    
                # Old format (sims before 2022): p2p id is in the infection column
                p2p0 = peer0[peer0.index("[") + 1:peer0.index(":")]
                if p2p0 in uuid_to_id:
                    id0 = uuid_to_id[p2p0]
                    if id0 in user_index:
                        n0 = user_index[id0]
                        if not (n0, n1) in ilist:                        
                            ilist += [(n0, n1)]
                        elif print_data_warnings:
                            print("Duplicated infection", id0, id1)  
                elif print_data_warnings:
                    print("Cannot find peer", p2p0)                        
            
    return ilist 

def get_node_state(user_index, events, state0, uuid_to_id=None, data_format=2, print_data_warnings=False):
    if state0 == None:
        state = [0] * len(user_index)
    else:
        state = state0
    
    inf = events[events["type"] == "infection"]
    infMap = pd.Series(inf.inf.values, index=inf.user_id).to_dict()
    for kid in infMap:
        src = infMap[kid]
        idx = user_index[kid]
        if "CASE0" in src:
            state[idx] = 1
        if "SOURCE" in src:
            state[idx] = 1
        if "PEER" in src:
            state[idx] = 2

            id0_str = src[5:].split(":")[0]
            if data_format == 0:
                if id0_str in uuid_to_id:
                    id0 = uuid_to_id[id0_str]
                else:
                    id0 = -1
                    if print_data_warnings:
                        print("Cannot find ID of infecting peer given its UUID", id0_str)
            else:
                id0 = int(id0_str)

            if -1 < id0:
                idx0 = user_index[id0]       
                if state[idx0] == 0:
                    state[idx0] = 1
                    if print_data_warnings:
                        print("Infecting peer did not have correct state", idx0)

    out = events[events["type"] == "outcome"]
    outMap = pd.Series(out.out.values, index=out.user_id).to_dict()
    for kid in outMap:
        out = outMap[kid]
        idx = user_index[kid]
        if out == "DEAD":
            state[idx] = 3
        if out == "RECOVERED":
            state[idx] = 4
        if out == "VACCINATED":
            state[idx] = 5
    
    return state

def get_user_name(users, uid, data_format=2):
    col = 'generated_id' if data_format == 1 else 'random_id'
    return users[users['id'] == uid][col].values[0]

def get_all_infection_events(events, users, tmin, data_format=2, print_data_warnings=False):
    etype = []
    user0 = []
    user1 = []
    etime = []
    
    infections = events[(events["type"] == "infection")]
    
    ilist = []
    itimes = {}
    infected = infections.user_id.values
    peers = infections.inf.values
    timestamp = infections.time.values
    for id1, peer0, ts in zip(infected, peers, timestamp):      
        tinf_minutes = round((ts-tmin)/60, 2)
        if "PEER" in peer0:
            id0 = None
            if 0 < data_format:
                # New schema
                id0 = int(peer0[peer0.index("[") + 1:peer0.index(":")])                
            else:    
                # Old format (sims before 2022): p2p id is in the infection column
                p2p0 = peer0[peer0.index("[") + 1:peer0.index(":")]
                if p2p0 in uuid_to_id:
                    id0 = uuid_to_id[p2p0]
            if id0: 
                etype += ['person-to-person infection']
                user0 += [get_user_name(users, id0, data_format)]
                user1 += [get_user_name(users, id1, data_format)]
                etime += [tinf_minutes]                
            if print_data_warnings:
                print('Cannot infecting peer of', id1)
        else: 
            etype += ['index case infection']
            user0 += [np.nan]
            user1 += [get_user_name(users, id1, data_format)]
            etime += [tinf_minutes]

    return pd.DataFrame({'Event': etype, 'Infector': user0, 'Infected': user1, 'Time': etime})

def get_all_illness_and_outcome_events(events, users, tmin, data_format=2):
    etype = []
    euser = []
    etime = []
    
    illness = events[(events["type"] == "illness")]
    userid = illness.user_id.values
    info = illness.inf.values
    timestamp = illness.time.values
    for uid, info, ts in zip(userid, info, timestamp):
        till_minutes = round((ts-tmin)/60, 2)
        if info == 'symptomatic':
            etype += ['Started showing symptoms']
        elif info == 'asymptomatic':
            etype += ['Stopped showing symptoms']
        else:
            etype += ['Unknown']            
        euser += [get_user_name(users, uid, data_format)]
        etime += [till_minutes]

    outcomes = events[(events["type"] == "outcome")]
    userid = outcomes.user_id.values
    userout = outcomes.out.values
    timestamp = outcomes.time.values
    for uid, info, ts in zip(userid, userout, timestamp):
        out_minutes = round((ts-tmin)/60, 2)
        etype += [info.capitalize()] 
        euser += [get_user_name(users, uid, data_format)]
        etime += [out_minutes]

    if data_format == 1:
        # Add vaccination events from OO Prod sims
        vaccine = events[(events["type"] == "vaccine")]
        userid = vaccine.user_id.values
        vaxinfo = vaccine.inf.values
        timestamp = vaccine.time.values
        for uid, info, ts in zip(userid, vaxinfo, timestamp):
            vax_minutes = round((ts-tmin)/60, 2)
            if info == 'immunized':
                edesc = 'Vaccinated and immunized'
            elif info == 'notImmunized':
                edesc = 'Vaccinated and NOT immunized'
            else:
                edesc = 'Unknown'            
            etype += [edesc] 
            euser += [get_user_name(users, uid, data_format)]
            etime += [vax_minutes]        
    
    return pd.DataFrame({'Event': etype, 'user': euser, 'Time': etime})

# Some utilities

# https://stackoverflow.com/a/48938464
def hour_rounder(t):
    # Rounds to nearest hour by adding a timedelta hour if minute >= 30
    return (t.replace(second=0, microsecond=0, minute=0, hour=t.hour) + timedelta(hours=t.minute//30))

# Create the network, skipping edges between nodes that spend less than min_contact_time
# in contact during the entire sim
def create_contact_network(user_index, contacts, state_values, state_name="state", minw=0):
    nodes = [i for i in range(0, len(user_index))]
    edges = []
    weights = []    
    if 0 < len(contacts):
        for p in contacts:
            n0 = p[0]
            n1 = p[1]
            w = contacts[p]            
            if minw < w:
                edges += [(n0, n1)]
                weights += [w]
    
    g = nx.Graph()
    g.add_nodes_from(nodes)
    g.add_weighted_edges_from([(edges[i][0], edges[i][1], weights[i]) for i in range(len(edges))])
    state_dict = dict(zip(g.nodes(), state_values))    
    nx.set_node_attributes(g, state_dict, state_name)
    return g

def remove_nodes_with_less_edges(G, k):
    nodes_to_remove = [node for node, degree in dict(G.degree()).items() if degree < k]
    G.remove_nodes_from(nodes_to_remove)
    return nodes_to_remove

def get_modifiers_list(events):
    modifiers = events[(events["type"] == "modifier")]
    return modifiers

def get_score_events(events):
    scores = events[(events["type"] == "score")]
    return scores

def get_modifier_count(events, name):
    return len(events[(events["modifier"] == name)])    

def get_info_count(events, name):
    return len(events[(events["inf"] == name)])

def get_info_users(events, name):
    sel_events = events[(events["inf"] == name)]
    return sel_events["user_id"].tolist()

def get_node_index(users, user_index, node_labels):
    idxs = [user_index[x] for x in users]
    return [node_labels.index(x) for x in idxs if x in node_labels]    