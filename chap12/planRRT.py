import numpy as np
from message_types.msg_waypoints import msg_waypoints
from parameters import planner_parameters as PLN

class Node():
    def __init__(self, ned=[0,0,0], parent_id=0, cost=0, connectToGoalFlag=0, 
                 isGoal=0):
        self.ned = np.array(ned)
        self.parent_id = parent_id
        self.cost = cost
        self.connectToGoalFlag = connectToGoalFlag
        self.isGoal = isGoal
  
class planRRT():
    def __init__(self, map):
        self.waypoints = msg_waypoints()
        self.segmentLength = 300 # standard length of path segments
        self.dx = 5
        self.min_distance = 50
    def planPath(self, wpp_start, wpp_end, map):

        # desired down position is down position of end node
        pd = wpp_end.item(2)

        # specify start and end nodes from wpp_start and wpp_end
        # format: N, E, D, cost, parentIndex, connectsToGoalFlag,
        start_node = np.array([wpp_start.item(0), wpp_start.item(1), pd, 0, 0, 0])
        end_node = np.array([wpp_end.item(0), wpp_end.item(1), pd, 0, 0, 0])

        # establish tree starting with the start node
        tree = start_node

        # check to see if start_node connects directly to end_node
        if ((np.linalg.norm(start_node[0:3] - end_node[0:3]) < self.segmentLength ) and not self.collision(start_node, end_node, map)):
            self.waypoints.ned = end_node[0:3]
        else:
            numPaths = 0
            while numPaths < 3:
                tree, flag = self.extendTree(tree, end_node, map, pd)
                numPaths = numPaths + flag


        # find path with minimum cost to end_node
        path = self.findMinimumPath(tree, end_node, start_node)
        return self.smoothPath(path, map)

    def generateRandomNode(self, map, pd):
        pn =  np.random.uniform(0, PLN.city_width)
        pe = np.random.uniform(0,PLN.city_width)
        return np.array([pn,pe,pd])


    def collision(self, start_node, end_node, map):
        points = self.pointsAlongPath(start_node,end_node, self.segmentLength)
        radius = np.sqrt(2)*map.building_width + self.min_distance   
        for i in points:
            for j in map.num_city_blocks:
                b_n = map.building_north[j]
                b_e = map.building_east[j]
                p_n = points[i][0]
                p_e = points[i][1]
                cllsn_dis = np.sqrt((p_n-b_n)**2 + (p_e-b_e)**2)
                if cllsn_dis <= radius:
                    return True
        return False

    def pointsAlongPath(self, start_node, end_node, Del):
        D = Del/5
        vec = end_node - start_node
        norm = np.linalg.norm(end_node-start_node)
        q = vec/norm
        points = np.array([])
        for i in range(1,D):
            point = np.array([start_node + q*i])
            # points.append(point)
            np.append(points, point)
        return points

    # def downAtNE(map, n, e):

    def v_plus(self, p,tree):
        tr_p = tree[:,0:3]
        lengths = []
        for i in tr_p:
            lengths.append(np.linalg.norm(p-tr_p[i]))
        p_index = np.argmin(lengths)
        parent = np.array([tr_p[p_index]])
        q = (p-parent[0:3])/np.linalg.norm(p-parent[0:3])        
        v_p = parent + q*self.segmentLength
        v_p_cost = self.segmentLength + tree[p_index][3]
        v_plus = np.array([ v_p[0], v_p[1], v_p[2], v_p_cost, p_index, 0])      
        return v_plus
    
    def extendTree(self,tree, end_node, map, pd):           
        p = self.generateRandomNode(map, pd)        
        v_plus = self.v_plus(p, tree)
        D = np.linalg.norm(end_node-v_plus)
        if not self.collision(p, v_plus, map):
            np.append(tree,v_plus)

        while D <= self.segmentLength:
            p = self.generateRandomNode(map, pd)        
            v_plus = self.v_plus(p, tree)
            D = np.linalg.norm(end_node-v_plus)
            if not self.collision(p, v_plus, map):
                np.append(tree,v_plus)            
        tree[-1][5] = 1 
        np.append(tree, end_node)
        return tree,1

    def findMinimumPath(self, tree, end_node, start_node):
        last_nodes = tree[:,-1]==1
        min_value= last_nodes[0][3]        
        for i in last_nodes:
            if last_nodes[i][3] < min_value:
                min_value = last_nodes[i][3]
                index = i
        path = np.array([end_node])
        point = last_nodes[i]

        while point[3]!=0:
            np.apeend(path, point )
            point = tree[point[4]]
        
        np.append(path, start_node)
        path_NED = path[:,0:3]
        return path_NED

    def smoothPath(self, path, map):
        path = np.flip(path,0)
        sm_path = np.array([])
        current_point = path[0]
        considering_point = path[1]
        i = 1
        while considering_point == path[-1]:
            while not self.collision(current_point,considering_point, map):
                i = i+1
                considering_point = path[i]
                
            np.append(sm_path, current_point)
            current_point = considering_point
        np.append(sm_path, path[-1])



