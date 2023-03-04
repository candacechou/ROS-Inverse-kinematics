#! /usr/bin/env python3
import math
import numpy as np
from numpy.linalg import inv, multi_dot,pinv
"""
    # {Chia-Hsuan CHou}
    # {chchou@kth.se}
"""

def scara_IK(point):
    x = point[0]
    y = point[1]
    z = point[2]
    q = [0.0, 0.0, 0.0]
    I_0 = 0.07
    I_1 = 0.3
    I_2 = 0.35
    w = math.atan2(y,x-I_0)
    line = math.sqrt((x-I_0)*(x-I_0)+y*y)
    foo = ((x-I_0)*(x-I_0) + y*y -I_1*I_1-I_2*I_2)/(2*I_1*I_2)
    #theta_1 = math.asin(foo)
    #theta_2 = math.asin(I_2*foo/I_1)
    q[0] = w-math.acos(((x-I_0)*(x-I_0) + y*y +I_1*I_1-I_2*I_2)/(2*I_1*line))
    q[1] = math.acos(foo)
    q[2] = z
    

    return q

def kuka_IK(point, R, joint_positions):
    x = point[0]
    y = point[1]
    z = point[2]
    q = joint_positions
    L = 0.4
    M = 0.39
    a = [0,0,0,0,0,0,0]
    d = [0,0,L,0,M,0,0]
    alpha = [np.pi/2,
            -np.pi/2,
            -np.pi/2,
            np.pi/2,
            np.pi/2,
            -np.pi/2,
            0]
    error_x = np.array([[1],[1],[1],[1],[1],[1]])
    error_to_x = 0.01
    error_to_w = 0.01
    error_w = [1,1,1,1,1,1,1]
    x = 2
    while(x!=1):
        x = 1
        TB = np.eye(4)
        TB[2][3] = 0.311
        T0 = trans_matrix(a[0],alpha[0],d[0],q[0])
        T1 = trans_matrix(a[1],alpha[1],d[1],q[1])
        T2 = trans_matrix(a[2],alpha[2],d[2],q[2])
        T3 = trans_matrix(a[3],alpha[3],d[3],q[3])
        T4 = trans_matrix(a[4],alpha[4],d[4],q[4])
        T5 = trans_matrix(a[5],alpha[5],d[5],q[5])
        T6 = trans_matrix(a[6],alpha[6],d[6],q[6])
        T7 = np.eye(4)
        T7[2][3] = 0.078
        J1 = T0
        J2 = np.linalg.multi_dot([T0 ,T1])
        J3 = np.linalg.multi_dot([T0 , T1 , T2])
        J4 = np.linalg.multi_dot([T0 , T1 , T2 , T3])
        J5 = np.linalg.multi_dot([T0 , T1 , T2 , T3 , T4])
        J6 = np.linalg.multi_dot([T0 , T1 , T2 , T3 , T4 , T5])
        J7 = np.linalg.multi_dot([T0 , T1 , T2 , T3 , T4 , T5 , T6])
        J_7_back = np.linalg.multi_dot([J7,T7])
        Tx = np.eye(4)
        p = J_7_back[0:3,3].T
        z1 = Tx[0:3,2]
        p1 = Tx[0:3,3]
        J_1 = np.cross(z1,(p-p1))
    
        J_1 = np.hstack((J_1 ,z1))
    
    #### jacobian 2 
        p2 = J1[0:3,3]
        z2 = J1[0:3,2]
        J_2 = np.cross(z2,(p-p2))
  
        J_2 = np.hstack((J_2,z2))
    #print(J_2)   
    #### jacobian 3 
        p3 = J2[0:3,3]
        z3 = J2[0:3,2]
        J_3 = np.cross(z3,(p-p3))
    #print(J_3)
        J_3 = np.hstack((J_3,z3))   
    #### jacobian 4
        p4 = J3[0:3,3]
        z4 = J3[0:3,2]
        J_4 = np.cross(z4,(p-p4))
    #print(J_4)
        J_4 = np.hstack((J_4,z4))
    #print('z4')
    
    #### jacobian 5 
        p5 = J4[0:3,3]
        z5 = J4[0:3,2]
        J_5 = np.cross(z5,(p-p5))
    #print(J_5)
        J_5 = np.hstack((J_5,z5))
    #print('z5')
    
    #### jacobian 6
        p6 = J5[0:3,3]
        z6 = J5[0:3,2]
        J_6 = np.cross(z6,(p-p6))
        J_6 = np.hstack((J_6,z6))
    #print('z6')
    #print(J_6)
    #### jacobian 7
        p7 = J6[0:3,3]
        z7 = J6[0:3,2]
        J_7 = np.cross(z7,(p-p7))
        J_7 = np.hstack((J_7,z7))
    #print('z7')
    #print(J_7)
        J = np.vstack((J_1,J_2,J_3,J_4,J_5,J_6,J_7))
        J =J.T
        J_7_back = np.linalg.multi_dot([TB,J_7_back])
        x1 = -J_7_back[0][3] + point[0]
        y1 = -J_7_back[1][3] + point[1]
        z1 = -J_7_back[2][3] + point[2]
        nd = np.array([
            [R[0][0]],
            [R[1][0]],
            [R[2][0]]
        ])
        sd = np.array([
            [R[0][1]],
            [R[1][1]],
            [R[2][1]]
        ])
        ad = np.array([
            [R[0][2]],
            [R[1][2]],
            [R[2][2]]
        ])
    
        ne = J_7_back[:3,[0]]
        se = J_7_back[:3,[1]]
        ae = J_7_back[:3,[2]]
        Euler_error = 0.5*np.cross(ne.T, nd.T).T + np.cross(se.T, sd.T).T + np.cross(ae.T, ad.T).T
        point_error = np.array([[x1],[y1],[z1]])
        current = np.vstack((point_error,Euler_error))
        invjac = pinv(J)
        error_w = np.dot(invjac , current)
        q = np.array([q]).T
        q = q + error_w
        for i in range(len(current)):
            
            if abs(current[i]) > 0.01:
                x = x + 1
        
        q = q.T
        q = q.tolist()
        q = q[0]
    return q

    
def trans_matrix(a,alpha,d,theta):
    T =np.zeros((4,4))
    T[0][0] = np.cos(theta)
    T[0][1] = -np.sin(theta) * np.cos(alpha)
    T[0][2] = np.sin(theta) * np.sin(alpha)
    T[0][3] = a * np.cos(theta)
    T[1][0] = np.sin(theta)
    T[1][1] = np.cos(theta) * np.cos(alpha)
    T[1][2] = -np.cos(theta) * np.sin(alpha)
    T[1][3] = np.sin(alpha) * a
    T[2][0] = 0
    T[2][1] = np.sin(alpha)
    T[2][2] = np.cos(alpha)
    T[2][3] = d
    T[3][3] = 1
    return T