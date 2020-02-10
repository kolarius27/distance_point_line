from math import sqrt, acos, pi, inf
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

def points_collector():
    while True:
        try:
            pointA = (float(input("X-coordinate of point A: ")), float(input("Y-coordinate of point A: ")))
            break
        except ValueError:
            continue
    print("A = [", pointA[0], ",", pointA[1], "]")
    while True:
        try:
            pointB = (float(input("X-coordinate of point B: ")), float(input("Y-coordinate of point B: ")))
            break
        except ValueError:
            continue
    print("B = [", pointB[0], ",", pointB[1], "]")
    while True:
        try:
            pointP = (float(input("X-coordinate of point P: ")), float(input("Y-coordinate of point P: ")))
            break
        except ValueError:
            continue
    print("P = [", pointP[0], ",", pointP[1], "]")
    return pointA, pointB, pointP

def points_collector_3D():
    while True:
        try:
            pointA = (float(input("X-coordinate of point A: ")), float(input("Y-coordinate of point A: ")), float(input("Z-coordinate of point A: ")))
            break
        except ValueError:
            continue
    print("A = [", pointA[0], ",", pointA[1], ",", pointA[2], "]")
    while True:
        try:
            pointB = (float(input("X-coordinate of point B: ")), float(input("Y-coordinate of point B: ")), float(input("Z-coordinate of point B: ")))
            break
        except ValueError:
            continue
    print("B = [", pointB[0], ",", pointB[1], ",", pointB[2], "]")
    while True:
        try:
            pointP = (float(input("X-coordinate of point P: ")), float(input("Y-coordinate of point P: ")), float(input("Z-coordinate of point P: ")))
            break
        except ValueError:
            continue
    print("P = [", pointP[0], ",", pointP[1], ",", pointP[2], "]")
    return pointA, pointB, pointP

def distance_point_line_calculator(point_A, point_B, point_P):
    plt.plot([point_A[0], point_B[0]], [point_A[1], point_B[1]])
    plt.plot([point_A[0], point_B[0]], [point_A[1], point_B[1]], 'bo')
    plt.plot(point_P[0], point_P[1], 'ro')

    vector_AB = (point_B[0] - point_A[0], point_B[1] - point_A[1])
    vector_AP = (point_P[0] - point_A[0], point_P[1] - point_A[1])
    vector_BP = (point_P[0] - point_B[0], point_P[1] - point_B[1])
    vector_BA = (point_A[0] - point_B[0], point_A[1] - point_B[1])

    dist_AB = sqrt(vector_AB[0] * vector_AB[0] + vector_AB[1] * vector_AB[1])
    dist_AP = sqrt(vector_AP[0] * vector_AP[0] + vector_AP[1] * vector_AP[1])
    dist_BP = sqrt(vector_BP[0] * vector_BP[0] + vector_BP[1] * vector_BP[1])

    cos_alfa = (vector_AP[0] * vector_AB[0] + vector_AP[1] * vector_AB[1]) / (dist_AP * dist_AB)
    cos_beta = (vector_BP[0] * vector_BA[0] + vector_BP[1] * vector_BA[1]) / (dist_BP * dist_AB)

    if cos_alfa > 1:
        alfa = acos(1)
    elif cos_alfa < -1:
        alfa = acos(-1)
    else:
        alfa = acos(cos_alfa)
    if cos_beta > 1:
        beta = acos(1)
    elif cos_beta < -1:
        beta = acos(-1)
    else:
        beta = acos(cos_beta)

    if alfa > (pi / 2):
        dist = sqrt(vector_AP[0] * vector_AP[0] + vector_AP[1] * vector_AP[1])
        plt.plot([point_P[0], point_A[0]], [point_P[1], point_A[1]], '--r')
    elif beta > (pi / 2):
        dist = sqrt(vector_BP[0] * vector_BP[0] + vector_BP[1] * vector_BP[1])
        plt.plot([point_P[0], point_B[0]], [point_P[1], point_B[1]], '--r')
    else:
        dist = abs(vector_AB[0] * vector_AP[1] - vector_AB[1] * vector_AP[0]) / sqrt(
            vector_AB[0] * vector_AB[0] + vector_AB[1] * vector_AB[1])
        point_C = (point_P[0]-vector_AB[1]*cos_alfa*dist_AP/dist_AB, point_P[1]+vector_AB[0]*cos_alfa*dist_AP/dist_AB)
        plt.plot([point_P[0], point_C[0]], [point_P[1], point_C[1]], '--r')
    plt.title("Distance between point P and line segment AB = {}".format(round(dist, 2)))
    plt.grid()
    plt.show()
    return dist

def distance_point_line_calculator_3D(point_A, point_B, point_P):
    ax = plt.axes(projection='3d')
    ax.plot3D([point_A[0], point_B[0]], [point_A[1], point_B[1]], [point_A[2], point_B[2]])
    ax.plot3D([point_A[0], point_B[0]], [point_A[1], point_B[1]], [point_A[2], point_B[2]], 'bo')
    ax.plot3D([point_P[0]], [point_P[1]], [point_P[2]], 'ro')

    vector_AB = (point_B[0] - point_A[0], point_B[1] - point_A[1], point_B[2] - point_A[2])
    vector_AP = (point_P[0] - point_A[0], point_P[1] - point_A[1], point_B[2] - point_A[2])
    vector_BP = (point_P[0] - point_B[0], point_P[1] - point_B[1], point_B[2] - point_A[2])
    vector_BA = (point_A[0] - point_B[0], point_A[1] - point_B[1], point_B[2] - point_A[2])

    dist_AB = sqrt(vector_AB[0] * vector_AB[0] + vector_AB[1] * vector_AB[1] + vector_AB[2] * vector_AB[2])
    dist_AP = sqrt(vector_AP[0] * vector_AP[0] + vector_AP[1] * vector_AP[1] + vector_AP[2] * vector_AP[2])
    dist_BP = sqrt(vector_BP[0] * vector_BP[0] + vector_BP[1] * vector_BP[1] + vector_BP[2] * vector_BP[2])


    cos_alfa = (vector_AP[0] * vector_AB[0] + vector_AP[1] * vector_AB[1] + vector_AP[2] * vector_AB[2]) / (dist_AP * dist_AB)
    cos_beta = (vector_BP[0] * vector_BA[0] + vector_BP[1] * vector_BA[1] + vector_BP[2] * vector_BA[2]) / (dist_BP * dist_AB)

    if cos_alfa > 1:
        alfa = acos(1)
    elif cos_alfa < -1:
        alfa = acos(-1)
    else:
        alfa = acos(cos_alfa)
    if cos_beta > 1:
        beta = acos(1)
    elif cos_beta < -1:
        beta = acos(-1)
    else:
        beta = acos(cos_beta)

    if alfa > (pi / 2):
        dist = dist_AP
        ax.plot3D([point_A[0], point_P[0]], [point_A[1], point_P[1]], [point_A[2], point_P[2]], '--r')
    elif beta > (pi / 2):
        dist = dist_BP
        ax.plot3D([point_B[0], point_P[0]], [point_B[1], point_P[1]], [point_B[2], point_P[2]], '--r')
    else:
        vector_APxAB = np.cross(vector_AP, vector_AB)
        dist = sqrt(vector_APxAB[0] * vector_APxAB[0] + vector_APxAB[1] * vector_APxAB[1] +
                    vector_APxAB[2] * vector_APxAB[2]) / dist_AB
        coef = (dist_AP*dist_AP - dist_BP*dist_BP + dist_AB*dist_AB)/(2*dist_AB*dist_AB)
        point_C = (point_A[0] + vector_AB[0] * coef,
                   point_A[1] + vector_AB[1] * coef,
                   point_A[2] + vector_AB[2] * coef)
        plt.plot([point_P[0], point_C[0]], [point_P[1], point_C[1]], [point_P[2], point_C[2]], '--r')
    plt.title("Distance between point P and line segment AB = {}".format(round(dist, 2)))
    plt.grid()
    plt.show()
    return dist


dimension = input("2D or 3D? ")
while dimension not in ("2D", "3D"):
    print("Please type '2D' or '3D': ")
    dimension = input("2D or 3D? ")

if dimension == "2D":
    A_point, B_point, P_point = points_collector()
    distance = distance_point_line_calculator(A_point, B_point, P_point)
elif dimension == "3D":
    A_point, B_point, P_point = points_collector_3D()
    distance = distance_point_line_calculator_3D(A_point, B_point, P_point)

print("Distance between point P and line segment AB is ", distance, "≐", round(distance, 2), ".")

