import numpy as np

def generate_sphere_points(radius, num_points):
    points = []
    increment = np.pi * (3 - np.sqrt(5))
    
    for i in range(num_points):
        y = ((i * 2.0) / num_points) - 1.0 + (1.0 / num_points)
        r = np.sqrt(1 - y * y)
        phi = (i % num_points) * increment
        
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        
        points.append([x * radius, y * radius, z * radius])
    
    return np.array(points)

radius = 6371
num_points = 1000

sphere_points = generate_sphere_points(radius, num_points)
