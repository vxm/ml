import bpy
from random import random, seed
import numpy as np
import colorsys

range_limit = 10
point_count = 9650
n_clusters = 10
K_MEANS_ITERATIONS = 10


# Define a function to generate distinct colors
def generate_distinct_colors(n, alpha=1.0, saturation=1.0, value=1.0):
    """
    Generate n distinct colors by evenly spacing hue values.
    """
    colors = []
    for i in range(n):
        hue = i / n
        value = value
        r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
        colors.append((r, g, b, alpha))
    return colors


# Create a list of colors for your clusters (RGB and Alpha)
colors = generate_distinct_colors(n_clusters,
                                  random() * 0.2 + 0.7,
                                  random() * 0.1 + 0.7,
                                  random() * 0.3 + 0.4)

n_clusters = len(colors)

positions = [((random()-0.5) * range_limit,
              (random()-0.5) * range_limit,
              (random()-0.5) * range_limit)
             for _ in range(point_count)]


def k_means_clustering(positions, n_clusters):
    # Convert positions to a numpy array for easier manipulation
    points = np.array(positions)

    # Randomly initialize centroids
    centroids = points[np.random.choice(points.shape[0],
                                        n_clusters,
                                        replace=False), :]
    print(centroids)
    cluster_assignments = [0] * len(points)

    for _ in range(K_MEANS_ITERATIONS):
        clusters = {}
        for i in range(n_clusters):
            clusters[i] = []
        for p, point in enumerate(points):
            distances = np.linalg.norm(point - centroids, axis=1)
            cluster_assignment = np.argmin(distances)
            clusters[cluster_assignment].append(point)
            cluster_assignments[p] = cluster_assignment
        for i in range(n_clusters):
            centroids[i] = np.mean(clusters[i], axis=0)
            print(centroids[i])

    return centroids, clusters, cluster_assignments


def filter_points_by_distance(points, reference_point, max_distance):
    """
    Filter a list of points, returning only those
    within a certain distance from a reference point.
    Returns:
        list: A list of points within the specified
        maximum distance from the reference point.
    """
    filtered_points = []
    for point in points:
        # Calculate Euclidean distance from the reference point
        distance = np.linalg.norm(np.array(point) - np.array(reference_point))
        if distance < max_distance:
            filtered_points.append(point)
    return filtered_points


positions = filter_points_by_distance(positions, (0, 0, 0), range_limit/1.9)


def create_blender_material(name, color, metal=0.5, roughness=0.1):
    """
    Create a material with specified color, metallicity, and roughness.
    Parameters:
        name (str): The name of the new material.
        color (tuple): The base color of the material as
            a tuple of (R, G, B, Alpha), with each value from 0 to 1.
        metal (float, optional): The metallicity of the material,
            ranging from 0 to 1. Defaults to 0.5.
        roughness (float, optional): The roughness of the material,
            ranging from 0 to 1. Defaults to 0.5.
    Returns:
        bpy.types.Material: The created material.
    """
    # Create a new material
    material = bpy.data.materials.new(name=name)
    material.use_nodes = True
    bsdf = material.node_tree.nodes.get('Principled BSDF')
    # Set the material properties
    bsdf.inputs['Base Color'].default_value = color
    bsdf.inputs['Metallic'].default_value = metal
    bsdf.inputs['Roughness'].default_value = roughness
    return material


# Example usage with 2 clusters
centroids, clusters, cluster_assignments = k_means_clustering(positions, n_clusters)


# Create materials for each cluster
materials = [create_blender_material(f"ClusterMaterial_{i}", color, random(), random()) for i, color in enumerate(colors)]


# Create a sphere for each position, assigning a material based on cluster
for idx, pos in enumerate(positions):
    cluster_idx = cluster_assignments[idx]
    seed(cluster_idx)
    radius = (random() * 0.03) + 0.09
    bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=pos)
    obj = bpy.context.object  # Get the newly created object
    obj.data.materials.append(materials[cluster_idx])

print("cluster_assignments:", cluster_assignments)
