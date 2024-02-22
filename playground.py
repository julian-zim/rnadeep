from scipy.spatial.distance import pdist, squareform

# Define a custom distance metric function for strings
def custom_distance(str1, str2):
    # Example of a custom distance metric (Hamming distance)
    distance = sum(c1 != c2 for c1, c2 in zip(str1, str2))
    return distance

# Example dataset of strings
strings = ["hello", "world", "python", "science"]

# Compute pairwise distances using the custom metric
distances = pdist(strings, metric=custom_distance)

# Convert pairwise distances to a square matrix
distance_matrix = squareform(distances)

print("Distance matrix:")
print(distance_matrix)
