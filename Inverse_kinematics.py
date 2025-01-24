import numpy as np


def rotation_matrix(psi, theta, phi):
    # Compute the rotation matrix P_RB
    cos_psi = np.cos(psi)
    sin_psi = np.sin(psi)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)

    PRB = np.array([
        [cos_psi * cos_theta, -sin_psi * cos_phi + cos_psi * sin_theta * sin_phi,
         sin_psi * sin_phi + cos_psi * sin_theta * cos_phi],
        [sin_psi * cos_theta, cos_psi * cos_phi + sin_psi * sin_theta * sin_phi,
         -cos_psi * sin_phi + sin_psi * sin_theta * cos_phi],
        [-sin_theta, cos_theta * sin_phi, cos_theta * cos_phi]
    ])

    return PRB


def compute_li(T, psi, theta, phi, p_i, b_i):
    PRB = rotation_matrix(psi, theta, phi)
    # Compute l_i as T + PRB * p_i - b_i
    l_i = T + np.dot(PRB, p_i) - b_i
    return l_i


def compute_li_length(T, psi, theta, phi, p_i, b_i):
    l_i = compute_li(T, psi, theta, phi, p_i, b_i)
    # Compute the length (magnitude) of l_i
    length = np.linalg.norm(l_i)
    return length



#Leg vectors
#leg name -> leg[0] = platform & leg[1] = base
# List of l_i vectors
l_vectors = [
[np.array([537.69, 515.43, 0]), np.array([177.53, 723.37, 0])],#C
[np.array([-537.69, 515.43, 0]), np.array([-177.53, 723.37, 0])],#B
    [np.array([-715.23, 207.94, 0]), np.array([-715.23, -207.94, 0])],#A
    [np.array([-177.53, -723.37, 0]), np.array([-537.69, -515.43, 0])],  # F
    [np.array([177.53, -723.37, 0]), np.array([537.69, -515.43, 0])],  # E
    [np.array([715.23, 207.94, 0]), np.array([715.23, -207.94, 0])]#D
]
T = np.array([0, 0, 1079])  # Start height of the platform
psi, theta, phi = np.radians([0, -2, 0])  # Ya1w, Roll, Pitch in degrees

# Loop over each l vector
for i, l in enumerate(l_vectors):
    p_i = l[0]
    b_i = l[1]

    # Calculate the vector l_i
    l_i = compute_li(T, psi, theta, phi, p_i, b_i)
    # Calculate the length (magnitude) of l_i
    li_length = compute_li_length(T, psi, theta, phi, p_i, b_i)
    print(f"Length of l_i for l{i}:", li_length-1156.372420286821)