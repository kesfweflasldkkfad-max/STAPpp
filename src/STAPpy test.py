import math

import numpy as np


E = 1000.0
NU = 0.3


def constitutive_matrix_plane_stress(e_modulus: float, poisson_ratio: float) -> np.ndarray:
    factor = e_modulus / (1.0 - poisson_ratio**2)
    return factor * np.array(
        [
            [1.0, poisson_ratio, 0.0],
            [poisson_ratio, 1.0, 0.0],
            [0.0, 0.0, (1.0 - poisson_ratio) / 2.0],
        ],
        dtype=float,
    )


def exact_displacement(x: float, y: float) -> np.ndarray:
    return np.array([0.01 * x, -0.003 * y], dtype=float)


def exact_strain() -> np.ndarray:
    return np.array([0.01, -0.003, 0.0], dtype=float)


def shape_function_gradients(xi: float, eta: float) -> tuple[np.ndarray, np.ndarray]:
    dndxi = 0.25 * np.array(
        [-(1.0 - eta), 1.0 - eta, 1.0 + eta, -(1.0 + eta)],
        dtype=float,
    )
    dndeta = 0.25 * np.array(
        [-(1.0 - xi), -(1.0 + xi), 1.0 + xi, 1.0 - xi],
        dtype=float,
    )
    return dndxi, dndeta


def b_matrix(element_coords: np.ndarray, xi: float, eta: float) -> tuple[np.ndarray, float]:
    dndxi, dndeta = shape_function_gradients(xi, eta)
    natural_gradients = np.vstack((dndxi, dndeta))
    jacobian = natural_gradients @ element_coords
    det_jacobian = np.linalg.det(jacobian)

    if det_jacobian <= 0.0:
        raise ValueError(f"Element has a non-positive Jacobian determinant: {det_jacobian}")

    cartesian_gradients = np.linalg.solve(jacobian, natural_gradients)
    dndx = cartesian_gradients[0]
    dndy = cartesian_gradients[1]

    b = np.zeros((3, 8), dtype=float)
    for local_node in range(4):
        column = 2 * local_node
        b[0, column] = dndx[local_node]
        b[1, column + 1] = dndy[local_node]
        b[2, column] = dndy[local_node]
        b[2, column + 1] = dndx[local_node]

    return b, det_jacobian


def q4_element_stiffness(element_coords: np.ndarray, constitutive: np.ndarray) -> np.ndarray:
    stiffness = np.zeros((8, 8), dtype=float)
    gauss = 1.0 / math.sqrt(3.0)

    for xi in (-gauss, gauss):
        for eta in (-gauss, gauss):
            b, det_jacobian = b_matrix(element_coords, xi, eta)
            stiffness += b.T @ constitutive @ b * det_jacobian

    return stiffness


def assemble_global_stiffness(nodes: np.ndarray, elements: np.ndarray, constitutive: np.ndarray) -> np.ndarray:
    total_dofs = nodes.shape[0] * 2
    stiffness = np.zeros((total_dofs, total_dofs), dtype=float)

    for element in elements:
        coords = nodes[element]
        element_stiffness = q4_element_stiffness(coords, constitutive)
        dofs = []
        for node_id in element:
            dofs.extend((2 * node_id, 2 * node_id + 1))

        stiffness[np.ix_(dofs, dofs)] += element_stiffness

    return stiffness


def evaluate_stresses(
    nodes: np.ndarray, elements: np.ndarray, constitutive: np.ndarray, displacement: np.ndarray
) -> list[np.ndarray]:
    gauss = 1.0 / math.sqrt(3.0)
    stresses = []

    for element in elements:
        coords = nodes[element]
        dofs = []
        for node_id in element:
            dofs.extend((2 * node_id, 2 * node_id + 1))
        element_displacement = displacement[dofs]

        for xi in (-gauss, gauss):
            for eta in (-gauss, gauss):
                b, _ = b_matrix(coords, xi, eta)
                strain = b @ element_displacement
                stresses.append(constitutive @ strain)

    return stresses


def main() -> None:
    constitutive = constitutive_matrix_plane_stress(E, NU)

    nodes = np.array(
        [
            [0.0, 0.0],
            [1.2, 0.1],
            [2.0, 0.0],
            [-0.1, 0.9],
            [1.0, 1.0],
            [2.1, 1.1],
            [0.0, 2.0],
            [0.9, 1.8],
            [2.0, 2.0],
        ],
        dtype=float,
    )

    elements = np.array(
        [
            [0, 1, 4, 3],
            [1, 2, 5, 4],
            [3, 4, 7, 6],
            [4, 5, 8, 7],
        ],
        dtype=int,
    )

    boundary_nodes = np.array([0, 1, 2, 3, 5, 6, 7, 8], dtype=int)
    free_nodes = np.array([4], dtype=int)

    total_dofs = nodes.shape[0] * 2
    exact_nodal_displacement = np.zeros(total_dofs, dtype=float)
    for node_id, (x_coord, y_coord) in enumerate(nodes):
        exact_nodal_displacement[2 * node_id : 2 * node_id + 2] = exact_displacement(x_coord, y_coord)

    stiffness = assemble_global_stiffness(nodes, elements, constitutive)

    prescribed_dofs = np.sort(np.concatenate((2 * boundary_nodes, 2 * boundary_nodes + 1)))
    free_dofs = np.sort(np.concatenate((2 * free_nodes, 2 * free_nodes + 1)))

    displacement = np.zeros(total_dofs, dtype=float)
    displacement[prescribed_dofs] = exact_nodal_displacement[prescribed_dofs]

    stiffness_ff = stiffness[np.ix_(free_dofs, free_dofs)]
    stiffness_fc = stiffness[np.ix_(free_dofs, prescribed_dofs)]
    rhs_free = -stiffness_fc @ displacement[prescribed_dofs]
    displacement[free_dofs] = np.linalg.solve(stiffness_ff, rhs_free)

    solved_center = displacement[free_dofs]
    exact_center = exact_nodal_displacement[free_dofs]
    free_residual = stiffness_ff @ solved_center + stiffness_fc @ displacement[prescribed_dofs]

    expected_stress = constitutive @ exact_strain()
    stresses = np.array(evaluate_stresses(nodes, elements, constitutive, displacement))
    max_stress_error = np.max(np.abs(stresses - expected_stress))
    max_center_error = np.max(np.abs(solved_center - exact_center))
    max_free_residual = np.max(np.abs(free_residual))

    print("Q4 patch test on a distorted 2x2 mesh")
    print(f"E = {E:.1f}, nu = {NU:.1f}")
    print(f"Exact strain = {exact_strain()}")
    print(f"Expected stress = {expected_stress}")
    print(f"Solved center displacement = {solved_center}")
    print(f"Exact center displacement = {exact_center}")
    print(f"Max center displacement error = {max_center_error:.6e}")
    print(f"Max Gauss-point stress error = {max_stress_error:.6e}")
    print(f"Max free-DOF residual = {max_free_residual:.6e}")

    if max_center_error < 1.0e-12 and max_stress_error < 1.0e-11 and max_free_residual < 1.0e-11:
        print("Patch test passed.")
    else:
        raise SystemExit("Patch test failed.")


if __name__ == "__main__":
    main()