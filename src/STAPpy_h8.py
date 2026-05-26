import argparse
import json
import math
from pathlib import Path

import numpy as np


VTK_HEXAHEDRON = 12


def constitutive_matrix_3d(e_modulus: float, poisson_ratio: float) -> np.ndarray:
    lam = e_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
    mu = e_modulus / (2.0 * (1.0 + poisson_ratio))

    constitutive = np.zeros((6, 6), dtype=float)
    constitutive[:3, :3] = lam
    np.fill_diagonal(constitutive[:3, :3], lam + 2.0 * mu)
    np.fill_diagonal(constitutive[3:, 3:], mu)
    return constitutive


def local_coordinates_h8() -> np.ndarray:
    return np.array(
        [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ],
        dtype=float,
    )


def gauss_points_2x2x2() -> list[tuple[float, float, float, float]]:
    gauss = 1.0 / math.sqrt(3.0)
    return [
        (xi, eta, zeta, 1.0)
        for xi in (-gauss, gauss)
        for eta in (-gauss, gauss)
        for zeta in (-gauss, gauss)
    ]


def shape_functions_h8(xi: float, eta: float, zeta: float) -> np.ndarray:
    shape = np.zeros(8, dtype=float)
    for node_id, (xi_i, eta_i, zeta_i) in enumerate(local_coordinates_h8()):
        shape[node_id] = 0.125 * (1.0 + xi * xi_i) * (1.0 + eta * eta_i) * (1.0 + zeta * zeta_i)
    return shape


def shape_function_gradients_h8(xi: float, eta: float, zeta: float) -> np.ndarray:
    gradients = np.zeros((8, 3), dtype=float)

    for node_id, (xi_i, eta_i, zeta_i) in enumerate(local_coordinates_h8()):
        gradients[node_id, 0] = 0.125 * xi_i * (1.0 + eta * eta_i) * (1.0 + zeta * zeta_i)
        gradients[node_id, 1] = 0.125 * eta_i * (1.0 + xi * xi_i) * (1.0 + zeta * zeta_i)
        gradients[node_id, 2] = 0.125 * zeta_i * (1.0 + xi * xi_i) * (1.0 + eta * eta_i)

    return gradients


def b_matrix_h8(element_coords: np.ndarray, xi: float, eta: float, zeta: float) -> tuple[np.ndarray, float, np.ndarray]:
    natural_gradients = shape_function_gradients_h8(xi, eta, zeta)
    jacobian = element_coords.T @ natural_gradients
    det_jacobian = np.linalg.det(jacobian)

    if det_jacobian <= 0.0:
        raise ValueError(f"Element has a non-positive Jacobian determinant: {det_jacobian}")

    cartesian_gradients = natural_gradients @ np.linalg.inv(jacobian)
    b = np.zeros((6, 24), dtype=float)

    for local_node in range(8):
        dndx, dndy, dndz = cartesian_gradients[local_node]
        column = 3 * local_node

        b[0, column] = dndx
        b[1, column + 1] = dndy
        b[2, column + 2] = dndz
        b[3, column] = dndy
        b[3, column + 1] = dndx
        b[4, column + 1] = dndz
        b[4, column + 2] = dndy
        b[5, column] = dndz
        b[5, column + 2] = dndx

    return b, det_jacobian, shape_functions_h8(xi, eta, zeta)


def h8_element_stiffness(element_coords: np.ndarray, constitutive: np.ndarray) -> np.ndarray:
    stiffness = np.zeros((24, 24), dtype=float)

    for xi, eta, zeta, weight in gauss_points_2x2x2():
        b, det_jacobian, _ = b_matrix_h8(element_coords, xi, eta, zeta)
        stiffness += b.T @ constitutive @ b * det_jacobian * weight

    return stiffness


def h8_element_body_force(element_coords: np.ndarray, body_force: np.ndarray) -> np.ndarray:
    element_force = np.zeros(24, dtype=float)

    for xi, eta, zeta, weight in gauss_points_2x2x2():
        _, det_jacobian, shape = b_matrix_h8(element_coords, xi, eta, zeta)
        n_matrix = np.zeros((3, 24), dtype=float)
        for local_node in range(8):
            column = 3 * local_node
            n_matrix[0, column] = shape[local_node]
            n_matrix[1, column + 1] = shape[local_node]
            n_matrix[2, column + 2] = shape[local_node]
        element_force += n_matrix.T @ body_force * det_jacobian * weight

    return element_force


def element_dofs(element: np.ndarray) -> list[int]:
    dofs: list[int] = []
    for node_id in element:
        dofs.extend((3 * int(node_id), 3 * int(node_id) + 1, 3 * int(node_id) + 2))
    return dofs


def assemble_global_stiffness(nodes: np.ndarray, elements: np.ndarray, constitutive: np.ndarray) -> np.ndarray:
    total_dofs = nodes.shape[0] * 3
    stiffness = np.zeros((total_dofs, total_dofs), dtype=float)

    for element in elements:
        coords = nodes[element]
        element_stiffness = h8_element_stiffness(coords, constitutive)
        dofs = element_dofs(element)
        stiffness[np.ix_(dofs, dofs)] += element_stiffness

    return stiffness


def assemble_global_force(
    nodes: np.ndarray,
    elements: np.ndarray,
    nodal_forces: np.ndarray | None = None,
    body_force: np.ndarray | None = None,
) -> np.ndarray:
    total_dofs = nodes.shape[0] * 3
    force = np.zeros(total_dofs, dtype=float)

    if nodal_forces is not None:
        force += nodal_forces

    if body_force is not None and np.linalg.norm(body_force) > 0.0:
        for element in elements:
            coords = nodes[element]
            dofs = element_dofs(element)
            force[dofs] += h8_element_body_force(coords, body_force)

    return force


def solve_linear_system(
    stiffness: np.ndarray,
    force: np.ndarray,
    prescribed_dofs: np.ndarray,
    prescribed_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    total_dofs = stiffness.shape[0]
    all_dofs = np.arange(total_dofs, dtype=int)
    free_dofs = np.setdiff1d(all_dofs, prescribed_dofs)

    displacement = np.zeros(total_dofs, dtype=float)
    displacement[prescribed_dofs] = prescribed_values

    stiffness_ff = stiffness[np.ix_(free_dofs, free_dofs)]
    stiffness_fc = stiffness[np.ix_(free_dofs, prescribed_dofs)]
    rhs_free = force[free_dofs] - stiffness_fc @ displacement[prescribed_dofs]
    displacement[free_dofs] = np.linalg.solve(stiffness_ff, rhs_free)

    residual_free = stiffness_ff @ displacement[free_dofs] + stiffness_fc @ displacement[prescribed_dofs] - force[free_dofs]
    reactions = stiffness @ displacement - force
    return displacement, free_dofs, residual_free, reactions


def element_gauss_stress_data(
    nodes: np.ndarray, elements: np.ndarray, constitutive: np.ndarray, displacement: np.ndarray
) -> list[dict[str, np.ndarray]]:
    results: list[dict[str, np.ndarray]] = []

    for element_id, element in enumerate(elements):
        coords = nodes[element]
        dofs = element_dofs(element)
        element_displacement = displacement[dofs]

        for xi, eta, zeta, _ in gauss_points_2x2x2():
            b, _, shape = b_matrix_h8(coords, xi, eta, zeta)
            strain = b @ element_displacement
            stress = constitutive @ strain
            physical_point = shape @ coords
            results.append(
                {
                    "element_id": np.array([element_id], dtype=int),
                    "point": physical_point,
                    "strain": strain,
                    "stress": stress,
                }
            )

    return results


def average_element_stress(
    nodes: np.ndarray, elements: np.ndarray, constitutive: np.ndarray, displacement: np.ndarray
) -> np.ndarray:
    cell_stress = np.zeros((elements.shape[0], 6), dtype=float)

    for element_id, element in enumerate(elements):
        coords = nodes[element]
        dofs = element_dofs(element)
        element_displacement = displacement[dofs]
        stress_sum = np.zeros(6, dtype=float)

        for xi, eta, zeta, _ in gauss_points_2x2x2():
            b, _, _ = b_matrix_h8(coords, xi, eta, zeta)
            strain = b @ element_displacement
            stress_sum += constitutive @ strain

        cell_stress[element_id] = stress_sum / 8.0

    return cell_stress


def write_vtk_legacy(
    file_path: str | Path,
    nodes: np.ndarray,
    elements: np.ndarray,
    displacement: np.ndarray,
    cell_stress: np.ndarray,
) -> None:
    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    with file_path.open("w", encoding="ascii") as vtk_file:
        vtk_file.write("# vtk DataFile Version 3.0\n")
        vtk_file.write("STAPpy H8 results\n")
        vtk_file.write("ASCII\n")
        vtk_file.write("DATASET UNSTRUCTURED_GRID\n")

        vtk_file.write(f"POINTS {nodes.shape[0]} float\n")
        for point in nodes:
            vtk_file.write(f"{point[0]:.10e} {point[1]:.10e} {point[2]:.10e}\n")

        total_cell_entries = elements.shape[0] * 9
        vtk_file.write(f"CELLS {elements.shape[0]} {total_cell_entries}\n")
        for element in elements:
            vtk_file.write("8 " + " ".join(str(int(node_id)) for node_id in element) + "\n")

        vtk_file.write(f"CELL_TYPES {elements.shape[0]}\n")
        for _ in range(elements.shape[0]):
            vtk_file.write(f"{VTK_HEXAHEDRON}\n")

        displacements = displacement.reshape((-1, 3))
        displacement_magnitude = np.linalg.norm(displacements, axis=1)

        vtk_file.write(f"POINT_DATA {nodes.shape[0]}\n")
        vtk_file.write("VECTORS displacement float\n")
        for vector in displacements:
            vtk_file.write(f"{vector[0]:.10e} {vector[1]:.10e} {vector[2]:.10e}\n")

        vtk_file.write("SCALARS displacement_magnitude float 1\n")
        vtk_file.write("LOOKUP_TABLE default\n")
        for value in displacement_magnitude:
            vtk_file.write(f"{value:.10e}\n")

        vtk_file.write(f"CELL_DATA {elements.shape[0]}\n")
        stress_names = ("sigma_xx", "sigma_yy", "sigma_zz", "tau_xy", "tau_yz", "tau_xz")
        for component_id, stress_name in enumerate(stress_names):
            vtk_file.write(f"SCALARS {stress_name} float 1\n")
            vtk_file.write("LOOKUP_TABLE default\n")
            for value in cell_stress[:, component_id]:
                vtk_file.write(f"{value:.10e}\n")


def write_summary(
    file_path: str | Path,
    title: str,
    nodes: np.ndarray,
    displacement: np.ndarray,
    reactions: np.ndarray,
    cell_stress: np.ndarray,
) -> None:
    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    with file_path.open("w", encoding="utf-8") as summary_file:
        summary_file.write(f"{title}\n\n")
        summary_file.write("Nodal displacements\n")
        summary_file.write("node, x, y, z, ux, uy, uz, rx, ry, rz\n")
        reshaped_displacement = displacement.reshape((-1, 3))
        reshaped_reactions = reactions.reshape((-1, 3))
        for node_id, point in enumerate(nodes, start=1):
            disp = reshaped_displacement[node_id - 1]
            reac = reshaped_reactions[node_id - 1]
            summary_file.write(
                f"{node_id}, {point[0]:.10e}, {point[1]:.10e}, {point[2]:.10e}, "
                f"{disp[0]:.10e}, {disp[1]:.10e}, {disp[2]:.10e}, "
                f"{reac[0]:.10e}, {reac[1]:.10e}, {reac[2]:.10e}\n"
            )

        summary_file.write("\nElement average stress\n")
        summary_file.write("element, sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz\n")
        for element_id, stress in enumerate(cell_stress, start=1):
            summary_file.write(
                f"{element_id}, {stress[0]:.10e}, {stress[1]:.10e}, {stress[2]:.10e}, "
                f"{stress[3]:.10e}, {stress[4]:.10e}, {stress[5]:.10e}\n"
            )


def linear_exact_displacement(point: np.ndarray) -> np.ndarray:
    x_coord, y_coord, z_coord = point
    return np.array(
        [
            0.010 * x_coord + 0.002 * y_coord - 0.001 * z_coord,
            -0.003 * x_coord + 0.015 * y_coord + 0.004 * z_coord,
            0.002 * x_coord - 0.005 * y_coord + 0.020 * z_coord,
        ],
        dtype=float,
    )


def linear_exact_strain() -> np.ndarray:
    return np.array([0.010, 0.015, 0.020, -0.001, -0.001, 0.001], dtype=float)


def quadratic_exact_displacement(point: np.ndarray) -> np.ndarray:
    x_coord, y_coord, z_coord = point
    return np.array([0.02 * x_coord * x_coord, -0.015 * y_coord * y_coord, 0.01 * z_coord * z_coord], dtype=float)


def quadratic_exact_strain(point: np.ndarray) -> np.ndarray:
    x_coord, y_coord, z_coord = point
    return np.array([0.04 * x_coord, -0.03 * y_coord, 0.02 * z_coord, 0.0, 0.0, 0.0], dtype=float)


def quadratic_body_force(e_modulus: float, poisson_ratio: float) -> np.ndarray:
    lam = e_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
    mu = e_modulus / (2.0 * (1.0 + poisson_ratio))
    return -np.array(
        [
            0.04 * (lam + 2.0 * mu),
            -0.03 * (lam + 2.0 * mu),
            0.02 * (lam + 2.0 * mu),
        ],
        dtype=float,
    )


def structured_block_mesh(
    divisions: tuple[int, int, int], transform: np.ndarray, offset: np.ndarray
) -> tuple[np.ndarray, np.ndarray, dict[str, list[int]]]:
    nx_div, ny_div, nz_div = divisions
    nodes = []
    grid_to_node: dict[tuple[int, int, int], int] = {}

    for k in range(nz_div + 1):
        for j in range(ny_div + 1):
            for i in range(nx_div + 1):
                xi = i / nx_div if nx_div else 0.0
                eta = j / ny_div if ny_div else 0.0
                zeta = k / nz_div if nz_div else 0.0
                point = offset + transform @ np.array([xi, eta, zeta], dtype=float)
                node_id = len(nodes)
                nodes.append(point)
                grid_to_node[(i, j, k)] = node_id

    elements = []
    for k in range(nz_div):
        for j in range(ny_div):
            for i in range(nx_div):
                elements.append(
                    [
                        grid_to_node[(i, j, k)],
                        grid_to_node[(i + 1, j, k)],
                        grid_to_node[(i + 1, j + 1, k)],
                        grid_to_node[(i, j + 1, k)],
                        grid_to_node[(i, j, k + 1)],
                        grid_to_node[(i + 1, j, k + 1)],
                        grid_to_node[(i + 1, j + 1, k + 1)],
                        grid_to_node[(i, j + 1, k + 1)],
                    ]
                )

    sets = {
        "xmin": [],
        "xmax": [],
        "ymin": [],
        "ymax": [],
        "zmin": [],
        "zmax": [],
        "boundary": [],
        "interior": [],
    }
    for (i, j, k), node_id in grid_to_node.items():
        if i == 0:
            sets["xmin"].append(node_id)
        if i == nx_div:
            sets["xmax"].append(node_id)
        if j == 0:
            sets["ymin"].append(node_id)
        if j == ny_div:
            sets["ymax"].append(node_id)
        if k == 0:
            sets["zmin"].append(node_id)
        if k == nz_div:
            sets["zmax"].append(node_id)
        if i in (0, nx_div) or j in (0, ny_div) or k in (0, nz_div):
            sets["boundary"].append(node_id)
        else:
            sets["interior"].append(node_id)

    return np.array(nodes, dtype=float), np.array(elements, dtype=int), sets


def distorted_patch_mesh(divisions: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray, dict[str, list[int]]]:
    return structured_block_mesh(
        divisions,
        np.array([[1.30, 0.20, 0.10], [0.00, 1.10, 0.16], [0.08, 0.12, 0.90]], dtype=float),
        np.array([-0.6, -0.73, -0.25], dtype=float),
    )


def prescribed_from_node_values(node_values: list[dict[str, object]], node_count: int) -> tuple[np.ndarray, np.ndarray]:
    prescribed: dict[int, float] = {}

    for item in node_values:
        node_id = int(item["node"]) - 1
        values = item["values"]
        if not isinstance(values, list) or len(values) != 3:
            raise ValueError("Each prescribed displacement entry must have three values.")
        for local_dof, value in enumerate(values):
            if value is None:
                continue
            global_dof = 3 * node_id + local_dof
            prescribed[global_dof] = float(value)

    if not prescribed:
        return np.array([], dtype=int), np.array([], dtype=float)

    prescribed_dofs = np.array(sorted(prescribed.keys()), dtype=int)
    prescribed_values = np.array([prescribed[dof] for dof in prescribed_dofs], dtype=float)
    if np.any(prescribed_dofs < 0) or np.any(prescribed_dofs >= 3 * node_count):
        raise ValueError("Prescribed DOF index is out of range.")
    return prescribed_dofs, prescribed_values


def nodal_force_vector(node_values: list[dict[str, object]], node_count: int) -> np.ndarray:
    force = np.zeros(node_count * 3, dtype=float)
    for item in node_values:
        node_id = int(item["node"]) - 1
        values = item["values"]
        if not isinstance(values, list) or len(values) != 3:
            raise ValueError("Each nodal load entry must have three values.")
        for local_dof, value in enumerate(values):
            force[3 * node_id + local_dof] += float(value)
    return force


def solve_case(case_data: dict[str, object]) -> dict[str, object]:
    material = case_data["material"]
    nodes = np.array(case_data["nodes"], dtype=float)
    elements = np.array(case_data["elements"], dtype=int) - 1
    constitutive = constitutive_matrix_3d(float(material["E"]), float(material["nu"]))

    prescribed_dofs, prescribed_values = prescribed_from_node_values(
        case_data.get("prescribed_displacements", []), nodes.shape[0]
    )
    nodal_forces = nodal_force_vector(case_data.get("nodal_loads", []), nodes.shape[0])
    body_force = None
    if "body_force" in case_data:
        body_force = np.array(case_data["body_force"], dtype=float)

    stiffness = assemble_global_stiffness(nodes, elements, constitutive)
    force = assemble_global_force(nodes, elements, nodal_forces, body_force)
    displacement, free_dofs, residual_free, reactions = solve_linear_system(
        stiffness, force, prescribed_dofs, prescribed_values
    )

    cell_stress = average_element_stress(nodes, elements, constitutive, displacement)
    return {
        "title": case_data.get("title", "STAPpy H8 solve"),
        "nodes": nodes,
        "elements": elements,
        "displacement": displacement,
        "free_dofs": free_dofs,
        "residual_free": residual_free,
        "reactions": reactions,
        "cell_stress": cell_stress,
    }


def run_patch_test() -> dict[str, float]:
    e_modulus = 1000.0
    poisson_ratio = 0.3
    constitutive = constitutive_matrix_3d(e_modulus, poisson_ratio)
    nodes, elements, sets = distorted_patch_mesh((2, 2, 2))

    total_dofs = nodes.shape[0] * 3
    exact_nodal_displacement = np.zeros(total_dofs, dtype=float)
    for node_id, point in enumerate(nodes):
        exact_nodal_displacement[3 * node_id : 3 * node_id + 3] = linear_exact_displacement(point)

    boundary_nodes = np.array(sets["boundary"], dtype=int)
    interior_nodes = np.array(sets["interior"], dtype=int)
    prescribed_dofs = np.sort(np.concatenate((3 * boundary_nodes, 3 * boundary_nodes + 1, 3 * boundary_nodes + 2)))
    free_dofs = np.sort(np.concatenate((3 * interior_nodes, 3 * interior_nodes + 1, 3 * interior_nodes + 2)))

    stiffness = assemble_global_stiffness(nodes, elements, constitutive)
    force = np.zeros(total_dofs, dtype=float)
    displacement, solved_free_dofs, residual_free, reactions = solve_linear_system(
        stiffness, force, prescribed_dofs, exact_nodal_displacement[prescribed_dofs]
    )

    gauss_data = element_gauss_stress_data(nodes, elements, constitutive, displacement)
    expected_stress = constitutive @ linear_exact_strain()
    stresses = np.array([item["stress"] for item in gauss_data], dtype=float)

    max_center_error = float(np.max(np.abs(displacement[solved_free_dofs] - exact_nodal_displacement[free_dofs])))
    max_stress_error = float(np.max(np.abs(stresses - expected_stress)))
    max_free_residual = float(np.max(np.abs(residual_free)))
    max_reaction = float(np.max(np.abs(reactions[prescribed_dofs])))

    print("H8 patch test on a distorted 2x2x2 mesh")
    print(f"E = {e_modulus:.1f}, nu = {poisson_ratio:.1f}")
    print(f"Number of nodes = {nodes.shape[0]}, number of elements = {elements.shape[0]}")
    print(f"Exact strain = {linear_exact_strain()}")
    print(f"Expected stress = {expected_stress}")
    print(f"Solved center displacement = {displacement[solved_free_dofs]}")
    print(f"Exact center displacement = {exact_nodal_displacement[free_dofs]}")
    print(f"Max center displacement error = {max_center_error:.6e}")
    print(f"Max Gauss-point stress error = {max_stress_error:.6e}")
    print(f"Max free-DOF residual = {max_free_residual:.6e}")
    print(f"Max prescribed-DOF reaction = {max_reaction:.6e}")

    if max_center_error < 1.0e-12 and max_stress_error < 1.0e-11 and max_free_residual < 1.0e-11:
        print("H8 patch test passed.")
    else:
        raise SystemExit("H8 patch test failed.")

    return {
        "max_center_error": max_center_error,
        "max_stress_error": max_stress_error,
        "max_free_residual": max_free_residual,
        "max_reaction": max_reaction,
    }


def convergence_metrics(divisions: int, e_modulus: float = 1000.0, poisson_ratio: float = 0.3) -> dict[str, float]:
    constitutive = constitutive_matrix_3d(e_modulus, poisson_ratio)
    nodes, elements, sets = structured_block_mesh(
        (divisions, divisions, divisions),
        np.array([[1.0, 0.20, 0.10], [0.00, 0.90, 0.10], [0.05, 0.10, 0.80]], dtype=float),
        np.array([0.0, 0.0, 0.0], dtype=float),
    )

    total_dofs = nodes.shape[0] * 3
    exact_nodal_displacement = np.zeros(total_dofs, dtype=float)
    for node_id, point in enumerate(nodes):
        exact_nodal_displacement[3 * node_id : 3 * node_id + 3] = quadratic_exact_displacement(point)

    boundary_nodes = np.array(sets["boundary"], dtype=int)
    interior_nodes = np.array(sets["interior"], dtype=int)
    prescribed_dofs = np.sort(np.concatenate((3 * boundary_nodes, 3 * boundary_nodes + 1, 3 * boundary_nodes + 2)))
    interior_dofs = np.sort(np.concatenate((3 * interior_nodes, 3 * interior_nodes + 1, 3 * interior_nodes + 2)))

    stiffness = assemble_global_stiffness(nodes, elements, constitutive)
    force = assemble_global_force(nodes, elements, np.zeros(total_dofs, dtype=float), quadratic_body_force(e_modulus, poisson_ratio))
    displacement, _, residual_free, _ = solve_linear_system(
        stiffness, force, prescribed_dofs, exact_nodal_displacement[prescribed_dofs]
    )

    gauss_data = element_gauss_stress_data(nodes, elements, constitutive, displacement)
    gauss_point_stress_error = []
    for item in gauss_data:
        point = np.asarray(item["point"], dtype=float)
        exact_stress = constitutive @ quadratic_exact_strain(point)
        gauss_point_stress_error.append(np.max(np.abs(np.asarray(item["stress"], dtype=float) - exact_stress)))

    interior_error = 0.0
    if interior_dofs.size > 0:
        interior_error = float(np.max(np.abs(displacement[interior_dofs] - exact_nodal_displacement[interior_dofs])))

    max_free_residual = 0.0
    if residual_free.size > 0:
        max_free_residual = float(np.max(np.abs(residual_free)))

    return {
        "divisions": float(divisions),
        "elements": float(elements.shape[0]),
        "nodes": float(nodes.shape[0]),
        "max_interior_displacement_error": interior_error,
        "max_gauss_stress_error": float(np.max(gauss_point_stress_error)),
        "max_free_residual": max_free_residual,
    }


def run_convergence(levels: list[int], csv_path: str | None = None) -> list[dict[str, float]]:
    results = [convergence_metrics(level) for level in levels]

    print("H8 convergence study with quadratic manufactured solution")
    print("divisions,elements,nodes,max_interior_displacement_error,max_gauss_stress_error,max_free_residual")
    for item in results:
        print(
            f"{int(item['divisions'])},{int(item['elements'])},{int(item['nodes'])},"
            f"{item['max_interior_displacement_error']:.6e},{item['max_gauss_stress_error']:.6e},"
            f"{item['max_free_residual']:.6e}"
        )

    if csv_path is not None:
        csv_file = Path(csv_path)
        csv_file.parent.mkdir(parents=True, exist_ok=True)
        with csv_file.open("w", encoding="ascii") as output:
            output.write(
                "divisions,elements,nodes,max_interior_displacement_error,max_gauss_stress_error,max_free_residual\n"
            )
            for item in results:
                output.write(
                    f"{int(item['divisions'])},{int(item['elements'])},{int(item['nodes'])},"
                    f"{item['max_interior_displacement_error']:.10e},{item['max_gauss_stress_error']:.10e},"
                    f"{item['max_free_residual']:.10e}\n"
                )

    return results


def solve_from_input(input_path: str | Path) -> dict[str, object]:
    with Path(input_path).open("r", encoding="utf-8") as input_file:
        case_data = json.load(input_file)

    result = solve_case(case_data)
    output_data = case_data.get("output", {})
    if "vtk" in output_data:
        write_vtk_legacy(
            output_data["vtk"],
            result["nodes"],
            result["elements"],
            result["displacement"],
            result["cell_stress"],
        )
    if "summary" in output_data:
        write_summary(
            output_data["summary"],
            str(result["title"]),
            result["nodes"],
            result["displacement"],
            result["reactions"],
            result["cell_stress"],
        )

    max_displacement = float(np.max(np.abs(result["displacement"])))
    max_residual = float(np.max(np.abs(result["residual_free"]))) if result["residual_free"].size > 0 else 0.0
    print(str(result["title"]))
    print(f"Nodes = {result['nodes'].shape[0]}, elements = {result['elements'].shape[0]}")
    print(f"Max displacement magnitude component = {max_displacement:.6e}")
    print(f"Max free-DOF residual = {max_residual:.6e}")
    if "vtk" in output_data:
        print(f"VTK written to {output_data['vtk']}")
    if "summary" in output_data:
        print(f"Summary written to {output_data['summary']}")
    return result


def write_sample_input(output_path: str | Path) -> None:
    sample = {
        "title": "Single H8 cantilever block",
        "material": {"E": 210000000000.0, "nu": 0.3},
        "nodes": [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.2, 0.0],
            [0.0, 0.2, 0.0],
            [0.0, 0.0, 0.2],
            [1.0, 0.0, 0.2],
            [1.0, 0.2, 0.2],
            [0.0, 0.2, 0.2],
        ],
        "elements": [[1, 2, 3, 4, 5, 6, 7, 8]],
        "prescribed_displacements": [
            {"node": 1, "values": [0.0, 0.0, 0.0]},
            {"node": 4, "values": [0.0, 0.0, 0.0]},
            {"node": 5, "values": [0.0, 0.0, 0.0]},
            {"node": 8, "values": [0.0, 0.0, 0.0]},
        ],
        "nodal_loads": [
            {"node": 2, "values": [0.0, 0.0, -250.0]},
            {"node": 3, "values": [0.0, 0.0, -250.0]},
            {"node": 6, "values": [0.0, 0.0, -250.0]},
            {"node": 7, "values": [0.0, 0.0, -250.0]},
        ],
        "output": {
            "vtk": "outputs/h8_cantilever.vtk",
            "summary": "outputs/h8_cantilever_summary.txt",
        },
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(sample, indent=2), encoding="utf-8")
    print(f"Sample input written to {output_path}")


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="STAPpy H8 finite element solver")
    subparsers = parser.add_subparsers(dest="command")

    subparsers.add_parser("patch-test", help="Run the H8 patch test")

    solve_parser = subparsers.add_parser("solve", help="Solve an input JSON case")
    solve_parser.add_argument("input", help="Path to the case JSON file")

    convergence_parser = subparsers.add_parser("convergence", help="Run the H8 convergence study")
    convergence_parser.add_argument("--levels", nargs="+", type=int, default=[1, 2, 4], help="Mesh divisions per axis")
    convergence_parser.add_argument("--csv", default="outputs/h8_convergence.csv", help="CSV output path")

    sample_parser = subparsers.add_parser("write-sample", help="Write a sample JSON input file")
    sample_parser.add_argument("output", help="Path to the sample JSON file")

    return parser


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()

    if args.command in (None, "patch-test"):
        run_patch_test()
        return

    if args.command == "solve":
        solve_from_input(args.input)
        return

    if args.command == "convergence":
        run_convergence(args.levels, args.csv)
        return

    if args.command == "write-sample":
        write_sample_input(args.output)
        return

    parser.error(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()