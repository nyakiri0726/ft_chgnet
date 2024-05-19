import numpy as np
from ase import Atoms
from ase.build import fcc100, molecule
from ase.io import write
from ase.visualize import view
from scipy.spatial.transform import Rotation


def validate_positive_integer(value: int, name: str) -> None:
    """
    Validate that a value is a positive integer.

    Parameters:
        value (int): The value to validate.
        name (str): The name of the value.

    Raises:
        ValueError: If the value is not a positive integer.
    """
    if not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer.")


def validate_non_negative(value: float, name: str) -> None:
    """
    Validate that a value is non-negative.

    Parameters:
        value (float): The value to validate.
        name (str): The name of the value.

    Raises:
        ValueError: If the value is negative.
    """
    if value < 0:
        raise ValueError(f"{name} must be non-negative.")


def validate_atoms_object(obj: Atoms, name: str) -> None:
    """
    Validate that an object is an ASE Atoms object.

    Parameters:
        obj (Atoms): The object to validate.
        name (str): The name of the object.

    Raises:
        TypeError: If the object is not an ASE Atoms object.
    """
    if not isinstance(obj, Atoms):
        raise TypeError(f"{name} must be an ASE Atoms object.")


class RandomMoleculeStructure:
    def __init__(
        self,
        molecule_templates: list,
        num_molecules_per_type: list,
        cell_size: list,
        max_attempts: int = 10000,
        overlap_threshold: float = 2.0,
    ) -> None:
        """
        Initialize the RandomMoleculeStructure.

        Parameters:
            molecule_templates (list): List of ASE Atoms objects representing molecule templates.
            num_molecules_per_type (list): List of integers representing the number of molecules of each type.
            cell_size (list): List of three non-negative numbers representing the cell size in x, y, and z dimensions.
            max_attempts (int): Maximum number of attempts to find non-overlapping configurations. Default is 10000.
            overlap_threshold (float): Threshold for overlap detection. Default is 2.0.

        Raises:
            ValueError: If the input parameters are invalid.
        """
        self.molecule_templates = molecule_templates
        self.num_molecules_per_type = num_molecules_per_type
        self.cell_size = cell_size
        self.max_attempts = max_attempts
        self.overlap_threshold = overlap_threshold
        self.positions = np.empty((0, 3))

        # Input validation
        validate_positive_integer(max_attempts, "max_attempts")
        validate_non_negative(overlap_threshold, "overlap_threshold")

        for num, name in zip(
            num_molecules_per_type, ["num_molecules_per_type[{}]".format(i) for i in range(len(num_molecules_per_type))]
        ):
            validate_positive_integer(num, name)

        if len(cell_size) != 3 or not all(isinstance(x, (int, float)) and x >= 0 for x in cell_size):
            raise ValueError("cell_size must be a list/tuple of three non-negative numbers.")

    def random_rotation_matrix(self) -> np.ndarray:
        """
        Generate a random rotation matrix.

        Returns:
            np.ndarray: Random rotation matrix.
        """
        random_rotation = Rotation.from_rotvec(np.random.rand(3) * 2 * np.pi)
        return random_rotation.as_matrix()

    def check_overlap(self, new_pos: np.ndarray) -> bool:
        """
        Check if there is overlap between molecules.

        Parameters:
            new_pos (np.ndarray): New position to check for overlap.

        Returns:
            bool: True if overlap is detected, False otherwise.
        """
        if self.positions.size == 0:
            return False
        dists = np.linalg.norm(self.positions - new_pos, axis=1)
        return np.any(dists < self.overlap_threshold)

    def generate_structure(self) -> Atoms:
        """
        Generate the molecular structure.

        Returns:
            Atoms: raodom molecular structure.
        """
        for mol_index, mol_template in enumerate(self.molecule_templates):
            num_molecules = self.num_molecules_per_type[mol_index]
            if num_molecules == 1:  # Place at the center
                self.positions = np.vstack(
                    [self.positions, np.array([self.cell_size[0] / 2, self.cell_size[1] / 2, 0.0])]
                )
            else:
                attempts = 0
                while len(self.positions) < sum(self.num_molecules_per_type[: mol_index + 1]):
                    new_pos = np.random.rand(3) * self.cell_size
                    if not self.check_overlap(new_pos):
                        self.positions = np.vstack([self.positions, new_pos])
                    else:
                        attempts += 1
                        if attempts >= self.max_attempts:
                            raise RuntimeError(
                                "Failed to find non-overlapping configuration after {} attempts.".format(
                                    self.max_attempts
                                )
                            )

        generated_structure = Atoms(cell=self.cell_size)
        for i, pos in enumerate(self.positions):
            mol_index = 0
            while i >= sum(self.num_molecules_per_type[: mol_index + 1]):
                mol_index += 1
            mol_template = self.molecule_templates[mol_index]
            mol = mol_template.copy()
            rotation_matrix = self.random_rotation_matrix()
            mol.positions = np.dot(mol.positions, rotation_matrix.T)
            mol.translate(pos)
            generated_structure += mol

        generated_structure.set_pbc([True, True, True])
        return generated_structure


class MolecularStructureManager:
    def __init__(
        self,
        slab: Atoms,
        molecule_templates: list,
        num_molecules_per_type: list,
        max_attempts: int = 10000,
        overlap_threshold: float = 2.0,
        cell_size_z: float = 30.0,
    ) -> None:
        """
        Initialize the MolecularStructureManager.

        Parameters:
            slab (Atoms): ASE Atoms object representing the slab.
            molecule_templates (list): List of ASE Atoms objects representing molecule templates.
            num_molecules_per_type (list): List of integers representing the number of molecules of each type.
            max_attempts (int): Maximum number of attempts to find non-overlapping configurations. Default is 10000.
            overlap_threshold (float): Threshold for overlap detection. Default is 2.0.
            cell_size_z (float): Cell size of random structure. Default is 30.0
        """
        self.slab = slab
        self.combined_structure = None
        slab_cell = slab.get_cell()
        cell_size = [
            slab_cell[0, 0],
            slab_cell[1, 1],
            cell_size_z,
        ]  # Use slab cell dimensions for x and y, set z manually
        self.random_molecule_structure = RandomMoleculeStructure(
            molecule_templates, num_molecules_per_type, cell_size, max_attempts, overlap_threshold
        )

    def merge_structures(self, distance_from_slab: int = 5, vacuum: int = 40) -> None:
        """
        Merge the molecular structure with the slab.

        Parameters:
            distance_from_slab (int): Distance from the slab to place the molecular structure. Default is 5.
            vacuum (int): Vacuum space to add above the molecular structure. Default is 40.
        """
        generated_structure = self.random_molecule_structure.generate_structure()
        slab_cell = self.slab.get_cell()
        slab_center = (
            np.mean(self.slab.positions, axis=0)
            + generated_structure.cell[2] / 2
            + slab_cell[2] / 2
            + distance_from_slab
        )
        generated_structure.translate(slab_center - np.mean(generated_structure.positions, axis=0))
        self.combined_structure = self.slab + generated_structure
        self.combined_structure.set_cell(
            (slab_cell[0], slab_cell[1], slab_cell[2] + generated_structure.cell[2] + [0, 0, vacuum])
        )
        self.combined_structure.wrap()

    def save_structure(self, filename: str) -> None:
        """
        Save the combined structure to a file.

        Parameters:
            filename (str): Name of the file to save the structure to.

        Raises:
            RuntimeError: If the combined structure has not been generated.
        """
        if self.combined_structure is not None:
            write(filename, self.combined_structure)
        else:
            raise RuntimeError("Combined structure has not been generated. Call merge_structures() first.")

    def get_structure(self) -> Atoms:
        """
        Get the combined structure.

        Returns:
            Atoms: Combined structure.

        Raises:
            RuntimeError: If the combined structure has not been generated.
        """
        if self.combined_structure is not None:
            return self.combined_structure
        else:
            raise RuntimeError("Combined structure has not been generated. Call merge_structures() first.")


if __name__ == "__main__":
    # Example usage
    pt100_atoms = fcc100("Au", size=(10, 10, 5), vacuum=0.0)
    num_molecules_per_type = [1, 800, 200, 200]
    water_molecule = molecule("H2O")
    carbon_dioxide_molecule = molecule("CO2")
    methane_molecule = molecule("CH4")
    center_mol = molecule("CH3CH2NH2")
    molecule_templates = [center_mol, water_molecule, carbon_dioxide_molecule, methane_molecule]

    manager = MolecularStructureManager(pt100_atoms, molecule_templates, num_molecules_per_type)
    manager.merge_structures(distance_from_slab=5, vacuum=30)
    combined_structure = manager.get_structure()
    manager.save_structure("sample.xyz")
    view(combined_structure)
    print(len(combined_structure))
