# Set matplotlib backend to 'Agg' for non-GUI environments
import matplotlib

matplotlib.use("Agg")

import os
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Project imports
from utils import from_comp_dict_to_str, replace_line, convert_species_name

# Define molecular weights for common elements (in g/mol)
C_MW = 12.01099968
H_MW = 1.00800002
O_MW = 15.99899960
N_MW = 14.00699997
K_MW = 39.10200000

PATH_TO_BIOSMOKE = "/root/OpenSMOKEppSolvers/build/OpenSMOKEpp_BioSMOKE.sh"


class ThermoCompass:
    """Class to handle thermodynamics refitting and simulation runs."""

    def __init__(self, kinetic_model_path: str = None):
        """Initialize ThermoCompass with paths to kinetic models."""
        self.gas_species_list = []  # List to hold gas species names
        self.solid_species_list = []  # List to hold solid species names

        self.gas_kinetics_tree = None  # XML tree for gas kinetics
        self.solid_kinetics_tree = None  # XML tree for solid kinetics

        self.loaded_kinetics = False

        if kinetic_model_path is not None:
            self.load_kinetic_model(kinetic_model_path)

    def load_kinetic_model(self, model_path: str):
        """Load kinetic model from the specified folder."""
        print(f"Loading kinetic model from: {model_path}...")

        # Load the gas species kinetic model
        self.gas_kinetics_tree = ET.parse(f"{model_path}/kinetics-compiled/kinetics.xml").getroot()
        self.solid_kinetics_tree = ET.parse(
            f"{model_path}/kinetics-compiled/kinetics.solid.xml"
        ).getroot()

        self.loaded_kinetics = True

        # Extract species and solid species lists
        self.species_list = self.solid_kinetics_tree.find("NamesOfSpecies").text.split()

        # Extract solid species list only
        self.solid_species_list = self.solid_kinetics_tree.find("SolidSpecies").text.split()

        # Derive gas species list: all species minus solid species
        self.gas_species_list = [
            species for species in self.species_list if species not in self.solid_species_list
        ]

        # Extract elements names
        elements_name = self.solid_kinetics_tree.find("NamesOfElements").text.split()

        # Load the actual atomic matrix
        vector_atomic_matrix = np.array(
            self.solid_kinetics_tree.find("AtomicComposition").text.split(), dtype=float
        )
        vector_atomic_matrix = vector_atomic_matrix.reshape(
            (len(self.species_list), len(elements_name))
        )  # Reshape to 2D array

        # Store into a dataframe
        self.atomic_matrix = pd.DataFrame(
            vector_atomic_matrix, index=self.species_list, columns=elements_name
        )

        # Caluclate the Molecular Weights
        molecular_weights_vector = (
            self.atomic_matrix.get("C", pd.Series(0)) * C_MW
            + self.atomic_matrix.get("H", pd.Series(0)) * H_MW
            + self.atomic_matrix.get("O", pd.Series(0)) * O_MW
            + self.atomic_matrix.get("N", pd.Series(0)) * N_MW
            + self.atomic_matrix.get("K", pd.Series(0)) * K_MW
        )

        # Append Molecular Weights to the atomic matrix
        self.atomic_matrix["Molecular_Weight"] = molecular_weights_vector

        print("Kinetic model loaded successfully.\n")

    def plot_results(self, models_to_plot: list = None, title: str = "Simulation Results"):
        """Plot the simulation results."""
        print(f"\nPlotting simulation results: {title}...")
        if models_to_plot is None:  # plot all models by default
            models_to_plot = [
                d for d in os.listdir("kinetics/") if os.path.isdir(os.path.join("kinetics/", d))
            ]

        # Ax[0]: Solid Mass over time, Ax[1]: Heat Released over time
        fig, ax = plt.subplots(1, 2, figsize=(14, 6))
        for model in models_to_plot:
            output_folder = f"simulations/output-{model}"
            results_file = f"{output_folder}/Output.out"

            if not os.path.exists(results_file):
                print(f"Results file not found for model: {model}. Skipping plot.")
                continue

            # Load results
            results_df = pd.read_csv(results_file, delim_whitespace=True)

            # Plot total solid mass fraction over time
            time = results_df["t[s](1)"]
            total_solid_mass = results_df["Ms/Ms0[-](3)"]
            heat_released = results_df["Qr[W/m3](5)"]
            integral_heat_released = (
                np.trapz(heat_released * total_solid_mass, time) * 1e-3
            )  # in kJ

            ax[0].plot(time, total_solid_mass, label=model)
            ax[1].plot(
                time,
                heat_released * total_solid_mass,
                label=f"{model} (Integral: {integral_heat_released:.2e} kJ)",
            )

        ax[0].set_xlabel("Time (s)")
        ax[0].set_ylabel("Total Solid Mass")
        ax[0].set_title("Solid Mass Over Time")
        ax[0].legend()
        ax[0].grid()

        ax[1].set_xlabel("Time (s)")
        ax[1].set_ylabel("Heat Released (W)")
        ax[1].set_title("Heat Released Over Time")
        ax[1].legend()
        ax[1].grid()

        plt.suptitle(title)
        plt.savefig("plots/plot.svg")
        print("Plotting completed.\n")

    def run_simulation(
        self,
        inputfile_path: str = "simulations/input.dic",
        solid_composition: dict = None,
        plot_results: bool = False,
        models_to_run: list = None,
        keep_results: bool = False,
    ):
        """Run the same simulation for the kinetic models."""
        if models_to_run is None:  # run all models by default
            models_to_run = [
                d for d in os.listdir("kinetics/") if os.path.isdir(os.path.join("kinetics/", d))
            ]

        # set default solid composition if not provided
        if solid_composition is None:
            solid_composition = {"CELL": 1.0}
        else:
            if sum(solid_composition.values()) != 1.0:
                raise ValueError("The solid composition fractions must sum to 1.0")

        # Run simulations for each model
        for model in models_to_run:
            print(f"Running simulation for model: {model}...")
            model_path = f"kinetics/{model}"
            # replace solid composition in the input file

            solid_composition_modified = solid_composition.copy()
            # Older model have differen solid species names
            if model in ["Polimi_1402", "Polimi_1805"] and "MOIST" in solid_composition_modified:
                solid_composition_modified["ACQUA"] = solid_composition_modified.pop("MOIST")
            if model == "Polimi_1402":
                if "XYHW" in solid_composition_modified:
                    solid_composition_modified["HCE"] = solid_composition_modified.pop("XYHW")
                if "GMSW" in solid_composition_modified:
                    solid_composition_modified["HCE"] = solid_composition_modified.pop("GMSW")
                if "XYGR" in solid_composition_modified:
                    solid_composition_modified["HCE"] = solid_composition_modified.pop("XYGR")

            replace_line(
                inputfile_path,
                "@MassFractions",
                f"\t@MassFractions \t\t{from_comp_dict_to_str(solid_composition_modified)};",
            )

            # Replace the kinetic model path in the input file
            replace_line(
                inputfile_path,
                "@KineticsFolder",
                f"\t@KineticsFolder \t\t{model_path}/kinetics-compiled;",
            )

            # Replace the output folder in the input file
            replace_line(
                inputfile_path,
                "@OutputFolder",
                f"\t@OutputFolder \t\tsimulations/output-{model};",
            )

            # Run the simulation using OpenSMOKEpp BioSMOKE solver
            os.system(f"{PATH_TO_BIOSMOKE} --input {inputfile_path} > /dev/null")

        if plot_results:
            self.plot_results(
                models_to_plot=models_to_run, title=from_comp_dict_to_str(solid_composition)
            )

        if not keep_results:
            # Clean up output folders
            for model in models_to_run:
                output_folder = f"simulations/output-{model}"
                if os.path.exists(output_folder):
                    os.system(f"rm -r {output_folder}")
            print("Cleaned up output folders.\n")

    def load_nasa_coefficents(self, kinetic_model_path):
        """Load NASA coefficients from the kinetic model."""
        print("Loading NASA coefficients...")

        # Check if kinetic model is loaded
        if self.loaded_kinetics is False:
            print("Kinetic model not loaded. Loading now...")
            self.load_kinetic_model(kinetic_model_path)

        nasa_index = [
            "a1_HT",
            "a2_HT",
            "a3_HT",
            "a4_HT",
            "a5_HT",
            "a6_HT",
            "a7_HT",
            "a1_LT",
            "a2_LT",
            "a3_LT",
            "a4_LT",
            "a5_LT",
            "a6_LT",
            "a7_LT",
            "T_min",
            "T_max",
            "T_mid",
            "MW",
        ]

        # Extract NASA coefficients for all the species
        nasa_data = np.array(
            self.solid_kinetics_tree.find("Thermodynamics/NASA-coefficients").text.split(),
            dtype=float,
        )

        # generate dataframe
        self.nasa_data = pd.DataFrame(
            nasa_data.reshape((len(self.species_list), 18)),
            index=self.species_list,
            columns=nasa_index,
        )

        print("NASA Coefficients loaded successfully.\n")

    def calculate_entalpy(self, species_name: str, temperature: float):
        """Calculate the enthalpy of a species at a given temperature using NASA polynomials in kJ/mol.
        Handles conversion of species names if needed."""
        # Check if NASA coefficients are loaded
        if not hasattr(self, "nasa_data"):
            raise ValueError("NASA coefficients not loaded. Please load them first.")

        # Check if species exists in the NASA data
        if species_name not in self.species_list:
            raise ValueError(f"Species {species_name} not found in model species.")

        coeffs = self.nasa_data.loc[species_name]

        # check if temperature is within the valid range
        if temperature < coeffs["T_min"] or temperature > coeffs["T_max"]:
            raise ValueError(
                f"Temperature {temperature} K is out of range for species {species_name} "
                f"({coeffs['T_min']} K - {coeffs['T_max']} K)."
            )

        if temperature < coeffs["T_mid"]:
            a = coeffs[["a1_LT", "a2_LT", "a3_LT", "a4_LT", "a5_LT", "a6_LT", "a7_LT"]].values
        else:
            a = coeffs[["a1_HT", "a2_HT", "a3_HT", "a4_HT", "a5_HT", "a6_HT", "a7_HT"]].values

        t = temperature
        H_RT = a[0] + a[1] * t / 2 + a[2] * t**2 / 3 + a[3] * t**3 / 4 + a[4] * t**4 / 5 + a[5] / t

        R = 8.314462618  # J/(mol·K)
        enthalpy = H_RT * R * temperature * 1e-3  # in kJ/mol

        return enthalpy

    def calculate_specific_heat(self, species_name: str, temperature: float):
        """Calculate the specific heat capacity of a species at a given temperature using NASA polynomials in kJ/(mol·K).
        Handles conversion of species names if needed."""

        # Check if NASA coefficients are loaded
        if not hasattr(self, "nasa_data"):
            raise ValueError("NASA coefficients not loaded. Please load them first.")

        # Check if species exists in the NASA data
        if species_name not in self.species_list:
            raise ValueError(f"Species {species_name} not found in model species.")

        coeffs = self.nasa_data.loc[species_name]

        # check if temperature is within the valid range
        if temperature < coeffs["T_min"] or temperature > coeffs["T_max"]:
            raise ValueError(
                f"Temperature {temperature} K is out of range for species {species_name} "
                f"({coeffs['T_min']} K - {coeffs['T_max']} K)."
            )

        if temperature < coeffs["T_mid"]:
            a = coeffs[["a1_LT", "a2_LT", "a3_LT", "a4_LT", "a5_LT", "a6_LT", "a7_LT"]].values
        else:
            a = coeffs[["a1_HT", "a2_HT", "a3_HT", "a4_HT", "a5_HT", "a6_HT", "a7_HT"]].values

        t = temperature
        Cp_R = a[0] + a[1] * t + a[2] * t**2 + a[3] * t**3 + a[4] * t**4

        R = 8.314462618  # J/(mol·K)
        specific_heat = Cp_R * R * 1e-3  # in kJ/(mol·K)

        return specific_heat

    def plot_entalpy(
        self,
        species_name: str,
        T_min: float = 300,
        T_max: float = 3000,
        num_points: int = 100,
        models_to_plot: list = None,
    ):
        """Plot enthalpy vs temperature for a given species."""
        if models_to_plot is None:  # plot all models by default
            models_to_plot = [
                d for d in os.listdir("kinetics/") if os.path.isdir(os.path.join("kinetics/", d))
            ]

        print(f"\nPlotting enthalpy for species: {species_name}...")
        temperatures = np.linspace(T_min, T_max, num_points)
        plt.figure(figsize=(8, 6))
        for model in models_to_plot:
            self.load_kinetic_model(f"kinetics/{model}")

            species_name_converted = convert_species_name(species_name, model)

            self.load_nasa_coefficents(kinetic_model_path=f"kinetics/{model}")

            enthalpies = [
                self.calculate_entalpy(species_name_converted, T) for T in temperatures
            ]  # in kJ/mol
            plt.plot(temperatures, enthalpies, label=f"{model} {species_name_converted}")

        plt.xlabel("Temperature (K)")
        plt.ylabel("Enthalpy (kJ/mol)")
        plt.title(f"Enthalpy vs Temperature for {species_name}")
        plt.legend()
        plt.grid()
        plt.savefig(f"plots/enthalpy_{species_name}.svg")
        print("Plotting completed.\n")

        def plot_specific_heat(
            self,
            species_name: str,
            T_min: float = 300,
            T_max: float = 3000,
            num_points: int = 100,
            models_to_plot: list = None,
        ):
            """Plot specific heat capacity vs temperature for a given species."""
            if models_to_plot is None:  # plot all models by default
                models_to_plot = [
                    d
                    for d in os.listdir("kinetics/")
                    if os.path.isdir(os.path.join("kinetics/", d))
                ]

            print(f"\nPlotting specific heat capacity for species: {species_name}...")
            temperatures = np.linspace(T_min, T_max, num_points)
            plt.figure(figsize=(8, 6))
            for model in models_to_plot:
                self.load_kinetic_model(f"kinetics/{model}")
                self.load_nasa_coefficents(kinetic_model_path=f"kinetics/{model}")

                # Convert species name if needed
                species_name_converted = convert_species_name(species_name, model)

                specific_heats = [
                    self.calculate_specific_heat(species_name_converted, T) for T in temperatures
                ]  # in kJ/(mol·K)

                plt.plot(temperatures, specific_heats, label=f"{model} {species_name_converted}")

            plt.xlabel("Temperature (K)")
            plt.ylabel("Specific Heat Capacity (kJ/(mol·K))")
            plt.title(f"Specific Heat Capacity vs Temperature for {species_name}")
            plt.legend()
            plt.grid()
            plt.savefig(f"plots/specific_heat_{species_name}.svg")
            print("Plotting completed.\n")


if __name__ == "__main__":
    tc = ThermoCompass()
    # tc.load_atomic_matrix()

    tc.plot_entalpy(species_name="CELL")

    exit()

    # Composition provided in the latest kinetic model solid species names
    # Conversions handeled:
    # if model == "Polimi_1402":
    #   MOIST -> ACQUA, XYHW/GMSW/XYGR -> HCE
    # if model == "Polimi_1805":
    #   MOIST -> ACQUA
    tc.run_simulation(
        keep_results=False,
        # models_to_run=["Polimi_2601", "Polimi_2401"],
        solid_composition={"CHAR": 0.9, "ASH": 0.1},
        plot_results=True,
    )
