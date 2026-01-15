import ThermoCompass

tc = ThermoCompass.ThermoCompass(new_kinetic_model_path="kinetics/Polimi_2601") # path relative to the main folder
tc.load_atomic_matrix() 