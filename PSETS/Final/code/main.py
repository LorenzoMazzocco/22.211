from lib.models import *


# ##################################################################
# #                   POINT 1 - COMPARE SALTS                      #
# ##################################################################

# ##################### PERTURBATION FUNCTION ######################

# def perturb_mod_density(hh, case_salt, case_filepath):
#     # reference case
#     model("output/1_choose_salts/{}/k_ref.csv".format(case_filepath), coolant=case_salt)

#     # perturb
#     for h in hh:
#         model("output/1_choose_salts/{}/mod_rho_+{:.2f}.csv".format(case_filepath,h), coolant=case_salt, mod_density=1.65+h)
#         model("output/1_choose_salts/{}/mod_rho_-{:.2f}.csv".format(case_filepath,h), coolant=case_salt, mod_density=1.65-h)


# hh = np.array([0.2, 0.1, 0.05])*1.65 # perturbation of 20%, 10% and 5%

# # Case 1 - NaF
# perturb_mod_density(hh=hh, case_salt="NaF", case_filepath="1_NaF")
# # Case 2 - FLiBe Natural
# perturb_mod_density(hh=hh, case_salt="FLiBe_nat", case_filepath="2_FLiBe_nat")
# # Case 3 - FLiBe Enr 1
# perturb_mod_density(hh=hh, case_salt="FLiBe_enr1", case_filepath="3_FLiBe_enr1")
# # Case 4 - FLiBe Enr 2
# perturb_mod_density(hh=hh, case_salt="FLiBe_enr2", case_filepath="4_FLiBe_enr2")



##################################################################
#                   POINT 2 - COMPARE SALTS                      #
##################################################################

# produce reference case data, in this case tallies are on and we perform a PHOTON TRANSPORT SIMULATION (this doesn't help with normalization BUT will help a lot with heat deposition in different regions, heating-local deposits all the photon heat in the place they were generated)
model("output/2_NaF/k_ref.csv", coolant="NaF", statepoint_filepath="output/2_NaF/statepoint_reference.h5", activate_tallies=True, simulate_photons=False)

#model("output/2_NaF/k_ref.csv", coolant="NaF", statepoint_filepath="output/2_NaF/statepoint_photons", activate_tallies=True, simulate_photons=True, inactive=, batches=, particles=)
