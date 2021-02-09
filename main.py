"""Spectrum Sharing Simulation"""

__author__ = 'Awab Asher'

# To vary the densities: go to simulation parameters and change it from there!
# To get graphs at each density: the simulation parameters were changed manually!
# All simulations were not done at once!
# Labelling of figures was done manually at each density!


from network import CellularNetwork

from simulation_parameters import (
    AREA_OF_SIMULATION,
    UE_DENSITY,
    BS_DENSITY_A,
    BS_DENSITY_B,
    ITERATIONS,
    TRANSMIT_POWER,
    CARRIER_FREQUENCY
)

from constants import (
    exclusive_bandwidth,
    exclusive_bandwidth_iii,
    exclusive_bandwidth_iv,
    exclusive_bandwidth_v,
    pooled_bandwidth,
    noise_power_exclusive,
    noise_power_pooled
)


from functions import (
    compute_rx_power,
    calculate_snr,
    create_ecdf,
    calculate_prelog_term,
    calculate_sinr_exclusive,
    calculate_throughput,
    calculate_sinr_pooled_ii,
    calculate_sinr_pooled_iii,
    calculate_sinr_pooled_iv,
    calculate_sinr_pooled_v,
    limit_throughput_snr,
    limit_throughput_sinr,
    compute_distance_matrix,
    compute_distance_matrix_m,
    compute_pathloss_matrix_fspl,
    compute_pathloss_matrix_nyc,
    cell_associate,
    compute_count_for_bs,
    compute_total_ue_for_bs,
    compute_distance_ue_bs
)

from pathloss import compute_pathloss_fspl

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 12, })


# Creates a cellular network object!
network = CellularNetwork(ue_density=UE_DENSITY,
                          bs_density_a=BS_DENSITY_A,
                          bs_density_b=BS_DENSITY_B,
                          area=AREA_OF_SIMULATION,
                          iterations=ITERATIONS)

# Runs the methods associated with the Cellular Network Class!
networkclass_methods = ["deploy_users",
                        "deploy_users_",
                        "deploy_basestations",
                        "deploy_basestations_",
                        "_deploy_basestations_",
                        "_deploy_basestations__",
                        "__deploy_basestations__"]

for method in networkclass_methods:
    getattr(network, method)()


users_A = network.UEs_A     # List of user objects of operator A!
users_B = network.UEs_B     # List of user objects of operator B!
basestations_A = network.BSs_A      # List of basestation objects of operator A!
basestations_B = network.BSs_B      # List of basestation objects of operator B!
basestations_C = network.BSs_C      # List of basestation objects of operator C!
basestations_D = network.BSs_D      # List of basestation objects of operator D!
basestations_E = network.BSs_E      # List of basestation objects of operator E!

# Computes the distance matrix for users in operator A in km!
distance_matrix_A = compute_distance_matrix(users=network.UEs_A, basestations=network.BSs_A)

# Computes the distance matrix for users in operator A in m!
distance_matrix_A_m = compute_distance_matrix_m(distance_matrix_km=distance_matrix_A)

# Computes the fspl matrix
pathloss_matrix_fspl_A = compute_pathloss_matrix_fspl(_dist_matrix_m=distance_matrix_A_m)

# Computes the nyc pathloss matrix for users in operator A!
nyc_matrix_A = compute_pathloss_matrix_nyc(_dist_matrix_m=distance_matrix_A_m)

# Inter-operator interference for users in operator A!
# Computes the distance matrix in km for users in operator A with basestations of operator B!
distance_matrix_A_B = compute_distance_matrix(users=network.UEs_A, basestations=network.BSs_B)

# Computes the distance matrix in m for users in operator A with basestations of operator B!
distance_matrix_A_B_m = compute_distance_matrix_m(distance_matrix_km=distance_matrix_A_B)

# Computes the pathloss matrix for users in operator A with basestations of operator B!
pathloss_matrix_fspl_A_B = compute_pathloss_matrix_fspl(_dist_matrix_m=distance_matrix_A_B_m)


# Computes the distance matrix in km for users in operator A with basestations of operator C!
distance_matrix_A_C = compute_distance_matrix(users=network.UEs_A, basestations=network.BSs_C)

# Computes the distance matrix in m for users in operator A with basestations of operator C!
distance_matrix_A_C_m = compute_distance_matrix_m(distance_matrix_km=distance_matrix_A_C)

# Computes the pathloss matrix for users in operator A with basestations of operator C!
pathloss_matrix_fspl_A_C = compute_pathloss_matrix_fspl(_dist_matrix_m=distance_matrix_A_C_m)


# Computes the distance matrix in km for users in operator A with basestations of operator D!
distance_matrix_A_D = compute_distance_matrix(users=network.UEs_A, basestations=network.BSs_D)

# Computes the distance matrix in m for users in operator A with basestations of operator D!
distance_matrix_A_D_m = compute_distance_matrix_m(distance_matrix_km=distance_matrix_A_D)

# Computes the pathloss matrix for users in operator A with basestations of operator D!
pathloss_matrix_fspl_A_D = compute_pathloss_matrix_fspl(_dist_matrix_m=distance_matrix_A_D_m)


# Computes the distance matrix in km for users in operator A with basestations of operator E!
distance_matrix_A_E = compute_distance_matrix(users=network.UEs_A, basestations=network.BSs_E)

# Computes the distance matrix in m for users in operator A with basestations of operator E!
distance_matrix_A_E_m = compute_distance_matrix_m(distance_matrix_km=distance_matrix_A_E)

# Computes the pathloss matrix for users in operator A with basestations of operator E!
pathloss_matrix_fspl_A_E = compute_pathloss_matrix_fspl(_dist_matrix_m=distance_matrix_A_E_m)


# Computes the cell association methodology used for users in operator A!
cell_associate_A = cell_associate(pathloss_matrix=pathloss_matrix_fspl_A,
                                  users=network.UEs_A,
                                  basestations=network.BSs_A)


count_for_bs_A = compute_count_for_bs(pathloss_matrix=pathloss_matrix_fspl_A, basestations=network.BSs_A)
count_for_ue_A = np.array(compute_total_ue_for_bs(_cell_associate=cell_associate_A, _bs_count=count_for_bs_A))

# Computes the distances between respective users and basestations belonging to operator A!
distances_A = np.array(compute_distance_ue_bs(_cell_associate=cell_associate_A))


# Pathlosses for each UE in operator A are calculated!
pathlosses_A = np.fromiter((compute_pathloss_fspl(distance=_,
                                                  frequency=CARRIER_FREQUENCY) for _ in distances_A),
                           distances_A.dtype)


# Received power for each UE in operator A is calculated in dB!
rxpower_A = compute_rx_power(P_tx=TRANSMIT_POWER,
                             G_tx=0,
                             G_rx=0,
                             PL=pathlosses_A)


# SNR for an exclusive licensed model in dB for operator A!
snr_exclusive_A = calculate_snr(rx_power=rxpower_A,
                                noise_power=noise_power_exclusive)


# SNR for a fully pooled case model in dB for operator A!
snr_pooled_A = calculate_snr(rx_power=rxpower_A,
                             noise_power=noise_power_pooled)


# Prelog term for an exclusive licensed model for operator A!
prelog_exclusive_A = calculate_prelog_term(bandwidth=exclusive_bandwidth,
                                           number=count_for_ue_A)

# Prelog term for a fully pooled case for operator A!
prelog_pooled_A = calculate_prelog_term(bandwidth=pooled_bandwidth,
                                        number=count_for_ue_A)


# SNR throughput for an exclusive licensed model for operator A!
snr_throughput_exclusive_A = calculate_throughput(dB_values=snr_exclusive_A,
                                                  prelog_term=prelog_exclusive_A)

# SNR throughput for a fully pooled case for operator A!
snr_throughput_pooled_A = calculate_throughput(dB_values=snr_pooled_A,
                                               prelog_term=prelog_pooled_A)


# SINR for an exclusive licensed model in dB for operator A!
sinr_dB_exclusive_A = calculate_sinr_exclusive(pathloss_matrix=pathloss_matrix_fspl_A,
                                               rx_power=rxpower_A,
                                               noise_power=noise_power_exclusive)


# SINR for a fully pooled case in dB for operator A (2 operators)!
sinr_dB_pooled_A = calculate_sinr_pooled_ii(pathloss_matrix_A=pathloss_matrix_fspl_A,
                                            pathloss_matrix_B=pathloss_matrix_fspl_A_B,
                                            rx_power=rxpower_A,
                                            noise_power=noise_power_pooled)


# SINR for a fully pooled case in dB for operator A (3 operators)!
sinr_dB_pooled_A_iii = calculate_sinr_pooled_iii(pathloss_matrix_A=pathloss_matrix_fspl_A,
                                                 pathloss_matrix_B=pathloss_matrix_fspl_A_B,
                                                 pathloss_matrix_C=pathloss_matrix_fspl_A_C,
                                                 rx_power=rxpower_A,
                                                 noise_power=noise_power_pooled)


# SINR for a fully pooled case in dB for operator A (4 operators)!
sinr_dB_pooled_A_iv = calculate_sinr_pooled_iv(pathloss_matrix_A=pathloss_matrix_fspl_A,
                                               pathloss_matrix_B=pathloss_matrix_fspl_A_B,
                                               pathloss_matrix_C=pathloss_matrix_fspl_A_C,
                                               pathloss_matrix_D=pathloss_matrix_fspl_A_D,
                                               rx_power=rxpower_A,
                                               noise_power=noise_power_pooled)

# SINR for a fully pooled case in dB for operator A (5 operators)!
sinr_dB_pooled_A_v = calculate_sinr_pooled_v(pathloss_matrix_A=pathloss_matrix_fspl_A,
                                             pathloss_matrix_B=pathloss_matrix_fspl_A_B,
                                             pathloss_matrix_C=pathloss_matrix_fspl_A_C,
                                             pathloss_matrix_D=pathloss_matrix_fspl_A_D,
                                             pathloss_matrix_E=pathloss_matrix_fspl_A_E,
                                             rx_power=rxpower_A,
                                             noise_power=noise_power_pooled)


# SINR throughput for an exclusive licensed model for operator A!
sinr_throughput_exclusive_A = calculate_throughput(dB_values=sinr_dB_exclusive_A,
                                                   prelog_term=prelog_exclusive_A)


# Prelog term for an exclusive licensed model for operator A!
prelog_exclusive_A_iii = calculate_prelog_term(bandwidth=exclusive_bandwidth_iii,
                                               number=count_for_ue_A)

prelog_exclusive_A_iv = calculate_prelog_term(bandwidth=exclusive_bandwidth_iv,
                                              number=count_for_ue_A)

prelog_exclusive_A_v = calculate_prelog_term(bandwidth=exclusive_bandwidth_v,
                                             number=count_for_ue_A)


sinr_throughput_exclusive_A_iii = calculate_throughput(dB_values=sinr_dB_exclusive_A,
                                                       prelog_term=prelog_exclusive_A_iii)

sinr_throughput_exclusive_A_iv = calculate_throughput(dB_values=sinr_dB_exclusive_A,
                                                      prelog_term=prelog_exclusive_A_iv)

sinr_throughput_exclusive_A_v = calculate_throughput(dB_values=sinr_dB_exclusive_A,
                                                     prelog_term=prelog_exclusive_A_v)

# SINR throughput for a fully pooled case for operator A (2 operators)!
sinr_throughput_pooled_A = calculate_throughput(dB_values=sinr_dB_pooled_A,
                                                prelog_term=prelog_pooled_A)


# SINR throughput for a fully pooled case for operator A (3 operators)!
sinr_throughput_pooled_A_iii = calculate_throughput(dB_values=sinr_dB_pooled_A_iii,
                                                    prelog_term=prelog_pooled_A)

# SINR throughput for a fully pooled case for operator A (4 operators)!
sinr_throughput_pooled_A_iv = calculate_throughput(dB_values=sinr_dB_pooled_A_iv,
                                                   prelog_term=prelog_pooled_A)

# SINR throughput for a fully pooled case for operator A (5 operators)!
sinr_throughput_pooled_A_v = calculate_throughput(dB_values=sinr_dB_pooled_A_v,
                                                  prelog_term=prelog_pooled_A)




# Calculations for operator 2!

# # Computes the distance matrix for users in operator B in km!
distance_matrix_B = compute_distance_matrix(users=network.UEs_B, basestations=network.BSs_B)
#
# # Computes the distance matrix for users in operator A in m!
distance_matrix_B_m = compute_distance_matrix_m(distance_matrix_km=distance_matrix_B)
#
# # Computes the fspl matrix
pathloss_matrix_fspl_B = compute_pathloss_matrix_fspl(_dist_matrix_m=distance_matrix_B_m)
#
# # Computes the nyc pathloss matrix for users in operator A!
nyc_matrix_B = compute_pathloss_matrix_nyc(_dist_matrix_m=distance_matrix_B_m)
#
# # Inter-operator interference for users in operator B!
# Computes the distance matrix in km for users in operator B with basestations of operator A!
distance_matrix_B_A = compute_distance_matrix(users=network.UEs_B, basestations=network.BSs_A)
#
# Computes the distance matrix in m for users in operator A with basestations of operator B!
distance_matrix_B_A_m = compute_distance_matrix_m(distance_matrix_km=distance_matrix_B_A)
#
# Computes the pathloss matrix for users in operator B with basestations of operator A!
pathloss_matrix_fspl_B_A = compute_pathloss_matrix_fspl(_dist_matrix_m=distance_matrix_B_A_m)
#
# # Computes the cell association methodology used for users in operator B!
cell_associate_B = cell_associate(pathloss_matrix=pathloss_matrix_fspl_B,
                                  users=network.UEs_B,
                                  basestations=network.BSs_B)
#
#
count_for_bs_B = compute_count_for_bs(pathloss_matrix=pathloss_matrix_fspl_B, basestations=network.BSs_B)
count_for_ue_B = np.array(compute_total_ue_for_bs(_cell_associate=cell_associate_B, _bs_count=count_for_bs_B))
#
# Computes the distances between respective users and basestations belonging to operator B!
distances_B = np.array(compute_distance_ue_bs(_cell_associate=cell_associate_B))
#
#
# Pathlosses for each UE in operator B are calculated!
pathlosses_B = np.fromiter((compute_pathloss_fspl(distance=_,
                                                  frequency=CARRIER_FREQUENCY) for _ in distances_B),
                           distances_B.dtype)
#
#
# Received power for each UE in operator B is calculated in dB!
rxpower_B = compute_rx_power(P_tx=TRANSMIT_POWER,
                             G_tx=0,
                             G_rx=0,
                             PL=pathlosses_B)
#
#
# SNR for an exclusive licensed model in dB for operator B!
snr_exclusive_B = calculate_snr(rx_power=rxpower_B,
                                noise_power=noise_power_exclusive)
#
#
# SNR for a fully pooled case model in dB for operator B!
snr_pooled_B = calculate_snr(rx_power=rxpower_B,
                             noise_power=noise_power_pooled)
#
#
# Prelog term for an exclusive licensed model for operator B!
prelog_exclusive_B = calculate_prelog_term(bandwidth=exclusive_bandwidth,
                                           number=count_for_ue_B)
#
# Prelog term for a fully pooled case for operator B!
prelog_pooled_B = calculate_prelog_term(bandwidth=pooled_bandwidth,
                                        number=count_for_ue_B)
#
#
# SNR throughput for an exclusive licensed model for operator B!
snr_throughput_exclusive_B = calculate_throughput(dB_values=snr_exclusive_B,
                                                  prelog_term=prelog_exclusive_B)
#
# SNR throughput for a fully pooled case for operator B!
snr_throughput_pooled_B = calculate_throughput(dB_values=snr_pooled_B,
                                               prelog_term=prelog_pooled_B)
#
#
# SINR for an exclusive licensed model in dB for operator B!
sinr_dB_exclusive_B = calculate_sinr_exclusive(pathloss_matrix=pathloss_matrix_fspl_B,
                                               rx_power=rxpower_B,
                                               noise_power=noise_power_exclusive)
#
#
# SINR for a fully pooled case in dB for operator B!
sinr_dB_pooled_B = calculate_sinr_pooled_ii(pathloss_matrix_A=pathloss_matrix_fspl_B,
                                            pathloss_matrix_B=pathloss_matrix_fspl_B_A,
                                            rx_power=rxpower_B,
                                            noise_power=noise_power_pooled)
#
#
# SINR throughput for an exclusive licensed model for operator B!
sinr_throughput_exclusive_B = calculate_throughput(dB_values=sinr_dB_exclusive_B,
                                                   prelog_term=prelog_exclusive_B)
#
# SINR throughput for a fully pooled case for operator B!
sinr_throughput_pooled_B = calculate_throughput(dB_values=sinr_dB_pooled_B,
                                                prelog_term=prelog_pooled_B)


# Simple algorithm test to see if a change in cell association would lead to an increase in performance!
# Concatenate pathloss matrix for operators A and B!


if __name__ == "__main__":


    # Create ECDF distributions for SNR metrics!
    x_snr_exclusive, y_snr_exclusive = create_ecdf(snr_exclusive_A)
    x_snr_pooled, y_snr_pooled = create_ecdf(snr_pooled_A)
    x_snr_throughput_exclusive, y_snr_throughput_exclusive = create_ecdf(limit_throughput_snr(snr_throughput_exclusive_A))
    x_snr_throughput_pooled, y_snr_throughput_pooled = create_ecdf(limit_throughput_snr(snr_throughput_pooled_A))

    # Create ECDF distributions for SINR metrics!
    x_sinr_exclusive, y_sinr_exclusive = create_ecdf(sinr_dB_exclusive_A)
    x_sinr_pooled, y_sinr_pooled = create_ecdf(sinr_dB_pooled_A)
    x_sinr_pooled_iii, y_sinr_pooled_iii = create_ecdf(sinr_dB_pooled_A_iii)
    x_sinr_pooled_iv, y_sinr_pooled_iv = create_ecdf(sinr_dB_pooled_A_iv)
    x_sinr_pooled_v, y_sinr_pooled_v = create_ecdf(sinr_dB_pooled_A_v)
    x_sinr_throughput_exclusive, y_sinr_throughput_exclusive = create_ecdf(limit_throughput_sinr(sinr_throughput_exclusive_A))
    x_sinr_throughput_pooled, y_sinr_throughput_pooled = create_ecdf(limit_throughput_sinr(sinr_throughput_pooled_A))
    x_sinr_throughput_pooled_iii, y_sinr_throughput_pooled_iii = create_ecdf(limit_throughput_sinr(sinr_throughput_pooled_A_iii))
    x_sinr_throughput_pooled_iv, y_sinr_throughput_pooled_iv = create_ecdf(
        limit_throughput_sinr(sinr_throughput_pooled_A_iv))
    x_sinr_throughput_pooled_v, y_sinr_throughput_pooled_v = create_ecdf(
        limit_throughput_sinr(sinr_throughput_pooled_A_v))

    fig1, axs1 = plt.subplots(1, figsize=(5, 5))
    axs1.plot(x_snr_exclusive, y_snr_exclusive, '#2F32FF', label='Exclusive Licensed')
    axs1.plot(x_snr_pooled, y_snr_pooled, '#FF0404', label='Fully Pooled')
    axs1.grid(color='k', linestyle='-', linewidth=0.1)
    axs1.margins(0.05)
    axs1.set_title(r'\bf{SNR at BS density of 25 BSs/$km^2$}', fontsize=12)
    axs1.legend(loc='lower right', frameon=True, prop={'size': 9})
    axs1.set_xlabel(r'\bf{SNR [dB]}', fontsize=12)
    axs1.set_ylabel(r'\bf{Empirical CDF}', fontsize=12)
    axs1.set_xlim(-10, 65)
    axs1.set_rasterized(True)
    fig1.savefig("snr.pdf", bbox_inches='tight', dpi=400)

    fig2, axs2 = plt.subplots(1, figsize=(5, 5))

    axs2.plot(x_snr_throughput_exclusive, y_snr_throughput_exclusive, '#2F32FF', label='Exclusive Licensed')
    axs2.plot(x_snr_throughput_pooled, y_snr_throughput_pooled, '#FF0404', label='Fully Pooled')
    axs2.grid(color='k', linestyle='-', linewidth=0.1)
    axs2.margins(0.05)
    axs2.set_title(r'\bf{User throughput at BS density of 25 BSs/$km^2$}', fontsize=12)
    axs2.legend(loc='lower right', frameon=True, prop={'size': 9})
    axs2.set_xlabel(r'\bf{User throughput [bps]}', fontsize=12)
    axs2.set_ylabel(r'\bf{Empirical CDF}', fontsize=12)
    axs2.set_rasterized(True)
    fig2.savefig("throughput_snr.pdf", bbox_inches='tight', dpi=400)

    fig3, axs3 = plt.subplots(1, figsize=(5, 5))

    axs3.plot(x_sinr_exclusive, y_sinr_exclusive, '#2F32FF', label='Exclusive Licensed')
    axs3.plot(x_sinr_pooled, y_sinr_pooled, '#FF0404', label='Fully Pooled')
    axs3.grid(color='k', linestyle='-', linewidth=0.1)
    axs3.margins(0.05)
    axs3.set_title(r'\bf{SINR at BS density of operator A equalling to 25 BSs/$km^2$}', fontsize=12)
    axs3.legend(loc='lower right', frameon=True, prop={'size': 9})
    axs3.set_xlabel(r'\bf{SINR [dB]}', fontsize=12)
    axs3.set_ylabel(r'\bf{Empirical CDF}', fontsize=12)
    axs3.set_xlim(-60, 60)
    axs3.set_rasterized(True)
    fig3.savefig("sinr.pdf", bbox_inches='tight', dpi=400)

    fig4, axs4 = plt.subplots(1, figsize=(5, 5))

    axs4.plot(x_sinr_throughput_exclusive, y_sinr_throughput_exclusive, '#2F32FF', label='Exclusive Licensed')
    axs4.plot(x_sinr_throughput_pooled, y_sinr_throughput_pooled, '#FF0404', label='Fully Pooled')
    axs4.grid(color='k', linestyle='-', linewidth=0.1)
    axs4.margins(0.05)
    axs4.set_title(r'\bf{BS density of operator A equalling to 25 BSs/$km^2$}', fontsize=12)
    axs4.legend(loc='lower right', frameon=True, prop={'size': 9})
    # r'\bf{phase field} $\phi$'
    axs4.set_xlabel(r'\bf{User throughput [bps]}', fontsize=12)
    axs4.set_ylabel(r'\bf{Empirical CDF}', fontsize=12)
    axs4.set_rasterized(True)
    # fig4.savefig('throughput.eps', format='eps', dpi=1000)
    fig4.savefig("throughput_sinr.pdf", bbox_inches='tight', dpi=400)
    # fig4.savefig("throughput_sinr.eps", bbox_inches='tight', dpi=100)

    fig5, axs5 = plt.subplots(1, figsize=(5, 5))

    # axs5.plot(x_sinr_exclusive, y_sinr_exclusive, '#2F32FF', label='Exclusive Licensed')
    # axs5.plot(x_sinr_pooled_iii, y_sinr_pooled_iii, '#FF0404', label='Fully Pooled')
    axs5.plot(x_sinr_throughput_exclusive, y_sinr_throughput_exclusive, '#2F32FF', label='Exclusive Licensed')
    axs5.plot(x_sinr_throughput_pooled_iii, y_sinr_throughput_pooled_iii, '#FF0404', label='Fully Pooled')
    axs5.grid(color='k', linestyle='-', linewidth=0.1)
    axs5.margins(0.05)
    axs5.set_title(r'\bf{User throughput at BS density of 25 BSs/$km^2$}', fontsize=12)
    axs5.legend(loc='lower right', frameon=True, prop={'size': 9})
    # r'\bf{phase field} $\phi$'
    axs5.set_xlabel(r'\bf{User throughput [bps]}', fontsize=12)
    axs5.set_ylabel(r'\bf{Empirical CDF}', fontsize=12)
    axs5.set_rasterized(True)
    # fig4.savefig('throughput.eps', format='eps', dpi=1000)
    fig5.savefig("sinr_iii_.pdf", bbox_inches='tight', dpi=400)
    # fig4.savefig("throughput_sinr.eps", bbox_inches='tight', dpi=100)

    fig6, axs6 = plt.subplots(1, figsize=(5, 5))

    # axs6.plot(x_sinr_exclusive, y_sinr_exclusive, '#2F32FF', label='Exclusive Licensed')
    # axs6.plot(x_sinr_pooled_iv, y_sinr_pooled_iv, '#FF0404', label='Fully Pooled')
    axs6.plot(x_sinr_throughput_exclusive, y_sinr_throughput_exclusive, '#2F32FF', label='Exclusive Licensed')
    axs6.plot(x_sinr_throughput_pooled_iv, y_sinr_throughput_pooled_iv, '#FF0404', label='Fully Pooled')
    axs6.grid(color='k', linestyle='-', linewidth=0.1)
    axs6.margins(0.05)
    axs6.set_title(r'\bf{User throughput at BS density of 25 BSs/$km^2$}', fontsize=12)
    axs6.legend(loc='lower right', frameon=True, prop={'size': 9})
    # r'\bf{phase field} $\phi$'
    axs6.set_xlabel(r'\bf{User throughput [bps]}', fontsize=12)
    axs6.set_ylabel(r'\bf{Empirical CDF}', fontsize=12)
    axs6.set_rasterized(True)
    # fig4.savefig('throughput.eps', format='eps', dpi=1000)
    fig6.savefig("sinr_iv_.pdf", bbox_inches='tight', dpi=400)
    # fig4.savefig("throughput_sinr.eps", bbox_inches='tight', dpi=100)

    fig7, axs7 = plt.subplots(1, figsize=(5, 5))

    # axs7.plot(x_sinr_exclusive, y_sinr_exclusive, '#2F32FF', label='Exclusive Licensed')
    # axs7.plot(x_sinr_pooled_v, y_sinr_pooled_v, '#FF0404', label='Fully Pooled')
    axs7.plot(x_sinr_throughput_exclusive, y_sinr_throughput_exclusive, '#2F32FF', label='Exclusive Licensed')
    axs7.plot(x_sinr_throughput_pooled_v, y_sinr_throughput_pooled_v, '#FF0404', label='Fully Pooled')
    axs7.grid(color='k', linestyle='-', linewidth=0.1)
    axs7.margins(0.05)
    axs7.set_title(r'\bf{User throughput at BS density of 25 BSs/$km^2$}', fontsize=12)
    axs7.legend(loc='lower right', frameon=True, prop={'size': 9})
    # r'\bf{phase field} $\phi$'
    axs7.set_xlabel(r'\bf{User throughput [bps]}', fontsize=12)
    axs7.set_ylabel(r'\bf{Empirical CDF}', fontsize=12)
    axs7.set_rasterized(True)
    # fig4.savefig('throughput.eps', format='eps', dpi=1000)
    fig7.savefig("sinr_v_.pdf", bbox_inches='tight', dpi=400)
    # fig4.savefig("throughput_sinr.eps", bbox_inches='tight', dpi=100)

    plt.show()


