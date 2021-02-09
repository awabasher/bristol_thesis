__author__ = 'Awab Asher'

import numpy as np
import statsmodels.api as sm
from scipy.spatial import distance
from collections import Counter

from simulation_parameters import(
    TRANSMIT_POWER,
    CARRIER_FREQUENCY
)

from pathloss import (
    compute_pathloss_fspl,
    compute_pathloss_nyc_28_los
)


def convert_km_to_m(km_distance):
    """Function that converts distance in km to m!"""
    _ = km_distance * 10**3
    return _


def compute_distance_matrix(users, basestations):
    """Distances between all users and basestations is calculated.

    Args:

        users: (obj) list of users!
        basestations: (obj) list of basestations!

    Returns:
        (list of) numpy arrays containing the distance between a user and all basestations in km!.
    """

    coords_list_ue = [getattr(ele, 'coordinates') for ele in users]
    coords_list_bs = [getattr(ele, 'coordinates') for ele in basestations]
    distance_matrix = []

    count = 0
    for _ in coords_list_ue:
        element = [coords_list_ue[count]]
        coords = element + coords_list_bs
        dist = distance.cdist(coords, coords, 'euclidean')
        new_dist = np.delete(dist[0], 0)
        distance_matrix.append(new_dist)
        count += 1

    return np.array(distance_matrix)


def compute_distance_matrix_m(distance_matrix_km):
    """Distance matrix is converted from km to m!

        Args:
            distance_matrix_km: (list of) numpy arrays containing the distance between a UE and all BSs in km!

        Returns:
            (list of) numpy arrays containing the distance between a user and all basestations in m!.
    """
    distance_matrix_m = [convert_km_to_m(_) for _ in distance_matrix_km]

    return distance_matrix_m


def compute_pathloss_matrix_fspl(_dist_matrix_m):
    """Pathlosses for all users is calculated.

    Args:

        _dist_matrix_m: (obj) list of users!

    Returns:
        (list of) numpy arrays containing the distance between a user and all basestations in km!.
    """
    pl_matrix = [compute_pathloss_fspl(distance=_, frequency=CARRIER_FREQUENCY) for _ in _dist_matrix_m]

    return pl_matrix


def compute_pathloss_matrix_nyc(_dist_matrix_m):
    """Pathlosses for all users is calculated.

        Args:

            _dist_matrix_m: (obj) list of users!

        Returns:
            (list of) numpy arrays containing the distance between a user and all basestations in km!.
    """
    pl_matrix = [compute_pathloss_nyc_28_los(_) for _ in _dist_matrix_m]

    return pl_matrix


def cell_associate(pathloss_matrix, users, basestations):
    """Associate a user with a basestation that provides the minimum pathloss.

        Args:

            pathloss_matrix: (list) of numpy arrays
            users: (obj) list of users!
            basestations: (obj) list of basestations!

        Returns:
            (list of) tuples containing the UE object and the BS it is associated to.
    """

    index_list_min_pl = []  # List that contains tuple of (index(min(pathloss)), pathloss) for eacb UE!
    list_bs = []  # List of basestation objects associated with each UE in order!
    for _ in pathloss_matrix:
        index_list_min_pl.append(min(enumerate(_), key=(lambda x: x[1])))

    for _ in index_list_min_pl:
        index = _[0]
        list_bs.append(basestations[index])

    cell_associate_list = list(zip(users, list_bs))  # List that contains tuple: (ue_object, bs_associated)!

    return cell_associate_list


def compute_count_for_bs(pathloss_matrix, basestations):
    """Computes the number of UEs associated with a BS object for BS.

        Args:

            pathloss_matrix: (list) of numpy arrays
            basestations: (obj) list of basestations!

        Returns:
            (list of) tuples containing the BS object and the number of UEs it is associated to.
    """

    index_list_min_pl = []  # List that contains tuple of (index(min(pathloss)), pathloss) for eacb UE!
    list_bs = []  # List of basestation objects associated with each UE in order!
    for _ in pathloss_matrix:
        index_list_min_pl.append(min(enumerate(_), key=(lambda x: x[1])))

    for _ in index_list_min_pl:
        index = _[0]
        list_bs.append(basestations[index])

    bs_cnt = Counter(list_bs)

    bs_count = list(zip(bs_cnt.keys(),
                        bs_cnt.values()))  # List of tuples that contain BS objects and the number of UEs they are associated to!

    return bs_count


def compute_total_ue_for_bs(_cell_associate, _bs_count):
    """Computes the number of UEs associated with a BS object for UE.

        Args:

            _cell_associate: (list) that contains tuple (ue_object, bs_associated)!
            _bs_count: (list) that contains tuple of BS objects and the number of UEs they are associated to!

        Returns:
            (list of) ints that correspond to the count of total number of UEs associated with a BS.
    """
    other_ue_count_bs = []  # List of the total number of UEs associated with a BS for each UE!

    for x in _cell_associate:
        for y in _bs_count:
            if x[1] == y[0]:
                other_ue_count_bs.append(y[1])

    return other_ue_count_bs


def compute_distance_ue_bs(_cell_associate):
    """Computes the distance between the user and the basestation it is associated to.

            Args:

                _cell_associate: (list) that contains tuple (ue_object, bs_associated)!

            Returns:
                (list of) distances in metres between each UE and BS.
    """
    distance_in_km = [distance.euclidean(_[0].coordinates, _[1].coordinates) for _ in _cell_associate]
    distance_in_m = [convert_km_to_m(_) for _ in distance_in_km]

    return distance_in_m


def compute_rx_power(P_tx, G_tx, G_rx, PL):
    """ Link Budget!
    Args:
        P_tx: transmit output power (dBm)
        G_tx: transmitter antenna gain (dBi)
        G_rx: receiver antenna gain (dBi)
        PL: pathloss in (dB)!

    Returns:
        P_rx:  received power (dB)!!!
    """
    p_rx = P_tx + G_tx + G_rx - PL - 30
    return p_rx     # returns P_rx in dB!


def calculate_snr(rx_power, noise_power):
    """ Function that calculates SNR in dB!

    Args:
        rx_power: (numpy array) received power in dB!
        noise_power: noise power in dB!

    Returns:
        snr: (numpy array) Signal-to-Noise ratio in dB!
    """
    snr = rx_power - noise_power  # rx_power and noise_power should be in dB!
    return snr  # This gives SNR in dB!


def calculate_log_2_factor(value):
    """ Function that calculates the logarithmic term!"""
    result = np.log2(1 + value)
    return result


def calculate_prelog_term(bandwidth, number):
    """ Function that calculates the prelog term!

        Args:
            bandwidth: (int) bandwidth allocated to basestation which is can be exclusive or pooled!
            number: (numpy array) (N) which is the number of total users associated to every basestation for the basestation with
            which the UE of interest is associated to!

        Returns:
            prelog_term: (float) prelog_term is the ratio of the bandwidth and number!
    """
    prelog_term = bandwidth / number
    return prelog_term


def create_ecdf(metric_values):
    """Function to create ecdf (Empirical Cumulative Distribution Function!).

            Args:
                metric_values: (numpy array) of values that needs ecdf!

            Returns:
                x: (list of or numpy array) x values for ecdf!
                y: (list of or numpy array) y values for ecdf!
    """
    ecdf = sm.distributions.ECDF(metric_values)
    x = np.linspace(min(metric_values), max(metric_values))
    y = ecdf(x)

    return x, y


def convert_dB_to_W(dB_value):
    """ Function that converts dB values into Watts!"""
    _ = 10 ** (dB_value / 10)
    return _


def convert_W_to_dB(W_value):
    """ Function that converts power values in Watts into dB values!"""
    _ = 10 * np.log10(W_value)
    return _


def calculate_sinr_exclusive(pathloss_matrix, rx_power, noise_power):
    """Function to calculate SINR in dB!

        Args:
            pathloss_matrix: (list of) numpy arrays!
            rx_power: (numpy array) received power in dB for each UE!
            noise_power: (constant term) noise power in dB!

        Returns:
            sinr: (numpy array) Signal-to-Interference-plus-Noise ratio in dB!
    """
    pathloss_matrix_ = np.array(pathloss_matrix)
    with np.nditer(pathloss_matrix_, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix = pathloss_matrix_

    with np.nditer(rxpower_matrix, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_W = rxpower_matrix
    all_rx_power_W = np.sum(rx_power_matrix_W, axis=1)
    rx_power_W = np.fromiter((convert_dB_to_W(_) for _ in rx_power), rx_power.dtype)
    inter_cell_interference = all_rx_power_W - rx_power_W
    sum_factor = inter_cell_interference + convert_dB_to_W(noise_power)
    sinr_watts = rx_power_W / sum_factor
    sinr_dB = np.fromiter((convert_W_to_dB(_) for _ in sinr_watts), sinr_watts.dtype)

    return sinr_dB



def calculate_sinr_pooled_ii(pathloss_matrix_A, pathloss_matrix_B, rx_power, noise_power):
    """Function to calculate SINR in dB!

        Args:
            pathloss_matrix_A: (list of) numpy arrays!
            pathloss_matrix_B: (list of) numpy arrays!
            rx_power: (numpy array) received power in dB for each UE!
            noise_power: (constant term) noise power in dB!

        Returns:
            sinr: (numpy array) Signal-to-Interference-plus-Noise ratio in dB!
    """
    pathloss_matrix_A = np.array(pathloss_matrix_A)
    with np.nditer(pathloss_matrix_A, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_A = pathloss_matrix_A

    pathloss_matrix_B = np.array(pathloss_matrix_B)
    with np.nditer(pathloss_matrix_B, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_B = pathloss_matrix_B

    with np.nditer(rxpower_matrix_A, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_A_W = rxpower_matrix_A

    with np.nditer(rxpower_matrix_B, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_B_W = rxpower_matrix_B

    all_rx_power_A_W = np.sum(rx_power_matrix_A_W, axis=1)
    all_rx_power_B_W = np.sum(rx_power_matrix_B_W, axis=1)
    rx_power_W = np.fromiter((convert_dB_to_W(_) for _ in rx_power), rx_power.dtype)
    inter_cell_interference = all_rx_power_A_W + all_rx_power_B_W - rx_power_W
    sum_factor = inter_cell_interference + convert_dB_to_W(noise_power)
    sinr_watts = rx_power_W / sum_factor
    sinr_dB = np.fromiter((convert_W_to_dB(_) for _ in sinr_watts), sinr_watts.dtype)

    return sinr_dB



def calculate_sinr_pooled_iii(pathloss_matrix_A,
                              pathloss_matrix_B,
                              pathloss_matrix_C,
                              rx_power,
                              noise_power):
    """Function to calculate SINR in dB!

        Args:
            pathloss_matrix_A: (list of) numpy arrays!
            pathloss_matrix_B: (list of) numpy arrays!
            pathloss_matrix_C: (list of) numpy arrays!
            rx_power: (numpy array) received power in dB for each UE!
            noise_power: (constant term) noise power in dB!

        Returns:
            sinr: (numpy array) Signal-to-Interference-plus-Noise ratio in dB!
    """
    pathloss_matrix_A = np.array(pathloss_matrix_A)
    with np.nditer(pathloss_matrix_A, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_A = pathloss_matrix_A

    pathloss_matrix_B = np.array(pathloss_matrix_B)
    with np.nditer(pathloss_matrix_B, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_B = pathloss_matrix_B

    pathloss_matrix_C = np.array(pathloss_matrix_C)
    with np.nditer(pathloss_matrix_C, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_C = pathloss_matrix_C

    with np.nditer(rxpower_matrix_A, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_A_W = rxpower_matrix_A

    with np.nditer(rxpower_matrix_B, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_B_W = rxpower_matrix_B

    with np.nditer(rxpower_matrix_C, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_C_W = rxpower_matrix_C

    all_rx_power_A_W = np.sum(rx_power_matrix_A_W, axis=1)
    all_rx_power_B_W = np.sum(rx_power_matrix_B_W, axis=1)
    all_rx_power_C_W = np.sum(rx_power_matrix_C_W, axis=1)
    rx_power_W = np.fromiter((convert_dB_to_W(_) for _ in rx_power), rx_power.dtype)
    inter_cell_interference = all_rx_power_A_W + all_rx_power_B_W + all_rx_power_C_W - rx_power_W
    sum_factor = inter_cell_interference + convert_dB_to_W(noise_power)
    sinr_watts = rx_power_W / sum_factor
    sinr_dB = np.fromiter((convert_W_to_dB(_) for _ in sinr_watts), sinr_watts.dtype)

    return sinr_dB

def calculate_sinr_pooled_iv(pathloss_matrix_A,
                             pathloss_matrix_B,
                             pathloss_matrix_C,
                             pathloss_matrix_D,
                             rx_power,
                             noise_power):
    """Function to calculate SINR in dB!

        Args:
            pathloss_matrix_A: (list of) numpy arrays!
            pathloss_matrix_B: (list of) numpy arrays!
            pathloss_matrix_C: (list of) numpy arrays!
            pathloss_matrix_D: (list of) numpy arrays!
            rx_power: (numpy array) received power in dB for each UE!
            noise_power: (constant term) noise power in dB!

        Returns:
            sinr: (numpy array) Signal-to-Interference-plus-Noise ratio in dB!
    """
    pathloss_matrix_A = np.array(pathloss_matrix_A)
    with np.nditer(pathloss_matrix_A, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_A = pathloss_matrix_A

    pathloss_matrix_B = np.array(pathloss_matrix_B)
    with np.nditer(pathloss_matrix_B, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_B = pathloss_matrix_B

    pathloss_matrix_C = np.array(pathloss_matrix_C)
    with np.nditer(pathloss_matrix_C, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_C = pathloss_matrix_C

    pathloss_matrix_D = np.array(pathloss_matrix_D)
    with np.nditer(pathloss_matrix_D, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_D = pathloss_matrix_D

    with np.nditer(rxpower_matrix_A, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_A_W = rxpower_matrix_A

    with np.nditer(rxpower_matrix_B, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_B_W = rxpower_matrix_B

    with np.nditer(rxpower_matrix_C, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_C_W = rxpower_matrix_C

    with np.nditer(rxpower_matrix_D, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_D_W = rxpower_matrix_D

    all_rx_power_A_W = np.sum(rx_power_matrix_A_W, axis=1)
    all_rx_power_B_W = np.sum(rx_power_matrix_B_W, axis=1)
    all_rx_power_C_W = np.sum(rx_power_matrix_C_W, axis=1)
    all_rx_power_D_W = np.sum(rx_power_matrix_D_W, axis=1)
    rx_power_W = np.fromiter((convert_dB_to_W(_) for _ in rx_power), rx_power.dtype)
    inter_cell_interference = all_rx_power_A_W + all_rx_power_B_W + all_rx_power_C_W + all_rx_power_D_W - rx_power_W
    sum_factor = inter_cell_interference + convert_dB_to_W(noise_power)
    sinr_watts = rx_power_W / sum_factor
    sinr_dB = np.fromiter((convert_W_to_dB(_) for _ in sinr_watts), sinr_watts.dtype)

    return sinr_dB


def calculate_sinr_pooled_v(pathloss_matrix_A,
                            pathloss_matrix_B,
                            pathloss_matrix_C,
                            pathloss_matrix_D,
                            pathloss_matrix_E,
                            rx_power,
                            noise_power):
    """Function to calculate SINR in dB!

        Args:
            pathloss_matrix_A: (list of) numpy arrays!
            pathloss_matrix_B: (list of) numpy arrays!
            pathloss_matrix_C: (list of) numpy arrays!
            pathloss_matrix_D: (list of) numpy arrays!
            pathloss_matrix_E: (list of) numpy arrays!
            rx_power: (numpy array) received power in dB for each UE!
            noise_power: (constant term) noise power in dB!

        Returns:
            sinr: (numpy array) Signal-to-Interference-plus-Noise ratio in dB!
    """
    pathloss_matrix_A = np.array(pathloss_matrix_A)
    with np.nditer(pathloss_matrix_A, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_A = pathloss_matrix_A

    pathloss_matrix_B = np.array(pathloss_matrix_B)
    with np.nditer(pathloss_matrix_B, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_B = pathloss_matrix_B

    pathloss_matrix_C = np.array(pathloss_matrix_C)
    with np.nditer(pathloss_matrix_C, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_C = pathloss_matrix_C

    pathloss_matrix_D = np.array(pathloss_matrix_D)
    with np.nditer(pathloss_matrix_D, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_D = pathloss_matrix_D

    pathloss_matrix_E = np.array(pathloss_matrix_E)
    with np.nditer(pathloss_matrix_E, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = compute_rx_power(P_tx=TRANSMIT_POWER, G_tx=0, G_rx=0, PL=x)

    rxpower_matrix_E = pathloss_matrix_E

    with np.nditer(rxpower_matrix_A, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_A_W = rxpower_matrix_A

    with np.nditer(rxpower_matrix_B, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_B_W = rxpower_matrix_B


    with np.nditer(rxpower_matrix_C, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_C_W = rxpower_matrix_C

    with np.nditer(rxpower_matrix_D, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_D_W = rxpower_matrix_D


    with np.nditer(rxpower_matrix_E, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = convert_dB_to_W(x)

    rx_power_matrix_E_W = rxpower_matrix_E


    all_rx_power_A_W = np.sum(rx_power_matrix_A_W, axis=1)
    all_rx_power_B_W = np.sum(rx_power_matrix_B_W, axis=1)
    all_rx_power_C_W = np.sum(rx_power_matrix_C_W, axis=1)
    all_rx_power_D_W = np.sum(rx_power_matrix_D_W, axis=1)
    all_rx_power_E_W = np.sum(rx_power_matrix_E_W, axis=1)

    rx_power_W = np.fromiter((convert_dB_to_W(_) for _ in rx_power), rx_power.dtype)
    inter_cell_interference = all_rx_power_A_W + all_rx_power_B_W + all_rx_power_C_W + all_rx_power_D_W + all_rx_power_E_W - rx_power_W
    sum_factor = inter_cell_interference + convert_dB_to_W(noise_power)
    sinr_watts = rx_power_W / sum_factor
    sinr_dB = np.fromiter((convert_W_to_dB(_) for _ in sinr_watts), sinr_watts.dtype)

    return sinr_dB


def limit_throughput_snr(throughput):

    _ = list(filter(lambda x: x < 1.5 * 10 ** 7, throughput))

    return _


def limit_throughput_sinr(throughput):
    """Function to limit throughput for SINR!

        Args:
            throughput: (limit) User throughput in bps!

        Returns:
            _: (list) Limited User throughput in bps!
    """

    _ = list(filter(lambda x: x < 5 * 10 ** 6, throughput))

    return _


def calculate_throughput(dB_values, prelog_term):
    """Function to calculate user throughput in bps!

        Args:
            dB_values: (numpy array) SNR or SINR values in dB!
            prelog_term: (numpy array) Prelog term!

        Returns:
            throughput: (numpy array) User throughput in bps!
    """
    # Converts dB values into Watts!
    watt_values = np.fromiter((convert_dB_to_W(_) for _ in dB_values), dB_values.dtype)

    # Calculates the log2 term!
    log_2_term = np.fromiter((calculate_log_2_factor(_) for _ in watt_values), watt_values.dtype)

    # Calculates the throughput in bps!
    throughput = prelog_term * log_2_term

    return throughput

