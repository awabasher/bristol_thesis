__author__ = 'Awab Asher'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.constants
import numpy as np
import math
plt.rc('text', usetex=True)


def compute_pathloss_fspl(distance, frequency):
    """Freespace pathloss between all users and basestations is calculated.

        Args:

            distance: (float) list of numpy arrays (Expected). This function works well for other data types too including 1-D arrays!
            frequency: (int) carrier frequency in Hz!

        Returns:
            (list of) numpy arrays containing the freespace pathloss between a user and all basestations.
    """
    f_Hz = frequency * 10**9
    c = scipy.constants.c   # Speed of light!
    fspl = 10 * np.log10(((4 * math.pi * distance * f_Hz) / c ) ** 2)
    return fspl # returns fspl in dB!


def compute_pathloss_nyc_28_nlos(distance):
    """Empirical NYC model is used to compute pathloss between all users and basestations is calculated. This model supports mmWave bands!

            Args:

                distance: (float) list of numpy array (Expected). This function works well for other data types too including 1-D arrays!

            Returns:
                (list of) numpy arrays containing the empirical NYC pathloss between a user and all basestations.
    """
    alpha = 72.0    # NLOS conditions at 28 GHz!
    beta = 2.92     # NLOS conditions at 28 GHz!
    mu = 0.
    sigma = 8.7

    log_normal_shadowing = np.random.lognormal(mu, sigma, size=None)

    path_loss = alpha + ((10 * beta) * np.log10(distance)) + log_normal_shadowing

    return path_loss    # returns pathloss in dB!


def compute_pathloss_nyc_28_los(distance):
    """Empirical NYC model is used to compute pathloss between all users and basestations is calculated. This model supports mmWave bands!

            Args:

                distance: (float) list of numpy array (Expected). This function works well for other data types too including 1-D arrays!

            Returns:
                (list of) numpy arrays containing the empirical NYC pathloss between a user and all basestations.
    """
    alpha = 61.4    # LOS conditions at 28 GHz!
    beta = 2.       # LOS conditions at 28 GHz!
    mu = 0.
    sigma = 5.8

    log_normal_shadowing = np.random.lognormal(mu, sigma, size=None)

    path_loss = alpha + ((10 * beta) * np.log10(distance)) + log_normal_shadowing

    return path_loss    # returns pathloss in dB!


def compute_pathloss_nyc_73_nlos(distance):
    """Empirical NYC model is used to compute pathloss between all users and basestations is calculated. This model supports mmWave bands!

            Args:

                distance: (float) list of numpy array (Expected). This function works well for other data types too including 1-D arrays!

            Returns:
                (list of) numpy arrays containing the empirical NYC pathloss between a user and all basestations.
    """
    alpha = 86.6    # NLOS conditions at 73 GHz!
    beta = 2.45     # NLOS conditions at 73 GHz!
    mu = 0.
    sigma = 8.0

    log_normal_shadowing = np.random.lognormal(mu, sigma, size=None)

    path_loss = alpha + ((10 * beta) * np.log10(distance)) + log_normal_shadowing

    return path_loss    # returns pathloss in dB!


def compute_pathloss_nyc_73_los(distance):
    """Empirical NYC model is used to compute pathloss between all users and basestations is calculated. This model supports mmWave bands!

            Args:

                distance: (float) list of numpy array (Expected). This function works well for other data types too including 1-D arrays!

            Returns:
                (list of) numpy arrays containing the empirical NYC pathloss between a user and all basestations.
    """
    alpha = 69.8    # LOS conditions at 73 GHz!
    beta = 2.   # LOS conditions at 73 GHz!
    mu = 0.
    sigma = 5.8

    log_normal_shadowing = np.random.lognormal(mu, sigma, size=None)

    path_loss = alpha + ((10 * beta) * np.log10(distance)) + log_normal_shadowing

    return path_loss    # returns pathloss in dB!


def compute_antenna_gain(phi):
    """Computes the antenna gain given the horizontal angle between user and basestation (3GPP spec)!

            Args:
                phi: (float) horizontal angle between user and basestation in degrees!

            Returns:
                (float) horizontal angle gain!
    """
    Am = 30     # Front-back ratio in dB!
    horizontal_beamwidth = 65   # Horizontal 3 dB beamwidth in degrees!
    gain = -min(12.*pow((phi / horizontal_beamwidth), 2), Am)
    return gain


# print(compute_antenna_gain(70))


if __name__ == "__main__":

    frequencies = [0.5, 1, 6, 28, 50, 73]

    distances = np.linspace(0.1, 200, num=200)  # distances in meters!
    distances_nyc = np.linspace(10, 200, num=200)  # distances in meters!

    fspl_at_500_mHz = compute_pathloss_fspl(distances, frequencies[0])
    fspl_at_1_gHz = compute_pathloss_fspl(distances, frequencies[1])
    fspl_at_6_gHz = compute_pathloss_fspl(distances, frequencies[2])
    fspl_at_28_gHz = compute_pathloss_fspl(distances, frequencies[3])
    fspl_at_50_gHz = compute_pathloss_fspl(distances, frequencies[4])
    fspl_at_73_gHz = compute_pathloss_fspl(distances, frequencies[5])

    empirical_nyc_28_nlos = compute_pathloss_nyc_28_nlos(distances_nyc)
    empirical_nyc_28_los = compute_pathloss_nyc_28_los(distances_nyc)
    empirical_nyc_73_nlos = compute_pathloss_nyc_73_nlos(distances_nyc)
    empirical_nyc_73_los = compute_pathloss_nyc_73_los(distances_nyc)

    fig1, axs = plt.subplots(1, figsize=(5, 5))
    axs.semilogx(distances, fspl_at_500_mHz, '#2F32FF', label='500 MHz')
    axs.semilogx(distances, fspl_at_1_gHz, '#FF0404', label='1 GHz')
    axs.semilogx(distances, fspl_at_6_gHz, '#00FFD8', label='6 GHz')
    axs.semilogx(distances, fspl_at_28_gHz, '#FFFF00', label='28 GHz')
    axs.semilogx(distances, fspl_at_50_gHz, '#46FF00', label='50 GHz')
    axs.semilogx(distances, fspl_at_73_gHz, '#FF4200', label='73 GHz')
    axs.set_xlim(0.1, 200)
    axs.grid(color='k', linestyle='-', linewidth=0.08)
    axs.set_title(r'\bf{Path Loss as a function of carrier frequency}', fontsize=12)
    axs.set_xlabel(r'\bf{Distance (m)}', fontsize=12)
    axs.set_ylabel(r'\bf{Path Loss in dB}', fontsize=12)
    axs.legend(loc='lower right', frameon=True, prop={'size': 7})
    axs.set_rasterized(True)
    fig1.savefig('fspl.pdf', bbox_inches = 'tight', dpi = 400)

    fig, (axs1, axs2) = plt.subplots(1, 2, figsize=(5, 5))
    axs1.semilogx(distances, empirical_nyc_28_nlos, '#2F32FF', label='NLOS')
    axs1.semilogx(distances, empirical_nyc_28_los, '#46FF00', label='LOS')
    axs1.grid(True, which="both", ls="--", color='0.85')
    axs1.set_title('28 GHz')
    axs1.set_xlabel('Distance (m)')
    axs1.set_ylabel('Path Loss in dB')
    axs1.set_ylim(80, 150)
    axs1.set_xlim(10, 200)
    axs1.legend(loc='lower right', frameon=True, prop={'size': 9})

    axs2.semilogx(distances, empirical_nyc_73_nlos, '#2F32FF', label='NLOS')
    axs2.semilogx(distances, empirical_nyc_73_los, '#46FF00', label='LOS')
    axs2.grid(True, which="both", ls="--", color='0.85')
    axs2.set_title('73 GHz')
    axs2.set_xlabel('Distance (m)')
    axs2.set_ylim(80, 150)
    axs2.set_xlim(10, 200)
    axs2.legend(loc='lower right', frameon=True, prop={'size': 9})

    plt.subplots_adjust(hspace=1.5)
    # fig.savefig('empirical_nyc.eps', format='eps', dpi=1000)
    # fig.savefig('empirical_nyc.png')

    plt.show()