__author__ = 'Awab Asher'

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from threading import Thread

from user import User
from basestation import Basestation
from simulation_parameters import TRANSMIT_POWER


plt.style.use('seaborn-whitegrid')
plt.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 12, })


class CellularNetwork(object):
    """Representation of a cellular network.

        Attributes:
            ue_density: (int) user density in UEs/km^2.
            bs_density_a: (int) basestation density in BSs/km^2.
            bs_density_b: (int) basestation density in BSs/km^2.
            area: (float) simulation area in km^2.
            iterations: number of iterations for Monte Carlo simulation.
            UEs_A: (list) of all user objects of operator A.
            UEs_B: (list) of all user objects of operator B.
            BSs_A: (list) of all basestation objects of operator A.
            BSs_B: (list) of all basestation objects of operator B.
    """

    def __init__(self, ue_density, bs_density_a, bs_density_b, area, iterations):
        self.ue_density = ue_density
        self.bs_density_a = bs_density_a
        self.bs_density_b = bs_density_b
        self.area = area
        self.iterations = iterations
        self.UEs_A = []
        self.UEs_B = []
        self.BSs_A = []
        self.BSs_B = []
        self.BSs_C = []
        self.BSs_D = []
        self.BSs_E = []


    def deploy_users(self):
        """Deploy users according to a Poisson Point Process."""
        # Poisson Point Process simulation adapted from: https://hpaulkeeler.com/poisson-point-process-simulation/
        x_min, x_max = 0, self.area
        y_min, y_max = 0, 1
        x_delta = x_max - x_min
        y_delta = y_max - y_min
        area_total = x_delta * y_delta  # rectangular area dimensions
        poisson_points = np.random.poisson(self.ue_density * area_total)  # Poisson number of points
        xcoords = []
        ycoords = []
        for _ in range(self.iterations): # This loop makes it a Monte Carlo Simulation: Repeated Sampling!
            x = x_delta * np.random.uniform(0, 1, poisson_points) + x_min  # x coordinates of Poisson points
            y = y_delta * np.random.uniform(0, 1, poisson_points) + y_min  # y coordinates of Poisson points
            xcoords.append(x)
            ycoords.append(y)

        coordinates = [list(zip(x, y)) for x, y in zip(xcoords, ycoords)]
        for count in coordinates:
            for cnt in range(len(count)):
                user = User(x_coord=count[cnt][0],
                            y_coord=count[cnt][1],
                            coordinates=count[cnt])
                self.UEs_A.append(user)


    def deploy_users_(self):
        """Deploy users according to a Poisson Point Process."""
        # Poisson Point Process simulation adapted from: https://hpaulkeeler.com/poisson-point-process-simulation/
        x_min, x_max = 0, self.area
        y_min, y_max = 0, 1
        x_delta = x_max - x_min
        y_delta = y_max - y_min
        area_total = x_delta * y_delta  # rectangular area dimensions
        poisson_points = np.random.poisson(self.ue_density * area_total)  # Poisson number of points
        xcoords = []
        ycoords = []
        for _ in range(self.iterations): # This loop makes it a Monte Carlo Simulation: Repeated Sampling!
            x = x_delta * np.random.uniform(0, 1, poisson_points) + x_min  # x coordinates of Poisson points
            y = y_delta * np.random.uniform(0, 1, poisson_points) + y_min  # y coordinates of Poisson points
            xcoords.append(x)
            ycoords.append(y)

        coordinates = [list(zip(x, y)) for x, y in zip(xcoords, ycoords)]
        for count in coordinates:
            for cnt in range(len(count)):
                user = User(x_coord=count[cnt][0],
                            y_coord=count[cnt][1],
                            coordinates=count[cnt])
                self.UEs_B.append(user)


    def deploy_basestations(self):
        """Deploy basestations according to a Poisson Point Process."""
        # Poisson Point Process simulation adapted from: https://hpaulkeeler.com/poisson-point-process-simulation/
        x_min, x_max = 0, self.area
        y_min, y_max = 0, 1
        x_delta = x_max - x_min
        y_delta = y_max - y_min
        area_total = x_delta * y_delta  # rectangular area dimensions
        poisson_points = np.random.poisson(self.bs_density_a * area_total)  # Poisson number of points
        xcoords = []
        ycoords = []
        for _ in range(1):
            x = x_delta * np.random.uniform(0, 1, poisson_points) + x_min  # x coordinates of Poisson points
            y = y_delta * np.random.uniform(0, 1, poisson_points) + y_min  # y coordinates of Poisson points
            xcoords.append(x)
            ycoords.append(y)

        coordinates = [list(zip(x, y)) for x, y in zip(xcoords, ycoords)]
        for count in coordinates:
            for cnt in range(len(count)):
                basestation = Basestation(x_coord=count[cnt][0],
                                          y_coord=count[cnt][1],
                                          coordinates=count[cnt],
                                          transmit_power=TRANSMIT_POWER)
                self.BSs_A.append(basestation)

    def deploy_basestations_(self):
        """Deploy basestations according to a Poisson Point Process."""
        # Poisson Point Process simulation adapted from: https://hpaulkeeler.com/poisson-point-process-simulation/
        x_min, x_max = 0, self.area
        y_min, y_max = 0, 1
        x_delta = x_max - x_min
        y_delta = y_max - y_min
        area_total = x_delta * y_delta  # rectangular area dimensions
        poisson_points = np.random.poisson(self.bs_density_b * area_total)  # Poisson number of points
        xcoords = []
        ycoords = []
        for _ in range(1):
            x = x_delta * np.random.uniform(0, 1, poisson_points) + x_min  # x coordinates of Poisson points
            y = y_delta * np.random.uniform(0, 1, poisson_points) + y_min  # y coordinates of Poisson points
            xcoords.append(x)
            ycoords.append(y)

        coordinates = [list(zip(x, y)) for x, y in zip(xcoords, ycoords)]
        for count in coordinates:
            for cnt in range(len(count)):
                basestation = Basestation(x_coord=count[cnt][0],
                                          y_coord=count[cnt][1],
                                          coordinates=count[cnt],
                                          transmit_power=TRANSMIT_POWER)
                self.BSs_B.append(basestation)

    def _deploy_basestations_(self):
        """Deploy basestations according to a Poisson Point Process."""
        # Poisson Point Process simulation adapted from: https://hpaulkeeler.com/poisson-point-process-simulation/
        x_min, x_max = 0, self.area
        y_min, y_max = 0, 1
        x_delta = x_max - x_min
        y_delta = y_max - y_min
        area_total = x_delta * y_delta  # rectangular area dimensions
        poisson_points = np.random.poisson(self.bs_density_b * area_total)  # Poisson number of points
        xcoords = []
        ycoords = []
        for _ in range(1):
            x = x_delta * np.random.uniform(0, 1, poisson_points) + x_min  # x coordinates of Poisson points
            y = y_delta * np.random.uniform(0, 1, poisson_points) + y_min  # y coordinates of Poisson points
            xcoords.append(x)
            ycoords.append(y)

        coordinates = [list(zip(x, y)) for x, y in zip(xcoords, ycoords)]
        for count in coordinates:
            for cnt in range(len(count)):
                basestation = Basestation(x_coord=count[cnt][0],
                                          y_coord=count[cnt][1],
                                          coordinates=count[cnt],
                                          transmit_power=TRANSMIT_POWER)
                self.BSs_C.append(basestation)

    def _deploy_basestations__(self):
        """Deploy basestations according to a Poisson Point Process."""
        # Poisson Point Process simulation adapted from: https://hpaulkeeler.com/poisson-point-process-simulation/
        x_min, x_max = 0, self.area
        y_min, y_max = 0, 1
        x_delta = x_max - x_min
        y_delta = y_max - y_min
        area_total = x_delta * y_delta  # rectangular area dimensions
        poisson_points = np.random.poisson(self.bs_density_b * area_total)  # Poisson number of points
        xcoords = []
        ycoords = []
        for _ in range(1):
            x = x_delta * np.random.uniform(0, 1, poisson_points) + x_min  # x coordinates of Poisson points
            y = y_delta * np.random.uniform(0, 1, poisson_points) + y_min  # y coordinates of Poisson points
            xcoords.append(x)
            ycoords.append(y)

        coordinates = [list(zip(x, y)) for x, y in zip(xcoords, ycoords)]
        for count in coordinates:
            for cnt in range(len(count)):
                basestation = Basestation(x_coord=count[cnt][0],
                                          y_coord=count[cnt][1],
                                          coordinates=count[cnt],
                                          transmit_power=TRANSMIT_POWER)
                self.BSs_D.append(basestation)


    def __deploy_basestations__(self):
        """Deploy basestations according to a Poisson Point Process."""
        # Poisson Point Process simulation adapted from: https://hpaulkeeler.com/poisson-point-process-simulation/
        x_min, x_max = 0, self.area
        y_min, y_max = 0, 1
        x_delta = x_max - x_min
        y_delta = y_max - y_min
        area_total = x_delta * y_delta  # rectangular area dimensions
        poisson_points = np.random.poisson(self.bs_density_b * area_total)  # Poisson number of points
        xcoords = []
        ycoords = []
        for _ in range(1):
            x = x_delta * np.random.uniform(0, 1, poisson_points) + x_min  # x coordinates of Poisson points
            y = y_delta * np.random.uniform(0, 1, poisson_points) + y_min  # y coordinates of Poisson points
            xcoords.append(x)
            ycoords.append(y)

        coordinates = [list(zip(x, y)) for x, y in zip(xcoords, ycoords)]
        for count in coordinates:
            for cnt in range(len(count)):
                basestation = Basestation(x_coord=count[cnt][0],
                                          y_coord=count[cnt][1],
                                          coordinates=count[cnt],
                                          transmit_power=TRANSMIT_POWER)
                self.BSs_E.append(basestation)


    def runall(self):
        if __name__ == '__main__':
            Thread(target=self.deploy_users).start()
            Thread(target=self.deploy_basestations).start()
            Thread(target=self.deploy_users_).start()
            Thread(target=self.deploy_basestations_).start()


def poisson_simulation(density, area):
    """"Args:
            density: (int) UE or BS density of an operator!
            area:  (int) area of simulation!

        Returns:
            (list) of x-coordinates of the Poisson Point Process!
            (list) of y-coordinates of the Poisson Point Process!
    """
    # Poisson Point Process simulation adapted from: https://hpaulkeeler.com/poisson-point-process-simulation/
    x_min, x_max = 0, area
    y_min, y_max = 0, 1
    x_delta = x_max - x_min
    y_delta = y_max - y_min
    area_total = x_delta * y_delta  # rectangular area dimensions
    poisson_points = np.random.poisson(density * area_total)  # Poisson number of points
    xcoords, ycoords = [], []
    x = x_delta * np.random.uniform(0, 1, poisson_points) + x_min  # x coordinates of Poisson points
    y = y_delta * np.random.uniform(0, 1, poisson_points) + y_min  # y coordinates of Poisson points
    xcoords.append(x)
    ycoords.append(y)
    return xcoords, ycoords


def plot_poisson_operators(ue_density, bs_density, area):
    """"Args:
            ue_density: (int) UE density of both operators!
            bs_density: (int) BS density of both operators!
            area: (int) area of simulation!

        Returns:
            (scatter plot) of the Poisson Point Process!
    """
    x_ue_coords_A, y_ue_coords_A = poisson_simulation(ue_density, area)
    x_ue_coords_B, y_ue_coords_B = poisson_simulation(ue_density, area)
    x_bs_coords_A, y_bs_coords_A = poisson_simulation(bs_density, area)
    x_bs_coords_B, y_bs_coords_B = poisson_simulation(bs_density, area)
    x_bs_coords_C, y_bs_coords_C = poisson_simulation(bs_density, area)
    x_bs_coords_D, y_bs_coords_D = poisson_simulation(bs_density, area)
    x_bs_coords_E, y_bs_coords_E = poisson_simulation(bs_density, area)

    fig1, axes1 = plt.subplots(1, figsize=(5, 5))
    axes1.plot(x_ue_coords_A, y_ue_coords_A, c='#3F3C3D', marker='+', linestyle='none', label='+')
    axes1.plot(x_ue_coords_B, y_ue_coords_B, c='#000FFA', marker='+', linestyle='none', label='+')
    axes1.scatter(x_bs_coords_A, y_bs_coords_A, edgecolor='#00FA09', facecolor='none', alpha=0.9, label='o')
    axes1.scatter(x_bs_coords_B, y_bs_coords_B, edgecolor='#F40606', facecolor='none', alpha=0.9, label='o')
    # axes1.scatter(x_bs_coords_C, y_bs_coords_C, edgecolor='#F40606', facecolor='none', alpha=0.9, label='o')
    # axes1.scatter(x_bs_coords_D, y_bs_coords_D, edgecolor='#06F4E6', facecolor='none', alpha=0.9, label='o')
    # axes1.scatter(x_bs_coords_E, y_bs_coords_E, edgecolor='#F4E906', facecolor='none', alpha=0.9, label='o')
    axes1.margins(0.05)
    axes1.set_xticks([])
    axes1.set_yticks([])
    # axes1.legend(loc='lower right', frameon=True, prop={'size': 9})
    axes1.grid(False)
    # axes1.set_title(r'\bf{Three Network Operators at BS density of 25 BSs/$km^2$}', fontsize=12)
    axes1.set_rasterized(True)
    fig1.savefig("network.pdf", bbox_inches='tight', dpi=400)

    plt.show()


if __name__ == "__main__":
    plot_poisson_operators(ue_density=500, bs_density=25, area=0.5)
