__author__ = 'Awab Asher'
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
plt.rc('text', usetex=True)

from simulation_parameters import (
    AREA_OF_SIMULATION
)


def plot_poisson(density, area):
    """"Args:
            density: (int) density to determine intensity of the Poisson Point Process!
            area: (float) simulation area!
        Returns:
            (scatter plot) of the Poisson Point Process for two operators having same densities at different BS densities!
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
    return x, y


def plot_two_operators_at_different_densities(ue_density, start_density, end_density, area):
    """"Args:
            ue_density: (int) UE density of both operators!
            start_density: (int) lower threshold of BS density!
            end_density: (int) upper threshold of BS density!
            area: (float) simulation area!

        Returns:
            (scatter plot) of the Poisson Point Process for two operators having same densities at different BS densities!
    """
    num = int((end_density / start_density))
    basestation_densities = np.linspace(start_density, end_density, num=num)

    coords_ue_x, coords_ue_y = plot_poisson(ue_density, area)
    x_bs_25, y_bs_25 = plot_poisson(basestation_densities[0], area)
    x_bs_50, y_bs_50 = plot_poisson(basestation_densities[1], area)
    x_bs_75, y_bs_75 = plot_poisson(basestation_densities[2], area)
    x_bs_100, y_bs_100 = plot_poisson(basestation_densities[3], area)
    x_bs_125, y_bs_125 = plot_poisson(basestation_densities[4], area)
    x_bs_150, y_bs_150 = plot_poisson(basestation_densities[5], area)
    x_bs_175, y_bs_175 = plot_poisson(basestation_densities[6], area)
    x_bs_200, y_bs_200 = plot_poisson(basestation_densities[7], area)

    x_bs_25_, y_bs_25_ = plot_poisson(basestation_densities[0], area)
    x_bs_50_, y_bs_50_ = plot_poisson(basestation_densities[1], area)
    x_bs_75_, y_bs_75_ = plot_poisson(basestation_densities[2], area)
    x_bs_100_, y_bs_100_ = plot_poisson(basestation_densities[3], area)
    x_bs_125_, y_bs_125_ = plot_poisson(basestation_densities[4], area)
    x_bs_150_, y_bs_150_ = plot_poisson(basestation_densities[5], area)
    x_bs_175_, y_bs_175_ = plot_poisson(basestation_densities[6], area)
    x_bs_200_, y_bs_200_ = plot_poisson(basestation_densities[7], area)

    # Non-homogeneous intensities!
    fig1, axes = plt.subplots(2, 4, figsize=(13, 6))
    axes[0, 0].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none')
    axes[0, 0].scatter(x_bs_25, y_bs_25, edgecolor='#00FA09', facecolor='none', alpha=0.9)
    axes[0, 0].scatter(x_bs_25_, y_bs_25_, edgecolor='#000FFA', facecolor='none', alpha=0.9)
    axes[0, 0].margins(0.05)
    axes[0, 0].set_title(r'\bf{25 BSs/$km^2$}', fontsize=12)
    axes[0, 0].set_xticks([])
    axes[0, 0].set_yticks([])

    axes[0, 1].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none')
    axes[0, 1].scatter(x_bs_25, y_bs_25, edgecolor='#00FA09', facecolor='none', alpha=0.9)
    axes[0, 1].scatter(x_bs_50, y_bs_50, edgecolor='#000FFA', facecolor='none', alpha=0.9)
    axes[0, 1].margins(0.05)
    axes[0, 1].set_title(r'\bf{50 BSs/$km^2$}', fontsize=12)
    axes[0, 1].set_xticks([])
    axes[0, 1].set_yticks([])

    axes[0, 2].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none')
    axes[0, 2].scatter(x_bs_25, y_bs_25, edgecolor='#00FA09', facecolor='none', alpha=0.9)
    axes[0, 2].scatter(x_bs_75, y_bs_75, edgecolor='#000FFA', facecolor='none', alpha=0.9)
    axes[0, 2].margins(0.05)
    axes[0, 2].set_title(r'\bf{75 BSs/$km^2$}', fontsize=12)
    axes[0, 2].set_xticks([])
    axes[0, 2].set_yticks([])

    axes[0, 3].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none')
    axes[0, 3].scatter(x_bs_25, y_bs_25, edgecolor='#00FA09', facecolor='none', alpha=0.9)
    axes[0, 3].scatter(x_bs_100, y_bs_100, edgecolor='#000FFA', facecolor='none', alpha=0.9)
    axes[0, 3].margins(0.05)
    axes[0, 3].set_title(r'\bf{100 BSs/$km^2$}', fontsize=12)
    axes[0, 3].set_xticks([])
    axes[0, 3].set_yticks([])

    axes[1, 0].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none')
    axes[1, 0].scatter(x_bs_25, y_bs_25, edgecolor='#00FA09', facecolor='none', alpha=0.9)
    axes[1, 0].scatter(x_bs_125, y_bs_125, edgecolor='#000FFA', facecolor='none', alpha=0.9)
    axes[1, 0].margins(0.05)
    axes[1, 0].set_title(r'\bf{125 BSs/$km^2$}', fontsize=12)
    axes[1, 0].set_xticks([])
    axes[1, 0].set_yticks([])

    axes[1, 1].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none')
    axes[1, 1].scatter(x_bs_25, y_bs_25, edgecolor='#00FA09', facecolor='none', alpha=0.9)
    axes[1, 1].scatter(x_bs_150, y_bs_150, edgecolor='#000FFA', facecolor='none', alpha=0.9)
    axes[1, 1].margins(0.05)
    axes[1, 1].set_title(r'\bf{150 BSs/$km^2$}', fontsize=12)
    axes[1, 1].set_xticks([])
    axes[1, 1].set_yticks([])

    axes[1, 2].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none')
    axes[1, 2].scatter(x_bs_25, y_bs_25, edgecolor='#00FA09', facecolor='none', alpha=0.9)
    axes[1, 2].scatter(x_bs_175, y_bs_175, edgecolor='#000FFA', facecolor='none', alpha=0.9)
    axes[1, 2].margins(0.05)
    axes[1, 2].set_title(r'\bf{175 BSs/$km^2$}', fontsize=12)
    axes[1, 2].set_xticks([])
    axes[1, 2].set_yticks([])

    axes[1, 3].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none')
    axes[1, 3].scatter(x_bs_25, y_bs_25, edgecolor='#00FA09', facecolor='none', alpha=0.9)
    axes[1, 3].scatter(x_bs_200, y_bs_200, edgecolor='#000FFA', facecolor='none', alpha=0.9)
    axes[1, 3].margins(0.05)
    axes[1, 3].set_title(r'\bf{200 BSs/$km^2$}', fontsize=12)
    axes[1, 3].set_xticks([])
    axes[1, 3].set_yticks([])

    fig1.savefig("poisson_nonhomogeneous.pdf", bbox_inches='tight', dpi=400)

    # Homogeneous intensities!
    fig, axs = plt.subplots(2, 4, figsize=(13, 6))
    axs[0, 0].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none', label='+')
    axs[0, 0].scatter(x_bs_25, y_bs_25, edgecolor='#00FA09', facecolor='none', alpha=0.9, label='o')
    axs[0, 0].scatter(x_bs_25_, y_bs_25_, edgecolor='#000FFA', facecolor='none', alpha=0.9, label='o')
    axs[0, 0].margins(0.05)
    axs[0, 0].set_title(r'\bf{25 BSs/$km^2$}', fontsize=12)
    axs[0, 0].set_xticks([])
    axs[0, 0].set_yticks([])

    axs[0, 1].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none', label='+')
    axs[0, 1].scatter(x_bs_50, y_bs_50, edgecolor='#00FA09', facecolor='none', alpha=0.9, label='o')
    axs[0, 1].scatter(x_bs_50_, y_bs_50_, edgecolor='#000FFA', facecolor='none', alpha=0.9, label='o')
    axs[0, 1].margins(0.05)
    axs[0, 1].set_title(r'\bf{50 BSs/$km^2$}', fontsize=12)
    axs[0, 1].set_xticks([])
    axs[0, 1].set_yticks([])

    axs[0, 2].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none', label='+')
    axs[0, 2].scatter(x_bs_75, y_bs_75, edgecolor='#00FA09', facecolor='none', alpha=0.9, label='o')
    axs[0, 2].scatter(x_bs_75_, y_bs_75_, edgecolor='#000FFA', facecolor='none', alpha=0.9, label='o')
    axs[0, 2].margins(0.05)
    axs[0, 2].set_title(r'\bf{75 BSs/$km^2$}', fontsize=12)
    axs[0, 2].set_xticks([])
    axs[0, 2].set_yticks([])

    axs[0, 3].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none', label='o')
    axs[0, 3].scatter(x_bs_100, y_bs_100, edgecolor='#00FA09', facecolor='none', alpha=0.9, label='o')
    axs[0, 3].scatter(x_bs_100_, y_bs_100_, edgecolor='#000FFA', facecolor='none', alpha=0.9, label='o')
    axs[0, 3].margins(0.05)
    axs[0, 3].set_title(r'\bf{100 BSs/$km^2$}', fontsize=12)
    axs[0, 3].set_xticks([])
    axs[0, 3].set_yticks([])

    axs[1, 0].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none', label='+')
    axs[1, 0].scatter(x_bs_100, y_bs_100, edgecolor='#00FA09', facecolor='none', alpha=0.9, label='o')
    axs[1, 0].scatter(x_bs_125_, y_bs_125_, edgecolor='#000FFA', facecolor='none', alpha=0.9, label='o')
    axs[1, 0].margins(0.05)
    axs[1, 0].set_title(r'\bf{125 BSs/$km^2$}', fontsize=12)
    axs[1, 0].set_xticks([])
    axs[1, 0].set_yticks([])

    axs[1, 1].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none', label='+')
    axs[1, 1].scatter(x_bs_150, y_bs_150, edgecolor='#00FA09', facecolor='none', alpha=0.9, label='o')
    axs[1, 1].scatter(x_bs_150_, y_bs_150_, edgecolor='#000FFA', facecolor='none', alpha=0.9, label='o')
    axs[1, 1].margins(0.05)
    axs[1, 1].set_title(r'\bf{150 BSs/$km^2$}', fontsize=12)
    axs[1, 1].set_xticks([])
    axs[1, 1].set_yticks([])

    axs[1, 2].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none', label='+')
    axs[1, 2].scatter(x_bs_175, y_bs_175, edgecolor='#00FA09', facecolor='none', alpha=0.9, label='o')
    axs[1, 2].scatter(x_bs_175_, y_bs_175_, edgecolor='#000FFA', facecolor='none', alpha=0.9, label='o')
    axs[1, 2].margins(0.05)
    axs[1, 2].set_title(r'\bf{175 BSs/$km^2$}', fontsize=12)
    axs[1, 2].set_xticks([])
    axs[1, 2].set_yticks([])

    axs[1, 3].plot(coords_ue_x, coords_ue_y, c='#3F3C3D', marker='+', linestyle='none', label='+')
    axs[1, 3].scatter(x_bs_200, y_bs_200, edgecolor='#00FA09', facecolor='none', alpha=0.9, label='o')
    axs[1, 3].scatter(x_bs_200_, y_bs_200_, edgecolor='#000FFA', facecolor='none', alpha=0.9, label='o')
    axs[1, 3].margins(0.05)
    axs[1, 3].set_title(r'\bf{200 BSs/$km^2$}', fontsize=12)
    axs[1, 3].set_xticks([])
    axs[1, 3].set_yticks([])
    # fig.savefig('poisson.eps', format='eps', dpi=1000)
    # fig.savefig('poisson.svg', format='svg', dpi=1200)
    # fig.savefig("poisson_homogeneous.pdf", bbox_inches='tight', dpi=400)

    plt.show()


# plot_two_operators_at_different_densities(500, 25, 200, area=AREA_OF_SIMULATION)



















