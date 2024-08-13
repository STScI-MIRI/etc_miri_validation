# import the necessary packages
import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table, vstack


def plot_comp_all(outpath, results):
    """
    Function to compare the SNR and background for all stars and all filters in one plot

    Parameters
    ----------
    outpath : string
        Path to store the output plot

    results : astropy Table
        Data to plot

    Returns
    -------
    Plots with results for all stars and all filters:
    - for SNR
    - for background
    """
    # plot the predicted SNR of the ETC vs. the measured SNR
    plt.scatter(results["SNR_data"], results["SNR_ETC"], color="k")

    # plot the one-to-one line
    x = np.arange(np.min(results["SNR_data"]), np.max(results["SNR_data"]), 1)
    plt.plot(x, x, ls=":")

    # finalize and save the figure
    plt.xlabel("SNR data", fontsize=16)
    plt.ylabel("SNR ETC", fontsize=16)
    plt.savefig(outpath + "SNR_comp_all.pdf", bbox_inches="tight")

    # plot the predicted background of the ETC vs. the measured background
    plt.clf()
    plt.scatter(results["bkg_data"], results["bkg_ETC"], color="k")

    # plot the one-to-one line
    x = np.arange(np.min(results["bkg_data"]), np.max(results["bkg_data"]), 1)
    plt.plot(x, x, ls=":")

    # finalize and save the figure
    plt.xlabel("background data (MJy/sr)", fontsize=16)
    plt.ylabel("background ETC (MJy/sr)", fontsize=16)
    plt.savefig(outpath + "bkg_comp_all.pdf", bbox_inches="tight")


def plot_comp_filter(outpath, stars, results, filters):
    """
    Function to plot the difference in SNR and background vs. the filter

    Parameters
    ----------
    outpath : string
        Path to store the output plot

    stars : list
        Star names

    results : astropy Table
        Data to plot

    filters : list
        Filters

    Returns
    -------
    Plots with fractional difference vs. filter:
    - for SNR
    - for background
    """
    # create the figures
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    colors = ["blue", "orange", "green"]

    for i, star in enumerate(stars):
        star_mask = results["star"] == star
        # plot SNR difference vs. filter
        ax1.scatter(
            results[star_mask]["filter"],
            results[star_mask]["SNR_diff=(data-ETC)/data"],
            color=colors[i],
            marker="x",
            alpha=0.9,
            s=60,
            label=star,
        )
        # plot background difference vs. filter
        ax2.scatter(
            results[star_mask]["filter"],
            results[star_mask]["bkg_diff=(data-ETC)/data"],
            color=colors[i],
            marker="x",
            alpha=0.9,
            s=60,
            label=star,
        )

    # finalize and save the plots
    ax1.axhline(ls=":", color="k")
    ax1.set_xlabel("filter", fontsize=15)
    ax1.set_ylabel("SNR diff. = (data-ETC)/data", fontsize=15)
    ax1.legend()
    fig1.savefig(outpath + "SNR_filter.pdf", bbox_inches="tight")

    ax2.axhline(ls=":", color="k")
    ax2.set_xlabel("filter", fontsize=15)
    ax2.set_ylabel("bkg diff. = (data-ETC)/data", fontsize=15)
    ax2.legend()
    fig2.savefig(outpath + "bkg_filter.pdf", bbox_inches="tight")


def main():
    # define the path
    path = "/Users/mdecleir/Documents/MIRI_func/"

    # obtain the results
    resultsHD = Table.read(path + "HD 180609_comp.txt", format="ascii")
    resultsHD["star"] = "HD 180609"
    resultsBD = Table.read(path + "BD+60 1753_comp.txt", format="ascii")
    resultsBD["star"] = "BD+60 1753"
    results2M = Table.read(path + "2MASS J17430448+6655015_comp.txt", format="ascii")
    results2M["star"] = "2MASS J17430448+6655015"

    # combine all tables
    all_results = vstack([resultsHD, resultsBD, results2M])

    # compare SNR and background for all measurements in one plot
    plot_comp_all(path, all_results)

    # define the stars
    stars = ["BD+60 1753", "HD 180609", "2MASS J17430448+6655015"]

    # define the filters
    filters = [
        "F560W",
        "F770W",
        "F1000W",
        "F1130W",
        "F1280W",
        "F1500W",
        "F1800W",
        "F2100W",
        "F2550W",
    ]

    # compare SNR per filter
    plot_comp_filter(path, stars, all_results, filters)


if __name__ == "__main__":
    main()
