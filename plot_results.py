# import the necessary packages
import matplotlib.pyplot as plt
import numpy as np

from astropy import visualization
from astropy.table import Table, vstack
from astropy.time import Time


def plot_comp_all(outpath, stars, results):
    """
    Function to compare the SNR and background for all stars and all filters in one plot

    Parameters
    ----------
    outpath : string
        Path to store the output plots

    stars : list
        Star names

    results : astropy Table
        Data to plot

    Returns
    -------
    Plots with results for all stars and all filters:
    - for SNR
    - for background
    """
    # create a figure
    fig = plt.figure(figsize=(6, 7))
    colors = plt.get_cmap("tab10")
    markers = ["o", "X", "v"]

    for i, star in enumerate(stars):
        star_mask = results["star"] == star
        # plot the predicted SNR of the ETC vs. the measured SNR
        plt.scatter(
            results[star_mask]["SNR_data"],
            results[star_mask]["SNR_ETC"],
            color=colors(i % 3),
            marker=markers[i],
            s=15,
            alpha=0.9,
            label=star,
        )

        # print general results
        print(
            star,
            ": Average SNR relative difference = |(data-ETC)|/data =",
            np.mean(np.abs(results[star_mask]["SNR_diff=(data-ETC)/data"])),
        )

    # plot the one-to-one line
    x = np.arange(np.min(results["SNR_data"]), np.max(results["SNR_data"]), 1)
    plt.plot(x, x, ls="-", c="k", alpha=0.9, label="S/N ETC = S/N data")
    plt.plot(
        x, x * 1.1, ls=":", c="purple", alpha=0.7, label="S/N ETC = 1.1 * S/N data"
    )
    plt.plot(x, x * 0.9, ls=":", c="red", alpha=0.7, label="S/N ETC = 0.9 * S/N data")
    print(
        "Percentage of data points above the one-to-one line for SNR",
        np.sum(results["SNR_data"] < results["SNR_ETC"]) / len(results["SNR_data"]),
    )

    # finalize and save the figure
    plt.xlabel("S/N data", fontsize=16)
    plt.ylabel("S/N ETC", fontsize=16)
    plt.legend()
    plt.savefig(outpath + "SNR_comp_all.pdf", bbox_inches="tight")

    # plot the predicted background of the ETC vs. the measured background
    fig = plt.figure(figsize=(6, 7))
    for i, star in enumerate(stars):
        star_mask = results["star"] == star
        plt.scatter(
            results[star_mask]["bkg_data"],
            results[star_mask]["bkg_ETC"],
            color=colors(i % 3),
            marker=markers[i],
            s=15,
            alpha=0.9,
            label=star,
        )

        # print general results
        print(
            star,
            ": Average background relative difference = |(data-ETC)|/data =",
            np.mean(np.abs(results[star_mask]["bkg_diff=(data-ETC)/data"])),
        )

    # plot the one-to-one line
    x = np.arange(np.min(results["bkg_data"]), np.max(results["bkg_data"]), 1)
    plt.plot(x, x, ls="-", c="k", alpha=0.9, label="bkg ETC = bkg data")
    plt.plot(
        x, x * 1.1, ls=":", c="purple", alpha=0.7, label="bkg ETC = 1.1 * bkg data"
    )
    plt.plot(x, x * 0.9, ls=":", c="red", alpha=0.7, label="bkg ETC = 0.9 * bkg data")
    print(
        "Percentage of data points under the one-to-one line for background",
        np.sum(results["bkg_data"] > results["bkg_ETC"]) / len(results["bkg_data"]),
    )

    # finalize and save the figure
    plt.xlabel("background data (MJy/sr)", fontsize=16)
    plt.ylabel("background ETC (MJy/sr)", fontsize=16)
    plt.legend()
    plt.savefig(outpath + "bkg_comp_all.pdf", bbox_inches="tight")


def plot_comp_filter(outpath, stars, results, filters):
    """
    Function to plot the difference in SNR and background vs. the filter

    Parameters
    ----------
    outpath : string
        Path to store the output plots

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
    colors = plt.get_cmap("tab10")
    markers = ["o", "X", "v"]

    for i, star in enumerate(stars):
        star_mask = results["star"] == star
        # plot SNR difference vs. filter
        ax1.scatter(
            results[star_mask]["filter"],
            results[star_mask]["SNR_diff=(data-ETC)/data"],
            color=colors(i % 3),
            marker=markers[i],
            s=15,
            alpha=0.9,
            label=star,
        )
        # plot background difference vs. filter
        ax2.scatter(
            results[star_mask]["filter"],
            results[star_mask]["bkg_diff=(data-ETC)/data"],
            color=colors(i % 3),
            marker=markers[i],
            s=15,
            alpha=0.9,
            label=star,
        )

    # finalize and save the plots
    ax1.axhline(ls=":", color="k")
    ax1.set_xlabel("filter", fontsize=15)
    ax1.set_ylabel("S/N diff. = (data-ETC)/data", fontsize=15)
    ax1.legend()
    fig1.savefig(outpath + "SNR_filter.pdf", bbox_inches="tight")

    ax2.axhline(ls=":", color="k")
    ax2.set_xlabel("filter", fontsize=15)
    ax2.set_ylabel("bkg diff. = (data-ETC)/data", fontsize=15)
    ax2.legend()
    fig2.savefig(outpath + "bkg_filter.pdf", bbox_inches="tight")


def plot_comp_time(outpath, stars, results):
    """
    Function to plot the difference in SNR and background vs. time

    Parameters
    ----------
    outpath : string
        Path to store the output plots

    stars : list
        Star names

    results : astropy Table
        Data to plot

    Returns
    -------
    Plots with fractional difference vs. time:
    - for SNR
    - for background
    """
    # create the figures
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    colors = plt.get_cmap("tab10")
    markers = ["o", "X", "v"]

    for i, star in enumerate(stars):
        star_mask = results["star"] == star
        dates = Time(results[star_mask]["date"], format="iso")
        with visualization.time_support(format="iso"):
            # plot SNR difference vs. time
            ax1.scatter(
                dates,
                results[star_mask]["SNR_diff=(data-ETC)/data"],
                color=colors(i % 3),
                marker=markers[i],
                s=15,
                alpha=0.9,
                label=star,
            )
            # plot background difference vs. time
            ax2.scatter(
                dates,
                results[star_mask]["bkg_diff=(data-ETC)/data"],
                color=colors(i % 3),
                marker=markers[i],
                s=15,
                alpha=0.9,
                label=star,
            )

    # finalize and save the plots
    xtick_labels = Time(
        [
            "2022-05-01",
            "2022-08-01",
            "2022-11-01",
            "2023-02-01",
            "2023-05-01",
            "2023-08-01",
            "2023-11-01",
            "2024-02-01",
            "2024-05-01",
            "2024-08-01",
        ],
        format="iso",
        out_subfmt="date",
    )
    ax1.axhline(ls=":", color="k")
    ax1.set_xlabel("date", fontsize=15)
    ax1.set_ylabel("S/N diff. = (data-ETC)/data", fontsize=15)
    ax1.set_xticks(xtick_labels)
    ax1.set_xticklabels(xtick_labels)
    ax1.tick_params(axis="x", labelrotation=27)
    ax1.legend()
    fig1.savefig(outpath + "SNR_time.pdf", bbox_inches="tight")

    ax2.axhline(ls=":", color="k")
    ax2.set_xlabel("date", fontsize=15)
    ax2.set_ylabel("bkg diff. = (data-ETC)/data", fontsize=15)
    ax2.set_xticks(xtick_labels)
    ax2.set_xticklabels(xtick_labels)
    ax2.tick_params(axis="x", labelrotation=27)
    ax2.legend()
    fig2.savefig(outpath + "bkg_time.pdf", bbox_inches="tight")


def main():
    # define the path
    path = "/Users/mdecleir/Documents/MIRI_func/"

    # obtain the results
    resultsBD = Table.read(path + "BD+60 1753_comp.txt", format="ascii")
    resultsBD["star"] = "BD+60 1753"
    resultsHD = Table.read(path + "HD 180609_comp.txt", format="ascii")
    resultsHD["star"] = "HD 180609"
    results2M = Table.read(path + "2MASS J17430448+6655015_comp.txt", format="ascii")
    results2M["star"] = "2MASS J17430448+6655015"

    # combine all tables
    all_results = vstack([resultsHD, resultsBD, results2M])

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

    # compare SNR and background for all measurements in one plot
    plot_comp_all(path, stars, all_results)

    # compare SNR and background difference per filter
    plot_comp_filter(path, stars, all_results, filters)

    # compare SNR and background difference over time
    plot_comp_time(path, stars, all_results)


if __name__ == "__main__":
    main()
