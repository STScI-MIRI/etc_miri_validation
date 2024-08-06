# import the necessary packages
import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table, vstack


def plot_comp_all(outpath, results):
    # plot the predicted SNR of the ETC vs the measured SNR
    plt.scatter(results["SNR_data"], results["SNR_ETC"], color="k")

    # plot the one-to-one line
    x = np.arange(np.min(results["SNR_data"]), np.max(results["SNR_data"]), 1)
    plt.plot(x, x, ls=":")

    # finalize and save the figure
    plt.xlabel("SNR ETC", fontsize=16)
    plt.ylabel("SNR data", fontsize=16)
    plt.savefig(outpath + "SNR_comp_all.pdf", bbox_inches="tight")

    # plot the predicted background of the ETC vs the measured background
    plt.clf()
    plt.scatter(results["bkg_data"], results["bkg_ETC"], color="k")

    # plot the one-to-one line
    x = np.arange(np.min(results["bkg_data"]), np.max(results["bkg_data"]), 1)
    plt.plot(x, x, ls=":")

    # finalize and save the figure
    plt.xlabel("background ETC (MJy/sr)", fontsize=16)
    plt.ylabel("background data (MJy/sr)", fontsize=16)
    plt.savefig(outpath + "bkg_comp_all.pdf", bbox_inches="tight")


def main():
    # define the path
    path = "/Users/mdecleir/Documents/MIRI_func/"

    # obtain the results
    resultsHD = Table.read(path + "HD 180609_comp.txt", format="ascii")
    resultsBD = Table.read(path + "BD+60 1753_comp.txt", format="ascii")
    results2M = Table.read(path + "2MASS J17430448+6655015_comp.txt", format="ascii")

    # combine all tables
    all_results = vstack([resultsHD, resultsBD, results2M])

    # compare SNR and background for all measurements in one plot
    plot_comp_all(path, all_results)


if __name__ == "__main__":
    main()
