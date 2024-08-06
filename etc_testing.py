# import the necessary packages
import json

import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table
from astropy.time import Time
from astropy import visualization
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pandeia.engine.perform_calculation import perform_calculation


def run_etc(datapath, outpath, star, filters, input_dic):
    """
    Function to run the ETC via Pandeia

    Parameters
    ----------
    datapath : string
        Path to the data files

    outpath : string
        Path to store the output table

    star : string
        Name of the star

    filters : list
        Filters

    input_dic : dict
        Input configuration

    Returns
    -------
    result_tab : astropy.table.table.Table
        Results of the ETC
    """
    # create a Table for the results
    result_tab = Table(
        names=("filter", "date", "SNR", "background"),
        dtype=("str", "str", "float64", "float64"),
    )

    # run over all filters
    for filter in filters:
        # obtain the table with the observation strategy
        filename = filter + "_eefrac0.7_phot.fits"
        phot_tab = Table.read(datapath + filename)
        observations = phot_tab[phot_tab["name"] == star]

        # skip this filter if there are no observations
        if len(observations) == 0:
            continue

        # run over all observations
        for obs in observations:
            # obtain the date of the observation
            time = Time(obs["timemid"], format="mjd")
            date = time.to_value("iso", subfmt="date")

            # change the filter
            input_dic["configuration"]["instrument"]["filter"] = filter.lower()

            # change the background
            back_table = Table.read(
                outpath
                + "etc_workbook_download/"
                + star
                + "/backgrounds/backgrounds_"
                + date
                + ".fits"
            )
            input_dic["background"] = [
                back_table["wavelength"],
                back_table["background"],
            ]

            # change the subarray and readout pattern
            input_dic["configuration"]["detector"]["subarray"] = obs["subarray"].lower()
            input_dic["configuration"]["detector"]["readout_pattern"] = obs[
                "readpattern"
            ].lower()

            # change the number of groups and integrations
            input_dic["configuration"]["detector"]["ngroup"] = int(obs["ngroups"])
            input_dic["configuration"]["detector"]["nint"] = int(obs["nints"])

            # change the aperture and sky annulus radii
            input_dic["strategy"]["aperture_size"] = obs["aprad"] * 0.11
            input_dic["strategy"]["sky_annulus"] = [
                obs["annrad1"] * 0.11,
                obs["annrad2"] * 0.11,
            ]
            print(
                date,
                input_dic["configuration"]["instrument"]["filter"],
                input_dic["configuration"]["detector"]["ngroup"],
                input_dic["configuration"]["detector"]["nint"],
                input_dic["configuration"]["detector"]["subarray"],
                input_dic["configuration"]["detector"]["readout_pattern"],
            )

            # run the ETC
            report = perform_calculation(input_dic)

            # add the results to the table
            result_tab.add_row(
                (
                    report["scalar"]["filter"],
                    date,
                    report["scalar"]["sn"],
                    report["scalar"]["background"],
                )
            )

    return result_tab


def comp_snr(datapath, outpath, star, filters, etc_res):
    """
    Function to compare the predicted and the measured SNR

    Parameters
    ----------
    datapath : string
        Path to the data files

    outpath : string
        Path to store the output table and plot

    star : string
        Name of the star

    filters : list
        Filters

    etc_res : astropy.table
        ETC predictions

    Returns
    -------
    comp_tab : astropy table
        Comparison between predicted and measured SNR
    """
    # create a figure
    fig, ax = plt.subplots()
    colors = plt.get_cmap("tab10")

    # create a table for the results
    result_tab = Table(
        names=("filter", "date", "SNR_ETC", "SNR_data", "SNR_diff=(data-ETC)/data"),
        dtype=("str", "str", "float64", "float64", "float64"),
    )

    # run over all filters
    for i, filter in enumerate(filters):
        # obtain the table with the photometry
        if filter == "F2550W":
            filename = filter + "_bkgsub_eefrac0.7_phot.fits"
        else:
            filename = filter + "_eefrac0.7_phot.fits"
        phot_tab = Table.read(datapath + filename)
        observations = phot_tab[phot_tab["name"] == star]

        # skip this filter if there are no observations
        if len(observations) == 0:
            continue

        # run over all observations
        for obs in observations:
            # obtain the date of the observation
            time = Time(obs["timemid"], format="mjd")
            date = time.to_value("iso", subfmt="date")

            # obtain the SNR of the observation
            snr_data = obs["aperture_sum_bkgsub"] / obs["aperture_sum_bkgsub_err"]

            # obtain the predicted SNR
            filt_mask = etc_res["filter"] == filter.lower()
            date_mask = etc_res["date"] == date
            snr_etc = etc_res[filt_mask & date_mask]["SNR"][0]

            # calculate the difference between the measured and predicted SNR
            diff = (snr_data - snr_etc) / snr_data

            # save the results to the table
            result_tab.add_row(
                (filter, date, "%.2f" % snr_etc, "%.2f" % snr_data, "%.2f" % diff)
            )

            # plot the measured and predicted SNR vs. time
            with visualization.time_support(format="iso"):
                plt.scatter(
                    time,
                    snr_data,
                    marker="o",
                    color=colors(i % 10),
                    facecolor="none",
                    lw=1.5,
                    s=40,
                    alpha=0.9,
                )
                plt.scatter(
                    time,
                    snr_etc,
                    marker="x",
                    color=colors(i % 10),
                    lw=1.5,
                    s=40,
                    alpha=0.9,
                )

    # finalize and save the figure
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
    xticks = xtick_labels.to_value("mjd")
    plt.xticks(xticks, xtick_labels, rotation=27)
    plt.xlabel("date", fontsize=16)
    plt.ylabel("SNR", fontsize=16)

    handle1 = Line2D(
        [],
        [],
        lw=0,
        color="grey",
        marker="o",
        fillstyle="none",
        markersize=7,
        alpha=0.9,
    )
    handle2 = Line2D(
        [],
        [],
        lw=0,
        color="grey",
        marker="x",
        markersize=7,
        alpha=0.9,
    )
    plt.figlegend(
        handles=[handle1, handle2],
        labels=["data", "ETC"],
        bbox_to_anchor=(0.9, 0.25),
    )
    handles = []
    for j in range(len(filters)):
        handles.append(Patch(facecolor=colors(j % 10), alpha=0.9))
    plt.figlegend(handles=handles, labels=filters, bbox_to_anchor=(0.9, 0.87))

    plt.savefig(outpath + star + "_SNR.pdf", bbox_inches="tight")

    return result_tab


def comp_back(datapath, outpath, star, filters, comp_tab, etc_res):
    """
    Function to compare the predicted and the measured background

    Parameters
    ----------
    datapath : string
        Path to the data files

    outpath : string
        Path to store the output table and plot

    star : string
        Name of the star

    filters : list
        Filters

    comp_tab : astropy table
        Comparison between predicted and measured SNR

    etc_res : astropy.table
        ETC predictions

    Returns
    -------
    comp_tab : astropy table
        Comparison between predicted and measured SNR and background
    """
    # create a figure
    fig, ax = plt.subplots()
    colors = plt.get_cmap("tab10")

    # open the file with the time dependent calibration factors
    cal_tab = Table.read(datapath + "jwst_miri_photom_coeff.dat", format="ascii")

    # add columns to the comparison table
    comp_tab["bkg_ETC"] = np.full(len(comp_tab), np.nan)
    comp_tab["bkg_data"] = np.full(len(comp_tab), np.nan)
    comp_tab["bkg_diff=(data-ETC)/data"] = np.full(len(comp_tab), np.nan)

    # run over all filters
    for i, filter in enumerate(filters):
        # obtain the table with the photometry
        filename = filter + "_eefrac0.7_phot.fits"
        phot_tab = Table.read(datapath + filename)
        observations = phot_tab[phot_tab["name"] == star]

        # skip this filter if there are no observations
        if len(observations) == 0:
            continue

        # run over all observations
        for obs in observations:
            # obtain the date of the observation
            time = Time(obs["timemid"], format="mjd")
            date = time.to_value("iso", subfmt="date")

            # obtain the background of the observations
            bkg_DNs = obs["mean_bkg"]

            # obtain the subarray of the observations
            sub = obs["subarray"]

            # calculate the calibration factor (Gordon+2024)
            # CF(t) = {A + B exp[−(t−to)/τ]} D(SA)
            # A, B and τ in table 4, Gordon+2024
            # D(SA) in table 3, Gordon+2024
            # to = 59720 d
            if sub == "FULL" or sub == "BRIGHTSKY":
                D = 1
            elif sub == "SUB256":
                D = 0.985
            elif sub == "SUB128":
                D = 0.977
            elif sub == "SUB64":
                D = 0.969
            else:
                print("Unknown subarray, please check.")
            cal_fil = cal_tab[cal_tab["filter"] == filter]
            CF = (
                cal_fil["photmjysr"]
                + cal_fil["amplitude"] * np.exp(-(time.value - 59720) / cal_fil["tau"])
            ) * D

            # convert the units from DN/s/pix to MJy/sr
            bkg_MJysr = bkg_DNs * CF[0]

            # obtain the predicted background
            filt_mask = etc_res["filter"] == filter.lower()
            date_mask = etc_res["date"] == date
            bkg_etc = etc_res[filt_mask & date_mask]["background"][0]

            # calculate the difference between the measured and predicted background
            diff = (bkg_MJysr - bkg_etc) / bkg_MJysr

            # save the results to the table
            filt_mask = comp_tab["filter"] == filter
            date_mask = comp_tab["date"] == date
            comp_tab["bkg_ETC"][filt_mask & date_mask] = "%.2f" % bkg_etc
            comp_tab["bkg_data"][filt_mask & date_mask] = "%.2f" % bkg_MJysr
            comp_tab["bkg_diff=(data-ETC)/data"][filt_mask & date_mask] = "%.2f" % diff

            # plot the measured and predicted background vs. time
            with visualization.time_support(format="iso"):
                plt.scatter(
                    time,
                    bkg_MJysr,
                    marker="o",
                    color=colors(i % 10),
                    facecolor="none",
                    lw=1.5,
                    s=40,
                    alpha=0.9,
                )
                plt.scatter(
                    time,
                    bkg_etc,
                    marker="x",
                    color=colors(i % 10),
                    lw=1.5,
                    s=40,
                    alpha=0.9,
                )

    # finalize and save the figure
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
    xticks = xtick_labels.to_value("mjd")
    plt.xticks(xticks, xtick_labels, rotation=27)
    plt.xlabel("date", fontsize=16)
    plt.ylabel("background (MJy/sr)", fontsize=16)

    handle1 = Line2D(
        [],
        [],
        lw=0,
        color="grey",
        marker="o",
        fillstyle="none",
        markersize=7,
        alpha=0.9,
    )
    handle2 = Line2D(
        [],
        [],
        lw=0,
        color="grey",
        marker="x",
        markersize=7,
        alpha=0.9,
    )
    plt.figlegend(
        handles=[handle1, handle2],
        labels=["data", "ETC"],
        bbox_to_anchor=(0.9, 0.25),
    )
    handles = []
    for j in range(len(filters)):
        handles.append(Patch(facecolor=colors(j % 10), alpha=0.9))
    plt.figlegend(handles=handles, labels=filters, bbox_to_anchor=(0.9, 0.87))

    plt.savefig(outpath + star + "_bkg.pdf", bbox_inches="tight")

    return comp_tab


def main():
    # define the path
    path = "/Users/mdecleir/Documents/MIRI_func/"

    # define the star
    star = "BD+60 1753"
    # star = "HD 180609"
    # star = "2MASS J17430448+6655015"

    # open the original ETC input file (for filter F560W)
    with open(path + "etc_workbook_download/" + star + "/input.json", "r") as infile:
        input_dic = json.loads(infile.read())

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

    # run the ETC for all filters
    etc_results = run_etc(path + "Absfluxcal/", path, star, filters, input_dic)

    # compare the predicted SNR from the ETC to the observed SNR
    comp_table = comp_snr(path + "Absfluxcal/", path, star, filters, etc_results)

    # compare the predicted background from the ETC to the observed background
    comp_table = comp_back(
        path + "Absfluxcal/", path, star, filters, comp_table, etc_results
    )

    # write the results table to a file
    comp_table.write(
        path + star + "_comp.txt",
        format="ascii",
        overwrite=True,
    )


if __name__ == "__main__":
    main()
