# import the necessary packages
import json

import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table
from astropy.time import Time
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
        names=("filter", "SNR", "background"),
        dtype=("str", "float64", "float64"),
    )

    for filter in filters:
        # obtain the table with the observation strategy
        filename = filter + "_bkgsub_eefrac0.7_phot.fits"
        phot_tab = Table.read(datapath + filename)
        star_mask = phot_tab["name"] == star

        # skip this filter if there are no observations
        if np.sum(star_mask) == 0:
            continue

        # obtain the MJD times of the observations and remove commissioning data
        times = phot_tab[star_mask]["timemid"]
        time_mask = times > 59762.2

        # obtain the subarray
        subs = phot_tab[star_mask][time_mask]["subarray"]
        if not np.all(subs == subs[0]):
            print(
                "Subarray for this star is not the same between different observations, please check."
            )
        else:
            sub = subs[0].lower()

        # obtain the number of groups
        ngroups = phot_tab[star_mask][time_mask]["ngroups"]
        if not np.all(ngroups == ngroups[0]):
            print(
                "Number of groups for this star is not the same between different observations, please check."
            )
        else:
            ngroup = int(ngroups[0])

        # obtain the aperture radius and the sky annulus radii
        ap_rads = phot_tab[star_mask][time_mask]["aprad"] * 0.11
        in_ann_rads = phot_tab[star_mask][time_mask]["annrad1"] * 0.11
        out_ann_rads = phot_tab[star_mask][time_mask]["annrad2"] * 0.11
        if (
            not np.all(ap_rads == ap_rads[0])
            or not np.all(in_ann_rads == in_ann_rads[0])
            or not np.all(out_ann_rads == out_ann_rads[0])
        ):
            print(
                "Aperture or sky annulus radii for this star are not the same between different observations, please check."
            )
        else:
            aprad = ap_rads[0]
            in_ann_rad = in_ann_rads[0]
            out_ann_rad = out_ann_rads[0]

        # change the filter
        input_dic["configuration"]["instrument"]["filter"] = filter.lower()

        # change the subarray
        input_dic["configuration"]["detector"]["subarray"] = sub

        # change the number of groups and integrations
        input_dic["configuration"]["detector"]["ngroup"] = ngroup
        if star == "BD+60 1753" and filter == "F2550W":
            input_dic["configuration"]["detector"]["nint"] = 10
        elif star == "HD 180609" and filter == "F560W":
            input_dic["configuration"]["detector"]["nint"] = 2
        elif star == "2MASS J17430448+6655015" and filter == "F1500W":
            input_dic["configuration"]["detector"]["nint"] = 2
        else:
            input_dic["configuration"]["detector"]["nint"] = 1

        # change the aperture and sky annulus radii
        input_dic["strategy"]["aperture_size"] = aprad
        input_dic["strategy"]["sky_annulus"] = [in_ann_rad, out_ann_rad]
        print(
            filter,
            ngroup,
            input_dic["configuration"]["detector"]["nint"],
            input_dic["configuration"]["detector"]["subarray"],
        )

        # run the ETC
        report = perform_calculation(input_dic)

        # add the results to the table
        result_tab.add_row(
            (
                report["scalar"]["filter"],
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

    # create a table for the results
    result_tab = Table(
        names=("filter", "SNR_ETC", "SNR_data", "SNR_diff=(data-ETC)/data"),
        dtype=("str", "float64", "float64", "float64"),
    )

    for i, filter in enumerate(filters):
        # obtain the table with the photometry
        if filter == "F2550W":
            filename = filter + "_bkgsub_eefrac0.7_phot.fits"
        else:
            filename = filter + "_eefrac0.7_phot.fits"
        phot_tab = Table.read(datapath + filename)
        star_mask = phot_tab["name"] == star

        # skip this filter if there are no observations
        if np.sum(star_mask) == 0:
            continue

        # obtain the MJD times of the observations and remove commissioning data
        times = phot_tab[star_mask]["timemid"]
        time_mask = times > 59762.2

        # obtain the SNRs of the observations and calculate the mean
        snrs = (
            phot_tab[star_mask]["aperture_sum_bkgsub"][time_mask]
            / phot_tab[star_mask]["aperture_sum_bkgsub_err"][time_mask]
        )
        mean_snr = np.mean(snrs)

        # obtain the predicted SNR
        etc_snr = etc_res[etc_res["filter"] == filter.lower()]["SNR"][0]

        # calculate the difference between the measured and predicted SNR
        diff = (mean_snr - etc_snr) / mean_snr

        # save the results to the table
        result_tab.add_row((filter, "%.2f" % etc_snr, "%.2f" % mean_snr, "%.2f" % diff))

        # plot the measured and predicted SNR vs. time
        sc = plt.scatter(
            times[time_mask],
            snrs,
            label=filter + ": " + "{:.2f}".format(diff),
        )
        color = sc.get_facecolor()
        plt.axhline(mean_snr, color=color, ls=":")
        plt.axhline(etc_snr, color=color)

    # finalize and save the figure
    xticks = [59700, 59800, 59900, 60000, 60100, 60200, 60300, 60400]
    xtick_times = Time(xticks, format="mjd")
    xtick_labels = xtick_times.to_value("iso", subfmt="date")
    plt.xticks(xticks, xtick_labels, rotation=27)
    plt.axvline(60128, color="k")
    plt.text(0.1, 0.7, "cycle 1", fontsize=18, transform=ax.transAxes)
    plt.text(0.7, 0.7, "cycle 2", fontsize=18, transform=ax.transAxes)
    plt.xlabel("date")
    plt.ylabel("SNR")
    plt.legend(loc="upper center", ncols=3)
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

    # open the file with the time dependent calibration factors
    cal_tab = Table.read(datapath + "jwst_miri_photom_coeff.dat", format="ascii")

    # add columns to the comparison table
    comp_tab["bkg_ETC"] = np.full(len(comp_tab), np.nan)
    comp_tab["bkg_data"] = np.full(len(comp_tab), np.nan)
    comp_tab["bkg_diff=(data-ETC)/data"] = np.full(len(comp_tab), np.nan)

    for i, filter in enumerate(filters):
        # obtain the table with the photometry
        filename = filter + "_eefrac0.7_phot.fits"
        phot_tab = Table.read(datapath + filename)
        star_mask = phot_tab["name"] == star

        # skip this filter if there are no observations
        if np.sum(star_mask) == 0:
            continue

        # obtain the MJD times of the observations and remove commissioning data
        times = phot_tab[star_mask]["timemid"]
        time_mask = times > 59762.2

        # obtain the background of the observations
        bkgs_DNs = phot_tab[star_mask][time_mask]["mean_bkg"]

        # obtain the subarray of the observations
        subs = phot_tab[star_mask][time_mask]["subarray"]
        if not np.all(subs == subs[0]):
            print(
                "Subarray for this star is not the same between different observations, please check."
            )
        else:
            sub = subs[0]

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
            + cal_fil["amplitude"]
            * np.exp(-(times[time_mask] - 59720) / cal_fil["tau"])
        ) * D

        # convert the units from DN/s/pix to MJy/sr and calculate the mean background
        bkgs_MJysr = bkgs_DNs * CF
        mean_bkg = np.mean(bkgs_MJysr)

        # obtain the predicted background
        etc_bkg = etc_res[etc_res["filter"] == filter.lower()]["background"][0]

        # calculate the difference between the measured and predicted background
        diff = (mean_bkg - etc_bkg) / mean_bkg

        # save the results to the table
        comp_tab["bkg_ETC"][comp_tab["filter"] == filter] = "%.2f" % etc_bkg
        comp_tab["bkg_data"][comp_tab["filter"] == filter] = "%.2f" % mean_bkg
        comp_tab["bkg_diff=(data-ETC)/data"][comp_tab["filter"] == filter] = (
            "%.2f" % diff
        )

        # plot the measured and predicted background vs. time
        sc = plt.scatter(
            times[time_mask],
            bkgs_MJysr,
            label=filter + ": " + "{:.2f}".format(diff),
        )
        color = sc.get_facecolor()
        plt.axhline(mean_bkg, color=color, ls=":")
        plt.axhline(etc_bkg, color=color)

    # finalize and save the figure
    xticks = [59700, 59800, 59900, 60000, 60100, 60200, 60300, 60400]
    xtick_times = Time(xticks, format="mjd")
    xtick_labels = xtick_times.to_value("iso", subfmt="date")
    plt.xticks(xticks, xtick_labels, rotation=20)
    plt.axvline(60128, color="k")
    plt.text(0.1, 0.7, "cycle 1", fontsize=18, transform=ax.transAxes)
    plt.text(0.7, 0.7, "cycle 2", fontsize=18, transform=ax.transAxes)
    plt.xlabel("date")
    plt.ylabel("background (MJy/sr)")
    plt.legend(loc="upper center", ncols=3)
    plt.savefig(outpath + star + "_bkg.pdf", bbox_inches="tight")

    return comp_tab


def main():
    # define the path
    path = "/Users/mdecleir/Documents/MIRI_func/"

    # define the star
    # star = "BD+60 1753"
    # star = "HD 180609"
    star = "2MASS J17430448+6655015"

    # open the original ETC input file (for BD+60 1753, filter F560W)
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
