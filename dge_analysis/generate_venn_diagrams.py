#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn3_unweighted

"""
3-circle Venn diagram generation functions.
"""


#-------# Functions Definitions #----------------------------------------------#

def generate_venn3_diagram(
        set1,
        set2,
        set3,
        plot_formats,
        title,
        file_name,
        ):
    """
    Generate a Venn diagram from 3 given sets and their respective names.
    """

    set_array = [
            set1[0],
            set2[0],
            set3[0],
            ]

    set_names = [
            set1[1],
            set2[1],
            set3[1],
            ]

    venn3(set_array[0:3], set_names[0:3])

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{file_name}.{format}")

    # Clears current figure
    plt.clf()

    venn3_unweighted(set_array[0:3], set_names[0:3])

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{file_name}_unweighted.{format}")

    # Clears current figure
    plt.clf()


def generate_venn3_diagram_with_regulation_labels(
        set1,
        set2,
        set3,
        up_regulation_labels,
        down_regulation_labels,
        plot_formats,
        title,
        file_name,
        ):
    """
    Generates a venn diagram of 3 sets with custom labels.
    Labels are represented binary numbers that represent each intersection:
        - First bit equals to the first set.
        - Second bit equals to the second set.
        - Third bit equals to the third set.
    """

    set_list = [
            set1[0],
            set2[0],
            set3[0],
            ]

    set_names_list = [
            set1[1],
            set2[1],
            set3[1],
            ]

    # Preparing the lables.
    label_100 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['100'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['100'])}"
            )
    label_110 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['110'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['110'])}"
            )
    label_101 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['101'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['101'])}"
            )
    label_111 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['111'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['111'])}"
            )
    label_010 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['010'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['010'])}"
            )
    label_011 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['011'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['011'])}"
            )
    label_001 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['001'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['001'])}"
            )

    venn = venn3(set_list[0:3], set_names_list[0:3])
    # Setting the labels
    venn.get_label_by_id('100').set_text(label_100)
    venn.get_label_by_id('110').set_text(label_110)
    venn.get_label_by_id('111').set_text(label_111)
    venn.get_label_by_id('011').set_text(label_011)
    venn.get_label_by_id('001').set_text(label_001)
    venn.get_label_by_id('101').set_text(label_101)
    venn.get_label_by_id('010').set_text(label_010)

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{file_name}.{format}")

    # Clears current figure
    plt.clf()

    venn_unw = venn3_unweighted(set_list[0:3], set_names_list[0:3])
    venn_unw.get_label_by_id('100').set_text(label_100)
    venn_unw.get_label_by_id('110').set_text(label_110)
    venn_unw.get_label_by_id('111').set_text(label_111)
    venn_unw.get_label_by_id('011').set_text(label_011)
    venn_unw.get_label_by_id('001').set_text(label_001)
    venn_unw.get_label_by_id('101').set_text(label_101)
    venn_unw.get_label_by_id('010').set_text(label_010)

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{file_name}_unweighted.{format}")

    # Clears current figure
    plt.clf()

