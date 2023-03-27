#!/usr/bin/env python3

"""3-circle Venn diagram generation functions. Uses matplotlib and
matplotlib-venn
"""

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_unweighted, venn3, venn3_unweighted
from venn import venn, draw_venn, generate_petal_labels, generate_colors


def generate_venn2_diagram(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        plot_formats,
        title,
        path,
        ):
    """Generate a Venn diagram from 2 given sets and their respective names.
    """

    plt.style.use(['default'])

    set_array = [mutant1_gene_set, mutant2_gene_set]

    set_names = [mutant1_name, mutant2_name]

    venn2(set_array[0:2], set_names[0:2])

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{path}.{format}", dpi=300)

        if format == "png":
            plt.savefig(f"{path}_transparent-bg.{format}", dpi=300, transparent=True)

    plt.close()

    venn2_unweighted(set_array[0:3], set_names[0:3])

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{path}_unweighted.{format}", dpi=300)

        if format == "png":
            plt.savefig(
                    f"{path}_unweighted_transparent-bg.{format}",
                    dpi=300,
                    transparent=True
                    )

    plt.close()


def generate_venn2_diagram_with_regulation_labels(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        up_regulation_labels,
        down_regulation_labels,
        plot_formats,
        title,
        file_path,
        ):
    """Generates a venn diagram of 2 sets with custom labels.
    Labels are represented binary numbers that represent each intersection:
        - First bit equals to the first set.
        - Second bit equals to the second set.
    """

    plt.style.use(['default'])

    set_array = [mutant1_gene_set, mutant2_gene_set]

    set_names = [mutant1_name, mutant2_name]

    # Preparing the lables.
    label_10 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['10'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['10'])}"
            )
    label_11 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['11'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['11'])}"
            )
    label_01 = (
            r"$\uparrow$"
                +
            f"{len(up_regulation_labels['01'])}\n"
                +
            r"$\downarrow$"
                +
            f"{len(down_regulation_labels['01'])}"
            )

    venn = venn2(set_array, set_names)
    # Setting the labels
    venn.get_label_by_id('10').set_text(label_10)
    venn.get_label_by_id('11').set_text(label_11)
    venn.get_label_by_id('01').set_text(label_01)

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{file_path}.{format}", dpi=300)

        if format == "png":
            plt.savefig(
                    f"{file_path}_transparent-bg.{format}",
                    dpi=300,
                    transparent=True
                    )

    # Clears current figure
    plt.close()

    venn_unw = venn2_unweighted(set_array, set_names)
    venn_unw.get_label_by_id('01').set_text(label_01)
    venn_unw.get_label_by_id('10').set_text(label_10)
    venn_unw.get_label_by_id('11').set_text(label_11)

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{file_path}_unweighted.{format}", dpi=300)

        if format == "png":
            plt.savefig(
                    f"{file_path}_unweighted_transparent-bg.{format}",
                    dpi=300,
                    transparent=True
                    )

    plt.close()

def generate_venn3_diagram(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        plot_formats,
        title,
        path,
        ):
    """Generate a Venn diagram from 3 given sets and their respective names.
    """

    plt.style.use(['default'])

    set_array = [mutant1_gene_set, mutant2_gene_set, mutant3_gene_set]

    set_names = [mutant1_name, mutant2_name, mutant3_name]

    venn3(set_array[0:3], set_names[0:3])

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{path}.{format}", dpi=300)

        if format == "png":
            plt.savefig(f"{path}_transparent-bg.{format}", dpi=300, transparent=True)

    # Clears current figure
    plt.clf()

    venn3_unweighted(set_array[0:3], set_names[0:3])

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{path}_unweighted.{format}", dpi=300)

        if format == "png":
            plt.savefig(
                    f"{path}_unweighted_transparent-bg.{format}",
                    dpi=300,
                    transparent=True
                    )

    plt.close()


def generate_venn3_diagram_with_regulation_labels(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        up_regulation_labels,
        down_regulation_labels,
        plot_formats,
        title,
        file_path,
        ):
    """Generates a venn diagram of 3 sets with custom labels.
    Labels are represented binary numbers that represent each intersection:
        - First bit equals to the first set.
        - Second bit equals to the second set.
        - Third bit equals to the third set.
    """

    plt.style.use(['default'])

    set_array = [mutant1_gene_set, mutant2_gene_set, mutant3_gene_set]

    set_names = [mutant1_name, mutant2_name, mutant3_name]

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

    venn = venn3(set_array, set_names)
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
        plt.savefig(f"{file_path}_labels.{format}", dpi=300)

        if format == "png":
            plt.savefig(
                    f"{file_path}_transparent-bg.{format}",
                    dpi=300,
                    transparent=True
                    )

    # Clears current figure
    plt.close()

    venn_unw = venn3_unweighted(set_array, set_names)
    venn_unw.get_label_by_id('100').set_text(label_100)
    venn_unw.get_label_by_id('110').set_text(label_110)
    venn_unw.get_label_by_id('111').set_text(label_111)
    venn_unw.get_label_by_id('011').set_text(label_011)
    venn_unw.get_label_by_id('001').set_text(label_001)
    venn_unw.get_label_by_id('101').set_text(label_101)
    venn_unw.get_label_by_id('010').set_text(label_010)

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{file_path}.{format}", dpi=300)

        if format == "png":
            plt.savefig(
                    f"{file_path}_unweigted_transparent-bg.{format}",
                    dpi=300,
                    transparent=True
                    )

    plt.close()


def generate_venn4_diagram(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        mutant4_gene_set,
        mutant4_name,
        plot_formats,
        title,
        path,
        ):
    """Generate a Venn diagram from 4 given sets and their respective names.
    """

    plt.style.use(['default'])

    dataset_dict = {
            mutant1_name : mutant1_gene_set,
            mutant2_name : mutant2_gene_set,
            mutant3_name : mutant3_gene_set,
            mutant4_name : mutant4_gene_set,
            }
    figure = venn(
            dataset_dict,
            cmap='plasma',
            hint_hidden=False,
            figsize=(8, 8),
            fontsize=10,
            legend_loc="best",
            ax=None,
            )


    plt.gcf()

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{path}.{format}", dpi=300)

        if format == "png":
            plt.savefig(f"{path}_transparent-bg.{format}", dpi=300, transparent=True)

    plt.close()


def generate_venn4_diagram_with_regulation_labels(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        mutant4_gene_set,
        mutant4_name,
        up_regulation_labels,
        down_regulation_labels,
        plot_formats,
        title,
        file_path,
        ):
    """Generates a venn diagram of 3 sets with custom labels.
    Labels are represented binary numbers that represent each intersection:
        - First bit equals to the first set.
        - Second bit equals to the second set.
        - Third bit equals to the third set.
    """

    plt.style.use(['default'])

    # Preparing the lables.
    labels = {
            '1000' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['1000'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['1000'])}"
                ),
            '1100' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['1100'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['1100'])}"
                ),
            '1110' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['1110'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['1110'])}"
                ),
            '1111' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['1111'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['1111'])}"
                ),
            '0111': (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['0111'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['0111'])}"
                ),
            '0011' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['0011'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['0011'])}"
                ),
            '0001' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['0001'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['0001'])}"
                ),
            '0100' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['0100'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['0100'])}"
                ),
            '0110' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['0110'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['0110'])}"
                ),
            '0010' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['0010'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['0010'])}"
                ),
            '1001' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['1001'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['1001'])}"
                ),
            '1101' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['1101'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['1101'])}"
                ),
            '1010' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['1010'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['1010'])}"
                ),
            '0101' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['0101'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['1010'])}"
                ),
            '1011' : (
                r"$\uparrow$"
                    +
                f"{len(up_regulation_labels['1011'])}\n"
                    +
                r"$\downarrow$"
                    +
                f"{len(down_regulation_labels['1011'])}"
                ),
            }

    dataset_dict = {
            mutant1_name : mutant1_gene_set,
            mutant2_name : mutant2_gene_set,
            mutant3_name : mutant3_gene_set,
            mutant4_name : mutant4_gene_set,
            }

    petal_labels = generate_petal_labels(dataset_dict.values())

    for k in petal_labels:
        petal_labels[k] = labels[k]

    draw_venn(
            petal_labels=petal_labels,
            dataset_labels=dataset_dict.keys(),
            colors=generate_colors(cmap='plasma', n_colors=4),
            hint_hidden=False,
            figsize=(8, 8),
            fontsize=10,
            legend_loc="best",
            ax=None,
            )

    plt.gcf()

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"{file_path}.{format}", dpi=300)

        if format == "png":
            plt.savefig(
                    f"{file_path}_transparent-bg.{format}",
                    dpi=300,
                    transparent=True
                    )

    # Clears current figure
    plt.close()

