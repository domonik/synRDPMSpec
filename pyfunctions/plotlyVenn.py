import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

import plotly as py
from plotly.offline import iplot
import plotly.graph_objs as go
import plotly.io as pio

import scipy


def venn_to_plotly(L_sets, L_labels=None, title=None, L_color = None):
    # get number of sets
    n_sets = len(L_sets)

    # choose and create matplotlib venn diagramm
    if n_sets == 2:
        if L_labels and len(L_labels) == n_sets:
            v = venn2(L_sets, L_labels)
        else:
            v = venn2(L_sets)
    elif n_sets == 3:
        if L_labels and len(L_labels) == n_sets:
            v = venn3(L_sets, L_labels)
        else:
            v = venn3(L_sets)
    # supress output of venn diagramm
    plt.close()

    # Create empty lists to hold shapes and annotations
    L_shapes = []
    L_annotation = []

    # Define color list for sets
    # check for other colors: https://css-tricks.com/snippets/css/named-colors-and-hex-equivalents/
    if L_color is None:
        L_color = ['FireBrick', 'DodgerBlue', 'DimGrey']

    # Create empty list to make hold of min and max values of set shapes
    L_x_max = []
    L_y_max = []
    L_x_min = []
    L_y_min = []
    final_shapes = []
    for i in range(0, n_sets):
        # create circle shape for current set

        shape = go.layout.Shape(
            type="circle",
            xref="x",
            yref="y",
            x0=v.centers[i][0] - v.radii[i],
            y0=v.centers[i][1] - v.radii[i],
            x1=v.centers[i][0] + v.radii[i],
            y1=v.centers[i][1] + v.radii[i],
            fillcolor=L_color[i],
            line_color=L_color[i],
            opacity=0.75
        )
        final_shapes.append(
            go.layout.Shape(
                type="circle",
                xref="x",
                yref="y",
                x0=v.centers[i][0] - v.radii[i],
                y0=v.centers[i][1] - v.radii[i],
                x1=v.centers[i][0] + v.radii[i],
                y1=v.centers[i][1] + v.radii[i],
                fillcolor="rgba(255,0,0,0)",
                line_color="black",
                opacity=1
            )
        )

        L_shapes.append(shape)

        # create set label for current set

        anno_set_label = go.layout.Annotation(
            xref="x",
            yref="y",
            x=v.set_labels[i].get_position()[0],
            y=v.set_labels[i].get_position()[1],
            text=v.set_labels[i].get_text(),
            showarrow=False
        )

        L_annotation.append(anno_set_label)

        # get min and max values of current set shape
        L_x_max.append(v.centers[i][0] + v.radii[i])
        L_x_min.append(v.centers[i][0] - v.radii[i])
        L_y_max.append(v.centers[i][1] + v.radii[i])
        L_y_min.append(v.centers[i][1] - v.radii[i])

    # determine number of subsets
    n_subsets = sum([scipy.special.binom(n_sets, i + 1) for i in range(0, n_sets)])
    for i in range(0, int(n_subsets)):
        # create subset label (number of common elements for current subset

        anno_subset_label = go.layout.Annotation(
            xref="x",
            yref="y",
            x=v.subset_labels[i].get_position()[0],
            y=v.subset_labels[i].get_position()[1],
            text=v.subset_labels[i].get_text(),
            showarrow=False
        )

        L_annotation.append(anno_subset_label)

    # define off_set for the figure range
    off_set = 0.2

    # get min and max for x and y dimension to set the figure range
    x_max = max(L_x_max) + off_set
    x_min = min(L_x_min) - off_set
    y_max = max(L_y_max) + off_set
    y_min = min(L_y_min) - off_set

    # create plotly figure
    p_fig = go.Figure()

    # set xaxes range and hide ticks and ticklabels
    p_fig.update_xaxes(
        range=[x_min, x_max],
        showticklabels=False,
        ticklen=0
    )

    # set yaxes range and hide ticks and ticklabels
    p_fig.update_yaxes(
        range=[y_min, y_max],
        scaleanchor="x",
        scaleratio=1,
        showticklabels=False,
        ticklen=0
    )

    # set figure properties and add shapes and annotations
    p_fig.update_layout(
        plot_bgcolor='white',
        margin=dict(b=0, l=10, pad=0, r=10, t=40),
        width=800,
        height=400,
        shapes=L_shapes + final_shapes,
        annotations=L_annotation,
        title=dict(text=title, x=0.5, xanchor='center')
    )

    return p_fig