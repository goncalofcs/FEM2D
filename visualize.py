import matplotlib.pyplot as plt

class Visualize:
    def __init__(self, model, scale=150):
        self.model = model
        self.scale = scale

    def plot_configuration(self, subplot, displacements, title, show_legend=False):
        for index, element in enumerate(self.model.elements_list):
            node_i = element.node_i
            node_j = element.node_j

            # original coordinates
            x_orig = [
                self.model.nodes[node_i][0],
                self.model.nodes[node_j][0]
            ]
            y_orig = [
                self.model.nodes[node_i][1],
                self.model.nodes[node_j][1]
            ]
            subplot.plot(
                x_orig, y_orig,
                color="black",
                linewidth=1.5,
                marker="o",
                markersize=4,
                label="Original" if index == 0 else ""
            )

            # deformed coordinates
            x_def = [
                x_orig[0] + displacements[2 * node_i] * self.scale,
                x_orig[1] + displacements[2 * node_j] * self.scale
            ]
            y_def = [
                y_orig[0] + displacements[2 * node_i + 1] * self.scale,
                y_orig[1] + displacements[2 * node_j + 1] * self.scale
            ]
            subplot.plot(
                x_def, y_def,
                color="red",
                linestyle="--",
                linewidth=1.5,
                marker="s",
                markersize=3.5,
                label=f"Deformada ({self.scale}x)" if index == 0 else ""
            )

            # deformed elements annotation
            subplot.annotate(
                str(element.element_id + 1),
                ((x_def[0] + x_def[1]) / 2, (y_def[0] + y_def[1]) / 2),
                xytext=(0, -10),
                textcoords="offset points",
                color="red",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.8)
            )

        # deformed nodes annotation
        for node_id, (x, y) in self.model.nodes.items():
            x_def = x + displacements[2 * node_id] * self.scale
            y_def = y + displacements[2 * node_id + 1] * self.scale
            subplot.annotate(
                str(node_id + 1),
                (x_def, y_def),
                xytext=(6, -11),
                textcoords="offset points",
                color="red",
                fontsize=8,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85)
            )

        subplot.set_title(title)
        if show_legend:
            subplot.legend()
        subplot.axis("off")