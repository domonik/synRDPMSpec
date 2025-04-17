from RAPDOR.datastructures import RAPDORData
import plotly.graph_objects as go
import numpy as np


def main(file):
    rapdor_data = RAPDORData.from_file(file)
    fig = go.Figure()
    z = np.nansum(rapdor_data.array, axis=0)
    fig.add_trace(
        go.Heatmap(
            x=np.arange(1, 21),
            y=np.arange(6),
            z=z
        )
    )
    p = 0

if __name__ == '__main__':
    file = "Pipeline/RAPDORAnalysis/GradRData.json"
    main(file)