import plotly.graph_objects as go
def create_scatter(pos_3d,
                   mode='markers',
                   marker_size=None,
                   name=None,
                   opacity=None,
                   color=None,
                   colorscale=None,
                   ):
    """Creates Scatter3D objects for use in plotly.
    Args:
        pos_3d: [N, 3] array containing N points with
            euclidean coordinates.
        mode: How to display points.
            Use 'markers' for scatter.
            Use 'lines' for lines connecting consecutive points.
            Use 'lines+markers' for scatter and lines.
        marker_size: Size of markers.
        name: Label of plotting layer to be displayed in legend.
        opacity: Transparency of points.
    """
    x, y, z = np.split(pos_3d, 3, axis=-1)
    args_dict = {
        'x': x[:, 0],
        'y': y[:, 0],
        'z': z[:, 0],
        'mode': mode,
        'marker': {}
    }
    if marker_size:
        args_dict['marker']['size'] = marker_size
    if name:
        args_dict['name'] = name
    if opacity:
        args_dict['opacity'] = opacity
    if color:
        args_dict['marker']['color'] = color
    if colorscale:
        args_dict['marker']['colorscale'] = colorscale
    return go.Scatter3d(**args_dict)