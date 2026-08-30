import os
from datetime import date, timedelta

import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_leaflet as dl

from TROPOMI_toolkit import (
    bbox_from_coordinates,
    download_TROPOMI_CH4_L2_data,
    generate_results,
    iter_screening_dates,
    Load_CH4,
    screening_plumes,
)

DEFAULT_THRESHOLD_DELTA = 15
DEFAULT_MIN_PIXELS = 1

app = dash.Dash(__name__, external_stylesheets=["https://codepen.io/chriddyp/pen/bWLwgP.css"])

app.layout = html.Div([
    html.Div([
        html.H3(
            children="TROPOMI Methane Plume Daily Screening Toolkit - V1.2",
            style={"textAlign": "center", "color": "#415F4A"},
        ),
        html.H5("1. Select Date Range", style={"textAlign": "left", "color": "#415F4A"}),
        dcc.DatePickerRange(
            id="date_picker",
            min_date_allowed=date(2019, 1, 1),
            max_date_allowed=date.today() + timedelta(days=3),
            initial_visible_month=date.today(),
        ),
        html.Br(),
        html.Div(id="date-picker-range"),
        html.H5(
            "3. Download Level-2 TROPOMI Methane Observations",
            style={"textAlign": "left", "color": "#415F4A"},
        ),
        html.P("Click to download TROPOMI Level-2 nc files to ./TROPOMI_data"),
        html.Button("Download", id="download_button", n_clicks=0),
        dcc.Loading(
            id="loading-download",
            type="circle",
            children=html.Div(id="download-log"),
        ),
        html.H5("4. Screening", style={"textAlign": "left", "color": "#415F4A"}),
        html.Div([
            dcc.Input(
                id="thda",
                type="number",
                placeholder="Threshold delta (default = 15)",
                style={"font-size": "11px", "width": "200px", "height": "40px", "padding": "10px"},
            )
        ]),
        html.Div([
            dcc.Input(
                id="min_pix",
                type="number",
                placeholder="Minimum pixel count (default = 1)",
                style={"font-size": "11px", "width": "200px", "height": "40px", "padding": "10px"},
            )
        ]),
        html.Div([html.Br()]),
        html.Button("Start Screening", id="screening-val", n_clicks=0),
    ], style={"width": "25%", "display": "inline-block", "verticalAlign": "top", "padding": "10px"}),
    html.Div([
        html.H5("2. Select Region on Map", style={"textAlign": "left", "color": "#415F4A"}),
        html.P(
            "Use the rectangle button on the map toolbar (below the zoom controls) to draw the region.",
            style={"marginBottom": "6px"},
        ),
        dl.Map([
            dl.TileLayer(id="basemap", url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"),
            dl.FeatureGroup([
                dl.EditControl(
                    id="edit_control",
                    position="topleft",
                    draw={
                        "polygon": False,
                        "polyline": False,
                        "rectangle": True,
                        "circle": False,
                        "circlemarker": False,
                        "marker": False,
                    },
                    edit={"remove": True},
                )
            ], id="feature_group"),
        ], id="map", style={"width": "95%", "height": "500px"}, center=[51.0447, -114.0719], zoom=5),
        html.Div(id="polygon_data", style={"whiteSpace": "pre-wrap"}),
        dcc.Loading(
            id="loading-screening",
            type="circle",
            children=html.Div(id="plot-log"),
        ),
    ], style={"width": "70%", "display": "inline-block", "padding": "10px"}),
])


def _first_polygon_ring(geojson):
    if not geojson:
        return None
    for feature in geojson.get("features") or []:
        geometry = feature.get("geometry") or {}
        if geometry.get("type") == "Polygon":
            coordinates = geometry.get("coordinates") or []
            if coordinates:
                return coordinates[0]
    return None


def _result_view(message, figure_path=None):
    children = [html.P(message)]
    if figure_path:
        children.append(html.P(f"Saved to {figure_path}"))
        children.append(html.Img(
            src=app.get_asset_url(os.path.basename(figure_path)),
            style={"width": "100%", "maxWidth": "720px"},
        ))
    return html.Div(children)


@app.callback(
    Output("date-picker-range", "children"),
    Input("date_picker", "start_date"),
    Input("date_picker", "end_date"),
)
def update_date_range(start_date, end_date):
    if start_date is None and end_date is None:
        return "Select a date range of plume screening"
    parts = ["You have selected:"]
    if start_date is not None:
        parts.append(date.fromisoformat(start_date).strftime("%B %d, %Y"))
    if end_date is not None:
        parts.append("to")
        parts.append(date.fromisoformat(end_date).strftime("%B %d, %Y"))
    elif start_date is not None:
        parts.append("(single day)")
    return " ".join(parts)


@app.callback(
    Output("polygon_data", "children"),
    Input("edit_control", "geojson"),
)
def update_polygon_data(geojson):
    ring = _first_polygon_ring(geojson)
    if ring is None:
        return "Draw a rectangle to define the region."
    try:
        minlat, maxlat, minlon, maxlon = bbox_from_coordinates(ring)
        return (
            f"Region selected: lon [{minlon:.2f}, {maxlon:.2f}], "
            f"lat [{minlat:.2f}, {maxlat:.2f}]"
        )
    except Exception as exc:
        return f"Error processing geometry: {exc}"


@app.callback(
    Output("download-log", "children"),
    Input("download_button", "n_clicks"),
    State("date_picker", "start_date"),
    State("date_picker", "end_date"),
    prevent_initial_call=True,
)
def download_data(n_clicks, start_date, end_date):
    if not n_clicks:
        raise PreventUpdate
    if start_date is None:
        return "Please select a start date."
    if end_date is None:
        end_date = start_date
    try:
        return download_TROPOMI_CH4_L2_data(start_date, end_date)
    except Exception as exc:
        return f"Download failed: {exc}"


@app.callback(
    Output("plot-log", "children"),
    Input("screening-val", "n_clicks"),
    State("date_picker", "start_date"),
    State("date_picker", "end_date"),
    State("edit_control", "geojson"),
    State("thda", "value"),
    State("min_pix", "value"),
    prevent_initial_call=True,
)
def run_analysis(n_clicks, start_date, end_date, geojson, thda, min_pix):
    if not n_clicks:
        raise PreventUpdate
    if start_date is None:
        return "Please select a start date."
    if end_date is None:
        end_date = start_date

    ring = _first_polygon_ring(geojson)
    if ring is None:
        return "Please define the region."

    threshold_delta = DEFAULT_THRESHOLD_DELTA if thda is None else int(thda)
    min_pixels = DEFAULT_MIN_PIXELS if min_pix is None else int(min_pix)

    try:
        minlat, maxlat, minlon, maxlon = bbox_from_coordinates(ring)
        last_path = None
        screened_days = 0
        plume_days = 0
        for screening_date in iter_screening_dates(start_date, end_date):
            grid_lon, grid_lat, fch4, fwind, fpressure = Load_CH4(
                minlat, maxlat, minlon, maxlon, screening_date, qa_pass=0.5
            )
            screened_days += 1
            detected = screening_plumes(
                fch4, fwind, fpressure, grid_lon, grid_lat, threshold_delta, min_pixels
            )
            last_path = generate_results(
                grid_lon,
                grid_lat,
                fch4,
                *detected,
                screening_date.strftime("%Y%m%d"),
            )
            if detected[0]:
                plume_days += 1

        if screened_days == 0:
            return "No dates to screen."
        if plume_days == 0:
            return _result_view(
                "No ultra-emitter event detected in your region.",
                last_path,
            )
        return _result_view(
            f"Suspect plumes found on {plume_days} of {screened_days} day(s). "
            "Maps and CSVs are in ./assets.",
            last_path,
        )
    except Exception as exc:
        return f"Screening failed: {exc}"


if __name__ == "__main__":
    app.run(debug=False)
