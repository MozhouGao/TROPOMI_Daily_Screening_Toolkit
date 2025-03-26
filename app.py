import dash
from dash import dcc, html
import dash_leaflet as dl
import dash_leaflet.express as dlx
from dash.dependencies import Input, Output, State
import json
import pandas as pd
import io
#import dash_extensions as de
from datetime import date, timedelta, datetime
from TROPOMI_toolkit import download_TROPOMI_CH4_L2_data, Load_CH4, screening_plumes, generate_results

app = dash.Dash(__name__, external_stylesheets=["https://codepen.io/chriddyp/pen/bWLwgP.css"])

app.layout = html.Div([
    html.Div([
        html.H3(children = "TROPOMI Methane Plume Daily Screening Toolkit - V1.1",
        style = {'textAlign': 'center','color': '#415F4A'}),
        
        html.H5(children = "1. Select Date Range",
               style = {'textAlign': 'left','color': '#415F4A'}),
        dcc.DatePickerRange(
            id='date_picker',
            min_date_allowed=date(2019, 1, 1),
            max_date_allowed=date.today() + timedelta(days=3),
            initial_visible_month=date.today()
        ),
        html.Br(),
        html.Div(id='date-picker-range'),
        
        html.H5(children = "3. Download Level-2 TROPOMI Methane Observations",
               style = {'textAlign': 'left','color': '#415F4A'}),
        html.P("Click to download TROPOMI Level-2 nc files to the local path: ~/TROPOMI_data"),
        html.Div([
                dcc.Loading(id='loading',
                            children=[html.Div([html.Div(id="loading-output")])],
                            type = 'circle')
                ]),
        html.Button("Download", id="download_button", n_clicks=0),
        html.Div(id = 'download-log'),

        html.H5(children = "4. Screening",
               style = {'textAlign': 'left','color': '#415F4A'}),
        html.Div(
            [dcc.Input(id='thda', type='number', placeholder='Threshold delta (default = 15)',
                      style={'font-size': '11px', 
                             'width': '200px',   # Set width
                             'height': '40px',
                            'padding': '10px'}   # Set height}
            )],
        ),
        html.Div(
            [dcc.Input(id='min_pix', type='number', placeholder='Minimum pixel count (default = 1)',
                      style={'font-size': '11px',
                             'width': '200px',  
                             'height': '40px',
                             'padding': '10px'})]
        ),
        html.Div([html.Br()]),
        html.Button("Start Screening", id="screening-val", n_clicks=0),
    ], style={"width": "25%", "display": "inline-block", "verticalAlign": "top", "padding": "10px"}),
    
    html.Div([
        html.H5(children = "2. Select Region on Map",
               style = {'textAlign': 'left','color': '#415F4A'}),
        dl.Map([
            dl.TileLayer(id="basemap", url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"),
            dl.FeatureGroup([dl.EditControl(
                id="edit_control", 
                draw={"polygon": False, "polyline": False, "rectangle": True, "circle": False, "marker": False},
                edit={"remove": True}
            )], id="feature_group")
        ], id="map", style={"width": "95%", "height": "500px"}, center = [51.0447, -114.0719], zoom = 5),
        html.Div(id="polygon_data", style={"whiteSpace": "pre-wrap"}),
        html.Div(id = 'plot-log'),
    ], style={"width": "70%", "display": "inline-block", "padding": "10px"}),
    

])

#### select date range 
@app.callback(
    Output('date-picker-range', 'children'),
    Input('date_picker','start_date'),
    Input('date_picker','end_date')
)

def update_date_range(start_date,end_date):

    string_prefix = 'You have selected: '
    if start_date is not None:
        start_date_object = date.fromisoformat(start_date)
        start_date_string = start_date_object.strftime('%B %d, %Y')
        string_prefix = string_prefix + start_date_string 
    if end_date is not None:
        end_date_object = date.fromisoformat(end_date)
        end_date_string = end_date_object.strftime('%B %d, %Y')
        string_prefix = string_prefix + ' to ' + end_date_string
    if len(string_prefix) == len('You have selected: '):
        return 'Select a date range of plume screening'
    else:
        return string_prefix

#### define region 
@app.callback(
    Output("polygon_data", "children"),
    Input("edit_control", "geojson"),
)
def update_polygon_data(geojson):
    if geojson is None:
        return "Draw a rectangle to define the region."
    
    try:
        features = geojson.get("features", [])
        for feature in features:
            if feature["geometry"]["type"] == "Polygon":
                coordinates = feature["geometry"]["coordinates"][0]
                return f"Region Selected: {coordinates[:-1]}"
    except Exception as e:
        return f"Error processing geometry: {str(e)}"
    
    return "Draw a rectangle to define the region."

@app.callback(
    Output("download-log","children"),
    Input('date_picker','start_date'),
    Input('date_picker','end_date'),
    Input('download_button', 'n_clicks')
    )

def download_data(start_date,end_date,n_clicks):
    download_prefix = "" 
    if n_clicks > 0: 
        if start_date and end_date is not None: 
            start_date_object = date.fromisoformat(start_date)
            input_start_date = start_date_object.strftime('%Y%m%d')
            end_date_object = date.fromisoformat(end_date)
            input_end_date = end_date_object.strftime('%Y%m%d')
            download_prefix = download_TROPOMI_CH4_L2_data(start_date, end_date)
                
        else: 
            download_prefix = "Please select date range"
    return download_prefix

@app.callback(
    Output("plot-log","children"),
    Input('date_picker','start_date'),
    Input('date_picker','end_date'),
    Input("edit_control", "geojson"),
    Input("thda","value"),
    Input("min_pix","value"),
    Input('screening-val', 'n_clicks')
    )

def run_analysis(start_date, end_date, geojson, thda, min_pix,n_clicks):
    if n_clicks > 0:
        if start_date is None and end_date is None: 
            return r"Please select date range"
        else: 
            if geojson is None: 
                return r"Please define the region"
            else: 
                if thda is None or min_pix is None: 
                    return r"Please input the threshold"
                else:
                    start_date_object = date.fromisoformat(start_date)
                    start_date_string = start_date_object.strftime('%Y-%m-%d')
                    end_date_object = date.fromisoformat(end_date)
                    end_date_string = end_date_object.strftime('%Y-%m-%d')
                    features = geojson.get("features", [])
                    for feature in features:
                        if feature["geometry"]["type"] == "Polygon":
                            coordinates = feature["geometry"]["coordinates"][0]

                            upper_left_x, upper_left_y = coordinates[1]
                            lower_right_x,lower_right_y = coordinates[3]
                            
                            cdate = start_date_object
                            fch4_list = [] 
                            fwind_list = []
                            fpressure_list = []
                            grid_lons_list = [] 
                            grid_lats_list = [] 
                            Dates = [] 
                            while cdate.day < end_date_object.day: 
                            
                                grid_lon,grid_lat,fch4,fwind,fpressure = Load_CH4(lower_right_y, upper_left_y, upper_left_x, lower_right_x, 
                                                              cdate, qa_pass = 0.5)
                                grid_lons_list.append(grid_lon)
                                grid_lats_list.append(grid_lat)
                                fch4_list.append(fch4)
                                fwind_list.append(fwind)
                                fpressure_list.append(fpressure)
                                Dates.append(cdate.strftime('%Y%m%d'))
                                cdate += timedelta(days = 1)

                            thda = int(thda)
                            min_pix = int(min_pix)
                    
                            path_list = [] 
                    
                            for ele in zip (fch4_list,fwind_list,fpressure_list, grid_lons_list,grid_lats_list,Dates):
                    
                                detected_plumes, detected_plume_wind, detected_plume_pressure,detected_plumes_lons, detected_plumes_lats = screening_plumes(ele[0],ele[1],ele[2],ele[3],ele[4],thda,min_pix)
                                if len(detected_plumes) > 0: 
                                    figure_path = generate_results(ele[3],ele[4],ele[0],detected_plumes,detected_plume_wind, 
                                                           detected_plume_pressure, detected_plumes_lons,detected_plumes_lats,
                                                           ele[5])
                            
                                    path_list.append(figure_path)
                                    
                            if len(path_list)>0:
                                return rf"You can find the result:{path_list[-1]}"
                            else: 
                                return r"No ultra-emitter event detected in your region"

if __name__ == "__main__":
    app.run(debug=False)


