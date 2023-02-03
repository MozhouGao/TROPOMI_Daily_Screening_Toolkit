#import plotly.graph_objects as go
import os 
import time
from sentinelsat import SentinelAPI
from dash import Dash, Input, Output, State, html, dcc, dash_table, callback
from datetime import date,timedelta
from TROPOMI_toolkit import download_TROPOMI_CH4_L2_data,Load_CH4,screening_plumes,create_figures

app = Dash(__name__)

app.layout = html.Div([
    html.H1(
        children='TROPOMI Plume Screener',
        style={
            'textAlign': 'center',
            'color': '#415F4A',
            "font-family":"cursive"
        }
    ),
    html.H3(children="1. Select Date",
            style={
                'textAlign': 'left',
                'color': '#415F4A',
                "font-family":"cursive"
            } ),
    dcc.DatePickerRange(
    id='my-date-picker-range',
    min_date_allowed=date(2019, 1, 1),
    max_date_allowed=date.today() + timedelta(days=3),
    initial_visible_month=date(2022, 9, 1),
    ),
    html.Br(),
    html.Br(),
    html.Div(id='date-picker-range'),
    
    html.Div([html.Br()]),
    html.H3(children="2. Define Regions",
            style={
                'textAlign': 'left',
                'color': '#415F4A',
                "font-family":"cursive"
            }),
    
    html.Div(
        [html.I("Latitudes",
                style ={"margin-left": "470px"}),
         html.Br(),
         dcc.Input(id="north_lat",
                   type='number',
                   style ={"margin-left": "420px",
                           "margin-bottom":"20px"
                       }),
         html.Br(),
         html.I("Longitudes",
                style ={"margin-left": "200px"
                    }),
         dcc.Input(id='west_long',
                   type='number',
                   style = {
                       "margin-left": "10px"
                       }),
        dcc.Input(id='east_long',
                  type='number',
                  style = {
                      "margin-left": "100px"
                      }),
        html.Br(),
        dcc.Input(id = "south_lat",
                  type='number',
                  style ={"margin-left": "420px",
                          "margin-top":"20px"
                      })]),
    html.Br(),
    html.Button('Define polygon', id='poly-button', n_clicks=0),
    html.Br(),
    html.Br(),
    html.Div(id='output-longitudes'),
    html.Div([html.Br()]),
    html.H3(children="3. Download Level-2 TROPOMI methane observations",
            style={
                'textAlign': 'left',
                'color': '#415F4A',
                "font-family":"cursive"
            }),
    html.Div([dcc.Loading(
                   id="loading",
                   children=[html.Div([html.Div(id="loading-output")])],
                   type="circle")]),
    html.Button("Download",id="submit-val",n_clicks=0),
    html.Div([html.Br()]),
    html.Div(id='download-log'),
    html.Div([html.Br()]),
    html.H3(children="4. Start screening",
            style={
                'textAlign': 'left',
                'color': '#415F4A',
                "font-family":"cursive"
            }),

    html.Div(
        [html.I("Threshold delta"),
          html.Br(),
        dcc.Input(id='thda',type='number')]),

    html.Div(
        [html.I("Minimum pixel count"),
          html.Br(),
        dcc.Input(id='min_pix', type='number',)]),
    html.Div([html.Br()]),
    html.Button("Start screening",id="screening-val",n_clicks=0),
    html.Div([html.Br()]),
    html.Div([html.Img(id='plot1')]),
    ])


@app.callback(
    Output('date-picker-range', 'children'),
    Input('my-date-picker-range', 'start_date'),
    Input('my-date-picker-range', 'end_date'))    

def update_output1(start_date,end_date):

    string_prefix = 'You have selected: '
    if start_date is not None:
        start_date_object = date.fromisoformat(start_date)
        start_date_string = start_date_object.strftime('%B %d, %Y')
        string_prefix = string_prefix + "Start Date:" + start_date_string 
    if end_date is not None:
        end_date_object = date.fromisoformat(end_date)
        end_date_string = end_date_object.strftime('%B %d, %Y')
        string_prefix = string_prefix + '|| End Date: ' + end_date_string
    if len(string_prefix) == len('You have selected: '):
        return 'Select a date to see it displayed here'
    else:
        return string_prefix
    
@app.callback(
    Output("output-longitudes","children"),
    Input("west_long","value"),
    Input("east_long","value"),
    Input("south_lat","value"),
    Input("north_lat","value"),
    Input('poly-button', 'n_clicks'))

def update_output2(west_long,east_long,south_lat,north_lat,n_clicks):
    region_prefix = "You have entered: "
    if west_long and east_long is not None: 
        west_long = float(west_long)
        east_long = float(east_long)
        if -180 < west_long < 180 and -180 < east_long < 180:  
            region_prefix = region_prefix + "Longitude range: " + "{} - {} degrees".format(west_long,east_long)
        else: 
            region_prefix = "longitude is ranging from -180 to 180 degrees"
    if south_lat and north_lat is not None: 
        south_lat = float(south_lat)
        north_lat = float(north_lat)
        if -90 < south_lat < 90 and -90 < north_lat < 90:   
            region_prefix = region_prefix + "|| Latitude range: " + "{} - {} degrees".format(south_lat,north_lat)
        else: 
            region_prefix = "latitude is ranging from -90 to 90 degrees"
    if len(region_prefix) == len('You have entered: '):
        return 'Please enter longitudes and latitudes to define your polygon'
    else:
        if n_clicks == 0:
            return "please click button to finish the polygon"
        else:
            return region_prefix 


@app.callback(
    Output("download-log","children"),
    Input('my-date-picker-range', 'start_date'),
    Input('my-date-picker-range', 'end_date'),
    Input("west_long","value"),
    Input("east_long","value"),
    Input("south_lat","value"),
    Input("north_lat","value"),
    Input('submit-val', 'n_clicks')
    )

def download_data(start_date,end_date,west_long,east_long,
                  south_lat,north_lat,n_clicks):
    
    
    if n_clicks == 0:
        download_prefix = "Click to download" 
    elif n_clicks > 0: 
        if start_date and end_date is not None: 
            if west_long and east_long and south_lat and north_lat is not None:
                start_date_object = date.fromisoformat(start_date)
                input_start_date = start_date_object.strftime('%Y%m%d')
                end_date_object = date.fromisoformat(end_date)
                input_end_date = end_date_object.strftime('%Y%m%d')
                
                
                download_prefix = download_TROPOMI_CH4_L2_data(input_start_date, input_end_date, 
                                                               west_long, east_long, 
                                                               south_lat, north_lat)
                
            else: 
                download_prefix = "Please define region"

    else: 
        download_prefix = download_prefix
    return download_prefix

@app.callback(
    Output("plot1","src"),
    Input('my-date-picker-range', 'start_date'),
    Input('my-date-picker-range', 'end_date'),
    Input("west_long","value"),
    Input("east_long","value"),
    Input("south_lat","value"),
    Input("north_lat","value"),
    Input("thda","value"),
    Input("min_pix","value"),
    Input('screening-val', 'n_clicks')
    )

def screening(start_date,end_date, west_long,east_long,
                south_lat,north_lat,thda,min_pix,n_clicks):
    if n_clicks > 0:
        if west_long and east_long and south_lat and north_lat is not None:
            start_date_object = date.fromisoformat(start_date)
            input_start_date = start_date_object.strftime('%Y%m%d')
            end_date_object = date.fromisoformat(end_date)
            input_end_date = end_date_object.strftime('%Y%m%d')
            
            if west_long and east_long and south_lat and north_lat is not None:
                west_long = float(west_long)
                east_long = float(east_long)
                south_lat = float(south_lat)
                north_lat = float(north_lat)
                grid_lon,grid_lat,fch4 = Load_CH4(south_lat, north_lat, west_long, east_long, 
                                                  start_date_object,end_date_object, qa_pass = 0.5)
                
                if thda and min_pix is not None:
                    thda = int(thda)
                    min_pix = int(min_pix)
                    detected_plumes, detected_plumes_lons, detected_plumes_lats = screening_plumes(fch4,
                                                                                           grid_lon,
                                                                                           grid_lat,
                                                                                           thda,
                                                                                           min_pix)
                    if len(detected_plumes) > 0: 
                        figure_path = create_figures(grid_lon,grid_lat,fch4,detected_plumes,
                                                     detected_plumes_lons,detected_plumes_lats,
                                                     input_start_date,input_end_date)
                        return figure_path
                    else: 
                        return r"assets/pic.JPG"
                
            
if __name__ == "__main__":
     app.run_server(debug=False)