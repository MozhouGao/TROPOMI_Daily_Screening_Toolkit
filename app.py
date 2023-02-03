#import plotly.graph_objects as go
import os 
import time
from sentinelsat import SentinelAPI
from dash import Dash, Input, Output, State, html, dcc, dash_table, callback
from datetime import date,timedelta
from TROPOMI_toolkit import download_TROPOMI_CH4_L2_data,Load_CH4,screening_plumes,generate_results

app = Dash(__name__)

app.layout = html.Div([
    html.H1(
        children='TROPOMI Daily Screening Toolkit - V1.0',
        style={
            'textAlign': 'center',
            'color': '#415F4A',
            "font-family":"cursive"
        }
    ),
    
    html.Div([
    html.P('This toolkit is developed to automatically screen the suspect methane (CH4) plumes over a user defined region based on the public accessible satellite observations (i.e., methane dry air mixing ratio from TROPOMI). Users are required to input several parameters to complete the screening process, including: the screening date(s), region boundaries, screening criteria (i.e., threshold enhancement, number of valid plume pixels. The output of each screening run is a XCH4 concentration map by highlighting the regions with high probability of detecting suspect methane plume(s). A list of the potential source locations is also available from the screening result. The noteworthy point is that this toolkit is designed neither for pinpointing the emission source at the facility or component level, nor for screening of small methane leaks (generally <25 tons/hour). For source attributions of the detected suspect plumes, the follow-up targeted fine-scale observations over the regions with high probability of detecting suspect methane plume(s) are required.'),
            ]),
    
    html.H3(children="1. Select the date",
            style={
                'textAlign': 'left',
                'color': '#415F4A',
                "font-family":"cursive"
            } ),
    html.P("Please specify an individual date or a period for daily screening Multiple days screening may cause longer data processing time. Note: If a period is selected, only the screening result of the last day will be displayed on the webpage. For full list of the daily screening results, check the local path: ~/TROPOMI_Daily_Screening_Toolkit/assets."),
    
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
    html.H3(children="2. Define regions",
            style={
                'textAlign': 'left',
                'color': '#415F4A',
                "font-family":"cursive"
            }),
    html.P("Currently the toolkit only supports screening over a rectangular region. Click to confirm. Longitude range: -180 ~ 180, latitude range: -90 ~ 90."),
    
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
    html.P("Click to download the data files to the local path: ~/TROPOMI_Daily_Screening_Toolkit/TROPOMI_data."),
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
    html.P("Enter the Threshold delta (ΔXCH4,thr; defaut = 15) and Minimum pixel count (n; defaut = 1). Then click to kick off the daily plume screening. The screening time may vary with region size and number of days. Please do NOT hit on multiple times. Thanks for your patience."),
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
            end_date_object = date.fromisoformat(end_date)
            
            if west_long and east_long and south_lat and north_lat is not None:
                west_long = float(west_long)
                east_long = float(east_long)
                south_lat = float(south_lat)
                north_lat = float(north_lat)
                
                cdate = start_date_object
                fch4_list = [] 
                fwind_list = []
                fpressure_list = []
                grid_lons_list = [] 
                grid_lats_list = [] 
                Dates = [] 
                while cdate.day < end_date_object.day: 
                
                    grid_lon,grid_lat,fch4,fwind,fpressure = Load_CH4(south_lat, north_lat, west_long, east_long, 
                                                  cdate, qa_pass = 0.5)
                    grid_lons_list.append(grid_lon)
                    grid_lats_list.append(grid_lat)
                    fch4_list.append(fch4)
                    fwind_list.append(fwind)
                    fpressure_list.append(fpressure)
                    Dates.append(cdate.strftime('%Y%m%d'))
                    cdate += timedelta(days = 1)
                
                if thda and min_pix is not None:
                    thda = int(thda)
                    min_pix = int(min_pix)
                    
                    path_list = [] 
                    
                    for ele in zip (fch4_list,fwind_list,fpressure_list, grid_lons_list,grid_lats_list,Dates):
                    
                        detected_plumes, detected_plume_wind, detected_plume_pressure,detected_plumes_lons, detected_plumes_lats = screening_plumes(ele[0],
                                                                                                                                                   ele[1],
                                                                                                                                                   ele[2],
                                                                                                                                                   ele[3],
                                                                                                                                                   ele[4],
                                                                                                                                                   thda,
                                                                                                                                                   min_pix)
                        if len(detected_plumes) > 0: 
                            figure_path = generate_results(ele[3],ele[4],ele[0],detected_plumes,detected_plume_wind, 
                                                           detected_plume_pressure, detected_plumes_lons,detected_plumes_lats,
                                                           ele[5])
                            
                            path_list.append(figure_path)
                    
                    if len(path_list)>0:
                        return path_list[-1]
                    
                    else: 
                        return r"assets/pic.JPG"
                
            
if __name__ == "__main__":
     app.run_server(debug=False)