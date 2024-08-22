import requests
import requests_cache
from math import cos, tan, atan, degrees, exp
from pprint import pprint
from requests import Session
from requests_ratelimiter import LimiterAdapter
import time
import streamlit as st
import os

RAPID_API_KEY = os.getenv('RAPID_API_KEY')

def visualize_grid(grid):
  """
  Displays the grid in a Streamlit format with colors

  Args:
      grid: The Grid object containing the plant data
  """
  rows, cols = grid_size, grid_size
  grid_representation = [[None for _ in range(cols)] for _ in range(rows)]

  state_colors = {
      state_list[0]: "green",       #state_list[0]: alive and healthy
      state_list[1]: "yellow",      #state_list[1]: lighting up...
      state_list[2]: "red",         #state_list[2]: on fire
      state_list[3]: "gray"         #state_list[3]: burnt / smoldering
  }

  for row in range(rows):
    for col in range(cols):
      plant = grid.get_cell(row, col)
      color = state_colors[plant.state]
      grid_representation[row][col] = f'<span style="color:{color};">{plant.state}</span>'

  for row in grid_representation:
    st.write(" ".join(row), unsafe_allow_html=True)

st.title("Fire Spread Simulator")
st.markdown("*this simulation takes inital coordinates, an area and a Δt as input for a display of a forest fires approximate spreading time and pattern.  The grid represents a pine tree in each cell, all 1 meter apart from eachother.  It takes into consideration the the weather(temperature, wind speed/direction, air humidity, etc.), the difference in elevation between each tree  and the average humidity level of a pine tree.*")

start_coordinate_lat = st.number_input("Insert starting fire latitude in degrees, as precisely as possible", value=None, format="%0.10f")

start_coordinate_long = st.number_input("Insert starting fire longitude in degrees, as precisely as possible", value=None, format="%0.10f")

grid_size = st.number_input("input dimentions of visualized area (square), maximum is 33", min_value=1, max_value=33, value=None)

Δt = st.number_input("input time Δt for each update (time in simulation), under a minute is reccommended for more accuracy", min_value=1, value=None)

rows = grid_size
cols = grid_size

coordinate_lat = start_coordinate_lat  #max is 90/-90
coordinate_long = start_coordinate_long   #max is 180/-180

state_list = [  "(1)",
                "(2)",
                "(3)",
                "(4)"
                      ]

if start_coordinate_lat and start_coordinate_long and grid_size and Δt != None:
    Δt = int(Δt)

    def get_weather_data (coord_lat, coord_long):
        """
        Returns current weather data of specific coordinates
            
            Args:
                    coord_lat(int, float); latitude coordinates
                    coord_long(int, float); longitude coordinates
            returns:
                    weather_response.json(); all weather data of location in json format
        """
        weather = Session()
        adapter = LimiterAdapter(per_hour=999)
        weather.mount('https://', adapter)

        try:
            url = "https://weatherapi-com.p.rapidapi.com/current.json"
            querystring = {"q":str(coord_lat)+","+str(coord_long)}
            headers = {
                "x-rapidapi-key": RAPID_API_KEY,
                "x-rapidapi-host": "weatherapi-com.p.rapidapi.com"
            }
            weather_response = weather.get(url, headers=headers, params=querystring) 

            return weather_response.json()
        except requests.exceptions.HTTPError as err:
            print(f"HTTP error occurred: {err}")
    #request rate is WAY too high. Idk what to do. ALso it HAS to be high or the entire process is too slow
    def get_elevation_data(coord_lat, coord_long):
        """
        Returns elevation data of specific coordinates
            
            Args:
                    coord_lat(int, float); latitude coordinates
                    coord_long(int, float); longitude coordinates
            returns:
                    elevation_response.json(); all elevation data of location in json format
        """
        data = Session()
        adapter = LimiterAdapter(per_hour=30)
        data.mount('https://', adapter)

        try:
            
            url = "https://elevation-from-latitude-and-longitude.p.rapidapi.com/"
            querystring = {
                "lat": str(coord_lat),
                "lng": str(coord_long)
            }
            headers = {
                "x-rapidapi-key": RAPID_API_KEY,
                "x-rapidapi-host": "elevation-from-latitude-and-longitude.p.rapidapi.com"
            }
            elevation_response = data.get(url, headers=headers, params=querystring)
            elevation_response.raise_for_status()

            return elevation_response.json()
        except requests.exceptions.HTTPError as err:
            print(f"HTTP error occurred: {err}")

    def get_static_variables(strt_coordinate_lat, strt_coordinate_long, coordinat_lat, coordinat_long):
        """
        Returns variables that don't change over time (elevation)

            Args:
                    strt_coordinate_lat(int, float); latitude of 'starting' cell
                    strt_coordinate_long(int, float); longitude of 'starting' cell
                    coordinat_lat(int, float); latitude of other/next interested cell
                    coordinat_long(int, float); longitude of other/next interested cell
            Returns:
                    g(int); value that determines whether difference in elevation is positive or negative
                    θ(float); slope angle
            Raises:
                    ValueError; if there is insufficient information for given coordinates or if subscription limit is reached
        """
        all_start_elevation_data = get_elevation_data(strt_coordinate_lat, strt_coordinate_long)
        all_elevation_data = get_elevation_data(coordinat_lat, coordinat_long)

        i_old = all_start_elevation_data["elevation"]
        i_new = all_elevation_data["elevation"]
        r = 1                                           #distance between centers of cells
        i = i_new - i_old
        θ = degrees(atan(i/r))

        if i > 0:
            g = 1       
        else:
            g = -1
        return g,θ

    def get_dynamic_variables(strt_coordinate_lat, strt_coordinate_long, k = None):
        """
        Returns variables that vary as time passes

            Args:
                    strt_coordinate_lat(int, float); latitude of interested cell
                    strt_coordinate_long(int, float); longitude of interested cell
            Returns:
                    v(float); wind speed in kph
                    h(int); air humidity %
                    T(float); temperature in °C
                    β(int); angle betw. fire spread direction and wind direction (North is 0°, East is 90°)
            Raises:
                    ValueError; if there is insufficient information for given coordinates
        """

        all_weather_data = get_weather_data(strt_coordinate_lat, strt_coordinate_long)
        try:
            v = all_weather_data['current']['wind_kph']
            h = all_weather_data['current']['humidity']
            temp = all_weather_data['current']['temp_c']
            η = all_weather_data['current']['wind_degree']
            if k != None:
                β = η - k * 45
            else:
                β = 0
            return v,h,temp,β
        except:
            raise TypeError("insufficient weather data for these coordinates")
    v, h, temp, β = get_dynamic_variables(start_coordinate_lat, start_coordinate_long)

    def get_R(v,h,temp,β,g,θ):
        """
        Returns the forest fire spreading speed (R) in m/min for the specified cell data
            
            Args:
                    v(float); wind speed in kph
                    h(int); air humidity %
                    T(float); temperature in °C
                    β(int); angle betw. fire spread and wind direction (North is 0°, East is 90°)
                    g(int); value that specifies whether difference in elevation is positive or negative
                    θ(float); slope angle
            Returns:
                    R(float); forest fire spreading speed in m/min
        """
        m = 110                                 #average wood moisture % for pine trees
        w = int((v/0.836)**(2/3))
        Ro = 0.03*temp + 0.05*w + 0.01*(100-h) - 0.3
        r1 = 1.0372*exp(-0.057*m)*Ro
        Kβ = exp(0.1783*v*cos(β))
        Kθ = exp(3.553*g*tan(1.2*θ))
        #R_grass = Ro*Kθ*Kβ*1.6
        r = r1*Kθ*Kβ                       #missing components, this is a big approximation

        return r, r1

    class Plant:
        """
        A single plant, has many of it's own values, they're updated later on.

            __init__():
                        temp; temperature
                        w_speed; wind_speed
                        humid; water air %
                        lat; latitude
                        long; longitude
        Attributes:
            self.lat; latitude of plant
            self.long; longitude of plant
            self.state; state of plant
            self.temperature; temperature of plant
            self.wind_speed; wind speed around plant
            self.humidity; air-humidity around plant
            self.k; k value of plant (defined later on)
            self.β; β value of plant
            self.g; g value of plant
            self.θ; θ value of plant
            self.R; R value of plant

            self.t_until_on_fire; total time passed necessary for going aflame (seconds)
            self.t_capture_burning; serves as a record for the time t at which the plant started burning
            self.t_until_burnt; (12 hour) average pine tree burning time (seconds)
            self.t_capture_fire; serves as a record for the time t at which the plant caught fire
        """
        def __init__(self,temp,w_speed,humid,lat,long,initial_state=state_list[0]):
            self.lat = lat
            self.long = long
            self.state = initial_state
            self.temperature = temp
            self.wind_speed = w_speed
            self.humidity = humid
            self.k = None
            result = get_dynamic_variables(self.lat, self.long)
            self.β = result[-1]
            self.g, self.θ = get_static_variables(start_coordinate_lat, start_coordinate_long, self.lat, self.long)
            r, r1 = get_R(self.wind_speed, self.humidity, self.temperature, self.β, self.g, self.θ)
            self.R = r
            self.R1 = r1                            #initial spreading speed

            self.t_until_caught_fire = 60/self.R
            self.t_capture_catching_fire = False
            self.t_until_on_fire = 20*60            #general time it takes for a pine tree to start burning heavily  
            self.t_capture_burning = False
            self.t_until_burnt = 12 * 60 * 60       #general time it takes for a pine tree to finish burning and enter a smoldering phase
            self.t_capture_fire = False

        def __repr__(self):
            return f"Plant({self.state})"
    class Grid:
        """
        A square 2D grid (list of lists) where each cell is a Class: Plant type object

            __init__:

                    grid_size; dimentions of grid

            get_cell: gets a specific cell

                    row; row of interested cell
                    col; column of interested cell

            find_row_index: finds row index of known cell

            find_col_index: finds column index of known cell

        Attr.:
                self.grid; the grid where each cell is a plant object
        """
        def __init__(self, grid_size):
            self.grid = [[Plant(temp,v,h,coordinate_lat,coordinate_long) for _ in range(grid_size)] for _ in range(grid_size)]
            
        def get_cell(self, row, col):
            return self.grid[row][col]

        def find_row_index(self, plant):
            for row_index, row in enumerate(self.grid):
                if plant in row:
                    return row_index
            return None
                
        def find_col_index(self, plant):
            for row in self.grid:
                if plant in row:
                    return row.index(plant)
            return None
        
        def __getitem__(self, idx):
            idx_row = idx
            idx_col = idx
            return self.grid[idx_row][idx_col]
        
        def __repr__(self):
            return f"Grid({self.grid})"

    def correct_data_for_each_cell(starting_coordinate_lat, starting_coordinate_long):
        """
        Returns grid with corrected data (corrects coordinates as if centers of cells were 1 meter apart from eachother --> corrected weather and elevation data)

            Args:
                    rows(int); number of rows in grid, defaults to previous input
                    cols(int); number of columns in grid, defaults to previous input
            Returns:
                    grid(Class: Grid); grid with correct data for each cell (each 1 meter apart)
        """
        grid = Grid(grid_size)
        for row in range(rows):
            for col in range(cols):
                cell = grid.get_cell(row, col)

                lat = grid.get_cell(rows-1, cols-1).lat + (rows - row) * 1/111111
                long = grid.get_cell(rows-1, cols-1).long + (cols - col) * 1/111320
                cell.lat = lat
                cell.long = long
                k = cell.k

                v, h, temp, β = get_dynamic_variables(cell.lat, cell.long, k)
                g, θ = get_static_variables(starting_coordinate_lat,starting_coordinate_long,cell.lat,cell.long)
                starting_coordinate_lat = cell.lat
                starting_coordinate_long = cell.long
                r = get_R(v,h,temp,β,g,θ)
        
        yield grid
    for i in correct_data_for_each_cell(start_coordinate_lat, start_coordinate_long):
        grid = i

    def get_neighbors(grid, row, col):
        """
        Returns neighbors for a given cell in a grid, accounts for missing neighbors (grid border)

            Args:
                    grid(Class: Grid); the grid that is being referenced, possibly the corrected grid
                    row(int); the row index of the interested cell within the grid
                    col(int); the column index of the interested cell within the grid
            Returns:
                    A list of neighboring cells
        """
        rows = grid_size
        cols = grid_size
        offsets = [(-1, -1), (-1, 0), (-1, 1),
                (0, -1),           (0, 1),
                (1, -1),  (1, 0),  (1, 1)]
        neighbor_list = []

        for dr, dc in offsets:
            new_row, new_col = row + dr, col + dc
            if 0 <= new_row < rows and 0 <= new_col < cols:
                neighbor_list.append(grid.grid[new_row][new_col])
        return neighbor_list

    def change_neighbors_state(grid, row, col, t, new_state=state_list[1]):
        """
        Changes the state of neighbors to new_state if neighbors state is in state 0 and center cell is in center_state 

            Args:
                    grid(Class: Grid); the grid
                    row(int): The row index of the center cell
                    col(int): The column index of the center cell
                    center_state: The state index of the center cell
                    neighbor_state: The state index of the neighbor to be changed
                    new_state: The new state index for the neighbor
            Returns:
                    updated grid
        """
        for plant in get_neighbors(grid, row, col):
            neighbor_row = grid.find_row_index(plant)
            neighbor_col = grid.find_col_index(plant)
            neighbor = grid.grid[neighbor_row][neighbor_col]
            if neighbor.state == state_list[0]:
                if neighbor.t_capture_catching_fire == False:
                    neighbor.t_until_caught_fire = t + plant.t_until_caught_fire
                    neighbor.t_capture_catching_fire = True
                if neighbor.t_until_caught_fire <= t:
                    neighbor.state = new_state

    def assign_k_to_neighbors(grid, row, col):
        #goes: top-left, top, top-right, right, ...
        """
        Returns value of k (0-7) depending on direction of fire travel. if heading up k = 0, if right k = 2, if bottom k = 4,...

            Args:
                    grid(Class: Grid); interested grid
                    row(int); row of center cell
            Returns:
                    k(int); value that's part of the β equation
        """
        x = 0
        for plant in get_neighbors(grid, row, col):
            neighbor_row = grid.find_row_index(plant)
            neighbor_col = grid.find_col_index(plant)
            cell = grid.grid[neighbor_row][neighbor_col]
            cell.k = 7
            if x > 0:
                cell.k = -1
            cell.k += x
            x += 1

    def are_all_plants_burnt(grid, target_value=state_list[3]):
        """
        Checks if all elements in a grid are not equal to a certain state. It's for ending the loop if all cells are smoldering.
        
        Args:
            grid; the grid
            target_value; The target state to check for.

        Returns:
            True if all elements are equal to the target state, False otherwise.
        """
        for row in range(rows):
            for col in range(cols):
                if grid.get_cell(row, col).state != target_value:
                    return False
        return True

    def expansion(corrected_grid, Δt):
        """
        "flame" expansion function, updates a time variable and simulates the flame spreading through a grid (based on real coordinates data)

        Args:
            corrected_grid; the grid from the "correct_data_for_each_cell" function.
            Δt; choose how fast time passes in the simulation (how many calculated seconds pass per-loop)
        returns:
            a generator that updates as time passes(?)
        """
        t = 0 
        #first_flaming_cell
        center = corrected_grid.grid[rows//2][cols//2]
        center.state = state_list[1]
        center.t_until_caught_fire = 60 / center.R1

        flag = are_all_plants_burnt(corrected_grid)
        while flag == False:

            for _ in range(Δt):

                for row in range(rows):
                    for col in range(cols):
                        plant = corrected_grid.get_cell(row, col)
                        print(f"Row: {row}, Col: {col}, State: {plant.state}, time: {t}, t_until_on_fire: {plant.t_until_on_fire}, t_until_burnt: {plant.t_until_burnt}")

                        if plant.state == state_list[0]:
                            continue

                        elif plant.state == state_list[1]:
                            if plant.t_capture_burning == False:
                                plant.t_until_on_fire = t + plant.t_until_on_fire
                                plant.t_capture_burning = True

                            assign_k_to_neighbors(corrected_grid, corrected_grid.find_row_index(plant), corrected_grid.find_col_index(plant))
                            change_neighbors_state(corrected_grid, corrected_grid.find_row_index(plant), corrected_grid.find_col_index(plant), t)

                            if plant.t_until_on_fire <= t:
                                plant.state = state_list[2]
                            continue

                        elif plant.state == state_list[2]:
                            if plant.t_capture_fire == False:
                                plant.t_until_burnt = t + plant.t_until_burnt
                                plant.t_capture_fire = True

                            if plant.t_until_burnt <= t:
                                plant.state = state_list[3]

                            continue

                        else:
                            continue
            
                t += 1
            flag = are_all_plants_burnt(corrected_grid)
            yield corrected_grid, t

    def main():
    
        """Main function for the Streamlit app"""
        
        st.markdown("**(1):** alive and healthy")
        st.markdown("**(2):** lighting up...")
        st.markdown("**(3):** on fire")
        st.markdown("**(4):** burnt/smoldering") 

        for updated_grid, new_t in expansion(grid, Δt):
            st.subheader(f"**Time:** {round(new_t)} seconds")
            st.markdown(f"➡{round(new_t/60)} minutes")
            st.markdown(f"➡{round(new_t/3600, 1)} hours")
            st.markdown(f"➡{round(new_t/86400, 1)} days")
            visualize_grid(updated_grid)

    if __name__ == "__main__":
        main()