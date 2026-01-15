# Author: Ryan Glenn
# Date: 02/23/23
# Purpose: City class file

from cs1lib import *
from random import randint


class City:
    def __init__(self, country_code, name, region, pop, lat, long):  # creates instance variables
        self.country_code = country_code
        self.name = name
        self.region = region
        self.pop = int(pop)
        self.lat = float(lat)
        self.long = float(long)

    def __str__(self):  # returns string that has city's name, population, latitude, and longitude
        return self.name + "," + str(self.pop) + "," + str(self.lat) + "," + str(self.long)

    def draw(self, cx, cy):
        pixels_per_long = 2  # scale constant for longitude
        pixels_per_lat = 2  # scale constant for latitude

        # scale self.lat and self.long to the size of the image
        px = self.long * pixels_per_long
        py = self.lat * pixels_per_lat

        # adjusting center of graphics of map so that center is where equator and prime meridian meet
        center_x = cx + 180 * pixels_per_long
        center_y = cy + 90 * pixels_per_lat

        disable_stroke()
        set_fill_color(randint(0, 1), randint(0, 1), randint(0, 1))  # random color
        draw_circle(center_x + px, center_y - py, 5)  # draw city in its real-life location on the map
