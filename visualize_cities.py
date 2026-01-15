# Author: Ryan Glenn
# Date: 02/23/23
# Purpose: Visualizes most populous 50 cities

from city import City
from sort_cities import city_list
from cs1lib import *


infile = open("cities_population.txt", "r")

for line in infile:  # reads through cities_population text file
    line = line.strip().split(",")  # gets rid of spaces on the ends of the string and separates elements by commas

    city_list.append(City(None, line[0], None, line[1], line[2], line[3]))

infile.close()  # closes file


i = 0


def draw():
    global i
    if i == 0:
        img = load_image("world.png")
        draw_image(img, 0, 0)  # draws image of world map onto graphics

    if i <= 50:
        city_list[i].draw(0, 0)  # draws cities' locations for first 50 cities in sorted population
        i = i + 1


start_graphics(draw, framerate=5, width=720, height=360)  # sets dimensions for graphics window and draws
clear()
