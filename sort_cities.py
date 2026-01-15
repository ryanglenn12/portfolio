# Author: Ryan Glenn
# Date: 02/27/23
# Purpose: sorting list of City objects in alphabetical order, by population, and by latitude

from quicksort import *
from city import City


def compare_population(city1, city2):
    if city2.pop <= city1.pop:
        return True  # allows for cities' populations to be sorted in decreasing order


def compare_names(a, b):  # allows for city names to be sorted in alphabetical order
    return a.name.lower() <= b.name.lower()


def compare_lat(city1, city2):
    if city1.lat <= city2.lat:
        return True  # allows for cities' latitudes to be sorted in increasing order


city_list = []
infile = open("world_cities.txt", "r")  # opens file that file pointer will read through
for line in infile:
    line = line.strip().split(",")  # gets rid of spaces on the ends of the string and separates elements by commas

    city_list.append(City(line[0], line[1], line[2], line[3], line[4], line[5]))

infile.close()  # closes file


def outfile(compare, text_file):  # define function that sets parameters for compare function and certain text file
    sort(city_list, compare)

    outfile = open(text_file, "w")  # opens new file that file pointer will write out

    for City in city_list:  # for each city in the list, write out the string of the city on a new line
        if City == city_list[len(city_list)-1]:
            outfile.write(str(City))  # for last line, a space is not included
        else:
            outfile.write(str(City) + "\n")

    outfile.close()  # closes file


outfile(compare_names, "cities_alpha.txt")
outfile(compare_lat, "cities_latitude.txt")
outfile(compare_population, "cities_population.txt")  # call for each outfile needed


