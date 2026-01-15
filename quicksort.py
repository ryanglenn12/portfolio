# Author: Ryan Glenn
# Date: 02/27/23
# Purpose: sorting list of City objects via quicksort

def compare_func(a, b):  # function that compares each sublist item with pivot
    return a <= b


def partition(the_list, p, r, compare_func):
    i = p - 1
    j = p
    pivot = the_list[r]  # sets pivot as the value of index r

    while j < r:
        if compare_func(the_list[j], pivot):  # if pivot is greater than the value at index j
            i = i + 1
            the_list[i], the_list[j] = the_list[j], the_list[i]  # switches values of i and j

        j = j + 1  # increment j regardless

    the_list[r], the_list[i + 1] = the_list[i + 1], the_list[r]  # switches values of i + 1 and r

    return i + 1  # returns index at which pivot is placed


a = [3, 9, 8, 4, 1, 7, 5]  # test for partition
print(partition(a, 0, 6, compare_func))


def quicksort(the_list, p, r, compare_func):
    if r - p >= 1:
        q = partition(the_list, p, r, compare_func)  # returns index at which pivot occurs with q

        quicksort(the_list, p, q - 1, compare_func)  # deals with first half
        quicksort(the_list, q + 1, r, compare_func)  # deals with second half


def sort(the_list, compare_func):  # calls quicksort with necessary optional parameters
    quicksort(the_list, p=0, r=len(the_list) - 1, compare_func=compare_func)

