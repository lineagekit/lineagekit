import os
from scipy.stats import poisson
import random


def get_file_path(input_request: str):
    while True:
        file_path = input(input_request)
        if not os.path.exists(file_path):
            print("The specified file does not exist, try again")
        elif not os.path.isfile(file_path):
            print("The specified path is not a file, try again")
        else:
            return file_path


def get_non_empty_string(input_request: str) -> str:
    while True:
        response = input(input_request)
        if not response:
            print("Specify a non-empty string")
        else:
            return response


def get_float_value(input_request: str) -> float:
    while True:
        response = input(input_request)
        if not response:
            print("Specify a non-empty string")
        try:
            return float(response)
        except ValueError:
            print("Specify a float")


def get_non_existing_path(input_request: str):
    while True:
        path = input(input_request)
        if not path:
            print("Specify a non-empty string")
        elif os.path.exists(path):
            print("The specified path exists, try again")
        else:
            return path


def get_natural_number_input(input_request: str):
    while True:
        try:
            result = int(input(input_request))
            if result < 1:
                print("Specify a positive value")
            else:
                return result
        except ValueError:
            print("You need to specify an integer")


def get_natural_number_input_in_bounds(input_request: str, lower_bound: int, upper_bound: int):
    if lower_bound >= upper_bound:
        raise ValueError("Lower bound cannot be greater than upper bound")
    while True:
        value = get_natural_number_input(input_request)
        if value < lower_bound or value > upper_bound:
            print("Value out of bounds, try again")
            continue
        return value


def random_subselect_poisson(input_list, percentage):
    n = len(input_list)
    math_expectation = n * percentage
    number_of_errors = min(poisson.rvs(math_expectation), n)
    return random.sample(input_list, number_of_errors)
