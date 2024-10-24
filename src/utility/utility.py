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


def random_subselect_poisson(input_list, percentage):
    n = len(input_list)
    math_expectation = n * percentage
    number_of_errors = min(poisson.rvs(math_expectation), n)
    return random.sample(input_list, number_of_errors)
