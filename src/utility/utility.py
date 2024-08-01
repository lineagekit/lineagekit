import os


def get_file_path(input_request: str):
    while True:
        file_path = input(input_request)
        if not os.path.exists(file_path):
            print("The specified file does not exist, try again")
        elif not os.path.isfile(file_path):
            print("The specified path is not a file, try again")
        else:
            return file_path
