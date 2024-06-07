import shlex
import subprocess
import importlib.resources as pkg_resources

from typing import Union, List


def get_executable_path(executable_name: str) -> str:
    # Use pkg_resources to find the binary within the package
    with pkg_resources.path('rhe', executable_name) as p:
        return str(p)


def run_genie(args: Union[str, List[str]]) -> str:

    if isinstance(args, str):
        args = shlex.split(args)  # Split the string into a list of arguments
    elif not isinstance(args, list):
        raise ValueError("The 'args' parameter should be either a string or a list of strings.")

    executable_path = get_executable_path('GENIE')
    result = subprocess.run([executable_path] + args, capture_output=True, text=True)
    return result.stdout
