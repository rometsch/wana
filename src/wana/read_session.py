""" Read a session xml file and parse the information. """
import argparse
from pprint import pprint

from session import Session


def main():
    args = parse_cli_args()
    s = Session(args.session_file)
    pprint(s.data)


def parse_cli_args():
    """ Parse command line arguments. 

    Returns
    -------
    Namespace
        A namespace containing the arguments as attributes.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "session_file", help="XML file describing the session.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
