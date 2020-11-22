""" Read a session xml file and parse the information. """
import argparse
import xml.etree.ElementTree as ET
from datetime import datetime
from pprint import pprint

from wana.xmldict import XmlDictConfig, XmlListConfig


class Session:

    def __init__(self, session_file):
        self.filename = session_file
        self.root = ET.parse(self.filename).getroot()
        self.data = XmlDictConfig(self.root)
        self.parse_tests()
        self.parse_meta()

    def parse_meta(self):
        self.start = datetime.fromisoformat(self.data["Start"])
        self.stop = datetime.fromisoformat(self.data["Stop"])
        self.duration = self.stop - self.start

    def parse_tests(self):
        if isinstance(self.data["TestList"]["Test"], XmlListConfig):
            self.tests = [x for x in self.data["TestList"]["Test"]]
        else:
            self.tests = [self.data["TestList"]["Test"]]
        self.tests.sort(key=lambda x: int(
            x["MoteList"]["Mote"][0]["Tag"]["Start"]))


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
