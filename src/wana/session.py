import xml.etree.ElementTree as ET
from datetime import datetime

from xmldict import XmlDictConfig


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
        self.tests = [x for x in self.data["TestList"]["Test"]]
        self.tests.sort(key=lambda x: int(
            x["MoteList"]["Mote"][0]["Tag"]["Start"]))
