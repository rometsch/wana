""" Read a session xml file and parse the information. """

from xml.etree import cElementTree as ElementTree
import xml.etree.ElementTree as ET
from pprint import pprint
from datetime import datetime


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


class XmlListConfig(list):
    def __init__(self, aList):
        for element in aList:
            if element:
                # treat like dict
                if len(element) == 1 or element[0].tag != element[1].tag:
                    self.append(XmlDictConfig(element))
                # treat like list
                elif element[0].tag == element[1].tag:
                    self.append(XmlListConfig(element))
            elif element.text:
                text = element.text.strip()
                if text:
                    self.append(text)


class XmlDictConfig(dict):
    '''
    Example usage:

    >>> tree = ElementTree.parse('your_file.xml')
    >>> root = tree.getroot()
    >>> xmldict = XmlDictConfig(root)

    Or, if you want to use an XML string:

    >>> root = ElementTree.XML(xml_string)
    >>> xmldict = XmlDictConfig(root)

    And then use xmldict for what it is... a dict.
    '''

    def __init__(self, parent_element):
        if parent_element.items():
            self.update(dict(parent_element.items()))
        for element in parent_element:
            if element:
                # treat like dict - we assume that if the first two tags
                # in a series are different, then they are all different.
                if len(element) == 1 or element[0].tag != element[1].tag:
                    aDict = XmlDictConfig(element)
                # treat like list - we assume that if the first two tags
                # in a series are the same, then the rest are the same.
                else:
                    # here, we put the list in dictionary; the key is the
                    # tag name the list elements all share in common, and
                    # the value is the list itself
                    aDict = {element[0].tag: XmlListConfig(element)}
                # if the tag has attributes, add those to the dict
                if element.items():
                    aDict.update(dict(element.items()))
                self.update({element.tag: aDict})
            # this assumes that if you've got an attribute in a tag,
            # you won't be having any text. This may or may not be a
            # good idea -- time will tell. It works for the way we are
            # currently doing XML configuration files...
            elif element.items():
                self.update({element.tag: dict(element.items())})
            # finally, if there are no child tags and no attributes, extract
            # the text
            else:
                self.update({element.tag: element.text})


if __name__ == "__main__":
    main()
