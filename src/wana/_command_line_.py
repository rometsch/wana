#!/usr/bin/env python3
# Caller for scripts in the simscripts package
import argparse
import sys

import wana.plot
import wana.session


def main():
    args = parse_command_line_args()
    name = args.script
    sys.argv = [sys.argv[0]] + args.args
    if name == "plot":
        wana.plot.main()
    elif name == "session":
        wana.plot.main()


def parse_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("script",
                        choices=["plot", "session"],
                        help="name of the script to run")
    parser.add_argument("args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
