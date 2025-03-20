""" Find sample set intersection

By: Robert Vogel
Affiliation: Palmer Lab at UCSD


"""
import importlib.resources
import subprocess


import hwas


def run(vcf: str) -> None:
    print(importlib.resources.read_text(hwas.R, "R/trait_intersect.R"))



