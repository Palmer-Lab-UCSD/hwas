""" Find sample set intersection

By: Robert Vogel
Affiliation: Palmer Lab at UCSD


"""
import importlib.resources
import subprocess
import configparser


def run(vcf: str) -> None:
    print(importlib.resources.files('hwas.R').joinpath('trait_intersect.R'))



