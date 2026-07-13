#!/usr/bin/env python3
from pathlib import Path
import yaml 

class Config:
    def __init__(self, file):
        file = Path(file)
        with open(file, "r") as f:
            self.__dict__.update(yaml.safe_load(f))

    # def validate(self):
        
# class Checks:
#     def __init__(self):
# validation