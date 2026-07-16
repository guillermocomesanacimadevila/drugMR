#!/usr/bin/env python3
import argparse
import polars as pl
import numpy as np 

# WE ARE NOT USING UKB PheWAS
# AS WALD RATIO H:0 TESTING -> WITH DELTA METHOD 
# SO Cov(B_JX, B_JY) != 0
# WE SHALL USE FinnGen then 


# Receive an input 
# Calc wald ratio
# Delta method 
# FDR correction here? - Maybe as a bool? 

def PheWAS(B_X, B_Y, SE_X, SE_Y):
    
    # wald ratio
    B_XY = B_Y / B_X
    return 