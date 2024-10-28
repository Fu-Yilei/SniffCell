#!/usr/bin/env python
import pysam, seaborn, re, scipy, os
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import pandas as pd
import warnings
from scipy.stats import zscore
warnings.filterwarnings("ignore", category=UserWarning, message=".*SmallSampleWarning.*")
warnings.filterwarnings("ignore", message="The index file is older than the data file")
from collections import Counter
import argparse

def parse_args(argv):
    parser = argparse.ArgumentParser(
                    prog='SniffMeth',
                    description='Sniffing CpG methylaiton changed around a (Mosaic SV)',
                    epilog='Version 0.1')
    parser = 
    