
%load_ext autoreload
%autoreload 2

import copy
import pyranges as pr
import pandas as pd
import numpy as np
import scipy.stats
import pyreadr
import pyBigWig
import swifter
import seaborn as sns
import torch
from sklearn.metrics import accuracy_score, roc_curve, roc_auc_score, precision_recall_curve, average_precision_score, f1_score
import matplotlib.pyplot as plt
import matplotlib
from statannotations.Annotator import Annotator
from glob import glob
from Bio import SeqIO


from helpers.plots import MotifMetrics, LoadedMotifMetrics, MetricsHandler

from helpers.motifs import motifs

# importing the sys module
import sys        
import os
 
project_dir = "./ML4RG-2023-project/"
sys.path.insert(0, project_dir)

results_dir = "./results/"

#sns.set_theme(context="poster")
sns.set_theme()
#sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
sns.set_context("talk", font_scale=1.7)#, rc={"font.size": 7})
#
#sns.set_context("poster", 
sns.set_style(style={'xtick.bottom': True,'ytick.left': True, 'axes.edgecolor': 'black'})
#sns.set_style("ticks")

model_paths = os.listdir(results_dir)
# remove the Random folder if it exists
if "random" in model_paths:
    model_paths.remove("random")
    model_paths.remove("entropy")
    model_paths.remove("dinuc_test")

print(model_paths)

# work around as plots.py doesn't use os.path.join to join paths...
model_names = list(map(lambda x: os.path.basename(x), model_paths))
model_map = {"full_species_aware": "Species aware all motifs"
             }

model_names = list(map(lambda x: model_map[x], model_names))
model_paths = list(map(lambda x: os.path.join(results_dir, x) + "/", model_paths))

print(model_paths)
print(model_names)
test_path = glob(os.path.join(data_dir ,"*.fa"))[0]

scer_mh = MetricsHandler(model_paths, model_names, test_path, 
                         motifs=motifs, seq_col="3UTR",
                         existing_probas=None,
                         n_random_kmers=6290,
                         )