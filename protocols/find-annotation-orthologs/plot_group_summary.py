#!/usr/bin/env python2

import json
import math
import sys

import matplotlib
# Force matplotlib not to use X11 backend, which produces exception when run
# over SSH.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_groups(groups, filename, make_y_scale_log, group_a_label, group_b_label):
  vals = []

  for group in groups:
    a_len = len(group['a'])
    b_len = len(group['b'])
    log = math.log(float(a_len) / float(b_len), 2)
    vals.append(log)

  bins = list(range(int(math.floor(min(vals))), int(math.ceil(max(vals))) + 1))
  plt.figure()
  plt.xticks(bins)
  plt.hist(vals, bins=bins, log=make_y_scale_log, facecolor='green', alpha=0.5)

  plt.title('Distribution of orthologous group gene counts')
  plt.xlabel('$log_2(%s / %s)$' % (group_a_label, group_b_label))
  plt.ylabel('Occurrences')
  plt.savefig(filename)

def main():
  group_a_label = sys.argv[1]
  group_b_label = sys.argv[2]
  output_linear = sys.argv[3]
  output_log    = sys.argv[4]

  ortho_groups = json.load(sys.stdin)['groups']
  plot_groups(ortho_groups, output_linear, False, group_a_label, group_b_label)
  plot_groups(ortho_groups, output_log,    True,  group_a_label, group_b_label)

if __name__ == '__main__':
  main()
