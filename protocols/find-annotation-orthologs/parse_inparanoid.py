#!/usr/bin/env python2
'''
This script parses the human-readable output file from InParanoid.
Unfortunately, only the human-readable output (as well as the HTML file
containing the same information) list all InParanoid results; the table and CSV
file both lack data.

Usage: parse_inparanoid.py Output.PRJEB506.munged.fa-PRJNA205202.munged.fa
'''

import math
import matplotlib
# Force matplotlib not to use X11 backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
from collections import defaultdict

# From http://stackoverflow.com/a/14620633/1691611
class AttrDict(dict):
  '''
  Permit access to dictionary keys as attributes (i.e., d['k'] can be accessed as d.k).
  '''
  def __init__(self, *args, **kwargs):
    super(AttrDict, self).__init__(*args, **kwargs)
    self.__dict__ = self

class NoTransitionError(Exception):
  '''
  Raised when no transition from state machine's current state for current input.
  '''
  pass

class StateMachine(object):
  def __init__(self, input_file):
    '''
    Initialize. Argument `input_file` corresponds to open file object.
    '''
    self._input_file = input_file
    self._state = self._file_header

  def run(self):
    '''
    Run state machine against results file, populating self._results.
    '''
    for line in self._input_file:
      # Strip trailing newline.
      self._current_line = line[:-1]
      self._state(self._current_line)

    self._results = {
      'groups': self._groups,
    }

  def results(self):
    '''
    Retrieve results after state machine run.
    '''
    return self._results

  def _switch(self, new_state, fire_for_current_line=False):
    '''
    Switch state machine's state. If fire_for_current_line is True, immediately
    run the new state with the current line. Otherwise, the state will not be
    run until the next input line.
    '''
    self._state = new_state
    if fire_for_current_line:
      self._state(self._current_line)


  ############
  # Miscellany
  ############
  def _is_bootstrap_line(self, line):
    return line.startswith('Bootstrap support for')

  def _is_separator(self, line):
    return line == 83*'_'
  
  def _parse_percentage(self, token):
    '''
    Convert '93.5%' to '0.935'.
    '''
    if not token.endswith('%'):
      raise Exception('Putative percentage "%s" does not end with "%"' % token)
    return float(token[:-1]) / 100


  ########
  # States
  ########
  def _file_header(self, line):
    if self._is_separator(line):
      self._groups = []
      self._switch(self._new_group, True)
    else:
      # Do nothing on file header lines.
      pass

  def _new_group(self, line):
    # Expect current line to be separator when in this state.
    if not self._is_separator(line):
      raise NoTransitionError()

    self._switch(self._group_header)
    self._current_group = AttrDict({'a': [], 'b': []})

  def _group_header(self, line):
    if line.startswith('Group of') or line.startswith('Score difference'):
      return
    self._switch(self._group_seq, True)

  def _group_seq(self, line):
    # Reached end of seqs within group. Now at list of bootstrap support values.
    if self._is_bootstrap_line(line):
      self._groups.append(self._current_group)
      self._switch(self._bootstrap_support, True)
      return

    fields = line.split('\t')
    fields = [f for f in fields if len(f) > 0] # Remove any empty fields.
    if len(fields) != 4:
      raise Exception('Unknown number of fields for group sequence: %s' % len(fields))

    seq_a = fields[0].strip()
    seq_b = fields[2].strip()

    if seq_a != '':
      score_a = self._parse_percentage(fields[1])
      self._current_group['a'].append((seq_a, score_a))
    if seq_b != '':
      score_b = self._parse_percentage(fields[3])
      self._current_group['b'].append((seq_b, score_b))

  def _bootstrap_support(self, line):
    if self._is_bootstrap_line(line):
      return
    self._switch(self._new_group, True)

def summarize_groups(groups):
  '''
  Suppose InParanoid was run on sequence sets A and B. This function returns a
  dictionary counting occurrences of the following orthologue groupp
  relationships:

    * 1_to_1: Instances where a single sequence in A corresponds to a single
      sequence in B
    * 1_to_n: Instances where a single sequence in A corresponds to multiple
      sequences in B
    * n_to_1: Instances where multiple sequences in A correspond to a single
      sequence in B
    * n_to_n: Instances where multiple sequences in A correspond to multiple
      sequences in B
  '''
  summary = defaultdict(int)

  for group in groups:
    a_len = len(group.a)
    b_len = len(group.b)

    if a_len == b_len == 1:
      summary['1_to_1'] += 1
    elif a_len > 1 and b_len > 1:
      summary['n_to_n'] += 1
    elif a_len > 1:
      summary['n_to_1'] += 1
    elif b_len > 1:
      summary['1_to_n'] += 1
    else:
      raise Exception('Unexpectedly empty group')


  return dict(summary)

def plot_groups(groups, filename, log_y_scale):
  vals = []

  for group in groups:
    a_len = len(group.a)
    b_len = len(group.b)
    log = math.log(float(a_len) / float(b_len), 2)
    vals.append(log)

  bins = list(range(int(math.floor(min(vals))), int(math.ceil(max(vals))) + 1))
  plt.figure()
  plt.xticks(bins)
  plt.hist(vals, bins=bins, log=log_y_scale, facecolor='green', alpha=0.5)

  plt.title('Distribution of orthologous group gene counts')
  plt.xlabel('$log_2(MHco3(ISE) / McMaster)$')
  plt.ylabel('Occurrences')
  plt.savefig(filename)

def main():
  for input_filename in sys.argv[1:]:
    with open(input_filename) as input_file:
      sm = StateMachine(input_file)
      sm.run()
    results = sm.results()

    summary = summarize_groups(results['groups'])
    plot_groups(results['groups'], 'groups_log.png',    True)
    plot_groups(results['groups'], 'groups_linear.png', False)

    print('%s:' % input_filename)
    for key in sorted(summary.keys()):
      print('  %s=%s' % (key, summary[key]))

if __name__ == '__main__':
  main()
