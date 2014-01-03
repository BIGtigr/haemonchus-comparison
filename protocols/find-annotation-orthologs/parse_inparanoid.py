#!/usr/bin/env python2

'''
This script parses the human-readable output file from InParanoid, then dumps
them to a JSON file Unfortunately, only the human-readable output (as well as
the HTML file containing the same information) list all InParanoid results; the
table and CSV file both lack data.

Usage: parse_paranoid.py <InParanoid results file>
'''

import sys
import json

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
    self._current_group = {'a': [], 'b': []}

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

def main():
  inparanoid_results_filename = sys.argv[1]

  with open(inparanoid_results_filename) as inparanoid_results_file:
    sm = StateMachine(inparanoid_results_file)
    sm.run()
  json.dump(sm.results(), sys.stdout)

if __name__ == '__main__':
  main()
