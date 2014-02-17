#!/usr/bin/env python2
import matplotlib
# Force matplotlib not to use X11 backend, which produces exception when run
# over SSH.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot(A, B, filename):
  plt.figure()
  plt.xticks(range(min(A), max(A) + 1))
  plt.yticks(range(min(B), max(B) + 1))
  plt.grid()
  plt.scatter(A, B, s=800, color='#085D82')
  plt.savefig(filename)

def main():
  A = list(range(1,11))
  B = list(range(1,11))
  plot(A, B, 'normal.png')

  B = [1, 2, 3, 8, 7, 6, 5, 4, 9, 10]
  plot(A, B, 'inversion.png')

  B = [(b < 5 and b) or (b + 5) for b in range(1, 11)]
  plot(A, B, 'insertion.png')

  B = [1, 2, 3, 5, 6, 7, 8, 9, 4, 10]
  # I like the appearance of the reversed ordering better.
  plot(B, A, 'translocation.png')

if __name__ == '__main__':
  main()
