# CEGMA's GFF file only lists KOG in attribute column, without providing it in
# key=value format. This script fixes that by prefixing the column with "kog=".
import sys

def main():
  for line in sys.stdin:
    fields = line.strip().split()
    fields[-1] = 'kog=' + fields[-1]
    print('\t'.join(fields))

main()
