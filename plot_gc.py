import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
import sys


def do(file):

    d = pd.read_csv(file, sep='\t', header=None)
    print(len(d))
    d = d[d.loc[:, 2]]

    xs, ys = [], []
    for x in range(10):
        print(x, sum(d[0] >= x))
        sb.distplot(d.loc[d[0] >= x, 1], hist=False, label=x, kde=True)
    plt.show()


def main():
    for f in sys.argv[1:]:
        do(f)


if __name__ == '__main__':
    main()
