import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# for d in pd.read_csv('NA12878_1.small.1.bam.raw', header=None, chunksize=10**8, sep='\t'):


def tel(t, g):
    return (t / g) * 332720800. / 1000. / 46.


def telr(d, l, h):
    return tel(
        sum(d.loc[:,
                  0] >= 7),
        sum((d.loc[:,
                   1] <= (h)) & (d.loc[:,
                                       1] >= (l)))
    )


def telr_new(d, l, h):
    return tel(
        sum((d.loc[:,
                   0] >= 7)
            & (d.loc[:,
                     1] <= (h))
            & (d.loc[:,
                     1] >= (l))),
        sum((d.loc[:,
                   1] <= (h))
            & (d.loc[:,
                     1] >= (l)))
    )


def do(file):

    d = pd.read_csv(file, sep='\t', header=None)

    dd = d[d.loc[:, 2]]

    print(sum(d.loc[:, 0] >= 7))
    print(sum(dd.loc[:, 0] >= 7))

    print(sum((d.loc[:, 1] <= 0.52) & (d.loc[:, 1] >= 0.48)))
    print(sum((dd.loc[:, 1] <= 0.52) & (dd.loc[:, 1] >= 0.48)))

    print(
        sum(
            (d.loc[:,
                   0] >= 7) & (d.loc[:,
                                     1] <= 0.52) & (d.loc[:,
                                                          1] >= 0.48)
        )
    )
    print(
        sum(
            (dd.loc[:,
                    0] >= 7) & (dd.loc[:,
                                       1] <= 0.52) & (dd.loc[:,
                                                             1] >= 0.48)
        )
    )

    g = pd.DataFrame(columns=['d1', 'dd1', 'd2', 'dd2'])
    for r in sorted(
        [0.001,
         0.0025,
         0.005,
         0.025] + list(np.linspace(0,
                                   1,
                                   10)) + list(np.linspace(0,
                                                           0.1,
                                                           10))
    ):
        temp = telr(d, 0.5 - r, 0.5 + r), telr(dd, 0.5 - r, 0.5 + r), telr_new(d, 0.5 - r, 0.5 + r), telr_new(dd, 0.5 - r, 0.5 + r)
        g.loc[r] = temp
        # print(*temp, sep='\t')
    print(g)
    g.plot.line()
    plt.show()


def do2(file):

    d = pd.read_csv(file, sep='\t', header=None)

    dd = d[d.loc[:, 2]]

    g = pd.DataFrame(columns=['d1', 'dd1', 'd2', 'dd2'])
    for r in sorted(
        list(np.linspace(.49,
                         .6,
                         12)) + list(np.linspace(.6,
                                                 1.,
                                                 5))
    ):
        temp = telr(d, 0.48, r), telr(dd, 0.48, r), telr_new(d, 0.48, r), telr_new(dd, 0.48, r)
        g.loc[r] = temp
        # print(*temp, sep='\t')
    print(g)
    g.plot.line()
    plt.show()


def main():
    for f in sys.argv[1:]:
        do(f)
        # do2(f)


if __name__ == '__main__':
    main()
