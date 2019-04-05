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
                  'tel_count'] >= 7),
        sum((d.loc[:,
                   'gc'] <= (h)) & (d.loc[:,
                                          'gc'] >= (l)))
    )


def telr_new(d, l, h):
    return tel(
        sum(
            (d.loc[:,
                   'tel_count'] >= 7)
            & (d.loc[:,
                     'gc'] <= (h))
            & (d.loc[:,
                     'gc'] >= (l))
        ),
        sum((d.loc[:,
                   'gc'] <= (h))
            & (d.loc[:,
                     'gc'] >= (l)))
    )


def yuleq(d, excl):
    def Q(OR):
        (OR - 1) / (OR + 1)

    gb = d.groupby(d.columns).size()
    print(gb)
    ret = pd.DataFrame(columns=d.columns)
    for c1 in d.columns:
        for c2 in d.columns:
            a = gb.xs(True, level=c1).xs(True, level=c2)
            b = gb.xs(True, level=c1).xs(False, level=c2)
            c = gb.xs(False, level=c1).xs(True, level=c2)
            d = gb.xs(False, level=c1).xs(False, level=c2)
            ret.loc[c1, c2] = Q((a * d) / (b * c))
    return ret


def three_binary_test(d):
    print(len(d))

    d['good_length'] = d.length > 7 * 6
    d['is_telomere_old'] = (d.tel_count >= 7)
    d['is_telomere'] = (d.tel_count >= 7) & (d.gc >= 0.48) & (d.gc < 0.52)
    d['is_gc'] = (d.gc >= 0.48) & (d.gc < 0.52)

    # print(d.pivot_table(index='is_telomere', columns='is_dup'))
    # print(d.loc[:, ['is_telomere', 'is_dup']].pivot_table(index='is_telomere', columns='is_dup', aggfunc=len), "\n")
    print(pd.crosstab(d.is_telomere, d.is_dup, margins=True, normalize=True), "\n")
    print(pd.crosstab(d.is_telomere_old, d.is_dup, margins=True, normalize=True), "\n")

    # dd = d.groupby(['is_dup', 'is_primary', 'good_length'])
    # dd = d[d.is_primary].groupby(['is_dup', 'good_length'])
    dd = d[d.is_primary].groupby('is_dup')
    cs = dd.sum().loc[:, ['is_telomere', 'is_telomere_old', 'is_gc']]

    print(cs)
    print()
    print("old\n", tel(cs.is_telomere_old, cs.is_gc))
    print("new\n", tel(cs.is_telomere, cs.is_gc))
    print()
    print(tel(sum(d.is_telomere_old), sum(d.is_gc)))
    print(tel(sum(d.is_telomere), sum(d.is_gc)))
    # print(dd.aggregate(lambda s: tel(sum(s.is_telomere), sum(s.is_gc))))

    # need a corr that accounts for weighting
    # print(d.loc[:, ['is_dup', 'is_primary', 'good_length', 'tel_count']].corr())
    # print(d.loc[:, ['is_dup', 'is_primary', 'good_length', 'is_telomere']].corr())
    # print(yuleq(d.loc[:, ['is_dup', 'is_primary', 'good_length', 'is_telomere']]))


def do(d):
    dd = d[d.loc[:, 'is_dup']]

    print(sum(d.loc[:, 'tel_count'] >= 7))
    print(sum(dd.loc[:, 'tel_count'] >= 7))

    print(sum((d.loc[:, 'gc'] <= 0.52) & (d.loc[:, 'gc'] >= 0.48)))
    print(sum((dd.loc[:, 'gc'] <= 0.52) & (dd.loc[:, 'gc'] >= 0.48)))

    print(
        sum(
            (d.loc[:,
                   'tel_count'] >= 7) & (d.loc[:,
                                               'gc'] <= 0.52) &
            (d.loc[:,
                   'gc'] >= 0.48)
        )
    )
    print(
        sum(
            (dd.loc[:,
                    'tel_count'] >= 7) & (dd.loc[:,
                                                 'gc'] <= 0.52) &
            (dd.loc[:,
                    'gc'] >= 0.48)
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
        d = pd.read_csv(
            f,
            sep='\t',
            header=None,
            names=["tel_count",
                   "gc",
                   "is_dup",
                   "is_primary",
                   "length"]
        )
        # do(d)
        # do2(d)
        three_binary_test(d)


if __name__ == '__main__':
    main()
