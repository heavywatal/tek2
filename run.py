#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import wtl.options as wopt

program = 'tek'


def iter(rest):
    const = [program] + rest
    params = wopt.OrderedDict()
    params['lambda'] = ['1e-4', '1e-3']
    params['xi'] = ['1e-5', '1e-4', '1e-3']
    params['nu'] = ['0', '1e-6', '1e-4']
    suffix = '_{}_{}'.format(wopt.now(), wopt.getpid())
    for i, x in enumerate(wopt.product(params)):
        label = wopt.join(x) + suffix + '_{:02}'.format(i)
        yield const + x + ['--outdir=' + label]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-j', '--jobs', type=int, default=wopt.cpu_count())
    (args, rest) = parser.parse_known_args()

    wopt.map_async(iter(rest),
                   args.jobs, args.dry_run)
    print('End of ' + __file__)
